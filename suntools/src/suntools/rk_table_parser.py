#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------

"""
rk_table_parser.py
==================

Parse Runge-Kutta Butcher tables from SUNDIALS / ARKODE X-macro definition
files (``arkode_butcher_erk.def`` or ``arkode_butcher_dirk.def``) into
:class:`~suntools.rk_butcher_table.ButcherTable` objects.

Each method is an ``ARK_BUTCHER_TABLE(NAME, { ... })`` block whose body is C
source filling a table struct.  Accommodated quirks:

  * ``SUN_RCONST(x)`` literal wrappers (treated as the identity).
  * Fractions such as ``SUN_RCONST(1.0)/SUN_RCONST(3.0)`` and ``+ - * /`` exprs.
  * Local named constants, ``[const] sunrealtype gamma = ...;``.
  * Chained assignments, ``B->A[2][0] = B->A[2][1] = SUN_RCONST(0.5);``.
  * Scientific notation and unary minus.
  * SUNDIALS math helpers, e.g. ``SUNRsqrt(...)``.
  * RHS references to earlier entries, ``B->A[4][0] = B->b[0];``.
  * ``B->d = NULL;`` meaning "no embedding".

Unset entries default to zero.  Methods that just ``return NULL`` are skipped.
"""

from __future__ import annotations

import math
import re

import numpy as np

from .rk_butcher_table import ButcherTable

# ---------------------------------------------------------------------------
# Parsing the .def file
# ---------------------------------------------------------------------------


def _strip_c_comments(text: str) -> str:
    """Remove C comments so they cannot confuse the statement/regex parsing downstream.

    Replaces /* ... */ blocks (DOTALL so they may span lines) and // ... line comments
    with a single space each, preserving surrounding token boundaries.
    """
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.DOTALL)
    text = re.sub(r"//[^\n]*", " ", text)
    return text


def _split_top_level_methods(text: str):
    """Yield (name, body) for each ARK_BUTCHER_TABLE(NAME, { ... }) block.

    A regex alone cannot isolate the body because the C initializer contains nested
    braces; instead we locate each macro header, then brace-match from the first "{" to
    its partner by tracking nesting depth.  ``body`` is the text strictly between the
    outermost braces.
    """
    for m in re.finditer(r"ARK_BUTCHER_TABLE\s*\(\s*([A-Za-z0-9_]+)\s*,", text):
        name = m.group(1)
        i = text.find("{", m.end())  # first brace after the macro header
        if i == -1:
            continue
        # Walk forward, +1 on "{" and -1 on "}", until depth returns to zero: that "}"
        # closes the block the opening brace started.
        depth, j = 0, i
        while j < len(text):
            if text[j] == "{":
                depth += 1
            elif text[j] == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        yield name, text[i + 1 : j]


# SUNDIALS scalar helpers seen in the .def files, mapped to Python callables.
_SUN_FUNCS = {
    "SUN_RCONST": (lambda x: x),
    "SUNRsqrt": math.sqrt,
    "SUNRabs": abs,
    "SUNRexp": math.exp,
    "SUNRpowerR": math.pow,
    "SUNRsin": math.sin,
    "SUNRcos": math.cos,
}

# Matches a tableau entry such as `B->b[0]` or `B->A[4][0]`.  Group 1 is the field letter
# (A/b/c/d); group 2 is the whole "[i]" or "[i][j]" index string.  The same pattern serves
# both sides of an assignment: anchored, it identifies an assignment target; unanchored, it
# finds references on a right-hand side.
_ENTRY_RE = re.compile(r"B->([Abcd])\s*((?:\[\d+\])+)")
_TARGET_RE = re.compile(_ENTRY_RE.pattern + r"\s*$")
_INDEX_RE = re.compile(r"\[(\d+)\]")


def _indices(index_string: str) -> tuple[int, ...]:
    """Turn an "[i]" / "[i][j]" index string into a tuple of ints."""
    return tuple(int(t) for t in _INDEX_RE.findall(index_string))


def _resolve_refs(expr: str, arrays: dict) -> str:
    """Splice the current value of each `B->A[i][j]`/`B->b[i]`/... reference into *expr*.

    Stiffly-accurate tableaux define entries by reference (e.g. `B->A[s-1][j] = B->b[j]`).
    Each reference is looked up in the arrays filled so far and replaced by a Python float
    literal so the expression can be eval'd.  This is order-dependent: a reference to an
    entry assigned *later* in the file resolves to that entry's current (default 0) value.
    """

    def repl(m):
        return repr(float(arrays[m.group(1)][_indices(m.group(2))]))

    return _ENTRY_RE.sub(repl, expr)


def _evaluate(expr: str, env: dict, arrays: dict | None = None) -> float:
    """Evaluate a C-style numeric expression from the .def file to a float.

    Handles arithmetic (+ - * /), fraction expressions, scientific notation, the SUN*
    math helpers, and locally-defined named constants (passed in via ``env``).  Table
    references are resolved to literals first (see _resolve_refs).

    Uses eval on a locked-down namespace: __builtins__ is emptied so no Python built-ins
    (open, __import__, etc.) are reachable, and only the whitelisted SUN* helpers plus the
    file's own named constants are exposed.  This is safe for the trusted,
    machine-generated SUNDIALS .def files it targets but is not a general-purpose sandbox.
    """
    expr = expr.strip()
    if arrays is not None and "B->" in expr:
        expr = _resolve_refs(expr, arrays)
    namespace = {"__builtins__": {}}
    namespace.update(_SUN_FUNCS)
    namespace.update(env)
    return float(eval(expr, namespace))  # noqa: S307 - sandboxed namespace


def _parse_body(name: str, body: str) -> ButcherTable | None:
    """Turn one macro body (C source filling a table struct) into a ButcherTable.

    Strategy: find the allocation call (which gives the stage count and whether an
    embedding was requested), pre-fill A/b/c/d with zeros so unset entries default to 0,
    then walk the body statement by statement (split on ';') updating those arrays.
    Returns None for stub bodies that never allocate a table (e.g. `return NULL`).
    """
    body = _strip_c_comments(body)
    # The Alloc call fixes the stage count and the SUNTRUE/SUNFALSE embedding flag.
    alloc = re.search(r"ARKodeButcherTable_Alloc\s*\(\s*(\d+)\s*,\s*(SUNTRUE|SUNFALSE)\s*\)", body)
    if not alloc:
        return None
    stages = int(alloc.group(1))

    # Unset entries stay zero; `d` holds the embedding weights (SUNDIALS field B->d).
    arrays = {
        "A": np.zeros((stages, stages)),
        "b": np.zeros(stages),
        "c": np.zeros(stages),
        "d": np.zeros(stages),
    }
    q = p = None  # method / embedding orders (B->q, B->p)
    has_embedding = alloc.group(2) == "SUNTRUE"
    env: dict = {}  # local `sunrealtype` constants defined in the body

    for raw in body.split(";"):
        stmt = raw.strip()
        # Skip blanks, the trailing `return B`, and the already-consumed Alloc call.
        if not stmt or stmt.startswith("return") or "ARKodeButcherTable_Alloc" in stmt:
            continue

        # Local scalar constant, e.g. `const sunrealtype gamma = ...;` -> remember it in
        # env so later expressions can reference it by name.
        cdef = re.match(r"(?:const\s+)?sunrealtype\s+([A-Za-z_]\w*)\s*=\s*(.+)$", stmt)
        if cdef:
            env[cdef.group(1)] = _evaluate(cdef.group(2), env, arrays)
            continue

        if "=" not in stmt:
            continue

        # Chained assignment `B->A[2][0] = B->A[2][1] = 0.5` splits into several targets
        # sharing one right-hand side (the last '='-separated piece).
        parts = [t.strip() for t in stmt.split("=")]
        rhs, targets = parts[-1], parts[:-1]

        # `B->d = NULL;` explicitly disables the embedding.
        if any(t == "B->d" for t in targets) and rhs.upper() == "NULL":
            has_embedding = False
            continue

        # Order fields are integers, not tableau entries.
        if all(t in ("B->q", "B->p") for t in targets):
            val = int(round(_evaluate(rhs, env, arrays)))
            for t in targets:
                if t == "B->q":
                    q = val
                else:
                    p = val
            continue

        # Otherwise an array assignment: evaluate the RHS once, then route the value to
        # each indexed target using the field letter and indices the regex captured.
        value = _evaluate(rhs, env, arrays)
        for t in targets:
            m = _TARGET_RE.match(t)
            if not m:
                continue
            field = m.group(1)
            arrays[field][_indices(m.group(2))] = value
            if field == "d":
                has_embedding = True  # writing any d[i] implies an embedding is present

    return ButcherTable(
        name,
        stages,
        arrays["A"],
        arrays["b"],
        arrays["c"],
        b_embedded=arrays["d"] if has_embedding else None,
        method_order=q,
        embedding_order=p,
    )


def parse_butcher_tables(path: str) -> dict[str, ButcherTable]:
    """Parse every Butcher table defined in *path*.

    Returns a dict mapping method name -> :class:`ButcherTable`, in file order.
    """
    with open(path, "r") as fh:
        text = fh.read()
    tables: dict[str, ButcherTable] = {}
    for name, body in _split_top_level_methods(text):
        table = _parse_body(name, body)
        if table is not None:
            tables[name] = table
    return tables
