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

"""Runge-Kutta Butcher table .def file parser."""

from __future__ import annotations

import math
import re

import numpy as np

from .rk_butcher_table import ButcherTable

# ---------------------------------------------------------------------------
# Parsing the .def file
# ---------------------------------------------------------------------------


def _strip_c_comments(text: str) -> str:
    """Remove C comments so they cannot confuse the brace/statement parsing downstream."""
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.DOTALL)  # DOTALL: block comments span lines
    text = re.sub(r"//[^\n]*", " ", text)
    return text


def _split_top_level_methods(text: str):
    """Yield (name, body) for each ARK_BUTCHER_TABLE(NAME, { ... }) block."""
    # *text* must already be comment-stripped; the depth counter below cannot tell a brace in code
    # from one in comments.
    for m in re.finditer(r"ARK_BUTCHER_TABLE\s*\(\s*([A-Za-z0-9_]+)\s*,", text):
        name = m.group(1)
        i = text.find("{", m.end())  # first brace after the macro header
        if i == -1:
            continue
        # A regex cannot isolate the body because the C initializer contains nested braces.
        # Instead track nesting depth from the opening brace; the "}" that brings depth back to
        # zero is the one that closes the block.
        depth, j = 0, i
        while j < len(text):
            if text[j] == "{":
                depth += 1
            elif text[j] == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        yield name, text[i + 1 : j]  # the text strictly between the outermost braces


# Map SUNDIALS functions/macros to Python.
_SUN_FUNCS = {
    "SUN_RCONST": (lambda x: x),
    "SUNRsqrt": math.sqrt,
    "SUNRabs": abs,
    "SUNRexp": math.exp,
    "SUNRpowerR": math.pow,
    "SUNRsin": math.sin,
    "SUNRcos": math.cos,
}

# Matches a table entry such as `B->b[0]` or `B->A[4][0]`. Group 1 is the field letter (A/b/c/d);
# group 2 is the whole index string, "[i]" or "[i][j]". The same pattern works for both sides of an
# assignment: anchored, it identifies an assignment target; unanchored, it finds references on a
# right-hand side.
_ENTRY_RE = re.compile(r"B->([Abcd])\s*((?:\[\d+\])+)")
_TARGET_RE = re.compile(_ENTRY_RE.pattern + r"\s*$")
_INDEX_RE = re.compile(r"\[(\d+)\]")


def _indices(index_string: str) -> tuple[int, ...]:
    """Turn an index string "[i]" or "[i][j]" into a tuple of ints."""
    return tuple(int(t) for t in _INDEX_RE.findall(index_string))


def _resolve_refs(expr: str, arrays: dict) -> str:
    """Splice the current value of each `B->A[i][j]`/`B->b[i]`/... reference into *expr*."""

    # Resolve entries defined by reference (e.g. `B->A[s-1][j] = B->b[j]`). Substituting a float
    # literal for each reference makes the expression eval-able. Order-dependent, a reference reads
    # whatever the entry holds at that point in the body, 0 if it is assigned later in the file.
    def repl(m):
        return repr(float(arrays[m.group(1)][_indices(m.group(2))]))

    return _ENTRY_RE.sub(repl, expr)


def _evaluate(expr: str, env: dict, arrays: dict | None = None) -> float:
    """Evaluate a numeric expression from the .def file. Covers arithmetic operations, fractions,
    scientific notation, SUN* helpers, and named constants.
    """
    # Collapse newlines and indentation: an expression split across lines would otherwise reach
    # eval as an indented continuation line and raise IndentationError.
    expr = " ".join(expr.split())
    if arrays is not None and "B->" in expr:
        expr = _resolve_refs(expr, arrays)  # table refs -> float literals, so eval can see them

    # __builtins__ is emptied so no Python built-ins (open, __import__, ...) are reachable, leaving
    # only the whitelisted SUN* helpers and the file's named constants in scope.  Adequate for the
    # trusted SUNDIALS .def files this targets; NOT a general-purpose sandbox.
    namespace = {"__builtins__": {}}
    namespace.update(_SUN_FUNCS)
    namespace.update(env)
    return float(eval(expr, namespace))  # noqa: S307 - sandboxed namespace


def _parse_body(name: str, body: str) -> ButcherTable | None:
    """Turn one macro body (C source filling a table struct) into a ButcherTable."""
    # *body* arrives already comment-stripped. The Alloc call fixes the stage count and the
    # embedding flag. Its absence means a stub body (e.g. `return NULL`) with no table to build.
    alloc = re.search(r"ARKodeButcherTable_Alloc\s*\(\s*(\d+)\s*,\s*(SUNTRUE|SUNFALSE)\s*\)", body)
    if not alloc:
        return None
    stages = int(alloc.group(1))

    # Pre-fill with zeros so entries (unassigned values default to 0); `d` holds the embedding.
    arrays = {
        "A": np.zeros((stages, stages)),
        "b": np.zeros(stages),
        "c": np.zeros(stages),
        "d": np.zeros(stages),
    }
    q = p = None  # method / embedding orders
    has_embedding = alloc.group(2) == "SUNTRUE"
    env: dict = {}  # local constants defined in the body

    # Walk the body one statement at a time, in file order, so that references to earlier entries
    # resolve the way they would during C execution.
    for raw in body.split(";"):
        stmt = raw.strip()
        # Skip blanks, the trailing `return B`, and the already-consumed Alloc call.
        if not stmt or stmt.startswith("return") or "ARKodeButcherTable_Alloc" in stmt:
            continue

        # Remember local scalar constants, e.g. `const sunrealtype gamma = ...;` in env so later
        # expressions can reference them by name. Use DOTALL so a definition split across lines is
        # captured as a whole. Without it `(.+)$` cannot cross the newline, the match fails, the
        # statement falls through to the array-assignment path, and the constant is silently never
        # defined (surfaces much later as a NameError from inside eval).
        cdef = re.match(
            r"(?:const\s+)?sunrealtype\s+([A-Za-z_]\w*)\s*=\s*(.+)$", stmt, flags=re.DOTALL
        )
        if cdef:
            env[cdef.group(1)] = _evaluate(cdef.group(2), env, arrays)
            continue

        if "=" not in stmt:
            continue

        # Chained assignment `B->A[2][0] = B->A[2][1] = 0.5` splits into several targets sharing
        # one right-hand side (the last '='-separated piece).
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

        # Otherwise an array assignment: evaluate the RHS once, then route the value to each
        # indexed target using the field letter and indices the regex captured.
        value = _evaluate(rhs, env, arrays)
        for t in targets:
            m = _TARGET_RE.match(t)
            if not m:
                continue  # target shape the regex does not cover, silently dropped
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
    """Parse every Butcher table in *path*; returns name -> ButcherTable, in file order."""
    with open(path, "r") as fh:
        text = _strip_c_comments(fh.read())  # strip first; _split_top_level_methods requires it
    tables: dict[str, ButcherTable] = {}
    for name, body in _split_top_level_methods(text):
        table = _parse_body(name, body)
        if table is not None:  # None == stub body with no Alloc call
            tables[name] = table
    return tables
