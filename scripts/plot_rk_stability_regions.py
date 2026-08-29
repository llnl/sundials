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

"""Plot absolute-stability regions from a SUNDIALS .def file."""

from __future__ import annotations

import argparse
import os
from dataclasses import fields

import matplotlib

matplotlib.use("Agg")  # batch rendering: no display needed, and much faster

import numpy as np  # noqa: E402

import sys

repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
suntools_src = os.path.join(repo_root, "suntools", "src")
if suntools_src not in sys.path:
    sys.path.insert(0, suntools_src)

from suntools.rk_butcher_table import ButcherTable  # noqa: E402
from suntools.rk_stability_plotting import PlotOptions, plot_stability_region  # noqa: E402
from suntools.rk_table_parser import parse_butcher_tables  # noqa: E402

# Field names match the argparse destinations so the two stay in step without the option list being
# written out a third time.
_OPTION_FIELDS = {f.name for f in fields(PlotOptions)}


# ---------------------------------------------------------------------------
# Print Utilities
# ---------------------------------------------------------------------------


def _poly_str(poly: np.poly1d, var: str = "z") -> str:
    return str(np.poly1d(poly.c, variable=var))


def describe(table: ButcherTable) -> None:
    """Print the tableau and stability function(s) for one method."""
    print(table)
    print()
    phi = table.stability_function
    p, q = phi.p, phi.q
    print("Method stability function phi(z) = p(z) / q(z):")
    print()
    print("p(z):")
    print(_poly_str(p))
    print()
    print("q(z):")
    print(_poly_str(q))

    phi = table.embedded_stability_function
    if phi is not None:
        print(f"\nEmbedding stability function phi_hat(z) = p_hat(z) / q(z):")
        print()
        print("p_hat(z):")
        print(_poly_str(phi.p))
        print()
        print("q(z):")
        print(_poly_str(q))
    print()


# ---------------------------------------------------------------------------
# Load, List, and Plotting Functions
# ---------------------------------------------------------------------------


def load_tables(paths: list[str]) -> dict[str, ButcherTable]:
    tables: dict[str, ButcherTable] = {}
    sources: dict[str, str] = {}
    for path in paths:
        parsed = parse_butcher_tables(path)
        if not parsed:
            raise SystemExit(f"No Butcher tables found in {path!r}")
        for name, table in parsed.items():
            if name in tables:
                raise SystemExit(
                    f"Duplicate method name {name!r} found in {path!r} "
                    f"(already loaded from {sources[name]!r})"
                )
            tables[name] = table
            sources[name] = path
    return tables


def list_methods(tables: dict[str, ButcherTable], source: str) -> None:
    """Print a one-line summary of every method in the file."""
    print(f"{len(tables)} methods defined in {source}:\n")
    for name, table in tables.items():
        order = f"q={table.method_order}" if table.method_order is not None else "q=?"
        print(f"  {name:<40s} s={table.stages:<3d} {order:<5s} {table.kind()}")


def plot_one(
    tables: dict[str, ButcherTable], method: str, outfile: str | None, options: PlotOptions
) -> None:
    """Describe and plot a single named method."""
    if method not in tables:
        raise SystemExit(f"Method {method!r} not found. Use --list to see available methods.")
    table = tables[method]
    describe(table)
    outfile = outfile or f"{method}_stab_region.png"
    plot_stability_region(table, options, filename=outfile, close=True)
    print(f"Saved stability region plot to: {outfile}")


def plot_all(
    tables: dict[str, ButcherTable], source: str, outdir: str, options: PlotOptions
) -> None:
    """Write one plot per method in the file into *outdir*."""
    os.makedirs(outdir, exist_ok=True)
    print(f"Plotting all {len(tables)} methods from {source} into {outdir}/\n")
    for i, (name, table) in enumerate(tables.items(), 1):
        outfile = os.path.join(outdir, f"{name}_stab_region.png")
        plot_stability_region(table, options, filename=outfile, close=True)
        print(f"  [{i:2d}/{len(tables)}] {name:<40s} -> {outfile}")
    print(f"\nDone. {len(tables)} plots written to {outdir}/")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Parse a SUNDIALS Runge-Kutta Butcher table file "
        "and plot absolute-stability regions."
    )
    ap.add_argument("deffile", nargs="+", help="path(s) to the .def file(s)")
    ap.add_argument(
        "-m",
        "--method",
        nargs="+",
        default=None,
        help="one or more exact method names; if omitted, ALL methods in the file are plotted",
    )
    ap.add_argument(
        "-o",
        "--out",
        nargs="+",
        default=None,
        help="one or more output image paths (requires --method)",
    )
    ap.add_argument(
        "--outdir",
        default="stability_plots",
        help="output directory for all-methods mode (default: stability_plots)",
    )
    ap.add_argument(
        "-l", "--list", action="store_true", help="list all methods in the file and exit"
    )

    # Everything below sets a PlotOptions field; the `dest` names must match its fields.
    ap.add_argument(
        "--bounds",
        type=float,
        nargs=4,
        default=None,
        metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
        help="manual plot bounds (overrides auto-framing); grown, never shrunk, to the "
        "aspect of the fixed canvas so the axes stay equal-aspect without letterboxing",
    )
    ap.add_argument(
        "--grid",
        dest="grid_points",
        type=int,
        default=PlotOptions.grid_points,
        help=f"grid resolution per axis (default: {PlotOptions.grid_points})",
    )
    ap.add_argument(
        "--figsize",
        dest="figure_size",
        type=float,
        nargs=2,
        default=PlotOptions.figure_size,
        metavar=("W", "H"),
        help="output canvas in inches; every plot is written at this size "
        f"(default: {PlotOptions.figure_size[0]:g} {PlotOptions.figure_size[1]:g})",
    )
    ap.add_argument(
        "--dpi",
        type=int,
        default=PlotOptions.dpi,
        help=f"output resolution; pixel size is figsize x dpi (default: {PlotOptions.dpi})",
    )
    ap.add_argument(
        "--font-size",
        dest="font_size",
        type=float,
        default=PlotOptions.font_size,
        help=f"font size in points (default: {PlotOptions.font_size})",
    )
    ap.add_argument(
        "--zeros",
        dest="show_zeros",
        action="store_true",
        help="mark the zeros of phi (numerator roots)",
    )
    ap.add_argument(
        "--poles",
        dest="show_poles",
        action="store_true",
        help="mark the poles of phi (denominator roots; implicit methods only)",
    )
    ap.add_argument(
        "--embedding",
        dest="show_embedding",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="overlay the embedded method's stability region",
    )
    ap.add_argument(
        "--suppress-islands",
        dest="suppress_tiny_islands",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="hide tiny disconnected stable islands around zeros of phi (render as stray dots)",
    )
    shade_group = ap.add_mutually_exclusive_group()
    shade_group.add_argument(
        "--shade", action="store_true", help="shade the method stable region by |phi(z)|"
    )
    shade_group.add_argument(
        "--shade-embedding", action="store_true", help="shade the embedding stable region |phi(z)|"
    )
    return ap


def main(argv=None):
    ap = build_parser()
    args = ap.parse_args(argv)
    if args.out is not None and args.method is None:
        ap.error("--out applies to methods; pass --method as well")

    tables = load_tables(args.deffile)
    source = ", ".join(os.path.basename(path) for path in args.deffile)
    if args.list:
        list_methods(tables, source)
        return

    # PlotOptions field names mirror the argparse destinations, so the option list is
    # written out once and flows through unchanged.
    options = PlotOptions(**{k: v for k, v in vars(args).items() if k in _OPTION_FIELDS})

    if args.method is not None:
        methods = args.method
        outs = args.out
        for method in methods:
            if method not in tables:
                ap.error(f"Method {method!r} not found. Use --list to see available methods.")
        if outs is not None and len(outs) != len(methods):
            ap.error("--out requires exactly one path per --method value")

        for idx, method in enumerate(methods):
            outfile = outs[idx] if outs is not None else None
            plot_one(tables, method, outfile, options)
    else:
        plot_all(tables, source, args.outdir, options)


if __name__ == "__main__":
    main()
