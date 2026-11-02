#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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

import suntools.rk as rk  # noqa: E402

# Field names match the argparse destinations so the two stay in step without the option list being
# written out a third time.
_OPTION_FIELDS = {f.name for f in fields(rk.PlotOptions)}


# ---------------------------------------------------------------------------
# Print Utilities
# ---------------------------------------------------------------------------


def _poly_str(poly: np.poly1d, var: str = "z") -> str:
    return str(np.poly1d(poly.c, variable=var))


def describe(table: rk.ButcherTable) -> None:
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
        print("\nEmbedding stability function phi_hat(z) = p_hat(z) / q(z):")
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


def load_tables(paths: list[str]) -> dict[str, rk.ButcherTable]:
    tables: dict[str, rk.ButcherTable] = {}
    sources: dict[str, str] = {}
    for path in paths:
        parsed = rk.parse_butcher_tables(path)
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


def list_methods(tables: dict[str, rk.ButcherTable], source: str) -> None:
    """Print a one-line summary of every method in the file."""
    print(f"{len(tables)} methods defined in {source}:\n")
    for name, table in tables.items():
        order = f"q={table.method_order}" if table.method_order is not None else "q=?"
        print(f"  {name:<40s} s={table.stages:<3d} {order:<5s} {table.kind()}")


# Column presets for the metrics tables: all explicit tables, all implicit tables, or (mixed) every
# column. Names are rk.StabilityMetrics fields.
_EXPLICIT_COLUMNS = (
    "name",
    "role",
    "order",
    "stage_evaluations",
    "real_extent",
    "imag_extent",
    "real_extent_per_evaluation",
    "imag_extent_per_evaluation",
    "ssp_coefficient",
    "effective_ssp_coefficient",
    "threshold_factor",
)
_IMPLICIT_COLUMNS = (
    "name",
    "role",
    "order",
    "stage_evaluations",
    "a_alpha",
    "r_infinity",
    "a_stable",
    "l_stable",
    "real_extent",
    "imag_extent",
    "ssp_coefficient",
    "effective_ssp_coefficient",
)
_ALL_COLUMNS = tuple(f.name for f in fields(rk.StabilityMetrics))

# Column headings of the printed table.
_TEXT_HEADINGS = {
    "name": "method",
    "role": "role",
    "kind": "kind",
    "order": "q",
    "stages": "stages",
    "stage_evaluations": "evals",
    "real_extent": "real",
    "imag_extent": "imag",
    "real_extent_per_evaluation": "real/eval",
    "imag_extent_per_evaluation": "imag/eval",
    "a_alpha": "A(alpha)",
    "r_infinity": "|R(inf)|",
    "a_stable": "A",
    "l_stable": "L",
    "ssp_coefficient": "C",
    "effective_ssp_coefficient": "C/eval",
    "threshold_factor": "threshold",
}
_LEFT_ALIGNED = ("name", "role", "kind")
_METRICS_EXTENSIONS = {".csv": "csv", ".rst": "rst", ".inc": "rst"}


def _metrics_columns(tables: list[rk.ButcherTable]) -> tuple[str, ...]:
    """The column preset for a selection: explicit-only, implicit-only, or all columns."""
    if all(t.is_explicit for t in tables):
        return _EXPLICIT_COLUMNS
    if not any(t.is_explicit for t in tables):
        return _IMPLICIT_COLUMNS
    return _ALL_COLUMNS


def _text_cell(value) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return "Y" if value else "N"
    if isinstance(value, float):
        return "inf" if np.isinf(value) else f"{value:.4g}"
    return str(value)


def _format_metrics_text(metrics: list[rk.StabilityMetrics], columns: tuple[str, ...]) -> str:
    """Fixed-width table of the selected columns, with a header and a rule."""
    rows = [[_TEXT_HEADINGS[c] for c in columns]]
    rows += [[_text_cell(getattr(m, c)) for c in columns] for m in metrics]
    widths = [max(len(row[i]) for row in rows) for i in range(len(columns))]

    def line(row):
        cells = [
            cell.ljust(width) if name in _LEFT_ALIGNED else cell.rjust(width)
            for cell, width, name in zip(row, widths, columns, strict=True)
        ]
        return "  ".join(cells).rstrip()

    return "\n".join(
        [line(rows[0]), "  ".join("-" * w for w in widths)] + [line(r) for r in rows[1:]]
    )


def _select_metrics(
    tables: dict[str, rk.ButcherTable], methods: list[str] | None, include_embedding: bool
) -> tuple[list[rk.StabilityMetrics], list[rk.ButcherTable]]:
    """Metrics rows in file order (method, then its embedding) and the tables they came from."""
    selected = [tables[name] for name in (methods if methods is not None else tables)]
    rows = []
    for table in selected:
        rows.append(rk.stability_metrics(table))
        if include_embedding and table.has_embedding:
            rows.append(rk.stability_metrics(table, embedding=True))
    return rows, selected


def write_metrics(
    path: str, rows: list[rk.StabilityMetrics], selected: list[rk.ButcherTable]
) -> None:
    """Write the metrics to *path* (``-`` is CSV on stdout); the extension picks the format."""
    if path == "-":
        sys.stdout.write(rk.format_metrics_csv(rows))
        return
    if _METRICS_EXTENSIONS[os.path.splitext(path)[1].lower()] == "csv":
        text = rk.format_metrics_csv(rows)
    else:
        text = (
            ".. Generated by scripts/plot_rk_stability_regions.py --metrics; do not edit.\n\n"
            + rk.format_metrics_rst(rows, columns=_metrics_columns(selected))
        )
    with open(path, "w", newline="") as f:
        f.write(text)
    print(f"Wrote metrics for {len(rows)} entries to: {path}")


def _star_path(base: str, suffix: str) -> str:
    """Insert *suffix* before the extension of *base*."""
    stem, ext = os.path.splitext(base)
    return f"{stem}{suffix}{ext}"


def _plot_stars(
    table: rk.ButcherTable,
    options: rk.PlotOptions,
    method_file: str,
    embedding_file: str,
    order_star: bool,
    order_star_embedding: bool,
    indent: str = "",
) -> None:
    """Write the requested order-star plots for one method."""
    if order_star:
        rk.plot_order_star(table, options, filename=method_file, close=True)
        print(f"{indent}Saved order star plot to: {method_file}")
    if order_star_embedding:
        if table.has_embedding:
            rk.plot_order_star(table, options, filename=embedding_file, close=True, embedding=True)
            print(f"{indent}Saved embedding order star plot to: {embedding_file}")
        else:
            print(f"{indent}Skipping embedding order star for {table.name}: no embedding")


def plot_one(
    tables: dict[str, rk.ButcherTable],
    method: str,
    outfile: str | None,
    options: rk.PlotOptions,
    order_star: bool = False,
    order_star_embedding: bool = False,
) -> None:
    """Describe and plot a single named method."""
    if method not in tables:
        raise SystemExit(f"Method {method!r} not found. Use --list to see available methods.")
    table = tables[method]
    describe(table)
    base = outfile or f"{method}_stab_region.png"
    rk.plot_stability_region(table, options, filename=base, close=True)
    print(f"Saved stability region plot to: {base}")
    if outfile is None:
        star_file, embedding_file = (
            f"{method}_order_star.png",
            f"{method}_order_star_embedding.png",
        )
    else:
        star_file = _star_path(outfile, "_order_star")
        embedding_file = _star_path(outfile, "_order_star_embedding")
    _plot_stars(table, options, star_file, embedding_file, order_star, order_star_embedding)


def plot_all(
    tables: dict[str, rk.ButcherTable],
    source: str,
    outdir: str,
    options: rk.PlotOptions,
    order_star: bool = False,
    order_star_embedding: bool = False,
) -> None:
    """Write one plot per method in the file into *outdir*."""
    os.makedirs(outdir, exist_ok=True)
    print(f"Plotting all {len(tables)} methods from {source} into {outdir}/\n")
    for i, (name, table) in enumerate(tables.items(), 1):
        outfile = os.path.join(outdir, f"{name}_stab_region.png")
        rk.plot_stability_region(table, options, filename=outfile, close=True)
        print(f"  [{i:2d}/{len(tables)}] {name:<40s} -> {outfile}")
        _plot_stars(
            table,
            options,
            os.path.join(outdir, f"{name}_order_star.png"),
            os.path.join(outdir, f"{name}_order_star_embedding.png"),
            order_star,
            order_star_embedding,
            indent="           ",
        )
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
    ap.add_argument(
        "--metrics",
        metavar="PATH",
        default=None,
        help="write the stability metrics of the selected methods (all, or those given by "
        "--method) and exit; the extension picks the format: .csv, or .rst/.inc for an RST "
        "table; '-' writes CSV to stdout. Embedding rows are included unless --no-embedding",
    )
    ap.add_argument(
        "--print-metrics",
        action="store_true",
        help="print a table of the stability metrics of the selected methods and exit",
    )
    ap.add_argument(
        "--order-star",
        action="store_true",
        help="also write the order star of each method next to its stability region plot "
        "(<name>_order_star.png)",
    )
    ap.add_argument(
        "--order-star-embedding",
        action="store_true",
        help="also write the order star of each embedding (<name>_order_star_embedding.png); "
        "methods without an embedding are skipped",
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
        default=rk.PlotOptions.grid_points,
        help=f"grid resolution per axis (default: {rk.PlotOptions.grid_points})",
    )
    ap.add_argument(
        "--figsize",
        dest="figure_size",
        type=float,
        nargs=2,
        default=rk.PlotOptions.figure_size,
        metavar=("W", "H"),
        help="output canvas in inches; every plot is written at this size "
        f"(default: {rk.PlotOptions.figure_size[0]:g} {rk.PlotOptions.figure_size[1]:g})",
    )
    ap.add_argument(
        "--dpi",
        type=int,
        default=rk.PlotOptions.dpi,
        help=f"output resolution; pixel size is figsize x dpi (default: {rk.PlotOptions.dpi})",
    )
    ap.add_argument(
        "--font-size",
        dest="font_size",
        type=float,
        default=rk.PlotOptions.font_size,
        help=f"font size in points (default: {rk.PlotOptions.font_size})",
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
        "--axis-extent",
        dest="show_axis_extent",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="label the real/imaginary stability interval endpoints",
    )
    ap.add_argument(
        "--suppress-islands",
        dest="suppress_tiny_islands",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="hide tiny disconnected stable islands around zeros of phi (render as stray dots)",
    )
    ap.add_argument(
        "--highlight-imag-axis-intervals",
        action="store_true",
        help="highlight finite stable intervals on the imaginary axis; "
        "ignored with --no-axis-extent",
    )
    ap.add_argument(
        "--label-detached-imag-axis-intervals",
        action="store_true",
        help="also label detached finite stable intervals on the imaginary axis; "
        "ignored with --no-axis-extent",
    )
    ap.add_argument(
        "--annotate-metrics",
        dest="annotate_metrics",
        action="store_true",
        help="draw the A(alpha) rays and a box of stability metrics on each region plot",
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

    metrics_mode = args.metrics is not None or args.print_metrics
    if metrics_mode:
        if args.out is not None or args.order_star or args.order_star_embedding:
            ap.error("--metrics and --print-metrics do not plot; drop --out and --order-star*")
        if args.metrics not in (None, "-") and (
            os.path.splitext(args.metrics)[1].lower() not in _METRICS_EXTENSIONS
        ):
            ap.error(f"--metrics: unknown extension in {args.metrics!r}; use .csv, .rst or .inc")

    tables = load_tables(args.deffile)
    source = ", ".join(os.path.basename(path) for path in args.deffile)
    if args.list:
        list_methods(tables, source)
        return

    if metrics_mode:
        if args.method is not None:
            for method in args.method:
                if method not in tables:
                    ap.error(f"Method {method!r} not found. Use --list to see available methods.")
        rows, selected = _select_metrics(tables, args.method, args.show_embedding)
        if args.metrics is not None:
            write_metrics(args.metrics, rows, selected)
        if args.print_metrics:
            print(_format_metrics_text(rows, _metrics_columns(selected)))
        return

    # PlotOptions field names mirror the argparse destinations, so the option list is
    # written out once and flows through unchanged.
    options = rk.PlotOptions(**{k: v for k, v in vars(args).items() if k in _OPTION_FIELDS})

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
            plot_one(tables, method, outfile, options, args.order_star, args.order_star_embedding)
    else:
        plot_all(tables, source, args.outdir, options, args.order_star, args.order_star_embedding)


if __name__ == "__main__":
    main()
