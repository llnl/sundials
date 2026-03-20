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

import argparse
import sys

from suntools import logs as sunlogs


def _cmd_parse_logs(args: argparse.Namespace) -> int:
    selected = sunlogs._split_filter_list(args.filter)
    valid = {"integrator", "nonlinear", "linear"}
    unknown = sorted(selected - valid)
    if unknown:
        sys.stderr.write(f"error: unknown filter categories: {', '.join(unknown)}\n")
        return 2

    if args.input == "-":
        if sys.stdin.isatty():
            sys.stderr.write("error: no stdin provided\n")
            return 2
        lines = sys.stdin.readlines()
    else:
        with open(args.input, "r") as fp:
            lines = fp.readlines()

    for out_line in sunlogs._iter_filtered_lines(lines, selected, invert=args.invert):
        sys.stdout.write(out_line)

    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="suntools")
    subparsers = parser.add_subparsers(dest="command", required=True)

    parse_logs = subparsers.add_parser(
        "parse_logs", help="Filter SUNLogger logs by category (reads stdin by default)."
    )
    parse_logs.add_argument(
        "-f",
        "--filter",
        default="integrator,nonlinear,linear",
        help='Comma-separated list: "integrator,nonlinear,linear".',
    )
    parse_logs.add_argument(
        "-i",
        "--invert",
        action="store_true",
        help="Invert the filter (exclude selected categories).",
    )
    parse_logs.add_argument(
        "input", nargs="?", default="-", help='Input logfile path (default: "-" for stdin).'
    )
    parse_logs.set_defaults(func=_cmd_parse_logs)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
