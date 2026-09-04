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
from typing import List, Optional

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


def _cmd_tune(args: argparse.Namespace) -> int:
    try:
        from suntools.tune.cli import run_from_args
    except ModuleNotFoundError as err:
        if err.name in ("pydantic", "yaml", "deephyper", "gptune", "ytopt"):
            sys.stderr.write(
                "error: suntools tune requires the suntools project dependencies "
                "(pydantic, PyYAML, DeepHyper, GPTune, and Ytopt)\n"
            )
            return 2
        raise

    return run_from_args(args)


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

    tune = subparsers.add_parser(
        "tune",
        help="Tune a SUNDIALS executable by appending SetOptions KEY VALUE pairs.",
    )
    tune.add_argument(
        "--config",
        help="YAML tune configuration. When set, command-line tune fields are ignored.",
    )
    tune.add_argument(
        "--params",
        action="append",
        nargs=2,
        metavar=("KEY", "SPEC"),
        help=(
            "Tunable SetOptions parameter. SPEC is LOW:HIGH, LOW:HIGH:log, "
            "int:LOW:HIGH, or choice:v1,v2,v3. Choice values may contain "
            "whitespace-separated SetOptions values. May be repeated."
        ),
    )
    tune.add_argument(
        "--max-evals",
        type=int,
        default=40,
        help="Maximum number of objective evaluations.",
    )
    tune.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of worker threads for backend evaluations.",
    )
    tune.add_argument(
        "--output-dir",
        default="suntools-tune",
        help="Directory for tune results.",
    )
    tune.add_argument(
        "--backend",
        default="deephyper",
        help=(
            'Optimization backend name: "deephyper", "gptune", or "ytopt" '
            '(default: "deephyper").'
        ),
    )
    tune.add_argument(
        "--backend-option",
        action="append",
        help="Backend option in KEY=VALUE form. May be repeated.",
    )
    tune.add_argument(
        "--cwd",
        default=".",
        help="Working directory for the executable.",
    )
    tune.add_argument(
        "--env",
        action="append",
        help="Environment override in KEY=VALUE form. May be repeated.",
    )
    tune.add_argument(
        "--metric",
        default="wall_time",
        help='Objective metric name (default: "wall_time").',
    )
    tune.add_argument(
        "--direction",
        choices=("minimize", "maximize"),
        default="minimize",
        help="Objective direction.",
    )
    tune.add_argument(
        "--objective-source",
        default=None,
        help='Source for regex objectives: "stdout", "stderr", or a file path.',
    )
    tune.add_argument(
        "--objective-regex",
        default=None,
        help="Regular expression used to extract a numeric objective.",
    )
    tune.add_argument(
        "--objective-group",
        default=1,
        help="Regex group index or name used for the objective value.",
    )
    tune.add_argument(
        "executable",
        nargs=argparse.REMAINDER,
        help="Executable command followed by its arguments.",
    )
    tune.set_defaults(func=_cmd_tune)

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
