#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
This script can be used to parse Caliper data if Caliper was run with the
hatchet-region-profile config. The expected format is json-split.
See http://software.llnl.gov/Caliper/BuiltinConfigurations.html and
https://hatchet.readthedocs.io/en/latest/analysis_examples.html#json-file.

The script requires numpy, pandas and hatchet:
  pip install hatchet
  pip install pandas
  pip install numpy

It can also use [Flamegraph](https://github.com/brendangregg/FlameGraph) to
generate a nice stack-trace-like plot, if Flamegraph is available in your path.
"""

import argparse
import hatchet as ht
import numpy as np
import pandas as pd
import subprocess
import tempfile


def run_flamegraph(txt_file, fg_args):
    """
    Call Flamegraph.
    """
    subprocess.run(["flamegraph.pl"] + fg_args + [txt_file])


def Flamegraph(gf, fg_args):
    with tempfile.NamedTemporaryFile() as folded_stack:
        folded_stack.write(str.encode(gf.to_flamegraph()))
        folded_stack.seek(0)
        run_flamegraph(folded_stack.name, fg_args)


parser = argparse.ArgumentParser(
    description="Utility for parsing/plotting Caliper data with Hatchet."
)
subparsers = parser.add_subparsers(help="sub-command help", dest="command")

fg_parser = subparsers.add_parser("flamegraph")
fg_parser.add_argument("--fg-bin", help="path to Flamegraph binary")
fg_parser.add_argument(
    "json_file",
    help="path to the Caliper json-split data file (generated when the hatchet-region-profile Caliper config is enabled)",
)
fg_parser.add_argument("fg_args", help="Flamegraph arguments", nargs=argparse.REMAINDER)

ht_parser = subparsers.add_parser("hatchet")
ht_parser.add_argument(
    "json_file",
    help="path to the Caliper json-split data file (generated when the hatchet-region-profile Caliper config is enabled)",
)

args = parser.parse_args()
# print(args)

# Use hatchet's ``from_caliper_json`` API with the resulting json-split.
# The result is stored into Hatchet's GraphFrame.
gf = ht.GraphFrame.from_caliper_json(args.json_file)

if args.command == "flamegraph":
    Flamegraph(gf, args.fg_args)
elif args.command == "hatchet":
    print(gf.tree())
