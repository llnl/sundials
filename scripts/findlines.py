#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------
# Find unique lines from a file containing a specific string and write them to
# a separate file.
#
# Example usage:
#   Find unique compiler warnings in make.log and write them to warnings.txt:
#   $ ./findlines.py make.log warning warnings.txt
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# main routine
# ------------------------------------------------------------------------------
def main():

    import argparse
    import sys, os
    import shlex

    parser = argparse.ArgumentParser(
        description="Find unique lines in <readfile> containing <key> "
        + "and write them to <outfile>",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("readfile", type=str, help="File to search")

    parser.add_argument("key", type=str, help="String to find in lines of <readfile>")

    parser.add_argument(
        "outfile", type=str, help="Output file for lines in <readfile> containing <key>"
    )

    # parse command line args
    args = parser.parse_args()

    # check if readfile exists
    if not os.path.isfile(args.readfile):
        print("Error:", args.readfile, "not found")
        return -1

    # check that outfile does not already exist
    if os.path.isfile(args.outfile):
        print("Error:", args.outfile, "already exists")
        return -1

    # create set to hold lines already seen
    lines_seen = set()

    with open(args.readfile, "r") as fr:
        with open(args.outfile, "w") as fw:
            # loop over every line in readfile
            for line in fr:
                # check if line contains the key
                if args.key in line:
                    # if not a duplicate line then
                    # add the line to the set and
                    # output it to the results file
                    if line not in lines_seen:
                        lines_seen.add(line)
                        print(line, file=fw)


# ------------------------------------------------------------------------------
# run the main routine
# ------------------------------------------------------------------------------
if __name__ == "__main__":
    import sys

    sys.exit(main())
