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
# git diff --numstat | awk '{ added += $1; removed += $2 } END { print "+" added " -" removed }'
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
    from subprocess import call

    parser = argparse.ArgumentParser(
        description="Count the number lines changed in a certain path",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("key1", type=str, help="First git identifier")

    parser.add_argument("key2", type=str, help="Second git identifier")

    parser.add_argument(
        "patterns", type=str, nargs="+", help="File path patterns to look for in diff"
    )

    # parse command line args
    args = parser.parse_args()

    # git stats from git
    cmd = "git diff --numstat " + args.key1 + " " + args.key2 + " > tmp.txt"
    print(cmd)
    call(cmd, shell=True)

    totaladded = 0
    totalremoved = 0
    netchange = 0

    with open("tmp.txt", "r") as fr:
        for line in fr:
            split_line = shlex.split(line)
            added = split_line[0]
            removed = split_line[1]
            filename = split_line[2]
            if any(x in filename for x in args.patterns):
                if (added == "-") or (removed == "-"):
                    continue
                else:
                    print(added, removed, filename)
                    totaladded += int(added)
                    totalremoved += int(removed)

    netchange = totaladded - totalremoved
    print(totaladded, totalremoved, netchange)


# ------------------------------------------------------------------------------
# run the main routine
# ------------------------------------------------------------------------------
if __name__ == "__main__":
    import sys

    sys.exit(main())
