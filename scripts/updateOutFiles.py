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
# Script to update example output files for failed tests
#
# Example usage:
#   $ ./updateOutFiles.py ../build ../.
# -----------------------------------------------------------------------------


class colors:
    SUCCESS = "\033[92m"
    WARNING = "\033[93m"
    ERROR = "\033[91m"
    END = "\033[0m"


def print_success(msg):
    print(f"{colors.SUCCESS}{msg}{colors.END}")


def print_warning(msg):
    print(f"{colors.WARNING}{msg}{colors.END}")


def print_error(msg):
    print(f"{colors.ERROR}{msg}{colors.END}")


# -----------------------------------------------------------------------------
# main routine
# -----------------------------------------------------------------------------
def main():

    import argparse
    import os, shutil

    parser = argparse.ArgumentParser(description="Update output files")

    parser.add_argument("source", type=str, help="Full path of build directory to read files from")
    parser.add_argument(
        "destination", type=str, help="Full path of sundials location to write files to"
    )
    parser.add_argument("--all", "-a", action="store_true", help="Update all output files")
    parser.add_argument(
        "--copy", "-c", action="store_true", help="Copy file to destination if not found"
    )
    parser.add_argument("--verbose", "-v", action="count", default=0, help="Enable verbose output")

    # parse command line args
    args = parser.parse_args()

    # check inputs
    if not os.path.isdir(args.source):
        print_error(f"Error: could not find {args.source}")
        return -1

    if not os.path.isdir(args.destination):
        print_error(f"Error: could not find {args.destination}")
        return -1

    # check that the output directory exists
    output = os.path.join(args.source, "Testing", "output")
    if not os.path.isdir(output):
        print_error(f"Error: could not find {output}")
        return -1

    # create a list of output files for all tests run or just the failed tests
    output_files = []

    if args.all:
        # get names of all .out files
        for f in os.listdir(output):
            if f.endswith(".out"):
                output_files.append(f)
    else:
        failed = os.path.join(args.source, "Testing", "Temporary", "LastTestsFailed.log")

        # extract test names from list and append .out
        with open(failed, "r") as f:
            for line in f:
                output_files.append(line.split(":")[1].rstrip() + ".out")

    # if failed tests were found update the output files
    if output_files:
        num_output_files = len(output_files)
        paths = [
            os.path.join(args.destination, "examples"),
            os.path.join(args.destination, "test", "unit_tests"),
            args.destination,
        ]
        for idx, out_file in enumerate(output_files):
            if args.verbose > 0:
                print(f"[{idx + 1} of {num_output_files}] Test {out_file[:-4]}")
            # Some tests do not have an output file (e.g., unit tests)
            if not os.path.isfile(os.path.join(output, out_file)):
                print_warning(f"  Warning: did not find the output file {out_file}")
                continue
            if args.verbose > 1:
                print(f"  Searching for answer file {out_file}")
            found = False
            for p in paths:
                if not os.path.isdir(p):
                    continue
                if args.verbose == 2:
                    print(f"  Looking in {p}")
                for root, _, files in os.walk(p):
                    if args.verbose == 3:
                        print(f"  Looking in {root}")
                    if out_file in files:
                        if args.verbose == 1:
                            print(f"  Found file in {root}")
                        if args.verbose > 1:
                            print("  Found file")
                        shutil.copy(os.path.join(output, out_file), os.path.join(root, out_file))
                        found = True
                        if args.verbose > 0:
                            print_success(f"  Answer file updated")
                        break
                if found:
                    break
            if not found:
                print_warning(f"  Warning: did not find the answer file {out_file}")
                if args.copy:
                    print(f"  Copying {out_file} to {args.destination}")
                    shutil.copy(
                        os.path.join(output, out_file), os.path.join(args.destination, out_file)
                    )


# -----------------------------------------------------------------------------
# run the main routine
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    import sys

    sys.exit(main())
