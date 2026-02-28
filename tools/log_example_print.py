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
# Example script that prints a parsed log file
# -----------------------------------------------------------------------------


def main():

    import argparse

    from suntools import logs as sunlog

    parser = argparse.ArgumentParser(description="Plots")

    parser.add_argument("logfile", type=str, help="Log file to print")

    # parse command line args
    args = parser.parse_args()

    # parse log and print
    sunlog.print_log(sunlog.log_file_to_list(args.logfile))


# run the main routine
if __name__ == "__main__":
    import sys

    sys.exit(main())
