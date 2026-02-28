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
# Test all installed Makefiles or CMakeLists.txt files under a given directory
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# main routine
# -----------------------------------------------------------------------------


def main():

    import argparse
    import fnmatch
    import os
    import re
    import subprocess

    parser = argparse.ArgumentParser(
        description="Find and build installed examples",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("directory", type=str, help="Directory to search for build files")

    parser.add_argument("--cmake", action="store_true", help="CMake build")

    parser.add_argument("--test", action="store_true", help="Test builds")

    parser.add_argument("--clean", action="store_true", help="Clean builds")

    parser.add_argument(
        "--regex", type=str, help="Regular expression for filtering example directories"
    )

    parser.add_argument("-v", "--verbose", action="count", default=0, help="Verbose output")

    parser.add_argument("--failfast", action="store_true", help="Stop on first failure")

    # parse command line args
    args = parser.parse_args()

    # get current directory
    cwd = os.getcwd()

    # which build file to search for
    if args.cmake:
        filename = "CMakeLists.txt"
    else:
        filename = "Makefile"

    # get list of build files
    buildfiles = []
    for root, dirs, files in os.walk(args.directory):
        for fn in fnmatch.filter(files, filename):
            buildfiles.append(os.path.join(root, fn))

    # output list of build files
    if args.verbose > 0:
        print(f"Total files: {len(buildfiles)}")
    if args.verbose > 2:
        for bf in buildfiles:
            print(bf)

    # filter files
    if args.regex:
        regex = re.compile(args.regex)
        buildfiles = [bf for bf in buildfiles if re.search(regex, bf)]
        if args.verbose > 0:
            print(f"Total files (filtered): {len(buildfiles)}")
        if args.verbose > 2:
            for bf in buildfiles:
                print(bf)

    # total files to test
    total = len(buildfiles)

    # build examples
    fail = False
    errors = []
    for i, bf in enumerate(buildfiles):

        print(f"[{i+1} of {total}] Building: {bf}")

        # move to example directory
        os.chdir(os.path.dirname(bf))

        # clean and move on
        if args.clean:
            ret = subprocess.call(
                "make clean", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )
            # return to original directory
            os.chdir(cwd)
            continue

        # configure cmake if necessary
        configfail = False
        if args.cmake:
            if os.path.isfile("Makefile"):
                os.remove("Makefile")
            ret = subprocess.call(
                "cmake -DCMAKE_VERBOSE_MAKEFILE=ON .",
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            if args.verbose > 0:
                print(f"  Config return: {ret}")
            if ret != 0:
                errors.append(os.path.dirname(bf))
                configfail = True

        # make examples
        buildfail = False
        if not configfail:
            ret = subprocess.call(
                "make", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )
            if args.verbose > 0:
                print(f"  Build return: {ret}")
            if ret != 0:
                errors.append(os.path.dirname(bf))
                buildfail = True

        # test examples
        testfail = False
        if not configfail and not buildfail and args.test:
            ret = subprocess.call(
                "make test", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )
            if args.verbose > 0:
                print(f"  Test return: {ret}")
            if ret != 0:
                errors.append(os.path.dirname(bf))
                testfail = True

        # return to original directory
        os.chdir(cwd)

        # set test status
        if configfail or buildfail or testfail:
            fail = True

        # stop on first error
        if args.failfast and fail:
            break

    # print list of failed builds
    if fail:
        print("The following examples failed:")
        for err in errors:
            print(err)
        sys.exit(1)
    else:
        print("All builds successful.")


# -----------------------------------------------------------------------------
# run the main routine
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    sys.exit(main())
