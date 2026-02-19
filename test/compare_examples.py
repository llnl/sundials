#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): Yu Pan @ LLNL
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

import os
import subprocess
import sys
import re
import glob
import argparse
import multiprocessing as mp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import hatchet as ht
import thicket as tt


def main():
    parser = argparse.ArgumentParser(
        description="Compare Sundials performance results against previous results"
    )

    parser.add_argument(
        "--release",
        dest="release",
        action="store_true",
        help="indicate if the current run to process is a release",
    )

    parser.add_argument(
        "--calidir",
        dest="caliDir",
        type=str,
        help="path to directory containing caliper files",
        default="/usr/workspace/sundials/califiles",
    )

    parser.add_argument(
        "--releasedir",
        dest="releaseDir",
        type=str,
        help="path to directory containing release caliper files",
        default="/usr/workspace/sundials/califiles/Release",
    )

    parser.add_argument(
        "--outpath",
        dest="outPath",
        type=str,
        help="path to directory to write results to",
        default="/dev/null",
    )

    parser.add_argument(
        "--threshold",
        dest="threshold",
        type=float,
        help="the percentage threshold in performance difference that indicates a regression",
        default=2.0,
    )

    args = parser.parse_args()

    release = args.release
    releaseDir = args.releaseDir
    caliDir = args.caliDir
    outPath = args.outPath
    threshold = args.threshold

    # Get the latest test run
    runDirs = glob.glob("%s/Testing/*" % caliDir, recursive=True)
    runDirs.sort(key=os.path.getmtime, reverse=True)
    runDir = runDirs[0]

    runFile = glob.glob(runDir)[0]
    th_temp = tt.Thicket.from_caliperreader(runFile)
    cluster = th_temp.metadata["cluster"]
    # get machine from the file
    if release:
        # Compare against the last release
        versionDirs = glob.glob("%s/%s/*" % (releaseDir, cluster))
        versionDirs.sort(key=os.path.getmtime, reverse=True)
        versionDir = versionDirs[1]
    else:
        # Compare against the release the run is a part of
        version = th_temp.metadata["sundials_version"].values[0]
        versionDir = "%s/%s/%s" % (releaseDir, cluster, version)

    # Gather files to process
    runFiles = glob.glob("%s/*.cali" % (runDir))

    if not os.path.exists(outPath):
        os.makedirs(outPath)
    outFile = open("%s/output.out" % outPath, "w")

    # Compare test results against past runs. If a test performs below a threshold, output test name to outFile.
    with mp.Pool() as pool:
        for res in pool.starmap(
            compare_against_release, [(versionDir, i, threshold) for i in runFiles]
        ):
            if res:
                outFile.write(res + "\n")
    outFile.close()

    outFile = open("%s/example_output.out" % outPath, "r")
    try:
        outLines = outFile.readlines()
    finally:
        outFile.close()

    if len(outLines) == 0:
        return -1
    return 0


def compare_against_release(releaseDir, file, threshold):
    th = tt.Thicket.from_caliperreader(file)

    testName = th.metadata["env.TEST_NAME"].values[0]

    # Gather release run
    releaseFile = glob.glob("%s/Testing/*/%s.*.cali" % (releaseDir, testName), recursive=True)
    th_release = tt.Thicket.from_caliperreader(releaseFile)

    metrics = ["Max time/rank"]
    tt.mean(th_release, columns=metrics)
    tt.mean(th, columns=metrics)

    ratio = (
        th.statsframe.dataframe["Max time/rank_mean"]
        / th_release.statsframe.dataframe["Max time/rank_mean"]
    )
    print(ratio[0])
    tolerance = threshold / 100
    if 1 - ratio[0] < tolerance:
        return testName


if __name__ == "__main__":
    main()
