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
import glob
import argparse
import multiprocessing as mp

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
        "--jobid", dest="jobID", type=int, help="job id of the current run to identify .cali files"
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
    jobID = args.jobID
    threshold = args.threshold

    # Get available benchmarks
    benchFiles = glob.glob("%s/Benchmarking/*/*" % caliDir)

    if not os.path.exists(outPath):
        os.makedirs(outPath)
    outFile = open("%s/benchmark_output.out" % outPath, "w")

    # thread per file
    with mp.Pool() as pool:
        for res in pool.starmap(
            process_benchmark, [(jobID, release, releaseDir, i, threshold) for i in benchFiles]
        ):
            if res:
                outFile.write(res + "\n")
    outFile.close()

    outFile = open("%s/benchmark_output.out" % outPath, "r")
    try:
        outLines = outFile.readlines()
    finally:
        outFile.close()

    if len(outLines) == 0:
        return -1
    return 0


def process_benchmark(jobID, isRelease, releaseDir, benchmarkDir, threshold):
    # Get the current benchmark run
    benchmarkFiles = glob.glob("%s/*.cali" % benchmarkDir)
    # Don't compare if the run didn't include this benchmark
    if len(benchmarkFiles) == 0:
        return

    th_files = tt.Thicket.from_caliperreader(benchmarkFiles)
    curFilter = lambda x: x["job_id"] == jobID
    th_current = th_files.filter_metadata(curFilter)

    # Get the release caliper file
    cluster = th_current.metadata["cluster"].values[0]
    if isRelease:
        # Get the last release
        versionDirs = glob.glob("%s/%s/*" % (releaseDir, cluster))
        versionDirs.sort(key=os.path.getmtime, reverse=True)
        versionDir = versionDirs[1]
    else:
        # Get the release the run is a part of
        version = th_current.metadata["sundials_version"].values[0]
        versionDir = "%s/%s/%s" % (releaseDir, cluster, version)
    benchmarkName = th_current.metadata["env.TEST_NAME"].values[0]
    releaseFile = glob.glob(
        "%s/Benchmarking/*/%s/*.cali" % (versionDir, benchmarkName), recursive=True
    )
    th_compare = tt.Thicket.from_caliperreader(releaseFile)
    metrics = ["Max time/rank"]
    tt.mean(th_current, columns=metrics)
    tt.mean(th_compare, columns=metrics)

    ratio = (
        th_current.statsframe.dataframe["Max time/rank_mean"]
        / th_compare.statsframe.dataframe["Max time/rank_mean"]
    )

    tolerance = threshold / 100
    if 1 - ratio[0] < tolerance:
        return benchmarkName


if __name__ == "__main__":
    main()
