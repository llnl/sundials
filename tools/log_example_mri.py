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
# Example script demonstrating how to use Python functions to extract and plot
# logs produced by the SUNLogger with an MRI method.
# -----------------------------------------------------------------------------


def main():

    import argparse
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    import numpy as np

    from suntools import logs as sunlog

    parser = argparse.ArgumentParser(description="Plots")

    parser.add_argument("logfiles", type=str, nargs="+", help="Log file to plot")

    parser.add_argument("--scatter", action="store_true", help="Use scatter plot for step sizes")

    parser.add_argument("--stats", action="store_true", help="Print step statistics")

    parser.add_argument("--labels", type=str, nargs="+", help="Labels for plot legend")

    parser.add_argument("--legend-loc", type=str, default="best", help="Legend location")

    parser.add_argument("--logx", action="store_true", help="Use log scale for x-axis")

    parser.add_argument("--logy", action="store_true", help="Use log scale for y-axis")

    parser.add_argument("--step-number", action="store_true", help="Plot value vs step number")

    parser.add_argument("--timescale", type=float, help="Time scaling factor")

    parser.add_argument("--stepscale", type=float, help="Step size scaling factor")

    parser.add_argument("--xlabel", type=str, help="X-axis label")

    parser.add_argument("--ylabel", type=str, help="Y-axis label")

    parser.add_argument(
        "--step-range",
        type=int,
        nargs=2,
        default=None,
        metavar=("LOWER_BOUND", "UPPER_BOUND"),
        help="Step range to plot",
    )

    parser.add_argument(
        "--time-range",
        type=float,
        nargs=2,
        default=None,
        metavar=("LOWER_BOUND", "UPPER_BOUND"),
        help="Time range to plot",
    )

    parser.add_argument(
        "--save",
        type=str,
        nargs="?",
        const="fig.pdf",
        default=None,
        metavar="FILE_NAME",
        help="Save figure to file",
    )

    # parse command line args
    args = parser.parse_args()

    # Enforce len(labels) == len(logfiles)
    if args.labels is not None and len(args.labels) != len(args.logfiles):
        parser.error(
            f"--labels expects {len(args.logfiles)} values " f"(you passed {len(args.labels)})"
        )

    # Delay making the plot axes until we know how many times levels to plot
    axes = None

    colors = plt.get_cmap("tab10")

    for log_idx, logfile in enumerate(args.logfiles):

        # parse log file
        log = sunlog.log_file_to_list(logfile)

        # extract step size history for all step attempts
        steps_a, times_a, vals_a = sunlog.get_history(
            log, "h", step_range=args.step_range, time_range=args.time_range, group_by_level=True
        )

        # extract the step size history for all failed steps
        steps_f, times_f, vals_f = sunlog.get_history(
            log,
            "h",
            "failed",
            step_range=args.step_range,
            time_range=args.time_range,
            group_by_level=True,
        )

        if args.stats:
            # extract step size history for all successful steps
            steps_s, times_s, vals_s = sunlog.get_history(
                log,
                "h",
                "success",
                step_range=args.step_range,
                time_range=args.time_range,
                group_by_level=True,
            )

        if args.step_number:
            x_a = steps_a
            x_f = steps_f
        else:
            x_a = times_a
            x_f = times_f

            if args.timescale:
                for level_idx in range(len(x_a)):
                    x_a[level_idx] = args.timescale * np.array(x_a[level_idx])
                    x_f[level_idx] = args.timescale * np.array(x_f[level_idx])

        if args.stepscale:
            for level_idx in range(len(vals_a)):
                vals_a[level_idx] = args.stepscale * np.array(vals_a[level_idx])
                vals_f[level_idx] = args.stepscale * np.array(vals_f[level_idx])

        # now that we've read the first log, make as many subplots as time levels
        if log_idx == 0:
            _, axes = plt.subplots(len(x_a), sharex=True)

        if args.stats:
            print()
            print(f"File: {args.logfiles[log_idx]}")

        for level_idx in range(len(x_a)):

            # plot step attempts
            if args.scatter:
                axes[level_idx].scatter(
                    x_a[level_idx],
                    vals_a[level_idx],
                    color=colors(log_idx),
                    marker=".",
                    zorder=0.1,
                    label="attempts",
                )
            else:
                axes[level_idx].plot(
                    x_a[level_idx],
                    vals_a[level_idx],
                    color=colors(log_idx),
                    marker=".",
                    zorder=0.1,
                    label="attempts",
                )

            # plot failed steps
            axes[level_idx].scatter(
                x_f[level_idx],
                vals_f[level_idx],
                facecolors=colors(log_idx),
                edgecolors="red",
                linewidth=1,
                marker="X",
                zorder=0.2,
                label="failed",
            )

            axes[level_idx].grid(alpha=0.3, linestyle="--")

            if args.ylabel:
                axes[level_idx].set_ylabel(args.ylabel)
            else:
                axes[level_idx].set_ylabel("step size")

            if len(x_a) == 2:
                if level_idx == 0:
                    axes[level_idx].set_title(f"Slow Time Scale")
                else:
                    axes[level_idx].set_title(f"Fast Time Scale")
            else:
                axes[level_idx].set_title(f"Level {level_idx}")

            if args.logx:
                axes[level_idx].set_xscale("log")
            if args.logy:
                axes[level_idx].set_yscale("log")

            if args.stats:
                nattempted = len(vals_a[level_idx])
                npassed = len(vals_s[level_idx])
                nfailed = len(vals_f[level_idx])
                print()
                print(f"Level {level_idx}")
                print(f"Step attempts: {nattempted:,}")
                print(f"Successful steps: {npassed:,} ({100 * npassed/nattempted:.2f}%)")
                print(f"Failed steps: {nfailed:,} ({100 * nfailed/nattempted:.2f}%)")
                print(f"Max step size: {np.max(vals_s[level_idx])}")
                print(f"Min step size: {np.min(vals_s[level_idx])}")
                print(f"Avg step size: {np.average(vals_s[level_idx])}")
                if level_idx > 0:
                    print(
                        f"Avg level {level_idx} attempts per level {level_idx - 1} attempts: {len(vals_a[level_idx])/len(vals_a[level_idx - 1]):.2f}"
                    )

        if args.xlabel:
            axes[len(x_a) - 1].set_xlabel(args.xlabel)
        else:
            if args.step_number:
                axes[len(x_a) - 1].set_xlabel("step")
            else:
                axes[len(x_a) - 1].set_xlabel("time")

    # number of logfiles
    nlogs = len(args.logfiles)

    # create legend handles with separate entries for line color and marker type
    handles = []
    if args.labels or nlogs > 1:
        handles = [Patch(color=colors(i)) for i in range(nlogs)]
    handles.extend(
        [
            Line2D([], [], marker=".", mec="black", mfc="black", ls=""),
            Line2D([], [], marker="X", mec="red", mfc="white", ls=""),
        ]
    )

    # create corresponding legend labels for line color and marker type
    labels = []
    if args.labels is not None:
        labels = args.labels
    elif nlogs > 1:
        labels = [f"file {i}" for i in range(nlogs)]
    labels.extend(["Attempted", "Failed"])

    # add the combined legend to the first subplot
    if args.legend_loc == "outside":
        axes[0].legend(
            handles, labels, bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0
        )
    else:
        axes[0].legend(handles, labels, loc=args.legend_loc)

    if args.save:
        plt.savefig(args.save, bbox_inches="tight")
    else:
        plt.tight_layout()
        plt.show()


# run the main routine
if __name__ == "__main__":
    import sys

    sys.exit(main())
