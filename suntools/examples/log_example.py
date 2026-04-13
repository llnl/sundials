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
# logs produced by the SUNLogger for adaptive integrators.
# -----------------------------------------------------------------------------


def main():

    import argparse
    import matplotlib.pyplot as plt
    import matplotlib.ticker as tik

    from suntools import logs as sunlog

    parser = argparse.ArgumentParser(description="Plots")

    parser.add_argument("logfiles", type=str, nargs="+", help="Log files to plot")

    parser.add_argument(
        "--val",
        type=str,
        default="h",
        choices=["h", "q", "dsm"],
        help="Value to plot (default: %(default)s)",
    )

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

    parser.add_argument("--step-number", action="store_true", help="Plot value vs step number")

    parser.add_argument("--scatter", action="store_true", help="Create scatter plot")

    parser.add_argument("--logx", action="store_true", help="Use log scale for x-axis")

    parser.add_argument("--logy", action="store_true", help="Use log scale for y-axis")

    parser.add_argument("--labels", type=str, nargs="+", help="Plot labels")

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

    fig, ax = plt.subplots()
    colors = plt.get_cmap("tab10")

    for idx, lf in enumerate(args.logfiles):

        # parse log file
        log = sunlog.log_file_to_list(lf)

        # get successful step data
        steps_s, times_s, vals_s = sunlog.get_history(
            log, args.val, "success", step_range=args.step_range, time_range=args.time_range
        )

        # get data for error test failures
        steps_etf, times_etf, vals_etf = sunlog.get_history(
            log,
            args.val,
            "failed error test",
            step_range=args.step_range,
            time_range=args.time_range,
        )

        # get data for solver failures
        steps_sf, times_sf, vals_sf = sunlog.get_history(
            log, args.val, "failed solve", step_range=args.step_range, time_range=args.time_range
        )

        # plot log data
        if args.step_number:
            x_s = steps_s
            x_etf = steps_etf
            x_sf = steps_sf
        else:
            x_s = times_s
            x_etf = times_etf
            x_sf = times_sf

        if len(args.logfiles) == 1:
            s_color = "green"
            etf_color = "red"
            sf_color = "darkorange"
        else:
            s_color = colors(idx)
            etf_color = s_color
            sf_color = s_color

        if args.labels:
            s_label = f"{args.labels[idx]} successful"
            etf_label = f"{args.labels[idx]} error test failed"
            sf_label = f"{args.labels[idx]} solver failed"
        else:
            s_label = "successful"
            etf_label = "error test failed"
            sf_label = "solver failed"

        # plot successful data
        if args.scatter:
            ax.scatter(x_s, vals_s, color=s_color, marker="o", label=s_label, zorder=0.1)
        else:
            ax.plot(x_s, vals_s, color=s_color, marker=".", label=s_label, zorder=0.1)

        # always add failures as scatter plot
        ax.scatter(x_etf, vals_etf, color=etf_color, marker="x", label=etf_label, zorder=0.2)

        ax.scatter(x_sf, vals_sf, color=sf_color, marker="d", label=sf_label, zorder=0.2)

    if args.logx:
        ax.set_xscale("log")
    if args.logy:
        ax.set_yscale("log")

    if args.step_number:
        ax.set_xlabel("step")
    else:
        ax.set_xlabel("time")

    if args.val == "h":
        ax.set_ylabel("step size")
    elif args.val == "q":
        ax.set_ylabel("order")
        ax.yaxis.set_major_locator(tik.MaxNLocator(integer=True))
    elif args.val == "dsm":
        ax.set_ylabel("LTE estimate")

    ax.legend(loc="best")

    ax.grid(alpha=0.3, linestyle="--")

    if args.save:
        plt.savefig(args.save, bbox_inches="tight")
    else:
        plt.show()


# run the main routine
if __name__ == "__main__":
    import sys

    sys.exit(main())
