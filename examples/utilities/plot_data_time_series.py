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
# Script to plot and compare time series data from multiple files
#
# The input data files must contain an M x (N + 1) matrix of values where the
# first column of each row is the output time and the remaining values are the
# N output quantities. Additionally the data files must contain header comment
# starting with "vars:" that lists the name for each data column i.e.,
#
# # vars: time quantity_0 quantity_1 . . . quantity_N-1
# t_0   q_0,0   q_0,1   q_0,2   . . . q_0,N-1
# t_1   q_1,0   q_1,1   q_1,2   . . . q_1,N-1
#   .     .       .       .       .     .
#   .     .       .       .       .     .
#   .     .       .       .       .     .
# t_M-1 q_M-1,0 q_M-1,1 q_M-1,2 . . . q_M-1,N-1
#
# where t_m is the m-th output time and q_m,n is the n-th quantity at the m-th
# output time.
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# main routine
# -----------------------------------------------------------------------------
def main():

    import argparse
    import matplotlib.pyplot as plt
    import numpy as np
    import shlex

    parser = argparse.ArgumentParser(description="""Plot data files""")

    parser.add_argument("quantity", type=str, help="""Quantity to plot""")

    parser.add_argument("datafiles", type=str, nargs="+", help="""Data files to plot""")

    # Plot display options

    parser.add_argument("--save", action="store_true", help="""Save figure to file""")

    parser.add_argument(
        "--labels", type=str, nargs="+", help="""Data file labels for plot legend"""
    )

    parser.add_argument("--title", type=str, help="""Plot title""")

    parser.add_argument("--xlabel", type=str, help="""x-axis label""")

    parser.add_argument("--ylabel", type=str, help="""y-axis label""")

    parser.add_argument("--grid", action="store_true", help="""Add grid to plot""")

    # Axis scaling
    logscale = parser.add_mutually_exclusive_group()

    logscale.add_argument("--logx", action="store_true", help="""Plot with log scale x-axis""")

    logscale.add_argument("--logy", action="store_true", help="""Plot with log scale y-axis""")

    logscale.add_argument("--loglog", action="store_true", help="""Use log scale x and y axes""")

    # Debugging options

    parser.add_argument("--debug", action="store_true", help="Enable debugging")

    # Parse command line args
    args = parser.parse_args()

    # Create figure and axes
    fig, ax = plt.subplots()

    for i, datafile in enumerate(args.datafiles):

        quantities = None
        with open(datafile) as fn:
            # read the file line by line
            for line in fn:
                # skip empty lines
                if not line.strip():
                    continue
                # exit after reading initial comment lines
                if "#" not in line:
                    break
                # split line into a list
                text = shlex.split(line)
                # extract quantity names
                if "vars:" in line:
                    quantities = text[2:]
                    continue

        # Check inputs
        if quantities is None:
            print("ERROR: quantity names not provided")
            sys.exit()

        if args.quantity not in quantities:
            print("ERROR: quantity not found")
            print(f"  Possible values: {quantities}")
            sys.exit()

        # Get data column index
        idx = quantities.index(args.quantity)

        # Load data
        data = np.loadtxt(datafile, dtype=np.double)

        if args.debug:
            print(np.shape(data))
            print(data)

        # Extract t and q data
        tdata = data[:, 0]  # first column has t values
        qdata = data[:, idx]  # remaining columns have q values

        # line colors: matplotlib.org/stable/tutorials/colors/colormaps.html
        # and colorbrewer2.org)
        if len(args.datafiles) < 22:
            colors = [
                "#d62728",
                "#1f77b4",
                "#2ca02c",
                "#9467bd",
                "#ff7f0e",
                "#8c564b",
                "#e377c2",
                "#7f7f7f",
                "#bcbd22",
                "#17becf",
                "#000000",
                "#ff9896",
                "#aec7e8",
                "#98df8a",
                "#c5b0d5",
                "#ffbb78",
                "#c49c94",
                "#f7b6d2",
                "#c7c7c7",
                "#dbdb8d",
                "#9edae5",
            ]
        else:
            print("ERROR: ncols > ncolors")
            sys.exit()

        # Set plot label for legend
        if args.labels:
            label = args.labels[i]
        else:
            label = None

        if args.logx or args.logy or args.loglog:
            ax.plot(tdata, np.abs(qdata), label=label, color=colors[i])
        else:
            ax.plot(tdata, qdata, label=label, color=colors[i])

    # Change axis scale
    if args.logx:
        ax.set_xscale("log")
    elif args.logy:
        ax.set_yscale("log")
    elif args.loglog:
        ax.set_xscale("log")
        ax.set_yscale("log")

    # Add title
    if args.title:
        plt.title(args.title)

    # Add x-axis label
    if args.xlabel:
        plt.xlabel(args.xlabel)
    else:
        plt.xlabel("time")

    # Add y-axis label
    if args.ylabel:
        plt.ylabel(args.ylabel)
    else:
        plt.ylabel(args.quantity.replace("_", " "))

    # Add legend
    if args.labels:
        ax.legend(bbox_to_anchor=(1, 0.95), loc="upper right")

    # Add grid
    if args.grid:
        plt.grid()

    # Save plot to file
    if args.save:
        plt.savefig("fig.pdf")
    else:
        plt.show()


# -----------------------------------------------------------------------------
# run the main routine
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    sys.exit(main())
