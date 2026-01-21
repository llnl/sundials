#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser()

parser.add_argument("data_file", type=str, help="File to plot")

parser.add_argument("--title", type=str, default=None, help="Plot title")

parser.add_argument(
    "--save",
    type=str,
    nargs="?",
    const="fig.pdf",
    default=None,
    metavar="FILE_NAME",
    help="Save figure to file",
)

args = parser.parse_args()

# Read the data file
try:
    data = np.loadtxt(
        args.data_file,
        dtype=[
            ("t", float),
            ("y1", float),
            ("y2", float),
            ("lambda1", complex),
            ("lambda2", complex),
            ("stiffness ratio", float),
        ],
    )
except FileNotFoundError:
    print(f"Error: File '{args.data_file}' not found.")
    sys.exit(1)
except Exception as e:
    print(f"Error reading file '{args.data_file}': {e}")
    sys.exit(1)

fig1, axes = plt.subplots(2, figsize=(12, 8))

axes[0].plot(data["t"], data["y1"], label="y1")
axes[0].plot(data["t"], data["y2"], label="y2")
axes[0].legend(loc="best")
axes[0].set_xlabel("time")
if args.title:
    axes[0].set_title(f"{args.title} Solution")
else:
    axes[0].set_title(f"Solution")
axes[0].grid(True, which="both", linestyle=":", alpha=0.5)

axes[1].plot(data["t"], data["stiffness ratio"])
axes[1].set_yscale("log")
axes[1].set_xlabel("time")
if args.title:
    axes[1].set_title(f"{args.title} Stiffness Ratio")
else:
    axes[1].set_title(f"Stiffness Ratio")
axes[1].grid(True, which="major", linestyle=":", alpha=0.5)

plt.tight_layout()

if args.save:
    plt.savefig(args.save, bbox_inches="tight")
else:
    plt.show()
