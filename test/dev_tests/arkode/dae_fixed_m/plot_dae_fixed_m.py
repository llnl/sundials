#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

from matplotlib.lines import Line2D

# Location of suntools directory
sys.path.append(os.path.join(os.environ['SUNDIALS_REPO'], "tools"))
from suntools import logs as sunlog

parser = argparse.ArgumentParser()

parser.add_argument("data_file_1", type=str, help="Data file to plot")
parser.add_argument("data_file_2", type=str, help="Data file to plot")

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

# Read the data files
try:
    data1 = np.loadtxt(
        args.data_file_1,
        dtype=[
            ("t", float),
            ("y1", float),
            ("y2", float),
            ("y3", float),
        ],
    )
except FileNotFoundError:
    print(f"Error: File '{args.data_file}' not found.")
    sys.exit(1)
except Exception as e:
    print(f"Error reading file '{args.data_file}': {e}")
    sys.exit(1)

try:
    data2 = np.loadtxt(
        args.data_file_2,
        dtype=[
            ("t", float),
            ("y1", float),
            ("y2", float),
            ("y3", float),
        ],
    )
except FileNotFoundError:
    print(f"Error: File '{args.data_file}' not found.")
    sys.exit(1)
except Exception as e:
    print(f"Error reading file '{args.data_file}': {e}")
    sys.exit(1)

fig1, axes = plt.subplots(2, sharex=True, figsize=(12, 8))

axes[0].plot(data1["t"], data1["y1"], linewidth=2, marker='.', label="y1")
axes[0].plot(data1["t"], data1["y2"], linewidth=2, marker='.', label="y2")
axes[0].plot(data1["t"], data1["y3"], linewidth=2, marker='.', label="y3")

axes[0].plot(data2["t"], data1["y1"], linewidth=2, marker=".", linestyle="--", label="y1")
axes[0].plot(data2["t"], data1["y2"], linewidth=2, marker=".", linestyle="--", label="y2")
axes[0].plot(data2["t"], data1["y3"], linewidth=2, marker=".", linestyle="--", label="y3")

axes[0].legend(loc="best")
axes[0].set_xlabel("time")
axes[0].tick_params(labelbottom=True)
axes[0].set_title(f"Solution")
axes[0].grid(True, which="both", linestyle=":", alpha=0.5)

axes[1].plot(data1["t"], data1["y1"] - data2["y1"], label="y1")
axes[1].plot(data1["t"], data1["y2"] - data2["y2"], label="y2")
axes[1].plot(data1["t"], data1["y3"] - data2["y3"], label="y3")

axes[1].legend(loc="best")
axes[1].set_xlabel("time")
axes[1].set_title(f"Solution Difference")
axes[1].grid(True, which="both", linestyle=":", alpha=0.5)

plt.tight_layout()

if args.save:
    plt.savefig(args.save, bbox_inches="tight")
else:
    plt.show()
