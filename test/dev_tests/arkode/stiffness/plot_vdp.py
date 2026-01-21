#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Location of suntools directory
sys.path.append(os.path.join(os.environ['SUNDIALS_REPO'], "tools"))
from suntools import logs as sunlog

parser = argparse.ArgumentParser()

parser.add_argument("data_file", type=str, help="File to plot")

parser.add_argument("log_file", type=str, help="Log file to plot")

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
            ("loc err est1", float),
            ("loc err est2", float),
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

# parse log file
log = sunlog.log_file_to_list(args.log_file)

# get step data
_, times_s, vals_s = sunlog.get_history(log, 'h', "success")
_, times_etf, vals_etf = sunlog.get_history(log, 'h', "failed error test")
_, times_sf, vals_sf = sunlog.get_history(log, 'h', "failed solve")

fig1, axes = plt.subplots(5, sharex=True, figsize=(18, 12))

axes[0].plot(data["t"], data["y1"], label="y1", linewidth=2)
axes[0].plot(data["t"], data["y2"], label="y2", linewidth=2, alpha=0.7)
axes[0].set_yscale('symlog') # symmetric logarithmic
axes[0].legend(loc="best")
axes[0].set_xlabel("time")
axes[0].tick_params(labelbottom=True)
if args.title:
    axes[0].set_title(f"{args.title} Solution")
else:
    axes[0].set_title(f"Solution")
axes[0].grid(True, which="both", linestyle=":", alpha=0.5)

axes[1].plot(data["t"], np.abs(data["lambda1"]), label=r"$|\lambda_1|$", linewidth=2)
axes[1].plot(data["t"], np.abs(data["lambda2"]), label=r"$|\lambda_2|$", linewidth=2, alpha=0.7)
axes[1].set_yscale("log")
axes[1].legend(loc="best")
axes[1].set_xlabel("time")
axes[1].tick_params(labelbottom=True)
if args.title:
    axes[1].set_title(f"{args.title} Eigenvalue Magnitudes")
else:
    axes[1].set_title(f"Eigenvalue Magnitudes")
axes[1].grid(True, which="major", linestyle=":", alpha=0.5)

axes[2].plot(data["t"], data["stiffness ratio"])
axes[2].set_yscale("log")
axes[2].set_xlabel("time")
axes[2].tick_params(labelbottom=True)
if args.title:
    axes[2].set_title(f"{args.title} Stiffness Ratio")
else:
    axes[2].set_title(f"Stiffness Ratio")
axes[2].grid(True, which="major", linestyle=":", alpha=0.5)

clip = np.finfo(np.float64).eps
axes[3].plot(data["t"], np.clip(np.abs(data["loc err est1"]), clip, None), label="y1", linewidth=2)
axes[3].plot(data["t"], np.clip(np.abs(data["loc err est2"]), clip, None), label="y2", linewidth=2, alpha=0.7)
axes[3].set_yscale("log")
axes[3].legend(loc="best")
axes[3].set_xlabel("time")
if args.title:
    axes[3].set_title(f"{args.title} Local Error Estimates")
else:
    axes[3].set_title(f"Local Error Estimates")
axes[3].grid(True, which="major", linestyle=":", alpha=0.5)

axes[4].scatter(times_s, vals_s, label="success")
axes[4].scatter(times_etf, vals_etf, label="error test fail")
axes[4].scatter(times_sf, vals_sf, label="solver fail") # need to get solver fail in successful step
#axes[4].set_yscale("log")
axes[4].legend(loc="best")
axes[4].set_xlabel("time")
if args.title:
    axes[4].set_title(f"{args.title} Step History")
else:
    axes[4].set_title(f"Step History")
axes[4].grid(True, which="major", linestyle=":", alpha=0.5)

plt.tight_layout()

if args.save:
    plt.savefig(args.save, bbox_inches="tight")
else:
    plt.show()
