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
            ("y3", float),
            ("loc err est1", float),
            ("loc err est2", float),
            ("loc err est3", float),
            ("lambda1", complex),
            ("lambda2", complex),
            ("lambda3", complex),
            ("stiffness_ratio", float),
        ],
    )
except FileNotFoundError:
    print(f"Error: File '{args.data_file}' not found.")
    sys.exit(1)
except Exception as e:
    print(f"Error reading file '{args.data_file}': {e}")
    sys.exit(1)

stiffness_ratio = data["stiffness_ratio"][0]
lambda_1 = data["lambda1"][0]
lambda_3 = data["lambda3"][0]
mag_lambda_1 = np.abs(lambda_1)
mag_lambda_3 = np.abs(lambda_3)
max_lambda = np.maximum(mag_lambda_1, mag_lambda_3)
min_lambda = np.minimum(mag_lambda_1, mag_lambda_3)

fig1, axes = plt.subplots(2, figsize=(12, 8))

axes[0].plot(data["t"], data["y1"], label="y1")
axes[0].plot(data["t"], data["y2"], label="y2")
axes[0].plot(data["t"], data["y3"], label="y3")
axes[0].legend(loc="best")
axes[0].set_xlabel("time")
# title = (
#     "Solution ("
#     + r"$\lambda_{{1,2}}$"
#     + f" = {lambda_1.real}"
#     + r"$\pm$"
#     + f"{lambda_1.imag}, "
#     + r"$\lambda_{{3}}$"
#     + f" = {lambda_3.real}, "
#     + f"Stiffness Ratio = {max_lambda:.2f} / {min_lambda:.2f} = {stiffness_ratio:.2f})"
# )
title = (
    "Solution ("
    + r"$|\lambda_{{1,2}}|$ = "
    + f"{mag_lambda_1:.2f}, "
    + r"$|\lambda_{{3}}|$ = "
    + f"{mag_lambda_3:.2f}, "
    + f"Stiffness Ratio = {stiffness_ratio:.2f})"
)
if args.title:
    axes[0].set_title(f"{args.title} " + title)
else:
    axes[0].set_title(title)
axes[0].grid(True, which="both", linestyle=":", alpha=0.5)

clip = np.finfo(np.float64).eps
axes[1].plot(data["t"], np.clip(np.abs(data["loc err est1"]), clip, None), label="y1")
axes[1].plot(data["t"], np.clip(np.abs(data["loc err est2"]), clip, None), label="y2")
axes[1].plot(data["t"], np.clip(np.abs(data["loc err est3"]), clip, None), label="y3")
axes[1].set_yscale("log")
axes[1].legend(loc="best")
axes[1].set_xlabel("time")
if args.title:
    axes[1].set_title(f"{args.title} Local Error Estimates")
else:
    axes[1].set_title(f"Local Error Estimates")
axes[1].grid(True, which="major", linestyle=":", alpha=0.5)

plt.tight_layout()

if args.save:
    plt.savefig(args.save, bbox_inches="tight")
else:
    plt.show()
