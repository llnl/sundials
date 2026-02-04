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

parser.add_argument("data_file", type=str, help="Data file to plot")

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
            ("y3", float),
            ("loc err est1", float),
            ("loc err est2", float),
            ("loc err est3", float),
            ("com loc err est1", float),
            ("com loc err est2", float),
            ("com loc err est3", float),
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

# parse log file
log = sunlog.log_file_to_list(args.log_file)

# get step data
_, times_s, vals_s = sunlog.get_history(log, 'h', "success")
_, times_etf, vals_etf = sunlog.get_history(log, 'h', "failed error test")
_, times_sf, vals_sf = sunlog.get_history(log, 'h', "failed solve")

_, times_1a, vals_1a = sunlog.get_history(log, 'compute-solution.err(:)')
_, times_1s, vals_1s = sunlog.get_history(log, 'compute-solution.err(:)', "success")
_, times_1f, vals_1f = sunlog.get_history(log, 'compute-solution.err(:)', "failed")

_, times_2a, vals_2a = sunlog.get_history(log, 'compute-solution.err2(:)')
_, times_2s, vals_2s = sunlog.get_history(log, 'compute-solution.err2(:)', "success")
_, times_2f, vals_2f = sunlog.get_history(log, 'compute-solution.err2(:)', "failed")

stiffness_ratio = data["stiffness_ratio"][0]
lambda_1 = data["lambda1"][0]
lambda_3 = data["lambda3"][0]
mag_lambda_1 = np.abs(lambda_1)
mag_lambda_3 = np.abs(lambda_3)
max_lambda = np.maximum(mag_lambda_1, mag_lambda_3)
min_lambda = np.minimum(mag_lambda_1, mag_lambda_3)

fig1, axes = plt.subplots(3, sharex=True, figsize=(12, 8))

axes[0].plot(data["t"], data["y1"], label="y1")
axes[0].plot(data["t"], data["y2"], label="y2")
axes[0].plot(data["t"], data["y3"], label="y3")
axes[0].legend(loc="best")
axes[0].set_xlabel("time")
axes[0].tick_params(labelbottom=True)
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
# axes[1].plot(data["t"], np.clip(np.abs(data["loc err est1"]), clip, None), label="y1")
# axes[1].plot(data["t"], np.clip(np.abs(data["loc err est2"]), clip, None), label="y2")
# axes[1].plot(data["t"], np.clip(np.abs(data["loc err est3"]), clip, None), label="y3")
# axes[1].plot(data["t"], np.clip(np.abs(data["com loc err est1"]), clip, None), linestyle="--", label="x1")
# axes[1].plot(data["t"], np.clip(np.abs(data["com loc err est2"]), clip, None), linestyle="--", label="x2")
# axes[1].plot(data["t"], np.clip(np.abs(data["com loc err est3"]), clip, None), linestyle="--", label="x3")

# data = np.array(vals_1s)
# for i in range(data.shape[1]):
#     axes[1].scatter(times_1s, np.abs(data[:, i]), marker='o', color=colors[i], label=f'Component {i}')
# data = np.array(vals_1f)
# for i in range(data.shape[1]):
#     axes[1].scatter(times_1f, np.abs(data[:, i]), marker='x', color=colors[i], label=f'Component {i}')

# data = np.array(vals_2s)
# for i in range(data.shape[1]):
#     axes[1].scatter(times_2s, np.abs(data[:, i]), marker='o', facecolors='none', edgecolors=colors[i], label=f'Component {i}')
# data = np.array(vals_2f)
# for i in range(data.shape[1]):
#     axes[1].scatter(times_2f, np.abs(data[:, i]), marker='X', facecolors='none', edgecolors=colors[i], label=f'Component {i}')

# colors = plt.cm.tab10.colors

# # Hairer and Wanner do companion / primary:
# # > 1 less stiff
# # < 1 more stiff
# # I'm going to flip it so larger -> stiffer
# # > 1 less stiff
# # < 1 more stiff
# data_success = np.array(vals_1s) / np.array(vals_2s)
# if len(vals_etf) > 0:
#     data_fail = np.array(vals_1f) / np.array(vals_2f)

# for i in range(data_success.shape[1]):
#     axes[1].scatter(times_1s, np.abs(data_success[:, i]), marker='o', facecolors='none', edgecolors=colors[i], label=f'y{i+1}')
# if len(vals_etf) > 0:
#     for i in range(data_fail.shape[1]):
#         axes[1].scatter(times_1f, np.abs(data_fail[:, i]), marker='X', facecolors='none', edgecolors=colors[i], label=f'y{i+1}')

# axes[1].axhline(1, color="gray", linestyle='--', alpha=0.7)

# legend_elements = []

# # Add component colors
# for i in range(data_success.shape[1]):
#     legend_elements.append(Line2D([0], [0], marker='s', color='w', markerfacecolor=colors[i],
#                                   label=f'y{i+1}'))

# # Create custom legend
# legend_elements.extend([
#     Line2D([0], [0], marker='o', color='w', markerfacecolor='none', markeredgecolor='gray', label='Successful step'),
#     Line2D([0], [0], marker='X', color='w', markerfacecolor='none', markeredgecolor='gray', label='Successful step')])

# axes[1].set_yscale("log")
# axes[1].legend(handles=legend_elements, loc="best", ncol=2)
# axes[1].set_xlabel("time")
# axes[1].grid(True, which="major", linestyle=":", alpha=0.5)

# title = "Absolute Ratio of Error Estimates " + r"($|e_p/e_s|$)"
# if args.title:
#     title = args.title + " " + title
# axes[1].set_title(title)
# axes[1].grid(True, which="major", linestyle=":", alpha=0.5)


colors = plt.cm.tab10.colors

# Hairer and Wanner do companion / primary:
# > 1 less stiff
# < 1 more stiff
# I'm going to flip it so larger -> stiffer
# > 1 less stiff
# < 1 more stiff
data_1 = np.array(vals_1a)
data_2 = np.array(vals_2a)

for i in range(data_1.shape[1]):
    axes[1].scatter(times_1a, np.abs(data_1[:, i]), marker='o', facecolors='none', edgecolors=colors[i], label=f'p y{i+1} ')
for i in range(data_2.shape[1]):
    axes[1].scatter(times_2a, np.abs(data_2[:, i]), marker='X', facecolors='none', edgecolors=colors[i], label=f's y{i+1}')

axes[1].axhline(1, color="gray", linestyle='--', alpha=0.7)

# legend_elements = []

# # Add component colors
# for i in range(data_success.shape[1]):
#     legend_elements.append(Line2D([0], [0], marker='s', color='w', markerfacecolor=colors[i],
#                                   label=f'y{i+1}'))

# # Create custom legend
# legend_elements.extend([
#     Line2D([0], [0], marker='o', color='w', markerfacecolor='none', markeredgecolor='gray', label='Successful step'),
#     Line2D([0], [0], marker='X', color='w', markerfacecolor='none', markeredgecolor='gray', label='Successful step')])

axes[1].set_yscale("log")
#axes[1].legend(handles=legend_elements, loc="best", ncol=2)
axes[1].legend(loc="best", ncol=2)
axes[1].set_xlabel("time")
axes[1].grid(True, which="major", linestyle=":", alpha=0.5)

title = "Absolute Ratio of Error Estimates " + r"($|e_p/e_s|$)"
if args.title:
    title = args.title + " " + title
axes[1].set_title(title)
axes[1].grid(True, which="major", linestyle=":", alpha=0.5)


axes[2].scatter(times_s, vals_s, label="success")
axes[2].scatter(times_etf, vals_etf, label="error test fail")
axes[2].scatter(times_sf, vals_sf, label="solver fail") # need to get solver fail in successful step
#axes[2].set_yscale("log")
axes[2].legend(loc="best")
axes[2].set_xlabel("time")
if args.title:
    axes[2].set_title(f"{args.title} Step History")
else:
    axes[2].set_title(f"Step History")
axes[2].grid(True, which="major", linestyle=":", alpha=0.5)

plt.tight_layout()

if args.save:
    plt.savefig(args.save, bbox_inches="tight")
else:
    plt.show()
