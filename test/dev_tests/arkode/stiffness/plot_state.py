#!/usr/bin/env python3

# ./plot_state.py data.text --labels '$i_2$' '$v_{C_3}$' '$i_4$' '$v_{C_4}$'

import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser()

parser.add_argument("data_file", type=str, help="File to plot")

parser.add_argument("--labels", type=str, nargs='+', help="File to plot")

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
    with open(args.data_file, "r") as f:
        header = f.readline().strip().split()
    data = np.loadtxt(args.data_file, skiprows=1)
except FileNotFoundError:
    print(f"Error: File '{args.data_file}' not found.")
    sys.exit(1)
except Exception as e:
    print(f"Error reading file '{args.data_file}': {e}")
    sys.exit(1)

time = data[:,0]
state = data[:,1:]

if args.labels:
    labels = args.labels
else:
    labels = header[1:]

_, ax = plt.subplots()

ax.plot(time, state, linewidth=2, label=labels)
ax.set_xlabel("time")
ax.tick_params(labelbottom=True)
ax.grid(True, which="both", linestyle=":", alpha=0.5)
ax.legend(loc="best")

plt.tight_layout()

if args.save:
    plt.savefig(args.save, bbox_inches="tight")
else:
    plt.show()
