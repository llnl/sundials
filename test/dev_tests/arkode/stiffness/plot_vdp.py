#!/usr/bin/env python3
"""
Plot EMT test
"""

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
    data = np.loadtxt(args.data_file,
                      dtype=[('t', float), ('y1', float), ('y2', float)])
except FileNotFoundError:
    print(f"Error: File '{args.data_file}' not found.")
    sys.exit(1)
except Exception as e:
    print(f"Error reading file '{args.data_file}': {e}")
    sys.exit(1)

fig1, ax1 = plt.subplots(figsize=(12, 8))

ax1.plot(data['t'], data['y1'], label='y1')
ax1.plot(data['t'], data['y2'], label='y2')
ax1.legend(loc='best')
ax1.set_xlabel('time')
ax1.set_title(args.title)
ax1.grid(True, which='both', linestyle=':', alpha=0.5)

plt.tight_layout()

if args.save:
    plt.savefig(args.save, bbox_inches="tight")
else:
    plt.show()
