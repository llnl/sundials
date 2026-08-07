#!/usr/bin/env python3
# ----------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ----------------------------------------------------------------
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
# ----------------------------------------------------------------
# matplotlib-based plotting script for ark_kepler.c
# ----------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(
    description="Script for plotting the energy, angular momentum, and phase space solution for ark_kepler.c"
)
parser.add_argument("output_times", help="file with the output times")
parser.add_argument("solution", help="file with the solution")
parser.add_argument("conserved_quantities", help="file with conserved quantities")

args = parser.parse_args()

t = np.loadtxt(args.output_times, dtype=np.float64)
y = np.loadtxt(args.solution, dtype=np.float64)
y = np.reshape(y, (y.shape[0] // 4, 4))

plt.figure(dpi=200)
plt.plot(y[:, 0], y[:, 1])
plt.savefig("ark_kepler_phase.png")
plt.close()

conserved = np.loadtxt(args.conserved_quantities, delimiter=",", dtype=np.float64)
energy = conserved[:, 0]
energy_0 = conserved[0, 0]
L = conserved[:, 1]
L_0 = conserved[0, 1]

plt.figure(dpi=200)
plt.title("Energy")
plt.plot(t, np.abs(energy))
plt.ylabel("H(t,p,q)")
plt.xlabel("<---  t  --->")
plt.xscale("log")
plt.savefig("ark_kepler_energy.png")
plt.close()

plt.figure(dpi=200)
plt.title("Momentum")
plt.plot(t, L)
plt.ylabel("L(t,p,q)")
plt.xlabel("<---  t  --->")
plt.xscale("log")
plt.savefig("ark_kepler_momentum.png")
plt.close()

#
#  Time plot.
#
plt.figure(dpi=200)
plt.plot(t, y[:, 0], linewidth=2)
plt.plot(t, y[:, 1], linewidth=2)
plt.plot(t, y[:, 2], linewidth=2)
plt.plot(t, y[:, 3], linewidth=2)
plt.grid(True)
plt.legend(["P", "P'", "Q", "Q'"])
plt.xlabel("<---  t  --->")
plt.ylabel("<---  y(1:4)  --->")
plt.title("Solution in Time")
plt.savefig("ark_kepler_plot.png")
plt.close()

#
#  Phase plot.
#
plt.figure(dpi=200)
plt.plot(y[:, 0], y[:, 1], linewidth=0.1)
plt.grid(True)
plt.xlabel("<---  y1  --->")
plt.ylabel("<---  y2  --->")
plt.title("Phase Plot")
plt.savefig("ark_kepler_phase.png")
plt.close()
