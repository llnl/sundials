#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ UMBC
# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------
# matplotlib-based plotting script for the serial ark_sod_lsrk example
# ------------------------------------------------------------------------------

# imports
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# data file name
datafile = "sod.out"

# return with an error if the file does not exist
if not os.path.isfile(datafile):
    msg = "Error: file " + datafile + " does not exist"
    sys.exit(msg)

# read solution file, storing each line as a string in a list
with open(datafile, "r") as file:
    lines = file.readlines()

    # extract header information
    title = lines.pop(0)
    nvar = int((lines.pop(0).split())[2])
    varnames = lines.pop(0)
    nt = int((lines.pop(0).split())[2])
    nx = int((lines.pop(0).split())[2])
    xl = float((lines.pop(0).split())[2])
    xr = float((lines.pop(0).split())[2])

    # allocate solution data as 2D Python arrays
    t = np.zeros((nt), dtype=float)
    rho = np.zeros((nt, nx), dtype=float)
    mx = np.zeros((nt, nx), dtype=float)
    my = np.zeros((nt, nx), dtype=float)
    mz = np.zeros((nt, nx), dtype=float)
    et = np.zeros((nt, nx), dtype=float)
    x = np.linspace(xl, xr, nx)

    # store remaining data into numpy arrays
    for it in range(nt):
        line = (lines.pop(0)).split()
        t[it] = line.pop(0)
        for ix in range(nx):
            rho[it, ix] = line.pop(0)
            mx[it, ix] = line.pop(0)
            my[it, ix] = line.pop(0)
            mz[it, ix] = line.pop(0)
            et[it, ix] = line.pop(0)

    gamma = 1.4
    u = mx / rho
    p = (gamma - 1.0) * (et - (mx * mx + my * my + mz * mz) / (2.0 * rho))


# utility routines for computing the analytical solution
def fsecant(p4, p1, p5, rho1, rho5, gamma):
    """
    Utility routine for exact_Riemann function below
    """
    z = p4 / p5 - 1.0
    c1 = np.sqrt(gamma * p1 / rho1)
    c5 = np.sqrt(gamma * p5 / rho5)
    gm1 = gamma - 1.0
    gp1 = gamma + 1.0
    g2 = 2.0 * gamma
    fact = gm1 / g2 * (c5 / c1) * z / np.sqrt(1.0 + gp1 / g2 * z)
    fact = (1.0 - fact) ** (g2 / gm1)
    return p1 * fact - p4


def exact_Riemann(t, x, xI):
    """
    Exact 1D Riemann problem solver (retrieves domain from EulerData structure),
    based on Fortran code at http://cococubed.asu.edu/codes/riemann/exact_riemann.f

    Inputs: (t,x) location for desired solution,
            xI location of discontinuity at t=0,
            gamma parameter for gas equation of state
    Outputs: rho, u, p (density, velocity, and pressure at (t,x))
    """

    # begin solution
    rho1 = 1.0
    p1 = 1.0
    u1 = 0.0
    rho5 = 0.125
    p5 = 0.1
    u5 = 0.0

    # solve for post-shock pressure by secant method initial guesses
    p40 = p1
    p41 = p5
    f0 = fsecant(p40, p1, p5, rho1, rho5, gamma)
    itmax = 50
    eps = 1.0e-14
    for iter in range(itmax):
        f1 = fsecant(p41, p1, p5, rho1, rho5, gamma)
        if f1 == f0:
            break
        p4 = p41 - (p41 - p40) * f1 / (f1 - f0)
        if (np.abs(p4 - p41) / np.abs(p41)) < eps:
            break
        p40 = p41
        p41 = p4
        f0 = f1
        if iter == itmax - 1:
            raise ValueError("exact_Riemann iteration failed to converge")

    # compute post-shock density and velocity
    z = p4 / p5 - 1.0
    c5 = np.sqrt(gamma * p5 / rho5)

    gm1 = gamma - 1.0
    gp1 = gamma + 1.0

    fact = np.sqrt(1.0 + 0.5 * gp1 * z / gamma)

    u4 = c5 * z / (gamma * fact)
    rho4 = rho5 * (1.0 + 0.5 * gp1 * z / gamma) / (1.0 + 0.5 * gm1 * z / gamma)

    # shock speed
    w = c5 * fact

    # compute values at foot of rarefaction
    p3 = p4
    u3 = u4
    rho3 = rho1 * (p3 / p1) ** (1.0 / gamma)

    # compute positions of waves
    c1 = np.sqrt(gamma * p1 / rho1)
    c3 = np.sqrt(gamma * p3 / rho3)

    xsh = xI + w * t
    xcd = xI + u3 * t
    xft = xI + (u3 - c3) * t
    xhd = xI - c1 * t

    # compute solution as a function of position
    if x < xhd:
        rho = rho1
        p = p1
        u = u1
    elif x < xft:
        u = 2.0 / gp1 * (c1 + (x - xI) / t)
        fact = 1.0 - 0.5 * gm1 * u / c1
        rho = rho1 * fact ** (2.0 / gm1)
        p = p1 * fact ** (2.0 * gamma / gm1)
    elif x < xcd:
        rho = rho3
        p = p3
        u = u3
    elif x < xsh:
        rho = rho4
        p = p4
        u = u4
    else:
        rho = rho5
        p = p5
        u = u5

    # return with success
    return rho, u, p


# generate analytical solutions over same mesh and times as loaded from data file
rhotrue = np.zeros((nt, nx), dtype=float)
utrue = np.zeros((nt, nx), dtype=float)
ptrue = np.zeros((nt, nx), dtype=float)
for it in range(nt):
    for ix in range(nx):
        rhotrue[it, ix], utrue[it, ix], ptrue[it, ix] = exact_Riemann(t[it], x[ix], 0.5)

#   plot defaults: increase default font size, increase plot width, enable LaTeX rendering
plt.rc("font", size=15)
plt.rcParams["figure.figsize"] = [7.2, 4.8]
plt.rcParams["text.usetex"] = True
plt.rcParams["figure.constrained_layout.use"] = True

#   subplots with time snapshots of the density, x-velocity, and pressure
fig = plt.figure(figsize=(10, 5))
gs = GridSpec(3, 3, figure=fig)
ax00 = fig.add_subplot(gs[0, 0])  # left column
ax10 = fig.add_subplot(gs[1, 0])
ax20 = fig.add_subplot(gs[2, 0])
ax01 = fig.add_subplot(gs[0, 1])  # middle column
ax11 = fig.add_subplot(gs[1, 1])
ax21 = fig.add_subplot(gs[2, 1])
ax02 = fig.add_subplot(gs[0, 2])  # right column
ax12 = fig.add_subplot(gs[1, 2])
ax22 = fig.add_subplot(gs[2, 2])
it = 0
ax00.plot(x, rho[it, :], "-b", x, rhotrue[it, :], ":k")
ax10.plot(x, u[it, :], "-b", x, utrue[it, :], ":k")
ax20.plot(x, p[it, :], "-b", x, ptrue[it, :], ":k")
ax00.set_title(r"$t =$ " + repr(t[it]).zfill(3))
ax00.set_ylabel(r"$\rho$")
ax10.set_ylabel(r"$v_x$")
ax20.set_ylabel(r"$p$")
ax20.set_xlabel(r"$x$")
it = nt // 2
ax01.plot(x, rho[it, :], "-b", x, rhotrue[it, :], ":k")
ax11.plot(x, u[it, :], "-b", x, utrue[it, :], ":k")
ax21.plot(x, p[it, :], "-b", x, ptrue[it, :], ":k")
ax01.set_title(r"$t =$ " + repr(t[it]).zfill(3))
ax21.set_xlabel(r"$x$")
it = nt - 1
ax02.plot(x, rho[it, :], "-b", x, rhotrue[it, :], ":k")
ax12.plot(x, u[it, :], "-b", x, utrue[it, :], ":k")
ax22.plot(x, p[it, :], "-b", x, ptrue[it, :], ":k")
ax02.set_title(r"$t =$ " + repr(t[it]).zfill(3))
ax22.set_xlabel(r"$x$")
plt.savefig("sod_frames.png")

plt.show()

##### end of script #####
