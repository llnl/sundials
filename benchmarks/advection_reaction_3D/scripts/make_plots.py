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
# matplotlib-based plotting script for the advection_reaction_3D benchmark codes
# ------------------------------------------------------------------------------

# imports
from os.path import exists
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------


# utility functions
def parallel_coords(rank):
    if rank == 0:
        return [0, 0, 0]
    if rank == 1:
        return [0, 0, 1]
    if rank == 2:
        return [0, 1, 0]
    if rank == 3:
        return [0, 1, 1]
    if rank == 4:
        return [1, 0, 0]
    if rank == 5:
        return [1, 0, 1]
    if rank == 6:
        return [1, 1, 0]
    if rank == 7:
        return [1, 1, 1]


def xslice(u, it, ix):
    return u[it, ix, :, :]


def yslice(u, it, iy):
    return u[it, :, iy, :]


def zslice(u, it, iz):
    return u[it, :, :, iz]


def xproj(u, it):
    return np.average(u[it, :, :, :], axis=0)


def yproj(u, it):
    return np.average(u[it, :, :, :], axis=1)


def zproj(u, it):
    return np.average(u[it, :, :, :], axis=2)


def myplot(axis, X, Y, Z, xlabel="none", ylabel="none"):
    frame = axis.contourf(X, Y, Z)
    plt.colorbar(frame, ax=axis)
    if xlabel != "none":
        axis.set_xlabel(xlabel)
    if ylabel != "none":
        axis.set_ylabel(ylabel)


# read time mesh
times = np.loadtxt("t.000000.txt")
nt = times.size

# read spatial mesh
mesh = np.loadtxt("mesh.txt", dtype=float)
x = mesh[0, :]
y = mesh[1, :]
z = mesh[2, :]
nx = x.size
ny = y.size
nz = z.size

# ensure that the run used exactly 1 or 8 MPI ranks
for i in range(9):
    if exists("u.00000" + str(i) + ".txt") and not exists("u.00000" + str(i + 1) + ".txt"):
        nprocs = i + 1
if (nprocs != 1) and (nprocs != 8):
    print("make_plots.py error: run must have used either 1 or 8 MPI ranks")
    exit()

# load data for run
if nprocs == 1:
    u = np.zeros((nt, nx, ny, nz), dtype=float)
    v = np.zeros((nt, nx, ny, nz), dtype=float)
    w = np.zeros((nt, nx, ny, nz), dtype=float)
    udata = np.loadtxt("u.000000.txt")
    vdata = np.loadtxt("v.000000.txt")
    wdata = np.loadtxt("w.000000.txt")
    if nt != udata.shape[0]:
        print("make_plots.py error: mesh and data have incompatible sizes")
        exit()
    if nx * ny * nz != udata.shape[1]:
        print("make_plots.py error: mesh and data have incompatible sizes")
        exit()
    for it in range(nt):
        u[it, :, :, :] = np.reshape(udata[it, :], (nx, ny, nz), order="C")
        v[it, :, :, :] = np.reshape(vdata[it, :], (nx, ny, nz), order="C")
        w[it, :, :, :] = np.reshape(wdata[it, :], (nx, ny, nz), order="C")
else:
    u = np.zeros((nt, nx, ny, nz), dtype=float)
    v = np.zeros((nt, nx, ny, nz), dtype=float)
    w = np.zeros((nt, nx, ny, nz), dtype=float)
    nxl = nx // 2
    nyl = ny // 2
    nzl = nz // 2
    for ip in range(8):
        udata = np.loadtxt("u.00000" + str(ip) + ".txt")
        vdata = np.loadtxt("v.00000" + str(ip) + ".txt")
        wdata = np.loadtxt("w.00000" + str(ip) + ".txt")
        if nt != udata.shape[0]:
            print("make_plots.py error: mesh and data have incompatible sizes")
            exit()
        if nxl * nyl * nzl != udata.shape[1]:
            print("make_plots.py error: mesh and data have incompatible sizes")
            exit()
        coords = parallel_coords(ip)
        ilo = coords[0] * nxl
        ihi = (coords[0] + 1) * nxl
        jlo = coords[1] * nyl
        jhi = (coords[1] + 1) * nyl
        klo = coords[2] * nzl
        khi = (coords[2] + 1) * nzl
        for it in range(nt):
            u[it, ilo:ihi, jlo:jhi, klo:khi] = np.reshape(udata[it, :], (nxl, nyl, nzl), order="C")
            v[it, ilo:ihi, jlo:jhi, klo:khi] = np.reshape(vdata[it, :], (nxl, nyl, nzl), order="C")
            w[it, ilo:ihi, jlo:jhi, klo:khi] = np.reshape(wdata[it, :], (nxl, nyl, nzl), order="C")


# set meshgrid objects
xy0, xy1 = np.meshgrid(x, y)
yz0, yz1 = np.meshgrid(y, z)
xz0, xz1 = np.meshgrid(x, z)

# generate plots
sliceidx = 25
tslice = [0, 5, 10]
figsize = (9, 7)

#    xy slices at various times
plt.figure(1)
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(
    3, 3, sharex=True, sharey=True, figsize=figsize
)
myplot(ax1, xy0, xy1, zslice(u, tslice[0], sliceidx), ylabel="u")
myplot(ax2, xy0, xy1, zslice(u, tslice[1], sliceidx))
myplot(ax3, xy0, xy1, zslice(u, tslice[2], sliceidx))
myplot(ax4, xy0, xy1, zslice(v, tslice[0], sliceidx), ylabel="v")
myplot(ax5, xy0, xy1, zslice(v, tslice[1], sliceidx))
myplot(ax6, xy0, xy1, zslice(v, tslice[2], sliceidx))
myplot(ax7, xy0, xy1, zslice(w, tslice[0], sliceidx), ylabel="w", xlabel="t = " + str(times[0]))
myplot(ax8, xy0, xy1, zslice(w, tslice[1], sliceidx), xlabel="t = " + str(times[1]))
myplot(ax9, xy0, xy1, zslice(w, tslice[2], sliceidx), xlabel="t = " + str(times[2]))
plt.savefig("xy-slices.png")

#    yz slices at various times
plt.figure(2)
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(
    3, 3, sharex=True, sharey=True, figsize=figsize
)
myplot(ax1, yz0, yz1, xslice(u, tslice[0], sliceidx), ylabel="u")
myplot(ax2, yz0, yz1, xslice(u, tslice[1], sliceidx))
myplot(ax3, yz0, yz1, xslice(u, tslice[2], sliceidx))
myplot(ax4, yz0, yz1, xslice(v, tslice[0], sliceidx), ylabel="v")
myplot(ax5, yz0, yz1, xslice(v, tslice[1], sliceidx))
myplot(ax6, yz0, yz1, xslice(v, tslice[2], sliceidx))
myplot(ax7, yz0, yz1, xslice(w, tslice[0], sliceidx), ylabel="w", xlabel="t = " + str(times[0]))
myplot(ax8, yz0, yz1, xslice(w, tslice[1], sliceidx), xlabel="t = " + str(times[1]))
myplot(ax9, yz0, yz1, xslice(w, tslice[2], sliceidx), xlabel="t = " + str(times[2]))
plt.savefig("yz-slices.png")

#    xz slices at various times
plt.figure(3)
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(
    3, 3, sharex=True, sharey=True, figsize=figsize
)
myplot(ax1, xz0, xz1, yslice(u, tslice[0], sliceidx), ylabel="u")
myplot(ax2, xz0, xz1, yslice(u, tslice[1], sliceidx))
myplot(ax3, xz0, xz1, yslice(u, tslice[2], sliceidx))
myplot(ax4, xz0, xz1, yslice(v, tslice[0], sliceidx), ylabel="v")
myplot(ax5, xz0, xz1, yslice(v, tslice[1], sliceidx))
myplot(ax6, xz0, xz1, yslice(v, tslice[2], sliceidx))
myplot(ax7, xz0, xz1, yslice(w, tslice[0], sliceidx), ylabel="w", xlabel="t = " + str(times[0]))
myplot(ax8, xz0, xz1, yslice(w, tslice[1], sliceidx), xlabel="t = " + str(times[1]))
myplot(ax9, xz0, xz1, yslice(w, tslice[2], sliceidx), xlabel="t = " + str(times[2]))
plt.savefig("xz-slices.png")

#    xy projection at various times
plt.figure(4)
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(
    3, 3, sharex=True, sharey=True, figsize=figsize
)
myplot(ax1, xy0, xy1, zproj(u, tslice[0]), ylabel="u")
myplot(ax2, xy0, xy1, zproj(u, tslice[1]))
myplot(ax3, xy0, xy1, zproj(u, tslice[2]))
myplot(ax4, xy0, xy1, zproj(v, tslice[0]), ylabel="v")
myplot(ax5, xy0, xy1, zproj(v, tslice[1]))
myplot(ax6, xy0, xy1, zproj(v, tslice[2]))
myplot(ax7, xy0, xy1, zproj(w, tslice[0]), ylabel="w", xlabel="t = " + str(times[0]))
myplot(ax8, xy0, xy1, zproj(w, tslice[1]), xlabel="t = " + str(times[1]))
myplot(ax9, xy0, xy1, zproj(w, tslice[2]), xlabel="t = " + str(times[2]))
plt.savefig("xy-projections.png")

#    yz projection at various times
fig = plt.figure(5)
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(
    3, 3, sharex=True, sharey=True, figsize=figsize
)
myplot(ax1, yz0, yz1, xproj(u, tslice[0]), ylabel="u")
myplot(ax2, yz0, yz1, xproj(u, tslice[1]))
myplot(ax3, yz0, yz1, xproj(u, tslice[2]))
myplot(ax4, yz0, yz1, xproj(v, tslice[0]), ylabel="v")
myplot(ax5, yz0, yz1, xproj(v, tslice[1]))
myplot(ax6, yz0, yz1, xproj(v, tslice[2]))
myplot(ax7, yz0, yz1, xproj(w, tslice[0]), ylabel="w", xlabel="t = " + str(times[0]))
myplot(ax8, yz0, yz1, xproj(w, tslice[1]), xlabel="t = " + str(times[1]))
myplot(ax9, yz0, yz1, xproj(w, tslice[2]), xlabel="t = " + str(times[2]))
plt.savefig("yz-projections.png")

#    xz projection at various times
fig = plt.figure(6)
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(
    3, 3, sharex=True, sharey=True, figsize=figsize
)
myplot(ax1, xz0, xz1, yproj(u, tslice[0]), ylabel="u")
myplot(ax2, xz0, xz1, yproj(u, tslice[1]))
myplot(ax3, xz0, xz1, yproj(u, tslice[2]))
myplot(ax4, xz0, xz1, yproj(v, tslice[0]), ylabel="v")
myplot(ax5, xz0, xz1, yproj(v, tslice[1]))
myplot(ax6, xz0, xz1, yproj(v, tslice[2]))
myplot(ax7, xz0, xz1, yproj(w, tslice[0]), ylabel="w", xlabel="t = " + str(times[0]))
myplot(ax8, xz0, xz1, yproj(w, tslice[1]), xlabel="t = " + str(times[1]))
myplot(ax9, xz0, xz1, yproj(w, tslice[2]), xlabel="t = " + str(times[2]))
plt.savefig("xz-projections.png")

# plt.show()
plt.close()

##### end of script #####
