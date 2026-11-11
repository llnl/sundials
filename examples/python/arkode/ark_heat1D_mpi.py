#!/usr/bin/env python3
# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------
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
# -----------------------------------------------------------------
# This example solves the 1D heat equation
#   u_t = k u_xx,  u(0,t) = u(1,t) = 0,
# distributed over MPI ranks with an MPIPlusX N_Vector over a serial local
# vector, mpi4py halo exchange on NumPy views, and ARKStep with SPGMR
# (classical Gram-Schmidt, which needs fewer global reductions).
#
# Requires a sundials4py build with SUNDIALS_ENABLE_MPI=ON and mpi4py.
# Run with:
#   mpiexec -n 4 python -m mpi4py ark_heat1D_mpi.py
# -----------------------------------------------------------------

import sys

import numpy as np

try:
    from mpi4py import MPI
except ImportError:
    MPI = None

from sundials4py.arkode import *
from sundials4py.core import *

NX = 201  # global number of interior grid points
K = 0.5
TF = 0.1


def mpi_available():
    return MPI is not None and "N_VMake_MPIPlusX" in globals()


class Heat1D:
    def __init__(self, comm, n_local, offset):
        self.comm = comm
        self.n_local = n_local
        self.offset = offset
        self.dx = 1.0 / (NX + 1)
        rank, size = comm.Get_rank(), comm.Get_size()
        self.left = rank - 1 if rank > 0 else MPI.PROC_NULL
        self.right = rank + 1 if rank < size - 1 else MPI.PROC_NULL

    def x(self):
        return (self.offset + np.arange(1, self.n_local + 1)) * self.dx

    def f(self, t, yvec, fvec, user_data):
        u = N_VGetNumpyArray(yvec)  # local slice of the distributed vector
        f = N_VGetNumpyArray(fvec)
        ext = np.zeros(self.n_local + 2, dtype=u.dtype)
        ext[1:-1] = u
        from_left = np.zeros(1, dtype=u.dtype)
        from_right = np.zeros(1, dtype=u.dtype)
        last = np.ascontiguousarray(u[-1:])
        first = np.ascontiguousarray(u[:1])
        self.comm.Sendrecv(last, dest=self.right, recvbuf=from_left, source=self.left)
        self.comm.Sendrecv(first, dest=self.left, recvbuf=from_right, source=self.right)
        ext[0] = from_left[0]
        ext[-1] = from_right[0]
        f[:] = K * (ext[:-2] - 2.0 * ext[1:-1] + ext[2:]) / self.dx**2
        return 0


def main():
    comm = MPI.COMM_WORLD
    rank, size = comm.Get_rank(), comm.Get_size()
    if size > NX:
        raise RuntimeError(
            f"this example needs at most {NX} MPI ranks (one grid point per "
            f"rank), but was launched with {size}"
        )

    status, sunctx = SUNContext_Create(comm)
    assert status == SUN_SUCCESS

    base, rem = divmod(NX, size)
    n_local = base + (1 if rank < rem else 0)
    offset = rank * base + min(rank, rem)

    local = N_VNew_Serial(n_local, sunctx)
    y = N_VMake_MPIPlusX(comm, local, sunctx)
    assert N_VGetLength(y) == NX

    problem = Heat1D(comm, n_local, offset)
    x = problem.x()
    N_VGetNumpyArray(y)[:] = np.sin(np.pi * x)

    ark = ARKStepCreate(None, problem.f, 0.0, y, sunctx)
    assert ARKodeSStolerances(ark.get(), 1.0e-6, 1.0e-10) == ARK_SUCCESS
    LS = SUNLinSol_SPGMR(y, 0, 0, sunctx)
    assert SUNLinSol_SPGMRSetGSType(LS, SUN_CLASSICAL_GS) == SUN_SUCCESS
    assert ARKodeSetLinearSolver(ark.get(), LS, None) == ARK_SUCCESS

    status, _ = ARKodeEvolve(ark.get(), TF, y, ARK_NORMAL)
    assert status == ARK_SUCCESS

    # exact solution of the first mode: exp(-k pi^2 t) sin(pi x)
    exact = np.exp(-K * np.pi**2 * TF) * np.sin(np.pi * x)
    err = N_VGetNumpyArray(y) - exact
    max_err = comm.allreduce(np.max(np.abs(err)), op=MPI.MAX)
    status, nst = ARKodeGetNumSteps(ark.get())
    assert status == ARK_SUCCESS
    status, nli = ARKodeGetNumLinIters(ark.get())
    assert status == ARK_SUCCESS

    if rank == 0:
        print(f"1D heat equation on {size} rank(s), {NX} points")
        print(f"   steps = {nst}, linear iterations = {nli}")
        print(f"   max error = {max_err:.3e}")
    assert max_err < 1.0e-3


def test_ark_heat1d_mpi():
    if not mpi_available():
        return
    main()


if __name__ == "__main__":
    if not mpi_available():
        print("This example needs mpi4py and an MPI-enabled sundials4py build")
        sys.exit(1)
    main()
