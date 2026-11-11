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

import numpy as np
import pytest

pytest.importorskip("mpi4py")

import sundials4py.core as sun  # noqa: E402
from sundials4py.arkode import *  # noqa: E402,F403
from sundials4py.cvodes import *  # noqa: E402,F403

if not hasattr(sun, "N_VMake_MPIPlusX"):
    pytest.skip("sundials4py was built without MPI support", allow_module_level=True)

from mpi4py import MPI  # noqa: E402

NX = 40  # global number of grid points
DIFF = 0.5
TF = 0.05
RTOL = 1.0e-8
ATOL = 1.0e-10


def initial_condition(x):
    return np.sin(np.pi * x) + 0.5 * np.sin(3 * np.pi * x)


class Diffusion1D:
    """u_t = DIFF u_xx on (0, 1) with u = 0 at both ends. If comm is None the
    whole grid lives in this process and no halo exchange happens."""

    def __init__(self, comm, n_local, offset):
        self.comm = comm
        self.n_local = n_local
        self.offset = offset
        self.dx = 1.0 / (NX + 1)
        rank = comm.Get_rank() if comm else 0
        size = comm.Get_size() if comm else 1
        self.left = rank - 1 if rank > 0 else MPI.PROC_NULL
        self.right = rank + 1 if rank < size - 1 else MPI.PROC_NULL

    def x(self):
        return (self.offset + np.arange(1, self.n_local + 1)) * self.dx

    def f(self, t, yvec, fvec, _):
        u = sun.N_VGetNumpyArray(yvec)
        f = sun.N_VGetNumpyArray(fvec)
        ext = np.zeros(self.n_local + 2, dtype=u.dtype)
        ext[1:-1] = u
        if self.comm is not None:
            # halo exchange of the boundary values (zeros at the domain ends)
            lo = np.ascontiguousarray(u[:1])
            hi = np.ascontiguousarray(u[-1:])
            from_left = np.zeros(1, dtype=u.dtype)
            from_right = np.zeros(1, dtype=u.dtype)
            self.comm.Sendrecv(hi, dest=self.right, recvbuf=from_left, source=self.left)
            self.comm.Sendrecv(lo, dest=self.left, recvbuf=from_right, source=self.right)
            ext[0] = from_left[0]
            ext[-1] = from_right[0]
        f[:] = DIFF * (ext[:-2] - 2.0 * ext[1:-1] + ext[2:]) / self.dx**2
        return 0


def decompose(comm):
    size = comm.Get_size()
    rank = comm.Get_rank()
    base, rem = divmod(NX, size)
    n_local = base + (1 if rank < rem else 0)
    offset = rank * base + min(rank, rem)
    return n_local, offset


def run(integrator, y, problem, sunctx):
    """Integrate with ARKStep+SPGMR or CVODES+SPGMR; returns counts."""
    ls = sun.SUNLinSol_SPGMR(y, 0, 0, sunctx)
    assert sun.SUNLinSol_SPGMRSetGSType(ls, sun.SUN_CLASSICAL_GS) == sun.SUN_SUCCESS
    if integrator == "arkstep":
        ark = ARKStepCreate(None, problem.f, 0.0, y, sunctx)
        assert ARKodeSStolerances(ark.get(), RTOL, ATOL) == ARK_SUCCESS
        assert ARKodeSetLinearSolver(ark.get(), ls, None) == ARK_SUCCESS
        status, _ = ARKodeEvolve(ark.get(), TF, y, ARK_NORMAL)
        assert status == ARK_SUCCESS
        _, nst = ARKodeGetNumSteps(ark.get())
        _, nli = ARKodeGetNumLinIters(ark.get())
    else:
        cv = CVodeCreate(CV_BDF, sunctx)
        assert CVodeInit(cv.get(), problem.f, 0.0, y) == CV_SUCCESS
        assert CVodeSStolerances(cv.get(), RTOL, ATOL) == CV_SUCCESS
        assert CVodeSetLinearSolver(cv.get(), ls, None) == CV_SUCCESS
        status, _ = CVode(cv.get(), TF, y, CV_NORMAL)
        assert status == CV_SUCCESS
        _, nst = CVodeGetNumSteps(cv.get())
        _, nli = CVodeGetNumLinIters(cv.get())
    return nst, nli


@pytest.mark.parametrize("integrator", ["arkstep", "cvodes"])
def test_diffusion_matches_single_rank_reference(integrator):
    comm = MPI.COMM_WORLD
    status, sunctx = sun.SUNContext_Create(comm)
    assert status == sun.SUN_SUCCESS

    # Distributed solve with MPIPlusX over a serial local vector
    n_local, offset = decompose(comm)
    local = sun.N_VNew_Serial(n_local, sunctx)
    y = sun.N_VMake_MPIPlusX(comm, local, sunctx)
    problem = Diffusion1D(comm, n_local, offset)
    sun.N_VGetNumpyArray(y)[:] = initial_condition(problem.x())
    counts = run(integrator, y, problem, sunctx)
    u_dist = np.concatenate(comm.allgather(np.array(sun.N_VGetNumpyArray(y))))

    # Reference: the whole problem, redundantly on every rank, serial vector
    status, ref_ctx = sun.SUNContext_Create(None)
    yref = sun.N_VNew_Serial(NX, ref_ctx)
    ref_problem = Diffusion1D(None, NX, 0)
    sun.N_VGetNumpyArray(yref)[:] = initial_condition(ref_problem.x())
    run(integrator, yref, ref_problem, ref_ctx)
    u_ref = np.array(sun.N_VGetNumpyArray(yref))

    np.testing.assert_allclose(u_dist, u_ref, rtol=1.0e-5, atol=1.0e-6)

    # Every rank must take the same steps and linear iterations
    all_counts = comm.allgather(counts)
    assert all(c == all_counts[0] for c in all_counts)
