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

import gc

import numpy as np
import pytest

pytest.importorskip("mpi4py")

import sundials4py.core as sun  # noqa: E402

if not hasattr(sun, "N_VMake_MPIPlusX"):
    pytest.skip("sundials4py was built without MPI support", allow_module_level=True)

from mpi4py import MPI  # noqa: E402


@pytest.fixture
def comm():
    return MPI.COMM_WORLD


@pytest.fixture
def sunctx(comm):
    status, sunctx = sun.SUNContext_Create(comm)
    assert status == sun.SUN_SUCCESS
    return sunctx


def make_vector(comm, sunctx):
    """Returns (y, local, n_local, offset); local lengths differ by rank."""
    n_local = 3 + comm.Get_rank()
    offset = comm.exscan(n_local) or 0
    local = sun.N_VNew_Serial(n_local, sunctx)
    y = sun.N_VMake_MPIPlusX(comm, local, sunctx)
    # fill with the 1-based global index
    arr = sun.N_VGetNumpyArray(y)
    arr[:] = offset + np.arange(1, n_local + 1)
    return y, local, n_local, offset


def check_same_on_all_ranks(comm, value):
    values = comm.allgather(value)
    assert all(v == values[0] for v in values)


def test_global_length(comm, sunctx):
    y, local, n_local, _ = make_vector(comm, sunctx)
    assert sun.N_VGetLength(y) == comm.allreduce(n_local)
    assert sun.N_VGetLocalLength(y) == n_local
    assert sun.N_VGetLocalLength_MPIPlusX(y) == n_local


def test_reductions(comm, sunctx):
    y, local, n_local, _ = make_vector(comm, sunctx)
    n = comm.allreduce(n_local)
    ones = sun.N_VClone(y)
    sun.N_VConst(1.0, ones)

    dot = sun.N_VDotProd(y, y)
    assert dot == pytest.approx(n * (n + 1) * (2 * n + 1) / 6)
    check_same_on_all_ranks(comm, dot)

    l1 = sun.N_VL1Norm(y)
    assert l1 == pytest.approx(n * (n + 1) / 2)
    check_same_on_all_ranks(comm, l1)

    mx = sun.N_VMaxNorm(y)
    assert mx == pytest.approx(n)
    check_same_on_all_ranks(comm, mx)

    mn = sun.N_VMin(y)
    assert mn == pytest.approx(1.0)
    check_same_on_all_ranks(comm, mn)

    wrms = sun.N_VWrmsNorm(y, ones)
    assert wrms == pytest.approx(np.sqrt(n * (n + 1) * (2 * n + 1) / 6 / n))
    check_same_on_all_ranks(comm, wrms)

    denom = sun.N_VClone(y)
    sun.N_VConst(2.0, denom)
    minq = sun.N_VMinQuotient(y, denom)
    assert minq == pytest.approx(0.5)
    check_same_on_all_ranks(comm, minq)


def test_array_view_is_local(comm, sunctx):
    y, local, n_local, offset = make_vector(comm, sunctx)

    arr = sun.N_VGetNumpyArray(y)
    assert len(arr) == n_local
    assert len(sun.N_VGetArrayPointer(y)) == n_local

    # writes through the view are visible in the local vector
    arr[:] = -1.0
    np.testing.assert_allclose(sun.N_VGetNumpyArray(local), -1.0)


def test_local_vector_survives_dropping_reference(comm, sunctx):
    y, local, n_local, offset = make_vector(comm, sunctx)
    del local
    gc.collect()
    n = comm.allreduce(n_local)
    assert sun.N_VMaxNorm(y) == pytest.approx(n)
    assert len(sun.N_VGetNumpyArray(y)) == n_local


def test_local_vector_wrapper_keeps_parent_alive(comm, sunctx):
    y, local, n_local, offset = make_vector(comm, sunctx)
    del local
    gc.collect()

    loc = sun.N_VGetLocalVector_MPIPlusX(y)
    del y
    gc.collect()

    # the local vector (and the parent that owns it) must still be valid
    assert sun.N_VGetLength(loc) == n_local
    np.testing.assert_allclose(sun.N_VGetNumpyArray(loc), offset + np.arange(1, n_local + 1))


def test_get_local_vector_requires_mpiplusx(sunctx):
    v = sun.N_VNew_Serial(2, sunctx)
    with pytest.raises(Exception):
        sun.N_VGetLocalVector_MPIPlusX(v)


def test_clone_and_destroy(comm, sunctx):
    y, local, n_local, offset = make_vector(comm, sunctx)
    z = sun.N_VClone(y)
    sun.N_VConst(2.0, z)
    sun.N_VLinearSum(1.0, y, 1.0, z, z)
    n = comm.allreduce(n_local)
    assert sun.N_VMaxNorm(z) == pytest.approx(n + 2.0)
    del z
    gc.collect()
    assert sun.N_VMaxNorm(y) == pytest.approx(n)


def test_null_comm_makes_no_communication(sunctx):
    local = sun.N_VNew_Serial(4, sunctx)
    y = sun.N_VMake_MPIPlusX(None, local, sunctx)
    assert sun.N_VGetLength(y) == 4
    assert sun.N_VGetCommunicator(y) == sun.SUN_COMM_NULL


def test_comm_null_makes_no_communication(sunctx):
    # MPI.COMM_NULL converts to the null communicator, like None
    local = sun.N_VNew_Serial(4, sunctx)
    y = sun.N_VMake_MPIPlusX(MPI.COMM_NULL, local, sunctx)
    assert sun.N_VGetLength(y) == 4
    assert sun.N_VGetCommunicator(y) == sun.SUN_COMM_NULL


def test_sub_communicator(comm):
    sub = comm.Split(comm.Get_rank() % 2)
    status, subctx = sun.SUNContext_Create(sub)
    assert status == sun.SUN_SUCCESS
    n_local = 2 + comm.Get_rank()
    local = sun.N_VNew_Serial(n_local, subctx)
    y = sun.N_VMake_MPIPlusX(sub, local, subctx)

    assert sun.N_VGetLength(y) == sub.allreduce(n_local)
    # the vector holds a duplicate of sub, so it is congruent but not identical
    assert sun.N_VGetCommunicator(y).Compare(sub) == MPI.CONGRUENT
    sun.N_VConst(1.0, y)
    assert sun.N_VL1Norm(y) == pytest.approx(sub.allreduce(n_local))

    # destroy everything that uses sub before freeing it
    del y, local, subctx
    gc.collect()
    sub.Free()
