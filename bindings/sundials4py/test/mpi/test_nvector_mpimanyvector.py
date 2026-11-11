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


def test_make_with_serial_subvectors(comm, sunctx):
    a = sun.N_VNew_Serial(2, sunctx)
    b = sun.N_VNew_Serial(3, sunctx)
    sun.N_VConst(1.0, a)
    sun.N_VConst(2.0, b)
    mv = sun.N_VMake_MPIManyVector(comm, 2, [a, b], sunctx)

    assert sun.N_VGetNumSubvectors_MPIManyVector(mv) == 2
    assert sun.N_VGetLength(mv) == comm.allreduce(5)
    assert sun.N_VGetSubvectorLocalLength_MPIManyVector(mv, 1) == 3
    assert sun.N_VL1Norm(mv) == pytest.approx(comm.allreduce(2 * 1.0 + 3 * 2.0))

    sub = sun.N_VGetSubvector_MPIManyVector(mv, 1)
    np.testing.assert_allclose(sun.N_VGetNumpyArray(sub), 2.0)


def test_new_over_mpiplusx_subvectors(comm, sunctx):
    subs = []
    for n in (2, 4):
        local = sun.N_VNew_Serial(n, sunctx)
        subs.append(sun.N_VMake_MPIPlusX(comm, local, sunctx))
    mv = sun.N_VNew_MPIManyVector(2, subs, sunctx)
    assert sun.N_VGetCommunicator(mv).Compare(comm) == MPI.CONGRUENT
    assert sun.N_VGetLength(mv) == comm.allreduce(6)


def test_new_with_only_serial_subvectors_raises(sunctx):
    subs = [sun.N_VNew_Serial(2, sunctx), sun.N_VNew_Serial(2, sunctx)]
    with pytest.raises(Exception):
        sun.N_VNew_MPIManyVector(2, subs, sunctx)


def test_invalid_subvector_lists(comm, sunctx):
    a = sun.N_VNew_Serial(2, sunctx)
    with pytest.raises(Exception):
        sun.N_VMake_MPIManyVector(comm, 2, [a], sunctx)
    with pytest.raises(Exception):
        sun.N_VMake_MPIManyVector(comm, 0, [], sunctx)
    with pytest.raises(Exception):
        sun.N_VMake_MPIManyVector(comm, 2, [a, None], sunctx)


def test_subvector_index_out_of_range(comm, sunctx):
    a = sun.N_VNew_Serial(2, sunctx)
    b = sun.N_VNew_Serial(2, sunctx)
    mv = sun.N_VMake_MPIManyVector(comm, 2, [a, b], sunctx)
    with pytest.raises(IndexError):
        sun.N_VGetSubvector_MPIManyVector(mv, 2)
    with pytest.raises(IndexError):
        sun.N_VGetSubvector_MPIManyVector(mv, -1)


def test_array_access_raises(comm, sunctx):
    a = sun.N_VNew_Serial(2, sunctx)
    b = sun.N_VNew_Serial(2, sunctx)
    mv = sun.N_VMake_MPIManyVector(comm, 2, [a, b], sunctx)
    with pytest.raises(Exception, match="N_VGetSubvector_MPIManyVector"):
        sun.N_VGetNumpyArray(mv)


def test_mutating_caller_list_does_not_free_subvectors(comm, sunctx):
    subs = [sun.N_VNew_Serial(2, sunctx), sun.N_VNew_Serial(3, sunctx)]
    for s in subs:
        sun.N_VConst(1.0, s)
    mv = sun.N_VMake_MPIManyVector(comm, 2, subs, sunctx)
    subs.clear()
    gc.collect()
    assert sun.N_VL1Norm(mv) == pytest.approx(comm.allreduce(5.0))
    assert sun.N_VGetSubvectorLocalLength_MPIManyVector(mv, 1) == 3
