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


def test_context_accepts_comm_variants(comm):
    for c in (comm, MPI.COMM_NULL, sun.SUN_COMM_NULL, None):
        status, sunctx = sun.SUNContext_Create(c)
        assert status == sun.SUN_SUCCESS
        assert sunctx is not None


@pytest.mark.parametrize("bad_comm", [1, True, 2.5, "comm", object()])
def test_context_rejects_invalid_comm(bad_comm):
    with pytest.raises(TypeError):
        sun.SUNContext_Create(bad_comm)


def test_serial_vector_communicator_is_null(sunctx):
    v = sun.N_VNew_Serial(4, sunctx)
    assert sun.N_VGetCommunicator(v) == sun.SUN_COMM_NULL


def test_mpiplusx_communicator_is_congruent(comm, sunctx):
    local = sun.N_VNew_Serial(4, sunctx)
    y = sun.N_VMake_MPIPlusX(comm, local, sunctx)
    vcomm = sun.N_VGetCommunicator(y)
    assert vcomm.Compare(comm) == MPI.CONGRUENT


def test_logger_and_profiler_with_comm(comm):
    status, logger = sun.SUNLogger_Create(comm, 0)
    assert status == sun.SUN_SUCCESS
    assert logger is not None

    status, logger = sun.SUNLogger_CreateFromEnv(comm)
    assert status == sun.SUN_SUCCESS

    status, profiler = sun.SUNProfiler_Create(comm, "test")
    assert status == sun.SUN_SUCCESS
    assert profiler is not None


def test_context_keeps_comm_alive(comm):
    # The context stores the communicator handle without duplicating it, so
    # it must keep the Python object alive.
    c = comm.Dup()
    status, sunctx = sun.SUNContext_Create(c)
    assert status == sun.SUN_SUCCESS
    del c
    gc.collect()

    # Using the context must still be safe, including creating MPI vectors.
    local = sun.N_VNew_Serial(3, sunctx)
    y = sun.N_VMake_MPIPlusX(comm, local, sunctx)
    assert sun.N_VGetLength(y) == comm.allreduce(3)
    assert sun.SUNContext_GetLastError(sunctx) == sun.SUN_SUCCESS
