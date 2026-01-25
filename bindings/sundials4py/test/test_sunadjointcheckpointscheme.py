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
# Unit/smoke tests for SUNAdjointCheckpointScheme module
# -----------------------------------------------------------------

import pytest
from fixtures import *
from sundials4py.core import *


def make_fixed_scheme(sunctx):
    io_mode = SUNDATAIOMODE_INMEM
    mem_helper = SUNMemoryHelper_Sys(sunctx)
    interval = 1
    estimate = 1
    keep = 0
    status, scheme = SUNAdjointCheckpointScheme_Create_Fixed(
        io_mode, mem_helper, interval, estimate, keep, sunctx
    )
    # must return mem_helper or it will get cleaned up
    return status, scheme, mem_helper


def test_needs_saving(sunctx):
    scheme_status, scheme, mem_helper = make_fixed_scheme(sunctx)
    step_num = 0
    stage_num = 0
    t = 0.0
    status, result = SUNAdjointCheckpointScheme_NeedsSaving(scheme, step_num, stage_num, t)
    assert status == SUN_SUCCESS
    assert isinstance(result, int)


def test_insert_vector(sunctx, nvec):
    scheme_status, scheme, mem_helper = make_fixed_scheme(sunctx)
    step_num = 0
    stage_num = 0
    t = 0.0
    status = SUNAdjointCheckpointScheme_InsertVector(scheme, step_num, stage_num, t, nvec)
    assert status == SUN_SUCCESS


def test_load_vector(sunctx, nvec):
    scheme_status, scheme, mem_helper = make_fixed_scheme(sunctx)

    step_num = 0
    stage_num = 0
    t = 0.0
    status = SUNAdjointCheckpointScheme_InsertVector(scheme, step_num, stage_num, t, nvec)
    assert status == SUN_SUCCESS

    peek = False
    status, vec, tout = SUNAdjointCheckpointScheme_LoadVector(
        scheme, step_num, stage_num, peek, nvec
    )
    assert status == SUN_SUCCESS


def test_enable_dense(sunctx):
    scheme_status, scheme, mem_helper = make_fixed_scheme(sunctx)
    status = SUNAdjointCheckpointScheme_EnableDense(scheme, True)
    assert status == SUN_SUCCESS
