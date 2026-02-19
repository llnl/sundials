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
# Unit/smoke tests for SUNAdjointStepper module
# -----------------------------------------------------------------

import pytest
from fixtures import *
from sundials4py.core import *


def make_adjoint_stepper(sunctx, sunstepper, nvec):
    mem_helper = SUNMemoryHelper_Sys(sunctx)
    status, scheme = SUNAdjointCheckpointScheme_Create_Fixed(0, mem_helper, 1, 1, 0, sunctx)
    own_stepper = False
    nsteps = 1
    t0 = 0.0
    y0 = nvec
    status, adj_stepper = SUNAdjointStepper_Create(
        sunstepper, own_stepper, sunstepper, own_stepper, nsteps, t0, y0, scheme, sunctx
    )
    return adj_stepper, scheme, mem_helper


def test_create_adjoint_stepper(sunctx, nvec, sunstepper):
    adj_stepper, scheme, mem_helper = make_adjoint_stepper(sunctx, sunstepper, nvec)
    assert adj_stepper is not None


def test_adjointstepper_reinit(sunctx, nvec, sunstepper):
    adj_stepper, scheme, mem_helper = make_adjoint_stepper(sunctx, sunstepper, nvec)
    t0 = 0.0
    y0 = nvec
    tf = 1.0
    err = SUNAdjointStepper_ReInit(adj_stepper, t0, y0, tf, nvec)
    assert err == SUN_SUCCESS


def test_adjointstepper_evolve(sunctx, nvec, sunstepper):
    adj_stepper, scheme, mem_helper = make_adjoint_stepper(sunctx, sunstepper, nvec)
    tout = 1.0
    sens = nvec
    err, tret = SUNAdjointStepper_Evolve(adj_stepper, tout, sens)
    assert err == SUN_SUCCESS
    assert isinstance(tret, float)


def test_adjointstepper_onestep(sunctx, nvec, sunstepper):
    adj_stepper, scheme, mem_helper = make_adjoint_stepper(sunctx, sunstepper, nvec)
    tout = 1.0
    sens = nvec
    err, tret = SUNAdjointStepper_OneStep(adj_stepper, tout, sens)
    assert err == SUN_SUCCESS
    assert isinstance(tret, float)


def test_adjointstepper_recomputefwd(sunctx, nvec, sunstepper):
    adj_stepper, scheme, mem_helper = make_adjoint_stepper(sunctx, sunstepper, nvec)
    start_idx = 0
    t0 = 0.0
    y0 = nvec
    tf = 1.0
    err = SUNAdjointStepper_RecomputeFwd(adj_stepper, start_idx, t0, y0, tf)
    assert err == SUN_SUCCESS


def test_adjointstepper_getnumsteps(sunctx, nvec, sunstepper):
    adj_stepper, scheme, mem_helper = make_adjoint_stepper(sunctx, sunstepper, nvec)
    err, num_steps = SUNAdjointStepper_GetNumSteps(adj_stepper)
    assert err == SUN_SUCCESS
    assert isinstance(num_steps, int)


def test_adjointstepper_getnumrecompute(sunctx, nvec, sunstepper):
    adj_stepper, scheme, mem_helper = make_adjoint_stepper(sunctx, sunstepper, nvec)
    err, num_recompute = SUNAdjointStepper_GetNumRecompute(adj_stepper)
    assert err == SUN_SUCCESS
    assert isinstance(num_recompute, int)
