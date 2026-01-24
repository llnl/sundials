# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2025, Lawrence Livermore National Security,
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
# Unit/smoke tests for SUNStepper module
# -----------------------------------------------------------------

import pytest
from fixtures import *
from sundials4py.core import *


def make_stepper(sunctx):
    # Create am empty stepper
    status, s = SUNStepper_Create(sunctx)
    assert status == SUN_SUCCESS
    return s


def test_create_stepper(sunctx):
    s = make_stepper(sunctx)
    assert s is not None


def test_stepper_evolve(sunctx, nvec):
    s = make_stepper(sunctx)
    vret = nvec
    err, tret = SUNStepper_Evolve(s, 1.0, vret)
    assert err == SUN_ERR_NOT_IMPLEMENTED
    assert isinstance(tret, float)


def test_stepper_one_step(sunctx, nvec):
    s = make_stepper(sunctx)
    vret = nvec
    err, tret = SUNStepper_OneStep(s, 1.0, vret)
    assert err == SUN_ERR_NOT_IMPLEMENTED
    assert isinstance(tret, float)


def test_stepper_reset(sunctx, nvec):
    s = make_stepper(sunctx)
    err = SUNStepper_Reset(s, 0.0, nvec)
    assert err == SUN_ERR_NOT_IMPLEMENTED


def test_stepper_set_evolve_fn(sunctx, nvec):
    s = make_stepper(sunctx)
    called = {"flag": False}

    def evolve_fn(stepper, tout, vret, tret):
        called["flag"] = True
        return 0

    err = SUNStepper_SetEvolveFn(s, evolve_fn)
    assert err == SUN_SUCCESS
    # Call evolve to trigger callback
    SUNStepper_Evolve(s, 1.0, nvec)
    assert called["flag"]


def test_stepper_set_one_step_fn(sunctx, nvec):
    s = make_stepper(sunctx)
    called = {"flag": False}

    def one_step_fn(stepper, tout, vret, tret):
        called["flag"] = True
        return 0

    err = SUNStepper_SetOneStepFn(s, one_step_fn)
    assert err == SUN_SUCCESS
    SUNStepper_OneStep(s, 1.0, nvec)
    assert called["flag"]


def test_stepper_set_full_rhs_fn(sunctx, nvec):
    s = make_stepper(sunctx)
    called = {"flag": False}

    def full_rhs_fn(stepper, t, v, f, mode):
        called["flag"] = True
        return 0

    err = SUNStepper_SetFullRhsFn(s, full_rhs_fn)
    assert err == SUN_SUCCESS
    # Call with dummy args
    SUNStepper_FullRhs(s, 0.0, nvec, nvec, 0)
    assert called["flag"]


def test_stepper_set_reinit_fn(sunctx, nvec):
    s = make_stepper(sunctx)
    called = {"flag": False}

    def reinit_fn(stepper, t, y):
        called["flag"] = True
        return 0

    err = SUNStepper_SetReInitFn(s, reinit_fn)
    assert err == SUN_SUCCESS
    SUNStepper_ReInit(s, 0.0, nvec)
    assert called["flag"]


def test_stepper_set_reset_fn(sunctx, nvec):
    s = make_stepper(sunctx)
    called = {"flag": False}

    def reset_fn(stepper, t, y):
        called["flag"] = True
        return 0

    err = SUNStepper_SetResetFn(s, reset_fn)
    assert err == SUN_SUCCESS
    SUNStepper_Reset(s, 0.0, nvec)
    assert called["flag"]


def test_stepper_set_reset_ckpt_idx_fn(sunctx):
    s = make_stepper(sunctx)
    called = {"flag": False}

    def reset_ckpt_idx_fn(stepper, idx):
        called["flag"] = True
        return 0

    err = SUNStepper_SetResetCheckpointIndexFn(s, reset_ckpt_idx_fn)
    assert err == SUN_SUCCESS
    SUNStepper_ResetCheckpointIndex(s, 1)
    assert called["flag"]


def test_stepper_set_stop_time_fn(sunctx):
    s = make_stepper(sunctx)
    called = {"flag": False}

    def stop_time_fn(stepper, tstop):
        called["flag"] = True
        return 0

    err = SUNStepper_SetStopTimeFn(s, stop_time_fn)
    assert err == SUN_SUCCESS
    SUNStepper_SetStopTime(s, 2.0)
    assert called["flag"]


def test_stepper_set_step_direction_fn(sunctx):
    s = make_stepper(sunctx)
    called = {"flag": False}

    def step_direction_fn(stepper, direction):
        called["flag"] = True
        return 0

    err = SUNStepper_SetStepDirectionFn(s, step_direction_fn)
    assert err == SUN_SUCCESS
    SUNStepper_SetStepDirection(s, 1.0)
    assert called["flag"]


def test_stepper_set_forcing_fn(sunctx, nvec):
    s = make_stepper(sunctx)
    called = {"flag": False}

    def forcing_fn(stepper, tshift, tscale, forcing, nforcing):
        assert type(forcing) is list
        assert len(forcing) == nforcing
        called["flag"] = True
        return 0

    err = SUNStepper_SetForcingFn(s, forcing_fn)
    assert err == SUN_SUCCESS
    SUNStepper_SetForcing(s, 0.0, 1.0, [nvec], 1)
    assert called["flag"]


def test_stepper_set_get_num_steps_fn(sunctx):
    s = make_stepper(sunctx)
    called = {"flag": False}

    def get_num_steps_fn(stepper):
        called["flag"] = True
        nst = 1
        return 0, nst

    err = SUNStepper_SetGetNumStepsFn(s, get_num_steps_fn)
    assert err == SUN_SUCCESS
    status, nst = SUNStepper_GetNumSteps(s)
    assert called["flag"]
    assert isinstance(nst, int)
    assert nst == 1
