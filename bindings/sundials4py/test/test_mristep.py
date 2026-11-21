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


import pytest
import weakref
import numpy as np
from numpy.testing import assert_allclose
from sundials4py.core import *
from sundials4py.arkode import *
from problems import AnalyticMultiscaleODE
from fixtures import sunctx


class RecordingMRIController(CustomSUNMRIController):
    """MRI controller that leaves the proposed step and tolerance unchanged."""

    def __init__(self, sunctx):
        self.calls = 0
        super().__init__(sunctx)

    def estimate_step_tol(self, H, tolfac, P, DSM, dsm):
        self.calls += 1
        return SUN_SUCCESS, 1e-3, tolfac


def test_multirate(sunctx):
    ode_problem = AnalyticMultiscaleODE()

    t0, tf = AnalyticMultiscaleODE.T0, 0.01

    def fslow(t, y, ydot, _):
        return ode_problem.f_linear(t, y, ydot, None)

    def ffast(t, y, ydot, _):
        # test MRIStepInnerStepper_GetForcingData
        status, tshift, tscale, forcing, nforcing = MRIStepInnerStepper_GetForcingData(
            ode_problem.inner_stepper
        )
        assert status == ARK_SUCCESS
        assert len(forcing) == nforcing

        return ode_problem.f_nonlinear(t, y, ydot, None)

    y = N_VNew_Serial(1, sunctx)
    y0 = N_VClone(y)

    ode_problem.set_init_cond(y)
    ode_problem.set_init_cond(y0)

    # create fast integrator
    inner_ark = ERKStepCreate(ffast, t0, y, sunctx)
    status = ARKodeSetFixedStep(inner_ark.get(), 5e-3)
    assert status == ARK_SUCCESS

    status, inner_stepper = ARKodeCreateMRIStepInnerStepper(inner_ark.get())
    assert status == ARK_SUCCESS

    # store inner_stepper in ode_problem so we can access it in ffast
    ode_problem.inner_stepper = inner_stepper

    # create slow integrator
    ark = MRIStepCreate(fslow, None, t0, y, inner_stepper, sunctx)
    status = ARKodeSetFixedStep(ark.get(), 1e-3)
    assert status == ARK_SUCCESS

    tout = tf
    status, tret = ARKodeEvolve(ark.get(), tout, y, ARK_NORMAL)
    assert status == ARK_SUCCESS

    sol = N_VClone(y)
    ode_problem.solution(y, sol, tret)
    # we use a fixed atol here since we use a fixed step size
    assert_allclose(N_VGetArrayPointer(sol), N_VGetArrayPointer(y), atol=1e-2)


def test_custom_mri_controller_through_mristep(sunctx):
    """A package-driven MRIStep solve invokes estimate_step_tol."""
    ode_problem = AnalyticMultiscaleODE()
    t0, tf = AnalyticMultiscaleODE.T0, 0.01

    y = N_VNew_Serial(1, sunctx)
    y0 = N_VClone(y)
    ode_problem.set_init_cond(y)
    ode_problem.set_init_cond(y0)

    def fslow(t, state, derivative, _):
        return ode_problem.f_linear(t, state, derivative, None)

    def ffast(t, state, derivative, _):
        return ode_problem.f_nonlinear(t, state, derivative, None)

    inner_ark = ERKStepCreate(ffast, t0, y, sunctx)
    assert ARKodeSetFixedStep(inner_ark.get(), 5e-3) == ARK_SUCCESS
    assert ARKodeSetAccumulatedErrorType(inner_ark.get(), ARK_ACCUMERROR_SUM) == ARK_SUCCESS
    status, inner_stepper = ARKodeCreateMRIStepInnerStepper(inner_ark.get())
    assert status == ARK_SUCCESS

    ark = MRIStepCreate(fslow, None, t0, y, inner_stepper, sunctx)
    assert ARKodeSStolerances(ark.get(), 1e-4, 1e-8) == ARK_SUCCESS
    assert ARKodeSetInitStep(ark.get(), 1e-3) == ARK_SUCCESS
    controller = RecordingMRIController(sunctx)
    assert ARKodeSetAdaptController(ark.get(), controller) == ARK_SUCCESS

    status, tret = ARKodeEvolve(ark.get(), tf, y, ARK_NORMAL)
    assert status == ARK_SUCCESS
    assert tret == pytest.approx(tf)

    expected = N_VClone(y)
    ode_problem.solution(y0, expected, tret)
    assert_allclose(N_VGetArrayPointer(y), N_VGetArrayPointer(expected), atol=1e-2)
    assert controller.calls > 0
