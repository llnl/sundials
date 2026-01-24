#!/usr/bin/env python3
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


import pytest
import numpy as np
from numpy.testing import assert_allclose
from fixtures import *
from sundials4py.core import *
from sundials4py.arkode import *
from problems import AnalyticMultiscaleODE


def test_forcingstep(sunctx):
    ode_problem = AnalyticMultiscaleODE()
    t0, tf = AnalyticMultiscaleODE.T0, 0.01

    y = N_VNew_Serial(1, sunctx)
    y0 = N_VClone(y)
    ode_problem.set_init_cond(y)
    ode_problem.set_init_cond(y0)

    linear_ark = ERKStepCreate(ode_problem.f_linear, t0, y, sunctx)
    status = ARKodeSetFixedStep(linear_ark.get(), 5e-3)
    assert status == ARK_SUCCESS

    nonlinear_ark = ARKStepCreate(ode_problem.f_nonlinear, None, t0, y, sunctx)
    status = ARKodeSetFixedStep(nonlinear_ark.get(), 1e-3)
    assert status == ARK_SUCCESS

    status, linear_stepper = ARKodeCreateSUNStepper(linear_ark.get())
    status, nonlinear_stepper = ARKodeCreateSUNStepper(nonlinear_ark.get())

    ark = ForcingStepCreate(linear_stepper, nonlinear_stepper, t0, y, sunctx)

    status = ARKodeSetFixedStep(ark.get(), 1e-2)
    assert status == ARK_SUCCESS

    tout = tf
    status, tret = ARKodeEvolve(ark.get(), tout, y, ARK_NORMAL)
    assert status == ARK_SUCCESS

    sol = N_VClone(y)
    ode_problem.solution(y0, sol, tf)
    # we use a fixed atol here since we use a fixed step size
    assert_allclose(N_VGetArrayPointer(sol), N_VGetArrayPointer(y), atol=1e-2)
