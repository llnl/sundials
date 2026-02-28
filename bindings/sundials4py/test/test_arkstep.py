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
import numpy as np
from numpy.testing import assert_allclose
from fixtures import *
from sundials4py.core import *
from sundials4py.arkode import *
from problems import AnalyticODE, AnalyticMultiscaleODE


def test_explicit(sunctx):
    y = N_VNew_Serial(1, sunctx)

    ode_problem = AnalyticODE()

    ode_problem.set_init_cond(y)

    ark = ARKStepCreate(ode_problem.f, None, 0, y, sunctx)

    status = ARKodeSStolerances(ark.get(), SUNREALTYPE_RTOL, SUNREALTYPE_ATOL)
    assert status == ARK_SUCCESS

    nrtfn = 2

    def rootfn(t, y, gout, _):
        # just a smoke test of the root finding callback
        gout[:] = 1.0
        assert len(gout) == nrtfn
        return 0

    status = ARKodeRootInit(ark.get(), nrtfn, rootfn)
    assert status == ARK_SUCCESS

    tout = 10.0
    status, tret = ARKodeEvolve(ark.get(), tout, y, ARK_NORMAL)
    assert status == ARK_SUCCESS

    status, num_steps = ARKodeGetNumSteps(ark.get())
    assert status == ARK_SUCCESS
    assert num_steps > 0

    sol = N_VClone(y)
    ode_problem.solution(y, sol, tret)

    assert_allclose(N_VGetArrayPointer(sol), N_VGetArrayPointer(y), atol=100 * SUNREALTYPE_RTOL)


def test_implicit(sunctx):
    y = N_VNew_Serial(1, sunctx)
    ls = SUNLinSol_SPGMR(y, 0, 0, sunctx)

    ode_problem = AnalyticODE()

    ode_problem.set_init_cond(y)

    ark = ARKStepCreate(None, ode_problem.f, 0, y, sunctx)

    status = ARKodeSStolerances(ark.get(), SUNREALTYPE_RTOL, SUNREALTYPE_ATOL)
    assert status == ARK_SUCCESS

    status = ARKodeSetLinearSolver(ark.get(), ls, None)
    assert status == ARK_SUCCESS

    tout = 10.0
    status, tret = ARKodeEvolve(ark.get(), tout, y, ARK_NORMAL)
    assert status == ARK_SUCCESS

    status, num_steps = ARKodeGetNumSteps(ark.get())
    assert status == ARK_SUCCESS
    assert num_steps > 0

    sol = N_VClone(y)
    ode_problem.solution(y, sol, tret)

    assert_allclose(N_VGetArrayPointer(sol), N_VGetArrayPointer(y), atol=100 * SUNREALTYPE_RTOL)


def test_implicit_with_dense_ls_and_jac(sunctx):
    y = N_VNew_Serial(1, sunctx)

    ode_problem = AnalyticODE()
    ode_problem.set_init_cond(y)

    A = SUNDenseMatrix(1, 1, sunctx)
    ls = SUNLinSol_Dense(y, A, sunctx)

    ark = ARKStepCreate(None, ode_problem.f, 0.0, y, sunctx)

    status = ARKodeSStolerances(ark.get(), SUNREALTYPE_RTOL, SUNREALTYPE_ATOL)
    assert status == ARK_SUCCESS

    status = ARKodeSetLinearSolver(ark.get(), ls, A)
    assert status == ARK_SUCCESS

    status = ARKodeSetJacFn(ark.get(), ode_problem.jac_fn)
    assert status == ARK_SUCCESS

    tout = 10.0
    status, tret = ARKodeEvolve(ark.get(), tout, y, ARK_NORMAL)
    assert status == ARK_SUCCESS

    status, num_steps = ARKodeGetNumSteps(ark.get())
    assert status == ARK_SUCCESS
    assert num_steps > 0

    sol = N_VClone(y)
    ode_problem.solution(y, sol, tret)

    assert_allclose(N_VGetArrayPointer(sol), N_VGetArrayPointer(y), atol=100 * SUNREALTYPE_RTOL)


def test_imex(sunctx):
    y = N_VNew_Serial(1, sunctx)
    ls = SUNLinSol_SPGMR(y, 0, 0, sunctx)

    ode_problem = AnalyticMultiscaleODE()

    ode_problem.set_init_cond(y)

    ark = ARKStepCreate(ode_problem.f_nonlinear, ode_problem.f_linear, 0, y, sunctx)

    status = ARKodeSStolerances(ark.get(), SUNREALTYPE_RTOL, SUNREALTYPE_ATOL)
    assert status == ARK_SUCCESS

    status = ARKodeSetLinearSolver(ark.get(), ls, None)
    assert status == ARK_SUCCESS

    tout = 10.0
    status, tret = ARKodeEvolve(ark.get(), tout, y, ARK_NORMAL)
    assert status == ARK_SUCCESS

    status, num_steps = ARKodeGetNumSteps(ark.get())
    assert status == ARK_SUCCESS
    assert num_steps > 0

    sol = N_VClone(y)
    ode_problem.solution(y, sol, tret)

    assert_allclose(N_VGetArrayPointer(sol), N_VGetArrayPointer(y), atol=100 * SUNREALTYPE_RTOL)
