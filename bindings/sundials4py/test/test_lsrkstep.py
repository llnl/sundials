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
from problems import AnalyticODE


def test_lsrkstep_sts(sunctx):
    y = N_VNew_Serial(1, sunctx)
    ode_problem = AnalyticODE()
    ode_problem.set_init_cond(y)

    def dom_eig(t, yvec, fnvec, _, tempv1, tempv2, tempv3):
        return ode_problem.dom_eig(t, yvec, fnvec, tempv1, tempv2, tempv3)

    lsrk = LSRKStepCreateSTS(ode_problem.f, 0, y, sunctx)
    status = LSRKStepSetDomEigFn(lsrk.get(), dom_eig)
    assert status == ARK_SUCCESS

    status = ARKodeSStolerances(lsrk.get(), SUNREALTYPE_RTOL, SUNREALTYPE_ATOL)
    assert status == ARK_SUCCESS

    status = ARKodeSetMaxNumSteps(lsrk.get(), 100000)
    assert status == ARK_SUCCESS

    tout = 10.0
    status, tret = ARKodeEvolve(lsrk.get(), tout, y, ARK_NORMAL)
    assert status == ARK_SUCCESS

    sol = N_VClone(y)
    ode_problem.solution(y, sol, tret)
    assert_allclose(N_VGetArrayPointer(sol), N_VGetArrayPointer(y), atol=100 * SUNREALTYPE_RTOL)


def test_lsrkstep_ssp(sunctx):
    y = N_VNew_Serial(1, sunctx)
    ode_problem = AnalyticODE()
    ode_problem.set_init_cond(y)

    def dom_eig(t, yvec, fnvec, _, tempv1, tempv2, tempv3):
        return ode_problem.dom_eig(t, yvec, fnvec, tempv1, tempv2, tempv3)

    lsrk = LSRKStepCreateSSP(ode_problem.f, 0, y, sunctx)
    status = LSRKStepSetDomEigFn(lsrk.get(), dom_eig)
    assert status == ARK_SUCCESS

    status = ARKodeSStolerances(lsrk.get(), SUNREALTYPE_RTOL, SUNREALTYPE_ATOL)
    assert status == ARK_SUCCESS

    status = ARKodeSetMaxNumSteps(lsrk.get(), 100000)
    assert status == ARK_SUCCESS

    tout = 10.0
    status, tret = ARKodeEvolve(lsrk.get(), tout, y, ARK_NORMAL)
    assert status == ARK_SUCCESS

    sol = N_VClone(y)
    ode_problem.solution(y, sol, tret)
    assert_allclose(N_VGetArrayPointer(sol), N_VGetArrayPointer(y), atol=100 * SUNREALTYPE_RTOL)
