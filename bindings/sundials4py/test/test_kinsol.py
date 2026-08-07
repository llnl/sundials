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
from sundials4py.kinsol import *
from problems import AnalyticNonlinearSys


def test_kinsol(sunctx):
    NEQ = 3
    m_aa = 2
    tol = SUNREALTYPE_ATOL
    kin_view = KINCreate(sunctx)
    problem = AnalyticNonlinearSys(None)
    u = N_VNew_Serial(NEQ, sunctx)

    def depth_fn(iter, u_val, g_val, f_val, df, R_mat, depth, _, remove_indices):
        if iter < 2:
            new_depth = 1
        else:
            new_depth = depth
        return 0, new_depth

    def damping_fn(iter, u_val, g_val, qt_fn, depth, _):
        # smoke test of damping
        damping_factor = 1.0
        return 0, damping_factor

    kin_status = KINSetMAA(kin_view.get(), m_aa)
    assert kin_status == KIN_SUCCESS
    kin_status = KINInit(kin_view.get(), problem.fixed_point_fn, u)
    assert kin_status == KIN_SUCCESS
    kin_status = KINSetFuncNormTol(kin_view.get(), tol)
    assert kin_status == KIN_SUCCESS
    kin_status = KINSetDepthFn(kin_view.get(), depth_fn)
    assert kin_status == KIN_SUCCESS
    kin_status = KINSetDampingFn(kin_view.get(), damping_fn)
    assert kin_status == KIN_SUCCESS

    # initial guess
    u_data = N_VGetArrayPointer(u)
    u_data[:] = [0.1, 0.1, -0.1]

    # no scaling used
    scale = N_VNew_Serial(NEQ, sunctx)
    N_VConst(1.0, scale)

    kin_status = KINSol(kin_view.get(), u, KIN_FP, scale, scale)
    assert kin_status == KIN_SUCCESS

    u_expected = N_VNew_Serial(NEQ, sunctx)
    u_expected_data = N_VGetArrayPointer(u_expected)
    problem.solution(u_expected)
    assert_allclose(u_data, u_expected_data, atol=10 * SUNREALTYPE_ATOL)


def test_kinsol_newton_dense(sunctx):
    NEQ = AnalyticNonlinearSys.NEQ
    tol = SUNREALTYPE_ATOL
    kin_view = KINCreate(sunctx)
    problem = AnalyticNonlinearSys(None)
    u = N_VNew_Serial(NEQ, sunctx)

    def sys_function(u, f, _):
        # Form the nonlinear residual f(u) = g(u) - u
        problem.fixed_point_fn(u, f, None)
        N_VLinearSum(1.0, f, -1.0, u, f)
        return 0

    kin_status = KINInit(kin_view.get(), sys_function, u)
    assert kin_status == KIN_SUCCESS
    kin_status = KINSetFuncNormTol(kin_view.get(), tol)
    assert kin_status == KIN_SUCCESS

    # Attach a dense linear solver
    J = SUNDenseMatrix(NEQ, NEQ, sunctx)
    LS = SUNLinSol_Dense(u, J, sunctx)
    kin_status = KINSetLinearSolver(kin_view.get(), LS, J)
    assert kin_status == KIN_SUCCESS

    # initial guess
    u_data = N_VGetArrayPointer(u)
    u_data[:] = [0.1, 0.1, -0.1]

    # no scaling used
    scale = N_VNew_Serial(NEQ, sunctx)
    N_VConst(1.0, scale)

    kin_status = KINSol(kin_view.get(), u, KIN_NONE, scale, scale)
    assert kin_status == KIN_SUCCESS

    u_expected = N_VNew_Serial(NEQ, sunctx)
    u_expected_data = N_VGetArrayPointer(u_expected)
    problem.solution(u_expected)
    assert_allclose(u_data, u_expected_data, atol=10 * SUNREALTYPE_ATOL)
