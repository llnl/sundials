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
# Unit/smoke tests for SUNNonlinearSolver module
# -----------------------------------------------------------------

import pytest
import numpy as np
from numpy.testing import assert_allclose
from fixtures import *
from sundials4py.core import *


@pytest.mark.parametrize("solver_type", ["fixedpoint", "newton"])
def test_create_solver(solver_type, sunctx, nvec):
    if solver_type == "fixedpoint":
        nls = SUNNonlinSol_FixedPoint(nvec, 5, sunctx)
    elif solver_type == "newton":
        nls = SUNNonlinSol_Newton(nvec, sunctx)
    else:
        raise ValueError("Unknown solver type")

    assert nls is not None


@pytest.mark.parametrize("solver_type", ["fixedpoint", "newton"])
def make_solver(solver_type, sunctx, nvec):
    if solver_type == "fixedpoint":
        nls = SUNNonlinSol_FixedPoint(nvec, 5, sunctx)
    elif solver_type == "newton":
        nls = SUNNonlinSol_Newton(nvec, sunctx)
    else:
        raise ValueError("Unknown solver type")

    def ctest(nls, u, delta, tol, ewt, _):
        return 0

    def sys(u, g, _):
        return 0

    def lsolve(u, _):
        return 0

    # Sys and CTest must be set before calling certain functions when error checking is enabled
    SUNNonlinSolSetConvTestFn(nls, ctest)
    SUNNonlinSolSetSysFn(nls, sys)

    # LSolve must also be set before calling certain functions when error checking is enabled
    if solver_type == "newton":
        SUNNonlinSolSetLSolveFn(nls, lsolve)

    return nls


@pytest.mark.parametrize(
    "solver_type, expected_type",
    [("newton", SUNNONLINEARSOLVER_ROOTFIND), ("fixedpoint", SUNNONLINEARSOLVER_FIXEDPOINT)],
)
def test_gettype(solver_type, expected_type, sunctx, nvec):
    nls = make_solver(solver_type, sunctx, nvec)
    typ = SUNNonlinSolGetType(nls)
    assert typ is expected_type


@pytest.mark.parametrize("solver_type", ["fixedpoint", "newton"])
def test_initialize(solver_type, sunctx, nvec):
    nls = make_solver(solver_type, sunctx, nvec)

    ret = SUNNonlinSolInitialize(nls)
    assert ret == SUN_SUCCESS


@pytest.mark.parametrize("solver_type,max_iters", [("newton", 5), ("fixedpoint", 10)])
def test_set_max_iters_and_get_num_iters(solver_type, max_iters, sunctx, nvec):
    nls = make_solver(solver_type, sunctx, nvec)
    ret = SUNNonlinSolSetMaxIters(nls, max_iters)
    assert ret == SUN_SUCCESS
    err, niters = SUNNonlinSolGetNumIters(nls)
    assert err == SUN_SUCCESS
    assert isinstance(niters, int)


@pytest.mark.parametrize("solver_type", ["fixedpoint", "newton"])
def test_get_cur_iter(solver_type, sunctx, nvec):
    nls = make_solver(solver_type, sunctx, nvec)
    err, cur_iter = SUNNonlinSolGetCurIter(nls)
    assert err == SUN_SUCCESS
    assert isinstance(cur_iter, int)


@pytest.mark.parametrize("solver_type", ["fixedpoint", "newton"])
def test_get_num_conv_fails(solver_type, sunctx, nvec):
    nls = make_solver(solver_type, sunctx, nvec)
    err, nconvfails = SUNNonlinSolGetNumConvFails(nls)
    assert err == SUN_SUCCESS
    assert isinstance(nconvfails, int)


def test_fixedpoint_setup_and_solve(sunctx):
    from problems import AnalyticNonlinearSys

    NEQ = AnalyticNonlinearSys.NEQ
    ucor = N_VNew_Serial(NEQ, sunctx)
    u0 = N_VNew_Serial(NEQ, sunctx)
    w = N_VNew_Serial(NEQ, sunctx)
    ucur = N_VNew_Serial(NEQ, sunctx)

    # Initial guess
    udata = N_VGetArrayPointer(u0)
    udata[:] = [0.1, 0.1, -0.1]

    # Initial correction
    N_VConst(0.0, ucor)

    # Set the weights
    N_VConst(1.0, w)

    # Create the problem
    with AnalyticNonlinearSys(u0) as problem:

        # Create the solver
        nls = SUNNonlinSol_FixedPoint(u0, 2, sunctx)

        # System function
        def g_fn(u, g, _):
            return problem.corrector_fp_fn(u, g)

        # Convergence test
        def conv_test(nls, u, delta, tol, ewt, _):
            return problem.conv_test(nls, u, delta, tol, ewt)

        ret = SUNNonlinSolSetSysFn(nls, g_fn)
        assert ret == SUN_SUCCESS

        ret = SUNNonlinSolSetConvTestFn(nls, conv_test)
        assert ret == SUN_SUCCESS

        ret = SUNNonlinSolSetMaxIters(nls, 500)

        ret = SUNNonlinSolSetup(nls, u0)
        assert ret == SUN_SUCCESS

        tol = SUNREALTYPE_ATOL
        ret = SUNNonlinSolSolve(nls, u0, ucor, w, tol, 0)
        assert ret == SUN_SUCCESS

        # Update the initial guess with the correction
        N_VLinearSum(1.0, u0, 1.0, ucor, ucur)

        # Compare to analytic solution
        utrue = N_VNew_Serial(NEQ, sunctx)
        problem.solution(utrue)
        # 1e-2 seems to be the best we can do here
        assert_allclose(N_VGetArrayPointer(ucur), N_VGetArrayPointer(utrue), atol=1e-2)
