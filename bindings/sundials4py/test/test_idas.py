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
from sundials4py.core import *
from sundials4py.idas import *
from problems import AnalyticDAE
from fixtures import *


# This test requires more precision than float32 offers to converge
@pytest.mark.skipif(
    sunrealtype == np.float32, reason="Test not supported for sunrealtype=np.float32"
)
def test_idas_ivp(sunctx):
    dae_problem = AnalyticDAE()

    ida = IDACreate(sunctx)
    yy = N_VNew_Serial(2, sunctx)
    yp = N_VNew_Serial(2, sunctx)

    # y and y' initial conditions
    dae_problem.set_init_cond(yy, yp, dae_problem.T0)

    ls = SUNLinSol_SPGMR(yy, SUN_PREC_LEFT, 0, sunctx)

    def psolve(t, yy, yp, rr, r, z, cj, delta, _):
        return dae_problem.psolve(t, yy, yp, rr, r, z, cj, delta)

    status = IDAInit(ida.get(), dae_problem.res, 0.0, yy, yp)
    assert status == IDA_SUCCESS

    status = IDASStolerances(ida.get(), 100 * SUNREALTYPE_RTOL, SUNREALTYPE_ATOL)
    assert status == IDA_SUCCESS

    status = IDASetLinearSolver(ida.get(), ls, None)
    assert status == IDA_SUCCESS

    status = IDASetPreconditioner(ida.get(), None, psolve)
    assert status == IDA_SUCCESS

    status = IDASetMaxNumSteps(ida.get(), 100000)
    assert status == IDA_SUCCESS

    nrtfn = 2

    def rootfn(t, yy, yp, gout, _):
        # just a smoke test of the root finding callback
        gout[:] = 1.0
        assert len(gout) == nrtfn
        return 0

    status = IDARootInit(ida.get(), nrtfn, rootfn)
    assert status == IDA_SUCCESS

    tout = dae_problem.TF
    status, tret = IDASolve(ida.get(), tout, yy, yp, IDA_NORMAL)
    assert status == IDA_SUCCESS

    status, num_steps = IDAGetNumSteps(ida.get())
    assert status == IDA_SUCCESS
    print("Number of steps: ", num_steps)

    sol_yy = N_VClone(yy)
    sol_yp = N_VClone(yp)

    dae_problem.solution(sol_yy, sol_yp, tret)
    assert_allclose(N_VGetArrayPointer(sol_yy), N_VGetArrayPointer(yy), rtol=1e-2)
    assert_allclose(N_VGetArrayPointer(sol_yp), N_VGetArrayPointer(yp), rtol=1e-2)


def test_idas_fsa(sunctx):
    dae_problem = AnalyticDAE()

    ida = IDACreate(sunctx)
    yy = N_VNew_Serial(2, sunctx)
    yp = N_VNew_Serial(2, sunctx)

    # y and y' initial conditions
    dae_problem.set_init_cond(yy, yp, dae_problem.T0)

    ls = SUNLinSol_SPGMR(yy, SUN_PREC_LEFT, 0, sunctx)

    def psolve(t, yy, yp, rr, r, z, cj, delta, _):
        return dae_problem.psolve(t, yy, yp, rr, r, z, cj, delta)

    status = IDAInit(ida.get(), dae_problem.res, 0.0, yy, yp)
    assert status == IDA_SUCCESS

    status = IDASStolerances(ida.get(), 100 * SUNREALTYPE_RTOL, SUNREALTYPE_ATOL)
    assert status == IDA_SUCCESS

    status = IDASetLinearSolver(ida.get(), ls, None)
    assert status == IDA_SUCCESS

    status = IDASetPreconditioner(ida.get(), None, psolve)
    assert status == IDA_SUCCESS

    status = IDASetMaxNumSteps(ida.get(), 100000)
    assert status == IDA_SUCCESS

    # Sensitivity setup
    Ns = 1  # number of sensitivities (wrt yy0)
    ism = IDA_SIMULTANEOUS
    yyS0 = [N_VClone(yy) for _ in range(Ns)]
    ypS0 = [N_VClone(yp) for _ in range(Ns)]
    for v in yyS0:
        N_VConst(0.0, v)
    for v in ypS0:
        N_VConst(0.0, v)

    def resS(Ns, t, yy, yp, resval, yS, ypS, resvalS, _, tmp1, tmp2, tmp3):
        # Sensitivity residuals: d(res)/d(yy) * yyS + d(res)/d(yp) * ypS
        # For smoke test, just zero out
        for i in range(Ns):
            N_VConst(0.0, resvalS[i])
        return 0

    status = IDASensInit(ida.get(), Ns, ism, resS, yyS0, ypS0)
    assert status == IDA_SUCCESS

    status = IDASensSStolerances(
        ida.get(), 100 * SUNREALTYPE_RTOL, np.array([SUNREALTYPE_ATOL], dtype=sunrealtype)
    )
    assert status == IDA_SUCCESS

    tout = dae_problem.TF
    yySout = [N_VClone(yy) for _ in range(Ns)]
    status, tret = IDASolve(ida.get(), tout, yy, yp, IDA_NORMAL)
    assert status == IDA_SUCCESS

    # Get sensitivities (smoke test: just check call and shape)
    status, tret_sens = IDAGetSens(ida.get(), yySout)
    assert status == IDA_SUCCESS
    assert len(yySout) == Ns


def test_idas_adjoint(sunctx):
    dae_problem = AnalyticDAE()

    ida = IDACreate(sunctx)
    yy = N_VNew_Serial(2, sunctx)
    yp = N_VNew_Serial(2, sunctx)

    # y and y' initial conditions
    dae_problem.set_init_cond(yy, yp, dae_problem.T0)

    ls = SUNLinSol_SPGMR(yy, SUN_PREC_LEFT, 0, sunctx)

    def psolve(t, yy, yp, rr, r, z, cj, delta, _):
        return dae_problem.psolve(t, yy, yp, rr, r, z, cj, delta)

    status = IDAInit(ida.get(), dae_problem.res, 0.0, yy, yp)
    assert status == IDA_SUCCESS

    status = IDASStolerances(ida.get(), 100 * SUNREALTYPE_RTOL, SUNREALTYPE_ATOL)
    assert status == IDA_SUCCESS

    status = IDASetLinearSolver(ida.get(), ls, None)
    assert status == IDA_SUCCESS

    status = IDASetPreconditioner(ida.get(), None, psolve)
    assert status == IDA_SUCCESS

    status = IDASetMaxNumSteps(ida.get(), 100000)
    assert status == IDA_SUCCESS

    # Adjoint (backward) problem setup
    steps = 5
    interp = IDA_HERMITE
    status = IDAAdjInit(ida.get(), steps, interp)
    assert status == IDA_SUCCESS

    # Integrate forward
    tout = dae_problem.TF
    status, tret, ncheck = IDASolveF(ida.get(), tout, yy, yp, IDA_NORMAL)
    assert status == IDA_SUCCESS

    # Create backward problem
    yyB = N_VClone(yy)
    ypB = N_VClone(yp)
    N_VConst(0.0, yyB)
    N_VConst(0.0, ypB)
    lsB = SUNLinSol_SPGMR(yyB, 0, 0, sunctx)

    status, whichB = IDACreateB(ida.get())
    assert status == IDA_SUCCESS

    def resB(t, y, yp, yB, ypB, resvalB, _):
        N_VConst(0.0, resvalB)
        return 0

    status = IDAInitB(ida.get(), whichB, resB, tout, yyB, ypB)
    assert status == IDA_SUCCESS

    status = IDASStolerancesB(ida.get(), whichB, 100 * SUNREALTYPE_RTOL, SUNREALTYPE_ATOL)
    assert status == IDA_SUCCESS

    status = IDASetLinearSolverB(ida.get(), whichB, lsB, None)
    assert status == IDA_SUCCESS

    # Integrate backward
    tB0 = dae_problem.T0
    status = IDASolveB(ida.get(), tB0, IDA_NORMAL)
    assert status >= IDA_SUCCESS

    # Get sensitivities
    yB0 = N_VClone(yyB)
    ypB0 = N_VClone(ypB)
    status, yyB0 = IDAGetB(ida.get(), whichB, yB0, ypB0)
    assert status == IDA_SUCCESS
