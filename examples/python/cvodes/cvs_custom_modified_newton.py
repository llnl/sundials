#!/usr/bin/env python3
# -----------------------------------------------------------------
# Programmer(s): SUNDIALS Developers
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------
"""CVODES with a generic modified-Newton solver implemented in Python.

Unlike the problem-specific quadratic example next to this file, this solver
uses the generic package callbacks: ``sys_fn``, ``lsetup_fn``, ``lsolve_fn``,
and ``conv_test_fn``.  It therefore demonstrates the interface needed by a
solver that can be reused across nonlinear systems.
"""

import numpy as np
from sundials4py.core import *
from sundials4py.cvodes import *


T0 = 0.0
TF = 0.2
RTOL = 1.0e-7
ATOL = 1.0e-10


class LinearODE:
    """The scalar test problem ``y' = -10 y`` with ``y(0) = 1``."""

    def rhs(self, t, y, ydot, _):
        N_VGetNumpyArray(ydot)[0] = -10.0 * N_VGetNumpyArray(y)[0]
        return SUN_SUCCESS

    def jacobian(self, t, y, fy, matrix, _, tmp1, tmp2, tmp3):
        SUNDenseMatrix_Data(matrix)[0, 0] = -10.0
        return SUN_SUCCESS

    def solution(self, t):
        return np.exp(-10.0 * t)


class ModifiedNewtonSolver(CustomSUNNonlinearSolver):
    """A reusable modified-Newton iteration driven by package callbacks."""

    def __init__(self, template, sunctx):
        self.delta = N_VClone(template)
        self.sys_fn = None
        self.lsetup_fn = None
        self.lsolve_fn = None
        self.conv_test_fn = None
        self.max_iters = 4
        self.cur_iter = 0
        self.num_iters = 0
        self.num_conv_fails = 0
        self.total_num_iters = 0
        self.total_num_conv_fails = 0
        self.calls = {"solve": 0, "lsetup": 0, "lsolve": 0, "conv_test": 0}
        super().__init__(sunctx, SUNNONLINEARSOLVER_ROOTFIND)

    def initialize(self):
        return SUN_SUCCESS

    def set_sys_fn(self, fn):
        self.sys_fn = fn
        return SUN_SUCCESS

    def set_lsetup_fn(self, fn):
        self.lsetup_fn = fn
        return SUN_SUCCESS

    def set_lsolve_fn(self, fn):
        self.lsolve_fn = fn
        return SUN_SUCCESS

    def set_conv_test_fn(self, fn):
        self.conv_test_fn = fn
        return SUN_SUCCESS

    def set_max_iters(self, max_iters):
        self.max_iters = max_iters
        return SUN_SUCCESS

    def get_cur_iter(self):
        return SUN_SUCCESS, self.cur_iter

    def get_num_iters(self):
        return SUN_SUCCESS, self.num_iters

    def get_num_conv_fails(self):
        return SUN_SUCCESS, self.num_conv_fails

    def solve(self, y0, y, w, tol, call_lsetup):
        self.calls["solve"] += 1
        # The SUNNonlinearSolver getters report statistics for the most recent
        # solve. CVODE accumulates those values into its integrator statistics.
        self.num_iters = 0
        self.num_conv_fails = 0
        N_VScale(1.0, y0, y)
        jbad = bool(call_lsetup)

        if jbad:
            status, _ = self.lsetup_fn(jbad)
            self.calls["lsetup"] += 1
            if status != SUN_SUCCESS:
                if status == SUN_NLS_CONV_RECVR:
                    self.num_conv_fails += 1
                    self.total_num_conv_fails += 1
                return status

        iteration_started = False
        try:
            for self.cur_iter in range(self.max_iters):
                iteration_started = True
                status = self.sys_fn(y, self.delta)
                if status != SUN_SUCCESS:
                    if status == SUN_NLS_CONV_RECVR:
                        self.num_conv_fails += 1
                    return status

                status = self.lsolve_fn(self.delta)
                self.calls["lsolve"] += 1
                if status != SUN_SUCCESS:
                    if status == SUN_NLS_CONV_RECVR:
                        self.num_conv_fails += 1
                    return status

                N_VLinearSum(1.0, y, -1.0, self.delta, y)
                status = self.conv_test_fn(y, self.delta, tol, w)
                self.calls["conv_test"] += 1
                if status == SUN_SUCCESS:
                    return SUN_SUCCESS
                if status != SUN_NLS_CONTINUE:
                    if status == SUN_NLS_CONV_RECVR:
                        self.num_conv_fails += 1
                    return status

            self.num_conv_fails += 1
            return SUN_NLS_CONV_RECVR
        finally:
            if iteration_started:
                self.num_iters = self.cur_iter + 1
                self.total_num_iters += self.num_iters
            self.total_num_conv_fails += self.num_conv_fails


def main():
    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert status == SUN_SUCCESS

    problem = LinearODE()
    y = N_VNew_Serial(1, sunctx)
    N_VGetNumpyArray(y)[0] = 1.0

    cvode = CVodeCreate(CV_BDF, sunctx)
    assert CVodeInit(cvode.get(), problem.rhs, T0, y) == CV_SUCCESS
    assert CVodeSStolerances(cvode.get(), RTOL, ATOL) == CV_SUCCESS

    matrix = SUNDenseMatrix(1, 1, sunctx)
    linear_solver = SUNLinSol_Dense(y, matrix, sunctx)
    assert CVodeSetLinearSolver(cvode.get(), linear_solver, matrix) == CV_SUCCESS
    assert CVodeSetJacFn(cvode.get(), problem.jacobian) == CV_SUCCESS

    nonlinear_solver = ModifiedNewtonSolver(y, sunctx)
    assert CVodeSetNonlinearSolver(cvode.get(), nonlinear_solver) == CV_SUCCESS

    status, tret = CVode(cvode.get(), TF, y, CV_NORMAL)
    assert status == CV_SUCCESS
    computed = float(N_VGetNumpyArray(y)[0])
    exact = float(problem.solution(tret))
    status, nni = CVodeGetNumNonlinSolvIters(cvode.get())
    assert status == CV_SUCCESS
    calls = dict(nonlinear_solver.calls)
    calls["num_iters"] = nonlinear_solver.total_num_iters
    calls["num_conv_fails"] = nonlinear_solver.total_num_conv_fails
    calls["cvode_num_iters"] = nni
    return computed, exact, calls


def test_cvs_custom_modified_newton():
    computed, exact, calls = main()
    assert np.isclose(computed, exact, rtol=1.0e-5, atol=1.0e-7)
    assert calls["solve"] > 0
    assert calls["lsetup"] > 0
    assert calls["lsolve"] > 0
    assert calls["conv_test"] > 0
    assert calls["num_iters"] >= calls["solve"]
    assert calls["cvode_num_iters"] == calls["num_iters"]


if __name__ == "__main__":
    main()
