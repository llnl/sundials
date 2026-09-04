# -----------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ UMBC
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
# End-to-end tests that drive Python-implemented SUNDIALS objects from real
# integrators.
#
# The per-family tests check each operation in isolation; these check that a
# package -- which calls the operations in an order the binding does not choose,
# and which supplies the memory pointer the callback adapters need -- reaches the
# same Python code and produces the right answer.
# -----------------------------------------------------------------

import numpy as np
import pytest
from fixtures import *
from numpy.testing import assert_allclose
from problems import AnalyticNonlinearSys, AnalyticODE
from sundials4py.arkode import *
from sundials4py.core import *
from sundials4py.cvodes import *
from sundials4py.kinsol import *


class DenseArrayMatrix(CustomSUNMatrix):
    # A dense SUNMatrix whose storage is an ordinary NumPy array. Everything a
    # package asks of a matrix is expressed in NumPy terms.
    def __init__(self, rows, cols, sunctx, data=None):
        self.data = np.zeros((rows, cols), dtype=sunrealtype) if data is None else data
        super().__init__(sunctx)

    def clone(self):
        return DenseArrayMatrix(*self.data.shape, self.sunctx)

    def zero(self):
        self.data[:] = 0.0
        return SUN_SUCCESS

    def copy(self, dst):
        dst.data[:] = self.data
        return SUN_SUCCESS

    def scaleadd(self, c, other):
        self.data[:] = c * self.data + other.data
        return SUN_SUCCESS

    def scaleaddi(self, c):
        self.data[:] = c * self.data
        self.data[np.diag_indices_from(self.data)] += 1.0
        return SUN_SUCCESS

    def matvec(self, x, y):
        N_VGetArrayPointer(y)[:] = self.data @ N_VGetArrayPointer(x)
        return SUN_SUCCESS


class NumpyLinearSolver(CustomSUNLinearSolver):
    # A direct linear solver that factors the matrix with NumPy.
    #
    # setup() and solve() receive the SUNMatrix as an opaque native handle, not
    # as the Python object behind it, so the matrix to read is supplied at
    # construction. This is the normal pattern for a Python implementation: the
    # solver and the matrix it is paired with are written together.
    def __init__(self, matrix, sunctx):
        self.matrix = matrix
        self.factor = None
        self.calls = {"initialize": 0, "setup": 0, "solve": 0}
        super().__init__(sunctx, SUNLINEARSOLVER_DIRECT)

    def initialize(self):
        self.calls["initialize"] += 1
        return SUN_SUCCESS

    def setup(self, A):
        self.calls["setup"] += 1
        # Copy, because the package is free to overwrite the matrix afterwards.
        self.factor = self.matrix.data.copy()
        return SUN_SUCCESS

    def solve(self, A, x, b, tol):
        self.calls["solve"] += 1
        try:
            N_VGetArrayPointer(x)[:] = np.linalg.solve(self.factor, N_VGetArrayPointer(b))
        except np.linalg.LinAlgError:
            # A singular matrix is a recoverable failure, not a Python error.
            return 1
        return SUN_SUCCESS


class NewtonNonlinearSolver(CustomSUNNonlinearSolver):
    # A modified-Newton rootfinding solver written entirely in Python.
    #
    # Every callback it uses -- the residual, the linear-solver setup and solve,
    # and the convergence test -- arrives from the integrator through a setter as
    # an ordinary Python callable, and each is valid only inside solve(). This is
    # the fullest exercise of the active-memory scope, because the memory pointer
    # SUNDIALS threads through those callbacks is the integrator's own.
    def __init__(self, template, sunctx):
        self.delta = N_VClone(template)
        self.sys_fn = None
        self.lsetup_fn = None
        self.lsolve_fn = None
        self.conv_test_fn = None
        self.max_iters = 3
        self.cur_iter = 0
        self.num_iters = 0
        self.num_conv_fails = 0
        super().__init__(sunctx, SUNNONLINEARSOLVER_ROOTFIND)

    def initialize(self):
        return SUN_SUCCESS

    def set_sys_fn(self, sys_fn):
        self.sys_fn = sys_fn
        return SUN_SUCCESS

    def set_lsetup_fn(self, lsetup_fn):
        self.lsetup_fn = lsetup_fn
        return SUN_SUCCESS

    def set_lsolve_fn(self, lsolve_fn):
        self.lsolve_fn = lsolve_fn
        return SUN_SUCCESS

    def set_conv_test_fn(self, conv_test_fn):
        self.conv_test_fn = conv_test_fn
        return SUN_SUCCESS

    def set_max_iters(self, maxiters):
        self.max_iters = maxiters
        return SUN_SUCCESS

    def get_num_iters(self):
        return SUN_SUCCESS, self.num_iters

    def get_cur_iter(self):
        # CVODE reads this during the linear-solver setup decision, so it must
        # reflect the iteration currently in progress rather than a total.
        return SUN_SUCCESS, self.cur_iter

    def get_num_conv_fails(self):
        return SUN_SUCCESS, self.num_conv_fails

    def solve(self, y0, y, w, tol, call_lsetup):
        # y0 is the predicted correction and y the correction being solved for;
        # this is the loop SUNNonlinSol_Newton implements in C.
        N_VScale(1.0, y0, y)
        jbad = bool(call_lsetup)
        jcur = False

        while True:
            if jbad:
                status, jcur = self.lsetup_fn(jbad)
                if status != 0:
                    return status
                jbad = False

            status = SUN_NLS_CONTINUE
            for self.cur_iter in range(self.max_iters):
                self.num_iters += 1

                status = self.sys_fn(y, self.delta)
                if status != 0:
                    break

                status = self.lsolve_fn(self.delta)
                if status != 0:
                    break

                N_VLinearSum(1.0, y, -1.0, self.delta, y)

                status = self.conv_test_fn(y, self.delta, tol, w)
                if status == SUN_SUCCESS:
                    return SUN_SUCCESS
                if status != SUN_NLS_CONTINUE:
                    break

            self.num_conv_fails += 1

            # An unrecoverable failure ends the solve. A recoverable one may be
            # retried, but only once the Jacobian has been refreshed -- retrying
            # with the same stale Jacobian would loop forever.
            if status < 0:
                return status
            if jcur:
                return SUN_NLS_CONV_RECVR
            jbad = True


class ProportionalController(CustomSUNHController):
    # A deliberately simple I-controller, so that the test can assert both that
    # ARKODE consulted it and that the steps it asked for were honored.
    def __init__(self, sunctx, safety=0.9):
        self.safety = safety
        self.bias = 1.0
        self.calls = {"estimate_step": 0, "reset": 0, "set_defaults": 0, "update_h": 0}
        super().__init__(sunctx)

    def estimate_step(self, h, p, dsm):
        self.calls["estimate_step"] += 1
        e = max(self.bias * dsm, 1.0e-10)
        return SUN_SUCCESS, self.safety * h * e ** (-1.0 / (p + 1))

    def reset(self):
        self.calls["reset"] += 1
        return SUN_SUCCESS

    def set_defaults(self):
        self.calls["set_defaults"] += 1
        self.bias = 1.0
        return SUN_SUCCESS

    def set_error_bias(self, bias):
        self.bias = bias
        return SUN_SUCCESS

    def update_h(self, h, dsm):
        self.calls["update_h"] += 1
        return SUN_SUCCESS


def test_custom_matrix_and_linear_solver_through_kinsol(sunctx):
    # Purpose:
    # Custom matrix and linear solver through kinsol.
    # KINSOL Newton, with both the Jacobian matrix and the linear solver
    # implemented in Python.
    NEQ = AnalyticNonlinearSys.NEQ
    problem = AnalyticNonlinearSys(None)
    kin = KINCreate(sunctx)
    u = N_VNew_Serial(NEQ, sunctx)

    def sys_function(u, f, _):
        # Residual form f(u) = g(u) - u of the fixed-point problem.
        problem.fixed_point_fn(u, f, None)
        N_VLinearSum(1.0, f, -1.0, u, f)
        return 0

    J = DenseArrayMatrix(NEQ, NEQ, sunctx)
    LS = NumpyLinearSolver(J, sunctx)

    def jac_fn(uvec, fuvec, Jmat, _, tmp1, tmp2):
        # Analytic Jacobian of f(u) = g(u) - u, written into the Python matrix.
        # Jmat is the same matrix as J, but arrives as an opaque native handle.
        x, y, z = N_VGetArrayPointer(uvec)
        r = np.sqrt(x * x + np.sin(z) + 1.06)
        e = np.exp(-x * (y - 1.0))
        J.data[:] = [
            [
                -1.0,
                -(1.0 / 3.0) * z * np.sin((y - 1.0) * z),
                -(1.0 / 3.0) * (y - 1.0) * np.sin((y - 1.0) * z),
            ],
            [x / (9.0 * r), -1.0, np.cos(z) / (18.0 * r)],
            [(y - 1.0) * e / 20.0, x * e / 20.0, -1.0],
        ]
        return 0

    assert KINInit(kin.get(), sys_function, u) == KIN_SUCCESS
    assert KINSetFuncNormTol(kin.get(), SUNREALTYPE_ATOL) == KIN_SUCCESS
    assert KINSetLinearSolver(kin.get(), LS, J) == KIN_SUCCESS
    assert KINSetJacFn(kin.get(), jac_fn) == KIN_SUCCESS

    N_VGetArrayPointer(u)[:] = [0.1, 0.1, -0.1]
    scale = N_VNew_Serial(NEQ, sunctx)
    N_VConst(1.0, scale)

    assert KINSol(kin.get(), u, KIN_NONE, scale, scale) == KIN_SUCCESS

    expected = N_VNew_Serial(NEQ, sunctx)
    problem.solution(expected)
    assert_allclose(
        N_VGetArrayPointer(u), N_VGetArrayPointer(expected), atol=100 * SUNREALTYPE_ATOL
    )

    # The package really did drive the Python implementations.
    assert LS.calls["initialize"] >= 1
    assert LS.calls["setup"] >= 1
    assert LS.calls["solve"] >= 1


def test_custom_nonlinear_solver_through_cvode(sunctx):
    # Purpose:
    # Custom nonlinear solver through cvode.
    # CVODE BDF with a Python nonlinear solver and a native dense linear solver,
    # so the only custom object in the loop is the one under test.
    problem = AnalyticODE(lamb=-10.0)
    y = N_VNew_Serial(1, sunctx)
    problem.set_init_cond(y)

    cvode = CVodeCreate(CV_BDF, sunctx)
    assert CVodeInit(cvode.get(), problem.f, 0.0, y) == CV_SUCCESS
    assert CVodeSStolerances(cvode.get(), 1.0e-8, 1.0e-10) == CV_SUCCESS

    A = SUNDenseMatrix(1, 1, sunctx)
    LS = SUNLinSol_Dense(y, A, sunctx)
    assert CVodeSetLinearSolver(cvode.get(), LS, A) == CV_SUCCESS
    assert CVodeSetJacFn(cvode.get(), problem.jac_fn) == CV_SUCCESS

    NLS = NewtonNonlinearSolver(y, sunctx)
    assert CVodeSetNonlinearSolver(cvode.get(), NLS) == CV_SUCCESS

    tf = 1.0
    status, tret = CVode(cvode.get(), tf, y, CV_NORMAL)
    assert status == CV_SUCCESS
    assert tret == pytest.approx(tf)

    expected = N_VNew_Serial(1, sunctx)
    problem.solution(None, expected, tf)
    assert_allclose(N_VGetArrayPointer(y), N_VGetArrayPointer(expected), rtol=1.0e-4, atol=1.0e-6)

    # CVODE drove the Python solver rather than falling back to its own.
    assert NLS.num_iters > 0


def test_custom_adapt_controller_through_arkode(sunctx):
    # Purpose:
    # Custom adapt controller through arkode.
    # ARKODE ERK with a Python time step controller.
    problem = AnalyticODE(lamb=-10.0)
    y = N_VNew_Serial(1, sunctx)
    problem.set_init_cond(y)

    ark = ERKStepCreate(problem.f, 0.0, y, sunctx)
    controller = ProportionalController(sunctx)

    assert ARKodeSStolerances(ark.get(), 1.0e-6, 1.0e-10) == ARK_SUCCESS
    assert ARKodeSetAdaptController(ark.get(), controller) == ARK_SUCCESS

    tf = 1.0
    status, tret = ARKodeEvolve(ark.get(), tf, y, ARK_NORMAL)
    assert status == ARK_SUCCESS
    assert tret == pytest.approx(tf)

    expected = N_VNew_Serial(1, sunctx)
    problem.solution(None, expected, tf)
    assert_allclose(N_VGetArrayPointer(y), N_VGetArrayPointer(expected), rtol=1.0e-3, atol=1.0e-5)

    # Adaptivity ran through Python, over more than one step.
    assert controller.calls["estimate_step"] > 1
    status, nsteps = ARKodeGetNumSteps(ark.get())
    assert status == ARK_SUCCESS
    assert nsteps > 1


def test_custom_adapt_controller_survives_repeated_evolutions(sunctx):
    # Purpose:
    # Custom adapt controller survives repeated evolutions.
    # 11.3: an H controller must keep working across many ARKODE steps and a
    # reset, which is where a stale weak reference or a released handle shows up.
    problem = AnalyticODE(lamb=-10.0)
    y = N_VNew_Serial(1, sunctx)
    problem.set_init_cond(y)

    ark = ERKStepCreate(problem.f, 0.0, y, sunctx)
    controller = ProportionalController(sunctx)
    assert ARKodeSStolerances(ark.get(), 1.0e-6, 1.0e-10) == ARK_SUCCESS
    assert ARKodeSetAdaptController(ark.get(), controller) == ARK_SUCCESS

    t = 0.0
    for tout in (0.25, 0.5, 0.75, 1.0):
        status, t = ARKodeEvolve(ark.get(), tout, y, ARK_NORMAL)
        assert status == ARK_SUCCESS
        assert t == pytest.approx(tout)

    expected = N_VNew_Serial(1, sunctx)
    problem.solution(None, expected, 1.0)
    assert_allclose(N_VGetArrayPointer(y), N_VGetArrayPointer(expected), rtol=1.0e-3, atol=1.0e-5)
    assert controller.calls["estimate_step"] > 4
