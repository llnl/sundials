#!/usr/bin/env python3
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
# TEMPLATE: implementing a SUNNonlinearSolver in Python.
#
# This file is a starting point, not a finished example. The scaffolding -- the
# problem, the CVODE setup, the output -- is complete; the parts marked
# `TODO(you):` are where your nonlinear solver goes. Fill them in, delete the
# pytest.skip() at the bottom, and run the file.
#
# The test problem is the scalar ODE
#
#    dy/dt = lambda*y + 1/(1+t^2) - lambda*atan(t),    y(0) = 0
#
# whose exact solution is y(t) = atan(t). CV_BDF is implicit, so CVODE solves a
# nonlinear system at every step; taking lambda strongly negative makes the
# problem stiff, which means those solves actually need several iterations rather
# than being finished by the predictor. Nothing in the solver below is specific to
# a scalar problem.
#
# HOW A CUSTOM NONLINEAR SOLVER FITS TOGETHER
#
# You subclass CustomSUNNonlinearSolver and override solve(), which is the only
# required method. Everything else is optional: override a method and SUNDIALS
# will call it, leave it alone and SUNDIALS behaves as though the operation were
# absent, exactly as for a solver written in C.
#
# CVODE does not hand your solve() the residual function directly. Instead it
# calls the setter methods -- set_sys_fn(), set_lsetup_fn(), set_lsolve_fn(),
# set_conv_test_fn() -- once, up front, passing an ordinary Python callable for
# each. You store those callables and use them inside solve(). The opaque
# integrator memory pointer that the C API threads through these callbacks is
# supplied for you, so the callables take only the arguments you care about.
#
# One rule follows from that: the callables given to set_sys_fn(),
# set_lsetup_fn() and set_lsolve_fn() are valid ONLY while SUNDIALS is inside
# your setup() or solve(). Calling one at any other time raises RuntimeError
# rather than reading a stale pointer. The callables from set_conv_test_fn() and
# the norm/rate setters carry their own data and may be called at any time.
# -----------------------------------------------------------------

import numpy as np
import pytest
from sundials4py.core import *
from sundials4py.cvodes import *

# Problem constants
NEQ = 1
LAMBDA = -100.0
T0 = 0.0
TF = 10.0
DTOUT = 1.0
RTOL = 1.0e-6
ATOL = 1.0e-10


class AnalyticODE:
    """dy/dt = lambda*y + 1/(1+t^2) - lambda*atan(t), with y(t) = atan(t)."""

    def __init__(self, lamb=LAMBDA):
        self.lamb = lamb

    def rhs(self, t, yvec, ydotvec, user_data):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)
        ydot[0] = self.lamb * y[0] + 1.0 / (1.0 + t * t) - self.lamb * np.arctan(t)
        return 0

    def jac(self, t, yvec, fyvec, J, user_data, tmp1, tmp2, tmp3):
        # df/dy = lambda for this scalar problem.
        SUNDenseMatrix_Data(J)[0, 0] = self.lamb
        return 0

    def solution(self, t):
        return np.arctan(t)


class MyNonlinearSolver(CustomSUNNonlinearSolver):
    """A SUNNonlinearSolver implemented in Python.

    Reference implementations worth reading alongside this template:
      src/sunnonlinsol/newton/sunnonlinsol_newton.c
      src/sunnonlinsol/fixedpoint/sunnonlinsol_fixedpoint.c
    """

    def __init__(self, template_vector, sunctx):
        # TODO(you): allocate whatever workspace your iteration needs. Cloning
        # the template vector is the usual way to get correctly sized scratch
        # space without assuming a particular N_Vector implementation.
        self.delta = N_VClone(template_vector)

        # Callables CVODE will install through the setters below. They start out
        # as None so that solve() can tell what the integrator actually provided.
        self.sys_fn = None
        self.lsetup_fn = None
        self.lsolve_fn = None
        self.conv_test_fn = None

        # Statistics SUNDIALS may ask for. CVODE reads cur_iter while deciding
        # whether to refresh the Jacobian, so keep it current.
        self.max_iters = 3
        self.cur_iter = 0
        self.num_iters = 0
        self.num_conv_fails = 0

        # ROOTFIND means "solve F(y) = 0"; FIXEDPOINT means "solve y = G(y)".
        # The choice tells CVODE which residual convention to hand you, so it
        # must match the iteration you write in solve().
        super().__init__(sunctx, SUNNONLINEARSOLVER_ROOTFIND)

    # -- optional: one-time setup ------------------------------------------

    def initialize(self):
        # Called once, after every callback has been set. Override this if your
        # iteration has state to reset before the first solve.
        return SUN_SUCCESS

    # -- the callback setters CVODE calls during CVodeSetNonlinearSolver ----

    def set_sys_fn(self, sys_fn):
        # sys_fn(y, F) -> status. Evaluates the nonlinear residual F at y.
        self.sys_fn = sys_fn
        return SUN_SUCCESS

    def set_lsetup_fn(self, lsetup_fn):
        # lsetup_fn(jbad) -> (status, jcur). Forms/refreshes the Newton matrix.
        # jbad says the caller believes the current one is stale; jcur reports
        # whether the returned one is freshly computed.
        self.lsetup_fn = lsetup_fn
        return SUN_SUCCESS

    def set_lsolve_fn(self, lsolve_fn):
        # lsolve_fn(b) -> status. Solves the Newton linear system in place: b
        # arrives as the right-hand side and is overwritten with the solution.
        self.lsolve_fn = lsolve_fn
        return SUN_SUCCESS

    def set_conv_test_fn(self, conv_test_fn):
        # conv_test_fn(y, delta, tol, w) -> SUN_SUCCESS when converged,
        # SUN_NLS_CONTINUE to keep iterating, or a failure code. Use the
        # integrator's test rather than inventing your own; it is what makes the
        # solver's accuracy consistent with the step size controller's.
        self.conv_test_fn = conv_test_fn
        return SUN_SUCCESS

    def set_max_iters(self, maxiters):
        self.max_iters = maxiters
        return SUN_SUCCESS

    # -- the required operation ---------------------------------------------

    def solve(self, y0, y, w, tol, call_lsetup):
        """Solve the nonlinear system.

        y0           predicted value, the initial iterate (do not modify)
        y            output: the solution (start by copying y0 into it)
        w            error weight vector, for the convergence test
        tol          convergence tolerance, for the convergence test
        call_lsetup  the integrator wants the Newton matrix refreshed first

        Return SUN_SUCCESS on convergence. Return SUN_NLS_CONV_RECVR for a
        failure the integrator can recover from by shrinking its step -- that is
        the right answer for "did not converge in max_iters". Return a negative
        code only for a failure no step size will fix.
        """
        N_VScale(1.0, y0, y)

        # TODO(you): write your iteration here.
        #
        # A modified-Newton iteration, for orientation:
        #
        #   if call_lsetup:
        #       status, jcur = self.lsetup_fn(True)
        #       if status != 0:
        #           return status
        #
        #   for self.cur_iter in range(self.max_iters):
        #       self.num_iters += 1
        #
        #       # residual at the current iterate
        #       status = self.sys_fn(y, self.delta)
        #       if status != 0:
        #           return status
        #
        #       # Newton update: solve J*delta = F in place
        #       status = self.lsolve_fn(self.delta)
        #       if status != 0:
        #           return status
        #
        #       N_VLinearSum(1.0, y, -1.0, self.delta, y)
        #
        #       status = self.conv_test_fn(y, self.delta, tol, w)
        #       if status == SUN_SUCCESS:
        #           return SUN_SUCCESS
        #       if status != SUN_NLS_CONTINUE:
        #           break
        #
        #   self.num_conv_fails += 1
        #   return SUN_NLS_CONV_RECVR
        #
        # Note that sys_fn, lsetup_fn and lsolve_fn are used only inside this
        # method. Stashing one and calling it later raises RuntimeError.
        raise NotImplementedError("TODO(you): implement MyNonlinearSolver.solve")

    # -- optional: statistics ----------------------------------------------
    #
    # Each of these returns a (status, value) pair, matching how sundials4py
    # binds C functions with output pointers.

    def get_num_iters(self):
        return SUN_SUCCESS, self.num_iters

    def get_cur_iter(self):
        return SUN_SUCCESS, self.cur_iter

    def get_num_conv_fails(self):
        return SUN_SUCCESS, self.num_conv_fails


def main():
    print("\nAnalytic ODE test problem:")
    print(f"   lambda = {LAMBDA}")
    print(f"   rtol = {RTOL}, atol = {ATOL}")
    print("Solution method: CVODE BDF with a Python nonlinear solver\n")

    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert status == SUN_SUCCESS

    problem = AnalyticODE()

    y = N_VNew_Serial(NEQ, sunctx)
    N_VGetArrayPointer(y)[0] = problem.solution(T0)

    cvode = CVodeCreate(CV_BDF, sunctx)
    assert CVodeInit(cvode.get(), problem.rhs, T0, y) == CV_SUCCESS
    assert CVodeSStolerances(cvode.get(), RTOL, ATOL) == CV_SUCCESS

    # The nonlinear solver needs a linear solver to form its Newton updates.
    # Here that is a native dense one, so the only Python-implemented object in
    # the loop is the nonlinear solver itself. Once your solver works, try
    # kin_custom_linsol.py's approach to replace this with your own.
    A = SUNDenseMatrix(NEQ, NEQ, sunctx)
    LS = SUNLinSol_Dense(y, A, sunctx)
    assert CVodeSetLinearSolver(cvode.get(), LS, A) == CV_SUCCESS
    assert CVodeSetJacFn(cvode.get(), problem.jac) == CV_SUCCESS

    # Attach the custom solver. `NLS` must stay referenced for as long as CVODE
    # uses it: the native handle holds only a weak reference back to the Python
    # object, so letting it go out of scope here would leave CVODE holding a
    # dangling pointer.
    NLS = MyNonlinearSolver(y, sunctx)
    assert CVodeSetNonlinearSolver(cvode.get(), NLS) == CV_SUCCESS

    print(f"{'t':>10}  {'y':>14}  {'error':>12}")
    print("-" * 40)

    t = T0
    tout = T0 + DTOUT
    while t < TF - 1.0e-12:
        status, t = CVode(cvode.get(), tout, y, CV_NORMAL)
        assert status == CV_SUCCESS, f"CVode returned {status}"

        computed = N_VGetArrayPointer(y)[0]
        exact = problem.solution(t)
        print(f"{t:10.4f}  {computed:14.8e}  {abs(computed - exact):12.4e}")

        tout = min(tout + DTOUT, TF)

    status, nst = CVodeGetNumSteps(cvode.get())
    assert status == CV_SUCCESS
    status, nfe = CVodeGetNumRhsEvals(cvode.get())
    assert status == CV_SUCCESS

    print("\nFinal Statistics..\n")
    print(f"nst      = {nst:6d}    nfe     = {nfe:6d}")
    print(f"nni      = {NLS.num_iters:6d}    ncfn    = {NLS.num_conv_fails:6d}")


def test_cvs_custom_nonlinsol():
    # This is a template, so there is nothing complete for CI to verify yet.
    # Delete this skip once you have filled in the TODO(you) sections.
    pytest.skip("template example: fill in the TODO(you) sections first")
    main()


if __name__ == "__main__":
    main()
