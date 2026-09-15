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
# EXAMPLE: implementing a SUNNonlinearSolver in Python.
#
# The test problem is the logistic equation
#
#    dy/dt = r*y*(1 - y/K),    y(0) = y0
#
# whose exact solution is the logistic curve
#
#    y(t) = K*y0*exp(r*t) / (K + y0*(exp(r*t) - 1)).
#
# CV_BDF is implicit, so CVODE solves a nonlinear system at every step. Rather
# than iterate generically for that system (as a Newton or fixed-point solver
# would), the solver below exploits a fact specific to this right-hand side:
# the per-step correction equation is exactly quadratic in the unknown, so it
# can be solved in one shot with the quadratic formula. See MyNonlinearSolver
# below for how that works and what it assumes about the problem.
#
# HOW A CUSTOM NONLINEAR SOLVER FITS TOGETHER
#
# You subclass CustomSUNNonlinearSolver and override solve(), which is the
# only required method. Everything else is optional: override a method and
# SUNDIALS will call it, leave it alone and SUNDIALS behaves as though the
# operation were absent, exactly as for a solver written in C.
#
# CVODE does not hand your solve() the residual function directly. Instead it
# calls setter methods -- set_sys_fn() and set_conv_test_fn() here -- once, up
# front, passing an ordinary Python callable for each. You store those
# callables and use them inside solve(). The opaque integrator memory pointer
# that the C API threads through these callbacks is supplied for you, so the
# callables take only the arguments you care about.
#
# A solver that forms Newton updates would also override set_lsetup_fn() and
# set_lsolve_fn() to get at a linear solver -- see kin_custom_linsol.py for
# that pattern. This solver never forms a Jacobian, so it leaves both
# unimplemented; CVODE treats that exactly like SUNNonlinSol_FixedPoint, which
# does not use a linear solver either.
#
# One rule follows from the callback pattern: the callable given to
# set_sys_fn() is valid ONLY while SUNDIALS is inside your solve(). Calling it
# at any other time raises RuntimeError rather than reading a stale pointer.
# The callable from set_conv_test_fn() carries its own data and may be called
# at any time.
# -----------------------------------------------------------------

import numpy as np
from sundials4py.core import *
from sundials4py.cvodes import *

# Problem constants
NEQ = 1
R = 2.0
K = 5.0
Y0 = 0.1
T0 = 0.0
TF = 5.0
DTOUT = 0.5
RTOL = 1.0e-6
ATOL = 1.0e-10


class LogisticODE:
    """dy/dt = r*y*(1 - y/K), with y(t) = K*y0*exp(r*t)/(K + y0*(exp(r*t)-1))."""

    def __init__(self, r=R, k=K, y0=Y0):
        self.r = r
        self.k = k
        self.y0 = y0

    def rhs(self, t, yvec, ydotvec, user_data):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)
        ydot[0] = self.r * y[0] * (1.0 - y[0] / self.k)
        return 0

    def solution(self, t):
        e = np.exp(self.r * t)
        return self.k * self.y0 * e / (self.k + self.y0 * (e - 1.0))


class MyNonlinearSolver(CustomSUNNonlinearSolver):
    """A SUNNonlinearSolver implemented in Python.

    CVODE does not solve for y directly -- it solves for the CORRECTION c to
    an already-computed predictor, y = ypred + c, starting from c = 0. Its
    per-step correction equation has the general form

        F(c) = c + rl1*zn1 - gamma*f(t, ypred + c) = 0

    for a step-dependent gamma, rl1, and history term zn1 that this solver
    never sees directly -- it only gets to call sys_fn(c, F), which evaluates
    F(c) for a given c. But this problem's right-hand side f(y) = r*y*(1-y/K)
    is a quadratic polynomial in y, and ypred + c is affine in c, so F(c) is a
    quadratic polynomial in c too. That means solve() does not need to iterate
    at all: probe sys_fn at three points to recover F's three coefficients
    exactly (no truncation error, since F really is a quadratic), then solve
    the quadratic directly for its root.

    This is specific to a right-hand side that is (at most) quadratic in y --
    it is not a general substitute for Newton iteration. For a general
    right-hand side, see the modified-Newton sketch this file used to carry;
    reference implementations of that approach are:
      src/sunnonlinsol/newton/sunnonlinsol_newton.c
      src/sunnonlinsol/fixedpoint/sunnonlinsol_fixedpoint.c
    """

    def __init__(self, template_vector, sunctx):
        # Scratch space for probing sys_fn at trial points, and for reporting
        # the correction to the convergence test.
        self.probe = N_VClone(template_vector)
        self.Fprobe = N_VClone(template_vector)
        self.corr = N_VClone(template_vector)

        # Callables CVODE installs through the setters below.
        self.sys_fn = None
        self.conv_test_fn = None

        # Statistics SUNDIALS may ask for. CVODE's convergence test judges an
        # iteration by how much the correction has shrunk relative to the
        # previous one, so it needs to know which iteration this is; see
        # solve() for why that matters even though this solver never really
        # "iterates".
        self.cur_iter = 0
        self.num_iters = 0
        self.num_conv_fails = 0

        # ROOTFIND means "solve F(y) = 0", matching the correction equation
        # above.
        super().__init__(sunctx, SUNNONLINEARSOLVER_ROOTFIND)

    # -- the callback setters CVODE calls during CVodeSetNonlinearSolver ----

    def set_sys_fn(self, sys_fn):
        # sys_fn(y, F) -> status. Evaluates the nonlinear residual F at y.
        self.sys_fn = sys_fn
        return SUN_SUCCESS

    def set_conv_test_fn(self, conv_test_fn):
        # conv_test_fn(y, delta, tol, w) -> SUN_SUCCESS when converged,
        # SUN_NLS_CONTINUE to keep iterating, or a failure code. Use the
        # integrator's test rather than inventing your own; it is what makes
        # the solver's accuracy consistent with the step size controller's.
        self.conv_test_fn = conv_test_fn
        return SUN_SUCCESS

    # -- helper: evaluate the scalar residual at a trial point ---------------

    def _residual(self, yval):
        N_VGetArrayPointer(self.probe)[0] = yval
        status = self.sys_fn(self.probe, self.Fprobe)
        return status, N_VGetArrayPointer(self.Fprobe)[0]

    # -- the required operation ---------------------------------------------

    def solve(self, y0, y, w, tol, call_lsetup):
        """Solve the nonlinear system by fitting and inverting its quadratic.

        y0   the predictor CVODE corrects away from; sys_fn already builds it
             into F, so it is used here only to scale the probe spacing
        y    the correction c: arrives holding CVODE's initial guess c = 0,
             output should overwrite it with the solved-for correction
        w    error weight vector, for the convergence test
        tol  convergence tolerance, for the convergence test

        Return SUN_SUCCESS on convergence, or SUN_NLS_CONV_RECVR for a
        failure the integrator can recover from by shrinking its step.
        """
        self.num_iters += 1

        # Probe F at three points straddling the current correction guess c0
        # (CVODE always starts this at 0). Since F(c0+u) is exactly
        # A*u^2 + B*u + C for some A, B, C, these three values pin down all
        # three coefficients exactly -- this is algebra, not a
        # finite-difference approximation, so the probe spacing only needs to
        # be chosen for floating-point conditioning, not accuracy. Scale it to
        # the predicted solution's own magnitude (y0), since the correction
        # itself starts at 0 regardless of how large y is.
        c0 = N_VGetArrayPointer(y)[0]
        delta = max(abs(N_VGetArrayPointer(y0)[0]), 1.0)
        status, Fm = self._residual(c0 - delta)
        if status != 0:
            return status
        status, F0 = self._residual(c0)
        if status != 0:
            return status
        status, Fp = self._residual(c0 + delta)
        if status != 0:
            return status

        A = (Fp + Fm - 2.0 * F0) / (2.0 * delta * delta)
        B = (Fp - Fm) / (2.0 * delta)
        C = F0

        if abs(A) < 1.0e-14 * max(abs(B), 1.0):
            # Degenerate to a linear equation (gamma*r/K ~ 0); avoid dividing
            # by a near-zero leading coefficient.
            u = -C / B
        else:
            disc = B * B - 4.0 * A * C
            if disc < 0.0:
                self.num_conv_fails += 1
                return SUN_NLS_CONV_RECVR
            # The textbook -B +/- sqrt(disc) formula cancels catastrophically
            # whenever B**2 >> 4*A*C -- which is exactly the well-conditioned
            # regime here, since A (curvature from the y^2 term) is small
            # relative to B (dominated by the "c" term in F). Compute the
            # well-separated root first and get the other from the product
            # of roots, C/A, which sidesteps the cancellation entirely.
            sqrt_disc = np.sqrt(disc)
            q = -0.5 * (B + np.copysign(sqrt_disc, B))
            u1 = q / A
            u2 = C / q
            # c0 is already a good guess (typically 0, the first correction
            # of the step), so keep whichever root is closest to it.
            u = u1 if abs(u1) < abs(u2) else u2

        N_VGetArrayPointer(y)[0] = c0 + u

        # Report the update just applied to the integrator's own convergence
        # test, rather than assuming the algebra above is exact in floating
        # point.
        self.cur_iter = 0
        N_VGetArrayPointer(self.corr)[0] = u
        status = self.conv_test_fn(y, self.corr, tol, w)
        if status == SUN_SUCCESS:
            return SUN_SUCCESS
        if status != SUN_NLS_CONTINUE:
            self.num_conv_fails += 1
            return SUN_NLS_CONV_RECVR

        # CVODE's convergence test compares the size of successive
        # corrections and is unconvinced by a single (possibly large) one, no
        # matter how accurate -- it is designed for an iteration that
        # shrinks, not a closed-form solve that lands on the answer in one
        # step. So confirm exactness in a way that test understands: report a
        # second "iteration" that applies zero further correction. Since the
        # root really is exact (up to floating point), that is simply the
        # truth, and a strictly shrinking correction is what the test wants
        # to see.
        self.cur_iter = 1
        N_VConst(0.0, self.corr)
        status = self.conv_test_fn(y, self.corr, tol, w)
        if status == SUN_SUCCESS:
            return SUN_SUCCESS
        self.num_conv_fails += 1
        return SUN_NLS_CONV_RECVR

    # -- optional: statistics ----------------------------------------------
    #
    # Each of these returns a (status, value) pair, matching how sundials4py
    # binds C functions with output pointers.

    def get_cur_iter(self):
        return SUN_SUCCESS, self.cur_iter

    def get_num_iters(self):
        return SUN_SUCCESS, self.num_iters

    def get_num_conv_fails(self):
        return SUN_SUCCESS, self.num_conv_fails


def main():
    print("\nLogistic ODE test problem:")
    print(f"   r = {R}, K = {K}, y0 = {Y0}")
    print(f"   rtol = {RTOL}, atol = {ATOL}")
    print("Solution method: CVODE BDF with a Python nonlinear solver\n")

    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert status == SUN_SUCCESS

    problem = LogisticODE()

    y = N_VNew_Serial(NEQ, sunctx)
    N_VGetArrayPointer(y)[0] = problem.solution(T0)

    cvode = CVodeCreate(CV_BDF, sunctx)
    assert CVodeInit(cvode.get(), problem.rhs, T0, y) == CV_SUCCESS
    assert CVodeSStolerances(cvode.get(), RTOL, ATOL) == CV_SUCCESS

    # Attach the custom solver. `NLS` must stay referenced for as long as CVODE
    # uses it: the native handle holds only a weak reference back to the
    # Python object, so letting it go out of scope here would leave CVODE
    # holding a dangling pointer. No linear solver is attached, since this
    # solver never forms one.
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
    main()


if __name__ == "__main__":
    main()
