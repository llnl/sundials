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
# IMPLEMENTING A CUSTOM NONLINEAR SOLVER
#
# Create a subclass of CustomSUNNonlinearSolver and override solve(). All other
# functions are optional. CVODE provides the nonlinear residual function or
# fixed-point function as an ordinary Python callable through set_sys_fn(). Store
# this callable and use it inside solve(). It is valid only while SUNDIALS is
# inside solve(); calling it at any other time raises RuntimeError.
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
        y = N_VGetNumpyArray(yvec)
        ydot = N_VGetNumpyArray(ydotvec)
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
        # Scratch space for probing sys_fn at trial points.
        self.probe = N_VClone(template_vector)
        self.Fprobe = N_VClone(template_vector)

        # Callable CVODE installs through the setter below.
        self.sys_fn = None

        # ROOTFIND means "solve F(y) = 0", matching the correction equation
        # above.
        super().__init__(sunctx, SUNNONLINEARSOLVER_ROOTFIND)

    # -- the callback setters CVODE calls during CVodeSetNonlinearSolver ----

    def set_sys_fn(self, sys_fn):
        # sys_fn(y, F) -> status. Evaluates the nonlinear residual F at y.
        self.sys_fn = sys_fn
        return SUN_SUCCESS

    # -- helper: evaluate the scalar residual at a trial point ---------------

    def _residual(self, yval):
        N_VGetNumpyArray(self.probe)[0] = yval
        status = self.sys_fn(self.probe, self.Fprobe)
        return status, N_VGetNumpyArray(self.Fprobe)[0]

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
        # Probe F at three points straddling the current correction guess c0
        # (CVODE always starts this at 0). Since F(c0+u) is exactly
        # A*u^2 + B*u + C for some A, B, C, these three values pin down all
        # three coefficients exactly -- this is algebra, not a
        # finite-difference approximation, so the probe spacing only needs to
        # be chosen for floating-point conditioning, not accuracy. Scale it to
        # the predicted solution's own magnitude (y0), since the correction
        # itself starts at 0 regardless of how large y is.
        c0 = N_VGetNumpyArray(y)[0]
        delta = max(abs(N_VGetNumpyArray(y0)[0]), 1.0)
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

        N_VGetNumpyArray(y)[0] = c0 + u
        return SUN_SUCCESS


def main():
    print("\nLogistic ODE test problem:")
    print(f"   r = {R}, K = {K}, y0 = {Y0}")
    print(f"   rtol = {RTOL}, atol = {ATOL}")
    print("Solution method: CVODE BDF with a Python nonlinear solver\n")

    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert status == SUN_SUCCESS

    problem = LogisticODE()

    y = N_VNew_Serial(NEQ, sunctx)
    N_VGetNumpyArray(y)[0] = problem.solution(T0)

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

        computed = N_VGetNumpyArray(y)[0]
        exact = problem.solution(t)
        print(f"{t:10.4f}  {computed:14.8e}  {abs(computed - exact):12.4e}")

        tout = min(tout + DTOUT, TF)

    status, nst = CVodeGetNumSteps(cvode.get())
    assert status == CV_SUCCESS
    status, nfe = CVodeGetNumRhsEvals(cvode.get())
    assert status == CV_SUCCESS

    print("\nFinal Statistics..\n")
    print(f"nst      = {nst:6d}    nfe     = {nfe:6d}")


def test_cvs_custom_nonlinsol():
    main()


if __name__ == "__main__":
    main()
