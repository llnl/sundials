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
# EXAMPLE: implementing a SUNAdaptController in Python.
#
# The test problem is the scalar ODE
#
#    dy/dt = (t+1)*exp(-y),    y(0) = 0
#
# whose exact solution is y(t) = log(t^2/2 + t + 1). It is non-stiff, so it runs
# with an explicit Runge-Kutta method. The solution's curvature falls off as t
# grows, so a controller should increase the step sizes over the run.
#
# WHAT A CONTROLLER IS ASKED TO DO
#
# After every step attempt, ARKODE forms a scaled error measure `dsm` from the
# embedded error estimate: dsm ~ 1 means the step just met the requested
# tolerance, dsm > 1 means it failed the error test and will be retried, and
# dsm < 1 means it was more accurate than needed. ARKODE passes this estimate to
# the controller which returns the step size to try next.
#
# Create a subclass of CustomSUNHController and override estimate_step(). All
# other functions are optional.
#
# For a multirate (MRIStep) controller, subclass CustomSUNMRIController instead
# and override estimate_step_tol(H, tolfac, P, DSM, dsm), which returns both a
# slow step size and a fast-solver tolerance factor.
# -----------------------------------------------------------------

import numpy as np
from sundials4py.arkode import *
from sundials4py.core import *

# Problem constants
NEQ = 1
T0 = 0.0
TF = 10.0
NOUT = 10
RTOL = 1.0e-6
ATOL = 1.0e-10


class NonlinearODE:
    """dy/dt = (t+1)*exp(-y), with y(t) = log(t^2/2 + t + 1)."""

    def rhs(self, t, yvec, ydotvec, user_data):
        y = N_VGetNumpyArray(yvec)
        ydot = N_VGetNumpyArray(ydotvec)
        ydot[0] = (t + 1.0) * np.exp(-y[0])
        return 0

    def solution(self, t):
        return np.log(0.5 * t * t + t + 1.0)


class MyController(CustomSUNHController):
    """A SUNAdaptController implemented in Python."""

    def __init__(self, sunctx, safety=0.9):
        # Step size safety and error bias factors for an I-controller.
        self.safety = safety
        self.bias = 1.0

        # The base constructor takes only the context, and must be called after
        # your own state is in place: it makes the object convertible to a native
        # handle, so the operations below have to be ready to run.
        super().__init__(sunctx)

    # -- the required operation ---------------------------------------------

    def estimate_step(self, h, p, dsm):
        """Return the step size to attempt next.

        h    the step size just attempted
        p    the order of the error estimate driving the adaptivity
        dsm  scaled error measure for that attempt: <1 accurate, >1 rejected

        Return a (status, hnew) pair -- SUN_SUCCESS and the new step size --
        matching how sundials4py binds C functions with output pointers. ARKODE
        applies its own step size bounds and change ratio limits to whatever you
        return, so the controller does not need to enforce them itself.
        """
        # Guard the error measure away from zero: an exact step would
        # otherwise ask for an infinite increase.
        e = max(self.bias * dsm, 1.0e-10)

        # Asymptotically the error scales like h^(p+1), so scaling h by
        # e^(-1/(p+1)) targets dsm == 1. The safety factor keeps the next
        # attempt on the accurate side of that target.
        hnew = self.safety * h * e ** (-1.0 / (p + 1))
        return SUN_SUCCESS, hnew

    # -- optional operations -------------------------------------------------

    def set_defaults(self):
        # Restore the parameters to their default values, undoing any tuning.
        self.bias = 1.0
        return SUN_SUCCESS

    def set_error_bias(self, bias):
        # An error multiplier letting a user aim below the requested tolerance.
        # Store it and use it in estimate_step(); ARKODE will not apply it for you.
        self.bias = bias
        return SUN_SUCCESS


def main():
    print("\nNonlinear ODE test problem:")
    print(f"   rtol = {RTOL}, atol = {ATOL}")
    print("Solution method: ARKODE ERK with a Python time step controller\n")

    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert status == SUN_SUCCESS

    problem = NonlinearODE()

    y = N_VNew_Serial(NEQ, sunctx)
    N_VGetNumpyArray(y)[0] = problem.solution(T0)

    ark = ERKStepCreate(problem.rhs, T0, y, sunctx)
    assert ARKodeSStolerances(ark.get(), RTOL, ATOL) == ARK_SUCCESS

    # Attach the custom controller. `controller` must stay referenced for as long
    # as ARKODE uses it: the native handle holds only a weak reference back to
    # the Python object, so letting it go out of scope here would leave ARKODE
    # holding a dangling pointer.
    controller = MyController(sunctx)
    assert ARKodeSetAdaptController(ark.get(), controller) == ARK_SUCCESS

    print(f"{'t':>10}  {'y':>14}  {'error':>12}  {'last h':>12}")
    print("-" * 54)

    dtout = (TF - T0) / NOUT
    tout = T0 + dtout
    for _ in range(NOUT):
        status, t = ARKodeEvolve(ark.get(), tout, y, ARK_NORMAL)
        assert status == ARK_SUCCESS, f"ARKodeEvolve returned {status}"

        status, hlast = ARKodeGetLastStep(ark.get())
        assert status == ARK_SUCCESS

        computed = N_VGetNumpyArray(y)[0]
        exact = problem.solution(t)
        print(f"{t:10.4f}  {computed:14.8e}  {abs(computed - exact):12.4e}  {hlast:12.4e}")

        tout = min(tout + dtout, TF)

    status, nst = ARKodeGetNumSteps(ark.get())
    assert status == ARK_SUCCESS
    status, nst_a = ARKodeGetNumStepAttempts(ark.get())
    assert status == ARK_SUCCESS
    status, nfe = ARKodeGetNumRhsEvals(ark.get(), 0)
    assert status == ARK_SUCCESS
    status, netf = ARKodeGetNumErrTestFails(ark.get())
    assert status == ARK_SUCCESS

    print("\nFinal Statistics..\n")
    print(f"nst      = {nst:6d}    nst_a   = {nst_a:6d}")
    print(f"nfe      = {nfe:6d}    netf    = {netf:6d}")


def test_ark_custom_adaptcontroller():
    main()


if __name__ == "__main__":
    main()
