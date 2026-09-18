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
# This example is based on examples/python/arkode/ark_heat1D.py,
# but it uses ARKODE's LSRKStep module with a super-time-stepping
# (STS) method for time evolution.
#
# The following test simulates a simple 1D heat equation,
#    u_t = k*u_xx + f
# for t in [0, 1], x in [0, 1], with initial conditions
#    u(0,x) =  0
# Dirichlet boundary conditions, i.e.
#    u_t(t,0) = u_t(t,1) = 0,
# and a point-source heating term,
#    f = 0.01 for x=0.5.
#
# The spatial derivatives are computed using second-order centered
# differences, with the data distributed over N points on a uniform
# spatial grid.
#
# LSRKStep's STS methods are designed for explicit integration of
# problems with dominant negative-real-axis eigenvalues, such as
# diffusion operators.  The example therefore provides a simple
# dominant eigenvalue estimate for the finite-difference Laplacian.
# -----------------------------------------------------------------

import numpy as np
import sundials4py.arkode as ark
import sundials4py.core as sun


def exact_semidiscrete_solution(n, k, t):
    dx = 1.0 / (n - 1)
    i = np.arange(1, n - 1)
    m = np.arange(1, n - 1)

    phi = np.sqrt(2.0 / (n - 1)) * np.sin(np.outer(i, m) * np.pi / (n - 1))
    lambdas = -4.0 * k / dx**2 * np.sin(0.5 * m * np.pi / (n - 1)) ** 2

    source = np.zeros(n - 2)
    source[(n // 2) - 1] = 0.01 / dx

    u = np.zeros(n)
    u[1:-1] = phi @ (((np.exp(lambdas * t) - 1.0) / lambdas) * (phi.T @ source))
    return u


class Heat1DProblem:
    def __init__(self, n=101, k=0.01):
        self.n = n
        self.k = k
        self.dx = 1.0 / (n - 1)
        self.isource = n // 2

    def set_init_cond(self, yvec):
        y = sun.N_VGetNumpyArray(yvec)
        y[:] = 0.0

    def f(self, t, yvec, ydotvec, user_data):
        y = sun.N_VGetNumpyArray(yvec)
        ydot = sun.N_VGetNumpyArray(ydotvec)

        # Initialize the RHS to zero.  This also imposes stationary
        # boundary conditions because the endpoints are never overwritten.
        ydot[:] = 0.0

        c1 = self.k / self.dx / self.dx
        c2 = -2.0 * self.k / self.dx / self.dx
        ydot[1:-1] = c1 * y[:-2] + c2 * y[1:-1] + c1 * y[2:]
        ydot[self.isource] += 0.01 / self.dx
        return 0

    def dom_eig(self, t, yvec, fnvec, user_data, tempv1, tempv2, tempv3):
        # The most negative eigenvalue of the 1D centered-difference
        # diffusion operator is bounded in magnitude by 4*k/dx^2.
        lambdaR = -4.0 * self.k / (self.dx * self.dx)
        lambdaI = 0.0
        return 0, lambdaR, lambdaI


def solve_heat1d_lsrk():
    n = 101
    k = 0.01
    tf = 1.0
    nt = 10
    reltol = 1e-6
    abstol = 1e-10

    status, sunctx = sun.SUNContext_Create(sun.SUN_COMM_NULL)
    assert status == sun.SUN_SUCCESS

    y = sun.N_VNew_Serial(n, sunctx)
    assert y is not None
    yarr = sun.N_VGetNumpyArray(y)

    problem = Heat1DProblem(n=n, k=k)
    problem.set_init_cond(y)

    # Create an explicit LSRKStep integrator using an STS method.  Unlike
    # ARKStep, this does not need a nonlinear solver or linear solver.
    stepper = ark.LSRKStepCreateSTS(problem.f, 0.0, y, sunctx)
    assert stepper is not None

    status = ark.ARKodeSStolerances(stepper.get(), reltol, abstol)
    assert status == ark.ARK_SUCCESS

    status = ark.ARKodeSetMaxNumSteps(stepper.get(), 100000)
    assert status == ark.ARK_SUCCESS

    status = ark.LSRKStepSetSTSMethod(stepper.get(), ark.ARKODE_LSRK_RKC_2)
    assert status == ark.ARK_SUCCESS

    status = ark.LSRKStepSetDomEigFn(stepper.get(), problem.dom_eig)
    assert status == ark.ARK_SUCCESS

    # These optional settings follow the C++ LSRK STS example: reuse the
    # eigenvalue estimate for several steps, allow many internal stages,
    # and add a small safety factor to the eigenvalue bound.
    status = ark.LSRKStepSetDomEigFrequency(stepper.get(), 25)
    assert status == ark.ARK_SUCCESS

    status = ark.LSRKStepSetMaxNumStages(stepper.get(), 1000)
    assert status == ark.ARK_SUCCESS

    status = ark.LSRKStepSetDomEigSafetyFactor(stepper.get(), 1.01)
    assert status == ark.ARK_SUCCESS

    status = ark.ARKodeSetAdaptControllerByName(stepper.get(), "I")
    assert status == ark.ARK_SUCCESS

    status = ark.ARKodeSetStopTime(stepper.get(), tf)
    assert status == ark.ARK_SUCCESS

    t = 0.0
    tout = tf / nt

    print("\n1D heat equation test problem (LSRKStep STS):")
    print(f"    N = {n},  k = {k}")
    print(f"    reltol = {reltol:.1e},  abstol = {abstol:.1e}\n")
    print("        t      ||u||_rms")
    print("   -------------------------")
    print(f"  {t:10.6f}  {0.0:10.6f}")

    for _ in range(nt):
        status, t = ark.ARKodeEvolve(stepper.get(), tout, y, ark.ARK_NORMAL)
        if status < 0:
            raise RuntimeError(f"ARKodeEvolve failed with status {status}")

        print(f"  {t:10.6f}  {np.sqrt(np.dot(yarr, yarr) / n):10.6f}")
        tout = min(tout + tf / nt, tf)

    print("   -------------------------")
    uexact = exact_semidiscrete_solution(n, k, tf)
    max_error = np.max(np.abs(yarr - uexact))
    print(f"\nFinal max error vs exact semi-discrete solution = {max_error:.6e}")
    np.testing.assert_allclose(yarr, uexact, rtol=1e-4, atol=1e-8)

    status, nst = ark.ARKodeGetNumSteps(stepper.get())
    assert status == ark.ARK_SUCCESS
    status, nst_a = ark.ARKodeGetNumStepAttempts(stepper.get())
    assert status == ark.ARK_SUCCESS
    status, nfe = ark.ARKodeGetNumRhsEvals(stepper.get(), 0)
    assert status == ark.ARK_SUCCESS
    status, netf = ark.ARKodeGetNumErrTestFails(stepper.get())
    assert status == ark.ARK_SUCCESS
    status, ndomeig = ark.LSRKStepGetNumDomEigUpdates(stepper.get())
    assert status == ark.ARK_SUCCESS

    print("\nFinal Solver Statistics:")
    print(f"   Internal solver steps = {nst} (attempted = {nst_a})")
    print(f"   Total RHS evals = {nfe}")
    print(f"   Dominant eigenvalue updates = {ndomeig}")
    print(f"   Total number of error test failures = {netf}")

    return yarr.copy()


def main():
    solve_heat1d_lsrk()


# This function allows pytest to discover the example as a test.
def test_ark_heat1D_lsrk():
    main()


if __name__ == "__main__":
    main()
