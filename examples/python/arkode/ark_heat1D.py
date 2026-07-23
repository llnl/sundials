#!/usr/bin/env python3
# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos
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
# This example is a copy of examples/arkode/C_serial/ark_heat1D.c,
# but ported to Python to use sundials4py.
#
# The following test simulates a simple 1D heat equation,
#    u_t = k*u_xx + f
# for t in [0, 10], x in [0, 1], with initial conditions
#    u(0,x) =  0
# Dirichlet boundary conditions, i.e.
#    u_t(t,0) = u_t(t,1) = 0,
# and a point-source heating term,
#    f = 0.01 for x=0.5.
#
# The spatial derivatives are computed using second-order
# centered differences, with the data distributed over N points
# on a uniform spatial grid.
#
# The final solution is checked against the exact solution of this
# semi-discrete ODE system.  With zero Dirichlet boundaries, the
# interior finite-difference Laplacian is diagonalized by the discrete
# sine basis.  For zero initial data and a constant point source b, the
# interior solution is
#
#    u(t) = Phi * diag((exp(lambda_m*t) - 1) / lambda_m) * Phi^T * b,
#
# where Phi contains the orthonormal sine modes and lambda_m are the
# corresponding finite-difference Laplacian eigenvalues.  This validates
# the time integration error without introducing error from a continuous
# PDE approximation.
#
# This program solves the problem with either an ERK or DIRK
# method.  For the DIRK method, we use a Newton iteration with
# the SUNLinSol_PCG linear solver, and a user-supplied Jacobian-vector
# product routine.
#
# 100 outputs are printed at equal intervals, and run statistics
# are printed at the end.
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
        y = sun.N_VGetArrayPointer(yvec)
        y[:] = 0.0

    def f(self, t, yvec, ydotvec, user_data):
        y = sun.N_VGetArrayPointer(yvec)
        ydot = sun.N_VGetArrayPointer(ydotvec)

        ydot[:] = 0.0
        c1 = self.k / self.dx / self.dx
        c2 = -2.0 * self.k / self.dx / self.dx
        ydot[1:-1] = c1 * y[:-2] + c2 * y[1:-1] + c1 * y[2:]
        ydot[0] = 0.0
        ydot[-1] = 0.0
        ydot[self.isource] += 0.01 / self.dx
        return 0


def solve_heat1d():
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
    yarr = sun.N_VGetArrayPointer(y)

    problem = Heat1DProblem(n=n, k=k)
    problem.set_init_cond(y)

    stepper = ark.ARKStepCreate(problem.f, None, 0.0, y, sunctx)
    assert stepper is not None

    status = ark.ARKodeSStolerances(stepper.get(), reltol, abstol)
    assert status == ark.ARK_SUCCESS

    status = ark.ARKodeSetMaxNumSteps(stepper.get(), 100000)
    assert status == ark.ARK_SUCCESS

    t = 0.0
    tout = tf / nt

    print("        t      ||u||_rms")
    print("   -------------------------")
    print(f"  {t:10.6f}  {0.0:10.6f}")

    for _ in range(nt):
        status, t = ark.ARKodeEvolve(stepper.get(), tout, y, ark.ARK_NORMAL)
        if status != ark.ARK_SUCCESS:
            raise RuntimeError(f"ARKodeEvolve failed with status {status}")

        print(f"  {t:10.6f}  {np.sqrt(np.dot(yarr, yarr) / n):10.6f}")
        tout = min(tout + tf / nt, tf)

    uexact = exact_semidiscrete_solution(n, k, tf)
    max_error = np.max(np.abs(yarr - uexact))
    print(f"\nFinal max error vs exact semi-discrete solution = {max_error:.6e}")
    np.testing.assert_allclose(yarr, uexact, rtol=1e-4, atol=1e-8)

    return yarr.copy()


def main():
    solve_heat1d()


# This function allows pytest to discover the example as a test
def test_ark_heat1D():
    main()


if __name__ == "__main__":
    main()
