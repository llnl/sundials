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
# This program solves the problem with either an ERK or DIRK
# method.  For the DIRK method, we use a Newton iteration with
# the SUNLinSol_PCG linear solver, and a user-supplied Jacobian-vector
# product routine.
#
# 100 outputs are printed at equal intervals, and run statistics
# are printed at the end.
# -----------------------------------------------------------------

import numpy as np
from sundials4py.core import *
from sundials4py.arkode import *


class Heat1DProblem:
    def __init__(self, N, k):
        self.N = N
        self.k = k
        self.dx = 1.0 / (N - 1)

    def set_init_cond(self, yvec):
        y = N_VGetArrayPointer(yvec)
        y[:] = 0.0
        return 0

    def f(self, t, yvec, ydotvec, user_data):
        N, k, dx = self.N, self.k, self.dx
        Y = N_VGetArrayPointer(yvec)
        Ydot = N_VGetArrayPointer(ydotvec)
        Ydot[:] = 0.0
        c1 = k / dx / dx
        c2 = -2.0 * k / dx / dx
        # Vectorized Laplacian
        Ydot[1:-1] = c1 * Y[:-2] + c2 * Y[1:-1] + c1 * Y[2:]
        # Dirichlet BCs
        Ydot[0] = 0.0
        Ydot[-1] = 0.0
        # Point source
        isource = N // 2
        Ydot[isource] += 0.01 / dx
        return 0

    def jtv(self, vvec, Jvvec, t, yvec, fyvec, tmpvec, user_data):
        N, k, dx = self.N, self.k, self.dx
        V = N_VGetArrayPointer(vvec)
        JV = N_VGetArrayPointer(Jvvec)
        JV[:] = 0.0
        c1 = k / dx / dx
        c2 = -2.0 * k / dx / dx
        # Vectorized tridiagonal product
        JV[1:-1] = c1 * V[:-2] + c2 * V[1:-1] + c1 * V[2:]
        JV[0] = 0.0
        JV[-1] = 0.0
        return 0


def main():
    # Problem parameters
    N = 201
    k = 0.5
    T0 = 0.0
    Tf = 1.0
    Nt = 10
    reltol = 1e-6
    abstol = 1e-10

    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert status == SUN_SUCCESS
    assert sunctx is not None

    y = N_VNew_Serial(N, sunctx)
    assert y is not None

    problem = Heat1DProblem(N, k)
    problem.set_init_cond(y)

    # Call ARKStepCreate to initialize the ARK timestepper module and
    # specify the right-hand side function in y'=f(t,y), the initial time
    # T0, and the initial dependent variable vector y.  Note: since this
    # problem is fully implicit, we set f_E to None and f_I to f. */
    ark = ARKStepCreate(None, problem.f, T0, y, sunctx)  # f_E (explicit)
    assert ark is not None

    # Set routines
    status = ARKodeSStolerances(ark.get(), reltol, abstol)
    assert status == ARK_SUCCESS

    status = ARKodeSetMaxNumSteps(ark.get(), 10000)
    assert status == ARK_SUCCESS

    status = ARKodeSetPredictorMethod(ark.get(), 1)
    assert status == ARK_SUCCESS

    # PCG linear solver with no preconditioning, with up to N iterations
    LS = SUNLinSol_PCG(y, SUN_PREC_NONE, N, sunctx)
    status = ARKodeSetLinearSolver(ark.get(), LS, None)
    assert status == ARK_SUCCESS

    status = ARKodeSetJacTimes(ark.get(), None, problem.jtv)
    assert status == ARK_SUCCESS

    status = ARKodeSetLinear(ark.get(), 0)
    assert status == ARK_SUCCESS

    # Output mesh
    with open("heat_mesh.txt", "w") as FID:
        for i in range(N):
            FID.write(f"  {problem.dx * i:.16e}\n")

    yarr = N_VGetArrayPointer(y)
    with open("heat1D.txt", "w") as UFID:
        # Output initial condition
        UFID.write(" ".join(f"{val:.16e}" for val in yarr) + "\n")

        # Main time-stepping loop
        t = T0
        dTout = (Tf - T0) / Nt
        tout = T0 + dTout
        print("        t      ||u||_rms")
        print("   -------------------------")
        print(f"  {t:10.6f}  {np.sqrt(np.dot(yarr, yarr) / N):10.6f}")
        for iout in range(Nt):
            status, tret = ARKodeEvolve(ark.get(), tout, y, ARK_NORMAL)
            yarr = N_VGetArrayPointer(y)
            print(f"  {tret:10.6f}  {np.sqrt(np.dot(yarr, yarr) / N):10.6f}")
            UFID.write(" ".join(f"{val:.16e}" for val in yarr) + "\n")
            if status == ARK_SUCCESS:
                tout += dTout
                tout = min(tout, Tf)
            else:
                print("Solver failure, stopping integration")
                break
        print("   -------------------------")

    # Print statistics
    status, nst = ARKodeGetNumSteps(ark.get())
    assert status == ARK_SUCCESS
    status, nst_a = ARKodeGetNumStepAttempts(ark.get())
    assert status == ARK_SUCCESS
    status, nfe = ARKodeGetNumRhsEvals(ark.get(), 0)
    assert status == ARK_SUCCESS
    status, nfi = ARKodeGetNumRhsEvals(ark.get(), 1)
    assert status == ARK_SUCCESS
    status, nsetups = ARKodeGetNumLinSolvSetups(ark.get())
    assert status == ARK_SUCCESS
    status, nli = ARKodeGetNumLinIters(ark.get())
    assert status == ARK_SUCCESS
    status, nJv = ARKodeGetNumJtimesEvals(ark.get())
    assert status == ARK_SUCCESS
    status, nlcf = ARKodeGetNumLinConvFails(ark.get())
    assert status == ARK_SUCCESS
    status, nni = ARKodeGetNumNonlinSolvIters(ark.get())
    assert status == ARK_SUCCESS
    status, ncfn = ARKodeGetNumNonlinSolvConvFails(ark.get())
    assert status == ARK_SUCCESS
    status, netf = ARKodeGetNumErrTestFails(ark.get())
    assert status == ARK_SUCCESS

    print("\nFinal Solver Statistics:")
    print(f"   Internal solver steps = {nst} (attempted = {nst_a})")
    print(f"   Total RHS evals:  Fe = {nfe},  Fi = {nfi}")
    print(f"   Total linear solver setups = {nsetups}")
    print(f"   Total linear iterations = {nli}")
    print(f"   Total number of Jacobian-vector products = {nJv}")
    print(f"   Total number of linear solver convergence failures = {nlcf}")
    print(f"   Total number of Newton iterations = {nni}")
    print(f"   Total number of nonlinear solver convergence failures = {ncfn}")
    print(f"   Total number of error test failures = {netf}")


# This function allows pytest to discover the example as a test
def test_ark_heat1D():
    main()


if __name__ == "__main__":
    main()
