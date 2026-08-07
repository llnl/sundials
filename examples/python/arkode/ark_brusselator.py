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
# This is a direct port of the C example,
#   examples/arkode/C_serial/ark_brusselator.c,
# specifically with the parameters from Test 2.
#
# The following test simulates a brusselator problem from chemical
# kinetics.  This is an ODE system with 3 components, Y = [u,v,w],
# satisfying the equations,
#    du/dt = a - (w+1)*u + v*u^2
#    dv/dt = w*u - v*u^2
#    dw/dt = (b-w)/ep - w*u
# for t in the interval [0.0, 10.0], with initial conditions
# Y0 = [u0,v0,w0].
#
#    u0=1.2,  v0=3.1,  w0=3,  a=1,  b=3.5,  ep=5.0e-6
#    Here, w experiences a fast initial transient, jumping 0.5
#    within a few steps.  All values proceed smoothly until
#    around t=6.5, when both u and v undergo a sharp transition,
#    with u increasing from around 0.5 to 5 and v decreasing
#    from around 6 to 1 in less than 0.5 time units.  After this
#    transition, both u and v continue to evolve somewhat
#    rapidly for another 1.4 time units, and finish off smoothly.
#
# This program solves the problem with the DIRK method, using a
# Newton iteration with the SUNDIALS dense linear solver, and a
# user-supplied Jacobian routine.
#
# 100 outputs are printed at equal intervals, and run statistics
# are printed at the end.
# -----------------------------------------------------------------

import numpy as np
from sundials4py.core import *
from sundials4py.arkode import *


# Brusselator ODE problem class
class BrusselatorODE:
    def __init__(self, u0, v0, w0, a, b, ep):
        self.u0 = u0
        self.v0 = v0
        self.w0 = w0
        self.a = a
        self.b = b
        self.ep = ep

    def set_init_cond(self, yvec):
        y = N_VGetArrayPointer(yvec)
        y[0] = self.u0
        y[1] = self.v0
        y[2] = self.w0
        return 0

    def f(self, t, yvec, ydotvec, user_data):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)
        a, b, ep = self.a, self.b, self.ep
        u, v, w = y[0], y[1], y[2]
        ydot[0] = a - (w + 1.0) * u + v * u * u
        ydot[1] = w * u - v * u * u
        ydot[2] = (b - w) / ep - w * u
        return 0

    def jac(self, t, yvec, fyvec, J, tmp1, tmp2, tmp3, user_data):
        y = N_VGetArrayPointer(yvec)
        a, b, ep = self.a, self.b, self.ep
        u, v, w = y[0], y[1], y[2]
        Jdata = SUNDenseMatrix_Data(J)
        Jdata[0, 0] = -(w + 1.0) + 2.0 * u * v
        Jdata[0, 1] = u * u
        Jdata[0, 2] = -u
        Jdata[1, 0] = w - 2.0 * u * v
        Jdata[1, 1] = -u * u
        Jdata[1, 2] = u
        Jdata[2, 0] = -w
        Jdata[2, 1] = 0.0
        Jdata[2, 2] = -1.0 / ep - u
        return 0


def main():
    # Problem parameters for Test 2
    u0 = 1.2
    v0 = 3.1
    w0 = 3.0
    a = 1.0
    b = 3.5
    ep = 5.0e-6
    T0 = 0.0
    Tf = 10.0
    dTout = 1.0
    NEQ = 3
    Nt = int(np.ceil(Tf / dTout))
    reltol = 1.0e-6
    abstol = 1.0e-10

    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    y = N_VNew_Serial(NEQ, sunctx)

    # Create ODE problem instance and set initial conditions
    ode_problem = BrusselatorODE(u0, v0, w0, a, b, ep)
    ode_problem.set_init_cond(y)

    ark = ARKStepCreate(None, ode_problem.f, T0, y, sunctx)  # f_E (explicit)  # f_I (implicit)
    assert ark is not None

    status = ARKodeSStolerances(ark.get(), reltol, abstol)
    assert status == ARK_SUCCESS

    status = ARKodeSetInterpolantType(ark.get(), ARK_INTERP_LAGRANGE)
    assert status == ARK_SUCCESS

    status = ARKodeSetDeduceImplicitRhs(ark.get(), 1)
    assert status == ARK_SUCCESS

    A = SUNDenseMatrix(NEQ, NEQ, sunctx)
    assert A is not None

    LS = SUNLinSol_Dense(y, A, sunctx)
    assert LS is not None

    status = ARKodeSetLinearSolver(ark.get(), LS, A)
    assert status == ARK_SUCCESS

    status = ARKodeSetJacFn(ark.get(), ode_problem.jac)
    assert status == ARK_SUCCESS

    # Signal that this problem does not explicitly depend on time
    status = ARKodeSetAutonomous(ark.get(), 1)
    assert status == ARK_SUCCESS

    # Initial problem output
    yarr = N_VGetArrayPointer(y)
    print("\nBrusselator ODE test problem:")
    print(f"    initial conditions:  u0 = {u0},  v0 = {v0},  w0 = {w0}")
    print(f"    problem parameters:  a = {a},  b = {b},  ep = {ep}")
    print(f"    reltol = {reltol},  abstol = {abstol}\n")
    print("        t           u           v           w")
    print("   -------------------------------------------")
    print(f"  {T0:10.6f}  {yarr[0]:10.6f}  {yarr[1]:10.6f}  {yarr[2]:10.6f}")

    # Main time-stepping loop: calls ARKodeEvolve to perform the integration,
    # then prints results.  Stops when the final time has been reached. The
    # solution is written out to a file.
    with open("solution.txt", "w") as UFID:
        UFID.write("# t u v w\n")
        UFID.write(f" {T0:.16e} {yarr[0]:.16e} {yarr[1]:.16e} {yarr[2]:.16e}\n")
        tout = T0 + dTout
        for iout in range(Nt):
            status, tret = ARKodeEvolve(ark.get(), tout, y, ARK_NORMAL)
            yarr = N_VGetArrayPointer(y)
            print(f"  {tret:10.6f}  {yarr[0]:10.6f}  {yarr[1]:10.6f}  {yarr[2]:10.6f}")
            UFID.write(f" {tret:.16e} {yarr[0]:.16e} {yarr[1]:.16e} {yarr[2]:.16e}\n")
            if status == ARK_SUCCESS:
                tout += dTout
                tout = min(tout, Tf)
            else:
                print("Solver failure, stopping integration")
                break
        print("   -------------------------------------------")

    # Print statistics (check status code)
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
    status, netf = ARKodeGetNumErrTestFails(ark.get())
    assert status == ARK_SUCCESS
    status, ncfn = ARKodeGetNumStepSolveFails(ark.get())
    assert status == ARK_SUCCESS
    status, nni = ARKodeGetNumNonlinSolvIters(ark.get())
    assert status == ARK_SUCCESS
    status, nnf = ARKodeGetNumNonlinSolvConvFails(ark.get())
    assert status == ARK_SUCCESS
    status, nje = ARKodeGetNumJacEvals(ark.get())
    assert status == ARK_SUCCESS
    status, nfeLS = ARKodeGetNumLinRhsEvals(ark.get())
    assert status == ARK_SUCCESS

    print("\nFinal Solver Statistics:")
    print(f"   Internal solver steps = {nst} (attempted = {nst_a})")
    print(f"   Total RHS evals:  Fe = {nfe},  Fi = {nfi}")
    print(f"   Total linear solver setups = {nsetups}")
    print(f"   Total RHS evals for setting up the linear system = {nfeLS}")
    print(f"   Total number of Jacobian evaluations = {nje}")
    print(f"   Total number of Newton iterations = {nni}")
    print(f"   Total number of nonlinear solver convergence failures = {nnf}")
    print(f"   Total number of error test failures = {netf}")
    print(f"   Total number of failed steps from solver failure = {ncfn}")


# This function allows pytest to discover the example as a test
def test_ark_brusselator():
    main()


if __name__ == "__main__":
    main()
