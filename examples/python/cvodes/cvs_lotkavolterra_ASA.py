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

import numpy as np
from sundials4py.core import *
from sundials4py.cvodes import *


class LotkaVolterraODE:
    def __init__(self, p):
        self.p = np.array(p, dtype=sunrealtype)
        self.NEQ = 2
        self.NP = 4

    def set_init_cond(self, yvec):
        # Set initial condition u0 = [1.0, 1.0]
        y = N_VGetArrayPointer(yvec)
        y[0] = 1.0
        y[1] = 1.0
        return 0

    def f(self, t, yvec, ydotvec, user_data):
        # Lotka-Volterra ODE right-hand side
        p = self.p
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)
        ydot[0] = p[0] * y[0] - p[1] * y[0] * y[1]
        ydot[1] = -p[2] * y[1] + p[3] * y[0] * y[1]
        return 0

    def vjp(self, vvec, Jvvec, t, yvec):
        # Jacobian-vector product v^T (df/du)
        p = self.p
        v = N_VGetArrayPointer(vvec)
        Jv = N_VGetArrayPointer(Jvvec)
        y = N_VGetArrayPointer(yvec)
        Jv[0] = (p[0] - p[1] * y[1]) * v[0] + p[3] * y[1] * v[1]
        Jv[1] = -p[1] * y[0] * v[0] + (-p[2] + p[3] * y[0]) * v[1]
        return 0

    def parameter_vjp(self, vvec, Jvvec, t, yvec):
        # Parameter Jacobian-vector product v^T (df/dp)
        v = N_VGetArrayPointer(vvec)
        Jv = N_VGetArrayPointer(Jvvec)
        y = N_VGetArrayPointer(yvec)
        # Derivatives w.r.t. each parameter
        Jv[0] = y[0] * v[0]
        Jv[1] = -y[0] * y[1] * v[0]
        Jv[2] = -y[1] * v[1]
        Jv[3] = y[0] * y[1] * v[1]
        return 0

    def dgdu(self, yvec):
        # Gradient of the cost function w.r.t. u
        y = N_VGetArrayPointer(yvec)
        # g(u) = 0.5 * ||1 - u||^2, so grad = u - 1
        return np.array([-1.0 + y[0], -1.0 + y[1]], dtype=sunrealtype)

    def adjoint_rhs(self, t, yvec, lvec, ldotvec, user_datas):
        # Adjoint ODE right-hand side: -mu^T (df/du)
        self.vjp(lvec, ldotvec, t, yvec)
        ldot = N_VGetArrayPointer(ldotvec)
        ldot *= -1.0
        return 0

    def quad_rhs(self, t, yvec, muvec, qBdotvec, user_data):
        # Quadrature right-hand side: mu^T (df/dp)
        self.parameter_vjp(muvec, qBdotvec, t, yvec)
        return 0


def main():
    # Problem parameters
    p = [1.5, 1.0, 3.0, 1.0]
    T0 = 0.0
    Tf = 10.0
    reltol = 1e-10
    abstol = 1e-14
    steps = 5
    ode = LotkaVolterraODE(p)
    NEQ = ode.NEQ
    NP = ode.NP

    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert status == SUN_SUCCESS

    y = N_VNew_Serial(NEQ, sunctx)
    assert y is not None

    # Set initial condition
    ode.set_init_cond(y)

    # Create CVODE solver and set up problem
    cvode = CVodeCreate(CV_BDF, sunctx)
    assert cvode is not None

    # Initialize CVODE with ODE RHS
    status = CVodeInit(cvode.get(), ode.f, T0, y)
    assert status == CV_SUCCESS

    # Set tolerances
    status = CVodeSStolerances(cvode.get(), reltol, abstol)
    assert status == CV_SUCCESS

    # Set max steps
    status = CVodeSetMaxNumSteps(cvode.get(), 100000)
    assert status == CV_SUCCESS

    # Set linear solver
    ls = SUNLinSol_SPGMR(y, 0, 3, sunctx)
    assert ls is not None
    status = CVodeSetLinearSolver(cvode.get(), ls, None)
    assert status == CV_SUCCESS

    # Enable adjoint sensitivity analysis
    status = CVodeAdjInit(cvode.get(), steps, CV_HERMITE)
    assert status == CV_SUCCESS

    # Output problem setup
    print("\nLotka-Volterra ODE test problem (CVODE, ASA):")
    print(f"    initial conditions:  y0 = [1.0, 1.0]")
    print(f"    parameters:  p = {p}")
    print(f"    reltol = {reltol}, abstol = {abstol}\n")
    print("        t           x           y")
    print("   ---------------------------------")
    yarr = N_VGetArrayPointer(y)
    print(f"  {T0:10.6f}  {yarr[0]:10.6f}  {yarr[1]:10.6f}")

    # Forward integration
    tout = Tf
    t = 0.0
    status, tret, ncheck = CVodeF(cvode.get(), tout, y, CV_NORMAL)
    assert status == CV_SUCCESS

    yarr = N_VGetArrayPointer(y)
    print(f"  {tout:10.6f}  {yarr[0]:10.6f}  {yarr[1]:10.6f}")
    print("   ---------------------------------")

    # Adjoint terminal condition
    uB = N_VNew_Serial(NEQ, sunctx)
    assert uB is not None
    arr_uB = N_VGetArrayPointer(uB)
    arr_uB[:] = ode.dgdu(y)
    qB = N_VNew_Serial(NP, sunctx)
    assert qB is not None
    N_VConst(0.0, qB)
    print("Adjoint terminal condition:")
    print(arr_uB)
    print(N_VGetArrayPointer(qB))

    # Create the CVODES object for the backward problem
    status, which = CVodeCreateB(cvode.get(), CV_BDF)
    assert status == CV_SUCCESS

    # Initialize the CVODES solver for the backward problem
    status = CVodeInitB(cvode.get(), which, ode.adjoint_rhs, Tf, uB)
    assert status == CV_SUCCESS

    # Set the tolerances for the backward problem
    status = CVodeSStolerancesB(cvode.get(), which, reltol, abstol)
    assert status == CV_SUCCESS

    # Create the linear solver for the backward problem
    lsb = SUNLinSol_SPGMR(uB, 0, 3, sunctx)
    assert lsb is not None
    status = CVodeSetLinearSolverB(cvode.get(), which, lsb, None)
    assert status == CV_SUCCESS

    # Call CVodeQuadInitB to allocate internal memory and initialize backward
    # quadrature integration. This gives the sensitivities w.r.t. the parameters.
    status = CVodeQuadInitB(cvode.get(), which, ode.quad_rhs, qB)
    assert status == CV_SUCCESS

    # Call CVodeSetQuadErrCon to specify whether or not the quadrature variables
    # are to be used in the step size control mechanism within CVODES. Call
    # CVodeQuadSStolerances or CVodeQuadSVtolerances to specify the integration
    # tolerances for the quadrature variables.
    status = CVodeSetQuadErrConB(cvode.get(), which, True)
    assert status == CV_SUCCESS

    # Call CVodeQuadSStolerancesB to specify the scalar relative and absolute tolerances
    # for the backward problem.
    status = CVodeQuadSStolerancesB(cvode.get(), which, reltol, abstol)
    assert status == CV_SUCCESS

    # Integrate the adjoint ODE
    status = CVodeB(cvode.get(), T0, CV_NORMAL)
    assert status >= 0

    # Get the final adjoint solution
    status, t = CVodeGetB(cvode.get(), which, uB)
    assert status == CV_SUCCESS

    # Call CVodeGetQuadB to get the quadrature solution vector after a
    # successful return from CVodeB.
    status, t = CVodeGetQuadB(cvode.get(), which, qB)
    assert status == CV_SUCCESS

    # dg/dp = -qB
    arr_qB = N_VGetArrayPointer(qB)
    arr_qB *= -1.0
    print(f"Adjoint Solution at t = {t}:")
    print(N_VGetArrayPointer(uB))
    print(arr_qB)


# This function allows pytest to discover the example as a test
def test_cvs_lotkavolterra_ASA():
    main()


if __name__ == "__main__":
    main()
