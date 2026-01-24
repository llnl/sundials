#!/usr/bin/env python3
# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025, Lawrence Livermore National Security,
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
# Python port of the SUNDIALS slider-crank DAE example in IDAS
# (idasSlCrank_dns.c).
#
# Simulation of a slider-crank mechanism modeled with 3 generalized
# coordinates: crank angle, connecting bar angle, and slider location.
# The mechanism moves under the action of a constant horizontal force
# applied to the connecting rod and a spring-damper connecting the crank
# and connecting rod.
#
# The equations of motion are formulated as a system of stabilized
# index-2 DAEs (Gear-Gupta-Leimkuhler formulation).
# IDAS also computes the average kinetic energy as the quadrature:
#   G = int_t0^tend g(t,y,p) dt,
# where
#   g(t,y,p) = 0.5*J1*v1^2 + 0.5*J2*v3^2 + 0.5*m2*v2^2
# -----------------------------------------------------------------

import numpy as np
from sundials4py.core import *
from sundials4py.idas import *


class SliderCrankDAE:
    def __init__(self, a=0.5, J1=1.0, m2=1.0, m1=1.0, J2=2.0, l0=1.0, F=1.0, k=1.0, c=1.0):
        self.a = a
        self.J1 = J1
        self.m2 = m2
        self.m1 = m1
        self.J2 = J2
        self.l0 = l0
        self.F = F
        self.k = k
        self.c = c

    def set_initial_conditions(self, yyvec, ypvec):
        pi = 4.0 * np.arctan(1.0)
        a = self.a
        J1 = self.J1
        m2 = self.m2
        J2 = self.J2
        q = pi / 2.0
        p = np.arcsin(-a)
        x = np.cos(p)
        yy = N_VGetArrayPointer(yyvec)
        yp = N_VGetArrayPointer(ypvec)
        yy[:] = 0.0
        yp[:] = 0.0
        yy[0] = q
        yy[1] = x
        yy[2] = p
        Q = self.force(yy)
        yp[3] = Q[0] / J1
        yp[4] = Q[1] / m2
        yp[5] = Q[2] / J2
        return 0

    def force(self, yy):
        a, k, c, l0, F = self.a, self.k, self.c, self.l0, self.F
        q, x, p = yy[0], yy[1], yy[2]
        qd, xd, pd = yy[3], yy[4], yy[5]
        s1, c1 = np.sin(q), np.cos(q)
        s2, c2 = np.sin(p), np.cos(p)
        s21 = s2 * c1 - c2 * s1
        c21 = c2 * c1 + s2 * s1
        l2 = x * x - x * (c2 + a * c1) + (1.0 + a * a) / 4.0 + a * c21 / 2.0
        ell = np.sqrt(l2)
        ld = (
            2.0 * x * xd
            - xd * (c2 + a * c1)
            + x * (s2 * pd + a * s1 * qd)
            - a * s21 * (pd - qd) / 2.0
        ) / (2.0 * ell)
        f = k * (ell - l0) + c * ld
        fl = f / ell
        Q = np.zeros(3)
        Q[0] = -fl * a * (s21 / 2.0 + x * s1) / 2.0
        Q[1] = fl * (c2 / 2.0 - x + a * c1 / 2.0) + F
        Q[2] = -fl * (x * s2 - a * s21 / 2.0) / 2.0 - F * s2
        return Q

    def residual(self, t, yyvec, ypvec, rvec, user_data):
        a, J1, m2, J2 = self.a, self.J1, self.m2, self.J2
        yy = N_VGetArrayPointer(yyvec)
        yp = N_VGetArrayPointer(ypvec)
        rr = N_VGetArrayPointer(rvec)
        q, x, p = yy[0], yy[1], yy[2]
        qd, xd, pd = yy[3], yy[4], yy[5]
        lam1, lam2 = yy[6], yy[7]
        mu1, mu2 = yy[8], yy[9]
        s1, c1 = np.sin(q), np.cos(q)
        s2, c2 = np.sin(p), np.cos(p)
        Q = self.force(yy)
        rr[0] = yp[0] - qd + a * s1 * mu1 - a * c1 * mu2
        rr[1] = yp[1] - xd + mu1
        rr[2] = yp[2] - pd + s2 * mu1 - c2 * mu2
        rr[3] = J1 * yp[3] - Q[0] + a * s1 * lam1 - a * c1 * lam2
        rr[4] = m2 * yp[4] - Q[1] + lam1
        rr[5] = J2 * yp[5] - Q[2] + s2 * lam1 - c2 * lam2
        rr[6] = x - c2 - a * c1
        rr[7] = -s2 - a * s1
        rr[8] = a * s1 * qd + xd + s2 * pd
        rr[9] = -a * c1 * qd - c2 * pd
        return 0

    def rhsQ(self, t, yyvec, ypvec, qdotvec, user_data):
        J1, m2, J2 = self.J1, self.m2, self.J2
        yy = N_VGetArrayPointer(yyvec)
        qdot = N_VGetArrayPointer(qdotvec)
        v1, v2, v3 = yy[3], yy[4], yy[5]
        qdot[0] = 0.5 * (J1 * v1 * v1 + m2 * v2 * v2 + J2 * v3 * v3)
        return 0


def main():
    # Problem parameters
    RTOLF = 1e-6
    ATOLF = 1e-7
    RTOLQ = 1e-6
    ATOLQ = 1e-8
    TBEGIN = 0.0
    TEND = 10.0
    NEQ = 10
    NOUT = 25

    # Create model
    dae = SliderCrankDAE()

    # SUNDIALS context
    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert status == SUN_SUCCESS

    id = N_VNew_Serial(NEQ, sunctx)
    assert id is not None
    yy = N_VClone(id)
    assert yy is not None
    yp = N_VClone(id)
    assert yp is not None
    q = N_VNew_Serial(1, sunctx)
    assert q is not None

    # Consistent IC
    id_arr = N_VGetArrayPointer(id)
    id_arr[:6] = 1.0
    id_arr[6:] = 0.0
    dae.set_initial_conditions(yy, yp)

    # IDAS initialization
    ida = IDACreate(sunctx)
    assert ida is not None
    status = IDAInit(ida.get(), dae.residual, TBEGIN, yy, yp)
    assert status == IDA_SUCCESS
    status = IDASStolerances(ida.get(), RTOLF, ATOLF)
    assert status == IDA_SUCCESS
    status = IDASetId(ida.get(), id)
    assert status == IDA_SUCCESS
    status = IDASetSuppressAlg(ida.get(), True)
    assert status == IDA_SUCCESS
    status = IDASetMaxNumSteps(ida.get(), 20000)
    assert status == IDA_SUCCESS

    # Create dense SUNMatrix to use with dense linear solver
    A = SUNDenseMatrix(NEQ, NEQ, sunctx)
    assert A is not None

    # Create dense linear solver
    LS = SUNLinSol_Dense(yy, A, sunctx)
    assert LS is not None

    # Attach the matrix and linear solver
    status = IDASetLinearSolver(ida.get(), LS, A)
    assert status == IDA_SUCCESS

    #  Setup quadrature
    N_VConst(0.0, q)
    status = IDAQuadInit(ida.get(), dae.rhsQ, q)
    assert status == IDA_SUCCESS
    status = IDASetQuadErrCon(ida.get(), True)
    assert status == IDA_SUCCESS
    status = IDAQuadSStolerances(ida.get(), RTOLQ, ATOLQ)
    assert status == IDA_SUCCESS

    # Output header
    print("\nidasSlCrank_dns: Slider-Crank DAE serial example problem for IDAS")
    print("Linear solver: DENSE, Jacobian is computed by IDAS.")
    print(f"Tolerance parameters:  rtol = {RTOLF}   atol = {ATOLF}")
    print("---------------------------------------------------------------------")
    print("  t         y1          y2           y3      | nst  k      h")
    print("---------------------------------------------------------------------")

    # Time stepping loop (C example style)
    yarr = N_VGetArrayPointer(yy)
    tout = TEND / NOUT
    tret = 0.0
    while True:
        status, nst = IDAGetNumSteps(ida.get())
        status, kused = IDAGetLastOrder(ida.get())
        status, hused = IDAGetLastStep(ida.get())
        print(
            f"{tret:5.2f} {yarr[0]:12.4e} {yarr[1]:12.4e} {yarr[2]:12.4e} | {nst:3d}  {kused:1d} {hused:12.4e}"
        )

        status, tret = IDASolve(ida.get(), tout, yy, yp, IDA_NORMAL)
        assert status >= 0

        tout += TEND / NOUT
        if tret > TEND:
            status, nst = IDAGetNumSteps(ida.get())
            status, kused = IDAGetLastOrder(ida.get())
            status, hused = IDAGetLastStep(ida.get())
            print(
                f"{tret:5.2f} {yarr[0]:12.4e} {yarr[1]:12.4e} {yarr[2]:12.4e} | {nst:3d}  {kused:1d} {hused:12.4e}"
            )
            break

    # Final statistics
    status, nst = IDAGetNumSteps(ida.get())
    status, nre = IDAGetNumResEvals(ida.get())
    status, nje = IDAGetNumJacEvals(ida.get())
    status, nni = IDAGetNumNonlinSolvIters(ida.get())
    status, netf = IDAGetNumErrTestFails(ida.get())
    status, nnf = IDAGetNumNonlinSolvConvFails(ida.get())
    status, ncfn = IDAGetNumStepSolveFails(ida.get())
    status, nreLS = IDAGetNumLinResEvals(ida.get())

    print("\nFinal Run Statistics: \n")
    print(f"Number of steps                    = {nst}")
    print(f"Number of residual evaluations     = {nre + nreLS}")
    print(f"Number of Jacobian evaluations     = {nje}")
    print(f"Number of nonlinear iterations     = {nni}")
    print(f"Number of error test failures      = {netf}")
    print(f"Number of nonlinear conv. failures = {nnf}")
    print(f"Number of step solver failures     = {ncfn}")

    status, tret = IDAGetQuad(ida.get(), q)
    print("--------------------------------------------")
    print(f"  G = {N_VGetArrayPointer(q)[0]}")
    print("--------------------------------------------\n")


# This function allows pytest to discover the example as a test
def test_idaSlCrank_dns():
    main()


if __name__ == "__main__":
    main()
