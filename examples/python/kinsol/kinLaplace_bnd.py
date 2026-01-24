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
# This is a direct port of the C example,
#   examples/kinsol/serial/kinLaplace_bnd.c
#
# This example solves a 2D elliptic PDE
#
#    d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u - 2.0
#
# subject to homogeneous Dirichlet boundary conditions.
# The PDE is discretized on a uniform NX+2 by NY+2 grid with
# central differencing, and with boundary values eliminated,
# leaving a system of size NEQ = NX*NY.
# The nonlinear system is solved by KINSOL using the SUNBAND linear
# solver.
# -----------------------------------------------------------------

import numpy as np
from sundials4py.core import *
from sundials4py.kinsol import *

# Problem Constants
NX = 31
NY = 31
NEQ = NX * NY
SKIP = 3
FTOL = 1e-12
ZERO = 0.0
ONE = 1.0
TWO = 2.0


class Laplace2D:
    def __init__(self, NX, NY):
        self.NX = NX
        self.NY = NY
        self.NEQ = NX * NY
        self.dx = ONE / (NX + 1)
        self.dy = ONE / (NY + 1)
        self.hdc = ONE / (self.dx * self.dx)
        self.vdc = ONE / (self.dy * self.dy)

    def set_init_cond(self, yvec):
        y = N_VGetArrayPointer(yvec)
        # Initial guess: zero everywhere
        y[:] = ZERO
        return 0

    def func(self, uvec, fvec, user_data):
        u = N_VGetArrayPointer(uvec)
        f = N_VGetArrayPointer(fvec)

        NX, NY = self.NX, self.NY
        hdc, vdc = self.hdc, self.vdc

        # 2D views of the interior data (zero-copy)
        u2d = u.reshape((NX, NY))
        f2d = f.reshape((NX, NY))

        for i in range(NX):
            for j in range(NY):
                uij = u2d[i, j]

                # homogeneous Dirichlet boundaries (0 outside domain)
                udn = u2d[i, j - 1] if j > 0 else 0.0
                uup = u2d[i, j + 1] if j < NY - 1 else 0.0
                ult = u2d[i - 1, j] if i > 0 else 0.0
                urt = u2d[i + 1, j] if i < NX - 1 else 0.0

                hdiff = hdc * (ult - TWO * uij + urt)
                vdiff = vdc * (uup - TWO * uij + udn)

                f2d[i, j] = hdiff + vdiff + uij - uij * uij * uij + 2.0

        return 0

    def print_output(self, yvec):
        y = N_VGetArrayPointer(yvec).reshape((self.NX, self.NY))
        print("            ", end="")
        for i in range(1, self.NX + 1, SKIP):
            x = i * self.dx
            print(f"{x:<8.5f} ", end="")
        print("\n")
        for j in range(1, self.NY + 1, SKIP):
            yval = j * self.dy
            print(f"{yval:<8.5f}    ", end="")
            for i in range(1, self.NX + 1, SKIP):
                print(f"{y[i-1,j-1]:<8.5f} ", end="")
            print()


def main():
    print("\n2D elliptic PDE on unit square")
    print("   d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u + 2.0")
    print(" + homogeneous Dirichlet boundary conditions\n")
    print("Solution method: Modified Newton with band linear solver")
    print(f"Problem size: {NX} x {NY} = {NEQ}\n")

    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert status == SUN_SUCCESS

    y = N_VNew_Serial(NEQ, sunctx)
    assert y is not None

    scale = N_VNew_Serial(NEQ, sunctx)
    assert scale is not None

    # Create problem instance and set initial guess
    problem = Laplace2D(NX, NY)
    problem.set_init_cond(y)

    # Initialize and allocate memory for KINSOL
    kin = KINCreate(sunctx)
    status = KINInit(kin.get(), problem.func, y)
    assert status == KIN_SUCCESS

    # Set function norm tolerance
    status = KINSetFuncNormTol(kin.get(), FTOL)
    assert status == KIN_SUCCESS

    # Create band matrix and linear solver
    J = SUNBandMatrix(NEQ, NX, NX, sunctx)
    assert J is not None

    LS = SUNLinSol_Band(y, J, sunctx)
    assert LS is not None

    status = KINSetLinearSolver(kin.get(), LS, J)
    assert status == KIN_SUCCESS

    # Set Modified Newton parameters
    status = KINSetMaxSetupCalls(kin.get(), 100)
    assert status == KIN_SUCCESS
    status = KINSetMaxSubSetupCalls(kin.get(), 1)
    assert status == KIN_SUCCESS

    # No scaling used
    N_VConst(ONE, scale)

    # Call KINSol to solve problem
    status = KINSol(kin.get(), y, KIN_LINESEARCH, scale, scale)
    assert status == KIN_SUCCESS

    # Get scaled norm of the system function
    status, fnorm = KINGetFuncNorm(kin.get())
    assert status == KIN_SUCCESS
    print(f"\nComputed solution (||F|| = {fnorm}):\n")
    problem.print_output(y)

    # Print final statistics (faithful to C PrintFinalStats)
    status, nni = KINGetNumNonlinSolvIters(kin.get())
    assert status == KIN_SUCCESS

    status, nfe = KINGetNumFuncEvals(kin.get())
    assert status == KIN_SUCCESS

    status, nbcfails = KINGetNumBetaCondFails(kin.get())
    assert status == KIN_SUCCESS

    status, nbacktr = KINGetNumBacktrackOps(kin.get())
    assert status == KIN_SUCCESS

    status, nje = KINGetNumJacEvals(kin.get())
    assert status == KIN_SUCCESS

    status, nfeD = KINGetNumLinFuncEvals(kin.get())
    assert status == KIN_SUCCESS

    print("\nFinal Statistics.. \n")
    print(f"nni      = {nni:6d}    nfe     = {nfe:6d}")
    print(f"nbcfails = {nbcfails:6d}    nbacktr = {nbacktr:6d}")
    print(f"nje      = {nje:6d}    nfeB    = {nfeD:6d}")


def test_kinLaplace_bnd():
    main()


if __name__ == "__main__":
    main()
