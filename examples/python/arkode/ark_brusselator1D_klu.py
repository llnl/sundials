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
#   examples/arkode/C_klu/ark_brusselator1D_klu.c
# to use sundials4py.
# -----------------------------------------------------------------

import sys
import numpy as np
from sundials4py.core import *
from sundials4py.arkode import *


def IDX(x, v):
    return 3 * x + v


class Brusselator1DProblem:
    def __init__(self, N, a, b, du, dv, dw, ep, sunctx):
        self.N = N
        self.dx = 1.0 / (N - 1)
        self.a = a
        self.b = b
        self.du = du
        self.dv = dv
        self.dw = dw
        self.ep = ep
        self.sunctx = sunctx
        self.R = None

    def set_init_cond(self, yvec):
        y = N_VGetArrayPointer(yvec).reshape((self.N, 3))
        x = np.linspace(0.0, 1.0, self.N)
        s = 0.1 * np.sin(np.pi * x)
        y[:, 0] = self.a + s
        y[:, 1] = self.b / self.a + s
        y[:, 2] = self.b + s
        return 0

    def f(self, t, yvec, ydotvec, user_data):
        y = N_VGetArrayPointer(yvec).reshape((self.N, 3))
        ydot = N_VGetArrayPointer(ydotvec).reshape((self.N, 3))
        ydot[:, :] = 0.0

        u = y[:, 0]
        v = y[:, 1]
        w = y[:, 2]

        uconst = self.du / self.dx / self.dx
        vconst = self.dv / self.dx / self.dx
        wconst = self.dw / self.dx / self.dx

        ydot[1:-1, 0] = (
            (u[0:-2] - 2.0 * u[1:-1] + u[2:]) * uconst
            + self.a
            - (w[1:-1] + 1.0) * u[1:-1]
            + v[1:-1] * u[1:-1] * u[1:-1]
        )
        ydot[1:-1, 1] = (
            (v[0:-2] - 2.0 * v[1:-1] + v[2:]) * vconst
            + w[1:-1] * u[1:-1]
            - v[1:-1] * u[1:-1] * u[1:-1]
        )

        ydot[1:-1, 2] = (
            (w[0:-2] - 2.0 * w[1:-1] + w[2:]) * wconst
            + (self.b - w[1:-1]) / self.ep
            - w[1:-1] * u[1:-1]
        )

        # enforce stationary boundary conditions
        ydot[0, :] = np.zeros((1, 3))

        return 0

    def laplace_matrix(self, A):
        N = self.N
        colptrs = SUNSparseMatrix_IndexPointers(A)
        rowvals = SUNSparseMatrix_IndexValues(A)
        data = SUNSparseMatrix_Data(A)
        SUNMatZero(A)

        nz = 0
        colptrs[IDX(0, 0)] = nz
        colptrs[IDX(0, 1)] = nz
        colptrs[IDX(0, 2)] = nz

        uconst = self.du / self.dx / self.dx
        vconst = self.dv / self.dx / self.dx
        wconst = self.dw / self.dx / self.dx
        uconst2 = -2.0 * uconst
        vconst2 = -2.0 * vconst
        wconst2 = -2.0 * wconst

        for i in range(1, N - 1):
            colptrs[IDX(i, 0)] = nz
            if i > 1:
                data[nz] = uconst
                rowvals[nz] = IDX(i - 1, 0)
                nz += 1
            data[nz] = uconst2
            rowvals[nz] = IDX(i, 0)
            nz += 1
            if i < N - 2:
                data[nz] = uconst
                rowvals[nz] = IDX(i + 1, 0)
                nz += 1

            colptrs[IDX(i, 1)] = nz
            if i > 1:
                data[nz] = vconst
                rowvals[nz] = IDX(i - 1, 1)
                nz += 1
            data[nz] = vconst2
            rowvals[nz] = IDX(i, 1)
            nz += 1
            if i < N - 2:
                data[nz] = vconst
                rowvals[nz] = IDX(i + 1, 1)
                nz += 1

            colptrs[IDX(i, 2)] = nz
            if i > 1:
                data[nz] = wconst
                rowvals[nz] = IDX(i - 1, 2)
                nz += 1
            data[nz] = wconst2
            rowvals[nz] = IDX(i, 2)
            nz += 1
            if i < N - 2:
                data[nz] = wconst
                rowvals[nz] = IDX(i + 1, 2)
                nz += 1

        colptrs[IDX(N - 1, 0)] = nz
        colptrs[IDX(N - 1, 1)] = nz
        colptrs[IDX(N - 1, 2)] = nz
        colptrs[IDX(N - 1, 2) + 1] = nz
        return 0

    def reaction_jac(self, yvec, A):
        N = self.N
        y = N_VGetArrayPointer(yvec).reshape((self.N, 3))
        colptrs = SUNSparseMatrix_IndexPointers(A)
        rowvals = SUNSparseMatrix_IndexValues(A)
        data = SUNSparseMatrix_Data(A)
        SUNMatZero(A)

        nz = 0
        colptrs[IDX(0, 0)] = 0
        colptrs[IDX(0, 1)] = 0
        colptrs[IDX(0, 2)] = 0

        for i in range(1, N - 1):
            u = y[i, 0]
            v = y[i, 1]
            w = y[i, 2]

            colptrs[IDX(i, 0)] = nz
            rowvals[nz] = IDX(i, 0)
            data[nz] = 2.0 * u * v - w - 1.0
            nz += 1
            rowvals[nz] = IDX(i, 1)
            data[nz] = w - 2.0 * u * v
            nz += 1
            rowvals[nz] = IDX(i, 2)
            data[nz] = -w
            nz += 1

            colptrs[IDX(i, 1)] = nz
            rowvals[nz] = IDX(i, 0)
            data[nz] = u * u
            nz += 1
            rowvals[nz] = IDX(i, 1)
            data[nz] = -u * u
            nz += 1

            colptrs[IDX(i, 2)] = nz
            rowvals[nz] = IDX(i, 0)
            data[nz] = -u
            nz += 1
            rowvals[nz] = IDX(i, 1)
            data[nz] = u
            nz += 1
            rowvals[nz] = IDX(i, 2)
            data[nz] = -1.0 / self.ep - u
            nz += 1

        colptrs[IDX(N - 1, 0)] = nz
        colptrs[IDX(N - 1, 1)] = nz
        colptrs[IDX(N - 1, 2)] = nz
        colptrs[IDX(N - 1, 2) + 1] = nz
        return 0

    def jac(self, t, yvec, fyvec, J, tmp1, tmp2, tmp3, user_data):
        status = self.laplace_matrix(J)
        if status != 0:
            return status
        if self.R is None:
            self.R = SUNSparseMatrix(
                SUNSparseMatrix_Rows(J),
                SUNSparseMatrix_Columns(J),
                SUNSparseMatrix_NNZ(J),
                SUN_CSC_MAT,
                self.sunctx,
            )
        status = self.reaction_jac(yvec, self.R)
        if status != 0:
            return status
        return SUNMatScaleAdd(1.0, J, self.R)


def main():
    if not hasattr(sys.modules[__name__], "SUNLinSol_KLU"):
        raise RuntimeError("SUNLinSol_KLU is unavailable in this build")

    T0 = 0.0
    Tf = 10.0
    Nt = 10
    N = 201
    a = 0.6
    b = 2.0
    du = 0.025
    dv = 0.025
    dw = 0.025
    ep = 1.0e-5
    reltol = 1.0e-6
    abstol = 1.0e-10

    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert status == SUN_SUCCESS
    NEQ = 3 * N
    NNZ = 5 * NEQ

    y = N_VNew_Serial(NEQ, sunctx)
    assert y is not None
    umask = N_VClone(y)
    vmask = N_VClone(y)
    wmask = N_VClone(y)

    problem = Brusselator1DProblem(N, a, b, du, dv, dw, ep, sunctx)
    problem.set_init_cond(y)

    N_VConst(0.0, umask)
    N_VConst(0.0, vmask)
    N_VConst(0.0, wmask)
    umask_data = N_VGetArrayPointer(umask)
    vmask_data = N_VGetArrayPointer(vmask)
    wmask_data = N_VGetArrayPointer(wmask)
    umask_data.reshape((N, 3))[:, 0] = 1.0
    vmask_data.reshape((N, 3))[:, 1] = 1.0
    wmask_data.reshape((N, 3))[:, 2] = 1.0

    ark = ARKStepCreate(None, problem.f, T0, y, sunctx)
    assert ark is not None
    status = ARKodeSStolerances(ark.get(), reltol, abstol)
    assert status == ARK_SUCCESS

    A = SUNSparseMatrix(NEQ, NEQ, NNZ, SUN_CSC_MAT, sunctx)
    assert A is not None
    LS = SUNLinSol_KLU(y, A, sunctx)
    assert LS is not None

    status = ARKodeSetLinearSolver(ark.get(), LS, A)
    assert status == ARK_SUCCESS
    status = ARKodeSetJacFn(ark.get(), problem.jac)
    assert status == ARK_SUCCESS
    status = ARKodeSetAutonomous(ark.get(), 1)
    assert status == ARK_SUCCESS

    print("\n1D Brusselator PDE test problem (KLU solver):")
    print(f"    N = {N},  NEQ = {NEQ}")
    print(f"    problem parameters:  a = {a},  b = {b},  ep = {ep}")
    print(f"    diffusion coefficients:  du = {du},  dv = {dv},  dw = {dw}")
    print(f"    reltol = {reltol:.1e},  abstol = {abstol:.1e}\n")

    np.savetxt("bruss_mesh.txt", np.linspace(0.0, 1.0, N), fmt="%.16e")
    ydata = N_VGetArrayPointer(y).reshape((N, 3))
    with (
        open("bruss_u.txt", "w") as ufid,
        open("bruss_v.txt", "w") as vfid,
        open("bruss_w.txt", "w") as wfid,
    ):
        np.savetxt(ufid, ydata[:, 0][None, :], fmt="%.16e")
        np.savetxt(vfid, ydata[:, 1][None, :], fmt="%.16e")
        np.savetxt(wfid, ydata[:, 2][None, :], fmt="%.16e")

        t = T0
        dt_out = (Tf - T0) / Nt
        tout = T0 + dt_out
        print("        t      ||u||_rms   ||v||_rms   ||w||_rms")
        print("   ----------------------------------------------")
        for _ in range(Nt):
            status, t = ARKodeEvolve(ark.get(), tout, y, ARK_NORMAL)
            ydata = N_VGetArrayPointer(y).reshape((N, 3))
            u = N_VWL2Norm(y, umask)
            u = np.sqrt(u * u / N)
            v = N_VWL2Norm(y, vmask)
            v = np.sqrt(v * v / N)
            w = N_VWL2Norm(y, wmask)
            w = np.sqrt(w * w / N)
            print(f"  {t:10.6f}  {u:10.6f}  {v:10.6f}  {w:10.6f}")

            if status >= 0:
                tout = min(tout + dt_out, Tf)
            else:
                print("Solver failure, stopping integration")
                break

            np.savetxt(ufid, ydata[:, 0][None, :], fmt="%.16e")
            np.savetxt(vfid, ydata[:, 1][None, :], fmt="%.16e")
            np.savetxt(wfid, ydata[:, 2][None, :], fmt="%.16e")
        print("   ----------------------------------------------")

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
    status, nni = ARKodeGetNumNonlinSolvIters(ark.get())
    assert status == ARK_SUCCESS
    status, ncfn = ARKodeGetNumNonlinSolvConvFails(ark.get())
    assert status == ARK_SUCCESS
    status, nje = ARKodeGetNumJacEvals(ark.get())
    assert status == ARK_SUCCESS

    print("\nFinal Solver Statistics:")
    print(f"   Internal solver steps = {nst} (attempted = {nst_a})")
    print(f"   Total RHS evals:  Fe = {nfe},  Fi = {nfi}")
    print(f"   Total linear solver setups = {nsetups}")
    print(f"   Total number of Jacobian evaluations = {nje}")
    print(f"   Total number of nonlinear iterations = {nni}")
    print(f"   Total number of nonlinear solver convergence failures = {ncfn}")
    print(f"   Total number of error test failures = {netf}")


def test_ark_brusselator1d_klu():
    if not hasattr(sys.modules[__name__], "SUNLinSol_KLU"):
        return
    main()


if __name__ == "__main__":
    main()
