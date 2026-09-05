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
#   examples/arkode/C_superlu-mt/ark_brusselator1D_FEM_slu.c
# to use sundials4py.
# -----------------------------------------------------------------

import sys

import numpy as np
from sundials4py.arkode import *
from sundials4py.core import *


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
        self.x = np.linspace(0.0, 1.0, N)
        self.R = None

    def set_init_cond(self, yvec):
        y = N_VGetArrayPointer(yvec).reshape((self.N, 3))
        s = 0.1 * np.sin(np.pi * self.x)
        y[:, 0] = self.a + s
        y[:, 1] = self.b / self.a + s
        y[:, 2] = self.b + s
        return 0

    def f_diff(self, yvec, ydotvec):
        y = N_VGetArrayPointer(yvec).reshape((self.N, 3))
        ydot = N_VGetArrayPointer(ydotvec).reshape((self.N, 3))
        ydot[:, :] = 0.0

        u = y[:, 0]
        v = y[:, 1]
        w = y[:, 2]

        uconst = self.du / self.dx / self.dx
        vconst = self.dv / self.dx / self.dx
        wconst = self.dw / self.dx / self.dx

        ydot[1:-1, 0] = (u[:-2] - 2.0 * u[1:-1] + u[2:]) * uconst
        ydot[1:-1, 1] = (v[:-2] - 2.0 * v[1:-1] + v[2:]) * vconst
        ydot[1:-1, 2] = (w[:-2] - 2.0 * w[1:-1] + w[2:]) * wconst
        return 0

    def f_rx(self, yvec, ydotvec):
        y = N_VGetArrayPointer(yvec).reshape((self.N, 3))
        ydot = N_VGetArrayPointer(ydotvec).reshape((self.N, 3))

        u = y[:, 0]
        v = y[:, 1]
        w = y[:, 2]

        ydot[1:-1, 0] += self.a - (w[1:-1] + 1.0) * u[1:-1] + v[1:-1] * u[1:-1] * u[1:-1]
        ydot[1:-1, 1] += w[1:-1] * u[1:-1] - v[1:-1] * u[1:-1] * u[1:-1]
        ydot[1:-1, 2] += (self.b - w[1:-1]) / self.ep - w[1:-1] * u[1:-1]
        return 0

    def f(self, t, yvec, ydotvec, _):
        N_VConst(0.0, ydotvec)
        status = self.f_rx(yvec, ydotvec)
        if status != 0:
            return status
        return self.f_diff(yvec, ydotvec)

    def laplace_matrix(self, A):
        N = self.N
        rowptrs = SUNSparseMatrix_IndexPointers(A)
        colinds = SUNSparseMatrix_IndexValues(A)
        data = SUNSparseMatrix_Data(A)
        SUNMatZero(A)

        nz = 0
        rowptrs[IDX(0, 0)] = nz
        rowptrs[IDX(0, 1)] = nz
        rowptrs[IDX(0, 2)] = nz

        uconst = self.du / self.dx / self.dx
        vconst = self.dv / self.dx / self.dx
        wconst = self.dw / self.dx / self.dx
        uconst2 = -2.0 * uconst
        vconst2 = -2.0 * vconst
        wconst2 = -2.0 * wconst

        for i in range(1, N - 1):
            rowptrs[IDX(i, 0)] = nz
            if i > 1:
                data[nz] = uconst
                colinds[nz] = IDX(i - 1, 0)
                nz += 1
            data[nz] = uconst2
            colinds[nz] = IDX(i, 0)
            nz += 1
            if i < N - 2:
                data[nz] = uconst
                colinds[nz] = IDX(i + 1, 0)
                nz += 1

            rowptrs[IDX(i, 1)] = nz
            if i > 1:
                data[nz] = vconst
                colinds[nz] = IDX(i - 1, 1)
                nz += 1
            data[nz] = vconst2
            colinds[nz] = IDX(i, 1)
            nz += 1
            if i < N - 2:
                data[nz] = vconst
                colinds[nz] = IDX(i + 1, 1)
                nz += 1

            rowptrs[IDX(i, 2)] = nz
            if i > 1:
                data[nz] = wconst
                colinds[nz] = IDX(i - 1, 2)
                nz += 1
            data[nz] = wconst2
            colinds[nz] = IDX(i, 2)
            nz += 1
            if i < N - 2:
                data[nz] = wconst
                colinds[nz] = IDX(i + 1, 2)
                nz += 1

        rowptrs[IDX(N - 1, 0)] = nz
        rowptrs[IDX(N - 1, 1)] = nz
        rowptrs[IDX(N - 1, 2)] = nz
        rowptrs[IDX(N - 1, 2) + 1] = nz
        return 0

    def reaction_jac(self, yvec, A):
        y = N_VGetArrayPointer(yvec).reshape((self.N, 3))
        rowptrs = SUNSparseMatrix_IndexPointers(A)
        colinds = SUNSparseMatrix_IndexValues(A)
        data = SUNSparseMatrix_Data(A)
        SUNMatZero(A)

        nz = 0
        rowptrs[IDX(0, 0)] = 0
        rowptrs[IDX(0, 1)] = 0
        rowptrs[IDX(0, 2)] = 0

        for i in range(1, self.N - 1):
            u = y[i, 0]
            v = y[i, 1]
            w = y[i, 2]

            rowptrs[IDX(i, 0)] = nz
            colinds[nz] = IDX(i, 0)
            data[nz] = 2.0 * u * v - w - 1.0
            nz += 1
            colinds[nz] = IDX(i, 1)
            data[nz] = w - 2.0 * u * v
            nz += 1
            colinds[nz] = IDX(i, 2)
            data[nz] = -w
            nz += 1

            rowptrs[IDX(i, 1)] = nz
            colinds[nz] = IDX(i, 0)
            data[nz] = u * u
            nz += 1
            colinds[nz] = IDX(i, 1)
            data[nz] = -u * u
            nz += 1

            rowptrs[IDX(i, 2)] = nz
            colinds[nz] = IDX(i, 0)
            data[nz] = -u
            nz += 1
            colinds[nz] = IDX(i, 1)
            data[nz] = u
            nz += 1
            colinds[nz] = IDX(i, 2)
            data[nz] = -1.0 / self.ep - u
            nz += 1

        rowptrs[IDX(self.N - 1, 0)] = nz
        rowptrs[IDX(self.N - 1, 1)] = nz
        rowptrs[IDX(self.N - 1, 2)] = nz
        rowptrs[IDX(self.N - 1, 2) + 1] = nz
        return 0

    def jac(self, t, yvec, fyvec, J, tmp1, tmp2, tmp3, _):
        status = self.laplace_matrix(J)
        if status != 0:
            return status
        if self.R is None:
            self.R = SUNSparseMatrix(
                SUNSparseMatrix_Rows(J),
                SUNSparseMatrix_Columns(J),
                SUNSparseMatrix_NNZ(J),
                SUN_CSR_MAT,
                self.sunctx,
            )
        status = self.reaction_jac(yvec, self.R)
        if status != 0:
            return status
        return SUNMatScaleAdd(1.0, J, self.R)

    def mass_matrix(self, t, M, _user_data, tmp1, tmp2, tmp3):
        N = self.N
        rowptrs = SUNSparseMatrix_IndexPointers(M)
        colinds = SUNSparseMatrix_IndexValues(M)
        data = SUNSparseMatrix_Data(M)
        SUNMatZero(M)

        nz = 0
        x = self.x

        for i in range(N):
            left = i != 0
            right = i != N - 1

            xl = x[i - 1] if left else 0.0
            xc = x[i]
            xr = x[i + 1] if right else 0.0

            Ml = 0.0
            Mc = 0.0
            Mr = 0.0

            if left:
                chi_l1 = self._chiL(xl, xc, self._x1(xl, xc))
                chi_l2 = self._chiL(xl, xc, self._x2(xl, xc))
                chi_l3 = self._chiL(xl, xc, self._x3(xl, xc))
                chi_r1 = self._chiR(xl, xc, self._x1(xl, xc))
                chi_r2 = self._chiR(xl, xc, self._x2(xl, xc))
                chi_r3 = self._chiR(xl, xc, self._x3(xl, xc))
                Ml += self._quad(chi_l1 * chi_r1, chi_l2 * chi_r2, chi_l3 * chi_r3, xl, xc)
                Mc += self._quad(chi_r1 * chi_r1, chi_r2 * chi_r2, chi_r3 * chi_r3, xl, xc)

            if right:
                chi_l1 = self._chiL(xc, xr, self._x1(xc, xr))
                chi_l2 = self._chiL(xc, xr, self._x2(xc, xr))
                chi_l3 = self._chiL(xc, xr, self._x3(xc, xr))
                chi_r1 = self._chiR(xc, xr, self._x1(xc, xr))
                chi_r2 = self._chiR(xc, xr, self._x2(xc, xr))
                chi_r3 = self._chiR(xc, xr, self._x3(xc, xr))
                Mc += self._quad(chi_l1 * chi_l1, chi_l2 * chi_l2, chi_l3 * chi_l3, xc, xr)
                Mr += self._quad(chi_l1 * chi_r1, chi_l2 * chi_r2, chi_l3 * chi_r3, xc, xr)

            rowptrs[IDX(i, 0)] = nz
            if left:
                data[nz] = Ml
                colinds[nz] = IDX(i - 1, 0)
                nz += 1
            data[nz] = Mc
            colinds[nz] = IDX(i, 0)
            nz += 1
            if right:
                data[nz] = Mr
                colinds[nz] = IDX(i + 1, 0)
                nz += 1

            rowptrs[IDX(i, 1)] = nz
            if left:
                data[nz] = Ml
                colinds[nz] = IDX(i - 1, 1)
                nz += 1
            data[nz] = Mc
            colinds[nz] = IDX(i, 1)
            nz += 1
            if right:
                data[nz] = Mr
                colinds[nz] = IDX(i + 1, 1)
                nz += 1

            rowptrs[IDX(i, 2)] = nz
            if left:
                data[nz] = Ml
                colinds[nz] = IDX(i - 1, 2)
                nz += 1
            data[nz] = Mc
            colinds[nz] = IDX(i, 2)
            nz += 1
            if right:
                data[nz] = Mr
                colinds[nz] = IDX(i + 1, 2)
                nz += 1

        rowptrs[IDX(N - 1, 2) + 1] = nz
        return 0

    @staticmethod
    def _x1(xl, xr):
        return 0.5 * (xl + xr) - 0.5 * (xr - xl) * 0.7745966692414834

    @staticmethod
    def _x2(xl, xr):
        return 0.5 * (xl + xr)

    @staticmethod
    def _x3(xl, xr):
        return 0.5 * (xl + xr) + 0.5 * (xr - xl) * 0.7745966692414834

    @staticmethod
    def _chiL(xl, xr, x):
        return (xr - x) / (xr - xl)

    @staticmethod
    def _chiR(xl, xr, x):
        return (x - xl) / (xr - xl)

    @staticmethod
    def _quad(f1, f2, f3, xl, xr):
        return (
            0.5
            * (xr - xl)
            * (0.5555555555555556 * f1 + 0.8888888888888888 * f2 + 0.5555555555555556 * f3)
        )


def main():
    if not hasattr(sys.modules[__name__], "SUNLinSol_SuperLUMT"):
        raise RuntimeError("SUNLinSol_SuperLUMT is unavailable in this build")

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
    num_threads = 1

    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert status == SUN_SUCCESS

    NEQ = 3 * N
    NNZ = 15 * NEQ

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
    umask_data = N_VGetArrayPointer(umask).reshape((N, 3))
    vmask_data = N_VGetArrayPointer(vmask).reshape((N, 3))
    wmask_data = N_VGetArrayPointer(wmask).reshape((N, 3))
    umask_data[:, 0] = 1.0
    vmask_data[:, 1] = 1.0
    wmask_data[:, 2] = 1.0

    ark = ARKStepCreate(None, problem.f, T0, y, sunctx)
    assert ark is not None
    assert ARKodeSStolerances(ark.get(), reltol, abstol) == ARK_SUCCESS
    assert ARKodeResStolerance(ark.get(), abstol) == ARK_SUCCESS
    assert ARKodeSetAutonomous(ark.get(), 1) == ARK_SUCCESS

    A = SUNSparseMatrix(NEQ, NEQ, NNZ, SUN_CSR_MAT, sunctx)
    assert A is not None
    LS = SUNLinSol_SuperLUMT(y, A, num_threads, sunctx)
    assert LS is not None
    assert SUNLinSol_SuperLUMTSetOrdering(LS, 3) == SUN_SUCCESS

    M = SUNMatClone(A)
    assert M is not None
    MLS = SUNLinSol_SuperLUMT(y, M, num_threads, sunctx)
    assert MLS is not None

    assert ARKodeSetLinearSolver(ark.get(), LS, A) == ARK_SUCCESS
    assert ARKodeSetJacFn(ark.get(), problem.jac) == ARK_SUCCESS
    assert ARKodeSetMassLinearSolver(ark.get(), MLS, M, 0) == ARK_SUCCESS
    assert ARKodeSetMassFn(ark.get(), problem.mass_matrix) == ARK_SUCCESS

    print("\n1D FEM Brusselator PDE test problem (SuperLU_MT solver):")
    print(f"    N = {N},  NEQ = {NEQ}")
    print(f"    num_threads = {num_threads}")
    print(f"    problem parameters:  a = {a},  b = {b},  ep = {ep}")
    print(f"    diffusion coefficients:  du = {du},  dv = {dv},  dw = {dw}")
    print(f"    reltol = {reltol:.1e},  abstol = {abstol:.1e}\n")

    print("        t      ||u||_rms   ||v||_rms   ||w||_rms")
    print("   ----------------------------------------------")
    t = T0
    dt_out = (Tf - T0) / Nt
    tout = T0 + dt_out
    for _ in range(Nt):
        status, t = ARKodeEvolve(ark.get(), tout, y, ARK_NORMAL)
        u = np.sqrt((N_VWL2Norm(y, umask) ** 2) / N)
        v = np.sqrt((N_VWL2Norm(y, vmask) ** 2) / N)
        w = np.sqrt((N_VWL2Norm(y, wmask) ** 2) / N)
        print(f"  {t:10.6f}  {u:10.6f}  {v:10.6f}  {w:10.6f}")

        if status >= 0:
            tout = min(tout + dt_out, Tf)
        else:
            print("Solver failure, stopping integration")
            break

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
    status, nmset = ARKodeGetNumMassSetups(ark.get())
    assert status == ARK_SUCCESS
    status, nms = ARKodeGetNumMassSolves(ark.get())
    assert status == ARK_SUCCESS
    status, nMv = ARKodeGetNumMassMult(ark.get())
    assert status == ARK_SUCCESS

    print("\nFinal Solver Statistics:")
    print(f"   Internal solver steps = {nst} (attempted = {nst_a})")
    print(f"   Total RHS evals:  Fe = {nfe},  Fi = {nfi}")
    print(f"   Total mass matrix setups = {nmset}")
    print(f"   Total mass matrix solves = {nms}")
    print(f"   Total mass times evals = {nMv}")
    print(f"   Total linear solver setups = {nsetups}")
    print(f"   Total number of Jacobian evaluations = {nje}")
    print(f"   Total number of nonlinear iterations = {nni}")
    print(f"   Total number of nonlinear solver convergence failures = {ncfn}")
    print(f"   Total number of error test failures = {netf}")


def test_ark_brusselator1d_superlumt():
    if not hasattr(sys.modules[__name__], "SUNLinSol_SuperLUMT"):
        return
    main()


if __name__ == "__main__":
    main()
