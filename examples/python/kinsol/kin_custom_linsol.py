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
# EXAMPLE: implementing a SUNMatrix *and* a SUNLinearSolver in Python.
#
# The test problem is the 1D nonlinear boundary value problem
#
#    -u'' + u^3 = f(x),    0 < x < 1,    u(0) = u(1) = 0,
#
# with f chosen (by substitution) so that the exact solution is
# u(x) = sin(pi*x). Discretizing on a uniform grid of N interior points with
# the standard second-order central difference for u'' turns this into a
# size-N nonlinear algebraic system in the nodal values. Its Jacobian is
# tridiagonal: the -u'' term contributes the constant off-diagonals, and the
# u^3 term only ever touches the diagonal. The tridiagonal system can be solved
# directly in O(N) work (Thomas algorithm) with no pivoting, since the diagonal
# here is always strictly larger in magnitude than the two off-diagonal entries
# in its row.
#
# WHY THESE TWO OBJECTS COME AS A PAIR
#
# A matrix and the linear solver that factors it are written together,
# because only they agree on how the entries are stored. That has a concrete
# consequence in Python: the matrix passed to the Jacobian function and to the
# solver's setup(A) and solve(A, x, b, tol) is the original CustomSUNMatrix
# subclass, so both can access its storage directly.
#
# WHICH OPERATIONS YOU MUST WRITE
#
# KINSOL only requires zero(); the other matrix operations below make this
# implementation reusable by packages with additional requirements.
#
# CustomSUNLinearSolver requires only solve(). Everything else -- initialize(),
# setup(), the statistics getters -- is optional: override a method and
# SUNDIALS will call it, leave it alone and SUNDIALS behaves as though the
# operation were absent, exactly as for a solver written in C.
# -----------------------------------------------------------------

import numpy as np
from sundials4py.core import *
from sundials4py.kinsol import *

# Problem constants
NEQ = 9
FTOL = 1.0e-10
STOL = 1.0e-10


class NonlinearBVP:
    """-u'' + u^3 = f(x) on N interior points, u(0) = u(1) = 0.

    f is chosen so that u(x) = sin(pi*x) solves the continuous problem; the
    discrete solution matches it to O(h^2), the discretization's own error.
    """

    def __init__(self, n=NEQ):
        self.n = n
        self.h = 1.0 / (n + 1)
        self.x = np.array([(i + 1) * self.h for i in range(n)])

        s = np.sin(np.pi * self.x)
        self.forcing = (np.pi**2) * s + s**3

    def residual(self, uvec, fvec, user_data):
        u = N_VGetNumpyArray(uvec)
        r = N_VGetNumpyArray(fvec)

        u_left = np.empty_like(u)
        u_left[0] = 0.0  # Dirichlet boundary value at x=0
        u_left[1:] = u[:-1]

        u_right = np.empty_like(u)
        u_right[:-1] = u[1:]
        u_right[-1] = 0.0  # Dirichlet boundary value at x=1

        h2 = self.h * self.h
        r[:] = (2.0 * u - u_left - u_right) / h2 + u**3 - self.forcing
        return 0

    def jacobian(self, uvec, fuvec, matrix, user_data, tmp1, tmp2):
        u = N_VGetNumpyArray(uvec)
        h2 = self.h * self.h

        diag = 2.0 / h2 + 3.0 * u**2
        sub = np.full(self.n, -1.0 / h2)
        sub[0] = 0.0  # row 0 has no u_{-1} neighbor
        sup = np.full(self.n, -1.0 / h2)
        sup[-1] = 0.0  # last row has no u_{n} neighbor

        matrix.set_entries(sub, diag, sup)
        return 0

    def solution(self):
        return np.sin(np.pi * self.x)


class MyMatrix(CustomSUNMatrix):
    """A tridiagonal SUNMatrix implemented in Python."""

    def __init__(self, n, sunctx):
        # Three arrays, all length n for uniform indexing: sub[i] is the
        # (i, i-1) entry (sub[0] unused), diag[i] is (i, i), and super[i] is
        # the (i, i+1) entry (super[n-1] unused).
        self.n = n
        self.sub = np.zeros(n, dtype=sunrealtype)
        self.diag = np.zeros(n, dtype=sunrealtype)
        self.super = np.zeros(n, dtype=sunrealtype)

        # The base constructor takes only the context. Note that it must be
        # called *after* your own state is in place: it makes the object
        # convertible to a native handle, and the operations below assume
        # that state exists.
        super().__init__(sunctx)

    def set_entries(self, sub, diag, super_):
        # Not a SUNMatrix operation -- just how this example's Jacobian
        # function writes into the matrix. Your own problem code can use
        # whatever accessor suits your storage.
        self.sub[:] = sub
        self.diag[:] = diag
        self.super[:] = super_

    # -- required operations -------------------------------------------------

    def clone(self):
        # Return a NEW, EMPTY matrix with the same shape and structure -- not
        # a copy of the entries. SUNDIALS may keep the result after your
        # reference to it is gone; that is expected and handled by the
        # binding, which holds a strong reference to the implementation for
        # as long as SUNDIALS owns the clone.
        return MyMatrix(self.n, self.sunctx)

    def zero(self):
        self.sub[:] = 0.0
        self.diag[:] = 0.0
        self.super[:] = 0.0
        return SUN_SUCCESS

    def copy(self, dst):
        # Copy self into dst.
        dst.sub[:] = self.sub
        dst.diag[:] = self.diag
        dst.super[:] = self.super
        return SUN_SUCCESS

    def scaleadd(self, c, other):
        # In place: self <- c*self + other.
        self.sub[:] = c * self.sub + other.sub
        self.diag[:] = c * self.diag + other.diag
        self.super[:] = c * self.super + other.super
        return SUN_SUCCESS

    def scaleaddi(self, c):
        # In place: self <- c*self + I. This is the operation an implicit ODE
        # integrator uses to form I - gamma*J, so it must add to the
        # diagonal rather than overwrite it.
        self.sub[:] *= c
        self.diag[:] = c * self.diag + 1.0
        self.super[:] *= c
        return SUN_SUCCESS

    def matvec(self, x, y):
        # y <- self*x. Unlike copy()/scaleadd(), x and y are N_Vectors, so
        # use the N_V* API (or N_VGetNumpyArray for a serial vector) rather
        # than assuming a representation.
        xv = N_VGetNumpyArray(x)
        yv = N_VGetNumpyArray(y)
        yv[:] = self.diag * xv
        yv[1:] += self.sub[1:] * xv[:-1]
        yv[:-1] += self.super[:-1] * xv[1:]
        return SUN_SUCCESS


class MyLinearSolver(CustomSUNLinearSolver):
    """A SUNLinearSolver implemented in Python, using the Thomas algorithm."""

    def __init__(self, matrix, sunctx):
        self.n = matrix.n

        # Populated by setup(): the modified diagonal and superdiagonal, and
        # the multipliers used to eliminate the subdiagonal. Precomputing
        # these here is what makes solve() itself cheap.
        self.diag_star = np.zeros(self.n)
        self.super_star = np.zeros(self.n)
        self.mult = np.zeros(self.n)

        # SUNLINEARSOLVER_DIRECT for solvers that perform an "exact" solve,
        # SUNLINEARSOLVER_ITERATIVE for a Krylov method,
        # SUNLINEARSOLVER_MATRIX_ITERATIVE for an iterative method that still
        # needs the matrix. The type tells the calling package how to use
        # you: a direct solver is asked for an exact solve and its `tol` is
        # ignored, while an iterative one is expected to honor `tol` and to
        # report num_iters() and res_norm().
        super().__init__(sunctx, SUNLINEARSOLVER_DIRECT)

    def setup(self, A):
        # Factor the matrix: the forward elimination pass of the Thomas
        # algorithm, which turns the tridiagonal system into a bidiagonal one
        # that solve() can finish cheaply. Called whenever the package
        # believes the matrix has changed.
        #
        n = self.n
        sub = A.sub
        diag = A.diag
        sup = A.super

        self.diag_star[0] = diag[0]
        self.super_star[0] = sup[0]
        for i in range(1, n):
            self.mult[i] = sub[i] / self.diag_star[i - 1]
            self.diag_star[i] = diag[i] - self.mult[i] * self.super_star[i - 1]
            self.super_star[i] = sup[i]
        return SUN_SUCCESS

    # -- the required operation ---------------------------------------------

    def solve(self, A, x, b, tol):
        """Solve A*x = b by the Thomas algorithm's back-substitution pass.

        A    the same matrix passed to setup()
        x    output N_Vector for the solution
        b    right-hand side N_Vector
        tol  requested residual tolerance; meaningful only for iterative
             types

        Return SUN_SUCCESS, or a POSITIVE code for a failure the caller can
        recover from (a singular matrix, or an iteration that did not
        converge -- KINSOL responds by shrinking its step), or a NEGATIVE
        code for a failure no retry will fix.
        """
        n = self.n
        d = N_VGetNumpyArray(b).copy()
        for i in range(1, n):
            d[i] -= self.mult[i] * d[i - 1]

        xv = N_VGetNumpyArray(x)
        xv[n - 1] = d[n - 1] / self.diag_star[n - 1]
        for i in range(n - 2, -1, -1):
            xv[i] = (d[i] - self.super_star[i] * xv[i + 1]) / self.diag_star[i]

        return SUN_SUCCESS


def main():
    print("\nNonlinear BVP test problem:")
    print(f"   -u'' + u^3 = f(x), 0 < x < 1, u(0) = u(1) = 0")
    print(f"   neq = {NEQ}")
    print(f"   ftol = {FTOL}, stol = {STOL}")
    print("Solution method: KINSOL Newton with a Python matrix and linear solver\n")

    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert status == SUN_SUCCESS

    J = MyMatrix(NEQ, sunctx)
    LS = MyLinearSolver(J, sunctx)
    problem = NonlinearBVP()

    u = N_VNew_Serial(NEQ, sunctx)
    scale = N_VNew_Serial(NEQ, sunctx)
    N_VConst(1.0, scale)

    kin = KINCreate(sunctx)
    assert KINInit(kin.get(), problem.residual, u) == KIN_SUCCESS
    assert KINSetFuncNormTol(kin.get(), FTOL) == KIN_SUCCESS
    assert KINSetScaledStepTol(kin.get(), STOL) == KIN_SUCCESS

    # Attach the custom objects. `J` and `LS` must stay referenced for as
    # long as KINSOL uses them: their native handles hold only a weak
    # reference back to the Python objects, so letting either go out of
    # scope here would leave KINSOL holding a dangling pointer.
    assert KINSetLinearSolver(kin.get(), LS, J) == KIN_SUCCESS
    assert KINSetJacFn(kin.get(), problem.jacobian) == KIN_SUCCESS

    # Initial guess: identically zero, away from the solution.
    N_VConst(0.0, u)

    status = KINSol(kin.get(), u, KIN_LINESEARCH, scale, scale)
    assert status == KIN_SUCCESS, f"KINSol returned {status}"

    computed = N_VGetNumpyArray(u)
    exact = problem.solution()
    print(f"{'x':>8}  {'computed':>14}  {'exact':>14}  {'error':>12}")
    print("-" * 54)
    for i in range(NEQ):
        err = abs(computed[i] - exact[i])
        print(f"{problem.x[i]:8.4f}  {computed[i]:14.8e}  {exact[i]:14.8e}  {err:12.4e}")

    status, nni = KINGetNumNonlinSolvIters(kin.get())
    assert status == KIN_SUCCESS
    status, nfe = KINGetNumFuncEvals(kin.get())
    assert status == KIN_SUCCESS
    status, nje = KINGetNumJacEvals(kin.get())
    assert status == KIN_SUCCESS

    print("\nFinal Statistics..\n")
    print(f"nni      = {nni:6d}    nfe     = {nfe:6d}")
    print(f"nje      = {nje:6d}")


def test_kin_custom_linsol():
    main()


if __name__ == "__main__":
    main()
