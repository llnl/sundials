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
# TEMPLATE: implementing a SUNMatrix *and* a SUNLinearSolver in Python.
#
# This file is a starting point, not a finished example. The scaffolding -- the
# problem, its Jacobian, the KINSOL setup, the output -- is complete; the parts
# marked `TODO(you):` are where your matrix and linear solver go. Fill them in,
# delete the pytest.skip() at the bottom, and run the file.
#
# The test problem is the nonlinear system
#
#    x - (1/3)*cos((y-1)*z) - 1/6                = 0
#    y - (1/9)*sqrt(x^2 + sin(z) + 1.06) - 0.9   = 0
#    z + (1/20)*exp(-x*(y-1)) + (10*pi - 3)/60   = 0
#
# whose solution is (x, y, z) = (1/2, 1, -pi/6). It is small enough to check by
# hand and nonlinear enough that KINSOL's Newton iteration needs several
# Jacobian factorizations to converge, so both objects below get exercised.
#
# WHY THESE TWO OBJECTS COME AS A PAIR
#
# A matrix and the linear solver that factors it are written together, because
# only they agree on how the entries are stored. That has a concrete consequence
# in Python: when KINSOL calls your solver's setup(A) and solve(A, x, b, tol), the
# `A` it passes is an *opaque native SUNMatrix handle*, not the Python object you
# created. So the solver reads the matrix through a reference it was given at
# construction -- see MyLinearSolver.__init__ below -- and ignores the handle. The
# same holds for the Jacobian function: it writes into the Python matrix object it
# closes over rather than into the handle KINSOL passes it.
#
# WHICH OPERATIONS YOU MUST WRITE
#
# CustomSUNMatrix requires all six of clone(), zero(), copy(), scaleadd(),
# scaleaddi(), and matvec(); the binding refuses to build a native handle without
# them, since a package has no way to ask which are missing. KINSOL itself only
# calls zero() and clone(), but an ODE package forming I - gamma*J needs
# scaleaddi(), and an iterative solver needs matvec(), so implementing all six
# makes the matrix reusable.
#
# CustomSUNLinearSolver requires only solve(). Everything else -- initialize(),
# setup(), the statistics getters -- is optional: override a method and SUNDIALS
# will call it, leave it alone and SUNDIALS behaves as though the operation were
# absent, exactly as for a solver written in C.
# -----------------------------------------------------------------

import numpy as np
import pytest
from sundials4py.core import *
from sundials4py.kinsol import *

# Problem constants
NEQ = 3
FTOL = 1.0e-10
STOL = 1.0e-10


class NonlinearSystem:
    """The 3-equation system described above, with its analytic Jacobian."""

    def __init__(self, matrix):
        # The Jacobian is written into this Python matrix object; see the note on
        # opaque handles in the file header.
        self.matrix = matrix

    def residual(self, uvec, fvec, user_data):
        x, y, z = N_VGetArrayPointer(uvec)
        f = N_VGetArrayPointer(fvec)
        f[0] = x - (1.0 / 3.0) * np.cos((y - 1.0) * z) - (1.0 / 6.0)
        f[1] = y - (1.0 / 9.0) * np.sqrt(x * x + np.sin(z) + 1.06) - 0.9
        f[2] = z + (1.0 / 20.0) * np.exp(-x * (y - 1.0)) + (10.0 * np.pi - 3.0) / 60.0
        return 0

    def jacobian(self, uvec, fuvec, Jhandle, user_data, tmp1, tmp2):
        # Jhandle is the same matrix as self.matrix, but arrives as an opaque
        # native handle, so the entries are written through the Python object.
        x, y, z = N_VGetArrayPointer(uvec)
        r = np.sqrt(x * x + np.sin(z) + 1.06)
        e = np.exp(-x * (y - 1.0))
        s = np.sin((y - 1.0) * z)

        self.matrix.set_entries(
            [
                [1.0, (1.0 / 3.0) * z * s, (1.0 / 3.0) * (y - 1.0) * s],
                [-x / (9.0 * r), 1.0, -np.cos(z) / (18.0 * r)],
                [-(y - 1.0) * e / 20.0, -x * e / 20.0, 1.0],
            ]
        )
        return 0

    def solution(self):
        return np.array([0.5, 1.0, -np.pi / 6.0])


class MyMatrix(CustomSUNMatrix):
    """A SUNMatrix implemented in Python.

    Reference implementations worth reading alongside this template:
      src/sunmatrix/dense/sunmatrix_dense.c
      src/sunmatrix/sparse/sunmatrix_sparse.c
    """

    def __init__(self, rows, cols, sunctx):
        # TODO(you): choose your storage. A dense NumPy array is used below
        # because it keeps the operations to one line each; a banded, sparse, or
        # GPU-resident layout works just as well, since nothing outside this
        # class ever looks at the representation.
        self.shape = (rows, cols)
        self.data = np.zeros((rows, cols), dtype=sunrealtype)

        # The base constructor takes only the context. Note that it must be
        # called *after* your own state is in place: it makes the object
        # convertible to a native handle, and the operations below assume that
        # state exists.
        super().__init__(sunctx)

    def set_entries(self, entries):
        # Not a SUNMatrix operation -- just how this example's Jacobian function
        # writes into the matrix. Your own problem code can use whatever
        # accessor suits your storage.
        self.data[:] = entries

    # -- required operations -------------------------------------------------

    def clone(self):
        # Return a NEW, EMPTY matrix with the same shape and structure -- not a
        # copy of the entries. SUNDIALS may keep the result after your reference
        # to it is gone; that is expected and handled by the binding, which holds
        # a strong reference to the implementation for as long as SUNDIALS owns
        # the clone.
        #
        # TODO(you): return an empty matrix laid out like self.
        return MyMatrix(*self.shape, self.sunctx)

    def zero(self):
        # Set every entry to zero. Return SUN_SUCCESS, or a nonzero error code.
        # TODO(you):
        #   self.data[:] = 0.0
        raise NotImplementedError("TODO(you): implement MyMatrix.zero")

    def copy(self, dst):
        # Copy self into dst, which is the PYTHON matrix object -- clone() built
        # it, so it is your own type and you may use your own attributes.
        # TODO(you):
        #   dst.data[:] = self.data
        raise NotImplementedError("TODO(you): implement MyMatrix.copy")

    def scaleadd(self, c, other):
        # In place: self <- c*self + other. `other` is a Python matrix object.
        # TODO(you):
        #   self.data[:] = c * self.data + other.data
        raise NotImplementedError("TODO(you): implement MyMatrix.scaleadd")

    def scaleaddi(self, c):
        # In place: self <- c*self + I. This is the operation an implicit ODE
        # integrator uses to form I - gamma*J, so it must add to the diagonal
        # rather than overwrite it.
        # TODO(you):
        #   self.data[:] = c * self.data
        #   self.data[np.diag_indices_from(self.data)] += 1.0
        raise NotImplementedError("TODO(you): implement MyMatrix.scaleaddi")

    def matvec(self, x, y):
        # y <- self*x. Unlike copy()/scaleadd(), x and y are N_Vectors, so use
        # the N_V* API (or N_VGetArrayPointer for a serial vector) rather than
        # assuming a representation.
        # TODO(you):
        #   N_VGetArrayPointer(y)[:] = self.data @ N_VGetArrayPointer(x)
        raise NotImplementedError("TODO(you): implement MyMatrix.matvec")


class MyLinearSolver(CustomSUNLinearSolver):
    """A SUNLinearSolver implemented in Python.

    Reference implementations worth reading alongside this template:
      src/sunlinsol/dense/sunlinsol_dense.c        (direct)
      src/sunlinsol/spgmr/sunlinsol_spgmr.c        (iterative, matrix-free)
    """

    def __init__(self, matrix, sunctx):
        # The matrix this solver factors. It is supplied here rather than read
        # from the handle setup()/solve() receive, because that handle is opaque.
        self.matrix = matrix

        # TODO(you): allocate whatever the factorization needs -- pivots, a
        # workspace copy, a preconditioner.
        self.factor = None

        # SUNLINEARSOLVER_DIRECT for a factorization,
        # SUNLINEARSOLVER_ITERATIVE for a Krylov method,
        # SUNLINEARSOLVER_MATRIX_ITERATIVE for an iterative method that still
        # needs the matrix. The type tells the calling package how to use you:
        # a direct solver is asked for an exact solve and its `tol` is ignored,
        # while an iterative one is expected to honor `tol` and to report
        # num_iters() and res_norm().
        super().__init__(sunctx, SUNLINEARSOLVER_DIRECT)

    # -- optional: one-time setup -------------------------------------------

    def initialize(self):
        # Called once before the first setup()/solve(). Override this if the
        # solver has state that must be established after construction but
        # before use -- for example, checking that the matrix shape is one you
        # can handle.
        rows, cols = self.matrix.shape
        if rows != cols:
            return SUN_ERR_ARG_DIMSMISMATCH
        return SUN_SUCCESS

    def setup(self, A):
        # Factor the matrix. Called whenever the package believes the matrix has
        # changed; a well-written solver does the expensive work here and leaves
        # solve() cheap.
        #
        # `A` is an opaque native handle for the same matrix as self.matrix --
        # read the entries through self.matrix instead.
        #
        # TODO(you): factor self.matrix into self.factor. Take a copy if your
        # factorization is destructive: the package is free to modify the matrix
        # after this returns, and is free to call solve() many times before the
        # next setup().
        #
        #   self.factor = self.matrix.data.copy()
        raise NotImplementedError("TODO(you): implement MyLinearSolver.setup")

    # -- the required operation ---------------------------------------------

    def solve(self, A, x, b, tol):
        """Solve A*x = b.

        A    opaque native handle for the matrix (read self.matrix instead)
        x    output N_Vector for the solution
        b    right-hand side N_Vector
        tol  requested residual tolerance; meaningful only for iterative types

        Return SUN_SUCCESS, or a POSITIVE code for a failure the caller can
        recover from (a singular matrix, or an iteration that did not converge --
        KINSOL responds by shrinking its step), or a NEGATIVE code for a failure
        no retry will fix.
        """
        # TODO(you): compute x from b using the factorization built in setup().
        #
        #   try:
        #       N_VGetArrayPointer(x)[:] = np.linalg.solve(
        #           self.factor, N_VGetArrayPointer(b)
        #       )
        #   except np.linalg.LinAlgError:
        #       # Recoverable: report it, do not raise. An exception here
        #       # propagates out through the C layer as an error, which denies
        #       # KINSOL the chance to retry with a fresh Jacobian.
        #       return 1
        #
        #   return SUN_SUCCESS
        raise NotImplementedError("TODO(you): implement MyLinearSolver.solve")

    # -- optional: the rest of the interface ---------------------------------
    #
    # An iterative solver would also override:
    #
    #   set_atimes(atimes)                  matrix-vector product callback
    #   set_preconditioner(psetup, psolve)  preconditioner callbacks
    #   set_scaling_vectors(s1, s2)         left/right scaling
    #   set_zero_guess(onoff)               x arrives as zero, skip the product
    #   num_iters()  -> int                 iterations in the last solve
    #   res_norm()   -> float               final residual norm
    #   resid()      -> N_Vector            the residual vector itself
    #
    # These are omitted here because a direct solver has nothing useful to
    # report for any of them, and an unimplemented optional method is simply
    # absent from the native operation table.


def main():
    print("\nNonlinear system test problem:")
    print(f"   neq = {NEQ}")
    print(f"   ftol = {FTOL}, stol = {STOL}")
    print("Solution method: KINSOL Newton with a Python matrix and linear solver\n")

    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert status == SUN_SUCCESS

    # The matrix, the solver that factors it, and the problem that fills it in
    # all refer to the same Python matrix object.
    J = MyMatrix(NEQ, NEQ, sunctx)
    LS = MyLinearSolver(J, sunctx)
    problem = NonlinearSystem(J)

    u = N_VNew_Serial(NEQ, sunctx)
    scale = N_VNew_Serial(NEQ, sunctx)
    N_VConst(1.0, scale)

    kin = KINCreate(sunctx)
    assert KINInit(kin.get(), problem.residual, u) == KIN_SUCCESS
    assert KINSetFuncNormTol(kin.get(), FTOL) == KIN_SUCCESS
    assert KINSetScaledStepTol(kin.get(), STOL) == KIN_SUCCESS

    # Attach the custom objects. `J` and `LS` must stay referenced for as long as
    # KINSOL uses them: their native handles hold only a weak reference back to
    # the Python objects, so letting either go out of scope here would leave
    # KINSOL holding a dangling pointer.
    assert KINSetLinearSolver(kin.get(), LS, J) == KIN_SUCCESS
    assert KINSetJacFn(kin.get(), problem.jacobian) == KIN_SUCCESS

    # Initial guess, deliberately away from the solution.
    N_VGetArrayPointer(u)[:] = [0.1, 0.1, -0.1]

    status = KINSol(kin.get(), u, KIN_LINESEARCH, scale, scale)
    assert status == KIN_SUCCESS, f"KINSol returned {status}"

    computed = N_VGetArrayPointer(u)
    exact = problem.solution()
    print(f"{'':>3}  {'computed':>16}  {'exact':>16}  {'error':>12}")
    print("-" * 54)
    for i, name in enumerate("xyz"):
        err = abs(computed[i] - exact[i])
        print(f"{name:>3}  {computed[i]:16.8e}  {exact[i]:16.8e}  {err:12.4e}")

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
    # This is a template, so there is nothing complete for CI to verify yet.
    # Delete this skip once you have filled in the TODO(you) sections.
    pytest.skip("template example: fill in the TODO(you) sections first")
    main()


if __name__ == "__main__":
    main()
