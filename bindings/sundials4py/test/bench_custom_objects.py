# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------
"""Opt-in performance measurements for Python custom SUNDIALS objects.

These tests intentionally report measurements instead of asserting a timing
threshold.  Run them with ``pytest -m bench -s``.  The normal test command
excludes this module through the ``bench`` marker in ``pyproject.toml``.
"""

import platform
import sys
import time

import numpy as np
import pytest
from fixtures import *
from sundials4py.core import *
from sundials4py.cvodes import *
from sundials4py.kinsol import *
from test_custom_nvector import ArrayVector
from test_custom_object_integration import DenseArrayMatrix, NumpyLinearSolver


pytestmark = pytest.mark.bench


class NoArgumentSolver(CustomSUNLinearSolver):
    """Minimal solver used to isolate a Python method with no arguments."""

    def __init__(self, sunctx):
        self.calls = 0
        super().__init__(sunctx, SUNLINEARSOLVER_DIRECT)

    def initialize(self):
        self.calls += 1
        return SUN_SUCCESS

    def solve(self, A, x, b, tol):
        return SUN_SUCCESS


def _timed(function, loops, samples=7):
    for _ in range(2):
        function()
    elapsed = []
    for _ in range(samples):
        start = time.perf_counter_ns()
        for _ in range(loops):
            function()
        elapsed.append((time.perf_counter_ns() - start) / loops)
    return float(np.median(elapsed))


def test_custom_dispatch_benchmark(sunctx):
    print(f"Python {platform.python_version()} ({sys.implementation.name}); {platform.platform()}")
    solver = NoArgumentSolver(sunctx)
    SUNLinSolInitialize(solver)
    per_call = _timed(lambda: SUNLinSolInitialize(solver), loops=5000)
    print(f"custom no-argument operation: {per_call:.1f} ns/call")
    assert solver.calls > 2


@pytest.mark.parametrize("length, loops", [(10, 1000), (1000, 100), (1_000_000, 3)])
def test_vector_operation_benchmark(sunctx, length, loops):
    native_x = N_VNew_Serial(length, sunctx)
    native_y = N_VNew_Serial(length, sunctx)
    native_z = N_VNew_Serial(length, sunctx)
    N_VConst(1.0, native_x)
    N_VConst(2.0, native_y)
    custom_x = ArrayVector(np.ones(length), sunctx)
    custom_y = ArrayVector(np.full(length, 2.0), sunctx)
    custom_z = ArrayVector(np.zeros(length), sunctx)

    native_sum = _timed(lambda: N_VLinearSum(1.0, native_x, 1.0, native_y, native_z), loops)
    custom_sum = _timed(lambda: N_VLinearSum(1.0, custom_x, 1.0, custom_y, custom_z), loops)
    native_norm = _timed(lambda: N_VWrmsNorm(native_x, native_y), loops)
    custom_norm = _timed(lambda: N_VWrmsNorm(custom_x, custom_y), loops)
    print(
        f"length={length}: N_VLinearSum native={native_sum:.1f} ns "
        f"custom={custom_sum:.1f} ns overhead={custom_sum - native_sum:.1f} ns; "
        f"N_VWrmsNorm native={native_norm:.1f} ns custom={custom_norm:.1f} ns "
        f"overhead={custom_norm - native_norm:.1f} ns"
    )


def _cvode_solve(sunctx, custom):
    y = ArrayVector([1.0], sunctx) if custom else N_VNew_Serial(1, sunctx)
    if not custom:
        N_VConst(1.0, y)

    def rhs(t, state, derivative, _):
        if custom:
            derivative.data[:] = -state.data
        else:
            N_VGetArrayPointer(derivative)[:] = -N_VGetArrayPointer(state)
        return SUN_SUCCESS

    cvode = CVodeCreate(CV_BDF, sunctx)
    assert CVodeInit(cvode.get(), rhs, 0.0, y) == CV_SUCCESS
    assert CVodeSStolerances(cvode.get(), 1.0e-6, 1.0e-10) == CV_SUCCESS
    ls = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, sunctx)
    assert CVodeSetLinearSolver(cvode.get(), ls, None) == CV_SUCCESS
    status, tret = CVode(cvode.get(), 1.0, y, CV_NORMAL)
    assert status == CV_SUCCESS
    assert tret == pytest.approx(1.0)


def test_cvode_native_and_custom_vector_benchmark(sunctx):
    native = _timed(lambda: _cvode_solve(sunctx, False), loops=3, samples=3)
    custom = _timed(lambda: _cvode_solve(sunctx, True), loops=3, samples=3)
    print(f"CVODE BDF native={native / 1e3:.1f} us custom={custom / 1e3:.1f} us")


def _kinsol_solve(sunctx, custom):
    u = N_VNew_Serial(1, sunctx)
    scale = N_VNew_Serial(1, sunctx)
    N_VConst(1.0, scale)

    def residual(state, result, _):
        N_VGetArrayPointer(result)[0] = N_VGetArrayPointer(state)[0] - 1.0
        return SUN_SUCCESS

    kin = KINCreate(sunctx)
    assert KINInit(kin.get(), residual, u) == KIN_SUCCESS
    if custom:
        matrix = DenseArrayMatrix(1, 1, sunctx)
        solver = NumpyLinearSolver(matrix, sunctx)

        def jacobian(state, result, matrix, _, tmp1, tmp2):
            matrix.data[0, 0] = 1.0
            return SUN_SUCCESS

        assert KINSetLinearSolver(kin.get(), solver, matrix) == KIN_SUCCESS
        assert KINSetJacFn(kin.get(), jacobian) == KIN_SUCCESS
    else:
        matrix = SUNDenseMatrix(1, 1, sunctx)
        solver = SUNLinSol_Dense(u, matrix, sunctx)
        assert KINSetLinearSolver(kin.get(), solver, matrix) == KIN_SUCCESS
    assert KINSol(kin.get(), u, KIN_NONE, scale, scale) == KIN_SUCCESS


def test_kinsol_custom_and_native_solver_benchmark(sunctx):
    native = _timed(lambda: _kinsol_solve(sunctx, False), loops=3, samples=3)
    custom = _timed(lambda: _kinsol_solve(sunctx, True), loops=3, samples=3)
    print(f"KINSOL dense={native / 1e3:.1f} us custom matrix/solver={custom / 1e3:.1f} us")
