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

from fixtures import *
from sundials4py.core import *


class FailingMatrix(CustomSUNMatrix):
    def __init__(self, sunctx):
        self.fail_zero = True
        self.fail_clone = True
        super().__init__(sunctx)

    def clone(self):
        if self.fail_clone:
            self.fail_clone = False
            raise LookupError("matrix clone boom")
        return type(self)(self.sunctx)

    def zero(self):
        if self.fail_zero:
            self.fail_zero = False
            raise ValueError("matrix zero boom")
        return SUN_SUCCESS

    def copy(self, dst):
        return SUN_SUCCESS

    def scaleadd(self, c, other):
        return SUN_SUCCESS

    def scaleaddi(self, c):
        return SUN_SUCCESS

    def matvec(self, x, y):
        return SUN_SUCCESS


class FailingLinearSolver(CustomSUNLinearSolver):
    def __init__(self, sunctx):
        self.fail_res_norm = True
        super().__init__(sunctx, SUNLINEARSOLVER_ITERATIVE)

    def solve(self, A, x, b, tol):
        return SUN_SUCCESS

    def res_norm(self):
        if self.fail_res_norm:
            self.fail_res_norm = False
            raise ArithmeticError("linear residual norm boom")
        return 0.5


class FailingNonlinearSolver(CustomSUNNonlinearSolver):
    def __init__(self, sunctx):
        self.fail_get_num_iters = True
        super().__init__(sunctx, SUNNONLINEARSOLVER_ROOTFIND)

    def solve(self, y0, y, w, tol, call_lsetup):
        return SUN_SUCCESS

    def get_num_iters(self):
        if self.fail_get_num_iters:
            self.fail_get_num_iters = False
            raise RuntimeError("nonlinear getter boom")
        return SUN_SUCCESS, 7


class FailingController(CustomSUNHController):
    def __init__(self, sunctx):
        self.fail_estimate = True
        super().__init__(sunctx)

    def estimate_step(self, h, p, dsm):
        if self.fail_estimate:
            self.fail_estimate = False
            raise ZeroDivisionError("controller estimate boom")
        return SUN_SUCCESS, h / 2


def record_errors(sunctx):
    errors = []

    def handler(line, operation, file, message, code, data, context):
        errors.append((operation, message, code))

    assert SUNContext_PushErrHandler(sunctx, handler) == SUN_SUCCESS
    return errors


def assert_reported(sunctx, errors, operation, exception):
    assert SUNContext_PeekLastError(sunctx) == SUN_ERR_EXT_FAIL
    assert errors[-1][0] == operation
    assert errors[-1][2] == SUN_ERR_EXT_FAIL
    assert operation in errors[-1][1]
    assert exception in errors[-1][1]


def test_status_exception_is_reported_and_python_state_is_cleared(sunctx):
    # Purpose:
    # Status exception is reported and python state is cleared.
    errors = record_errors(sunctx)
    A = FailingMatrix(sunctx)

    assert SUNMatZero(A) == SUN_ERR_EXT_FAIL
    assert_reported(sunctx, errors, "CustomSUNMatrix.zero", "ValueError: matrix zero boom")
    assert SUNMatZero(A) == SUN_SUCCESS


def test_pointer_exception_is_reported_and_python_state_is_cleared(sunctx):
    # Purpose:
    # Pointer exception is reported and python state is cleared.
    errors = record_errors(sunctx)
    A = FailingMatrix(sunctx)

    assert SUNMatClone(A) is None
    assert_reported(sunctx, errors, "CustomSUNMatrix.clone", "LookupError: matrix clone boom")
    assert SUNMatClone(A) is not None


def test_real_scalar_exception_is_reported_and_python_state_is_cleared(sunctx):
    # Purpose:
    # Real scalar exception is reported and python state is cleared.
    errors = record_errors(sunctx)
    LS = FailingLinearSolver(sunctx)

    assert SUNLinSolResNorm(LS) == 0.0
    assert_reported(
        sunctx, errors, "custom_linsol_resnorm", "ArithmeticError: linear residual norm boom"
    )
    assert SUNLinSolResNorm(LS) == 0.5


def test_nonlinear_output_exception_is_reported_and_python_state_is_cleared(sunctx):
    # Purpose:
    # Nonlinear output exception is reported and python state is cleared.
    errors = record_errors(sunctx)
    NLS = FailingNonlinearSolver(sunctx)

    status, value = SUNNonlinSolGetNumIters(NLS)
    assert status == SUN_ERR_EXT_FAIL
    assert value == 0
    assert_reported(sunctx, errors, "call_tuple_getter", "RuntimeError: nonlinear getter boom")
    assert SUNNonlinSolGetNumIters(NLS) == (SUN_SUCCESS, 7)


def test_controller_output_exception_is_reported_and_python_state_is_cleared(sunctx):
    # Purpose:
    # Controller output exception is reported and python state is cleared.
    errors = record_errors(sunctx)
    C = FailingController(sunctx)

    status, hnew = SUNAdaptController_EstimateStep(C, 2.0, 1, 1.0)
    assert status == SUN_ERR_EXT_FAIL
    assert hnew == 0.0
    assert_reported(
        sunctx,
        errors,
        "custom_controller_estimatestep",
        "ZeroDivisionError: controller estimate boom",
    )
    assert SUNAdaptController_EstimateStep(C, 2.0, 1, 1.0) == (SUN_SUCCESS, 1.0)


def test_error_handler_exception_does_not_cross_the_c_boundary(sunctx):
    # Purpose:
    # Error handler exception does not cross the c boundary.
    def failing_handler(line, operation, file, message, code, data, context):
        raise RuntimeError("error handler also failed")

    assert SUNContext_PushErrHandler(sunctx, failing_handler) == SUN_SUCCESS
    A = FailingMatrix(sunctx)

    assert SUNMatZero(A) == SUN_ERR_EXT_FAIL
    assert SUNContext_PeekLastError(sunctx) == SUN_ERR_EXT_FAIL
    assert SUNMatZero(A) == SUN_SUCCESS
