# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2026, Lawrence Livermore National Security and the SUNDIALS contributors.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------
"""Regression tests for exception state around Python-facing bindings."""

import pytest
from fixtures import *
from sundials4py.core import *
from sundials4py.cvodes import *
from sundials4py.kinsol import *
from test_custom_nvector import ArrayVector, FailingLinearSumVector


class MissingCompare(ArrayVector):
    """A vector whose ``compare`` method is deliberately not overridden."""

    compare = CustomNVector.compare


class NestedFailingVector(ArrayVector):
    """Fail in a vector operation only after an inner solver call returns."""

    def __init__(self, values, sunctx):
        self.fail_linear_sum = False
        super().__init__(values, sunctx)

    def linear_sum(self, a, x, b, y):
        if self.fail_linear_sum:
            raise LookupError("nested linear sum boom")
        super().linear_sum(a, x, b, y)


def test_failed_entry_point_conversion_does_not_poison_thread(sunctx):
    """Argument conversion failures must not disable later custom operations."""
    healthy = ArrayVector([1.0], sunctx)

    def rhs(t, y, ydot, _):
        ydot.data[:] = -y.data
        return SUN_SUCCESS

    cvode = CVodeCreate(CV_BDF, sunctx)
    assert CVodeInit(cvode.get(), rhs, 0.0, healthy) == CV_SUCCESS
    with pytest.raises(TypeError):
        CVode(cvode.get(), 1.0, MissingCompare([1.0], sunctx), CV_NORMAL)

    x = ArrayVector([1.0, 2.0], sunctx)
    z = ArrayVector([0.0, 0.0], sunctx)
    N_VLinearSum(2.0, x, 1.0, x, z)
    assert list(z.data) == [3.0, 6.0]


def test_exception_out_of_solver_does_not_leave_capture_enabled(sunctx):
    """A callback exception must not poison a following direct binding call."""

    def rhs_raises(t, y, ydot, _):
        raise ValueError("rhs boom")

    y = ArrayVector([1.0], sunctx)
    cvode = CVodeCreate(CV_BDF, sunctx)
    assert CVodeInit(cvode.get(), rhs_raises, 0.0, y) == CV_SUCCESS
    assert CVodeSStolerances(cvode.get(), 1.0e-6, 1.0e-8) == CV_SUCCESS
    with pytest.raises(ValueError, match="rhs boom"):
        CVode(cvode.get(), 1.0, y, CV_NORMAL)

    bad = FailingLinearSumVector([1.0], sunctx)
    with pytest.raises(LookupError, match="linear sum boom"):
        N_VLinearSum(1.0, bad, 1.0, bad, bad)

    x = ArrayVector([1.0], sunctx)
    z = ArrayVector([0.0], sunctx)
    N_VLinearSum(1.0, x, 1.0, x, z)
    assert z.data[0] == 2.0


def test_native_abstol_with_custom_state_is_rejected_loudly(sunctx):
    """Auxiliary vectors of a custom-vector solve must not fail silently."""
    y = ArrayVector([1.0], sunctx)
    abstol = N_VNew_Serial(1, sunctx)
    N_VConst(1.0e-8, abstol)
    cvode = CVodeCreate(CV_BDF, sunctx)
    assert CVodeInit(cvode.get(), lambda t, y, yd, _: SUN_SUCCESS, 0.0, y) == CV_SUCCESS
    with pytest.raises(TypeError, match="different vector type"):
        CVodeSVtolerances(cvode.get(), 1.0e-6, abstol)


def test_nested_solver_does_not_disable_outer_capture(sunctx):
    """An inner KINSOL call must preserve the enclosing CVODE scope."""
    y = NestedFailingVector([1.0], sunctx)
    inner = N_VNew_Serial(1, sunctx)
    inner_scale = N_VNew_Serial(1, sunctx)
    N_VGetArrayPointer(inner)[0] = 1.0
    N_VConst(1.0, inner_scale)

    kin = KINCreate(sunctx)

    def residual(u, f, _):
        N_VLinearSum(1.0, u, -0.5, u, f)
        return SUN_SUCCESS

    assert KINInit(kin.get(), residual, inner) == KIN_SUCCESS

    def rhs(t, state, derivative, _):
        assert KINSol(kin.get(), inner, KIN_FP, inner_scale, inner_scale) == KIN_SUCCESS
        y.fail_linear_sum = True
        derivative.data[:] = -state.data
        return SUN_SUCCESS

    cvode = CVodeCreate(CV_BDF, sunctx)
    assert CVodeInit(cvode.get(), rhs, 0.0, y) == CV_SUCCESS
    assert CVodeSStolerances(cvode.get(), 1.0e-7, 1.0e-10) == CV_SUCCESS
    with pytest.raises(LookupError, match="nested linear sum boom"):
        CVode(cvode.get(), 0.1, y, CV_NORMAL)
