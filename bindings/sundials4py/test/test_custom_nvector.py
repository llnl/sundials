# -----------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ UMBC
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

import numpy as np
import pytest
from fixtures import *
from numpy.testing import assert_allclose
from sundials4py.arkode import *
from sundials4py.core import *
from sundials4py.cvodes import *
from sundials4py.idas import *
from sundials4py.kinsol import *
from sundials4py.test import N_VTestHandleAddress, N_VTestStrongClone


class ArrayVector(CustomNVector):
    def __init__(self, values, sunctx):
        self.data = np.asarray(values, dtype=sunrealtype).copy()
        super().__init__(sunctx)

    def clone(self):
        return type(self)(np.zeros_like(self.data), self.sunctx)

    def get_array_pointer(self):
        return self.data.ctypes.data

    def get_length(self):
        return self.data.size

    def linear_sum(self, a, x, b, y):
        self.data[:] = a * x.data + b * y.data

    def const(self, c):
        self.data[:] = c

    def prod(self, x, y):
        self.data[:] = x.data * y.data

    def div(self, x, y):
        self.data[:] = x.data / y.data

    def scale(self, c, x):
        self.data[:] = c * x.data

    def abs(self, x):
        self.data[:] = np.abs(x.data)

    def inv(self, x):
        self.data[:] = 1.0 / x.data

    def add_const(self, x, b):
        self.data[:] = x.data + b

    def dot_prod(self, y):
        return np.dot(self.data, y.data)

    def max_norm(self):
        return np.max(np.abs(self.data))

    def wrms_norm(self, w):
        return np.sqrt(np.mean((self.data * w.data) ** 2))

    def wrms_norm_mask(self, w, mask):
        return np.sqrt(np.mean(((self.data * w.data) * (mask.data != 0)) ** 2))

    def min(self):
        return np.min(self.data)

    def wl2_norm(self, w):
        return np.sqrt(np.sum((self.data * w.data) ** 2))

    def l1_norm(self):
        return np.sum(np.abs(self.data))

    def compare(self, c, x):
        self.data[:] = np.abs(x.data) >= c

    def inv_test(self, x):
        if np.any(x.data == 0):
            return SUNFALSE
        self.data[:] = 1.0 / x.data
        return SUNTRUE

    def constr_mask(self, constraints, mask):
        mask.data[:] = 0
        return SUNTRUE

    def min_quotient(self, denom):
        nz = denom.data != 0
        return np.min(self.data[nz] / denom.data[nz])

    def linear_combination(self, coefficients, vectors):
        self.data[:] = sum(c * v.data for c, v in zip(coefficients, vectors))
        return SUN_SUCCESS

    def dot_prod_multi(self, vectors):
        return SUN_SUCCESS, [np.dot(self.data, v.data) for v in vectors]

    def dot_prod_local(self, other):
        return np.dot(self.data, other.data)


class SlotsVector(CustomNVector):
    """Complete custom vector whose payload is stored in ``__slots__``."""

    __slots__ = ("data",)

    def __init__(self, values, sunctx):
        self.data = np.asarray(values, dtype=sunrealtype).copy()
        super().__init__(sunctx)

    clone = ArrayVector.clone
    get_array_pointer = ArrayVector.get_array_pointer
    get_length = ArrayVector.get_length
    linear_sum = ArrayVector.linear_sum
    const = ArrayVector.const
    prod = ArrayVector.prod
    div = ArrayVector.div
    scale = ArrayVector.scale
    abs = ArrayVector.abs
    inv = ArrayVector.inv
    add_const = ArrayVector.add_const
    dot_prod = ArrayVector.dot_prod
    max_norm = ArrayVector.max_norm
    wrms_norm = ArrayVector.wrms_norm
    wrms_norm_mask = ArrayVector.wrms_norm_mask
    min = ArrayVector.min
    wl2_norm = ArrayVector.wl2_norm
    l1_norm = ArrayVector.l1_norm
    compare = ArrayVector.compare
    inv_test = ArrayVector.inv_test
    constr_mask = ArrayVector.constr_mask
    min_quotient = ArrayVector.min_quotient


class IncompleteVector(CustomNVector):
    def get_length(self):
        return 1


class OtherArrayVector(ArrayVector):
    pass


class CompatibleArrayVector(ArrayVector):
    def is_compatible(self, other):
        return isinstance(other, OtherArrayVector)


class WrongContextVector(ArrayVector):
    def __init__(self, values, sunctx, clone_context):
        self.clone_context = clone_context
        super().__init__(values, sunctx)

    def clone(self):
        return ArrayVector(np.zeros_like(self.data), self.clone_context)


class FailingLinearSumVector(ArrayVector):
    def linear_sum(self, a, x, b, y):
        raise LookupError("linear sum boom")


class FailingMaxNormVector(ArrayVector):
    def max_norm(self):
        raise ArithmeticError("max norm boom")


class FailingArrayReductionVector(ArrayVector):
    def wrms_norm_vector_array(self, vectors, weights):
        raise ArithmeticError("array norm boom")

    def dot_prod_multi(self, vectors):
        raise ArithmeticError("array dot boom")


def test_custom_nvector_dispatch_and_clone(sunctx):
    x = ArrayVector([1.0, -2.0, 3.0], sunctx)
    y = ArrayVector([4.0, 5.0, 6.0], sunctx)
    z = ArrayVector([0.0, 0.0, 0.0], sunctx)

    assert N_VGetVectorID(x) == SUNDIALS_NVEC_CUSTOM
    assert N_VGetLength(x) == 3
    N_VLinearSum(2.0, x, -1.0, y, z)
    assert_allclose(z.data, [-2.0, -9.0, 0.0])
    assert N_VDotProd(x, y) == pytest.approx(12.0)
    assert N_VMaxNorm(x) == pytest.approx(3.0)

    cloned = N_VClone(x)
    assert isinstance(cloned, ArrayVector)
    assert cloned is not x
    assert_allclose(cloned.data, 0.0)

    cloned_empty = N_VCloneEmpty(x)
    assert isinstance(cloned_empty, ArrayVector)
    assert cloned_empty is not x


def test_custom_nvector_clone_reuses_native_shell(sunctx):
    x = ArrayVector([1.0, 2.0], sunctx)
    cloned, keepalive = N_VTestStrongClone(x)

    first = N_VTestHandleAddress(cloned)
    second = N_VTestHandleAddress(cloned)
    assert first == second


def test_custom_nvector_slots_materialize_and_dispatch(sunctx):
    x = SlotsVector([1.0, 2.0], sunctx)
    y = SlotsVector([3.0, 4.0], sunctx)
    z = SlotsVector([0.0, 0.0], sunctx)

    assert N_VGetVectorID(x) == SUNDIALS_NVEC_CUSTOM
    N_VLinearSum(2.0, x, -1.0, y, z)
    assert_allclose(z.data, [-1.0, 0.0])


def test_custom_nvector_numpy_view(sunctx):
    x = ArrayVector([1.0, 2.0, 3.0], sunctx)
    view = N_VGetNumpyArray(x)
    view[:] = [4.0, 5.0, 6.0]
    assert_allclose(x.data, [4.0, 5.0, 6.0])


def test_custom_nvector_optional_fused_and_local_operations(sunctx):
    x = ArrayVector([1.0, 2.0], sunctx)
    y = ArrayVector([3.0, 4.0], sunctx)
    z = ArrayVector([0.0, 0.0], sunctx)

    assert N_VLinearCombination(2, np.array([2.0, -1.0]), [x, y], z) == SUN_SUCCESS
    assert_allclose(z.data, [-1.0, 0.0])
    dots = np.zeros(2, dtype=sunrealtype)
    assert N_VDotProdMulti(2, x, [y, z], dots) == SUN_SUCCESS
    assert_allclose(dots, [11.0, -1.0])
    assert N_VDotProdLocal(x, y) == pytest.approx(11.0)


def test_vector_array_reductions_write_failure_values(sunctx):
    x = FailingArrayReductionVector([1.0, 2.0], sunctx)
    w = FailingArrayReductionVector([1.0, 1.0], sunctx)

    norms = np.zeros(1, dtype=sunrealtype)
    with pytest.raises(ArithmeticError, match="array norm boom"):
        N_VWrmsNormVectorArray(1, [x], [w], norms)
    assert np.isposinf(norms[0])

    dots = np.zeros(1, dtype=sunrealtype)
    with pytest.raises(ArithmeticError, match="array dot boom"):
        N_VDotProdMulti(1, x, [x], dots)
    assert np.isnan(dots[0])


def test_custom_nvector_required_methods_are_validated(sunctx):
    x = IncompleteVector(sunctx)
    with pytest.raises(TypeError, match="missing required operations"):
        x.validate()
    with pytest.raises(TypeError, match="N_VGetVectorID"):
        N_VGetVectorID(x)
    assert not x._is_materialized()


def test_custom_nvector_required_error_lists_all_missing_operations(sunctx):
    errors = []
    SUNContext_PushErrHandler(
        sunctx, lambda line, function, file, message, code, data, context: errors.append(message)
    )
    with pytest.raises(TypeError):
        N_VGetVectorID(IncompleteVector(sunctx))
    assert "missing required operations" in errors[-1]
    assert "clone" in errors[-1]
    assert "constr_mask" in errors[-1]


def test_custom_nvector_rejects_incompatible_operands(sunctx):
    errors = []
    SUNContext_PushErrHandler(
        sunctx, lambda line, function, file, message, code, data, context: errors.append(message)
    )
    x = ArrayVector([1.0], sunctx)
    y = OtherArrayVector([2.0], sunctx)
    with pytest.raises(TypeError, match="different vector type"):
        N_VDotProd(x, y)
    assert "different vector type" in errors[-1]


def test_custom_nvector_allows_explicitly_compatible_operands(sunctx):
    x = CompatibleArrayVector([1.0], sunctx)
    y = OtherArrayVector([2.0], sunctx)
    assert N_VDotProd(x, y) == pytest.approx(2.0)


def test_custom_nvector_clone_rejects_different_context(sunctx):
    errors = []
    SUNContext_PushErrHandler(
        sunctx, lambda line, function, file, message, code, data, context: errors.append(message)
    )
    _, other_context = SUNContext_Create(SUN_COMM_NULL)
    x = WrongContextVector([1.0], sunctx, other_context)
    with pytest.raises(ValueError, match="same SUNContext"):
        N_VClone(x)
    assert "same SUNContext" in errors[-1]


def test_custom_nvector_integrates_with_arkode(sunctx):
    y = ArrayVector([1.0], sunctx)

    def rhs(t, state, derivative, _):
        assert isinstance(state, ArrayVector)
        assert isinstance(derivative, ArrayVector)
        derivative.data[:] = -state.data
        return SUN_SUCCESS

    arkode = ERKStepCreate(rhs, 0.0, y, sunctx)
    assert ARKodeSStolerances(arkode.get(), 1e-7, 1e-10) == ARK_SUCCESS
    status, tret = ARKodeEvolve(arkode.get(), 1.0, y, ARK_NORMAL)
    assert status == ARK_SUCCESS
    assert tret == pytest.approx(1.0)
    assert y.data[0] == pytest.approx(np.exp(-1.0), rel=2e-6)


def test_cvode_quadrature_init_without_nvspace(sunctx):
    """CVodeQuadInit accepts a custom vector without the deprecated N_VSpace."""
    y = ArrayVector([1.0], sunctx)
    yq = ArrayVector([0.0], sunctx)
    cvode = CVodeCreate(CV_BDF, sunctx)

    assert CVodeInit(cvode.get(), lambda t, y, yd, _: SUN_SUCCESS, 0.0, y) == CV_SUCCESS
    assert CVodeQuadInit(cvode.get(), lambda t, y, yqdot, _: SUN_SUCCESS, yq) == CV_SUCCESS


def test_ida_quadrature_init_without_nvspace(sunctx):
    """IDAQuadInit accepts a custom vector without the deprecated N_VSpace."""
    y = ArrayVector([1.0], sunctx)
    yp = ArrayVector([0.0], sunctx)
    yq = ArrayVector([0.0], sunctx)
    ida = IDACreate(sunctx)

    assert IDAInit(ida.get(), lambda t, y, yp, res, _: SUN_SUCCESS, 0.0, y, yp) == IDA_SUCCESS
    assert IDAQuadInit(ida.get(), lambda t, y, yp, yqdot, _: SUN_SUCCESS, yq) == IDA_SUCCESS


def test_void_operation_exception_reaches_cvode_caller(sunctx):
    y = FailingLinearSumVector([1.0], sunctx)

    def rhs(t, state, derivative, _):
        derivative.data[:] = -state.data
        return SUN_SUCCESS

    cvode = CVodeCreate(CV_BDF, sunctx)
    assert CVodeInit(cvode.get(), rhs, 0.0, y) == CV_SUCCESS
    assert CVodeSStolerances(cvode.get(), 1e-7, 1e-10) == CV_SUCCESS
    linear_solver = SUNLinSol_SPGMR(y, 0, 0, sunctx)
    assert CVodeSetLinearSolver(cvode.get(), linear_solver, None) == CV_SUCCESS
    with pytest.raises(LookupError, match="linear sum boom"):
        CVode(cvode.get(), 1.0, y, CV_NORMAL)

    # A later top-level solve starts with a clean pending-exception slot.
    healthy = ArrayVector([1.0], sunctx)
    cvode = CVodeCreate(CV_BDF, sunctx)
    assert CVodeInit(cvode.get(), rhs, 0.0, healthy) == CV_SUCCESS
    assert CVodeSStolerances(cvode.get(), 1e-7, 1e-10) == CV_SUCCESS
    linear_solver = SUNLinSol_SPGMR(healthy, 0, 0, sunctx)
    assert CVodeSetLinearSolver(cvode.get(), linear_solver, None) == CV_SUCCESS
    status, _ = CVode(cvode.get(), 0.1, healthy, CV_NORMAL)
    assert status == CV_SUCCESS


def test_reduction_exception_reaches_kinsol_caller(sunctx):
    y = FailingMaxNormVector([1.0], sunctx)
    scale = FailingMaxNormVector([1.0], sunctx)

    def residual(state, result, _):
        result.data[:] = state.data
        return SUN_SUCCESS

    kinsol = KINCreate(sunctx)
    assert KINInit(kinsol.get(), residual, y) == KIN_SUCCESS
    assert KINSetNumMaxIters(kinsol.get(), 1) == KIN_SUCCESS
    with pytest.raises(ArithmeticError, match="max norm boom"):
        KINSol(kinsol.get(), y, KIN_FP, scale, scale)
