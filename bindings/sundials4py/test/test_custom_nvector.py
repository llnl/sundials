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


class IncompleteVector(CustomNVector):
    def get_length(self):
        return 1


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


def test_custom_nvector_required_methods_are_validated(sunctx):
    x = IncompleteVector(sunctx)
    with pytest.raises(TypeError, match="N_VGetVectorID"):
        N_VGetVectorID(x)
    assert not x._is_materialized()


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
