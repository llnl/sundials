# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2025, Lawrence Livermore National Security,
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

import pytest
import numpy as np
from numpy.testing import assert_allclose
from fixtures import *
from sundials4py.core import *


def test_create_manyvector(sunctx):
    y = N_VNew_Serial(5, sunctx)
    z = N_VNew_Serial(5, sunctx)

    N_VConst(1.0, y)
    N_VConst(2.0, z)

    yz = N_VNew_ManyVector(2, [y, z], sunctx)

    yarr = N_VGetArrayPointer(N_VGetSubvector_ManyVector(yz, 0))
    assert_allclose(N_VGetArrayPointer(y), 1.0)

    zarr = N_VGetArrayPointer(N_VGetSubvector_ManyVector(yz, 1))
    assert_allclose(N_VGetArrayPointer(z), 2.0)

    N_VConst(3.0, yz)
    assert_allclose(3.0, yarr)
    assert_allclose(3.0, zarr)


@pytest.mark.parametrize("vector_type", ["serial"])
def test_create_nvector(vector_type, sunctx):
    if vector_type == "serial":
        nvec = N_VNew_Serial(5, sunctx)
    else:
        raise ValueError("Unknown vector type")
    assert nvec is not None

    arr = N_VGetArrayPointer(nvec)
    assert arr.shape[0] == 5

    arr[:] = np.array([5.0, 4.0, 3.0, 2.0, 1.0], dtype=sunrealtype)
    assert_allclose(N_VGetArrayPointer(nvec), [5.0, 4.0, 3.0, 2.0, 1.0])

    N_VConst(2.0, nvec)
    assert_allclose(arr, 2.0)


@pytest.mark.parametrize("vector_type", ["serial"])
def test_make_nvector(vector_type, sunctx):
    arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=sunrealtype)
    if vector_type == "serial":
        nvec = N_VMake_Serial(5, arr, sunctx)
    else:
        raise ValueError("Unknown vector type")
    assert nvec is not None

    assert_allclose(N_VGetArrayPointer(nvec), arr)

    N_VConst(2.0, nvec)
    assert_allclose(arr, 2.0)

    arr[:] = np.array([5.0, 4.0, 3.0, 2.0, 1.0], dtype=sunrealtype)
    assert_allclose(N_VGetArrayPointer(nvec), [5.0, 4.0, 3.0, 2.0, 1.0])


# Test an operation that involves vector arrays
@pytest.mark.parametrize("vector_type", ["serial"])
def test_nvlinearcombination(vector_type, sunctx):
    if vector_type == "serial":
        nvec1 = N_VNew_Serial(5, sunctx)
        nvec2 = N_VNew_Serial(5, sunctx)
    else:
        raise ValueError("Unknown vector type")

    arr1 = N_VGetArrayPointer(nvec1)
    arr1[:] = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=sunrealtype)

    arr2 = N_VGetArrayPointer(nvec2)
    arr2[:] = np.array([10.0, 20.0, 30.0, 40.0, 50.0], dtype=sunrealtype)

    c = np.array([1.0, 0.1], dtype=sunrealtype)
    X = [nvec1, nvec2]

    z = N_VNew_Serial(5, sunctx)
    N_VConst(0.0, z)

    N_VLinearCombination(2, c, X, z)

    assert_allclose(N_VGetArrayPointer(z), [2.0, 4.0, 6.0, 8.0, 10.0])


def test_nvscaleaddmultivectorarray_serial(sunctx):
    nvec = 2
    nsum = 2
    length = 3

    # c_1d shape (nsum,)
    c_1d = np.array([2.0, 3.0], dtype=sunrealtype)

    # X_1d shape (nvec,)
    X_1d = [N_VNew_Serial(length, sunctx) for _ in range(nvec)]

    for i, x in enumerate(X_1d):
        N_VConst(float(i + 1), x)

    # Y_2d shape (nsum, nvec)
    Y_2d = [[N_VNew_Serial(length, sunctx) for _ in range(nvec)] for _ in range(nsum)]
    for s in range(nsum):
        for v in range(nvec):
            N_VConst(float((s + 1) * 10 + v), Y_2d[s][v])

    # Z_2d shape (nsum, nvec)
    Z_2d = [[N_VNew_Serial(length, sunctx) for _ in range(nvec)] for _ in range(nsum)]

    err = N_VScaleAddMultiVectorArray(nvec, nsum, c_1d, X_1d, Y_2d, Z_2d)
    assert err == SUN_SUCCESS

    # Check Z_2d[s][v] = c_1d[s] * X_1d[v] + Y_2d[s][v]
    for s in range(nsum):
        for v in range(nvec):
            expected = c_1d[s] * N_VGetArrayPointer(X_1d[v]) + N_VGetArrayPointer(Y_2d[s][v])
            actual = N_VGetArrayPointer(Z_2d[s][v])
            assert_allclose(actual, expected)


def test_nvlinearcombinationvectorarray_serial(sunctx):
    nvec = 2
    nsum = 2
    length = 3

    # c_1d shape (nsum,)
    c_1d = np.array([2.0, 3.0], dtype=sunrealtype)

    # X_2d shape (nsum, nvec)
    X_2d = []
    for s in range(nsum):
        row = []
        for v in range(nvec):
            x = N_VNew_Serial(length, sunctx)
            N_VConst(float((s + 1) * 10 + v), x)
            row.append(x)
        X_2d.append(row)

    # Z_1d shape (nvec,)
    Z_1d = [N_VNew_Serial(length, sunctx) for _ in range(nvec)]

    err = N_VLinearCombinationVectorArray(nvec, nsum, c_1d, X_2d, Z_1d)
    assert err == SUN_SUCCESS

    # Check Z_1d[v] = sum_s c_1d[s] * X_2d[s][v]
    for v in range(nvec):
        expected = sum(c_1d[s] * N_VGetArrayPointer(X_2d[s][v]) for s in range(nsum))
        actual = N_VGetArrayPointer(Z_1d[v])
        assert_allclose(actual, expected)
