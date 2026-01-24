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
from numpy.testing import assert_allclose, assert_array_equal
from fixtures import *
from sundials4py.core import *


def test_create_dense_matrix(sunctx):
    rows, cols = 3, 2
    A = SUNDenseMatrix(rows, cols, sunctx)
    assert A is not None
    assert SUNMatGetID(A) == SUNMATRIX_DENSE
    # Ensure the shape is being translated correctly
    dataA_shape = np.shape(SUNDenseMatrix_Data(A))
    assert dataA_shape[0] == rows and dataA_shape[1] == cols


def test_clone_matrix(sunctx):
    A = SUNDenseMatrix(2, 2, sunctx)
    B = SUNMatClone(A)
    assert B is not None
    assert np.shape(SUNDenseMatrix_Data(A)) == np.shape(SUNDenseMatrix_Data(B))


def test_zero_matrix(sunctx):
    A = SUNDenseMatrix(2, 2, sunctx)
    ret = SUNMatZero(A)
    assert ret == SUN_SUCCESS
    dataA = SUNDenseMatrix_Data(A)
    assert_array_equal(dataA, np.zeros_like(dataA))


def test_copy_matrix(sunctx):
    A = SUNDenseMatrix(2, 2, sunctx)
    B = SUNDenseMatrix(2, 2, sunctx)
    # Set some values in A
    dataA = SUNDenseMatrix_Data(A)
    dataA[0, 0] = 1.0
    ret = SUNMatCopy(A, B)
    assert ret == SUN_SUCCESS
    dataB = SUNDenseMatrix_Data(B)
    assert dataB[0, 0] == 1.0


def test_scale_add_matrix(sunctx):
    A = SUNDenseMatrix(2, 2, sunctx)
    B = SUNDenseMatrix(2, 2, sunctx)
    dataA = SUNDenseMatrix_Data(A)
    dataB = SUNDenseMatrix_Data(B)
    dataA[:, :] = 1.0
    dataB[:, :] = 2.0
    ret = SUNMatScaleAdd(3.0, A, B)
    assert ret == SUN_SUCCESS
    # A should now be 3*A + B = 3*1 + 2 = 5
    assert_allclose(dataA, 5.0)


def test_scale_add_identity(sunctx):
    A = SUNDenseMatrix(2, 2, sunctx)
    dataA = SUNDenseMatrix_Data(A)
    dataA[:, :] = 2.0
    ret = SUNMatScaleAddI(3.0, A)
    assert ret == SUN_SUCCESS
    # A should now be 3*A + I
    expected = np.eye(2) + 6.0
    assert_allclose(dataA, expected)


def test_matvec(sunctx, nvec):
    A = SUNDenseMatrix(2, 2, sunctx)
    dataA = SUNDenseMatrix_Data(A)
    dataA[:, :] = [[1.0, 2.0], [3.0, 4.0]]
    x = N_VNew_Serial(2, sunctx)
    y = N_VNew_Serial(2, sunctx)
    N_VConst(1.0, x)
    ret = SUNMatMatvec(A, x, y)
    assert ret == SUN_SUCCESS
    arr = N_VGetArrayPointer(y)
    assert_allclose(arr, [3.0, 7.0])
