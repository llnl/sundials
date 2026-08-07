# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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

import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from fixtures import *
from sundials4py.core import *


def test_create_band_matrix(sunctx):
    rows, mu, ml = 4, 1, 1
    A = SUNBandMatrix(rows, mu, ml, sunctx)
    assert A is not None
    assert SUNMatGetID(A) == SUNMATRIX_BAND
    # Ensure the shape is being translated correctly
    dataA = SUNBandMatrix_Data(A)
    ldata = SUNBandMatrix_LData(A)
    assert dataA.shape[0] == ldata


def test_clone_matrix(sunctx):
    rows, mu, ml = 4, 1, 1
    A = SUNBandMatrix(rows, mu, ml, sunctx)
    B = SUNMatClone(A)
    assert B is not None
    assert SUNBandMatrix_LData(A) == SUNBandMatrix_LData(B)


def test_zero_matrix(sunctx):
    rows, mu, ml = 4, 1, 1
    A = SUNBandMatrix(rows, mu, ml, sunctx)
    assert A is not None
    ret = SUNMatZero(A)
    assert ret == SUN_SUCCESS
    dataA = SUNBandMatrix_Data(A)
    assert_array_equal(dataA, np.zeros_like(dataA))


def test_copy_matrix(sunctx):
    rows, mu, ml = 4, 1, 1
    A = SUNBandMatrix(rows, mu, ml, sunctx)
    B = SUNBandMatrix(rows, mu, ml, sunctx)
    smu = SUNBandMatrix_StoredUpperBandwidth(A)
    dataA = SUNBandMatrix_Data(A)
    dataA[smu - mu] = 1.0
    ret = SUNMatCopy(A, B)
    assert ret == SUN_SUCCESS
    dataB = SUNBandMatrix_Data(B)
    assert dataB[smu - mu] == 1.0


def test_scale_add_matrix(sunctx):
    rows, mu, ml = 4, 1, 1
    A = SUNBandMatrix(rows, mu, ml, sunctx)
    B = SUNBandMatrix(rows, mu, ml, sunctx)
    smu = SUNBandMatrix_StoredUpperBandwidth(A)
    dataA = SUNBandMatrix_Data(A)
    dataB = SUNBandMatrix_Data(B)
    dataA[smu - mu : smu + ml] = 1.0  # column 0 set to 1.0
    dataB[smu - mu : smu + ml] = 2.0
    ret = SUNMatScaleAdd(3.0, A, B)
    assert ret == SUN_SUCCESS
    # A should now be 3*A + B = 3*1 + 2 = 5
    assert_allclose(dataA[smu - mu : smu + ml], 5.0)


def test_scale_add_identity(sunctx):
    rows, mu, ml = 4, 1, 1
    A = SUNBandMatrix(rows, mu, ml, sunctx)
    ldim = SUNBandMatrix_LDim(A)
    smu = SUNBandMatrix_StoredUpperBandwidth(A)
    dataA = SUNBandMatrix_Data(A)
    ret = SUNMatScaleAddI(0.0, A)
    assert ret == SUN_SUCCESS
    # A should now be I
    diag = np.array([dataA[smu + i * ldim] for i in range(rows)], dtype=sunrealtype)
    assert_allclose(diag, 1.0)


def test_matvec(sunctx):
    rows, mu, ml = 4, 1, 1
    A = SUNBandMatrix(rows, mu, ml, sunctx)
    x = N_VNew_Serial(rows, sunctx)
    y = N_VNew_Serial(rows, sunctx)

    N_VConst(1.0, x)

    # Fill band matrix data for a simple 4x4 banded matrix
    # [3 2 0 0]
    # [1 3 2 0]
    # [0 1 3 2]
    # [0 0 1 3]
    dataA = SUNBandMatrix_Data(A)
    ldim = SUNBandMatrix_LDim(A)
    smu = SUNBandMatrix_StoredUpperBandwidth(A)
    for j in range(rows):
        # Diagonal
        dataA[smu + j * ldim] = 3.0
        # Lower diagonal
        if j > 0:
            dataA[smu - 1 + j * ldim] = 2.0
        # Upper diagonal
        if j < rows - 1:
            dataA[smu + 1 + j * ldim] = 1.0

    ret = SUNMatMatvec(A, x, y)
    assert ret == SUN_SUCCESS

    assert_allclose(N_VGetArrayPointer(y), [5.0, 6.0, 6.0, 4.0])
