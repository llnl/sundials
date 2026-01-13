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


def test_create_sparse_matrix(sunctx):
    rows, cols, nnz = 3, 3, 4
    A = SUNSparseMatrix(rows, cols, nnz, SUN_CSR_MAT, sunctx)
    assert A is not None
    assert SUNMatGetID(A) == SUNMATRIX_SPARSE
    # Check shape of data
    dataA = SUNSparseMatrix_Data(A)
    assert dataA.shape[0] == nnz
    idxA = SUNSparseMatrix_IndexValues(A)
    assert idxA.shape[0] == nnz
    ptrA = SUNSparseMatrix_IndexPointers(A)
    assert ptrA.shape[0] == 4


def test_clone_matrix(sunctx):
    rows, cols, nnz = 3, 3, 4
    A = SUNSparseMatrix(rows, cols, nnz, SUN_CSR_MAT, sunctx)
    B = SUNMatClone(A)
    assert B is not None
    assert SUNSparseMatrix_NNZ(A) == SUNSparseMatrix_NNZ(B)
    assert SUNSparseMatrix_NP(A) == SUNSparseMatrix_NP(B)


def test_zero_matrix(sunctx):
    rows, cols, nnz = 3, 3, 4
    A = SUNSparseMatrix(rows, cols, nnz, SUN_CSR_MAT, sunctx)
    ret = SUNMatZero(A)
    assert ret == SUN_SUCCESS
    dataA = SUNSparseMatrix_Data(A)
    assert_array_equal(dataA, np.zeros_like(dataA))


def test_copy_matrix(sunctx):
    rows, cols, nnz = 3, 3, 4
    A = SUNSparseMatrix(rows, cols, nnz, SUN_CSR_MAT, sunctx)
    B = SUNSparseMatrix(rows, cols, nnz, SUN_CSR_MAT, sunctx)
    dataA = SUNSparseMatrix_Data(A)
    idx_vals = SUNSparseMatrix_IndexValues(A)
    idx_ptrs = SUNSparseMatrix_IndexPointers(A)
    # CSR: row 0: col 0 (1.0), row 1: col 1 (1.0), row 2: col 2 (1.0), row 2: col 0 (2.0)
    dataA[:] = [1.0, 1.0, 1.0, 2.0]
    idx_vals[:] = [0, 1, 2, 0]
    idx_ptrs[:] = [0, 1, 2, 4]
    ret = SUNMatCopy(A, B)
    assert ret == SUN_SUCCESS
    dataB = SUNSparseMatrix_Data(B)
    idx_valsB = SUNSparseMatrix_IndexValues(B)
    idx_ptrsB = SUNSparseMatrix_IndexPointers(B)
    assert_allclose(dataB, dataA)
    assert_allclose(idx_valsB, idx_vals)
    assert_allclose(idx_ptrsB, idx_ptrs)


def test_scale_add_matrix(sunctx):
    rows, cols, nnz = 3, 3, 4
    A = SUNSparseMatrix(rows, cols, nnz, SUN_CSR_MAT, sunctx)
    B = SUNSparseMatrix(rows, cols, nnz, SUN_CSR_MAT, sunctx)
    dataA = SUNSparseMatrix_Data(A)
    dataB = SUNSparseMatrix_Data(B)
    idx_vals = SUNSparseMatrix_IndexValues(A)
    idx_ptrs = SUNSparseMatrix_IndexPointers(A)
    idx_valsB = SUNSparseMatrix_IndexValues(B)
    idx_ptrsB = SUNSparseMatrix_IndexPointers(B)
    # CSR: row 0: col 0 (1.0), row 1: col 1 (1.0), row 2: col 2 (1.0), row 2: col 0 (2.0)
    dataA[:] = [1.0, 1.0, 1.0, 2.0]
    dataB[:] = [2.0, 2.0, 2.0, 4.0]
    idx_vals[:] = [0, 1, 2, 0]
    idx_ptrs[:] = [0, 1, 2, 4]
    idx_valsB[:] = [0, 1, 2, 0]
    idx_ptrsB[:] = [0, 1, 2, 4]
    ret = SUNMatScaleAdd(3.0, A, B)
    assert ret == SUN_SUCCESS
    # 3*A + B = [3+2, 3+2, 3+2, 6+4] = [5, 5, 5, 10]
    assert_allclose(dataA, [5.0, 5.0, 5.0, 10.0])


def test_scale_add_identity(sunctx):
    rows, cols, nnz = 3, 3, 4
    A = SUNSparseMatrix(rows, cols, nnz, SUN_CSR_MAT, sunctx)
    dataA = SUNSparseMatrix_Data(A)
    idx_vals = SUNSparseMatrix_IndexValues(A)
    idx_ptrs = SUNSparseMatrix_IndexPointers(A)
    dataA[:] = [2.0, 2.0, 2.0, 2.0]
    idx_vals[:] = [0, 1, 2, 0]
    idx_ptrs[:] = [0, 1, 2, 4]
    ret = SUNMatScaleAddI(3.0, A)
    assert ret == SUN_SUCCESS
    # Diagonal elements should be 3*2+1=7, off-diagonal should be 3*2=6
    assert_allclose(dataA, [7.0, 7.0, 7.0, 6.0])


def test_matvec(sunctx):
    rows, cols, nnz = 3, 3, 4
    A = SUNSparseMatrix(rows, cols, nnz, SUN_CSR_MAT, sunctx)
    x = N_VNew_Serial(cols, sunctx)
    y = N_VNew_Serial(rows, sunctx)
    N_VConst(1.0, x)
    # Fill a simple sparse matrix: 3x3 identity with one extra off-diagonal
    dataA = SUNSparseMatrix_Data(A)
    idx_vals = SUNSparseMatrix_IndexValues(A)
    idx_ptrs = SUNSparseMatrix_IndexPointers(A)
    # CSR: row 0: col 0 (1.0), row 1: col 1 (1.0), row 2: col 2 (1.0), row 2: col 0 (2.0)
    dataA[:] = [1.0, 1.0, 1.0, 2.0]
    idx_vals[:] = [0, 1, 2, 0]
    idx_ptrs[:] = [0, 1, 2, 4]
    ret = SUNMatMatvec(A, x, y)
    assert ret == SUN_SUCCESS
    # y[0] = 1.0*x[0] = 1.0
    # y[1] = 1.0*x[1] = 1.0
    # y[2] = 1.0*x[2] + 2.0*x[0] = 1.0 + 2.0 = 3.0
    assert_allclose(N_VGetArrayPointer(y), [1.0, 1.0, 3.0])
