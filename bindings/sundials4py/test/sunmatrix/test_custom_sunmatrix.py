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

import numpy as np
import pytest
from fixtures import *
from numpy.testing import assert_allclose
from sundials4py.core import *


class DiagonalMatrix(CustomSUNMatrix):
    # Minimal concrete matrix used to verify that SUNMatrix operations dispatch
    # from the native vtable back into Python methods.
    def __init__(self, diagonal, sunctx):
        self.diagonal = np.asarray(diagonal, dtype=sunrealtype).copy()
        self.calls = {"clone": 0, "zero": 0, "copy": 0, "scaleadd": 0, "scaleaddi": 0, "matvec": 0}
        super().__init__(sunctx)

    def clone(self):
        self.calls["clone"] += 1
        return type(self)(self.diagonal.copy(), self.sunctx)

    def zero(self):
        self.calls["zero"] += 1
        self.diagonal[:] = 0
        return SUN_SUCCESS

    def copy(self, dst):
        self.calls["copy"] += 1
        dst.diagonal[:] = self.diagonal
        return SUN_SUCCESS

    def scaleadd(self, c, other):
        self.calls["scaleadd"] += 1
        self.diagonal[:] = c * self.diagonal + other.diagonal
        return SUN_SUCCESS

    def scaleaddi(self, c):
        self.calls["scaleaddi"] += 1
        self.diagonal[:] = c * self.diagonal + 1
        return SUN_SUCCESS

    def matvec(self, x, y):
        self.calls["matvec"] += 1
        N_VGetArrayPointer(y)[:] = self.diagonal * N_VGetArrayPointer(x)
        return SUN_SUCCESS


class SetupDiagonalMatrix(DiagonalMatrix):
    # Optional methods should appear in the native ops table only when a
    # subclass overrides the base method.
    def __init__(self, diagonal, sunctx):
        self.matvecsetup_calls = 0
        super().__init__(diagonal, sunctx)

    def matvecsetup(self):
        self.matvecsetup_calls += 1
        return SUN_SUCCESS


class HermitianDiagonalMatrix(DiagonalMatrix):
    def __init__(self, diagonal, sunctx):
        self.hermitian_calls = 0
        super().__init__(diagonal, sunctx)

    def hermitian_transpose_matvec(self, x, y):
        self.hermitian_calls += 1
        N_VGetArrayPointer(y)[:] = self.diagonal * N_VGetArrayPointer(x)
        return SUN_SUCCESS


class IncompleteMatrix(CustomSUNMatrix):
    # Missing required overrides must fail before a native handle is cached.
    def __init__(self, sunctx):
        super().__init__(sunctx)


def test_custom_sunmatrix_converts_for_handwritten_bindings(sunctx):
    # Purpose:
    # Custom sunmatrix converts for handwritten bindings.
    A = DiagonalMatrix([1.0, 2.0], sunctx)

    assert A._materialization_count() == 0
    assert SUNMatGetID(A) == SUNMATRIX_CUSTOM
    assert A._materialization_count() == 1

    assert SUNMatZero(A) == SUN_SUCCESS
    assert A._materialization_count() == 1
    assert A.calls["zero"] == 1
    assert_allclose(A.diagonal, [0.0, 0.0])


def test_custom_sunmatrix_required_methods_are_validated(sunctx):
    # Purpose:
    # Custom sunmatrix required methods are validated.
    A = IncompleteMatrix(sunctx)

    with pytest.raises(TypeError, match="SUNMatGetID"):
        SUNMatGetID(A)

    assert A._materialization_count() == 0


def test_custom_sunmatrix_copy_scaleadd_and_scaleaddi(sunctx):
    # Purpose:
    # Custom sunmatrix copy scaleadd and scaleaddi.
    A = DiagonalMatrix([1.0, 2.0], sunctx)
    B = DiagonalMatrix([3.0, 4.0], sunctx)

    assert SUNMatCopy(A, B) == SUN_SUCCESS
    assert_allclose(B.diagonal, [1.0, 2.0])
    assert A.calls["copy"] == 1

    assert SUNMatScaleAdd(10.0, A, B) == SUN_SUCCESS
    assert_allclose(A.diagonal, [11.0, 22.0])
    assert A.calls["scaleadd"] == 1

    assert SUNMatScaleAddI(2.0, A) == SUN_SUCCESS
    assert_allclose(A.diagonal, [23.0, 45.0])
    assert A.calls["scaleaddi"] == 1


def test_custom_sunmatrix_matvec(sunctx):
    # Purpose:
    # Custom sunmatrix matvec.
    A = DiagonalMatrix([2.0, 3.0], sunctx)
    x = N_VNew_Serial(2, sunctx)
    y = N_VNew_Serial(2, sunctx)
    N_VGetArrayPointer(x)[:] = [4.0, 5.0]

    assert SUNMatMatvec(A, x, y) == SUN_SUCCESS
    assert_allclose(N_VGetArrayPointer(y), [8.0, 15.0])
    assert A.calls["matvec"] == 1


def test_custom_sunmatrix_clone_uses_python_clone(sunctx):
    # Purpose:
    # Custom sunmatrix clone uses python clone.
    A = DiagonalMatrix([1.0, 2.0], sunctx)
    B = SUNMatClone(A)

    assert B is not None
    assert A.calls["clone"] == 1
    assert SUNMatGetID(B) == SUNMATRIX_CUSTOM
    assert SUNMatScaleAddI(3.0, B) == SUN_SUCCESS


def test_custom_sunmatrix_optional_matvecsetup_detection(sunctx):
    # Purpose:
    # Custom sunmatrix optional matvecsetup detection.
    A = DiagonalMatrix([1.0, 2.0], sunctx)
    B = SetupDiagonalMatrix([1.0, 2.0], sunctx)

    assert SUNMatMatvecSetup(A) == SUN_SUCCESS
    assert SUNMatMatvecSetup(B) == SUN_SUCCESS
    assert B.matvecsetup_calls == 1


def test_optional_vtable_is_fixed_when_handle_is_materialized(sunctx, monkeypatch):
    # Purpose:
    # Optional vtable is fixed when handle is materialized.
    A = DiagonalMatrix([1.0, 2.0], sunctx)
    B = SetupDiagonalMatrix([1.0, 2.0], sunctx)
    assert SUNMatGetID(A) == SUNMATRIX_CUSTOM
    assert SUNMatGetID(B) == SUNMATRIX_CUSTOM

    calls = {"added": 0, "replacement": 0}

    def added(self):
        calls["added"] += 1
        return SUN_SUCCESS

    def replacement(self):
        calls["replacement"] += 1
        return SUN_SUCCESS

    monkeypatch.setattr(DiagonalMatrix, "matvecsetup", added)
    monkeypatch.setattr(SetupDiagonalMatrix, "matvecsetup", replacement)

    # A had no optional slot at materialization, so adding a method later does
    # not install one. B retains its slot, whose trampoline performs normal
    # dynamic Python method lookup and therefore sees the replacement.
    assert SUNMatMatvecSetup(A) == SUN_SUCCESS
    assert SUNMatMatvecSetup(B) == SUN_SUCCESS
    assert calls == {"added": 0, "replacement": 1}


def test_custom_sunmatrix_optional_hermitian_matvec_detection(sunctx):
    # Purpose:
    # Custom sunmatrix optional hermitian matvec detection.
    A = DiagonalMatrix([1.0, 2.0], sunctx)
    B = HermitianDiagonalMatrix([1.0, 2.0], sunctx)
    x = N_VNew_Serial(2, sunctx)
    y = N_VNew_Serial(2, sunctx)
    N_VGetArrayPointer(x)[:] = [3.0, 4.0]

    assert SUNMatHermitianTransposeVec(A, x, y) == SUN_ERR_NOT_IMPLEMENTED
    assert SUNMatHermitianTransposeVec(B, x, y) == SUN_SUCCESS
    assert_allclose(N_VGetArrayPointer(y), [3.0, 8.0])
    assert B.hermitian_calls == 1


def test_custom_sunmatrix_rejects_mismatched_custom_operand_types(sunctx):
    # Purpose:
    # Custom sunmatrix rejects mismatched custom operand types.
    A = DiagonalMatrix([1.0, 2.0], sunctx)
    B = SetupDiagonalMatrix([3.0, 4.0], sunctx)

    assert SUNMatCopy(A, B) == SUN_ERR_EXT_FAIL


def test_native_sunmatrix_conversion_still_works(sunctx):
    # Purpose:
    # Native sunmatrix conversion still works.
    A = SUNDenseMatrix(2, 2, sunctx)

    assert SUNMatGetID(A) == SUNMATRIX_DENSE
    assert SUNMatZero(A) == SUN_SUCCESS


def test_custom_sunmatrix_converts_for_generated_optional_bindings(sunctx):
    # Purpose:
    # Custom sunmatrix converts for generated optional bindings.
    x = N_VNew_Serial(2, sunctx)
    LS = SUNLinSol_PCG(x, SUN_PREC_NONE, 0, sunctx)
    A = DiagonalMatrix([1.0, 2.0], sunctx)

    assert SUNLinSolSetup(LS, A) == SUN_SUCCESS
    assert A._materialization_count() == 1


def test_generated_optional_bindings_still_accept_none(sunctx):
    # Purpose:
    # Generated optional bindings still accept none.
    x = N_VNew_Serial(2, sunctx)
    LS = SUNLinSol_PCG(x, SUN_PREC_NONE, 0, sunctx)

    assert SUNLinSolSetup(LS, None) == SUN_SUCCESS
