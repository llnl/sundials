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

import gc
import sys

import pytest
from fixtures import *
from numpy.testing import assert_allclose
from sundials4py.core import *


class CopyLinearSolver(CustomSUNLinearSolver):
    # Simple concrete solver: solve copies b into x so the tests can focus on
    # transparent conversion and vtable routing rather than numerical details.
    def __init__(self, sunctx, solver_type=SUNLINEARSOLVER_ITERATIVE):
        self.calls = {
            "initialize": 0,
            "setup": 0,
            "solve": 0,
            "set_scaling_vectors": 0,
            "set_zero_guess": 0,
        }
        self.zero_guess = False
        super().__init__(sunctx, solver_type)

    def initialize(self):
        self.calls["initialize"] += 1
        return SUN_SUCCESS

    def setup(self, A):
        self.calls["setup"] += 1
        self.last_setup_matrix = A
        return SUN_SUCCESS

    def solve(self, A, x, b, tol):
        self.calls["solve"] += 1
        self.last_solve_matrix = A
        self.last_tol = tol
        N_VGetArrayPointer(x)[:] = N_VGetArrayPointer(b)
        return SUN_SUCCESS

    def set_scaling_vectors(self, s1, s2):
        self.calls["set_scaling_vectors"] += 1
        self.last_scaling_vectors = (s1, s2)
        return SUN_SUCCESS

    def set_zero_guess(self, onoff):
        self.calls["set_zero_guess"] += 1
        self.zero_guess = bool(onoff)
        return SUN_SUCCESS

    def num_iters(self):
        return 3

    def res_norm(self):
        return 0.25


class CallbackLinearSolver(CopyLinearSolver):
    # Stores wrapped native callbacks so tests can confirm they are delivered as
    # ordinary Python callables to custom solver implementations.
    def __init__(self, sunctx):
        self.atimes = None
        self.psetup = None
        self.psolve = None
        super().__init__(sunctx)

    def set_atimes(self, atimes):
        self.atimes = atimes
        return SUN_SUCCESS

    def set_preconditioner(self, psetup, psolve):
        self.psetup = psetup
        self.psolve = psolve
        return SUN_SUCCESS


class IncompleteLinearSolver(CustomSUNLinearSolver):
    # Omitting solve() exercises the required-method validation path.
    def __init__(self, sunctx):
        super().__init__(sunctx, SUNLINEARSOLVER_DIRECT)


class ResidualLinearSolver(CopyLinearSolver):
    def __init__(self, sunctx):
        self.residual_refs = []
        super().__init__(sunctx)

    def resid(self):
        residual = N_VNew_Serial(1, self.sunctx)
        self.residual_refs.append(residual)
        return residual


def refcount(obj):
    # Keep pytest's assertion rewriting from retaining the object expression
    # while sys.getrefcount() measures it.
    return sys.getrefcount(obj)


def test_custom_sunlinearsolver_type_id_and_lazy_materialization(sunctx):
    # Purpose:
    # Custom sunlinearsolver type id and lazy materialization.
    LS = CopyLinearSolver(sunctx, SUNLINEARSOLVER_MATRIX_ITERATIVE)

    assert LS._materialization_count() == 0
    assert SUNLinSolGetType(LS) == SUNLINEARSOLVER_MATRIX_ITERATIVE
    assert SUNLinSolGetID(LS) == SUNLINEARSOLVER_CUSTOM
    assert LS._materialization_count() == 1


def test_custom_sunlinearsolver_required_solve_is_validated(sunctx):
    # Purpose:
    # Custom sunlinearsolver required solve is validated.
    LS = IncompleteLinearSolver(sunctx)

    with pytest.raises(TypeError, match="SUNLinSolGetType"):
        SUNLinSolGetType(LS)

    assert LS._materialization_count() == 0


def test_custom_sunlinearsolver_initialize_setup_and_solve(sunctx):
    # Purpose:
    # Custom sunlinearsolver initialize setup and solve.
    LS = CopyLinearSolver(sunctx)
    x = N_VNew_Serial(2, sunctx)
    b = N_VNew_Serial(2, sunctx)
    N_VGetArrayPointer(b)[:] = [2.0, 5.0]

    assert SUNLinSolInitialize(LS) == SUN_SUCCESS
    assert SUNLinSolSetup(LS, None) == SUN_SUCCESS
    assert SUNLinSolSolve(LS, None, x, b, 1.0e-8) == SUN_SUCCESS

    assert LS.calls["initialize"] == 1
    assert LS.calls["setup"] == 1
    assert LS.calls["solve"] == 1
    assert LS.last_setup_matrix is None
    assert LS.last_solve_matrix is None
    assert LS.last_tol == 1.0e-8
    assert_allclose(N_VGetArrayPointer(x), [2.0, 5.0])


def test_custom_sunlinearsolver_optional_methods(sunctx):
    # Purpose:
    # Custom sunlinearsolver optional methods.
    LS = CopyLinearSolver(sunctx)
    s1 = N_VNew_Serial(2, sunctx)
    s2 = N_VNew_Serial(2, sunctx)

    assert SUNLinSolSetScalingVectors(LS, s1, s2) == SUN_SUCCESS
    assert SUNLinSolSetZeroGuess(LS, True) == SUN_SUCCESS
    assert SUNLinSolNumIters(LS) == 3
    assert SUNLinSolResNorm(LS) == 0.25
    assert LS.calls["set_scaling_vectors"] == 1
    assert LS.calls["set_zero_guess"] == 1
    assert LS.zero_guess


def test_custom_sunlinearsolver_set_atimes_receives_python_adapter(sunctx):
    # Purpose:
    # Custom sunlinearsolver set atimes receives python adapter.
    LS = CallbackLinearSolver(sunctx)
    x = N_VNew_Serial(2, sunctx)
    y = N_VNew_Serial(2, sunctx)
    N_VGetArrayPointer(x)[:] = [1.0, 4.0]

    called = {"atimes": 0}

    def atimes(_, x, y):
        called["atimes"] += 1
        N_VGetArrayPointer(y)[:] = 2 * N_VGetArrayPointer(x)
        return SUN_SUCCESS

    assert SUNLinSolSetATimes(LS, atimes) == SUN_SUCCESS
    assert LS.atimes is not None
    assert LS.atimes(x, y) == SUN_SUCCESS
    assert called["atimes"] == 1
    assert_allclose(N_VGetArrayPointer(y), [2.0, 8.0])

    assert SUNLinSolSetATimes(LS, None) == SUN_SUCCESS
    assert LS.atimes is None


def test_custom_sunlinearsolver_atimes_adapters_are_revoked(sunctx):
    # Purpose:
    # Custom sunlinearsolver atimes adapters are revoked.
    LS = CallbackLinearSolver(sunctx)
    x = N_VNew_Serial(1, sunctx)
    y = N_VNew_Serial(1, sunctx)

    def atimes1(_, x, y):
        return 1

    def atimes2(_, x, y):
        return 2

    assert SUNLinSolSetATimes(LS, atimes1) == SUN_SUCCESS
    old = LS.atimes
    assert SUNLinSolSetATimes(LS, atimes2) == SUN_SUCCESS
    with pytest.raises(RuntimeError, match="no longer valid"):
        old(x, y)
    current = LS.atimes
    assert current(x, y) == 2

    assert SUNLinSolSetATimes(LS, None) == SUN_SUCCESS
    with pytest.raises(RuntimeError, match="no longer valid"):
        current(x, y)


def test_custom_sunlinearsolver_set_preconditioner_receives_python_adapters(sunctx):
    # Purpose:
    # Custom sunlinearsolver set preconditioner receives python adapters.
    LS = CallbackLinearSolver(sunctx)
    r = N_VNew_Serial(2, sunctx)
    z = N_VNew_Serial(2, sunctx)
    N_VGetArrayPointer(r)[:] = [3.0, 5.0]
    called = {"psetup": 0, "psolve": 0}

    def psetup(_):
        called["psetup"] += 1
        return SUN_SUCCESS

    def psolve(_, r, z, tol, lr):
        called["psolve"] += 1
        N_VGetArrayPointer(z)[:] = N_VGetArrayPointer(r) / tol
        return lr

    assert SUNLinSolSetPreconditioner(LS, psetup, psolve) == SUN_SUCCESS
    assert LS.psetup() == SUN_SUCCESS
    assert LS.psolve(r, z, 2.0, 1) == 1
    assert called == {"psetup": 1, "psolve": 1}
    assert_allclose(N_VGetArrayPointer(z), [1.5, 2.5])

    assert SUNLinSolSetPreconditioner(LS, None, None) == SUN_SUCCESS
    assert LS.psetup is None
    assert LS.psolve is None


@pytest.mark.parametrize("with_setup", [False, True])
@pytest.mark.parametrize("with_solve", [False, True])
def test_custom_sunlinearsolver_preconditioner_nullability(sunctx, with_setup, with_solve):
    # Purpose:
    # Custom sunlinearsolver preconditioner nullability.
    LS = CallbackLinearSolver(sunctx)

    def psetup(_):
        return SUN_SUCCESS

    def psolve(_, r, z, tol, lr):
        return SUN_SUCCESS

    assert (
        SUNLinSolSetPreconditioner(
            LS, psetup if with_setup else None, psolve if with_solve else None
        )
        == SUN_SUCCESS
    )
    assert callable(LS.psetup) is with_setup
    assert callable(LS.psolve) is with_solve


def test_custom_sunlinearsolver_preconditioner_adapters_are_revoked(sunctx):
    # Purpose:
    # Custom sunlinearsolver preconditioner adapters are revoked.
    LS = CallbackLinearSolver(sunctx)
    r = N_VNew_Serial(1, sunctx)
    z = N_VNew_Serial(1, sunctx)

    def psetup(_):
        return SUN_SUCCESS

    def psolve(_, r, z, tol, lr):
        return SUN_SUCCESS

    assert SUNLinSolSetPreconditioner(LS, psetup, psolve) == SUN_SUCCESS
    old_setup, old_solve = LS.psetup, LS.psolve
    assert SUNLinSolSetPreconditioner(LS, psetup, psolve) == SUN_SUCCESS
    with pytest.raises(RuntimeError, match="no longer valid"):
        old_setup()
    with pytest.raises(RuntimeError, match="no longer valid"):
        old_solve(r, z, 1.0, 1)

    current_setup, current_solve = LS.psetup, LS.psolve
    assert SUNLinSolSetPreconditioner(LS, None, None) == SUN_SUCCESS
    with pytest.raises(RuntimeError, match="no longer valid"):
        current_setup()
    with pytest.raises(RuntimeError, match="no longer valid"):
        current_solve(r, z, 1.0, 1)


def test_custom_sunlinearsolver_resid_retains_and_replaces_owner(sunctx):
    # Purpose:
    # Custom sunlinearsolver resid retains and replaces owner.
    LS = ResidualLinearSolver(sunctx)

    first = SUNLinSolResid(LS)
    del first
    gc.collect()
    first_owned_count = refcount(LS.residual_refs[0])

    second = SUNLinSolResid(LS)
    del second
    gc.collect()
    first_released_count = refcount(LS.residual_refs[0])
    assert first_released_count == first_owned_count - 1
    second_owned_count = refcount(LS.residual_refs[1])

    refs = LS.residual_refs
    del LS
    gc.collect()
    second_released_count = refcount(refs[1])
    assert second_released_count == second_owned_count - 1


def test_native_sunlinearsolver_conversion_still_works(sunctx):
    # Purpose:
    # Native sunlinearsolver conversion still works.
    x = N_VNew_Serial(2, sunctx)
    LS = SUNLinSol_PCG(x, SUN_PREC_NONE, 0, sunctx)

    assert SUNLinSolGetID(LS) == SUNLINEARSOLVER_PCG
    assert SUNLinSolSetup(LS, None) == SUN_SUCCESS
