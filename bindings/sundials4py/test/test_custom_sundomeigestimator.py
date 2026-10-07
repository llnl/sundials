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

import pytest
from fixtures import *
from sundials4py.arkode import *
from sundials4py.core import *


class Estimator(CustomSUNDomEigEstimator):
    def __init__(self, sunctx):
        self.max_iters = None
        self.atimes = None
        self.rhs = None
        self.rhs_point = None
        self.estimate_calls = 0
        super().__init__(sunctx)

    def set_atimes(self, fn):
        self.atimes = fn
        return SUN_SUCCESS

    def set_rhs(self, fn):
        self.rhs = fn
        return SUN_SUCCESS

    def set_max_iters(self, max_iters):
        self.max_iters = max_iters
        return SUN_SUCCESS

    def set_rhs_linearization_point(self, t, vector):
        self.rhs_point = (t, vector)
        return SUN_SUCCESS

    def estimate(self):
        self.estimate_calls += 1
        return SUN_SUCCESS, 4.0, -0.5

    def get_num_iters(self):
        return SUN_SUCCESS, 3


class IncompleteEstimator(CustomSUNDomEigEstimator):
    pass


def test_custom_domeigestimator_dispatch(sunctx):
    estimator = Estimator(sunctx)
    vector = N_VNew_Serial(1, sunctx)

    assert SUNDomEigEstimator_SetMaxIters(estimator, 12) == SUN_SUCCESS
    assert SUNDomEigEstimator_SetRhsLinearizationPoint(estimator, 2.0, vector) == SUN_SUCCESS
    assert estimator.max_iters == 12
    assert estimator.rhs_point == (2.0, vector)
    assert SUNDomEigEstimator_Estimate(estimator) == (SUN_SUCCESS, 4.0, -0.5)
    assert SUNDomEigEstimator_GetNumIters(estimator) == (SUN_SUCCESS, 3)


def test_custom_domeigestimator_callback_adapters(sunctx, nvec):
    estimator = Estimator(sunctx)
    calls = []

    def atimes(_, x, y):
        calls.append((x, y))
        return SUN_SUCCESS

    assert SUNDomEigEstimator_SetATimes(estimator, atimes) == SUN_SUCCESS
    assert estimator.atimes(nvec, nvec) == SUN_SUCCESS
    assert len(calls) == 1

    old_atimes = estimator.atimes
    assert SUNDomEigEstimator_SetATimes(estimator, atimes) == SUN_SUCCESS
    with pytest.raises(RuntimeError, match="no longer valid"):
        old_atimes(nvec, nvec)

    def rhs(t, y, ydot, _):
        calls.append((t, y, ydot))
        return SUN_SUCCESS

    assert SUNDomEigEstimator_SetRhs(estimator, rhs) == SUN_SUCCESS
    assert estimator.rhs(1.25, nvec, nvec) == SUN_SUCCESS
    assert calls[-1][0] == pytest.approx(1.25)


def test_custom_domeigestimator_required_method_is_validated(sunctx):
    estimator = IncompleteEstimator(sunctx)
    with pytest.raises(TypeError, match="SUNDomEigEstimator_Estimate"):
        SUNDomEigEstimator_Estimate(estimator)
    assert not estimator._is_materialized()


def test_custom_domeigestimator_integrates_with_lsrkstep(sunctx):
    class NegativeRealEstimator(Estimator):
        def estimate(self):
            self.estimate_calls += 1
            return SUN_SUCCESS, -10.0, 0.0

    y = N_VNew_Serial(1, sunctx)
    N_VGetNumpyArray(y)[:] = 1.0

    def rhs(t, state, derivative, _):
        N_VGetNumpyArray(derivative)[:] = -10.0 * N_VGetNumpyArray(state)
        return SUN_SUCCESS

    estimator = NegativeRealEstimator(sunctx)
    arkode = LSRKStepCreateSTS(rhs, 0.0, y, sunctx)
    assert LSRKStepSetDomEigEstimator(arkode.get(), estimator) == ARK_SUCCESS
    assert ARKodeSStolerances(arkode.get(), 1e-5, 1e-10) == ARK_SUCCESS
    status, tret = ARKodeEvolve(arkode.get(), 0.1, y, ARK_NORMAL)

    assert status == ARK_SUCCESS
    assert tret == pytest.approx(0.1)
    assert N_VGetNumpyArray(y)[0] == pytest.approx(0.367879, rel=3e-4)
    assert estimator.estimate_calls > 0
