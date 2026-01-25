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
# Unit/smoke tests for SUNDomEigEstimator module
# -----------------------------------------------------------------

import pytest
from fixtures import *
from sundials4py.core import *


def make_estimator(estimator_type, sunctx):
    if estimator_type == "power":
        nvec = N_VNew_Serial(5, sunctx)
        N_VConst(1.0, nvec)
        e = SUNDomEigEstimator_Power(nvec, 10, 1.0, sunctx)

        def atimes(_, v, z):
            # dummy atimes for smoke testing
            return 0

        SUNDomEigEstimator_SetATimes(e, atimes)

        return e, nvec
    else:
        raise ValueError("Unknown estimator type")


@pytest.mark.parametrize("estimator_type", ["power"])
def test_create_estimator(estimator_type, sunctx):
    est, nvec = make_estimator(estimator_type, sunctx)
    assert est is not None


@pytest.mark.parametrize("estimator_type", ["power"])
def test_set_max_iters(estimator_type, sunctx):
    est, nvec = make_estimator(estimator_type, sunctx)
    status = SUNDomEigEstimator_SetMaxIters(est, 10)
    assert status == SUN_SUCCESS


@pytest.mark.parametrize("estimator_type", ["power"])
def test_set_num_preprocess_iters(estimator_type, sunctx):
    est, nvec = make_estimator(estimator_type, sunctx)
    status = SUNDomEigEstimator_SetNumPreprocessIters(est, 2)
    assert status == SUN_SUCCESS


@pytest.mark.parametrize("estimator_type", ["power"])
def test_set_reltol(estimator_type, sunctx):
    est, nvec = make_estimator(estimator_type, sunctx)
    status = SUNDomEigEstimator_SetRelTol(est, 1e-6)
    assert status == SUN_SUCCESS


@pytest.mark.parametrize("estimator_type", ["power"])
def test_set_initial_guess(estimator_type, sunctx):
    est, nvec = make_estimator(estimator_type, sunctx)
    status = SUNDomEigEstimator_SetInitialGuess(est, nvec)
    assert status == SUN_SUCCESS


@pytest.mark.parametrize("estimator_type", ["power"])
def test_initialize(estimator_type, sunctx):
    est, nvec = make_estimator(estimator_type, sunctx)
    status = SUNDomEigEstimator_Initialize(est)
    assert status == SUN_SUCCESS


@pytest.mark.parametrize("estimator_type", ["power"])
def test_estimate(estimator_type, sunctx):
    est, nvec = make_estimator(estimator_type, sunctx)
    err, lambdaR, lambdaI = SUNDomEigEstimator_Estimate(est)
    assert err == SUN_SUCCESS
    assert isinstance(err, int)
    assert isinstance(lambdaR, float)
    assert isinstance(lambdaI, float)


@pytest.mark.parametrize("estimator_type", ["power"])
def test_get_res(estimator_type, sunctx):
    est, nvec = make_estimator(estimator_type, sunctx)
    err, res = SUNDomEigEstimator_GetRes(est)
    assert err == SUN_SUCCESS
    assert isinstance(err, int)
    assert isinstance(res, float)


@pytest.mark.parametrize("estimator_type", ["power"])
def test_get_num_iters(estimator_type, sunctx):
    est, nvec = make_estimator(estimator_type, sunctx)
    err, num_iters = SUNDomEigEstimator_GetNumIters(est)
    assert err == SUN_SUCCESS
    assert isinstance(err, int)
    assert isinstance(num_iters, int)


@pytest.mark.parametrize("estimator_type", ["power"])
def test_get_num_atimes_calls(estimator_type, sunctx):
    est, nvec = make_estimator(estimator_type, sunctx)
    err, num_atimes = SUNDomEigEstimator_GetNumATimesCalls(est)
    assert err == SUN_SUCCESS
    assert isinstance(err, int)
    assert isinstance(num_atimes, int)
