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
from fixtures import *
from sundials4py.core import *


def test_create_dense(sunctx, nvec):
    A = SUNDenseMatrix(N_VGetLength(nvec), N_VGetLength(nvec), sunctx)
    LS = SUNLinSol_Dense(nvec, A, sunctx)
    assert LS is not None


def test_create_band(sunctx, nvec):
    A = SUNBandMatrix(N_VGetLength(nvec), 1, 1, sunctx)
    LS = SUNLinSol_Band(nvec, A, sunctx)
    assert LS is not None


def test_create_spgmr(sunctx, nvec):
    LS = SUNLinSol_SPGMR(nvec, SUN_PREC_NONE, 0, sunctx)
    assert LS is not None


def test_create_pcg(sunctx, nvec):
    LS = SUNLinSol_PCG(nvec, SUN_PREC_NONE, 0, sunctx)
    assert LS is not None


def test_create_spbcgs(sunctx, nvec):
    LS = SUNLinSol_SPBCGS(nvec, SUN_PREC_NONE, 0, sunctx)
    assert LS is not None


def test_create_sptfqmr(sunctx, nvec):
    LS = SUNLinSol_SPTFQMR(nvec, SUN_PREC_NONE, 0, sunctx)
    assert LS is not None


def test_get_type_and_id(sunctx, nvec):
    A = SUNDenseMatrix(N_VGetLength(nvec), N_VGetLength(nvec), sunctx)
    LS = SUNLinSol_Dense(nvec, A, sunctx)
    typ = SUNLinSolGetType(LS)
    id = SUNLinSolGetID(LS)
    assert typ == SUNLINEARSOLVER_DIRECT
    assert id == SUNLINEARSOLVER_DENSE


def test_initialize_setup(sunctx, nvec):
    A = SUNDenseMatrix(N_VGetLength(nvec), N_VGetLength(nvec), sunctx)
    LS = SUNLinSol_Dense(nvec, A, sunctx)
    ret_init = SUNLinSolInitialize(LS)
    ret_setup = SUNLinSolSetup(LS, A)
    assert isinstance(ret_init, int)
    assert isinstance(ret_setup, int)


def test_num_iters_resnorm_lastflag(sunctx, nvec):
    LS = SUNLinSol_SPGMR(nvec, 0, 0, sunctx)
    niters = SUNLinSolNumIters(LS)
    resnorm = SUNLinSolResNorm(LS)
    lastflag = SUNLinSolLastFlag(LS)
    assert isinstance(niters, int)
    assert isinstance(resnorm, float)
    assert isinstance(lastflag, int)


def test_sunlinsol_set_atimes(sunctx):
    x = N_VNew_Serial(1, sunctx)
    y = N_VNew_Serial(1, sunctx)

    LS = SUNLinSol_SPGMR(x, 0, 0, sunctx)
    assert LS is not None

    # Define a dummy ATimes function
    called = {"flag": False}

    def atimes(LS, x, y):
        called["flag"] = True
        return 0

    # Set the ATimes function
    ret = SUNLinSolSetATimes(LS, atimes)
    assert ret == SUN_SUCCESS

    SUNLinSolInitialize(LS)
    SUNLinSolSetup(LS, None)
    SUNLinSolSolve(LS, None, x, y, 1e-2)
    assert called["flag"]


def test_sunlinsol_set_preconditioner(sunctx):
    x = N_VNew_Serial(1, sunctx)
    y = N_VNew_Serial(1, sunctx)

    LS = SUNLinSol_SPGMR(x, SUN_PREC_LEFT, 0, sunctx)
    assert LS is not None

    def atimes(_, x, y):
        return 0

    # Define a dummy preconditioner functions
    called = {"psetup": False, "psolve": False}

    def psetup(_):
        called["psetup"] = True
        return 0

    def psolve(_, r, z, tol, lr):
        called["psolve"] = True
        return 0

    # Set the ATimes function
    ret = SUNLinSolSetATimes(LS, atimes)
    assert ret == SUN_SUCCESS

    ret = SUNLinSolSetPreconditioner(LS, psetup, psolve)
    assert ret == SUN_SUCCESS

    SUNLinSolInitialize(LS)
    SUNLinSolSetup(LS, None)
    SUNLinSolSolve(LS, None, x, y, 1e-2)
    assert called["psetup"]
    assert called["psolve"]
