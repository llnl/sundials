#!/usr/bin/env python3
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
# Integrate y' = i * omega * y, y(0) = 1 + 0.5i, with CVODE.
# Build sundials4py with SUNDIALS_SCALAR_TYPE=COMPLEX before running.
# -----------------------------------------------------------------

import cmath

import numpy as np
import pytest

import sundials4py.core as sun
import sundials4py.cvodes as cv


def uses_complex_scalars():
    scalar_type = getattr(sun, "sunscalartype", sun.sunrealtype)
    return np.issubdtype(np.dtype(scalar_type), np.complexfloating)


def integrate():
    omega = 2.0
    initial_value = 1.0 + 0.5j
    initial_time = 0.0
    output_times = (0.25, 0.5, 0.75, 1.0)

    status, context = sun.SUNContext_Create(sun.SUN_COMM_NULL)
    assert status == sun.SUN_SUCCESS

    y = sun.N_VNew_Serial(1, context)
    if y is None:
        raise RuntimeError("N_VNew_Serial failed")
    sun.N_VGetNumpyArray(y)[0] = initial_value

    def rhs(_t, yvec, ydotvec, _user_data):
        sun.N_VGetNumpyArray(ydotvec)[0] = 1j * omega * sun.N_VGetNumpyArray(yvec)[0]
        return 0

    def jac(_t, _yvec, _fyvec, matrix, _user_data, _tmp1, _tmp2, _tmp3):
        sun.SUNDenseMatrix_Data(matrix)[0, 0] = 1j * omega
        return 0

    solver = cv.CVodeCreate(cv.CV_BDF, context)
    if solver is None:
        raise RuntimeError("CVodeCreate failed")
    assert cv.CVodeInit(solver.get(), rhs, initial_time, y) == cv.CV_SUCCESS
    assert cv.CVodeSStolerances(solver.get(), 1e-10, 1e-12) == cv.CV_SUCCESS

    matrix = sun.SUNDenseMatrix(1, 1, context)
    if matrix is None:
        raise RuntimeError("SUNDenseMatrix failed")
    linear_solver = sun.SUNLinSol_Dense(y, matrix, context)
    if linear_solver is None:
        raise RuntimeError("SUNLinSol_Dense failed")
    assert cv.CVodeSetLinearSolver(solver.get(), linear_solver, matrix) == cv.CV_SUCCESS
    assert cv.CVodeSetJacFn(solver.get(), jac) == cv.CV_SUCCESS

    results = []
    for tout in output_times:
        status, t = cv.CVode(solver.get(), tout, y, cv.CV_NORMAL)
        if status != cv.CV_SUCCESS:
            raise RuntimeError(f"CVode failed with status {status}")
        numerical = sun.N_VGetNumpyArray(y)[0]
        exact = initial_value * cmath.exp(1j * omega * t)
        results.append((t, numerical, exact))

    return results


def main():
    if not uses_complex_scalars():
        raise RuntimeError("This example requires SUNDIALS_SCALAR_TYPE=COMPLEX")

    print(f"scalar dtype: {np.dtype(getattr(sun, 'sunscalartype', sun.sunrealtype))}")
    print("    t                 numerical y                   exact y       abs error")
    for t, numerical, exact in integrate():
        print(f"{t:5.2f}  {numerical:>28}  {exact:>28}  {abs(numerical - exact):.2e}")


def test_cvs_complex_oscillator():
    if not uses_complex_scalars():
        pytest.skip("requires SUNDIALS_SCALAR_TYPE=COMPLEX")

    results = integrate()
    np.testing.assert_allclose(results[-1][1], results[-1][2], rtol=1e-8, atol=1e-10)


if __name__ == "__main__":
    main()
