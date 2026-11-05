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
# Integrate y' = i * omega * y, y(0) = 1 + 0.5i, with ERKStep.
# Build sundials4py with SUNDIALS_SCALAR_TYPE=COMPLEX before running.
# -----------------------------------------------------------------

import cmath
import numpy as np

import sundials4py.arkode as ark
import sundials4py.core as sun


def main():
    if not np.issubdtype(np.dtype(sun.sunscalartype), np.complexfloating):
        raise RuntimeError("This example requires SUNDIALS_SCALAR_TYPE=COMPLEX")

    omega = 2.0
    initial_value = 1.0 + 0.5j

    def rhs(_t, y, ydot, _user_data):
        sun.N_VScale(1j * omega, y, ydot)
        return 0

    status, context = sun.SUNContext_Create(sun.SUN_COMM_NULL)
    if status != sun.SUN_SUCCESS:
        raise RuntimeError("SUNContext_Create failed")

    initial_data = np.array([initial_value], dtype=sun.sunscalartype)
    y = sun.N_VMake_Serial(initial_data.size, initial_data, context)
    if y is None:
        raise RuntimeError("N_VMake_Serial failed")

    solver = ark.ERKStepCreate(rhs, 0.0, y, context)
    if solver is None:
        raise RuntimeError("ERKStepCreate failed")
    if ark.ARKodeSStolerances(solver.get(), 1e-10, 1e-12) != ark.ARK_SUCCESS:
        raise RuntimeError("ARKodeSStolerances failed")

    print(f"scalar dtype: {np.dtype(sun.sunscalartype)}")
    print("    t                 numerical y                   exact y       abs error")
    for tout in (0.25, 0.5, 0.75, 1.0):
        status, t = ark.ARKodeEvolve(solver.get(), tout, y, ark.ARK_NORMAL)
        if status != ark.ARK_SUCCESS:
            raise RuntimeError(f"ARKodeEvolve failed with status {status}")
        numerical = sun.N_VGetNumpyArray(y)[0]
        exact = initial_value * cmath.exp(1j * omega * t)
        print(f"{t:5.2f}  {numerical:>28}  {exact:>28}  {abs(numerical - exact):.2e}")


if __name__ == "__main__":
    main()
