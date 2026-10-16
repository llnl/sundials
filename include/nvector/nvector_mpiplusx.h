/* -----------------------------------------------------------------
 * Programmer(s): Cody Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the MPI+X implementation of the
 * NVECTOR module. The MPIPlusX NVECTOR is really just an extension
 * of the ManyVector.
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_NVECTOR_MPIPLUSX_H
#define SUNDIALS_NVECTOR_MPIPLUSX_H

#include <mpi.h>
#include <nvector/nvector_mpimanyvector.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef N_VectorContent_MPIManyVector N_VectorContent_MPIPlusX;

SUNDIALS_EXPORT
N_Vector N_VMake_MPIPlusX(MPI_Comm comm, N_Vector X, SUNContext sunctx);

SUNDIALS_EXPORT
N_Vector N_VGetLocalVector_MPIPlusX(N_Vector v);

SUNDIALS_EXPORT
SUNErrCode N_VEnableFusedOps_MPIPlusX(N_Vector v, sunbooleantype tf);

#ifdef __cplusplus
}
#endif

#endif
