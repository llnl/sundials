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
 * Temporary header file for deprecated functions.
 * TODO:(SBR) remove in version 8.1.0
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_NVECTOR_NVECTOR_MPIPLUSX_DEPRECATED_H
#define SUNDIALS_NVECTOR_NVECTOR_MPIPLUSX_DEPRECATED_H

#include <nvector/nvector_mpiplusx.h>

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VGetVectorID instead; will be removed in version 8.1.0")
N_Vector_ID N_VGetVectorID_MPIPlusX(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VGetArrayPointer instead; will be removed in version 8.1.0")
sunrealtype* N_VGetArrayPointer_MPIPlusX(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VSetArrayPointer instead; will be removed in version 8.1.0")
void N_VSetArrayPointer_MPIPlusX(sunrealtype* vdata, N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VGetLocalLength instead; will be removed in version 8.1.0")
sunindextype N_VGetLocalLength_MPIPlusX(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VPrint instead; will be removed in version 8.1.0")
void N_VPrint_MPIPlusX(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VPrintFile instead; will be removed in version 8.1.0")
void N_VPrintFile_MPIPlusX(N_Vector x, FILE* outfile);

#ifdef __cplusplus
}
#endif

#endif
