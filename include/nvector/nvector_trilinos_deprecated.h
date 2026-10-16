/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
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

#ifndef SUNDIALS_NVECTOR_NVECTOR_TRILINOS_DEPRECATED_H
#define SUNDIALS_NVECTOR_NVECTOR_TRILINOS_DEPRECATED_H

#include <nvector/nvector_trilinos.h>

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VGetVectorID instead; will be removed in version 8.1.0")
N_Vector_ID N_VGetVectorID_Trilinos(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VCloneEmpty instead; will be removed in version 8.1.0")
N_Vector N_VCloneEmpty_Trilinos(N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VClone instead; will be removed in version 8.1.0")
N_Vector N_VClone_Trilinos(N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDestroy instead; will be removed in version 8.1.0")
void N_VDestroy_Trilinos(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VGetCommunicator instead; will be removed in version 8.1.0")
SUNComm N_VGetCommunicator_Trilinos(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VGetLength instead; will be removed in version 8.1.0")
sunindextype N_VGetLength_Trilinos(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VLinearSum instead; will be removed in version 8.1.0")
void N_VLinearSum_Trilinos(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                           N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VConst instead; will be removed in version 8.1.0")
void N_VConst_Trilinos(sunrealtype c, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VProd instead; will be removed in version 8.1.0")
void N_VProd_Trilinos(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDiv instead; will be removed in version 8.1.0")
void N_VDiv_Trilinos(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VScale instead; will be removed in version 8.1.0")
void N_VScale_Trilinos(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VAbs instead; will be removed in version 8.1.0")
void N_VAbs_Trilinos(N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VInv instead; will be removed in version 8.1.0")
void N_VInv_Trilinos(N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VAddConst instead; will be removed in version 8.1.0")
void N_VAddConst_Trilinos(N_Vector x, sunrealtype b, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDotProd instead; will be removed in version 8.1.0")
sunrealtype N_VDotProd_Trilinos(N_Vector x, N_Vector y);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMaxNorm instead; will be removed in version 8.1.0")
sunrealtype N_VMaxNorm_Trilinos(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWrmsNorm instead; will be removed in version 8.1.0")
sunrealtype N_VWrmsNorm_Trilinos(N_Vector x, N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWrmsNormMask instead; will be removed in version 8.1.0")
sunrealtype N_VWrmsNormMask_Trilinos(N_Vector x, N_Vector w, N_Vector id);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMin instead; will be removed in version 8.1.0")
sunrealtype N_VMin_Trilinos(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWL2Norm instead; will be removed in version 8.1.0")
sunrealtype N_VWL2Norm_Trilinos(N_Vector x, N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VL1Norm instead; will be removed in version 8.1.0")
sunrealtype N_VL1Norm_Trilinos(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VCompare instead; will be removed in version 8.1.0")
void N_VCompare_Trilinos(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VInvTest instead; will be removed in version 8.1.0")
sunbooleantype N_VInvTest_Trilinos(N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VConstrMask instead; will be removed in version 8.1.0")
sunbooleantype N_VConstrMask_Trilinos(N_Vector c, N_Vector x, N_Vector m);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMinQuotient instead; will be removed in version 8.1.0")
sunrealtype N_VMinQuotient_Trilinos(N_Vector num, N_Vector denom);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDotProdLocal instead; will be removed in version 8.1.0")
sunrealtype N_VDotProdLocal_Trilinos(N_Vector x, N_Vector y);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMaxNormLocal instead; will be removed in version 8.1.0")
sunrealtype N_VMaxNormLocal_Trilinos(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMinLocal instead; will be removed in version 8.1.0")
sunrealtype N_VMinLocal_Trilinos(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VL1NormLocal instead; will be removed in version 8.1.0")
sunrealtype N_VL1NormLocal_Trilinos(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWSqrSumLocal instead; will be removed in version 8.1.0")
sunrealtype N_VWSqrSumLocal_Trilinos(N_Vector x, N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWSqrSumMaskLocal instead; will be removed in version 8.1.0")
sunrealtype N_VWSqrSumMaskLocal_Trilinos(N_Vector x, N_Vector w, N_Vector id);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VInvTestLocal instead; will be removed in version 8.1.0")
sunbooleantype N_VInvTestLocal_Trilinos(N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VConstrMaskLocal instead; will be removed in version 8.1.0")
sunbooleantype N_VConstrMaskLocal_Trilinos(N_Vector c, N_Vector x, N_Vector m);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMinQuotientLocal instead; will be removed in version 8.1.0")
sunrealtype N_VMinQuotientLocal_Trilinos(N_Vector num, N_Vector denom);

#ifdef __cplusplus
}
#endif

#endif
