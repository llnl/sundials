/* -------------------------------------------------------------------
 * Programmer(s): David J. Gardner and Shelby Lockhart @ LLNL
 * -------------------------------------------------------------------
 * Acknowledgements: This NVECTOR module is based on the NVECTOR
 *                   Serial module by Scott D. Cohen, Alan C.
 *                   Hindmarsh, Radu Serban, and Aaron Collier
 *                   @ LLNL
 * -------------------------------------------------------------------
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

#ifndef SUNDIALS_NVECTOR_NVECTOR_OPENMPDEV_DEPRECATED_H
#define SUNDIALS_NVECTOR_NVECTOR_OPENMPDEV_DEPRECATED_H

#include <nvector/nvector_openmpdev.h>

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VGetLength instead; will be removed in version 8.1.0")
sunindextype N_VGetLength_OpenMPDEV(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VGetArrayPointer instead; will be removed in version 8.1.0")
sunrealtype* N_VGetHostArrayPointer_OpenMPDEV(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VGetDeviceArrayPointer instead; will be removed in version 8.1.0")
sunrealtype* N_VGetDeviceArrayPointer_OpenMPDEV(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VPrint instead; will be removed in version 8.1.0")
void N_VPrint_OpenMPDEV(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VPrintFile instead; will be removed in version 8.1.0")
void N_VPrintFile_OpenMPDEV(N_Vector v, FILE* outfile);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VGetVectorID instead; will be removed in version 8.1.0")
N_Vector_ID N_VGetVectorID_OpenMPDEV(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VCloneEmpty instead; will be removed in version 8.1.0")
N_Vector N_VCloneEmpty_OpenMPDEV(N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VClone instead; will be removed in version 8.1.0")
N_Vector N_VClone_OpenMPDEV(N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDestroy instead; will be removed in version 8.1.0")
void N_VDestroy_OpenMPDEV(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VLinearSum instead; will be removed in version 8.1.0")
void N_VLinearSum_OpenMPDEV(sunrealtype a, N_Vector x, sunrealtype b,
                            N_Vector y, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VConst instead; will be removed in version 8.1.0")
void N_VConst_OpenMPDEV(sunrealtype c, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VProd instead; will be removed in version 8.1.0")
void N_VProd_OpenMPDEV(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDiv instead; will be removed in version 8.1.0")
void N_VDiv_OpenMPDEV(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VScale instead; will be removed in version 8.1.0")
void N_VScale_OpenMPDEV(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VAbs instead; will be removed in version 8.1.0")
void N_VAbs_OpenMPDEV(N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VInv instead; will be removed in version 8.1.0")
void N_VInv_OpenMPDEV(N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VAddConst instead; will be removed in version 8.1.0")
void N_VAddConst_OpenMPDEV(N_Vector x, sunrealtype b, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDotProd instead; will be removed in version 8.1.0")
sunrealtype N_VDotProd_OpenMPDEV(N_Vector x, N_Vector y);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMaxNorm instead; will be removed in version 8.1.0")
sunrealtype N_VMaxNorm_OpenMPDEV(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWrmsNorm instead; will be removed in version 8.1.0")
sunrealtype N_VWrmsNorm_OpenMPDEV(N_Vector x, N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWrmsNormMask instead; will be removed in version 8.1.0")
sunrealtype N_VWrmsNormMask_OpenMPDEV(N_Vector x, N_Vector w, N_Vector id);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMin instead; will be removed in version 8.1.0")
sunrealtype N_VMin_OpenMPDEV(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWL2Norm instead; will be removed in version 8.1.0")
sunrealtype N_VWL2Norm_OpenMPDEV(N_Vector x, N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VL1Norm instead; will be removed in version 8.1.0")
sunrealtype N_VL1Norm_OpenMPDEV(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VCompare instead; will be removed in version 8.1.0")
void N_VCompare_OpenMPDEV(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VInvTest instead; will be removed in version 8.1.0")
sunbooleantype N_VInvTest_OpenMPDEV(N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VConstrMask instead; will be removed in version 8.1.0")
sunbooleantype N_VConstrMask_OpenMPDEV(N_Vector c, N_Vector x, N_Vector m);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMinQuotient instead; will be removed in version 8.1.0")
sunrealtype N_VMinQuotient_OpenMPDEV(N_Vector num, N_Vector denom);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VLinearCombination instead; will be removed in version 8.1.0")
SUNErrCode N_VLinearCombination_OpenMPDEV(int nvec, sunrealtype* c, N_Vector* V,
                                          N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VScaleAddMulti instead; will be removed in version 8.1.0")
SUNErrCode N_VScaleAddMulti_OpenMPDEV(int nvec, sunrealtype* a, N_Vector x,
                                      N_Vector* Y, N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDotProdMulti instead; will be removed in version 8.1.0")
SUNErrCode N_VDotProdMulti_OpenMPDEV(int nvec, N_Vector x, N_Vector* Y,
                                     sunrealtype* dotprods);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VLinearSumVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VLinearSumVectorArray_OpenMPDEV(int nvec, sunrealtype a,
                                             N_Vector* X, sunrealtype b,
                                             N_Vector* Y, N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VScaleVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VScaleVectorArray_OpenMPDEV(int nvec, sunrealtype* c, N_Vector* X,
                                         N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VConstVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VConstVectorArray_OpenMPDEV(int nvecs, sunrealtype c, N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWrmsNormVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VWrmsNormVectorArray_OpenMPDEV(int nvecs, N_Vector* X, N_Vector* W,
                                            sunrealtype* nrm);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWrmsNormMaskVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VWrmsNormMaskVectorArray_OpenMPDEV(int nvecs, N_Vector* X,
                                                N_Vector* W, N_Vector id,
                                                sunrealtype* nrm);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VScaleAddMultiVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VScaleAddMultiVectorArray_OpenMPDEV(int nvec, int nsum,
                                                 sunrealtype* a, N_Vector* X,
                                                 N_Vector** Y, N_Vector** Z);

SUNDIALS_DEPRECATED_EXPORT_MSG("use N_VLinearCombinationVectorArray instead; "
                               "will be removed in version 8.1.0")
SUNErrCode N_VLinearCombinationVectorArray_OpenMPDEV(int nvec, int nsum,
                                                     sunrealtype* c,
                                                     N_Vector** X, N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWSqrSumLocal instead; will be removed in version 8.1.0")
sunrealtype N_VWSqrSumLocal_OpenMPDEV(N_Vector x, N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWSqrSumMaskLocal instead; will be removed in version 8.1.0")
sunrealtype N_VWSqrSumMaskLocal_OpenMPDEV(N_Vector x, N_Vector w, N_Vector id);

#ifdef __cplusplus
}
#endif

#endif
