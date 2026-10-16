/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
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

#ifndef SUNDIALS_NVECTOR_NVECTOR_MPIMANYVECTOR_DEPRECATED_H
#define SUNDIALS_NVECTOR_NVECTOR_MPIMANYVECTOR_DEPRECATED_H

#include <nvector/nvector_mpimanyvector.h>

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VGetVectorID instead; will be removed in version 8.1.0")
N_Vector_ID N_VGetVectorID_MPIManyVector(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VPrint instead; will be removed in version 8.1.0")
void N_VPrint_MPIManyVector(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VPrintFile instead; will be removed in version 8.1.0")
void N_VPrintFile_MPIManyVector(N_Vector v, FILE* outfile);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VCloneEmpty instead; will be removed in version 8.1.0")
N_Vector N_VCloneEmpty_MPIManyVector(N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VClone instead; will be removed in version 8.1.0")
N_Vector N_VClone_MPIManyVector(N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDestroy instead; will be removed in version 8.1.0")
void N_VDestroy_MPIManyVector(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VGetCommunicator instead; will be removed in version 8.1.0")
MPI_Comm N_VGetCommunicator_MPIManyVector(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VGetLength instead; will be removed in version 8.1.0")
sunindextype N_VGetLength_MPIManyVector(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VLinearSum instead; will be removed in version 8.1.0")
void N_VLinearSum_MPIManyVector(sunrealtype a, N_Vector x, sunrealtype b,
                                N_Vector y, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VConst instead; will be removed in version 8.1.0")
void N_VConst_MPIManyVector(sunrealtype c, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VProd instead; will be removed in version 8.1.0")
void N_VProd_MPIManyVector(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDiv instead; will be removed in version 8.1.0")
void N_VDiv_MPIManyVector(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VScale instead; will be removed in version 8.1.0")
void N_VScale_MPIManyVector(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VAbs instead; will be removed in version 8.1.0")
void N_VAbs_MPIManyVector(N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VInv instead; will be removed in version 8.1.0")
void N_VInv_MPIManyVector(N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VAddConst instead; will be removed in version 8.1.0")
void N_VAddConst_MPIManyVector(N_Vector x, sunrealtype b, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDotProd instead; will be removed in version 8.1.0")
sunrealtype N_VDotProd_MPIManyVector(N_Vector x, N_Vector y);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMaxNorm instead; will be removed in version 8.1.0")
sunrealtype N_VMaxNorm_MPIManyVector(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWrmsNorm instead; will be removed in version 8.1.0")
sunrealtype N_VWrmsNorm_MPIManyVector(N_Vector x, N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWrmsNormMask instead; will be removed in version 8.1.0")
sunrealtype N_VWrmsNormMask_MPIManyVector(N_Vector x, N_Vector w, N_Vector id);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMin instead; will be removed in version 8.1.0")
sunrealtype N_VMin_MPIManyVector(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWL2Norm instead; will be removed in version 8.1.0")
sunrealtype N_VWL2Norm_MPIManyVector(N_Vector x, N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VL1Norm instead; will be removed in version 8.1.0")
sunrealtype N_VL1Norm_MPIManyVector(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VCompare instead; will be removed in version 8.1.0")
void N_VCompare_MPIManyVector(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VInvTest instead; will be removed in version 8.1.0")
sunbooleantype N_VInvTest_MPIManyVector(N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VConstrMask instead; will be removed in version 8.1.0")
sunbooleantype N_VConstrMask_MPIManyVector(N_Vector c, N_Vector x, N_Vector m);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMinQuotient instead; will be removed in version 8.1.0")
sunrealtype N_VMinQuotient_MPIManyVector(N_Vector num, N_Vector denom);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VLinearCombination instead; will be removed in version 8.1.0")
SUNErrCode N_VLinearCombination_MPIManyVector(int nvec, sunrealtype* c,
                                              N_Vector* V, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VScaleAddMulti instead; will be removed in version 8.1.0")
SUNErrCode N_VScaleAddMulti_MPIManyVector(int nvec, sunrealtype* a, N_Vector x,
                                          N_Vector* Y, N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDotProdMulti instead; will be removed in version 8.1.0")
SUNErrCode N_VDotProdMulti_MPIManyVector(int nvec, N_Vector x, N_Vector* Y,
                                         sunrealtype* dotprods);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VLinearSumVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VLinearSumVectorArray_MPIManyVector(int nvec, sunrealtype a,
                                                 N_Vector* X, sunrealtype b,
                                                 N_Vector* Y, N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VScaleVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VScaleVectorArray_MPIManyVector(int nvec, sunrealtype* c,
                                             N_Vector* X, N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VConstVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VConstVectorArray_MPIManyVector(int nvecs, sunrealtype c,
                                             N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDotProdMultiLocal instead; will be removed in version 8.1.0")
SUNErrCode N_VDotProdMultiLocal_MPIManyVector(int nvec, N_Vector x, N_Vector* Y,
                                              sunrealtype* dotprods);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDotProdMultiAllReduce instead; will be removed in version 8.1.0")
SUNErrCode N_VDotProdMultiAllReduce_MPIManyVector(int nvec_total, N_Vector x,
                                                  sunrealtype* sum);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWrmsNormVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VWrmsNormVectorArray_MPIManyVector(int nvecs, N_Vector* X,
                                                N_Vector* W, sunrealtype* nrm);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWrmsNormMaskVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VWrmsNormMaskVectorArray_MPIManyVector(int nvec, N_Vector* X,
                                                    N_Vector* W, N_Vector id,
                                                    sunrealtype* nrm);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDotProdLocal instead; will be removed in version 8.1.0")
sunrealtype N_VDotProdLocal_MPIManyVector(N_Vector x, N_Vector y);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMaxNormLocal instead; will be removed in version 8.1.0")
sunrealtype N_VMaxNormLocal_MPIManyVector(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMinLocal instead; will be removed in version 8.1.0")
sunrealtype N_VMinLocal_MPIManyVector(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VL1NormLocal instead; will be removed in version 8.1.0")
sunrealtype N_VL1NormLocal_MPIManyVector(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWSqrSumLocal instead; will be removed in version 8.1.0")
sunrealtype N_VWSqrSumLocal_MPIManyVector(N_Vector x, N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWSqrSumMaskLocal instead; will be removed in version 8.1.0")
sunrealtype N_VWSqrSumMaskLocal_MPIManyVector(N_Vector x, N_Vector w,
                                              N_Vector id);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VInvTestLocal instead; will be removed in version 8.1.0")
sunbooleantype N_VInvTestLocal_MPIManyVector(N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VConstrMaskLocal instead; will be removed in version 8.1.0")
sunbooleantype N_VConstrMaskLocal_MPIManyVector(N_Vector c, N_Vector x,
                                                N_Vector m);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMinQuotientLocal instead; will be removed in version 8.1.0")
sunrealtype N_VMinQuotientLocal_MPIManyVector(N_Vector num, N_Vector denom);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VBufSize instead; will be removed in version 8.1.0")
SUNErrCode N_VBufSize_MPIManyVector(N_Vector x, sunindextype* size);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VBufPack instead; will be removed in version 8.1.0")
SUNErrCode N_VBufPack_MPIManyVector(N_Vector x, void* buf);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VBufUnpack instead; will be removed in version 8.1.0")
SUNErrCode N_VBufUnpack_MPIManyVector(N_Vector x, void* buf);

#ifdef __cplusplus
}
#endif

#endif
