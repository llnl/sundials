/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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

#ifndef SUNDIALS_NVECTOR_NVECTOR_SYCL_DEPRECATED_H
#define SUNDIALS_NVECTOR_NVECTOR_SYCL_DEPRECATED_H

#include <nvector/nvector_sycl.h>

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VSetArrayPointer instead; will be removed in version 8.1.0")
void N_VSetHostArrayPointer_Sycl(sunrealtype* h_vdata_1d, N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VSetDeviceArrayPointer instead; will be removed in version 8.1.0")
void N_VSetDeviceArrayPointer_Sycl(sunrealtype* d_vdata_1d, N_Vector v);

SUNDIALS_DEPRECATED_MSG(
  "use N_VGetLength instead; will be removed in version 8.1.0")

static inline sunindextype N_VGetLength_Sycl(N_Vector x)
{
  N_VectorContent_Sycl content = (N_VectorContent_Sycl)x->content;

  return content->length;
}
SUNDIALS_DEPRECATED_MSG(
  "use N_VGetArrayPointer instead; will be removed in version 8.1.0")

static inline sunrealtype* N_VGetHostArrayPointer_Sycl(N_Vector x)
{
  N_VectorContent_Sycl content = (N_VectorContent_Sycl)x->content;

  return (content->host_data == NULL ? NULL
                                     : (sunrealtype*)content->host_data->ptr);
}
SUNDIALS_DEPRECATED_MSG(
  "use N_VGetDeviceArrayPointer instead; will be removed in version 8.1.0")

static inline sunrealtype* N_VGetDeviceArrayPointer_Sycl(N_Vector x)
{
  N_VectorContent_Sycl content = (N_VectorContent_Sycl)x->content;

  return (content->device_data == NULL ? NULL
                                       : (sunrealtype*)content->device_data->ptr);
}
SUNDIALS_DEPRECATED_MSG(
  "use N_VGetVectorID instead; will be removed in version 8.1.0")

static inline N_Vector_ID N_VGetVectorID_Sycl(N_Vector)
{
  return SUNDIALS_NVEC_SYCL;
}
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VCloneEmpty instead; will be removed in version 8.1.0")
N_Vector N_VCloneEmpty_Sycl(N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VClone instead; will be removed in version 8.1.0")
N_Vector N_VClone_Sycl(N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDestroy instead; will be removed in version 8.1.0")
void N_VDestroy_Sycl(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VLinearSum instead; will be removed in version 8.1.0")
void N_VLinearSum_Sycl(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                       N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VConst instead; will be removed in version 8.1.0")
void N_VConst_Sycl(sunrealtype c, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VProd instead; will be removed in version 8.1.0")
void N_VProd_Sycl(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDiv instead; will be removed in version 8.1.0")
void N_VDiv_Sycl(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VScale instead; will be removed in version 8.1.0")
void N_VScale_Sycl(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VAbs instead; will be removed in version 8.1.0")
void N_VAbs_Sycl(N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VInv instead; will be removed in version 8.1.0")
void N_VInv_Sycl(N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VAddConst instead; will be removed in version 8.1.0")
void N_VAddConst_Sycl(N_Vector x, sunrealtype b, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDotProd instead; will be removed in version 8.1.0")
sunrealtype N_VDotProd_Sycl(N_Vector x, N_Vector y);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMaxNorm instead; will be removed in version 8.1.0")
sunrealtype N_VMaxNorm_Sycl(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWrmsNorm instead; will be removed in version 8.1.0")
sunrealtype N_VWrmsNorm_Sycl(N_Vector x, N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWrmsNormMask instead; will be removed in version 8.1.0")
sunrealtype N_VWrmsNormMask_Sycl(N_Vector x, N_Vector w, N_Vector id);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMin instead; will be removed in version 8.1.0")
sunrealtype N_VMin_Sycl(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWL2Norm instead; will be removed in version 8.1.0")
sunrealtype N_VWL2Norm_Sycl(N_Vector x, N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VL1Norm instead; will be removed in version 8.1.0")
sunrealtype N_VL1Norm_Sycl(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VCompare instead; will be removed in version 8.1.0")
void N_VCompare_Sycl(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VInvTest instead; will be removed in version 8.1.0")
sunbooleantype N_VInvTest_Sycl(N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VConstrMask instead; will be removed in version 8.1.0")
sunbooleantype N_VConstrMask_Sycl(N_Vector c, N_Vector x, N_Vector m);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMinQuotient instead; will be removed in version 8.1.0")
sunrealtype N_VMinQuotient_Sycl(N_Vector num, N_Vector denom);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VLinearCombination instead; will be removed in version 8.1.0")
SUNErrCode N_VLinearCombination_Sycl(int nvec, sunrealtype* c, N_Vector* X,
                                     N_Vector Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VScaleAddMulti instead; will be removed in version 8.1.0")
SUNErrCode N_VScaleAddMulti_Sycl(int nvec, sunrealtype* c, N_Vector X,
                                 N_Vector* Y, N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VLinearSumVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VLinearSumVectorArray_Sycl(int nvec, sunrealtype a, N_Vector* X,
                                        sunrealtype b, N_Vector* Y, N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VScaleVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VScaleVectorArray_Sycl(int nvec, sunrealtype* c, N_Vector* X,
                                    N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VConstVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VConstVectorArray_Sycl(int nvec, sunrealtype c, N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VScaleAddMultiVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VScaleAddMultiVectorArray_Sycl(int nvec, int nsum, sunrealtype* a,
                                            N_Vector* X, N_Vector** Y,
                                            N_Vector** Z);

SUNDIALS_DEPRECATED_EXPORT_MSG("use N_VLinearCombinationVectorArray instead; "
                               "will be removed in version 8.1.0")
SUNErrCode N_VLinearCombinationVectorArray_Sycl(int nvec, int nsum,
                                                sunrealtype* c, N_Vector** X,
                                                N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWSqrSumLocal instead; will be removed in version 8.1.0")
sunrealtype N_VWSqrSumLocal_Sycl(N_Vector x, N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWSqrSumMaskLocal instead; will be removed in version 8.1.0")
sunrealtype N_VWSqrSumMaskLocal_Sycl(N_Vector x, N_Vector w, N_Vector id);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VBufSize instead; will be removed in version 8.1.0")
SUNErrCode N_VBufSize_Sycl(N_Vector x, sunindextype* size);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VBufPack instead; will be removed in version 8.1.0")
SUNErrCode N_VBufPack_Sycl(N_Vector x, void* buf);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VBufUnpack instead; will be removed in version 8.1.0")
SUNErrCode N_VBufUnpack_Sycl(N_Vector x, void* buf);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VPrint instead; will be removed in version 8.1.0")
void N_VPrint_Sycl(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VPrintFile instead; will be removed in version 8.1.0")
void N_VPrintFile_Sycl(N_Vector v, FILE* outfile);

#ifdef __cplusplus
}
#endif

#endif
