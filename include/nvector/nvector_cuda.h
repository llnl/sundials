/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles and Cody J. Balos @ LLNL
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
 * This is the header file for the CUDA implementation of the
 * NVECTOR module.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_CUDA_H
#define _NVECTOR_CUDA_H

#include <cuda_runtime.h>
#include <stdio.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_cuda_policies.hpp>
#include <sundials/sundials_nvector.h>
#include <sunmemory/sunmemory_cuda.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * CUDA implementation of N_Vector
 * -----------------------------------------------------------------
 */

struct _N_VectorContent_Cuda
{
  sunindextype length;
  sunbooleantype own_helper;
  SUNMemory host_data;
  SUNMemory device_data;
  SUNCudaExecPolicy* stream_exec_policy;
  SUNCudaExecPolicy* reduce_exec_policy;
  SUNMemoryHelper mem_helper;
  void* priv; /* 'private' data */
};

typedef struct _N_VectorContent_Cuda* N_VectorContent_Cuda;

/*
 * -----------------------------------------------------------------
 * NVECTOR_CUDA implementation specific functions
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_Cuda(SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNew_Cuda(sunindextype length, SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNewManaged_Cuda(sunindextype length,
                                            SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNewWithMemHelp_Cuda(sunindextype length,
                                                sunbooleantype use_managed_mem,
                                                SUNMemoryHelper helper,
                                                SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VMake_Cuda(sunindextype length,
                                      sunrealtype* h_vdata_1d,
                                      sunrealtype* d_vdata_1d, SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VMakeManaged_Cuda(sunindextype length,
                                             sunrealtype* vdata_1d,
                                             SUNContext sunctx);
SUNDIALS_EXPORT void N_VSetHostArrayPointer_Cuda(sunrealtype* h_vdata_1d,
                                                 N_Vector v);
SUNDIALS_EXPORT void N_VSetDeviceArrayPointer_Cuda(sunrealtype* d_vdata_1d,
                                                   N_Vector v);
SUNDIALS_EXPORT sunbooleantype N_VIsManagedMemory_Cuda(N_Vector x);
SUNDIALS_EXPORT
SUNErrCode N_VSetKernelExecPolicy_Cuda(N_Vector x,
                                       SUNCudaExecPolicy* stream_exec_policy,
                                       SUNCudaExecPolicy* reduce_exec_policy);
SUNDIALS_EXPORT void N_VCopyToDevice_Cuda(N_Vector v);
SUNDIALS_EXPORT void N_VCopyFromDevice_Cuda(N_Vector v);

static inline sunindextype N_VGetLength_Cuda(N_Vector x)
{
  N_VectorContent_Cuda content = (N_VectorContent_Cuda)x->content;
  return content->length;
}

static inline sunrealtype* N_VGetHostArrayPointer_Cuda(N_Vector x)
{
  N_VectorContent_Cuda content = (N_VectorContent_Cuda)x->content;
  return (content->host_data == NULL ? NULL
                                     : (sunrealtype*)content->host_data->ptr);
}

static inline sunrealtype* N_VGetDeviceArrayPointer_Cuda(N_Vector x)
{
  N_VectorContent_Cuda content = (N_VectorContent_Cuda)x->content;
  return (content->device_data == NULL ? NULL
                                       : (sunrealtype*)content->device_data->ptr);
}

/*
 * -----------------------------------------------------------------
 * NVECTOR API functions
 * -----------------------------------------------------------------
 */

static inline N_Vector_ID N_VGetVectorID_Cuda(N_Vector /*v*/)
{
  return SUNDIALS_NVEC_CUDA;
}

SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Cuda(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_Cuda(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_Cuda(N_Vector v);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
void N_VSpace_Cuda(N_Vector v, sunindextype* lrw, sunindextype* liw);

/* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum_Cuda(sunrealtype a, N_Vector x, sunrealtype b,
                                       N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_Cuda(sunrealtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_Cuda(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_Cuda(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_Cuda(sunrealtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_Cuda(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_Cuda(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_Cuda(N_Vector x, sunrealtype b, N_Vector z);
SUNDIALS_EXPORT sunrealtype N_VDotProd_Cuda(N_Vector x, N_Vector y);
SUNDIALS_EXPORT sunrealtype N_VMaxNorm_Cuda(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWrmsNorm_Cuda(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VWrmsNormMask_Cuda(N_Vector x, N_Vector w,
                                                 N_Vector id);
SUNDIALS_EXPORT sunrealtype N_VMin_Cuda(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWL2Norm_Cuda(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VL1Norm_Cuda(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_Cuda(sunrealtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VInvTest_Cuda(N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VConstrMask_Cuda(N_Vector c, N_Vector x,
                                                  N_Vector m);
SUNDIALS_EXPORT sunrealtype N_VMinQuotient_Cuda(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT SUNErrCode N_VLinearCombination_Cuda(int nvec, sunrealtype* c,
                                                     N_Vector* X, N_Vector Z);
SUNDIALS_EXPORT SUNErrCode N_VScaleAddMulti_Cuda(int nvec, sunrealtype* c,
                                                 N_Vector X, N_Vector* Y,
                                                 N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VDotProdMulti_Cuda(int nvec, N_Vector x, N_Vector* Y,
                                                sunrealtype* dotprods);

/* vector array operations */
SUNDIALS_EXPORT SUNErrCode N_VLinearSumVectorArray_Cuda(
  int nvec, sunrealtype a, N_Vector* X, sunrealtype b, N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VScaleVectorArray_Cuda(int nvec, sunrealtype* c,
                                                    N_Vector* X, N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VConstVectorArray_Cuda(int nvec, sunrealtype c,
                                                    N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VScaleAddMultiVectorArray_Cuda(
  int nvec, int nsum, sunrealtype* a, N_Vector* X, N_Vector** Y, N_Vector** Z);
SUNDIALS_EXPORT SUNErrCode N_VLinearCombinationVectorArray_Cuda(
  int nvec, int nsum, sunrealtype* c, N_Vector** X, N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VWrmsNormVectorArray_Cuda(int nvec, N_Vector* X,
                                                       N_Vector* W,
                                                       sunrealtype* nrm);
SUNDIALS_EXPORT SUNErrCode N_VWrmsNormMaskVectorArray_Cuda(int nvec, N_Vector* X,
                                                           N_Vector* W,
                                                           N_Vector id,
                                                           sunrealtype* nrm);

/* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT sunrealtype N_VWSqrSumLocal_Cuda(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VWSqrSumMaskLocal_Cuda(N_Vector x, N_Vector w,
                                                     N_Vector id);

/* OPTIONAL XBraid interface operations */
SUNDIALS_EXPORT SUNErrCode N_VBufSize_Cuda(N_Vector x, sunindextype* size);
SUNDIALS_EXPORT SUNErrCode N_VBufPack_Cuda(N_Vector x, void* buf);
SUNDIALS_EXPORT SUNErrCode N_VBufUnpack_Cuda(N_Vector x, void* buf);

/* OPTIONAL operations for debugging */
SUNDIALS_EXPORT void N_VPrint_Cuda(N_Vector v);
SUNDIALS_EXPORT void N_VPrintFile_Cuda(N_Vector v, FILE* outfile);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNErrCode N_VEnableFusedOps_Cuda(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableLinearCombination_Cuda(N_Vector v,
                                                           sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleAddMulti_Cuda(N_Vector v,
                                                       sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableDotProdMulti_Cuda(N_Vector v,
                                                      sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableLinearSumVectorArray_Cuda(N_Vector v,
                                                              sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleVectorArray_Cuda(N_Vector v,
                                                          sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableConstVectorArray_Cuda(N_Vector v,
                                                          sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableWrmsNormVectorArray_Cuda(N_Vector v,
                                                             sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormMaskVectorArray_Cuda(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMultiVectorArray_Cuda(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombinationVectorArray_Cuda(N_Vector v,
                                                      sunbooleantype tf);

#ifdef __cplusplus
}
#endif

#endif
