/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
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
 * This is the header file for the SYCL implementation of the
 * NVECTOR module.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_SYCL_H
#define _NVECTOR_SYCL_H

#include <stdio.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sycl_policies.hpp>
#include <sunmemory/sunmemory_sycl.h>
#include <sycl/sycl.hpp>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------------------
 * SYCL implementation of N_Vector
 * ----------------------------------------------------------------- */

struct _N_VectorContent_Sycl
{
  sunindextype length;
  sunbooleantype own_helper;
  SUNMemory host_data;
  SUNMemory device_data;
  SUNSyclExecPolicy* stream_exec_policy;
  SUNSyclExecPolicy* reduce_exec_policy;
  SUNMemoryHelper mem_helper;
  ::sycl::queue* queue;
  void* priv; /* 'private' data */
};

typedef struct _N_VectorContent_Sycl* N_VectorContent_Sycl;

/* -----------------------------------------------------------------
 * NVECTOR_SYCL implementation specific functions
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_Sycl(SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNew_Sycl(sunindextype length, ::sycl::queue* Q,
                                     SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNewManaged_Sycl(sunindextype length,
                                            ::sycl::queue* Q, SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNewWithMemHelp_Sycl(sunindextype length,
                                                sunbooleantype use_managed_mem,
                                                SUNMemoryHelper helper,
                                                ::sycl::queue* Q,
                                                SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VMake_Sycl(sunindextype length,
                                      sunrealtype* h_vdata_1d,
                                      sunrealtype* d_vdata_1d, ::sycl::queue* Q,
                                      SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VMakeManaged_Sycl(sunindextype length,
                                             sunrealtype* vdata_1d,
                                             ::sycl::queue* Q, SUNContext sunctx);

SUNDIALS_EXPORT void N_VSetHostArrayPointer_Sycl(sunrealtype* h_vdata_1d,
                                                 N_Vector v);
SUNDIALS_EXPORT void N_VSetDeviceArrayPointer_Sycl(sunrealtype* d_vdata_1d,
                                                   N_Vector v);
SUNDIALS_EXPORT sunbooleantype N_VIsManagedMemory_Sycl(N_Vector x);
SUNDIALS_EXPORT
SUNErrCode N_VSetKernelExecPolicy_Sycl(N_Vector x,
                                       SUNSyclExecPolicy* stream_exec_policy,
                                       SUNSyclExecPolicy* reduce_exec_policy);
SUNDIALS_EXPORT void N_VCopyToDevice_Sycl(N_Vector v);
SUNDIALS_EXPORT void N_VCopyFromDevice_Sycl(N_Vector v);

static inline sunindextype N_VGetLength_Sycl(N_Vector x)
{
  N_VectorContent_Sycl content = (N_VectorContent_Sycl)x->content;
  return content->length;
}

static inline sunrealtype* N_VGetHostArrayPointer_Sycl(N_Vector x)
{
  N_VectorContent_Sycl content = (N_VectorContent_Sycl)x->content;
  return (content->host_data == NULL ? NULL
                                     : (sunrealtype*)content->host_data->ptr);
}

static inline sunrealtype* N_VGetDeviceArrayPointer_Sycl(N_Vector x)
{
  N_VectorContent_Sycl content = (N_VectorContent_Sycl)x->content;
  return (content->device_data == NULL ? NULL
                                       : (sunrealtype*)content->device_data->ptr);
}

/* -----------------------------------------------------------------
 * NVECTOR API functions
 * ----------------------------------------------------------------- */

static inline N_Vector_ID N_VGetVectorID_Sycl(N_Vector)
{
  return SUNDIALS_NVEC_SYCL;
}

SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Sycl(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_Sycl(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_Sycl(N_Vector v);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
void N_VSpace_Sycl(N_Vector v, sunindextype* lrw, sunindextype* liw);

/* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum_Sycl(sunrealtype a, N_Vector x, sunrealtype b,
                                       N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_Sycl(sunrealtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_Sycl(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_Sycl(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_Sycl(sunrealtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_Sycl(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_Sycl(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_Sycl(N_Vector x, sunrealtype b, N_Vector z);
SUNDIALS_EXPORT sunrealtype N_VDotProd_Sycl(N_Vector x, N_Vector y);
SUNDIALS_EXPORT sunrealtype N_VMaxNorm_Sycl(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWrmsNorm_Sycl(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VWrmsNormMask_Sycl(N_Vector x, N_Vector w,
                                                 N_Vector id);
SUNDIALS_EXPORT sunrealtype N_VMin_Sycl(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWL2Norm_Sycl(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VL1Norm_Sycl(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_Sycl(sunrealtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VInvTest_Sycl(N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VConstrMask_Sycl(N_Vector c, N_Vector x,
                                                  N_Vector m);
SUNDIALS_EXPORT sunrealtype N_VMinQuotient_Sycl(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT SUNErrCode N_VLinearCombination_Sycl(int nvec, sunrealtype* c,
                                                     N_Vector* X, N_Vector Z);
SUNDIALS_EXPORT SUNErrCode N_VScaleAddMulti_Sycl(int nvec, sunrealtype* c,
                                                 N_Vector X, N_Vector* Y,
                                                 N_Vector* Z);

/* vector array operations */
SUNDIALS_EXPORT SUNErrCode N_VLinearSumVectorArray_Sycl(
  int nvec, sunrealtype a, N_Vector* X, sunrealtype b, N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VScaleVectorArray_Sycl(int nvec, sunrealtype* c,
                                                    N_Vector* X, N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VConstVectorArray_Sycl(int nvec, sunrealtype c,
                                                    N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VScaleAddMultiVectorArray_Sycl(
  int nvec, int nsum, sunrealtype* a, N_Vector* X, N_Vector** Y, N_Vector** Z);
SUNDIALS_EXPORT SUNErrCode N_VLinearCombinationVectorArray_Sycl(
  int nvec, int nsum, sunrealtype* c, N_Vector** X, N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VWrmsNormVectorArray_Sycl(int nvec, N_Vector* X,
                                                       N_Vector* W,
                                                       sunrealtype* nrm);
SUNDIALS_EXPORT SUNErrCode N_VWrmsNormMaskVectorArray_Sycl(int nvec, N_Vector* X,
                                                           N_Vector* W,
                                                           N_Vector id,
                                                           sunrealtype* nrm);

/* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT sunrealtype N_VWSqrSumLocal_Sycl(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VWSqrSumMaskLocal_Sycl(N_Vector x, N_Vector w,
                                                     N_Vector id);

/* OPTIONAL XBraid interface operations */
SUNDIALS_EXPORT SUNErrCode N_VBufSize_Sycl(N_Vector x, sunindextype* size);
SUNDIALS_EXPORT SUNErrCode N_VBufPack_Sycl(N_Vector x, void* buf);
SUNDIALS_EXPORT SUNErrCode N_VBufUnpack_Sycl(N_Vector x, void* buf);

/* OPTIONAL operations for debugging */
SUNDIALS_EXPORT void N_VPrint_Sycl(N_Vector v);
SUNDIALS_EXPORT void N_VPrintFile_Sycl(N_Vector v, FILE* outfile);

/* -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT SUNErrCode N_VEnableFusedOps_Sycl(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableLinearCombination_Sycl(N_Vector v,
                                                           sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleAddMulti_Sycl(N_Vector v,
                                                       sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableDotProdMulti_Sycl(N_Vector v,
                                                      sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableLinearSumVectorArray_Sycl(N_Vector v,
                                                              sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleVectorArray_Sycl(N_Vector v,
                                                          sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableConstVectorArray_Sycl(N_Vector v,
                                                          sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableWrmsNormVectorArray_Sycl(N_Vector v,
                                                             sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormMaskVectorArray_Sycl(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMultiVectorArray_Sycl(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombinationVectorArray_Sycl(N_Vector v,
                                                      sunbooleantype tf);

#ifdef __cplusplus
}
#endif

#endif
