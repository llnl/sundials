/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles, Cody J. Balos, Daniel McGreer @ LLNL
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
 * This is the header file for the RAJA implementation of the
 * NVECTOR module.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_RAJA_H
#define _NVECTOR_RAJA_H

#include <stdio.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_memory.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * RAJA implementation of N_Vector
 * -----------------------------------------------------------------
 */

/* RAJA implementation of the N_Vector 'content' structure
   contains the length of the vector, pointers to host and device
   arrays of 'sunrealtype' components, a flag indicating ownership of
   the data, and a private data pointer  */

struct _N_VectorContent_Raja
{
  sunindextype length;
  sunbooleantype own_helper;
  SUNMemory host_data;
  SUNMemory device_data;
  SUNMemoryHelper mem_helper;
  void* priv; /* 'private' data */
};

typedef struct _N_VectorContent_Raja* N_VectorContent_Raja;

/*
 * -----------------------------------------------------------------
 * NVECTOR_RAJA implementation specific functions
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_Raja(SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNew_Raja(sunindextype length, SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNewManaged_Raja(sunindextype length,
                                            SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNewWithMemHelp_Raja(sunindextype length,
                                                sunbooleantype use_managed_mem,
                                                SUNMemoryHelper helper,
                                                SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VMake_Raja(sunindextype length,
                                      sunrealtype* h_vdata_1d,
                                      sunrealtype* d_vdata_1d, SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VMakeManaged_Raja(sunindextype length,
                                             sunrealtype* vdata_1d,
                                             SUNContext sunctx);
SUNDIALS_EXPORT void N_VSetHostArrayPointer_Raja(sunrealtype* h_vdata_1d,
                                                 N_Vector v);
SUNDIALS_EXPORT void N_VSetDeviceArrayPointer_Raja(sunrealtype* d_vdata_1d,
                                                   N_Vector v);
SUNDIALS_EXPORT sunbooleantype N_VIsManagedMemory_Raja(N_Vector x);
SUNDIALS_EXPORT void N_VCopyToDevice_Raja(N_Vector v);
SUNDIALS_EXPORT void N_VCopyFromDevice_Raja(N_Vector v);

static inline sunindextype N_VGetLength_Raja(N_Vector x)
{
  N_VectorContent_Raja content = (N_VectorContent_Raja)x->content;
  return content->length;
}

static inline sunrealtype* N_VGetHostArrayPointer_Raja(N_Vector x)
{
  N_VectorContent_Raja content = (N_VectorContent_Raja)x->content;
  return (content->host_data == NULL ? NULL
                                     : (sunrealtype*)content->host_data->ptr);
}

static inline sunrealtype* N_VGetDeviceArrayPointer_Raja(N_Vector x)
{
  N_VectorContent_Raja content = (N_VectorContent_Raja)x->content;
  return (content->device_data == NULL ? NULL
                                       : (sunrealtype*)content->device_data->ptr);
}

/*
 * -----------------------------------------------------------------
 * NVECTOR API functions
 * -----------------------------------------------------------------
 */

static inline N_Vector_ID N_VGetVectorID_Raja(N_Vector v)
{
  return SUNDIALS_NVEC_RAJA;
}

SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Raja(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_Raja(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_Raja(N_Vector v);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
void N_VSpace_Raja(N_Vector v, sunindextype* lrw, sunindextype* liw);
SUNDIALS_EXPORT void N_VSetArrayPointer_Raja(sunrealtype* v_data, N_Vector v);

/* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum_Raja(sunrealtype a, N_Vector x, sunrealtype b,
                                       N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_Raja(sunrealtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_Raja(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_Raja(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_Raja(sunrealtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_Raja(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_Raja(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_Raja(N_Vector x, sunrealtype b, N_Vector z);
SUNDIALS_EXPORT sunrealtype N_VDotProd_Raja(N_Vector x, N_Vector y);
SUNDIALS_EXPORT sunrealtype N_VMaxNorm_Raja(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWrmsNorm_Raja(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VWrmsNormMask_Raja(N_Vector x, N_Vector w,
                                                 N_Vector id);
SUNDIALS_EXPORT sunrealtype N_VMin_Raja(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWL2Norm_Raja(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VL1Norm_Raja(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_Raja(sunrealtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VInvTest_Raja(N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VConstrMask_Raja(N_Vector c, N_Vector x,
                                                  N_Vector m);
SUNDIALS_EXPORT sunrealtype N_VMinQuotient_Raja(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT SUNErrCode N_VLinearCombination_Raja(int nvec, sunrealtype* c,
                                                     N_Vector* X, N_Vector z);
SUNDIALS_EXPORT SUNErrCode N_VScaleAddMulti_Raja(int nvec, sunrealtype* c,
                                                 N_Vector x, N_Vector* Y,
                                                 N_Vector* Z);

/* vector array operations */
SUNDIALS_EXPORT SUNErrCode N_VLinearSumVectorArray_Raja(
  int nvec, sunrealtype a, N_Vector* X, sunrealtype b, N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VScaleVectorArray_Raja(int nvec, sunrealtype* c,
                                                    N_Vector* X, N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VConstVectorArray_Raja(int nvec, sunrealtype c,
                                                    N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VScaleAddMultiVectorArray_Raja(
  int nvec, int nsum, sunrealtype* a, N_Vector* X, N_Vector** Y, N_Vector** Z);
SUNDIALS_EXPORT SUNErrCode N_VLinearCombinationVectorArray_Raja(
  int nvec, int nsum, sunrealtype* c, N_Vector** X, N_Vector* Z);

/* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT sunrealtype N_VWSqrSumLocal_Raja(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VWSqrSumMaskLocal_Raja(N_Vector x, N_Vector w,
                                                     N_Vector id);

/* OPTIONAL XBraid interface operations */
SUNDIALS_EXPORT SUNErrCode N_VBufSize_Raja(N_Vector x, sunindextype* size);
SUNDIALS_EXPORT SUNErrCode N_VBufPack_Raja(N_Vector x, void* buf);
SUNDIALS_EXPORT SUNErrCode N_VBufUnpack_Raja(N_Vector x, void* buf);

/* OPTIONAL operations for debugging */
SUNDIALS_EXPORT void N_VPrint_Raja(N_Vector v);
SUNDIALS_EXPORT void N_VPrintFile_Raja(N_Vector v, FILE* outfile);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNErrCode N_VEnableFusedOps_Raja(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableLinearCombination_Raja(N_Vector v,
                                                           sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleAddMulti_Raja(N_Vector v,
                                                       sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableLinearSumVectorArray_Raja(N_Vector v,
                                                              sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleVectorArray_Raja(N_Vector v,
                                                          sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableConstVectorArray_Raja(N_Vector v,
                                                          sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMultiVectorArray_Raja(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombinationVectorArray_Raja(N_Vector v,
                                                      sunbooleantype tf);

#ifdef __cplusplus
}
#endif

#endif
