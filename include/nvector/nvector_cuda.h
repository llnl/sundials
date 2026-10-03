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

#ifndef SUNDIALS_NVECTOR_CUDA_H
#define SUNDIALS_NVECTOR_CUDA_H

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
SUNDIALS_EXPORT sunbooleantype N_VIsManagedMemory_Cuda(N_Vector x);
SUNDIALS_EXPORT
SUNErrCode N_VSetKernelExecPolicy_Cuda(N_Vector x,
                                       SUNCudaExecPolicy* stream_exec_policy,
                                       SUNCudaExecPolicy* reduce_exec_policy);
SUNDIALS_EXPORT void N_VCopyToDevice_Cuda(N_Vector v);
SUNDIALS_EXPORT void N_VCopyFromDevice_Cuda(N_Vector v);

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
