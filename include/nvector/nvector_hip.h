/* -----------------------------------------------------------------
 * Programmer(s): Daniel McGreer and Cody J. Balos @ LLNL
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
 * This is the header file for the HIP implementation of the
 * NVECTOR module.
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_NVECTOR_HIP_H
#define SUNDIALS_NVECTOR_HIP_H

#include <hip/hip_runtime.h>
#include <stdio.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_hip_policies.hpp>
#include <sundials/sundials_nvector.h>
#include <sunmemory/sunmemory_hip.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * HIP implementation of N_Vector
 * -----------------------------------------------------------------
 */

struct N_VectorContent_Hip_
{
  sunindextype length;
  sunbooleantype own_helper;
  SUNMemory host_data;
  SUNMemory device_data;
  SUNHipExecPolicy* stream_exec_policy;
  SUNHipExecPolicy* reduce_exec_policy;
  SUNMemoryHelper mem_helper;
  void* priv; /* 'private' data */
};

typedef struct N_VectorContent_Hip_* N_VectorContent_Hip;

/*
 * -----------------------------------------------------------------
 * NVECTOR_HIP implementation specific functions
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_Hip(SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNew_Hip(sunindextype length, SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNewManaged_Hip(sunindextype length,
                                           SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNewWithMemHelp_Hip(sunindextype length,
                                               sunbooleantype use_managed_mem,
                                               SUNMemoryHelper helper,
                                               SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VMake_Hip(sunindextype length, sunrealtype* h_vdata_1d,
                                     sunrealtype* d_vdata_1d, SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VMakeManaged_Hip(sunindextype length,
                                            sunrealtype* vdata_1d,
                                            SUNContext sunctx);
SUNDIALS_EXPORT sunbooleantype N_VIsManagedMemory_Hip(N_Vector x);
SUNDIALS_EXPORT
SUNErrCode N_VSetKernelExecPolicy_Hip(N_Vector x,
                                      SUNHipExecPolicy* stream_exec_policy,
                                      SUNHipExecPolicy* reduce_exec_policy);
SUNDIALS_EXPORT void N_VCopyToDevice_Hip(N_Vector v);
SUNDIALS_EXPORT void N_VCopyFromDevice_Hip(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNErrCode N_VEnableFusedOps_Hip(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableLinearCombination_Hip(N_Vector v,
                                                          sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleAddMulti_Hip(N_Vector v,
                                                      sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableDotProdMulti_Hip(N_Vector v,
                                                     sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableLinearSumVectorArray_Hip(N_Vector v,
                                                             sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleVectorArray_Hip(N_Vector v,
                                                         sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableConstVectorArray_Hip(N_Vector v,
                                                         sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableWrmsNormVectorArray_Hip(N_Vector v,
                                                            sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableWrmsNormMaskVectorArray_Hip(N_Vector v,
                                                                sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMultiVectorArray_Hip(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombinationVectorArray_Hip(N_Vector v,
                                                     sunbooleantype tf);

#ifdef __cplusplus
}
#endif

#endif
