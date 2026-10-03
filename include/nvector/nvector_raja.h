/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles, Cody J. Balos, Daniel McGreer @ LLNL
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
 * This is the header file for the RAJA implementation of the
 * NVECTOR module.
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_NVECTOR_RAJA_H
#define SUNDIALS_NVECTOR_RAJA_H

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

struct N_VectorContent_Raja_
{
  sunindextype length;
  sunbooleantype own_helper;
  SUNMemory host_data;
  SUNMemory device_data;
  SUNMemoryHelper mem_helper;
  void* priv; /* 'private' data */
};

typedef struct N_VectorContent_Raja_* N_VectorContent_Raja;

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
SUNDIALS_EXPORT sunbooleantype N_VIsManagedMemory_Raja(N_Vector x);
SUNDIALS_EXPORT void N_VCopyToDevice_Raja(N_Vector v);
SUNDIALS_EXPORT void N_VCopyFromDevice_Raja(N_Vector v);

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
