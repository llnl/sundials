/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * SUNDIALS CUDA memory helper header file.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_CUDAMEMORY_H
#define _SUNDIALS_CUDAMEMORY_H

#include <cuda_runtime.h>
#include <sundials/sundials_memory.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Implementation specific functions */

SUNDIALS_EXPORT
SUNMemoryHelper SUNMemoryHelper_Cuda(SUNContext sunctx);

/* SUNMemoryHelper functions */

SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_Alloc_Cuda(SUNMemoryHelper helper, SUNMemory* memptr,
                                      size_t mem_size, SUNMemoryType mem_type,
                                      void* queue);

SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_AllocStrided_Cuda(SUNMemoryHelper helper,
                                             SUNMemory* memptr, size_t mem_size,
                                             size_t stride,
                                             SUNMemoryType mem_type, void* queue);

SUNDIALS_EXPORT SUNMemoryHelper SUNMemoryHelper_Clone_Cuda(SUNMemoryHelper helper);

SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_Dealloc_Cuda(SUNMemoryHelper helper, SUNMemory mem,
                                        void* queue);

SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_Copy_Cuda(SUNMemoryHelper helper, SUNMemory dst,
                                     SUNMemory src, size_t memory_size,
                                     void* queue);

SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_CopyAsync_Cuda(SUNMemoryHelper helper, SUNMemory dst,
                                          SUNMemory src, size_t memory_size,
                                          void* queue);

SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_Destroy_Cuda(SUNMemoryHelper helper);

SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_GetAllocStats_Cuda(SUNMemoryHelper helper,
                                              SUNMemoryType mem_type,
                                              unsigned long* num_allocations,
                                              unsigned long* num_deallocations,
                                              size_t* bytes_allocated,
                                              size_t* bytes_high_watermark);

#ifdef __cplusplus
}
#endif

#endif
