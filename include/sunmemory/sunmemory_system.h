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
 * SUNDIALS system memory helper header file.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_SYSMEMORY_H
#define _SUNDIALS_SYSMEMORY_H

#include <sundials/sundials_memory.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Implementation specific functions */

SUNDIALS_EXPORT
SUNMemoryHelper SUNMemoryHelper_Sys(SUNContext sunctx);

/* SUNMemoryHelper functions */

SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_Alloc_Sys(SUNMemoryHelper helper, SUNMemory* memptr,
                                     size_t mem_size, SUNMemoryType mem_type,
                                     void* queue);

SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_AllocStrided_Sys(SUNMemoryHelper helper,
                                            SUNMemory* memptr, size_t mem_size,
                                            size_t stride,
                                            SUNMemoryType mem_type, void* queue);

SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_Dealloc_Sys(SUNMemoryHelper helper, SUNMemory mem,
                                       void* queue);

SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_Copy_Sys(SUNMemoryHelper helper, SUNMemory dst,
                                    SUNMemory src, size_t memory_size,
                                    void* queue);

SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_GetAllocStats_Sys(SUNMemoryHelper helper,
                                             SUNMemoryType mem_type,
                                             unsigned long* num_allocations,
                                             unsigned long* num_deallocations,
                                             size_t* bytes_allocated,
                                             size_t* bytes_high_watermark);

SUNDIALS_EXPORT
SUNMemoryHelper SUNMemoryHelper_Clone_Sys(SUNMemoryHelper helper);

SUNDIALS_EXPORT
SUNErrCode SUNMemoryHelper_Destroy_Sys(SUNMemoryHelper helper);

#ifdef __cplusplus
}
#endif

#endif
