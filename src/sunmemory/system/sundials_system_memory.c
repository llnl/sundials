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
 * SUNDIALS memory helper implementation that uses the standard
 * system memory allocators.
 * ----------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_memory.h>
#include <sunmemory/sunmemory_system.h>

#include "sundials_debug.h"
#include "sundials_macros.h"

struct SUNMemoryHelper_Content_Sys_
{
  unsigned long num_allocations;
  unsigned long num_deallocations;
  size_t bytes_allocated;
  size_t bytes_high_watermark;
};

typedef struct SUNMemoryHelper_Content_Sys_ SUNMemoryHelper_Content_Sys;

#define SUNHELPER_CONTENT(h) ((SUNMemoryHelper_Content_Sys*)h->content)

SUNMemoryHelper SUNMemoryHelper_Sys(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNMemoryHelper helper;

  /* Allocate the helper */
  helper = SUNMemoryHelper_NewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Set the ops */
  helper->ops->alloc         = SUNMemoryHelper_Alloc_Sys;
  helper->ops->allocstrided  = SUNMemoryHelper_AllocStrided_Sys;
  helper->ops->dealloc       = SUNMemoryHelper_Dealloc_Sys;
  helper->ops->copy          = SUNMemoryHelper_Copy_Sys;
  helper->ops->getallocstats = SUNMemoryHelper_GetAllocStats_Sys;
  helper->ops->clone         = SUNMemoryHelper_Clone_Sys;
  helper->ops->destroy       = SUNMemoryHelper_Destroy_Sys;

  /* Attach content and ops */
  helper->content =
    (SUNMemoryHelper_Content_Sys*)malloc(sizeof(SUNMemoryHelper_Content_Sys));
  SUNAssertNull(helper->content, SUN_ERR_MALLOC_FAIL);

  SUNHELPER_CONTENT(helper)->num_allocations      = 0;
  SUNHELPER_CONTENT(helper)->num_deallocations    = 0;
  SUNHELPER_CONTENT(helper)->bytes_allocated      = 0;
  SUNHELPER_CONTENT(helper)->bytes_high_watermark = 0;

  return helper;
}

SUNErrCode SUNMemoryHelper_Alloc_Sys(SUNMemoryHelper helper, SUNMemory* memptr,
                                     size_t mem_size, SUNMemoryType mem_type,
                                     SUNDIALS_MAYBE_UNUSED void* queue)
{
  SUNFunctionBegin(helper->sunctx);

  SUNAssert(mem_type == SUNMEMTYPE_HOST || mem_type == SUNMEMTYPE_UVM,
            SUN_ERR_ARG_INCOMPATIBLE);

  SUNMemory mem = SUNMemoryNewEmpty(helper->sunctx);
  SUNCheckLastErr();

  mem->ptr   = NULL;
  mem->own   = SUNTRUE;
  mem->type  = mem_type;
  mem->bytes = mem_size;

  if (mem_type == SUNMEMTYPE_HOST || mem_type == SUNMEMTYPE_UVM)
  {
    mem->ptr = malloc(mem_size);
    SUNAssert(mem->ptr, SUN_ERR_MALLOC_FAIL);
    SUNHELPER_CONTENT(helper)->bytes_allocated += mem_size;
    SUNHELPER_CONTENT(helper)->num_allocations++;
    SUNHELPER_CONTENT(helper)->bytes_high_watermark =
      SUNMAX(SUNHELPER_CONTENT(helper)->bytes_allocated,
             SUNHELPER_CONTENT(helper)->bytes_high_watermark);
  }

  *memptr = mem;
  return SUN_SUCCESS;
}

SUNErrCode SUNMemoryHelper_AllocStrided_Sys(SUNMemoryHelper helper,
                                            SUNMemory* memptr, size_t mem_size,
                                            size_t stride,
                                            SUNMemoryType mem_type, void* queue)
{
  SUNFunctionBegin(helper->sunctx);

  SUNCheckCall(
    SUNMemoryHelper_Alloc_Sys(helper, memptr, mem_size, mem_type, queue));

  (*memptr)->stride = stride;

  return SUN_SUCCESS;
}

SUNErrCode SUNMemoryHelper_Dealloc_Sys(SUNMemoryHelper helper, SUNMemory mem,
                                       SUNDIALS_MAYBE_UNUSED void* queue)
{
  SUNFunctionBegin(helper->sunctx);

  if (mem == NULL) { return SUN_SUCCESS; }

  SUNAssert(mem->type == SUNMEMTYPE_HOST || mem->type == SUNMEMTYPE_UVM,
            SUN_ERR_ARG_INCOMPATIBLE);

  if (mem->ptr != NULL && mem->own)
  {
    if (mem->type == SUNMEMTYPE_HOST || mem->type == SUNMEMTYPE_UVM)
    {
      SUNHELPER_CONTENT(helper)->num_deallocations++;
      SUNHELPER_CONTENT(helper)->bytes_allocated -= mem->bytes;
      free(mem->ptr);
      mem->ptr = NULL;
    }
  }

  free(mem);
  return SUN_SUCCESS;
}

SUNErrCode SUNMemoryHelper_Copy_Sys(SUNMemoryHelper helper, SUNMemory dst,
                                    SUNMemory src, size_t memory_size,
                                    SUNDIALS_MAYBE_UNUSED void* queue)
{
  SUNFunctionBegin(helper->sunctx);
  SUNAssert(src->type == SUNMEMTYPE_HOST || src->type == SUNMEMTYPE_UVM,
            SUN_ERR_ARG_INCOMPATIBLE);
  SUNAssert(dst->type == SUNMEMTYPE_HOST || dst->type == SUNMEMTYPE_UVM,
            SUN_ERR_ARG_INCOMPATIBLE);
  memcpy(dst->ptr, src->ptr, memory_size);
  return SUN_SUCCESS;
}

SUNErrCode SUNMemoryHelper_GetAllocStats_Sys(
  SUNMemoryHelper helper, SUNDIALS_MAYBE_UNUSED SUNMemoryType mem_type,
  unsigned long* num_allocations, unsigned long* num_deallocations,
  size_t* bytes_allocated, size_t* bytes_high_watermark)
{
  SUNFunctionBegin(helper->sunctx);
  SUNAssert(mem_type == SUNMEMTYPE_HOST || mem_type == SUNMEMTYPE_UVM,
            SUN_ERR_ARG_INCOMPATIBLE);
  *num_allocations      = SUNHELPER_CONTENT(helper)->num_allocations;
  *num_deallocations    = SUNHELPER_CONTENT(helper)->num_deallocations;
  *bytes_allocated      = SUNHELPER_CONTENT(helper)->bytes_allocated;
  *bytes_high_watermark = SUNHELPER_CONTENT(helper)->bytes_high_watermark;
  return SUN_SUCCESS;
}

SUNMemoryHelper SUNMemoryHelper_Clone_Sys(SUNMemoryHelper helper)
{
  SUNFunctionBegin(helper->sunctx);
  SUNMemoryHelper hclone = SUNMemoryHelper_Sys(helper->sunctx);
  SUNCheckLastErrNull();
  return hclone;
}

SUNErrCode SUNMemoryHelper_Destroy_Sys(SUNMemoryHelper helper)
{
  if (helper)
  {
    if (helper->content) { free(helper->content); }
    if (helper->ops) { free(helper->ops); }
    free(helper);
  }
  return SUN_SUCCESS;
}
