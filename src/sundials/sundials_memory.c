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
 * SUNDIALS memory helper.
 * ----------------------------------------------------------------*/

#include "sundials/sundials_memory.h"

#include <string.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>

#include "sundials/sundials_errors.h"
#include "sundials/sundials_types.h"
#include "sundials_debug.h"

#if defined(SUNDIALS_BUILD_WITH_PROFILING)
static SUNProfiler getSUNProfiler(SUNMemoryHelper H)
{
  return (H->sunctx->profiler);
}
#endif

SUNMemory SUNMemoryNewEmpty(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNMemory mem = NULL;

  mem = (SUNMemory)malloc(sizeof(struct SUNMemory_));
  SUNAssertNull(mem, SUN_ERR_MALLOC_FAIL);

  mem->bytes  = 0;
  mem->stride = 1;

  return (mem);
}

SUNMemoryHelper SUNMemoryHelper_NewEmpty(SUNContext sunctx)
{
  if (sunctx == NULL) { return NULL; }

  SUNFunctionBegin(sunctx);
  SUNMemoryHelper helper = NULL;

  helper = (SUNMemoryHelper)malloc(sizeof(struct SUNMemoryHelper_));
  SUNAssertNull(helper, SUN_ERR_MALLOC_FAIL);

  helper->ops = (SUNMemoryHelper_Ops)malloc(sizeof(struct SUNMemoryHelper_Ops_));
  SUNAssertNull(helper->ops, SUN_ERR_MALLOC_FAIL);

  /* Set all ops to NULL */
  memset(helper->ops, 0, sizeof(struct SUNMemoryHelper_Ops_));
  helper->content = NULL;
  helper->queue   = NULL;
  helper->sunctx  = sunctx;

  return helper;
}

SUNErrCode SUNMemoryHelper_CopyOps(SUNMemoryHelper src, SUNMemoryHelper dst)
{
  memcpy(dst->ops, src->ops, sizeof(struct SUNMemoryHelper_Ops_));
  return SUN_SUCCESS;
}

sunbooleantype SUNMemoryHelper_ImplementsRequiredOps(SUNMemoryHelper helper)
{
  if (helper->ops->alloc == NULL || helper->ops->dealloc == NULL ||
      helper->ops->copy == NULL)
  {
    return SUNFALSE;
  }
  return SUNTRUE;
}

SUNMemory SUNMemoryHelper_Alias(SUNMemoryHelper helper, SUNMemory mem)
{
  SUNFunctionBegin(helper->sunctx);

  SUNMemory alias = SUNMemoryNewEmpty(helper->sunctx);
  SUNCheckLastErrNull();

  alias->ptr  = mem->ptr;
  alias->type = mem->type;
  alias->own  = SUNFALSE;

  return alias;
}

SUNMemory SUNMemoryHelper_Wrap(SUNMemoryHelper helper, void* ptr,
                               SUNMemoryType mem_type)
{
  SUNFunctionBegin(helper->sunctx);

  SUNCheckNull(mem_type == SUNMEMTYPE_HOST || mem_type == SUNMEMTYPE_PINNED ||
                 mem_type == SUNMEMTYPE_DEVICE || mem_type == SUNMEMTYPE_UVM,
               SUN_ERR_ARG_OUTOFRANGE);

  SUNMemory mem = SUNMemoryNewEmpty(helper->sunctx);
  SUNCheckLastErrNull();

  mem->ptr  = ptr;
  mem->own  = SUNFALSE;
  mem->type = mem_type;

  return mem;
}

SUNErrCode SUNMemoryHelper_GetAllocStats(SUNMemoryHelper helper,
                                         SUNMemoryType mem_type,
                                         unsigned long* num_allocations,
                                         unsigned long* num_deallocations,
                                         size_t* bytes_allocated,
                                         size_t* bytes_high_watermark)
{
  SUNErrCode ier = SUN_SUCCESS;
  SUNFunctionBegin(helper->sunctx);
  SUNAssert(helper->ops->getallocstats, SUN_ERR_NOT_IMPLEMENTED);
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(helper));
  ier = helper->ops->getallocstats(helper, mem_type, num_allocations,
                                   num_deallocations, bytes_allocated,
                                   bytes_high_watermark);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(helper));
  return (ier);
}

SUNErrCode SUNMemoryHelper_Alloc(SUNMemoryHelper helper, SUNMemory* memptr,
                                 size_t mem_size, SUNMemoryType mem_type,
                                 void* queue)
{
  SUNErrCode ier = SUN_SUCCESS;
  SUNFunctionBegin(helper->sunctx);
  SUNAssert(helper->ops->alloc, SUN_ERR_NOT_IMPLEMENTED);
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(helper));
  ier = helper->ops->alloc(helper, memptr, mem_size, mem_type, queue);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(helper));
  return ier;
}

SUNErrCode SUNMemoryHelper_AllocStrided(SUNMemoryHelper helper, SUNMemory* memptr,
                                        size_t mem_size, size_t stride,
                                        SUNMemoryType mem_type, void* queue)
{
  SUNErrCode ier = SUN_SUCCESS;
  SUNFunctionBegin(helper->sunctx);
  SUNAssert(helper->ops->allocstrided, SUN_ERR_NOT_IMPLEMENTED);
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(helper));
  ier = helper->ops->allocstrided(helper, memptr, mem_size, stride, mem_type,
                                  queue);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(helper));
  return ier;
}

SUNErrCode SUNMemoryHelper_Dealloc(SUNMemoryHelper helper, SUNMemory mem,
                                   void* queue)
{
  SUNErrCode ier = SUN_SUCCESS;
  SUNFunctionBegin(helper->sunctx);
  SUNAssert(helper->ops->dealloc, SUN_ERR_NOT_IMPLEMENTED);
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(helper));
  if (!mem) { ier = SUN_SUCCESS; }
  else { ier = helper->ops->dealloc(helper, mem, queue); }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(helper));
  return ier;
}

SUNErrCode SUNMemoryHelper_Copy(SUNMemoryHelper helper, SUNMemory dst,
                                SUNMemory src, size_t memory_size, void* queue)
{
  SUNErrCode ier = SUN_SUCCESS;
  SUNFunctionBegin(helper->sunctx);
  SUNAssert(helper->ops->copy, SUN_ERR_NOT_IMPLEMENTED);
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(helper));
  ier = helper->ops->copy(helper, dst, src, memory_size, queue);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(helper));
  return ier;
}

SUNErrCode SUNMemoryHelper_CopyAsync(SUNMemoryHelper helper, SUNMemory dst,
                                     SUNMemory src, size_t memory_size,
                                     void* queue)
{
  SUNErrCode ier = SUN_SUCCESS;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(helper));
  if (!helper->ops->copyasync)
  {
    ier = SUNMemoryHelper_Copy(helper, dst, src, memory_size, queue);
  }
  else { ier = helper->ops->copyasync(helper, dst, src, memory_size, queue); }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(helper));
  return ier;
}

SUNErrCode SUNMemoryHelper_Destroy(SUNMemoryHelper helper)
{
  SUNErrCode err = SUN_SUCCESS;
  SUNFunctionBegin(helper->sunctx);

  if (!helper) { return err; }

  if (helper->content)
  {
    SUNAssert(helper->ops->destroy, SUN_ERR_NOT_IMPLEMENTED);
  }

  if (helper->ops->destroy)
  {
    /* user helper defined destroy */
    err = helper->ops->destroy(helper);
  }
  else
  {
    /* default destroy */
    free(helper->ops);
    free(helper);
  }

  return err;
}

SUNMemoryHelper SUNMemoryHelper_Clone(SUNMemoryHelper helper)
{
  SUNFunctionBegin(helper->sunctx);
  if (!helper->ops->clone)
  {
    if (helper->content)
    {
      SUNAssertNull(helper->ops->clone, SUN_ERR_NOT_IMPLEMENTED);
      return (NULL);
    }
    else
    {
      SUNMemoryHelper hclone = SUNMemoryHelper_NewEmpty(helper->sunctx);
      if (hclone) { SUNMemoryHelper_CopyOps(helper, hclone); }
      return (hclone);
    }
  }
  else { return (helper->ops->clone(helper)); }
}

SUNErrCode SUNMemoryHelper_SetDefaultQueue(SUNMemoryHelper helper, void* queue)
{
  SUNFunctionBegin(helper->sunctx);

  helper->queue = queue;

  return SUN_SUCCESS;
}
