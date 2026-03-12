/* ----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ----------------------------------------------------------------
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
 * ----------------------------------------------------------------
 * SUNDIALS SYCL memory helper implementation.
 * ----------------------------------------------------------------*/

#include <cstdlib>
#include <sundials/sundials_math.h>
#include <sunmemory/sunmemory_sycl.h>

#include "sundials/priv/sundials_errors_impl.h"
#include "sundials/sundials_errors.h"
#include "sundials_debug.h"

struct SUNMemoryHelper_Content_Sycl_
{
  unsigned long num_allocations_host;
  unsigned long num_deallocations_host;
  unsigned long num_allocations_device;
  unsigned long num_deallocations_device;
  unsigned long num_allocations_pinned;
  unsigned long num_deallocations_pinned;
  unsigned long num_allocations_uvm;
  unsigned long num_deallocations_uvm;
  size_t bytes_allocated_host;
  size_t bytes_high_watermark_host;
  size_t bytes_allocated_device;
  size_t bytes_high_watermark_device;
  size_t bytes_allocated_pinned;
  size_t bytes_high_watermark_pinned;
  size_t bytes_allocated_uvm;
  size_t bytes_high_watermark_uvm;
};

typedef struct SUNMemoryHelper_Content_Sycl_ SUNMemoryHelper_Content_Sycl;

#define SUNHELPER_CONTENT(h) ((SUNMemoryHelper_Content_Sycl*)h->content)
#define GET_SYCL_QUEUE(h, q) (static_cast<::sycl::queue*>(q ? q : h->queue))

SUNMemoryHelper SUNMemoryHelper_Sycl(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  // Allocate the helper
  SUNMemoryHelper helper = SUNMemoryHelper_NewEmpty(sunctx);
  SUNCheckLastErrNull();

  // Set the ops
  helper->ops->alloc         = SUNMemoryHelper_Alloc_Sycl;
  helper->ops->dealloc       = SUNMemoryHelper_Dealloc_Sycl;
  helper->ops->clone         = SUNMemoryHelper_Clone_Sycl;
  helper->ops->getallocstats = SUNMemoryHelper_GetAllocStats_Sycl;
  helper->ops->copy          = SUNMemoryHelper_Copy_Sycl;
  helper->ops->copyasync     = SUNMemoryHelper_CopyAsync_Sycl;
  helper->ops->destroy       = SUNMemoryHelper_Destroy_Sycl;

  // Attach content
  helper->content =
    (SUNMemoryHelper_Content_Sycl*)malloc(sizeof(SUNMemoryHelper_Content_Sycl));
  SUNAssertNull(helper->content, SUN_ERR_MALLOC_FAIL);

  SUNHELPER_CONTENT(helper)->num_allocations_host        = 0;
  SUNHELPER_CONTENT(helper)->num_deallocations_host      = 0;
  SUNHELPER_CONTENT(helper)->bytes_allocated_host        = 0;
  SUNHELPER_CONTENT(helper)->bytes_high_watermark_host   = 0;
  SUNHELPER_CONTENT(helper)->num_allocations_device      = 0;
  SUNHELPER_CONTENT(helper)->num_deallocations_device    = 0;
  SUNHELPER_CONTENT(helper)->bytes_allocated_device      = 0;
  SUNHELPER_CONTENT(helper)->bytes_high_watermark_device = 0;
  SUNHELPER_CONTENT(helper)->num_allocations_pinned      = 0;
  SUNHELPER_CONTENT(helper)->num_deallocations_pinned    = 0;
  SUNHELPER_CONTENT(helper)->bytes_allocated_pinned      = 0;
  SUNHELPER_CONTENT(helper)->bytes_high_watermark_pinned = 0;
  SUNHELPER_CONTENT(helper)->num_allocations_uvm         = 0;
  SUNHELPER_CONTENT(helper)->num_deallocations_uvm       = 0;
  SUNHELPER_CONTENT(helper)->bytes_allocated_uvm         = 0;
  SUNHELPER_CONTENT(helper)->bytes_high_watermark_uvm    = 0;

  return helper;
}

SUNMemoryHelper SUNMemoryHelper_Clone_Sycl(SUNMemoryHelper helper)
{
  SUNFunctionBegin(helper->sunctx);
  SUNMemoryHelper hclone = SUNMemoryHelper_Sycl(helper->sunctx);
  SUNCheckLastErrNull();
  return hclone;
}

SUNErrCode SUNMemoryHelper_Alloc_Sycl(SUNMemoryHelper helper, SUNMemory* memptr,
                                      size_t mem_size, SUNMemoryType mem_type,
                                      void* queue)
{
  SUNFunctionBegin(helper->sunctx);

  // Check inputs
  ::sycl::queue* sycl_queue = GET_SYCL_QUEUE(helper, queue);
  SUNAssert(sycl_queue, SUN_ERR_ARG_CORRUPT);

  // Allocate the memory struct
  SUNMemory mem = SUNMemoryNewEmpty(helper->sunctx);
  SUNCheckLastErr();

  // Initialize the memory content
  mem->ptr   = nullptr;
  mem->own   = SUNTRUE;
  mem->type  = mem_type;
  mem->bytes = mem_size;

  // Allocate the data pointer
  if (mem_type == SUNMEMTYPE_HOST)
  {
    mem->ptr = malloc(mem_size);
    if (!(mem->ptr))
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in SUNMemoryHelper_Alloc_Sycl: malloc returned NULL\n");
      free(mem);
      return SUN_ERR_EXT_FAIL;
    }
    else
    {
      SUNHELPER_CONTENT(helper)->bytes_allocated_host += mem_size;
      SUNHELPER_CONTENT(helper)->num_allocations_host++;
      SUNHELPER_CONTENT(helper)->bytes_high_watermark_host =
        SUNMAX(SUNHELPER_CONTENT(helper)->bytes_allocated_host,
               SUNHELPER_CONTENT(helper)->bytes_high_watermark_host);
    }
  }
  else if (mem_type == SUNMEMTYPE_PINNED)
  {
    mem->ptr = ::sycl::malloc_host(mem_size, *sycl_queue);
    if (!(mem->ptr))
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in SUNMemoryHelper_Alloc_Sycl: malloc_host returned NULL\n");
      free(mem);
      return SUN_ERR_EXT_FAIL;
    }
    else
    {
      SUNHELPER_CONTENT(helper)->bytes_allocated_pinned += mem_size;
      SUNHELPER_CONTENT(helper)->num_allocations_pinned++;
      SUNHELPER_CONTENT(helper)->bytes_high_watermark_pinned =
        SUNMAX(SUNHELPER_CONTENT(helper)->bytes_allocated_pinned,
               SUNHELPER_CONTENT(helper)->bytes_high_watermark_pinned);
    }
  }
  else if (mem_type == SUNMEMTYPE_DEVICE)
  {
    mem->ptr = ::sycl::malloc_device(mem_size, *sycl_queue);
    if (!(mem->ptr))
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in SUNMemoryHelper_Alloc_Sycl: malloc_device returned NULL\n");
      free(mem);
      return SUN_ERR_EXT_FAIL;
    }
    else
    {
      SUNHELPER_CONTENT(helper)->bytes_allocated_device += mem_size;
      SUNHELPER_CONTENT(helper)->num_allocations_device++;
      SUNHELPER_CONTENT(helper)->bytes_high_watermark_device =
        SUNMAX(SUNHELPER_CONTENT(helper)->bytes_allocated_device,
               SUNHELPER_CONTENT(helper)->bytes_high_watermark_device);
    }
  }
  else if (mem_type == SUNMEMTYPE_UVM)
  {
    mem->ptr = ::sycl::malloc_shared(mem_size, *sycl_queue);
    if (!(mem->ptr))
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in SUNMemoryHelper_Alloc_Sycl: malloc_shared returned NULL\n");
      free(mem);
      return SUN_ERR_EXT_FAIL;
    }
    else
    {
      SUNHELPER_CONTENT(helper)->bytes_allocated_uvm += mem_size;
      SUNHELPER_CONTENT(helper)->num_allocations_uvm++;
      SUNHELPER_CONTENT(helper)->bytes_high_watermark_uvm =
        SUNMAX(SUNHELPER_CONTENT(helper)->bytes_allocated_uvm,
               SUNHELPER_CONTENT(helper)->bytes_high_watermark_uvm);
    }
  }
  else
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in SUNMemoryHelper_Alloc_Sycl: unknown memory type\n");
    free(mem);
    return SUN_ERR_ARG_OUTOFRANGE;
  }

  *memptr = mem;
  return SUN_SUCCESS;
}

SUNErrCode SUNMemoryHelper_AllocStrided_Sycl(SUNMemoryHelper helper,
                                             SUNMemory* memptr, size_t mem_size,
                                             size_t stride,
                                             SUNMemoryType mem_type, void* queue)
{
  SUNFunctionBegin(helper->sunctx);

  SUNCheckCall(
    SUNMemoryHelper_Alloc_Sycl(helper, memptr, mem_size, mem_type, queue));

  (*memptr)->stride = stride;

  return SUN_SUCCESS;
}

SUNErrCode SUNMemoryHelper_Dealloc_Sycl(SUNMemoryHelper helper, SUNMemory mem,
                                        void* queue)
{
  SUNFunctionBegin(helper->sunctx);

  if (!mem) { return SUN_SUCCESS; }

  if (mem->ptr && mem->own)
  {
    ::sycl::queue* sycl_queue = GET_SYCL_QUEUE(helper, queue);
    SUNAssert(sycl_queue, SUN_ERR_ARG_CORRUPT);

    if (mem->type == SUNMEMTYPE_HOST)
    {
      SUNHELPER_CONTENT(helper)->num_deallocations_host++;
      SUNHELPER_CONTENT(helper)->bytes_allocated_host -= mem->bytes;
      free(mem->ptr);
      mem->ptr = nullptr;
    }
    else if (mem->type == SUNMEMTYPE_PINNED || mem->type == SUNMEMTYPE_DEVICE ||
             mem->type == SUNMEMTYPE_UVM)
    {
      if (mem->type == SUNMEMTYPE_PINNED)
      {
        SUNHELPER_CONTENT(helper)->num_deallocations_pinned++;
        SUNHELPER_CONTENT(helper)->bytes_allocated_pinned -= mem->bytes;
      }
      else if (mem->type == SUNMEMTYPE_DEVICE)
      {
        SUNHELPER_CONTENT(helper)->num_deallocations_device++;
        SUNHELPER_CONTENT(helper)->bytes_allocated_device -= mem->bytes;
      }
      else if (mem->type == SUNMEMTYPE_UVM)
      {
        SUNHELPER_CONTENT(helper)->num_deallocations_uvm++;
        SUNHELPER_CONTENT(helper)->bytes_allocated_uvm -= mem->bytes;
      }
      ::sycl::free(mem->ptr, *sycl_queue);
      mem->ptr = nullptr;
    }
    else
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in SUNMemoryHelper_Dealloc_Sycl: unknown memory type\n");
      return SUN_ERR_OUTOFRANGE;
    }
  }

  free(mem);
  return SUN_SUCCESS;
}

SUNErrCode SUNMemoryHelper_Copy_Sycl(SUNMemoryHelper helper, SUNMemory dst,
                                     SUNMemory src, size_t memory_size,
                                     void* queue)
{
  SUNFunctionBegin(helper->sunctx);
  ::sycl::queue* sycl_queue = GET_SYCL_QUEUE(helper, queue);
  SUNAssert(sycl_queue, SUN_ERR_ARG_CORRUPT);

  if (SUNMemoryHelper_CopyAsync_Sycl(helper, dst, src, memory_size, queue))
  {
    return SUN_ERR_EXT_FAIL;
  }
  sycl_queue->wait_and_throw();
  return SUN_SUCCESS;
}

SUNErrCode SUNMemoryHelper_CopyAsync_Sycl(SUNMemoryHelper helper, SUNMemory dst,
                                          SUNMemory src, size_t memory_size,
                                          void* queue)
{
  SUNFunctionBegin(helper->sunctx);
  ::sycl::queue* sycl_queue = GET_SYCL_QUEUE(helper, queue);
  SUNAssert(sycl_queue, SUN_ERR_ARG_CORRUPT);

  if (src->type == SUNMEMTYPE_HOST && dst->type == SUNMEMTYPE_HOST)
  {
    memcpy(dst->ptr, src->ptr, memory_size);
  }
  else { sycl_queue->memcpy(dst->ptr, src->ptr, memory_size); }
  return SUN_SUCCESS;
}

SUNErrCode SUNMemoryHelper_Destroy_Sycl(SUNMemoryHelper helper)
{
  if (helper)
  {
    if (helper->content) { free(helper->content); }
    if (helper->ops) { free(helper->ops); }
    free(helper);
  }
  return SUN_SUCCESS;
}

SUNErrCode SUNMemoryHelper_GetAllocStats_Sycl(SUNMemoryHelper helper,
                                              SUNMemoryType mem_type,
                                              unsigned long* num_allocations,
                                              unsigned long* num_deallocations,
                                              size_t* bytes_allocated,
                                              size_t* bytes_high_watermark)
{
  if (mem_type == SUNMEMTYPE_HOST)
  {
    *num_allocations   = SUNHELPER_CONTENT(helper)->num_allocations_host;
    *num_deallocations = SUNHELPER_CONTENT(helper)->num_deallocations_host;
    *bytes_allocated   = SUNHELPER_CONTENT(helper)->bytes_allocated_host;
    *bytes_high_watermark = SUNHELPER_CONTENT(helper)->bytes_high_watermark_host;
  }
  else if (mem_type == SUNMEMTYPE_PINNED)
  {
    *num_allocations   = SUNHELPER_CONTENT(helper)->num_allocations_pinned;
    *num_deallocations = SUNHELPER_CONTENT(helper)->num_deallocations_pinned;
    *bytes_allocated   = SUNHELPER_CONTENT(helper)->bytes_allocated_pinned;
    *bytes_high_watermark = SUNHELPER_CONTENT(helper)->bytes_high_watermark_pinned;
  }
  else if (mem_type == SUNMEMTYPE_DEVICE)
  {
    *num_allocations   = SUNHELPER_CONTENT(helper)->num_allocations_device;
    *num_deallocations = SUNHELPER_CONTENT(helper)->num_deallocations_device;
    *bytes_allocated   = SUNHELPER_CONTENT(helper)->bytes_allocated_device;
    *bytes_high_watermark = SUNHELPER_CONTENT(helper)->bytes_high_watermark_device;
  }
  else if (mem_type == SUNMEMTYPE_UVM)
  {
    *num_allocations      = SUNHELPER_CONTENT(helper)->num_allocations_uvm;
    *num_deallocations    = SUNHELPER_CONTENT(helper)->num_deallocations_uvm;
    *bytes_allocated      = SUNHELPER_CONTENT(helper)->bytes_allocated_uvm;
    *bytes_high_watermark = SUNHELPER_CONTENT(helper)->bytes_high_watermark_uvm;
  }
  else { return SUN_ERR_ARG_OUTOFRANGE; }
  return SUN_SUCCESS;
}
