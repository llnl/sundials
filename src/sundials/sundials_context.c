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
 * SUNDIALS context class. A context object holds data that all
 * SUNDIALS objects in a simulation share.
 * ----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/priv/sundials_context_impl.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_logger.h>
#include <sundials/sundials_profiler.h>
#include <sundials/sundials_types.h>

#include "sundials/sundials_allocator.h"
#include "sundials/sundials_memory.h"
#include "sundials_adiak_metadata.h"
#include "sundials_macros.h"
#include "sundials_allocator_system.h"

/* Forward declaration of function used to destroy any data allocated for Python */
#if defined(SUNDIALS_ENABLE_PYTHON)
void SUNContextFunctionTable_Destroy(void* ptr);
#endif

SUNErrCode SUNContext_Create(SUNComm comm, SUNContext* sunctx_out)
{
  SUNErrCode err              = SUN_SUCCESS;
  SUNProfiler profiler        = NULL;
  SUNLogger logger            = NULL;
  SUNContext sunctx           = NULL;
  SUNErrHandler eh            = NULL;
  SUNAllocator host_allocator = NULL;

  /* Initialize output */
  *sunctx_out = NULL;
  sunctx      = (SUNContext)malloc(sizeof(struct SUNContext_));

  /* Cannot assert or log since the SUNContext is not yet created */
  if (!sunctx) { return SUN_ERR_MALLOC_FAIL; }

  /* Now we can assert and log errors */

  SUNFunctionBegin(sunctx);

#ifdef SUNDIALS_ADIAK_ENABLED
  adiak_init(&comm);
  sunAdiakCollectMetadata();
#endif

  do {
#if SUNDIALS_LOGGING_LEVEL > 0
#if SUNDIALS_MPI_ENABLED
    err = SUNLogger_CreateFromEnv(comm, &logger);
    SUNCheckCallNoRet(err);
    if (err) { break; }
#else
    err = SUNLogger_CreateFromEnv(SUN_COMM_NULL, &logger);
    SUNCheckCallNoRet(err);
    if (err) { break; }
#endif
#else
    err = SUNLogger_Create(SUN_COMM_NULL, 0, &logger);
    SUNCheckCallNoRet(err);
    if (err) { break; }
    err = SUNLogger_SetErrorFilename(logger, "");
    SUNCheckCallNoRet(err);
    if (err) { break; }
    err = SUNLogger_SetWarningFilename(logger, "");
    SUNCheckCallNoRet(err);
    if (err) { break; }
    err = SUNLogger_SetInfoFilename(logger, "");
    SUNCheckCallNoRet(err);
    if (err) { break; }
    err = SUNLogger_SetDebugFilename(logger, "");
    SUNCheckCallNoRet(err);
    if (err) { break; }
#endif

#if defined(SUNDIALS_ENABLE_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
    err = SUNProfiler_Create(comm, "SUNContext Default", &profiler);
    SUNCheckCallNoRet(err);
    if (err) { break; }
#endif

    err = SUNErrHandler_Create(SUNLogErrHandlerFn, NULL, &eh);
    SUNCheckCallNoRet(err);
    if (err) { break; }

    err = SUNAllocator_Create_System(&host_allocator);
    SUNCheckCallNoRet(err);
    if (err) { break; }

    sunctx->python             = NULL;
    sunctx->logger             = logger;
    sunctx->own_logger         = logger != NULL;
    sunctx->profiler           = profiler;
    sunctx->own_profiler       = profiler != NULL;
    sunctx->last_err           = SUN_SUCCESS;
    sunctx->err_handler        = eh;
    sunctx->host_allocator     = host_allocator;
    sunctx->own_host_allocator = SUNTRUE;
    sunctx->comm               = comm;
  }
  while (0);

  if (err)
  {
#if defined(SUNDIALS_ENABLE_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
    SUNCheckCallNoRet(SUNProfiler_Free(&profiler));
#endif
    SUNCheckCallNoRet(SUNLogger_Destroy(&logger));
    SUNCheckCallNoRet(SUNAllocator_Destroy(&host_allocator));
    free(sunctx);
  }
  else { *sunctx_out = sunctx; }

  return err;
}

SUNErrCode SUNContext_GetLastError(SUNContext sunctx)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);
  SUNErrCode err   = sunctx->last_err;
  sunctx->last_err = SUN_SUCCESS;
  return err;
}

SUNErrCode SUNContext_PeekLastError(SUNContext sunctx)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);
  return sunctx->last_err;
}

SUNErrCode SUNContext_PushErrHandler(SUNContext sunctx, SUNErrHandlerFn err_fn,
                                     void* err_user_data)
{
  if (!sunctx || !err_fn) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);
  SUNErrHandler new_err_handler = NULL;
  if (SUNErrHandler_Create(err_fn, err_user_data, &new_err_handler))
  {
    return SUN_ERR_CORRUPT;
  }
  new_err_handler->previous = sunctx->err_handler;
  sunctx->err_handler       = new_err_handler;
  return SUN_SUCCESS;
}

SUNErrCode SUNContext_PopErrHandler(SUNContext sunctx)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);
  if (sunctx->err_handler)
  {
    SUNErrHandler eh = sunctx->err_handler;
    if (sunctx->err_handler->previous)
    {
      sunctx->err_handler = sunctx->err_handler->previous;
    }
    else { sunctx->err_handler = NULL; }
    SUNErrHandler_Destroy(&eh);
  }
  return SUN_SUCCESS;
}

SUNErrCode SUNContext_ClearErrHandlers(SUNContext sunctx)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);
  while (sunctx->err_handler != NULL)
  {
    SUNCheckCall(SUNContext_PopErrHandler(sunctx));
  }
  return SUN_SUCCESS;
}

SUNErrCode SUNContext_GetProfiler(SUNContext sunctx, SUNProfiler* profiler)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);

#ifdef SUNDIALS_ENABLE_PROFILING
  /* get profiler */
  *profiler = sunctx->profiler;
#else
  *profiler = NULL;
#endif

  return SUN_SUCCESS;
}

SUNErrCode SUNContext_SetProfiler(SUNContext sunctx, SUNProfiler profiler)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);

#ifdef SUNDIALS_ENABLE_PROFILING
  /* free any existing profiler */
  if (sunctx->profiler && sunctx->own_profiler)
  {
    SUNCheckCall(SUNProfiler_Free(&(sunctx->profiler)));
    sunctx->profiler = NULL;
  }

  /* set profiler */
  sunctx->profiler     = profiler;
  sunctx->own_profiler = SUNFALSE;
#else
  /* silence warnings when profiling is disabled */
  ((void)profiler);
#endif

  return SUN_SUCCESS;
}

SUNErrCode SUNContext_GetLogger(SUNContext sunctx, SUNLogger* logger)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);

  /* get logger */
  *logger = sunctx->logger;
  return SUN_SUCCESS;
}

SUNErrCode SUNContext_SetLogger(SUNContext sunctx, SUNLogger logger)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }

  SUNFunctionBegin(sunctx);

  /* free any existing logger */
  if (sunctx->logger && sunctx->own_logger)
  {
    if (SUNLogger_Destroy(&(sunctx->logger))) { return SUN_ERR_DESTROY_FAIL; }
    sunctx->logger = NULL;
  }

  /* set logger */
  sunctx->logger     = logger;
  sunctx->own_logger = SUNFALSE;

  return SUN_SUCCESS;
}

SUNErrCode SUNContext_SetAllocator(SUNContext sunctx, SUNMemoryType type,
                                   SUNAllocator allocator)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }
  SUNFunctionBegin(sunctx);

  SUNAssert(allocator, SUN_ERR_SUNCTX_CORRUPT);

  switch (type)
  {
    case SUNMEMTYPE_HOST:
      sunctx->host_allocator = allocator;
      break;
    case SUNMEMTYPE_DEVICE:
      sunctx->device_allocator = allocator;
      break;
    case SUNMEMTYPE_UVM:
      sunctx->uvm_allocator = allocator;
      break;
    case SUNMEMTYPE_PINNED:
      sunctx->pinned_allocator = allocator;
      break;
    default:
      return SUN_ERR_ARG_OUTOFRANGE;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNContext_GetAllocator(SUNContext sunctx, SUNMemoryType type,
                                   SUNAllocator* allocator)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }
  SUNFunctionBegin(sunctx);

  SUNAssert(allocator, SUN_ERR_SUNCTX_CORRUPT);

  switch (type)
  {
    case SUNMEMTYPE_HOST:
      *allocator = sunctx->host_allocator;
      break;
    case SUNMEMTYPE_DEVICE:
      *allocator = sunctx->device_allocator;
      break;
    case SUNMEMTYPE_UVM:
      *allocator = sunctx->uvm_allocator;
      break;
    case SUNMEMTYPE_PINNED:
      *allocator = sunctx->pinned_allocator;
      break;
    default:
      return SUN_ERR_ARG_OUTOFRANGE;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNContext_PrintAllocatorStats(SUNContext sunctx, FILE* outfile,
                                          SUNOutputFormat fmt)
{
  if (!sunctx) { return SUN_ERR_SUNCTX_CORRUPT; }
  if (!outfile) { return SUN_ERR_ARG_CORRUPT; }

  SUNFunctionBegin(sunctx);

  SUNCheckCall(SUNAllocator_PrintStats(sunctx->host_allocator, outfile, fmt));

  if (sunctx->device_allocator)
  {
    SUNCheckCall(SUNAllocator_PrintStats(sunctx->host_allocator, outfile, fmt));
  }
  if (sunctx->uvm_allocator)
  {
    SUNCheckCall(SUNAllocator_PrintStats(sunctx->uvm_allocator, outfile, fmt));
  }
  if (sunctx->pinned_allocator)
  {
    SUNCheckCall(SUNAllocator_PrintStats(sunctx->pinned_allocator, outfile, fmt));
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNContext_Free(SUNContext* sunctx)
{
#ifdef SUNDIALS_ADIAK_ENABLED
  adiak_fini();
#endif

  if (!sunctx || !(*sunctx)) { return SUN_SUCCESS; }

#if defined(SUNDIALS_ENABLE_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
  /* Find out where we are printing to */
  FILE* fp                    = NULL;
  char* sunprofiler_print_env = getenv("SUNPROFILER_PRINT");
  fp                          = NULL;
  if (sunprofiler_print_env)
  {
    if (!strcmp(sunprofiler_print_env, "0")) { fp = NULL; }
    else if (!strcmp(sunprofiler_print_env, "1") ||
             !strcmp(sunprofiler_print_env, "TRUE") ||
             !strcmp(sunprofiler_print_env, "stdout"))
    {
      fp = stdout;
    }
    else { fp = fopen(sunprofiler_print_env, "a"); }
  }

  /* Enforce that the profiler is freed before finalizing,
     if it is not owned by the sunctx. */
  if ((*sunctx)->profiler)
  {
    if (fp) { SUNProfiler_Print((*sunctx)->profiler, fp); }
    if (fp) { fclose(fp); }
    if ((*sunctx)->own_profiler) { SUNProfiler_Free(&(*sunctx)->profiler); }
  }
#endif

  if ((*sunctx)->logger && (*sunctx)->own_logger)
  {
    SUNLogger_Destroy(&(*sunctx)->logger);
  }

  if ((*sunctx)->own_host_allocator)
  {
    SUNAllocator_Destroy(&(*sunctx)->host_allocator);
  }
  if ((*sunctx)->own_device_allocator)
  {
    SUNAllocator_Destroy(&(*sunctx)->device_allocator);
  }
  if ((*sunctx)->own_uvm_allocator)
  {
    SUNAllocator_Destroy(&(*sunctx)->uvm_allocator);
  }
  if ((*sunctx)->own_pinned_allocator)
  {
    SUNAllocator_Destroy(&(*sunctx)->pinned_allocator);
  }

  SUNContext_ClearErrHandlers(*sunctx);

#if defined(SUNDIALS_ENABLE_PYTHON)
  SUNContextFunctionTable_Destroy((*sunctx)->python);
#endif
  (*sunctx)->python = NULL;

  free(*sunctx);
  *sunctx = NULL;

  return SUN_SUCCESS;
}
