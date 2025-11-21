/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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

#include "sundials_adiak_metadata.h"
#include "sundials_macros.h"

SUNErrCode SUNContext_Create(SUNComm comm, SUNContext* sunctx_out)
{
  SUNErrCode err       = SUN_SUCCESS;
  SUNProfiler profiler = NULL;
  SUNLogger logger     = NULL;
  SUNContext sunctx    = NULL;
  SUNErrHandler eh     = NULL;

  *sunctx_out = NULL;
  sunctx      = (SUNContext)malloc(sizeof(struct SUNContext_));

  /* SUNContext_Create cannot assert or log since the SUNContext is not yet
   * created */
  if (!sunctx) { return SUN_ERR_MALLOC_FAIL; }

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

#if defined(SUNDIALS_BUILD_WITH_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
    err = SUNProfiler_Create(comm, "SUNContext Default", &profiler);
    SUNCheckCallNoRet(err);
    if (err) { break; }
#endif

    err = SUNErrHandler_Create(SUNLogErrHandlerFn, NULL, &eh);
    SUNCheckCallNoRet(err);
    if (err) { break; }

    sunctx->python       = NULL;
    sunctx->logger       = logger;
    sunctx->own_logger   = logger != NULL;
    sunctx->profiler     = profiler;
    sunctx->own_profiler = profiler != NULL;
    sunctx->last_err     = SUN_SUCCESS;
    sunctx->err_handler  = eh;
    sunctx->comm         = comm;
  }
  while (0);

  if (err)
  {
#if defined(SUNDIALS_BUILD_WITH_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
    SUNCheckCallNoRet(SUNProfiler_Free(&profiler));
#endif
    SUNCheckCallNoRet(SUNLogger_Destroy(&logger));
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

#ifdef SUNDIALS_BUILD_WITH_PROFILING
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

#ifdef SUNDIALS_BUILD_WITH_PROFILING
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

SUNErrCode SUNContext_Free(SUNContext* sunctx)
{
#ifdef SUNDIALS_ADIAK_ENABLED
  adiak_fini();
#endif

  if (!sunctx || !(*sunctx)) { return SUN_SUCCESS; }

#if defined(SUNDIALS_BUILD_WITH_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
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

  SUNContext_ClearErrHandlers(*sunctx);

  free((*sunctx)->python);
  (*sunctx)->python = NULL;

  free(*sunctx);
  *sunctx = NULL;

  return SUN_SUCCESS;
}
