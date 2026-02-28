/* -----------------------------------------------------------------
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
 * -----------------------------------------------------------------*/

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_logger.h>
#include <sundials/sundials_types.h>

#include "sundials_logger_impl.h"
#include "sundials_macros.h"
#include "sundials_utils.h"

SUNErrCode SUNErrHandler_Create(SUNErrHandlerFn eh_fn, void* eh_data,
                                SUNErrHandler* eh_out)
{
  SUNErrHandler eh = NULL;

  eh = (SUNErrHandler)malloc(sizeof(struct SUNErrHandler_));
  if (!eh) { return SUN_ERR_MALLOC_FAIL; }

  eh->previous = NULL;
  eh->call     = eh_fn;
  eh->data     = eh_data;

  *eh_out = eh;
  return SUN_SUCCESS;
}

void SUNErrHandler_Destroy(SUNErrHandler* eh)
{
  if (!eh || !(*eh)) { return; }
  free(*eh);
  *eh = NULL;
}

const char* SUNGetErrMsg(SUNErrCode code)
{
#define SUN_EXPAND_TO_CASES(name, description) \
  case name: return description; break;

  switch (code)
  {
    SUN_ERR_CODE_LIST(SUN_EXPAND_TO_CASES)
  default: return "unknown error";
  }

  return NULL;
}

void SUNLogErrHandlerFn(int line, const char* func, const char* file,
                        const char* msg, SUNErrCode err_code,
                        SUNDIALS_MAYBE_UNUSED void* err_user_data,
                        SUNContext sunctx)
{
  char* file_and_line = sunCombineFileAndLine(line, file);
  if (msg == NULL) { msg = SUNGetErrMsg(err_code); }
  SUNLogger_QueueMsg(sunctx->logger, SUN_LOGLEVEL_ERROR, file_and_line, func,
                     msg);
  free(file_and_line);
}

void SUNAbortErrHandlerFn(int line, const char* func, const char* file,
                          SUNDIALS_MAYBE_UNUSED const char* msg,
                          SUNDIALS_MAYBE_UNUSED SUNErrCode err_code,
                          SUNDIALS_MAYBE_UNUSED void* err_user_data,
                          SUNContext sunctx)
{
  /* Flush all buffered logging messages now before we abort */
  SUNLogger_Flush(sunctx->logger, SUN_LOGLEVEL_ALL);

  char* file_and_line = sunCombineFileAndLine(line, file);
  if (msg == NULL) { msg = SUNGetErrMsg(err_code); }
  SUNLogger_QueueMsg(sunctx->logger, SUN_LOGLEVEL_ERROR, file_and_line, func,
                     msg);
  free(file_and_line);
  /* It is convenient to have the exit message point to the message line,
     so we add 1 to the line number. As such, do not separate the following lines! */
  file_and_line = sunCombineFileAndLine(__LINE__ + 1, __FILE__);
  SUNLogger_QueueMsg(sunctx->logger, SUN_LOGLEVEL_ERROR, file_and_line, __func__,
                     "SUNAbortErrHandler: Calling abort now, use a different "
                     "error handler to avoid program termination.\n");
  free(file_and_line);
  abort();
}

void SUNGlobalFallbackErrHandler(int line, const char* func, const char* file,
                                 const char* msgfmt, SUNErrCode err_code, ...)
{
  va_list ap;
  char* log_msg       = NULL;
  char* file_and_line = NULL;

  file_and_line = sunCombineFileAndLine(__LINE__, __FILE__);
  va_start(ap, err_code);
  sunCreateLogMessage(SUN_LOGLEVEL_ERROR, 0, file_and_line,
                      __func__, "The SUNDIALS SUNContext was corrupt or NULL when an error occurred. As such, error messages have been printed to stderr.",
                      ap, &log_msg);
  va_end(ap);
  fprintf(stderr, "%s", log_msg);
  free(log_msg);
  free(file_and_line);

  file_and_line = sunCombineFileAndLine(line, file);
  if (msgfmt == NULL) { msgfmt = SUNGetErrMsg(err_code); }
  va_start(ap, err_code);
  sunCreateLogMessage(SUN_LOGLEVEL_ERROR, 0, file_and_line, func, msgfmt, ap,
                      &log_msg);
  va_end(ap);
  fprintf(stderr, "%s", log_msg);
  free(log_msg);
  free(file_and_line);
}
