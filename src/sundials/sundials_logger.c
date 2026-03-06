/* -----------------------------------------------------------------
 * Programmer: Cody J. Balos @ LLNL
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
 * -----------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/sundials_config.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_logger.h>

#include <sundials/priv/sundials_errors_impl.h>
#include "sundials_logger_impl.h"

#include "sundials/sundials_errors.h"
#include "sundials/sundials_types.h"

#if SUNDIALS_MPI_ENABLED
#include <mpi.h>
#endif

#include "sundials_hashmap_impl.h"
#include "sundials_macros.h"
#include "sundials_utils.h"

/* default number of files that we allocate space for */
#define SUN_DEFAULT_LOGFILE_HANDLES_ 8

void sunCreateLogMessage(SUNLogLevel lvl, int rank, const char* scope,
                         const char* label, const char* txt, va_list args,
                         char** log_msg)
{
  const char* prefix;
  char* formatted_txt;
  int msg_length;

  prefix        = NULL;
  formatted_txt = NULL;
  msg_length    = 0;
  *log_msg      = NULL;

  msg_length = sunvasnprintf(&formatted_txt, txt, args);
  if (msg_length < 0)
  {
    char* fileAndLine = sunCombineFileAndLine(__LINE__ + 1, __FILE__);
    fprintf(stderr, "[ERROR][rank %d][%s][%s] %s\n", rank, fileAndLine,
            __func__, "FATAL LOGGER ERROR: message size too large");
    free(fileAndLine);
  }

  if (lvl == SUN_LOGLEVEL_DEBUG) { prefix = "DEBUG"; }
  else if (lvl == SUN_LOGLEVEL_WARNING) { prefix = "WARNING"; }
  else if (lvl == SUN_LOGLEVEL_INFO) { prefix = "INFO"; }
  else if (lvl == SUN_LOGLEVEL_ERROR) { prefix = "ERROR"; }

  msg_length = snprintf(NULL, 0, "[%s][rank %d][%s][%s] %s\n", prefix, rank,
                        scope, label, formatted_txt);
  *log_msg   = (char*)malloc(msg_length + 1);
  snprintf(*log_msg, msg_length + 1, "[%s][rank %d][%s][%s] %s\n", prefix, rank,
           scope, label, formatted_txt);
  free(formatted_txt);
}

#if SUNDIALS_LOGGING_LEVEL > 0
static FILE* sunOpenLogFile(const char* fname, const char* mode)
{
  FILE* fp = NULL;

  if (fname)
  {
    if (!strcmp(fname, "stdout")) { fp = stdout; }
    else if (!strcmp(fname, "stderr")) { fp = stderr; }
    else { fp = fopen(fname, mode); }
  }

  return fp;
}
#endif

static void sunCloseLogFile(void* fp)
{
  if (fp && fp != stdout && fp != stderr) { fclose((FILE*)fp); }
}

static sunbooleantype sunLoggerIsOutputRank(SUNDIALS_MAYBE_UNUSED SUNLogger logger,
                                            int* rank_ref)
{
  sunbooleantype retval;

#if SUNDIALS_MPI_ENABLED
  int rank = 0;

  if (logger->comm != SUN_COMM_NULL)
  {
    MPI_Comm_rank(logger->comm, &rank);

    if (logger->output_rank < 0)
    {
      if (rank_ref) { *rank_ref = rank; }
      retval = SUNTRUE; /* output all ranks */
    }
    else
    {
      if (rank_ref) { *rank_ref = rank; }
      retval = logger->output_rank == rank;
    }
  }
  else { retval = SUNTRUE; /* output all ranks */ }
#else
  if (rank_ref) { *rank_ref = 0; }
  retval = SUNTRUE;
#endif

  return retval;
}

static SUNErrCode sunLoggerFreeKeyValue(SUNHashMapKeyValue* kv_ptr)
{
  if (!kv_ptr || !(*kv_ptr)) { return SUN_SUCCESS; }
  sunCloseLogFile((*kv_ptr)->value);
  free((*kv_ptr)->key);
  free(*kv_ptr);
  return SUN_SUCCESS;
}

SUNErrCode SUNLogger_Create(SUNComm comm, int output_rank, SUNLogger* logger_ptr)
{
  SUNLogger logger = NULL;

  *logger_ptr = logger = (SUNLogger)malloc(sizeof(struct SUNLogger_));
  if (logger == NULL) { return SUN_ERR_MALLOC_FAIL; }

  /* Attach the comm, duplicating it if MPI is used. */
#if SUNDIALS_MPI_ENABLED
  logger->comm = SUN_COMM_NULL;
  if (comm != SUN_COMM_NULL) { MPI_Comm_dup(comm, &logger->comm); }
#else
  logger->comm = SUN_COMM_NULL;
  if (comm != SUN_COMM_NULL)
  {
    free(logger);
    return SUN_ERR_ARG_CORRUPT;
  }
#endif
  logger->output_rank = output_rank;
  logger->content     = NULL;

  /* use default routines */
  logger->queuemsg = NULL;
  logger->flush    = NULL;
  logger->destroy  = NULL;

  /* set the output file handles */
  logger->filenames  = NULL;
  logger->error_fp   = stderr;
  logger->warning_fp = stdout;
  logger->debug_fp   = NULL;
  logger->info_fp    = NULL;
  if (sunLoggerIsOutputRank(logger, NULL))
  {
    /* We store the FILE* in a hash map so that we can ensure
       that we do not open a file twice if the same file is used
       for multiple output levels */
    SUNHashMap_New(SUN_DEFAULT_LOGFILE_HANDLES_, sunLoggerFreeKeyValue,
                   &logger->filenames);
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNLogger_CreateFromEnv(SUNComm comm, SUNLogger* logger_out)
{
  SUNErrCode err   = SUN_SUCCESS;
  SUNLogger logger = NULL;

  const char* output_rank_env   = getenv("SUNLOGGER_OUTPUT_RANK");
  int output_rank               = (output_rank_env) ? atoi(output_rank_env) : 0;
  const char* error_fname_env   = getenv("SUNLOGGER_ERROR_FILENAME");
  const char* warning_fname_env = getenv("SUNLOGGER_WARNING_FILENAME");
  const char* info_fname_env    = getenv("SUNLOGGER_INFO_FILENAME");
  const char* debug_fname_env   = getenv("SUNLOGGER_DEBUG_FILENAME");

  if (SUNLogger_Create(comm, output_rank, &logger))
  {
    err = SUN_ERR_CORRUPT;
    return err;
  }

  do {
    err = SUNLogger_SetErrorFilename(logger, error_fname_env);
    if (err) { break; }
    err = SUNLogger_SetWarningFilename(logger, warning_fname_env);
    if (err) { break; }
    err = SUNLogger_SetDebugFilename(logger, debug_fname_env);
    if (err) { break; }
    err = SUNLogger_SetInfoFilename(logger, info_fname_env);
  }
  while (0);

  if (err) { SUNLogger_Destroy(&logger); }
  else { *logger_out = logger; }

  return err;
}

SUNErrCode SUNLogger_SetErrorFilename(SUNLogger logger, const char* error_filename)
{
  if (!logger) { return SUN_ERR_ARG_CORRUPT; }

  if (!sunLoggerIsOutputRank(logger, NULL)) { return SUN_SUCCESS; }

  if (error_filename && strcmp(error_filename, ""))
  {
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_ERROR
    FILE* fp = NULL;
    if (!SUNHashMap_GetValue(logger->filenames, error_filename, (void*)&fp))
    {
      logger->error_fp = fp;
    }
    else
    {
      logger->error_fp = sunOpenLogFile(error_filename, "w+");
      if (logger->error_fp)
      {
        SUNHashMap_Insert(logger->filenames, error_filename,
                          (void*)logger->error_fp);
      }
      else { return SUN_ERR_FILE_OPEN; }
    }
#endif
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNLogger_SetWarningFilename(SUNLogger logger,
                                        const char* warning_filename)
{
  if (!logger) { return SUN_ERR_ARG_CORRUPT; }

  if (!sunLoggerIsOutputRank(logger, NULL)) { return SUN_SUCCESS; }

  if (warning_filename && strcmp(warning_filename, ""))
  {
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_WARNING
    FILE* fp = NULL;
    if (!SUNHashMap_GetValue(logger->filenames, warning_filename, (void*)&fp))
    {
      logger->warning_fp = fp;
    }
    else
    {
      logger->warning_fp = sunOpenLogFile(warning_filename, "w+");
      if (logger->warning_fp)
      {
        SUNHashMap_Insert(logger->filenames, warning_filename,
                          (void*)logger->warning_fp);
      }
      else { return SUN_ERR_FILE_OPEN; }
    }
#endif
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNLogger_SetInfoFilename(SUNLogger logger, const char* info_filename)
{
  if (!logger) { return SUN_ERR_ARG_CORRUPT; }

  if (!sunLoggerIsOutputRank(logger, NULL)) { return SUN_SUCCESS; }

  if (info_filename && strcmp(info_filename, ""))
  {
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
    FILE* fp = NULL;
    if (!SUNHashMap_GetValue(logger->filenames, info_filename, (void*)&fp))
    {
      logger->info_fp = fp;
    }
    else
    {
      logger->info_fp = sunOpenLogFile(info_filename, "w+");
      if (logger->info_fp)
      {
        SUNHashMap_Insert(logger->filenames, info_filename,
                          (void*)logger->info_fp);
      }
      else { return SUN_ERR_FILE_OPEN; }
    }
#endif
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNLogger_SetDebugFilename(SUNLogger logger, const char* debug_filename)
{
  if (!logger) { return SUN_ERR_ARG_CORRUPT; }

  if (!sunLoggerIsOutputRank(logger, NULL)) { return SUN_SUCCESS; }

  if (debug_filename && strcmp(debug_filename, ""))
  {
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_DEBUG
    FILE* fp = NULL;
    if (!SUNHashMap_GetValue(logger->filenames, debug_filename, (void*)&fp))
    {
      logger->debug_fp = fp;
    }
    else
    {
      logger->debug_fp = sunOpenLogFile(debug_filename, "w+");
      if (logger->debug_fp)
      {
        SUNHashMap_Insert(logger->filenames, debug_filename,
                          (void*)logger->debug_fp);
      }
      else { return SUN_ERR_FILE_OPEN; }
    }
#endif
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNLogger_QueueMsg(SUNLogger logger, SUNLogLevel lvl,
                              const char* scope, const char* label,
                              const char* msg_txt, ...)
{
  SUNErrCode retval = SUN_SUCCESS;

#if SUNDIALS_LOGGING_LEVEL > 0
  {
    if (!logger)
    {
      retval = SUN_ERR_ARG_CORRUPT;
      return retval;
    }

    if (logger->queuemsg)
    {
      va_list args;
      va_start(args, msg_txt);
      retval = logger->queuemsg(logger, lvl, scope, label, msg_txt, args);
      va_end(args);
    }
    else
    {
      /* Default implementation */
      int rank = 0;
      if (sunLoggerIsOutputRank(logger, &rank))
      {
        char* log_msg = NULL;
        va_list args;
        va_start(args, msg_txt);
        sunCreateLogMessage(lvl, rank, scope, label, msg_txt, args, &log_msg);
        va_end(args);

        switch (lvl)
        {
        case (SUN_LOGLEVEL_DEBUG):
          if (logger->debug_fp) { fprintf(logger->debug_fp, "%s", log_msg); }
          break;
        case (SUN_LOGLEVEL_WARNING):
          if (logger->warning_fp)
          {
            fprintf(logger->warning_fp, "%s", log_msg);
          }
          break;
        case (SUN_LOGLEVEL_INFO):
          if (logger->info_fp) { fprintf(logger->info_fp, "%s", log_msg); }
          break;
        case (SUN_LOGLEVEL_ERROR):
          if (logger->error_fp) { fprintf(logger->error_fp, "%s", log_msg); }
          break;
        default: retval = SUN_ERR_UNREACHABLE;
        }

        free(log_msg);
      }
    }
  }
#else
  /* silence warnings when all logging is disabled */
  ((void)logger);
  ((void)lvl);
  ((void)scope);
  ((void)label);
  ((void)msg_txt);
#endif

  return retval;
}

SUNErrCode SUNLogger_Flush(SUNLogger logger, SUNLogLevel lvl)
{
  SUNErrCode retval = SUN_SUCCESS;

  if (!logger)
  {
    retval = SUN_ERR_ARG_CORRUPT;
    return retval;
  }

#if SUNDIALS_LOGGING_LEVEL > 0
  if (logger->flush) { retval = logger->flush(logger, lvl); }
  else
  {
    /* Default implementation */
    if (sunLoggerIsOutputRank(logger, NULL))
    {
      switch (lvl)
      {
      case (SUN_LOGLEVEL_DEBUG):
        if (logger->debug_fp) { fflush(logger->debug_fp); }
        break;
      case (SUN_LOGLEVEL_WARNING):
        if (logger->warning_fp) { fflush(logger->warning_fp); }
        break;
      case (SUN_LOGLEVEL_INFO):
        if (logger->info_fp) { fflush(logger->info_fp); }
        break;
      case (SUN_LOGLEVEL_ERROR):
        if (logger->error_fp) { fflush(logger->error_fp); }
        break;
      case (SUN_LOGLEVEL_ALL):
        if (logger->debug_fp) { fflush(logger->debug_fp); }
        if (logger->warning_fp) { fflush(logger->warning_fp); }
        if (logger->info_fp) { fflush(logger->info_fp); }
        if (logger->error_fp) { fflush(logger->error_fp); }
        break;
      default: retval = SUN_ERR_UNREACHABLE;
      }
    }
  }
#else
  /* silence warnings when all logging is disabled */
  ((void)lvl);
#endif

  return retval;
}

SUNErrCode SUNLogger_GetOutputRank(SUNLogger logger, int* output_rank)
{
  if (!logger) { return SUN_ERR_ARG_CORRUPT; }
  *output_rank = logger->output_rank;
  return SUN_SUCCESS;
}

SUNErrCode SUNLogger_Destroy(SUNLogger* logger_ptr)
{
  int retval       = 0;
  SUNLogger logger = NULL;

  if (!logger_ptr) { return SUN_SUCCESS; }

  logger = *logger_ptr;

  if (logger && logger->destroy) { retval = logger->destroy(logger_ptr); }
  else if (logger)
  {
    /* Default implementation */

    if (sunLoggerIsOutputRank(logger, NULL))
    {
      SUNHashMap_Destroy(&logger->filenames);
    }

#if SUNDIALS_MPI_ENABLED
    if (logger->comm != SUN_COMM_NULL) { MPI_Comm_free(&logger->comm); }
#endif

    free(logger);
    logger = NULL;
  }

  return retval;
}
