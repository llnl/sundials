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
 * SUNDIALS logging class implementation.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_LOGGER_IMPL_H
#define _SUNDIALS_LOGGER_IMPL_H

#include <stdarg.h>
#include <sundials/priv/sundials_logger_macros.h>
#include <sundials/sundials_logger.h>
#include <sundials/sundials_types.h>

#include "sundials_hashmap_impl.h"
#include "sundials_utils.h"

#define SUNDIALS_LOGGING_ERROR   1
#define SUNDIALS_LOGGING_WARNING 2
#define SUNDIALS_LOGGING_INFO    3
#define SUNDIALS_LOGGING_DEBUG   4
#if SUNDIALS_LOGGING_LEVEL > SUNDIALS_LOGGING_DEBUG
#define SUNDIALS_LOGGING_EXTRA_DEBUG
#endif

/*
  In the variadic logging macros below, the message text (msg_txt) is not
  explicitly included as a macro parameter and instead inserted by
  __VA_ARGS__. This allows us to omit adding an empty string for the optional
  arguments when the message does not include format specifiers e.g.,

  SUNLogInfo(logger, "label", "message");

  instead of

  SUNLogInfo(logger, "label", "message", "");

  Without this workaround, an orphaned comma is placed end of the argument list
  if the empty string is not included. Note the C23 standard adds __VA_OPT__
  which will allow for explicitly including msg_txt while removing the trailing
  comma when no additional arguments are needed.
*/

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
#define SUNLogInfo(logger, label, /* msg_txt, */...)             \
  SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_INFO, __func__, label, \
                     /* msg_txt, */ __VA_ARGS__)
#define SUNLogInfoIf(condition, logger, label, /* msg_txt, */...)    \
  do {                                                               \
    if ((condition))                                                 \
    {                                                                \
      SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_INFO, __func__, label, \
                         /* msg_txt, */ __VA_ARGS__);                \
    }                                                                \
  }                                                                  \
  while (0)
#else
#define SUNLogInfo(logger, label, /* msg_txt, */...)
#define SUNLogInfoIf(condition, logger, label, /* msg_txt, */...)
#endif

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_DEBUG
#define SUNLogDebug(logger, label, /* msg_txt, */...)             \
  SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_DEBUG, __func__, label, \
                     /* msg_txt, */ __VA_ARGS__)
#define SUNLogDebugIf(condition, logger, label, /* msg_txt, */...)    \
  do {                                                                \
    if ((condition))                                                  \
    {                                                                 \
      SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_DEBUG, __func__, label, \
                         /* msg_txt, */ __VA_ARGS__);                 \
    }                                                                 \
  }                                                                   \
  while (0)
#else
#define SUNLogDebug(logger, label, /* msg_txt, */...)
#define SUNLogDebugIf(condition, logger, label, /* msg_txt, */...)
#endif

#ifdef SUNDIALS_LOGGING_EXTRA_DEBUG
#define SUNLogExtraDebug(logger, label, /* msg_txt, */...)        \
  SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_DEBUG, __func__, label, \
                     /* msg_txt, */ __VA_ARGS__)
#define SUNLogExtraDebugIf(condition, logger, label, /* msg_txt, */...) \
  do {                                                                  \
    if ((condition))                                                    \
    {                                                                   \
      SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_DEBUG, __func__, label,   \
                         /* msg_txt, */ __VA_ARGS__);                   \
    }                                                                   \
  }                                                                     \
  while (0)
#define SUNLogExtraDebugVec(logger, label, vec, /*msg_txt, */...)   \
  do {                                                              \
    SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_DEBUG, __func__, label, \
                       /* msg_txt, */ __VA_ARGS__);                 \
    N_VPrintFile(vec, logger->debug_fp);                            \
  }                                                                 \
  while (0)
#define SUNLogExtraDebugVecIf(condition, logger, label, vec, /* msg_txt, */...) \
  do {                                                                          \
    if ((condition))                                                            \
    {                                                                           \
      SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_DEBUG, __func__, label,           \
                         /* msg_txt, */ __VA_ARGS__);                           \
      N_VPrintFile(vec, logger->debug_fp);                                      \
    }                                                                           \
  }                                                                             \
  while (0)
#define SUNLogExtraDebugVecArray(logger, label, nvecs, vecs, msg_txt)          \
  do {                                                                         \
    for (int vi = 0; vi < (nvecs); ++vi)                                       \
    {                                                                          \
      SUNLogger_QueueMsg(logger, SUN_LOGLEVEL_DEBUG, __func__, label, msg_txt, \
                         vi);                                                  \
      N_VPrintFile(vecs[vi], logger->debug_fp);                                \
    }                                                                          \
  }                                                                            \
  while (0)
#else
#define SUNLogExtraDebug(logger, label, /* msg_txt, */...)
#define SUNLogExtraDebugIf(condition, logger, label, /* msg_txt, */...)
#define SUNLogExtraDebugVec(logger, label, vec, /* msg_txt, */...)
#define SUNLogExtraDebugVecIf(condition, logger, label, vec, /* msg_txt, */...)
#define SUNLogExtraDebugVecArray(logger, label, nvecs, vecs, msg_txt)
#endif

struct SUNLogger_
{
  /* MPI information */
  SUNComm comm;
  int output_rank;

  /* Output files */
  FILE* debug_fp;
  FILE* warning_fp;
  FILE* info_fp;
  FILE* error_fp;

  /* Hashmap used to store filename, FILE* pairs */
  SUNHashMap filenames;

  /* Slic-style format string */
  const char* format;

  /* Content for custom implementations */
  void* content;

  /* Overridable operations */
  SUNErrCode (*queuemsg)(SUNLogger logger, SUNLogLevel lvl, const char* scope,
                         const char* label, const char* msg_txt, va_list args);
  SUNErrCode (*flush)(SUNLogger logger, SUNLogLevel lvl);
  SUNErrCode (*destroy)(SUNLogger* logger);
};

/*
  This function creates a log message string in the correct format.
  It allocates the log_msg parameter, which must be freed by the caller.
  The format of the log message is:

    [ERROR][rank <rank>][<scope>][<label>] <formatted txt>

  :param lvl: the logging level (ERROR, WARNING, INFO, DEBUG)
  :param rank: the MPI rank of the caller
  :param scope: the scope part of the log message (see above)
  :param label: the label part of the log message (see above)
  :param txt: descriptive text for the log message in the form of a format string
  :param args: format string substitutions
  :param log_msg: on output this is an allocated string containing the log message

  :return: void
*/
void sunCreateLogMessage(SUNLogLevel lvl, int rank, const char* scope,
                         const char* label, const char* txt, va_list args,
                         char** log_msg);

#endif /* _SUNDIALS_LOGGER_IMPL_H */
