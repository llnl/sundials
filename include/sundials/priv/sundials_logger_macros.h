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
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_LOGGER_MACROS_H
#define _SUNDIALS_LOGGER_MACROS_H

#include <stdarg.h>

#include <sundials/sundials_logger.h>
#include <sundials/sundials_types.h>

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

#endif /* _SUNDIALS_LOGGER_MACROS_H */
