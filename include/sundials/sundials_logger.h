/* -----------------------------------------------------------------
 * Programmer: Cody J. Balos @ LLNL
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
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_LOGGER_H
#define _SUNDIALS_LOGGER_H

#include <stdio.h>

#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

enum SUNLogLevel
{
  SUN_LOGLEVEL_ALL     = -1,
  SUN_LOGLEVEL_NONE    = 0,
  SUN_LOGLEVEL_ERROR   = 1,
  SUN_LOGLEVEL_WARNING = 2,
  SUN_LOGLEVEL_INFO    = 3,
  SUN_LOGLEVEL_DEBUG   = 4
};

#ifndef SWIG
typedef enum SUNLogLevel SUNLogLevel;
#endif

SUNDIALS_EXPORT
SUNErrCode SUNLogger_Create(SUNComm comm, int output_rank, SUNLogger* logger);

SUNDIALS_EXPORT
SUNErrCode SUNLogger_CreateFromEnv(SUNComm comm, SUNLogger* logger);

SUNDIALS_EXPORT
SUNErrCode SUNLogger_SetErrorFilename(SUNLogger logger,
                                      const char* error_filename);

SUNDIALS_EXPORT
SUNErrCode SUNLogger_SetWarningFilename(SUNLogger logger,
                                        const char* warning_filename);

SUNDIALS_EXPORT
SUNErrCode SUNLogger_SetDebugFilename(SUNLogger logger,
                                      const char* debug_filename);

SUNDIALS_EXPORT
SUNErrCode SUNLogger_SetInfoFilename(SUNLogger logger, const char* info_filename);

SUNDIALS_EXPORT
SUNErrCode SUNLogger_QueueMsg(SUNLogger logger, SUNLogLevel lvl,
                              const char* scope, const char* label,
                              const char* msg_txt, ...);

SUNDIALS_EXPORT
SUNErrCode SUNLogger_Flush(SUNLogger logger, SUNLogLevel lvl);

SUNDIALS_EXPORT
SUNErrCode SUNLogger_GetOutputRank(SUNLogger logger, int* output_rank);

SUNDIALS_EXPORT
SUNErrCode SUNLogger_Destroy(SUNLogger* logger);

#ifdef __cplusplus /* wrapper to enable C++ usage */
}
#endif
#endif /* SUNDIALS_LOGGER_H_ */
