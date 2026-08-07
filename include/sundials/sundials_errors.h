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
 * -----------------------------------------------------------------
 * Contains all error handling interfaces for SUNDIALS that do
 * not depend on MPI.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_ERRORS_H
#define _SUNDIALS_ERRORS_H

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_types.h>

/* ----------------------------------------------------------------------------
 * Error code definitions
 * ---------------------------------------------------------------------------*/

/* SUN_ERR_CODE_LIST is an X macro that can be expanded in various ways */
#define SUN_ERR_CODE_LIST(ENTRY)                                               \
  ENTRY(SUN_ERR_ARG_CORRUPT, "argument provided is NULL or corrupted")         \
  ENTRY(SUN_ERR_ARG_INCOMPATIBLE, "argument provided is not compatible")       \
  ENTRY(SUN_ERR_ARG_OUTOFRANGE, "argument is out of the valid range")          \
  ENTRY(SUN_ERR_ARG_WRONGTYPE, "argument provided is not the right type")      \
  ENTRY(SUN_ERR_ARG_DIMSMISMATCH, "argument dimensions do not agree")          \
                                                                               \
  ENTRY(SUN_ERR_GENERIC, "an error occurred")                                  \
  ENTRY(SUN_ERR_CORRUPT, "value is NULL or corrupt")                           \
  ENTRY(SUN_ERR_OUTOFRANGE, "Value is out of the expected range")              \
  ENTRY(SUN_ERR_FILE_OPEN, "Unable to open file")                              \
  ENTRY(SUN_ERR_OP_FAIL, "an operation failed")                                \
  ENTRY(SUN_ERR_MEM_FAIL, "a memory operation failed")                         \
  ENTRY(SUN_ERR_MALLOC_FAIL, "malloc returned NULL")                           \
  ENTRY(SUN_ERR_EXT_FAIL, "a failure occurred in an external library")         \
  ENTRY(SUN_ERR_DESTROY_FAIL, "a destroy function returned an error")          \
  ENTRY(SUN_ERR_NOT_IMPLEMENTED,                                               \
        "operation is not implemented: function pointer is NULL")              \
  ENTRY(SUN_ERR_USER_FCN_FAIL, "the user provided callback function failed")   \
                                                                               \
  ENTRY(SUN_ERR_DATANODE_NODENOTFOUND, "the data node could not be found")     \
                                                                               \
  ENTRY(SUN_ERR_PROFILER_MAPFULL,                                              \
        "the number of profiler entries exceeded SUNPROFILER_MAX_ENTRIES")     \
  ENTRY(SUN_ERR_PROFILER_MAPGET, "unknown error getting SUNProfiler timer")    \
  ENTRY(SUN_ERR_PROFILER_MAPINSERT,                                            \
        "unknown error inserting SUNProfiler timer")                           \
  ENTRY(SUN_ERR_PROFILER_MAPKEYNOTFOUND, "timer was not found in SUNProfiler") \
  ENTRY(SUN_ERR_PROFILER_MAPSORT, "error sorting SUNProfiler map")             \
                                                                               \
  ENTRY(SUN_ERR_ADJOINT_STEPPERFAILED,                                         \
        "SUNStepper stopped without successfully reaching the requested "      \
        "output time when solving the adjoint system")                         \
  ENTRY(SUN_ERR_ADJOINT_STEPPERINVALIDSTOP,                                    \
        "SUNStepper stopped with a flag not supported by the "                 \
        "adjoint solver")                                                      \
                                                                               \
  ENTRY(SUN_ERR_CHECKPOINT_NOT_FOUND,                                          \
        "the requested checkpoint was not found")                              \
  ENTRY(SUN_ERR_CHECKPOINT_MISMATCH, "the expected time for the checkpoint "   \
                                     "and the stored time do not match")       \
                                                                               \
  ENTRY(SUN_ERR_SUNCTX_CORRUPT, "SUNContext is NULL or corrupt")               \
                                                                               \
  ENTRY(SUN_ERR_MPI_FAIL,                                                      \
        "an MPI call returned something other than MPI_SUCCESS")               \
                                                                               \
  ENTRY(SUN_ERR_UNREACHABLE,                                                   \
        "Reached code that should be unreachable: open an issue at: "          \
        "https://github.com/LLNL/sundials")                                    \
  ENTRY(SUN_ERR_UNKNOWN, "Unknown error occurred: open an issue at "           \
                         "https://github.com/LLNL/sundials")

/* Expand SUN_ERR_CODE_LIST to enum */
#define SUN_EXPAND_TO_ENUM(name, description) name,

/* SUNErrorCode range is [-10000, -1000] to avoid conflicts with package error
   codes, and old/deprecated codes for matrix and (non)linear solvers. */

/* clang-format off */
enum SUNErrCode_
{
  SUN_ERR_MINIMUM                                       = -10000,
  SUN_ERR_CODE_LIST(SUN_EXPAND_TO_ENUM)
  SUN_ERR_MAXIMUM                                       = -1000,
  SUN_SUCCESS                                           = 0
};

/* clang-format on */

/* ----------------------------------------------------------------------------
 * Error handler definitions
 * ---------------------------------------------------------------------------*/

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

SUNDIALS_EXPORT
void SUNLogErrHandlerFn(int line, const char* func, const char* file,
                        const char* msg, SUNErrCode err_code,
                        void* err_user_data, SUNContext sunctx);

SUNDIALS_EXPORT
void SUNAbortErrHandlerFn(int line, const char* func, const char* file,
                          const char* msg, SUNErrCode err_code,
                          void* err_user_data, SUNContext sunctx);

/* ----------------------------------------------------------------------------
 * Error functions
 * ---------------------------------------------------------------------------*/

/* Turn error code into error message */
SUNDIALS_EXPORT
const char* SUNGetErrMsg(SUNErrCode code);

#ifdef __cplusplus /* wrapper to enable C++ usage */
} /* extern "C" */
#endif
#endif /* _SUNDIALS_ERRORS_H */
