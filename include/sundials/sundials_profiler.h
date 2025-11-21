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

#ifndef _SUNDIALS_PROFILER_H
#define _SUNDIALS_PROFILER_H

#include <stdio.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>

#if defined(SUNDIALS_BUILD_WITH_PROFILING) && defined(SUNDIALS_CALIPER_ENABLED)
#include "caliper/cali.h"
#endif

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

SUNDIALS_EXPORT
SUNErrCode SUNProfiler_Create(SUNComm comm, const char* title, SUNProfiler* p);

SUNDIALS_EXPORT
SUNErrCode SUNProfiler_Free(SUNProfiler* p);

SUNDIALS_EXPORT
SUNErrCode SUNProfiler_Begin(SUNProfiler p, const char* name);

SUNDIALS_EXPORT
SUNErrCode SUNProfiler_End(SUNProfiler p, const char* name);

SUNDIALS_EXPORT
SUNErrCode SUNProfiler_GetTimerResolution(SUNProfiler p, double* resolution);

SUNDIALS_EXPORT
SUNErrCode SUNProfiler_GetElapsedTime(SUNProfiler p, const char* name,
                                      double* time);

SUNDIALS_EXPORT
SUNErrCode SUNProfiler_Print(SUNProfiler p, FILE* fp);

SUNDIALS_EXPORT
SUNErrCode SUNProfiler_Reset(SUNProfiler p);

#if defined(SUNDIALS_BUILD_WITH_PROFILING) && defined(SUNDIALS_CALIPER_ENABLED)

#define SUNDIALS_MARK_FUNCTION_BEGIN(profobj) CALI_MARK_FUNCTION_BEGIN

#define SUNDIALS_MARK_FUNCTION_END(profobj) CALI_MARK_FUNCTION_END

#define SUNDIALS_WRAP_STATEMENT(profobj, name, stmt) \
  CALI_WRAP_STATEMENT(name, stmt)

#define SUNDIALS_MARK_BEGIN(profobj, name) CALI_MARK_BEGIN(name)

#define SUNDIALS_MARK_END(profobj, name) CALI_MARK_END(name)

#elif defined(SUNDIALS_BUILD_WITH_PROFILING)

#define SUNDIALS_MARK_FUNCTION_BEGIN(profobj) \
  SUNProfiler_Begin(profobj, __func__)

#define SUNDIALS_MARK_FUNCTION_END(profobj) SUNProfiler_End(profobj, __func__)

#define SUNDIALS_WRAP_STATEMENT(profobj, name, stmt) \
  SUNProfiler_Begin(profobj, (name));                \
  stmt;                                              \
  SUNProfiler_End(profobj, (name));

#define SUNDIALS_MARK_BEGIN(profobj, name) SUNProfiler_Begin(profobj, (name))

#define SUNDIALS_MARK_END(profobj, name) SUNProfiler_End(profobj, (name))

#else

#define SUNDIALS_MARK_FUNCTION_BEGIN(profobj)

#define SUNDIALS_MARK_FUNCTION_END(profobj)

#define SUNDIALS_WRAP_STATEMENT(profobj, name, stmt)

#define SUNDIALS_MARK_BEGIN(profobj, name)

#define SUNDIALS_MARK_END(profobj, name)

#endif

#ifdef __cplusplus
}

#endif
#endif /* SUNDIALS_PROFILER_H */
