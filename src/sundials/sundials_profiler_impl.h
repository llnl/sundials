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

#ifndef SUNDIALS_PROFILER_IMPL_H_
#define SUNDIALS_PROFILER_IMPL_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_profiler.h>
#include <sundials/sundials_types.h>

#if SUNDIALS_MPI_ENABLED
#include <mpi.h>
#endif

#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
#include <stddef.h>
#include <time.h>
#include <unistd.h>
#elif defined(WIN32) || defined(_WIN32)
#include <windows.h>
#else
#error SUNProfiler needs POSIX or Windows timers
#endif

#include "sundials_debug.h"
#include "sundials_hashmap_impl.h"
#include "sundials_macros.h"

#define SUNDIALS_ROOT_TIMER ((const char*)"From profiler epoch")

#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
typedef struct timespec sunTimespec;
#else
typedef struct _sunTimespec
{
  long int tv_sec;
  long int tv_nsec;
} sunTimespec;
#endif

/*
  sunTimerStruct.
  A private structure holding timing information.
 */

struct _sunTimerStruct
{
  sunTimespec* tic;
  sunTimespec* toc;
  double average;
  double maximum;
  double elapsed;
  long count;
};

typedef struct _sunTimerStruct sunTimerStruct;

struct SUNProfiler_
{
  SUNComm comm;
  char* title;
  SUNHashMap map;
  sunTimerStruct* overhead;
  double sundials_time;
};

#endif
