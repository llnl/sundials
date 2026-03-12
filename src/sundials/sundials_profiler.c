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

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_profiler.h>
#include <sundials/sundials_types.h>

#include "sundials_debug.h"
#include "sundials_hashmap_impl.h"
#include "sundials_macros.h"
#include "sundials_profiler_impl.h"

/* Private functions */
#if SUNDIALS_MPI_ENABLED
static SUNErrCode sunCollectTimers(SUNProfiler p);
#endif
static void sunPrintTimer(SUNHashMapKeyValue kv, FILE* fp, void* pvoid);
static int sunCompareTimes(const void* l, const void* r);
static int sunclock_gettime_monotonic(sunTimespec* tp);

static sunTimerStruct* sunTimerStructNew(void)
{
  sunTimerStruct* ts = (sunTimerStruct*)malloc(sizeof(sunTimerStruct));
  ts->tic            = (sunTimespec*)malloc(sizeof(sunTimespec));
  ts->toc            = (sunTimespec*)malloc(sizeof(sunTimespec));
  ts->tic->tv_sec    = 0;
  ts->tic->tv_nsec   = 0;
  ts->elapsed        = 0.0;
  ts->average        = 0.0;
  ts->maximum        = 0.0;
  ts->count          = 0;
  return ts;
}

static void sunTimerStructFree(void* TS)
{
  sunTimerStruct* ts = (sunTimerStruct*)TS;
  if (ts)
  {
    if (ts->tic) { free(ts->tic); }
    if (ts->toc) { free(ts->toc); }
    free(ts);
  }
}

static SUNErrCode sunProfilerDestroyKeyValue(SUNHashMapKeyValue* kv_ptr)
{
  if (!kv_ptr || !(*kv_ptr)) { return SUN_SUCCESS; }
  sunTimerStructFree((*kv_ptr)->value);
  free((*kv_ptr)->key);
  free(*kv_ptr);
  return SUN_SUCCESS;
}

static void sunStartTiming(sunTimerStruct* entry)
{
  sunclock_gettime_monotonic(entry->tic);
}

static void sunStopTiming(sunTimerStruct* entry)
{
  long s_difference  = 0;
  long ns_difference = 0;

  sunclock_gettime_monotonic(entry->toc);

  s_difference  = entry->toc->tv_sec - entry->tic->tv_sec;
  ns_difference = entry->toc->tv_nsec - entry->tic->tv_nsec;
  if (ns_difference < 0)
  {
    s_difference--;
    ns_difference = 1000000000 + entry->toc->tv_nsec - entry->tic->tv_nsec;
  }

  entry->elapsed += ((double)s_difference) + ((double)ns_difference) * 1e-9;
  entry->average = entry->elapsed;
  entry->maximum = entry->elapsed;
}

static void sunResetTiming(sunTimerStruct* entry)
{
  entry->tic->tv_sec  = 0;
  entry->tic->tv_nsec = 0;
  entry->toc->tv_sec  = 0;
  entry->toc->tv_nsec = 0;
  entry->elapsed      = 0.0;
  entry->average      = 0.0;
  entry->maximum      = 0.0;
  entry->count        = 0;
}

SUNErrCode SUNProfiler_Create(SUNComm comm, const char* title, SUNProfiler* p)
{
  SUNProfiler profiler;
  int max_entries;
  char* max_entries_env;

  *p = profiler = (SUNProfiler)malloc(sizeof(struct SUNProfiler_));

  if (profiler == NULL) { return SUN_SUCCESS; }

  profiler->overhead = sunTimerStructNew();
  if (!profiler->overhead)
  {
    free(profiler);
    *p = profiler = NULL;
    return SUN_ERR_MALLOC_FAIL;
  }

  sunStartTiming(profiler->overhead);

  /* Check to see if max entries env variable was set, and use if it was. */
  max_entries     = 2560;
  max_entries_env = getenv("SUNPROFILER_MAX_ENTRIES");
  if (max_entries_env) { max_entries = atoi(max_entries_env); }
  if (max_entries <= 0) { max_entries = 2560; }

  /* Create the hashmap used to store the timers */
  if (SUNHashMap_New(max_entries, sunProfilerDestroyKeyValue, &profiler->map))
  {
    sunTimerStructFree((void*)profiler->overhead);
    free(profiler);
    *p = profiler = NULL;
    return SUN_ERR_MALLOC_FAIL;
  }

  /* Attach the comm, duplicating it if MPI is used. */
#if SUNDIALS_MPI_ENABLED
  profiler->comm = SUN_COMM_NULL;
  if (comm != SUN_COMM_NULL) { MPI_Comm_dup(comm, &profiler->comm); }
#else
  if (comm != SUN_COMM_NULL)
  {
    free(profiler);
    return -1;
  }
  profiler->comm = SUN_COMM_NULL;
#endif

  /* Copy the title of the profiler (note strlen does not include terminating
     null character hence the +1) */
  profiler->title = malloc((strlen(title) + 1) * sizeof(char));
  strcpy(profiler->title, title);

  /* Initialize the overall timer to 0. */
  profiler->sundials_time = 0.0;

  SUNDIALS_MARK_BEGIN(profiler, SUNDIALS_ROOT_TIMER);
  sunStopTiming(profiler->overhead);

  return SUN_SUCCESS;
}

SUNErrCode SUNProfiler_Free(SUNProfiler* p)
{
  if (!p || !(*p)) { return SUN_SUCCESS; }

  SUNDIALS_MARK_END(*p, SUNDIALS_ROOT_TIMER);

  if (*p)
  {
    SUNHashMap_Destroy(&(*p)->map);
    sunTimerStructFree((void*)(*p)->overhead);
#if SUNDIALS_MPI_ENABLED
    if ((*p)->comm != SUN_COMM_NULL) { MPI_Comm_free(&(*p)->comm); }
#endif
    free((*p)->title);
    free(*p);
  }
  *p = NULL;

  return SUN_SUCCESS;
}

SUNErrCode SUNProfiler_Begin(SUNProfiler p, const char* name)
{
  int64_t ier;
  sunTimerStruct* timer = NULL;

  if (!p) { return SUN_ERR_ARG_CORRUPT; }

  sunStartTiming(p->overhead);

  if (SUNHashMap_GetValue(p->map, name, (void**)&timer))
  {
    timer = sunTimerStructNew();
    ier   = SUNHashMap_Insert(p->map, name, (void*)timer);
    if (ier)
    {
      sunTimerStructFree(timer);
      sunStopTiming(p->overhead);
      if (ier == SUNHASHMAP_ERROR) { return SUN_ERR_PROFILER_MAPINSERT; }
      if (ier == SUNHASHMAP_DUPLICATE) { return SUN_ERR_PROFILER_MAPFULL; }
    }
  }

  timer->count++;
  sunStartTiming(timer);

  sunStopTiming(p->overhead);
  return SUN_SUCCESS;
}

SUNErrCode SUNProfiler_End(SUNProfiler p, const char* name)
{
  int64_t ier;
  sunTimerStruct* timer;

  if (!p) { return SUN_ERR_ARG_CORRUPT; }

  sunStartTiming(p->overhead);

  ier = SUNHashMap_GetValue(p->map, name, (void**)&timer);
  if (ier)
  {
    sunStopTiming(p->overhead);
    if (ier == SUNHASHMAP_ERROR) { return SUN_ERR_PROFILER_MAPGET; }
    if (ier == SUNHASHMAP_KEYNOTFOUND)
    {
      return SUN_ERR_PROFILER_MAPKEYNOTFOUND;
    }
  }

  sunStopTiming(timer);

  sunStopTiming(p->overhead);
  return SUN_SUCCESS;
}

SUNErrCode SUNProfiler_GetTimerResolution(SUNProfiler p, double* resolution)
{
  if (!p) { return SUN_ERR_ARG_CORRUPT; }

#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
  sunTimespec spec;
  clock_getres(CLOCK_MONOTONIC, &spec);
  *resolution = 1e-9 * ((double)spec.tv_nsec);

  return SUN_SUCCESS;
#elif (defined(WIN32) || defined(_WIN32))
  static LARGE_INTEGER ticks_per_sec;

  if (!ticks_per_sec.QuadPart)
  {
    QueryPerformanceFrequency(&ticks_per_sec);
    if (!ticks_per_sec.QuadPart) { return -1; }
  }

  *resolution = (double)ticks_per_sec.QuadPart;

  return SUN_SUCCESS;
#else
#error SUNProfiler needs POSIX or Windows timers
#endif
}

SUNErrCode SUNProfiler_GetElapsedTime(SUNProfiler p, const char* name,
                                      double* time)
{
  sunTimerStruct* timer;

  if (!p) { return SUN_ERR_ARG_CORRUPT; }

  if (SUNHashMap_GetValue(p->map, name, (void**)&timer)) { return (-1); }

  *time = timer->elapsed;

  return SUN_SUCCESS;
}

SUNErrCode SUNProfiler_Reset(SUNProfiler p)
{
  if (!p) { return SUN_ERR_ARG_CORRUPT; }

  /* Reset the overhead timer */
  sunResetTiming(p->overhead);
  sunStartTiming(p->overhead);

  /* Reset all timers */
  for (int64_t i = 0; i < SUNHashMap_Capacity(p->map); i++)
  {
    SUNHashMapKeyValue* kvp =
      SUNStlVector_SUNHashMapKeyValue_At(p->map->buckets, i);

    if (!kvp || !(*kvp)) { continue; }
    sunTimerStruct* timer = (*kvp)->value;
    if (timer) { sunResetTiming(timer); }
  }

  /* Reset the overall timer. */
  p->sundials_time = 0.0;

  SUNDIALS_MARK_BEGIN(p, SUNDIALS_ROOT_TIMER);
  sunStopTiming(p->overhead);

  return SUN_SUCCESS;
}

SUNErrCode SUNProfiler_Print(SUNProfiler p, FILE* fp)
{
  int64_t ier                = 0;
  int rank                   = 0;
  sunTimerStruct* timer      = NULL;
  SUNHashMapKeyValue* sorted = NULL;

  if (!p) { return SUN_ERR_ARG_CORRUPT; }

  sunStartTiming(p->overhead);

  /* Get the total SUNDIALS time up to this point */
  SUNDIALS_MARK_END(p, SUNDIALS_ROOT_TIMER);
  SUNDIALS_MARK_BEGIN(p, SUNDIALS_ROOT_TIMER);

  ier = SUNHashMap_GetValue(p->map, SUNDIALS_ROOT_TIMER, (void**)&timer);
  if (ier == SUNHASHMAP_ERROR) { return SUN_ERR_PROFILER_MAPGET; }
  if (ier == SUNHASHMAP_KEYNOTFOUND) { return SUN_ERR_PROFILER_MAPKEYNOTFOUND; }
  p->sundials_time = timer->elapsed;

#if SUNDIALS_MPI_ENABLED
  if (p->comm != SUN_COMM_NULL)
  {
    MPI_Comm_rank(p->comm, &rank);
    /* Find the max and average time across all ranks */
    sunCollectTimers(p);
  }
#endif

  if (rank == 0)
  {
    double resolution;
    /* Sort the timers in descending order */
    if (SUNHashMap_Sort(p->map, &sorted, sunCompareTimes))
    {
      return SUN_ERR_PROFILER_MAPSORT;
    }
    SUNProfiler_GetTimerResolution(p, &resolution);
    fprintf(fp, "\n============================================================"
                "====================================================\n");
    fprintf(fp, "SUNDIALS GIT VERSION: %s\n", SUNDIALS_GIT_VERSION);
    fprintf(fp, "SUNDIALS PROFILER: %s\n", p->title);
    fprintf(fp, "TIMER RESOLUTION: %gs\n", resolution);
    fprintf(fp, "%-40s\t %% time (inclusive) \t max/rank \t average/rank \t count \n",
            "RESULTS:");
    fprintf(fp, "=============================================================="
                "==================================================\n");

#if SUNDIALS_MPI_ENABLED
    if (p->comm == SUN_COMM_NULL)
    {
      printf(
        "WARNING: no MPI communicator provided, times shown are for rank 0\n");
    }
#endif

    /* Print all the other timers out */
    for (int64_t i = 0; i < SUNHashMap_Capacity(p->map); i++)
    {
      if (sorted[i]) { sunPrintTimer(sorted[i], fp, (void*)p); }
    }
    free(sorted);
  }

  sunStopTiming(p->overhead);

  if (rank == 0)
  {
    /* Print out the total time and the profiler overhead */
    fprintf(fp, "%-40s\t %6.2f%% \t         %.6fs \t -- \t\t -- \n",
            "Est. profiler overhead", p->overhead->elapsed / p->sundials_time,
            p->overhead->elapsed);

    /* End of output */
    fprintf(fp, "\n");
  }

  return SUN_SUCCESS;
}

#if SUNDIALS_MPI_ENABLED
static void sunTimerStructReduceMaxAndSum(void* a, void* b, int* len,
                                          SUNDIALS_MAYBE_UNUSED MPI_Datatype* dType)
{
  sunTimerStruct* a_ts = (sunTimerStruct*)a;
  sunTimerStruct* b_ts = (sunTimerStruct*)b;
  int i;
  for (i = 0; i < *len; ++i)
  {
    b_ts[i].average += a_ts[i].elapsed;
    b_ts[i].maximum = SUNMAX(a_ts[i].maximum, b_ts[i].maximum);
  }
}

/* Find the max and average time across all ranks */
SUNErrCode sunCollectTimers(SUNProfiler p)
{
  int rank, nranks;

  MPI_Comm comm = p->comm;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nranks);

  sunTimerStruct** values = NULL;

  /* MPI restricts us to int, but the hashmap allows int64_t.
     We add a check here to make sure that the capacity does
     not exceed an int, although it is unlikely we ever will. */
  if (SUNHashMap_Capacity(p->map) > INT_MAX)
  {
    return SUN_ERR_PROFILER_MAPFULL;
  }
  int map_size = (int)SUNHashMap_Capacity(p->map);

  /* Extract the elapsed times from the hash map */
  SUNHashMap_Values(p->map, (void***)&values, sizeof(sunTimerStruct));
  sunTimerStruct* reduced =
    (sunTimerStruct*)malloc(map_size * sizeof(sunTimerStruct));
  if (!reduced) { return SUN_ERR_MALLOC_FAIL; }
  for (int i = 0; i < map_size; ++i) { reduced[i] = *values[i]; }

  /* Register MPI datatype for sunTimerStruct */
  MPI_Datatype tmp_type, MPI_sunTimerStruct;
  const int block_lens[2]     = {5, 1};
  const MPI_Datatype types[2] = {MPI_DOUBLE, MPI_LONG};
  const MPI_Aint displ[2]     = {offsetof(sunTimerStruct, tic),
                                 offsetof(sunTimerStruct, count)};
  MPI_Aint lb, extent;

  MPI_Type_create_struct(2, block_lens, displ, types, &tmp_type);
  MPI_Type_get_extent(tmp_type, &lb, &extent);
  extent = sizeof(sunTimerStruct);
  MPI_Type_create_resized(tmp_type, lb, extent, &MPI_sunTimerStruct);
  MPI_Type_commit(&MPI_sunTimerStruct);

  /* Register max and sum MPI reduction operations for our datatype */
  MPI_Op MPI_sunTimerStruct_MAXANDSUM;
  MPI_Op_create(sunTimerStructReduceMaxAndSum, 1, &MPI_sunTimerStruct_MAXANDSUM);

  /* Compute max and average time across all ranks */
  if (rank == 0)
  {
    MPI_Reduce(MPI_IN_PLACE, reduced, map_size, MPI_sunTimerStruct,
               MPI_sunTimerStruct_MAXANDSUM, 0, comm);
  }
  else
  {
    MPI_Reduce(reduced, reduced, map_size, MPI_sunTimerStruct,
               MPI_sunTimerStruct_MAXANDSUM, 0, comm);
  }

  /* Cleanup custom MPI datatype and operations */
  MPI_Type_free(&tmp_type);
  MPI_Type_free(&MPI_sunTimerStruct);
  MPI_Op_free(&MPI_sunTimerStruct_MAXANDSUM);

  /* Update the values that are in this rank's hash map. */
  for (int i = 0; i < map_size; ++i)
  {
    values[i]->average = reduced[i].average / (double)nranks;
    values[i]->maximum = reduced[i].maximum;
  }

  free(reduced);
  free(values);

  return SUN_SUCCESS;
}
#endif

/* Print out the: timer name, percentage of exec time (based on the max),
   max across ranks, average across ranks, and the timer counter. */
void sunPrintTimer(SUNHashMapKeyValue kv, FILE* fp, void* pvoid)
{
  SUNProfiler p      = (SUNProfiler)pvoid;
  sunTimerStruct* ts = (sunTimerStruct*)kv->value;
  double maximum     = ts->maximum;
  double average     = ts->average;
  double percent = strcmp((const char*)kv->key, (const char*)SUNDIALS_ROOT_TIMER)
                     ? maximum / p->sundials_time * 100
                     : 100;
  fprintf(fp, "%-40s\t %6.2f%% \t         %.6fs \t %.6fs \t %ld\n", kv->key,
          percent, maximum, average, ts->count);
}

/* Comparator for qsort that compares key-value pairs
   based on the maximum time in the sunTimerStruct. */
int sunCompareTimes(const void* l, const void* r)
{
  double left_max;
  double right_max;

  const SUNHashMapKeyValue left  = *((const SUNHashMapKeyValue*)l);
  const SUNHashMapKeyValue right = *((const SUNHashMapKeyValue*)r);

  if (left == NULL && right == NULL) { return 0; }
  if (left == NULL) { return 1; }
  if (right == NULL) { return -1; }

  left_max  = ((sunTimerStruct*)left->value)->maximum;
  right_max = ((sunTimerStruct*)right->value)->maximum;

  if (left_max < right_max) { return 1; }
  if (left_max > right_max) { return -1; }

  return 0;
}

int sunclock_gettime_monotonic(sunTimespec* ts)
{
#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
  return clock_gettime(CLOCK_MONOTONIC, ts);
#elif (defined(WIN32) || defined(_WIN32))
  static LARGE_INTEGER ticks_per_sec;
  LARGE_INTEGER ticks;

  if (!ticks_per_sec.QuadPart)
  {
    QueryPerformanceFrequency(&ticks_per_sec);
    if (!ticks_per_sec.QuadPart) { return -1; }
  }

  QueryPerformanceCounter(&ticks);

  /* QueryPerformanceCounter is ticks in microseconds */

  ts->tv_sec  = (long)(ticks.QuadPart / ticks_per_sec.QuadPart);
  ts->tv_nsec = (long)(((ticks.QuadPart % ticks_per_sec.QuadPart) * 1000000) /
                       ticks_per_sec.QuadPart);

  return 0;
#else
#error SUNProfiler needs POSIX or Windows timers
#endif
}
