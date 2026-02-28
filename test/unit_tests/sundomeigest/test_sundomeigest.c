/*
 * -----------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ SMU
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
 * These test functions are designed to check a SUNDomEigEstimator
 * module implementation.
 * -----------------------------------------------------------------
 */

#include "test_sundomeigest.h"

#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
#include <time.h>
#include <unistd.h>
#endif

/* private functions */
static double get_time(void);

int print_time = 0;

#define PRINT_TIME(format, time) \
  if (print_time) printf(format, time)

/* ----------------------------------------------------------------------
 * SUNDomEigEstimator_SetATimes Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstimator_SetATimes(SUNDomEigEstimator DEE, void* ATdata,
                                      SUNATimesFn ATimes, int myid)
{
  int failure;
  double start_time, stop_time;

  /* try calling SetATimes routine: should pass/fail based on expected input */
  start_time = get_time();
  failure    = SUNDomEigEstimator_SetATimes(DEE, ATdata, ATimes);
  stop_time  = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstimator_SetATimes returned %d on "
           "Proc %d "
           "\n",
           failure, myid);
    PRINT_TIME("    SUNDomEigEstimator_SetATimes Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstimator_SetATimes \n");
    PRINT_TIME("    SUNDomEigEstimator_SetATimes Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEstimator_SetMaxIters Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstimator_SetMaxIters(SUNDomEigEstimator DEE,
                                        long int max_iters, int myid)
{
  int failure;
  double start_time, stop_time;

  /* try calling SetMaxIters routine: should pass/fail based on expected input */
  start_time = get_time();
  failure    = SUNDomEigEstimator_SetMaxIters(DEE, max_iters);
  stop_time  = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstimator_SetMaxIters returned %d on "
           "Proc "
           "%d \n",
           failure, myid);
    PRINT_TIME("    SUNDomEigEstimator_SetMaxIters Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstimator_SetMaxIters \n");
    PRINT_TIME("    SUNDomEigEstimator_SetMaxIters Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

int Test_SUNDomEigEstimator_SetNumPreprocessIters(SUNDomEigEstimator DEE,
                                                  int num_warmups, int myid)
{
  int failure;
  double start_time, stop_time;

  /* try calling SUNDomEigEstimator_SetNumPreprocessIters routine: should pass/fail based on expected input */
  start_time = get_time();
  failure    = SUNDomEigEstimator_SetNumPreprocessIters(DEE, num_warmups);
  stop_time  = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstimator_SetNumPreprocessIters check, "
           "Proc "
           "%d \n",
           myid);
    PRINT_TIME("    SUNDomEigEstimator_SetNumPreprocessIters Time: %22.15e \n "
               "\n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstimator_SetNumPreprocessIters \n");
    PRINT_TIME("    SUNDomEigEstimator_SetNumPreprocessIters Time: %22.15e \n "
               "\n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEstimator_SetRelTol Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstimator_SetRelTol(SUNDomEigEstimator DEE, sunrealtype tol,
                                      int myid)
{
  int failure;
  double start_time, stop_time;

  /* try calling SUNDomEigEstimator_SetRelTol routine: should pass/fail based on expected input */
  start_time = get_time();
  failure    = SUNDomEigEstimator_SetRelTol(DEE, tol);
  stop_time  = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstimator_SetRelTol check, Proc %d \n",
           myid);
    PRINT_TIME("    SUNDomEigEstimator_SetRelTol Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstimator_SetRelTol \n");
    PRINT_TIME("    SUNDomEigEstimator_SetRelTol Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEstimator_SetInitialGuess Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstimator_SetInitialGuess(SUNDomEigEstimator DEE, N_Vector q,
                                            int myid)
{
  int failure;
  double start_time, stop_time;

  /* try calling SUNDomEigEstimator_SetInitialGuess routine: should pass/fail based on expected input */
  start_time = get_time();
  failure    = SUNDomEigEstimator_SetInitialGuess(DEE, q);
  stop_time  = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstimator_SetInitialGuess check, Proc "
           "%d \n",
           myid);
    PRINT_TIME("    SUNDomEigEstimator_SetInitialGuess Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstimator_SetInitialGuess \n");
    PRINT_TIME("    SUNDomEigEstimator_SetInitialGuess Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEstimator_Initialize Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstimator_Initialize(SUNDomEigEstimator DEE, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEstimator_Initialize(DEE);
  stop_time  = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstimator_Initialize check, Proc %d \n",
           myid);
    PRINT_TIME("    SUNDomEigEstimator_Initialize Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstimator_Initialize \n");
    PRINT_TIME("    SUNDomEigEstimator_Initialize Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEstimator_Estimate Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstimator_Estimate(SUNDomEigEstimator DEE, sunrealtype* lambdaR,
                                     sunrealtype* lambdaI, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEstimator_Estimate(DEE, lambdaR, lambdaI);
  stop_time  = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstimator_Estimate check, Proc %d \n",
           myid);
    PRINT_TIME("    SUNDomEigEstimator_Estimate Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstimator_Estimate \n");
    PRINT_TIME("    SUNDomEigEstimator_Estimate Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEstimator_GetRes Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstimator_GetRes(SUNDomEigEstimator DEE, sunrealtype* cur_res,
                                   int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEstimator_GetRes(DEE, cur_res);
  stop_time  = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstimator_GetRes check, Proc %d \n",
           myid);
    PRINT_TIME("    SUNDomEigEstimator_GetRes Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstimator_GetRes \n");
    PRINT_TIME("    SUNDomEigEstimator_GetRes Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEstimator_GetNumIters Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstimator_GetNumIters(SUNDomEigEstimator DEE,
                                        long int* curniter, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEstimator_GetNumIters(DEE, curniter);
  stop_time  = get_time();

  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstimator_GetNumIters check, Proc %d "
           "\n",
           myid);
    PRINT_TIME("    SUNDomEigEstimator_GetNumIters Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstimator_GetNumIters \n");
    PRINT_TIME("    SUNDomEigEstimator_GetNumIters Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEstimator_GetNumATimesCalls Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstimator_GetNumATimesCalls(SUNDomEigEstimator DEE,
                                              long int* num_ATimes, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();
  failure    = SUNDomEigEstimator_GetNumATimesCalls(DEE, num_ATimes);
  stop_time  = get_time();
  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstimator_GetNumATimesCalls check, "
           "Proc %d "
           "\n",
           myid);
    PRINT_TIME("    SUNDomEigEstimator_GetNumATimesCalls Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstimator_GetNumATimesCalls \n");
    PRINT_TIME("    SUNDomEigEstimator_GetNumATimesCalls Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ----------------------------------------------------------------------
 * SUNDomEigEstimator_Write Test
 * --------------------------------------------------------------------*/
int Test_SUNDomEigEstimator_Write(SUNDomEigEstimator DEE, int myid)
{
  int failure;
  double start_time, stop_time;

  start_time = get_time();

  printf("    SUNDomEigEstimator_Write: \n");
  printf("    --------------------------------------------\n");
  failure = SUNDomEigEstimator_Write(DEE, stdout);
  printf("    --------------------------------------------\n");
  stop_time = get_time();
  if (failure)
  {
    printf(">>> FAILED test -- SUNDomEigEstimator_Write check, Proc %d \n", myid);
    PRINT_TIME("    SUNDomEigEstimator_Write Time: %22.15e \n \n",
               stop_time - start_time);
    return (1);
  }
  else if (myid == 0)
  {
    printf("    PASSED test -- SUNDomEigEstimator_Write \n");
    PRINT_TIME("    SUNDomEigEstimator_Write Time: %22.15e \n \n",
               stop_time - start_time);
  }

  return (0);
}

/* ======================================================================
 * Private functions
 * ====================================================================*/

#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
time_t base_time_tv_sec = 0; /* Base time; makes time values returned
                                by get_time easier to read when
                                printed since they will be zero
                                based.
                              */
#endif

void SetTiming(int onoff)
{
  print_time = onoff;

#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
  struct timespec spec;
  clock_gettime(CLOCK_MONOTONIC, &spec);
  base_time_tv_sec = spec.tv_sec;
#endif
}

/* ----------------------------------------------------------------------
 * Timer
 * --------------------------------------------------------------------*/
static double get_time(void)
{
#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
  struct timespec spec;
  clock_gettime(CLOCK_MONOTONIC, &spec);
  double time = (double)(spec.tv_sec - base_time_tv_sec) +
                ((double)(spec.tv_nsec) / 1E9);
#else
  double time = 0;
#endif
  return time;
}
