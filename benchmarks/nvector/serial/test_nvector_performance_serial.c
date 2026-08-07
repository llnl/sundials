/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * This is the testing routine to evaluate the performance of the
 * serial NVECTOR module implementation.
 * -----------------------------------------------------------------*/

#include <nvector/nvector_serial.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "test_nvector_performance.h"

/* private functions */
static int InitializeClearCache(int cachesize);
static int FinalizeClearCache();

/* private data for clearing cache */
static sunindextype N;    /* data length */
static sunrealtype* data; /* host data   */

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  SUNContext ctx = NULL; /* SUNDIALS context */
  N_Vector X     = NULL; /* test vector      */
  sunindextype veclen;   /* vector length    */

  int print_timing; /* output timings     */
  int ntests;       /* number of tests    */
  int nvecs;        /* number of tests    */
  int nsums;        /* number of sums     */
  int cachesize;    /* size of cache (MB) */
  int flag;         /* return flag        */

  printf("\nStart Tests\n");
  printf("Vector Name: Serial\n");

  /* check input and set vector length */
  if (argc < 7)
  {
    printf("ERROR: SIX (6) arguments required: ");
    printf("<vector length> <number of vectors> <number of sums> <number of "
           "tests> ");
    printf("<cache size (MB)> <print timing>\n");
    return (-1);
  }

  veclen = (sunindextype)atol(argv[1]);
  if (veclen <= 0)
  {
    printf("ERROR: length of vector must be a positive integer \n");
    return (-1);
  }

  nvecs = (int)atol(argv[2]);
  if (nvecs < 1) { printf("WARNING: Fused operation test disabled\n"); }

  nsums = (int)atol(argv[3]);
  if (nsums < 1) { printf("WARNING: Some fused operation tests disabled\n"); }

  ntests = (int)atol(argv[4]);
  if (ntests <= 0)
  {
    printf("ERROR: number of tests must be a positive integer \n");
    return (-1);
  }

  cachesize = (int)atol(argv[5]);
  if (cachesize < 0)
  {
    printf("ERROR: cache size (MB) must be a non-negative integer \n");
    return (-1);
  }
  InitializeClearCache(cachesize);

  print_timing = atoi(argv[6]);
  SetTiming(print_timing, 0);

  printf("\nRunning with: \n");
  printf("  vector length         %ld \n", (long int)veclen);
  printf("  max number of vectors %d  \n", nvecs);
  printf("  max number of sums    %d  \n", nsums);
  printf("  number of tests       %d  \n", ntests);
  printf("  timing on/off         %d  \n", print_timing);

  flag = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (flag) { return flag; }

  /* Create vectors */
  X = N_VNew_Serial(veclen, ctx);

  /* run tests */
  if (print_timing) { printf("\n\n standard operations:\n"); }
  if (print_timing) { PrintTableHeader(1); }
  flag = Test_N_VLinearSum(X, veclen, ntests);
  flag = Test_N_VConst(X, veclen, ntests);
  flag = Test_N_VProd(X, veclen, ntests);
  flag = Test_N_VDiv(X, veclen, ntests);
  flag = Test_N_VScale(X, veclen, ntests);
  flag = Test_N_VAbs(X, veclen, ntests);
  flag = Test_N_VInv(X, veclen, ntests);
  flag = Test_N_VAddConst(X, veclen, ntests);
  flag = Test_N_VDotProd(X, veclen, ntests);
  flag = Test_N_VMaxNorm(X, veclen, ntests);
  flag = Test_N_VWrmsNorm(X, veclen, ntests);
  flag = Test_N_VWrmsNormMask(X, veclen, ntests);
  flag = Test_N_VMin(X, veclen, ntests);
  flag = Test_N_VWL2Norm(X, veclen, ntests);
  flag = Test_N_VL1Norm(X, veclen, ntests);
  flag = Test_N_VCompare(X, veclen, ntests);
  flag = Test_N_VInvTest(X, veclen, ntests);
  flag = Test_N_VConstrMask(X, veclen, ntests);
  flag = Test_N_VMinQuotient(X, veclen, ntests);

  if (nvecs > 0)
  {
    if (print_timing) { printf("\n\n fused operations 1: nvecs= %d\n", nvecs); }
    if (print_timing) { PrintTableHeader(2); }
    flag = Test_N_VLinearCombination(X, veclen, nvecs, ntests);
    flag = Test_N_VScaleAddMulti(X, veclen, nvecs, ntests);
    flag = Test_N_VDotProdMulti(X, veclen, nvecs, ntests);
    flag = Test_N_VLinearSumVectorArray(X, veclen, nvecs, ntests);
    flag = Test_N_VScaleVectorArray(X, veclen, nvecs, ntests);
    flag = Test_N_VConstVectorArray(X, veclen, nvecs, ntests);
    flag = Test_N_VWrmsNormVectorArray(X, veclen, nvecs, ntests);
    flag = Test_N_VWrmsNormMaskVectorArray(X, veclen, nvecs, ntests);

    if (nsums > 0)
    {
      if (print_timing)
      {
        printf("\n\n fused operations 2: nvecs= %d nsums= %d\n", nvecs, nsums);
      }
      if (print_timing) { PrintTableHeader(2); }
      flag = Test_N_VScaleAddMultiVectorArray(X, veclen, nvecs, nsums, ntests);
      flag = Test_N_VLinearCombinationVectorArray(X, veclen, nvecs, nsums,
                                                  ntests);
    }
  }

  /* Free vectors */
  N_VDestroy(X);

  FinalizeClearCache();

  flag = SUNContext_Free(&ctx);
  if (flag) { return flag; }

  printf("\nFinished Tests\n");

  return (flag);
}

/* ----------------------------------------------------------------------
 * Functions required by testing routines to fill vector data
 * --------------------------------------------------------------------*/

/* random data between lower and upper */
void N_VRand(N_Vector Xvec, sunindextype Xlen, sunrealtype lower,
             sunrealtype upper)
{
  sunrealtype* Xdata;

  Xdata = N_VGetArrayPointer(Xvec);
  rand_realtype(Xdata, Xlen, lower, upper);
}

/* series of 0 and 1 */
void N_VRandZeroOne(N_Vector Xvec, sunindextype Xlen)
{
  sunrealtype* Xdata;

  Xdata = N_VGetArrayPointer(Xvec);
  rand_realtype_zero_one(Xdata, Xlen);
}

/* random values for constraint array */
void N_VRandConstraints(N_Vector Xvec, sunindextype Xlen)
{
  sunrealtype* Xdata;

  Xdata = N_VGetArrayPointer(Xvec);
  rand_realtype_constraints(Xdata, Xlen);
}

/* ----------------------------------------------------------------------
 * Functions required for MPI or GPU testing
 * --------------------------------------------------------------------*/

void collect_times(N_Vector X, double* times, int ntimes)
{
  /* not running with MPI, just return */
  return;
}

void sync_device(N_Vector x)
{
  /* not running on GPU, just return */
  return;
}

/* ----------------------------------------------------------------------
 * Functions required for clearing cache
 * --------------------------------------------------------------------*/

static int InitializeClearCache(int cachesize)
{
  size_t nbytes; /* cache size in bytes */

  if (!cachesize)
  {
    N    = 0;
    data = NULL;
    return 0;
  }

  /* determine size of vector to clear cache, N = ceil(2 * nbytes/sunrealtype) */
  nbytes = (size_t)(2 * cachesize * 1024 * 1024);
  N = (sunindextype)((nbytes + sizeof(sunrealtype) - 1) / sizeof(sunrealtype));

  /* allocate data and fill random values */
  data = (sunrealtype*)malloc(N * sizeof(sunrealtype));
  rand_realtype(data, N, SUN_RCONST(-1.0), SUN_RCONST(1.0));

  return (0);
}

static int FinalizeClearCache()
{
  if (data) { free(data); }
  return (0);
}

void ClearCache()
{
  if (data)
  {
    sunrealtype sum;
    sunindextype i;

    sum = SUN_RCONST(0.0);
    for (i = 0; i < N; i++) { sum += data[i]; }
    (void)sum;
  }

  return;
}
