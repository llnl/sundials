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
 * This is the testing routine to check the performance of the
 * NVECTOR Kokkos module implementation.
 * -----------------------------------------------------------------*/

#include <nvector/nvector_kokkos.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_core.hpp>

#include "test_nvector_performance.h"

#if defined(USE_CUDA)
using ExecSpace = Kokkos::DefaultExecutionSpace;
// using ExecSpace = Kokkos::Cuda;
#elif defined(USE_HIP)
#if KOKKOS_VERSION / 10000 > 3
using ExecSpace = Kokkos::HIP;
#else
using ExecSpace = Kokkos::Experimental::HIP;
#endif
#elif defined(USE_OPENMP)
using ExecSpace = Kokkos::OpenMP;
#else
using ExecSpace = Kokkos::Serial;
#endif

using VecType  = sundials::kokkos::Vector<ExecSpace>;
using SizeType = VecType::size_type;

/* private functions */
static int InitializeClearCache(size_t cachesize);
static int FinalizeClearCache();

/* private data for clearing cache */
struct Cache
{
  Kokkos::View<sunrealtype*, ExecSpace::memory_space> device;
  Kokkos::View<sunrealtype*, ExecSpace::memory_space>::HostMirror host;
};

static Cache* cache;

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  sundials::Context sunctx;
  int flag; /* return flag     */

  /* CLI args */
  sunindextype veclen; /* vector length   */
  int print_timing;    /* output timings  */
  int ntests;          /* number of tests */
  int nvecs;           /* number of tests */
  int nsums;           /* number of sums  */
  size_t cachesize;    /* cache size      */

  Kokkos::initialize(argc, argv);
  {
    printf("\nStart Tests\n");
    printf("Vector Name: Kokkos\n");

    /* check input and set vector length */
    if (argc < 6)
    {
      printf("ERROR: SIX (6) arguments required: ");
      printf("<vector length> <number of vectors> <number of sums> <number of "
             "tests> ");
      printf("<cachesize> <print timing>\n");
      return (-1);
    }

#ifdef SUNDIALS_INT64_T
    veclen = atol(argv[1]);
#else
    veclen = atoi(argv[1]);
#endif
    if (veclen <= 0)
    {
      printf("ERROR: length of vector must be a positive integer \n");
      return (-1);
    }

    nvecs = atoi(argv[2]);
    if (nvecs < 1) { printf("WARNING: Fused operation tests disabled\n"); }

    nsums = atoi(argv[3]);
    if (nsums < 1) { printf("WARNING: Some fused operation tests disabled\n"); }

    ntests = atoi(argv[4]);
    if (ntests <= 0)
    {
      printf("ERROR: number of tests must be a positive integer \n");
      return (-1);
    }

    cachesize = static_cast<size_t>(atol(argv[5]));
    InitializeClearCache(cachesize);

    print_timing = atoi(argv[6]);
    SetTiming(print_timing, 0);

    printf("\nRunning with: \n");
    printf("  vector length         %ld \n", (long int)veclen);
    printf("  max number of vectors %d  \n", nvecs);
    printf("  max number of sums    %d  \n", nsums);
    printf("  number of tests       %d  \n", ntests);
    printf("  timing on/off         %d  \n", print_timing);

    /* Create vectors */
    VecType X{static_cast<SizeType>(veclen), sunctx};

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
      if (print_timing)
      {
        printf("\n\n fused operations 1: nvecs= %d\n", nvecs);
      }
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

    sync_device(X);
    FinalizeClearCache();

    printf("\nFinished Tests\n");
  }
  Kokkos::finalize();

  return flag;
}

/* ----------------------------------------------------------------------
 * Functions required by testing routines to fill vector data
 * --------------------------------------------------------------------*/

/* random data between lower and upper */
void N_VRand(N_Vector Xvec, sunindextype Xlen, sunrealtype lower,
             sunrealtype upper)
{
  auto X{sundials::kokkos::GetVec<VecType>(Xvec)};
  rand_realtype(X->HostView().data(), Xlen, lower, upper);
  sundials::kokkos::CopyToDevice<VecType>(*X);
}

/* series of 0 and 1 */
void N_VRandZeroOne(N_Vector Xvec, sunindextype Xlen)
{
  auto X{sundials::kokkos::GetVec<VecType>(Xvec)};
  rand_realtype_zero_one(X->HostView().data(), Xlen);
  sundials::kokkos::CopyToDevice<VecType>(*X);
}

/* random values for constraint array */
void N_VRandConstraints(N_Vector Xvec, sunindextype Xlen)
{
  auto X{sundials::kokkos::GetVec<VecType>(Xvec)};
  rand_realtype_constraints(X->HostView().data(), Xlen);
  sundials::kokkos::CopyToDevice<VecType>(*X);
}

/* ----------------------------------------------------------------------
 * Functions required for MPI or GPU testing
 * --------------------------------------------------------------------*/

void collect_times(N_Vector Xvec, double* times, int ntimes)
{
  /* not running with MPI, just return */
  return;
}

void sync_device(N_Vector x)
{
  Kokkos::fence();
  return;
}

/* ----------------------------------------------------------------------
 * Functions required for clearing cache
 * --------------------------------------------------------------------*/

int InitializeClearCache(size_t cachesize)
{
  if (!cachesize) { return 0; }

  cache = new Cache();

  // determine size of vector to clear cache, N = ceil(2 * nbytes/sunrealtype)
  size_t nbytes  = static_cast<size_t>(2 * cachesize * 1024 * 1024);
  sunindextype N = static_cast<sunindextype>(
    (nbytes + sizeof(sunrealtype) - 1) / sizeof(sunrealtype));

  cache->device =
    Kokkos::View<sunrealtype*, ExecSpace::memory_space>("device cache", N);
  cache->host = Kokkos::create_mirror_view(cache->device);

  rand_realtype(cache->host.data(), N, sunrealtype{-1.0}, sunrealtype{1.0});
  Kokkos::deep_copy(cache->device, cache->host);

  return 0;
}

int FinalizeClearCache()
{
  delete cache;
  return (0);
}

void ClearCache()
{
  if (cache->device.data())
  {
    auto device_cache{cache->device};
    sunrealtype sum{0.0};
    Kokkos::parallel_reduce(
      "ClearCache",
      Kokkos::RangePolicy<ExecSpace>(0, static_cast<SizeType>(
                                          device_cache.extent(0))),
      KOKKOS_LAMBDA(const SizeType i, sunrealtype& accum) {
        accum += device_cache(i);
      },
      sum);
  }
  return;
}
