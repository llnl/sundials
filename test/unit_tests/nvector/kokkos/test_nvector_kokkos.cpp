/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel McGreer and Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 * This is the testing routine for the NVector implementation using Kokkos.
 * ---------------------------------------------------------------------------*/

#include <nvector/nvector_kokkos.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "test_nvector.h"

#if defined(USE_CUDA)
using ExecSpace = Kokkos::Cuda;
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

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  int fails{0};        /* counter for test failures */
  sunindextype length; /* vector length             */
  int print_timing;    /* turn timing on/off        */

  Test_Init(SUN_COMM_NULL);

  /* check input and set vector length */
  if (argc < 3)
  {
    printf("ERROR: TWO (2) Inputs required: vector length, print timing \n");
    return (-1);
  }

  length = (sunindextype)atol(argv[1]);
  if (length <= 0)
  {
    printf("ERROR: length of vector must be a positive integer \n");
    return (-1);
  }

  print_timing = atoi(argv[2]);
  SetTiming(print_timing, 0);

  printf("Testing KOKKOS N_Vector \n");

  printf("Vector length %ld \n\n", (long int)length);

  Kokkos::initialize(argc, argv);
  {
    VecType X{static_cast<SizeType>(length), sunctx};

    /* Check vector ID */
    fails += Test_N_VGetVectorID(X, SUNDIALS_NVEC_KOKKOS, 0);

    /* Test clone functions */
    fails += Test_N_VClone(X, length, 0);
    fails += Test_N_VCloneVectorArray(5, X, length, 0);

    /* Check vector length */
    fails += Test_N_VGetLength(X, 0);

    /* Check vector communicator */
    fails += Test_N_VGetCommunicator(X, SUN_COMM_NULL, 0);

    /* Clone additional vectors for testing */
    VecType Y{X};
    VecType Z{X};

    /* Standard vector operation tests */
    printf("\nTesting standard vector operations:\n\n");

    fails += Test_N_VAbs(X, Z, length, 0);
    fails += Test_N_VAddConst(X, Z, length, 0);
    fails += Test_N_VCompare(X, Z, length, 0);
    fails += Test_N_VConst(X, length, 0);
    fails += Test_N_VConstrMask(X, Y, Z, length, 0);
    fails += Test_N_VDiv(X, Y, Z, length, 0);
    fails += Test_N_VDotProd(X, Y, length, 0);
    fails += Test_N_VInv(X, Z, length, 0);
    fails += Test_N_VInvTest(X, Z, length, 0);
    fails += Test_N_VL1Norm(X, length, 0);
    fails += Test_N_VLinearSum(X, Y, Z, length, 0);
    fails += Test_N_VMaxNorm(X, length, 0);
    fails += Test_N_VMin(X, length, 0);
    fails += Test_N_VMinQuotient(X, Y, length, 0);
    fails += Test_N_VProd(X, Y, Z, length, 0);
    fails += Test_N_VScale(X, Z, length, 0);
    fails += Test_N_VWL2Norm(X, Y, length, 0);
    fails += Test_N_VWrmsNorm(X, Y, length, 0);
    fails += Test_N_VWrmsNormMask(X, Y, Z, length, 0);

    /* Fused and vector array operations tests (disabled) */
    printf("\nTesting fused and vector array operations (disabled):\n\n");

    /* create vector and test vector array operations */
    VecType U{X};

    /* fused operations */
    fails += Test_N_VLinearCombination(U, length, 0);
    fails += Test_N_VScaleAddMulti(U, length, 0);
    fails += Test_N_VDotProdMulti(U, length, 0);

    /* vector array operations */
    fails += Test_N_VLinearSumVectorArray(U, length, 0);
    fails += Test_N_VScaleVectorArray(U, length, 0);
    fails += Test_N_VConstVectorArray(U, length, 0);
    fails += Test_N_VWrmsNormVectorArray(U, length, 0);
    fails += Test_N_VWrmsNormMaskVectorArray(U, length, 0);
    fails += Test_N_VScaleAddMultiVectorArray(U, length, 0);
    fails += Test_N_VLinearCombinationVectorArray(U, length, 0);

    /* local reduction operations */
    printf("\nTesting local reduction operations:\n\n");

    fails += Test_N_VDotProdLocal(X, Y, length, 0);
    fails += Test_N_VMaxNormLocal(X, length, 0);
    fails += Test_N_VMinLocal(X, length, 0);
    fails += Test_N_VL1NormLocal(X, length, 0);
    fails += Test_N_VWSqrSumLocal(X, Y, length, 0);
    fails += Test_N_VWSqrSumMaskLocal(X, Y, Z, length, 0);
    fails += Test_N_VInvTestLocal(X, Z, length, 0);
    fails += Test_N_VConstrMaskLocal(X, Y, Z, length, 0);
    fails += Test_N_VMinQuotientLocal(X, Y, length, 0);
  }
  Kokkos::finalize();

  /* Print result */
  if (fails) { printf("FAIL: NVector module failed %i tests \n\n", fails); }
  else { printf("SUCCESS: NVector module passed all tests \n\n"); }

  Test_Finalize();

  return (fails);
}

/* ----------------------------------------------------------------------
 * Implementation specific utility functions for vector tests
 * --------------------------------------------------------------------*/

int check_ans(sunrealtype ans, N_Vector X, sunindextype local_length)
{
  int failure{0};
  auto Xvec{static_cast<VecType*>(X->content)};
  auto Xdata{Xvec->HostView()};

  sundials::kokkos::CopyFromDevice<VecType>(*Xvec);
  for (sunindextype i = 0; i < local_length; i++)
  {
    failure += SUNRCompare(Xdata[i], ans);
  }

  return (failure > ZERO) ? (1) : (0);
}

sunbooleantype has_data(N_Vector X)
{
  /* check if vector data is non-null */
  return SUNTRUE;
}

void set_element(N_Vector X, sunindextype i, sunrealtype val)
{
  /* set i-th element of data array */
  set_element_range(X, i, i, val);
}

void set_element_range(N_Vector X, sunindextype is, sunindextype ie,
                       sunrealtype val)
{
  auto Xvec{static_cast<VecType*>(X->content)};
  auto Xdata{Xvec->HostView()};

  /* set elements [is,ie] of the data array */
  sundials::kokkos::CopyFromDevice<VecType>(X);
  for (sunindextype i = is; i <= ie; i++) { Xdata[i] = val; }
  sundials::kokkos::CopyToDevice<VecType>(X);
}

sunrealtype get_element(N_Vector X, sunindextype i)
{
  /* get i-th element of data array */
  auto Xvec{static_cast<VecType*>(X->content)};
  auto Xdata{Xvec->HostView()};
  sundials::kokkos::CopyFromDevice<VecType>(X);
  return Xdata[i];
}

double max_time(N_Vector X, double time)
{
  /* not running in parallel, just return input time */
  return time;
}

void sync_device(N_Vector x)
{
  /* sync with GPU */
  Kokkos::fence();
  return;
}
