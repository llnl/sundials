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
 *
 * Three construction scenarios are exercised:
 *
 *   A. Allocating constructor     -- Vector(length, sunctx)
 *   B. Managed view constructor   -- Vector(managed_view, sunctx)
 *   C. Unmanaged view constructor -- Vector(unmanaged_view, sunctx)
 *
 * Additional vectors are obtained via N_VClone, which always returns a managed
 * Vector, to mirror the pattern in SUNDIALS packages where the user supplies
 * an initial vector and the integrator clones it to create internal workspace.
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

using VecType          = sundials::kokkos::Vector<ExecSpace>;
using UnmanagedVecType = sundials::kokkos::Vector<ExecSpace,
                                                   ExecSpace::memory_space,
                                                   Kokkos::MemoryUnmanaged>;
using SizeType         = VecType::size_type;

/* -----------------------------------------------------------------------------
 * Run the vector operation tests
 * ---------------------------------------------------------------------------*/
static int run_all_tests(N_Vector X, N_Vector Y, N_Vector Z, sunindextype length)
{
  int fails{0};

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

  printf("\nTesting fused and vector array operations (disabled):\n\n");

  fails += Test_N_VLinearCombination(X, length, 0);
  fails += Test_N_VScaleAddMulti(X, length, 0);
  fails += Test_N_VDotProdMulti(X, length, 0);
  fails += Test_N_VLinearSumVectorArray(X, length, 0);
  fails += Test_N_VScaleVectorArray(X, length, 0);
  fails += Test_N_VConstVectorArray(X, length, 0);
  fails += Test_N_VWrmsNormVectorArray(X, length, 0);
  fails += Test_N_VWrmsNormMaskVectorArray(X, length, 0);
  fails += Test_N_VScaleAddMultiVectorArray(X, length, 0);
  fails += Test_N_VLinearCombinationVectorArray(X, length, 0);

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

  return fails;
}

/* -----------------------------------------------------------------------------
 * Main NVector Testing Routine
 * ---------------------------------------------------------------------------*/
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

  Kokkos::initialize(argc, argv);
  {
    /* -------------------------------------------------------------------------
     * Scenario A: allocating constructor
     *
     * X owns its memory. Y and Z are managed clones of X. All three vectors
     * share the same VecType specialisation so there is no managed/unmanaged
     * mixing.
     * -----------------------------------------------------------------------*/
    printf("====================================================\n");
    printf("Testing KOKKOS N_Vector: allocating constructor\n");
    printf("Vector length %ld\n", (long int)length);
    printf("====================================================\n");

    {
      VecType X{static_cast<SizeType>(length), sunctx};

      fails += Test_N_VGetVectorID(X, SUNDIALS_NVEC_KOKKOS, 0);
      fails += Test_N_VClone(X, length, 0);
      fails += Test_N_VCloneVectorArray(5, X, length, 0);
      fails += Test_N_VGetLength(X, 0);
      fails += Test_N_VGetCommunicator(X, SUN_COMM_NULL, 0);

      N_Vector Y_nv = N_VClone(X);
      N_Vector Z_nv = N_VClone(X);

      fails += run_all_tests(X, Y_nv, Z_nv, length);

      N_VDestroy(Y_nv);
      N_VDestroy(Z_nv);
    }

    /* -------------------------------------------------------------------------
     * Scenario B: managed view constructor
     *
     * X is constructed from a pre-existing managed Kokkos view. The view is
     * shared between the caller and the Vector (reference counted). Y and Z are
     * managed clones. All three vectors are the same VecType specialisation; no
     * managed/unmanaged mixing occurs.
     *
     * This exercises the view constructor path and verifies that constructing
     * from a managed view does not allocate a second buffer.
     * -----------------------------------------------------------------------*/
    printf("\n====================================================\n");
    printf("Testing KOKKOS N_Vector: managed view constructor\n");
    printf("Vector length %ld\n", (long int)length);
    printf("====================================================\n");

    {
      VecType::view_type managed_view("managed_view",
                                      static_cast<SizeType>(length));
      VecType X{managed_view, sunctx};

      /* Verify zero-copy: X wraps the same allocation, no copy occurred */
      if (X.View().data() != managed_view.data())
      {
        printf("FAILED: managed view constructor: data pointer mismatch -- "
               "Vector did not wrap the provided view\n");
        fails++;
      }
      else
      {
        printf("PASSED: managed view constructor: data pointer matches "
               "provided view\n");
      }

      fails += Test_N_VGetVectorID(X, SUNDIALS_NVEC_KOKKOS, 0);
      fails += Test_N_VClone(X, length, 0);
      fails += Test_N_VCloneVectorArray(5, X, length, 0);
      fails += Test_N_VGetLength(X, 0);
      fails += Test_N_VGetCommunicator(X, SUN_COMM_NULL, 0);

      N_Vector Y_nv = N_VClone(X);
      N_Vector Z_nv = N_VClone(X);

      fails += run_all_tests(X, Y_nv, Z_nv, length);

      N_VDestroy(Y_nv);
      N_VDestroy(Z_nv);
    }

#if KOKKOS_VERSION / 10000 > 4
#warning "Unmanaged views are not currently supported with Kokkos 5+"
#else
    /* -------------------------------------------------------------------------
     * Scenario C: unmanaged view constructor
     *
     * X is constructed from an unmanaged view wrapping a raw pointer. Y and Z
     * are managed clones produced by N_VClone, which always returns a managed
     * VecType regardless of the type of the vector being cloned.
     *
     * The user supplies X and the integrator clones it to create managed
     * internal workspace vectors. Operations such as N_VScale(c, X, Z) then
     * dispatch through Z->ops (managed), operating on X (unmanaged) and Z
     * (managed) simultaneously.
     *
     * This is the exact configuration that triggered a Kokkos 5 abort in
     * SharedAllocationRecord::increment before the to_managed fix was applied
     * to the Vector view constructor.
     * -----------------------------------------------------------------------*/
    printf("\n====================================================\n");
    printf("Testing KOKKOS N_Vector: unmanaged view constructor\n");
    printf("Vector length %ld\n", (long int)length);
    printf("====================================================\n");

    /* Raw pointer allocated in the correct memory space for the execution
     * backend. Lifetime spans the op tests so that it outlives the Vector
     * constructed over it. */
    sunrealtype* unmanaged_ptr =
      static_cast<sunrealtype*>(Kokkos::kokkos_malloc<ExecSpace::memory_space>(
        "unmanaged_storage",
        static_cast<size_t>(length) * sizeof(sunrealtype)));

    {
      UnmanagedVecType::view_type unmanaged_view(unmanaged_ptr,
                                                 static_cast<SizeType>(length));
      UnmanagedVecType X{unmanaged_view, sunctx};

      /* Verify zero-copy: X wraps the raw pointer without copying */
      if (X.View().data() != unmanaged_ptr)
      {
        printf("FAILED: unmanaged view constructor: data pointer mismatch -- "
               "Vector did not wrap the provided view\n");
        fails++;
      }
      else
      {
        printf("PASSED: unmanaged view constructor: data pointer matches "
               "unmanaged storage\n");
      }

      fails += Test_N_VGetVectorID(X, SUNDIALS_NVEC_KOKKOS, 0);

      /* N_VClone on an unmanaged vector must return a managed clone that
       * owns separate memory -- it must not alias the raw pointer. */
      {
        N_Vector clone = N_VClone(X);
        auto clone_vec{sundials::kokkos::GetVec<VecType>(clone)};
        if (clone_vec->View().data() == unmanaged_ptr)
        {
          printf("FAILED: N_VClone of unmanaged vector aliased the original "
                 "unmanaged storage -- clone must own independent memory\n");
          fails++;
        }
        else
        {
          printf("PASSED: N_VClone of unmanaged vector owns independent "
                 "memory\n");
        }
        N_VDestroy(clone);
      }

      fails += Test_N_VCloneVectorArray(5, X, length, 0);
      fails += Test_N_VGetLength(X, 0);
      fails += Test_N_VGetCommunicator(X, SUN_COMM_NULL, 0);

      /* Y and Z are managed clones of the unmanaged X. Vector operations will
       * mix unmanaged (X) and managed views (Y, Z), exercising the type-mixing
       * code path that failed in Kokkos 5. */
      N_Vector Y_nv = N_VClone(X);
      N_Vector Z_nv = N_VClone(X);

      fails += run_all_tests(X, Y_nv, Z_nv, length);

      N_VDestroy(Y_nv);
      N_VDestroy(Z_nv);
    }
    Kokkos::kokkos_free<ExecSpace::memory_space>(unmanaged_ptr);
#endif
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
