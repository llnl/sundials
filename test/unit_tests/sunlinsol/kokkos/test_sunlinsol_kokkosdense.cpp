/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * This is the testing routine for the dense SUNLinearSolver using Kokkos.
 * ---------------------------------------------------------------------------*/

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <nvector/nvector_kokkos.hpp>
#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_kokkosdense.hpp>
#include <sunmatrix/sunmatrix_kokkosdense.hpp>

#include "test_sunlinsol.h"

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
using MatType  = sundials::kokkos::DenseMatrix<ExecSpace>;
using LSType   = sundials::kokkos::DenseLinearSolver<ExecSpace>;
using SizeType = VecType::size_type;

/* -----------------------------------------------------------------------------
 * SUNLinearSolver Testing
 * ---------------------------------------------------------------------------*/

int main(int argc, char* argv[])
{
  // Create SUNDIALS context
  sundials::Context sunctx;

  // Counter for test failures
  int fails = 0;

  // Check input and set matrix dimensions
  if (argc < 4)
  {
    std::cerr << "ERROR: THREE (3) Inputs required:\n"
              << " (1) matrix cols\n"
              << " (2) number of blocks\n"
              << " (3) print timing\n";
    return -1;
  }

  SizeType cols{static_cast<SizeType>(atol(argv[1]))};
  if (cols <= 0)
  {
    std::cerr << "ERROR: number of matrix columns must be a positive integer\n";
    return -1;
  }
  SizeType rows{cols};

  SizeType nblocks{static_cast<SizeType>(atol(argv[2]))};
  if (nblocks <= 0)
  {
    std::cerr << "ERROR: number of blocks must be a positive integer\n";
    return -1;
  }

  int print_timing = atoi(argv[3]);
  SetTiming(print_timing);

  std::cout << "\nKokkos dense linear solver test:\n"
            << "  Size:   " << rows << " by " << cols << "\n"
            << "  Blocks: " << nblocks << "\n";

  Kokkos::initialize(argc, argv);
  {
    auto exec_instance = ExecSpace();

    // Create matrices and vectors
    MatType A{nblocks, rows, cols, exec_instance, sunctx};
    VecType x{cols * nblocks, sunctx};
    VecType b{cols * nblocks, sunctx};

    using RandPoolType = Kokkos::Random_XorShift64_Pool<ExecSpace>;
    using GenType      = RandPoolType::generator_type;

    RandPoolType rand_pool(5374857);

    // Fill A matrix with uniform random data 1 + rand[0,1]
    auto A_data = A.View();

    Kokkos::parallel_for(
      "fill_A",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>({0, 0, 0},
                                                        {nblocks, rows, cols}),
      KOKKOS_LAMBDA(const sunindextype i, const sunindextype j,
                    const sunindextype k) {
        auto rgen       = rand_pool.get_state();
        auto rval       = Kokkos::rand<GenType, sunrealtype>::draw(rgen, ONE);
        A_data(i, j, k) = ONE + rval;
        rand_pool.free_state(rgen);
      });

    // Fill x vector with uniform random data in 1 + rand[0,1]
    auto x_data = x.View();

    Kokkos::parallel_for(
      "fill_x", Kokkos::RangePolicy<ExecSpace>(0, cols * nblocks),
      KOKKOS_LAMBDA(const sunindextype j) {
        auto rgen = rand_pool.get_state();
        auto rval = Kokkos::rand<GenType, sunrealtype>::draw(rgen, ONE);
        x_data(j) = ONE + rval;
        rand_pool.free_state(rgen);
      });

    // Create right-hand side vector for linear solve
    fails += SUNMatMatvecSetup(A);
    fails += SUNMatMatvec(A, x, b);
    if (fails)
    {
      std::cerr << "FAIL: SUNLinSol SUNMatMatvec failure\n";
      return 1;
    }

    // Create dense linear solver
    LSType LS{sunctx};

    // Run Tests
    fails += Test_SUNLinSolInitialize(LS, 0);
    fails += Test_SUNLinSolSetup(LS, A, 0);
    fails += Test_SUNLinSolSolve(LS, A, x, b, SUN_RCONST(1.0e-10), SUNTRUE, 0);
    fails += Test_SUNLinSolGetType(LS, SUNLINEARSOLVER_DIRECT, 0);
    fails += Test_SUNLinSolGetID(LS, SUNLINEARSOLVER_KOKKOSDENSE, 0);

    // Print result
    if (fails)
      std::cout << "FAIL: SUNLinSol module failed " << fails << " tests\n\n";
    else std::cout << "SUCCESS: SUNLinSol module passed all tests\n\n";
  }
  Kokkos::finalize();

  return fails;
}

/* -----------------------------------------------------------------------------
 * Check functions
 * ---------------------------------------------------------------------------*/

KOKKOS_FUNCTION
static int CompareTol(sunrealtype a, sunrealtype b, sunrealtype tol)
{
  if (a == b) { return 0; }
  if (std::isnan(a) || std::isnan(b)) { return 1; }
  if (std::isinf(a) || std::isinf(b)) { return 1; }
  sunrealtype diff = std::abs(a - b);
  sunrealtype norm = std::min(std::abs(a + b),
                              std::numeric_limits<sunrealtype>::max());
  return diff >=
         std::max(10 * std::numeric_limits<sunrealtype>::epsilon(), tol * norm);
}

int check_vector(N_Vector expected, N_Vector computed, sunrealtype tol)
{
  int failure = 0;

  auto e_vec = sundials::kokkos::GetVec<VecType>(expected);
  auto c_vec = sundials::kokkos::GetVec<VecType>(computed);

  auto e_data = e_vec->View();
  auto c_data = c_vec->View();
  auto length = e_vec->Length();

  // check vector data
  Kokkos::parallel_reduce(
    "check_vector", Kokkos::RangePolicy<ExecSpace>(0, length),
    KOKKOS_LAMBDA(const sunindextype i, int& l_fail) {
      l_fail += CompareTol(c_data(i), e_data(i), tol);
    },
    failure);

  if (failure > ZERO) { return 1; }
  else { return 0; }
}

void sync_device() { Kokkos::fence(); }
