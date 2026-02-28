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
 * This is the testing routine for the dense SUNMatrix using Kokkos.
 * ---------------------------------------------------------------------------*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <nvector/nvector_kokkos.hpp>
#include <sundials/sundials_core.hpp>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_kokkosdense.hpp>

#include "test_sunmatrix.h"

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
using SizeType = VecType::size_type;

/* -----------------------------------------------------------------------------
 * SUNMatrix Testing
 * ---------------------------------------------------------------------------*/

int main(int argc, char* argv[])
{
  // Create SUNDIALS context
  sundials::Context sunctx;

  // counter for test failures
  int fails = 0;

  // check input and set vector length
  if (argc < 5)
  {
    std::cerr << "ERROR: FOUR (4) Input required:\n"
              << " (1) matrix rows\n"
              << " (2) matrix cols\n"
              << " (3) number of matrix blocks\n"
              << " (4) print timing\n";
    return -1;
  }

  SizeType matrows{static_cast<SizeType>(atol(argv[1]))};
  if (matrows <= 0)
  {
    std::cerr << "ERROR: number of rows must be a positive integer\n";
    return -1;
  }

  SizeType matcols{static_cast<SizeType>(atol(argv[2]))};
  if (matcols <= 0)
  {
    std::cerr << "ERROR: number of cols must be a positive integer\n";
    return -1;
  }

  SizeType nblocks{static_cast<SizeType>(atol(argv[3]))};
  if (nblocks <= 0)
  {
    std::cerr << "ERROR: number of blocks must be a positive integer\n";
    return -1;
  }

  int print_timing = atoi(argv[4]);
  SetTiming(print_timing);

  bool square = (matrows == matcols);
  std::cout << "\nKokkos dense matrix test:\n"
            << "  Size:   " << matrows << " by " << matcols << "\n"
            << "  Blocks: " << nblocks << "\n";

  Kokkos::initialize(argc, argv);
  {
    auto exec_instance = ExecSpace();

    // Create vectors and matrices
    VecType x{matcols * nblocks, sunctx};
    VecType y{matrows * nblocks, sunctx};
    MatType A{nblocks, matrows, matcols, exec_instance, sunctx};
    MatType I{nblocks, matrows, matcols, exec_instance, sunctx};

    // Fill matrices and vectors
    auto A_data = A.View();

    Kokkos::parallel_for(
      "fill_A",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>({0, 0, 0},
                                                        {nblocks, matrows,
                                                         matcols}),
      KOKKOS_LAMBDA(const SizeType i, const SizeType j, const SizeType k) {
        A_data(i, j, k) = static_cast<sunrealtype>((k + 1) * (j + k));
      });

    if (square)
    {
      auto I_data = I.View();
      Kokkos::parallel_for(
        "fill_I",
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>({0, 0, 0},
                                                          {nblocks, matrows,
                                                           matcols}),
        KOKKOS_LAMBDA(const SizeType i, const SizeType j, const SizeType k) {
          I_data(i, j, k) = (j == k) ? ONE : ZERO;
        });
    }

    auto x_data = x.View();

    Kokkos::parallel_for(
      "fill_x", Kokkos::RangePolicy<ExecSpace>(0, matcols * nblocks),
      KOKKOS_LAMBDA(const SizeType j) {
        x_data(j) = ONE / static_cast<sunrealtype>((j % matcols) + 1);
      });

    auto y_data = y.View();

    Kokkos::parallel_for(
      "fill_y", Kokkos::RangePolicy<ExecSpace>(0, matrows * nblocks),
      KOKKOS_LAMBDA(const SizeType j) {
        auto m    = j % matrows;
        auto n    = m + matcols - 1;
        y_data(j) = HALF * static_cast<sunrealtype>((n + 1 - m) * (n + m));
      });

    // SUNMatrix Tests
    fails += Test_SUNMatGetID(A, SUNMATRIX_KOKKOSDENSE, 0);
    fails += Test_SUNMatClone(A, 0);
    fails += Test_SUNMatCopy(A, 0);
    fails += Test_SUNMatZero(A, 0);
    if (square)
    {
      fails += Test_SUNMatScaleAdd(A, I, 0);
      fails += Test_SUNMatScaleAddI(A, I, 0);
    }
    fails += Test_SUNMatMatvec(A, x, y, 0);

    // Print result
    if (fails)
      std::cout << "FAIL: SUNMatrix module failed " << fails << " tests\n\n";
    else std::cout << "SUCCESS: SUNMatrix module passed all tests\n\n";
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

extern "C" int check_matrix(SUNMatrix A, SUNMatrix B, sunrealtype tol)
{
  int failure = 0;

  auto A_mat{sundials::kokkos::GetDenseMat<MatType>(A)};
  auto B_mat{sundials::kokkos::GetDenseMat<MatType>(B)};

  const auto nblocks = A_mat->Blocks();
  const auto matrows = A_mat->BlockRows();
  const auto matcols = A_mat->BlockCols();

  auto A_data{A_mat->View()};
  auto B_data{B_mat->View()};

  // compare data
  Kokkos::parallel_reduce(
    "check_matrix",
    Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>({0, 0, 0},
                                                      {nblocks, matrows, matcols}),
    KOKKOS_LAMBDA(const SizeType i, const SizeType j, const SizeType k,
                  int& l_fail) {
      l_fail += CompareTol(A_data(i, j, k), B_data(i, j, k), tol);
    },
    failure);

  if (failure > ZERO) { return 1; }
  else { return 0; }
}

extern "C" int check_matrix_entry(SUNMatrix A, sunrealtype val, sunrealtype tol)
{
  int failure = 0;

  auto A_mat{sundials::kokkos::GetDenseMat<MatType>(A)};

  const auto nblocks = A_mat->Blocks();
  const auto matrows = A_mat->BlockRows();
  const auto matcols = A_mat->BlockCols();

  auto A_data{A_mat->View()};

  // compare data
  Kokkos::parallel_reduce(
    "check_matrix_entry",
    Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>({0, 0, 0},
                                                      {nblocks, matrows, matcols}),
    KOKKOS_LAMBDA(const SizeType i, const SizeType j, const SizeType k,
                  int& l_fail) {
      l_fail += CompareTol(A_data(i, j, k), val, tol);
    },
    failure);

  if (failure > ZERO) { return 1; }
  else { return 0; }
}

extern "C" int check_vector(N_Vector actual, N_Vector expected, sunrealtype tol)
{
  int failure = 0;

  auto a_vec{sundials::kokkos::GetVec<VecType>(actual)};
  auto e_vec{sundials::kokkos::GetVec<VecType>(expected)};
  auto length = a_vec->Length();

  auto a_data = a_vec->View();
  auto e_data = e_vec->View();

  // check vector data
  Kokkos::parallel_reduce(
    "check_vector", Kokkos::RangePolicy<ExecSpace>(0, length),
    KOKKOS_LAMBDA(const SizeType i, int& l_fail) {
      l_fail += CompareTol(a_data(i), e_data(i), tol);
    },
    failure);

  if (failure > ZERO) { return 1; }
  else { return 0; }
}

extern "C" sunbooleantype has_data(SUNMatrix A) { return SUNTRUE; }

extern "C" sunbooleantype is_square(SUNMatrix A)
{
  auto A_mat{sundials::kokkos::GetDenseMat<MatType>(A)};

  const auto matrows = A_mat->Rows();
  const auto matcols = A_mat->Cols();

  if (matrows == matcols) { return SUNTRUE; }
  else { return SUNFALSE; }
}

extern "C" void sync_device(SUNMatrix A) { Kokkos::fence(); }
