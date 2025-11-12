/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * -----------------------------------------------------------------
 * This is the testing routine to check the SUNLinSol Dense module
 * implementation.
 * ----------------------------------------------------------------- */

#include <algorithm>
#include <cctype>
#include <map>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_ginkgobatch.hpp>

#include "test_sunlinsol.h"

#if (GKO_VERSION_MAJOR == 1) && (GKO_VERSION_MINOR < 9)
#error GINKGO >= 1.9.0 required
#endif

#if defined(USE_HIP)
#include <nvector/nvector_hip.h>
#define HIP_OR_CUDA_OR_SYCL(a, b, c) a
constexpr auto N_VNew = N_VNew_Hip;
#elif defined(USE_CUDA)
#include <nvector/nvector_cuda.h>
#define HIP_OR_CUDA_OR_SYCL(a, b, c) b
constexpr auto N_VNew = N_VNew_Cuda;
#elif defined(USE_SYCL)
#include <nvector/nvector_sycl.h>
#define HIP_OR_CUDA_OR_SYCL(a, b, c) c
constexpr auto N_VNew = N_VNew_Sycl;
#elif defined(USE_OMP)
#include <nvector/nvector_openmp.h>
#define HIP_OR_CUDA_OR_SYCL(a, b, c)
auto N_VNew = [](sunindextype length, SUNContext sunctx)
{
  auto omp_num_threads_var{std::getenv("OMP_NUM_THREADS")};
  int num_threads{1};
  if (omp_num_threads_var) { num_threads = std::atoi(omp_num_threads_var); }
  return N_VNew_OpenMP(length, num_threads, sunctx);
};
#else
#include <nvector/nvector_serial.h>
#define HIP_OR_CUDA_OR_SYCL(a, b, c)
constexpr auto N_VNew = N_VNew_Serial;
#endif

using namespace sundials::ginkgo;

const std::unordered_map<std::string, int> methods{{"bicgstab", 0}};
const std::unordered_map<std::string, int> matrix_types{{"csr", 0}, {"dense", 1}};

constexpr sunrealtype solve_tolerance =
  1000 * std::numeric_limits<sunrealtype>::epsilon();

static void fill_matrix_data(gko::matrix_data<sunrealtype, sunindextype>& data)
{
  auto num_rows = data.size[0];
  for (gko::size_type row = 0; row < num_rows; ++row)
  {
    if (row > 0)
    {
      data.nonzeros.emplace_back(row - 1, row, sunrealtype{-1.0});
    }
    data.nonzeros.emplace_back(row, row, sunrealtype{5.0});
    if (row < num_rows - 1)
    {
      data.nonzeros.emplace_back(row, row + 1, sunrealtype{-1.0});
    }
  }
}

/* -------------------------------------------------------------------------- *
 * Global Executor for sync_device                                            *
 * -------------------------------------------------------------------------- */
// sycl only provides synchronize on queue
std::shared_ptr<const gko::Executor> global_exec;

/* ----------------------------------------------------------------------
  * SUNLinSol_Ginkgo Testing Routine
  * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  int argi{0};
  int fails{0}; /* counter for test failures */

  sundials::Context sunctx;

#if defined(USE_HIP)
  auto gko_exec{gko::HipExecutor::create(0, gko::OmpExecutor::create())};
#elif defined(USE_CUDA)
  auto gko_exec{gko::CudaExecutor::create(0, gko::OmpExecutor::create())};
#elif defined(USE_SYCL)
  auto gko_exec{gko::DpcppExecutor::create(0, gko::ReferenceExecutor::create())};
#elif defined(USE_OMP)
  auto gko_exec{gko::OmpExecutor::create()};
#else
  auto gko_exec{gko::ReferenceExecutor::create()};
#endif

  // For sync_device
  global_exec = gko_exec;

  /* ------------ *
   * Check inputs *
   * ------------ */

  if (argc < 8)
  {
    std::cerr << "ERROR: SEVEN (7) inputs required:\n"
              << "  1) batch method\n"
              << "  2) matrix type (csr or dense)\n"
              << "  3) number of matrix columns in a batch\n"
              << "  4) number of matrix batches\n"
              << "  5) condition number of a batch\n"
              << "  6) max iterations per batch\n"
              << "  7) print timing\n";
    return 1;
  }

  std::string method{argv[++argi]};
  std::transform(method.begin(), method.end(), method.begin(),
                 [](unsigned char c) { return std::tolower(c); });

  if (!methods.count(method))
  {
    std::cerr << "ERROR: method must be one of ";
    for (const auto& m : methods) { std::cout << m.first << ", "; }
    std::cout << std::endl;
    return 1;
  }

  std::string matrix_type{argv[++argi]};
  std::transform(matrix_type.begin(), matrix_type.end(), matrix_type.begin(),
                 [](unsigned char c) { return std::tolower(c); });

#ifdef SUNDIALS_INT64_T
  if (matrix_type == "csr")
  {
    std::cerr << "ERROR: the CSR matrix type is not compatible with 64-bit "
                 "index types\n";
    return 1;
  }
#endif

  if (!matrix_types.count(matrix_type))
  {
    std::cerr << "ERROR: matrix type must be one of ";
    for (const auto& m : matrix_types) { std::cout << m.first << ", "; }
    std::cout << std::endl;
    return 1;
  }

  const auto matcols{static_cast<sunindextype>(atoll(argv[++argi]))};
  if (matcols <= 0)
  {
    std::cerr << "ERROR: number of matrix columns must be a positive integer\n";
    return 1;
  }
  const auto matrows{matcols};

  const auto num_batches{static_cast<sunindextype>(atoll(argv[++argi]))};
  if (num_batches <= 0)
  {
    std::cerr << "ERROR: number of batches must be a positive integer\n";
    return 1;
  }

  const auto matcond{SUNStrToReal(argv[++argi])};
  if (matcond < 0)
  {
    std::cerr << "ERROR: matrix condition number must be positive or 0 "
                 "(poisson test)\n";
    return 1;
  }

  const auto max_iters{static_cast<unsigned long>(atoll(argv[++argi]))};
  if (max_iters <= 0)
  {
    std::cerr << "ERROR: max iterations must be a positive integer\n";
    return 1;
  }

  const int print_timing{atoi(argv[++argi])};
  SetTiming(print_timing);

  std::cout << "Ginkgo linear solver test:\n"
            << "  method     = " << method.c_str() << "\n"
            << "  matrix     = " << matrix_type.c_str() << "\n"
            << "  size       = " << matrows << " x " << matcols << " x "
            << num_batches << " (batches)\n"
            << "  cond       = " << matcond << "\n"
            << "  max iters. = " << max_iters << "\n\n";

  /* ------------------------------- *
   * Create solution and RHS vectors *
   * ------------------------------- */

#if defined(USE_SYCL)
  N_Vector x{N_VNew(num_batches * matcols, gko_exec->get_queue(), sunctx)};
#else
  N_Vector x{N_VNew(num_batches * matcols, sunctx)};
#endif
  N_Vector b{N_VClone(x)};

  /* Fill x with random data */
  std::default_random_engine engine;
  std::uniform_real_distribution<sunrealtype> distribution_real(8, 10);
  sunrealtype* xdata = N_VGetArrayPointer(x);
  for (sunindextype batchi = 0; batchi < num_batches; batchi++)
  {
    for (sunindextype i = 0; i < matcols; i++)
    {
      xdata[batchi * matcols + i] = distribution_real(engine);
    }
  }
  HIP_OR_CUDA_OR_SYCL(N_VCopyToDevice_Hip(x), N_VCopyToDevice_Cuda(x),
                      N_VCopyToDevice_Sycl(x));

  auto matrix_dim = gko::dim<2>(matrows, matcols);

  auto gko_matdata =
    matcond > 0
      ? gko::matrix_data<sunrealtype,
                         sunindextype>::cond(matrows,
                                             gko::remove_complex<sunrealtype>{
                                               matcond},
                                             distribution_real, engine)
      : gko::matrix_data<sunrealtype, sunindextype>(matrix_dim);

  if (matcond <= 0) { fill_matrix_data(gko_matdata); }

  std::unique_ptr<sundials::ConvertibleTo<SUNMatrix>> A, B;

  if (matrix_type == "csr")
  {
#ifdef SUNDIALS_INT32_T
    using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
    using GkoBatchMatrixType = gko::batch::matrix::Csr<sunrealtype, sunindextype>;

    gko_matdata.remove_zeros();

    auto gko_matrix{GkoMatrixType::create(gko_exec, matrix_dim)};
    gko_matrix->read(gko_matdata);

    auto num_nnz     = gko_matrix->get_num_stored_elements();
    auto common_size = gko_matrix->get_size();
    auto gko_batch_matrix{
      GkoBatchMatrixType::create(gko_exec,
                                 gko::batch_dim<2>(num_batches, common_size),
                                 num_nnz)};
    for (sunindextype blk = 0; blk < num_batches; ++blk)
    {
      gko_batch_matrix->create_view_for_item(blk)->read(gko_matdata);
    }

    A = std::make_unique<
      sundials::ginkgo::BatchMatrix<GkoBatchMatrixType>>(std::move(
                                                           gko_batch_matrix),
                                                         sunctx);
#else
    return 1;
#endif
  }
  else if (matrix_type == "dense")
  {
    using GkoMatrixType      = gko::matrix::Dense<sunrealtype>;
    using GkoBatchMatrixType = gko::batch::matrix::Dense<sunrealtype>;

    auto gko_matrix{GkoMatrixType::create(gko_exec, matrix_dim)};
    gko_matrix->read(gko_matdata);

    auto common_size = gko_matrix->get_size();
    auto gko_batch_matrix{
      GkoBatchMatrixType::create(gko_exec,
                                 gko::batch_dim<2>(num_batches, common_size))};
    for (sunindextype blk = 0; blk < num_batches; ++blk)
    {
      gko_batch_matrix->create_view_for_item(blk)->read(gko_matdata);
    }

    A = std::make_unique<
      sundials::ginkgo::BatchMatrix<GkoBatchMatrixType>>(std::move(
                                                           gko_batch_matrix),
                                                         sunctx);
  }

  /* Create scaling vectors */
  N_Vector s1{N_VClone(x)};
  N_Vector s2{N_VClone(x)};
  N_VConst(sunrealtype{2.0}, s1);
  N_VConst(sunrealtype{3.0}, s2);

  /* Create right-hand side vector for linear solve */
  fails += SUNMatMatvecSetup(A->Convert());
  fails += SUNMatMatvec(A->Convert(), x, b);
  if (fails)
  {
    printf("FAIL: SUNLinSol SUNMatMatvec failure\n");

    N_VDestroy(x);
    N_VDestroy(b);
    N_VDestroy(s1);
    N_VDestroy(s2);

    return (1);
  }

  /*
   * Create linear solver.
   */

  auto precond_factory = gko::share(
    gko::batch::preconditioner::Jacobi<sunrealtype>::build().on(gko_exec));

  std::unique_ptr<sundials::ConvertibleTo<SUNLinearSolver>> LS;
  if (method == "bicgstab")
  {
    using GkoSolverType = gko::batch::solver::Bicgstab<sunrealtype>;
    if (matrix_type == "csr")
    {
#ifdef SUNDIALS_INT32_T
      using GkoBatchMatrixType = gko::batch::matrix::Csr<sunrealtype>;
      using SUNGkoLinearSolverType =
        BatchLinearSolver<GkoSolverType, GkoBatchMatrixType>;
      LS =
        std::make_unique<SUNGkoLinearSolverType>(gko_exec,
                                                 gko::batch::stop::tolerance_type::absolute,
                                                 precond_factory, max_iters,
                                                 num_batches, sunctx);
#else
      return 1;
#endif
    }
    else if (matrix_type == "dense")
    {
      using GkoBatchMatrixType = gko::batch::matrix::Dense<sunrealtype>;
      using SUNGkoLinearSolverType =
        BatchLinearSolver<GkoSolverType, GkoBatchMatrixType>;
      LS =
        std::make_unique<SUNGkoLinearSolverType>(gko_exec,
                                                 gko::batch::stop::tolerance_type::absolute,
                                                 nullptr, max_iters,
                                                 num_batches, sunctx);
    }
  }

  /* Run Tests */
  fails += Test_SUNLinSolGetID(LS->Convert(), SUNLINEARSOLVER_GINKGOBATCH, 0);
  fails += Test_SUNLinSolGetType(LS->Convert(),
                                 SUNLINEARSOLVER_MATRIX_ITERATIVE, 0);
  fails += Test_SUNLinSolInitialize(LS->Convert(), 0);
  fails += Test_SUNLinSolSetup(LS->Convert(), A->Convert(), 0);
  fails += Test_SUNLinSolSolve(LS->Convert(), A->Convert(), x, b,
                               solve_tolerance, SUNTRUE, 0);

  /* Try with scaling */
  fails += Test_SUNLinSolSetScalingVectors(LS->Convert(), s1, s2, 0);
  fails += Test_SUNLinSolSetup(LS->Convert(), A->Convert(), 0);
  fails += Test_SUNLinSolSolve(LS->Convert(), A->Convert(), x, b,
                               solve_tolerance, SUNTRUE, 0);

  /* Print result */
  if (fails)
  {
    std::cerr << "FAIL: SUNLinSol module failed " << fails << " tests\n\n";
  }
  else { std::cout << "\nSUCCESS: SUNLinSol module passed all tests\n\n"; }

  // clear global_exec
  global_exec = nullptr;

  /* Free solver, matrix and vectors */
  N_VDestroy(x);
  N_VDestroy(b);
  N_VDestroy(s1);
  N_VDestroy(s2);

  return (fails);
}

/* ----------------------------------------------------------------------
  * Implementation-specific 'check' routines
  * --------------------------------------------------------------------*/

extern "C" int check_vector(N_Vector expected, N_Vector actual, sunrealtype tol)
{
  int failure = 0;
  sunrealtype *xdata, *ydata;
  sunindextype xldata, yldata;
  sunindextype i;

  /* copy vectors to host */
  HIP_OR_CUDA_OR_SYCL(N_VCopyFromDevice_Hip(actual),
                      N_VCopyFromDevice_Cuda(actual),
                      N_VCopyFromDevice_Sycl(actual));
  HIP_OR_CUDA_OR_SYCL(N_VCopyFromDevice_Hip(expected),
                      N_VCopyFromDevice_Cuda(expected),
                      N_VCopyFromDevice_Sycl(expected));

  /* get vector data */
  xdata = N_VGetArrayPointer(actual);
  ydata = N_VGetArrayPointer(expected);

  /* check data lengths */
  xldata = N_VGetLength(actual);
  yldata = N_VGetLength(expected);

  if (xldata != yldata)
  {
    printf(">>> ERROR: check_vector: Different data array lengths \n");
    return (1);
  }

  /* check vector data */
  for (i = 0; i < xldata; i++)
    failure += SUNRCompareTol(xdata[i], ydata[i], tol);

  if (failure > ZERO)
  {
    printf("Check_vector failures:\n");
    for (i = 0; i < xldata; i++)
      if (SUNRCompareTol(xdata[i], ydata[i], tol) != 0)
        printf("  actual[%ld] = %g != %e (err = %g)\n", (long int)i, xdata[i],
               ydata[i], SUNRabs(xdata[i] - ydata[i]));
  }

  return failure > 0;
}

extern "C" void sync_device(void)
{
  HIP_OR_CUDA_OR_SYCL(hipDeviceSynchronize(), cudaDeviceSynchronize(),
                      global_exec->synchronize());
}
