/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 * Test the Ginkgo SUNLinearSolver implementation
 * ---------------------------------------------------------------------------*/

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <ginkgo/ginkgo.hpp>
#include <map>
#include <memory>
#include <random>
#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_ginkgo.hpp>
#include <sunmatrix/sunmatrix_ginkgo.hpp>
#include <unordered_map>

#include "test_sunlinsol.h"

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
  int num_threads = omp_num_threads_var ? std::atoi(omp_num_threads_var) : 1;
  return N_VNew_OpenMP(length, num_threads, sunctx);
};
#else
#include <nvector/nvector_serial.h>
#define HIP_OR_CUDA_OR_SYCL(a, b, c)
constexpr auto N_VNew = N_VNew_Serial;
#endif

/* -------------------------------------------------------------------------- *
 * Matrix and solver options                                                  *
 * -------------------------------------------------------------------------- */

// "multigrid" does not support the combined stopping criteria we use
// "cbgmres" does not support setting stopping criteria
const std::unordered_map<std::string, int> methods{{"bicg", 0}, {"bicgstab", 1},
                                                   {"cg", 2},   {"cgs", 3},
                                                   {"fcg", 4},  {"gmres", 5},
                                                   {"idr", 6}};

const std::unordered_map<std::string, int> matrix_types{{"csr", 0}, {"dense", 1}};

/* -------------------------------------------------------------------------- *
 * Matrix fill functions                                                      *
 * -------------------------------------------------------------------------- */

#if defined(USE_CUDA) || defined(USE_HIP)
__global__ void fill_kernel(sunindextype mat_rows, sunindextype mat_cols,
                            sunindextype* row_ptrs, sunindextype* col_idxs,
                            sunrealtype* mat_data)
{
  const sunindextype row = blockIdx.x * blockDim.x + threadIdx.x;
  const sunindextype nnz = 3 * mat_rows - 2;

  if (row == 0)
  {
    // first row
    mat_data[0] = 2;
    mat_data[1] = -1;
    col_idxs[0] = 0;
    col_idxs[1] = 1;
    row_ptrs[0] = 0;
  }
  else if (row == mat_rows - 1)
  {
    // last row
    mat_data[nnz - 2]      = -1;
    mat_data[nnz - 1]      = 2;
    col_idxs[nnz - 2]      = mat_rows - 2;
    col_idxs[nnz - 1]      = mat_rows - 1;
    row_ptrs[mat_rows - 1] = nnz - 2;
    row_ptrs[mat_rows]     = nnz;
  }
  else if (row < mat_rows)
  {
    // other rows
    sunindextype idx  = 3 * row - 1;
    mat_data[idx]     = -1;
    mat_data[idx + 1] = 2;
    mat_data[idx + 2] = -1;
    col_idxs[idx]     = row - 1;
    col_idxs[idx + 1] = row;
    col_idxs[idx + 2] = row + 1;
    row_ptrs[row]     = idx;
  }
}

__global__ void fill_kernel(sunindextype mat_rows, sunindextype mat_cols,
                            sunrealtype* mat_data)
{
  const sunindextype row = blockIdx.x * blockDim.x + threadIdx.x;

  if (row == 0)
  {
    // first row
    mat_data[0] = 2;
    mat_data[1] = -1;
  }
  else if (row == mat_rows - 1)
  {
    // last row
    mat_data[mat_cols * mat_rows - 2] = -1;
    mat_data[mat_cols * mat_rows - 1] = 2;
  }
  else if (row < mat_rows)
  {
    // other rows
    sunindextype idx  = mat_cols * row + row;
    mat_data[idx - 1] = -1;
    mat_data[idx]     = 2;
    mat_data[idx + 1] = -1;
  }
}
#endif

static void fill_matrix(
  std::shared_ptr<gko::matrix::Csr<sunrealtype, sunindextype>> matrix)
{
  sunindextype mat_rows  = static_cast<sunindextype>(matrix->get_size()[0]);
  sunindextype mat_cols  = static_cast<sunindextype>(matrix->get_size()[1]);
  sunindextype* row_ptrs = matrix->get_row_ptrs();
  sunindextype* col_idxs = matrix->get_col_idxs();
  sunrealtype* mat_data  = matrix->get_values();

#if defined(USE_CUDA) || defined(USE_HIP)
  unsigned threads_per_block = 256;
  unsigned num_blocks = (mat_rows + threads_per_block - 1) / threads_per_block;

  fill_kernel<<<num_blocks, threads_per_block>>>(mat_rows, mat_cols, row_ptrs,
                                                 col_idxs, mat_data);
  HIP_OR_CUDA_OR_SYCL(hipDeviceSynchronize(), cudaDeviceSynchronize(), );
#elif defined(USE_SYCL)
  std::dynamic_pointer_cast<const gko::DpcppExecutor>(matrix->get_executor())
    ->get_queue()
    ->submit(
      [&](sycl::handler& cgh)
      {
        cgh.parallel_for(mat_rows,
                         [=](sycl::id<1> id)
                         {
                           const sunindextype row = id[0];
                           // copied from fill_kernel for csr`
                           const sunindextype nnz = 3 * mat_rows - 2;

                           if (row == 0)
                           {
                             // first row
                             mat_data[0] = 2;
                             mat_data[1] = -1;
                             col_idxs[0] = 0;
                             col_idxs[1] = 1;
                             row_ptrs[0] = 0;
                           }
                           else if (row == mat_rows - 1)
                           {
                             // last row
                             mat_data[nnz - 2]      = -1;
                             mat_data[nnz - 1]      = 2;
                             col_idxs[nnz - 2]      = mat_rows - 2;
                             col_idxs[nnz - 1]      = mat_rows - 1;
                             row_ptrs[mat_rows - 1] = nnz - 2;
                             row_ptrs[mat_rows]     = nnz;
                           }
                           else if (row < mat_rows)
                           {
                             // other rows
                             sunindextype idx  = 3 * row - 1;
                             mat_data[idx]     = -1;
                             mat_data[idx + 1] = 2;
                             mat_data[idx + 2] = -1;
                             col_idxs[idx]     = row - 1;
                             col_idxs[idx + 1] = row;
                             col_idxs[idx + 2] = row + 1;
                             row_ptrs[row]     = idx;
                           }
                         });
      });
  matrix->get_executor()->synchronize();
#else
  // Matrix entries
  const sunrealtype vals[] = {-1, 2, -1};

  // Fill matrix
  int idx     = 0;
  row_ptrs[0] = idx;
  for (auto row = 0; row < mat_rows; ++row)
  {
    for (auto diag_offset : {-1, 0, 1})
    {
      auto col = row + diag_offset;
      if (0 <= col && col < mat_cols)
      {
        mat_data[idx] = vals[diag_offset + 1];
        col_idxs[idx] = col;
        ++idx;
      }
    }
    row_ptrs[row + 1] = idx;
  }
#endif
}

static void fill_matrix(std::shared_ptr<gko::matrix::Dense<sunrealtype>> matrix)
{
  sunindextype mat_rows = static_cast<sunindextype>(matrix->get_size()[0]);
  sunindextype mat_cols = static_cast<sunindextype>(matrix->get_size()[1]);
  sunrealtype* mat_data = matrix->get_values();

#if defined(USE_CUDA) || defined(USE_HIP)
  unsigned threads_per_block = 256;
  unsigned num_blocks = (mat_rows + threads_per_block - 1) / threads_per_block;

  fill_kernel<<<num_blocks, threads_per_block>>>(mat_rows, mat_cols, mat_data);
  HIP_OR_CUDA_OR_SYCL(hipDeviceSynchronize(), cudaDeviceSynchronize(), );
#elif defined(USE_SYCL)
  std::dynamic_pointer_cast<const gko::DpcppExecutor>(matrix->get_executor())
    ->get_queue()
    ->submit(
      [&](sycl::handler& cgh)
      {
        cgh.parallel_for(mat_rows,
                         [=](sycl::id<1> id)
                         {
                           const sunindextype row = id[0];
                           // copied from fill_kernel for dense
                           if (row == 0)
                           {
                             // first row
                             mat_data[0] = 2;
                             mat_data[1] = -1;
                           }
                           else if (row == mat_rows - 1)
                           {
                             // last row
                             mat_data[mat_cols * mat_rows - 2] = -1;
                             mat_data[mat_cols * mat_rows - 1] = 2;
                           }
                           else if (row < mat_rows)
                           {
                             // other rows
                             sunindextype idx  = mat_cols * row + row;
                             mat_data[idx - 1] = -1;
                             mat_data[idx]     = 2;
                             mat_data[idx + 1] = -1;
                           }
                         });
      });
  matrix->get_executor()->synchronize();
#else
  // Matrix entries
  const sunrealtype vals[] = {-1, 2, -1};

  // Fill matrix
  for (auto row = 0; row < mat_rows; ++row)
  {
    for (auto diag_offset : {-1, 0, 1})
    {
      auto col = row + diag_offset;
      if (0 <= col && col < mat_cols)
      {
        // Data stored in row-major format
        auto idx      = row * mat_cols + col;
        mat_data[idx] = vals[diag_offset + 1];
      }
    }
  }
#endif
}

/* -------------------------------------------------------------------------- *
 * Extra test functions                                                       *
 * -------------------------------------------------------------------------- */

template<class GkoSolverType, class GkoMatrixType>
void Test_Move(std::unique_ptr<typename GkoSolverType::Factory>&& gko_solver_factory,
               sundials::Context& sunctx)
{
  // Move constructor
  sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType>
    solver{std::move(gko_solver_factory), sunctx};
  sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType> solver2{
    std::move(solver)};
  assert(solver2.GkoFactory());
  assert(solver2.GkoExec());
  assert(SUNLinSolNumIters(solver2) == 0);

  // Move assignment
  sundials::ginkgo::LinearSolver<GkoSolverType, GkoMatrixType> solver3;
  solver3 = std::move(solver2);
  assert(solver3.GkoFactory());
  assert(solver3.GkoExec());
  assert(SUNLinSolNumIters(solver3) == 0);

  std::cout << "    PASSED test -- Test_Move\n";
}

/* -------------------------------------------------------------------------- *
 * Global Executor for sync_device                                            *
 * -------------------------------------------------------------------------- */
// sycl only provides synchronize on queue
std::shared_ptr<const gko::Executor> global_exec;

/* -------------------------------------------------------------------------- *
 * SUNLinSol_Ginkgo Testing Routine                                           *
 * -------------------------------------------------------------------------- */

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

  if (argc < 7)
  {
    std::cerr << "ERROR: SIX (6) inputs required:\n"
              << "  1) method\n"
              << "  2) matrix type\n"
              << "  3) number of matrix columns\n"
              << "  4) condition number\n"
              << "  5) max iterations\n"
              << "  6) print timing\n";
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
            << "  size       = " << matrows << " x " << matcols << "\n"
            << "  cond       = " << matcond << "\n"
            << "  max iters. = " << max_iters << "\n\n";

  /* ------------------------------- *
   * Create solution and RHS vectors *
   * ------------------------------- */

#if defined(USE_SYCL)
  N_Vector x{N_VNew(matcols, gko_exec->get_queue(), sunctx)};
#else
  N_Vector x{N_VNew(matcols, sunctx)};
#endif
  N_Vector b{N_VClone(x)};

  /* Fill x with random data */
  std::default_random_engine engine;
  std::uniform_real_distribution<sunrealtype> distribution_real(8, 10);

  auto xdata{N_VGetArrayPointer(x)};
  for (sunindextype i = 0; i < matcols; i++)
  {
    xdata[i] = distribution_real(engine);
  }
  HIP_OR_CUDA_OR_SYCL(N_VCopyToDevice_Hip(x), N_VCopyToDevice_Cuda(x),
                      N_VCopyToDevice_Sycl(x));

  /* -------------------- *
   * Create system matrix *
   * -------------------- */

  /* Wrap ginkgo matrix for SUNDIALS,
     Matrix is overloaded to a SUNMatrix. */

  std::unique_ptr<sundials::ConvertibleTo<SUNMatrix>> A;

  auto matrix_dim{gko::dim<2>(matrows, matcols)};

  if (matrix_type == "csr")
  {
    using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
    auto matrix_nnz{3 * matrows - 2};
    auto gko_matrix =
      gko::share(GkoMatrixType::create(gko_exec, matrix_dim, matrix_nnz));

    if (matcond > 0)
    {
      auto gko_matdata{gko::matrix_data<
        sunrealtype, sunindextype>::cond(matrows,
                                         gko::remove_complex<sunrealtype>{matcond},
                                         distribution_real, engine)};
      gko_matdata.remove_zeros();
      gko_matrix->read(gko_matdata);
    }
    else { fill_matrix(gko_matrix); }
    A = std::make_unique<sundials::ginkgo::Matrix<GkoMatrixType>>(std::move(
                                                                    gko_matrix),
                                                                  sunctx);
  }
  else if (matrix_type == "dense")
  {
    using GkoMatrixType = gko::matrix::Dense<sunrealtype>;
    auto gko_matrix = gko::share(GkoMatrixType::create(gko_exec, matrix_dim));
    if (matcond > 0)
    {
      auto gko_matdata{gko::matrix_data<
        sunrealtype, sunindextype>::cond(matrows,
                                         gko::remove_complex<sunrealtype>{matcond},
                                         distribution_real, engine)};
      gko_matdata.remove_zeros();
      gko_matrix->read(gko_matdata);
    }
    else
    {
      gko_matrix->fill(0.0);
      fill_matrix(gko_matrix);
    }
    A = std::make_unique<sundials::ginkgo::Matrix<GkoMatrixType>>(std::move(
                                                                    gko_matrix),
                                                                  sunctx);
  }

  /* Create right-hand side vector for linear solve */
  fails += SUNMatMatvecSetup(A->get());
  fails += SUNMatMatvec(A->get(), x, b);
  if (fails)
  {
    std::cerr << "FAIL: SUNLinSol SUNMatMatvec failure\n";
    N_VDestroy(x);
    N_VDestroy(b);
    return 1;
  }

  /* -------------------- *
   * Create linear solver *
   * -------------------- */

  /* Use default stopping criteria */
  auto crit{sundials::ginkgo::DefaultStop::build().with_max_iters(max_iters).on(
    gko_exec)};

  /* Use a Jacobi preconditioner */
  auto precon{
    gko::preconditioner::Jacobi<sunrealtype, sunindextype>::build().on(gko_exec)};

  /* Wrap ginkgo matrix for SUNDIALS,
     Matrix is overloaded to a SUNLinearSolver. */
  std::unique_ptr<sundials::ConvertibleTo<SUNLinearSolver>> LS;

  if (method == "bicg")
  {
    using GkoSolverType = gko::solver::Bicg<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                              .with_criteria(std::move(crit))
                              .with_preconditioner(std::move(precon))
                              .on(gko_exec)};
    if (matrix_type == "csr")
    {
      using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(),
                                              sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<
        GkoSolverType, GkoMatrixType>>(std::move(gko_solver_factory), sunctx);
    }
    else if (matrix_type == "dense")
    {
      using GkoMatrixType = gko::matrix::Dense<sunrealtype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(),
                                              sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<
        GkoSolverType, GkoMatrixType>>(std::move(gko_solver_factory), sunctx);
    }
  }
  else if (method == "bicgstab")
  {
    using GkoSolverType = gko::solver::Bicgstab<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                              .with_criteria(std::move(crit))
                              .with_preconditioner(std::move(precon))
                              .on(gko_exec)};
    if (matrix_type == "csr")
    {
      using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(),
                                              sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<
        GkoSolverType, GkoMatrixType>>(std::move(gko_solver_factory), sunctx);
    }
    else if (matrix_type == "dense")
    {
      using GkoMatrixType = gko::matrix::Dense<sunrealtype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(),
                                              sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<
        GkoSolverType, GkoMatrixType>>(std::move(gko_solver_factory), sunctx);
    }
  }
  else if (method == "cg")
  {
    using GkoSolverType = gko::solver::Cg<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                              .with_criteria(std::move(crit))
                              .with_preconditioner(std::move(precon))
                              .on(gko_exec)};
    if (matrix_type == "csr")
    {
      using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(),
                                              sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<
        GkoSolverType, GkoMatrixType>>(std::move(gko_solver_factory), sunctx);
    }
    else if (matrix_type == "dense")
    {
      using GkoMatrixType = gko::matrix::Dense<sunrealtype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(),
                                              sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<
        GkoSolverType, GkoMatrixType>>(std::move(gko_solver_factory), sunctx);
    }
  }
  else if (method == "cgs")
  {
    using GkoSolverType = gko::solver::Cgs<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                              .with_criteria(std::move(crit))
                              .with_preconditioner(std::move(precon))
                              .on(gko_exec)};
    if (matrix_type == "csr")
    {
      using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(),
                                              sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<
        GkoSolverType, GkoMatrixType>>(std::move(gko_solver_factory), sunctx);
    }
    else if (matrix_type == "dense")
    {
      using GkoMatrixType = gko::matrix::Dense<sunrealtype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(),
                                              sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<
        GkoSolverType, GkoMatrixType>>(std::move(gko_solver_factory), sunctx);
    }
  }
  else if (method == "fcg")
  {
    using GkoSolverType = gko::solver::Fcg<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                              .with_criteria(std::move(crit))
                              .with_preconditioner(std::move(precon))
                              .on(gko_exec)};
    if (matrix_type == "csr")
    {
      using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(),
                                              sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<
        GkoSolverType, GkoMatrixType>>(std::move(gko_solver_factory), sunctx);
    }
    else if (matrix_type == "dense")
    {
      using GkoMatrixType = gko::matrix::Dense<sunrealtype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(),
                                              sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<
        GkoSolverType, GkoMatrixType>>(std::move(gko_solver_factory), sunctx);
    }
  }
  else if (method == "gmres")
  {
    using GkoSolverType = gko::solver::Gmres<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                              .with_criteria(std::move(crit))
                              .with_preconditioner(std::move(precon))
                              .on(gko_exec)};
    if (matrix_type == "csr")
    {
      using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(),
                                              sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<
        GkoSolverType, GkoMatrixType>>(std::move(gko_solver_factory), sunctx);
    }
    else if (matrix_type == "dense")
    {
      using GkoMatrixType = gko::matrix::Dense<sunrealtype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(),
                                              sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<
        GkoSolverType, GkoMatrixType>>(std::move(gko_solver_factory), sunctx);
    }
  }
  else if (method == "idr")
  {
    using GkoSolverType = gko::solver::Idr<sunrealtype>;
    auto gko_solver_factory{GkoSolverType::build()
                              .with_criteria(std::move(crit))
                              .with_preconditioner(std::move(precon))
                              .on(gko_exec)};
    if (matrix_type == "csr")
    {
      using GkoMatrixType = gko::matrix::Csr<sunrealtype, sunindextype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(),
                                              sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<
        GkoSolverType, GkoMatrixType>>(std::move(gko_solver_factory), sunctx);
    }
    else if (matrix_type == "dense")
    {
      using GkoMatrixType = gko::matrix::Dense<sunrealtype>;
      Test_Move<GkoSolverType, GkoMatrixType>(gko_solver_factory->clone(),
                                              sunctx);
      LS = std::make_unique<sundials::ginkgo::LinearSolver<
        GkoSolverType, GkoMatrixType>>(std::move(gko_solver_factory), sunctx);
    }
  }

  /* Run Tests */
  fails += Test_SUNLinSolGetID(LS->get(), SUNLINEARSOLVER_GINKGO, 0);
  fails += Test_SUNLinSolGetType(LS->get(), SUNLINEARSOLVER_MATRIX_ITERATIVE, 0);
  fails += Test_SUNLinSolInitialize(LS->get(), 0);
  fails += Test_SUNLinSolSetup(LS->get(), A->get(), 0);
  fails += Test_SUNLinSolSolve(LS->get(), A->get(), x, b,
                               1e4 * SUN_UNIT_ROUNDOFF, SUNTRUE, 0);

  /* Print result */
  if (fails)
  {
    std::cerr << "FAIL: SUNLinSol module failed " << fails << " tests\n\n";
  }
  else { std::cout << "\nSUCCESS: SUNLinSol module passed all tests\n\n"; }

  /* Print solve information */
  std::cout << "Number of linear solver iterations: "
            << static_cast<long int>(SUNLinSolNumIters(LS->get())) << std::endl;
  std::cout << "Final residual norm: " << SUNLinSolResNorm(LS->get())
            << std::endl;

  // clear global_exec
  global_exec = nullptr;

  /* Free solver, matrix and vectors */
  N_VDestroy(x);
  N_VDestroy(b);

  return fails;
}

/* -------------------------------------------------------------------------- *
 * Implementation-specific 'check' routines                                   *
 * -------------------------------------------------------------------------- */

extern "C" int check_vector(N_Vector expected, N_Vector actual,
                            sunrealtype check_tol)
{
  int failure{0};

  /* copy vectors to host */
  HIP_OR_CUDA_OR_SYCL(N_VCopyFromDevice_Hip(actual),
                      N_VCopyFromDevice_Cuda(actual),
                      N_VCopyFromDevice_Sycl(actual));
  HIP_OR_CUDA_OR_SYCL(N_VCopyFromDevice_Hip(expected),
                      N_VCopyFromDevice_Cuda(expected),
                      N_VCopyFromDevice_Sycl(expected));

  /* get vector data */
  auto xdata{N_VGetArrayPointer(actual)};
  auto ydata{N_VGetArrayPointer(expected)};

  /* check data lengths */
  auto xldata{N_VGetLength(actual)};
  auto yldata{N_VGetLength(expected)};

  if (xldata != yldata)
  {
    std::cerr << ">>> ERROR: check_vector: Different data array lengths\n";
    return 1;
  }

  /* check vector data */
  for (sunindextype i = 0; i < xldata; i++)
  {
    failure += SUNRCompareTol(xdata[i], ydata[i], check_tol);
  }

  if (failure > ZERO)
  {
    std::cerr << "check_vector failures:\n";
    for (sunindextype i = 0; i < xldata; i++)
    {
      if (SUNRCompareTol(xdata[i], ydata[i], check_tol) != 0)
      {
        std::cerr << "  x[" << i << "] = " << xdata[i] << " != " << ydata[i]
                  << " (err = " << abs(xdata[i] - ydata[i]) << ")\n";
      }
    }
  }

  return failure > 0;
}

extern "C" void sync_device()
{
  HIP_OR_CUDA_OR_SYCL(hipDeviceSynchronize(), cudaDeviceSynchronize(),
                      global_exec->synchronize());
}
