/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * Test the Ginkgo SUNMatrix implementation
 * ---------------------------------------------------------------------------*/

#include <cstring>
#include <iostream>
#include <memory>
#include <sundials/sundials_core.hpp>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_ginkgo.hpp>

#include "test_sunmatrix.h"

#if defined(USE_HIP)
#define REF_OR_OMP_OR_HIP_OR_CUDA_OR_SYCL(a, b, c, d, e) c
#elif defined(USE_CUDA)
#define REF_OR_OMP_OR_HIP_OR_CUDA_OR_SYCL(a, b, c, d, e) d
#elif defined(USE_SYCL)
#define REF_OR_OMP_OR_HIP_OR_CUDA_OR_SYCL(a, b, c, d, e) e
#elif defined(USE_OMP)
#define REF_OR_OMP_OR_HIP_OR_CUDA_OR_SYCL(a, b, c, d, e) b
#else
#define REF_OR_OMP_OR_HIP_OR_CUDA_OR_SYCL(a, b, c, d, e) a
#endif

#if defined(USE_CUDA)
#include <nvector/nvector_cuda.h>
#elif defined(USE_HIP)
#include <nvector/nvector_hip.h>
#elif defined(USE_OMP)
#include <nvector/nvector_openmp.h>
#elif defined(USE_SYCL)
#include <nvector/nvector_sycl.h>
#else
#include <nvector/nvector_serial.h>
#endif

using GkoDenseMat = gko::matrix::Dense<sunrealtype>;
using GkoCsrMat   = gko::matrix::Csr<sunrealtype, sunindextype>;
using GkoVecType  = GkoDenseMat;

bool using_csr_matrix_type   = false;
bool using_dense_matrix_type = false;

template<class GkoMatType>
int Test_CopyAndMove(std::unique_ptr<GkoMatType>&& gko_mat, SUNContext sunctx)
{
  // Copy constructor
  sundials::ginkgo::Matrix<GkoMatType> mat{std::move(gko_mat), sunctx};
  sundials::ginkgo::Matrix<GkoMatType> mat2{mat};
  SUNMatZero(mat2);
  check_matrix_entry(mat2, sunrealtype{0.0}, SUN_UNIT_ROUNDOFF);

  // Move constructor
  sundials::ginkgo::Matrix<GkoMatType> mat3{std::move(mat2)};
  SUNMatZero(mat3);
  check_matrix_entry(mat3, sunrealtype{0.0}, SUN_UNIT_ROUNDOFF);

  // Copy assignment
  sundials::ginkgo::Matrix<GkoMatType> mat4;
  mat4 = mat3;
  SUNMatZero(mat4);
  check_matrix_entry(mat4, sunrealtype{0.0}, SUN_UNIT_ROUNDOFF);

  // Move assignment
  sundials::ginkgo::Matrix<GkoMatType> mat5;
  mat5 = std::move(mat4);
  SUNMatZero(mat5);
  check_matrix_entry(mat5, sunrealtype{0.0}, SUN_UNIT_ROUNDOFF);

  std::cout << "    PASSED test -- Test_CopyAndMove\n";

  return 0;
}

/* ----------------------------------------------------------------------
 * Main SUNMatrix Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  int fails{0}; /* counter for test failures */

  /* Create SUNDIALS context before calling any other SUNDIALS function*/
  sundials::Context sunctx;

#if defined(USE_CUDA)
  auto gko_exec{gko::CudaExecutor::create(0, gko::OmpExecutor::create())};
#elif defined(USE_HIP)
  auto gko_exec{gko::HipExecutor::create(0, gko::OmpExecutor::create())};
#elif defined(USE_SYCL)
  auto gko_exec{gko::DpcppExecutor::create(0, gko::ReferenceExecutor::create())};
#elif defined(USE_OMP)
  auto gko_exec{gko::OmpExecutor::create()};
#else
  auto gko_exec{gko::ReferenceExecutor::create()};
#endif

  /* check input and set vector length */
  if (argc < 4)
  {
    std::cerr << "ERROR: THREE (3) Input required: matrix rows, matrix cols, "
                 "format (0 = csr, 1 = dense)\n";
    return 1;
  }

  int argi{0};

  auto matrows{static_cast<sunindextype>(atol(argv[++argi]))};
  if (matrows <= 0)
  {
    std::cerr << "ERROR: number of rows must be a positive integer \n";
    return 1;
  }

  auto matcols{static_cast<sunindextype>(atol(argv[++argi]))};
  if (matcols <= 0)
  {
    std::cerr << "ERROR: number of cols must be a positive integer \n";
    return 1;
  }

  auto format{static_cast<int>(atoi(argv[++argi]))};
  if (format != 0 && format != 1)
  {
    std::cerr << "ERROR: format must be 0 (csr) or 1 (dense) \n";
    return 1;
  }

  if (format == 0) { using_csr_matrix_type = true; }
  else if (format == 1) { using_dense_matrix_type = true; }

#if defined(USE_OMP)
  int num_threads{1};
  auto omp_num_threads_var{std::getenv("OMP_NUM_THREADS")};
  if (omp_num_threads_var) { num_threads = std::atoi(omp_num_threads_var); }
#endif

  SetTiming(0);

  int square{matrows == matcols ? 1 : 0};
  std::cout << "\n SUNMATRIX_GINKGO test: size " << matrows << " x " << matcols
            << ", format ";
  if (using_csr_matrix_type) { std::cout << "csr\n"; }
  else if (using_dense_matrix_type) { std::cout << "dense\n"; }

  /* Create vectors and matrices */
  std::default_random_engine generator;
  std::uniform_real_distribution<sunrealtype>
    distribution{0.0, static_cast<sunrealtype>(matrows)};

  N_Vector x{
    REF_OR_OMP_OR_HIP_OR_CUDA_OR_SYCL(N_VNew_Serial(matcols, sunctx),
                                      N_VNew_OpenMP(matcols, num_threads, sunctx),
                                      N_VNew_Hip(matcols, sunctx),
                                      N_VNew_Cuda(matcols, sunctx),
                                      N_VNew_Sycl(matcols, gko_exec->get_queue(),
                                                  sunctx))};
  N_Vector y{
    REF_OR_OMP_OR_HIP_OR_CUDA_OR_SYCL(N_VNew_Serial(matrows, sunctx),
                                      N_VNew_OpenMP(matrows, num_threads, sunctx),
                                      N_VNew_Hip(matrows, sunctx),
                                      N_VNew_Cuda(matrows, sunctx),
                                      N_VNew_Sycl(matrows, gko_exec->get_queue(),
                                                  sunctx))};

  auto matrix_dim{gko::dim<2>(matrows, matcols)};
  auto gko_matdata{gko::matrix_data<sunrealtype, sunindextype>(matrix_dim,
                                                               distribution,
                                                               generator)};

  /* Wrap ginkgo matrices for SUNDIALS.
     sundials::ginkgo::Matrix is overloaded to a SUNMatrix. */

  std::unique_ptr<sundials::ConvertibleTo<SUNMatrix>> A;
  std::unique_ptr<sundials::ConvertibleTo<SUNMatrix>> I;

  auto xdata{N_VGetArrayPointer(x)};
  for (sunindextype i = 0; i < matcols; i++)
  {
    xdata[i] = distribution(generator);
  }
  REF_OR_OMP_OR_HIP_OR_CUDA_OR_SYCL(, , N_VCopyToDevice_Hip(x),
                                    N_VCopyToDevice_Cuda(x),
                                    N_VCopyToDevice_Sycl(x));

  /* Compute true solution */
  SUNMatrix Aref{SUNDenseMatrix(matrows, matcols, sunctx)};
  if (using_csr_matrix_type)
  {
    auto gko_matrix{GkoCsrMat::create(gko_exec, matrix_dim)};
    gko_matrix->read(gko_matdata);
    auto gko_ident{GkoCsrMat::create(gko_exec, matrix_dim)};
    if (square)
    {
      gko_ident->read(
        gko::matrix_data<sunrealtype, sunindextype>::diag(matrix_dim, 1.0));
    }

    auto Arowptrs{gko_matrix->get_const_row_ptrs()};
    auto Acolidxs{gko_matrix->get_const_col_idxs()};
    auto Avalues{gko_matrix->get_const_values()};
    for (sunindextype irow = 0;
         irow < static_cast<sunindextype>(gko_matrix->get_size()[0]); irow++)
    {
      for (auto inz = gko_exec->copy_val_to_host(Arowptrs + irow);
           inz < gko_exec->copy_val_to_host(Arowptrs + irow + 1); inz++)
      {
        SM_ELEMENT_D(Aref, irow, gko_exec->copy_val_to_host(Acolidxs + inz)) =
          gko_exec->copy_val_to_host(Avalues + inz);
      }
    }

    fails += Test_CopyAndMove(gko_matrix->clone(), sunctx);

    A = std::make_unique<sundials::ginkgo::Matrix<GkoCsrMat>>(std::move(gko_matrix),
                                                              sunctx);
    I = std::make_unique<sundials::ginkgo::Matrix<GkoCsrMat>>(std::move(gko_ident),
                                                              sunctx);
  }
  else if (using_dense_matrix_type)
  {
    auto gko_matrix{GkoDenseMat::create(gko_exec, matrix_dim)};
    gko_matrix->read(gko_matdata);
    auto gko_ident{GkoDenseMat::create(gko_exec, matrix_dim)};
    if (square)
    {
      gko_ident->read(
        gko::matrix_data<sunrealtype, sunindextype>::diag(matrix_dim, 1.0));
    }

    auto Avalues{gko_matrix->get_const_values()};
    for (sunindextype j = 0; j < matcols; j++)
    {
      for (sunindextype i = 0; i < matrows; i++)
      {
        SM_ELEMENT_D(Aref, i, j) =
          gko_exec->copy_val_to_host(Avalues + i * gko_matrix->get_stride() + j);
      }
    }

    fails += Test_CopyAndMove(gko_matrix->clone(), sunctx);

    A = std::make_unique<sundials::ginkgo::Matrix<GkoDenseMat>>(std::move(
                                                                  gko_matrix),
                                                                sunctx);
    I = std::make_unique<sundials::ginkgo::Matrix<GkoDenseMat>>(std::move(
                                                                  gko_ident),
                                                                sunctx);
  }
  SUNMatMatvec_Dense(Aref, x, y);
  SUNMatDestroy(Aref);
  REF_OR_OMP_OR_HIP_OR_CUDA_OR_SYCL(, , N_VCopyToDevice_Hip(y),
                                    N_VCopyToDevice_Cuda(y),
                                    N_VCopyToDevice_Sycl(y));
  /* SUNMatrix Tests */
  fails += Test_SUNMatGetID(*A, SUNMATRIX_GINKGO, 0);
  fails += Test_SUNMatClone(*A, 0);
  fails += Test_SUNMatCopy(*A, 0);
  fails += Test_SUNMatZero(*A, 0);
  if (square)
  {
#if !defined(USE_OMP)
    // TODO(CJB): ScaleAdd with a dense matrix is not supported in develop
    // branch, possibly supported on the batch-develop branch. CSR matrix with
    // OMP executor currently fails on the develop branch.
    if (!using_dense_matrix_type) { fails += Test_SUNMatScaleAdd(*A, *I, 0); }
#endif
    fails += Test_SUNMatScaleAddI(*A, *I, 0);
  }
  fails += Test_SUNMatMatvec(*A, x, y, 0);

  /* Print result */
  if (fails)
  {
    std::cerr << " FAIL: SUNMatrix module failed " << fails << " tests \n\n";
  }
  else { std::cout << " SUCCESS: SUNMatrix module passed all tests \n\n"; }

  /* Free vectors */
  N_VDestroy(x);
  N_VDestroy(y);

  return fails;
}

/* ----------------------------------------------------------------------
 * Check matrix
 * --------------------------------------------------------------------*/
static int check_matrix_csr(SUNMatrix A, SUNMatrix B, sunrealtype tol)
{
  int failure{0};
  auto Amat{
    static_cast<sundials::ginkgo::Matrix<GkoCsrMat>*>(A->content)->GkoMtx()};
  auto Bmat{
    static_cast<sundials::ginkgo::Matrix<GkoCsrMat>*>(B->content)->GkoMtx()};
  auto Amat_ref = Amat->clone(Amat->get_executor()->get_master());
  auto Bmat_ref = Bmat->clone(Bmat->get_executor()->get_master());
  auto Arowptrs{Amat_ref->get_const_row_ptrs()};
  auto Avalues{Amat_ref->get_const_values()};
  auto Bvalues{Bmat_ref->get_const_values()};

  /* check lengths */
  if (Amat_ref->get_size() != Bmat_ref->get_size())
  {
    std::cerr << ">>> ERROR: check_matrix: Different data array lengths \n";
    return 1;
  }

  /* compare data */
  for (sunindextype irow = 0;
       irow < static_cast<sunindextype>(Amat_ref->get_size()[0]); irow++)
  {
    for (sunindextype inz = Arowptrs[irow]; inz < Arowptrs[irow + 1]; inz++)
    {
      failure += SUNRCompareTol(Avalues[inz], Bvalues[inz], tol);
    }
  }

  return failure > 0;
}

static int check_matrix_dense(SUNMatrix A, SUNMatrix B, sunrealtype tol)
{
  int failure{0};
  auto Amat{
    static_cast<sundials::ginkgo::Matrix<GkoDenseMat>*>(A->content)->GkoMtx()};
  auto Bmat{
    static_cast<sundials::ginkgo::Matrix<GkoDenseMat>*>(B->content)->GkoMtx()};
  auto Amat_ref = Amat->clone(Amat->get_executor()->get_master());
  auto Bmat_ref = Bmat->clone(Bmat->get_executor()->get_master());
  auto rows{Amat->get_size()[0]};
  auto cols{Amat->get_size()[1]};

  /* check lengths */
  if (Amat->get_size() != Bmat->get_size())
  {
    std::cerr << ">>> ERROR: check_matrix: Different data array lengths \n";
    return 1;
  }

  /* compare data */
  for (sunindextype i = 0; i < static_cast<sunindextype>(rows); i++)
  {
    for (sunindextype j = 0; j < static_cast<sunindextype>(cols); j++)
    {
      failure += SUNRCompareTol(Amat_ref->at(i, j), Bmat_ref->at(i, j), tol);
    }
  }

  return failure > 0;
}

extern "C" int check_matrix(SUNMatrix A, SUNMatrix B, sunrealtype tol)
{
  if (using_csr_matrix_type) { return check_matrix_csr(A, B, tol); }
  else if (using_dense_matrix_type) { return check_matrix_dense(A, B, tol); }
  else { return 1; }
}

static int check_matrix_entry_csr(SUNMatrix A, sunrealtype val, sunrealtype tol)
{
  int failure{0};
  auto Amat{
    static_cast<sundials::ginkgo::Matrix<GkoCsrMat>*>(A->content)->GkoMtx()};
  auto Amat_ref = Amat->clone(Amat->get_executor()->get_master());
  auto Arowptrs{Amat_ref->get_const_row_ptrs()};
  auto Avalues{Amat_ref->get_const_values()};

  /* compare data */
  for (sunindextype irow = 0;
       irow < static_cast<sunindextype>(Amat_ref->get_size()[0]); irow++)
  {
    for (sunindextype inz = Arowptrs[irow]; inz < Arowptrs[irow + 1]; inz++)
    {
      int check = SUNRCompareTol(Avalues[inz], val, tol);
      if (check)
      {
        std::cerr << "  actual = " << Avalues[inz] << " != " << val
                  << " (err = " << SUNRabs(Avalues[inz] - val) << ")\n";
        failure += check;
      }
    }
  }

  return failure > 0;
}

static int check_matrix_entry_dense(SUNMatrix A, sunrealtype val, sunrealtype tol)
{
  int failure{0};
  auto Amat{
    static_cast<sundials::ginkgo::Matrix<GkoDenseMat>*>(A->content)->GkoMtx()};
  auto rows{Amat->get_size()[0]};
  auto cols{Amat->get_size()[1]};

  auto Amat_ref = Amat->clone(Amat->get_executor()->get_master());
  /* compare data */
  for (sunindextype i = 0; i < static_cast<sunindextype>(rows); i++)
  {
    for (sunindextype j = 0; j < static_cast<sunindextype>(cols); j++)
    {
      int check = SUNRCompareTol(Amat_ref->at(i, j), val, tol);
      if (check)
      {
        std::cerr << "  actual[" << i << "," << j
                  << "] = " << Amat_ref->at(i, j) << " != " << val
                  << " (err = " << SUNRabs(Amat_ref->at(i, j) - val) << ")\n";
        failure += check;
      }
    }
  }

  return failure > 0;
}

extern "C" int check_matrix_entry(SUNMatrix A, sunrealtype val, sunrealtype tol)
{
  if (using_csr_matrix_type) { return check_matrix_entry_csr(A, val, tol); }
  else if (using_dense_matrix_type)
  {
    return check_matrix_entry_dense(A, val, tol);
  }
  else { return 1; }
}

extern "C" int check_vector(N_Vector expected, N_Vector computed, sunrealtype tol)
{
  int failure{0};

  /* copy vectors to host */
  REF_OR_OMP_OR_HIP_OR_CUDA_OR_SYCL(, , N_VCopyFromDevice_Hip(computed),
                                    N_VCopyFromDevice_Cuda(computed),
                                    N_VCopyFromDevice_Sycl(computed));
  REF_OR_OMP_OR_HIP_OR_CUDA_OR_SYCL(, , N_VCopyFromDevice_Hip(expected),
                                    N_VCopyFromDevice_Cuda(expected),
                                    N_VCopyFromDevice_Sycl(expected));

  /* get vector data */
  auto xdata{N_VGetArrayPointer(computed)};
  auto ydata{N_VGetArrayPointer(expected)};

  /* check data lengths */
  auto xldata{N_VGetLength(computed)};
  auto yldata{N_VGetLength(expected)};

  if (xldata != yldata)
  {
    std::cerr << "ERROR check_vector: different vector lengths\n";
    return 1;
  }

  /* check vector data */
  for (sunindextype i = 0; i < xldata; i++)
  {
    failure += SUNRCompareTol(xdata[i], ydata[i], tol);
  }

  if (failure > ZERO)
  {
    std::cerr << "Check_vector failures:\n";
    for (sunindextype i = 0; i < xldata; i++)
    {
      if (SUNRCompareTol(xdata[i], ydata[i], tol) != 0)
      {
        std::cerr << "  computed[" << i << "] = " << xdata[i]
                  << " != " << ydata[i]
                  << " (err = " << SUNRabs(xdata[i] - ydata[i]) << ")\n";
      }
    }
  }

  return failure > 0;
}

extern "C" sunbooleantype has_data(SUNMatrix A)
{
  if (using_csr_matrix_type)
  {
    auto Amat{
      static_cast<sundials::ginkgo::Matrix<GkoCsrMat>*>(A->content)->GkoMtx()};
    return !(Amat->get_values() == nullptr || Amat->get_size()[0] == 0 ||
             Amat->get_size()[1] == 0);
  }
  else if (using_dense_matrix_type)
  {
    auto Amat{
      static_cast<sundials::ginkgo::Matrix<GkoDenseMat>*>(A->content)->GkoMtx()};
    return !(Amat->get_values() == nullptr || Amat->get_size()[0] == 0 ||
             Amat->get_size()[1] == 0);
  }
  else { return SUNFALSE; }
}

extern "C" sunbooleantype is_square(SUNMatrix A)
{
  if (using_csr_matrix_type)
  {
    auto Amat{
      static_cast<sundials::ginkgo::Matrix<GkoCsrMat>*>(A->content)->GkoMtx()};
    return Amat->get_size()[0] == Amat->get_size()[1];
  }
  else if (using_dense_matrix_type)
  {
    auto Amat{
      static_cast<sundials::ginkgo::Matrix<GkoDenseMat>*>(A->content)->GkoMtx()};
    return Amat->get_size()[0] == Amat->get_size()[1];
  }

  // unknown matrix type
  return SUNFALSE;
}

extern "C" void sync_device(SUNMatrix A)
{
  if (using_csr_matrix_type)
  {
    static_cast<sundials::ginkgo::Matrix<GkoCsrMat>*>(A->content)
      ->GkoExec()
      ->synchronize();
  }
  else if (using_dense_matrix_type)
  {
    static_cast<sundials::ginkgo::Matrix<GkoDenseMat>*>(A->content)
      ->GkoExec()
      ->synchronize();
  }
}
