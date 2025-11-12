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
 * Test the Ginkgo SUNMatrix batch implementation
 * ---------------------------------------------------------------------------*/

#include <iostream>
#include <memory>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_ginkgobatch.hpp>

#include "test_sunmatrix.h"

#if defined(USE_HIP)
#define REF_OR_OMP_OR_HIP_OR_CUDA(a, b, c, d) c
#elif defined(USE_CUDA)
#define REF_OR_OMP_OR_HIP_OR_CUDA(a, b, c, d) d
#elif defined(USE_OMP)
#define REF_OR_OMP_OR_HIP_OR_CUDA(a, b, c, d) b
#else
#define REF_OR_OMP_OR_HIP_OR_CUDA(a, b, c, d) a
#endif

#if defined(USE_CUDA)
#include <nvector/nvector_cuda.h>
#elif defined(USE_HIP)
#include <nvector/nvector_hip.h>
#elif defined(USE_OMP)
#include <nvector/nvector_openmp.h>
#else
#include <nvector/nvector_serial.h>
#endif

using namespace sundials;
using namespace sundials::ginkgo;

using GkoDenseMat = gko::matrix::Dense<sunrealtype>;
#ifdef SUNDIALS_INT32_T
using GkoCsrMat = gko::matrix::Csr<sunrealtype, sunindextype>;
#endif
using GkoVecType = GkoDenseMat;

bool using_csr_matrix_type   = false;
bool using_dense_matrix_type = false;

template<class GkoMatType>
int Test_CopyAndMove(std::unique_ptr<GkoMatType>&& gko_mat, SUNContext sunctx)
{
  // Copy constructor
  sundials::ginkgo::BatchMatrix<GkoMatType> mat{std::move(gko_mat), sunctx};
  sundials::ginkgo::BatchMatrix<GkoMatType> mat2{mat};
  SUNMatScaleAddI(sunrealtype{1.0}, mat2);
  SUNMatScaleAddI(sunrealtype{-1.0}, mat2);
  check_matrix_entry(mat2, sunrealtype{0.0}, SUN_UNIT_ROUNDOFF);

  // Move constructor
  sundials::ginkgo::BatchMatrix<GkoMatType> mat3{std::move(mat2)};
  SUNMatScaleAddI(sunrealtype{1.0}, mat3);
  SUNMatScaleAddI(sunrealtype{-1.0}, mat3);
  check_matrix_entry(mat3, sunrealtype{0.0}, SUN_UNIT_ROUNDOFF);

  // Copy assignment
  sundials::ginkgo::BatchMatrix<GkoMatType> mat4;
  mat4 = mat3;
  SUNMatScaleAddI(sunrealtype{1.0}, mat4);
  SUNMatScaleAddI(sunrealtype{-1.0}, mat4);
  check_matrix_entry(mat4, sunrealtype{0.0}, SUN_UNIT_ROUNDOFF);

  // Move assignment
  sundials::ginkgo::BatchMatrix<GkoMatType> mat5;
  mat5 = std::move(mat4);
  SUNMatScaleAddI(sunrealtype{1.0}, mat5);
  SUNMatScaleAddI(sunrealtype{-1.0}, mat5);
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
  Context sunctx;

  auto gko_exec{
    REF_OR_OMP_OR_HIP_OR_CUDA(gko::ReferenceExecutor::create(),
                              gko::OmpExecutor::create(),
                              gko::HipExecutor::create(0,
                                                       gko::OmpExecutor::create()),
                              gko::CudaExecutor::create(0,
                                                        gko::OmpExecutor::create()))};

  /* check input and set vector length */
  if (argc < 5)
  {
    std::cerr << "ERROR: FOUR (4) Input required: matrix rows, matrix cols, "
                 "number of batches (batch entries), format "
                 "(0 = csr, 1 = dense)\n";
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

  auto num_batches{static_cast<gko::size_type>(atol(argv[++argi]))};

  auto format{static_cast<int>(atoi(argv[++argi]))};
#ifdef SUNDIALS_INT32_T
  if (format != 0 && format != 1)
  {
    std::cerr << "ERROR: format must be 0 (csr) or 1 (dense) \n";
    return 1;
  }
#else
  if (format != 1)
  {
    std::cerr << "ERROR: format must be 1 (dense) with 64-bit index types\n";
    return 1;
  }
#endif

  if (format == 0) { using_csr_matrix_type = true; }
  else if (format == 1) { using_dense_matrix_type = true; }

#if defined(USE_OMP)
  int num_threads{1};
  auto omp_num_threads_var{std::getenv("OMP_NUM_THREADS")};
  if (omp_num_threads_var) { num_threads = std::atoi(omp_num_threads_var); }
#endif

  SetTiming(0);

  int square{matrows == matcols ? 1 : 0};
  std::cout << "\n SUNMATRIX_GINKGOBATCH test: " << num_batches
            << " batches of size " << matrows << " x " << matcols << ", format ";
  if (using_csr_matrix_type) { std::cout << "csr\n"; }
  else if (using_dense_matrix_type) { std::cout << "dense\n"; }

  /* Create vectors and matrices */
  std::default_random_engine engine;
  std::uniform_real_distribution<sunrealtype> distribution{0.0, sunrealtype{1.0}};

  sunindextype tot_cols = static_cast<sunindextype>(num_batches) * matcols;
  sunindextype tot_rows = static_cast<sunindextype>(num_batches) * matrows;

  N_Vector x{
    REF_OR_OMP_OR_HIP_OR_CUDA(N_VNew_Serial(tot_cols, sunctx),
                              N_VNew_OpenMP(tot_cols, num_threads, sunctx),
                              N_VNew_Hip(tot_cols, sunctx),
                              N_VNew_Cuda(tot_cols, sunctx))};
  N_Vector y{
    REF_OR_OMP_OR_HIP_OR_CUDA(N_VNew_Serial(tot_rows, sunctx),
                              N_VNew_OpenMP(tot_rows, num_threads, sunctx),
                              N_VNew_Hip(tot_rows, sunctx),
                              N_VNew_Cuda(tot_rows, sunctx))};

  auto xdata{N_VGetArrayPointer(x)};
  for (sunindextype i = 0; i < tot_cols; i++)
  {
    xdata[i] = distribution(engine);
  }
  REF_OR_OMP_OR_HIP_OR_CUDA(, , N_VCopyToDevice_Hip(x), N_VCopyToDevice_Cuda(x));

  auto matrix_dim{gko::dim<2>(matrows, matcols)};
  auto batch_mat_size{gko::batch_dim<2>(num_batches, matrix_dim)};
  auto batch_vec_size{gko::batch_dim<2>(num_batches, gko::dim<2>(matrows, 1))};
  auto gko_matdata{
    gko::matrix_data<sunrealtype, sunindextype>(matrix_dim, distribution, engine)};

  /* Wrap ginkgo matrices for SUNDIALS.
     sundials::ginkgo::Matrix is overloaded to a SUNMatrix. */
  std::unique_ptr<ConvertibleTo<SUNMatrix>> A;
  std::unique_ptr<ConvertibleTo<SUNMatrix>> I;

  /* Compute true solution */
  SUNMatrix Aref{SUNDenseMatrix(tot_rows, tot_cols, sunctx)};
  if (using_csr_matrix_type)
  {
#ifdef SUNDIALS_INT32_T
    auto gko_matrix{GkoCsrMat::create(gko_exec, matrix_dim)};
    gko_matrix->read(gko_matdata);
    auto num_nnz     = gko_matrix->get_num_stored_elements();
    auto common_size = gko_matrix->get_size();
    auto gko_batch_matrix{
      GkoBatchCsrMat::create(gko_exec, gko::batch_dim<2>(num_batches, common_size),
                             num_nnz)};
    for (gko::size_type b = 0; b < num_batches; ++b)
    {
      gko_batch_matrix->create_view_for_item(b)->read(gko_matdata);
    }

    auto gko_ident{GkoCsrMat::create(gko_exec, matrix_dim)};

    auto ident_data =
      gko::matrix_data<sunrealtype, sunindextype>::diag(batch_mat_size
                                                          .get_common_size(),
                                                        sunrealtype{1.0});
    gko_ident->read(ident_data);
    auto gko_batch_ident{
      GkoBatchCsrMat::create(gko_exec, gko::batch_dim<2>(num_batches, common_size),
                             common_size[0])};
    for (gko::size_type b = 0; b < num_batches; ++b)
    {
      gko_batch_ident->create_view_for_item(b)->read(ident_data);
    }

    auto Arowptrs{gko_batch_matrix->get_const_row_ptrs()};
    auto Acolidxs{gko_batch_matrix->get_const_col_idxs()};
    auto Avalues{gko_batch_matrix->get_const_values()};
    for (gko::size_type batchi = 0; batchi < num_batches; batchi++)
    {
      for (sunindextype irow = 0;
           irow < static_cast<sunindextype>(gko_matrix->get_size()[0]); irow++)
      {
        for (sunindextype inz = gko_exec->copy_val_to_host(Arowptrs + irow);
             inz < gko_exec->copy_val_to_host(Arowptrs + irow + 1); inz++)
        {
          SM_ELEMENT_D(Aref, batchi * matrows + irow,
                       batchi * matcols +
                         gko_exec->copy_val_to_host(Acolidxs + inz)) =
            gko_exec->copy_val_to_host(Avalues + inz);
        }
      }
    }

    fails += Test_CopyAndMove(gko_batch_matrix->clone(), sunctx);

    A = std::make_unique<BatchMatrix<GkoBatchCsrMat>>(std::move(gko_batch_matrix),
                                                      sunctx);
    I = std::make_unique<BatchMatrix<GkoBatchCsrMat>>(std::move(gko_batch_ident),
                                                      sunctx);
#else
    return 1;
#endif
  }
  else if (using_dense_matrix_type)
  {
    auto gko_matrix{GkoDenseMat::create(gko_exec, matrix_dim)};
    gko_matrix->read(gko_matdata);

    auto common_size = gko_matrix->get_size();
    auto gko_batch_matrix{
      GkoBatchDenseMat::create(gko_exec,
                               gko::batch_dim<2>(num_batches, common_size))};
    for (gko::size_type b = 0; b < num_batches; ++b)
    {
      gko_batch_matrix->create_view_for_item(b)->read(gko_matdata);
    }

    auto gko_ident{GkoDenseMat::create(gko_exec, matrix_dim)};
    auto ident_data =
      gko::matrix_data<sunrealtype, sunindextype>::diag(batch_mat_size
                                                          .get_common_size(),
                                                        sunrealtype{1.0});
    gko_ident->read(ident_data);
    auto gko_batch_ident{
      GkoBatchDenseMat::create(gko_exec,
                               gko::batch_dim<2>(num_batches, common_size))};
    for (gko::size_type b = 0; b < num_batches; ++b)
    {
      gko_batch_ident->create_view_for_item(b)->read(ident_data);
    }

    auto Avalues{gko_batch_matrix->get_const_values()};
    for (gko::size_type batchi = 0; batchi < num_batches; batchi++)
    {
      for (sunindextype i = 0; i < matrows; i++)
      {
        for (sunindextype j = 0; j < matcols; j++)
        {
          SM_ELEMENT_D(Aref, batchi * matrows + i, batchi * matcols + j) =
            gko_exec->copy_val_to_host(
              Avalues + gko_batch_matrix->get_cumulative_offset(batchi) +
              i * gko_batch_matrix->get_size().get_common_size()[1] + j);
        }
      }
    }

    fails += Test_CopyAndMove(gko_batch_matrix->clone(), sunctx);

    A = std::make_unique<BatchMatrix<GkoBatchDenseMat>>(std::move(gko_batch_matrix),
                                                        sunctx);
    I = std::make_unique<BatchMatrix<GkoBatchDenseMat>>(std::move(gko_batch_ident),
                                                        sunctx);
  }
  SUNMatMatvec_Dense(Aref, x, y);
  SUNMatDestroy(Aref);

  REF_OR_OMP_OR_HIP_OR_CUDA(, , N_VCopyToDevice_Hip(y), N_VCopyToDevice_Cuda(y));

  /* SUNMatrix Tests */
  fails += Test_SUNMatGetID(*A, SUNMATRIX_GINKGOBATCH, 0);
  fails += Test_SUNMatClone(*A, 0);
  fails += Test_SUNMatCopy(*A, 0);
  if (square)
  {
    if (using_dense_matrix_type) { fails += Test_SUNMatScaleAdd(*A, *I, 0); }
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
#ifdef SUNDIALS_INT32_T
  int failure{0};
  auto Amat{static_cast<BatchMatrix<GkoBatchCsrMat>*>(A->content)->GkoMtx()};
  auto Bmat{static_cast<BatchMatrix<GkoBatchCsrMat>*>(B->content)->GkoMtx()};
  auto Amat_ref = Amat->clone(Amat->get_executor()->get_master());
  auto Bmat_ref = Bmat->clone(Bmat->get_executor()->get_master());
  const auto Avalues{Amat_ref->get_const_values()};
  const auto Bvalues{Bmat_ref->get_const_values()};

  /* check lengths */
  if (Amat->get_size() != Bmat->get_size())
  {
    std::cerr << ">>> ERROR: check_matrix: Different data array lengths \n";
    return 1;
  }

  /* compare data */
  for (sunindextype ivalue = 0;
       ivalue < static_cast<sunindextype>(Amat_ref->get_num_stored_elements());
       ivalue++)
  {
    failure += SUNRCompareTol(Avalues[ivalue], Bvalues[ivalue], tol);
  }

  return failure > 0;
#else
  return 1;
#endif
}

static int check_matrix_dense(SUNMatrix A, SUNMatrix B, sunrealtype tol)
{
  int failure{0};
  auto Amat{static_cast<BatchMatrix<GkoBatchDenseMat>*>(A->content)->GkoMtx()};
  auto Bmat{static_cast<BatchMatrix<GkoBatchDenseMat>*>(B->content)->GkoMtx()};
  auto Amat_ref = Amat->clone(Amat->get_executor()->get_master());
  auto Bmat_ref = Bmat->clone(Bmat->get_executor()->get_master());
  auto batches{Amat->get_num_batch_items()};
  auto rows{Amat->get_size().get_common_size()[0]};
  auto cols{Amat->get_size().get_common_size()[1]};

  /* check lengths */
  if (Amat->get_size() != Bmat->get_size())
  {
    std::cerr << ">>> ERROR: check_matrix: Different data array lengths \n";
    return 1;
  }

  /* compare data */
  for (sunindextype batchi = 0; batchi < static_cast<sunindextype>(batches);
       batchi++)
  {
    for (sunindextype i = 0; i < static_cast<sunindextype>(rows); i++)
    {
      for (sunindextype j = 0; j < static_cast<sunindextype>(cols); j++)
      {
        failure += SUNRCompareTol(Amat_ref->at(batchi, i, j),
                                  Bmat_ref->at(batchi, i, j), tol);
      }
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
#ifdef SUNDIALS_INT32_T
  int failure{0};
  auto Amat{static_cast<BatchMatrix<GkoBatchCsrMat>*>(A->content)->GkoMtx()};
  auto Amat_ref = Amat->clone(Amat->get_executor()->get_master());
  auto Avalues{Amat_ref->get_const_values()};

  /* compare data */
  for (sunindextype ivalue = 0;
       ivalue < static_cast<sunindextype>(Amat_ref->get_num_stored_elements());
       ivalue++)
  {
    failure += SUNRCompareTol(Avalues[ivalue], val, tol);
  }

  return failure > 0;
#else
  return 1;
#endif
}

static int check_matrix_entry_dense(SUNMatrix A, sunrealtype val, sunrealtype tol)
{
  int failure{0};
  auto Amat{static_cast<BatchMatrix<GkoBatchDenseMat>*>(A->content)->GkoMtx()};
  auto batches{Amat->get_num_batch_items()};
  auto rows{Amat->get_size().get_common_size()[0]};
  auto cols{Amat->get_size().get_common_size()[1]};
  auto Amat_ref = Amat->clone(Amat->get_executor()->get_master());

  /* compare data */
  for (sunindextype batchi = 0; batchi < static_cast<sunindextype>(batches);
       batchi++)
  {
    for (sunindextype i = 0; i < static_cast<sunindextype>(rows); i++)
    {
      for (sunindextype j = 0; j < static_cast<sunindextype>(cols); j++)
      {
        failure += SUNRCompareTol(Amat_ref->at(batchi, i, j), val, tol);
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
  REF_OR_OMP_OR_HIP_OR_CUDA(, , N_VCopyFromDevice_Hip(computed),
                            N_VCopyFromDevice_Cuda(computed));
  REF_OR_OMP_OR_HIP_OR_CUDA(, , N_VCopyFromDevice_Hip(expected),
                            N_VCopyFromDevice_Cuda(expected));

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
    failure += SUNRCompareTol(xdata[i], ydata[i], tol);

  if (failure > ZERO)
  {
    std::cerr << "Check_vector failures:\n";
    for (sunindextype i = 0; i < xldata; i++)
      if (SUNRCompareTol(xdata[i], ydata[i], tol) != 0)
        std::cerr << "  computed[" << i << "] = " << xdata[i]
                  << " != " << ydata[i]
                  << " (err = " << SUNRabs(xdata[i] - ydata[i]) << ")\n";
  }

  return failure > 0;
}

extern "C" sunbooleantype has_data(SUNMatrix A)
{
  if (using_csr_matrix_type)
  {
#ifdef SUNDIALS_INT32_T
    auto Amat{static_cast<BatchMatrix<GkoBatchCsrMat>*>(A->content)->GkoMtx()};
    return !(Amat->get_values() == NULL || Amat->get_num_batch_items() == 0);
#else
    return SUNFALSE;
#endif
  }
  else if (using_dense_matrix_type)
  {
    auto Amat{static_cast<BatchMatrix<GkoBatchDenseMat>*>(A->content)->GkoMtx()};
    return !(Amat->get_values() == NULL || Amat->get_num_batch_items() == 0);
  }
  else { return SUNFALSE; }
}

extern "C" sunbooleantype is_square(SUNMatrix A)
{
  if (using_csr_matrix_type)
  {
#ifdef SUNDIALS_INT32_T
    auto Amat{static_cast<BatchMatrix<GkoBatchCsrMat>*>(A->content)->GkoMtx()};
    if (Amat->get_size().get_common_size()[0] !=
        Amat->get_size().get_common_size()[1])
    {
      return SUNFALSE;
    }
    return SUNTRUE;
#else
    return SUNFALSE;
#endif
  }
  else if (using_dense_matrix_type)
  {
    auto Amat{static_cast<BatchMatrix<GkoBatchDenseMat>*>(A->content)->GkoMtx()};
    if (Amat->get_size().get_common_size()[0] !=
        Amat->get_size().get_common_size()[1])
    {
      return SUNFALSE;
    }
    return SUNTRUE;
  }

  // unknown matrix type
  return SUNFALSE;
}

extern "C" void sync_device(SUNMatrix A)
{
  REF_OR_OMP_OR_HIP_OR_CUDA(, , hipDeviceSynchronize(), cudaDeviceSynchronize());
}
