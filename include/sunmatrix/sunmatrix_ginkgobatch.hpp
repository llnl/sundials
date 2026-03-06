/* -----------------------------------------------------------------
 * Programmer: Cody J. Balos @ LLNL
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
 * ----------------------------------------------------------------*/

#ifndef _SUNMATRIX_GINKGOBATCH_HPP
#define _SUNMATRIX_GINKGOBATCH_HPP

#include <memory>

#include <ginkgo/core/base/batch_multi_vector.hpp>
#include <ginkgo/ginkgo.hpp>

#include <sundials/sundials_config.h>
#include <sundials/sundials_matrix.hpp>

#if (GKO_VERSION_MAJOR < 1) || (GKO_VERSION_MAJOR == 1 && GKO_VERSION_MINOR < 9)
#error "Ginkgo 1.9.0 or later is required."
#endif

namespace sundials {
namespace ginkgo {

using GkoBatchDenseMat = gko::batch::matrix::Dense<sunrealtype>;
#ifdef SUNDIALS_INT32_T
// Currently (Ginkgo v1.10) only the dense matrix supports 64-bit index types
using GkoBatchCsrMat = gko::batch::matrix::Csr<sunrealtype, sunindextype>;
using GkoBatchEllMat = gko::batch::matrix::Ell<sunrealtype, sunindextype>;
#endif
using GkoBatchVecType = gko::batch::MultiVector<sunrealtype>;

// Forward declare BatchMatrix class
template<class GkoBatchMatType>
class BatchMatrix;

namespace impl {

//
// Prototypes for non-class methods that operate on Matrix
//

template<class GkoBatchMatType>
void Matvec(BatchMatrix<GkoBatchMatType>& A, GkoBatchVecType* x,
            GkoBatchVecType* y);

template<class GkoBatchMatType>
void Matvec(BatchMatrix<GkoBatchMatType>& A, N_Vector x, N_Vector y);

void ScaleAdd(const sunrealtype c, BatchMatrix<GkoBatchDenseMat>& A,
              BatchMatrix<GkoBatchDenseMat>& B);

#ifdef SUNDIALS_INT32_T
void ScaleAdd(const sunrealtype c, BatchMatrix<GkoBatchCsrMat>& A,
              BatchMatrix<GkoBatchCsrMat>& B);

void ScaleAdd(const sunrealtype c, BatchMatrix<GkoBatchEllMat>& A,
              BatchMatrix<GkoBatchEllMat>& B);
#endif

template<class GkoBatchMatType>
void ScaleAddI(const sunrealtype c, BatchMatrix<GkoBatchMatType>& A);

template<class GkoBatchMatType>
void Copy(BatchMatrix<GkoBatchMatType>& A, BatchMatrix<GkoBatchMatType>& B);

//
// Methods that operate on SUNMatrix
//
template<class GkoBatchMatType>
SUNMatrix_ID SUNMatGetID_GinkgoBatch(SUNMatrix A)
{
  return SUNMATRIX_GINKGOBATCH;
}

template<class GkoBatchMatType>
SUNMatrix SUNMatClone_GinkgoBatch(SUNMatrix A)
{
  auto Amat{static_cast<BatchMatrix<GkoBatchMatType>*>(A->content)};
  auto new_mat{new BatchMatrix<GkoBatchMatType>(*Amat)};
  return new_mat->get();
}

template<class GkoBatchMatType>
void SUNMatDestroy_GinkgoBatch(SUNMatrix A)
{
  auto Amat{static_cast<BatchMatrix<GkoBatchMatType>*>(A->content)};
  delete Amat; // NOLINT
  return;
}

template<class GkoBatchMatType>
SUNErrCode SUNMatCopy_GinkgoBatch(SUNMatrix A, SUNMatrix B)
{
  auto Amat{static_cast<BatchMatrix<GkoBatchMatType>*>(A->content)};
  Copy(*Amat, *static_cast<ginkgo::BatchMatrix<GkoBatchMatType>*>(B->content));
  return SUN_SUCCESS;
}

template<class GkoBatchMatType>
SUNErrCode SUNMatScaleAdd_GinkgoBatch(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  auto Amat{static_cast<BatchMatrix<GkoBatchMatType>*>(A->content)};
  ScaleAdd(c, *Amat,
           *static_cast<ginkgo::BatchMatrix<GkoBatchMatType>*>(B->content));
  return SUN_SUCCESS;
}

template<class GkoBatchMatType>
SUNErrCode SUNMatScaleAddI_GinkgoBatch(sunrealtype c, SUNMatrix A)
{
  auto Amat{static_cast<BatchMatrix<GkoBatchMatType>*>(A->content)};
  ScaleAddI(c, *Amat);
  return SUN_SUCCESS;
}

template<class GkoBatchMatType>
SUNErrCode SUNMatMatvec_GinkgoBatch(SUNMatrix A, N_Vector x, N_Vector y)
{
  auto Amat{static_cast<BatchMatrix<GkoBatchMatType>*>(A->content)};
  Matvec(*Amat, x, y);
  return SUN_SUCCESS;
}

} // namespace impl

template<class GkoBatchMatType>
class BatchMatrix : public sundials::impl::BaseMatrix,
                    public sundials::ConvertibleTo<SUNMatrix>
{
public:
  // Default constructor means the matrix must be copied or moved to
  BatchMatrix() = default;

  // We do not have implementations of these two constructors for general
  // GkoBatchMatType
  BatchMatrix(gko::size_type num_batches, sunindextype M, sunindextype N,
              std::shared_ptr<const gko::Executor> gko_exec, SUNContext sunctx);
  BatchMatrix(gko::size_type num_batches, sunindextype M, sunindextype N,
              sunindextype num_nonzeros,
              std::shared_ptr<const gko::Executor> gko_exec, SUNContext sunctx);

  BatchMatrix(std::shared_ptr<GkoBatchMatType> gko_mat, SUNContext sunctx)
    : sundials::impl::BaseMatrix(sunctx), gkomtx_(gko_mat)
  {
    initSUNMatrix();
  }

  // Move constructor
  BatchMatrix(BatchMatrix&& that_matrix) noexcept
    : sundials::impl::BaseMatrix(std::forward<BatchMatrix>(that_matrix)),
      gkomtx_(std::move(that_matrix.gkomtx_))
  {}

  // Copy constructor clones the gko::matrix and SUNMatrix
  BatchMatrix(const BatchMatrix& that_matrix)
    : sundials::impl::BaseMatrix(that_matrix),
      gkomtx_(gko::clone(that_matrix.gkomtx_))
  {}

  // Move assignment
  BatchMatrix& operator=(BatchMatrix&& rhs) noexcept
  {
    gkomtx_ = std::move(rhs.gkomtx_);
    sundials::impl::BaseMatrix::operator=(std::forward<BatchMatrix>(rhs));
    return *this;
  }

  // Copy assignment clones the gko::matrix and SUNMatrix
  BatchMatrix& operator=(const BatchMatrix& rhs)
  {
    gkomtx_ = gko::clone(rhs.gkomtx_);
    sundials::impl::BaseMatrix::operator=(rhs);
    return *this;
  }

  // Default destructor is fine since all members are RAII
  virtual ~BatchMatrix() = default;

  // Getters
  std::shared_ptr<GkoBatchMatType> GkoMtx() const { return gkomtx_; }

  std::shared_ptr<const gko::Executor> GkoExec() const
  {
    return GkoMtx()->get_executor();
  }

  const gko::batch_dim<2>& GkoSize() const { return GkoMtx()->get_size(); }

  gko::size_type NumBatches() const { return GkoSize().get_num_batch_items(); }

  using sundials::impl::BaseMatrix::sunctx;

  // Override the ConvertibleTo methods
  operator SUNMatrix() noexcept override { return object_.get(); }

  operator SUNMatrix() const noexcept override { return object_.get(); }

  SUNMatrix get() noexcept override { return object_.get(); }

  SUNMatrix get() const noexcept override { return object_.get(); }

private:
  std::shared_ptr<GkoBatchMatType> gkomtx_;

  void initSUNMatrix()
  {
    this->object_->content = this;

    this->object_->ops->getid = impl::SUNMatGetID_GinkgoBatch<GkoBatchMatType>;
    this->object_->ops->clone = impl::SUNMatClone_GinkgoBatch<GkoBatchMatType>;
    this->object_->ops->copy  = impl::SUNMatCopy_GinkgoBatch<GkoBatchMatType>;
    // Ginkgo does not provide what we need for ScaleAdd except with BatchDense
    this->object_->ops->scaleadd =
      impl::SUNMatScaleAdd_GinkgoBatch<GkoBatchMatType>;
    this->object_->ops->scaleaddi =
      impl::SUNMatScaleAddI_GinkgoBatch<GkoBatchMatType>;
    this->object_->ops->matvec = impl::SUNMatMatvec_GinkgoBatch<GkoBatchMatType>;
    this->object_->ops->destroy = impl::SUNMatDestroy_GinkgoBatch<GkoBatchMatType>;
  }
};

//
// Class method specializations for specific types of Ginkgo matrices.
//

template<>
inline BatchMatrix<GkoBatchDenseMat>::BatchMatrix(
  gko::size_type num_batches, sunindextype M, sunindextype N,
  std::shared_ptr<const gko::Executor> gko_exec, SUNContext sunctx)
  : sundials::impl::BaseMatrix(sunctx),
    gkomtx_(
      GkoBatchDenseMat::create(gko_exec, gko::batch_dim<2>(num_batches,
                                                           gko::dim<2>(M, N))))
{
  initSUNMatrix();
}

#ifdef SUNDIALS_INT32_T
template<>
inline BatchMatrix<GkoBatchCsrMat>::BatchMatrix(
  gko::size_type num_batches, sunindextype M, sunindextype N,
  sunindextype num_nonzeros, std::shared_ptr<const gko::Executor> gko_exec,
  SUNContext sunctx)
  : sundials::impl::BaseMatrix(sunctx),
    gkomtx_(
      GkoBatchCsrMat::create(gko_exec,
                             gko::batch_dim<2>(num_batches, gko::dim<2>(M, N)),
                             num_nonzeros))
{
  initSUNMatrix();
}

template<>
inline BatchMatrix<GkoBatchEllMat>::BatchMatrix(
  gko::size_type num_batches, sunindextype M, sunindextype N,
  sunindextype num_nonzeros, std::shared_ptr<const gko::Executor> gko_exec,
  SUNContext sunctx)
  : sundials::impl::BaseMatrix(sunctx),
    gkomtx_(
      GkoBatchEllMat::create(gko_exec,
                             gko::batch_dim<2>(num_batches, gko::dim<2>(M, N)),
                             num_nonzeros))
{
  initSUNMatrix();
}
#endif

// =============================================================================
// Everything in the implementation (impl) namespace is private and should not
// be referred to directly in user code.
// =============================================================================

namespace impl {

//
// Non-class methods
//

inline std::unique_ptr<GkoBatchVecType> WrapBatchVector(
  std::shared_ptr<const gko::Executor> gko_exec, gko::size_type num_batches,
  N_Vector x)
{
  auto x_arr{(x->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(x)
                                               : N_VGetArrayPointer(x)};
  const auto xvec_len{N_VGetLength(x)};
  auto batch_xvec_size{
    gko::batch_dim<2>(num_batches, gko::dim<2>(xvec_len / num_batches, 1))};
  auto xvec_view{gko::array<sunrealtype>::view(gko_exec, xvec_len, x_arr)};
  return GkoBatchVecType::create(gko_exec, batch_xvec_size, std::move(xvec_view));
}

inline std::unique_ptr<const GkoBatchVecType> WrapConstBatchVector(
  std::shared_ptr<const gko::Executor> gko_exec, gko::size_type num_batches,
  N_Vector x)
{
  auto x_arr{(x->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(x)
                                               : N_VGetArrayPointer(x)};
  const auto xvec_len{N_VGetLength(x)};
  auto batch_xvec_size{
    gko::batch_dim<2>(num_batches, gko::dim<2>(xvec_len / num_batches, 1))};
  auto xvec_view{gko::array<sunrealtype>::const_view(gko_exec, xvec_len, x_arr)};
  return GkoBatchVecType::create_const(gko_exec, batch_xvec_size,
                                       std::move(xvec_view));
}

template<class GkoBatchMatType>
void Matvec(BatchMatrix<GkoBatchMatType>& A, GkoBatchVecType* x,
            GkoBatchVecType* y)
{
  A.GkoMtx()->apply(x, y);
}

template<class GkoBatchMatType>
void Matvec(BatchMatrix<GkoBatchMatType>& A, N_Vector x, N_Vector y)
{
  if (x != y)
  {
    auto x_vec = WrapConstBatchVector(A.GkoExec(), A.NumBatches(), x);
    auto y_vec = WrapBatchVector(A.GkoExec(), A.NumBatches(), y);

    // y = Ax
    A.GkoMtx()->apply(x_vec.get(), y_vec.get());
  }
  else
  {
    auto x_vec = WrapBatchVector(A.GkoExec(), A.NumBatches(), x);

    // x = Ax
    A.GkoMtx()->apply(x_vec.get(), x_vec.get());
  }
}

void ScaleAdd(const sunrealtype c, BatchMatrix<GkoBatchDenseMat>& A,
              BatchMatrix<GkoBatchDenseMat>& B)
{
  auto cmat = GkoBatchVecType::create(A.GkoExec(),
                                      gko::batch_dim<2>(A.NumBatches(),
                                                        gko::dim<2>(1, 1)));
  cmat->fill(c);
  // A = B + cA
  A.GkoMtx()->scale_add(cmat.get(), B.GkoMtx().get());
}

#ifdef SUNDIALS_INT32_T
void ScaleAdd(const sunrealtype c, BatchMatrix<GkoBatchCsrMat>& A,
              BatchMatrix<GkoBatchCsrMat>& B)
{
  // NOTE: This is not implemented by Ginkgo for BatchCsr yet
  throw("scale_add not implemented for gko::batch::matrix::Csr");
}

void ScaleAdd(const sunrealtype c, BatchMatrix<GkoBatchEllMat>& A,
              BatchMatrix<GkoBatchEllMat>& B)
{
  // NOTE: This is not implemented by Ginkgo for BatchEll yet
  throw("scale_add not implemented for gko::batch::matrix::Ell");
}
#endif

template<class GkoBatchMatType>
void ScaleAddI(const sunrealtype c, BatchMatrix<GkoBatchMatType>& A)
{
  auto one = GkoBatchVecType::create(A.GkoExec(),
                                     gko::batch_dim<2>(A.NumBatches(),
                                                       gko::dim<2>(1, 1)));
  one->fill(gko::one<sunrealtype>());
  auto cmat = GkoBatchVecType::create(A.GkoExec(),
                                      gko::batch_dim<2>(A.NumBatches(),
                                                        gko::dim<2>(1, 1)));
  cmat->fill(c);
  // A = 1*I + c*A = cA + I
  A.GkoMtx()->add_scaled_identity(one.get(), cmat.get());
}

template<class GkoBatchMatType>
void Copy(BatchMatrix<GkoBatchMatType>& A, BatchMatrix<GkoBatchMatType>& B)
{
  B.GkoMtx()->copy_from(A.GkoMtx().get());
}

template<class GkoBatchMatType>
void Print(BatchMatrix<GkoBatchMatType>& A, std::ostream& ost = std::cout)
{
  gko::write(ost, A.GkoMtx().get());
}

} // namespace impl

} // namespace ginkgo
} // namespace sundials

#endif
