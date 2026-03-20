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
 * SUNMatrix interface to Ginkgo matrices
 * ---------------------------------------------------------------------------*/

#ifndef _SUNMATRIX_GINKGO_HPP
#define _SUNMATRIX_GINKGO_HPP

#include <ginkgo/ginkgo.hpp>
#include <memory>
#include <sundials/sundials_base.hpp>
#include <sundials/sundials_matrix.hpp>
#include <utility>

#include "sundials/sundials_errors.h"

namespace sundials {
namespace ginkgo {

using GkoDenseMat = gko::matrix::Dense<sunrealtype>;
using GkoCsrMat   = gko::matrix::Csr<sunrealtype, sunindextype>;
using GkoVecType  = GkoDenseMat;

// Forward declaration of regular Matrix class
template<typename GkoMatType>
class Matrix;

// =============================================================================
// Everything in the implementation (impl) namespace is private and should not
// be referred to directly in user code.
// =============================================================================

namespace impl {

//
// Prototypes for non-class methods that operate on Matrix
//

inline std::unique_ptr<GkoVecType> WrapVector(
  std::shared_ptr<const gko::Executor> gko_exec, N_Vector x);

inline std::unique_ptr<const GkoVecType> WrapConstVector(
  std::shared_ptr<const gko::Executor> gko_exec, N_Vector x);

template<typename GkoMatType>
void Print(Matrix<GkoMatType>& A, std::ostream& ost = std::cout);

template<typename GkoMatType>
void Matvec(Matrix<GkoMatType>& A, GkoVecType* x, GkoVecType* y);

template<typename GkoMatType>
void Matvec(Matrix<GkoMatType>& A, N_Vector x, N_Vector y);

template<typename GkoMatType>
void ScaleAdd(const sunrealtype c, Matrix<GkoMatType>& A, Matrix<GkoMatType>& B);

template<typename GkoMatType>
void ScaleAddI(const sunrealtype c, Matrix<GkoMatType>& A);

template<typename GkoMatType>
void Zero(Matrix<GkoMatType>& A);

template<>
inline void Zero(Matrix<GkoDenseMat>& A);

template<typename GkoMatType>
void Copy(Matrix<GkoMatType>& A, Matrix<GkoMatType>& B);

//
// Methods that operate on SUNMatrix
//

template<typename GkoMatType>
SUNMatrix_ID SUNMatGetID_Ginkgo(SUNMatrix A)
{
  return SUNMATRIX_GINKGO;
}

template<typename GkoMatType>
SUNMatrix SUNMatClone_Ginkgo(SUNMatrix A)
{
  auto A_mat{static_cast<Matrix<GkoMatType>*>(A->content)};
  auto new_mat{new Matrix<GkoMatType>(*A_mat)}; // NOLINT
  return new_mat->get();
}

template<typename GkoMatType>
void SUNMatDestroy_Ginkgo(SUNMatrix A)
{
  auto A_mat{static_cast<Matrix<GkoMatType>*>(A->content)};
  delete A_mat; // NOLINT
  return;
}

template<typename GkoMatType>
SUNErrCode SUNMatZero_Ginkgo(SUNMatrix A)
{
  auto A_mat{static_cast<Matrix<GkoMatType>*>(A->content)};
  impl::Zero(*A_mat);
  return SUN_SUCCESS;
}

template<typename GkoMatType>
SUNErrCode SUNMatCopy_Ginkgo(SUNMatrix A, SUNMatrix B)
{
  auto A_mat{static_cast<Matrix<GkoMatType>*>(A->content)};
  auto B_mat{static_cast<Matrix<GkoMatType>*>(B->content)};
  impl::Copy(*A_mat, *B_mat);
  return SUN_SUCCESS;
}

template<typename GkoMatType>
SUNErrCode SUNMatScaleAdd_Ginkgo(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  auto A_mat{static_cast<Matrix<GkoMatType>*>(A->content)};
  auto B_mat{static_cast<Matrix<GkoMatType>*>(B->content)};
  impl::ScaleAdd(c, *A_mat, *B_mat);
  return SUN_SUCCESS;
}

template<typename GkoMatType>
SUNErrCode SUNMatScaleAddI_Ginkgo(sunrealtype c, SUNMatrix A)
{
  auto A_mat{static_cast<Matrix<GkoMatType>*>(A->content)};
  impl::ScaleAddI(c, *A_mat);
  return SUN_SUCCESS;
}

template<typename GkoMatType>
SUNErrCode SUNMatMatvec_Ginkgo(SUNMatrix A, N_Vector x, N_Vector y)
{
  auto A_mat{static_cast<Matrix<GkoMatType>*>(A->content)};
  impl::Matvec(*A_mat, x, y);
  return SUN_SUCCESS;
}

} // namespace impl

// =============================================================================
// Public namespace
// =============================================================================

/// Class that wraps a Ginkgo matrix and allows it to convert to a fully
/// functioning `SUNMatrix`.
template<typename GkoMatType>
class Matrix : public sundials::impl::BaseMatrix,
               public sundials::ConvertibleTo<SUNMatrix>
{
public:
  /// Default constructor - means the matrix must be copied or moved to
  Matrix() = default;

  /// Constructs a Matrix from an existing Ginkgo matrix object.
  /// \param gko_mat A Ginkgo matrix object
  /// \param sunctx The SUNDIALS simulation context object
  Matrix(std::shared_ptr<GkoMatType> gko_mat, SUNContext sunctx)
    : sundials::impl::BaseMatrix(sunctx), gkomtx_(gko_mat)
  {
    initSUNMatrix();
  }

  /// Move constructor
  Matrix(Matrix&& that_matrix) noexcept
    : sundials::impl::BaseMatrix(std::forward<Matrix>(that_matrix)),
      gkomtx_(std::move(that_matrix.gkomtx_))
  {}

  /// Copy constructor clones the ``gko::matrix`` and ``SUNMatrix``
  Matrix(const Matrix& that_matrix)
    : sundials::impl::BaseMatrix(that_matrix),
      gkomtx_(gko::clone(that_matrix.gkomtx_))
  {}

  /// Move assignment
  Matrix& operator=(Matrix&& rhs) noexcept
  {
    gkomtx_ = std::move(rhs.gkomtx_);

    sundials::impl::BaseMatrix::operator=(std::forward<Matrix>(rhs));

    return *this;
  }

  /// Copy assignment clones the gko::matrix and SUNMatrix
  Matrix& operator=(const Matrix& rhs)
  {
    gkomtx_ = gko::clone(rhs.gkomtx_);

    sundials::impl::BaseMatrix::operator=(rhs);

    return *this;
  }

  /// Default destructor
  // fine since all members are RAII
  virtual ~Matrix() = default;

  /// Get the underlying Ginkgo matrix object
  std::shared_ptr<GkoMatType> GkoMtx() const { return gkomtx_; }

  /// Get the ``gko::Executor`` associated with the Ginkgo matrix
  std::shared_ptr<const gko::Executor> GkoExec() const
  {
    return GkoMtx()->get_executor();
  }

  /// Get the size, i.e. ``gko::dim``, for the Ginkgo matrix
  const gko::dim<2>& GkoSize() const { return GkoMtx()->get_size(); }

  using sundials::impl::BaseMatrix::sunctx;

  // Override the ConvertibleTo methods

  /// Implicit conversion to a :c:type:`SUNMatrix`
  operator SUNMatrix() noexcept override { return object_.get(); }

  /// Implicit conversion to a :c:type:`SUNMatrix`
  operator SUNMatrix() const noexcept override { return object_.get(); }

  /// Explicit conversion to a :c:type:`SUNMatrix`
  SUNMatrix get() noexcept override { return object_.get(); }

  /// Explicit conversion to a :c:type:`SUNMatrix`
  SUNMatrix get() const noexcept override { return object_.get(); }

private:
  std::shared_ptr<GkoMatType> gkomtx_;

  void initSUNMatrix()
  {
    this->object_->content = this;

    this->object_->ops->getid     = impl::SUNMatGetID_Ginkgo<GkoMatType>;
    this->object_->ops->clone     = impl::SUNMatClone_Ginkgo<GkoMatType>;
    this->object_->ops->zero      = impl::SUNMatZero_Ginkgo<GkoMatType>;
    this->object_->ops->copy      = impl::SUNMatCopy_Ginkgo<GkoMatType>;
    this->object_->ops->scaleadd  = impl::SUNMatScaleAdd_Ginkgo<GkoMatType>;
    this->object_->ops->scaleaddi = impl::SUNMatScaleAddI_Ginkgo<GkoMatType>;
    this->object_->ops->matvec    = impl::SUNMatMatvec_Ginkgo<GkoMatType>;
    this->object_->ops->destroy   = impl::SUNMatDestroy_Ginkgo<GkoMatType>;
  }
};

// =============================================================================
// Everything in the implementation (impl) namespace is private and should not
// be referred to directly in user code.
// =============================================================================

namespace impl {

//
// Non-class methods
//

inline std::unique_ptr<GkoVecType> WrapVector(
  std::shared_ptr<const gko::Executor> gko_exec, N_Vector x)
{
  sunrealtype* x_arr{(x->ops->nvgetdevicearraypointer)
                       ? N_VGetDeviceArrayPointer(x)
                       : N_VGetArrayPointer(x)};
  const sunindextype x_len{N_VGetLength(x)};
  return GkoVecType::create(gko_exec, gko::dim<2>(x_len, 1),
                            gko::array<sunrealtype>::view(gko_exec, x_len, x_arr),
                            1);
}

inline std::unique_ptr<const GkoVecType> WrapConstVector(
  std::shared_ptr<const gko::Executor> gko_exec, N_Vector x)
{
  sunrealtype* x_arr{(x->ops->nvgetdevicearraypointer)
                       ? N_VGetDeviceArrayPointer(x)
                       : N_VGetArrayPointer(x)};
  const sunindextype x_len{N_VGetLength(x)};
  return GkoVecType::create_const(gko_exec, gko::dim<2>(x_len, 1),
                                  gko::array<sunrealtype>::const_view(gko_exec,
                                                                      x_len,
                                                                      x_arr),
                                  1);
}

template<typename GkoMatType>
void Print(Matrix<GkoMatType>& A, std::ostream& ost)
{
  gko::write(ost, A.GkoMtx().get());
}

template<typename GkoMatType>
void Matvec(Matrix<GkoMatType>& A, GkoVecType* x, GkoVecType* y)
{
  A.GkoMtx()->apply(x, y);
}

template<typename GkoMatType>
void Matvec(Matrix<GkoMatType>& A, N_Vector x, N_Vector y)
{
  if (x != y)
  {
    auto x_vec{WrapConstVector(A.GkoExec(), x)};
    auto y_vec{WrapVector(A.GkoExec(), y)};

    // y = Ax
    A.GkoMtx()->apply(x_vec.get(), y_vec.get());
  }
  else
  {
    auto x_vec{WrapVector(A.GkoExec(), x)};

    // x = Ax
    A.GkoMtx()->apply(x_vec.get(), x_vec.get());
  }
}

template<typename GkoMatType>
void ScaleAdd(const sunrealtype c, Matrix<GkoMatType>& A, Matrix<GkoMatType>& B)
{
  // Identity constructor always creates a square matrix
  const auto I{
    gko::matrix::Identity<sunrealtype>::create(A.GkoExec(), A.GkoSize()[1])};
  const auto one{gko::initialize<GkoDenseMat>({1.0}, A.GkoExec())};
  const auto cmat{gko::initialize<GkoDenseMat>({c}, A.GkoExec())};
  // A = B + cA
  B.GkoMtx()->apply(one.get(), I.get(), cmat.get(), A.GkoMtx().get());
}

template<typename GkoMatType>
void ScaleAddI(const sunrealtype c, Matrix<GkoMatType>& A)
{
  const auto one{gko::initialize<GkoDenseMat>({1.0}, A.GkoExec())};
  const auto cmat{gko::initialize<GkoDenseMat>({c}, A.GkoExec())};
  // A = 1*I + c*A = cA + I
  A.GkoMtx()->add_scaled_identity(one.get(), cmat.get());
}

template<typename GkoMatType>
void Zero(Matrix<GkoMatType>& A)
{
  A.GkoMtx()->scale(gko::initialize<GkoDenseMat>({0.0}, A.GkoExec()).get());
}

template<>
inline void Zero(Matrix<GkoDenseMat>& A)
{
  A.GkoMtx()->fill(0.0);
}

template<typename GkoMatType>
void Copy(Matrix<GkoMatType>& A, Matrix<GkoMatType>& B)
{
  B.GkoMtx()->copy_from(A.GkoMtx().get());
}

} // namespace impl

} // namespace ginkgo
} // namespace sundials

#endif
