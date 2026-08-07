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
 * This is the header file for a SUNLinearSolver using Kokkoks Kernels
 * ---------------------------------------------------------------------------*/

#ifndef _SUNLINSOL_KOKKOSDENSE_HPP
#define _SUNLINSOL_KOKKOSDENSE_HPP

#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_Trsv_Decl.hpp>
#include <Kokkos_Core.hpp>
#include <nvector/nvector_kokkos.hpp>
#include <sundials/sundials_base.hpp>
#include <sundials/sundials_linearsolver.hpp>
#include <sunmatrix/sunmatrix_kokkosdense.hpp>

namespace sundials {
namespace kokkos {

// Forward declaration of DenseLinearSolver class
template<class ExecutionSpace = Kokkos::DefaultExecutionSpace,
         class MemorySpace    = typename ExecutionSpace::memory_space>
class DenseLinearSolver;

// =============================================================================
// Everything in the implementation (impl) namespace is private and should not
// be referred to directly in user code.
// =============================================================================

namespace impl {

static SUNLinearSolver_Type SUNLinSolGetType_KokkosDense(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_DIRECT;
}

static SUNLinearSolver_ID SUNLinSolGetID_KokkosDense(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_KOKKOSDENSE;
}

template<class MatrixType, class LinearSolverType>
int SUNLinSolSetup_KokkosDense(SUNLinearSolver S, SUNMatrix A)
{
  // Access matrix data
  auto A_mat{sundials::kokkos::GetDenseMat<MatrixType>(A)};

  auto A_exec = A_mat->ExecSpace();
  auto A_data = A_mat->View();

  const auto blocks = A_mat->Blocks();

  // Compute LU factorization of A (no pivoting)
  using team_policy = typename LinearSolverType::team_policy;
  using member_type = typename LinearSolverType::member_type;

  Kokkos::parallel_for(
    "sunlinsol_lu",
    team_policy(A_exec, static_cast<int>(blocks), Kokkos::AUTO, Kokkos::AUTO),
    KOKKOS_LAMBDA(const member_type& team_member) {
      const auto idx = team_member.league_rank();
      auto A_subdata = Kokkos::subview(A_data, idx, Kokkos::ALL(), Kokkos::ALL());
      KokkosBatched::TeamLU<
        member_type, KokkosBatched::Algo::LU::Unblocked>::invoke(team_member,
                                                                 A_subdata);
    });

  return SUN_SUCCESS;
}

template<class VectorType, class MatrixType, class LinearSolverType>
int SUNLinSolSolve_KokkosDense(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                               N_Vector b, sunrealtype tol)
{
  // Copy b into x
  N_VScale(SUN_RCONST(1.0), b, x);

  // Access matrix and vector data
  auto A_mat{sundials::kokkos::GetDenseMat<MatrixType>(A)};
  auto x_vec{sundials::kokkos::GetVec<VectorType>(x)};

  auto A_exec = A_mat->ExecSpace();
  auto A_data = A_mat->View();
  auto x_data = x_vec->View();

  const auto blocks = A_mat->Blocks();
  const auto rows   = A_mat->BlockRows();

  // Solve the linear system
  using team_policy = typename LinearSolverType::team_policy;
  using member_type = typename LinearSolverType::member_type;
  using size_type   = typename VectorType::size_type;

  Kokkos::parallel_for(
    "sunlinsol_trsv",
    team_policy(A_exec, static_cast<int>(blocks), Kokkos::AUTO, Kokkos::AUTO),
    KOKKOS_LAMBDA(const member_type& team_member) {
      const auto idx = team_member.league_rank();
      auto A_subdata = Kokkos::subview(A_data, idx, Kokkos::ALL(), Kokkos::ALL());
      auto x_subdata =
        Kokkos::subview(x_data,
                        Kokkos::pair<size_type, size_type>(idx * rows,
                                                           (idx + 1) * rows));
      // Lower triangular solve
      KokkosBatched::TeamVectorTrsv<
        member_type, KokkosBatched::Uplo::Lower,
        KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::Unit,
        KokkosBatched::Algo::Trsv::Unblocked>::invoke(team_member,
                                                      SUN_RCONST(1.0),
                                                      A_subdata, x_subdata);
      // Upper triangular solve
      KokkosBatched::TeamVectorTrsv<
        member_type, KokkosBatched::Uplo::Upper,
        KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::NonUnit,
        KokkosBatched::Algo::Trsv::Unblocked>::invoke(team_member,
                                                      SUN_RCONST(1.0),
                                                      A_subdata, x_subdata);
    });

  return SUN_SUCCESS;
}

template<class LinearSolverType>
SUNErrCode SUNLinSolFree_KokkosDense(SUNLinearSolver S)
{
  auto S_ls{static_cast<LinearSolverType*>(S->content)};
  delete S_ls; // NOLINT
  return SUN_SUCCESS;
}

} // namespace impl

// =============================================================================
// Public namespace
// =============================================================================

// -----------------------------------------------------------------------------
// Kokkos dense linear solver class, convertible to a SUNLinearSolver
// -----------------------------------------------------------------------------

template<class ExecutionSpace, class MemorySpace>
class DenseLinearSolver : public sundials::impl::BaseLinearSolver,
                          public sundials::ConvertibleTo<SUNLinearSolver>
{
public:
  using exec_space   = ExecutionSpace;
  using memory_space = MemorySpace;
  using team_policy  = typename Kokkos::TeamPolicy<exec_space>;
  using member_type  = typename Kokkos::TeamPolicy<exec_space>::member_type;

  // Default constructor - means the linear solver must be copied or moved to
  DenseLinearSolver() = default;

  DenseLinearSolver(SUNContext sunctx)
    : sundials::impl::BaseLinearSolver(sunctx)
  {
    initSUNLinearSolver();
  }

  // Move constructor
  DenseLinearSolver(DenseLinearSolver&& that_solver) noexcept
    : sundials::impl::BaseLinearSolver(
        std::forward<DenseLinearSolver>(that_solver))
  {}

  // Copy constructor
  DenseLinearSolver(const DenseLinearSolver& that_solver)
    : sundials::impl::BaseLinearSolver(that_solver)
  {}

  // Move assignment
  DenseLinearSolver& operator=(DenseLinearSolver&& rhs) noexcept
  {
    sundials::impl::BaseLinearSolver::operator=(
      std::forward<DenseLinearSolver>(rhs));
    return *this;
  }

  // Copy assignment
  DenseLinearSolver& operator=(const DenseLinearSolver& rhs)
  {
    sundials::impl::BaseLinearSolver::operator=(rhs);
    return *this;
  }

  // Default destructor since all members are RAII
  virtual ~DenseLinearSolver() = default;

  // Override the ConvertibleTo methods

  // Implicit conversion to a SUNLinearSolver
  operator SUNLinearSolver() noexcept override { return object_.get(); }

  // Implicit conversion to SUNLinearSolver
  operator SUNLinearSolver() const noexcept override { return object_.get(); }

  // Explicit conversion to a SUNLinearSolver
  SUNLinearSolver get() noexcept override { return object_.get(); }

  // Explicit conversion to a SUNLinearSolver
  SUNLinearSolver get() const noexcept override { return object_.get(); }

private:
  void initSUNLinearSolver()
  {
    using vec_type = Vector<exec_space, memory_space>;
    using mat_type = DenseMatrix<exec_space, memory_space>;
    using ls_type  = DenseLinearSolver<exec_space, memory_space>;

    this->object_->content = this;

    this->object_->ops->gettype = impl::SUNLinSolGetType_KokkosDense;
    this->object_->ops->getid   = impl::SUNLinSolGetID_KokkosDense;
    this->object_->ops->setup =
      impl::SUNLinSolSetup_KokkosDense<mat_type, ls_type>;
    this->object_->ops->solve =
      impl::SUNLinSolSolve_KokkosDense<vec_type, mat_type, ls_type>;
    this->object_->ops->free = impl::SUNLinSolFree_KokkosDense<ls_type>;
  }
};

} // namespace kokkos
} // namespace sundials

#endif
