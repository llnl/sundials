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
 * SUNLinearSolver interface to Ginkgo linear solvers
 * ---------------------------------------------------------------------------*/

#ifndef _SUNLINSOL_GINKGO_HPP
#define _SUNLINSOL_GINKGO_HPP

#include <cstring>
#include <ginkgo/ginkgo.hpp>
#include <memory>
#include <sundials/sundials_linearsolver.hpp>
#include <sundials/sundials_matrix.hpp>
#include <sunmatrix/sunmatrix_ginkgo.hpp>

namespace sundials {
namespace ginkgo {

template<class GkoSolverType, class GkoMatrixType>
class LinearSolver;

// =============================================================================
// Everything in the implementation (impl) namespace is private and should not
// be referred to directly in user code.
// =============================================================================

namespace impl {

inline SUNLinearSolver_Type SUNLinSolGetType_Ginkgo(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_MATRIX_ITERATIVE;
}

inline SUNLinearSolver_ID SUNLinSolGetID_Ginkgo(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_GINKGO;
}

template<class GkoSolverType, class GkoMatrixType>
int SUNLinSolSetup_Ginkgo(SUNLinearSolver S, SUNMatrix A)
{
  auto solver{
    static_cast<LinearSolver<GkoSolverType, GkoMatrixType>*>(S->content)};
  solver->Setup(static_cast<Matrix<GkoMatrixType>*>(A->content));
  return SUN_SUCCESS;
}

template<class GkoSolverType, class GkoMatrixType>
int SUNLinSolSolve_Ginkgo(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                          N_Vector b, sunrealtype tol)
{
  auto solver{
    static_cast<LinearSolver<GkoSolverType, GkoMatrixType>*>(S->content)};
  solver->Solve(b, x, tol);
  return SUN_SUCCESS;
}

template<class GkoSolverType, class GkoMatrixType>
SUNErrCode SUNLinSolFree_Ginkgo(SUNLinearSolver S)
{
  auto solver{
    static_cast<LinearSolver<GkoSolverType, GkoMatrixType>*>(S->content)};
  delete solver; // NOLINT
  return SUN_SUCCESS;
}

template<class GkoSolverType, class GkoMatrixType>
int SUNLinSolNumIters_Ginkgo(SUNLinearSolver S)
{
  auto solver{
    static_cast<LinearSolver<GkoSolverType, GkoMatrixType>*>(S->content)};
  return solver->NumIters();
}

template<class GkoSolverType, class GkoMatrixType>
sunrealtype SUNLinSolResNorm_Ginkgo(SUNLinearSolver S)
{
  auto solver{
    static_cast<LinearSolver<GkoSolverType, GkoMatrixType>*>(S->content)};
  return solver->ResNorm();
}

} // namespace impl

/// Custom gko::stop::Criterion that does the normal SUNDIALS stopping checks.
/// This checks if:
/// 1. Was the absolute residual tolerance met?
/// 2. Was the max iteration count reached?
class DefaultStop
  : public gko::EnablePolymorphicObject<DefaultStop, gko::stop::Criterion>
{
  friend class gko::EnablePolymorphicObject<DefaultStop, gko::stop::Criterion>;

public:
  GKO_CREATE_FACTORY_PARAMETERS(parameters, Factory)
  {
    sunrealtype GKO_FACTORY_PARAMETER_SCALAR(reduction_factor,
                                             SUN_UNIT_ROUNDOFF); // NOLINT(cppcoreguidelines-avoid-magic-numbers)
    gko::uint64 GKO_FACTORY_PARAMETER_SCALAR(max_iters, 5);
  };

  GKO_ENABLE_CRITERION_FACTORY(DefaultStop, parameters, Factory);
  GKO_ENABLE_BUILD_METHOD(Factory);

  gko::uint64 get_max_iters() const { return parameters_.max_iters; }

  sunrealtype get_reduction_factor() const
  {
    return parameters_.reduction_factor;
  }

protected:
  bool check_impl(gko::uint8 stoppingId, bool setFinalized,
                  gko::array<gko::stopping_status>* stop_status,
                  bool* one_changed, const Updater&) override;

  explicit DefaultStop(std::shared_ptr<const gko::Executor> exec)
    : EnablePolymorphicObject<DefaultStop, gko::stop::Criterion>(std::move(exec))
  {}

  explicit DefaultStop(const Factory* factory,
                       const gko::stop::CriterionArgs& args)
    : EnablePolymorphicObject<DefaultStop, gko::stop::Criterion>(
        factory->get_executor()),
      parameters_{factory->get_parameters()}
  {
    criteria_.push_back(gko::stop::ResidualNorm<sunrealtype>::build()
                          .with_reduction_factor(parameters_.reduction_factor)
                          .with_baseline(gko::stop::mode::absolute)
                          .on(factory->get_executor())
                          ->generate(args));
    criteria_.push_back(gko::stop::Iteration::build()
                          .with_max_iters(parameters_.max_iters)
                          .on(factory->get_executor())
                          ->generate(args));
  }

private:
  std::vector<std::unique_ptr<Criterion>> criteria_{};
};

inline bool DefaultStop::check_impl(gko::uint8 stoppingId, bool setFinalized,
                                    gko::array<gko::stopping_status>* stop_status,
                                    bool* one_changed, const Updater& updater)
{
  bool one_converged = false;
  gko::uint8 ids{1};
  *one_changed = false;
  for (auto& c : criteria_)
  {
    bool local_one_changed = false;
    one_converged |= c->check(ids, setFinalized, stop_status,
                              &local_one_changed, updater);
    *one_changed |= local_one_changed;
    if (one_converged) { break; }
    ids++;
  }
  return one_converged;
}

/// Class that wraps a Ginkgo solver (factory) and is convertible to a fully
/// functioning ``SUNLinearSolver``.
template<class GkoSolverType, class GkoMatrixType>
class LinearSolver : public ConvertibleTo<SUNLinearSolver>
{
public:
  /// Default constructor - means the solver must be moved to
  LinearSolver() = default;

  /// Constructs a new LinearSolver from a Ginkgo solver factory
  /// \param gko_solver_factory The Ginkgo solver factory (typically
  /// `gko::matrix::<type>::Factory`) \param sunctx The SUNDIALS simulation
  /// context (:c:type:`SUNContext`)
  LinearSolver(std::shared_ptr<typename GkoSolverType::Factory> gko_solver_factory,
               SUNContext sunctx)
    : gko_solver_factory_(gko_solver_factory),
      gko_solver_(nullptr),
      sunlinsol_(std::make_unique<_generic_SUNLinearSolver>()),
      sunlinsol_ops_(std::make_unique<_generic_SUNLinearSolver_Ops>()),
      res_norm_(sunrealtype{0.0})
  {
    sunlinsol_->content = this;
    sunlinsol_->ops     = sunlinsol_ops_.get();
    sunlinsol_->sunctx  = sunctx;

    sunlinsol_->ops->gettype = impl::SUNLinSolGetType_Ginkgo;
    sunlinsol_->ops->getid   = impl::SUNLinSolGetID_Ginkgo;
    sunlinsol_->ops->setup =
      impl::SUNLinSolSetup_Ginkgo<GkoSolverType, GkoMatrixType>;
    sunlinsol_->ops->solve =
      impl::SUNLinSolSolve_Ginkgo<GkoSolverType, GkoMatrixType>;
    sunlinsol_->ops->numiters =
      impl::SUNLinSolNumIters_Ginkgo<GkoSolverType, GkoMatrixType>;
    sunlinsol_->ops->resnorm =
      impl::SUNLinSolResNorm_Ginkgo<GkoSolverType, GkoMatrixType>;
    sunlinsol_->ops->free =
      impl::SUNLinSolFree_Ginkgo<GkoSolverType, GkoMatrixType>;
  }

  // Copy constructor is deleted
  LinearSolver(const LinearSolver& that_solver) = delete;

  /// Move constructor
  LinearSolver(LinearSolver&& that_solver) noexcept
    : gko_solver_factory_(std::move(that_solver.gko_solver_factory_)),
      gko_solver_(std::move(that_solver.gko_solver_)),
      sunlinsol_(std::move(that_solver.sunlinsol_)),
      sunlinsol_ops_(std::move(that_solver.sunlinsol_ops_)),
      iter_count_(that_solver.iter_count_),
      res_norm_(that_solver.res_norm_)
  {
    sunlinsol_->content = this;
    sunlinsol_->ops     = sunlinsol_ops_.get();
  }

  // Don't allow copy assignment
  LinearSolver& operator=(const LinearSolver& rhs) = delete;

  /// Move assignment
  LinearSolver& operator=(LinearSolver&& rhs)
  {
    gko_solver_factory_ = std::move(rhs.gko_solver_factory_);
    gko_solver_         = std::move(rhs.gko_solver_);
    sunlinsol_          = std::move(rhs.sunlinsol_);
    sunlinsol_ops_      = std::move(rhs.sunlinsol_ops_);
    iter_count_         = rhs.iter_count_;
    res_norm_           = rhs.res_norm_;
    sunlinsol_->content = this;
    sunlinsol_->ops     = sunlinsol_ops_.get();
    return *this;
  }

  /// Default destructor
  ~LinearSolver() override = default;

  /// Implicit conversion to a :c:type:`SUNLinearSolver`
  operator SUNLinearSolver() noexcept override { return sunlinsol_.get(); }

  /// Implicit conversion to a :c:type:`SUNLinearSolver`
  operator SUNLinearSolver() const noexcept override
  {
    return sunlinsol_.get();
  }

  /// Explicit conversion to a :c:type:`SUNLinearSolver`
  SUNLinearSolver get() noexcept override { return sunlinsol_.get(); }

  /// Explicit conversion to a :c:type:`SUNLinearSolver`
  SUNLinearSolver get() const noexcept override { return sunlinsol_.get(); }

  /// Get the ``gko::Executor`` associated with the Ginkgo solver
  std::shared_ptr<const gko::Executor> GkoExec() const
  {
    return gko_solver_factory_->get_executor();
  }

  /// Get the underlying Ginkgo solver factory
  std::shared_ptr<typename GkoSolverType::Factory> GkoFactory()
  {
    return gko_solver_factory_;
  }

  /// Get the underlying Ginkgo solver
  /// \note This will be `nullptr` until the linear solver setup phase.
  GkoSolverType* GkoSolver() { return gko_solver_.get(); }

  /// Get the number of linear solver iterations in the most recent solve.
  int NumIters() const { return iter_count_; }

  /// Get the residual norm of the solution at the end of the last solve.
  /// The type of residual norm depends on the Ginkgo stopping criteria
  /// used with the solver. With the \ref DefaultStop "DefaultStop" criteria
  /// this would be the absolute residual 2-norm.
  sunrealtype ResNorm() const { return res_norm_; }

  /// Setup the linear system
  /// \param A the linear system matrix
  GkoSolverType* Setup(Matrix<GkoMatrixType>* A)
  {
    gko_solver_ = gko_solver_factory_->generate(A->GkoMtx());
    return gko_solver_.get();
  }

  /// Solve the linear system Ax = b to the specified tolerance.
  /// \param b the right-hand side vector
  /// \param x the solution vector
  /// \param tol the tolerance to solve the system to
  gko::LinOp* Solve(N_Vector b, N_Vector x, sunrealtype tol)
  {
    auto logger{gko::share(gko::log::Convergence<sunrealtype>::create())};

    // Ginkgo provides a lot of options for stopping criterion,
    // so we make it possible to use it, but default to using
    // our normal iterative linear solver criterion.
    // If the criterion on the solver is of type DefaultStop,
    // then we will override the reduction_factor (tolerance).
    auto crit{dynamic_cast<const DefaultStop::Factory*>(
      GkoSolver()->get_stop_criterion_factory().get())};
    if (crit != nullptr)
    {
      auto new_crit = DefaultStop::build()
                        .with_reduction_factor(tol)
                        .with_max_iters(crit->get_parameters().max_iters)
                        .on(GkoExec());
      new_crit->add_logger(logger);
      GkoSolver()->set_stop_criterion_factory(std::move(new_crit));
    }

    gko::LinOp* result{nullptr};
    if (x != b)
    {
      auto x_vec = impl::WrapVector(GkoExec(), x);
      auto b_vec = impl::WrapVector(GkoExec(), b);

      // x = A^{-1} b
      result = GkoSolver()->apply(b_vec.get(), x_vec.get());
    }
    else
    {
      auto x_vec = impl::WrapVector(GkoExec(), x);

      // x = A^{-1} x
      result = GkoSolver()->apply(x_vec.get(), x_vec.get());
    }

    iter_count_ = static_cast<int>(logger->get_num_iterations());
    GkoExec()->get_master()->copy_from(GkoExec(), 1,
                                       gko::as<GkoDenseMat>(
                                         logger->get_residual_norm())
                                         ->get_const_values(),
                                       &res_norm_);

    return result;
  }

private:
  std::shared_ptr<typename GkoSolverType::Factory> gko_solver_factory_;
  std::unique_ptr<GkoSolverType> gko_solver_;
  std::unique_ptr<_generic_SUNLinearSolver> sunlinsol_;
  std::unique_ptr<_generic_SUNLinearSolver_Ops> sunlinsol_ops_;
  int iter_count_{};
  sunrealtype res_norm_{};
};

} // namespace ginkgo
} // namespace sundials

#endif // SUNLINSOL_GINKGO_HPP
