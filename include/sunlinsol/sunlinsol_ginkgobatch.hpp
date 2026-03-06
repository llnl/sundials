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
 * ---------------------------------------------------------------------------*/

#include <cmath>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/priv/sundials_logger_macros.h>
#include <sundials/sundials_core.hpp>

#include <sunmatrix/sunmatrix_ginkgobatch.hpp>

#if (GKO_VERSION_MAJOR < 1) || (GKO_VERSION_MAJOR == 1 && GKO_VERSION_MINOR < 9)
#error "Ginkgo 1.9.0 or later is required."
#endif

#ifndef _SUNLINSOL_GINKGOBATCH_HPP
#define _SUNLINSOL_GINKGOBATCH_HPP

namespace sundials {
namespace ginkgo {

template<class GkoBatchSolverType, class GkoBatchMatType>
class BatchLinearSolver;

namespace impl {

inline SUNLinearSolver_Type SUNLinSolGetType_GinkgoBatch(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_MATRIX_ITERATIVE;
}

inline SUNLinearSolver_ID SUNLinSolGetID_GinkgoBatch(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_GINKGOBATCH;
}

template<class GkoBatchLinearSolverType>
SUNErrCode SUNLinSolInitialize_GinkgoBatch(SUNLinearSolver S)
{
  return SUN_SUCCESS;
}

template<class GkoBatchLinearSolverType>
SUNErrCode SUNLinSolSetScalingVectors_GinkgoBatch(SUNLinearSolver S,
                                                  N_Vector s1, N_Vector s2)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  SUNFunctionBegin(S->sunctx);
  SUNAssert(s1 || s2, SUN_ERR_ARG_INCOMPATIBLE);
  return solver->SetScalingVectors(s1, s2);
}

template<class GkoBatchLinearSolverType, class GkoBatchMatType>
int SUNLinSolSetup_GinkgoBatch(SUNLinearSolver S, SUNMatrix A)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  return solver->Setup(static_cast<BatchMatrix<GkoBatchMatType>*>(A->content));
}

template<class GkoBatchLinearSolverType, class GkoBatchMatType>
int SUNLinSolSolve_GinkgoBatch(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                               N_Vector b, sunrealtype tol)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  return solver->Solve(static_cast<BatchMatrix<GkoBatchMatType>*>(A->content),
                       b, x, tol);
}

template<class GkoBatchLinearSolverType>
SUNErrCode SUNLinSolFree_GinkgoBatch(SUNLinearSolver S)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  delete solver; // NOLINT(cppcoreguidelines-owning-memory)
  return SUN_SUCCESS;
}

template<class GkoBatchLinearSolverType>
int SUNLinSolNumIters_GinkgoBatch(SUNLinearSolver S)
{
  auto solver = static_cast<GkoBatchLinearSolverType*>(S->content);
  return static_cast<int>(std::round(solver->AvgNumIters()));
}

inline gko::array<sunrealtype> WrapBatchScalingArray(
  std::shared_ptr<const gko::Executor> gko_exec, gko::size_type num_batches,
  N_Vector x)
{
  auto x_arr{(x->ops->nvgetdevicearraypointer) ? N_VGetDeviceArrayPointer(x)
                                               : N_VGetArrayPointer(x)};
  auto xvec_len{N_VGetLength(x)};
  return gko::array<sunrealtype>(gko_exec,
                                 gko::array<sunrealtype>::view(gko_exec,
                                                               xvec_len, x_arr));
}

inline std::unique_ptr<GkoBatchVecType> WrapAsMultiVector(
  gko::array<sunrealtype>& arr, gko::size_type num_batches)
{
  return GkoBatchVecType::create(arr.get_executor(),
                                 gko::batch_dim<2>(num_batches,
                                                   gko::dim<2>(arr.get_size() /
                                                                 num_batches,
                                                               1)),
                                 arr);
}

} // namespace impl

template<class GkoBatchSolverType, class GkoBatchMatType>
class BatchLinearSolver : public sundials::impl::BaseLinearSolver,
                          public ConvertibleTo<SUNLinearSolver>
{
public:
  constexpr static int NO_SCALING     = 0;
  constexpr static int LAGGED_SCALING = 1;
  constexpr static int SOLVE_SCALING  = 2;

  BatchLinearSolver(std::shared_ptr<const gko::Executor> gko_exec,
                    gko::batch::stop::tolerance_type tolerance_type,
                    std::shared_ptr<gko::batch::BatchLinOpFactory> precon_factory,
                    int max_iters, gko::size_type num_batches, SUNContext sunctx)
    : sundials::impl::BaseLinearSolver(sunctx),
      gko_exec_(std::move(gko_exec)),
      tolerance_type_(tolerance_type),
      precon_factory_(std::move(precon_factory)),
      solver_factory_(nullptr),
      solver_(nullptr),
      col_scale_vec_(gko_exec_),
      row_scale_vec_(gko_exec_),
      logger_(nullptr),
      num_batches_(num_batches),
      max_iters_(max_iters),
      res_norm_array_(gko_exec_->get_master(), num_batches_),
      iter_count_array_(gko_exec_->get_master(), num_batches_),
      avg_iter_count_(sunrealtype{0.0}),
      sum_of_avg_iters_(sunrealtype{0.0}),
      stddev_iter_count_(sunrealtype{0.0}),
      s1_(nullptr),
      s2inv_(nullptr),
      scaling_mode_(LAGGED_SCALING),
      scaling_initialized_(false),
      do_setup_(true)
  {
    initSUNLinSol(sunctx);
  }

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BatchLinearSolver(std::shared_ptr<const gko::Executor> gko_exec,
                    gko::size_type num_batches, SUNContext sunctx)
    : BatchLinearSolver(std::move(gko_exec),
                        gko::batch::stop::tolerance_type::absolute, nullptr,
                        default_max_iters_, num_batches, sunctx)
  {}

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BatchLinearSolver(std::shared_ptr<const gko::Executor> gko_exec,
                    gko::batch::stop::tolerance_type tolerance_type,
                    gko::size_type num_batches, SUNContext sunctx)
    : BatchLinearSolver(std::move(gko_exec), tolerance_type, nullptr,
                        default_max_iters_, num_batches, sunctx)
  {}

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BatchLinearSolver(std::shared_ptr<const gko::Executor> gko_exec,
                    std::shared_ptr<gko::batch::BatchLinOpFactory> precon_factory,
                    gko::size_type num_batches, SUNContext sunctx)
    : BatchLinearSolver(std::move(gko_exec),
                        gko::batch::stop::tolerance_type::absolute,
                        precon_factory, default_max_iters_, num_batches, sunctx)
  {}

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BatchLinearSolver(std::shared_ptr<const gko::Executor> gko_exec,
                    int max_iters, gko::size_type num_batches, SUNContext sunctx)
    : BatchLinearSolver(std::move(gko_exec),
                        gko::batch::stop::tolerance_type::absolute, nullptr,
                        max_iters, num_batches, sunctx)
  {}

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BatchLinearSolver(std::shared_ptr<const gko::Executor> gko_exec,
                    gko::batch::stop::tolerance_type tolerance_type,
                    int max_iters, gko::size_type num_batches, SUNContext sunctx)
    : BatchLinearSolver(std::move(gko_exec), tolerance_type, nullptr, max_iters,
                        num_batches, sunctx)
  {}

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BatchLinearSolver(std::shared_ptr<const gko::Executor> gko_exec,
                    std::shared_ptr<gko::batch::BatchLinOpFactory> precon_factory,
                    int max_iters, gko::size_type num_batches, SUNContext sunctx)
    : BatchLinearSolver(std::move(gko_exec),
                        gko::batch::stop::tolerance_type::absolute,
                        precon_factory, max_iters, num_batches, sunctx)
  {}

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  BatchLinearSolver(std::shared_ptr<const gko::Executor> gko_exec,
                    gko::batch::stop::tolerance_type tolerance_type,
                    std::shared_ptr<gko::batch::BatchLinOpFactory> precon_factory,
                    gko::size_type num_batches, SUNContext sunctx)
    : BatchLinearSolver(std::move(gko_exec), tolerance_type, precon_factory,
                        default_max_iters_, num_batches, sunctx)
  {}

  /// Get the ``gko::Executor`` associated with the Ginkgo matrix
  std::shared_ptr<const gko::Executor> GkoExec() const { return gko_exec_; }

  // Override the ConvertibleTo methods

  /// Implicit conversion to a :c:type:`SUNLinearSolver`
  operator SUNLinearSolver() noexcept override { return object_.get(); }

  /// Implicit conversion to a :c:type:`SUNLinearSolver`
  operator SUNLinearSolver() const noexcept override { return object_.get(); }

  /// Explicit conversion to a :c:type:`SUNLinearSolver`
  SUNLinearSolver get() noexcept override { return object_.get(); }

  /// Explicit conversion to a :c:type:`SUNLinearSolver`
  SUNLinearSolver get() const noexcept override { return object_.get(); }

  /// Get the underlying Ginkgo solver
  /// \note This will be `nullptr` until the linear solver setup phase.
  GkoBatchSolverType* GkoSolver() { return solver_.get(); }

  /// Average number of iterations across the batches during the last solve.
  sunrealtype AvgNumIters() const { return avg_iter_count_; }

  /// Standard deviation of the number of iterations across the batches during the last solve.
  sunrealtype StddevNumIters() const { return stddev_iter_count_; }

  /// Running sum of the average number of iterations in this solvers lifetime.
  sunrealtype SumAvgNumIters() const { return sum_of_avg_iters_; }

  SUNErrCode SetScalingMode(int scaling_mode = LAGGED_SCALING)
  {
    SUNFunctionBegin(sunCtx());
    SUNAssert(scaling_mode == NO_SCALING || scaling_mode == LAGGED_SCALING ||
                scaling_mode == SOLVE_SCALING,
              SUN_ERR_ARG_INCOMPATIBLE);
    scaling_mode_ = scaling_mode;
    return SUN_SUCCESS;
  }

  /// Sets the left and right scaling vectors to be used.
  SUNErrCode SetScalingVectors(N_Vector s1, N_Vector s2)
  {
    SUNFunctionBegin(sunCtx());

    if (s1 && s2)
    {
      s1_ = s1;
      col_scale_vec_ =
        std::move(impl::WrapBatchScalingArray(GkoExec(), num_batches_, s1_));

      if (!s2inv_.get())
      {
        s2inv_ = std::move(sundials::experimental::N_VectorView(N_VClone(s2)));
      }

      // SUNLinearSolver API wants s2inv_
      N_VInv(s2, s2inv_);
      SUNCheckLastErr();

      row_scale_vec_ =
        std::move(impl::WrapBatchScalingArray(GkoExec(), num_batches_, s2inv_));

      scaling_initialized_ = true;
    }
    else
    {
      scaling_mode_        = NO_SCALING;
      scaling_initialized_ = false;
    }

    return SUN_SUCCESS;
  }

  int Setup(BatchMatrix<GkoBatchMatType>* A)
  {
    if (num_batches_ != A->NumBatches()) { return SUN_ERR_ARG_OUTOFRANGE; }

    do_setup_ = true;

    return SUN_SUCCESS;
  }

  int Solve(BatchMatrix<GkoBatchMatType>* A, N_Vector b, N_Vector x,
            sunrealtype tol)
  {
    SUNLogInfo(sunLogger(), "linear-solver", "solver = ginkgo (batched)");

    // Applying scaling every solve requires us to redo all the Ginkgo solver setup
    // This is because the scaling happens in memory, and the Ginkgo solver needs
    // to know about the matrix values when it is constructed.
    if (scaling_mode_ == SOLVE_SCALING) { do_setup_ = true; }

    // SetScaling was called and scaling is used
    bool apply_scaling = scaling_initialized_ && scaling_mode_ != NO_SCALING;

    if (do_setup_)
    {
      // Create the solver factory
      solver_factory_ = GkoBatchSolverType::build()             //
                          .with_max_iterations(max_iters_)      //
                          .with_tolerance(SUN_UNIT_ROUNDOFF)    //
                          .with_tolerance_type(tolerance_type_) //
                          .with_preconditioner(precon_factory_) //
                          .on(GkoExec());

      if (apply_scaling)
      {
        // This will scale the matrix in place.
        // So if a Ginkgo preconditioner is used, it is applied to the scaled matrix:
        //    \tilde{A} = P_1^{-1} S_1 A S_2^{-1} P_2^{-1},
        // but in memory, we just have
        //    \tilde{A} = S_1 A S_2^{-1}.
        A->GkoMtx()->scale(row_scale_vec_, col_scale_vec_);
      }

      solver_ = solver_factory_->generate(A->GkoMtx());

      do_setup_ = false;
    }

    solver_->reset_tolerance(tol);

    if (!logger_)
    {
      logger_ =
        gko::share(gko::batch::log::BatchConvergence<sunrealtype>::create());
    }
    solver_->add_logger(logger_);

    std::unique_ptr<GkoBatchVecType> x_vec{
      impl::WrapBatchVector(GkoExec(), num_batches_, x)};
    std::unique_ptr<GkoBatchVecType> b_vec{
      impl::WrapBatchVector(GkoExec(), num_batches_, b)};

    if (apply_scaling)
    {
      // \tilde{b} = S_1 b
      b_vec->scale(impl::WrapAsMultiVector(col_scale_vec_, num_batches_));
    }

    // \tilde{x} = \tilde{A}^{-1} \tilde{b}
    [[maybe_unused]] gko::batch::BatchLinOp* result =
      solver_->apply(b_vec.get(), x_vec.get());

    if (apply_scaling)
    {
      // x = S_2^{-1} \tilde{x}
      x_vec->scale(impl::WrapAsMultiVector(row_scale_vec_, num_batches_));
    }

    if (apply_scaling && scaling_mode_ == SOLVE_SCALING)
    {
      // We must undo the scaling of the matrix, so first invert the scaling vectors
      N_VInv(s2inv_, s2inv_);
      N_VInv(s1_, s1_);

      // Unscale the matrix (row_scale_vec_ and col_scale_vec_ point to the same data as s2inv_ and s1_)
      A->GkoMtx()->scale(row_scale_vec_, col_scale_vec_);

      // Fix the scaling vectors
      N_VInv(s2inv_, s2inv_);
      N_VInv(s1_, s1_);
    }

    // Check if any batch entry did not reach the tolerance.
    bool at_least_one_did_not_converge{false};
    sunrealtype max_res_norm{0.0};
    sunrealtype min_res_norm{SUN_BIG_REAL};
    gko_exec_->get_master()->copy_from(gko_exec_, num_batches_,
                                       logger_->get_residual_norm().get_const_data(),
                                       res_norm_array_.get_data());
    auto res_norm = res_norm_array_.get_data();

    for (gko::size_type i = 0; i < num_batches_; i++)
    {
      max_res_norm =
        std::max(max_res_norm,
                 res_norm[i]); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      min_res_norm =
        std::min(min_res_norm,
                 res_norm[i]); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    }
    if (max_res_norm > tol) { at_least_one_did_not_converge = true; }

    // Compute the average number of iterations across all batch entries
    // as well as the maximum and minimum iteration counts.
    gko_exec_->get_master()->copy_from(gko_exec_, num_batches_,
                                       logger_->get_num_iterations().get_const_data(),
                                       iter_count_array_.get_data());
    auto iter_count = iter_count_array_.get_data();

    int max_iter_count{0};
    int min_iter_count{max_iters_};
    avg_iter_count_ = sunrealtype{0.0};
    for (gko::size_type i = 0; i < num_batches_; i++)
    {
      avg_iter_count_ +=
        iter_count[i]; // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      max_iter_count =
        std::max(max_iter_count,
                 iter_count[i]); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      min_iter_count =
        std::min(min_iter_count,
                 iter_count[i]); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    }
    avg_iter_count_ /= static_cast<sunrealtype>(num_batches_);
    sum_of_avg_iters_ += avg_iter_count_;

    // Compute the std. dev. in iteration count across all batch entries.
    // This helps us understand how varied (in difficulty to solve) the entries are.
    stddev_iter_count_ = sunrealtype{0.0};
    for (gko::size_type i = 0; i < num_batches_; i++)
    {
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      stddev_iter_count_ += std::pow(std::abs(iter_count[i] - avg_iter_count_),
                                     2);
    }
    stddev_iter_count_ =
      std::sqrt(stddev_iter_count_ / static_cast<sunrealtype>(num_batches_));

    int retval{0};
    if (at_least_one_did_not_converge) { retval = SUNLS_CONV_FAIL; }
    else { retval = SUN_SUCCESS; }

    SUNLogInfo(sunLogger(), "linear-solver",
               "avg. iter count = " SUN_FORMAT_G
               ", stddev. iter count = " SUN_FORMAT_G ", "
               "max iter count = %d, "
               "min iter count = %d, max res. norm = " SUN_FORMAT_G
               ", min res. "
               "norm = " SUN_FORMAT_G ", tol = " SUN_FORMAT_G ","
               " return code = %d",
               avg_iter_count_, stddev_iter_count_, max_iter_count,
               min_iter_count, max_res_norm, min_res_norm, tol, retval);

    return retval;
  }

private:
  std::shared_ptr<const gko::Executor> gko_exec_;
  gko::batch::stop::tolerance_type tolerance_type_;
  std::shared_ptr<gko::batch::BatchLinOpFactory> precon_factory_;
  std::unique_ptr<typename GkoBatchSolverType::Factory> solver_factory_;
  std::unique_ptr<GkoBatchSolverType> solver_;
  gko::array<sunrealtype> col_scale_vec_;
  gko::array<sunrealtype> row_scale_vec_;
  std::shared_ptr<gko::batch::log::BatchConvergence<sunrealtype>> logger_;
  gko::size_type num_batches_;
  int max_iters_;
  gko::array<sunrealtype> res_norm_array_;
  gko::array<int> iter_count_array_;
  sunrealtype avg_iter_count_;
  sunrealtype sum_of_avg_iters_;
  sunrealtype stddev_iter_count_;
  N_Vector s1_;
  sundials::experimental::N_VectorView s2inv_;
  int scaling_mode_;
  bool scaling_initialized_;
  bool do_setup_;

  void initSUNLinSol(SUNContext sunctx)
  {
    // Attach function pointers for SUNLinearSolver
    using this_type = BatchLinearSolver<GkoBatchSolverType, GkoBatchMatType>;

    object_->content = this;
    object_->ops     = object_ops_.get();
    object_->sunctx  = sunctx;

    object_->ops->gettype = impl::SUNLinSolGetType_GinkgoBatch;
    object_->ops->getid   = impl::SUNLinSolGetID_GinkgoBatch;
    object_->ops->setscalingvectors =
      impl::SUNLinSolSetScalingVectors_GinkgoBatch<this_type>;
    object_->ops->initialize = impl::SUNLinSolInitialize_GinkgoBatch<this_type>;
    object_->ops->setup =
      impl::SUNLinSolSetup_GinkgoBatch<this_type, GkoBatchMatType>;
    object_->ops->solve =
      impl::SUNLinSolSolve_GinkgoBatch<this_type, GkoBatchMatType>;
    object_->ops->free     = impl::SUNLinSolFree_GinkgoBatch<this_type>;
    object_->ops->numiters = impl::SUNLinSolNumIters_GinkgoBatch<this_type>;
  }

  SUNLogger sunLogger()
  {
    SUNLogger log{nullptr};
    SUNContext_GetLogger(object_->sunctx, &log);
    return log;
  }

  SUNContext sunCtx() { return object_->sunctx; }

  static constexpr int default_max_iters_ = 5;
  static constexpr int default_restart_   = 0;
};

} // namespace ginkgo
} // namespace sundials

#endif // SUNLINSOL_GINKGOBATCH_HPP
