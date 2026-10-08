/*-----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *----------------------------------------------------------------------------*/

#include "sundials4py.hpp"

namespace sundials4py {

std::shared_ptr<std::remove_pointer_t<SUNDomEigEstimator>> CustomSUNDomEigEstimator::make_handle(
  nb::handle impl)
{ return make_handle(impl, sunctx_owner_); }

std::shared_ptr<std::remove_pointer_t<SUNDomEigEstimator>> CustomSUNDomEigEstimator::make_handle(
  nb::handle impl, std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner)
{
  validate_required_methods(impl);

  const bool has_set_atimes = method_overridden(impl, "set_atimes");
  const bool has_set_rhs    = method_overridden(impl, "set_rhs");
  const bool has_set_rhs_point =
    method_overridden(impl, "set_rhs_linearization_point");
  const bool has_set_options   = method_overridden(impl, "set_options");
  const bool has_set_max_iters = method_overridden(impl, "set_max_iters");
  const bool has_set_num_preprocess =
    method_overridden(impl, "set_num_preprocess_iters");
  const bool has_set_rel_tol       = method_overridden(impl, "set_rel_tol");
  const bool has_set_initial_guess = method_overridden(impl,
                                                       "set_initial_guess");
  const bool has_initialize        = method_overridden(impl, "initialize");
  const bool has_get_res           = method_overridden(impl, "get_res");
  const bool has_get_num_iters     = method_overridden(impl, "get_num_iters");
  const bool has_get_num_rhs_evals = method_overridden(impl,
                                                       "get_num_rhs_evals");
  const bool has_get_num_atimes    = method_overridden(impl,
                                                       "get_num_atimes_calls");
  const bool has_write             = method_overridden(impl, "write");

  std::unique_ptr<Content> content(new Content);
  content->weak_impl    = nb::weakref(impl);
  content->sunctx_owner = std::move(sunctx_owner);

  SUNDomEigEstimator dee =
    SUNDomEigEstimator_NewEmpty(content->sunctx_owner.get());
  if (!dee) { throw error_returned("SUNDomEigEstimator_NewEmpty failed"); }

  dee->content       = content.release();
  dee->ops->estimate = custom_estimate;
  dee->ops->destroy  = custom_destroy;
  if (has_set_atimes) { dee->ops->setatimes = custom_set_atimes; }
  if (has_set_rhs) { dee->ops->setrhs = custom_set_rhs; }
  if (has_set_rhs_point)
  {
    dee->ops->setrhslinearizationpoint = custom_set_rhs_linearization_point;
  }
  if (has_set_options) { dee->ops->setoptions = custom_set_options; }
  if (has_set_max_iters) { dee->ops->setmaxiters = custom_set_max_iters; }
  if (has_set_num_preprocess)
  {
    dee->ops->setnumpreprocessiters = custom_set_num_preprocess_iters;
  }
  if (has_set_rel_tol) { dee->ops->setreltol = custom_set_rel_tol; }
  if (has_set_initial_guess)
  {
    dee->ops->setinitialguess = custom_set_initial_guess;
  }
  if (has_initialize) { dee->ops->initialize = custom_initialize; }
  if (has_get_res) { dee->ops->getres = custom_get_res; }
  if (has_get_num_iters) { dee->ops->getnumiters = custom_get_num_iters; }
  if (has_get_num_rhs_evals)
  {
    dee->ops->getnumrhsevals = custom_get_num_rhs_evals;
  }
  if (has_get_num_atimes)
  {
    dee->ops->getnumatimescalls = custom_get_num_atimes_calls;
  }
  if (has_write) { dee->ops->write = custom_write; }

  return sundials::experimental::our_make_shared<
    std::remove_pointer_t<SUNDomEigEstimator>,
    sundials::experimental::SUNDomEigEstimatorDeleter>(dee);
}

SUNErrCode CustomSUNDomEigEstimator::custom_set_atimes(SUNDomEigEstimator dee,
                                                       void* data, SUNATimesFn fn)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    Content* content = get_content(dee);
    if (!content) { throw nb::type_error("invalid custom estimator content"); }
    auto state = NativeCallbackRegistry::prepare(fn, data);
    PreparedCallback prepared(state);
    nb::object callback = nb::none();
    if (state)
    {
      callback = nb::cpp_function(
        [state](N_Vector x, N_Vector y) -> int
        {
          require_valid_callback(state.get(), "ATimes");
          return state->fn(state->data, x, y);
        },
        nb::arg("x"), nb::arg("y"));
    }
    SUNErrCode status = static_cast<SUNErrCode>(
      nb::cast<int>(get_impl(dee).attr("set_atimes")(callback)));
    if (status == SUN_SUCCESS)
    {
      prepared.commit(content->callbacks, NativeCallbackSlot::atimes);
    }
    return status;
  }
  SUNDIALS4PY_CATCH_AND_REPORT(dee ? dee->sunctx : nullptr,
                               "CustomSUNDomEigEstimator.set_atimes",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNDomEigEstimator::custom_set_rhs(SUNDomEigEstimator dee,
                                                    void* data, SUNRhsFn fn)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    Content* content = get_content(dee);
    if (!content) { throw nb::type_error("invalid custom estimator content"); }
    auto state = NativeCallbackRegistry::prepare(fn, data);
    PreparedCallback prepared(state);
    nb::object callback = nb::none();
    if (state)
    {
      callback = nb::cpp_function(
        [state](sunrealtype t, N_Vector y, N_Vector ydot) -> int
        {
          require_valid_callback(state.get(), "RHS");
          return state->fn(t, y, ydot, state->data);
        },
        nb::arg("t"), nb::arg("y"), nb::arg("ydot"));
    }
    SUNErrCode status = static_cast<SUNErrCode>(
      nb::cast<int>(get_impl(dee).attr("set_rhs")(callback)));
    if (status == SUN_SUCCESS)
    {
      prepared.commit(content->callbacks, NativeCallbackSlot::rhs);
    }
    return status;
  }
  SUNDIALS4PY_CATCH_AND_REPORT(dee ? dee->sunctx : nullptr,
                               "CustomSUNDomEigEstimator.set_rhs",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNDomEigEstimator::custom_set_rhs_linearization_point(
  SUNDomEigEstimator dee, sunrealtype t, N_Vector v)
{
  return call_status(dee, "set_rhs_linearization_point",
                     "CustomSUNDomEigEstimator.set_rhs_linearization_point", t,
                     nb::cast(v, nb::rv_policy::reference));
}

SUNErrCode CustomSUNDomEigEstimator::custom_set_options(SUNDomEigEstimator dee,
                                                        const char* id,
                                                        const char* file_name,
                                                        int argc, char* argv[])
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    std::vector<std::string> args;
    args.reserve(static_cast<size_t>(argc));
    for (int i = 0; i < argc; ++i) { args.emplace_back(argv[i]); }
    return static_cast<SUNErrCode>(nb::cast<int>(get_impl(dee).attr(
      "set_options")(id ? id : "", file_name ? file_name : "", args)));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(dee ? dee->sunctx : nullptr,
                               "CustomSUNDomEigEstimator.set_options",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNDomEigEstimator::custom_set_max_iters(SUNDomEigEstimator dee,
                                                          long int value)
{
  return call_status(dee, "set_max_iters",
                     "CustomSUNDomEigEstimator.set_max_iters", value);
}

SUNErrCode CustomSUNDomEigEstimator::custom_set_num_preprocess_iters(
  SUNDomEigEstimator dee, int value)
{
  return call_status(dee, "set_num_preprocess_iters",
                     "CustomSUNDomEigEstimator.set_num_preprocess_iters", value);
}

SUNErrCode CustomSUNDomEigEstimator::custom_set_rel_tol(SUNDomEigEstimator dee,
                                                        sunrealtype value)
{
  return call_status(dee, "set_rel_tol", "CustomSUNDomEigEstimator.set_rel_tol",
                     value);
}

SUNErrCode CustomSUNDomEigEstimator::custom_set_initial_guess(SUNDomEigEstimator dee,
                                                              N_Vector q)
{
  return call_status(dee, "set_initial_guess",
                     "CustomSUNDomEigEstimator.set_initial_guess",
                     nb::cast(q, nb::rv_policy::reference));
}

SUNErrCode CustomSUNDomEigEstimator::custom_initialize(SUNDomEigEstimator dee)
{
  return call_status(dee, "initialize", "CustomSUNDomEigEstimator.initialize");
}

SUNErrCode CustomSUNDomEigEstimator::custom_estimate(SUNDomEigEstimator dee,
                                                     sunrealtype* lambda_real,
                                                     sunrealtype* lambda_imag)
{
  *lambda_real = SUN_RCONST(0.0);
  *lambda_imag = SUN_RCONST(0.0);
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    auto result = nb::cast<std::tuple<int, sunrealtype, sunrealtype>>(
      get_impl(dee).attr("estimate")());
    *lambda_real = std::get<1>(result);
    *lambda_imag = std::get<2>(result);
    return static_cast<SUNErrCode>(std::get<0>(result));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(dee ? dee->sunctx : nullptr,
                               "CustomSUNDomEigEstimator.estimate",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNDomEigEstimator::custom_get_res(SUNDomEigEstimator dee,
                                                    sunrealtype* value)
{
  return call_getter(dee, "get_res", value, "CustomSUNDomEigEstimator.get_res");
}

SUNErrCode CustomSUNDomEigEstimator::custom_get_num_iters(SUNDomEigEstimator dee,
                                                          long int* value)
{
  return call_getter(dee, "get_num_iters", value,
                     "CustomSUNDomEigEstimator.get_num_iters");
}

SUNErrCode CustomSUNDomEigEstimator::custom_get_num_rhs_evals(SUNDomEigEstimator dee,
                                                              long int* value)
{
  return call_getter(dee, "get_num_rhs_evals", value,
                     "CustomSUNDomEigEstimator.get_num_rhs_evals");
}

SUNErrCode CustomSUNDomEigEstimator::custom_get_num_atimes_calls(
  SUNDomEigEstimator dee, long int* value)
{
  return call_getter(dee, "get_num_atimes_calls", value,
                     "CustomSUNDomEigEstimator.get_num_atimes_calls");
}

SUNErrCode CustomSUNDomEigEstimator::custom_write(SUNDomEigEstimator dee,
                                                  FILE* outfile)
{
  return call_status(dee, "write", "CustomSUNDomEigEstimator.write",
                     reinterpret_cast<std::uintptr_t>(outfile));
}

SUNErrCode CustomSUNDomEigEstimator::custom_destroy(SUNDomEigEstimator* dee_ptr)
{
  if (!dee_ptr || !*dee_ptr) { return SUN_SUCCESS; }
  SUNDomEigEstimator dee = *dee_ptr;
  Content* content       = get_content(dee);
  custom_content_destroy(content);
  dee->content = nullptr;
  SUNDomEigEstimator_FreeEmpty(dee);
  *dee_ptr = nullptr;
  return SUN_SUCCESS;
}

} // namespace sundials4py
