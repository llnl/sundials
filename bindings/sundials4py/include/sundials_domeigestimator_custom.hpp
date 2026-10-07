/*-----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *-----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *-----------------------------------------------------------------------------
 * CustomSUNDomEigEstimator: Python implementations of SUNDomEigEstimator.
 *---------------------------------------------------------------------------*/

#ifndef SUNDIALS4PY_DOMEIGESTIMATOR_CUSTOM_HPP
#define SUNDIALS4PY_DOMEIGESTIMATOR_CUSTOM_HPP

#include <cstdint>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <sundials/sundials_classview.hpp>
#include <sundials/sundials_domeigestimator.h>
#include <sundials/sundials_domeigestimator.hpp>
#include <sundials/sundials_errors.h>

#include "sundials4py_core_types.hpp"
#include "sundials4py_custom_object.hpp"

namespace sundials4py {

class CustomSUNDomEigEstimator
{
public:
  enum class HandleState
  {
    unmaterialized,
    materializing,
    materialized
  };

  explicit CustomSUNDomEigEstimator(
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx)
    : sunctx_owner_(std::move(sunctx))
  {}

  std::shared_ptr<std::remove_pointer_t<SUNDomEigEstimator>> _get_sundials_handle(
    nb::handle self)
  {
    if (!sunctx_owner_)
    {
      throw nb::type_error(
        "CustomSUNDomEigEstimator base constructor was not initialized");
    }
    if (state_ == HandleState::materializing)
    {
      throw std::runtime_error(
        "reentrant CustomSUNDomEigEstimator native handle materialization");
    }
    if (state_ == HandleState::materialized) { return estimator_; }

    state_ = HandleState::materializing;
    try
    {
      estimator_ = make_handle(self, sunctx_owner_);
      state_     = HandleState::materialized;
    }
    catch (...)
    {
      estimator_.reset();
      state_ = HandleState::unmaterialized;
      throw;
    }
    return estimator_;
  }

  bool _is_materialized() const { return state_ == HandleState::materialized; }

  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx() const
  {
    return sunctx_owner_;
  }

  static SUNErrCode base_method_status(const char* name)
  {
    PyErr_SetString(PyExc_NotImplementedError, name);
    nb::raise_python_error();
    return SUN_ERR_EXT_FAIL;
  }

private:
  struct Content
  {
    SUNDIALS4PY_CUSTOM_CONTENT_MEMBERS;
  };

  static constexpr const char* label = "CustomSUNDomEigEstimator";

  static std::shared_ptr<std::remove_pointer_t<SUNDomEigEstimator>> make_handle(
    nb::handle impl,
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner)
  {
    if (!method_overridden(impl, "estimate"))
    {
      throw nb::type_error(
        "CustomSUNDomEigEstimator subclass must override estimate()");
    }

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

  static Content* get_content(SUNDomEigEstimator dee)
  {
    if (!dee || !dee->ops || dee->ops->destroy != custom_destroy || !dee->content)
    {
      return nullptr;
    }
    return static_cast<Content*>(dee->content);
  }

public:
  static nb::object _python_object_for(SUNDomEigEstimator dee) noexcept
  {
    try
    {
      Content* content = get_content(dee);
      if (!content) { return nb::object(); }
      nb::object impl = content->weak_impl();
      return impl.is_none() ? nb::object() : impl;
    }
    catch (...)
    {
      return nb::object();
    }
  }

private:
  static nb::object get_impl(SUNDomEigEstimator dee)
  {
    return custom_content_impl(get_content(dee), label);
  }

  static bool method_overridden(nb::handle impl, const char* name)
  {
    return custom_method_overridden<CustomSUNDomEigEstimator>(impl, name);
  }

  template<typename... Args>
  static SUNErrCode call_status(SUNDomEigEstimator dee, const char* name,
                                const char* operation, Args&&... args)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      return static_cast<SUNErrCode>(
        nb::cast<int>(get_impl(dee).attr(name)(std::forward<Args>(args)...)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(dee ? dee->sunctx : nullptr, operation,
                                 SUN_ERR_EXT_FAIL)
  }

  template<typename Value>
  static SUNErrCode call_getter(SUNDomEigEstimator dee, const char* name,
                                Value* value, const char* operation)
  {
    *value = Value{};
    try
    {
      nb::gil_scoped_acquire gil;
      auto result = nb::cast<std::tuple<int, Value>>(get_impl(dee).attr(name)());
      *value = std::get<1>(result);
      return static_cast<SUNErrCode>(std::get<0>(result));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(dee ? dee->sunctx : nullptr, operation,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_set_atimes(SUNDomEigEstimator dee, void* data,
                                      SUNATimesFn fn)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      Content* content = get_content(dee);
      auto state = content->callbacks.install(NativeCallbackSlot::atimes, fn,
                                              data);
      nb::object callback = nb::none();
      if (state)
      {
        callback = nb::module_::import_("functools")
                     .attr("partial")(nb::cpp_function(
                       [state](N_Vector x, N_Vector y) -> int
                       {
                         require_valid_callback(state.get(), "ATimes");
                         return state->fn(state->data, x, y);
                       },
                       nb::arg("x"), nb::arg("y")));
      }
      return static_cast<SUNErrCode>(
        nb::cast<int>(get_impl(dee).attr("set_atimes")(callback)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(dee ? dee->sunctx : nullptr,
                                 "CustomSUNDomEigEstimator.set_atimes",
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_set_rhs(SUNDomEigEstimator dee, void* data, SUNRhsFn fn)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      Content* content = get_content(dee);
      auto state = content->callbacks.install(NativeCallbackSlot::rhs, fn, data);
      nb::object callback = nb::none();
      if (state)
      {
        callback = nb::module_::import_("functools")
                     .attr("partial")(nb::cpp_function(
                       [state](sunrealtype t, N_Vector y, N_Vector ydot) -> int
                       {
                         require_valid_callback(state.get(), "RHS");
                         return state->fn(t, y, ydot, state->data);
                       },
                       nb::arg("t"), nb::arg("y"), nb::arg("ydot")));
      }
      return static_cast<SUNErrCode>(
        nb::cast<int>(get_impl(dee).attr("set_rhs")(callback)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(dee ? dee->sunctx : nullptr,
                                 "CustomSUNDomEigEstimator.set_rhs",
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_set_rhs_linearization_point(SUNDomEigEstimator dee,
                                                       sunrealtype t, N_Vector v)
  {
    return call_status(dee, "set_rhs_linearization_point",
                       "CustomSUNDomEigEstimator.set_rhs_linearization_point",
                       t, nb::cast(v, nb::rv_policy::reference));
  }

  static SUNErrCode custom_set_options(SUNDomEigEstimator dee, const char* id,
                                       const char* file_name, int argc,
                                       char* argv[])
  {
    std::vector<std::string> args;
    args.reserve(argc);
    for (int i = 0; i < argc; ++i) { args.emplace_back(argv[i]); }
    return call_status(dee, "set_options", "CustomSUNDomEigEstimator.set_options",
                       id ? id : "", file_name ? file_name : "", args);
  }

  static SUNErrCode custom_set_max_iters(SUNDomEigEstimator dee, long int value)
  {
    return call_status(dee, "set_max_iters",
                       "CustomSUNDomEigEstimator.set_max_iters", value);
  }

  static SUNErrCode custom_set_num_preprocess_iters(SUNDomEigEstimator dee,
                                                    int value)
  {
    return call_status(dee, "set_num_preprocess_iters",
                       "CustomSUNDomEigEstimator.set_num_preprocess_iters",
                       value);
  }

  static SUNErrCode custom_set_rel_tol(SUNDomEigEstimator dee, sunrealtype value)
  {
    return call_status(dee, "set_rel_tol",
                       "CustomSUNDomEigEstimator.set_rel_tol", value);
  }

  static SUNErrCode custom_set_initial_guess(SUNDomEigEstimator dee, N_Vector q)
  {
    return call_status(dee, "set_initial_guess",
                       "CustomSUNDomEigEstimator.set_initial_guess",
                       nb::cast(q, nb::rv_policy::reference));
  }

  static SUNErrCode custom_initialize(SUNDomEigEstimator dee)
  {
    return call_status(dee, "initialize", "CustomSUNDomEigEstimator.initialize");
  }

  static SUNErrCode custom_estimate(SUNDomEigEstimator dee,
                                    sunrealtype* lambda_real,
                                    sunrealtype* lambda_imag)
  {
    *lambda_real = SUN_RCONST(0.0);
    *lambda_imag = SUN_RCONST(0.0);
    try
    {
      nb::gil_scoped_acquire gil;
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

  static SUNErrCode custom_get_res(SUNDomEigEstimator dee, sunrealtype* value)
  {
    return call_getter(dee, "get_res", value, "CustomSUNDomEigEstimator.get_res");
  }

  static SUNErrCode custom_get_num_iters(SUNDomEigEstimator dee, long int* value)
  {
    return call_getter(dee, "get_num_iters", value,
                       "CustomSUNDomEigEstimator.get_num_iters");
  }

  static SUNErrCode custom_get_num_rhs_evals(SUNDomEigEstimator dee,
                                             long int* value)
  {
    return call_getter(dee, "get_num_rhs_evals", value,
                       "CustomSUNDomEigEstimator.get_num_rhs_evals");
  }

  static SUNErrCode custom_get_num_atimes_calls(SUNDomEigEstimator dee,
                                                long int* value)
  {
    return call_getter(dee, "get_num_atimes_calls", value,
                       "CustomSUNDomEigEstimator.get_num_atimes_calls");
  }

  static SUNErrCode custom_write(SUNDomEigEstimator dee, FILE* outfile)
  {
    return call_status(dee, "write", "CustomSUNDomEigEstimator.write",
                       reinterpret_cast<std::uintptr_t>(outfile));
  }

  static SUNErrCode custom_destroy(SUNDomEigEstimator* dee_ptr)
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

  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner_;
  std::shared_ptr<std::remove_pointer_t<SUNDomEigEstimator>> estimator_;
  HandleState state_{HandleState::unmaterialized};
};

} // namespace sundials4py

#endif // SUNDIALS4PY_DOMEIGESTIMATOR_CUSTOM_HPP
