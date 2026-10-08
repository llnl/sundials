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
  : public CustomObjectBase<CustomSUNDomEigEstimator,
                            std::remove_pointer_t<SUNDomEigEstimator>>
{
public:
  explicit CustomSUNDomEigEstimator(
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx)
    : CustomObjectBase(std::move(sunctx), "CustomSUNDomEigEstimator")
  {}

  void validate(nb::handle self) const { validate_required_methods(self); }

  static SUNErrCode base_method_status(const char* name)
  {
    PyErr_SetString(PyExc_NotImplementedError, name);
    nb::raise_python_error();
    return SUN_ERR_EXT_FAIL;
  }

private:
  struct Content : CustomContentBase
  {};

  static constexpr const char* label = "CustomSUNDomEigEstimator";

  friend class CustomObjectBase<CustomSUNDomEigEstimator,
                                std::remove_pointer_t<SUNDomEigEstimator>>;

  std::shared_ptr<std::remove_pointer_t<SUNDomEigEstimator>> make_handle(
    nb::handle impl);

  static std::shared_ptr<std::remove_pointer_t<SUNDomEigEstimator>> make_handle(
    nb::handle impl,
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner);

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
  { return custom_content_impl(get_content(dee), label); }

  static bool method_overridden(nb::handle impl, const char* name)
  { return custom_method_overridden<CustomSUNDomEigEstimator>(impl, name); }

  static void validate_required_methods(nb::handle impl)
  {
    if (!method_overridden(impl, "estimate"))
    {
      throw nb::type_error(
        "CustomSUNDomEigEstimator subclass must override estimate()");
    }
  }

  template<typename... Args>
  static SUNErrCode call_status(SUNDomEigEstimator dee, const char* name,
                                const char* operation, Args&&... args)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
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
      if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
      auto result = nb::cast<std::tuple<int, Value>>(get_impl(dee).attr(name)());
      *value = std::get<1>(result);
      return static_cast<SUNErrCode>(std::get<0>(result));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(dee ? dee->sunctx : nullptr, operation,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_set_atimes(SUNDomEigEstimator dee, void* data,
                                      SUNATimesFn fn);

  static SUNErrCode custom_set_rhs(SUNDomEigEstimator dee, void* data,
                                   SUNRhsFn fn);

  static SUNErrCode custom_set_rhs_linearization_point(SUNDomEigEstimator dee,
                                                       sunrealtype t, N_Vector v);

  static SUNErrCode custom_set_options(SUNDomEigEstimator dee, const char* id,
                                       const char* file_name, int argc,
                                       char* argv[]);

  static SUNErrCode custom_set_max_iters(SUNDomEigEstimator dee, long int value);

  static SUNErrCode custom_set_num_preprocess_iters(SUNDomEigEstimator dee,
                                                    int value);

  static SUNErrCode custom_set_rel_tol(SUNDomEigEstimator dee, sunrealtype value);

  static SUNErrCode custom_set_initial_guess(SUNDomEigEstimator dee, N_Vector q);

  static SUNErrCode custom_initialize(SUNDomEigEstimator dee);

  static SUNErrCode custom_estimate(SUNDomEigEstimator dee,
                                    sunrealtype* lambda_real,
                                    sunrealtype* lambda_imag);

  static SUNErrCode custom_get_res(SUNDomEigEstimator dee, sunrealtype* value);

  static SUNErrCode custom_get_num_iters(SUNDomEigEstimator dee, long int* value);

  static SUNErrCode custom_get_num_rhs_evals(SUNDomEigEstimator dee,
                                             long int* value);

  static SUNErrCode custom_get_num_atimes_calls(SUNDomEigEstimator dee,
                                                long int* value);

  static SUNErrCode custom_write(SUNDomEigEstimator dee, FILE* outfile);

  static SUNErrCode custom_destroy(SUNDomEigEstimator* dee_ptr);
};

} // namespace sundials4py

#endif // SUNDIALS4PY_DOMEIGESTIMATOR_CUSTOM_HPP
