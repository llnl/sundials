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
 * CustomNVector: Python implementations of the required N_Vector operations.
 *---------------------------------------------------------------------------*/

#ifndef SUNDIALS4PY_NVECTOR_CUSTOM_HPP
#define SUNDIALS4PY_NVECTOR_CUSTOM_HPP

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <sundials/sundials_classview.hpp>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_nvector.hpp>

#include "sundials4py_core_types.hpp"
#include "sundials4py_custom_object.hpp"

namespace sundials4py {

class CustomNVector
  : public CustomObjectBase<CustomNVector, std::remove_pointer_t<N_Vector>>
{
public:
  explicit CustomNVector(std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx)
    : CustomObjectBase(std::move(sunctx), "CustomNVector")
  {}

  void validate(nb::handle self) const { validate_required_methods(self); }

  static void base_method(const char* name)
  {
    PyErr_SetString(PyExc_NotImplementedError, name);
    nb::raise_python_error();
  }

  template<typename T>
  static T base_method_value(const char* name)
  {
    base_method(name);
    return T{};
  }

private:
  enum class Ownership
  {
    weak,
    strong
  };

  struct Content : CustomContentBase
  {
    /* Non-owning pointer used only to clear a SUNDIALS-created clone's
       back-reference before its strong Python reference is released. */
    CustomNVector* owner{nullptr};
    bool has_is_compatible{false};
  };

  static constexpr const char* label = "CustomNVector";

  friend class CustomObjectBase<CustomNVector, std::remove_pointer_t<N_Vector>>;

  std::shared_ptr<std::remove_pointer_t<N_Vector>> make_handle(nb::handle impl);

  static N_Vector create_raw_handle(
    nb::handle impl,
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
    Ownership ownership);

  static N_Vector clone_raw_handle(
    nb::handle impl, N_Vector source,
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner);

  static std::shared_ptr<std::remove_pointer_t<N_Vector>> make_handle(
    nb::handle impl,
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
    Ownership ownership);

  static Content* get_content(N_Vector v)
  {
    if (!v || !v->ops || v->ops->nvdestroy != custom_destroy || !v->content)
    {
      return nullptr;
    }
    return static_cast<Content*>(v->content);
  }

public:
  static nb::object clone_python(nb::handle impl, const char* operation)
  {
    auto* source = nb::cast<CustomNVector*>(impl);
    // The native vtable uses clone() as the fallback when clone_empty() is not
    // overridden. Keep the Python-facing shortcut identical to that behavior.
    const char* method = operation;
    if (std::string(operation) == "clone_empty" &&
        !method_overridden(impl, "clone_empty"))
    {
      method = "clone";
    }
    nb::object cloned_impl = impl.attr(method)();
    if (!nb::isinstance<CustomNVector>(cloned_impl))
    {
      throw nb::type_error((std::string("CustomNVector.") + operation +
                            "() must return a CustomNVector")
                             .c_str());
    }

    auto* cloned = nb::cast<CustomNVector*>(cloned_impl);
    if (!source->sunctx_owner_ || !cloned->sunctx_owner_)
    {
      throw nb::type_error(
        "CustomNVector clone did not initialize the base constructor");
    }
    require_same_context(source->sunctx_owner_.get(), cloned->sunctx_owner_,
                         (std::string("CustomNVector.") + operation).c_str());
    return cloned_impl;
  }

  static nb::object _python_object_for(N_Vector v) noexcept
  {
    try
    {
      Content* content = get_content(v);
      if (!content) { return nb::object(); }
      if (content->strong_impl.is_valid()) { return content->strong_impl; }
      nb::object impl = content->weak_impl();
      return impl.is_none() ? nb::object() : impl;
    }
    catch (...)
    {
      return nb::object();
    }
  }

private:
  static nb::object get_impl(N_Vector v)
  { return custom_content_impl(get_content(v), label); }

  static bool method_overridden(nb::handle impl, const char* name)
  { return custom_method_overridden<CustomNVector>(impl, name); }

  static void validate_required_methods(nb::handle impl)
  {
    static constexpr const char* required[] = {"clone",       "get_length",
                                               "linear_sum",  "const",
                                               "prod",        "div",
                                               "scale",       "abs",
                                               "inv",         "add_const",
                                               "dot_prod",    "max_norm",
                                               "wrms_norm",   "wrms_norm_mask",
                                               "min",         "wl2_norm",
                                               "l1_norm",     "compare",
                                               "inv_test",    "constr_mask",
                                               "min_quotient"};
    std::string missing;
    for (const char* name : required)
    {
      if (method_overridden(impl, name)) { continue; }
      if (!missing.empty()) { missing += ", "; }
      missing += name;
    }
    if (!missing.empty())
    {
      nb::object type_name =
        nb::borrow<nb::object>(reinterpret_cast<PyObject*>(Py_TYPE(impl.ptr())))
          .attr("__name__");
      throw nb::type_error((std::string("CustomNVector subclass ") +
                            nb::cast<std::string>(type_name) +
                            " is missing required operations: " + missing)
                             .c_str());
    }
  }

  static nb::object arg(N_Vector v)
  { return nb::cast(v, nb::rv_policy::reference); }

  static nb::object operand(nb::handle self, N_Vector owner, N_Vector v,
                            const char* operation)
  {
    nb::object value = arg(v);
    if (Py_TYPE(value.ptr()) == Py_TYPE(self.ptr())) { return value; }

    Content* self_content = get_content(owner);
    if (self_content && self_content->has_is_compatible &&
        nb::cast<bool>(self.attr("is_compatible")(value)))
    {
      return value;
    }
    nb::object value_type =
      nb::borrow<nb::object>(reinterpret_cast<PyObject*>(Py_TYPE(value.ptr())));
    nb::object self_type =
      nb::borrow<nb::object>(reinterpret_cast<PyObject*>(Py_TYPE(self.ptr())));
    throw nb::type_error(
      (std::string(operation) + ": operand is a different vector type (" +
       nb::cast<std::string>(value_type.attr("__name__")) + ") than " +
       nb::cast<std::string>(self_type.attr("__name__")))
        .c_str());
  }

  template<typename Value, typename Invoke>
  static Value call_value_with(N_Vector v, const char* operation, Value failure,
                               Invoke&& invoke)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      if (custom_exception_pending()) { return failure; }
      nb::object self = get_impl(v);
      return nb::cast<Value>(invoke(self));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(v ? v->sunctx : nullptr, operation, failure)
  }

  template<typename Invoke>
  static void call_void_with(N_Vector v, const char* operation, Invoke&& invoke)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      if (custom_exception_pending()) { return; }
      invoke(get_impl(v));
    }
    catch (const std::exception& error)
    {
      report_custom_exception(v ? v->sunctx : nullptr, operation, error,
                              __FILE__, __LINE__);
    }
    catch (...)
    {
      report_custom_unknown_exception(v ? v->sunctx : nullptr, operation,
                                      __FILE__, __LINE__);
    }
  }

  template<typename Value, typename... Args>
  static Value call_value(N_Vector v, const char* name, const char* operation,
                          Value failure, Args&&... args)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      if (custom_exception_pending()) { return failure; }
      return nb::cast<Value>(get_impl(v).attr(name)(std::forward<Args>(args)...));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(v ? v->sunctx : nullptr, operation, failure)
  }

  template<typename... Args>
  static void call_void(N_Vector v, const char* name, const char* operation,
                        Args&&... args)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      if (custom_exception_pending()) { return; }
      get_impl(v).attr(name)(std::forward<Args>(args)...);
    }
    catch (const std::exception& error)
    {
      report_custom_exception(v ? v->sunctx : nullptr, operation, error,
                              __FILE__, __LINE__);
    }
    catch (...)
    {
      report_custom_unknown_exception(v ? v->sunctx : nullptr, operation,
                                      __FILE__, __LINE__);
    }
  }

  template<typename... Args>
  static SUNErrCode call_status(N_Vector v, const char* name,
                                const char* operation, Args&&... args)
  {
    return static_cast<SUNErrCode>(
      call_value<int>(v, name, operation, static_cast<int>(SUN_ERR_EXT_FAIL),
                      std::forward<Args>(args)...));
  }

  template<typename Invoke>
  static SUNErrCode call_status_with(N_Vector v, const char* operation,
                                     Invoke&& invoke)
  {
    return static_cast<SUNErrCode>(
      call_value_with<int>(v, operation, static_cast<int>(SUN_ERR_EXT_FAIL),
                           std::forward<Invoke>(invoke)));
  }

  static N_Vector_ID custom_get_vector_id(N_Vector);

  static N_Vector clone_impl(N_Vector v, const char* name, const char* operation);

  static N_Vector custom_clone(N_Vector v);

  static N_Vector custom_clone_empty(N_Vector v);

  static void custom_destroy(N_Vector v);

  static void custom_space(N_Vector v, sunindextype* lrw, sunindextype* liw);

  static sunrealtype* custom_get_array_pointer(N_Vector v);

  static sunrealtype* custom_get_device_array_pointer(N_Vector v);

  static void custom_set_array_pointer(sunrealtype* ptr, N_Vector v);

  static void custom_set_device_array_pointer(sunrealtype* ptr, N_Vector v);

  static SUNComm custom_get_communicator(N_Vector v);

  static sunindextype custom_get_length(N_Vector v);

  static sunindextype custom_get_local_length(N_Vector v);

  static void custom_linear_sum(sunrealtype a, N_Vector x, sunrealtype b,
                                N_Vector y, N_Vector z);

  static void custom_const(sunrealtype c, N_Vector z);

  static void custom_prod(N_Vector x, N_Vector y, N_Vector z);

  static void custom_div(N_Vector x, N_Vector y, N_Vector z);

  static void custom_scale(sunrealtype c, N_Vector x, N_Vector z);

  static void custom_abs(N_Vector x, N_Vector z);

  static void custom_inv(N_Vector x, N_Vector z);

  static void custom_add_const(N_Vector x, sunrealtype b, N_Vector z);

  static sunrealtype custom_dot_prod(N_Vector x, N_Vector y);

  static sunrealtype custom_max_norm(N_Vector x);

  static sunrealtype custom_wrms_norm(N_Vector x, N_Vector w);

  static sunrealtype custom_wrms_norm_mask(N_Vector x, N_Vector w, N_Vector id);

  static sunrealtype custom_min(N_Vector x);

  static sunrealtype custom_wl2_norm(N_Vector x, N_Vector w);

  static sunrealtype custom_l1_norm(N_Vector x);

  static void custom_compare(sunrealtype c, N_Vector x, N_Vector z);

  static sunbooleantype custom_inv_test(N_Vector x, N_Vector z);

  static sunbooleantype custom_constr_mask(N_Vector c, N_Vector x, N_Vector m);

  static sunrealtype custom_min_quotient(N_Vector num, N_Vector denom);

  static nb::list vector_list(nb::handle self, N_Vector owner, N_Vector* vectors,
                              int count, const char* operation)
  {
    nb::list result;
    for (int i = 0; i < count; ++i)
    {
      result.append(operand(self, owner, vectors[i], operation));
    }
    return result;
  }

  static nb::list vector_lists(nb::handle self, N_Vector owner,
                               N_Vector** vectors, int outer_count,
                               int inner_count, const char* operation)
  {
    nb::list result;
    for (int i = 0; i < outer_count; ++i)
    {
      result.append(vector_list(self, owner, vectors[i], inner_count, operation));
    }
    return result;
  }

  static std::vector<sunrealtype> real_values(const sunrealtype* values, int count)
  { return std::vector<sunrealtype>(values, values + count); }

  static SUNErrCode custom_linear_combination(int nvec, sunrealtype* c,
                                              N_Vector* X, N_Vector z);

  static SUNErrCode custom_scale_add_multi(int nvec, sunrealtype* a, N_Vector x,
                                           N_Vector* Y, N_Vector* Z);

  static SUNErrCode custom_dot_prod_multi(int nvec, N_Vector x, N_Vector* Y,
                                          sunrealtype* dots);

  static SUNErrCode custom_linear_sum_vector_array(int nvec, sunrealtype a,
                                                   N_Vector* X, sunrealtype b,
                                                   N_Vector* Y, N_Vector* Z);

  static SUNErrCode custom_scale_vector_array(int nvec, sunrealtype* c,
                                              N_Vector* X, N_Vector* Z);

  static SUNErrCode custom_const_vector_array(int nvec, sunrealtype c,
                                              N_Vector* Z);

  template<typename Invoke>
  static SUNErrCode vector_array_reduction(N_Vector owner, int nvec,
                                           sunrealtype* values,
                                           sunrealtype failure,
                                           const char* operation, Invoke&& invoke)
  {
    auto fail = [&]
    {
      if (values && nvec > 0) { std::fill(values, values + nvec, failure); }
    };
    try
    {
      nb::gil_scoped_acquire gil;
      if (custom_exception_pending())
      {
        fail();
        return SUN_ERR_EXT_FAIL;
      }
      auto result = nb::cast<std::tuple<int, std::vector<sunrealtype>>>(
        invoke(get_impl(owner)));
      const auto& output = std::get<1>(result);
      if (static_cast<int>(output.size()) != nvec)
      {
        throw nb::value_error(
          "vector-array reduction returned the wrong result size");
      }
      std::copy(output.begin(), output.end(), values);
      return static_cast<SUNErrCode>(std::get<0>(result));
    }
    catch (const std::exception& error)
    {
      fail();
      report_custom_exception(owner ? owner->sunctx : nullptr, operation, error,
                              __FILE__, __LINE__);
      return SUN_ERR_EXT_FAIL;
    }
    catch (...)
    {
      fail();
      report_custom_unknown_exception(owner ? owner->sunctx : nullptr,
                                      operation, __FILE__, __LINE__);
      return SUN_ERR_EXT_FAIL;
    }
  }

  static SUNErrCode custom_wrms_norm_vector_array(int nvec, N_Vector* X,
                                                  N_Vector* W,
                                                  sunrealtype* norms);

  static SUNErrCode custom_wrms_norm_mask_vector_array(int nvec, N_Vector* X,
                                                       N_Vector* W, N_Vector id,
                                                       sunrealtype* norms);

  static SUNErrCode custom_scale_add_multi_vector_array(
    int nvec, int nsum, sunrealtype* a, N_Vector* X, N_Vector** Y, N_Vector** Z);

  static SUNErrCode custom_linear_combination_vector_array(int nvec, int nsum,
                                                           sunrealtype* c,
                                                           N_Vector** X,
                                                           N_Vector* Z);

  static sunrealtype custom_dot_prod_local(N_Vector x, N_Vector y);

  static sunrealtype custom_max_norm_local(N_Vector x);

  static sunrealtype custom_min_local(N_Vector x);

  static sunrealtype custom_l1_norm_local(N_Vector x);

  static sunbooleantype custom_inv_test_local(N_Vector x, N_Vector z);

  static sunbooleantype custom_constr_mask_local(N_Vector c, N_Vector x,
                                                 N_Vector m);

  static sunrealtype custom_min_quotient_local(N_Vector num, N_Vector denom);

  static sunrealtype custom_wsqrsum_local(N_Vector x, N_Vector w);

  static sunrealtype custom_wsqrsum_mask_local(N_Vector x, N_Vector w,
                                               N_Vector id);

  static SUNErrCode custom_dot_prod_multi_local(int nvec, N_Vector x,
                                                N_Vector* Y, sunrealtype* dots);

  static SUNErrCode custom_dot_prod_multi_all_reduce(int nvec, N_Vector x,
                                                     sunrealtype* dots);

  static SUNErrCode custom_buf_size(N_Vector x, sunindextype* size);

  static SUNErrCode custom_buf_pack(N_Vector x, void* buffer);

  static SUNErrCode custom_buf_unpack(N_Vector x, void* buffer);

  static void custom_print(N_Vector x);

  static void custom_print_file(N_Vector x, FILE* outfile);
};

} // namespace sundials4py

#endif // SUNDIALS4PY_NVECTOR_CUSTOM_HPP
