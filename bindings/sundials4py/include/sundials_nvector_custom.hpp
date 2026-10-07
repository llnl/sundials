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
#include <memory>
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
{
public:
  enum class HandleState
  {
    unmaterialized,
    materializing,
    materialized
  };

  explicit CustomNVector(std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx)
    : sunctx_owner_(std::move(sunctx))
  {}

  std::shared_ptr<std::remove_pointer_t<N_Vector>> _get_sundials_handle(
    nb::handle self)
  {
    if (!sunctx_owner_)
    {
      throw nb::type_error(
        "CustomNVector base constructor was not initialized");
    }
    if (state_ == HandleState::materializing)
    {
      throw std::runtime_error(
        "reentrant CustomNVector native handle materialization");
    }
    if (state_ == HandleState::materialized) { return vector_; }

    state_ = HandleState::materializing;
    try
    {
      vector_ = make_handle(self, sunctx_owner_, Ownership::weak);
      state_  = HandleState::materialized;
    }
    catch (...)
    {
      vector_.reset();
      state_ = HandleState::unmaterialized;
      throw;
    }
    return vector_;
  }

  bool _is_materialized() const { return state_ == HandleState::materialized; }

  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx() const
  {
    return sunctx_owner_;
  }

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

  struct Content
  {
    SUNDIALS4PY_CUSTOM_CONTENT_MEMBERS;
  };

  static constexpr const char* label = "CustomNVector";

  static N_Vector create_raw_handle(
    nb::handle impl,
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
    Ownership ownership)
  {
    validate_required_methods(impl);

    const bool has_clone_empty = method_overridden(impl, "clone_empty");
    const bool has_get_array   = method_overridden(impl, "get_array_pointer");
    const bool has_get_device  = method_overridden(impl,
                                                   "get_device_array_pointer");
    const bool has_set_array   = method_overridden(impl, "set_array_pointer");
    const bool has_set_device  = method_overridden(impl,
                                                   "set_device_array_pointer");
    const bool has_get_comm    = method_overridden(impl, "get_communicator");
    const bool has_get_local_length = method_overridden(impl,
                                                        "get_local_length");
#define SUNDIALS4PY_NVEC_HAS(NAME) \
  const bool has_##NAME = method_overridden(impl, #NAME)
    SUNDIALS4PY_NVEC_HAS(linear_combination);
    SUNDIALS4PY_NVEC_HAS(scale_add_multi);
    SUNDIALS4PY_NVEC_HAS(dot_prod_multi);
    SUNDIALS4PY_NVEC_HAS(linear_sum_vector_array);
    SUNDIALS4PY_NVEC_HAS(scale_vector_array);
    SUNDIALS4PY_NVEC_HAS(const_vector_array);
    SUNDIALS4PY_NVEC_HAS(wrms_norm_vector_array);
    SUNDIALS4PY_NVEC_HAS(wrms_norm_mask_vector_array);
    SUNDIALS4PY_NVEC_HAS(scale_add_multi_vector_array);
    SUNDIALS4PY_NVEC_HAS(linear_combination_vector_array);
    SUNDIALS4PY_NVEC_HAS(dot_prod_local);
    SUNDIALS4PY_NVEC_HAS(max_norm_local);
    SUNDIALS4PY_NVEC_HAS(min_local);
    SUNDIALS4PY_NVEC_HAS(l1_norm_local);
    SUNDIALS4PY_NVEC_HAS(inv_test_local);
    SUNDIALS4PY_NVEC_HAS(constr_mask_local);
    SUNDIALS4PY_NVEC_HAS(min_quotient_local);
    SUNDIALS4PY_NVEC_HAS(wsqrsum_local);
    SUNDIALS4PY_NVEC_HAS(wsqrsum_mask_local);
    SUNDIALS4PY_NVEC_HAS(dot_prod_multi_local);
    SUNDIALS4PY_NVEC_HAS(dot_prod_multi_all_reduce);
    SUNDIALS4PY_NVEC_HAS(buf_size);
    SUNDIALS4PY_NVEC_HAS(buf_pack);
    SUNDIALS4PY_NVEC_HAS(buf_unpack);
    SUNDIALS4PY_NVEC_HAS(print);
    SUNDIALS4PY_NVEC_HAS(print_file);
#undef SUNDIALS4PY_NVEC_HAS

    std::unique_ptr<Content> content(new Content);
    content->sunctx_owner = std::move(sunctx_owner);
    if (ownership == Ownership::weak)
    {
      content->weak_impl = nb::weakref(impl);
    }
    else { content->strong_impl = nb::borrow<nb::object>(impl); }

    N_Vector v = N_VNewEmpty(content->sunctx_owner.get());
    if (!v) { throw error_returned("N_VNewEmpty failed"); }

    v->content            = content.release();
    v->ops->nvgetvectorid = custom_get_vector_id;
    v->ops->nvclone       = custom_clone;
    v->ops->nvcloneempty  = has_clone_empty ? custom_clone_empty : custom_clone;
    v->ops->nvdestroy     = custom_destroy;
    v->ops->nvspace       = custom_space;
    v->ops->nvgetlength   = custom_get_length;
    v->ops->nvgetlocallength = has_get_local_length ? custom_get_local_length
                                                    : custom_get_length;
    v->ops->nvlinearsum      = custom_linear_sum;
    v->ops->nvconst          = custom_const;
    v->ops->nvprod           = custom_prod;
    v->ops->nvdiv            = custom_div;
    v->ops->nvscale          = custom_scale;
    v->ops->nvabs            = custom_abs;
    v->ops->nvinv            = custom_inv;
    v->ops->nvaddconst       = custom_add_const;
    v->ops->nvdotprod        = custom_dot_prod;
    v->ops->nvmaxnorm        = custom_max_norm;
    v->ops->nvwrmsnorm       = custom_wrms_norm;
    v->ops->nvwrmsnormmask   = custom_wrms_norm_mask;
    v->ops->nvmin            = custom_min;
    v->ops->nvwl2norm        = custom_wl2_norm;
    v->ops->nvl1norm         = custom_l1_norm;
    v->ops->nvcompare        = custom_compare;
    v->ops->nvinvtest        = custom_inv_test;
    v->ops->nvconstrmask     = custom_constr_mask;
    v->ops->nvminquotient    = custom_min_quotient;

    if (has_get_array) { v->ops->nvgetarraypointer = custom_get_array_pointer; }
    if (has_get_device)
    {
      v->ops->nvgetdevicearraypointer = custom_get_device_array_pointer;
    }
    if (has_set_array) { v->ops->nvsetarraypointer = custom_set_array_pointer; }
    if (has_set_device)
    {
      v->ops->nvsetdevicearraypointer = custom_set_device_array_pointer;
    }
    if (has_get_comm) { v->ops->nvgetcommunicator = custom_get_communicator; }

#define SUNDIALS4PY_NVEC_SET_OP(NAME, FIELD) \
  if (has_##NAME) { v->ops->FIELD = custom_##NAME; }
    SUNDIALS4PY_NVEC_SET_OP(linear_combination, nvlinearcombination)
    SUNDIALS4PY_NVEC_SET_OP(scale_add_multi, nvscaleaddmulti)
    SUNDIALS4PY_NVEC_SET_OP(dot_prod_multi, nvdotprodmulti)
    SUNDIALS4PY_NVEC_SET_OP(linear_sum_vector_array, nvlinearsumvectorarray)
    SUNDIALS4PY_NVEC_SET_OP(scale_vector_array, nvscalevectorarray)
    SUNDIALS4PY_NVEC_SET_OP(const_vector_array, nvconstvectorarray)
    SUNDIALS4PY_NVEC_SET_OP(wrms_norm_vector_array, nvwrmsnormvectorarray)
    SUNDIALS4PY_NVEC_SET_OP(wrms_norm_mask_vector_array, nvwrmsnormmaskvectorarray)
    SUNDIALS4PY_NVEC_SET_OP(scale_add_multi_vector_array,
                            nvscaleaddmultivectorarray)
    SUNDIALS4PY_NVEC_SET_OP(linear_combination_vector_array,
                            nvlinearcombinationvectorarray)
    SUNDIALS4PY_NVEC_SET_OP(dot_prod_local, nvdotprodlocal)
    SUNDIALS4PY_NVEC_SET_OP(max_norm_local, nvmaxnormlocal)
    SUNDIALS4PY_NVEC_SET_OP(min_local, nvminlocal)
    SUNDIALS4PY_NVEC_SET_OP(l1_norm_local, nvl1normlocal)
    SUNDIALS4PY_NVEC_SET_OP(inv_test_local, nvinvtestlocal)
    SUNDIALS4PY_NVEC_SET_OP(constr_mask_local, nvconstrmasklocal)
    SUNDIALS4PY_NVEC_SET_OP(min_quotient_local, nvminquotientlocal)
    SUNDIALS4PY_NVEC_SET_OP(wsqrsum_local, nvwsqrsumlocal)
    SUNDIALS4PY_NVEC_SET_OP(wsqrsum_mask_local, nvwsqrsummasklocal)
    SUNDIALS4PY_NVEC_SET_OP(dot_prod_multi_local, nvdotprodmultilocal)
    SUNDIALS4PY_NVEC_SET_OP(dot_prod_multi_all_reduce, nvdotprodmultiallreduce)
    SUNDIALS4PY_NVEC_SET_OP(buf_size, nvbufsize)
    SUNDIALS4PY_NVEC_SET_OP(buf_pack, nvbufpack)
    SUNDIALS4PY_NVEC_SET_OP(buf_unpack, nvbufunpack)
    SUNDIALS4PY_NVEC_SET_OP(print, nvprint)
    SUNDIALS4PY_NVEC_SET_OP(print_file, nvprintfile)
#undef SUNDIALS4PY_NVEC_SET_OP

    return v;
  }

  static std::shared_ptr<std::remove_pointer_t<N_Vector>> make_handle(
    nb::handle impl,
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
    Ownership ownership)
  {
    N_Vector v = create_raw_handle(impl, std::move(sunctx_owner), ownership);
    return sundials::experimental::our_make_shared<
      std::remove_pointer_t<N_Vector>, sundials::experimental::N_VectorDeleter>(v);
  }

  static Content* get_content(N_Vector v)
  {
    if (!v || !v->ops || v->ops->nvdestroy != custom_destroy || !v->content)
    {
      return nullptr;
    }
    return static_cast<Content*>(v->content);
  }

public:
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
  {
    return custom_content_impl(get_content(v), label);
  }

  static bool method_overridden(nb::handle impl, const char* name)
  {
    return custom_method_overridden<CustomNVector>(impl, name);
  }

  static void validate_required_methods(nb::handle impl)
  {
    if (!method_overridden(impl, "clone"))
    {
      throw nb::type_error("CustomNVector subclass must override clone()");
    }
  }

  static nb::object arg(N_Vector v)
  {
    return nb::cast(v, nb::rv_policy::reference);
  }

  template<typename Value, typename... Args>
  static Value call_value(N_Vector v, const char* name, const char* operation,
                          Value failure, Args&&... args)
  {
    try
    {
      nb::gil_scoped_acquire gil;
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

  static N_Vector_ID custom_get_vector_id(N_Vector)
  {
    return SUNDIALS_NVEC_CUSTOM;
  }

  static N_Vector clone_impl(N_Vector v, const char* name, const char* operation)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      nb::object cloned_impl = get_impl(v).attr(name)();
      if (!nb::isinstance<CustomNVector>(cloned_impl))
      {
        throw nb::type_error((std::string("CustomNVector.") + name +
                              "() must return a CustomNVector")
                               .c_str());
      }
      auto* cloned = nb::cast<CustomNVector*>(cloned_impl);
      if (!cloned->sunctx_owner_)
      {
        throw nb::type_error(
          "CustomNVector clone did not initialize the base constructor");
      }
      return create_raw_handle(cloned_impl, cloned->sunctx_owner_,
                               Ownership::strong);
    }
    SUNDIALS4PY_CATCH_AND_REPORT(v ? v->sunctx : nullptr, operation, nullptr)
  }

  static N_Vector custom_clone(N_Vector v)
  {
    return clone_impl(v, "clone", "CustomNVector.clone");
  }

  static N_Vector custom_clone_empty(N_Vector v)
  {
    return clone_impl(v, "clone_empty", "CustomNVector.clone_empty");
  }

  static void custom_destroy(N_Vector v)
  {
    if (!v) { return; }
    Content* content = get_content(v);
    custom_content_destroy(content);
    v->content = nullptr;
    N_VFreeEmpty(v);
  }

  static void custom_space(N_Vector v, sunindextype* lrw, sunindextype* liw)
  {
    *lrw = 0;
    *liw = 0;
    try
    {
      nb::gil_scoped_acquire gil;
      auto result = nb::cast<std::tuple<sunindextype, sunindextype>>(
        get_impl(v).attr("space")());
      *lrw = std::get<0>(result);
      *liw = std::get<1>(result);
    }
    catch (const std::exception& error)
    {
      report_custom_exception(v ? v->sunctx : nullptr, "CustomNVector.space",
                              error, __FILE__, __LINE__);
    }
    catch (...)
    {
      report_custom_unknown_exception(v ? v->sunctx : nullptr,
                                      "CustomNVector.space", __FILE__, __LINE__);
    }
  }

  static sunrealtype* custom_get_array_pointer(N_Vector v)
  {
    auto address = call_value<std::uintptr_t>(v, "get_array_pointer",
                                              "CustomNVector.get_array_pointer",
                                              0);
    return reinterpret_cast<sunrealtype*>(address);
  }

  static sunrealtype* custom_get_device_array_pointer(N_Vector v)
  {
    auto address =
      call_value<std::uintptr_t>(v, "get_device_array_pointer",
                                 "CustomNVector.get_device_array_pointer", 0);
    return reinterpret_cast<sunrealtype*>(address);
  }

  static void custom_set_array_pointer(sunrealtype* ptr, N_Vector v)
  {
    call_void(v, "set_array_pointer", "CustomNVector.set_array_pointer",
              reinterpret_cast<std::uintptr_t>(ptr));
  }

  static void custom_set_device_array_pointer(sunrealtype* ptr, N_Vector v)
  {
    call_void(v, "set_device_array_pointer",
              "CustomNVector.set_device_array_pointer",
              reinterpret_cast<std::uintptr_t>(ptr));
  }

  static SUNComm custom_get_communicator(N_Vector v)
  {
    return call_value<SUNComm>(v, "get_communicator",
                               "CustomNVector.get_communicator", SUN_COMM_NULL);
  }

  static sunindextype custom_get_length(N_Vector v)
  {
    return call_value<sunindextype>(v, "get_length", "CustomNVector.get_length",
                                    0);
  }

  static sunindextype custom_get_local_length(N_Vector v)
  {
    return call_value<sunindextype>(v, "get_local_length",
                                    "CustomNVector.get_local_length", 0);
  }

  static void custom_linear_sum(sunrealtype a, N_Vector x, sunrealtype b,
                                N_Vector y, N_Vector z)
  {
    call_void(z, "linear_sum", "CustomNVector.linear_sum", a, arg(x), b, arg(y));
  }

  static void custom_const(sunrealtype c, N_Vector z)
  {
    call_void(z, "const", "CustomNVector.const", c);
  }

  static void custom_prod(N_Vector x, N_Vector y, N_Vector z)
  {
    call_void(z, "prod", "CustomNVector.prod", arg(x), arg(y));
  }

  static void custom_div(N_Vector x, N_Vector y, N_Vector z)
  {
    call_void(z, "div", "CustomNVector.div", arg(x), arg(y));
  }

  static void custom_scale(sunrealtype c, N_Vector x, N_Vector z)
  {
    call_void(z, "scale", "CustomNVector.scale", c, arg(x));
  }

  static void custom_abs(N_Vector x, N_Vector z)
  {
    call_void(z, "abs", "CustomNVector.abs", arg(x));
  }

  static void custom_inv(N_Vector x, N_Vector z)
  {
    call_void(z, "inv", "CustomNVector.inv", arg(x));
  }

  static void custom_add_const(N_Vector x, sunrealtype b, N_Vector z)
  {
    call_void(z, "add_const", "CustomNVector.add_const", arg(x), b);
  }

  static sunrealtype custom_dot_prod(N_Vector x, N_Vector y)
  {
    return call_value<sunrealtype>(x, "dot_prod", "CustomNVector.dot_prod",
                                   SUN_RCONST(0.0), arg(y));
  }

  static sunrealtype custom_max_norm(N_Vector x)
  {
    return call_value<sunrealtype>(x, "max_norm", "CustomNVector.max_norm",
                                   SUN_RCONST(0.0));
  }

  static sunrealtype custom_wrms_norm(N_Vector x, N_Vector w)
  {
    return call_value<sunrealtype>(x, "wrms_norm", "CustomNVector.wrms_norm",
                                   SUN_RCONST(0.0), arg(w));
  }

  static sunrealtype custom_wrms_norm_mask(N_Vector x, N_Vector w, N_Vector id)
  {
    return call_value<sunrealtype>(x, "wrms_norm_mask",
                                   "CustomNVector.wrms_norm_mask",
                                   SUN_RCONST(0.0), arg(w), arg(id));
  }

  static sunrealtype custom_min(N_Vector x)
  {
    return call_value<sunrealtype>(x, "min", "CustomNVector.min",
                                   SUN_RCONST(0.0));
  }

  static sunrealtype custom_wl2_norm(N_Vector x, N_Vector w)
  {
    return call_value<sunrealtype>(x, "wl2_norm", "CustomNVector.wl2_norm",
                                   SUN_RCONST(0.0), arg(w));
  }

  static sunrealtype custom_l1_norm(N_Vector x)
  {
    return call_value<sunrealtype>(x, "l1_norm", "CustomNVector.l1_norm",
                                   SUN_RCONST(0.0));
  }

  static void custom_compare(sunrealtype c, N_Vector x, N_Vector z)
  {
    call_void(z, "compare", "CustomNVector.compare", c, arg(x));
  }

  static sunbooleantype custom_inv_test(N_Vector x, N_Vector z)
  {
    return call_value<sunbooleantype>(z, "inv_test", "CustomNVector.inv_test",
                                      SUNFALSE, arg(x));
  }

  static sunbooleantype custom_constr_mask(N_Vector c, N_Vector x, N_Vector m)
  {
    return call_value<sunbooleantype>(x, "constr_mask",
                                      "CustomNVector.constr_mask", SUNFALSE,
                                      arg(c), arg(m));
  }

  static sunrealtype custom_min_quotient(N_Vector num, N_Vector denom)
  {
    return call_value<sunrealtype>(num, "min_quotient",
                                   "CustomNVector.min_quotient",
                                   SUN_RCONST(0.0), arg(denom));
  }

  static nb::list vector_list(N_Vector* vectors, int count)
  {
    nb::list result;
    for (int i = 0; i < count; ++i) { result.append(arg(vectors[i])); }
    return result;
  }

  static nb::list vector_lists(N_Vector** vectors, int outer_count,
                               int inner_count)
  {
    nb::list result;
    for (int i = 0; i < outer_count; ++i)
    {
      result.append(vector_list(vectors[i], inner_count));
    }
    return result;
  }

  static std::vector<sunrealtype> real_values(const sunrealtype* values, int count)
  {
    return std::vector<sunrealtype>(values, values + count);
  }

  static SUNErrCode custom_linear_combination(int nvec, sunrealtype* c,
                                              N_Vector* X, N_Vector z)
  {
    return call_status(z, "linear_combination",
                       "CustomNVector.linear_combination", real_values(c, nvec),
                       vector_list(X, nvec));
  }

  static SUNErrCode custom_scale_add_multi(int nvec, sunrealtype* a, N_Vector x,
                                           N_Vector* Y, N_Vector* Z)
  {
    return call_status(x, "scale_add_multi", "CustomNVector.scale_add_multi",
                       real_values(a, nvec), vector_list(Y, nvec),
                       vector_list(Z, nvec));
  }

  static SUNErrCode custom_dot_prod_multi(int nvec, N_Vector x, N_Vector* Y,
                                          sunrealtype* dots)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      auto result = nb::cast<std::tuple<int, std::vector<sunrealtype>>>(
        get_impl(x).attr("dot_prod_multi")(vector_list(Y, nvec)));
      const auto& values = std::get<1>(result);
      if (static_cast<int>(values.size()) != nvec)
      {
        throw nb::value_error(
          "dot_prod_multi() returned the wrong result size");
      }
      std::copy(values.begin(), values.end(), dots);
      return static_cast<SUNErrCode>(std::get<0>(result));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(x ? x->sunctx : nullptr,
                                 "CustomNVector.dot_prod_multi", SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_linear_sum_vector_array(int nvec, sunrealtype a,
                                                   N_Vector* X, sunrealtype b,
                                                   N_Vector* Y, N_Vector* Z)
  {
    return call_status(Z[0], "linear_sum_vector_array",
                       "CustomNVector.linear_sum_vector_array", a,
                       vector_list(X, nvec), b, vector_list(Y, nvec),
                       vector_list(Z, nvec));
  }

  static SUNErrCode custom_scale_vector_array(int nvec, sunrealtype* c,
                                              N_Vector* X, N_Vector* Z)
  {
    return call_status(Z[0], "scale_vector_array",
                       "CustomNVector.scale_vector_array", real_values(c, nvec),
                       vector_list(X, nvec), vector_list(Z, nvec));
  }

  static SUNErrCode custom_const_vector_array(int nvec, sunrealtype c, N_Vector* Z)
  {
    return call_status(Z[0], "const_vector_array",
                       "CustomNVector.const_vector_array", c,
                       vector_list(Z, nvec));
  }

  template<typename Invoke>
  static SUNErrCode vector_array_reduction(N_Vector owner, int nvec,
                                           sunrealtype* values,
                                           const char* operation, Invoke&& invoke)
  {
    try
    {
      nb::gil_scoped_acquire gil;
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
    SUNDIALS4PY_CATCH_AND_REPORT(owner ? owner->sunctx : nullptr, operation,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_wrms_norm_vector_array(int nvec, N_Vector* X,
                                                  N_Vector* W, sunrealtype* norms)
  {
    return vector_array_reduction(X[0], nvec, norms,
                                  "CustomNVector.wrms_norm_vector_array",
                                  [&](nb::object impl)
                                  {
                                    return impl.attr(
                                      "wrms_norm_vector_array")(vector_list(X,
                                                                            nvec),
                                                                vector_list(W,
                                                                            nvec));
                                  });
  }

  static SUNErrCode custom_wrms_norm_mask_vector_array(int nvec, N_Vector* X,
                                                       N_Vector* W, N_Vector id,
                                                       sunrealtype* norms)
  {
    return vector_array_reduction(X[0], nvec, norms,
                                  "CustomNVector.wrms_norm_mask_vector_array",
                                  [&](nb::object impl)
                                  {
                                    return impl.attr(
                                      "wrms_norm_mask_vector_array")(vector_list(X,
                                                                                 nvec),
                                                                     vector_list(W,
                                                                                 nvec),
                                                                     arg(id));
                                  });
  }

  static SUNErrCode custom_scale_add_multi_vector_array(
    int nvec, int nsum, sunrealtype* a, N_Vector* X, N_Vector** Y, N_Vector** Z)
  {
    return call_status(X[0], "scale_add_multi_vector_array",
                       "CustomNVector.scale_add_multi_vector_array",
                       real_values(a, nsum), vector_list(X, nvec),
                       vector_lists(Y, nsum, nvec), vector_lists(Z, nsum, nvec));
  }

  static SUNErrCode custom_linear_combination_vector_array(int nvec, int nsum,
                                                           sunrealtype* c,
                                                           N_Vector** X,
                                                           N_Vector* Z)
  {
    return call_status(Z[0], "linear_combination_vector_array",
                       "CustomNVector.linear_combination_vector_array",
                       real_values(c, nsum), vector_lists(X, nsum, nvec),
                       vector_list(Z, nvec));
  }

  static sunrealtype custom_dot_prod_local(N_Vector x, N_Vector y)
  {
    return call_value<sunrealtype>(x, "dot_prod_local",
                                   "CustomNVector.dot_prod_local",
                                   SUN_RCONST(0.0), arg(y));
  }

  static sunrealtype custom_max_norm_local(N_Vector x)
  {
    return call_value<sunrealtype>(x, "max_norm_local",
                                   "CustomNVector.max_norm_local",
                                   SUN_RCONST(0.0));
  }

  static sunrealtype custom_min_local(N_Vector x)
  {
    return call_value<sunrealtype>(x, "min_local", "CustomNVector.min_local",
                                   SUN_RCONST(0.0));
  }

  static sunrealtype custom_l1_norm_local(N_Vector x)
  {
    return call_value<sunrealtype>(x, "l1_norm_local",
                                   "CustomNVector.l1_norm_local",
                                   SUN_RCONST(0.0));
  }

  static sunbooleantype custom_inv_test_local(N_Vector x, N_Vector z)
  {
    return call_value<sunbooleantype>(z, "inv_test_local",
                                      "CustomNVector.inv_test_local", SUNFALSE,
                                      arg(x));
  }

  static sunbooleantype custom_constr_mask_local(N_Vector c, N_Vector x,
                                                 N_Vector m)
  {
    return call_value<sunbooleantype>(x, "constr_mask_local",
                                      "CustomNVector.constr_mask_local",
                                      SUNFALSE, arg(c), arg(m));
  }

  static sunrealtype custom_min_quotient_local(N_Vector num, N_Vector denom)
  {
    return call_value<sunrealtype>(num, "min_quotient_local",
                                   "CustomNVector.min_quotient_local",
                                   SUN_RCONST(0.0), arg(denom));
  }

  static sunrealtype custom_wsqrsum_local(N_Vector x, N_Vector w)
  {
    return call_value<sunrealtype>(x, "wsqrsum_local",
                                   "CustomNVector.wsqrsum_local",
                                   SUN_RCONST(0.0), arg(w));
  }

  static sunrealtype custom_wsqrsum_mask_local(N_Vector x, N_Vector w, N_Vector id)
  {
    return call_value<sunrealtype>(x, "wsqrsum_mask_local",
                                   "CustomNVector.wsqrsum_mask_local",
                                   SUN_RCONST(0.0), arg(w), arg(id));
  }

  static SUNErrCode custom_dot_prod_multi_local(int nvec, N_Vector x,
                                                N_Vector* Y, sunrealtype* dots)
  {
    return vector_array_reduction(x, nvec, dots,
                                  "CustomNVector.dot_prod_multi_local",
                                  [&](nb::object impl) {
                                    return impl.attr("dot_prod_multi_local")(
                                      vector_list(Y, nvec));
                                  });
  }

  static SUNErrCode custom_dot_prod_multi_all_reduce(int nvec, N_Vector x,
                                                     sunrealtype* dots)
  {
    return vector_array_reduction(x, nvec, dots,
                                  "CustomNVector.dot_prod_multi_all_reduce",
                                  [&](nb::object impl) {
                                    return impl.attr(
                                      "dot_prod_multi_all_reduce")(
                                      real_values(dots, nvec));
                                  });
  }

  static SUNErrCode custom_buf_size(N_Vector x, sunindextype* size)
  {
    *size = 0;
    try
    {
      nb::gil_scoped_acquire gil;
      auto result =
        nb::cast<std::tuple<int, sunindextype>>(get_impl(x).attr("buf_size")());
      *size = std::get<1>(result);
      return static_cast<SUNErrCode>(std::get<0>(result));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(x ? x->sunctx : nullptr,
                                 "CustomNVector.buf_size", SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_buf_pack(N_Vector x, void* buffer)
  {
    return call_status(x, "buf_pack", "CustomNVector.buf_pack",
                       reinterpret_cast<std::uintptr_t>(buffer));
  }

  static SUNErrCode custom_buf_unpack(N_Vector x, void* buffer)
  {
    return call_status(x, "buf_unpack", "CustomNVector.buf_unpack",
                       reinterpret_cast<std::uintptr_t>(buffer));
  }

  static void custom_print(N_Vector x)
  {
    call_void(x, "print", "CustomNVector.print");
  }

  static void custom_print_file(N_Vector x, FILE* outfile)
  {
    call_void(x, "print_file", "CustomNVector.print_file",
              reinterpret_cast<std::uintptr_t>(outfile));
  }

  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner_;
  std::shared_ptr<std::remove_pointer_t<N_Vector>> vector_;
  HandleState state_{HandleState::unmaterialized};
};

} // namespace sundials4py

#endif // SUNDIALS4PY_NVECTOR_CUSTOM_HPP
