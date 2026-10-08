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

std::shared_ptr<std::remove_pointer_t<N_Vector>> CustomNVector::make_handle(
  nb::handle impl)
{ return make_handle(impl, sunctx_owner_, Ownership::weak); }

N_Vector CustomNVector::create_raw_handle(
  nb::handle impl,
  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
  Ownership ownership)
{
  validate_required_methods(impl);

  const bool has_is_compatible = method_overridden(impl, "is_compatible");
  const bool has_clone_empty   = method_overridden(impl, "clone_empty");
  const bool has_space         = method_overridden(impl, "space");
  const bool has_get_array     = method_overridden(impl, "get_array_pointer");
  const bool has_get_device    = method_overridden(impl,
                                                   "get_device_array_pointer");
  const bool has_set_array     = method_overridden(impl, "set_array_pointer");
  const bool has_set_device    = method_overridden(impl,
                                                   "set_device_array_pointer");
  const bool has_get_comm      = method_overridden(impl, "get_communicator");
  const bool has_get_local_length = method_overridden(impl, "get_local_length");
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
  CustomNVector* owner = nullptr;
  if (ownership == Ownership::strong)
  {
    owner          = nb::cast<CustomNVector*>(impl);
    content->owner = owner;
  }
  content->sunctx_owner      = std::move(sunctx_owner);
  content->has_is_compatible = has_is_compatible;
  if (ownership == Ownership::weak) { content->weak_impl = nb::weakref(impl); }
  else
  {
    content->strong_impl = nb::borrow<nb::object>(impl);
  }

  N_Vector v = N_VNewEmpty(content->sunctx_owner.get());
  if (!v) { throw error_returned("N_VNewEmpty failed"); }

  v->content = content.release();
  if (owner) { owner->raw_handle_ = v; }
  v->ops->nvgetvectorid = custom_get_vector_id;
  v->ops->nvclone       = custom_clone;
  v->ops->nvcloneempty  = has_clone_empty ? custom_clone_empty : custom_clone;
  v->ops->nvdestroy     = custom_destroy;
  if (has_space) { v->ops->nvspace = custom_space; }
  v->ops->nvgetlength      = custom_get_length;
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
  SUNDIALS4PY_NVEC_SET_OP(scale_add_multi_vector_array, nvscaleaddmultivectorarray)
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

N_Vector CustomNVector::clone_impl(N_Vector v, const char* name,
                                   const char* operation)
{
  try
  {
    nb::gil_scoped_acquire gil;
    nb::object impl        = get_impl(v);
    nb::object cloned_impl = impl.attr(name)();
    if (!nb::isinstance<CustomNVector>(cloned_impl))
    {
      throw nb::type_error(
        (std::string("CustomNVector.") + name + "() must return a CustomNVector")
          .c_str());
    }
    auto* cloned = nb::cast<CustomNVector*>(cloned_impl);
    if (!cloned->sunctx_owner_)
    {
      throw nb::type_error(
        "CustomNVector clone did not initialize the base constructor");
    }
    require_same_context(v->sunctx, cloned->sunctx_owner_, operation);
    if (Py_TYPE(cloned_impl.ptr()) == Py_TYPE(impl.ptr()))
    {
      return clone_raw_handle(cloned_impl, v, cloned->sunctx_owner_);
    }
    return create_raw_handle(cloned_impl, cloned->sunctx_owner_,
                             Ownership::strong);
  }
  SUNDIALS4PY_CATCH_AND_REPORT(v ? v->sunctx : nullptr, operation, nullptr)
}

N_Vector CustomNVector::clone_raw_handle(
  nb::handle impl, N_Vector source,
  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner)
{
  Content* source_content = get_content(source);
  if (!source_content)
  {
    throw nb::type_error("source is not a custom N_Vector");
  }

  auto* owner = nb::cast<CustomNVector*>(impl);
  std::unique_ptr<Content> content(new Content);
  content->owner             = owner;
  content->sunctx_owner      = std::move(sunctx_owner);
  content->has_is_compatible = source_content->has_is_compatible;
  content->strong_impl       = nb::borrow<nb::object>(impl);

  N_Vector v = N_VNewEmpty(content->sunctx_owner.get());
  if (!v) { throw error_returned("N_VNewEmpty failed"); }
  if (N_VCopyOps(source, v) != SUN_SUCCESS)
  {
    N_VFreeEmpty(v);
    throw error_returned("N_VCopyOps failed");
  }

  v->content         = content.release();
  owner->raw_handle_ = v;
  return v;
}

std::shared_ptr<std::remove_pointer_t<N_Vector>> CustomNVector::make_handle(
  nb::handle impl,
  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
  Ownership ownership)
{
  N_Vector v = create_raw_handle(impl, std::move(sunctx_owner), ownership);
  return sundials::experimental::our_make_shared<
    std::remove_pointer_t<N_Vector>, sundials::experimental::N_VectorDeleter>(v);
}

N_Vector_ID CustomNVector::custom_get_vector_id(N_Vector)
{ return SUNDIALS_NVEC_CUSTOM; }

N_Vector CustomNVector::custom_clone(N_Vector v)
{ return clone_impl(v, "clone", "CustomNVector.clone"); }

N_Vector CustomNVector::custom_clone_empty(N_Vector v)
{ return clone_impl(v, "clone_empty", "CustomNVector.clone_empty"); }

void CustomNVector::custom_destroy(N_Vector v)
{
  if (!v) { return; }
  Content* content = get_content(v);
  if (content && content->owner && content->owner->raw_handle_ == v)
  {
    content->owner->raw_handle_ = nullptr;
  }
  custom_content_destroy(content);
  v->content = nullptr;
  N_VFreeEmpty(v);
}

void CustomNVector::custom_space(N_Vector v, sunindextype* lrw, sunindextype* liw)
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

sunrealtype* CustomNVector::custom_get_array_pointer(N_Vector v)
{
  auto address = call_value<std::uintptr_t>(v, "get_array_pointer",
                                            "CustomNVector.get_array_pointer", 0);
  return reinterpret_cast<sunrealtype*>(address);
}

sunrealtype* CustomNVector::custom_get_device_array_pointer(N_Vector v)
{
  auto address =
    call_value<std::uintptr_t>(v, "get_device_array_pointer",
                               "CustomNVector.get_device_array_pointer", 0);
  return reinterpret_cast<sunrealtype*>(address);
}

void CustomNVector::custom_set_array_pointer(sunrealtype* ptr, N_Vector v)
{
  call_void(v, "set_array_pointer", "CustomNVector.set_array_pointer",
            reinterpret_cast<std::uintptr_t>(ptr));
}

void CustomNVector::custom_set_device_array_pointer(sunrealtype* ptr, N_Vector v)
{
  call_void(v, "set_device_array_pointer",
            "CustomNVector.set_device_array_pointer",
            reinterpret_cast<std::uintptr_t>(ptr));
}

SUNComm CustomNVector::custom_get_communicator(N_Vector v)
{
  return call_value<SUNComm>(v, "get_communicator",
                             "CustomNVector.get_communicator", SUN_COMM_NULL);
}

sunindextype CustomNVector::custom_get_length(N_Vector v)
{
  return call_value<sunindextype>(v, "get_length", "CustomNVector.get_length", 0);
}

sunindextype CustomNVector::custom_get_local_length(N_Vector v)
{
  return call_value<sunindextype>(v, "get_local_length",
                                  "CustomNVector.get_local_length", 0);
}

void CustomNVector::custom_linear_sum(sunrealtype a, N_Vector x, sunrealtype b,
                                      N_Vector y, N_Vector z)
{
  call_void_with(z, "CustomNVector.linear_sum",
                 [&](nb::handle self)
                 {
                   self.attr("linear_sum")(a, operand(self, z, x, "linear_sum"),
                                           b, operand(self, z, y, "linear_sum"));
                 });
}

void CustomNVector::custom_const(sunrealtype c, N_Vector z)
{ call_void(z, "const", "CustomNVector.const", c); }

void CustomNVector::custom_prod(N_Vector x, N_Vector y, N_Vector z)
{
  call_void_with(z, "CustomNVector.prod",
                 [&](nb::handle self)
                 {
                   self.attr("prod")(operand(self, z, x, "prod"),
                                     operand(self, z, y, "prod"));
                 });
}

void CustomNVector::custom_div(N_Vector x, N_Vector y, N_Vector z)
{
  call_void_with(z, "CustomNVector.div",
                 [&](nb::handle self)
                 {
                   self.attr("div")(operand(self, z, x, "div"),
                                    operand(self, z, y, "div"));
                 });
}

void CustomNVector::custom_scale(sunrealtype c, N_Vector x, N_Vector z)
{
  call_void_with(z, "CustomNVector.scale", [&](nb::handle self)
                 { self.attr("scale")(c, operand(self, z, x, "scale")); });
}

void CustomNVector::custom_abs(N_Vector x, N_Vector z)
{
  call_void_with(z, "CustomNVector.abs", [&](nb::handle self)
                 { self.attr("abs")(operand(self, z, x, "abs")); });
}

void CustomNVector::custom_inv(N_Vector x, N_Vector z)
{
  call_void_with(z, "CustomNVector.inv", [&](nb::handle self)
                 { self.attr("inv")(operand(self, z, x, "inv")); });
}

void CustomNVector::custom_add_const(N_Vector x, sunrealtype b, N_Vector z)
{
  call_void_with(z, "CustomNVector.add_const", [&](nb::handle self)
                 { self.attr("add_const")(operand(self, z, x, "add_const"), b); });
}

sunrealtype CustomNVector::custom_dot_prod(N_Vector x, N_Vector y)
{
  return call_value_with<sunrealtype>(x, "CustomNVector.dot_prod",
                                      std::numeric_limits<sunrealtype>::quiet_NaN(),
                                      [&](nb::handle self)
                                      {
                                        return self.attr("dot_prod")(
                                          operand(self, x, y, "dot_prod"));
                                      });
}

sunrealtype CustomNVector::custom_max_norm(N_Vector x)
{
  return call_value<sunrealtype>(x, "max_norm", "CustomNVector.max_norm",
                                 std::numeric_limits<sunrealtype>::infinity());
}

sunrealtype CustomNVector::custom_wrms_norm(N_Vector x, N_Vector w)
{
  return call_value_with<sunrealtype>(x, "CustomNVector.wrms_norm",
                                      std::numeric_limits<sunrealtype>::infinity(),
                                      [&](nb::handle self)
                                      {
                                        return self.attr("wrms_norm")(
                                          operand(self, x, w, "wrms_norm"));
                                      });
}

sunrealtype CustomNVector::custom_wrms_norm_mask(N_Vector x, N_Vector w,
                                                 N_Vector id)
{
  return call_value_with<
    sunrealtype>(x, "CustomNVector.wrms_norm_mask",
                 std::numeric_limits<sunrealtype>::infinity(),
                 [&](nb::handle self)
                 {
                   return self.attr(
                     "wrms_norm_mask")(operand(self, x, w, "wrms_norm_mask"),
                                       operand(self, x, id, "wrms_norm_mask"));
                 });
}

sunrealtype CustomNVector::custom_min(N_Vector x)
{
  return call_value<sunrealtype>(x, "min", "CustomNVector.min",
                                 -std::numeric_limits<sunrealtype>::infinity());
}

sunrealtype CustomNVector::custom_wl2_norm(N_Vector x, N_Vector w)
{
  return call_value_with<sunrealtype>(x, "CustomNVector.wl2_norm",
                                      std::numeric_limits<sunrealtype>::infinity(),
                                      [&](nb::handle self)
                                      {
                                        return self.attr("wl2_norm")(
                                          operand(self, x, w, "wl2_norm"));
                                      });
}

sunrealtype CustomNVector::custom_l1_norm(N_Vector x)
{
  return call_value<sunrealtype>(x, "l1_norm", "CustomNVector.l1_norm",
                                 std::numeric_limits<sunrealtype>::infinity());
}

void CustomNVector::custom_compare(sunrealtype c, N_Vector x, N_Vector z)
{
  call_void_with(z, "CustomNVector.compare", [&](nb::handle self)
                 { self.attr("compare")(c, operand(self, z, x, "compare")); });
}

sunbooleantype CustomNVector::custom_inv_test(N_Vector x, N_Vector z)
{
  return call_value_with<sunbooleantype>(z, "CustomNVector.inv_test", SUNFALSE,
                                         [&](nb::handle self)
                                         {
                                           return self.attr("inv_test")(
                                             operand(self, z, x, "inv_test"));
                                         });
}

sunbooleantype CustomNVector::custom_constr_mask(N_Vector c, N_Vector x,
                                                 N_Vector m)
{
  return call_value_with<
    sunbooleantype>(x, "CustomNVector.constr_mask", SUNFALSE,
                    [&](nb::handle self)
                    {
                      return self.attr(
                        "constr_mask")(operand(self, x, c, "constr_mask"),
                                       operand(self, x, m, "constr_mask"));
                    });
}

sunrealtype CustomNVector::custom_min_quotient(N_Vector num, N_Vector denom)
{
  return call_value_with<sunrealtype>(num, "CustomNVector.min_quotient",
                                      std::numeric_limits<sunrealtype>::quiet_NaN(),
                                      [&](nb::handle self)
                                      {
                                        return self.attr("min_quotient")(
                                          operand(self, num, denom,
                                                  "min_quotient"));
                                      });
}

SUNErrCode CustomNVector::custom_linear_combination(int nvec, sunrealtype* c,
                                                    N_Vector* X, N_Vector z)
{
  return call_status_with(z, "CustomNVector.linear_combination",
                          [&](nb::handle self)
                          {
                            return self.attr(
                              "linear_combination")(real_values(c, nvec),
                                                    vector_list(self, z, X,
                                                                nvec, "linear_combination"));
                          });
}

SUNErrCode CustomNVector::custom_scale_add_multi(int nvec, sunrealtype* a,
                                                 N_Vector x, N_Vector* Y,
                                                 N_Vector* Z)
{
  return call_status_with(x, "CustomNVector.scale_add_multi",
                          [&](nb::handle self)
                          {
                            return self.attr(
                              "scale_add_multi")(real_values(a, nvec),
                                                 vector_list(self, x, Y, nvec,
                                                             "scale_add_multi"),
                                                 vector_list(self, x, Z,
                                                             nvec, "scale_add_multi"));
                          });
}

SUNErrCode CustomNVector::custom_dot_prod_multi(int nvec, N_Vector x,
                                                N_Vector* Y, sunrealtype* dots)
{
  return vector_array_reduction(x, nvec, dots,
                                std::numeric_limits<sunrealtype>::quiet_NaN(),
                                "CustomNVector.dot_prod_multi",
                                [&](nb::object impl)
                                {
                                  return impl.attr("dot_prod_multi")(
                                    vector_list(impl, x, Y, nvec,
                                                "dot_prod_multi"));
                                });
}

SUNErrCode CustomNVector::custom_linear_sum_vector_array(
  int nvec, sunrealtype a, N_Vector* X, sunrealtype b, N_Vector* Y, N_Vector* Z)
{
  return call_status_with(Z[0], "CustomNVector.linear_sum_vector_array",
                          [&](nb::handle self)
                          {
                            return self.attr(
                              "linear_sum_vector_array")(a,
                                                         vector_list(self, Z[0], X,
                                                                     nvec, "linear_sum_vector_array"),
                                                         b,
                                                         vector_list(self, Z[0], Y,
                                                                     nvec, "linear_sum_vector_array"),
                                                         vector_list(self, Z[0], Z,
                                                                     nvec, "linear_sum_vector_array"));
                          });
}

SUNErrCode CustomNVector::custom_scale_vector_array(int nvec, sunrealtype* c,
                                                    N_Vector* X, N_Vector* Z)
{
  return call_status_with(Z[0], "CustomNVector.scale_vector_array",
                          [&](nb::handle self)
                          {
                            return self.attr(
                              "scale_vector_array")(real_values(c, nvec),
                                                    vector_list(self, Z[0], X,
                                                                nvec, "scale_vector_array"),
                                                    vector_list(self, Z[0], Z,
                                                                nvec, "scale_vector_array"));
                          });
}

SUNErrCode CustomNVector::custom_const_vector_array(int nvec, sunrealtype c,
                                                    N_Vector* Z)
{
  return call_status_with(Z[0], "CustomNVector.const_vector_array",
                          [&](nb::handle self)
                          {
                            return self.attr(
                              "const_vector_array")(c,
                                                    vector_list(self, Z[0], Z,
                                                                nvec, "const_vector_array"));
                          });
}

SUNErrCode CustomNVector::custom_wrms_norm_vector_array(int nvec, N_Vector* X,
                                                        N_Vector* W,
                                                        sunrealtype* norms)
{
  return vector_array_reduction(X[0], nvec, norms,
                                std::numeric_limits<sunrealtype>::infinity(),
                                "CustomNVector.wrms_norm_vector_array",
                                [&](nb::object impl)
                                {
                                  return impl.attr(
                                    "wrms_norm_vector_array")(vector_list(impl,
                                                                          X[0], X,
                                                                          nvec, "wrms_norm_vector_array"),
                                                              vector_list(impl,
                                                                          X[0], W,
                                                                          nvec, "wrms_norm_vector_array"));
                                });
}

SUNErrCode CustomNVector::custom_wrms_norm_mask_vector_array(
  int nvec, N_Vector* X, N_Vector* W, N_Vector id, sunrealtype* norms)
{
  return vector_array_reduction(X[0], nvec, norms,
                                std::numeric_limits<sunrealtype>::infinity(),
                                "CustomNVector.wrms_norm_mask_vector_array",
                                [&](nb::object impl)
                                {
                                  return impl.attr(
                                    "wrms_norm_mask_vector_"
                                    "array")(vector_list(impl, X[0], X,
                                                         nvec, "wrms_norm_mask_vector_array"),
                                             vector_list(impl, X[0], W,
                                                         nvec, "wrms_norm_mask_vector_array"),
                                             operand(impl, X[0],
                                                     id, "wrms_norm_mask_vector_array"));
                                });
}

SUNErrCode CustomNVector::custom_scale_add_multi_vector_array(
  int nvec, int nsum, sunrealtype* a, N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  return call_status_with(X[0], "CustomNVector.scale_add_multi_vector_array",
                          [&](nb::handle self)
                          {
                            return self.attr(
                              "scale_add_multi_vector_array")(real_values(a, nsum),
                                                              vector_list(self,
                                                                          X[0], X,
                                                                          nvec, "scale_add_multi_vector_array"),
                                                              vector_lists(self,
                                                                           X[0],
                                                                           Y, nsum,
                                                                           nvec, "scale_add_multi_vector_array"),
                                                              vector_lists(self,
                                                                           X[0],
                                                                           Z, nsum,
                                                                           nvec, "scale_add_multi_vector_array"));
                          });
}

SUNErrCode CustomNVector::custom_linear_combination_vector_array(
  int nvec, int nsum, sunrealtype* c, N_Vector** X, N_Vector* Z)
{
  return call_status_with(Z[0], "CustomNVector.linear_combination_vector_array",
                          [&](nb::handle self)
                          {
                            return self.attr(
                              "linear_combination_vector_"
                              "array")(real_values(c, nsum),
                                       vector_lists(self, Z[0], X, nsum,
                                                    nvec, "linear_combination_vector_array"),
                                       vector_list(self, Z[0], Z,
                                                   nvec, "linear_combination_vector_array"));
                          });
}

sunrealtype CustomNVector::custom_dot_prod_local(N_Vector x, N_Vector y)
{
  return call_value_with<sunrealtype>(x, "CustomNVector.dot_prod_local",
                                      std::numeric_limits<sunrealtype>::quiet_NaN(),
                                      [&](nb::handle self)
                                      {
                                        return self.attr("dot_prod_local")(
                                          operand(self, x, y, "dot_prod_local"));
                                      });
}

sunrealtype CustomNVector::custom_max_norm_local(N_Vector x)
{
  return call_value<sunrealtype>(x, "max_norm_local",
                                 "CustomNVector.max_norm_local",
                                 std::numeric_limits<sunrealtype>::infinity());
}

sunrealtype CustomNVector::custom_min_local(N_Vector x)
{
  return call_value<sunrealtype>(x, "min_local", "CustomNVector.min_local",
                                 -std::numeric_limits<sunrealtype>::infinity());
}

sunrealtype CustomNVector::custom_l1_norm_local(N_Vector x)
{
  return call_value<sunrealtype>(x, "l1_norm_local",
                                 "CustomNVector.l1_norm_local",
                                 std::numeric_limits<sunrealtype>::infinity());
}

sunbooleantype CustomNVector::custom_inv_test_local(N_Vector x, N_Vector z)
{
  return call_value_with<sunbooleantype>(z, "CustomNVector.inv_test_local",
                                         SUNFALSE,
                                         [&](nb::handle self)
                                         {
                                           return self.attr("inv_test_local")(
                                             operand(self, z, x,
                                                     "inv_test_local"));
                                         });
}

sunbooleantype CustomNVector::custom_constr_mask_local(N_Vector c, N_Vector x,
                                                       N_Vector m)
{
  return call_value_with<
    sunbooleantype>(x, "CustomNVector.constr_mask_local", SUNFALSE,
                    [&](nb::handle self)
                    {
                      return self.attr(
                        "constr_mask_local")(operand(self, x, c,
                                                     "constr_mask_local"),
                                             operand(self, x, m,
                                                     "constr_mask_local"));
                    });
}

sunrealtype CustomNVector::custom_min_quotient_local(N_Vector num, N_Vector denom)
{
  return call_value_with<sunrealtype>(num, "CustomNVector.min_quotient_local",
                                      std::numeric_limits<sunrealtype>::quiet_NaN(),
                                      [&](nb::handle self)
                                      {
                                        return self.attr("min_quotient_local")(
                                          operand(self, num, denom,
                                                  "min_quotient_local"));
                                      });
}

sunrealtype CustomNVector::custom_wsqrsum_local(N_Vector x, N_Vector w)
{
  return call_value_with<sunrealtype>(x, "CustomNVector.wsqrsum_local",
                                      std::numeric_limits<sunrealtype>::infinity(),
                                      [&](nb::handle self)
                                      {
                                        return self.attr("wsqrsum_local")(
                                          operand(self, x, w, "wsqrsum_local"));
                                      });
}

sunrealtype CustomNVector::custom_wsqrsum_mask_local(N_Vector x, N_Vector w,
                                                     N_Vector id)
{
  return call_value_with<
    sunrealtype>(x, "CustomNVector.wsqrsum_mask_local",
                 std::numeric_limits<sunrealtype>::infinity(),
                 [&](nb::handle self)
                 {
                   return self.attr(
                     "wsqrsum_mask_local")(operand(self, x, w,
                                                   "wsqrsum_mask_local"),
                                           operand(self, x, id,
                                                   "wsqrsum_mask_local"));
                 });
}

SUNErrCode CustomNVector::custom_dot_prod_multi_local(int nvec, N_Vector x,
                                                      N_Vector* Y,
                                                      sunrealtype* dots)
{
  return vector_array_reduction(x, nvec, dots,
                                std::numeric_limits<sunrealtype>::quiet_NaN(),
                                "CustomNVector.dot_prod_multi_local",
                                [&](nb::object impl)
                                {
                                  return impl.attr("dot_prod_multi_local")(
                                    vector_list(impl, x, Y, nvec,
                                                "dot_prod_multi_local"));
                                });
}

SUNErrCode CustomNVector::custom_dot_prod_multi_all_reduce(int nvec, N_Vector x,
                                                           sunrealtype* dots)
{
  return vector_array_reduction(x, nvec, dots,
                                std::numeric_limits<sunrealtype>::quiet_NaN(),
                                "CustomNVector.dot_prod_multi_all_reduce",
                                [&](nb::object impl)
                                {
                                  return impl.attr("dot_prod_multi_all_reduce")(
                                    real_values(dots, nvec));
                                });
}

SUNErrCode CustomNVector::custom_buf_size(N_Vector x, sunindextype* size)
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

SUNErrCode CustomNVector::custom_buf_pack(N_Vector x, void* buffer)
{
  return call_status(x, "buf_pack", "CustomNVector.buf_pack",
                     reinterpret_cast<std::uintptr_t>(buffer));
}

SUNErrCode CustomNVector::custom_buf_unpack(N_Vector x, void* buffer)
{
  return call_status(x, "buf_unpack", "CustomNVector.buf_unpack",
                     reinterpret_cast<std::uintptr_t>(buffer));
}

void CustomNVector::custom_print(N_Vector x)
{ call_void(x, "print", "CustomNVector.print"); }

void CustomNVector::custom_print_file(N_Vector x, FILE* outfile)
{
  call_void(x, "print_file", "CustomNVector.print_file",
            reinterpret_cast<std::uintptr_t>(outfile));
}

} // namespace sundials4py
