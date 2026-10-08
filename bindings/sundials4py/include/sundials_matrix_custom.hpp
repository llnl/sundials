/*------------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *------------------------------------------------------------------------------
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
 *------------------------------------------------------------------------------
 * CustomSUNMatrix: the base class Python code subclasses to implement a
 * SUNMatrix.
 *----------------------------------------------------------------------------*/

#ifndef SUNDIALS4PY_MATRIX_CUSTOM_HPP
#define SUNDIALS4PY_MATRIX_CUSTOM_HPP

#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include <sundials/sundials_classview.hpp>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_matrix.hpp>

#include "sundials4py_core_types.hpp"
#include "sundials4py_custom_object.hpp"

namespace sundials4py {

/*
 * Python-owned SUNMatrix implementation.
 *
 * A Python subclass remains an ordinary nanobind object until a generated or
 * hand-written wrapper asks for SUNMatrix. At that point the custom type caster
 * calls _get_sundials_handle(), which lazily builds a SUNMatrix shell, fills in
 * the SUNDIALS operation table, and stores a reference back to the Python
 * implementation.
 *
 * Matrices are the one family where SUNDIALS creates handles of its own accord:
 * SUNMatClone() must hand back a matrix that native code owns outright. Those
 * clones therefore hold a strong reference to their Python implementation, while
 * handles materialized for a user-held object hold only a weak one.
 */
class CustomSUNMatrix
  : public CustomObjectBase<CustomSUNMatrix, std::remove_pointer_t<SUNMatrix>>
{
public:
  explicit CustomSUNMatrix(std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx)
    : CustomObjectBase(std::move(sunctx), "CustomSUNMatrix")
  {}

  void validate(nb::handle self) const { validate_required_methods(self); }

  static SUNErrCode base_method_status(const char* name)
  {
    PyErr_SetString(PyExc_NotImplementedError, name);
    nb::raise_python_error();
    return SUN_ERR_EXT_FAIL;
  }

private:
  enum class Ownership
  {
    weak,
    strong
  };

  struct Content : CustomContentBase
  {
    bool has_is_compatible{false};
  };

  static constexpr const char* label = "CustomSUNMatrix";

  friend class CustomObjectBase<CustomSUNMatrix, std::remove_pointer_t<SUNMatrix>>;

  std::shared_ptr<std::remove_pointer_t<SUNMatrix>> make_handle(nb::handle impl);

  static SUNMatrix create_raw_handle(
    nb::handle impl,
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
    Ownership ownership);

  static std::shared_ptr<std::remove_pointer_t<SUNMatrix>> make_handle(
    nb::handle impl,
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
    Ownership ownership);

  static Content* get_content(SUNMatrix A)
  {
    if (!A || !A->ops || A->ops->destroy != custom_matrix_destroy || !A->content)
    {
      return nullptr;
    }
    return static_cast<Content*>(A->content);
  }

public:
  static nb::object clone_python(nb::handle impl, const char* operation)
  {
    auto* source           = nb::cast<CustomSUNMatrix*>(impl);
    nb::object cloned_impl = impl.attr(operation)();
    if (!nb::isinstance<CustomSUNMatrix>(cloned_impl))
    {
      throw nb::type_error((std::string("CustomSUNMatrix.") + operation +
                            "() must return a CustomSUNMatrix")
                             .c_str());
    }

    auto* cloned = nb::cast<CustomSUNMatrix*>(cloned_impl);
    if (!source->sunctx_owner_ || !cloned->sunctx_owner_)
    {
      throw nb::type_error(
        "CustomSUNMatrix clone did not initialize the base constructor");
    }
    require_same_context(source->sunctx_owner_.get(), cloned->sunctx_owner_,
                         (std::string("CustomSUNMatrix.") + operation).c_str());
    return cloned_impl;
  }

  /* Recover the Python implementation for borrowed custom handles. */
  static nb::object _python_object_for(SUNMatrix A) noexcept
  {
    try
    {
      Content* content = get_content(A);
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
  static nb::object get_impl(SUNMatrix A)
  {
    return custom_content_impl(get_content(A), label);
  }

  static bool method_overridden(nb::handle impl, const char* name)
  {
    return custom_method_overridden<CustomSUNMatrix>(impl, name);
  }

  static void validate_required_methods(nb::handle impl)
  {
    if (!method_overridden(impl, "zero"))
    {
      throw nb::type_error("CustomSUNMatrix subclass must override zero()");
    }
  }

  static SUNErrCode call_status_method(SUNMatrix A, const char* name,
                                       const char* operation)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
      return static_cast<SUNErrCode>(nb::cast<int>(get_impl(A).attr(name)()));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(A ? A->sunctx : nullptr, operation,
                                 SUN_ERR_EXT_FAIL)
  }

  static int call_status_method(SUNMatrix A, const char* name, N_Vector x,
                                N_Vector y, const char* operation)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
      return nb::cast<int>(
        get_impl(A).attr(name)(nb::cast(x, nb::rv_policy::reference),
                               nb::cast(y, nb::rv_policy::reference)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(A ? A->sunctx : nullptr, operation,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode check_same_custom_type(SUNMatrix A, SUNMatrix B,
                                           nb::object& impl_A, nb::object& impl_B)
  {
    // Binary matrix operations are forwarded to Python objects, so require the
    // same concrete Python type before exposing one object to another.
    impl_A = get_impl(A);
    impl_B = get_impl(B);
    if (Py_TYPE(impl_A.ptr()) != Py_TYPE(impl_B.ptr()) &&
        !(get_content(A)->has_is_compatible &&
          nb::cast<bool>(impl_A.attr("is_compatible")(impl_B))))
    {
      throw nb::type_error(
        "custom SUNMatrix operands must have the same Python type");
    }
    return SUN_SUCCESS;
  }

  static SUNMatrix_ID custom_matrix_getid(SUNMatrix);

  static SUNMatrix custom_matrix_clone(SUNMatrix A);

  static SUNMatrix clone_raw_handle(
    nb::handle impl, SUNMatrix source,
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner);

  static void custom_matrix_destroy(SUNMatrix A);

  static SUNErrCode custom_matrix_zero(SUNMatrix A);

  static SUNErrCode custom_matrix_space(SUNMatrix A, long int* lenrw,
                                        long int* leniw);

  static SUNErrCode custom_matrix_copy(SUNMatrix A, SUNMatrix B);

  static SUNErrCode custom_matrix_scaleadd(sunrealtype c, SUNMatrix A,
                                           SUNMatrix B);

  static SUNErrCode custom_matrix_scaleaddi(sunrealtype c, SUNMatrix A);

  static SUNErrCode custom_matrix_matvecsetup(SUNMatrix A);

  static SUNErrCode custom_matrix_matvec(SUNMatrix A, N_Vector x, N_Vector y);

  static SUNErrCode custom_matrix_hermitian_transpose_matvec(SUNMatrix A,
                                                             N_Vector x,
                                                             N_Vector y);
};

} // namespace sundials4py

#endif // SUNDIALS4PY_MATRIX_CUSTOM_HPP
