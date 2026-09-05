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

#ifndef _SUNDIALS4PY_MATRIX_CUSTOM_HPP
#define _SUNDIALS4PY_MATRIX_CUSTOM_HPP

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
{
public:
  enum class HandleState
  {
    unmaterialized,
    materializing,
    materialized
  };

  explicit CustomSUNMatrix(std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx)
    : sunctx_owner_(std::move(sunctx))
  {}

  std::shared_ptr<std::remove_pointer_t<SUNMatrix>> _get_sundials_handle(
    nb::handle self)
  {
    if (!sunctx_owner_)
    {
      throw nb::type_error(
        "CustomSUNMatrix base constructor was not initialized");
    }

    if (state_ == HandleState::materializing)
    {
      throw std::runtime_error(
        "reentrant CustomSUNMatrix native handle materialization");
    }

    if (state_ == HandleState::materialized) { return matrix_; }

    // Materialization is transactional: failures leave the Python object in the
    // unmaterialized state so a corrected subclass can be tried again.
    state_ = HandleState::materializing;
    try
    {
      matrix_ = make_handle(self, sunctx_owner_, Ownership::weak);
      materialization_count_++;
      state_ = HandleState::materialized;
    }
    catch (...)
    {
      matrix_.reset();
      state_ = HandleState::unmaterialized;
      throw;
    }

    return matrix_;
  }

  int _materialization_count() const { return materialization_count_; }

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
  enum class Ownership
  {
    weak,
    strong
  };

  struct Content
  {
    SUNDIALS4PY_CUSTOM_CONTENT_MEMBERS(0x53554e4d41545059ULL);
  };

  static constexpr const char* label = "CustomSUNMatrix";

  static SUNMatrix create_raw_handle(
    nb::handle impl,
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
    Ownership ownership)
  {
    validate_required_methods(impl);

    // Complete every potentially throwing Python operation before allocating
    // the native shell. Once A exists, the remaining assignments are noexcept.
    const bool has_matvecsetup = method_overridden(impl, "matvecsetup");
    const bool has_hermitian   = method_overridden(impl,
                                                   "hermitian_transpose_matvec");

    std::unique_ptr<Content> content(new Content);
    content->sunctx_owner = std::move(sunctx_owner);
    // User-created handles should not extend the Python object's lifetime.
    // Clones returned to SUNDIALS are C-owned, so they hold a strong reference.
    if (ownership == Ownership::weak)
    {
      content->weak_impl = nb::weakref(impl);
    }
    else { content->strong_impl = nb::borrow<nb::object>(impl); }

    SUNMatrix A = SUNMatNewEmpty(content->sunctx_owner.get());
    if (!A) { throw error_returned("SUNMatNewEmpty failed"); }

    // Required operations are always installed; optional operations are exposed
    // only when the Python subclass overrides the corresponding base method.
    A->content        = content.release();
    A->ops->getid     = custom_matrix_getid;
    A->ops->clone     = custom_matrix_clone;
    A->ops->destroy   = custom_matrix_destroy;
    A->ops->zero      = custom_matrix_zero;
    A->ops->copy      = custom_matrix_copy;
    A->ops->scaleadd  = custom_matrix_scaleadd;
    A->ops->scaleaddi = custom_matrix_scaleaddi;
    A->ops->matvec    = custom_matrix_matvec;

    if (has_matvecsetup) { A->ops->matvecsetup = custom_matrix_matvecsetup; }
    if (has_hermitian)
    {
      A->ops->mathermitiantransposevec = custom_matrix_hermitian_transpose_matvec;
    }

    return A;
  }

  static std::shared_ptr<std::remove_pointer_t<SUNMatrix>> make_handle(
    nb::handle impl,
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
    Ownership ownership)
  {
    SUNMatrix A = create_raw_handle(impl, std::move(sunctx_owner), ownership);
    return sundials::experimental::our_make_shared<
      std::remove_pointer_t<SUNMatrix>, sundials::experimental::SUNMatrixDeleter>(
      A);
  }

  static Content* get_content(SUNMatrix A)
  {
    return custom_content_cast<Content>(A);
  }

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
    static constexpr const char* required[] = {"clone",     "zero",
                                               "copy",      "scaleadd",
                                               "scaleaddi", "matvec"};

    for (const char* name : required)
    {
      if (!method_overridden(impl, name))
      {
        throw nb::type_error(
          (std::string("CustomSUNMatrix subclass must override ") + name + "()")
            .c_str());
      }
    }
  }

  static SUNErrCode call_status_method(SUNMatrix A, const char* name,
                                       const char* operation)
  {
    try
    {
      nb::gil_scoped_acquire gil;
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
    if (Py_TYPE(impl_A.ptr()) != Py_TYPE(impl_B.ptr()))
    {
      throw nb::type_error(
        "custom SUNMatrix operands must have the same Python type");
    }
    return SUN_SUCCESS;
  }

  static SUNMatrix_ID custom_matrix_getid(SUNMatrix)
  {
    return SUNMATRIX_CUSTOM;
  }

  static SUNMatrix custom_matrix_clone(SUNMatrix A)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      nb::object impl        = get_impl(A);
      nb::object cloned_impl = impl.attr("clone")();
      if (!nb::isinstance<CustomSUNMatrix>(cloned_impl))
      {
        throw nb::type_error(
          "CustomSUNMatrix.clone() must return a CustomSUNMatrix");
      }

      auto* cloned = nb::cast<CustomSUNMatrix*>(cloned_impl);
      if (!cloned->sunctx_owner_)
      {
        throw nb::type_error(
          "CustomSUNMatrix clone did not initialize the base constructor");
      }
      // Native SUNDIALS owns clones, so the returned handle must keep its
      // Python implementation alive until SUNMatDestroy is called.
      return create_raw_handle(cloned_impl, cloned->sunctx_owner_,
                               Ownership::strong);
    }
    SUNDIALS4PY_CATCH_AND_REPORT(A ? A->sunctx : nullptr,
                                 "CustomSUNMatrix.clone", nullptr)
  }

  static void custom_matrix_destroy(SUNMatrix A)
  {
    if (!A) { return; }
    Content* content = get_content(A);
    custom_content_destroy(content);
    A->content = nullptr;
    SUNMatFreeEmpty(A);
  }

  static SUNErrCode custom_matrix_zero(SUNMatrix A)
  {
    return call_status_method(A, "zero", "CustomSUNMatrix.zero");
  }

  static SUNErrCode custom_matrix_copy(SUNMatrix A, SUNMatrix B)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      nb::object impl_A;
      nb::object impl_B;
      check_same_custom_type(A, B, impl_A, impl_B);
      return static_cast<SUNErrCode>(nb::cast<int>(impl_A.attr("copy")(impl_B)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(A ? A->sunctx : nullptr,
                                 "CustomSUNMatrix.copy", SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_matrix_scaleadd(sunrealtype c, SUNMatrix A, SUNMatrix B)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      nb::object impl_A;
      nb::object impl_B;
      check_same_custom_type(A, B, impl_A, impl_B);
      return static_cast<SUNErrCode>(
        nb::cast<int>(impl_A.attr("scaleadd")(c, impl_B)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(A ? A->sunctx : nullptr,
                                 "CustomSUNMatrix.scaleadd", SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_matrix_scaleaddi(sunrealtype c, SUNMatrix A)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      return static_cast<SUNErrCode>(
        nb::cast<int>(get_impl(A).attr("scaleaddi")(c)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(A ? A->sunctx : nullptr,
                                 "CustomSUNMatrix.scaleaddi", SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_matrix_matvecsetup(SUNMatrix A)
  {
    return call_status_method(A, "matvecsetup", "CustomSUNMatrix.matvecsetup");
  }

  static SUNErrCode custom_matrix_matvec(SUNMatrix A, N_Vector x, N_Vector y)
  {
    return static_cast<SUNErrCode>(
      call_status_method(A, "matvec", x, y, "CustomSUNMatrix.matvec"));
  }

  static SUNErrCode custom_matrix_hermitian_transpose_matvec(SUNMatrix A,
                                                             N_Vector x,
                                                             N_Vector y)
  {
    return static_cast<SUNErrCode>(
      call_status_method(A, "hermitian_transpose_matvec", x, y,
                         "CustomSUNMatrix.hermitian_transpose_matvec"));
  }

  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner_;
  std::shared_ptr<std::remove_pointer_t<SUNMatrix>> matrix_;
  HandleState state_{HandleState::unmaterialized};
  int materialization_count_{0};
};

} // namespace sundials4py

#endif // _SUNDIALS4PY_MATRIX_CUSTOM_HPP
