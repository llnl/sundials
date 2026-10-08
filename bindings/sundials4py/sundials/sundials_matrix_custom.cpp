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

std::shared_ptr<std::remove_pointer_t<SUNMatrix>> CustomSUNMatrix::make_handle(
  nb::handle impl)
{ return make_handle(impl, sunctx_owner_, Ownership::weak); }

SUNMatrix CustomSUNMatrix::create_raw_handle(
  nb::handle impl,
  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
  Ownership ownership)
{
  validate_required_methods(impl);

  // Complete every potentially throwing Python operation before allocating
  // the native shell. Once A exists, the remaining assignments are noexcept.
  const bool has_clone       = method_overridden(impl, "clone");
  const bool has_copy        = method_overridden(impl, "copy");
  const bool has_scaleadd    = method_overridden(impl, "scaleadd");
  const bool has_scaleaddi   = method_overridden(impl, "scaleaddi");
  const bool has_matvec      = method_overridden(impl, "matvec");
  const bool has_matvecsetup = method_overridden(impl, "matvecsetup");
  const bool has_hermitian   = method_overridden(impl,
                                                 "hermitian_transpose_matvec");
  const bool has_is_compatible = method_overridden(impl, "is_compatible");

  std::unique_ptr<Content> content(new Content);
  CustomSUNMatrix* owner = nullptr;
  if (ownership == Ownership::strong)
  {
    owner          = nb::cast<CustomSUNMatrix*>(impl);
    content->owner = owner;
  }
  content->sunctx_owner      = std::move(sunctx_owner);
  content->has_is_compatible = has_is_compatible;
  // User-created handles should not extend the Python object's lifetime.
  // Clones returned to SUNDIALS are C-owned, so they hold a strong reference.
  if (ownership == Ownership::weak) { content->weak_impl = nb::weakref(impl); }
  else
  {
    content->strong_impl = nb::borrow<nb::object>(impl);
  }

  SUNMatrix A = SUNMatNewEmpty(content->sunctx_owner.get());
  if (!A) { throw error_returned("SUNMatNewEmpty failed"); }

  // getid, destroy, and zero are universal. Other operations are exposed only
  // when the subclass overrides them, allowing a minimal KINSOL/IDA matrix.
  A->content = content.release();
  if (owner) { owner->raw_handle_ = A; }
  A->ops->getid   = custom_matrix_getid;
  A->ops->destroy = custom_matrix_destroy;
  A->ops->zero    = custom_matrix_zero;

  if (has_clone) { A->ops->clone = custom_matrix_clone; }
  if (has_copy) { A->ops->copy = custom_matrix_copy; }
  if (has_scaleadd) { A->ops->scaleadd = custom_matrix_scaleadd; }
  if (has_scaleaddi) { A->ops->scaleaddi = custom_matrix_scaleaddi; }
  if (has_matvec) { A->ops->matvec = custom_matrix_matvec; }
  if (has_matvecsetup) { A->ops->matvecsetup = custom_matrix_matvecsetup; }
  if (has_hermitian)
  {
    A->ops->mathermitiantransposevec = custom_matrix_hermitian_transpose_matvec;
  }

  return A;
}

std::shared_ptr<std::remove_pointer_t<SUNMatrix>> CustomSUNMatrix::make_handle(
  nb::handle impl,
  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
  Ownership ownership)
{
  SUNMatrix A = create_raw_handle(impl, std::move(sunctx_owner), ownership);
  return sundials::experimental::our_make_shared<
    std::remove_pointer_t<SUNMatrix>, sundials::experimental::SUNMatrixDeleter>(A);
}

SUNMatrix_ID CustomSUNMatrix::custom_matrix_getid(SUNMatrix)
{ return SUNMATRIX_CUSTOM; }

SUNMatrix CustomSUNMatrix::custom_matrix_clone(SUNMatrix A)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return nullptr; }
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
    require_same_context(A->sunctx, cloned->sunctx_owner_,
                         "CustomSUNMatrix.clone");
    if (Py_TYPE(cloned_impl.ptr()) == Py_TYPE(impl.ptr()))
    {
      return clone_raw_handle(cloned_impl, A, cloned->sunctx_owner_);
    }
    // Native SUNDIALS owns clones, so the returned handle must keep its
    // Python implementation alive until SUNMatDestroy is called.
    return create_raw_handle(cloned_impl, cloned->sunctx_owner_,
                             Ownership::strong);
  }
  SUNDIALS4PY_CATCH_AND_REPORT(A ? A->sunctx : nullptr, "CustomSUNMatrix.clone",
                               nullptr)
}

SUNMatrix CustomSUNMatrix::clone_raw_handle(
  nb::handle impl, SUNMatrix source,
  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner)
{
  Content* source_content = get_content(source);
  if (!source_content)
  {
    throw nb::type_error("source is not a custom SUNMatrix");
  }

  auto* owner = nb::cast<CustomSUNMatrix*>(impl);
  std::unique_ptr<Content> content(new Content);
  content->owner             = owner;
  content->sunctx_owner      = std::move(sunctx_owner);
  content->has_is_compatible = source_content->has_is_compatible;
  content->strong_impl       = nb::borrow<nb::object>(impl);

  SUNMatrix A = SUNMatNewEmpty(content->sunctx_owner.get());
  if (!A) { throw error_returned("SUNMatNewEmpty failed"); }
  if (SUNMatCopyOps(source, A) != SUN_SUCCESS)
  {
    SUNMatFreeEmpty(A);
    throw error_returned("SUNMatCopyOps failed");
  }

  A->content         = content.release();
  owner->raw_handle_ = A;
  return A;
}

void CustomSUNMatrix::custom_matrix_destroy(SUNMatrix A)
{
  if (!A) { return; }
  Content* content = get_content(A);
  if (content && content->owner && content->owner->raw_handle_ == A)
  {
    content->owner->raw_handle_ = nullptr;
  }
  custom_content_destroy(content);
  A->content = nullptr;
  SUNMatFreeEmpty(A);
}

SUNErrCode CustomSUNMatrix::custom_matrix_zero(SUNMatrix A)
{ return call_status_method(A, "zero", "CustomSUNMatrix.zero"); }

SUNErrCode CustomSUNMatrix::custom_matrix_copy(SUNMatrix A, SUNMatrix B)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    nb::object impl_A;
    nb::object impl_B;
    check_same_custom_type(A, B, impl_A, impl_B);
    return static_cast<SUNErrCode>(nb::cast<int>(impl_A.attr("copy")(impl_B)));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(A ? A->sunctx : nullptr, "CustomSUNMatrix.copy",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNMatrix::custom_matrix_scaleadd(sunrealtype c, SUNMatrix A,
                                                   SUNMatrix B)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    nb::object impl_A;
    nb::object impl_B;
    check_same_custom_type(A, B, impl_A, impl_B);
    return static_cast<SUNErrCode>(
      nb::cast<int>(impl_A.attr("scaleadd")(c, impl_B)));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(A ? A->sunctx : nullptr,
                               "CustomSUNMatrix.scaleadd", SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNMatrix::custom_matrix_scaleaddi(sunrealtype c, SUNMatrix A)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    return static_cast<SUNErrCode>(
      nb::cast<int>(get_impl(A).attr("scaleaddi")(c)));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(A ? A->sunctx : nullptr,
                               "CustomSUNMatrix.scaleaddi", SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNMatrix::custom_matrix_matvecsetup(SUNMatrix A)
{ return call_status_method(A, "matvecsetup", "CustomSUNMatrix.matvecsetup"); }

SUNErrCode CustomSUNMatrix::custom_matrix_matvec(SUNMatrix A, N_Vector x,
                                                 N_Vector y)
{
  return static_cast<SUNErrCode>(
    call_status_method(A, "matvec", x, y, "CustomSUNMatrix.matvec"));
}

SUNErrCode CustomSUNMatrix::custom_matrix_hermitian_transpose_matvec(SUNMatrix A,
                                                                     N_Vector x,
                                                                     N_Vector y)
{
  return static_cast<SUNErrCode>(
    call_status_method(A, "hermitian_transpose_matvec", x, y,
                       "CustomSUNMatrix.hermitian_transpose_matvec"));
}

} // namespace sundials4py
