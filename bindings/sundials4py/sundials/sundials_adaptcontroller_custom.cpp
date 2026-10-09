/*-----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *----------------------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace sundials4py {

std::shared_ptr<std::remove_pointer_t<SUNAdaptController>> CustomSUNAdaptController::make_handle(
  nb::handle impl)
{
  return make_handle(impl, sunctx_owner_, controller_type_);
}

std::shared_ptr<std::remove_pointer_t<SUNAdaptController>> CustomSUNAdaptController::make_handle(
  nb::handle impl,
  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
  SUNAdaptController_Type controller_type)
{
  validate_required_methods(impl, controller_type);

  // Resolve Python overrides before allocating the native shell so a failed
  // lookup leaves no C allocation to unwind.
  const bool has_reset            = method_overridden(impl, "reset");
  const bool has_set_defaults     = method_overridden(impl, "set_defaults");
  const bool has_set_options      = method_overridden(impl, "set_options");
  const bool has_write            = method_overridden(impl, "write");
  const bool has_set_error_bias   = method_overridden(impl, "set_error_bias");
  const bool has_update_h         = method_overridden(impl, "update_h");
  const bool has_update_mri_h_tol = method_overridden(impl, "update_mri_h_tol");

  // Common optional operations are shared; the estimate operation is selected
  // by controller type because H and MRI controllers have different C APIs.
  // Keep the content under RAII until the native shell is ready to own it.
  std::unique_ptr<Content> content(new Content);
  content->weak_impl       = nb::weakref(impl);
  content->sunctx_owner    = std::move(sunctx_owner);
  content->controller_type = controller_type;

  SUNAdaptController C = SUNAdaptController_NewEmpty(content->sunctx_owner.get());
  if (!C) { throw error_returned("SUNAdaptController_NewEmpty failed"); }

  C->content      = content.release();
  C->ops->gettype = custom_controller_gettype;
  C->ops->destroy = custom_controller_destroy;
  if (controller_type == SUN_ADAPTCONTROLLER_H)
  {
    C->ops->estimatestep = custom_controller_estimatestep;
  }
  if (controller_type == SUN_ADAPTCONTROLLER_MRI_H_TOL)
  {
    C->ops->estimatesteptol = custom_controller_estimatesteptol;
  }

  if (has_reset) { C->ops->reset = custom_controller_reset; }
  if (has_set_options) { C->ops->setoptions = custom_controller_setoptions; }
  if (has_set_defaults) { C->ops->setdefaults = custom_controller_setdefaults; }
  if (has_write) { C->ops->write = custom_controller_write; }
  if (has_set_error_bias)
  {
    C->ops->seterrorbias = custom_controller_seterrorbias;
  }
  if (has_update_h) { C->ops->updateh = custom_controller_updateh; }
  if (has_update_mri_h_tol)
  {
    C->ops->updatemrihtol = custom_controller_updatemrihtol;
  }

  return sundials::experimental::our_make_shared<
    std::remove_pointer_t<SUNAdaptController>,
    sundials::experimental::SUNAdaptControllerDeleter>(C);
}

SUNAdaptController_Type CustomSUNAdaptController::custom_controller_gettype(
  SUNAdaptController C)
{
  Content* content = get_content(C);
  return content ? content->controller_type : SUN_ADAPTCONTROLLER_NONE;
}

SUNErrCode CustomSUNAdaptController::custom_controller_estimatestep(
  SUNAdaptController C, sunrealtype h, int p, sunrealtype dsm, sunrealtype* hnew)
{
  // Python returns both the status code and output value, matching the public
  // binding style for C routines with output pointers.
  *hnew = SUN_RCONST(0.0);
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    auto result = nb::cast<std::tuple<int, sunrealtype>>(
      get_impl(C).attr("estimate_step")(h, p, dsm));
    *hnew = std::get<1>(result);
    return static_cast<SUNErrCode>(std::get<0>(result));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(C ? C->sunctx : nullptr,
                               "CustomSUNAdaptController.estimate_step",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNAdaptController::custom_controller_estimatesteptol(
  SUNAdaptController C, sunrealtype H, sunrealtype tolfac, int P,
  sunrealtype DSM, sunrealtype dsm, sunrealtype* Hnew, sunrealtype* tolfacnew)
{
  // MRI controllers produce two output values, so the Python method returns a
  // three-tuple: (status, Hnew, tolfacnew).
  *Hnew      = SUN_RCONST(0.0);
  *tolfacnew = SUN_RCONST(0.0);
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    auto result = nb::cast<std::tuple<int, sunrealtype, sunrealtype>>(
      get_impl(C).attr("estimate_step_tol")(H, tolfac, P, DSM, dsm));
    *Hnew      = std::get<1>(result);
    *tolfacnew = std::get<2>(result);
    return static_cast<SUNErrCode>(std::get<0>(result));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(C ? C->sunctx : nullptr,
                               "CustomSUNAdaptController.estimate_step_tol",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNAdaptController::custom_controller_reset(SUNAdaptController C)
{
  return call_status(C, "reset", "CustomSUNAdaptController.reset");
}

SUNErrCode CustomSUNAdaptController::custom_controller_setdefaults(
  SUNAdaptController C)
{
  return call_status(C, "set_defaults", "CustomSUNAdaptController.set_defaults");
}

SUNErrCode CustomSUNAdaptController::custom_controller_setoptions(
  SUNAdaptController C, const char* id, const char* file_name, int argc,
  char* argv[])
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    std::vector<std::string> args;
    args.reserve(static_cast<size_t>(argc));
    for (int i = 0; i < argc; ++i) { args.emplace_back(argv[i]); }
    return static_cast<SUNErrCode>(nb::cast<int>(get_impl(C).attr(
      "set_options")(id ? id : "", file_name ? file_name : "", args)));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(C ? C->sunctx : nullptr,
                               "CustomSUNAdaptController.set_options",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNAdaptController::custom_controller_write(SUNAdaptController C,
                                                             FILE* outfile)
{
  return call_status(C, "write", "CustomSUNAdaptController.write",
                     reinterpret_cast<std::uintptr_t>(outfile));
}

SUNErrCode CustomSUNAdaptController::custom_controller_seterrorbias(
  SUNAdaptController C, sunrealtype bias)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    return static_cast<SUNErrCode>(
      nb::cast<int>(get_impl(C).attr("set_error_bias")(bias)));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(C ? C->sunctx : nullptr,
                               "CustomSUNAdaptController.set_error_bias",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNAdaptController::custom_controller_updateh(
  SUNAdaptController C, sunrealtype h, sunrealtype dsm)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    return static_cast<SUNErrCode>(
      nb::cast<int>(get_impl(C).attr("update_h")(h, dsm)));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(C ? C->sunctx : nullptr,
                               "CustomSUNAdaptController.update_h",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNAdaptController::custom_controller_updatemrihtol(
  SUNAdaptController C, sunrealtype H, sunrealtype tolfac, sunrealtype DSM,
  sunrealtype dsm)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    return static_cast<SUNErrCode>(
      nb::cast<int>(get_impl(C).attr("update_mri_h_tol")(H, tolfac, DSM, dsm)));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(C ? C->sunctx : nullptr,
                               "CustomSUNAdaptController.update_mri_h_tol",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNAdaptController::custom_controller_destroy(SUNAdaptController C)
{
  if (!C) { return SUN_SUCCESS; }
  Content* content = get_content(C);
  custom_content_destroy(content);
  C->content = nullptr;
  SUNAdaptController_DestroyEmpty(C);
  return SUN_SUCCESS;
}

} // namespace sundials4py
