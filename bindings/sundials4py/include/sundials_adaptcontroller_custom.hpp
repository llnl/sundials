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
 * CustomSUNAdaptController and its two concrete flavors, CustomSUNHController and
 * CustomSUNMRIController.
 *----------------------------------------------------------------------------*/

#ifndef _SUNDIALS4PY_ADAPTCONTROLLER_CUSTOM_HPP
#define _SUNDIALS4PY_ADAPTCONTROLLER_CUSTOM_HPP

#include <memory>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>

#include <sundials/sundials_adaptcontroller.h>
#include <sundials/sundials_adaptcontroller.hpp>
/* sundials_adaptcontroller.hpp only declares the deleter, so the error codes and
   our_make_shared are requested explicitly rather than relied on transitively. */
#include <sundials/sundials_classview.hpp>
#include <sundials/sundials_errors.h>

#include "sundials4py_core_types.hpp"
#include "sundials4py_custom_object.hpp"

namespace sundials4py {

/*
 * Python-owned SUNAdaptController implementation.
 *
 * H and MRI controllers share the same native content and operation plumbing;
 * the controller type determines which estimate callback is mandatory and which
 * SUNDIALS vtable slot is populated.
 */
class CustomSUNAdaptController
{
public:
  enum class HandleState
  {
    unmaterialized,
    materializing,
    materialized
  };

  CustomSUNAdaptController(std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx,
                           SUNAdaptController_Type controller_type)
    : sunctx_owner_(std::move(sunctx)), controller_type_(controller_type)
  {}

  virtual ~CustomSUNAdaptController() = default;

  std::shared_ptr<std::remove_pointer_t<SUNAdaptController>> _get_sundials_handle(
    nb::handle self)
  {
    if (!sunctx_owner_)
    {
      throw nb::type_error(
        "CustomSUNAdaptController base constructor was not initialized");
    }
    if (state_ == HandleState::materializing)
    {
      throw std::runtime_error(
        "reentrant CustomSUNAdaptController native handle materialization");
    }
    if (state_ == HandleState::materialized) { return controller_; }

    // Lazily create the native controller the first time a SUNDIALS wrapper
    // needs SUNAdaptController.
    state_ = HandleState::materializing;
    try
    {
      controller_ = make_handle(self, sunctx_owner_, controller_type_);
      materialization_count_++;
      state_ = HandleState::materialized;
    }
    catch (...)
    {
      controller_.reset();
      state_ = HandleState::unmaterialized;
      throw;
    }
    return controller_;
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
  struct Content
  {
    SUNDIALS4PY_CUSTOM_CONTENT_MEMBERS(0x53554e4144505931ULL);

    // Store the requested controller kind so gettype and the typed estimate
    // trampoline remain available from the opaque C handle.
    SUNAdaptController_Type controller_type{SUN_ADAPTCONTROLLER_NONE};
  };

  static constexpr const char* label = "CustomSUNAdaptController";

  static std::shared_ptr<std::remove_pointer_t<SUNAdaptController>> make_handle(
    nb::handle impl,
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
    SUNAdaptController_Type controller_type)
  {
    validate_required_methods(impl, controller_type);

    // Resolve Python overrides before allocating the native shell so a failed
    // lookup leaves no C allocation to unwind.
    const bool has_reset            = method_overridden(impl, "reset");
    const bool has_set_defaults     = method_overridden(impl, "set_defaults");
    const bool has_set_error_bias   = method_overridden(impl, "set_error_bias");
    const bool has_update_h         = method_overridden(impl, "update_h");
    const bool has_update_mri_h_tol = method_overridden(impl,
                                                        "update_mri_h_tol");

    // Common optional operations are shared; the estimate operation is selected
    // by controller type because H and MRI controllers have different C APIs.
    std::unique_ptr<Content> content(new Content);
    content->weak_impl       = nb::weakref(impl);
    content->sunctx_owner    = std::move(sunctx_owner);
    content->controller_type = controller_type;

    SUNAdaptController C =
      SUNAdaptController_NewEmpty(content->sunctx_owner.get());
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
    if (has_set_defaults)
    {
      C->ops->setdefaults = custom_controller_setdefaults;
    }
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

  static Content* get_content(SUNAdaptController C)
  {
    return custom_content_cast<Content>(C);
  }

  static nb::object get_impl(SUNAdaptController C)
  {
    return custom_content_impl(get_content(C), label);
  }

  static bool method_overridden(nb::handle impl, const char* name)
  {
    return custom_method_overridden<CustomSUNAdaptController>(impl, name);
  }

  static void validate_required_methods(nb::handle impl,
                                        SUNAdaptController_Type controller_type)
  {
    if (controller_type == SUN_ADAPTCONTROLLER_H &&
        !method_overridden(impl, "estimate_step"))
    {
      throw nb::type_error(
        "CustomSUNHController subclass must override estimate_step()");
    }
    if (controller_type == SUN_ADAPTCONTROLLER_MRI_H_TOL &&
        !method_overridden(impl, "estimate_step_tol"))
    {
      throw nb::type_error(
        "CustomSUNMRIController subclass must override estimate_step_tol()");
    }
  }

  static SUNAdaptController_Type custom_controller_gettype(SUNAdaptController C)
  {
    Content* content = get_content(C);
    return content ? content->controller_type : SUN_ADAPTCONTROLLER_NONE;
  }

  static SUNErrCode custom_controller_estimatestep(SUNAdaptController C,
                                                   sunrealtype h, int p,
                                                   sunrealtype dsm,
                                                   sunrealtype* hnew)
  {
    // Python returns both the status code and output value, matching the public
    // binding style for C routines with output pointers.
    *hnew = SUN_RCONST(0.0);
    try
    {
      nb::gil_scoped_acquire gil;
      auto result = nb::cast<std::tuple<int, sunrealtype>>(
        get_impl(C).attr("estimate_step")(h, p, dsm));
      *hnew = std::get<1>(result);
      return static_cast<SUNErrCode>(std::get<0>(result));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(C ? C->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_controller_estimatesteptol(
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
      auto result = nb::cast<std::tuple<int, sunrealtype, sunrealtype>>(
        get_impl(C).attr("estimate_step_tol")(H, tolfac, P, DSM, dsm));
      *Hnew      = std::get<1>(result);
      *tolfacnew = std::get<2>(result);
      return static_cast<SUNErrCode>(std::get<0>(result));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(C ? C->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode call_status(SUNAdaptController C, const char* name)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      return static_cast<SUNErrCode>(nb::cast<int>(get_impl(C).attr(name)()));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(C ? C->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_controller_reset(SUNAdaptController C)
  {
    return call_status(C, "reset");
  }

  static SUNErrCode custom_controller_setdefaults(SUNAdaptController C)
  {
    return call_status(C, "set_defaults");
  }

  static SUNErrCode custom_controller_seterrorbias(SUNAdaptController C,
                                                   sunrealtype bias)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      return static_cast<SUNErrCode>(
        nb::cast<int>(get_impl(C).attr("set_error_bias")(bias)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(C ? C->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_controller_updateh(SUNAdaptController C,
                                              sunrealtype h, sunrealtype dsm)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      return static_cast<SUNErrCode>(
        nb::cast<int>(get_impl(C).attr("update_h")(h, dsm)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(C ? C->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_controller_updatemrihtol(SUNAdaptController C,
                                                    sunrealtype H,
                                                    sunrealtype tolfac,
                                                    sunrealtype DSM,
                                                    sunrealtype dsm)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      return static_cast<SUNErrCode>(nb::cast<int>(
        get_impl(C).attr("update_mri_h_tol")(H, tolfac, DSM, dsm)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(C ? C->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_controller_destroy(SUNAdaptController C)
  {
    if (!C) { return SUN_SUCCESS; }
    Content* content = get_content(C);
    custom_content_destroy(content);
    C->content = nullptr;
    SUNAdaptController_DestroyEmpty(C);
    return SUN_SUCCESS;
  }

  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner_;
  SUNAdaptController_Type controller_type_;
  std::shared_ptr<std::remove_pointer_t<SUNAdaptController>> controller_;
  HandleState state_{HandleState::unmaterialized};
  int materialization_count_{0};
};

/* Time step (H) controller: must implement estimate_step(). */
class CustomSUNHController : public CustomSUNAdaptController
{
public:
  explicit CustomSUNHController(
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx)
    : CustomSUNAdaptController(std::move(sunctx), SUN_ADAPTCONTROLLER_H)
  {}
};

/* Multirate (H, tolerance) controller: must implement estimate_step_tol(). */
class CustomSUNMRIController : public CustomSUNAdaptController
{
public:
  explicit CustomSUNMRIController(
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx)
    : CustomSUNAdaptController(std::move(sunctx), SUN_ADAPTCONTROLLER_MRI_H_TOL)
  {}
};

} // namespace sundials4py

#endif // _SUNDIALS4PY_ADAPTCONTROLLER_CUSTOM_HPP
