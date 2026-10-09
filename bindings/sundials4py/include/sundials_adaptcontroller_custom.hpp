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

#ifndef SUNDIALS4PY_ADAPTCONTROLLER_CUSTOM_HPP
#define SUNDIALS4PY_ADAPTCONTROLLER_CUSTOM_HPP

#include <cstdint>
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
  : public CustomObjectBase<CustomSUNAdaptController,
                            std::remove_pointer_t<SUNAdaptController>>
{
public:
  CustomSUNAdaptController(std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx,
                           SUNAdaptController_Type controller_type)
    : CustomObjectBase(std::move(sunctx), "CustomSUNAdaptController"),
      controller_type_(controller_type)
  {}

  virtual ~CustomSUNAdaptController() = default;

  void validate(nb::handle self) const
  {
    validate_required_methods(self, controller_type_);
  }

  static SUNErrCode base_method_status(const char* name)
  {
    PyErr_SetString(PyExc_NotImplementedError, name);
    nb::raise_python_error();
    return SUN_ERR_EXT_FAIL;
  }

private:
  struct Content : CustomContentBase
  {
    // Store the requested controller kind so gettype and the typed estimate
    // trampoline remain available from the opaque C handle.
    SUNAdaptController_Type controller_type{SUN_ADAPTCONTROLLER_NONE};
  };

  static constexpr const char* label = "CustomSUNAdaptController";

  friend class CustomObjectBase<CustomSUNAdaptController,
                                std::remove_pointer_t<SUNAdaptController>>;

  std::shared_ptr<std::remove_pointer_t<SUNAdaptController>> make_handle(
    nb::handle impl);

  static std::shared_ptr<std::remove_pointer_t<SUNAdaptController>> make_handle(
    nb::handle impl,
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
    SUNAdaptController_Type controller_type);

  static Content* get_content(SUNAdaptController C)
  {
    if (!C || !C->ops || C->ops->destroy != custom_controller_destroy ||
        !C->content)
    {
      return nullptr;
    }
    return static_cast<Content*>(C->content);
  }

public:
  static nb::object _python_object_for(SUNAdaptController C) noexcept
  {
    try
    {
      Content* content = get_content(C);
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

  static SUNAdaptController_Type custom_controller_gettype(SUNAdaptController C);

  static SUNErrCode custom_controller_estimatestep(SUNAdaptController C,
                                                   sunrealtype h, int p,
                                                   sunrealtype dsm,
                                                   sunrealtype* hnew);

  static SUNErrCode custom_controller_estimatesteptol(
    SUNAdaptController C, sunrealtype H, sunrealtype tolfac, int P,
    sunrealtype DSM, sunrealtype dsm, sunrealtype* Hnew, sunrealtype* tolfacnew);

  template<typename... Args>
  static SUNErrCode call_status(SUNAdaptController C, const char* name,
                                const char* operation, Args&&... args)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
      return static_cast<SUNErrCode>(
        nb::cast<int>(get_impl(C).attr(name)(std::forward<Args>(args)...)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(C ? C->sunctx : nullptr, operation,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_controller_reset(SUNAdaptController C);

  static SUNErrCode custom_controller_setdefaults(SUNAdaptController C);

  static SUNErrCode custom_controller_setoptions(SUNAdaptController C,
                                                 const char* id,
                                                 const char* file_name,
                                                 int argc, char* argv[]);

  static SUNErrCode custom_controller_write(SUNAdaptController C, FILE* outfile);

  static SUNErrCode custom_controller_seterrorbias(SUNAdaptController C,
                                                   sunrealtype bias);

  static SUNErrCode custom_controller_updateh(SUNAdaptController C,
                                              sunrealtype h, sunrealtype dsm);

  static SUNErrCode custom_controller_updatemrihtol(SUNAdaptController C,
                                                    sunrealtype H,
                                                    sunrealtype tolfac,
                                                    sunrealtype DSM,
                                                    sunrealtype dsm);

  static SUNErrCode custom_controller_destroy(SUNAdaptController C);

  SUNAdaptController_Type controller_type_;
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

#endif // SUNDIALS4PY_ADAPTCONTROLLER_CUSTOM_HPP
