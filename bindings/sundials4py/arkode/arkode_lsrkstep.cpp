/*------------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 *----------------------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <arkode/arkode.hpp>
#include <arkode/arkode_lsrkstep.h>
#include <sundials/sundials_core.hpp>

#include "arkode_impl.h"
#include "arkode_usersupplied.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_arkode_lsrkstep(nb::module_& m)
{
#include "arkode_lsrkstep_generated.hpp"

  m.def(
    "LSRKStepCreateSTS",
    [](std::function<std::remove_pointer_t<ARKRhsFn>> rhs, sunrealtype t0,
       N_Vector y0, SUNContext sunctx)
    {
      if (!rhs) { throw sundials4py::illegal_value("rhs was None"); }

      void* ark_mem = LSRKStepCreateSTS(lsrkstep_f_wrapper, t0, y0, sunctx);
      if (ark_mem == nullptr)
      {
        throw sundials4py::error_returned("Failed to create LSRKStep memory");
      }

      auto fn_table = new arkode_user_supplied_fn_table;

      static_cast<ARKodeMem>(ark_mem)->python = fn_table;

      int ark_status = ARKodeSetUserData(ark_mem, ark_mem);
      if (ark_status != ARK_SUCCESS)
      {
        free(fn_table);
        throw sundials4py::error_returned(
          "Failed to set user data in ARKODE memory");
      }

      fn_table->lsrkstep_f = nb::cast(rhs);

      return std::make_shared<ARKodeView>(ark_mem);
    },
    nb::arg("rhs"), nb::arg("t0"), nb::arg("y0"), nb::arg("sunctx"),
    nb::keep_alive<0, 4>());

  m.def(
    "LSRKStepCreateSSP",
    [](std::function<std::remove_pointer_t<ARKRhsFn>> rhs, sunrealtype t0,
       N_Vector y0, SUNContext sunctx)
    {
      if (!rhs) { throw sundials4py::illegal_value("rhs was null"); }

      void* ark_mem = LSRKStepCreateSSP(lsrkstep_f_wrapper, t0, y0, sunctx);
      if (ark_mem == nullptr)
      {
        throw sundials4py::error_returned("Failed to create LSRKStep memory");
      }

      auto fn_table = new arkode_user_supplied_fn_table;

      static_cast<ARKodeMem>(ark_mem)->python = fn_table;

      int ark_status = ARKodeSetUserData(ark_mem, ark_mem);
      if (ark_status != ARK_SUCCESS)
      {
        free(fn_table);
        throw sundials4py::error_returned(
          "Failed to set user data in ARKODE memory");
      }

      fn_table->lsrkstep_f = nb::cast(rhs);

      return std::make_shared<ARKodeView>(ark_mem);
    },
    nb::arg("rhs"), nb::arg("t0"), nb::arg("y0"), nb::arg("sunctx"),
    nb::keep_alive<0, 4>());

  m.def(
    "LSRKStepSetDomEigFn",
    [](void* ark_mem, std::function<std::remove_pointer_t<ARKDomEigFn>> fn)
    {
      auto fn_table             = get_arkode_fn_table(ark_mem);
      fn_table->lsrkstep_domeig = nb::cast(fn);
      if (fn) { return LSRKStepSetDomEigFn(ark_mem, lsrkstep_domeig_wrapper); }
      else { return LSRKStepSetDomEigFn(ark_mem, nullptr); }
    },
    nb::arg("arkode_mem"), nb::arg("eig_fn"));
}

} // namespace sundials4py
