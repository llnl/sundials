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
#include <arkode/arkode_sprkstep.hpp>
#include <sundials/sundials_core.hpp>

#include "arkode_usersupplied.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_arkode_sprkstep(nb::module_& m)
{
#include "arkode_sprkstep_generated.hpp"

  /////////////////////////////////////////////////////////////////////////////
  // Manual attributes to capture static const int declarations that nanobind does not catch
  /////////////////////////////////////////////////////////////////////////////

  m.attr("SPRKSTEP_DEFAULT_1")  = static_cast<int>(SPRKSTEP_DEFAULT_1);
  m.attr("SPRKSTEP_DEFAULT_2")  = static_cast<int>(SPRKSTEP_DEFAULT_2);
  m.attr("SPRKSTEP_DEFAULT_3")  = static_cast<int>(SPRKSTEP_DEFAULT_3);
  m.attr("SPRKSTEP_DEFAULT_4")  = static_cast<int>(SPRKSTEP_DEFAULT_4);
  m.attr("SPRKSTEP_DEFAULT_5")  = static_cast<int>(SPRKSTEP_DEFAULT_5);
  m.attr("SPRKSTEP_DEFAULT_6")  = static_cast<int>(SPRKSTEP_DEFAULT_6);
  m.attr("SPRKSTEP_DEFAULT_8")  = static_cast<int>(SPRKSTEP_DEFAULT_8);
  m.attr("SPRKSTEP_DEFAULT_10") = static_cast<int>(SPRKSTEP_DEFAULT_10);

  /////////////////////////////////////////////////////////////////////////////
  // SPRKStep user-supplied function setters
  /////////////////////////////////////////////////////////////////////////////

  m.def(
    "SPRKStepCreate",
    [](std::function<std::remove_pointer_t<ARKRhsFn>> f1,
       std::function<std::remove_pointer_t<ARKRhsFn>> f2, sunrealtype t0,
       N_Vector y0, SUNContext sunctx)
    {
      if (!f1 && !f2)
      {
        throw sundials4py::illegal_value("f1 and f2 cannot be null");
      }

      void* ark_mem = SPRKStepCreate(sprkstep_f1_wrapper, sprkstep_f2_wrapper,
                                     t0, y0, sunctx);
      if (ark_mem == nullptr)
      {
        throw sundials4py::error_returned("Failed to create SPRKStep memory");
      }

      auto fn_table         = new arkode_user_supplied_fn_table;
      fn_table->sprkstep_f1 = nb::cast(f1);
      fn_table->sprkstep_f2 = nb::cast(f2);

      static_cast<ARKodeMem>(ark_mem)->python = fn_table;

      int ark_status = ARKodeSetUserData(ark_mem, ark_mem);
      if (ark_status != ARK_SUCCESS)
      {
        free(fn_table);
        throw sundials4py::error_returned(
          "Failed to set user data in SPRKStep memory");
      }

      return std::make_shared<ARKodeView>(ark_mem);
    },
    nb::arg("f1"), nb::arg("f2"), nb::arg("t0"), nb::arg("y0"),
    nb::arg("sunctx"), nb::keep_alive<0, 5>());
}

} // namespace sundials4py
