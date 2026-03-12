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
#include <arkode/arkode_erkstep.h>
#include <sundials/sundials_core.hpp>

#include "arkode_impl.h"
#include "arkode_usersupplied.hpp"
#include "sundials_adjointstepper_impl.h"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_arkode_erkstep(nb::module_& m)
{
#include "arkode_erkstep_generated.hpp"

  m.def(
    "ERKStepCreate",
    [](std::function<std::remove_pointer_t<ARKRhsFn>> rhs, sunrealtype t0,
       N_Vector y0, SUNContext sunctx)
    {
      if (!rhs) { throw sundials4py::illegal_value("rhs was None"); }

      void* ark_mem = ERKStepCreate(erkstep_f_wrapper, t0, y0, sunctx);
      if (ark_mem == nullptr)
      {
        throw sundials4py::error_returned("Failed to create ARKODE memory");
      }

      // Create the user-supplied function table to store the Python user functions
      auto fn_table = new arkode_user_supplied_fn_table;

      // Smuggle the user-supplied function table into callback wrappers through the user_data pointer
      static_cast<ARKodeMem>(ark_mem)->python = fn_table;
      int ark_status = ARKodeSetUserData(ark_mem, ark_mem);
      if (ark_status != ARK_SUCCESS)
      {
        free(fn_table);
        throw sundials4py::error_returned(
          "Failed to set user data in ARKODE memory");
      }

      // Finally, set the RHS function
      fn_table->erkstep_f = nb::cast(rhs);

      return std::make_shared<ARKodeView>(ark_mem);
    },
    nb::arg("rhs"), nb::arg("t0"), nb::arg("y0"), nb::arg("sunctx"),
    nb::keep_alive<0, 4>());

  m.def(
    "ERKStepCreateAdjointStepper",
    [](void* arkode_mem,
       std::function<std::remove_pointer_t<SUNAdjRhsFn>> adj_f, sunrealtype tf,
       N_Vector sf, SUNContext sunctx) -> std::tuple<int, SUNAdjointStepper>
    {
      if (!adj_f) { throw sundials4py::illegal_value("adj_f was null"); }

      SUNAdjointStepper adj_stepper = nullptr;
      int ark_status                = ERKStepCreateAdjointStepper(arkode_mem,
                                                                  erkstep_adjf_wrapper, tf, sf,
                                                                  sunctx, &adj_stepper);
      if (ark_status != ARK_SUCCESS)
      {
        throw sundials4py::error_returned("Failed to create adjoint stepper");
      }

      auto fn_table = get_arkode_fn_table(arkode_mem);

      fn_table->erkstep_adjf = nb::cast(adj_f);

      return std::make_tuple(ark_status, adj_stepper);
    },
    nb::arg("arkode_mem"), nb::arg("adj_f").none(), nb::arg("tf"),
    nb::arg("sf"), nb::arg("sunctx"));
}

} // namespace sundials4py
