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
#include <arkode/arkode_arkstep.hpp>
#include <sundials/sundials_core.hpp>

#include "arkode_impl.h"
#include "arkode_usersupplied.hpp"
#include "sundials_adjointstepper_impl.h"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_arkode_arkstep(nb::module_& m)
{
#include "arkode_arkstep_generated.hpp"

  m.def(
    "ARKStepCreate",
    [](std::function<std::remove_pointer_t<ARKRhsFn>> fe,
       std::function<std::remove_pointer_t<ARKRhsFn>> fi, sunrealtype t0,
       N_Vector y0, SUNContext sunctx)
    {
      auto fe_wrapper = fe ? arkstep_fe_wrapper : nullptr;
      auto fi_wrapper = fi ? arkstep_fi_wrapper : nullptr;

      void* ark_mem = ARKStepCreate(fe_wrapper, fi_wrapper, t0, y0, sunctx);
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

      // Finally, set the RHS functions
      if (fe) { fn_table->arkstep_fe = nb::cast(fe); }
      if (fi) { fn_table->arkstep_fi = nb::cast(fi); }

      return std::make_shared<ARKodeView>(ark_mem);
    }, // .none() must be added to functions that accept nullptr as a valid argument
    nb::arg("fe").none(), nb::arg("fi").none(), nb::arg("t0"), nb::arg("y0"),
    nb::arg("sunctx"), nb::keep_alive<0, 5>());

  m.def(
    "ARKStepCreateAdjointStepper",
    [](void* arkode_mem, std::function<std::remove_pointer_t<SUNAdjRhsFn>> adj_fe,
       std::function<std::remove_pointer_t<SUNAdjRhsFn>> adj_fi, sunrealtype tf,
       N_Vector sf, SUNContext sunctx) -> std::tuple<int, SUNAdjointStepper>
    {
      auto fe_wrapper = adj_fe ? arkstep_adjfe_wrapper : nullptr;
      auto fi_wrapper = adj_fi ? arkstep_adjfi_wrapper : nullptr;

      SUNAdjointStepper adj_stepper = nullptr;
      int ark_status = ARKStepCreateAdjointStepper(arkode_mem, fe_wrapper,
                                                   fi_wrapper, tf, sf, sunctx,
                                                   &adj_stepper);
      if (ark_status != ARK_SUCCESS)
      {
        throw sundials4py::error_returned(
          "Failed to create adjoint stepper in py-sundials memory");
      }

      // Finally, set the RHS functions
      void* user_data = nullptr;
      ark_status      = ARKodeGetUserData(arkode_mem, &user_data);
      if (ark_status != ARK_SUCCESS)
      {
        throw sundials4py::error_returned("Failed to extract ARKODE user data");
      }

      auto fn_table = get_arkode_fn_table(arkode_mem);

      if (adj_fe) { fn_table->arkstep_adjfe = nb::cast(adj_fe); }
      if (adj_fi) { fn_table->arkstep_adjfi = nb::cast(adj_fi); }

      return std::make_tuple(ark_status, adj_stepper);
    },
    nb::arg("arkode_mem"), nb::arg("adj_fe").none(), nb::arg("adj_fi").none(),
    nb::arg("tf"), nb::arg("sf"), nb::arg("sunctx"));
}

} // namespace sundials4py
