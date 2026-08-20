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

#include "arkode/arkode.h"
#include "sundials4py.hpp"

#include <arkode/arkode.hpp>
#include <arkode/arkode_mristep.hpp>
#include <sundials/sundials_core.hpp>

#include "arkode_usersupplied.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_arkode_mristep(nb::module_& m)
{
#include "arkode_mristep_generated.hpp"

  m.def(
    "MRIStepCreate",
    [](std::function<std::remove_pointer_t<ARKRhsFn>> fse,
       std::function<std::remove_pointer_t<ARKRhsFn>> fsi, sunrealtype t0,
       N_Vector y0, SUNStepper stepper, SUNContext sunctx)
    {
      auto fse_wrapper = fse ? mristep_fse_wrapper : nullptr;
      auto fsi_wrapper = fsi ? mristep_fsi_wrapper : nullptr;

      void* ark_mem = MRIStepCreate(fse_wrapper, fsi_wrapper, t0, y0, stepper,
                                    sunctx);
      if (ark_mem == nullptr)
      {
        throw sundials4py::error_returned("MRIStepCreate returned NULL");
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
      fn_table->mristep_fse = nb::cast(fse);
      fn_table->mristep_fsi = nb::cast(fsi);

      return std::make_shared<ARKodeView>(ark_mem);
    },
    nb::arg("fse").none(), nb::arg("fsi").none(), nb::arg("t0"), nb::arg("y0"),
    nb::arg("inner_stepper"), nb::arg("sunctx"), nb::keep_alive<0, 6>());
}

} // namespace sundials4py
