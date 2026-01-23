/*------------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *------------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2025, Lawrence Livermore National Security,
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

#include "arkode_mristep_impl.h"
#include "arkode_usersupplied.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_arkode_mristep(nb::module_& m)
{
#include "arkode_mristep_generated.hpp"

  // _MRIStepInnerStepper is a opaque/private class forward declared in a public header but
  // defined in a source file elsewhere. As such, we need to declare it here since its
  // not picked up in any header files by the generator.
  nb::class_<_MRIStepInnerStepper>(m, "_MRIStepInnerStepper");

  m.def(
    "MRIStepInnerStepper_Create",
    [](SUNContext sunctx)
    {
      MRIStepInnerStepper stepper = nullptr;

      int status = MRIStepInnerStepper_Create(sunctx, &stepper);
      if (status == ARK_SUCCESS)
      {
        auto fn_table   = new mristepinnerstepper_user_supplied_fn_table;
        stepper->python = static_cast<void*>(fn_table);
      }

      return std::make_tuple(status,
                             our_make_shared<std::remove_pointer_t<MRIStepInnerStepper>,
                                             MRIStepInnerStepperDeleter>(stepper));
    },
    nb::arg("sunctx"));

  m.def(
    "MRIStepInnerStepper_CreateFromSUNStepper",
    [](SUNStepper stepper)
    {
      MRIStepInnerStepper inner_stepper = nullptr;

      int status = MRIStepInnerStepper_CreateFromSUNStepper(stepper,
                                                            &inner_stepper);

      return std::make_tuple(status,
                             our_make_shared<std::remove_pointer_t<MRIStepInnerStepper>,
                                             MRIStepInnerStepperDeleter>(
                               inner_stepper));
    },
    nb::arg("stepper"));

  m.def(
    "MRIStepInnerStepper_GetForcingData",
    [](MRIStepInnerStepper stepper)
      -> std::tuple<int, sunrealtype, sunrealtype, std::vector<N_Vector>, int>
    {
      sunrealtype tshift   = 0.0;
      sunrealtype tscale   = 0.0;
      N_Vector* forcing_1d = nullptr;
      int nforcing         = 0;

      int status = MRIStepInnerStepper_GetForcingData(stepper, &tshift, &tscale,
                                                      &forcing_1d, &nforcing);
      std::vector<N_Vector> forcing(nforcing);
      if (status == ARK_SUCCESS)
      {
        for (int i = 0; i < nforcing; i++) { forcing[i] = forcing_1d[i]; }
      }
      return std::make_tuple(status, tshift, tscale, forcing, nforcing);
    },
    nb::rv_policy::reference);

  m.def(
    "ARKodeCreateMRIStepInnerStepper",
    [](void* inner_arkode_mem)
    {
      MRIStepInnerStepper stepper = nullptr;

      int status = ARKodeCreateMRIStepInnerStepper(inner_arkode_mem, &stepper);

      return std::make_tuple(status,
                             our_make_shared<std::remove_pointer_t<MRIStepInnerStepper>,
                                             MRIStepInnerStepperDeleter>(stepper));
    },
    nb::arg("inner_arkode_mem"));

  m.def(
    "MRIStepCreate",
    [](std::function<std::remove_pointer_t<ARKRhsFn>> fse,
       std::function<std::remove_pointer_t<ARKRhsFn>> fsi, sunrealtype t0,
       N_Vector y0, MRIStepInnerStepper stepper, SUNContext sunctx)
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

// The destroy functions gets called in our C code
extern "C" void mristepinnerstepper_user_supplied_fn_table_destroy(void* ptr)
{
  delete static_cast<mristepinnerstepper_user_supplied_fn_table*>(ptr);
}
