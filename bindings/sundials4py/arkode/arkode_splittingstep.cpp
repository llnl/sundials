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

#include "sundials4py.hpp"

#include <arkode/arkode.hpp>
#include <arkode/arkode_splittingstep.hpp>
#include <sundials/sundials_core.hpp>

#include "sundials_stepper_impl.h"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_arkode_splittingstep(nb::module_& m)
{
#include "arkode_splittingstep_generated.hpp"

  m.def(
    "SplittingStepCreate",
    [](std::vector<SUNStepper> steppers, int partitions, sunrealtype t0,
       N_Vector y0, SUNContext sunctx)
    {
      auto stepper = SplittingStepCreate(steppers.data(), partitions, t0, y0,
                                         sunctx);
      if (!stepper)
      {
        throw sundials4py::error_returned("SplittingStepCreate returned NULL");
      }
      return std::make_shared<ARKodeView>(stepper);
    },
    nb::arg("steppers"), nb::arg("partitions"), nb::arg("t0"), nb::arg("y0"),
    nb::arg("sunctx"), nb::keep_alive<0, 5>());

  m.def(
    "SplittingStepReInit",
    [](void* arkode_mem, std::vector<SUNStepper> steppers, int partitions,
       sunrealtype t0, N_Vector y0) -> int
    {
      return SplittingStepReInit(arkode_mem, steppers.data(), partitions, t0, y0);
    },
    nb::arg("arkode_mem"), nb::arg("steppers"), nb::arg("partitions"),
    nb::arg("t0"), nb::arg("y0"));
}

} // namespace sundials4py
