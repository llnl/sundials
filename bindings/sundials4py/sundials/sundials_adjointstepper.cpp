/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * This file is the entrypoint for the Python binding code for the
 * SUNDIALS SUNAdjointStepper class. It contains hand-written code for
 * functions that require special treatment, and includes the generated
 * code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <sundials/sundials_adjointstepper.hpp>
#include <sundials/sundials_classview.hpp>

#include "sundials_adjointcheckpointscheme_impl.h"
#include "sundials_adjointstepper_impl.h"
#include "sundials_stepper_impl.h"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_sunadjointstepper(nb::module_& m)
{
#include "sundials_adjointstepper_generated.hpp"

  nb::class_<SUNAdjointStepper_>(m, "SUNAdjointStepper_");
}

} // namespace sundials4py
