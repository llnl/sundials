/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_pcg.h>

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_sunlinsol_pcg(nb::module_& m)
{
#include "sunlinsol_pcg_generated.hpp"
}

} // namespace sundials4py
