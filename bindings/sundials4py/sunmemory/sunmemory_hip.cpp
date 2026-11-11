/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * All rights reserved.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <sundials/sundials_core.hpp>
#include <sunmemory/sunmemory_hip.h>

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_sunmemoryhelper_hip(nb::module_& m)
{
#include "sunmemory_hip_generated.hpp"
}

} // namespace sundials4py
