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
 * This file contains hand-written bindings for sundials_errors.h
 * -----------------------------------------------------------------*/

#include <sundials/sundials_errors.h>

#include "sundials_errors_generated.hpp"

/* Expand SUN_ERR_CODE_LIST to enum */
#define SUN_EXPAND_TO_NB_BINDING(name, description) \
  .value(#name, name, description)

auto pyEnumSUNErrCode_ = nb::enum_<SUNErrCode_>(m, "SUNErrCode",
                                                nb::is_arithmetic(), "")
                           .value("SUN_ERR_MINIMUM", SUN_ERR_MINIMUM, "")
                             SUN_ERR_CODE_LIST(SUN_EXPAND_TO_NB_BINDING)
                           .value("SUN_ERR_MAXIMUM", SUN_ERR_MAXIMUM, "")
                           .value("SUN_SUCCESS", SUN_SUCCESS, "")
                           .export_values();
