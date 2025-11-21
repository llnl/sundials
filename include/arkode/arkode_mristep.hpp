/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
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
 * -----------------------------------------------------------------------------
 * C++ specific ARKODE definitions.
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_ARKODE_MRISTEP_HPP
#define _SUNDIALS_ARKODE_MRISTEP_HPP

#include <arkode/arkode.hpp>
#include <arkode/arkode_mristep.h>

namespace sundials {
namespace experimental {

struct MRIStepCouplingDeleter
{
  void operator()(MRIStepCoupling t) { MRIStepCoupling_Free(t); }
};

struct MRIStepInnerStepperDeleter
{
  void operator()(MRIStepInnerStepper s) { MRIStepInnerStepper_Free(&s); }
};

} // namespace experimental
} // namespace sundials

#endif
