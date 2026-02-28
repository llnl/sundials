/* -----------------------------------------------------------------
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
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_ADJOINTSTEPPER_IMPL_H_
#define SUNDIALS_ADJOINTSTEPPER_IMPL_H_

#include <sundials/sundials_adjointcheckpointscheme.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_stepper.h>

struct SUNAdjointStepper_
{
  suncountertype nrecompute;
  suncountertype final_step_idx;

  SUNStepper adj_sunstepper, fwd_sunstepper;
  sunbooleantype own_adj_sunstepper, own_fwd_sunstepper;
  SUNAdjointCheckpointScheme checkpoint_scheme;

  void* user_data;
  void* content;
  SUNContext sunctx;

  sunrealtype tf;
};

#endif
