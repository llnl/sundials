/*---------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 *---------------------------------------------------------------
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
 *---------------------------------------------------------------
 * This header defines the step memory for ForcingStep.
 *--------------------------------------------------------------*/

#ifndef ARKODE_FORCINGSTEP_IMPL_H_
#define ARKODE_FORCINGSTEP_IMPL_H_

#include <sundials/sundials_stepper.h>

#define NUM_PARTITIONS 2

typedef struct ARKodeForcingStepMemRec
{
  SUNStepper stepper[NUM_PARTITIONS];
  long int n_stepper_evolves[NUM_PARTITIONS];
}* ARKodeForcingStepMem;

#endif
