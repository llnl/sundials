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

#ifndef _SUNDIALS_STEPPER_IMPL_H
#define _SUNDIALS_STEPPER_IMPL_H

#include <sundials/sundials_core.h>
#include <sundials/sundials_stepper.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SUNStepper_Ops_* SUNStepper_Ops;

struct SUNStepper_Ops_
{
  SUNStepperEvolveFn evolve;
  SUNStepperOneStepFn onestep;
  SUNStepperFullRhsFn fullrhs;
  SUNStepperReInitFn reinit;
  SUNStepperResetFn reset;
  SUNStepperResetCheckpointIndexFn resetcheckpointindex;
  SUNStepperSetStopTimeFn setstoptime;
  SUNStepperSetStepDirectionFn setstepdirection;
  SUNStepperSetForcingFn setforcing;
  SUNStepperGetNumStepsFn getnumsteps;
  SUNStepperDestroyFn destroy;
};

struct SUNStepper_
{
  /* python interface specific content */
  void* python;

  /* stepper specific content and operations */
  void* content;
  SUNStepper_Ops ops;

  /* stepper context */
  SUNContext sunctx;

  /* last stepper return flag */
  int last_flag;
};

#ifdef __cplusplus
}
#endif

#endif /* _SUNDIALS_STEPPER_IMPL_H */
