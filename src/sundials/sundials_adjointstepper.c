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

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_adjointcheckpointscheme.h>
#include <sundials/sundials_adjointstepper.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_stepper.h>
#include "sundials/sundials_types.h"
#include "sundials_adjointstepper_impl.h"
#include "sundials_macros.h"
#include "sundials_stepper_impl.h"
#include "sundials_utils.h"

SUNErrCode SUNAdjointStepper_Create(
  SUNStepper fwd_sunstepper, sunbooleantype own_fwd, SUNStepper adj_sunstepper,
  sunbooleantype own_adj, suncountertype final_step_idx, sunrealtype tf,
  SUNDIALS_MAYBE_UNUSED N_Vector sf, SUNAdjointCheckpointScheme checkpoint_scheme,
  SUNContext sunctx, SUNAdjointStepper* adj_stepper_ptr)
{
  SUNFunctionBegin(sunctx);

  SUNAdjointStepper adj_stepper = NULL;
  adj_stepper = (SUNAdjointStepper)malloc(sizeof(*adj_stepper));
  SUNAssert(adj_stepper, SUN_ERR_MALLOC_FAIL);

  adj_stepper->fwd_sunstepper     = fwd_sunstepper;
  adj_stepper->own_fwd_sunstepper = own_fwd;
  adj_stepper->adj_sunstepper     = adj_sunstepper;
  adj_stepper->own_adj_sunstepper = own_adj;
  adj_stepper->checkpoint_scheme  = checkpoint_scheme;

  adj_stepper->tf             = tf;
  adj_stepper->final_step_idx = final_step_idx;

  adj_stepper->nrecompute = 0;

  adj_stepper->user_data = NULL;
  adj_stepper->sunctx    = sunctx;

  *adj_stepper_ptr = adj_stepper;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_ReInit(SUNAdjointStepper self, sunrealtype t0,
                                    N_Vector y0, sunrealtype tf, N_Vector sf)
{
  SUNFunctionBegin(self->sunctx);
  self->tf         = tf;
  self->nrecompute = 0;
  SUNStepper_ReInit(self->adj_sunstepper, tf, sf);
  SUNStepper_ReInit(self->fwd_sunstepper, t0, y0);
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_Evolve(SUNAdjointStepper self, sunrealtype tout,
                                    N_Vector sens, sunrealtype* tret)

{
  SUNFunctionBegin(self->sunctx);
  SUNCheckCall(SUNStepper_Evolve(self->adj_sunstepper, tout, sens, tret));
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_OneStep(SUNAdjointStepper self, sunrealtype tout,
                                     N_Vector sens, sunrealtype* tret)

{
  SUNFunctionBegin(self->sunctx);
  SUNCheckCall(SUNStepper_OneStep(self->adj_sunstepper, tout, sens, tret));
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_RecomputeFwd(SUNAdjointStepper self,
                                          suncountertype start_idx,
                                          sunrealtype t0, N_Vector y0,
                                          sunrealtype tf)
{
  SUNFunctionBegin(self->sunctx);

  SUNErrCode retcode = SUN_SUCCESS;

  sunrealtype fwd_t      = t0;
  SUNStepper fwd_stepper = self->fwd_sunstepper;
  SUNCheckCall(SUNStepper_Reset(fwd_stepper, t0, y0));
  SUNCheckCall(SUNStepper_ResetCheckpointIndex(fwd_stepper, start_idx));

  SUNCheckCall(SUNAdjointCheckpointScheme_EnableDense(self->checkpoint_scheme, 1));

  SUNCheckCall(SUNStepper_SetStopTime(fwd_stepper, tf));

  suncountertype nst_before, nst_after;
  SUNCheckCall(SUNStepper_GetNumSteps(fwd_stepper, &nst_before));
  SUNCheckCall(SUNStepper_Evolve(fwd_stepper, tf, y0, &fwd_t));
  SUNCheckCall(SUNStepper_GetNumSteps(fwd_stepper, &nst_after));
  self->nrecompute += nst_after - nst_before;

  SUNCheckCall(SUNAdjointCheckpointScheme_EnableDense(self->checkpoint_scheme, 0));

  return retcode;
}

SUNErrCode SUNAdjointStepper_Destroy(SUNAdjointStepper* self_ptr)
{
  SUNAdjointStepper self = *self_ptr;
  if (self->own_fwd_sunstepper) { SUNStepper_Destroy(&self->fwd_sunstepper); }
  if (self->own_adj_sunstepper) { SUNStepper_Destroy(&self->adj_sunstepper); }
  free(self);
  *self_ptr = NULL;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_SetUserData(SUNAdjointStepper self, void* user_data)
{
  SUNFunctionBegin(self->sunctx);

  self->user_data = user_data;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_GetNumSteps(SUNAdjointStepper self,
                                         suncountertype* num_steps)
{
  SUNFunctionBegin(self->sunctx);
  return SUNStepper_GetNumSteps(self->adj_sunstepper, num_steps);
}

SUNErrCode SUNAdjointStepper_GetNumRecompute(SUNAdjointStepper self,
                                             suncountertype* num_recompute)
{
  SUNFunctionBegin(self->sunctx);
  *num_recompute = self->nrecompute;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointStepper_PrintAllStats(SUNAdjointStepper self,
                                           FILE* outfile, SUNOutputFormat fmt)
{
  SUNFunctionBegin(self->sunctx);
  suncountertype nst = 0;
  SUNCheckCall(SUNStepper_GetNumSteps(self->adj_sunstepper, &nst));
  sunfprintf_long(outfile, fmt, SUNTRUE, "Num backwards steps", nst);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Num recompute steps",
                  self->nrecompute);

  return SUN_SUCCESS;
}
