/*---------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
                  Cody J. Balos @ LLNL
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
 * This is the implementation file for ARKODE's interfacing with
 * SUNStepper
 *--------------------------------------------------------------*/

#include <arkode/arkode.h>
#include <sundials/sundials_stepper.h>
#include "arkode_impl.h"
#include "sundials/sundials_types.h"
#include "sundials_macros.h"

static SUNErrCode arkSUNStepperEvolveHelper(SUNStepper stepper,
                                            sunrealtype tout, N_Vector y,
                                            sunrealtype* tret, int mode)
{
  SUNFunctionBegin(stepper->sunctx);
  /* extract the ARKODE memory struct */
  void* arkode_mem;
  SUNCheckCall(SUNStepper_GetContent(stepper, &arkode_mem));

  /* evolve inner ODE */
  stepper->last_flag = ARKodeEvolve(arkode_mem, tout, y, tret, mode);
  if (stepper->last_flag < 0) { return SUN_ERR_OP_FAIL; }

  return SUN_SUCCESS;
}

static SUNErrCode arkSUNStepperEvolve(SUNStepper stepper, sunrealtype tout,
                                      N_Vector y, sunrealtype* tret)
{
  return arkSUNStepperEvolveHelper(stepper, tout, y, tret, ARK_NORMAL);
}

static SUNErrCode arkSUNStepperOneStep(SUNStepper stepper, sunrealtype tout,
                                       N_Vector y, sunrealtype* tret)
{
  return arkSUNStepperEvolveHelper(stepper, tout, y, tret, ARK_ONE_STEP);
}

/*------------------------------------------------------------------------------
  Implementation of SUNStepperFullRhsFn to compute the full inner
  (fast) ODE IVP RHS.
  ----------------------------------------------------------------------------*/

static SUNErrCode arkSUNStepperFullRhs(SUNStepper stepper, sunrealtype t,
                                       N_Vector y, N_Vector f,
                                       SUNFullRhsMode mode)
{
  SUNFunctionBegin(stepper->sunctx);
  /* extract the ARKODE memory struct */
  void* arkode_mem;
  SUNCheckCall(SUNStepper_GetContent(stepper, &arkode_mem));
  ARKodeMem ark_mem = (ARKodeMem)arkode_mem;

  int ark_mode;
  switch (mode)
  {
  case SUN_FULLRHS_START: ark_mode = ARK_FULLRHS_START; break;
  case SUN_FULLRHS_END: ark_mode = ARK_FULLRHS_END; break;
  case SUN_FULLRHS_OTHER: ark_mode = ARK_FULLRHS_OTHER; break;
  default: ark_mode = -1; break;
  }

  stepper->last_flag = ark_mem->step_fullrhs(ark_mem, t, y, f, ark_mode);
  if (stepper->last_flag != ARK_SUCCESS) { return SUN_ERR_OP_FAIL; }

  return SUN_SUCCESS;
}

/*------------------------------------------------------------------------------
  Implementation of SUNStepperResetFn to reset the stepper state.
  ----------------------------------------------------------------------------*/

static SUNErrCode arkSUNStepperReset(SUNStepper stepper, sunrealtype tR,
                                     N_Vector yR)
{
  SUNFunctionBegin(stepper->sunctx);

  /* extract the ARKODE memory struct */
  void* arkode_mem;
  SUNCheckCall(SUNStepper_GetContent(stepper, &arkode_mem));

  stepper->last_flag = ARKodeReset(arkode_mem, tR, yR);
  if (stepper->last_flag != ARK_SUCCESS) { return SUN_ERR_OP_FAIL; }

  return SUN_SUCCESS;
}

/*------------------------------------------------------------------------------
  Implementation of SUNStepperResetCheckpointIndexFn.
  ----------------------------------------------------------------------------*/

static SUNErrCode arkSUNStepperResetCheckpointIndex(SUNStepper stepper,
                                                    suncountertype ckptIdxR)
{
  SUNFunctionBegin(stepper->sunctx);

  /* extract the ARKODE memory struct */
  void* arkode_mem;
  SUNCheckCall(SUNStepper_GetContent(stepper, &arkode_mem));

  stepper->last_flag = ARKodeSetAdjointCheckpointIndex(arkode_mem, ckptIdxR);
  if (stepper->last_flag != ARK_SUCCESS) { return SUN_ERR_OP_FAIL; }

  return SUN_SUCCESS;
}

/*------------------------------------------------------------------------------
  Implementation of SUNStepperStopTimeFn to set the tstop time
  ----------------------------------------------------------------------------*/

static SUNErrCode arkSUNStepperSetStopTime(SUNStepper stepper, sunrealtype tstop)
{
  SUNFunctionBegin(stepper->sunctx);
  /* extract the ARKODE memory struct */
  void* arkode_mem;
  SUNCheckCall(SUNStepper_GetContent(stepper, &arkode_mem));

  stepper->last_flag = ARKodeSetStopTime(arkode_mem, tstop);
  if (stepper->last_flag != ARK_SUCCESS) { return SUN_ERR_OP_FAIL; }

  return SUN_SUCCESS;
}

static SUNErrCode arkSUNStepperSetStepDirection(SUNStepper stepper,
                                                sunrealtype stepdir)
{
  SUNFunctionBegin(stepper->sunctx);
  /* extract the ARKODE memory struct */
  void* arkode_mem;
  SUNCheckCall(SUNStepper_GetContent(stepper, &arkode_mem));

  stepper->last_flag = ARKodeSetStepDirection(arkode_mem, stepdir);
  if (stepper->last_flag != ARK_SUCCESS) { return SUN_ERR_OP_FAIL; }

  return SUN_SUCCESS;
}

static SUNErrCode arkSUNStepperSetForcing(SUNStepper stepper, sunrealtype tshift,
                                          sunrealtype tscale, N_Vector* forcing,
                                          int nforcing)
{
  SUNFunctionBegin(stepper->sunctx);
  /* extract the ARKODE memory struct */
  void* arkode_mem;
  SUNCheckCall(SUNStepper_GetContent(stepper, &arkode_mem));
  ARKodeMem ark_mem = (ARKodeMem)arkode_mem;

  stepper->last_flag = ark_mem->step_setforcing(ark_mem, tshift, tscale,
                                                forcing, nforcing);
  if (stepper->last_flag != ARK_SUCCESS) { return SUN_ERR_OP_FAIL; }

  return SUN_SUCCESS;
}

SUNErrCode arkSUNStepperSelfDestruct(SUNStepper stepper)
{
  /* This function is useful when we create a ARKodeMem/SUNStepper internally,
     and want it to be destroyed with the SUNStepper. */
  ARKodeMem ark_mem;

  SUNErrCode errcode = SUNStepper_GetContent(stepper, (void**)&ark_mem);
  if (errcode) { return errcode; }

  ARKodeFree((void**)&ark_mem);

  return SUN_SUCCESS;
}

static SUNErrCode arkSUNStepperGetNumSteps(SUNStepper stepper, suncountertype* nst)
{
  ARKodeMem ark_mem;

  SUNErrCode errcode = SUNStepper_GetContent(stepper, (void**)&ark_mem);
  if (errcode) { return errcode; }

  *nst = ark_mem->nst;

  return SUN_SUCCESS;
}

int ARKodeCreateSUNStepper(void* arkode_mem, SUNStepper* stepper)
{
  /* unpack ark_mem */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return ARK_MEM_NULL;
  }
  ARKodeMem ark_mem = (ARKodeMem)arkode_mem;

  SUNErrCode err = SUNStepper_Create(ark_mem->sunctx, stepper);
  if (err != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_SUNSTEPPER_ERR, __LINE__, __func__, __FILE__,
                    "Failed to create SUNStepper");
    return ARK_SUNSTEPPER_ERR;
  }

  err = SUNStepper_SetContent(*stepper, arkode_mem);
  if (err != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_SUNSTEPPER_ERR, __LINE__, __func__, __FILE__,
                    "Failed to set SUNStepper content");
    return ARK_SUNSTEPPER_ERR;
  }

  err = SUNStepper_SetEvolveFn(*stepper, arkSUNStepperEvolve);
  if (err != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_SUNSTEPPER_ERR, __LINE__, __func__, __FILE__,
                    "Failed to set SUNStepper evolve function");
    return ARK_SUNSTEPPER_ERR;
  }

  err = SUNStepper_SetOneStepFn(*stepper, arkSUNStepperOneStep);
  if (err != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_SUNSTEPPER_ERR, __LINE__, __func__, __FILE__,
                    "Failed to set SUNStepper one step function");
    return ARK_SUNSTEPPER_ERR;
  }

  err = SUNStepper_SetFullRhsFn(*stepper, arkSUNStepperFullRhs);
  if (err != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_SUNSTEPPER_ERR, __LINE__, __func__, __FILE__,
                    "Failed to set SUNStepper full RHS function");
    return ARK_SUNSTEPPER_ERR;
  }

  err = SUNStepper_SetResetFn(*stepper, arkSUNStepperReset);
  if (err != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_SUNSTEPPER_ERR, __LINE__, __func__, __FILE__,
                    "Failed to set SUNStepper reset function");
    return ARK_SUNSTEPPER_ERR;
  }

  err = SUNStepper_SetResetCheckpointIndexFn(*stepper,
                                             arkSUNStepperResetCheckpointIndex);
  if (err != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_SUNSTEPPER_ERR, __LINE__, __func__, __FILE__,
                    "Failed to set SUNStepper reset function");
    return ARK_SUNSTEPPER_ERR;
  }

  err = SUNStepper_SetStopTimeFn(*stepper, arkSUNStepperSetStopTime);
  if (err != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_SUNSTEPPER_ERR, __LINE__, __func__, __FILE__,
                    "Failed to set SUNStepper stop time function");
    return ARK_SUNSTEPPER_ERR;
  }

  err = SUNStepper_SetStepDirectionFn(*stepper, arkSUNStepperSetStepDirection);
  if (err != SUN_SUCCESS) { return ARK_SUNSTEPPER_ERR; }

  if (ark_mem->step_setforcing != NULL)
  {
    err = SUNStepper_SetForcingFn(*stepper, arkSUNStepperSetForcing);
    if (err != SUN_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_SUNSTEPPER_ERR, __LINE__, __func__, __FILE__,
                      "Failed to set SUNStepper forcing function");
      return ARK_SUNSTEPPER_ERR;
    }
  }

  err = SUNStepper_SetGetNumStepsFn(*stepper, arkSUNStepperGetNumSteps);
  if (err != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_SUNSTEPPER_ERR, __LINE__, __func__, __FILE__,
                    "Failed to set SUNStepper get number of steps function");
    return ARK_SUNSTEPPER_ERR;
  }

  return ARK_SUCCESS;
}
