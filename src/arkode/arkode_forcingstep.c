/*------------------------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 *------------------------------------------------------------------------------
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
 *------------------------------------------------------------------------------
 * This is the implementation file for ARKODE's forcing method
 *----------------------------------------------------------------------------*/

#include <arkode/arkode_forcingstep.h>
#include <sundials/sundials_nvector.h>

#include "arkode_forcingstep_impl.h"
#include "arkode_impl.h"
#include "sundials_utils.h"

/*------------------------------------------------------------------------------
  Shortcut routine to unpack step_mem structure from ark_mem. If missing it
  returns ARK_MEM_NULL.
  ----------------------------------------------------------------------------*/
static int forcingStep_AccessStepMem(ARKodeMem ark_mem, const char* fname,
                                     ARKodeForcingStepMem* step_mem)
{
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    "Time step module memory is NULL.");
    return ARK_MEM_NULL;
  }
  *step_mem = (ARKodeForcingStepMem)ark_mem->step_mem;
  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  Shortcut routine to unpack ark_mem and step_mem structures from void* pointer.
  If either is missing it returns ARK_MEM_NULL.
  ----------------------------------------------------------------------------*/
static int forcingStep_AccessARKODEStepMem(void* arkode_mem, const char* fname,
                                           ARKodeMem* ark_mem,
                                           ARKodeForcingStepMem* step_mem)
{
  /* access ARKodeMem structure */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_ARK_NO_MEM);
    return ARK_MEM_NULL;
  }
  *ark_mem = (ARKodeMem)arkode_mem;

  return forcingStep_AccessStepMem(*ark_mem, __func__, step_mem);
}

/*------------------------------------------------------------------------------
  This routine is called just prior to performing internal time steps (after
  all user "set" routines have been called) from within arkInitialSetup.
  ----------------------------------------------------------------------------*/
static int forcingStep_Init(ARKodeMem ark_mem,
                            SUNDIALS_MAYBE_UNUSED sunrealtype tout, int init_type)
{
  ARKodeForcingStepMem step_mem = NULL;
  int retval = forcingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  /* assume fixed outer step size */
  if (!ark_mem->fixedstep)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Adaptive outer time stepping is not currently supported");
    return ARK_ILL_INPUT;
  }

  if (ark_mem->interp_type == ARK_INTERP_HERMITE &&
      (step_mem->stepper[0]->ops->fullrhs == NULL ||
       step_mem->stepper[1]->ops->fullrhs == NULL))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The SUNSteppers must implement SUNStepper_FullRhs when "
                    "using Hermite interpolation");
    return ARK_ILL_INPUT;
  }

  /* immediately return if resize or reset */
  if (init_type == RESIZE_INIT || init_type == RESET_INIT)
  {
    return ARK_SUCCESS;
  }

  /* On first initialization, make the SUNStepper consistent with the current
   * state in case a user provided a different initial condition for the
   * ForcingStep integrator and SUNStepper. */
  SUNErrCode err = SUNStepper_Reset(step_mem->stepper[1], ark_mem->tn,
                                    ark_mem->yn);
  if (err != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_SUNSTEPPER_ERR, __LINE__, __func__, __FILE__,
                    "Resetting the second partition SUNStepper failed");
    return ARK_SUNSTEPPER_ERR;
  }

  ark_mem->interp_degree = 1;

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  This routine resets the ForcingStep integrator by resetting the partition
  integrators
  ----------------------------------------------------------------------------*/
static int forcingStep_Reset(ARKodeMem ark_mem, sunrealtype tR, N_Vector yR)
{
  ARKodeForcingStepMem step_mem = NULL;
  int retval = forcingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  SUNErrCode err = SUNStepper_Reset(step_mem->stepper[0], tR, yR);
  if (err != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_SUNSTEPPER_ERR, __LINE__, __func__, __FILE__,
                    "Resetting the first partition SUNStepper failed");
    return ARK_SUNSTEPPER_ERR;
  }

  err = SUNStepper_Reset(step_mem->stepper[1], tR, yR);
  if (err != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_SUNSTEPPER_ERR, __LINE__, __func__, __FILE__,
                    "Resetting the second partition SUNStepper failed");
    return ARK_SUNSTEPPER_ERR;
  }

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  This routine sets the step direction of the partition integrators and is
  called once the ForcingStep integrator has updated its step direction.
  ----------------------------------------------------------------------------*/
static int forcingStep_SetStepDirection(ARKodeMem ark_mem, sunrealtype stepdir)
{
  ARKodeForcingStepMem step_mem = NULL;
  int retval = forcingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  SUNErrCode err = SUNStepper_SetStepDirection(step_mem->stepper[0], stepdir);
  if (err != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_SUNSTEPPER_ERR, __LINE__, __func__,
                    __FILE__, "Setting the step direction for the first partition SUNStepper failed");
    return ARK_SUNSTEPPER_ERR;
  }

  err = SUNStepper_SetStepDirection(step_mem->stepper[1], stepdir);
  if (err != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_SUNSTEPPER_ERR, __LINE__, __func__,
                    __FILE__, "Setting the step direction for the second partition SUNStepper failed");
    return ARK_SUNSTEPPER_ERR;
  }

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  This is just a wrapper to call the user-supplied RHS function,
  f^1(t,y) + f^2(t,y).

  This will be called in one of three 'modes':

     ARK_FULLRHS_START -> called at the beginning of a simulation i.e., at
                          (tn, yn) = (t0, y0) or (tR, yR)

     ARK_FULLRHS_END   -> called at the end of a successful step i.e, at
                          (tcur, ycur) or the start of the subsequent step i.e.,
                          at (tn, yn) = (tcur, ycur) from the end of the last
                          step

     ARK_FULLRHS_OTHER -> called elsewhere (e.g. for dense output)

  The stepper for partition 1 has a state that is inconsistent with the
  ForcingStep integrator, so we cannot pass it the SUN_FULLRHS_END option. For
  partition 2, the state should be consistent, and we can use SUN_FULLRHS_END.
  ----------------------------------------------------------------------------*/
static int forcingStep_FullRHS(ARKodeMem ark_mem, sunrealtype t, N_Vector y,
                               N_Vector f, int mode)
{
  ARKodeForcingStepMem step_mem = NULL;
  int retval = forcingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  /* TODO(SBR): Possible optimization in FULLRHS_START mode. Currently that
   * mode is not forwarded to the SUNSteppers */
  SUNErrCode err = SUNStepper_FullRhs(step_mem->stepper[0], t, y,
                                      ark_mem->tempv1, SUN_FULLRHS_OTHER);
  if (err != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_RHSFUNC_FAILED, t);
    return ARK_RHSFUNC_FAIL;
  }

  err = SUNStepper_FullRhs(step_mem->stepper[1], t, y, f,
                           mode == ARK_FULLRHS_END ? SUN_FULLRHS_END
                                                   : SUN_FULLRHS_OTHER);
  if (err != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_RHSFUNC_FAILED, t);
    return ARK_RHSFUNC_FAIL;
  }
  N_VLinearSum(SUN_RCONST(1.0), f, SUN_RCONST(1.0), ark_mem->tempv1, f);

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  This routine performs a single step of the forcing method.
  ----------------------------------------------------------------------------*/
static int forcingStep_TakeStep(ARKodeMem ark_mem, sunrealtype* dsmPtr,
                                int* nflagPtr)
{
  ARKodeForcingStepMem step_mem = NULL;
  int retval = forcingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  *nflagPtr = ARK_SUCCESS; /* No algebraic solver */
  *dsmPtr   = ZERO;        /* No error estimate */

  SUNStepper s0    = step_mem->stepper[0];
  sunrealtype tout = ark_mem->tn + ark_mem->h;
  sunrealtype tret = ZERO;

  /* Evolve stepper 0 on its own */
  SUNLogInfo(ARK_LOGGER, "begin-partitions-list", "partition = 0");

  SUNErrCode err = SUNStepper_Reset(s0, ark_mem->tn, ark_mem->yn);
  if (err != SUN_SUCCESS)
  {
    SUNLogInfo(ARK_LOGGER, "end-partitions-list",
               "status = failed stepper reset, err = %i", err);
    return ARK_SUNSTEPPER_ERR;
  }

  err = SUNStepper_SetStopTime(s0, tout);
  if (err != SUN_SUCCESS)
  {
    SUNLogInfo(ARK_LOGGER, "end-partitions-list",
               "status = failed set stop time, err = %i", err);
    return ARK_SUNSTEPPER_ERR;
  }

  err = SUNStepper_Evolve(s0, tout, ark_mem->ycur, &tret);
  SUNLogExtraDebugVec(ARK_LOGGER, "partition state", ark_mem->ycur, "y_par(:) =");
  if (err != SUN_SUCCESS)
  {
    SUNLogInfo(ARK_LOGGER, "end-partitions-list",
               "status = failed evolve, err = %i", err);
    return ARK_SUNSTEPPER_ERR;
  }
  step_mem->n_stepper_evolves[0]++;

  SUNLogInfo(ARK_LOGGER, "end-partitions-list", "status = success");
  SUNLogInfo(ARK_LOGGER, "begin-partitions-list", "partition = 1");

  SUNStepper s1 = step_mem->stepper[1];
  /* A reset is not needed because steeper 1's state is consistent with the
   * forcing method */
  err = SUNStepper_SetStopTime(s1, tout);
  if (err != SUN_SUCCESS)
  {
    SUNLogInfo(ARK_LOGGER, "end-partitions-list",
               "status = failed set stop time, err = %i", err);
    return ARK_SUNSTEPPER_ERR;
  }

  /* Write tendency (ycur - yn)/h into stepper 1 forcing */
  sunrealtype hinv = SUN_RCONST(1.0) / ark_mem->h;
  N_VLinearSum(hinv, ark_mem->ycur, -hinv, ark_mem->yn, ark_mem->tempv1);
  err = SUNStepper_SetForcing(s1, ZERO, ZERO, &ark_mem->tempv1, 1);
  SUNLogExtraDebugVec(ARK_LOGGER, "forcing", ark_mem->tempv1, "forcing(:) =");
  if (err != SUN_SUCCESS)
  {
    SUNLogInfo(ARK_LOGGER, "end-partitions-list",
               "status = failed set forcing, err = %i", err);
    return ARK_SUNSTEPPER_ERR;
  }

  /* Evolve stepper 1 with the forcing */
  err = SUNStepper_Evolve(s1, tout, ark_mem->ycur, &tret);
  SUNLogExtraDebugVec(ARK_LOGGER, "partition state", ark_mem->ycur, "y_par(:) =");
  if (err != SUN_SUCCESS)
  {
    SUNLogInfo(ARK_LOGGER, "end-partitions-list",
               "status = failed evolve, err = %i", err);
    return ARK_SUNSTEPPER_ERR;
  }
  step_mem->n_stepper_evolves[1]++;

  /* Clear the forcing so it doesn't get included in a fullRhs call */
  err = SUNStepper_SetForcing(s1, ZERO, ZERO, NULL, 0);
  if (err != SUN_SUCCESS)
  {
    SUNLogInfo(ARK_LOGGER, "end-partitions-list",
               "status = failed set forcing, err = %i", err);
    return ARK_SUNSTEPPER_ERR;
  }

  SUNLogInfo(ARK_LOGGER, "end-partitions-list", "status = success");
  SUNLogExtraDebugVec(ARK_LOGGER, "current state", ark_mem->ycur, "y_cur(:) =");

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  Prints integrator statistics
  ----------------------------------------------------------------------------*/
static int forcingStep_PrintAllStats(ARKodeMem ark_mem, FILE* outfile,
                                     SUNOutputFormat fmt)
{
  ARKodeForcingStepMem step_mem = NULL;
  int retval = forcingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  sunfprintf_long(outfile, fmt, SUNFALSE, "Partition 1 evolves",
                  step_mem->n_stepper_evolves[0]);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Partition 2 evolves",
                  step_mem->n_stepper_evolves[1]);

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  Frees all ForcingStep memory.
  ----------------------------------------------------------------------------*/
static void forcingStep_Free(ARKodeMem ark_mem)
{
  if (ark_mem->step_mem != NULL) { free(ark_mem->step_mem); }
  ark_mem->step_mem = NULL;
}

/*------------------------------------------------------------------------------
  This routine outputs the memory from the ForcingStep structure to a specified
  file pointer (useful when debugging).
  ----------------------------------------------------------------------------*/
static void forcingStep_PrintMem(ARKodeMem ark_mem, FILE* outfile)
{
  ARKodeForcingStepMem step_mem = NULL;
  int retval = forcingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return; }

  /* output long integer quantities */
  for (int k = 0; k < NUM_PARTITIONS; k++)
  {
    fprintf(outfile, "ForcingStep: partition %i: n_stepper_evolves = %li\n", k,
            step_mem->n_stepper_evolves[k]);
  }
}

/*------------------------------------------------------------------------------
  This routine checks if all required SUNStepper operations are present. If any
  of them are missing it return SUNFALSE.
  ----------------------------------------------------------------------------*/
static sunbooleantype forcingStep_CheckSUNStepper(SUNStepper stepper,
                                                  sunbooleantype needs_forcing)
{
  SUNStepper_Ops ops = stepper->ops;
  return ops->evolve != NULL && ops->reset != NULL &&
         ops->setstoptime != NULL && (!needs_forcing || ops->setforcing != NULL);
}

/*------------------------------------------------------------------------------
  This routine validates arguments when (re)initializing a ForcingStep
  integrator
  ----------------------------------------------------------------------------*/
static int forcingStep_CheckArgs(ARKodeMem ark_mem, SUNStepper stepper1,
                                 SUNStepper stepper2, N_Vector y0)
{
  if (stepper1 == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "stepper1 = NULL illegal.");
    return ARK_ILL_INPUT;
  }
  if (!forcingStep_CheckSUNStepper(stepper1, SUNFALSE))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "stepper1 does not implement the required operations.");
    return ARK_ILL_INPUT;
  }

  if (stepper2 == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "stepper2 = NULL illegal.");
    return ARK_ILL_INPUT;
  }
  if (!forcingStep_CheckSUNStepper(stepper2, SUNTRUE))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "stepper2 does not implement the required operations.");
    return ARK_ILL_INPUT;
  }

  if (y0 == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_Y0);
    return ARK_ILL_INPUT;
  }

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  This routine initializes the step memory and resets the statistics
  ----------------------------------------------------------------------------*/
static void forcingStep_InitStepMem(ARKodeForcingStepMem step_mem,
                                    SUNStepper stepper1, SUNStepper stepper2)
{
  step_mem->stepper[0]           = stepper1;
  step_mem->stepper[1]           = stepper2;
  step_mem->n_stepper_evolves[0] = 0;
  step_mem->n_stepper_evolves[1] = 0;
}

/*------------------------------------------------------------------------------
  Creates the ForcingStep integrator
  ----------------------------------------------------------------------------*/
void* ForcingStepCreate(SUNStepper stepper1, SUNStepper stepper2,
                        sunrealtype t0, N_Vector y0, SUNContext sunctx)
{
  int retval = forcingStep_CheckArgs(NULL, stepper1, stepper2, y0);
  if (retval != ARK_SUCCESS) { return NULL; }

  if (sunctx == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_SUNCTX);
    return NULL;
  }

  /* Create ark_mem structure and set default values */
  ARKodeMem ark_mem = arkCreate(sunctx);
  if (ark_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return NULL;
  }

  ARKodeForcingStepMem step_mem = (ARKodeForcingStepMem)malloc(sizeof(*step_mem));
  if (step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }
  forcingStep_InitStepMem(step_mem, stepper1, stepper2);

  /* Attach step_mem structure and function pointers to ark_mem */
  ark_mem->step_init             = forcingStep_Init;
  ark_mem->step_fullrhs          = forcingStep_FullRHS;
  ark_mem->step_reset            = forcingStep_Reset;
  ark_mem->step_setstepdirection = forcingStep_SetStepDirection;
  ark_mem->step                  = forcingStep_TakeStep;
  ark_mem->step_printallstats    = forcingStep_PrintAllStats;
  ark_mem->step_free             = forcingStep_Free;
  ark_mem->step_printmem         = forcingStep_PrintMem;
  ark_mem->step_mem              = (void*)step_mem;

  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Error setting default solver options");
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to initialize main ARKODE infrastructure");
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }

  ARKodeSetInterpolantType(ark_mem, ARK_INTERP_LAGRANGE);

  return ark_mem;
}

/*------------------------------------------------------------------------------
  This routine re-initializes the ForcingStep module to solve a new problem of
  the same size as was previously solved. This routine should also be called
  when the problem dynamics or desired solvers have changed dramatically, so
  that the problem integration should resume as if started from scratch.

  Note all internal counters are set to 0 on re-initialization.
  ----------------------------------------------------------------------------*/
int ForcingStepReInit(void* arkode_mem, SUNStepper stepper1,
                      SUNStepper stepper2, sunrealtype t0, N_Vector y0)
{
  ARKodeMem ark_mem             = NULL;
  ARKodeForcingStepMem step_mem = NULL;

  int retval = forcingStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                               &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE)
  {
    arkProcessError(ark_mem, ARK_NO_MALLOC, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MALLOC);
    return ARK_NO_MALLOC;
  }

  retval = forcingStep_CheckArgs(ark_mem, stepper1, stepper2, y0);
  if (retval != ARK_SUCCESS) { return retval; }

  forcingStep_InitStepMem(step_mem, stepper1, stepper2);

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to initialize main ARKODE infrastructure");
    return retval;
  }

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  Accesses the number of times a given partition was evolved
  ----------------------------------------------------------------------------*/
int ForcingStepGetNumEvolves(void* arkode_mem, int partition, long int* evolves)
{
  ARKodeMem ark_mem             = NULL;
  ARKodeForcingStepMem step_mem = NULL;
  int retval = forcingStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                               &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (partition >= NUM_PARTITIONS)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The partition index is %i but there are only 2 partitions",
                    partition);
    return ARK_ILL_INPUT;
  }

  if (partition < 0)
  {
    *evolves = step_mem->n_stepper_evolves[0] + step_mem->n_stepper_evolves[1];
  }
  else { *evolves = step_mem->n_stepper_evolves[partition]; }

  return ARK_SUCCESS;
}
