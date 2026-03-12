/*------------------------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 *------------------------------------------------------------------------------
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
 *------------------------------------------------------------------------------
 * This is the implementation file for ARKODE's operator splitting module
 *----------------------------------------------------------------------------*/

#include <string.h>

#include <arkode/arkode_splittingstep.h>
#include <sundials/sundials_nvector.h>

#include "arkode_impl.h"
#include "arkode_splittingstep_impl.h"
#include "sundials_utils.h"

/*------------------------------------------------------------------------------
  Shortcut routine to unpack step_mem structure from ark_mem. If missing it
  returns ARK_MEM_NULL.
  ----------------------------------------------------------------------------*/
static int splittingStep_AccessStepMem(ARKodeMem ark_mem, const char* fname,
                                       ARKodeSplittingStepMem* step_mem)
{
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    "Time step module memory is NULL.");
    return ARK_MEM_NULL;
  }
  *step_mem = (ARKodeSplittingStepMem)ark_mem->step_mem;
  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  Shortcut routine to unpack ark_mem and step_mem structures from void* pointer.
  If either is missing it returns ARK_MEM_NULL.
  ----------------------------------------------------------------------------*/
static int splittingStep_AccessARKODEStepMem(void* arkode_mem, const char* fname,
                                             ARKodeMem* ark_mem,
                                             ARKodeSplittingStepMem* step_mem)
{
  /* access ARKodeMem structure */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_ARK_NO_MEM);
    return ARK_MEM_NULL;
  }
  *ark_mem = (ARKodeMem)arkode_mem;

  return splittingStep_AccessStepMem(*ark_mem, __func__, step_mem);
}

/*------------------------------------------------------------------------------
  This routine determines the splitting coefficients to use based on the desired
  accuracy.
  ----------------------------------------------------------------------------*/
static int splittingStep_SetCoefficients(ARKodeMem ark_mem,
                                         ARKodeSplittingStepMem step_mem)
{
  if (step_mem->coefficients != NULL) { return ARK_SUCCESS; }

  if (step_mem->order <= 1)
  {
    /* Lie-Trotter is the default (order < 1) */
    step_mem->coefficients =
      SplittingStepCoefficients_LieTrotter(step_mem->partitions);
  }
  else if (step_mem->order == 3)
  {
    step_mem->coefficients =
      SplittingStepCoefficients_ThirdOrderSuzuki(step_mem->partitions);
  }
  else if (step_mem->order % 2 == 0)
  {
    /* Triple jump only works for even order */
    step_mem->coefficients =
      SplittingStepCoefficients_TripleJump(step_mem->partitions, step_mem->order);
  }
  else
  {
    /* Bump the order up to be even but with a warning */
    int new_order = step_mem->order + 1;
    arkProcessError(ark_mem, ARK_WARNING, __LINE__, __func__, __FILE__,
                    "No splitting method at requested order, using q=%i.",
                    new_order);
    step_mem->coefficients =
      SplittingStepCoefficients_TripleJump(step_mem->partitions, new_order);
  }

  if (step_mem->coefficients == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Failed to allocate splitting coefficients");
    return ARK_MEM_FAIL;
  }

  return ARK_SUCCESS;
}

/*-----------------------------------------------------------------------------
  This routine is called just prior to performing internal time steps (after all
  user "set" routines have been called) from within arkInitialSetup.

  With initialization types FIRST_INIT this routine:
  - sets/checks the splitting coefficients to be used

  With other initialization types, this routine does nothing.
  ----------------------------------------------------------------------------*/
static int splittingStep_Init(ARKodeMem ark_mem,
                              SUNDIALS_MAYBE_UNUSED sunrealtype tout,
                              int init_type)
{
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (ark_mem->interp_type == ARK_INTERP_HERMITE)
  {
    for (int i = 0; i < step_mem->partitions; i++)
    {
      if (step_mem->steppers[i]->ops->fullrhs == NULL)
      {
        arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                        __FILE__, "steppers[%d] must implement SUNStepper_FullRhs when using Hermite interpolation",
                        i);
        return ARK_ILL_INPUT;
      }
    }
  }

  /* immediately return if resize or reset */
  if (init_type == RESIZE_INIT || init_type == RESET_INIT)
  {
    return ARK_SUCCESS;
  }

  /* assume fixed step size */
  if (!ark_mem->fixedstep)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Adaptive outer time stepping is not currently supported");
    return ARK_ILL_INPUT;
  }

  retval = splittingStep_SetCoefficients(ark_mem, step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  ark_mem->interp_degree =
    SUNMAX(1, SUNMIN(step_mem->coefficients->order - 1, ark_mem->interp_degree));

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  This is just a wrapper to call the user-supplied RHS function,
  f^1(t,y) + f^2(t,y) + ... + f^P(t,y).

  This will be called in one of three 'modes':

     ARK_FULLRHS_START -> called at the beginning of a simulation i.e., at
                          (tn, yn) = (t0, y0) or (tR, yR)

     ARK_FULLRHS_END   -> called at the end of a successful step i.e, at
                          (tcur, ycur) or the start of the subsequent step i.e.,
                          at (tn, yn) = (tcur, ycur) from the end of the last
                          step

     ARK_FULLRHS_OTHER -> called elsewhere (e.g. for dense output)

  In SplittingStep, we accumulate the RHS functions in ARK_FULLRHS_OTHER mode.
  Generally, inner steppers will not have the correct yn when this function is
  called and will not be able to reuse a function evaluation since their state
  resets at the next SUNStepper_Evolve call.
  ----------------------------------------------------------------------------*/
static int splittingStep_FullRHS(ARKodeMem ark_mem, sunrealtype t, N_Vector y,
                                 N_Vector f, SUNDIALS_MAYBE_UNUSED int mode)
{
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  for (int i = 0; i < step_mem->partitions; i++)
  {
    SUNErrCode err = SUNStepper_FullRhs(step_mem->steppers[i], t, y,
                                        i == 0 ? f : ark_mem->tempv1,
                                        SUN_FULLRHS_OTHER);
    if (err != SUN_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_RHSFUNC_FAILED, t);
      return ARK_RHSFUNC_FAIL;
    }
    if (i > 0) { N_VLinearSum(ONE, f, ONE, ark_mem->tempv1, f); }
  }

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  This routine performs a sequential operator splitting method
  ----------------------------------------------------------------------------*/
static int splittingStep_SequentialMethod(ARKodeMem ark_mem,
                                          ARKodeSplittingStepMem step_mem,
                                          int i, N_Vector y)
{
  SplittingStepCoefficients coefficients = step_mem->coefficients;

  for (int j = 0; j < coefficients->stages; j++)
  {
    SUNLogInfo(ARK_LOGGER, "begin-stages-list", "stage = %i", j);

    for (int k = 0; k < coefficients->partitions; k++)
    {
      sunrealtype beta_start = coefficients->beta[i][j][k];
      sunrealtype beta_end   = coefficients->beta[i][j + 1][k];

      if (beta_start == beta_end) { continue; }

      sunrealtype t_start = ark_mem->tn + beta_start * ark_mem->h;
      sunrealtype t_end   = ark_mem->tn + beta_end * ark_mem->h;

      SUNLogInfo(ARK_LOGGER, "begin-partitions-list",
                 "partition = %i, t_start = " SUN_FORMAT_G
                 ", t_end = " SUN_FORMAT_G,
                 k, t_start, t_end);

      SUNStepper stepper = step_mem->steppers[k];
      /* TODO(SBR): A potential future optimization is removing this reset and
       * a call to SUNStepper_SetStopTime later for methods that start a step
       * evolving the same partition the last step ended with (essentially a
       * FSAL property). Care is needed when a reset occurs, the step direction
       * changes, the coefficients change, etc. */
      SUNErrCode err = SUNStepper_Reset(stepper, t_start, y);
      if (err != SUN_SUCCESS)
      {
        SUNLogInfo(ARK_LOGGER, "end-partitions-list",
                   "status = failed stepper reset, err = %i", err);
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed partition, err = %i", err);
        return ARK_SUNSTEPPER_ERR;
      }

      err = SUNStepper_SetStepDirection(stepper, t_end - t_start);
      if (err != SUN_SUCCESS)
      {
        SUNLogInfo(ARK_LOGGER, "end-partitions-list",
                   "status = failed set direction, err = %i", err);
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed partition, err = %i", err);
        return ARK_SUNSTEPPER_ERR;
      }

      err = SUNStepper_SetStopTime(stepper, t_end);
      if (err != SUN_SUCCESS)
      {
        SUNLogInfo(ARK_LOGGER, "end-partitions-list",
                   "status = failed set stop time, err = %i", err);
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed partition, err = %i", err);
        return ARK_SUNSTEPPER_ERR;
      }

      sunrealtype tret = ZERO;
      err              = SUNStepper_Evolve(stepper, t_end, y, &tret);
      SUNLogExtraDebugVec(ARK_LOGGER, "partition state", y, "y_par(:) =");
      if (err != SUN_SUCCESS)
      {
        SUNLogInfo(ARK_LOGGER, "end-partitions-list",
                   "status = failed evolve, err = %i", err);
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed partition, err = %i", err);
        return ARK_SUNSTEPPER_ERR;
      }
      step_mem->n_stepper_evolves[k]++;

      SUNLogInfo(ARK_LOGGER, "end-partitions-list", "status = success");
    }
    SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");
  }

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  This routine performs a single step of the splitting method.
  ----------------------------------------------------------------------------*/
static int splittingStep_TakeStep(ARKodeMem ark_mem, sunrealtype* dsmPtr,
                                  int* nflagPtr)
{
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  *nflagPtr = ARK_SUCCESS; /* No algebraic solver */
  *dsmPtr   = ZERO;        /* No error estimate */

  SplittingStepCoefficients coefficients = step_mem->coefficients;

  SUNLogInfo(ARK_LOGGER, "begin-sequential-methods-list",
             "sequential method = 0");

  N_VScale(ONE, ark_mem->yn, ark_mem->ycur);
  retval = splittingStep_SequentialMethod(ark_mem, step_mem, 0, ark_mem->ycur);
  SUNLogExtraDebugVec(ARK_LOGGER, "sequential state", ark_mem->ycur,
                      "y_seq(:) =");
  if (retval != ARK_SUCCESS)
  {
    SUNLogInfo(ARK_LOGGER, "end-sequential-methods-list",
               "status = failed sequential method, retval = %i", retval);
    return retval;
  }

  if (coefficients->alpha[0] != ONE)
  {
    N_VScale(coefficients->alpha[0], ark_mem->ycur, ark_mem->ycur);
  }
  SUNLogExtraDebugVec(ARK_LOGGER, "current state", ark_mem->ycur, "y_cur(:) =");
  SUNLogInfo(ARK_LOGGER, "end-sequential-methods-list", "status = success");

  for (int i = 1; i < coefficients->sequential_methods; i++)
  {
    SUNLogInfo(ARK_LOGGER, "begin-sequential-methods-list",
               "sequential method = %i", i);

    N_VScale(ONE, ark_mem->yn, ark_mem->tempv1);
    retval = splittingStep_SequentialMethod(ark_mem, step_mem, i,
                                            ark_mem->tempv1);
    SUNLogExtraDebugVec(ARK_LOGGER, "sequential state", ark_mem->tempv1,
                        "y_seq(:) =");
    if (retval != ARK_SUCCESS)
    {
      SUNLogInfo(ARK_LOGGER, "end-sequential-methods-list",
                 "status = failed sequential method, retval = %i", retval);
      return retval;
    }
    N_VLinearSum(ONE, ark_mem->ycur, coefficients->alpha[i], ark_mem->tempv1,
                 ark_mem->ycur);

    SUNLogExtraDebugVec(ARK_LOGGER, "current state", ark_mem->ycur, "y_cur(:) =");
    SUNLogInfo(ARK_LOGGER, "end-sequential-methods-list", "status = success");
  }

  SUNLogExtraDebugVec(ARK_LOGGER, "current state", ark_mem->ycur, "y_cur(:) =");

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  Prints integrator statistics
  ----------------------------------------------------------------------------*/
static int splittingStep_PrintAllStats(ARKodeMem ark_mem, FILE* outfile,
                                       SUNOutputFormat fmt)
{
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  char name_buf[SUN_TABLE_WIDTH];
  for (int k = 0; k < step_mem->partitions; k++)
  {
    snprintf(name_buf, sizeof(name_buf), "Partition %i evolves", k + 1);
    sunfprintf_long(outfile, fmt, SUNFALSE, name_buf,
                    step_mem->n_stepper_evolves[k]);
  }

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  Outputs all solver parameters to the provided file pointer.
  ----------------------------------------------------------------------------*/
static int splittingStep_WriteParameters(ARKodeMem ark_mem, FILE* fp)
{
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  fprintf(fp, "SplittingStep time step module parameters:\n  Method order %i\n\n",
          step_mem->order);

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  Frees all SplittingStep memory.
  ----------------------------------------------------------------------------*/
static void splittingStep_Free(ARKodeMem ark_mem)
{
  ARKodeSplittingStepMem step_mem = (ARKodeSplittingStepMem)ark_mem->step_mem;
  if (step_mem != NULL)
  {
    if (step_mem->steppers != NULL) { free(step_mem->steppers); }
    if (step_mem->n_stepper_evolves != NULL)
    {
      free(step_mem->n_stepper_evolves);
    }
    SplittingStepCoefficients_Destroy(&step_mem->coefficients);
    free(step_mem);
  }
  ark_mem->step_mem = NULL;
}

/*------------------------------------------------------------------------------
  This routine outputs the memory from the SplittingStep structure to a
  specified file pointer (useful when debugging).
  ----------------------------------------------------------------------------*/
static void splittingStep_PrintMem(ARKodeMem ark_mem, FILE* outfile)
{
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return; }

  /* output integer quantities */
  fprintf(outfile, "SplittingStep: partitions = %i\n", step_mem->partitions);
  fprintf(outfile, "SplittingStep: order = %i\n", step_mem->order);

  /* output long integer quantities */
  for (int k = 0; k < step_mem->partitions; k++)
  {
    fprintf(outfile, "SplittingStep: partition %i: n_stepper_evolves = %li\n",
            k, step_mem->n_stepper_evolves[k]);
  }

  /* output sunrealtype quantities */
  fprintf(outfile, "SplittingStep: Coefficients:\n");
  SplittingStepCoefficients_Write(step_mem->coefficients, outfile);
}

/*------------------------------------------------------------------------------
  Specifies the method order
  ----------------------------------------------------------------------------*/
static int splittingStep_SetOrder(ARKodeMem ark_mem, int order)
{
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  /* set user-provided value, or default, depending on argument */
  step_mem->order = SUNMAX(1, order);

  SplittingStepCoefficients_Destroy(&step_mem->coefficients);

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  Routine to set SplittingStep options
  ----------------------------------------------------------------------------*/
static int splittingStep_SetOptions(ARKodeMem ark_mem, int* argidx, char* argv[],
                                    size_t offset, sunbooleantype* arg_used)
{
  /* The only SplittingStep-specific "Set" routine takes a custom set of
     coefficients; however, these may be specified by name, so here we'll support
     a key to specify the SplittingStepCoefficients by name,
     create the coefficients with that name, attach it to SplittingStep (who copies its
     values), and then frees the coefficients. */
  if (strcmp(argv[*argidx] + offset, "splitting_coefficients_name") == 0)
  {
    (*argidx)++;
    SplittingStepCoefficients Coefficients =
      SplittingStepCoefficients_LoadCoefficientsByName(argv[*argidx]);
    if (Coefficients == NULL)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "error setting key %s %s (invalid coefficients name)",
                      argv[(*argidx) - 1], argv[*argidx]);
      return ARK_ILL_INPUT;
    }
    int retval = SplittingStepSetCoefficients(ark_mem, Coefficients);
    SplittingStepCoefficients_Destroy(&Coefficients);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting key %s %s (SetCoefficients failed)",
                      argv[(*argidx) - 1], argv[*argidx]);
      return retval;
    }
    *arg_used = SUNTRUE;
    return ARK_SUCCESS;
  }

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  Resets all SplittingStep optional inputs to their default values. Does not
  change problem-defining function pointers or user_data pointer.
  ----------------------------------------------------------------------------*/
static int splittingStep_SetDefaults(ARKodeMem ark_mem)
{
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  return splittingStep_SetOrder(ark_mem, 0);
}

/*------------------------------------------------------------------------------
  This routine checks if all required SUNStepper operations are present. If any
  of them are missing it return SUNFALSE.
  ----------------------------------------------------------------------------*/
static sunbooleantype splittingStep_CheckSUNStepper(SUNStepper stepper)
{
  SUNStepper_Ops ops = stepper->ops;
  return ops->evolve != NULL && ops->reset != NULL &&
         ops->setstoptime != NULL && ops->setstepdirection != NULL;
}

/*------------------------------------------------------------------------------
  This routine validates arguments when (re)initializing a SplittingStep
  integrator
  ----------------------------------------------------------------------------*/
static int splittingStep_CheckArgs(ARKodeMem ark_mem, SUNStepper* steppers,
                                   int partitions, N_Vector y0)
{
  if (steppers == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "steppers = NULL illegal.");
    return ARK_ILL_INPUT;
  }

  if (partitions <= 1)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The number of partitions must be greater than one");
    return ARK_ILL_INPUT;
  }

  for (int i = 0; i < partitions; i++)
  {
    if (steppers[i] == NULL)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "steppers[%d] = NULL illegal.", i);
      return ARK_ILL_INPUT;
    }

    if (!splittingStep_CheckSUNStepper(steppers[i]))
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "stepper[%d] does not implement the required operations.",
                      i);
      return ARK_ILL_INPUT;
    }
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
static int splittingStep_InitStepMem(ARKodeMem ark_mem,
                                     ARKodeSplittingStepMem step_mem,
                                     SUNStepper* steppers, int partitions)
{
  if (step_mem->steppers != NULL) { free(step_mem->steppers); }
  step_mem->steppers = malloc(partitions * sizeof(*steppers));
  if (step_mem->steppers == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    return ARK_MEM_FAIL;
  }
  memcpy(step_mem->steppers, steppers, partitions * sizeof(*steppers));

  if (step_mem->n_stepper_evolves != NULL)
  {
    free(step_mem->n_stepper_evolves);
  }
  step_mem->n_stepper_evolves = calloc(partitions,
                                       sizeof(*step_mem->n_stepper_evolves));

  /* If the number of partitions changed, the coefficients are no longer
   * compatible and must be cleared. If a user previously called ARKodeSetOrder
   * that will still be respected at the next call to ARKodeEvolve */
  if (step_mem->partitions != partitions)
  {
    SplittingStepCoefficients_Destroy(&step_mem->coefficients);
  }
  step_mem->partitions = partitions;

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  Creates the SplittingStep integrator
  ---------------------------------------------------------------*/
void* SplittingStepCreate(SUNStepper* steppers, int partitions, sunrealtype t0,
                          N_Vector y0, SUNContext sunctx)
{
  int retval = splittingStep_CheckArgs(NULL, steppers, partitions, y0);
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

  ARKodeSplittingStepMem step_mem =
    (ARKodeSplittingStepMem)malloc(sizeof(*step_mem));
  if (step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }

  step_mem->partitions        = partitions;
  step_mem->order             = 0;
  step_mem->steppers          = NULL;
  step_mem->n_stepper_evolves = NULL;
  step_mem->coefficients      = NULL;
  retval = splittingStep_InitStepMem(ark_mem, step_mem, steppers, partitions);
  if (retval != ARK_SUCCESS)
  {
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }

  /* Attach step_mem structure and function pointers to ark_mem */
  ark_mem->step_init            = splittingStep_Init;
  ark_mem->step_fullrhs         = splittingStep_FullRHS;
  ark_mem->step                 = splittingStep_TakeStep;
  ark_mem->step_printallstats   = splittingStep_PrintAllStats;
  ark_mem->step_writeparameters = splittingStep_WriteParameters;
  ark_mem->step_free            = splittingStep_Free;
  ark_mem->step_printmem        = splittingStep_PrintMem;
  ark_mem->step_setoptions      = splittingStep_SetOptions;
  ark_mem->step_setdefaults     = splittingStep_SetDefaults;
  ark_mem->step_setorder        = splittingStep_SetOrder;
  ark_mem->step_mem             = (void*)step_mem;

  /* Set default values for ARKStep optional inputs */
  retval = splittingStep_SetDefaults(ark_mem);
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
  This routine re-initializes the SplittingStep module to solve a new problem of
  the same size as was previously solved. This routine should also be called
  when the problem dynamics or desired solvers have changed dramatically, so
  that the problem integration should resume as if started from scratch.

  Note all internal counters are set to 0 on re-initialization.
  ----------------------------------------------------------------------------*/
int SplittingStepReInit(void* arkode_mem, SUNStepper* steppers, int partitions,
                        sunrealtype t0, N_Vector y0)
{
  ARKodeMem ark_mem               = NULL;
  ARKodeSplittingStepMem step_mem = NULL;

  int retval = splittingStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                                 &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE)
  {
    arkProcessError(ark_mem, ARK_NO_MALLOC, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MALLOC);
    return ARK_NO_MALLOC;
  }

  retval = splittingStep_CheckArgs(ark_mem, steppers, partitions, y0);
  if (retval != ARK_SUCCESS) { return retval; }

  splittingStep_InitStepMem(ark_mem, step_mem, steppers, partitions);

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

/*---------------------------------------------------------------
  Sets the SplittingStep coefficients.
  ---------------------------------------------------------------*/
int SplittingStepSetCoefficients(void* arkode_mem,
                                 SplittingStepCoefficients coefficients)
{
  ARKodeMem ark_mem               = NULL;
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                                 &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (coefficients == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Splitting coefficients must be non-NULL");
    return ARK_ILL_INPUT;
  }

  if (step_mem->partitions != coefficients->partitions)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "The splitting method has %i partitions but the coefficients have %i.",
                    step_mem->partitions, coefficients->partitions);
    return ARK_ILL_INPUT;
  }

  SplittingStepCoefficients_Destroy(&step_mem->coefficients);
  step_mem->coefficients = SplittingStepCoefficients_Copy(coefficients);
  if (step_mem->coefficients == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Failed to copy splitting coefficients");
    return ARK_MEM_NULL;
  }

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  Accesses the number of times a given partition was evolved
  ----------------------------------------------------------------------------*/
int SplittingStepGetNumEvolves(void* arkode_mem, int partition, long int* evolves)
{
  ARKodeMem ark_mem               = NULL;
  ARKodeSplittingStepMem step_mem = NULL;
  int retval = splittingStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                                 &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (partition >= step_mem->partitions)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "The partition index is %i but there are only %i partitions",
                    partition, step_mem->partitions);
    return ARK_ILL_INPUT;
  }

  if (partition < 0)
  {
    *evolves = 0;
    for (int k = 0; k < step_mem->partitions; k++)
    {
      *evolves += step_mem->n_stepper_evolves[k];
    }
  }
  else { *evolves = step_mem->n_stepper_evolves[partition]; }

  return ARK_SUCCESS;
}
