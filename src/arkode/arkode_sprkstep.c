/*---------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *---------------------------------------------------------------
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
 *---------------------------------------------------------------
 * This is the implementation file for ARKODE's SPRK time stepper
 * module.
 *--------------------------------------------------------------*/

#include "arkode/arkode_sprkstep.h"

#include <arkode/arkode.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#include "arkode_impl.h"
#include "arkode_interp_impl.h"
#include "arkode_sprkstep_impl.h"

/*===============================================================
  Exported functions
  ===============================================================*/

void* SPRKStepCreate(ARKRhsFn f1, ARKRhsFn f2, sunrealtype t0, N_Vector y0,
                     SUNContext sunctx)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* Check that f1 and f2 are supplied */
  if (!f1)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_F);
    return (NULL);
  }

  if (!f2)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_F);
    return (NULL);
  }

  /* Check for legal input parameters */
  if (!y0)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_Y0);
    return (NULL);
  }

  if (!sunctx)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_SUNCTX);
    return (NULL);
  }

  /* Create ark_mem structure and set default values */
  ark_mem = arkCreate(sunctx);
  if (ark_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (NULL);
  }

  /* Allocate ARKodeSPRKStepMem structure, and initialize to zero */
  step_mem = NULL;
  step_mem = (ARKodeSPRKStepMem)malloc(sizeof(struct ARKodeSPRKStepMemRec));
  if (step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    return (NULL);
  }
  memset(step_mem, 0, sizeof(struct ARKodeSPRKStepMemRec));

  /* Allocate vectors in stepper mem */
  if (!arkAllocVec(ark_mem, y0, &(step_mem->sdata)))
  {
    ARKodeFree((void**)&ark_mem);
    return (NULL);
  }

  if (ark_mem->use_compensated_sums)
  {
    if (!arkAllocVec(ark_mem, y0, &(step_mem->yerr)))
    {
      ARKodeFree((void**)&ark_mem);
      return (NULL);
    }
    /* Zero yerr for compensated summation */
    N_VConst(ZERO, step_mem->yerr);
  }
  else { step_mem->yerr = NULL; }
  ark_mem->step_init                  = sprkStep_Init;
  ark_mem->step_fullrhs               = sprkStep_FullRHS;
  ark_mem->step                       = sprkStep_TakeStep;
  ark_mem->step_printallstats         = sprkStep_PrintAllStats;
  ark_mem->step_writeparameters       = sprkStep_WriteParameters;
  ark_mem->step_setusecompensatedsums = sprkStep_SetUseCompensatedSums;
  ark_mem->step_setoptions            = sprkStep_SetOptions;
  ark_mem->step_resize                = sprkStep_Resize;
  ark_mem->step_free                  = sprkStep_Free;
  ark_mem->step_setdefaults           = sprkStep_SetDefaults;
  ark_mem->step_setorder              = sprkStep_SetOrder;
  ark_mem->step_getnumrhsevals        = sprkStep_GetNumRhsEvals;
  ark_mem->step_mem                   = (void*)step_mem;

  /* Set default values for optional inputs */
  retval = sprkStep_SetDefaults((void*)ark_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Error setting default solver options");
    ARKodeFree((void**)&ark_mem);
    return (NULL);
  }

  /* Copy the input parameters into ARKODE state */
  step_mem->f1 = f1;
  step_mem->f2 = f2;

  /* Initialize the counters */
  step_mem->nf1    = 0;
  step_mem->nf2    = 0;
  step_mem->istage = 0;

  /* SPRKStep uses Lagrange interpolation by default, since Hermite is
     less compatible with these methods. */
  ARKodeSetInterpolantType(ark_mem, ARK_INTERP_LAGRANGE);

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to initialize main ARKODE infrastructure");
    ARKodeFree((void**)&ark_mem);
    return (NULL);
  }

  return ((void*)ark_mem);
}

/*---------------------------------------------------------------
  SPRKStepReInit:

  This routine re-initializes the SPRKStep module to solve a new
  problem of the same size as was previously solved. This routine
  should also be called when the problem dynamics or desired solvers
  have changed dramatically, so that the problem integration should
  resume as if started from scratch.

  Note all internal counters are set to 0 on re-initialization.
  ---------------------------------------------------------------*/
int SPRKStepReInit(void* arkode_mem, ARKRhsFn f1, ARKRhsFn f2, sunrealtype t0,
                   N_Vector y0)
{
  ARKodeMem ark_mem          = NULL;
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeMem and ARKodeSPRKStepMem structures */
  retval = sprkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE)
  {
    arkProcessError(ark_mem, ARK_NO_MALLOC, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MALLOC);
    return (ARK_NO_MALLOC);
  }

  /* Check that f1 and f2 are supplied */
  if (!f1 || !f2)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_F);
    return (ARK_ILL_INPUT);
  }

  /* Check that y0 is supplied */
  if (!y0)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_Y0);
    return (ARK_ILL_INPUT);
  }

  /* Copy the input parameters into ARKODE state */
  step_mem->f1 = f1;
  step_mem->f2 = f2;

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to reinitialize main ARKODE infrastructure");
    return (retval);
  }

  /* Initialize the counters */
  step_mem->nf1    = 0;
  step_mem->nf2    = 0;
  step_mem->istage = 0;

  /* Zero yerr for compensated summation */
  if (ark_mem->use_compensated_sums) { N_VConst(ZERO, step_mem->yerr); }

  return (ARK_SUCCESS);
}

/*===============================================================
  Interface routines supplied to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  sprkStep_Resize:

  This routine resizes the memory within the SPRKStep module.
  ---------------------------------------------------------------*/
int sprkStep_Resize(ARKodeMem ark_mem, N_Vector y0,
                    SUNDIALS_MAYBE_UNUSED sunrealtype hscale,
                    SUNDIALS_MAYBE_UNUSED sunrealtype t0, ARKVecResizeFn resize,
                    void* resize_data)
{
  ARKodeSPRKStepMem step_mem = NULL;
  sunindextype lrw1, liw1, lrw_diff, liw_diff;
  int retval;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Determine change in vector sizes */
  lrw1 = liw1 = 0;
  if (y0->ops->nvspace != NULL) { N_VSpace(y0, &lrw1, &liw1); }
  lrw_diff      = lrw1 - ark_mem->lrw1;
  liw_diff      = liw1 - ark_mem->liw1;
  ark_mem->lrw1 = lrw1;
  ark_mem->liw1 = liw1;

  /* Resize the local vectors */
  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, y0,
                    &step_mem->sdata))
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Unable to resize vector");
    return (ARK_MEM_FAIL);
  }

  if (step_mem->yerr)
  {
    if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, y0,
                      &step_mem->yerr))
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to resize vector");
      return (ARK_MEM_FAIL);
    }
    /* Zero yerr for compensated summation */
    N_VConst(ZERO, step_mem->yerr);
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  sprkStep_Reset:

  This routine resets the SPRKStep module state to solve the same
  problem from the given time with the input state (all counter
  values are retained).
  ---------------------------------------------------------------*/
int sprkStep_Reset(ARKodeMem ark_mem, SUNDIALS_MAYBE_UNUSED sunrealtype tR,
                   SUNDIALS_MAYBE_UNUSED N_Vector yR)
{
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  if (ark_mem->use_compensated_sums)
  {
    N_VConst(SUN_RCONST(0.0), step_mem->yerr);
  }
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  sprkStep_Free frees all SPRKStep memory.
  ---------------------------------------------------------------*/
void sprkStep_Free(ARKodeMem ark_mem)
{
  ARKodeSPRKStepMem step_mem = NULL;

  /* nothing to do if ark_mem is already NULL */
  if (ark_mem == NULL) { return; }

  /* conditional frees on non-NULL SPRKStep module */
  if (ark_mem->step_mem != NULL)
  {
    step_mem = (ARKodeSPRKStepMem)ark_mem->step_mem;

    if (step_mem->sdata != NULL)
    {
      arkFreeVec(ark_mem, &step_mem->sdata);
      step_mem->sdata = NULL;
    }

    if (step_mem->yerr != NULL)
    {
      arkFreeVec(ark_mem, &step_mem->yerr);
      step_mem->yerr = NULL;
    }

    ARKodeSPRKTable_Free(step_mem->method);

    free(ark_mem->step_mem);
    ark_mem->step_mem = NULL;
  }
}

/*---------------------------------------------------------------
  sprkStep_Init:

  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup.

  For all initialization types, this routine sets the relevant
  TakeStep routine based on the current problem configuration.

  With initialization type FIRST_INIT or RESIZE_INIT, this routine
  this routines loads the default method of the selected order
  if necessary.

  With initialization type RESET_INIT, this routine does nothing.
  ---------------------------------------------------------------*/
int sprkStep_Init(ARKodeMem ark_mem, SUNDIALS_MAYBE_UNUSED sunrealtype tout,
                  int init_type)
{
  ARKodeSPRKStepMem step_mem = NULL;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* immediately return if reset */
  if (init_type == RESET_INIT) { return (ARK_SUCCESS); }

  /* initializations/checks for (re-)initialization call */
  if (init_type == FIRST_INIT)
  {
    if (!step_mem->method)
    {
      switch (step_mem->q)
      {
      case 1:
        step_mem->method = ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_1);
        break;
      case 2:
        step_mem->method = ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_2);
        break;
      case 3:
        step_mem->method = ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_3);
        break;
      case 4:
        step_mem->method = ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_4);
        break;
      case 5:
        step_mem->method = ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_5);
        break;
      case 6:
        step_mem->method = ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_6);
        break;
      case 7:
      case 8:
        step_mem->method = ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_8);
        break;
      case 9:
      case 10:
        step_mem->method = ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_10);
        break;
      default:
        arkProcessError(ark_mem, ARK_WARNING, __LINE__, __func__, __FILE__,
                        "No SPRK method at requested order, using q=4.");
        step_mem->method = ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_4);
        break;
      }
    }
  }

  /* Override the interpolant degree (if needed), used in arkInitialSetup */
  if (step_mem->method->q > 1 &&
      ark_mem->interp_degree > (step_mem->method->q - 1))
  {
    /* Limit max degree to at most one less than the method global order */
    ark_mem->interp_degree = step_mem->method->q - 1;
  }
  else if (step_mem->method->q == 1 && ark_mem->interp_degree > 1)
  {
    /* Allow for linear interpolant with first order methods to ensure
       solution values are returned at the time interval end points */
    ark_mem->interp_degree = 1;
  }

  /* Zero yerr for compensated summation */
  if (ark_mem->use_compensated_sums) { N_VConst(ZERO, step_mem->yerr); }

  return (ARK_SUCCESS);
}

/* Utility to call f1 and increment the counter */
inline int sprkStep_f1(ARKodeSPRKStepMem step_mem, sunrealtype tcur,
                       N_Vector ycur, N_Vector f1, void* user_data)
{
  int retval = step_mem->f1(tcur, ycur, f1, user_data);
  step_mem->nf1++;
  return retval;
}

/* Utility to call f2 and increment the counter */
inline int sprkStep_f2(ARKodeSPRKStepMem step_mem, sunrealtype tcur,
                       N_Vector ycur, N_Vector f2, void* user_data)
{
  int retval = step_mem->f2(tcur, ycur, f2, user_data);
  step_mem->nf2++;
  return retval;
}

/*------------------------------------------------------------------------------
  sprkStep_FullRHS:

  This is just a wrapper to call the user-supplied RHS,
  f1(t,y) + f2(t,y).

  This will be called in one of three 'modes':

     ARK_FULLRHS_START -> called at the beginning of a simulation i.e., at
                          (tn, yn) = (t0, y0) or (tR, yR)

     ARK_FULLRHS_END   -> called at the end of a successful step i.e, at
                          (tcur, ycur) or the start of the subsequent step i.e.,
                          at (tn, yn) = (tcur, ycur) from the end of the last
                          step

     ARK_FULLRHS_OTHER -> called elsewhere (e.g. for dense output)

  Since RHS values are not stored in SPRKStep we evaluate the RHS functions for
  all modes.
  ----------------------------------------------------------------------------*/
int sprkStep_FullRHS(ARKodeMem ark_mem, sunrealtype t, N_Vector y, N_Vector f,
                     int mode)
{
  int retval                 = 0;
  ARKodeSPRKStepMem step_mem = NULL;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* perform RHS functions contingent on 'mode' argument */
  switch (mode)
  {
  case ARK_FULLRHS_START:
  case ARK_FULLRHS_END:
  case ARK_FULLRHS_OTHER:

    /* Since f1 and f2 do not have overlapping outputs and so the f vector is
       passed to both RHS functions. */

    retval = sprkStep_f1(step_mem, t, y, f, ark_mem->user_data);
    if (retval != 0)
    {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_RHSFUNC_FAILED, t);
      return (ARK_RHSFUNC_FAIL);
    }

    retval = sprkStep_f2(step_mem, t, y, f, ark_mem->user_data);
    if (retval != 0)
    {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_RHSFUNC_FAILED, t);
      return (ARK_RHSFUNC_FAIL);
    }

    break;

  default:
    /* return with RHS failure if unknown mode is passed */
    arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                    "Unknown full RHS mode");
    return (ARK_RHSFUNC_FAIL);
  }

  return (ARK_SUCCESS);
}

/* Standard formulation of SPRK.
   This requires only 2 vectors in principle, but we use three
   since we persist the stage data. Only the stage data vector
   belongs to SPRKStep, the other two are reused from the ARKODE core. */
int sprkStep_TakeStep(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  ARKodeSPRKStepMem step_mem = NULL;
  N_Vector prev_stage        = NULL;
  N_Vector curr_stage        = NULL;
  sunrealtype ci             = SUN_RCONST(0.0);
  sunrealtype chati          = SUN_RCONST(0.0);
  int is                     = 0;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  prev_stage = ark_mem->yn;
  curr_stage = ark_mem->ycur;
  for (is = 0; is < step_mem->method->stages; is++)
  {
    /* load/compute coefficients */
    sunrealtype ai    = step_mem->method->a[is];
    sunrealtype ahati = step_mem->method->ahat[is];

    ci += ai;
    chati += ahati;

    /* store current stage index */
    step_mem->istage = is;

    SUNLogInfo(ARK_LOGGER, "begin-stages-list",
               "stage = %i, t = " SUN_FORMAT_G ", that = " SUN_FORMAT_G, is,
               ark_mem->tn + ci * ark_mem->h, ark_mem->tn + chati * ark_mem->h);
    SUNLogExtraDebugVec(ARK_LOGGER, "stage", prev_stage, "z2_%i(:) =", is);

    /* evaluate p' with the previous velocity */
    if (SUNRabs(ahati) > TINY)
    {
      N_VConst(ZERO, step_mem->sdata); /* either have to do this or ask user to
                                          set other outputs to zero */

      retval = sprkStep_f1(step_mem, ark_mem->tn + chati * ark_mem->h,
                           prev_stage, step_mem->sdata, ark_mem->user_data);

      SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", step_mem->sdata,
                          "f1_%i(:) =", is);

      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed rhs eval, retval = %i", retval);
        return ARK_RHSFUNC_FAIL;
      }
    }

    /* position update */
    N_VLinearSum(ONE, prev_stage, ark_mem->h * ahati, step_mem->sdata,
                 curr_stage);

    /* set current stage time(s) */
    ark_mem->tcur = ark_mem->tn + chati * ark_mem->h;

    SUNLogExtraDebugVec(ARK_LOGGER, "stage", curr_stage, "z1_%i(:) =", is);

    /* evaluate q' with the current positions and update velocity */
    if (SUNRabs(ai) > TINY)
    {
      N_VConst(ZERO, step_mem->sdata); /* either have to do this or ask user to
                                        set other outputs to zero */

      retval = sprkStep_f2(step_mem, ark_mem->tn + ci * ark_mem->h, curr_stage,
                           step_mem->sdata, ark_mem->user_data);

      SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", step_mem->sdata,
                          "f2_%i(:) =", is);

      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed rhs eval, retval = %i", retval);
        return ARK_RHSFUNC_FAIL;
      }

      /* velocity update */
      N_VLinearSum(ONE, curr_stage, ark_mem->h * ai, step_mem->sdata, curr_stage);
    }

    /* apply user-supplied stage postprocessing function (if supplied) */
    if (ark_mem->ProcessStage != NULL)
    {
      retval = ark_mem->ProcessStage(ark_mem->tcur, ark_mem->ycur,
                                     ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed postprocess stage, retval = %i", retval);
        return (ARK_POSTPROCESS_STAGE_FAIL);
      }
    }

    /* keep track of the stage number */
    step_mem->istage++;

    prev_stage = curr_stage;

    SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");
  }

  *nflagPtr = 0;
  *dsmPtr   = 0;

  SUNLogExtraDebugVec(ARK_LOGGER, "updated solution", ark_mem->ycur, "ycur(:) =");

  return ARK_SUCCESS;
}

/* Increment SPRK algorithm with compensated summation.
   This algorithm requires 6 vectors, but 5 of them are reused
   from the ARKODE core. */
int sprkStep_TakeStep_Compensated(ARKodeMem ark_mem, sunrealtype* dsmPtr,
                                  int* nflagPtr)
{
  ARKodeSPRKStepMem step_mem = NULL;
  N_Vector delta_Yi          = NULL;
  N_Vector yn_plus_delta_Yi  = NULL;
  N_Vector diff              = NULL;
  sunrealtype ci             = SUN_RCONST(0.0);
  sunrealtype chati          = SUN_RCONST(0.0);
  int is                     = 0;
  int retval                 = 0;

  /* access ARKodeSPRKStepMem structure */
  retval = sprkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Vector shortcuts */
  delta_Yi         = ark_mem->tempv1;
  yn_plus_delta_Yi = ark_mem->tempv2;
  diff             = ark_mem->tempv3;

  /* [ \Delta P_0 ] = [ 0 ]
     [ \Delta Q_0 ] = [ 0 ] */
  N_VConst(ZERO, delta_Yi);

  /* loop over internal stages to the step */
  for (is = 0; is < step_mem->method->stages; is++)
  {
    /* load/compute coefficients */
    sunrealtype ai    = step_mem->method->a[is];
    sunrealtype ahati = step_mem->method->ahat[is];

    ci += ai;
    chati += ahati;

    /* store current stage index */
    step_mem->istage = is;

    SUNLogInfo(ARK_LOGGER, "begin-stages-list",
               "stage = %i, t = " SUN_FORMAT_G ", that = " SUN_FORMAT_G, is,
               ark_mem->tn + ci * ark_mem->h, ark_mem->tn + chati * ark_mem->h);

    /* [     ] + [            ]
       [ q_n ] + [ \Delta Q_i ] */
    N_VLinearSum(ONE, ark_mem->yn, ONE, delta_Yi, yn_plus_delta_Yi);

    SUNLogExtraDebugVec(ARK_LOGGER, "stage", yn_plus_delta_Yi, "z2_%i(:) =", is);

    if (SUNRabs(ahati) > TINY)
    {
      /* Evaluate p' with the previous velocity */
      N_VConst(ZERO, step_mem->sdata); /* either have to do this or ask user to
                                          set other outputs to zero */
      retval = sprkStep_f1(step_mem, ark_mem->tn + chati * ark_mem->h,
                           yn_plus_delta_Yi, step_mem->sdata, ark_mem->user_data);

      SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", step_mem->sdata,
                          "f1_%i(:) =", is);

      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed rhs eval, retval = %i", retval);
        return (ARK_RHSFUNC_FAIL);
      }

      /* Incremental position update:
         [ \Delta P_i ] = [ \Delta P_{i-1} ] + [ sdata ]
         [            ] = [                ] + [       ] */
      N_VLinearSum(ONE, delta_Yi, ark_mem->h * ahati, step_mem->sdata, delta_Yi);
    }

    /* [ p_n ] + [ \Delta P_i ]
       [     ] + [            ] */
    N_VLinearSum(ONE, ark_mem->yn, ONE, delta_Yi, yn_plus_delta_Yi);

    SUNLogExtraDebugVec(ARK_LOGGER, "stage", yn_plus_delta_Yi, "z1_%i(:) =", is);

    /* set current stage time(s) */
    ark_mem->tcur = ark_mem->tn + chati * ark_mem->h;

    if (SUNRabs(ai) > TINY)
    {
      /* Evaluate q' with the current positions */
      N_VConst(ZERO, step_mem->sdata); /* either have to do this or ask user to
                                          set other outputs to zero */
      retval = sprkStep_f2(step_mem, ark_mem->tn + ci * ark_mem->h,
                           yn_plus_delta_Yi, step_mem->sdata, ark_mem->user_data);

      SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", step_mem->sdata,
                          "f2_%i(:) =", is);

      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed rhs eval, retval = %i", retval);
        return (ARK_RHSFUNC_FAIL);
      }

      /* Incremental velocity update:
         [            ] = [                ] + [       ]
         [ \Delta Q_i ] = [ \Delta Q_{i-1} ] + [ sdata ] */
      N_VLinearSum(ONE, delta_Yi, ark_mem->h * ai, step_mem->sdata, delta_Yi);
    }

    /* if user-supplied stage postprocessing function, we error out since it
     * won't work with the increment form */
    if (ark_mem->ProcessStage != NULL)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed postprocess stage, retval = %i", retval);
      arkProcessError(ark_mem, ARK_POSTPROCESS_STAGE_FAIL, __LINE__, __func__,
                      __FILE__,
                      "Compensated summation is not compatible with stage "
                      "PostProcessing!\n");
      return (ARK_POSTPROCESS_STAGE_FAIL);
    }

    SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");
  }

  /*
    Now we compute the step solution via compensated summation.
     [ p_{n+1} ] = [ p_n ] + [ \Delta P_i ]
     [ q_{n+1} ] = [ q_n ] + [ \Delta Q_i ] */
  N_VLinearSum(ONE, delta_Yi, -ONE, step_mem->yerr, delta_Yi);
  N_VLinearSum(ONE, ark_mem->yn, ONE, delta_Yi, ark_mem->ycur);
  N_VLinearSum(ONE, ark_mem->ycur, -ONE, ark_mem->yn, diff);
  N_VLinearSum(ONE, diff, -ONE, delta_Yi, step_mem->yerr);

  *nflagPtr = 0;
  *dsmPtr   = SUN_RCONST(0.0);

  SUNLogExtraDebugVec(ARK_LOGGER, "updated solution", ark_mem->ycur, "ycur(:) =");

  return 0;
}

/*===============================================================
  Internal utility routines
  ===============================================================*/

/*---------------------------------------------------------------
  sprkStep_AccessARKODEStepMem:

  Shortcut routine to unpack ark_mem and step_mem structures from
  void* pointer.  If either is missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int sprkStep_AccessARKODEStepMem(void* arkode_mem, const char* fname,
                                 ARKodeMem* ark_mem, ARKodeSPRKStepMem* step_mem)
{
  /* access ARKodeMem structure */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *ark_mem = (ARKodeMem)arkode_mem;

  /* access ARKodeSPRKStepMem structure */
  if ((*ark_mem)->step_mem == NULL)
  {
    arkProcessError(*ark_mem, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_SPRKSTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *step_mem = (ARKodeSPRKStepMem)(*ark_mem)->step_mem;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  sprkStep_AccessStepMem:

  Shortcut routine to unpack step_mem structure from ark_mem.
  If missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int sprkStep_AccessStepMem(ARKodeMem ark_mem, const char* fname,
                           ARKodeSPRKStepMem* step_mem)
{
  /* access ARKodeSPRKStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_SPRKSTEP_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *step_mem = (ARKodeSPRKStepMem)ark_mem->step_mem;
  return (ARK_SUCCESS);
}

/*===============================================================
  EOF
  ===============================================================*/
