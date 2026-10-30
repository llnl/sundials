/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 * This is the implementation file for MRIStep's ExtSTS utility routines.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <arkode/arkode_lsrkstep.h>
#include "arkode_impl.h"
#include "arkode_lsrkstep_impl.h"
#include "arkode_mristep_impl.h"

/* MRIStep ExtSTS constructor routine */
void* MRIStepCreateExtSTS(ARKRhsFn fd, ARKRhsFn fe, ARKRhsFn fi, sunrealtype t0,
                          N_Vector y0, SUNContext sunctx)
{
  /* Create LSRKStep integrator and configure with ExtSTS defaults */
  void* sts_mem = LSRKStepCreateSTS(fd, t0, y0, sunctx);
  if (sts_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Failed to create LSRKStep memory for ExtSTS method.");
    return NULL;
  }

  /* Increase the maximum number of internal stages */
  int retval = LSRKStepSetMaxNumStages(sts_mem, 10000);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to set maximum number of STS stages.");
    ARKodeFree(&sts_mem);
    return NULL;
  }

  /* Disable temporal interpolation for inner STS method */
  retval = ARKodeSetInterpolantType(sts_mem, ARK_INTERP_NONE);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__,
                    __FILE__, "Failed to disable LRSKStep interpolation in ExtSTS method.");
    ARKodeFree(&sts_mem);
    return NULL;
  }

  /* Create the inner stepper and attach objects and function pointers */
  MRIStepInnerStepper inner_stepper = NULL;
  retval = MRIStepInnerStepper_Create(sunctx, &inner_stepper);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__,
                    __FILE__, "Failed to create MRIStep inner stepper for ExtSTS method.");
    ARKodeFree(&sts_mem);
    return NULL;
  }
  retval = MRIStepInnerStepper_SetContent(inner_stepper, sts_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to set MRIStep inner stepper content.");
    ARKodeFree(&sts_mem);
    MRIStepInnerStepper_Free(&inner_stepper);
    return NULL;
  }
  retval = MRIStepInnerStepper_SetEvolveFn(inner_stepper,
                                           extSTSInnerStepper_Evolve);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to set MRIStep inner stepper evolve function.");
    ARKodeFree(&sts_mem);
    MRIStepInnerStepper_Free(&inner_stepper);
    return NULL;
  }
  retval = MRIStepInnerStepper_SetFullRhsFn(inner_stepper,
                                            ark_MRIStepInnerFullRhs);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to set MRIStep inner stepper full RHS function.");
    ARKodeFree(&sts_mem);
    MRIStepInnerStepper_Free(&inner_stepper);
    return NULL;
  }
  retval = MRIStepInnerStepper_SetResetFn(inner_stepper, ark_MRIStepInnerReset);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to set MRIStep inner stepper reset function.");
    ARKodeFree(&sts_mem);
    MRIStepInnerStepper_Free(&inner_stepper);
    return NULL;
  }

  /* Create the MRIStep integrator, attaching the inner stepper for diffusion */
  void* arkode_mem = MRIStepCreate(fe, fi, t0, y0, inner_stepper, sunctx);
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Failed to create MRIStep integrator.");
    ARKodeFree(&sts_mem);
    MRIStepInnerStepper_Free(&inner_stepper);
    return NULL;
  }

  /* Store a pointer to the inner stepper content directly in MRIStep */
  ARKodeMem ark_mem         = NULL;
  ARKodeMRIStepMem step_mem = NULL;
  retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to access stepper memory from MRIStep.");
    ARKodeFree(&arkode_mem);
    ARKodeFree(&sts_mem);
    MRIStepInnerStepper_Free(&inner_stepper);
    return NULL;
  }
  step_mem->extsts_method = SUNTRUE;

  /* Configure the MRIStep integrator with ExtSTS defaults */
  /*   Select default ExtSTS method based on provided RHS functions */
  MRIStepCoupling MRIC = NULL;
  if (fe == NULL && fi != NULL) /* Implicit ExtSTS method */
  {
    MRIC = MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_GARK_GKC21);
  }
  else if (fe != NULL && fi == NULL) /* Explicit ExtSTS method */
  {
    MRIC = MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_EXP_GKC21);
  }
  else /* ImEx ExtSTS method */
  {
    MRIC = MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_GARK_GKC21);
  }
  if (MRIC == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__,
                    __FILE__, "Failed to create MRIStep coupling table for ExtSTS method.");
    ARKodeFree(&arkode_mem);
    ARKodeFree(&sts_mem);
    MRIStepInnerStepper_Free(&inner_stepper);
    return NULL;
  }
  retval = MRIStepSetCoupling(arkode_mem, MRIC);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to set MRIStep coupling table for ExtSTS method.");
    ARKodeFree(&arkode_mem);
    ARKodeFree(&sts_mem);
    MRIStepInnerStepper_Free(&inner_stepper);
    return NULL;
  }
  MRIStepCoupling_Free(MRIC);

  /* return with the constructed MRIStep object */
  return arkode_mem;
}

/* MRIStep ExtSTS reinitialization routine */
int MRIStepReInitExtSTS(void* arkode_mem, ARKRhsFn fd, ARKRhsFn fe, ARKRhsFn fi,
                        sunrealtype t0, N_Vector y0)
{
  /* access ARKodeMem and ARKodeMRIStepMem structures */
  ARKodeMem ark_mem         = NULL;
  ARKodeMRIStepMem step_mem = NULL;
  int retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                           &step_mem);
  if (retval) { return retval; }

  /* Reinitialize the MRIStep integrator */
  retval = MRIStepReInit(arkode_mem, fe, fi, t0, y0);
  if (retval)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to reinitialize MRIStep for ExtSTS method.");
    return retval;
  }

  /* Reinitialize the LSRKStep integrator */
  retval = LSRKStepReInitSTS(step_mem->stepper->content, fd, t0, y0);
  if (retval)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to reinitialize LSRKStep for ExtSTS method.");
    return retval;
  }

  return retval;
}

/* Accessor routine for the inner LSRKStep solver from an ExtSTS method */
int MRIStepGetSTS(void* arkode_mem, void** sts_mem)
{
  /* access ARKodeMem and ARKodeMRIStepMem structures */
  ARKodeMem ark_mem         = NULL;
  ARKodeMRIStepMem step_mem = NULL;
  int retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                           &step_mem);
  if (retval) { return retval; }

  /* return pointer to stored STS integrator */
  *sts_mem = step_mem->stepper->content;
  return ARK_SUCCESS;
}

/* ExtSTS temporal error estimation routines.  The STS error committed in each
   fast solve is estimated by the RKC/RKL embedding, measured in the MRIStep
   error weights, and combined over the K fast solves of a step attempt as a
   subinterval-weighted RMS,

     extsts_dsm = sqrt( sum_k |dt_k| ||e_k||^2 / sum_k |dt_k| ),

   which is then combined with the MRI embedding estimate dsm_S as the
   root-sum-square sqrt(dsm_S^2 + extsts_dsm^2). */

/* Enable/disable the STS error estimate and zero the accumulators at the start
   of a step attempt */
int extSTS_BeginStepError(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem,
                          sunbooleantype enable)
{
  ARKodeLSRKStepMem lsrkstep_mem = NULL;
  int retval = lsrkStep_AccessStepMem((ARKodeMem)step_mem->stepper->content,
                                      __func__, &lsrkstep_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to access LSRKStep memory for ExtSTS method.");
    return retval;
  }

  lsrkstep_mem->compute_err_est = enable;
  step_mem->extsts_dsm          = ZERO;
  step_mem->extsts_acc          = ZERO;
  step_mem->extsts_dt_acc       = ZERO;

  return ARK_SUCCESS;
}

/* Add the STS error estimate from one fast solve over [t0, tf].  The RKC/RKL
   embedding leaves the error vector in the inner tempv1; this is measured
   with the outer ewt so that the STS and MRI estimates share one norm. */
int extSTS_AccumulateStageError(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem,
                                sunrealtype t0, sunrealtype tf)
{
  const sunrealtype dt = SUNRabs(tf - t0);
  sunrealtype err;

  /* zero-length solves do not take an STS step, so tempv1 is not current */
  if (dt <= ZERO) { return ARK_SUCCESS; }

  err = N_VWrmsNorm(((ARKodeMem)step_mem->stepper->content)->tempv1,
                    ark_mem->ewt);
  step_mem->extsts_acc += dt * err * err;
  step_mem->extsts_dt_acc += dt;

  SUNLogDebug(ARK_LOGGER, "extsts-stage-error",
              "t0 = " SUN_FORMAT_G ", tf = " SUN_FORMAT_G
              ", err = " SUN_FORMAT_G,
              t0, tf, err);

  return ARK_SUCCESS;
}

/* Form the STS error estimate for the step attempt and combine it with the MRI
   embedding estimate stored in *dsm */
void extSTS_FinalizeStepError(SUNDIALS_MAYBE_UNUSED ARKodeMem ark_mem,
                              ARKodeMRIStepMem step_mem, sunrealtype* dsm)
{
  const sunrealtype dsm_slow = *dsm;

  step_mem->extsts_dsm =
    (step_mem->extsts_dt_acc > ZERO)
      ? SUNRsqrt(step_mem->extsts_acc / step_mem->extsts_dt_acc)
      : ZERO;
  *dsm =
    SUNRsqrt(dsm_slow * dsm_slow + step_mem->extsts_dsm * step_mem->extsts_dsm);

  SUNLogInfo(ARK_LOGGER, "extsts-error-estimate",
             "dsm_slow = " SUN_FORMAT_G ", dsm_sts = " SUN_FORMAT_G
             ", dsm = " SUN_FORMAT_G,
             dsm_slow, step_mem->extsts_dsm, *dsm);
}

/* Inner stepper utility routines */
int extSTSInnerStepper_Evolve(MRIStepInnerStepper sts_mem, sunrealtype t0,
                              sunrealtype tout, N_Vector y)
{
  /* Get the forcing data */
  ARKodeMem ark_mem              = (ARKodeMem)sts_mem->content;
  ARKodeLSRKStepMem lsrkstep_mem = NULL;
  sunrealtype tshift, tscale, dsm;
  N_Vector* forcing;
  const sunrealtype h = tout - t0;
  int nforcing, nflag, retval;

  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &lsrkstep_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  retval = MRIStepInnerStepper_GetForcingData(sts_mem, &tshift, &tscale,
                                              &forcing, &nforcing);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to retrieve forcing data for ExtSTS method.");
    return retval;
  }

  /* Reset LSRKStep to current state, using ARKodeReset so first-call
     setup invariants are preserved. */
  retval = ARKodeReset(ark_mem, t0, y);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to reset LSRKStep for ExtSTS method.");
    return retval;
  }

  /* Set step size to reach tout in a single fixed step */
  ark_mem->h = ark_mem->hin = h;
  ark_mem->fixedstep        = SUNTRUE;

  /* Run setup before enabling forcing since setup may allocate stepper data. */
  if (ark_mem->initsetup)
  {
    retval = arkInitialSetup(ark_mem, tout);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                      "Failed to initialize LSRKStep for ExtSTS method.");
      return retval;
    }
  }

  /* Set the inner forcing data. */
  retval = ark_mem->step_setforcing(ark_mem, tshift, tscale, forcing, nforcing);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to set LSRKStep forcing for ExtSTS method.");
    return retval;
  }

  /* The forcing wrapper changes the RHS, so recompute the start RHS. */
  ark_mem->fn_is_current = SUNFALSE;

  /* Call the user-supplied pre-step function (if supplied) */
  ark_mem->tcur = ark_mem->tn;
  if (ark_mem->ensure_ycur) { N_VScale(ONE, ark_mem->yn, ark_mem->ycur); }
  if (ark_mem->PreStepFn)
  {
    retval = ark_mem->PreStepFn(ark_mem->tcur, ark_mem->ycur, ark_mem->nst, 1,
                                ark_mem->user_data);
    if (retval != 0)
    {
      arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                      "Failure in user-supplied PreStepFn for ExtSTS method.");
      return retval;
    }
  }

  /* Take a single inner STS step */
  ark_mem->nst_attempts++;
  SUNLogInfo(ARK_LOGGER, "begin-step-attempt",
             "step = %li, tn = " SUN_FORMAT_G ", h = " SUN_FORMAT_G,
             ark_mem->nst + 1, ark_mem->tn, ark_mem->h);
  lsrkstep_mem->suppress_max_stage_limit_error = SUNTRUE;
  retval = ark_mem->step(ark_mem, &dsm, &nflag);
  lsrkstep_mem->suppress_max_stage_limit_error = SUNFALSE;

  if (retval != ARK_SUCCESS)
  {
    SUNLogInfo(ARK_LOGGER, "end-step-attempt",
               "status = failed step, kflag = %i", retval);
    if (retval == ARK_MAX_STAGE_LIMIT_FAIL) { return ARK_RETRY_STEP; }
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failure in LSRKStep step function for ExtSTS method.");
    return retval;
  }

  SUNLogInfo(ARK_LOGGER, "end-step-attempt",
             "status = success, dsm = " SUN_FORMAT_G, ZERO);

  /* Disable inner forcing */
  retval = ark_mem->step_setforcing(ark_mem, ZERO, ONE, NULL, 0);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to reset LSRKStep forcing for ExtSTS method.");
    return retval;
  }

  /* Complete successful steps to update stats and call the inner PostStepFn. */
  retval = arkCompleteStep(ark_mem, dsm);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to complete LSRKStep for ExtSTS method.");
  }

  return retval;
}

/*===============================================================
  EOF
  ===============================================================*/
