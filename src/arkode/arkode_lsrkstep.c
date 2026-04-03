/*---------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ UMBC
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
 * This is the implementation file for ARKODE's LSRK time stepper
 * module.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_math.h>

#include "arkode_impl.h"
#include "arkode_interp_impl.h"
#include "arkode_lsrkstep_impl.h"

/*===============================================================
  Exported functions
  ===============================================================*/

void* LSRKStepCreateSTS(ARKRhsFn rhs, sunrealtype t0, N_Vector y0,
                        SUNContext sunctx)
{
  ARKodeMem ark_mem;
  int retval;

  /* Create shared LSRKStep memory structure */
  ark_mem = lsrkStep_Create_Commons(rhs, t0, y0, sunctx);

  /* set default ARKODE_LSRK_RKC_2 method */
  retval = LSRKStepSetSTSMethod((void*)ark_mem, ARKODE_LSRK_RKC_2);
  if (retval != ARK_SUCCESS)
  {
    lsrkStep_Free(ark_mem);
    return NULL;
  }

  return (void*)ark_mem;
}

void* LSRKStepCreateSSP(ARKRhsFn rhs, sunrealtype t0, N_Vector y0,
                        SUNContext sunctx)
{
  ARKodeMem ark_mem;
  int retval;

  /* Create shared LSRKStep memory structure */
  ark_mem = lsrkStep_Create_Commons(rhs, t0, y0, sunctx);

  /* set default ARKODE_LSRK_SSP_S_2 method */
  retval = LSRKStepSetSSPMethod((void*)ark_mem, ARKODE_LSRK_SSP_S_2);
  if (retval != ARK_SUCCESS)
  {
    lsrkStep_Free(ark_mem);
    return NULL;
  }

  return (void*)ark_mem;
}

/*---------------------------------------------------------------
  LSRKStepReInitSTS:

  This routine re-initializes the LSRK STS module to solve a new
  problem of the same size as was previously solved. This routine
  should also be called when the problem dynamics or desired solvers
  have changed dramatically, so that the problem integration should
  resume as if started from scratch.

  Note all internal counters are set to 0 on re-initialization.
  ---------------------------------------------------------------*/
int LSRKStepReInitSTS(void* arkode_mem, ARKRhsFn rhs, sunrealtype t0, N_Vector y0)
{
  int retval;
  retval = lsrkStep_ReInit_Commons(arkode_mem, rhs, t0, y0);
  return retval;
}

/*---------------------------------------------------------------
  LSRKStepReInitSSP:

  This routine re-initializes the LSRK SSP module to solve a new
  problem of the same size as was previously solved. This routine
  should also be called when the problem dynamics or desired solvers
  have changed dramatically, so that the problem integration should
  resume as if started from scratch.

  Note all internal counters are set to 0 on re-initialization.
  ---------------------------------------------------------------*/
int LSRKStepReInitSSP(void* arkode_mem, ARKRhsFn rhs, sunrealtype t0, N_Vector y0)
{
  int retval;
  retval = lsrkStep_ReInit_Commons(arkode_mem, rhs, t0, y0);
  return retval;
}

/*===============================================================
  Interface routines supplied to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  lsrkStep_Create_Commons:

  A submodule for creating the common features of
  LSRKStepCreateSTS and LSRKStepCreateSSP.
  ---------------------------------------------------------------*/

void* lsrkStep_Create_Commons(ARKRhsFn rhs, sunrealtype t0, N_Vector y0,
                              SUNContext sunctx)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* Check that rhs is supplied */
  if (rhs == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_F);
    return NULL;
  }

  /* Check for legal input parameters */
  if (y0 == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_Y0);
    return NULL;
  }

  if (sunctx == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_SUNCTX);
    return NULL;
  }

  /* Create ark_mem structure and set default values */
  ark_mem = arkCreate(sunctx);
  if (ark_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return NULL;
  }

  /* Allocate ARKodeLSRKStepMem structure, and initialize to zero */
  step_mem = (ARKodeLSRKStepMem)calloc(1, sizeof(*step_mem));
  if (step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }

  /* Attach step_mem structure and function pointers to ark_mem */
  ark_mem->step_init              = lsrkStep_Init;
  ark_mem->step_fullrhs           = lsrkStep_FullRHS;
  ark_mem->step                   = lsrkStep_TakeStepRKC;
  ark_mem->step_printallstats     = lsrkStep_PrintAllStats;
  ark_mem->step_writeparameters   = lsrkStep_WriteParameters;
  ark_mem->step_free              = lsrkStep_Free;
  ark_mem->step_printmem          = lsrkStep_PrintMem;
  ark_mem->step_setoptions        = lsrkStep_SetOptions;
  ark_mem->step_setdefaults       = lsrkStep_SetDefaults;
  ark_mem->step_getnumrhsevals    = lsrkStep_GetNumRhsEvals;
  ark_mem->step_getestlocalerrors = lsrkStep_GetEstLocalErrors;
  ark_mem->step_getstageindex     = lsrkStep_GetStageIndex;
  ark_mem->step_mem               = (void*)step_mem;
  ark_mem->step_supports_adaptive = SUNTRUE;

  /* Set default values for optional inputs */
  retval = lsrkStep_SetDefaults((void*)ark_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Error setting default solver options");
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }

  /* Copy the input parameters into ARKODE state */
  step_mem->fe = rhs;

  /* Initialize spectral radius info */
  step_mem->lambdaR             = ZERO;
  step_mem->lambdaI             = ZERO;
  step_mem->spectral_radius     = ZERO;
  step_mem->spectral_radius_max = ZERO;
  step_mem->spectral_radius_min = ZERO;

  /* Initialize flags */
  step_mem->dom_eig_update     = SUNTRUE;
  step_mem->dom_eig_is_current = SUNFALSE;
  step_mem->is_SSP             = SUNFALSE;
  step_mem->init_warmup        = SUNTRUE;

  /* Set NULL for dom_eig_fn */
  step_mem->dom_eig_fn = NULL;

  /* Set NULL for DEE */
  step_mem->DEE = NULL;

  /* Initialize all the counters */
  step_mem->nfe               = 0;
  step_mem->nfeDQ             = 0;
  step_mem->stage_max         = 0;
  step_mem->dom_eig_num_evals = 0;
  step_mem->stage_max_limit   = STAGE_MAX_LIMIT_DEFAULT;
  step_mem->dom_eig_nst       = 0;
  step_mem->num_dee_iters     = 0;

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to initialize main ARKODE infrastructure");
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }

  /* Specify preferred interpolation type */
  ARKodeSetInterpolantType(ark_mem, ARK_INTERP_LAGRANGE);

  return (void*)ark_mem;
}

/*---------------------------------------------------------------
  lsrkStep_ReInit_Commons:

  A submodule designed to reinitialize the common features of
  LSRKStepCreateSTS and LSRKStepCreateSSP.
  ---------------------------------------------------------------*/

int lsrkStep_ReInit_Commons(void* arkode_mem, ARKRhsFn rhs, sunrealtype t0,
                            N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE)
  {
    arkProcessError(ark_mem, ARK_NO_MALLOC, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MALLOC);
    return ARK_NO_MALLOC;
  }

  /* Check that rhs is supplied */
  if (rhs == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_F);
    return ARK_ILL_INPUT;
  }

  /* Check for legal input parameters */
  if (y0 == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_Y0);
    return ARK_ILL_INPUT;
  }

  /* Copy the input parameters into ARKODE state */
  step_mem->fe = rhs;

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(arkode_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to initialize main ARKODE infrastructure");
    return retval;
  }

  /* Initialize all the counters, flags and stats */
  step_mem->nfe                 = 0;
  step_mem->nfeDQ               = 0;
  step_mem->dom_eig_num_evals   = 0;
  step_mem->stage_max           = 0;
  step_mem->lambdaR             = ZERO;
  step_mem->lambdaI             = ZERO;
  step_mem->spectral_radius     = ZERO;
  step_mem->spectral_radius_max = 0;
  step_mem->spectral_radius_min = 0;
  step_mem->dom_eig_nst         = 0;
  step_mem->num_dee_iters       = 0;
  step_mem->dom_eig_update      = SUNTRUE;
  step_mem->dom_eig_is_current  = SUNFALSE;
  step_mem->init_warmup         = SUNTRUE;

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  lsrkStep_Init:

  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup.

  With initialization types FIRST_INIT this routine performs
  setup and allocations needed for the method and sets
  the call_fullrhs flag.

  With other initialization types, this routine does nothing.
  ---------------------------------------------------------------*/
int lsrkStep_Init(ARKodeMem ark_mem, int init_type)
{
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  /* immediately return if resize or reset */
  if (init_type == RESIZE_INIT || init_type == RESET_INIT)
  {
    return ARK_SUCCESS;
  }
  /* enforce use of arkEwtSmallReal if using a fixed step size
     and an internal error weight function */
  if (ark_mem->fixedstep && !ark_mem->user_efun)
  {
    ark_mem->user_efun = SUNFALSE;
    ark_mem->efun      = arkEwtSetSmallReal;
    ark_mem->e_data    = ark_mem;
  }

  /* Check if user has provided dom_eig_fn or DEE */
  if (!step_mem->is_SSP && step_mem->dom_eig_fn == NULL && step_mem->DEE == NULL)
  {
    arkProcessError(ark_mem, ARK_DOMEIG_FAIL, __LINE__, __func__, __FILE__,
                    "STS methods require either a user provided dominant "
                    "eigenvalue function or a SUNDomEigEstimator");
    return ARK_DOMEIG_FAIL;
  }

  /* Initialize the DEE */
  if (step_mem->DEE != NULL)
  {
    retval = SUNDomEigEstimator_Initialize(step_mem->DEE);
    if (retval != SUN_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_DEE_FAIL, __LINE__, __func__, __FILE__,
                      "SUNDomEigEstimator_Initialize failed");
      return ARK_DEE_FAIL;
    }

    /* Set number of DEE preprocessing iterations for the initial estimate */
    retval = SUNDomEigEstimator_SetNumPreprocessIters(step_mem->DEE,
                                                      step_mem->num_init_warmups);
    if (retval != SUN_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_DEE_FAIL, __LINE__, __func__, __FILE__,
                      "SUNDomEigEstimator_SetNumPreprocessIters failed");
      return ARK_DEE_FAIL;
    }
  }

  /* Allocate reusable arrays for fused vector interface */
  if (step_mem->cvals == NULL)
  {
    step_mem->cvals = (sunrealtype*)calloc(step_mem->nfusedopvecs,
                                           sizeof(sunrealtype));
    if (step_mem->cvals == NULL) { return ARK_MEM_FAIL; }
    ark_mem->lrw += step_mem->nfusedopvecs;
  }
  if (step_mem->Xvecs == NULL)
  {
    step_mem->Xvecs = (N_Vector*)calloc(step_mem->nfusedopvecs, sizeof(N_Vector));
    if (step_mem->Xvecs == NULL) { return ARK_MEM_FAIL; }
    ark_mem->liw += step_mem->nfusedopvecs; /* pointers */
  }

  /* While LSRKStep does not currently call the full RHS function directly (later
     optimizations might) we do need the fn vector to always be allocated. Signaling
     to shared arkode module that full RHS evaluations are required will ensure
     fn is always allocated. */
  ark_mem->call_fullrhs = SUNTRUE;

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  lsrkStep_FullRHS:

  This is just a wrapper to call the user-supplied RHS function, f(t,y).

  This will be called in one of three 'modes':

     ARK_FULLRHS_START -> called at the beginning of a simulation i.e., at
                          (tn, yn) = (t0, y0) or (tR, yR)

     ARK_FULLRHS_END   -> called at the end of a successful step i.e, at
                          (tcur, ycur) or the start of the subsequent step i.e.,
                          at (tn, yn) = (tcur, ycur) from the end of the last
                          step

     ARK_FULLRHS_OTHER -> called elsewhere (e.g. for dense output)

  If this function is called in ARK_FULLRHS_START the RHS function is always
  evaluated.

  In ARK_FULLRHS_END mode we evaluate the RHS if an SSP method is being
  used otherwise we copy the RHS evaluation from the end of the STS step.

  ARK_FULLRHS_OTHER mode is only called for dense output in-between steps, or
  when estimating the initial time step size, so we strive to store the
  intermediate parts so that they do not interfere with the other two modes.
  ----------------------------------------------------------------------------*/
int lsrkStep_FullRHS(ARKodeMem ark_mem, sunrealtype t, N_Vector y, N_Vector f,
                     int mode)
{
  int retval;
  ARKodeLSRKStepMem step_mem;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  /* perform RHS functions contingent on 'mode' argument */
  switch (mode)
  {
  case ARK_FULLRHS_START:

    /* compute the RHS */
    if (!ark_mem->fn_is_current)
    {
      /* call the user-supplied pre-rhs function (if supplied) */
      if (ark_mem->PreRhsFn)
      {
        retval = ark_mem->PreRhsFn(t, y, ark_mem->user_data);
        if (retval != 0) { return (ARK_PRERHSFN_FAIL); }
      }
      retval = step_mem->fe(t, y, f, ark_mem->user_data);
      step_mem->nfe++;
      if (retval != 0)
      {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                        MSG_ARK_RHSFUNC_FAILED, t);
        return ARK_RHSFUNC_FAIL;
      }
    }

    break;

  case ARK_FULLRHS_END:
    /* No further action is needed if STS since the currently available STS methods
       evaluate the RHS at the end of each time step. If the stepper is an SSP, fn is
       updated and reused at the beginning of the step unless
       ark_mem->fn_is_current is changed by ARKODE. */
    if (step_mem->is_SSP)
    {
      /* call the user-supplied pre-rhs function (if supplied) */
      if (ark_mem->PreRhsFn)
      {
        retval = ark_mem->PreRhsFn(t, y, ark_mem->user_data);
        if (retval != 0) { return (ARK_PRERHSFN_FAIL); }
      }
      retval = step_mem->fe(t, y, ark_mem->fn, ark_mem->user_data);
      step_mem->nfe++;
      if (retval != 0)
      {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                        MSG_ARK_RHSFUNC_FAILED, t);
        return ARK_RHSFUNC_FAIL;
      }
      ark_mem->fn_is_current = SUNTRUE;
    }
    N_VScale(ONE, ark_mem->fn, f);

    break;

  case ARK_FULLRHS_OTHER:

    /* call the user-supplied pre-rhs function (if supplied) */
    if (ark_mem->PreRhsFn)
    {
      retval = ark_mem->PreRhsFn(t, y, ark_mem->user_data);
      if (retval != 0) { return (ARK_PRERHSFN_FAIL); }
    }

    /* call f */
    retval = step_mem->fe(t, y, f, ark_mem->user_data);
    step_mem->nfe++;
    if (retval != 0)
    {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_RHSFUNC_FAILED, t);
      return ARK_RHSFUNC_FAIL;
    }

    break;

  default:
    /* return with RHS failure if unknown mode is passed */
    arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                    "Unknown full RHS mode");
    return ARK_RHSFUNC_FAIL;
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  lsrkStep_TakeStepRKC:

  This routine serves the primary purpose of the LSRKStepRKC module:
  it performs a single RKC step (with embedding, if possible).

  The output variable dsmPtr should contain estimate of the
  weighted local error if an embedding is present; otherwise it
  should be 0.

  The variables (ark_mem->tcur, ark_mem->ycur) should
  contain the current time and solution at the end of this time step.

  The input/output variable nflagPtr is generally used in ARKODE
  to gauge the convergence of any algebraic solvers. However, since
  the STS step routines do not involve an algebraic solve, this variable
  instead serves to identify possible ARK_RETRY_STEP returns within this
  routine.

  The return value from this routine is:
            0 => step completed successfully
           >0 => step encountered recoverable failure;
                 reduce step and retry (if possible)
                 ARK_RETRY_STEP indicates that the required stage
                 number has reached the stage_max_limit with the
                 current value of h. The step is then returned to
                 adjust the step size.
           <0 => step encountered unrecoverable failure
  ---------------------------------------------------------------*/
int lsrkStep_TakeStepRKC(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  int retval;
  sunrealtype hmax, w0, w1, temp1, temp2, arg, bjm1, bjm2, mus, thjm1, thjm2,
    zjm1, zjm2, dzjm1, dzjm2, d2zjm1, d2zjm2, zj, dzj, d2zj, bj, ajm1, mu, nu,
    thj;
  const sunrealtype onep54 = SUN_RCONST(1.54), c13 = SUN_RCONST(13.0),
                    p8 = SUN_RCONST(0.8), p4 = SUN_RCONST(0.4);
  ARKodeLSRKStepMem step_mem;

  /* initialize algebraic solver convergence flag to success,
     temporal error estimate to zero */
  *nflagPtr = ARK_SUCCESS;
  *dsmPtr   = ZERO;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  sunrealtype* cvals = step_mem->cvals;
  N_Vector* Xvecs    = step_mem->Xvecs;
  N_Vector tmp1      = ark_mem->tempv1;
  N_Vector tmp2      = ark_mem->tempv2;

  /* Initialize the current stage index */
  step_mem->istage     = 0;
  step_mem->req_stages = step_mem->stage_max_limit;

  /* Compute dominant eigenvalue and update stats */
  if (step_mem->dom_eig_update)
  {
    retval = lsrkStep_ComputeNewDomEig(ark_mem, step_mem);
    if (retval != ARK_SUCCESS) { return retval; }
  }

  sunrealtype ss =
    SUNRceil(SUNRsqrt(onep54 * SUNRabs(ark_mem->h) * step_mem->spectral_radius));
  ss = SUNMAX(ss, SUN_RCONST(2.0));

  if (ss >= step_mem->stage_max_limit)
  {
    SUNLogInfo(ARK_LOGGER, "compute-num-stages",
               "spectral radius = " SUN_FORMAT_G ", num stages = " SUN_FORMAT_G
               ", max stages = %i, max stage limit = %i",
               step_mem->spectral_radius, ss, step_mem->stage_max,
               step_mem->stage_max_limit);

    if (!ark_mem->fixedstep)
    {
      hmax = ark_mem->hadapt_mem->safety * SUNSQR(step_mem->stage_max_limit) /
             (onep54 * step_mem->spectral_radius);
      ark_mem->eta = hmax / ark_mem->h;
      *nflagPtr    = ARK_RETRY_STEP;
      ark_mem->hadapt_mem->nst_exp++;
      return ARK_RETRY_STEP;
    }
    else
    {
      arkProcessError(ark_mem, ARK_MAX_STAGE_LIMIT_FAIL, __LINE__, __func__,
                      __FILE__,
                      "Unable to achieve stable results: Either reduce the "
                      "step size or increase the stage_max_limit");
      return ARK_MAX_STAGE_LIMIT_FAIL;
    }
  }

  step_mem->req_stages = (int)ss;
  step_mem->stage_max  = SUNMAX(step_mem->req_stages, step_mem->stage_max);

  SUNLogInfo(ARK_LOGGER, "compute-num-stages",
             "spectral radius = " SUN_FORMAT_G
             ", num stages = %i, max stages = %i, max stage limit = %i",
             step_mem->spectral_radius, step_mem->req_stages,
             step_mem->stage_max, step_mem->stage_max_limit);
  SUNLogInfo(ARK_LOGGER, "begin-stages-list",
             "stage = %i, tcur = " SUN_FORMAT_G, 0, ark_mem->tcur);
  SUNLogExtraDebugVec(ARK_LOGGER, "stage", ark_mem->yn, "z_0(:) =");

  /* Compute RHS function for the start of the step, if necessary. */
  if ((!ark_mem->fn_is_current && ark_mem->initsetup) ||
      (step_mem->step_nst != ark_mem->nst))
  {
    /* call the user-supplied pre-rhs function (if supplied) */
    if (ark_mem->PreRhsFn)
    {
      retval = ark_mem->PreRhsFn(ark_mem->tn, ark_mem->yn, ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed preprocess rhs, retval = %i", retval);
        return ARK_PRERHSFN_FAIL;
      }
    }

    /* call fe */
    retval = step_mem->fe(ark_mem->tn, ark_mem->yn, ark_mem->fn,
                          ark_mem->user_data);
    step_mem->nfe++;
    if (retval != ARK_SUCCESS)
    {
      SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->fn, "F_0(:) =");
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed rhs eval, retval = %i", retval);
      return (ARK_RHSFUNC_FAIL);
    }
    ark_mem->fn_is_current = SUNTRUE;
  }

  SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->fn, "F_0(:) =");
  SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

  /* Track the number of successful steps to determine if the previous step failed. */
  step_mem->step_nst = ark_mem->nst + 1;

  /* Initialize constants */
  w0    = (ONE + TWO / (c13 * SUNSQR((sunrealtype)(step_mem->req_stages))));
  temp1 = SUNSQR(w0) - ONE;
  temp2 = SUNRsqrt(temp1);
  arg   = step_mem->req_stages * SUNRlog(w0 + temp2);
  w1    = SUNRsinh(arg) * temp1 /
       (SUNRcosh(arg) * step_mem->req_stages * temp2 - w0 * SUNRsinh(arg));
  bjm1 = ONE / SUNSQR(TWO * w0);
  bjm2 = bjm1;
  mus  = w1 * bjm1;

  /* Begin stage 1 (store in tmp2) and initialize embedding */
  ark_mem->tcur    = ark_mem->tn + ark_mem->h * mus;
  step_mem->istage = 1;
  SUNLogInfo(ARK_LOGGER, "begin-stages-list",
             "stage = %i, tcur = " SUN_FORMAT_G, 1, ark_mem->tcur);
  N_VLinearSum(ONE, ark_mem->yn, ark_mem->h * mus, ark_mem->fn, tmp2);
  N_VScale(ONE, ark_mem->yn, tmp1);

  /* apply user-supplied stage postprocessing function (if supplied) */
  if (ark_mem->PostProcessStageFn)
  {
    retval = ark_mem->PostProcessStageFn(ark_mem->tcur, tmp2, ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed postprocess stage, retval = %i", retval);
      return ARK_POSTPROCESS_STAGE_FAIL;
    }
  }

  /* Initialize constants for stage loop */
  thjm2  = ZERO;
  thjm1  = mus;
  zjm1   = w0;
  zjm2   = ONE;
  dzjm1  = ONE;
  dzjm2  = ZERO;
  d2zjm1 = ZERO;
  d2zjm2 = ZERO;

  /* Evaluate stages j = 2,...,step_mem->req_stages */
  for (int j = 2; j <= step_mem->req_stages; j++)
  {
    /* Complete the previous stage (evaluate the RHS and store it in ycur) */

    /* call the user-supplied pre-RHS function (if supplied) */
    if (ark_mem->PreRhsFn)
    {
      retval = ark_mem->PreRhsFn(ark_mem->tcur, tmp2, ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed preprocess rhs, retval = %i", retval);
        return ARK_PRERHSFN_FAIL;
      }
    }

    retval = step_mem->fe(ark_mem->tcur, tmp2, ark_mem->ycur, ark_mem->user_data);
    step_mem->nfe++;

    SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->ycur,
                        "F_%i(:) =", j - 1);
    SUNLogInfoIf(retval != 0, ARK_LOGGER, "end-stages-list",
                 "status = failed rhs eval, retval = %i", retval);

    if (retval < 0) { return ARK_RHSFUNC_FAIL; }
    if (retval > 0) { return RHSFUNC_RECVR; }

    SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

    /* Begin stage j (store in ycur) */
    zj               = TWO * w0 * zjm1 - zjm2;
    dzj              = TWO * w0 * dzjm1 - dzjm2 + TWO * zjm1;
    d2zj             = TWO * w0 * d2zjm1 - d2zjm2 + FOUR * dzjm1;
    bj               = d2zj / SUNSQR(dzj);
    ajm1             = ONE - zjm1 * bjm1;
    mu               = TWO * w0 * bj / bjm1;
    nu               = -bj / bjm2;
    mus              = mu * w1 / w0;
    thj              = mu * thjm1 + nu * thjm2 + mus * (ONE - ajm1);
    ark_mem->tcur    = ark_mem->tn + ark_mem->h * thj;
    step_mem->istage = j;
    SUNLogInfo(ARK_LOGGER, "begin-stages-list",
               "stage = %i, tcur = " SUN_FORMAT_G, j, ark_mem->tcur);
    cvals[0] = mus * ark_mem->h;
    Xvecs[0] = ark_mem->ycur;
    cvals[1] = nu;
    Xvecs[1] = tmp1;
    cvals[2] = ONE - mu - nu;
    Xvecs[2] = ark_mem->yn;
    cvals[3] = mu;
    Xvecs[3] = tmp2;
    cvals[4] = -mus * ajm1 * ark_mem->h;
    Xvecs[4] = ark_mem->fn;
    retval   = N_VLinearCombination(5, cvals, Xvecs, ark_mem->ycur);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed vector op, retval = %i", retval);
      return ARK_VECTOROP_ERR;
    }

    /* apply user-supplied stage or step postprocessing function (if supplied) */
    if (j < step_mem->req_stages && ark_mem->PostProcessStageFn)
    {
      retval = ark_mem->PostProcessStageFn(ark_mem->tcur, ark_mem->ycur,
                                           ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed postprocess stage, retval = %i", retval);
        return ARK_POSTPROCESS_STAGE_FAIL;
      }
    }
    else if (j == step_mem->req_stages && ark_mem->PostProcessStepFn)
    {
      retval = ark_mem->PostProcessStepFn(ark_mem->tcur, ark_mem->ycur,
                                          ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed postprocess step, retval = %i", retval);
        return ARK_POSTPROCESS_STEP_FAIL;
      }
    }

    /* Shift the data for the next stage */
    if (j < step_mem->req_stages)
    {
      /* Swap tempv1 and tempv2 pointers to handle two-previous-stage logic */
      N_Vector temp = tmp1;
      tmp1          = tmp2;
      tmp2          = temp;

      N_VScale(ONE, ark_mem->ycur, tmp2);

      /* Update coefficients to handle the two-previous stage logic */
      thjm2  = thjm1;
      thjm1  = thj;
      bjm2   = bjm1;
      bjm1   = bj;
      zjm2   = zjm1;
      zjm1   = zj;
      dzjm2  = dzjm1;
      dzjm1  = dzj;
      d2zjm2 = d2zjm1;
      d2zjm1 = d2zj;
    }
  }

  SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");
  SUNLogExtraDebugVec(ARK_LOGGER, "updated solution", ark_mem->ycur, "ycur(:) =");
  SUNLogInfo(ARK_LOGGER, "begin-compute-embedding", "");

  /* final stage processing */
  ark_mem->tcur = ark_mem->tn + ark_mem->h;

  /* call the user-supplied pre-RHS function (if supplied) */
  if (ark_mem->PreRhsFn)
  {
    retval = ark_mem->PreRhsFn(ark_mem->tcur, ark_mem->ycur, ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-compute-embedding",
                 "status = failed preprocess rhs, retval = %i", retval);
      return ARK_PRERHSFN_FAIL;
    }
  }

  retval = step_mem->fe(ark_mem->tcur, ark_mem->ycur, ark_mem->tempv2,
                        ark_mem->user_data);
  step_mem->nfe++;

  SUNLogExtraDebugVec(ARK_LOGGER, "solution RHS", ark_mem->tempv2, "F_n(:) =");
  SUNLogInfoIf(retval != 0, ARK_LOGGER, "end-compute-embedding",
               "status = failed rhs eval, retval = %i", retval);

  if (retval < 0) { return ARK_RHSFUNC_FAIL; }
  if (retval > 0) { return RHSFUNC_RECVR; }

  /* Compute yerr (if step adaptivity enabled) */
  if (!ark_mem->fixedstep)
  {
    /* Estimate the local error and compute its weighted RMS norm */
    cvals[0] = p8;
    Xvecs[0] = ark_mem->yn;
    cvals[1] = -p8;
    Xvecs[1] = ark_mem->ycur;
    cvals[2] = p4 * ark_mem->h;
    Xvecs[2] = ark_mem->fn;
    cvals[3] = p4 * ark_mem->h;
    Xvecs[3] = ark_mem->tempv2;

    retval = N_VLinearCombination(4, cvals, Xvecs, ark_mem->tempv1);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-compute-embedding",
                 "status = failed vector op, retval = %i", retval);
      return ARK_VECTOROP_ERR;
    }
    *dsmPtr = N_VWrmsNorm(ark_mem->tempv1, ark_mem->ewt);
  }
  lsrkStep_DomEigUpdateLogic(ark_mem, step_mem, *dsmPtr, ark_mem->tempv2);

  SUNLogInfo(ARK_LOGGER, "end-compute-embedding", "status = success");

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  lsrkStep_TakeStepRKL:

  This routine serves the primary purpose of the LSRKStepRKL module:
  it performs a single RKL step (with embedding, if possible).

  The output variable dsmPtr should contain estimate of the
  weighted local error if an embedding is present; otherwise it
  should be 0.

  The variables (ark_mem->tcur, ark_mem->ycur) should
  contain the current time and solution at the end of this time step.

  The input/output variable nflagPtr is generally used in ARKODE
  to gauge the convergence of any algebraic solvers. However, since
  the STS step routines do not involve an algebraic solve, this variable
  instead serves to identify possible ARK_RETRY_STEP returns within this
  routine.

  The return value from this routine is:
            0 => step completed successfully
           >0 => step encountered recoverable failure;
                 reduce step and retry (if possible)
                 ARK_RETRY_STEP indicates that the required stage
                 number has reached the stage_max_limit with the
                 current value of h. The step is then returned to
                 adjust the step size.
           <0 => step encountered unrecoverable failure
  ---------------------------------------------------------------*/
int lsrkStep_TakeStepRKL(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  int retval;
  sunrealtype hmax, w1, bjm1, bjm2, mus, bj, ajm1, temj, cj, mu, nu;
  const sunrealtype p8 = SUN_RCONST(0.8), p4 = SUN_RCONST(0.4);
  ARKodeLSRKStepMem step_mem;

  /* initialize algebraic solver convergence flag to success,
     temporal error estimate to zero */
  *nflagPtr = ARK_SUCCESS;
  *dsmPtr   = ZERO;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  sunrealtype* cvals = step_mem->cvals;
  N_Vector* Xvecs    = step_mem->Xvecs;
  N_Vector tmp1      = ark_mem->tempv1;
  N_Vector tmp2      = ark_mem->tempv2;

  /* Initialize the current stage index */
  step_mem->istage     = 0;
  step_mem->req_stages = step_mem->stage_max_limit;

  /* Compute dominant eigenvalue and update stats */
  if (step_mem->dom_eig_update)
  {
    retval = lsrkStep_ComputeNewDomEig(ark_mem, step_mem);
    if (retval != ARK_SUCCESS) { return retval; }
  }

  sunrealtype ss =
    SUNRceil((SUNRsqrt(SUN_RCONST(9.0) + SUN_RCONST(8.0) * SUNRabs(ark_mem->h) *
                                           step_mem->spectral_radius) -
              ONE) /
             TWO);

  ss = SUNMAX(ss, SUN_RCONST(2.0));

  if (ss >= step_mem->stage_max_limit)
  {
    SUNLogInfo(ARK_LOGGER, "compute-num-stages",
               "spectral radius = " SUN_FORMAT_G ", num stages = " SUN_FORMAT_G
               ", max stages = %i, max stage limit = %i",
               step_mem->spectral_radius, ss, step_mem->stage_max,
               step_mem->stage_max_limit);

    if (!ark_mem->fixedstep)
    {
      hmax =
        ark_mem->hadapt_mem->safety *
        (SUNSQR(step_mem->stage_max_limit) + step_mem->stage_max_limit - TWO) /
        (TWO * step_mem->spectral_radius);
      ark_mem->eta = hmax / ark_mem->h;
      *nflagPtr    = ARK_RETRY_STEP;
      ark_mem->hadapt_mem->nst_exp++;
      return ARK_RETRY_STEP;
    }
    else
    {
      arkProcessError(ark_mem, ARK_MAX_STAGE_LIMIT_FAIL, __LINE__, __func__,
                      __FILE__,
                      "Unable to achieve stable results: Either reduce the "
                      "step size or increase the stage_max_limit");
      return ARK_MAX_STAGE_LIMIT_FAIL;
    }
  }

  step_mem->req_stages = (int)ss;
  step_mem->stage_max  = SUNMAX(step_mem->req_stages, step_mem->stage_max);

  SUNLogInfo(ARK_LOGGER, "compute-num-stages",
             "spectral radius = " SUN_FORMAT_G
             ", num stages = %i, max stages = %i, max stage limit = %i",
             step_mem->spectral_radius, step_mem->req_stages,
             step_mem->stage_max, step_mem->stage_max_limit);
  SUNLogInfo(ARK_LOGGER, "begin-stages-list",
             "stage = %i, tcur = " SUN_FORMAT_G, 0, ark_mem->tcur);
  SUNLogExtraDebugVec(ARK_LOGGER, "stage", ark_mem->yn, "z_0(:) =");

  /* Compute RHS function for the start of the step, if necessary. */
  if ((!ark_mem->fn_is_current && ark_mem->initsetup) ||
      (step_mem->step_nst != ark_mem->nst))
  {
    /* call the user-supplied pre-RHS function (if supplied) */
    if (ark_mem->PreRhsFn)
    {
      retval = ark_mem->PreRhsFn(ark_mem->tn, ark_mem->yn, ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed preprocess rhs, retval = %i", retval);
        return ARK_PRERHSFN_FAIL;
      }
    }
    retval = step_mem->fe(ark_mem->tn, ark_mem->yn, ark_mem->fn,
                          ark_mem->user_data);
    step_mem->nfe++;
    if (retval != ARK_SUCCESS)
    {
      SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->fn, "F_0(:) =");
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed rhs eval, retval = %i", retval);
      return (ARK_RHSFUNC_FAIL);
    }
    ark_mem->fn_is_current = SUNTRUE;
  }

  SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->fn, "F_0(:) =");
  SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

  /* Track the number of successful steps to determine if the previous step failed. */
  step_mem->step_nst = ark_mem->nst + 1;

  /* Initialize constants */
  w1   = FOUR / ((step_mem->req_stages + TWO) * (step_mem->req_stages - ONE));
  bjm2 = ONE / THREE;
  bjm1 = bjm2;
  mus  = w1 * bjm1;

  /* Begin stage 1 (store in tmp2) and initialize embedding */
  ark_mem->tcur    = ark_mem->tn + ark_mem->h * mus;
  step_mem->istage = 1;
  SUNLogInfo(ARK_LOGGER, "begin-stages-list",
             "stage = %i, tcur = " SUN_FORMAT_G, 1, ark_mem->tcur);
  N_VLinearSum(ONE, ark_mem->yn, ark_mem->h * mus, ark_mem->fn, tmp2);
  N_VScale(ONE, ark_mem->yn, tmp1);

  /* apply user-supplied stage postprocessing function (if supplied) */
  if (ark_mem->PostProcessStageFn)
  {
    retval = ark_mem->PostProcessStageFn(ark_mem->tcur, tmp2, ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed postprocess stage, retval = %i", retval);
      return ARK_POSTPROCESS_STAGE_FAIL;
    }
  }

  /* Evaluate stages j = 2,...,step_mem->req_stages */
  for (int j = 2; j <= step_mem->req_stages; j++)
  {
    /* Complete the previous stage (evaluate the RHS and store it in ycur) */

    /* call the user-supplied pre-RHS function (if supplied) */
    if (ark_mem->PreRhsFn)
    {
      retval = ark_mem->PreRhsFn(ark_mem->tcur, tmp2, ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed preprocess rhs, retval = %i", retval);
        return ARK_PRERHSFN_FAIL;
      }
    }

    retval = step_mem->fe(ark_mem->tcur, tmp2, ark_mem->ycur, ark_mem->user_data);
    step_mem->nfe++;

    SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->ycur,
                        "F_%i(:) =", j - 1);
    SUNLogInfoIf(retval != 0, ARK_LOGGER, "end-stages-list",
                 "status = failed rhs eval, retval = %i", retval);

    if (retval < 0) { return ARK_RHSFUNC_FAIL; }
    if (retval > 0) { return RHSFUNC_RECVR; }

    SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

    /* Begin stage j (store in ycur) */
    temj             = (j + TWO) * (j - ONE);
    bj               = temj / (TWO * j * (j + ONE));
    ajm1             = ONE - bjm1;
    mu               = (TWO * j - ONE) / j * (bj / bjm1);
    nu               = -(j - ONE) / j * (bj / bjm2);
    mus              = w1 * mu;
    cj               = temj * w1 / FOUR;
    ark_mem->tcur    = ark_mem->tn + ark_mem->h * cj;
    step_mem->istage = j;
    SUNLogInfo(ARK_LOGGER, "begin-stages-list",
               "stage = %i, tcur = " SUN_FORMAT_G, j, ark_mem->tcur);
    cvals[0] = mus * ark_mem->h;
    Xvecs[0] = ark_mem->ycur;
    cvals[1] = nu;
    Xvecs[1] = tmp1;
    cvals[2] = ONE - mu - nu;
    Xvecs[2] = ark_mem->yn;
    cvals[3] = mu;
    Xvecs[3] = tmp2;
    cvals[4] = -mus * ajm1 * ark_mem->h;
    Xvecs[4] = ark_mem->fn;
    retval   = N_VLinearCombination(5, cvals, Xvecs, ark_mem->ycur);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed vector op, retval = %i", retval);
      return ARK_VECTOROP_ERR;
    }

    /* apply user-supplied stage or step postprocessing function (if supplied) */
    if (j < step_mem->req_stages && ark_mem->PostProcessStageFn)
    {
      retval = ark_mem->PostProcessStageFn(ark_mem->tcur, ark_mem->ycur,
                                           ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed postprocess stage, retval = %i", retval);
        return ARK_POSTPROCESS_STAGE_FAIL;
      }
    }
    else if (j == step_mem->req_stages && ark_mem->PostProcessStepFn)
    {
      retval = ark_mem->PostProcessStepFn(ark_mem->tcur + ark_mem->h * cj,
                                          ark_mem->ycur, ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed postprocess step, retval = %i", retval);
        return ARK_POSTPROCESS_STEP_FAIL;
      }
    }

    /* Shift the data for the next stage */
    if (j < step_mem->req_stages)
    {
      /* To avoid two data copies we swap ARKODE's tempv1 and tempv2 pointers*/
      N_Vector temp = tmp1;
      tmp1          = tmp2;
      tmp2          = temp;

      N_VScale(ONE, ark_mem->ycur, tmp2);

      bjm2 = bjm1;
      bjm1 = bj;
    }
  }

  SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");
  SUNLogExtraDebugVec(ARK_LOGGER, "updated solution", ark_mem->ycur, "ycur(:) =");

  /* final stage processing */
  SUNLogInfo(ARK_LOGGER, "begin-compute-embedding", "");
  ark_mem->tcur = ark_mem->tn + ark_mem->h;

  /* call the user-supplied pre-RHS function (if supplied) */
  if (ark_mem->PreRhsFn)
  {
    retval = ark_mem->PreRhsFn(ark_mem->tcur, ark_mem->ycur, ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-compute-embedding",
                 "status = failed preprocess rhs, retval = %i", retval);
      return ARK_PRERHSFN_FAIL;
    }
  }
  retval = step_mem->fe(ark_mem->tcur, ark_mem->ycur, ark_mem->tempv2,
                        ark_mem->user_data);
  step_mem->nfe++;

  SUNLogExtraDebugVec(ARK_LOGGER, "solution RHS", ark_mem->tempv2, "F_n(:) =");
  SUNLogInfoIf(retval != 0, ARK_LOGGER, "end-compute-embedding",
               "status = failed rhs eval, retval = %i", retval);

  if (retval < 0) { return ARK_RHSFUNC_FAIL; }
  if (retval > 0) { return RHSFUNC_RECVR; }

  /* Compute yerr (if step adaptivity enabled) */
  if (!ark_mem->fixedstep)
  {
    /* Estimate the local error and compute its weighted RMS norm */
    cvals[0] = p8;
    Xvecs[0] = ark_mem->yn;
    cvals[1] = -p8;
    Xvecs[1] = ark_mem->ycur;
    cvals[2] = p4 * ark_mem->h;
    Xvecs[2] = ark_mem->fn;
    cvals[3] = p4 * ark_mem->h;
    Xvecs[3] = ark_mem->tempv2;
    retval   = N_VLinearCombination(4, cvals, Xvecs, ark_mem->tempv1);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-compute-embedding",
                 "status = failed vector op, retval = %i", retval);
      return ARK_VECTOROP_ERR;
    }
    *dsmPtr = N_VWrmsNorm(ark_mem->tempv1, ark_mem->ewt);
    lsrkStep_DomEigUpdateLogic(ark_mem, step_mem, *dsmPtr, ark_mem->tempv2);
  }
  else
  {
    lsrkStep_DomEigUpdateLogic(ark_mem, step_mem, *dsmPtr, ark_mem->tempv2);
  }

  SUNLogInfo(ARK_LOGGER, "end-compute-embedding", "status = success");

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  lsrkStep_TakeStepSSPs2:

  This routine serves the primary purpose of the LSRKStepSSPs2 module:
  it performs a single SSPs2 step (with embedding).

  The output variable dsmPtr should contain estimate of the
  weighted local error if an embedding is present; otherwise it
  should be 0.

  The variables (ark_mem->tcur, ark_mem->ycur) should
  contain the current time and solution at the end of this time step.

  The input/output variable nflagPtr is used to gauge convergence
  of any algebraic solvers within the step.  As this routine
  involves no algebraic solve, it is set to 0 (success).

  The return value from this routine is:
            0 => step completed successfully
           >0 => step encountered recoverable failure;
                 reduce step and retry (if possible)
           <0 => step encountered unrecoverable failure
  ---------------------------------------------------------------*/
int lsrkStep_TakeStepSSPs2(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  int retval;
  ARKodeLSRKStepMem step_mem;

  /* initialize algebraic solver convergence flag to success,
     temporal error estimate to zero */
  *nflagPtr = ARK_SUCCESS;
  *dsmPtr   = ZERO;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  sunrealtype* cvals = step_mem->cvals;
  N_Vector* Xvecs    = step_mem->Xvecs;

  /* Initialize the current stage index */
  step_mem->istage = 0;

  /* Initialize method coefficients */
  const sunrealtype rs      = (sunrealtype)step_mem->req_stages;
  const sunrealtype sm1inv  = ONE / (rs - ONE);
  const sunrealtype hsm1inv = ark_mem->h * sm1inv;
  const sunrealtype rsinv   = ONE / rs;
  const sunrealtype hrsinv  = ark_mem->h * rsinv;
  sunrealtype hbt1, hbt2, hbt3;

  /* Embedding coefficients differ when req_stages == 2 */
  if (step_mem->req_stages == 2)
  {
    // from https://doi.org/10.1016/j.cam.2022.114325 pg 5
    hbt1 = ark_mem->h * SUN_RCONST(0.694021459207626);
    hbt2 = ZERO;
    hbt3 = ark_mem->h - hbt1;
  }
  else
  {
    hbt1 = hrsinv * (ONE + rsinv);
    hbt2 = hrsinv;
    hbt3 = hrsinv * (ONE - rsinv);
  }

  /* Begin stage 0 */
  SUNLogInfo(ARK_LOGGER, "begin-stages-list",
             "stage = %i, tcur = " SUN_FORMAT_G, 0, ark_mem->tcur);
  SUNLogExtraDebugVec(ARK_LOGGER, "stage", ark_mem->yn, "z_0(:) =");

  /* The method is not FSAL. Therefore, fn ​is computed at the beginning
     of the step unless the previous step failed or ARKODE updated fn. */
  if (!ark_mem->fn_is_current)
  {
    /* call the user-supplied pre-RHS function (if supplied) */
    if (ark_mem->PreRhsFn)
    {
      retval = ark_mem->PreRhsFn(ark_mem->tn, ark_mem->yn, ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed preprocess rhs, retval = %i", retval);
        return ARK_PRERHSFN_FAIL;
      }
    }

    retval = step_mem->fe(ark_mem->tn, ark_mem->yn, ark_mem->fn,
                          ark_mem->user_data);
    step_mem->nfe++;
    if (retval != ARK_SUCCESS)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed rhs eval, retval = %i", retval);
      return (ARK_RHSFUNC_FAIL);
    }
    ark_mem->fn_is_current = SUNTRUE;
  }

  SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->fn, "F_0(:) =");
  SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

  /* Begin stage 1 and accumulate embedding into tempv1 */
  ark_mem->tcur    = ark_mem->tn + hsm1inv;
  step_mem->istage = 1;
  SUNLogInfo(ARK_LOGGER, "begin-stages-list",
             "stage = %i, tcur = " SUN_FORMAT_G, 1, ark_mem->tcur);
  N_VLinearSum(ONE, ark_mem->yn, hsm1inv, ark_mem->fn, ark_mem->ycur);
  if (!ark_mem->fixedstep)
  {
    N_VLinearSum(ONE, ark_mem->yn, hbt1, ark_mem->fn, ark_mem->tempv1);
  }

  /* apply user-supplied stage postprocessing function (if supplied) */
  if (ark_mem->PostProcessStageFn)
  {
    retval = ark_mem->PostProcessStageFn(ark_mem->tcur, ark_mem->ycur,
                                         ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed postprocess stage, retval = %i", retval);
      return ARK_POSTPROCESS_STAGE_FAIL;
    }
  }

  /* Evaluate stages j = 2,...,step_mem->req_stages - 1 */
  for (int j = 2; j < step_mem->req_stages; j++)
  {
    /* Complete the previous stage (evaluate the RHS and store it in tempv2) */

    /* apply user-supplied stage preprocessing function (if supplied) */
    if (ark_mem->PreRhsFn)
    {
      retval = ark_mem->PreRhsFn(ark_mem->tcur, ark_mem->ycur,
                                 ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed preprocess rhs, retval = %i", retval);
        return ARK_PRERHSFN_FAIL;
      }
    }

    retval = step_mem->fe(ark_mem->tcur, ark_mem->ycur, ark_mem->tempv2,
                          ark_mem->user_data);
    step_mem->nfe++;

    SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->tempv2,
                        "F_%i(:) =", j - 1);
    SUNLogInfoIf(retval != 0, ARK_LOGGER, "end-stages-list",
                 "status = failed rhs eval, retval = %i", retval);

    if (retval < 0) { return ARK_RHSFUNC_FAIL; }
    if (retval > 0) { return RHSFUNC_RECVR; }

    SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

    /* Begin stage j (update the state and embedding) */
    step_mem->istage = j;
    ark_mem->tcur    = ark_mem->tn + j * hsm1inv;
    SUNLogInfo(ARK_LOGGER, "begin-stages-list",
               "stage = %i, tcur = " SUN_FORMAT_G, j, ark_mem->tcur);
    N_VLinearSum(ONE, ark_mem->ycur, hsm1inv, ark_mem->tempv2, ark_mem->ycur);
    if (!ark_mem->fixedstep)
    {
      N_VLinearSum(ONE, ark_mem->tempv1, hbt2, ark_mem->tempv2, ark_mem->tempv1);
    }

    /* apply user-supplied stage postprocessing function (if supplied) */
    if (ark_mem->PostProcessStageFn)
    {
      retval = ark_mem->PostProcessStageFn(ark_mem->tcur, ark_mem->ycur,
                                           ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed postprocess stage, retval = %i", retval);
        return ARK_POSTPROCESS_STAGE_FAIL;
      }
    }
  }

  /* Complete the next-to-last stage by evaluating the RHS and storing it in tempv2 */
  if (ark_mem->PreRhsFn)
  {
    retval = ark_mem->PreRhsFn(ark_mem->tcur, ark_mem->ycur, ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed preprocess rhs, retval = %i", retval);
      return ARK_PRERHSFN_FAIL;
    }
  }
  retval = step_mem->fe(ark_mem->tcur, ark_mem->ycur, ark_mem->tempv2,
                        ark_mem->user_data);
  step_mem->nfe++;

  SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->tempv2,
                      "F_%i(:) =", step_mem->req_stages - 1);
  SUNLogInfoIf(retval != 0, ARK_LOGGER, "end-stages-list",
               "status = failed rhs eval, retval = %i", retval);

  if (retval < 0) { return ARK_RHSFUNC_FAIL; }
  if (retval > 0) { return RHSFUNC_RECVR; }

  SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

  /* Compute the step solution */
  ark_mem->tcur    = ark_mem->tn + ark_mem->h;
  step_mem->istage = step_mem->req_stages;
  SUNLogInfo(ARK_LOGGER, "begin-stages-list", "stage = %i, tcur = " SUN_FORMAT_G,
             step_mem->req_stages, ark_mem->tcur);
  cvals[0] = ONE / (sm1inv * rs);
  Xvecs[0] = ark_mem->ycur;
  cvals[1] = rsinv;
  Xvecs[1] = ark_mem->yn;
  cvals[2] = hrsinv;
  Xvecs[2] = ark_mem->tempv2;
  retval   = N_VLinearCombination(3, cvals, Xvecs, ark_mem->ycur);
  if (retval != 0)
  {
    SUNLogInfo(ARK_LOGGER, "end-stages-list",
               "status = failed vector op, retval = %i", retval);
    return ARK_VECTOROP_ERR;
  }

  /* apply user-supplied step postprocessing function (if supplied) */
  if (ark_mem->PostProcessStepFn)
  {
    retval = ark_mem->PostProcessStepFn(ark_mem->tcur, ark_mem->ycur,
                                        ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed postprocess step, retval = %i", retval);
      return ARK_POSTPROCESS_STEP_FAIL;
    }
  }

  SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");
  SUNLogExtraDebugVec(ARK_LOGGER, "updated solution", ark_mem->ycur, "ycur(:) =");

  /* Compute yerr (if step adaptivity enabled) */
  if (!ark_mem->fixedstep)
  {
    N_VLinearSum(ONE, ark_mem->tempv1, hbt3, ark_mem->tempv2, ark_mem->tempv1);
    SUNLogExtraDebugVec(ARK_LOGGER, "embedded solution", ark_mem->tempv1,
                        "y_embedded(:) =");
    N_VLinearSum(ONE, ark_mem->ycur, -ONE, ark_mem->tempv1, ark_mem->tempv1);
    *dsmPtr = N_VWrmsNorm(ark_mem->tempv1, ark_mem->ewt);
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  lsrkStep_TakeStepSSPs3:

  This routine serves the primary purpose of the LSRKStepSSPs3 module:
  it performs a single SSPs3 step (with embedding).

  The SSP3 method differs significantly when s = 4. Therefore, the case
  where num_of_stages = 4 is considered separately to avoid unnecessary
  boolean checks and improve computational efficiency.

  The output variable dsmPtr should contain estimate of the
  weighted local error if an embedding is present; otherwise it
  should be 0.

  The variables (ark_mem->tcur, ark_mem->ycur) should
  contain the current time and solution at the end of this time step.

  The input/output variable nflagPtr is used to gauge convergence
  of any algebraic solvers within the step.  As this routine
  involves no algebraic solve, it is set to 0 (success).

  The return value from this routine is:
            0 => step completed successfully
           >0 => step encountered recoverable failure;
                 reduce step and retry (if possible)
           <0 => step encountered unrecoverable failure
  ---------------------------------------------------------------*/
int lsrkStep_TakeStepSSPs3(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  int retval;
  ARKodeLSRKStepMem step_mem;

  /* initialize algebraic solver convergence flag to success,
     temporal error estimate to zero */
  *nflagPtr = ARK_SUCCESS;
  *dsmPtr   = ZERO;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  sunrealtype* cvals = step_mem->cvals;
  N_Vector* Xvecs    = step_mem->Xvecs;

  /* Initialize the current stage index */
  step_mem->istage = 0;

  /* Initialize method coefficients */
  const sunrealtype rs     = (sunrealtype)step_mem->req_stages;
  const sunrealtype rn     = SUNRsqrt(rs);
  const sunrealtype hrat   = ark_mem->h / (rs - rn);
  const sunrealtype hrsinv = ark_mem->h / rs;
  const int in             = (int)SUNRround(rn);

  /* Begin stage 0 */
  SUNLogInfo(ARK_LOGGER, "begin-stages-list",
             "stage = %i, tcur = " SUN_FORMAT_G, 0, ark_mem->tcur);
  SUNLogExtraDebugVec(ARK_LOGGER, "stage", ark_mem->yn, "z_0(:) =");

  /* The method is not FSAL. Therefore, fn ​is computed at the beginning
     of the step unless ARKODE updated fn. */
  if (!ark_mem->fn_is_current)
  {
    /* call the user-supplied pre-RHS function (if supplied) */
    if (ark_mem->PreRhsFn)
    {
      retval = ark_mem->PreRhsFn(ark_mem->tn, ark_mem->yn, ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed preprocess rhs, retval = %i", retval);
        return ARK_PRERHSFN_FAIL;
      }
    }

    retval = step_mem->fe(ark_mem->tn, ark_mem->yn, ark_mem->fn,
                          ark_mem->user_data);
    step_mem->nfe++;
    if (retval != ARK_SUCCESS)
    {
      SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->fn, "F_0(:) =");
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed rhs eval, retval = %i", retval);
      return (ARK_RHSFUNC_FAIL);
    }
    ark_mem->fn_is_current = SUNTRUE;
  }

  SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->fn, "F_0(:) =");
  SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

  /* Begin stage 1 and accumulate embedding into tempv1 */
  ark_mem->tcur    = ark_mem->tn + hrat;
  step_mem->istage = 1;
  SUNLogInfo(ARK_LOGGER, "begin-stages-list",
             "stage = %i, tcur = " SUN_FORMAT_G, 1, ark_mem->tcur);
  N_VLinearSum(ONE, ark_mem->yn, hrat, ark_mem->fn, ark_mem->ycur);
  if (!ark_mem->fixedstep)
  {
    N_VLinearSum(ONE, ark_mem->yn, hrsinv, ark_mem->fn, ark_mem->tempv1);
  }

  /* apply user-supplied stage postprocessing function (if supplied) */
  if (ark_mem->PostProcessStageFn)
  {
    retval = ark_mem->PostProcessStageFn(ark_mem->tcur, ark_mem->ycur,
                                         ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed postprocess stage, retval = %i", retval);
      return ARK_POSTPROCESS_STAGE_FAIL;
    }
  }

  /* Evaluate first stage group */
  for (int j = 2; j <= ((in - 1) * (in - 2) / 2); j++)
  {
    /* Complete the previous stage (evaluate the RHS and store it in tempv3) */

    /* apply user-supplied stage preprocessing function (if supplied) */
    if (ark_mem->PreRhsFn)
    {
      retval = ark_mem->PreRhsFn(ark_mem->tcur, ark_mem->ycur,
                                 ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed preprocess rhs, retval = %i", retval);
        return ARK_PRERHSFN_FAIL;
      }
    }

    retval = step_mem->fe(ark_mem->tcur, ark_mem->ycur, ark_mem->tempv3,
                          ark_mem->user_data);
    step_mem->nfe++;

    SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->tempv3,
                        "F_%i(:) =", j - 1);
    SUNLogInfoIf(retval != 0, ARK_LOGGER, "end-stages-list",
                 "status = failed rhs eval, retval = %i", retval);

    if (retval < 0) { return ARK_RHSFUNC_FAIL; }
    if (retval > 0) { return RHSFUNC_RECVR; }

    SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

    /* Begin stage j (update the state and embedding) */
    ark_mem->tcur    = ark_mem->tn + j * hrat;
    step_mem->istage = j;
    SUNLogInfo(ARK_LOGGER, "begin-stages-list",
               "stage = %i, tcur = " SUN_FORMAT_G, j, ark_mem->tcur);
    N_VLinearSum(ONE, ark_mem->ycur, hrat, ark_mem->tempv3, ark_mem->ycur);
    if (!ark_mem->fixedstep)
    {
      N_VLinearSum(ONE, ark_mem->tempv1, hrsinv, ark_mem->tempv3,
                   ark_mem->tempv1);
    }

    /* apply user-supplied stage postprocessing function (if supplied) */
    if (ark_mem->PostProcessStageFn)
    {
      retval = ark_mem->PostProcessStageFn(ark_mem->tcur, ark_mem->ycur,
                                           ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed postprocess stage, retval = %i", retval);
        return ARK_POSTPROCESS_STAGE_FAIL;
      }
    }
  }

  /* Copy ycur into tempv2 before looping over second stage group */
  N_VScale(ONE, ark_mem->ycur, ark_mem->tempv2);

  /* Evaluate second stage group */
  for (int j = ((in - 1) * (in - 2) / 2 + 1); j <= (in * (in + 1) / 2 - 1); j++)
  {
    /* Complete the previous stage (evaluate the RHS and store it in tempv3) */

    /* apply user-supplied stage preprocessing function (if supplied) */
    if (ark_mem->PreRhsFn)
    {
      retval = ark_mem->PreRhsFn(ark_mem->tcur, ark_mem->ycur,
                                 ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed preprocess rhs, retval = %i", retval);
        return ARK_PRERHSFN_FAIL;
      }
    }

    retval = step_mem->fe(ark_mem->tcur, ark_mem->ycur, ark_mem->tempv3,
                          ark_mem->user_data);
    step_mem->nfe++;

    SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->tempv3,
                        "F_%i(:) =", j - 1);
    SUNLogInfoIf(retval != 0, ARK_LOGGER, "end-stages-list",
                 "status = failed rhs eval, retval = %i", retval);

    if (retval < 0) { return ARK_RHSFUNC_FAIL; }
    if (retval > 0) { return RHSFUNC_RECVR; }

    SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

    /* Begin stage j (update the state and embedding) */
    ark_mem->tcur    = ark_mem->tn + j * hrat;
    step_mem->istage = j;
    SUNLogInfo(ARK_LOGGER, "begin-stages-list",
               "stage = %i, tcur = " SUN_FORMAT_G, j, ark_mem->tcur);
    N_VLinearSum(ONE, ark_mem->ycur, hrat, ark_mem->tempv3, ark_mem->ycur);
    if (!ark_mem->fixedstep)
    {
      N_VLinearSum(ONE, ark_mem->tempv1, hrsinv, ark_mem->tempv3,
                   ark_mem->tempv1);
    }

    /* apply user-supplied stage postprocessing function (if supplied) */
    if (ark_mem->PostProcessStageFn)
    {
      retval = ark_mem->PostProcessStageFn(ark_mem->tcur, ark_mem->ycur,
                                           ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed postprocess stage, retval = %i", retval);
        return ARK_POSTPROCESS_STAGE_FAIL;
      }
    }
  }

  /* apply user-supplied stage preprocessing function (if supplied) */
  if (ark_mem->PreRhsFn)
  {
    retval = ark_mem->PreRhsFn(ark_mem->tcur, ark_mem->ycur, ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed preprocess rhs, retval = %i", retval);
      return ARK_PRERHSFN_FAIL;
    }
  }

  retval = step_mem->fe(ark_mem->tcur, ark_mem->ycur, ark_mem->tempv3,
                        ark_mem->user_data);
  step_mem->nfe++;

  SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->tempv3,
                      "F_%i(:) =", in * (in + 1) / 2 - 1);
  SUNLogInfoIf(retval != 0, ARK_LOGGER, "end-stages-list",
               "status = failed rhs eval, retval = %i", retval);

  if (retval < 0) { return ARK_RHSFUNC_FAIL; }
  if (retval > 0) { return RHSFUNC_RECVR; }

  SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

  /* Begin the next stage before final stage group */
  ark_mem->tcur    = ark_mem->tn + (in * (in - 1) / 2) * hrat;
  step_mem->istage = in * (in + 1) / 2;
  SUNLogInfo(ARK_LOGGER, "begin-stages-list", "stage = %i, tcur = " SUN_FORMAT_G,
             (in * (in + 1) / 2), ark_mem->tcur);
  cvals[0] = (rn - ONE) / (TWO * rn - ONE);
  Xvecs[0] = ark_mem->ycur;
  cvals[1] = rn / (TWO * rn - ONE);
  Xvecs[1] = ark_mem->tempv2;
  cvals[2] = (rn - ONE) * hrat / (TWO * rn - ONE);
  Xvecs[2] = ark_mem->tempv3;
  retval   = N_VLinearCombination(3, cvals, Xvecs, ark_mem->ycur);
  if (retval != 0)
  {
    SUNLogInfo(ARK_LOGGER, "end-stages-list",
               "status = failed vector op, retval = %i", retval);
    return ARK_VECTOROP_ERR;
  }
  if (!ark_mem->fixedstep)
  {
    N_VLinearSum(ONE, ark_mem->tempv1, hrsinv, ark_mem->tempv3, ark_mem->tempv1);
  }

  /* apply user-supplied stage postprocessing function (if supplied) */
  if (ark_mem->PostProcessStageFn)
  {
    retval = ark_mem->PostProcessStageFn(ark_mem->tcur, ark_mem->ycur,
                                         ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed postprocess stage, retval = %i", retval);
      return ARK_POSTPROCESS_STAGE_FAIL;
    }
  }

  /* Evaluate final stage group */
  for (int j = (in * (in + 1) / 2 + 1); j <= step_mem->req_stages; j++)
  {
    /* Complete the previous stage (evaluate the RHS and store it in tempv3) */

    /* apply user-supplied stage preprocessing function (if supplied) */
    if (ark_mem->PreRhsFn)
    {
      retval = ark_mem->PreRhsFn(ark_mem->tcur, ark_mem->ycur,
                                 ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed preprocess rhs, retval = %i", retval);
        return ARK_PRERHSFN_FAIL;
      }
    }

    retval = step_mem->fe(ark_mem->tcur, ark_mem->ycur, ark_mem->tempv3,
                          ark_mem->user_data);
    step_mem->nfe++;

    SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->tempv3,
                        "F_%i(:) =", j - 1);
    SUNLogInfoIf(retval != 0, ARK_LOGGER, "end-stages-list",
                 "status = failed rhs eval, retval = %i", retval);

    if (retval < 0) { return ARK_RHSFUNC_FAIL; }
    if (retval > 0) { return RHSFUNC_RECVR; }

    SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

    /* Begin stage j (update the state and embedding) */
    ark_mem->tcur    = ark_mem->tn + (j - in) * hrat;
    step_mem->istage = j;
    SUNLogInfo(ARK_LOGGER, "begin-stages-list",
               "stage = %i, tcur = " SUN_FORMAT_G, j, ark_mem->tcur);
    N_VLinearSum(ONE, ark_mem->ycur, hrat, ark_mem->tempv3, ark_mem->ycur);
    if (!ark_mem->fixedstep)
    {
      N_VLinearSum(ONE, ark_mem->tempv1, hrsinv, ark_mem->tempv3,
                   ark_mem->tempv1);
    }

    /* apply user-supplied stage or step postprocessing function (if supplied) */
    if (j < step_mem->req_stages && ark_mem->PostProcessStageFn)
    {
      retval = ark_mem->PostProcessStageFn(ark_mem->tcur, ark_mem->ycur,
                                           ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed postprocess stage, retval = %i", retval);
        return ARK_POSTPROCESS_STAGE_FAIL;
      }
    }
    else if (j == step_mem->req_stages && ark_mem->PostProcessStepFn)
    {
      retval = ark_mem->PostProcessStepFn(ark_mem->tn + ark_mem->h,
                                          ark_mem->ycur, ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed postprocess step, retval = %i", retval);
        return ARK_POSTPROCESS_STEP_FAIL;
      }
    }
  }

  SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");
  SUNLogExtraDebugVec(ARK_LOGGER, "updated solution", ark_mem->ycur, "ycur(:) =");

  /* Compute yerr (if step adaptivity enabled) */
  SUNLogInfo(ARK_LOGGER, "begin-compute-embedding", "");
  if (!ark_mem->fixedstep)
  {
    SUNLogExtraDebugVec(ARK_LOGGER, "embedded solution", ark_mem->tempv1,
                        "y_embedded(:) =");
    N_VLinearSum(ONE, ark_mem->ycur, -ONE, ark_mem->tempv1, ark_mem->tempv1);
    *dsmPtr = N_VWrmsNorm(ark_mem->tempv1, ark_mem->ewt);
  }
  SUNLogInfo(ARK_LOGGER, "end-compute-embedding", "status = success");

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  lsrkStep_TakeStepSSP43:

  This routine serves the primary purpose of the LSRKStepSSP43 module:
  it performs a single SSP43 step (with embedding).

  The SSP3 method differs significantly when s = 4. Therefore, the case
  where num_of_stages = 4 is considered separately to avoid unnecessary
  boolean checks and improve computational efficiency.

  The output variable dsmPtr should contain estimate of the
  weighted local error if an embedding is present; otherwise it
  should be 0.

  The variables (ark_mem->tcur, ark_mem->ycur) should
  contain the current time and solution at the end of this time step.

  The input/output variable nflagPtr is used to gauge convergence
  of any algebraic solvers within the step.  As this routine
  involves no algebraic solve, it is set to 0 (success).

  The return value from this routine is:
            0 => step completed successfully
           >0 => step encountered recoverable failure;
                 reduce step and retry (if possible)
           <0 => step encountered unrecoverable failure
  ---------------------------------------------------------------*/
int lsrkStep_TakeStepSSP43(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  int retval;
  ARKodeLSRKStepMem step_mem;

  /* initialize algebraic solver convergence flag to success,
     temporal error estimate to zero */
  *nflagPtr = ARK_SUCCESS;
  *dsmPtr   = ZERO;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  sunrealtype* cvals = step_mem->cvals;
  N_Vector* Xvecs    = step_mem->Xvecs;

  /* Initialize the current stage index */
  step_mem->istage = 0;

  /* Initialize method coefficients */
  const sunrealtype rs     = SUN_RCONST(4.0);
  const sunrealtype hp5    = ark_mem->h * SUN_RCONST(0.5);
  const sunrealtype hrsinv = ark_mem->h / rs;

  /* Begin stage 0 */
  SUNLogInfo(ARK_LOGGER, "begin-stages-list",
             "stage = %i, tcur = " SUN_FORMAT_G, 0, ark_mem->tcur);
  SUNLogExtraDebugVec(ARK_LOGGER, "stage", ark_mem->yn, "z_0(:) =");

  /* The method is not FSAL. Therefore, fn ​is computed at the beginning
     of the step unless ARKODE updated fn. */
  if (!ark_mem->fn_is_current)
  {
    /* call the user-supplied pre-RHS function (if supplied) */
    if (ark_mem->PreRhsFn)
    {
      retval = ark_mem->PreRhsFn(ark_mem->tn, ark_mem->yn, ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed preprocess rhs, retval = %i", retval);
        return ARK_PRERHSFN_FAIL;
      }
    }

    retval = step_mem->fe(ark_mem->tn, ark_mem->yn, ark_mem->fn,
                          ark_mem->user_data);
    step_mem->nfe++;
    if (retval != ARK_SUCCESS)
    {
      SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->fn, "F_0(:) =");
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed rhs eval, retval = %i", retval);
      return (ARK_RHSFUNC_FAIL);
    }
    ark_mem->fn_is_current = SUNTRUE;
  }

  SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->fn, "F_0(:) =");
  SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

  /* Begin stage 1 and accumulate embedding into tempv1 */
  ark_mem->tcur    = ark_mem->tn + hp5;
  step_mem->istage = 1;
  SUNLogInfo(ARK_LOGGER, "begin-stages-list",
             "stage = %i, tcur = " SUN_FORMAT_G, 1, ark_mem->tcur);
  N_VLinearSum(ONE, ark_mem->yn, hp5, ark_mem->fn, ark_mem->ycur);
  if (!ark_mem->fixedstep)
  {
    N_VLinearSum(ONE, ark_mem->yn, hrsinv, ark_mem->fn, ark_mem->tempv1);
  }

  /* apply user-supplied stage postprocessing function (if supplied) */
  if (ark_mem->PostProcessStageFn)
  {
    retval = ark_mem->PostProcessStageFn(ark_mem->tcur, ark_mem->ycur,
                                         ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed postprocess stage, retval = %i", retval);
      return ARK_POSTPROCESS_STAGE_FAIL;
    }
  }

  /* call the user-supplied pre-RHS function (if supplied) */
  if (ark_mem->PreRhsFn)
  {
    retval = ark_mem->PreRhsFn(ark_mem->tcur, ark_mem->ycur, ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed preprocess rhs, retval = %i", retval);
      return ARK_PRERHSFN_FAIL;
    }
  }

  /* Evaluate stage RHS */
  retval = step_mem->fe(ark_mem->tcur, ark_mem->ycur, ark_mem->tempv3,
                        ark_mem->user_data);
  step_mem->nfe++;

  SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->tempv3, "F_1(:) =");
  SUNLogInfoIf(retval != 0, ARK_LOGGER, "end-stages-list",
               "status = failed rhs eval, retval = %i", retval);

  if (retval < 0) { return ARK_RHSFUNC_FAIL; }
  if (retval > 0) { return RHSFUNC_RECVR; }

  SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

  /* Begin stage 2 */
  ark_mem->tcur    = ark_mem->tn + ark_mem->h;
  step_mem->istage = 2;
  SUNLogInfo(ARK_LOGGER, "begin-stages-list",
             "stage = %i, tcur = " SUN_FORMAT_G, 2, ark_mem->tcur);
  N_VLinearSum(ONE, ark_mem->ycur, hp5, ark_mem->tempv3, ark_mem->ycur);
  if (!ark_mem->fixedstep)
  {
    N_VLinearSum(ONE, ark_mem->tempv1, hrsinv, ark_mem->tempv3, ark_mem->tempv1);
  }

  /* apply user-supplied stage postprocessing function (if supplied) */
  if (ark_mem->PostProcessStageFn)
  {
    retval = ark_mem->PostProcessStageFn(ark_mem->tcur, ark_mem->ycur,
                                         ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed postprocess stage, retval = %i", retval);
      return ARK_POSTPROCESS_STAGE_FAIL;
    }
  }

  /* Evaluate stage RHS */

  /* apply user-supplied stage preprocessing function (if supplied) */
  if (ark_mem->PreRhsFn)
  {
    retval = ark_mem->PreRhsFn(ark_mem->tcur, ark_mem->ycur, ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed preprocess rhs, retval = %i", retval);
      return ARK_PRERHSFN_FAIL;
    }
  }

  retval = step_mem->fe(ark_mem->tcur, ark_mem->ycur, ark_mem->tempv3,
                        ark_mem->user_data);
  step_mem->nfe++;

  SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->tempv3, "F_2(:) =");
  SUNLogInfoIf(retval != 0, ARK_LOGGER, "end-stages-list",
               "status = failed rhs eval, retval = %i", retval);

  if (retval < 0) { return ARK_RHSFUNC_FAIL; }
  if (retval > 0) { return RHSFUNC_RECVR; }

  SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

  /* Begin stage 3 */
  ark_mem->tcur    = ark_mem->tn + hp5;
  step_mem->istage = 3;
  SUNLogInfo(ARK_LOGGER, "begin-stages-list",
             "stage = %i, tcur = " SUN_FORMAT_G, 3, ark_mem->tcur);
  cvals[0] = ONE / THREE;
  Xvecs[0] = ark_mem->ycur;
  cvals[1] = TWO / THREE;
  Xvecs[1] = ark_mem->yn;
  cvals[2] = ONE / SIX * ark_mem->h;
  Xvecs[2] = ark_mem->tempv3;
  retval   = N_VLinearCombination(3, cvals, Xvecs, ark_mem->ycur);
  if (retval != 0)
  {
    SUNLogInfo(ARK_LOGGER, "end-stages-list",
               "status = failed vector op, retval = %i", retval);
    return ARK_VECTOROP_ERR;
  }
  if (!ark_mem->fixedstep)
  {
    N_VLinearSum(ONE, ark_mem->tempv1, hrsinv, ark_mem->tempv3, ark_mem->tempv1);
  }

  /* apply user-supplied stage postprocessing function (if supplied) */
  if (ark_mem->PostProcessStageFn)
  {
    retval = ark_mem->PostProcessStageFn(ark_mem->tcur, ark_mem->ycur,
                                         ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed postprocess stage, retval = %i", retval);
      return ARK_POSTPROCESS_STAGE_FAIL;
    }
  }

  /* Evaluate stage RHS */

  /* apply user-supplied stage preprocessing function (if supplied) */
  if (ark_mem->PreRhsFn)
  {
    retval = ark_mem->PreRhsFn(ark_mem->tcur, ark_mem->ycur, ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed preprocess rhs, retval = %i", retval);
      return ARK_PRERHSFN_FAIL;
    }
  }

  retval = step_mem->fe(ark_mem->tcur, ark_mem->ycur, ark_mem->tempv3,
                        ark_mem->user_data);
  step_mem->nfe++;

  SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->tempv3, "F_3(:) =");
  SUNLogInfoIf(retval != 0, ARK_LOGGER, "end-stages-list",
               "status = failed rhs eval, retval = %i", retval);

  if (retval < 0) { return ARK_RHSFUNC_FAIL; }
  if (retval > 0) { return RHSFUNC_RECVR; }

  SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

  /* Compute the time step solution and embedding */
  ark_mem->tcur    = ark_mem->tn + ark_mem->h;
  step_mem->istage = 4;
  SUNLogInfo(ARK_LOGGER, "begin-stages-list",
             "stage = %i, tcur = " SUN_FORMAT_G, 4, ark_mem->tcur);
  N_VLinearSum(ONE, ark_mem->ycur, hp5, ark_mem->tempv3, ark_mem->ycur);

  /* apply user-supplied step postprocessing function (if supplied) */
  if (ark_mem->PostProcessStepFn)
  {
    retval = ark_mem->PostProcessStepFn(ark_mem->tcur, ark_mem->ycur,
                                        ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed postprocess step, retval = %i", retval);
      return ARK_POSTPROCESS_STEP_FAIL;
    }
  }

  SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");
  SUNLogExtraDebugVec(ARK_LOGGER, "updated solution", ark_mem->ycur, "ycur(:) =");

  /* Compute yerr (if step adaptivity enabled) */
  if (!ark_mem->fixedstep)
  {
    N_VLinearSum(ONE, ark_mem->tempv1, hrsinv, ark_mem->tempv3, ark_mem->tempv1);
    SUNLogExtraDebugVec(ARK_LOGGER, "embedded solution", ark_mem->tempv1,
                        "y_embedded(:) =");

    N_VLinearSum(ONE, ark_mem->ycur, -ONE, ark_mem->tempv1, ark_mem->tempv1);
    *dsmPtr = N_VWrmsNorm(ark_mem->tempv1, ark_mem->ewt);
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  lsrkStep_TakeStepSSP104:

  This routine serves the primary purpose of the LSRKStepSSP104 module:
  it performs a single SSP104 step (with embedding).

  The output variable dsmPtr should contain estimate of the
  weighted local error if an embedding is present; otherwise it
  should be 0.

  The variables (ark_mem->tcur, ark_mem->ycur) should
  contain the current time and solution at the end of this time step.

  The input/output variable nflagPtr is used to gauge convergence
  of any algebraic solvers within the step.  As this routine
  involves no algebraic solve, it is set to 0 (success).

  The return value from this routine is:
            0 => step completed successfully
           >0 => step encountered recoverable failure;
                 reduce step and retry (if possible)
           <0 => step encountered unrecoverable failure
  ---------------------------------------------------------------*/
int lsrkStep_TakeStepSSP104(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  int retval;
  ARKodeLSRKStepMem step_mem;

  /* initialize algebraic solver convergence flag to success,
     temporal error estimate to zero */
  *nflagPtr = ARK_SUCCESS;
  *dsmPtr   = ZERO;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  sunrealtype* cvals = step_mem->cvals;
  N_Vector* Xvecs    = step_mem->Xvecs;

  /* Initialize the current stage index */
  step_mem->istage = 0;

  /* Initialize method coefficients */
  const sunrealtype hsixth = ark_mem->h / SIX;
  const sunrealtype hfifth = ark_mem->h / FIVE;

  /* Begin stage 0 */
  SUNLogInfo(ARK_LOGGER, "begin-stages-list",
             "stage = %i, tcur = " SUN_FORMAT_G, 0, ark_mem->tcur);
  SUNLogExtraDebugVec(ARK_LOGGER, "stage", ark_mem->yn, "z_0(:) =");

  /* The method is not FSAL. Therefore, fn ​is computed at the beginning
     of the step unless ARKODE updated fn. */
  if (!ark_mem->fn_is_current)
  {
    /* call the user-supplied pre-RHS function (if supplied) */
    if (ark_mem->PreRhsFn)
    {
      retval = ark_mem->PreRhsFn(ark_mem->tn, ark_mem->yn, ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed preprocess rhs, retval = %i", retval);
        return ARK_PRERHSFN_FAIL;
      }
    }

    retval = step_mem->fe(ark_mem->tn, ark_mem->yn, ark_mem->fn,
                          ark_mem->user_data);
    step_mem->nfe++;
    if (retval != ARK_SUCCESS)
    {
      SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->fn, "F_0(:) =");
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed rhs eval, retval = %i", retval);
      return (ARK_RHSFUNC_FAIL);
    }
    ark_mem->fn_is_current = SUNTRUE;
  }

  SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->fn, "F_0(:) =");
  SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

  /* Copy yn into tempv2 for use in later stages */
  N_VScale(ONE, ark_mem->yn, ark_mem->tempv2);

  /* Begin stage 1 and accumulate embedding into tempv1 */
  ark_mem->tcur    = ark_mem->tn + hsixth;
  step_mem->istage = 1;
  SUNLogInfo(ARK_LOGGER, "begin-stages-list",
             "stage = %i, tcur = " SUN_FORMAT_G, 1, ark_mem->tcur);
  N_VLinearSum(ONE, ark_mem->yn, hsixth, ark_mem->fn, ark_mem->ycur);
  if (!ark_mem->fixedstep)
  {
    N_VLinearSum(ONE, ark_mem->yn, hfifth, ark_mem->fn, ark_mem->tempv1);
  }

  /* Evaluate stages j = 2,...,5 */
  for (int j = 2; j <= 5; j++)
  {
    /* Complete the previous stage (postprocesses the stage, evaluate the RHS, and
       store it in tempv3) */

    /* apply user-supplied stage postprocessing function (if supplied) */
    if (ark_mem->PostProcessStageFn)
    {
      retval = ark_mem->PostProcessStageFn(ark_mem->tcur, ark_mem->ycur,
                                           ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed postprocess stage, retval = %i", retval);
        return ARK_POSTPROCESS_STAGE_FAIL;
      }
    }

    /* apply user-supplied stage preprocessing function (if supplied) */
    if (ark_mem->PreRhsFn)
    {
      retval = ark_mem->PreRhsFn(ark_mem->tcur, ark_mem->ycur,
                                 ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed preprocess rhs, retval = %i", retval);
        return ARK_PRERHSFN_FAIL;
      }
    }

    retval = step_mem->fe(ark_mem->tcur, ark_mem->ycur, ark_mem->tempv3,
                          ark_mem->user_data);
    step_mem->nfe++;

    SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->tempv3,
                        "F_%i(:) =", j - 1);
    SUNLogInfoIf(retval != 0, ARK_LOGGER, "end-stages-list",
                 "status = failed rhs eval, retval = %i", retval);

    if (retval < 0) { return ARK_RHSFUNC_FAIL; }
    if (retval > 0) { return RHSFUNC_RECVR; }

    SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

    /* Begin stage j (update the state and embedding) */
    if (j == 5) { ark_mem->tcur = ark_mem->tn + 2 * hsixth; }
    else { ark_mem->tcur = ark_mem->tn + j * hsixth; }
    step_mem->istage = j;
    SUNLogInfo(ARK_LOGGER, "begin-stages-list",
               "stage = %i, tcur = " SUN_FORMAT_G, j, ark_mem->tcur);
    N_VLinearSum(ONE, ark_mem->ycur, hsixth, ark_mem->tempv3, ark_mem->ycur);
    if (j == 4 && !ark_mem->fixedstep)
    {
      N_VLinearSum(ONE, ark_mem->tempv1, SUN_RCONST(0.3) * ark_mem->h,
                   ark_mem->tempv3, ark_mem->tempv1);
    }
  }

  /* no need to call RHS preprocessing here, since the stage does not require
     a RHS function evaluation */

  /* Finish stage 5 by preparing for the final stage group */
  N_VLinearSum(SUN_RCONST(1.0) / SUN_RCONST(25.0), ark_mem->tempv2,
               SUN_RCONST(9.0) / SUN_RCONST(25.0), ark_mem->ycur,
               ark_mem->tempv2);
  N_VLinearSum(SUN_RCONST(15.0), ark_mem->tempv2, SUN_RCONST(-5.0),
               ark_mem->ycur, ark_mem->ycur);

  /* apply user-supplied stage postprocessing function (if supplied) */
  if (ark_mem->PostProcessStageFn)
  {
    retval = ark_mem->PostProcessStageFn(ark_mem->tcur, ark_mem->ycur,
                                         ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed postprocess stage, retval = %i", retval);
      return ARK_POSTPROCESS_STAGE_FAIL;
    }
  }

  /* Evaluate stages j = 6,...,9 */
  for (int j = 6; j <= 9; j++)
  {
    /* Complete the previous stage (evaluate the RHS and store in tempv3) */

    /* apply user-supplied stage preprocessing function (if supplied) */
    if (ark_mem->PreRhsFn)
    {
      retval = ark_mem->PreRhsFn(ark_mem->tcur, ark_mem->ycur,
                                 ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed preprocess rhs, retval = %i", retval);
        return ARK_PRERHSFN_FAIL;
      }
    }

    retval = step_mem->fe(ark_mem->tcur, ark_mem->ycur, ark_mem->tempv3,
                          ark_mem->user_data);
    step_mem->nfe++;

    SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->tempv3,
                        "F_%i(:) =", j - 1);
    SUNLogInfoIf(retval != 0, ARK_LOGGER, "end-stages-list",
                 "status = failed rhs eval, retval = %i", retval);

    if (retval < 0) { return ARK_RHSFUNC_FAIL; }
    if (retval > 0) { return RHSFUNC_RECVR; }

    SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");

    /* Begin stage j (update the state and embedding) */
    ark_mem->tcur    = ark_mem->tn + (j - 3) * hsixth;
    step_mem->istage = j;
    SUNLogInfo(ARK_LOGGER, "begin-stages-list",
               "stage = %i, tcur = " SUN_FORMAT_G, j, ark_mem->tcur);
    N_VLinearSum(ONE, ark_mem->ycur, hsixth, ark_mem->tempv3, ark_mem->ycur);

    if (j == 7 && !ark_mem->fixedstep)
    {
      N_VLinearSum(ONE, ark_mem->tempv1, hfifth, ark_mem->tempv3,
                   ark_mem->tempv1);
    }
    if (j == 9 && !ark_mem->fixedstep)
    {
      N_VLinearSum(ONE, ark_mem->tempv1, SUN_RCONST(0.3) * ark_mem->h,
                   ark_mem->tempv3, ark_mem->tempv1);
    }

    /* apply user-supplied stage postprocessing function (if supplied) */
    if (ark_mem->PostProcessStageFn)
    {
      retval = ark_mem->PostProcessStageFn(ark_mem->tcur, ark_mem->ycur,
                                           ark_mem->user_data);
      if (retval != 0)
      {
        SUNLogInfo(ARK_LOGGER, "end-stages-list",
                   "status = failed postprocess stage, retval = %i", retval);
        return ARK_POSTPROCESS_STAGE_FAIL;
      }
    }
  }

  /* Complete the previous stage (evaluate the RHS and store it in tempv3) */

  /* apply user-supplied stage preprocessing function (if supplied) */
  if (ark_mem->PreRhsFn)
  {
    retval = ark_mem->PreRhsFn(ark_mem->tcur, ark_mem->ycur, ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed preprocess rhs, retval = %i", retval);
      return ARK_PRERHSFN_FAIL;
    }
  }

  retval = step_mem->fe(ark_mem->tcur, ark_mem->ycur, ark_mem->tempv3,
                        ark_mem->user_data);
  step_mem->nfe++;

  SUNLogExtraDebugVec(ARK_LOGGER, "stage RHS", ark_mem->tempv3, "F_9(:) =");
  SUNLogInfoIf(retval != 0, ARK_LOGGER, "end-stages-list",
               "status = failed rhs eval, retval = %i", retval);

  if (retval < 0) { return ARK_RHSFUNC_FAIL; }
  if (retval > 0) { return RHSFUNC_RECVR; }

  SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");
  SUNLogInfo(ARK_LOGGER, "begin-stages-list",
             "stage = %i, tcur = " SUN_FORMAT_G, 10, ark_mem->tcur);

  /* Compute the final time step solution */
  step_mem->istage = 10;
  cvals[0]         = SUN_RCONST(0.6);
  Xvecs[0]         = ark_mem->ycur;
  cvals[1]         = ONE;
  Xvecs[1]         = ark_mem->tempv2;
  cvals[2]         = SUN_RCONST(0.1) * ark_mem->h;
  Xvecs[2]         = ark_mem->tempv3;
  retval           = N_VLinearCombination(3, cvals, Xvecs, ark_mem->ycur);
  if (retval != 0)
  {
    SUNLogInfo(ARK_LOGGER, "end-stages-list",
               "status = failed vector op, retval = %i", retval);
    return ARK_VECTOROP_ERR;
  }

  /* apply user-supplied step postprocessing function (if supplied) */
  if (ark_mem->PostProcessStepFn)
  {
    retval = ark_mem->PostProcessStepFn(ark_mem->tcur, ark_mem->ycur,
                                        ark_mem->user_data);
    if (retval != 0)
    {
      SUNLogInfo(ARK_LOGGER, "end-stages-list",
                 "status = failed postprocess step, retval = %i", retval);
      return ARK_POSTPROCESS_STEP_FAIL;
    }
  }

  SUNLogInfo(ARK_LOGGER, "end-stages-list", "status = success");
  SUNLogExtraDebugVec(ARK_LOGGER, "updated solution", ark_mem->ycur, "ycur(:) =");

  /* Compute yerr (if step adaptivity enabled) */
  if (!ark_mem->fixedstep)
  {
    SUNLogExtraDebugVec(ARK_LOGGER, "embedded solution", ark_mem->tempv1,
                        "y_embedded(:) =");
    N_VLinearSum(ONE, ark_mem->ycur, -ONE, ark_mem->tempv1, ark_mem->tempv1);
    *dsmPtr = N_VWrmsNorm(ark_mem->tempv1, ark_mem->ewt);
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  lsrkStep_Free frees all LSRKStep memory.
  ---------------------------------------------------------------*/
void lsrkStep_Free(ARKodeMem ark_mem)
{
  ARKodeLSRKStepMem step_mem;

  /* nothing to do if ark_mem is already NULL */
  if (ark_mem == NULL) { return; }

  /* conditional frees on non-NULL LSRKStep module */
  if (ark_mem->step_mem != NULL)
  {
    step_mem = (ARKodeLSRKStepMem)ark_mem->step_mem;

    /* free the reusable arrays for fused vector interface */
    if (step_mem->cvals != NULL)
    {
      free(step_mem->cvals);
      step_mem->cvals = NULL;
      ark_mem->lrw -= step_mem->nfusedopvecs;
    }
    if (step_mem->Xvecs != NULL)
    {
      free(step_mem->Xvecs);
      step_mem->Xvecs = NULL;
      ark_mem->liw -= step_mem->nfusedopvecs;
    }

    /* free the time stepper module itself */
    free(ark_mem->step_mem);
    ark_mem->step_mem = NULL;
  }
}

/*---------------------------------------------------------------
  lsrkStep_PrintMem:

  This routine outputs the memory from the LSRKStep structure to
  a specified file pointer (useful when debugging).
  ---------------------------------------------------------------*/
void lsrkStep_PrintMem(ARKodeMem ark_mem, FILE* outfile)
{
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return; }

  /* print integrator memory to file */
  switch (step_mem->LSRKmethod)
  {
  case ARKODE_LSRK_RKC_2:
    fprintf(outfile, "LSRKStep RKC time step module memory:\n");
    break;
  case ARKODE_LSRK_RKL_2:
    fprintf(outfile, "LSRKStep RKL time step module memory:\n");
    break;
  case ARKODE_LSRK_SSP_S_2:
    fprintf(outfile, "LSRKStep SSP(s,2) time step module memory:\n");
    break;
  case ARKODE_LSRK_SSP_S_3:
    fprintf(outfile, "LSRKStep SSP(s,3) time step module memory:\n");
    break;
  case ARKODE_LSRK_SSP_10_4:
    fprintf(outfile, "LSRKStep SSP(10,4) time step module memory:\n");
    break;
  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid method option.");
    return;
  }

  fprintf(outfile, "LSRKStep: q                   = %i\n", step_mem->q);
  fprintf(outfile, "LSRKStep: p                   = %i\n", step_mem->p);
  fprintf(outfile, "LSRKStep: istage              = %i\n", step_mem->istage);
  fprintf(outfile, "LSRKStep: req_stages          = %i\n", step_mem->req_stages);

  if (step_mem->is_SSP)
  {
    fprintf(outfile, "LSRKStep: nfe                 = %li\n", step_mem->nfe);
  }
  else if (!step_mem->is_SSP)
  {
    /* output integer quantities */
    fprintf(outfile, "LSRKStep: dom_eig_nst           = %li\n",
            step_mem->dom_eig_nst);
    fprintf(outfile, "LSRKStep: stage_max             = %i\n",
            step_mem->stage_max);
    fprintf(outfile, "LSRKStep: stage_max_limit       = %i\n",
            step_mem->stage_max_limit);
    fprintf(outfile, "LSRKStep: dom_eig_freq          = %li\n",
            step_mem->dom_eig_freq);
    fprintf(outfile, "LSRKStep: num_init_warmups      = %i\n",
            step_mem->num_init_warmups);
    fprintf(outfile, "LSRKStep: num_warmups           = %i\n",
            step_mem->num_warmups);

    /* output long integer quantities */
    fprintf(outfile, "LSRKStep: nfe                   = %li\n", step_mem->nfe);
    if (step_mem->DEE != NULL)
    {
      fprintf(outfile, "LSRKStep: nfeDQ               = %li\n", step_mem->nfeDQ);
      fprintf(outfile, "LSRKStep: num_iters           = %li\n",
              step_mem->num_dee_iters);
    }
    fprintf(outfile, "LSRKStep: dom_eig_num_evals     = %li\n",
            step_mem->dom_eig_num_evals);

    /* output sunrealtype quantities */
    // TODO(SRB): temporary fix for complex numbers
    fprintf(outfile,
            "LSRKStep: dom_eig               = " SUN_FORMAT_G SUN_FORMAT_SG
            "i\n",
            step_mem->lambdaR, step_mem->lambdaI);
    fprintf(outfile, "LSRKStep: spectral_radius       = " SUN_FORMAT_G "\n",
            step_mem->spectral_radius);
    fprintf(outfile, "LSRKStep: spectral_radius_max   = " SUN_FORMAT_G "\n",
            step_mem->spectral_radius_max);
    fprintf(outfile, "LSRKStep: spectral_radius_min   = " SUN_FORMAT_G "\n",
            step_mem->spectral_radius_min);
    fprintf(outfile, "LSRKStep: dom_eig_safety        = " SUN_FORMAT_G "\n",
            step_mem->dom_eig_safety);

    /* output sunbooleantype quantities */
    fprintf(outfile, "LSRKStep: dom_eig_update        = %d\n",
            step_mem->dom_eig_update);
    fprintf(outfile, "LSRKStep: dom_eig_is_current    = %d\n",
            step_mem->dom_eig_is_current);

    if (step_mem->DEE != NULL)
    {
      retval = SUNDomEigEstimator_Write(step_mem->DEE, outfile);
      if (retval != SUN_SUCCESS)
      {
        arkProcessError(ark_mem, ARK_DEE_FAIL, __LINE__, __func__, __FILE__,
                        "SUNDomEigEstimator_Write failed");
        return;
      }
    }
  }
  else
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid method type.");
    return;
  }
}

/*===============================================================
  Internal utility routines
  ===============================================================*/

/*---------------------------------------------------------------
  lsrkStep_AccessARKODEStepMem:

  Shortcut routine to unpack both ark_mem and step_mem structures
  from void* pointer.  If either is missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int lsrkStep_AccessARKODEStepMem(void* arkode_mem, const char* fname,
                                 ARKodeMem* ark_mem, ARKodeLSRKStepMem* step_mem)
{
  /* access ARKodeMem structure */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_ARK_NO_MEM);
    return ARK_MEM_NULL;
  }
  *ark_mem = (ARKodeMem)arkode_mem;

  /* access ARKodeLSRKStepMem structure */
  if ((*ark_mem)->step_mem == NULL)
  {
    arkProcessError(*ark_mem, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_LSRKSTEP_NO_MEM);
    return ARK_MEM_NULL;
  }
  *step_mem = (ARKodeLSRKStepMem)(*ark_mem)->step_mem;
  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  lsrkStep_AccessStepMem:

  Shortcut routine to unpack the step_mem structure from
  ark_mem.  If missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int lsrkStep_AccessStepMem(ARKodeMem ark_mem, const char* fname,
                           ARKodeLSRKStepMem* step_mem)
{
  /* access ARKodeLSRKStepMem structure */
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_LSRKSTEP_NO_MEM);
    return ARK_MEM_NULL;
  }
  *step_mem = (ARKodeLSRKStepMem)ark_mem->step_mem;
  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  lsrkStep_DomEigUpdateLogic:

  This routine checks if the step is accepted or not and reassigns
  the dom_eig update flags accordingly.
  ---------------------------------------------------------------*/

void lsrkStep_DomEigUpdateLogic(ARKodeMem ark_mem, ARKodeLSRKStepMem step_mem,
                                sunrealtype dsm, N_Vector fnew)
{
  if (dsm <= ONE)
  {
    N_VScale(ONE, fnew, ark_mem->fn);
    ark_mem->fn_is_current = SUNTRUE;

    step_mem->dom_eig_is_current = (step_mem->const_Jac == SUNTRUE);

    step_mem->dom_eig_update = SUNFALSE;
    if (ark_mem->nst + 1 >= step_mem->dom_eig_nst + step_mem->dom_eig_freq)
    {
      step_mem->dom_eig_update = !step_mem->dom_eig_is_current;
    }
  }
  else { step_mem->dom_eig_update = !step_mem->dom_eig_is_current; }
}

/*---------------------------------------------------------------
  lsrkStep_ComputeNewDomEig:

  This routine computes new dom_eig and returns SUN_SUCCESS.
  ---------------------------------------------------------------*/

int lsrkStep_ComputeNewDomEig(ARKodeMem ark_mem, ARKodeLSRKStepMem step_mem)
{
  int retval = SUN_SUCCESS;

  if (step_mem->DEE != NULL)
  {
    retval = SUNDomEigEstimator_Estimate(step_mem->DEE, &step_mem->lambdaR,
                                         &step_mem->lambdaI);
    step_mem->dom_eig_num_evals++;
    if (retval != SUN_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_DEE_FAIL, __LINE__, __func__, __FILE__,
                      "SUNDomEigEstimator_Estimate failed");
      return ARK_DEE_FAIL;
    }

    long int num_iters;
    retval = SUNDomEigEstimator_GetNumIters(step_mem->DEE, &num_iters);
    if (retval != SUN_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_DEE_FAIL, __LINE__, __func__, __FILE__,
                      "SUNDomEigEstimator_GetNumIters failed");
      return ARK_DEE_FAIL;
    }
    step_mem->num_dee_iters += num_iters;

    /* After the first call to SUNDomEigEstimator_Estimate, the number of warmups is set to
       num_warmups, this allows the successive calls to
       SUNDomEigEstimator_Estimate to use a different number of warmups. */
    if (step_mem->init_warmup)
    {
      retval = SUNDomEigEstimator_SetNumPreprocessIters(step_mem->DEE,
                                                        step_mem->num_warmups);
      if (retval != SUN_SUCCESS)
      {
        arkProcessError(ark_mem, ARK_DEE_FAIL, __LINE__, __func__, __FILE__,
                        "SUNDomEigEstimator_SetNumPreprocessIters failed");
        return ARK_DEE_FAIL;
      }
      step_mem->init_warmup = SUNFALSE;
    }
  }
  else if (step_mem->dom_eig_fn != NULL)
  {
    retval = step_mem->dom_eig_fn(ark_mem->tn, ark_mem->yn, ark_mem->fn,
                                  &step_mem->lambdaR, &step_mem->lambdaI,
                                  ark_mem->user_data, ark_mem->tempv1,
                                  ark_mem->tempv2, ark_mem->tempv3);
    step_mem->dom_eig_num_evals++;
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_DOMEIG_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to estimate the dominant eigenvalue");
      return ARK_DOMEIG_FAIL;
    }
  }
  else
  {
    arkProcessError(ark_mem, ARK_DOMEIG_FAIL, __LINE__, __func__, __FILE__,
                    "Unable to estimate the dominant eigenvalue: Either a user "
                    "provided function or a SUNDomEigEstimator is required");
    return ARK_DOMEIG_FAIL;
  }

  if (step_mem->lambdaR * ark_mem->h > ZERO)
  {
    arkProcessError(NULL, ARK_DOMEIG_FAIL, __LINE__, __func__, __FILE__,
                    "lambdaR*h must be nonpositive");
    return ARK_DOMEIG_FAIL;
  }
  else if (step_mem->lambdaR == SUN_RCONST(0.0) &&
           SUNRabs(step_mem->lambdaI) > SUN_RCONST(0.0))
  {
    arkProcessError(NULL, ARK_DOMEIG_FAIL, __LINE__, __func__, __FILE__,
                    "DomEig cannot be purely imaginary");
    return ARK_DOMEIG_FAIL;
  }

  step_mem->lambdaR *= step_mem->dom_eig_safety;
  step_mem->lambdaI *= step_mem->dom_eig_safety;
  step_mem->spectral_radius =
    SUNRsqrt(SUNSQR(step_mem->lambdaR) + SUNSQR(step_mem->lambdaI));

  step_mem->dom_eig_is_current = SUNTRUE;
  step_mem->dom_eig_nst        = ark_mem->nst;

  step_mem->spectral_radius_max =
    (step_mem->spectral_radius > step_mem->spectral_radius_max)
      ? step_mem->spectral_radius
      : step_mem->spectral_radius_max;

  if (step_mem->spectral_radius < step_mem->spectral_radius_min ||
      ark_mem->nst == 0)
  {
    step_mem->spectral_radius_min = step_mem->spectral_radius;
  }

  step_mem->dom_eig_update = SUNFALSE;

  return retval;
}

/*---------------------------------------------------------------
  lsrkStep_DQJtimes:

  This routine generates a difference quotient approximation to
  the Jacobian-vector product f_y(t,y) * v. The approximation is
  Jv = [f(y + v*sig) - f(y)]/sig, where sig = 1 / ||v||_WRMS,
  i.e. the WRMS norm of v*sig is 1.
  ---------------------------------------------------------------*/
int lsrkStep_DQJtimes(void* arkode_mem, N_Vector v, N_Vector Jv)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;

  sunrealtype sig, siginv;
  int iter, retval;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  sunrealtype t = ark_mem->tn;
  N_Vector y    = ark_mem->yn;
  N_Vector work = ark_mem->tempv3;

  /* Compute RHS function, if necessary. */
  if ((!ark_mem->fn_is_current && ark_mem->initsetup) ||
      (step_mem->step_nst != ark_mem->nst))
  {
    /* call the user-supplied pre-RHS function (if supplied) */
    if (ark_mem->PreRhsFn)
    {
      retval = ark_mem->PreRhsFn(t, y, ark_mem->user_data);
      if (retval != 0) { return ARK_PRERHSFN_FAIL; }
    }

    retval = step_mem->fe(t, y, ark_mem->fn, ark_mem->user_data);
    step_mem->nfeDQ++;
    if (retval != ARK_SUCCESS)
    {
      SUNLogExtraDebugVec(ARK_LOGGER, "DomEig JvTimes RHS", ark_mem->fn,
                          "F_n(:) =");
      SUNLogInfo(ARK_LOGGER, "DomEig JvTimes",
                 "status = failed rhs eval, retval = %i", retval);
      return (ARK_RHSFUNC_FAIL);
    }
    ark_mem->fn_is_current = SUNTRUE;
  }

  /* Initialize perturbation to 1/||v|| */
  sig = ONE / N_VWrmsNorm(v, ark_mem->ewt);

  for (iter = 0; iter < MAX_DQITERS; iter++)
  {
    /* Set work = y + sig*v */
    N_VLinearSum(sig, v, ONE, y, work);

    /* call the user-supplied pre-RHS function (if supplied) */
    if (ark_mem->PreRhsFn)
    {
      retval = ark_mem->PreRhsFn(t, work, ark_mem->user_data);
      if (retval != 0) { return ARK_PRERHSFN_FAIL; }
    }
    /* Set Jv = f(tn, y+sig*v) */
    retval = step_mem->fe(t, work, Jv, ark_mem->user_data);
    step_mem->nfeDQ++;
    if (retval == 0) { break; }
    if (retval < 0) { return (-1); }

    /* If f failed recoverably, shrink sig and retry */
    sig *= SUN_RCONST(0.25);
  }

  /* If retval still isn't 0, return with a recoverable failure */
  if (retval > 0) { return (+1); }

  /* Replace Jv by (Jv - fn)/sig */
  siginv = ONE / sig;
  N_VLinearSum(siginv, Jv, -siginv, ark_mem->fn, Jv);

  return ARK_SUCCESS;
}

/*===============================================================
  EOF
  ===============================================================*/
