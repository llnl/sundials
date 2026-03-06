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
 * This is the implementation file for the optional input and
 * output functions for the ARKODE LSRKStep time stepper module.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "arkode_lsrkstep_impl.h"

#include "sundials_cli.h"
#include "sundials_utils.h"

/*===============================================================
  Exported optional input functions.
  ===============================================================*/

/*---------------------------------------------------------------
  LSRKStepSetSTSMethod sets method
    ARKODE_LSRK_RKC_2
    ARKODE_LSRK_RKL_2
  ---------------------------------------------------------------*/
int LSRKStepSetSTSMethod(void* arkode_mem, ARKODE_LSRKMethodType method)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  switch (method)
  {
  case ARKODE_LSRK_RKC_2:
    ark_mem->step          = lsrkStep_TakeStepRKC;
    step_mem->is_SSP       = SUNFALSE;
    step_mem->nfusedopvecs = 5;
    step_mem->q = ark_mem->hadapt_mem->q = 2;
    step_mem->p = ark_mem->hadapt_mem->p = 2;
    step_mem->step_nst                   = 0;
    break;
  case ARKODE_LSRK_RKL_2:
    ark_mem->step          = lsrkStep_TakeStepRKL;
    step_mem->is_SSP       = SUNFALSE;
    step_mem->nfusedopvecs = 5;
    step_mem->q = ark_mem->hadapt_mem->q = 2;
    step_mem->p = ark_mem->hadapt_mem->p = 2;
    step_mem->step_nst                   = 0;
    break;
  case ARKODE_LSRK_SSP_S_2:
  case ARKODE_LSRK_SSP_S_3:
  case ARKODE_LSRK_SSP_10_4:
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "Invalid method option: Call LSRKStepCreateSSP to create an SSP method first.");
    break;

  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid method option.");
    return ARK_ILL_INPUT;
  }

  step_mem->LSRKmethod = method;

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  LSRKStepSetSSPMethod sets method
    ARKODE_LSRK_SSP_S_2
    ARKODE_LSRK_SSP_S_3
    ARKODE_LSRK_SSP_10_4
  ---------------------------------------------------------------*/

int LSRKStepSetSSPMethod(void* arkode_mem, ARKODE_LSRKMethodType method)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  switch (method)
  {
  case ARKODE_LSRK_RKC_2:
  case ARKODE_LSRK_RKL_2:
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "Invalid method option: Call LSRKStepCreateSTS to create an STS method first.");
    break;
  case ARKODE_LSRK_SSP_S_2:
    ark_mem->step          = lsrkStep_TakeStepSSPs2;
    step_mem->is_SSP       = SUNTRUE;
    step_mem->req_stages   = 10;
    step_mem->nfusedopvecs = 3;
    step_mem->q = ark_mem->hadapt_mem->q = 2;
    step_mem->p = ark_mem->hadapt_mem->p = 1;
    break;
  case ARKODE_LSRK_SSP_S_3:
    ark_mem->step          = lsrkStep_TakeStepSSPs3;
    step_mem->is_SSP       = SUNTRUE;
    step_mem->req_stages   = 9;
    step_mem->nfusedopvecs = 3;
    step_mem->q = ark_mem->hadapt_mem->q = 3;
    step_mem->p = ark_mem->hadapt_mem->p = 2;
    break;
  case ARKODE_LSRK_SSP_10_4:
    ark_mem->step          = lsrkStep_TakeStepSSP104;
    step_mem->is_SSP       = SUNTRUE;
    step_mem->req_stages   = 10;
    step_mem->nfusedopvecs = 3;
    step_mem->q = ark_mem->hadapt_mem->q = 4;
    step_mem->p = ark_mem->hadapt_mem->p = 3;
    break;

  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid method option.");
    return ARK_ILL_INPUT;
  }

  step_mem->LSRKmethod = method;

  return ARK_SUCCESS;
}

int LSRKStepSetSTSMethodByName(void* arkode_mem, const char* emethod)
{
  if (strcmp(emethod, "ARKODE_LSRK_RKC_2") == 0)
  {
    return LSRKStepSetSTSMethod(arkode_mem, ARKODE_LSRK_RKC_2);
  }
  if (strcmp(emethod, "ARKODE_LSRK_RKL_2") == 0)
  {
    return LSRKStepSetSTSMethod(arkode_mem, ARKODE_LSRK_RKL_2);
  }
  if ((strcmp(emethod, "ARKODE_LSRK_SSP_S_2") == 0) ||
      (strcmp(emethod, "ARKODE_LSRK_SSP_S_3") == 0) ||
      (strcmp(emethod, "ARKODE_LSRK_SSP_10_4") == 0))
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "Invalid method option: Call LSRKStepCreateSTS to create an STS method first.");
  }

  arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                  "Unknown method type");

  return ARK_ILL_INPUT;
}

int LSRKStepSetSSPMethodByName(void* arkode_mem, const char* emethod)
{
  if ((strcmp(emethod, "ARKODE_LSRK_RKC_2") == 0) ||
      (strcmp(emethod, "ARKODE_LSRK_RKL_2") == 0))
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "Invalid method option: Call LSRKStepCreateSSP to create an SSP method first.");
  }
  if (strcmp(emethod, "ARKODE_LSRK_SSP_S_2") == 0)
  {
    return LSRKStepSetSSPMethod(arkode_mem, ARKODE_LSRK_SSP_S_2);
  }
  if (strcmp(emethod, "ARKODE_LSRK_SSP_S_3") == 0)
  {
    return LSRKStepSetSSPMethod(arkode_mem, ARKODE_LSRK_SSP_S_3);
  }
  if (strcmp(emethod, "ARKODE_LSRK_SSP_10_4") == 0)
  {
    return LSRKStepSetSSPMethod(arkode_mem, ARKODE_LSRK_SSP_10_4);
  }

  arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                  "Unknown method type");

  return ARK_ILL_INPUT;
}

/*---------------------------------------------------------------
  LSRKStepSetDomEigFn specifies the dom_eig function.
  Specifies the dominant eigenvalue approximation routine to be used for determining
  the number of stages that will be used by either the RKC or RKL methods.
  ---------------------------------------------------------------*/
int LSRKStepSetDomEigFn(void* arkode_mem, ARKDomEigFn dom_eig)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  /* set the dom_eig routine pointer, and update relevant flags */
  step_mem->dom_eig_fn = dom_eig;

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  LSRKStepSetDomEigFrequency sets dom_eig computation frequency -
  Dominated Eigenvalue is recomputed after "nsteps" successful steps.

  nsteps = 0 refers to constant dominant eigenvalue
  nsteps < 0 resets the default value 25 and sets nonconstant dominant eigenvalue
  ---------------------------------------------------------------*/
int LSRKStepSetDomEigFrequency(void* arkode_mem, long int nsteps)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (nsteps < 0)
  {
    step_mem->dom_eig_freq = DOM_EIG_FREQ_DEFAULT;
    step_mem->const_Jac    = SUNFALSE;
  }

  if (nsteps == 0)
  {
    step_mem->const_Jac    = SUNTRUE;
    step_mem->dom_eig_freq = 1;
  }
  else
  {
    step_mem->dom_eig_freq = nsteps;
    step_mem->const_Jac    = SUNFALSE;
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  LSRKStepSetMaxNumStages sets the maximum number of stages allowed.
  If the combination of the maximum number of stages and the current
  time step size in the LSRKStep module does not allow for a stable
  step, the step routine returns to ARKODE for an updated (refined)
  step size. The number of such returns is tracked in a counter,
  which can be accessed using ARKodeGetNumExpSteps.
  ---------------------------------------------------------------*/
int LSRKStepSetMaxNumStages(void* arkode_mem, int stage_max_limit)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (stage_max_limit < 2)
  {
    step_mem->stage_max_limit = STAGE_MAX_LIMIT_DEFAULT;
  }
  else { step_mem->stage_max_limit = stage_max_limit; }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  LSRKStepSetDomEigSafetyFactor sets the safety factor for the DomEigs.
  Specifies a safety factor to use for the result of the dominant eigenvalue estimation function.
  This value is used to scale the magnitude of the dominant eigenvalue, in the hope of ensuring
  a sufficient number of stages for the method to be stable.  This input is only used for RKC
  and RKL methods.

  Calling this function with dom_eig_safety < 0 resets the default value
  ---------------------------------------------------------------*/
int LSRKStepSetDomEigSafetyFactor(void* arkode_mem, sunrealtype dom_eig_safety)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (dom_eig_safety < SUN_RCONST(1.0))
  {
    step_mem->dom_eig_safety = DOM_EIG_SAFETY_DEFAULT;
  }
  else { step_mem->dom_eig_safety = dom_eig_safety; }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  LSRKStepSetNumDomEigEstInitPreprocessIters sets the number of the preprocessing
  iterations before the very first estimate call.
  ---------------------------------------------------------------*/
int LSRKStepSetNumDomEigEstInitPreprocessIters(void* arkode_mem, int num_iters)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  /* This value will be used in lsrkStep_Init to set the number of preprocessing
     iterations for the first dominant eigenvalue estimate. If num_iters < 0,
     then the DEE's default will be used. */
  step_mem->num_init_warmups = num_iters;

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  LSRKStepSetNumDomEigEstPreprocessIters sets the number of the preprocessing
  iterations before each estimate call after the initial estimate call.
  ---------------------------------------------------------------*/
int LSRKStepSetNumDomEigEstPreprocessIters(void* arkode_mem, int num_iters)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (num_iters < 0) { step_mem->num_warmups = DOM_EIG_NUM_WARMUPS_DEFAULT; }
  else { step_mem->num_warmups = num_iters; }

  /* Set the number of iterations immediately (if possible) to allow the user to
     can change the value at any time during the integration. This value will be
     overridden for the first estimate by the value set with InitPreprocessIters
     above then reset to the supplied value afterward.

     A (perhaps pathological) corner case can occur where this value is not
     applied if a user detaches the DEE, calls this function then attaches a new
     DEE (or reattached the same DEE), and reinit is not called before
     evolve. In that case, whatever value the DEE has will be used. This could
     be avoided (with additional overhead) by calling SetNumPreprocessIters
     before every Estimate call. */
  if (step_mem->DEE)
  {
    retval = SUNDomEigEstimator_SetNumPreprocessIters(step_mem->DEE,
                                                      step_mem->num_init_warmups);
    if (retval != SUN_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_DEE_FAIL, __LINE__, __func__, __FILE__,
                      "SUNDomeEigEstimator_SetNumPreprocessIters failed");
      return ARK_DEE_FAIL;
    }
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  LSRKStepSetNumSSPStages sets the number of stages in the following
  SSP methods:

      ARKODE_LSRK_SSP_S_2  -- num_of_stages must be greater than or equal to 2
      ARKODE_LSRK_SSP_S_3  -- num_of_stages must be a perfect square greater than or equal to 4
      ARKODE_LSRK_SSP_10_4 -- num_of_stages must be equal to 10 - no need to call!

   Sets the number of stages, s in SSP(s, p) methods. This input is only utilized by
   SSPRK methods. Thus, this set routine must be called after calling LSRKStepSetSSPMethod.

   Calling this function with num_of_stages =< 0 resets the default value.
  ---------------------------------------------------------------*/
int LSRKStepSetNumSSPStages(void* arkode_mem, int num_of_stages)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (!step_mem->is_SSP)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Call this function only for SSP methods: Use "
                    "LSRKStepSetSSPMethod to declare SSP method type first!");
    return ARK_ILL_INPUT;
  }

  if (num_of_stages <= 0)
  {
    switch (step_mem->LSRKmethod)
    {
    case ARKODE_LSRK_SSP_S_2: step_mem->req_stages = 10; break;

    case ARKODE_LSRK_SSP_S_3: step_mem->req_stages = 9; break;

    case ARKODE_LSRK_SSP_10_4: step_mem->req_stages = 10; break;

    default:
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                      __FILE__, "Call LSRKStepSetSSPMethod to declare SSP method type first!");
      return ARK_ILL_INPUT;
      break;
    }
    return ARK_SUCCESS;
  }
  else
  {
    switch (step_mem->LSRKmethod)
    {
    case ARKODE_LSRK_SSP_S_2:
      if (num_of_stages < 2)
      {
        arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                        "num_of_stages must be greater than or equal to 2, or "
                        "set it less than or equal to 0 to reset the default "
                        "value");
        return ARK_ILL_INPUT;
      }
      break;

    case ARKODE_LSRK_SSP_S_3:
      /* The SSP3 method differs significantly when s = 4. Therefore, the case
      where num_of_stages = 4 is considered separately to avoid unnecessary
      boolean checks and improve computational efficiency. */

      /* We check that num_of_stages is a perfect square. Note the call to sqrt
       * rather than SUNRsqrt which could cause loss of precision if
       * sunrealtype is float. sqrt cannot produce a number bigger than INT_MAX
       * here so there's no problem with the cast back to int */
      if (num_of_stages < 4 || SUNSQR((int)sqrt(num_of_stages)) != num_of_stages)
      {
        arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                        "num_of_stages must be a perfect square greater than "
                        "or equal to 4, or set it less than or equal to 0 to "
                        "reset the default value");
        return ARK_ILL_INPUT;
      }
      if (num_of_stages == 4) { ark_mem->step = lsrkStep_TakeStepSSP43; }
      break;

    case ARKODE_LSRK_SSP_10_4:
      if (num_of_stages != 10)
      {
        arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                        "SSP10_4 method has a prefixed num_of_stages = 10");
        return ARK_ILL_INPUT;
      }
      break;

    default:
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                      __FILE__, "Call LSRKStepSetSSPMethod to declare SSP method type first!");
      return ARK_ILL_INPUT;
      break;
    }
    step_mem->req_stages = num_of_stages;
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  LSRKStepSetDomEigEstimator:

  This routine sets the dominant eigenvalue estimator DEE.
  ---------------------------------------------------------------*/
int LSRKStepSetDomEigEstimator(void* arkode_mem, SUNDomEigEstimator DEE)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;

  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (DEE == NULL)
  {
    step_mem->DEE = DEE;
    return ARK_SUCCESS;
  }

  if (DEE->ops == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Null SUNDomEigEstimator operations structure");
    return ARK_ILL_INPUT;
  }

  if (DEE->ops->estimate == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Null SUNDomEigEstimator estimate operation");
    return ARK_ILL_INPUT;
  }

  /* Attach the DEE pointer to the step memory */
  step_mem->DEE = DEE;

  /* Set the ATimes function for the DEE with A_data = arkode_mem */
  retval = SUNDomEigEstimator_SetATimes(DEE, arkode_mem, lsrkStep_DQJtimes);
  if (retval != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_DEE_FAIL, __LINE__, __func__, __FILE__,
                    "SUNDomEigEstimator_SetATimes failed");
    return ARK_DEE_FAIL;
  }

  return ARK_SUCCESS;
}

/*===============================================================
  Exported optional output functions.
  ===============================================================*/

/*---------------------------------------------------------------
  LSRKStepGetNumDomEigUpdates:

  Returns the number of dominant eigenvalue updates
  ---------------------------------------------------------------*/
int LSRKStepGetNumDomEigUpdates(void* arkode_mem, long int* dom_eig_num_evals)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (dom_eig_num_evals == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "dom_eig_num_evals cannot be NULL");
    return ARK_ILL_INPUT;
  }

  /* get values from step_mem */
  *dom_eig_num_evals = step_mem->dom_eig_num_evals;

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  LSRKStepGetMaxNumStages:

  Returns the max number of stages used
  ---------------------------------------------------------------*/
int LSRKStepGetMaxNumStages(void* arkode_mem, int* stage_max)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (stage_max == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "stage_max cannot be NULL");
    return ARK_ILL_INPUT;
  }

  /* get values from step_mem */
  *stage_max = step_mem->stage_max;

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  LSRKStepGetNumDomEigEstRhsEvals:

  Returns the number of RHS evals in DQ Jacobian computations
  ---------------------------------------------------------------*/
int LSRKStepGetNumDomEigEstRhsEvals(void* arkode_mem, long int* nfeDQ)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (nfeDQ == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "nfeDQ cannot be NULL");
    return ARK_ILL_INPUT;
  }

  /* get values from step_mem */
  *nfeDQ = step_mem->nfeDQ;

  return ARK_SUCCESS;
}

int LSRKStepGetNumDomEigEstIters(void* arkode_mem, long int* num_iters)
{
  ARKodeMem ark_mem;
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeLSRKStepMem structures */
  retval = lsrkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (num_iters == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "num_iters cannot be NULL");
    return ARK_ILL_INPUT;
  }

  /* get values from step_mem */
  *num_iters = step_mem->num_dee_iters;

  return ARK_SUCCESS;
}

/*===============================================================
  Private functions attached to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  lsrkStep_SetOption:

  Provides string-based control over LSRKStep-specific "set" routines.
  ---------------------------------------------------------------*/
int lsrkStep_SetOptions(ARKodeMem ark_mem, int* argidx, char* argv[],
                        size_t offset, sunbooleantype* arg_used)
{
  /* Set lists of keys, and the corresponding set routines */
  static const struct sunKeyCharPair char_pairs[] =
    {{"sts_method_name", LSRKStepSetSTSMethodByName},
     {"ssp_method_name", LSRKStepSetSSPMethodByName}};
  static const int num_char_keys = sizeof(char_pairs) / sizeof(*char_pairs);

  static const struct sunKeyLongPair long_pairs[] = {
    {"dom_eig_frequency", LSRKStepSetDomEigFrequency}};
  static const int num_long_keys = sizeof(long_pairs) / sizeof(*long_pairs);

  static const struct sunKeyIntPair int_pairs[] =
    {{"max_num_stages", LSRKStepSetMaxNumStages},
     {"num_ssp_stages", LSRKStepSetNumSSPStages},
     {"num_dom_eig_est_init_preprocess_iters",
      LSRKStepSetNumDomEigEstInitPreprocessIters},
     {"num_dom_eig_est_preprocess_iters", LSRKStepSetNumDomEigEstPreprocessIters}};
  static const int num_int_keys = sizeof(int_pairs) / sizeof(*int_pairs);

  static const struct sunKeyRealPair real_pairs[] = {
    {"dom_eig_safety_factor", LSRKStepSetDomEigSafetyFactor}};
  static const int num_real_keys = sizeof(real_pairs) / sizeof(*real_pairs);

  /* check all "char" keys */
  int j, retval;
  retval = sunCheckAndSetCharArgs((void*)ark_mem, argidx, argv, offset,
                                  char_pairs, num_char_keys, arg_used, &j);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "error setting key: %s", char_pairs[j].key);
    return retval;
  }
  if (*arg_used) { return ARK_SUCCESS; }

  /* check all "long int" keys */
  retval = sunCheckAndSetLongArgs((void*)ark_mem, argidx, argv, offset,
                                  long_pairs, num_long_keys, arg_used, &j);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "error setting key: %s", long_pairs[j].key);
    return retval;
  }
  if (*arg_used) { return ARK_SUCCESS; }

  /* check all "int" keys */
  retval = sunCheckAndSetIntArgs((void*)ark_mem, argidx, argv, offset,
                                 int_pairs, num_int_keys, arg_used, &j);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "error setting key: %s", int_pairs[j].key);
    return retval;
  }
  if (*arg_used) { return ARK_SUCCESS; }

  /* check all "real" keys */
  retval = sunCheckAndSetRealArgs((void*)ark_mem, argidx, argv, offset,
                                  real_pairs, num_real_keys, arg_used, &j);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "error setting key: %s", real_pairs[j].key);
    return retval;
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  lsrkStep_SetDefaults:

  Resets all LSRKStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.
  ---------------------------------------------------------------*/
int lsrkStep_SetDefaults(ARKodeMem ark_mem)
{
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  /* Set default values for integrator optional inputs
     (overwrite some adaptivity params for LSRKStep use) */
  step_mem->req_stages = 0; /* no stages */

  /* Spectral info */
  step_mem->dom_eig_safety   = DOM_EIG_SAFETY_DEFAULT;
  step_mem->dom_eig_freq     = DOM_EIG_FREQ_DEFAULT;
  step_mem->const_Jac        = SUNFALSE;
  step_mem->num_init_warmups = DOM_EIG_NUM_INIT_WARMUPS_DEFAULT;
  step_mem->num_warmups      = DOM_EIG_NUM_WARMUPS_DEFAULT;

  /* Load the default SUNAdaptController */
  retval = arkReplaceAdaptController(ark_mem, NULL, SUNTRUE);
  if (retval) { return retval; }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  lsrkStep_PrintAllStats:

  Prints integrator statistics for STS methods
  ---------------------------------------------------------------*/
int lsrkStep_PrintAllStats(ARKodeMem ark_mem, FILE* outfile, SUNOutputFormat fmt)
{
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  sunfprintf_long(outfile, fmt, SUNFALSE, "RHS fn evals", step_mem->nfe);
  if (step_mem->is_SSP)
  {
    sunfprintf_long(outfile, fmt, SUNFALSE, "Number of stages used",
                    step_mem->req_stages);
  }
  else
  {
    sunfprintf_long(outfile, fmt, SUNFALSE, "Number of dom_eig updates",
                    step_mem->dom_eig_num_evals);
    if (step_mem->DEE != NULL)
    {
      sunfprintf_long(outfile, fmt, SUNFALSE, "Number of fe calls for DEE",
                      step_mem->nfeDQ);
      sunfprintf_long(outfile, fmt, SUNFALSE, "Number of iterations for DEE",
                      step_mem->num_dee_iters);
    }
    sunfprintf_long(outfile, fmt, SUNFALSE, "Max. num. of stages used",
                    step_mem->stage_max);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Max. num. of stages allowed",
                    step_mem->stage_max_limit);
    sunfprintf_real(outfile, fmt, SUNFALSE, "Max. spectral radius",
                    step_mem->spectral_radius_max);
    sunfprintf_real(outfile, fmt, SUNFALSE, "Min. spectral radius",
                    step_mem->spectral_radius_min);
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  lsrkStep_WriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int lsrkStep_WriteParameters(ARKodeMem ark_mem, FILE* fp)
{
  ARKodeLSRKStepMem step_mem;
  int retval;

  /* access ARKodeLSRKStepMem structure */
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  /* print integrator parameters to file */
  switch (step_mem->LSRKmethod)
  {
  case ARKODE_LSRK_RKC_2:
    fprintf(fp, "LSRKStep RKC time step module parameters:\n");
    break;
  case ARKODE_LSRK_RKL_2:
    fprintf(fp, "LSRKStep RKL time step module parameters:\n");
    break;
  case ARKODE_LSRK_SSP_S_2:
    fprintf(fp, "LSRKStep SSP(s,2) time step module parameters:\n");
    break;
  case ARKODE_LSRK_SSP_S_3:
    fprintf(fp, "LSRKStep SSP(s,3) time step module parameters:\n");
    break;
  case ARKODE_LSRK_SSP_10_4:
    fprintf(fp, "LSRKStep SSP(10,4) time step module parameters:\n");

    break;

  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid method option.");
    return ARK_ILL_INPUT;
  }

  fprintf(fp, "  Method order %i\n", step_mem->q);
  fprintf(fp, "  Embedding order %i\n", step_mem->p);

  switch (step_mem->is_SSP)
  {
  case SUNTRUE:
    fprintf(fp, "  Number of stages used = %i\n", step_mem->req_stages);
    break;
  case SUNFALSE:
    fprintf(fp, "  Maximum number of stages allowed = %i\n",
            step_mem->stage_max_limit);
    if (step_mem->DEE != NULL)
    {
      fprintf(fp, "  Number of fe calls for DEE = %li\n", step_mem->nfeDQ);
    }
    fprintf(fp, "  Current spectral radius = " SUN_FORMAT_G "\n",
            step_mem->spectral_radius);
    fprintf(fp, "  Safety factor for the dom eig = " SUN_FORMAT_G "\n",
            step_mem->dom_eig_safety);
    fprintf(fp, "  Max num of successful steps before new dom eig update = %li\n",
            step_mem->dom_eig_freq);
    fprintf(fp, "  Number of first preprocessing warmups = %i\n",
            step_mem->num_init_warmups);
    fprintf(fp, "  Number of subsequent preprocessing warmups = %i\n",
            step_mem->num_warmups);
    fprintf(fp, "  Flag to indicate Jacobian is constant = %d\n",
            step_mem->const_Jac);
    break;
  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid method type.");
    return ARK_ILL_INPUT;
  }

  fprintf(fp, "\n");

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  lsrkStep_GetNumRhsEvals:

  Returns the current number of RHS calls
  ---------------------------------------------------------------*/
int lsrkStep_GetNumRhsEvals(ARKodeMem ark_mem, int partition_index,
                            long int* rhs_evals)
{
  ARKodeLSRKStepMem step_mem = NULL;

  /* access ARKodeLSRKStepMem structure */
  int retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (rhs_evals == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "rhs_evals is NULL");
    return ARK_ILL_INPUT;
  }

  if (partition_index > 0)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid partition index");
    return ARK_ILL_INPUT;
  }

  *rhs_evals = step_mem->nfe;

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  lsrkStep_GetEstLocalErrors: Returns the current local truncation
  error estimate vector
  ---------------------------------------------------------------*/
int lsrkStep_GetEstLocalErrors(ARKodeMem ark_mem, N_Vector ele)
{
  int retval;
  ARKodeLSRKStepMem step_mem;
  retval = lsrkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  /* return an error if local truncation error is not computed */
  if (ark_mem->fixedstep) { return ARK_STEPPER_UNSUPPORTED; }

  /* otherwise, copy local truncation error vector to output */
  N_VScale(ONE, ark_mem->tempv1, ele);
  return ARK_SUCCESS;
}

/*===============================================================
  EOF
  ===============================================================*/
