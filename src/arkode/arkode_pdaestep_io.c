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
 * This is the implementation file for the optional input and
 * output functions for the ARKODE PDAEStep time stepper module.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <ida/ida.h>

#include "arkode_pdaestep_impl.h"

#include "sundials_utils.h"

/*===============================================================
  Exported optional input functions.
  ===============================================================*/

/*---------------------------------------------------------------
  pdaeStep_SetUserData:

  Passes user-data pointer to attached linear solver modules.
  ---------------------------------------------------------------*/
int pdaeStep_SetUserData(ARKodeMem ark_mem, void* user_data)
{
  ARKodePDAEStepMem step_mem;
  int retval;

  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set user data in ARKODE LS mem */
  if (step_mem->lmem != NULL)
  {
    retval = arkLSSetUserData(ark_mem, user_data);
    if (retval != ARKLS_SUCCESS) { return (retval); }
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  pdaeStep_SetOrder:

  Specifies the method order
  ---------------------------------------------------------------*/
int pdaeStep_SetOrder(ARKodeMem ark_mem, int ord)
{
  ARKodePDAEStepMem step_mem;
  int retval;

  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set user-provided value, or default, depending on argument */
  if (ord <= 0) { step_mem->order = DEFAULT_ORDER; }
  else { step_mem->order = ord; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  pdaeStep_SetDefaults:

  Resets all PDAEStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.  Also leaves alone any data
  structures/options related to the ARKODE infrastructure itself
  (e.g., root-finding and post-process step).
  ---------------------------------------------------------------*/
int pdaeStep_SetDefaults(ARKodeMem ark_mem)
{
  ARKodePDAEStepMem step_mem;
  int retval;

  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  //TODO(SBR): determine what's actually needed
  /* Set default values for integrator optional inputs */
  step_mem->order          = DEFAULT_ORDER; /* method order */
  step_mem->maxcor         = DEFAULT_MAX_COR;    /* max nonlinear iters/stage */
  step_mem->nlscoef        = NLSCOEF;   /* nonlinear tolerance coefficient */
  step_mem->crdown         = CRDOWN; /* nonlinear convergence estimate coeff. */
  step_mem->rdiv           = RDIV;   /* nonlinear divergence tolerance */
  step_mem->dgmax    = DGMAX; /* max step change before recomputing J or P */
  step_mem->msbp     = MSBP;  /* max steps between updates to J or P */
  step_mem->jcur     = SUNFALSE;
  step_mem->convfail = ARK_NO_FAILURES;

  /* Remove pre-existing nonlinear solver object */
  if (step_mem->NLS && step_mem->ownNLS) { SUNNonlinSolFree(step_mem->NLS); }
  step_mem->NLS = NULL;

  //TODO(SBR): call arkReplaceAdaptController if we support adaptive steps

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  pdaeStep_SetMaxNonlinIters:

  Specifies the maximum number of nonlinear iterations during
  one solve.  A non-positive input implies a reset to the
  default value.
  ---------------------------------------------------------------*/
int pdaeStep_SetMaxNonlinIters(ARKodeMem ark_mem, int maxcor)
{
  ARKodePDAEStepMem step_mem;
  int retval;

  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  SUNFunctionBegin(ark_mem->sunctx);

  /* Return error message if no NLS module is present */
  if (step_mem->NLS == NULL)
  {
    arkProcessError(ark_mem, ARK_NLS_OP_ERR, __LINE__, __func__, __FILE__,
                    "No SUNNonlinearSolver object is present");
    return (ARK_ILL_INPUT);
  }

  /* argument <= 0 sets default, otherwise set input */
  if (maxcor <= 0) { step_mem->maxcor = DEFAULT_MAX_COR; }
  else { step_mem->maxcor = maxcor; }

  /* send argument to NLS structure */
  retval = SUNNonlinSolSetMaxIters(step_mem->NLS, step_mem->maxcor);
  if (retval != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_NLS_OP_ERR, __LINE__, __func__, __FILE__,
                    "Error setting maxcor in SUNNonlinearSolver object");
    return (ARK_NLS_OP_ERR);
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  pdaeStep_SetNonlinConvCoef:

  Specifies the coefficient in the nonlinear solver convergence
  test.  A non-positive input implies a reset to the default value.
  ---------------------------------------------------------------*/
int pdaeStep_SetNonlinConvCoef(ARKodeMem ark_mem, sunrealtype nlscoef)
{
  ARKodePDAEStepMem step_mem;
  int retval;

  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* argument <= 0 sets default, otherwise set input */
  if (nlscoef <= ZERO) { step_mem->nlscoef = DEFAULT_NLSCOEF; }
  else { step_mem->nlscoef = nlscoef; }

  return (ARK_SUCCESS);
}


/*===============================================================
  Exported optional output functions.
  ===============================================================*/

/*---------------------------------------------------------------
  pdaeStep_GetNumRhsEvals:

  Returns the current number of calls
  ---------------------------------------------------------------*/
int pdaeStep_GetNumRhsEvals(ARKodeMem ark_mem, int partition_index,
                            long int* rhs_evals)
{
  ARKodePDAEStepMem step_mem = NULL;

  int retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (rhs_evals == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "rhs_evals is NULL");
    return ARK_ILL_INPUT;
  }

  if (partition_index > step_mem->partitions)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid partition index");
    return ARK_ILL_INPUT;
  }

  if (partition_index < step_mem->partitions) {
    IDAGetNumResEvals(step_mem->ida_mems[partition_index], rhs_evals);
  } else {
    *rhs_evals = step_mem->nh;
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  pdaeStep_GetNumLinSolvSetups:

  Returns the current number of calls to the lsetup routine
  ---------------------------------------------------------------*/
int pdaeStep_GetNumLinSolvSetups(ARKodeMem ark_mem, long int* nlinsetups)
{
  ARKodePDAEStepMem step_mem;
  int retval;

  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* get value from step_mem */
  *nlinsetups = step_mem->nsetups;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  pdaeStep_GetNumNonlinSolvIters:

  Returns the current number of nonlinear solver iterations
  ---------------------------------------------------------------*/
int arkStep_GetNumNonlinSolvIters(ARKodeMem ark_mem, long int* nniters)
{
  ARKodePDAEStepMem step_mem;
  int retval;

  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  *nniters = step_mem->nls_iters;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  pdaeStep_GetNumNonlinSolvConvFails:

  Returns the current number of nonlinear solver convergence fails
  ---------------------------------------------------------------*/
int pdaeStep_GetNumNonlinSolvConvFails(ARKodeMem ark_mem, long int* nnfails)
{
  ARKodePDAEStepMem step_mem;
  int retval;

  /* access ARKodePDAEStepMem structure */
  retval = pdaeStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set output from step_mem */
  *nnfails = step_mem->nls_fails;

  return (ARK_SUCCESS);
}
/*---------------------------------------------------------------
  arkStep_GetNonlinSolvStats:

  Returns nonlinear solver statistics
  ---------------------------------------------------------------*/
int arkStep_GetNonlinSolvStats(ARKodeMem ark_mem, long int* nniters,
                               long int* nnfails)
{
  ARKodePDAEStepMem step_mem;
  int retval;

  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  *nniters = step_mem->nls_iters;
  *nnfails = step_mem->nls_fails;

  return (ARK_SUCCESS);
}