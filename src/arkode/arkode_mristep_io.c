/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ UMBC
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
 * This is the implementation file for the optional input and output functions
 * for the ARKODE MRIStep time stepper module.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "arkode_mristep_impl.h"

#include "sundials_cli.h"
#include "sundials_utils.h"

/*===============================================================
  Exported optional input functions.
  ===============================================================*/

/*---------------------------------------------------------------
  MRIStepSetCoupling:

  Specifies to use a customized coupling structure for the slow
  portion of the system.
  ---------------------------------------------------------------*/
int MRIStepSetCoupling(void* arkode_mem, MRIStepCoupling MRIC)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  sunindextype Tlrw, Tliw;

  /* access ARKodeMem and ARKodeMRIStepMem structures */
  retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* check for illegal inputs */
  if (MRIC == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_MRISTEP_NO_COUPLING);
    return (ARK_ILL_INPUT);
  }

  /* clear any existing parameters and coupling structure */
  step_mem->stages = 0;
  step_mem->q      = 0;
  step_mem->p      = 0;
  MRIStepCoupling_Space(step_mem->MRIC, &Tliw, &Tlrw);
  MRIStepCoupling_Free(step_mem->MRIC);
  step_mem->MRIC = NULL;
  ark_mem->liw -= Tliw;
  ark_mem->lrw -= Tlrw;

  /* set the relevant parameters */
  step_mem->stages = MRIC->stages;
  step_mem->q      = MRIC->q;
  step_mem->p      = MRIC->p;

  /* copy the coupling structure in step memory */
  step_mem->MRIC = MRIStepCoupling_Copy(MRIC);
  if (step_mem->MRIC == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_MRISTEP_NO_COUPLING);
    return (ARK_MEM_NULL);
  }
  MRIStepCoupling_Space(step_mem->MRIC, &Tliw, &Tlrw);
  ark_mem->liw += Tliw;
  ark_mem->lrw += Tlrw;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  MRIStepSetPreInnerFn:

  Sets the user-supplied function called BEFORE the inner evolve
  ---------------------------------------------------------------*/
int MRIStepSetPreInnerFn(void* arkode_mem, MRIStepPreInnerFn prefn)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeMRIStepMem structures */
  retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Set pre inner evolve function */
  step_mem->pre_inner_evolve = prefn;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  MRIStepSetPostInnerFn:

  Sets the user-supplied function called AFTER the inner evolve
  ---------------------------------------------------------------*/
int MRIStepSetPostInnerFn(void* arkode_mem, MRIStepPostInnerFn postfn)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeMRIStepMem structures */
  retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Set pre inner evolve function */
  step_mem->post_inner_evolve = postfn;

  return (ARK_SUCCESS);
}

/*===============================================================
  Exported optional output functions.
  ===============================================================*/

/*---------------------------------------------------------------
  mriStep_GetNumRhsEvals:

  Returns the current number of RHS calls
  ---------------------------------------------------------------*/
int mriStep_GetNumRhsEvals(ARKodeMem ark_mem, int partition_index,
                           long int* rhs_evals)
{
  ARKodeMRIStepMem step_mem = NULL;

  /* access ARKodeMRIStepMem structure */
  int retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (rhs_evals == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "rhs_evals is NULL");
    return ARK_ILL_INPUT;
  }

  if (partition_index > 1)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid partition index");
    return ARK_ILL_INPUT;
  }

  switch (partition_index)
  {
  case 0: *rhs_evals = step_mem->nfse; break;
  case 1: *rhs_evals = step_mem->nfsi; break;
  default: *rhs_evals = step_mem->nfse + step_mem->nfsi; break;
  }

  return ARK_SUCCESS;
}

int MRIStepGetNumRhsEvals(void* arkode_mem, long int* nfse_evals,
                          long int* nfsi_evals)
{
  int retval = ARK_SUCCESS;

  retval = ARKodeGetNumRhsEvals(arkode_mem, 0, nfse_evals);
  if (retval != ARK_SUCCESS) { return retval; }

  retval = ARKodeGetNumRhsEvals(arkode_mem, 1, nfsi_evals);
  if (retval != ARK_SUCCESS) { return retval; }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  MRIStepGetCurrentCoupling:

  Sets pointer to the slow coupling structure currently in use.
  ---------------------------------------------------------------*/
int MRIStepGetCurrentCoupling(void* arkode_mem, MRIStepCoupling* MRIC)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeMRIStepMem structures */
  retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* get coupling structure from step_mem */
  *MRIC = step_mem->MRIC;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  MRIStepGetLastInnerStepFlag:

  Returns the last return value from the inner stepper.
  ---------------------------------------------------------------*/
int MRIStepGetLastInnerStepFlag(void* arkode_mem, int* flag)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeMRIStepMem structures */
  retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* get the last return value from the inner stepper */
  *flag = step_mem->stepper->last_flag;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  MRIStepGetNumInnerStepperFails:

  Returns the number of recoverable failures encountered by the
  inner stepper.
  ---------------------------------------------------------------*/
int MRIStepGetNumInnerStepperFails(void* arkode_mem, long int* inner_fails)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeMRIStepMem structures */
  retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set output from step_mem */
  *inner_fails = step_mem->inner_fails;

  return (ARK_SUCCESS);
}

/*===============================================================
  Private functions attached to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  mriStep_SetOption:

  Provides command-line control over MRIStep-specific "set" routines.
  ---------------------------------------------------------------*/
int mriStep_SetOptions(ARKodeMem ark_mem, int* argidx, char* argv[],
                       size_t offset, sunbooleantype* arg_used)
{
  /* The only MRIStep-specific "Set" routine takes a custom MRIStepCoupling
     table; however, these may be specified by name, so here we'll support
     a key to specify the MRIStepCoupling table name,
     create the table with that name, attach it to MRIStep (who copies its
     values), and then free the table. */
  if (strcmp(argv[*argidx] + offset, "coupling_table_name") == 0)
  {
    (*argidx)++;
    MRIStepCoupling Coupling = MRIStepCoupling_LoadTableByName(argv[*argidx]);
    if (Coupling == NULL)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "error setting key %s %s (invalid table name)",
                      argv[(*argidx) - 1], argv[*argidx]);
      return ARK_ILL_INPUT;
    }
    int retval = MRIStepSetCoupling(ark_mem, Coupling);
    MRIStepCoupling_Free(Coupling);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                      "error setting key %s %s (SetCoupling failed)",
                      argv[(*argidx) - 1], argv[*argidx]);
      return retval;
    }
    *arg_used = SUNTRUE;
    return ARK_SUCCESS;
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  mriStep_SetAdaptController:

  Specifies a temporal adaptivity controller for MRIStep to use.
  If a non-MRI controller is provided, this just passes that
  through to arkReplaceAdaptController.  However, if an MRI
  controller is provided, then this wraps that inside a
  "SUNAdaptController_MRIStep" wrapper, which will properly
  interact with the fast integration module.
  ---------------------------------------------------------------*/
int mriStep_SetAdaptController(ARKodeMem ark_mem, SUNAdaptController C)
{
  /* Retrieve the controller type */
  SUNAdaptController_Type ctype = SUNAdaptController_GetType(C);

  /* If this does not have MRI type, then just pass to ARKODE */
  if (ctype != SUN_ADAPTCONTROLLER_MRI_H_TOL)
  {
    return (arkReplaceAdaptController(ark_mem, C, SUNFALSE));
  }

  /* Create the mriStepControl wrapper, pass that to ARKODE, and give ownership
     of the wrapper to ARKODE */
  SUNAdaptController Cwrapper = SUNAdaptController_MRIStep(ark_mem, C);
  return (arkReplaceAdaptController(ark_mem, Cwrapper, SUNTRUE));
}

/*---------------------------------------------------------------
  mriStep_SetUserData:

  Passes user-data pointer to attached linear solver module.
  ---------------------------------------------------------------*/
int mriStep_SetUserData(ARKodeMem ark_mem, void* user_data)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set user data in ARKODELS mem */
  if (step_mem->lmem != NULL)
  {
    retval = arkLSSetUserData(ark_mem, user_data);
    if (retval != ARKLS_SUCCESS) { return (retval); }
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_SetDefaults:

  Resets all MRIStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.
  ---------------------------------------------------------------*/
int mriStep_SetDefaults(ARKodeMem ark_mem)
{
  ARKodeMRIStepMem step_mem;
  sunindextype Clenrw, Cleniw;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Set default values for integrator optional inputs */
  step_mem->q              = 3;        /* method order */
  step_mem->p              = 0;        /* embedding order */
  step_mem->predictor      = 0;        /* trivial predictor */
  step_mem->linear         = SUNFALSE; /* nonlinear problem */
  step_mem->linear_timedep = SUNTRUE;  /* dfs/dy depends on t */
  step_mem->deduce_rhs     = SUNFALSE; /* deduce fi on result of NLS */
  step_mem->maxcor         = MAXCOR;   /* max nonlinear iters/stage */
  step_mem->nlscoef        = NLSCOEF;  /* nonlinear tolerance coefficient */
  step_mem->crdown         = CRDOWN; /* nonlinear convergence estimate coeff. */
  step_mem->rdiv           = RDIV;   /* nonlinear divergence tolerance */
  step_mem->dgmax          = DGMAX;  /* max gamma change to recompute J or P */
  step_mem->msbp           = MSBP;   /* max steps between updating J or P */
  step_mem->stages         = 0;      /* no stages */
  step_mem->istage         = 0;      /* current stage index */
  step_mem->jcur           = SUNFALSE;
  step_mem->convfail       = ARK_NO_FAILURES;
  step_mem->stage_predict  = NULL; /* no user-supplied stage predictor */

  /* Remove pre-existing nonlinear solver object */
  if (step_mem->NLS && step_mem->ownNLS) { SUNNonlinSolFree(step_mem->NLS); }
  step_mem->NLS = NULL;

  /* Remove pre-existing coupling table */
  if (step_mem->MRIC)
  {
    MRIStepCoupling_Space(step_mem->MRIC, &Cleniw, &Clenrw);
    ark_mem->lrw -= Clenrw;
    ark_mem->liw -= Cleniw;
    MRIStepCoupling_Free(step_mem->MRIC);
  }
  step_mem->MRIC = NULL;

  /* Load the default SUNAdaptController */
  retval = arkReplaceAdaptController(ark_mem, NULL, SUNTRUE);
  if (retval) { return retval; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_SetLinear:

  Specifies that the implicit slow function, fsi(t,y), is linear
  in y, and to tighten the linear solver tolerances while taking
  only one Newton iteration.  DO NOT USE IN COMBINATION WITH THE
  FIXED-POINT SOLVER.  Automatically tightens DeltaGammaMax
  to ensure that step size changes cause Jacobian recomputation.

  The argument should be 1 or 0, where 1 indicates that the
  Jacobian of fs with respect to y depends on time, and
  0 indicates that it is not time dependent.  Alternately, when
  using an iterative linear solver this flag denotes time
  dependence of the preconditioner.
  ---------------------------------------------------------------*/
int mriStep_SetLinear(ARKodeMem ark_mem, int timedepend)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set parameters */
  step_mem->linear         = SUNTRUE;
  step_mem->linear_timedep = (timedepend == 1);
  step_mem->dgmax          = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_SetNonlinear:

  Specifies that the implicit slow function, fsi(t,y), is
  nonlinear in y.  Used to undo a previous call to
  mriStep_SetLinear.  Automatically loosens DeltaGammaMax back to
  default value.
  ---------------------------------------------------------------*/
int mriStep_SetNonlinear(ARKodeMem ark_mem)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set parameters */
  step_mem->linear         = SUNFALSE;
  step_mem->linear_timedep = SUNTRUE;
  step_mem->dgmax          = DGMAX;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_SetOrder:

  Specifies the method order
  ---------------------------------------------------------------*/
int mriStep_SetOrder(ARKodeMem ark_mem, int ord)
{
  int retval;
  ARKodeMRIStepMem step_mem;
  sunindextype Tlrw, Tliw;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval) { return (retval); }

  /* check for illegal inputs */
  if (ord <= 0) { step_mem->q = 3; }
  else { step_mem->q = ord; }

  /* Clear tables, the user is requesting a change in method or a reset to
     defaults. Tables will be set in InitialSetup. */
  step_mem->stages = 0;
  step_mem->p      = 0;
  MRIStepCoupling_Space(step_mem->MRIC, &Tliw, &Tlrw);
  MRIStepCoupling_Free(step_mem->MRIC);
  step_mem->MRIC = NULL;
  ark_mem->liw -= Tliw;
  ark_mem->lrw -= Tlrw;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_SetNonlinCRDown:

  Specifies the user-provided nonlinear convergence constant
  crdown.  Legal values are strictly positive; illegal values
  imply a reset to the default.
  ---------------------------------------------------------------*/
int mriStep_SetNonlinCRDown(ARKodeMem ark_mem, sunrealtype crdown)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* if argument legal set it, otherwise set default */
  if (crdown <= ZERO) { step_mem->crdown = CRDOWN; }
  else { step_mem->crdown = crdown; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_SetNonlinRDiv:

  Specifies the user-provided nonlinear convergence constant
  rdiv.  Legal values are strictly positive; illegal values
  imply a reset to the default.
  ---------------------------------------------------------------*/
int mriStep_SetNonlinRDiv(ARKodeMem ark_mem, sunrealtype rdiv)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* if argument legal set it, otherwise set default */
  if (rdiv <= ZERO) { step_mem->rdiv = RDIV; }
  else { step_mem->rdiv = rdiv; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_SetDeltaGammaMax:

  Specifies the user-provided linear setup decision constant
  dgmax.  Legal values are strictly positive; illegal values imply
  a reset to the default.
  ---------------------------------------------------------------*/
int mriStep_SetDeltaGammaMax(ARKodeMem ark_mem, sunrealtype dgmax)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* if argument legal set it, otherwise set default */
  if (dgmax <= ZERO) { step_mem->dgmax = DGMAX; }
  else { step_mem->dgmax = dgmax; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_SetLSetupFrequency:

  Specifies the user-provided linear setup decision constant
  msbp.  Positive values give the frequency for calling lsetup;
  negative values imply recomputation of lsetup at each nonlinear
  solve; a zero value implies a reset to the default.
  ---------------------------------------------------------------*/
int mriStep_SetLSetupFrequency(ARKodeMem ark_mem, int msbp)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* if argument legal set it, otherwise set default */
  if (msbp == 0) { step_mem->msbp = MSBP; }
  else { step_mem->msbp = msbp; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_SetPredictorMethod:

  Specifies the method to use for predicting implicit solutions.
  Non-default choices are {1,2,3,4}, all others will use default
  (trivial) predictor.
  ---------------------------------------------------------------*/
int mriStep_SetPredictorMethod(ARKodeMem ark_mem, int pred_method)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set parameter */
  step_mem->predictor = pred_method;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_SetMaxNonlinIters:

  Specifies the maximum number of nonlinear iterations during
  one solve.  A non-positive input implies a reset to the
  default value.
  ---------------------------------------------------------------*/
int mriStep_SetMaxNonlinIters(ARKodeMem ark_mem, int maxcor)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Return error message if no NLS module is present */
  if (step_mem->NLS == NULL)
  {
    arkProcessError(ark_mem, ARK_NLS_OP_ERR, __LINE__, __func__, __FILE__,
                    "No SUNNonlinearSolver object is present");
    return (ARK_ILL_INPUT);
  }

  /* argument <= 0 sets default, otherwise set input */
  if (maxcor <= 0) { step_mem->maxcor = MAXCOR; }
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
  mriStep_SetNonlinConvCoef:

  Specifies the coefficient in the nonlinear solver convergence
  test.  A non-positive input implies a reset to the default value.
  ---------------------------------------------------------------*/
int mriStep_SetNonlinConvCoef(ARKodeMem ark_mem, sunrealtype nlscoef)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* argument <= 0 sets default, otherwise set input */
  if (nlscoef <= ZERO) { step_mem->nlscoef = NLSCOEF; }
  else { step_mem->nlscoef = nlscoef; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_SetStagePredictFn:  Specifies a user-provided step
  predictor function having type ARKStagePredictFn.  A
  NULL input function disables calls to this routine.
  ---------------------------------------------------------------*/
int mriStep_SetStagePredictFn(ARKodeMem ark_mem, ARKStagePredictFn PredictStage)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure and set function pointer */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  step_mem->stage_predict = PredictStage;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_SetDeduceImplicitRhs:

  Specifies if an optimization is used to avoid an evaluation of
  fi after a nonlinear solve for an implicit stage.  If stage
  postprocessecing in enabled, this option is ignored, and fi is
  never deduced.

  An argument of SUNTRUE indicates that fi is deduced to compute
  fi(z_i), and SUNFALSE indicates that fi(z_i) is computed with
  an additional evaluation of fi.
  ---------------------------------------------------------------*/
int mriStep_SetDeduceImplicitRhs(ARKodeMem ark_mem, sunbooleantype deduce)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure and set function pointer */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  step_mem->deduce_rhs = deduce;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_GetCurrentGamma: Returns the current value of gamma
  ---------------------------------------------------------------*/
int mriStep_GetCurrentGamma(ARKodeMem ark_mem, sunrealtype* gamma)
{
  int retval;
  ARKodeMRIStepMem step_mem;
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }
  *gamma = step_mem->gamma;
  return (retval);
}

/*---------------------------------------------------------------
  mriStep_GetEstLocalErrors: Returns the current local truncation
  error estimate vector
  ---------------------------------------------------------------*/
int mriStep_GetEstLocalErrors(ARKodeMem ark_mem, N_Vector ele)
{
  int retval;
  ARKodeMRIStepMem step_mem;
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* return an error if local truncation error is not computed */
  if ((ark_mem->fixedstep && (ark_mem->AccumErrorType == ARK_ACCUMERROR_NONE)) ||
      (step_mem->p <= 0))
  {
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* otherwise, copy local truncation error vector to output */
  N_VScale(ONE, ark_mem->tempv1, ele);
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_GetNumLinSolvSetups:

  Returns the current number of calls to the lsetup routine
  ---------------------------------------------------------------*/
int mriStep_GetNumLinSolvSetups(ARKodeMem ark_mem, long int* nlinsetups)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* get value from step_mem */
  *nlinsetups = step_mem->nsetups;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_GetNumNonlinSolvIters:

  Returns the current number of nonlinear solver iterations
  ---------------------------------------------------------------*/
int mriStep_GetNumNonlinSolvIters(ARKodeMem ark_mem, long int* nniters)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  *nniters = step_mem->nls_iters;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_GetNumNonlinSolvConvFails:

  Returns the current number of nonlinear solver convergence fails
  ---------------------------------------------------------------*/
int mriStep_GetNumNonlinSolvConvFails(ARKodeMem ark_mem, long int* nnfails)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set output from step_mem */
  *nnfails = step_mem->nls_fails;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_GetNonlinSolvStats:

  Returns nonlinear solver statistics
  ---------------------------------------------------------------*/
int mriStep_GetNonlinSolvStats(ARKodeMem ark_mem, long int* nniters,
                               long int* nnfails)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  *nniters = step_mem->nls_iters;
  *nnfails = step_mem->nls_fails;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_PrintAllStats:

  Prints integrator statistics
  ---------------------------------------------------------------*/
int mriStep_PrintAllStats(ARKodeMem ark_mem, FILE* outfile, SUNOutputFormat fmt)
{
  ARKodeMRIStepMem step_mem;
  ARKLsMem arkls_mem;
  int retval;

  /* access ARKode MRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* function evaluations */
  sunfprintf_long(outfile, fmt, SUNTRUE, "Explicit slow RHS fn evals",
                  step_mem->nfse);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Implicit slow RHS fn evals",
                  step_mem->nfsi);

  /* inner stepper and nonlinear solver stats */
  sunfprintf_long(outfile, fmt, SUNFALSE, "Inner stepper failures",
                  step_mem->inner_fails);
  sunfprintf_long(outfile, fmt, SUNFALSE, "NLS iters", step_mem->nls_iters);
  sunfprintf_long(outfile, fmt, SUNFALSE, "NLS fails", step_mem->nls_fails);
  if (ark_mem->nst > 0)
  {
    sunfprintf_real(outfile, fmt, SUNFALSE, "NLS iters per step",
                    (sunrealtype)step_mem->nls_iters / (sunrealtype)ark_mem->nst);
  }

  /* linear solver stats */
  sunfprintf_long(outfile, fmt, SUNFALSE, "LS setups", step_mem->nsetups);
  if (ark_mem->step_getlinmem(ark_mem))
  {
    arkls_mem = (ARKLsMem)(ark_mem->step_getlinmem(ark_mem));
    sunfprintf_long(outfile, fmt, SUNFALSE, "Jac fn evals", arkls_mem->nje);
    sunfprintf_long(outfile, fmt, SUNFALSE, "LS RHS fn evals", arkls_mem->nfeDQ);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Prec setup evals", arkls_mem->npe);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Prec solves", arkls_mem->nps);
    sunfprintf_long(outfile, fmt, SUNFALSE, "LS iters", arkls_mem->nli);
    sunfprintf_long(outfile, fmt, SUNFALSE, "LS fails", arkls_mem->ncfl);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Jac-times setups",
                    arkls_mem->njtsetup);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Jac-times evals",
                    arkls_mem->njtimes);
    if (step_mem->nls_iters > 0)
    {
      sunfprintf_real(outfile, fmt, SUNFALSE, "LS iters per NLS iter",
                      (sunrealtype)arkls_mem->nli /
                        (sunrealtype)step_mem->nls_iters);
      sunfprintf_real(outfile, fmt, SUNFALSE, "Jac evals per NLS iter",
                      (sunrealtype)arkls_mem->nje /
                        (sunrealtype)step_mem->nls_iters);
      sunfprintf_real(outfile, fmt, SUNFALSE, "Prec evals per NLS iter",
                      (sunrealtype)arkls_mem->npe /
                        (sunrealtype)step_mem->nls_iters);
    }
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  mriStep_WriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int mriStep_WriteParameters(ARKodeMem ark_mem, FILE* fp)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* print integrator parameters to file */
  fprintf(fp, "MRIStep time step module parameters:\n");
  fprintf(fp, "  Method order %i\n", step_mem->q);
  if (step_mem->linear)
  {
    fprintf(fp, "  Linear implicit problem");
    if (step_mem->linear_timedep)
    {
      fprintf(fp, " (time-dependent Jacobian)\n");
    }
    else { fprintf(fp, " (time-independent Jacobian)\n"); }
  }
  if (step_mem->explicit_rhs && step_mem->implicit_rhs)
  {
    fprintf(fp, "  ImEx slow time scale\n");
  }
  else if (step_mem->implicit_rhs)
  {
    fprintf(fp, "  Implicit slow time scale\n");
  }
  else { fprintf(fp, "  Explicit slow time scale\n"); }

  if (step_mem->implicit_rhs)
  {
    fprintf(fp, "  Implicit predictor method = %i\n", step_mem->predictor);
    fprintf(fp, "  Implicit solver tolerance coefficient = " SUN_FORMAT_G "\n",
            step_mem->nlscoef);
    fprintf(fp, "  Maximum number of nonlinear corrections = %i\n",
            step_mem->maxcor);
    fprintf(fp, "  Nonlinear convergence rate constant = " SUN_FORMAT_G "\n",
            step_mem->crdown);
    fprintf(fp, "  Nonlinear divergence tolerance = " SUN_FORMAT_G "\n",
            step_mem->rdiv);
    fprintf(fp, "  Gamma factor LSetup tolerance = " SUN_FORMAT_G "\n",
            step_mem->dgmax);
    fprintf(fp, "  Number of steps between LSetup calls = %i\n", step_mem->msbp);
  }
  fprintf(fp, "\n");

  return (ARK_SUCCESS);
}

/*===============================================================
  Exported-but-deprecated user-callable functions.
  ===============================================================*/

int MRIStepResize(void* arkode_mem, N_Vector y0, sunrealtype t0,
                  ARKVecResizeFn resize, void* resize_data)
{
  return (ARKodeResize(arkode_mem, y0, ONE, t0, resize, resize_data));
}

int MRIStepReset(void* arkode_mem, sunrealtype tR, N_Vector yR)
{
  return (ARKodeReset(arkode_mem, tR, yR));
}

int MRIStepSStolerances(void* arkode_mem, sunrealtype reltol, sunrealtype abstol)
{
  return (ARKodeSStolerances(arkode_mem, reltol, abstol));
}

int MRIStepSVtolerances(void* arkode_mem, sunrealtype reltol, N_Vector abstol)
{
  return (ARKodeSVtolerances(arkode_mem, reltol, abstol));
}

int MRIStepWFtolerances(void* arkode_mem, ARKEwtFn efun)
{
  return (ARKodeWFtolerances(arkode_mem, efun));
}

int MRIStepSetLinearSolver(void* arkode_mem, SUNLinearSolver LS, SUNMatrix A)
{
  return (ARKodeSetLinearSolver(arkode_mem, LS, A));
}

int MRIStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g)
{
  return (ARKodeRootInit(arkode_mem, nrtfn, g));
}

int MRIStepSetDefaults(void* arkode_mem)
{
  return (ARKodeSetDefaults(arkode_mem));
}

int MRIStepSetOrder(void* arkode_mem, int ord)
{
  return (ARKodeSetOrder(arkode_mem, ord));
}

int MRIStepSetInterpolantType(void* arkode_mem, int itype)
{
  return (ARKodeSetInterpolantType(arkode_mem, itype));
}

int MRIStepSetInterpolantDegree(void* arkode_mem, int degree)
{
  return (ARKodeSetInterpolantDegree(arkode_mem, degree));
}

int MRIStepSetDenseOrder(void* arkode_mem, int dord)
{
  return (ARKodeSetInterpolantDegree(arkode_mem, dord));
}

int MRIStepSetNonlinearSolver(void* arkode_mem, SUNNonlinearSolver NLS)
{
  return (ARKodeSetNonlinearSolver(arkode_mem, NLS));
}

int MRIStepSetNlsRhsFn(void* arkode_mem, ARKRhsFn nls_fi)
{
  return (ARKodeSetNlsRhsFn(arkode_mem, nls_fi));
}

int MRIStepSetLinear(void* arkode_mem, int timedepend)
{
  return (ARKodeSetLinear(arkode_mem, timedepend));
}

int MRIStepSetNonlinear(void* arkode_mem)
{
  return (ARKodeSetNonlinear(arkode_mem));
}

int MRIStepSetMaxNumSteps(void* arkode_mem, long int mxsteps)
{
  return (ARKodeSetMaxNumSteps(arkode_mem, mxsteps));
}

int MRIStepSetNonlinCRDown(void* arkode_mem, sunrealtype crdown)
{
  return (ARKodeSetNonlinCRDown(arkode_mem, crdown));
}

int MRIStepSetNonlinRDiv(void* arkode_mem, sunrealtype rdiv)
{
  return (ARKodeSetNonlinRDiv(arkode_mem, rdiv));
}

int MRIStepSetDeltaGammaMax(void* arkode_mem, sunrealtype dgmax)
{
  return (ARKodeSetDeltaGammaMax(arkode_mem, dgmax));
}

int MRIStepSetLSetupFrequency(void* arkode_mem, int msbp)
{
  return (ARKodeSetLSetupFrequency(arkode_mem, msbp));
}

int MRIStepSetPredictorMethod(void* arkode_mem, int pred_method)
{
  return (ARKodeSetPredictorMethod(arkode_mem, pred_method));
}

int MRIStepSetMaxNonlinIters(void* arkode_mem, int maxcor)
{
  return (ARKodeSetMaxNonlinIters(arkode_mem, maxcor));
}

int MRIStepSetNonlinConvCoef(void* arkode_mem, sunrealtype nlscoef)
{
  return (ARKodeSetNonlinConvCoef(arkode_mem, nlscoef));
}

int MRIStepSetMaxHnilWarns(void* arkode_mem, int mxhnil)
{
  return (ARKodeSetMaxHnilWarns(arkode_mem, mxhnil));
}

int MRIStepSetInterpolateStopTime(void* arkode_mem, sunbooleantype interp)
{
  return (ARKodeSetInterpolateStopTime(arkode_mem, interp));
}

int MRIStepSetStopTime(void* arkode_mem, sunrealtype tstop)
{
  return (ARKodeSetStopTime(arkode_mem, tstop));
}

int MRIStepClearStopTime(void* arkode_mem)
{
  return (ARKodeClearStopTime(arkode_mem));
}

int MRIStepSetFixedStep(void* arkode_mem, sunrealtype hfixed)
{
  return (ARKodeSetFixedStep(arkode_mem, hfixed));
}

int MRIStepSetRootDirection(void* arkode_mem, int* rootdir)
{
  return (ARKodeSetRootDirection(arkode_mem, rootdir));
}

int MRIStepSetNoInactiveRootWarn(void* arkode_mem)
{
  return (ARKodeSetNoInactiveRootWarn(arkode_mem));
}

int MRIStepSetUserData(void* arkode_mem, void* user_data)
{
  return (ARKodeSetUserData(arkode_mem, user_data));
}

int MRIStepSetPostprocessStepFn(void* arkode_mem, ARKPostProcessFn ProcessStep)
{
  return (ARKodeSetPostprocessStepFn(arkode_mem, ProcessStep));
}

int MRIStepSetPostprocessStageFn(void* arkode_mem, ARKPostProcessFn ProcessStage)
{
  return (ARKodeSetPostprocessStageFn(arkode_mem, ProcessStage));
}

int MRIStepSetStagePredictFn(void* arkode_mem, ARKStagePredictFn PredictStage)
{
  return (ARKodeSetStagePredictFn(arkode_mem, PredictStage));
}

int MRIStepSetDeduceImplicitRhs(void* arkode_mem, sunbooleantype deduce)
{
  return (ARKodeSetDeduceImplicitRhs(arkode_mem, deduce));
}

int MRIStepSetJacFn(void* arkode_mem, ARKLsJacFn jac)
{
  return (ARKodeSetJacFn(arkode_mem, jac));
}

int MRIStepSetJacEvalFrequency(void* arkode_mem, long int msbj)
{
  return (ARKodeSetJacEvalFrequency(arkode_mem, msbj));
}

int MRIStepSetLinearSolutionScaling(void* arkode_mem, sunbooleantype onoff)
{
  return (ARKodeSetLinearSolutionScaling(arkode_mem, onoff));
}

int MRIStepSetEpsLin(void* arkode_mem, sunrealtype eplifac)
{
  return (ARKodeSetEpsLin(arkode_mem, eplifac));
}

int MRIStepSetLSNormFactor(void* arkode_mem, sunrealtype nrmfac)
{
  return (ARKodeSetLSNormFactor(arkode_mem, nrmfac));
}

int MRIStepSetPreconditioner(void* arkode_mem, ARKLsPrecSetupFn psetup,
                             ARKLsPrecSolveFn psolve)
{
  return (ARKodeSetPreconditioner(arkode_mem, psetup, psolve));
}

int MRIStepSetJacTimes(void* arkode_mem, ARKLsJacTimesSetupFn jtsetup,
                       ARKLsJacTimesVecFn jtimes)
{
  return (ARKodeSetJacTimes(arkode_mem, jtsetup, jtimes));
}

int MRIStepSetJacTimesRhsFn(void* arkode_mem, ARKRhsFn jtimesRhsFn)
{
  return (ARKodeSetJacTimesRhsFn(arkode_mem, jtimesRhsFn));
}

int MRIStepSetLinSysFn(void* arkode_mem, ARKLsLinSysFn linsys)
{
  return (ARKodeSetLinSysFn(arkode_mem, linsys));
}

int MRIStepEvolve(void* arkode_mem, sunrealtype tout, N_Vector yout,
                  sunrealtype* tret, int itask)
{
  return (ARKodeEvolve(arkode_mem, tout, yout, tret, itask));
}

int MRIStepGetDky(void* arkode_mem, sunrealtype t, int k, N_Vector dky)
{
  return (ARKodeGetDky(arkode_mem, t, k, dky));
}

int MRIStepComputeState(void* arkode_mem, N_Vector zcor, N_Vector z)
{
  return (ARKodeComputeState(arkode_mem, zcor, z));
}

int MRIStepGetNumLinSolvSetups(void* arkode_mem, long int* nlinsetups)
{
  return (ARKodeGetNumLinSolvSetups(arkode_mem, nlinsetups));
}

int MRIStepGetWorkSpace(void* arkode_mem, long int* lenrw, long int* leniw)
{
  return (ARKodeGetWorkSpace(arkode_mem, lenrw, leniw));
}

int MRIStepGetNumSteps(void* arkode_mem, long int* nssteps)
{
  return (ARKodeGetNumSteps(arkode_mem, nssteps));
}

int MRIStepGetLastStep(void* arkode_mem, sunrealtype* hlast)
{
  return (ARKodeGetLastStep(arkode_mem, hlast));
}

int MRIStepGetCurrentTime(void* arkode_mem, sunrealtype* tcur)
{
  return (ARKodeGetCurrentTime(arkode_mem, tcur));
}

int MRIStepGetCurrentState(void* arkode_mem, N_Vector* state)
{
  return (ARKodeGetCurrentState(arkode_mem, state));
}

int MRIStepGetCurrentGamma(void* arkode_mem, sunrealtype* gamma)
{
  return (ARKodeGetCurrentGamma(arkode_mem, gamma));
}

int MRIStepGetTolScaleFactor(void* arkode_mem, sunrealtype* tolsfact)
{
  return (ARKodeGetTolScaleFactor(arkode_mem, tolsfact));
}

int MRIStepGetErrWeights(void* arkode_mem, N_Vector eweight)
{
  return (ARKodeGetErrWeights(arkode_mem, eweight));
}

int MRIStepGetNumGEvals(void* arkode_mem, long int* ngevals)
{
  return (ARKodeGetNumGEvals(arkode_mem, ngevals));
}

int MRIStepGetRootInfo(void* arkode_mem, int* rootsfound)
{
  return (ARKodeGetRootInfo(arkode_mem, rootsfound));
}

int MRIStepGetUserData(void* arkode_mem, void** user_data)
{
  return (ARKodeGetUserData(arkode_mem, user_data));
}

int MRIStepPrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt)
{
  return (ARKodePrintAllStats(arkode_mem, outfile, fmt));
}

char* MRIStepGetReturnFlagName(long int flag)
{
  return (ARKodeGetReturnFlagName(flag));
}

int MRIStepWriteParameters(void* arkode_mem, FILE* fp)
{
  return (ARKodeWriteParameters(arkode_mem, fp));
}

int MRIStepWriteCoupling(void* arkode_mem, FILE* fp)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeMRIStepMem structures */
  retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* check that coupling structure is non-NULL (otherwise report error) */
  if (step_mem->MRIC == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "Coupling structure is NULL");
    return (ARK_MEM_NULL);
  }

  /* write coupling structure to specified file */
  fprintf(fp, "\nMRIStep coupling structure:\n");
  MRIStepCoupling_Write(step_mem->MRIC, fp);

  return (ARK_SUCCESS);
}

int MRIStepGetNonlinearSystemData(void* arkode_mem, sunrealtype* tcur,
                                  N_Vector* zpred, N_Vector* z, N_Vector* Fi,
                                  sunrealtype* gamma, N_Vector* sdata,
                                  void** user_data)
{
  return (ARKodeGetNonlinearSystemData(arkode_mem, tcur, zpred, z, Fi, gamma,
                                       sdata, user_data));
}

int MRIStepGetNumNonlinSolvIters(void* arkode_mem, long int* nniters)
{
  return (ARKodeGetNumNonlinSolvIters(arkode_mem, nniters));
}

int MRIStepGetNumNonlinSolvConvFails(void* arkode_mem, long int* nnfails)
{
  return (ARKodeGetNumNonlinSolvConvFails(arkode_mem, nnfails));
}

int MRIStepGetNonlinSolvStats(void* arkode_mem, long int* nniters,
                              long int* nnfails)
{
  return (ARKodeGetNonlinSolvStats(arkode_mem, nniters, nnfails));
}

int MRIStepGetNumStepSolveFails(void* arkode_mem, long int* nncfails)
{
  return (ARKodeGetNumStepSolveFails(arkode_mem, nncfails));
}

int MRIStepGetJac(void* arkode_mem, SUNMatrix* J)
{
  return (ARKodeGetJac(arkode_mem, J));
}

int MRIStepGetJacTime(void* arkode_mem, sunrealtype* t_J)
{
  return (ARKodeGetJacTime(arkode_mem, t_J));
}

int MRIStepGetJacNumSteps(void* arkode_mem, long* nst_J)
{
  return (ARKodeGetJacNumSteps(arkode_mem, nst_J));
}

int MRIStepGetLinWorkSpace(void* arkode_mem, long int* lenrwLS, long int* leniwLS)
{
  return (ARKodeGetLinWorkSpace(arkode_mem, lenrwLS, leniwLS));
}

int MRIStepGetNumJacEvals(void* arkode_mem, long int* njevals)
{
  return (ARKodeGetNumJacEvals(arkode_mem, njevals));
}

int MRIStepGetNumPrecEvals(void* arkode_mem, long int* npevals)
{
  return (ARKodeGetNumPrecEvals(arkode_mem, npevals));
}

int MRIStepGetNumPrecSolves(void* arkode_mem, long int* npsolves)
{
  return (ARKodeGetNumPrecSolves(arkode_mem, npsolves));
}

int MRIStepGetNumLinIters(void* arkode_mem, long int* nliters)
{
  return (ARKodeGetNumLinIters(arkode_mem, nliters));
}

int MRIStepGetNumLinConvFails(void* arkode_mem, long int* nlcfails)
{
  return (ARKodeGetNumLinConvFails(arkode_mem, nlcfails));
}

int MRIStepGetNumJTSetupEvals(void* arkode_mem, long int* njtsetups)
{
  return (ARKodeGetNumJTSetupEvals(arkode_mem, njtsetups));
}

int MRIStepGetNumJtimesEvals(void* arkode_mem, long int* njvevals)
{
  return (ARKodeGetNumJtimesEvals(arkode_mem, njvevals));
}

int MRIStepGetNumLinRhsEvals(void* arkode_mem, long int* nfevalsLS)
{
  return (ARKodeGetNumLinRhsEvals(arkode_mem, nfevalsLS));
}

int MRIStepGetLastLinFlag(void* arkode_mem, long int* flag)
{
  return (ARKodeGetLastLinFlag(arkode_mem, flag));
}

char* MRIStepGetLinReturnFlagName(long int flag)
{
  return (ARKodeGetLinReturnFlagName(flag));
}

void MRIStepFree(void** arkode_mem) { ARKodeFree(arkode_mem); }

void MRIStepPrintMem(void* arkode_mem, FILE* outfile)
{
  ARKodePrintMem(arkode_mem, outfile);
}

/*===============================================================
  EOF
  ===============================================================*/
