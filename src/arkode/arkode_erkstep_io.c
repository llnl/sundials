/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
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
 * output functions for the ARKODE ERKStep time stepper module.
 *
 * NOTE: many functions currently in arkode_io.c will move here,
 * with slightly different names.  The code transition will be
 * minimal, but the documentation changes will be significant.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "arkode_erkstep_impl.h"

#include "sundials_cli.h"
#include "sundials_utils.h"

/*===============================================================
  Exported optional input functions.
  ===============================================================*/

/*---------------------------------------------------------------
  ERKStepSetTable:

  Specifies to use a customized Butcher table for the explicit
  portion of the system.

  If d==NULL, then the method is automatically flagged as a
  fixed-step method; a user MUST also call either
  ERKStepSetFixedStep or ERKStepSetInitStep to set the desired
  time step size.
  ---------------------------------------------------------------*/
int ERKStepSetTable(void* arkode_mem, ARKodeButcherTable B)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  sunindextype Blrw, Bliw;
  int retval;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  retval = erkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* check for legal inputs */
  if (B == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }

  /* clear any existing parameters and Butcher tables */
  step_mem->stages = 0;
  step_mem->q      = 0;
  step_mem->p      = 0;

  ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
  ARKodeButcherTable_Free(step_mem->B);
  step_mem->B = NULL;
  ark_mem->liw -= Bliw;
  ark_mem->lrw -= Blrw;

  /* set the relevant parameters */
  step_mem->stages = B->stages;
  step_mem->q      = B->q;
  step_mem->p      = B->p;

  /* copy the table into step memory */
  step_mem->B = ARKodeButcherTable_Copy(B);
  if (step_mem->B == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }

  ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ERKStepSetTableNum:

  Specifies to use a pre-existing Butcher table for the problem,
  based on the integer flag passed to ARKodeButcherTable_LoadERK()
  within the file arkode_butcher_erk.c.
  ---------------------------------------------------------------*/
int ERKStepSetTableNum(void* arkode_mem, ARKODE_ERKTableID etable)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  sunindextype Blrw, Bliw;
  int retval;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  retval = erkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* check that argument specifies an explicit table */
  if (etable < ARKODE_MIN_ERK_NUM || etable > ARKODE_MAX_ERK_NUM)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "Illegal ERK table number");
    return (ARK_ILL_INPUT);
  }

  /* clear any existing parameters and Butcher tables */
  step_mem->stages = 0;
  step_mem->q      = 0;
  step_mem->p      = 0;

  ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
  ARKodeButcherTable_Free(step_mem->B);
  step_mem->B = NULL;
  ark_mem->liw -= Bliw;
  ark_mem->lrw -= Blrw;

  /* fill in table based on argument */
  step_mem->B = ARKodeButcherTable_LoadERK(etable);
  if (step_mem->B == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "Error setting table with that index");
    return (ARK_ILL_INPUT);
  }
  step_mem->stages = step_mem->B->stages;
  step_mem->q      = step_mem->B->q;
  step_mem->p      = step_mem->B->p;

  ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ERKStepSetTableName:

  Specifies to use a pre-existing Butcher table for the problem,
  based on the string passed to ARKodeButcherTable_LoadERKByNmae()
  within the file arkode_butcher_erk.c.
  ---------------------------------------------------------------*/
int ERKStepSetTableName(void* arkode_mem, const char* etable)
{
  return ERKStepSetTableNum(arkode_mem, arkButcherTableERKNameToID(etable));
}

/*===============================================================
  Exported optional output functions.
  ===============================================================*/

/*---------------------------------------------------------------
  erkStep_GetNumRhsEvals:

  Returns the current number of RHS calls
  ---------------------------------------------------------------*/
int erkStep_GetNumRhsEvals(ARKodeMem ark_mem, int partition_index,
                           long int* rhs_evals)
{
  ARKodeERKStepMem step_mem = NULL;

  /* access ARKodeERKStepMem structure */
  int retval = erkStep_AccessStepMem(ark_mem, __func__, &step_mem);
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

int ERKStepGetNumRhsEvals(void* arkode_mem, long int* fevals)
{
  return ARKodeGetNumRhsEvals(arkode_mem, 0, fevals);
}

/*---------------------------------------------------------------
  ERKStepGetCurrentButcherTable:

  Sets pointers to the Butcher table currently in use.
  ---------------------------------------------------------------*/
int ERKStepGetCurrentButcherTable(void* arkode_mem, ARKodeButcherTable* B)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  retval = erkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* get tables from step_mem */
  *B = step_mem->B;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ERKStepGetTimestepperStats:

  Returns integrator statistics
  ---------------------------------------------------------------*/
int ERKStepGetTimestepperStats(void* arkode_mem, long int* expsteps,
                               long int* accsteps, long int* attempts,
                               long int* fevals, long int* netfails)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  retval = erkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set expsteps and accsteps from adaptivity structure */
  *expsteps = ark_mem->hadapt_mem->nst_exp;
  *accsteps = ark_mem->hadapt_mem->nst_acc;

  /* set remaining outputs */
  *attempts = ark_mem->nst_attempts;
  *fevals   = step_mem->nfe;
  *netfails = ark_mem->netf;

  return (ARK_SUCCESS);
}

/*===============================================================
  Private functions attached to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  erkStep_SetOption:

  Provides string-based control over ERKStep-specific "set" routines.
  ---------------------------------------------------------------*/
int erkStep_SetOptions(ARKodeMem ark_mem, int* argidx, char* argv[],
                       size_t offset, sunbooleantype* arg_used)
{
  /* Set lists of keys, and the corresponding set routines */
  static const struct sunKeyCharPair char_pairs[] = {
    {"table_name", ERKStepSetTableName}};
  static const int num_char_keys = sizeof(char_pairs) / sizeof(*char_pairs);

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

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  erkStep_SetRelaxFn:

  Sets up the relaxation module using ERKStep's utility routines.
  ---------------------------------------------------------------*/
int erkStep_SetRelaxFn(ARKodeMem ark_mem, ARKRelaxFn rfn, ARKRelaxJacFn rjac)
{
  return (
    arkRelaxCreate(ark_mem, rfn, rjac, erkStep_RelaxDeltaE, erkStep_GetOrder));
}

/*---------------------------------------------------------------
  erkStep_SetDefaults:

  Resets all ERKStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.
  ---------------------------------------------------------------*/
int erkStep_SetDefaults(ARKodeMem ark_mem)
{
  ARKodeERKStepMem step_mem;
  sunindextype Blrw, Bliw;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Set default values for integrator optional inputs */
  step_mem->q      = Q_DEFAULT;                  /* method order */
  step_mem->p      = 0;                          /* embedding order */
  step_mem->stages = 0;                          /* no stages */
  ark_mem->hadapt_mem->etamxf = SUN_RCONST(0.3); /* max change on error-failed step */
  ark_mem->hadapt_mem->safety = SUN_RCONST(0.99); /* step adaptivity safety factor  */
  ark_mem->hadapt_mem->growth = SUN_RCONST(25.0); /* step adaptivity growth factor */

  /* Remove pre-existing Butcher table */
  if (step_mem->B)
  {
    ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
    ark_mem->liw -= Bliw;
    ark_mem->lrw -= Blrw;
    ARKodeButcherTable_Free(step_mem->B);
  }
  step_mem->B = NULL;

  /* Load the default SUNAdaptController */
  retval = arkReplaceAdaptController(ark_mem, NULL, SUNTRUE);
  if (retval) { return retval; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  erkStep_SetOrder:

  Specifies the method order
  ---------------------------------------------------------------*/
int erkStep_SetOrder(ARKodeMem ark_mem, int ord)
{
  ARKodeERKStepMem step_mem;
  sunindextype Blrw, Bliw;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set user-provided value, or default, depending on argument */
  if (ord <= 0) { step_mem->q = Q_DEFAULT; }
  else { step_mem->q = ord; }

  /* clear Butcher tables, since user is requesting a change in method
     or a reset to defaults.  Tables will be set in ARKInitialSetup. */
  step_mem->stages = 0;
  step_mem->p      = 0;

  ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
  ARKodeButcherTable_Free(step_mem->B);
  step_mem->B = NULL;
  ark_mem->liw -= Bliw;
  ark_mem->lrw -= Blrw;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  erkStep_GetEstLocalErrors: Returns the current local truncation
  error estimate vector
  ---------------------------------------------------------------*/
int erkStep_GetEstLocalErrors(ARKodeMem ark_mem, N_Vector ele)
{
  int retval;
  ARKodeERKStepMem step_mem;
  retval = erkStep_AccessStepMem(ark_mem, __func__, &step_mem);
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
  erkStep_GetStageIndex:

  Returns the current stage index and number of stages
  ---------------------------------------------------------------*/
int erkStep_GetStageIndex(ARKodeMem ark_mem, int* stage, int* max_stages)
{
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  *stage      = step_mem->istage;
  *max_stages = step_mem->stages;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  erkStep_PrintAllStats:

  Prints integrator statistics
  ---------------------------------------------------------------*/
int erkStep_PrintAllStats(ARKodeMem ark_mem, FILE* outfile, SUNOutputFormat fmt)
{
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  sunfprintf_long(outfile, fmt, SUNFALSE, "RHS fn evals", step_mem->nfe);

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  erkStep_WriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int erkStep_WriteParameters(ARKodeMem ark_mem, FILE* fp)
{
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* print integrator parameters to file */
  fprintf(fp, "ERKStep time step module parameters:\n");
  fprintf(fp, "  Method order %i\n", step_mem->q);
  fprintf(fp, "\n");

  return (ARK_SUCCESS);
}

/*===============================================================
  Exported-but-deprecated user-callable functions.
  ===============================================================*/

int ERKStepResize(void* arkode_mem, N_Vector y0, sunrealtype hscale,
                  sunrealtype t0, ARKVecResizeFn resize, void* resize_data)
{
  return (ARKodeResize(arkode_mem, y0, hscale, t0, resize, resize_data));
}

int ERKStepReset(void* arkode_mem, sunrealtype tR, N_Vector yR)
{
  return (ARKodeReset(arkode_mem, tR, yR));
}

int ERKStepSStolerances(void* arkode_mem, sunrealtype reltol, sunrealtype abstol)
{
  return (ARKodeSStolerances(arkode_mem, reltol, abstol));
}

int ERKStepSVtolerances(void* arkode_mem, sunrealtype reltol, N_Vector abstol)
{
  return (ARKodeSVtolerances(arkode_mem, reltol, abstol));
}

int ERKStepWFtolerances(void* arkode_mem, ARKEwtFn efun)
{
  return (ARKodeWFtolerances(arkode_mem, efun));
}

int ERKStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g)
{
  return (ARKodeRootInit(arkode_mem, nrtfn, g));
}

int ERKStepSetDefaults(void* arkode_mem)
{
  return (ARKodeSetDefaults(arkode_mem));
}

int ERKStepSetOrder(void* arkode_mem, int ord)
{
  return (ARKodeSetOrder(arkode_mem, ord));
}

int ERKStepSetInterpolantType(void* arkode_mem, int itype)
{
  return (ARKodeSetInterpolantType(arkode_mem, itype));
}

int ERKStepSetInterpolantDegree(void* arkode_mem, int degree)
{
  return (ARKodeSetInterpolantDegree(arkode_mem, degree));
}

int ERKStepSetDenseOrder(void* arkode_mem, int dord)
{
  return (ARKodeSetInterpolantDegree(arkode_mem, dord));
}

int ERKStepSetAdaptController(void* arkode_mem, SUNAdaptController C)
{
  return (ARKodeSetAdaptController(arkode_mem, C));
}

int ERKStepSetAdaptivityAdjustment(void* arkode_mem, int adjust)
{
  return (ARKodeSetAdaptivityAdjustment(arkode_mem, adjust));
}

int ERKStepSetCFLFraction(void* arkode_mem, sunrealtype cfl_frac)
{
  return (ARKodeSetCFLFraction(arkode_mem, cfl_frac));
}

int ERKStepSetSafetyFactor(void* arkode_mem, sunrealtype safety)
{
  return (ARKodeSetSafetyFactor(arkode_mem, safety));
}

int ERKStepSetErrorBias(void* arkode_mem, sunrealtype bias)
{
  return (ARKodeSetErrorBias(arkode_mem, bias));
}

int ERKStepSetMaxGrowth(void* arkode_mem, sunrealtype mx_growth)
{
  return (ARKodeSetMaxGrowth(arkode_mem, mx_growth));
}

int ERKStepSetMinReduction(void* arkode_mem, sunrealtype eta_min)
{
  return (ARKodeSetMinReduction(arkode_mem, eta_min));
}

int ERKStepSetFixedStepBounds(void* arkode_mem, sunrealtype lb, sunrealtype ub)
{
  return (ARKodeSetFixedStepBounds(arkode_mem, lb, ub));
}

int ERKStepSetAdaptivityMethod(void* arkode_mem, int imethod, int idefault,
                               int pq, sunrealtype adapt_params[3])
{
  return (arkSetAdaptivityMethod(arkode_mem, imethod, idefault, pq, adapt_params));
}

int ERKStepSetAdaptivityFn(void* arkode_mem, ARKAdaptFn hfun, void* h_data)
{
  return (arkSetAdaptivityFn(arkode_mem, hfun, h_data));
}

int ERKStepSetMaxFirstGrowth(void* arkode_mem, sunrealtype etamx1)
{
  return (ARKodeSetMaxFirstGrowth(arkode_mem, etamx1));
}

int ERKStepSetMaxEFailGrowth(void* arkode_mem, sunrealtype etamxf)
{
  return (ARKodeSetMaxEFailGrowth(arkode_mem, etamxf));
}

int ERKStepSetSmallNumEFails(void* arkode_mem, int small_nef)
{
  return (ARKodeSetSmallNumEFails(arkode_mem, small_nef));
}

int ERKStepSetStabilityFn(void* arkode_mem, ARKExpStabFn EStab, void* estab_data)
{
  return (ARKodeSetStabilityFn(arkode_mem, EStab, estab_data));
}

int ERKStepSetMaxErrTestFails(void* arkode_mem, int maxnef)
{
  return (ARKodeSetMaxErrTestFails(arkode_mem, maxnef));
}

int ERKStepSetConstraints(void* arkode_mem, N_Vector constraints)
{
  return (ARKodeSetConstraints(arkode_mem, constraints));
}

int ERKStepSetMaxNumSteps(void* arkode_mem, long int mxsteps)
{
  return (ARKodeSetMaxNumSteps(arkode_mem, mxsteps));
}

int ERKStepSetMaxHnilWarns(void* arkode_mem, int mxhnil)
{
  return (ARKodeSetMaxHnilWarns(arkode_mem, mxhnil));
}

int ERKStepSetInitStep(void* arkode_mem, sunrealtype hin)
{
  return (ARKodeSetInitStep(arkode_mem, hin));
}

int ERKStepSetMinStep(void* arkode_mem, sunrealtype hmin)
{
  return (ARKodeSetMinStep(arkode_mem, hmin));
}

int ERKStepSetMaxStep(void* arkode_mem, sunrealtype hmax)
{
  return (ARKodeSetMaxStep(arkode_mem, hmax));
}

int ERKStepSetInterpolateStopTime(void* arkode_mem, sunbooleantype interp)
{
  return (ARKodeSetInterpolateStopTime(arkode_mem, interp));
}

int ERKStepSetStopTime(void* arkode_mem, sunrealtype tstop)
{
  return (ARKodeSetStopTime(arkode_mem, tstop));
}

int ERKStepClearStopTime(void* arkode_mem)
{
  return (ARKodeClearStopTime(arkode_mem));
}

int ERKStepSetFixedStep(void* arkode_mem, sunrealtype hfixed)
{
  return (ARKodeSetFixedStep(arkode_mem, hfixed));
}

int ERKStepSetMaxNumConstrFails(void* arkode_mem, int maxfails)
{
  return (ARKodeSetMaxNumConstrFails(arkode_mem, maxfails));
}

int ERKStepSetRootDirection(void* arkode_mem, int* rootdir)
{
  return (ARKodeSetRootDirection(arkode_mem, rootdir));
}

int ERKStepSetNoInactiveRootWarn(void* arkode_mem)
{
  return (ARKodeSetNoInactiveRootWarn(arkode_mem));
}

int ERKStepSetUserData(void* arkode_mem, void* user_data)
{
  return (ARKodeSetUserData(arkode_mem, user_data));
}

int ERKStepSetPostprocessStepFn(void* arkode_mem, ARKPostProcessFn ProcessStep)
{
  return (ARKodeSetPostprocessStepFn(arkode_mem, ProcessStep));
}

int ERKStepSetPostprocessStageFn(void* arkode_mem, ARKPostProcessFn ProcessStage)
{
  return (ARKodeSetPostprocessStageFn(arkode_mem, ProcessStage));
}

int ERKStepEvolve(void* arkode_mem, sunrealtype tout, N_Vector yout,
                  sunrealtype* tret, int itask)
{
  return (ARKodeEvolve(arkode_mem, tout, yout, tret, itask));
}

int ERKStepGetDky(void* arkode_mem, sunrealtype t, int k, N_Vector dky)
{
  return (ARKodeGetDky(arkode_mem, t, k, dky));
}

int ERKStepGetNumExpSteps(void* arkode_mem, long int* nsteps)
{
  return (ARKodeGetNumExpSteps(arkode_mem, nsteps));
}

int ERKStepGetNumAccSteps(void* arkode_mem, long int* nsteps)
{
  return (ARKodeGetNumAccSteps(arkode_mem, nsteps));
}

int ERKStepGetNumStepAttempts(void* arkode_mem, long int* nstep_attempts)
{
  return (ARKodeGetNumStepAttempts(arkode_mem, nstep_attempts));
}

int ERKStepGetNumErrTestFails(void* arkode_mem, long int* netfails)
{
  return (ARKodeGetNumErrTestFails(arkode_mem, netfails));
}

int ERKStepGetEstLocalErrors(void* arkode_mem, N_Vector ele)
{
  return (ARKodeGetEstLocalErrors(arkode_mem, ele));
}

int ERKStepGetWorkSpace(void* arkode_mem, long int* lenrw, long int* leniw)
{
  return (ARKodeGetWorkSpace(arkode_mem, lenrw, leniw));
}

int ERKStepGetNumSteps(void* arkode_mem, long int* nsteps)
{
  return (ARKodeGetNumSteps(arkode_mem, nsteps));
}

int ERKStepGetActualInitStep(void* arkode_mem, sunrealtype* hinused)
{
  return (ARKodeGetActualInitStep(arkode_mem, hinused));
}

int ERKStepGetLastStep(void* arkode_mem, sunrealtype* hlast)
{
  return (ARKodeGetLastStep(arkode_mem, hlast));
}

int ERKStepGetCurrentStep(void* arkode_mem, sunrealtype* hcur)
{
  return (ARKodeGetCurrentStep(arkode_mem, hcur));
}

int ERKStepGetCurrentTime(void* arkode_mem, sunrealtype* tcur)
{
  return (ARKodeGetCurrentTime(arkode_mem, tcur));
}

int ERKStepGetTolScaleFactor(void* arkode_mem, sunrealtype* tolsfact)
{
  return (ARKodeGetTolScaleFactor(arkode_mem, tolsfact));
}

int ERKStepGetErrWeights(void* arkode_mem, N_Vector eweight)
{
  return (ARKodeGetErrWeights(arkode_mem, eweight));
}

int ERKStepGetNumGEvals(void* arkode_mem, long int* ngevals)
{
  return (ARKodeGetNumGEvals(arkode_mem, ngevals));
}

int ERKStepGetRootInfo(void* arkode_mem, int* rootsfound)
{
  return (ARKodeGetRootInfo(arkode_mem, rootsfound));
}

int ERKStepGetNumConstrFails(void* arkode_mem, long int* nconstrfails)
{
  return (ARKodeGetNumConstrFails(arkode_mem, nconstrfails));
}

int ERKStepGetUserData(void* arkode_mem, void** user_data)
{
  return (ARKodeGetUserData(arkode_mem, user_data));
}

int ERKStepPrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt)
{
  return (ARKodePrintAllStats(arkode_mem, outfile, fmt));
}

char* ERKStepGetReturnFlagName(long int flag)
{
  return (ARKodeGetReturnFlagName(flag));
}

int ERKStepWriteParameters(void* arkode_mem, FILE* fp)
{
  return (ARKodeWriteParameters(arkode_mem, fp));
}

int ERKStepWriteButcher(void* arkode_mem, FILE* fp)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;

  /* access ARKodeMem and ARKodeERKStepMem structures */
  retval = erkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* check that Butcher table is non-NULL (otherwise report error) */
  if (step_mem->B == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "Butcher table memory is NULL");
    return (ARK_MEM_NULL);
  }

  /* print Butcher table to file */
  fprintf(fp, "\nERKStep Butcher table (stages = %i):\n", step_mem->stages);
  ARKodeButcherTable_Write(step_mem->B, fp);
  fprintf(fp, "\n");

  return (ARK_SUCCESS);
}

int ERKStepGetStepStats(void* arkode_mem, long int* nsteps, sunrealtype* hinused,
                        sunrealtype* hlast, sunrealtype* hcur, sunrealtype* tcur)
{
  return (ARKodeGetStepStats(arkode_mem, nsteps, hinused, hlast, hcur, tcur));
}

void ERKStepFree(void** arkode_mem) { ARKodeFree(arkode_mem); }

void ERKStepPrintMem(void* arkode_mem, FILE* outfile)
{
  ARKodePrintMem(arkode_mem, outfile);
}

int ERKStepSetRelaxFn(void* arkode_mem, ARKRelaxFn rfn, ARKRelaxJacFn rjac)
{
  return (ARKodeSetRelaxFn(arkode_mem, rfn, rjac));
}

int ERKStepSetRelaxEtaFail(void* arkode_mem, sunrealtype eta_rf)
{
  return (ARKodeSetRelaxEtaFail(arkode_mem, eta_rf));
}

int ERKStepSetRelaxLowerBound(void* arkode_mem, sunrealtype lower)
{
  return (ARKodeSetRelaxLowerBound(arkode_mem, lower));
}

int ERKStepSetRelaxMaxFails(void* arkode_mem, int max_fails)
{
  return (ARKodeSetRelaxMaxFails(arkode_mem, max_fails));
}

int ERKStepSetRelaxMaxIters(void* arkode_mem, int max_iters)
{
  return (ARKodeSetRelaxMaxIters(arkode_mem, max_iters));
}

int ERKStepSetRelaxSolver(void* arkode_mem, ARKRelaxSolver solver)
{
  return (ARKodeSetRelaxSolver(arkode_mem, solver));
}

int ERKStepSetRelaxResTol(void* arkode_mem, sunrealtype res_tol)
{
  return (ARKodeSetRelaxResTol(arkode_mem, res_tol));
}

int ERKStepSetRelaxTol(void* arkode_mem, sunrealtype rel_tol, sunrealtype abs_tol)
{
  return (ARKodeSetRelaxTol(arkode_mem, rel_tol, abs_tol));
}

int ERKStepSetRelaxUpperBound(void* arkode_mem, sunrealtype upper)
{
  return (ARKodeSetRelaxUpperBound(arkode_mem, upper));
}

int ERKStepGetNumRelaxFnEvals(void* arkode_mem, long int* r_evals)
{
  return (ARKodeGetNumRelaxFnEvals(arkode_mem, r_evals));
}

int ERKStepGetNumRelaxJacEvals(void* arkode_mem, long int* J_evals)
{
  return (ARKodeGetNumRelaxJacEvals(arkode_mem, J_evals));
}

int ERKStepGetNumRelaxFails(void* arkode_mem, long int* relax_fails)
{
  return (ARKodeGetNumRelaxFails(arkode_mem, relax_fails));
}

int ERKStepGetNumRelaxBoundFails(void* arkode_mem, long int* fails)
{
  return (ARKodeGetNumRelaxBoundFails(arkode_mem, fails));
}

int ERKStepGetNumRelaxSolveFails(void* arkode_mem, long int* fails)
{
  return (ARKodeGetNumRelaxSolveFails(arkode_mem, fails));
}

int ERKStepGetNumRelaxSolveIters(void* arkode_mem, long int* iters)
{
  return (ARKodeGetNumRelaxSolveIters(arkode_mem, iters));
}

/*===============================================================
  EOF
  ===============================================================*/
