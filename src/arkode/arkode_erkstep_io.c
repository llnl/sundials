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
  ARKodeSetFixedStep or ARKodeSetInitStep to set the desired
  time step size.
  ---------------------------------------------------------------*/
int ERKStepSetTable(void* arkode_mem, ARKodeButcherTable B)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
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

  if (step_mem->B != NULL)
  {
    ARKodeButcherTable_Free(step_mem->B);
    step_mem->B = NULL;
  }

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

  if (step_mem->B != NULL)
  {
    ARKodeButcherTable_Free(step_mem->B);
    step_mem->B = NULL;
  }

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
  if (step_mem->B) { ARKodeButcherTable_Free(step_mem->B); }
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

  if (step_mem->B != NULL)
  {
    ARKodeButcherTable_Free(step_mem->B);
    step_mem->B = NULL;
  }

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
  EOF
  ===============================================================*/
