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
 * output functions for the ARKODE ARKStep time stepper module.
 *
 * NOTE: many functions currently in arkode_io.c will move here,
 * with slightly different names.  The code transition will be
 * minimal, but the documentation changes will be significant.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "arkode_arkstep_impl.h"

#include "sundials_cli.h"
#include "sundials_utils.h"

/*===============================================================
  Exported optional input functions.
  ===============================================================*/

/*---------------------------------------------------------------
  ARKStepSetExplicit:

  Specifies that the implicit portion of the problem is disabled,
  and to use an explicit RK method.
  ---------------------------------------------------------------*/
int ARKStepSetExplicit(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  retval = arkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* ensure that fe is defined */
  if (step_mem->fe == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_MISSING_FE);
    return (ARK_ILL_INPUT);
  }

  /* set the relevant parameters */
  step_mem->explicit = SUNTRUE;
  step_mem->implicit = SUNFALSE;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKStepSetImplicit:

  Specifies that the explicit portion of the problem is disabled,
  and to use an implicit RK method.
  ---------------------------------------------------------------*/
int ARKStepSetImplicit(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  retval = arkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* ensure that fi is defined */
  if (step_mem->fi == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_MISSING_FI);
    return (ARK_ILL_INPUT);
  }

  /* set the relevant parameters */
  step_mem->implicit = SUNTRUE;
  step_mem->explicit = SUNFALSE;

  /* re-attach internal error weight functions if necessary */
  if (!ark_mem->user_efun)
  {
    if (ark_mem->itol == ARK_SV && ark_mem->Vabstol != NULL)
    {
      retval = ARKodeSVtolerances(ark_mem, ark_mem->reltol, ark_mem->Vabstol);
    }
    else
    {
      retval = ARKodeSStolerances(ark_mem, ark_mem->reltol, ark_mem->Sabstol);
    }
    if (retval != ARK_SUCCESS) { return (retval); }
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKStepSetImEx:

  Specifies that the specifies that problem has both implicit and
  explicit parts, and to use an ARK method (this is the default).
  ---------------------------------------------------------------*/
int ARKStepSetImEx(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  retval = arkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* ensure that fe and fi are defined */
  if (step_mem->fe == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_MISSING_FE);
    return (ARK_ILL_INPUT);
  }
  if (step_mem->fi == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_MISSING_FI);
    return (ARK_ILL_INPUT);
  }

  /* set the relevant parameters */
  step_mem->explicit = SUNTRUE;
  step_mem->implicit = SUNTRUE;

  /* re-attach internal error weight functions if necessary */
  if (!ark_mem->user_efun)
  {
    if (ark_mem->itol == ARK_SV && ark_mem->Vabstol != NULL)
    {
      retval = ARKodeSVtolerances(ark_mem, ark_mem->reltol, ark_mem->Vabstol);
    }
    else
    {
      retval = ARKodeSStolerances(ark_mem, ark_mem->reltol, ark_mem->Sabstol);
    }
    if (retval != ARK_SUCCESS) { return (retval); }
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKStepSetTables:

  Specifies to use customized Butcher tables for the system.

  If Bi is NULL, then this sets the integrator in 'explicit' mode.

  If Be is NULL, then this sets the integrator in 'implicit' mode.

  Returns ARK_ILL_INPUT if both Butcher tables are not supplied.
  ---------------------------------------------------------------*/
int ARKStepSetTables(void* arkode_mem, int q, int p, ARKodeButcherTable Bi,
                     ARKodeButcherTable Be)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  retval = arkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* check for illegal inputs */
  if ((Bi == NULL) && (Be == NULL))
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "At least one complete table must be supplied");
    return (ARK_ILL_INPUT);
  }

  /* if both tables are set, check that they have the same number of stages */
  if ((Bi != NULL) && (Be != NULL))
  {
    if (Bi->stages != Be->stages)
    {
      arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                      "Both tables must have the same number of stages");
      return (ARK_ILL_INPUT);
    }
  }

  /* clear any existing parameters and Butcher tables */
  step_mem->stages = 0;
  step_mem->q      = 0;
  step_mem->p      = 0;

  if (step_mem->Be != NULL)
  {
    ARKodeButcherTable_Free(step_mem->Be);
    step_mem->Be = NULL;
  }

  if (step_mem->Bi != NULL)
  {
    ARKodeButcherTable_Free(step_mem->Bi);
    step_mem->Bi = NULL;
  }

  /*
   * determine mode (implicit/explicit/ImEx), and perform appropriate actions
   */

  /* explicit */
  if (Bi == NULL)
  {
    /* set the relevant parameters (use table q and p) */
    step_mem->stages = Be->stages;
    step_mem->q      = Be->q;
    step_mem->p      = Be->p;

    /* copy the table in step memory */
    step_mem->Be = ARKodeButcherTable_Copy(Be);
    if (step_mem->Be == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                      MSG_ARK_NO_MEM);
      return (ARK_MEM_NULL);
    }

    /* set method as purely explicit */
    retval = ARKStepSetExplicit(arkode_mem);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Error in ARKStepSetExplicit");
      return (retval);
    }

    /* implicit */
  }
  else if (Be == NULL)
  {
    /* set the relevant parameters (use table q and p) */
    step_mem->stages = Bi->stages;
    step_mem->q      = Bi->q;
    step_mem->p      = Bi->p;

    /* copy the table in step memory */
    step_mem->Bi = ARKodeButcherTable_Copy(Bi);
    if (step_mem->Bi == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                      MSG_ARK_NO_MEM);
      return (ARK_MEM_NULL);
    }

    /* set method as purely implicit */
    retval = ARKStepSetImplicit(arkode_mem);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Error in ARKStepSetImplicit");
      return (ARK_ILL_INPUT);
    }

    /* ImEx */
  }
  else
  {
    /* set the relevant parameters (use input q and p) */
    step_mem->stages = Bi->stages;
    step_mem->q      = q;
    step_mem->p      = p;

    /* copy the explicit table into step memory */
    step_mem->Be = ARKodeButcherTable_Copy(Be);
    if (step_mem->Be == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                      MSG_ARK_NO_MEM);
      return (ARK_MEM_NULL);
    }

    /* copy the implicit table into step memory */
    step_mem->Bi = ARKodeButcherTable_Copy(Bi);
    if (step_mem->Bi == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                      MSG_ARK_NO_MEM);
      return (ARK_MEM_NULL);
    }

    /* set method as ImEx */
    retval = ARKStepSetImEx(arkode_mem);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Error in ARKStepSetImEx");
      return (ARK_ILL_INPUT);
    }
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKStepSetTableNum:

  Specifies to use pre-existing Butcher tables for the system,
  based on the integer flags passed to
  ARKodeButcherTable_LoadERK() and ARKodeButcherTable_LoadDIRK()
  within the files arkode_butcher_erk.c and arkode_butcher_dirk.c
  (automatically calls ARKStepSetImEx).

  If either argument is negative (illegal), then this disables the
  corresponding table (e.g. itable = -1  ->  explicit)

  Note: this routine should NOT be used in conjunction with
  ARKodeSetOrder.
  ---------------------------------------------------------------*/
int ARKStepSetTableNum(void* arkode_mem, ARKODE_DIRKTableID itable,
                       ARKODE_ERKTableID etable)
{
  int flag, retval;
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  retval = arkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* clear any existing parameters and Butcher tables */
  step_mem->stages = 0;
  step_mem->q      = 0;
  step_mem->p      = 0;

  if (step_mem->Be != NULL)
  {
    ARKodeButcherTable_Free(step_mem->Be);
    step_mem->Be = NULL;
  }

  if (step_mem->Bi != NULL)
  {
    ARKodeButcherTable_Free(step_mem->Bi);
    step_mem->Bi = NULL;
  }

  /* determine mode (implicit/explicit/ImEx), and perform
     appropriate actions  */

  /*     illegal inputs */
  if ((itable < 0) && (etable < 0))
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "At least one valid table number must be supplied");
    return (ARK_ILL_INPUT);

    /* explicit */
  }
  else if (itable < 0)
  {
    /* check that argument specifies an explicit table */
    if (etable < ARKODE_MIN_ERK_NUM || etable > ARKODE_MAX_ERK_NUM)
    {
      arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                      "Illegal ERK table number");
      return (ARK_ILL_INPUT);
    }

    /* fill in table based on argument */
    step_mem->Be = ARKodeButcherTable_LoadERK(etable);
    if (step_mem->Be == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                      "Error setting explicit table with that index");
      return (ARK_ILL_INPUT);
    }
    step_mem->stages = step_mem->Be->stages;
    step_mem->q      = step_mem->Be->q;
    step_mem->p      = step_mem->Be->p;

    /* set method as purely explicit */
    flag = ARKStepSetExplicit(arkode_mem);
    if (flag != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Error in ARKStepSetExplicit");
      return (flag);
    }

    /* implicit */
  }
  else if (etable < 0)
  {
    /* check that argument specifies an implicit table */
    if (itable < ARKODE_MIN_DIRK_NUM || itable > ARKODE_MAX_DIRK_NUM)
    {
      arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                      "Illegal IRK table number");
      return (ARK_ILL_INPUT);
    }

    /* fill in table based on argument */
    step_mem->Bi = ARKodeButcherTable_LoadDIRK(itable);
    if (step_mem->Bi == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                      "Error setting table with that index");
      return (ARK_ILL_INPUT);
    }
    step_mem->stages = step_mem->Bi->stages;
    step_mem->q      = step_mem->Bi->q;
    step_mem->p      = step_mem->Bi->p;

    /* set method as purely implicit */
    flag = ARKStepSetImplicit(arkode_mem);
    if (flag != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Error in ARKStepSetImplicit");
      return (flag);
    }

    /* ImEx */
  }
  else
  {
    /* ensure that tables match */
    if (!((etable == ARKODE_ARK324L2SA_ERK_4_2_3) &&
          (itable == ARKODE_ARK324L2SA_DIRK_4_2_3)) &&
        !((etable == ARKODE_ARK436L2SA_ERK_6_3_4) &&
          (itable == ARKODE_ARK436L2SA_DIRK_6_3_4)) &&
        !((etable == ARKODE_ARK437L2SA_ERK_7_3_4) &&
          (itable == ARKODE_ARK437L2SA_DIRK_7_3_4)) &&
        !((etable == ARKODE_ARK548L2SA_ERK_8_4_5) &&
          (itable == ARKODE_ARK548L2SA_DIRK_8_4_5)) &&
        !((etable == ARKODE_ARK548L2SAb_ERK_8_4_5) &&
          (itable == ARKODE_ARK548L2SAb_DIRK_8_4_5)) &&
        !((etable == ARKODE_GKC21_ERK_3_1_2) &&
          (itable == ARKODE_GKC21_ESDIRK_3_1_2)) &&
        !((etable == ARKODE_ARS222_ERK_3_1_2) &&
          (itable == ARKODE_ARS222_ESDIRK_3_1_2)) &&
        /*New Embedded IMEX-SSP Methods*/
        !((etable == ARKODE_SSP_ERK_3_1_2) && (itable == ARKODE_SSP_DIRK_3_1_2)) &&
        !((etable == ARKODE_SSP_LSPUM_ERK_3_1_2) &&
          (itable == ARKODE_SSP_LSPUM_SDIRK_3_1_2)) &&
        !((etable == ARKODE_SSP_ERK_4_2_3) && (itable == ARKODE_ESDIRK_4_2_3)))
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Incompatible Butcher tables for ARK method");
      return (ARK_ILL_INPUT);
    }

    /* fill in tables based on arguments */
    step_mem->Bi = ARKodeButcherTable_LoadDIRK(itable);
    step_mem->Be = ARKodeButcherTable_LoadERK(etable);
    if (step_mem->Bi == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                      "Illegal IRK table number");
      return (ARK_ILL_INPUT);
    }
    if (step_mem->Be == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                      "Illegal ERK table number");
      return (ARK_ILL_INPUT);
    }
    step_mem->stages = step_mem->Bi->stages;
    step_mem->q      = step_mem->Bi->q;
    step_mem->p      = step_mem->Bi->p;

    /* set method as ImEx */
    if (ARKStepSetImEx(arkode_mem) != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      MSG_ARK_MISSING_F);
      return (ARK_ILL_INPUT);
    }
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKStepSetTableName:

  Specifies to use pre-existing Butcher tables for the system,
  based on the string passed to
  ARKodeButcherTable_LoadERKByName() and
  ARKodeButcherTable_LoadDIRKByName() within the files
  arkode_butcher_erk.c and arkode_butcher_dirk.c (automatically
  calls ARKStepSetImEx).

  If itable is "ARKODE_DIRK_NONE" or etable is "ARKODE_ERK_NONE",
  then this disables the corresponding table.
  ---------------------------------------------------------------*/
int ARKStepSetTableName(void* arkode_mem, const char* itable, const char* etable)
{
  return (ARKStepSetTableNum(arkode_mem, arkButcherTableDIRKNameToID(itable),
                             arkButcherTableERKNameToID(etable)));
}

/*===============================================================
  Exported optional output functions.
  ===============================================================*/

/*---------------------------------------------------------------
  arkStep_GetNumRhsEvals:

  Returns the current number of calls
  ---------------------------------------------------------------*/
int arkStep_GetNumRhsEvals(ARKodeMem ark_mem, int partition_index,
                           long int* rhs_evals)
{
  ARKodeARKStepMem step_mem = NULL;

  /* access ARKodeARKStepMem structure */
  int retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
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
  case 0: *rhs_evals = step_mem->nfe; break;
  case 1: *rhs_evals = step_mem->nfi; break;
  default: *rhs_evals = step_mem->nfe + step_mem->nfi; break;
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  ARKStepGetCurrentButcherTables:

  Sets pointers to the explicit and implicit Butcher tables
  currently in use.
  ---------------------------------------------------------------*/
int ARKStepGetCurrentButcherTables(void* arkode_mem, ARKodeButcherTable* Bi,
                                   ARKodeButcherTable* Be)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  retval = arkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* get tables from step_mem */
  *Bi = step_mem->Bi;
  *Be = step_mem->Be;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKStepGetTimestepperStats:

  Returns integrator statistics
  ---------------------------------------------------------------*/
int ARKStepGetTimestepperStats(void* arkode_mem, long int* expsteps,
                               long int* accsteps, long int* step_attempts,
                               long int* fe_evals, long int* fi_evals,
                               long int* nlinsetups, long int* netfails)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  retval = arkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set expsteps and accsteps from adaptivity structure */
  *expsteps = ark_mem->hadapt_mem->nst_exp;
  *accsteps = ark_mem->hadapt_mem->nst_acc;

  /* set remaining outputs */
  *step_attempts = ark_mem->nst_attempts;
  *fe_evals      = step_mem->nfe;
  *fi_evals      = step_mem->nfi;
  *nlinsetups    = step_mem->nsetups;
  *netfails      = ark_mem->netf;

  return (ARK_SUCCESS);
}

/*===============================================================
  Private functions attached to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  arkStep_SetOption:

  Provides string-based control over ARKStep-specific "set" routines.
  ---------------------------------------------------------------*/
int arkStep_SetOptions(ARKodeMem ark_mem, int* argidx, char* argv[],
                       size_t offset, sunbooleantype* arg_used)
{
  /* Set lists of keys, and the corresponding set routines */
  static const struct sunKeyTwoCharPair twochar_pairs[] = {
    {"table_names", ARKStepSetTableName}};
  static const int num_twochar_keys = sizeof(twochar_pairs) /
                                      sizeof(*twochar_pairs);

  static const struct sunKeyActionPair action_pairs[] =
    {{"explicit", ARKStepSetExplicit},
     {"implicit", ARKStepSetImplicit},
     {"imex", ARKStepSetImEx}};
  static const int num_action_keys = sizeof(action_pairs) / sizeof(*action_pairs);

  /* check all "twochar" keys */
  int j, retval;
  retval = sunCheckAndSetTwoCharArgs((void*)ark_mem, argidx, argv, offset,
                                     twochar_pairs, num_twochar_keys, arg_used,
                                     &j);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "error setting command-line argument: %s",
                    twochar_pairs[j].key);
    return retval;
  }
  if (*arg_used) { return ARK_SUCCESS; }

  /* check all action keys */
  retval = sunCheckAndSetActionArgs((void*)ark_mem, argidx, argv, offset,
                                    action_pairs, num_action_keys, arg_used, &j);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "error setting command-line argument: %s",
                    action_pairs[j].key);
    return retval;
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_SetRelaxFn:

  Sets up the relaxation module using ARKStep's utility routines.
  ---------------------------------------------------------------*/
int arkStep_SetRelaxFn(ARKodeMem ark_mem, ARKRelaxFn rfn, ARKRelaxJacFn rjac)
{
  return (
    arkRelaxCreate(ark_mem, rfn, rjac, arkStep_RelaxDeltaE, arkStep_GetOrder));
}

/*---------------------------------------------------------------
  arkStep_SetUserData:

  Passes user-data pointer to attached linear solver modules.
  ---------------------------------------------------------------*/
int arkStep_SetUserData(ARKodeMem ark_mem, void* user_data)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set user data in ARKODE LS mem */
  if (step_mem->lmem != NULL)
  {
    retval = arkLSSetUserData(ark_mem, user_data);
    if (retval != ARKLS_SUCCESS) { return (retval); }
  }

  /* set user data in ARKODE LSMass mem */
  if (step_mem->mass_mem != NULL)
  {
    retval = arkLSSetMassUserData(ark_mem, user_data);
    if (retval != ARKLS_SUCCESS) { return (retval); }
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_SetDefaults:

  Resets all ARKStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.  Also leaves alone any data
  structures/options related to the ARKODE infrastructure itself
  (e.g., root-finding and post-process step).
  ---------------------------------------------------------------*/
int arkStep_SetDefaults(ARKodeMem ark_mem)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Set default values for integrator optional inputs */
  step_mem->q              = Q_DEFAULT; /* method order */
  step_mem->p              = 0;         /* embedding order */
  step_mem->predictor      = 0;         /* trivial predictor */
  step_mem->linear         = SUNFALSE;  /* nonlinear problem */
  step_mem->linear_timedep = SUNTRUE;   /* dfi/dy depends on t */
  step_mem->autonomous     = SUNFALSE;  /* non-autonomous problem */
  step_mem->explicit       = SUNTRUE;   /* fe(t,y) will be used */
  step_mem->implicit       = SUNTRUE;   /* fi(t,y) will be used */
  step_mem->deduce_rhs     = SUNFALSE;  /* deduce fi on result of NLS */
  step_mem->maxcor         = MAXCOR;    /* max nonlinear iters/stage */
  step_mem->nlscoef        = NLSCOEF;   /* nonlinear tolerance coefficient */
  step_mem->crdown         = CRDOWN; /* nonlinear convergence estimate coeff. */
  step_mem->rdiv           = RDIV;   /* nonlinear divergence tolerance */
  step_mem->dgmax    = DGMAX; /* max step change before recomputing J or P */
  step_mem->msbp     = MSBP;  /* max steps between updates to J or P */
  step_mem->stages   = 0;     /* no stages */
  step_mem->istage   = 0;     /* current stage */
  step_mem->jcur     = SUNFALSE;
  step_mem->convfail = ARK_NO_FAILURES;
  step_mem->stage_predict = NULL; /* no user-supplied stage predictor */

  /* Remove pre-existing Butcher tables */
  if (step_mem->Be) { ARKodeButcherTable_Free(step_mem->Be); }
  step_mem->Be = NULL;
  if (step_mem->Bi) { ARKodeButcherTable_Free(step_mem->Bi); }
  step_mem->Bi = NULL;

  /* Remove pre-existing nonlinear solver object */
  if (step_mem->NLS && step_mem->ownNLS) { SUNNonlinSolFree(step_mem->NLS); }
  step_mem->NLS = NULL;

  /* Load the default SUNAdaptController */
  retval = arkReplaceAdaptController(ark_mem, NULL, SUNTRUE);
  if (retval) { return retval; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_SetOrder:

  Specifies the method order
  ---------------------------------------------------------------*/
int arkStep_SetOrder(ARKodeMem ark_mem, int ord)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set user-provided value, or default, depending on argument */
  if (ord <= 0) { step_mem->q = Q_DEFAULT; }
  else { step_mem->q = ord; }

  /* clear Butcher tables, since user is requesting a change in method
     or a reset to defaults.  Tables will be set in ARKInitialSetup. */
  step_mem->stages = 0;
  step_mem->istage = 0;
  step_mem->p      = 0;

  if (step_mem->Be != NULL)
  {
    ARKodeButcherTable_Free(step_mem->Be);
    step_mem->Be = NULL;
  }

  if (step_mem->Bi != NULL)
  {
    ARKodeButcherTable_Free(step_mem->Bi);
    step_mem->Bi = NULL;
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_SetLinear:

  Specifies that the implicit portion of the problem is linear,
  and to tighten the linear solver tolerances while taking only
  one Newton iteration.  DO NOT USE IN COMBINATION WITH THE
  FIXED-POINT SOLVER.  Automatically tightens DeltaGammaMax
  to ensure that step size changes cause Jacobian recomputation.

  The argument should be 1 or 0, where 1 indicates that the
  Jacobian of fi with respect to y depends on time, and
  0 indicates that it is not time dependent.  Alternately, when
  using an iterative linear solver this flag denotes time
  dependence of the preconditioner.
  ---------------------------------------------------------------*/
int arkStep_SetLinear(ARKodeMem ark_mem, int timedepend)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  if (timedepend && step_mem->autonomous)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Incompatible settings, the problem is autonomous but the "
                    "Jacobian is time dependent");
    return ARK_ILL_INPUT;
  }

  /* set parameters */
  step_mem->linear         = SUNTRUE;
  step_mem->linear_timedep = (timedepend == 1);
  step_mem->dgmax          = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_SetNonlinear:

  Specifies that the implicit portion of the problem is nonlinear.
  Used to undo a previous call to arkStep_SetLinear.  Automatically
  loosens DeltaGammaMax back to default value.
  ---------------------------------------------------------------*/
int arkStep_SetNonlinear(ARKodeMem ark_mem)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set parameters */
  step_mem->linear         = SUNFALSE;
  step_mem->linear_timedep = SUNTRUE;
  step_mem->dgmax          = DGMAX;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_SetAutonomous:

  Indicates if the problem is autonomous (True) or non-autonomous
  (False).
  ---------------------------------------------------------------*/
int arkStep_SetAutonomous(ARKodeMem ark_mem, sunbooleantype autonomous)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  step_mem->autonomous = autonomous;

  if (autonomous && step_mem->linear) { step_mem->linear_timedep = SUNFALSE; }

  /* Reattach the nonlinear system function e.g., switching to/from an
     autonomous problem with the trivial predictor requires swapping the
     nonlinear system function provided to the nonlinear solver */
  retval = arkStep_SetNlsSysFn(ark_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Setting nonlinear system function failed");
    return ARK_ILL_INPUT;
  }

  /* This will be better handled when the temp vector stack is added */
  if (autonomous)
  {
    /* Allocate tempv5 if needed */
    if (!arkAllocVec(ark_mem, ark_mem->yn, &ark_mem->tempv5))
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_MEM_FAIL);
      return ARK_MEM_FAIL;
    }
  }
  else
  {
    /* Free tempv5 if necessary */
    arkFreeVec(&ark_mem->tempv5);
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  arkStep_SetNonlinCRDown:

  Specifies the user-provided nonlinear convergence constant
  crdown.  Legal values are strictly positive; illegal values
  imply a reset to the default.
  ---------------------------------------------------------------*/
int arkStep_SetNonlinCRDown(ARKodeMem ark_mem, sunrealtype crdown)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* if argument legal set it, otherwise set default */
  if (crdown <= ZERO) { step_mem->crdown = CRDOWN; }
  else { step_mem->crdown = crdown; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_SetNonlinRDiv:

  Specifies the user-provided nonlinear convergence constant
  rdiv.  Legal values are strictly positive; illegal values
  imply a reset to the default.
  ---------------------------------------------------------------*/
int arkStep_SetNonlinRDiv(ARKodeMem ark_mem, sunrealtype rdiv)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* if argument legal set it, otherwise set default */
  if (rdiv <= ZERO) { step_mem->rdiv = RDIV; }
  else { step_mem->rdiv = rdiv; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_SetDeltaGammaMax:

  Specifies the user-provided linear setup decision constant
  dgmax.  Legal values are strictly positive; illegal values imply
  a reset to the default.
  ---------------------------------------------------------------*/
int arkStep_SetDeltaGammaMax(ARKodeMem ark_mem, sunrealtype dgmax)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* if argument legal set it, otherwise set default */
  if (dgmax <= ZERO) { step_mem->dgmax = DGMAX; }
  else { step_mem->dgmax = dgmax; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_SetLSetupFrequency:

  Specifies the user-provided linear setup decision constant
  msbp.  Positive values give the frequency for calling lsetup;
  negative values imply recomputation of lsetup at each nonlinear
  solve; a zero value implies a reset to the default.
  ---------------------------------------------------------------*/
int arkStep_SetLSetupFrequency(ARKodeMem ark_mem, int msbp)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* if argument legal set it, otherwise set default */
  if (msbp == 0) { step_mem->msbp = MSBP; }
  else { step_mem->msbp = msbp; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_SetPredictorMethod:

  Specifies the method to use for predicting implicit solutions.
  Non-default choices are {1,2,3,4}, all others will use default
  (trivial) predictor.
  ---------------------------------------------------------------*/
int arkStep_SetPredictorMethod(ARKodeMem ark_mem, int pred_method)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* set parameter */
  step_mem->predictor = pred_method;

  /* Reattach the nonlinear system function e.g., switching to/from the trivial
     predictor with an autonomous problem requires swapping the nonlinear system
     function provided to the nonlinear solver */
  retval = arkStep_SetNlsSysFn(ark_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Setting nonlinear system function failed");
    return ARK_ILL_INPUT;
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_SetMaxNonlinIters:

  Specifies the maximum number of nonlinear iterations during
  one solve.  A non-positive input implies a reset to the
  default value.
  ---------------------------------------------------------------*/
int arkStep_SetMaxNonlinIters(ARKodeMem ark_mem, int maxcor)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
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
  arkStep_SetNonlinConvCoef:

  Specifies the coefficient in the nonlinear solver convergence
  test.  A non-positive input implies a reset to the default value.
  ---------------------------------------------------------------*/
int arkStep_SetNonlinConvCoef(ARKodeMem ark_mem, sunrealtype nlscoef)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* argument <= 0 sets default, otherwise set input */
  if (nlscoef <= ZERO) { step_mem->nlscoef = NLSCOEF; }
  else { step_mem->nlscoef = nlscoef; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_SetStagePredictFn:  Specifies a user-provided step
  predictor function having type ARKStagePredictFn.  A
  NULL input function disables calls to this routine.
  ---------------------------------------------------------------*/
int arkStep_SetStagePredictFn(ARKodeMem ark_mem, ARKStagePredictFn PredictStage)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure and set function pointer */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  step_mem->stage_predict = PredictStage;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_SetDeduceImplicitRhs:

  Specifies if an optimization is used to avoid an evaluation of
  fi after a nonlinear solve for an implicit stage.  If stage
  postprocessecing in enabled, this option is ignored, and fi is
  never deduced.

  An argument of SUNTRUE indicates that fi is deduced to compute
  fi(z_i), and SUNFALSE indicates that fi(z_i) is computed with
  an additional evaluation of fi.
  ---------------------------------------------------------------*/
int arkStep_SetDeduceImplicitRhs(ARKodeMem ark_mem, sunbooleantype deduce)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure and set function pointer */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  step_mem->deduce_rhs = deduce;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_GetCurrentGamma: Returns the current value of gamma
  ---------------------------------------------------------------*/
int arkStep_GetCurrentGamma(ARKodeMem ark_mem, sunrealtype* gamma)
{
  int retval;
  ARKodeARKStepMem step_mem;
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }
  *gamma = step_mem->gamma;
  return (retval);
}

/*---------------------------------------------------------------
  arkStep_GetEstLocalErrors: Returns the current local truncation
  error estimate vector
  ---------------------------------------------------------------*/
int arkStep_GetEstLocalErrors(ARKodeMem ark_mem, N_Vector ele)
{
  int retval;
  ARKodeARKStepMem step_mem;
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
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
  arkStep_GetNumLinSolvSetups:

  Returns the current number of calls to the lsetup routine
  ---------------------------------------------------------------*/
int arkStep_GetNumLinSolvSetups(ARKodeMem ark_mem, long int* nlinsetups)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* get value from step_mem */
  *nlinsetups = step_mem->nsetups;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_GetNumNonlinSolvIters:

  Returns the current number of nonlinear solver iterations
  ---------------------------------------------------------------*/
int arkStep_GetNumNonlinSolvIters(ARKodeMem ark_mem, long int* nniters)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  *nniters = step_mem->nls_iters;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_GetNumNonlinSolvConvFails:

  Returns the current number of nonlinear solver convergence fails
  ---------------------------------------------------------------*/
int arkStep_GetNumNonlinSolvConvFails(ARKodeMem ark_mem, long int* nnfails)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
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
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  *nniters = step_mem->nls_iters;
  *nnfails = step_mem->nls_fails;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_GetStageIndex:

  Returns the current stage index and number of stages
  ---------------------------------------------------------------*/
int arkStep_GetStageIndex(ARKodeMem ark_mem, int* stage, int* max_stages)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  *stage      = step_mem->istage;
  *max_stages = step_mem->stages;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_PrintAllStats:

  Prints integrator statistics
  ---------------------------------------------------------------*/
int arkStep_PrintAllStats(ARKodeMem ark_mem, FILE* outfile, SUNOutputFormat fmt)
{
  ARKodeARKStepMem step_mem;
  ARKLsMem arkls_mem;
  ARKLsMassMem arklsm_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* function evaluations */
  sunfprintf_long(outfile, fmt, SUNFALSE, "Explicit RHS fn evals", step_mem->nfe);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Implicit RHS fn evals", step_mem->nfi);

  /* nonlinear solver stats */
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

  /* mass solve stats */
  if (ark_mem->step_getmassmem(ark_mem))
  {
    arklsm_mem = (ARKLsMassMem)(ark_mem->step_getmassmem(ark_mem));
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass setups", arklsm_mem->nmsetups);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass solves", arklsm_mem->nmsolves);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass Prec setup evals",
                    arklsm_mem->npe);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass Prec solves", arklsm_mem->nps);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass LS iters", arklsm_mem->nli);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass LS fails", arklsm_mem->ncfl);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass-times setups",
                    arklsm_mem->nmtsetup);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Mass-times evals",
                    arklsm_mem->nmtimes);
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStep_WriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int arkStep_WriteParameters(ARKodeMem ark_mem, FILE* fp)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* print integrator parameters to file */
  fprintf(fp, "ARKStep time step module parameters:\n");
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
  if (step_mem->explicit && step_mem->implicit)
  {
    fprintf(fp, "  ImEx integrator\n");
  }
  else if (step_mem->implicit) { fprintf(fp, "  Implicit integrator\n"); }
  else { fprintf(fp, "  Explicit integrator\n"); }

  if (step_mem->implicit)
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
  EOF
  ===============================================================*/
