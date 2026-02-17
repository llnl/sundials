/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
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
  sunindextype Blrw, Bliw;

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

  ARKodeButcherTable_Space(step_mem->Be, &Bliw, &Blrw);
  ARKodeButcherTable_Free(step_mem->Be);
  step_mem->Be = NULL;
  ark_mem->liw -= Bliw;
  ark_mem->lrw -= Blrw;

  ARKodeButcherTable_Space(step_mem->Bi, &Bliw, &Blrw);
  ARKodeButcherTable_Free(step_mem->Bi);
  step_mem->Bi = NULL;
  ark_mem->liw -= Bliw;
  ark_mem->lrw -= Blrw;

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

  /* note Butcher table space requirements */
  ARKodeButcherTable_Space(step_mem->Be, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;

  ARKodeButcherTable_Space(step_mem->Bi, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;

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
  sunindextype Blrw, Bliw;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  retval = arkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* clear any existing parameters and Butcher tables */
  step_mem->stages = 0;
  step_mem->q      = 0;
  step_mem->p      = 0;

  ARKodeButcherTable_Space(step_mem->Be, &Bliw, &Blrw);
  ARKodeButcherTable_Free(step_mem->Be);
  step_mem->Be = NULL;
  ark_mem->liw -= Bliw;
  ark_mem->lrw -= Blrw;

  ARKodeButcherTable_Space(step_mem->Bi, &Bliw, &Blrw);
  ARKodeButcherTable_Free(step_mem->Bi);
  step_mem->Bi = NULL;
  ark_mem->liw -= Bliw;
  ark_mem->lrw -= Blrw;

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
        !((etable == ARKODE_ARK2_ERK_3_1_2) && (itable == ARKODE_ARK2_DIRK_3_1_2)) &&
        /*New Embedded IMEX-SSP Methods*/
        !((etable == ARKODE_SSP_ERK_2_1_2) &&
          (itable == ARKODE_SSP_SDIRK_2_1_2)) &&
        !((etable == ARKODE_SSP_ERK_3_1_2) &&
          (itable == ARKODE_SSP_DIRK_3_1_2)) &&
        !((etable == ARKODE_SSP_LSPUM_ERK_3_1_2) &&
          (itable == ARKODE_SSP_LSPUM_SDIRK_3_1_2)) &&
        !((etable == ARKODE_SSP_ERK_4_2_3) &&
          (itable == ARKODE_SSP_ESDIRK_4_2_3)) &&
        !((etable == ARKODE_SSP_ERK_9_2_3) &&
          (itable == ARKODE_SSP_ESDIRK_9_2_3)))
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

    /* hard code the bias factor for  ARKODE_SSP_LSPUM_3_1_2 */
    if ((etable == ARKODE_SSP_LSPUM_ERK_3_1_2) &&
        (itable == ARKODE_SSP_LSPUM_SDIRK_3_1_2))
    {
        (void)SUNAdaptController_SetErrorBias(ark_mem->hadapt_mem->hcontroller,
                                              SUN_RCONST(2.0));
    }
    /* END hard code the bias factor for  ARKODE_SSP_LSPUM_3_1_2 */

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

int ARKStepGetNumRhsEvals(void* arkode_mem, long int* fe_evals, long int* fi_evals)
{
  int retval = ARK_SUCCESS;

  retval = ARKodeGetNumRhsEvals(arkode_mem, 0, fe_evals);
  if (retval != ARK_SUCCESS) { return retval; }

  retval = ARKodeGetNumRhsEvals(arkode_mem, 1, fi_evals);
  if (retval != ARK_SUCCESS) { return retval; }

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
  sunindextype Blrw, Bliw;
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
  if (step_mem->Be)
  {
    ARKodeButcherTable_Space(step_mem->Be, &Bliw, &Blrw);
    ark_mem->liw -= Bliw;
    ark_mem->lrw -= Blrw;
    ARKodeButcherTable_Free(step_mem->Be);
  }
  step_mem->Be = NULL;
  if (step_mem->Bi)
  {
    ARKodeButcherTable_Space(step_mem->Bi, &Bliw, &Blrw);
    ark_mem->liw -= Bliw;
    ark_mem->lrw -= Blrw;
    ARKodeButcherTable_Free(step_mem->Bi);
  }
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
  sunindextype Blrw, Bliw;
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

  ARKodeButcherTable_Space(step_mem->Be, &Bliw, &Blrw);
  ARKodeButcherTable_Free(step_mem->Be);
  step_mem->Be = NULL;
  ark_mem->liw -= Bliw;
  ark_mem->lrw -= Blrw;

  ARKodeButcherTable_Space(step_mem->Bi, &Bliw, &Blrw);
  ARKodeButcherTable_Free(step_mem->Bi);
  step_mem->Bi = NULL;
  ark_mem->liw -= Bliw;
  ark_mem->lrw -= Blrw;

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
    arkFreeVec(ark_mem, &ark_mem->tempv5);
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
  Exported-but-deprecated user-callable functions.
  ===============================================================*/

int ARKStepCreateMRIStepInnerStepper(void* inner_arkode_mem,
                                     MRIStepInnerStepper* stepper)
{
  return (ARKodeCreateMRIStepInnerStepper(inner_arkode_mem, stepper));
}

int ARKStepResize(void* arkode_mem, N_Vector y0, sunrealtype hscale,
                  sunrealtype t0, ARKVecResizeFn resize, void* resize_data)
{
  return (ARKodeResize(arkode_mem, y0, hscale, t0, resize, resize_data));
}

int ARKStepReset(void* arkode_mem, sunrealtype tR, N_Vector yR)
{
  return (ARKodeReset(arkode_mem, tR, yR));
}

int ARKStepSStolerances(void* arkode_mem, sunrealtype reltol, sunrealtype abstol)
{
  return (ARKodeSStolerances(arkode_mem, reltol, abstol));
}

int ARKStepSVtolerances(void* arkode_mem, sunrealtype reltol, N_Vector abstol)
{
  return (ARKodeSVtolerances(arkode_mem, reltol, abstol));
}

int ARKStepWFtolerances(void* arkode_mem, ARKEwtFn efun)
{
  return (ARKodeWFtolerances(arkode_mem, efun));
}

int ARKStepResStolerance(void* arkode_mem, sunrealtype rabstol)
{
  return (ARKodeResStolerance(arkode_mem, rabstol));
}

int ARKStepResVtolerance(void* arkode_mem, N_Vector rabstol)
{
  return (ARKodeResVtolerance(arkode_mem, rabstol));
}

int ARKStepResFtolerance(void* arkode_mem, ARKRwtFn rfun)
{
  return (ARKodeResFtolerance(arkode_mem, rfun));
}

int ARKStepSetLinearSolver(void* arkode_mem, SUNLinearSolver LS, SUNMatrix A)
{
  return (ARKodeSetLinearSolver(arkode_mem, LS, A));
}

int ARKStepSetMassLinearSolver(void* arkode_mem, SUNLinearSolver LS,
                               SUNMatrix M, sunbooleantype time_dep)
{
  return (ARKodeSetMassLinearSolver(arkode_mem, LS, M, time_dep));
}

int ARKStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g)
{
  return (ARKodeRootInit(arkode_mem, nrtfn, g));
}

int ARKStepSetDefaults(void* arkode_mem)
{
  return (ARKodeSetDefaults(arkode_mem));
}

int ARKStepSetOptimalParams(void* arkode_mem)
{
  /* TODO: do we need to do something here? This is deprecated with no
   * ARKodeSetOptimalParams to replace it */
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;
  long int lenrw, leniw;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  retval = arkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* access ARKodeHAdaptMem structure */
  if (ark_mem->hadapt_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARKADAPT_NO_MEM);
    return (ARK_MEM_NULL);
  }
  hadapt_mem = ark_mem->hadapt_mem;

  /* Remove current SUNAdaptController object */
  retval = SUNAdaptController_Space(hadapt_mem->hcontroller, &lenrw, &leniw);
  if (retval == SUN_SUCCESS)
  {
    ark_mem->liw -= leniw;
    ark_mem->lrw -= lenrw;
  }
  if (hadapt_mem->owncontroller)
  {
    retval = SUNAdaptController_Destroy(hadapt_mem->hcontroller);
    ark_mem->hadapt_mem->owncontroller = SUNFALSE;
    if (retval != SUN_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "SUNAdaptController_Destroy failure");
      return (ARK_MEM_FAIL);
    }
  }
  hadapt_mem->hcontroller = NULL;

  /* Choose values based on method, order */

  /*    explicit */
  if (step_mem->explicit && !step_mem->implicit)
  {
    hadapt_mem->hcontroller = SUNAdaptController_PI(ark_mem->sunctx);
    if (hadapt_mem->hcontroller == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "SUNAdaptController_PI allocation failure");
      return (ARK_MEM_FAIL);
    }
    (void)SUNAdaptController_SetErrorBias(hadapt_mem->hcontroller,
                                          SUN_RCONST(1.2));
    (void)SUNAdaptController_SetParams_PI(hadapt_mem->hcontroller,
                                          SUN_RCONST(0.8), -SUN_RCONST(0.31));
    hadapt_mem->safety = SUN_RCONST(0.99);
    hadapt_mem->growth = SUN_RCONST(25.0);
    hadapt_mem->etamxf = SUN_RCONST(0.3);
    hadapt_mem->pq     = PQ;

    /*    implicit */
  }
  else if (step_mem->implicit && !step_mem->explicit)
  {
    switch (step_mem->q)
    {
    case 2: /* just use standard defaults since better ones unknown */
      hadapt_mem->hcontroller = SUNAdaptController_PID(ark_mem->sunctx);
      if (hadapt_mem->hcontroller == NULL)
      {
        arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                        "SUNAdaptController_PID allocation failure");
        return (ARK_MEM_FAIL);
      }
      hadapt_mem->safety    = SAFETY;
      hadapt_mem->growth    = GROWTH;
      hadapt_mem->etamxf    = ETAMXF;
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      hadapt_mem->pq        = PQ;
      step_mem->nlscoef     = SUN_RCONST(0.001);
      step_mem->maxcor      = 5;
      step_mem->crdown      = CRDOWN;
      step_mem->rdiv        = RDIV;
      step_mem->dgmax       = DGMAX;
      step_mem->msbp        = MSBP;
      break;
    case 3:
      hadapt_mem->hcontroller = SUNAdaptController_I(ark_mem->sunctx);
      if (hadapt_mem->hcontroller == NULL)
      {
        arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                        "SUNAdaptController_I allocation failure");
        return (ARK_MEM_FAIL);
      }
      (void)SUNAdaptController_SetErrorBias(hadapt_mem->hcontroller,
                                            SUN_RCONST(1.9));
      hadapt_mem->safety    = SUN_RCONST(0.957);
      hadapt_mem->growth    = SUN_RCONST(17.6);
      hadapt_mem->etamxf    = SUN_RCONST(0.45);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      hadapt_mem->pq        = PQ;
      step_mem->nlscoef     = SUN_RCONST(0.22);
      step_mem->crdown      = SUN_RCONST(0.17);
      step_mem->rdiv        = SUN_RCONST(2.3);
      step_mem->dgmax       = SUN_RCONST(0.19);
      step_mem->msbp        = 60;
      break;
    case 4:
      hadapt_mem->hcontroller = SUNAdaptController_PID(ark_mem->sunctx);
      if (hadapt_mem->hcontroller == NULL)
      {
        arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                        "SUNAdaptController_PID allocation failure");
        return (ARK_MEM_FAIL);
      }
      (void)SUNAdaptController_SetErrorBias(hadapt_mem->hcontroller,
                                            SUN_RCONST(1.2));
      (void)SUNAdaptController_SetParams_PID(hadapt_mem->hcontroller,
                                             SUN_RCONST(0.535),
                                             -SUN_RCONST(0.209),
                                             SUN_RCONST(0.148));
      hadapt_mem->safety    = SUN_RCONST(0.988);
      hadapt_mem->growth    = SUN_RCONST(31.5);
      hadapt_mem->etamxf    = SUN_RCONST(0.33);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      hadapt_mem->pq        = PQ;
      step_mem->nlscoef     = SUN_RCONST(0.24);
      step_mem->crdown      = SUN_RCONST(0.26);
      step_mem->rdiv        = SUN_RCONST(2.3);
      step_mem->dgmax       = SUN_RCONST(0.16);
      step_mem->msbp        = 31;
      break;
    case 5:
      hadapt_mem->hcontroller = SUNAdaptController_PID(ark_mem->sunctx);
      if (hadapt_mem->hcontroller == NULL)
      {
        arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                        "SUNAdaptController_PID allocation failure");
        return (ARK_MEM_FAIL);
      }
      (void)SUNAdaptController_SetErrorBias(hadapt_mem->hcontroller,
                                            SUN_RCONST(3.3));
      (void)SUNAdaptController_SetParams_PID(hadapt_mem->hcontroller,
                                             SUN_RCONST(0.56), -SUN_RCONST(0.338),
                                             SUN_RCONST(0.14));
      hadapt_mem->safety    = SUN_RCONST(0.937);
      hadapt_mem->growth    = SUN_RCONST(22.0);
      hadapt_mem->etamxf    = SUN_RCONST(0.44);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      hadapt_mem->pq        = PQ;
      step_mem->nlscoef     = SUN_RCONST(0.25);
      step_mem->crdown      = SUN_RCONST(0.4);
      step_mem->rdiv        = SUN_RCONST(2.3);
      step_mem->dgmax       = SUN_RCONST(0.32);
      step_mem->msbp        = 31;
      break;
    }

    /*    imex */
  }
  else
  {
    switch (step_mem->q)
    {
    case 2: /* just use standard defaults since better ones unknown */
      hadapt_mem->hcontroller = SUNAdaptController_PID(ark_mem->sunctx);
      if (hadapt_mem->hcontroller == NULL)
      {
        arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                        "SUNAdaptController_PID allocation failure");
        return (ARK_MEM_FAIL);
      }
      hadapt_mem->safety    = SAFETY;
      hadapt_mem->growth    = GROWTH;
      hadapt_mem->etamxf    = ETAMXF;
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      hadapt_mem->pq        = PQ;
      step_mem->nlscoef     = SUN_RCONST(0.001);
      step_mem->maxcor      = 5;
      step_mem->crdown      = CRDOWN;
      step_mem->rdiv        = RDIV;
      step_mem->dgmax       = DGMAX;
      step_mem->msbp        = MSBP;
      break;
    case 3:
      hadapt_mem->hcontroller = SUNAdaptController_PID(ark_mem->sunctx);
      if (hadapt_mem->hcontroller == NULL)
      {
        arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                        "SUNAdaptController_PID allocation failure");
        return (ARK_MEM_FAIL);
      }
      (void)SUNAdaptController_SetErrorBias(hadapt_mem->hcontroller,
                                            SUN_RCONST(1.42));
      (void)SUNAdaptController_SetParams_PID(hadapt_mem->hcontroller,
                                             SUN_RCONST(0.54), -SUN_RCONST(0.36),
                                             SUN_RCONST(0.14));
      hadapt_mem->safety    = SUN_RCONST(0.965);
      hadapt_mem->growth    = SUN_RCONST(28.7);
      hadapt_mem->etamxf    = SUN_RCONST(0.46);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      hadapt_mem->pq        = PQ;
      step_mem->nlscoef     = SUN_RCONST(0.22);
      step_mem->crdown      = SUN_RCONST(0.17);
      step_mem->rdiv        = SUN_RCONST(2.3);
      step_mem->dgmax       = SUN_RCONST(0.19);
      step_mem->msbp        = 60;
      break;
    case 4:
      hadapt_mem->hcontroller = SUNAdaptController_PID(ark_mem->sunctx);
      if (hadapt_mem->hcontroller == NULL)
      {
        arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                        "SUNAdaptController_PID allocation failure");
        return (ARK_MEM_FAIL);
      }
      (void)SUNAdaptController_SetErrorBias(hadapt_mem->hcontroller,
                                            SUN_RCONST(1.35));
      (void)SUNAdaptController_SetParams_PID(hadapt_mem->hcontroller,
                                             SUN_RCONST(0.543),
                                             -SUN_RCONST(0.297),
                                             SUN_RCONST(0.14));
      hadapt_mem->safety    = SUN_RCONST(0.97);
      hadapt_mem->growth    = SUN_RCONST(25.0);
      hadapt_mem->etamxf    = SUN_RCONST(0.47);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      hadapt_mem->pq        = PQ;
      step_mem->nlscoef     = SUN_RCONST(0.24);
      step_mem->crdown      = SUN_RCONST(0.26);
      step_mem->rdiv        = SUN_RCONST(2.3);
      step_mem->dgmax       = SUN_RCONST(0.16);
      step_mem->msbp        = 31;
      break;
    case 5:
      hadapt_mem->hcontroller = SUNAdaptController_PI(ark_mem->sunctx);
      if (hadapt_mem->hcontroller == NULL)
      {
        arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                        "SUNAdaptController_PI allocation failure");
        return (ARK_MEM_FAIL);
      }
      (void)SUNAdaptController_SetErrorBias(hadapt_mem->hcontroller,
                                            SUN_RCONST(1.15));
      (void)SUNAdaptController_SetParams_PI(hadapt_mem->hcontroller,
                                            SUN_RCONST(0.8), -SUN_RCONST(0.35));
      hadapt_mem->safety    = SUN_RCONST(0.993);
      hadapt_mem->growth    = SUN_RCONST(28.5);
      hadapt_mem->etamxf    = SUN_RCONST(0.3);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      hadapt_mem->pq        = PQ;
      step_mem->nlscoef     = SUN_RCONST(0.25);
      step_mem->crdown      = SUN_RCONST(0.4);
      step_mem->rdiv        = SUN_RCONST(2.3);
      step_mem->dgmax       = SUN_RCONST(0.32);
      step_mem->msbp        = 31;
      break;
    }
    hadapt_mem->owncontroller = SUNTRUE;

    retval = SUNAdaptController_Space(hadapt_mem->hcontroller, &lenrw, &leniw);
    if (retval == SUN_SUCCESS)
    {
      ark_mem->liw += leniw;
      ark_mem->lrw += lenrw;
    }
  }
  return (ARK_SUCCESS);
}

int ARKStepSetOrder(void* arkode_mem, int ord)
{
  return (ARKodeSetOrder(arkode_mem, ord));
}

int ARKStepSetInterpolantType(void* arkode_mem, int itype)
{
  return (ARKodeSetInterpolantType(arkode_mem, itype));
}

int ARKStepSetInterpolantDegree(void* arkode_mem, int degree)
{
  return (ARKodeSetInterpolantDegree(arkode_mem, degree));
}

int ARKStepSetDenseOrder(void* arkode_mem, int dord)
{
  return (ARKodeSetInterpolantDegree(arkode_mem, dord));
}

int ARKStepSetNonlinearSolver(void* arkode_mem, SUNNonlinearSolver NLS)
{
  return (ARKodeSetNonlinearSolver(arkode_mem, NLS));
}

int ARKStepSetNlsRhsFn(void* arkode_mem, ARKRhsFn nls_fi)
{
  return (ARKodeSetNlsRhsFn(arkode_mem, nls_fi));
}

int ARKStepSetLinear(void* arkode_mem, int timedepend)
{
  return (ARKodeSetLinear(arkode_mem, timedepend));
}

int ARKStepSetNonlinear(void* arkode_mem)
{
  return (ARKodeSetNonlinear(arkode_mem));
}

int ARKStepSetDeduceImplicitRhs(void* arkode_mem, sunbooleantype deduce)
{
  return (ARKodeSetDeduceImplicitRhs(arkode_mem, deduce));
}

int ARKStepSetAdaptController(void* arkode_mem, SUNAdaptController C)
{
  return (ARKodeSetAdaptController(arkode_mem, C));
}

int ARKStepSetAdaptivityAdjustment(void* arkode_mem, int adjust)
{
  return (ARKodeSetAdaptivityAdjustment(arkode_mem, adjust));
}

int ARKStepSetCFLFraction(void* arkode_mem, sunrealtype cfl_frac)
{
  return (ARKodeSetCFLFraction(arkode_mem, cfl_frac));
}

int ARKStepSetSafetyFactor(void* arkode_mem, sunrealtype safety)
{
  return (ARKodeSetSafetyFactor(arkode_mem, safety));
}

int ARKStepSetErrorBias(void* arkode_mem, sunrealtype bias)
{
  return (ARKodeSetErrorBias(arkode_mem, bias));
}

int ARKStepSetMaxGrowth(void* arkode_mem, sunrealtype mx_growth)
{
  return (ARKodeSetMaxGrowth(arkode_mem, mx_growth));
}

int ARKStepSetMinReduction(void* arkode_mem, sunrealtype eta_min)
{
  return (ARKodeSetMinReduction(arkode_mem, eta_min));
}

int ARKStepSetFixedStepBounds(void* arkode_mem, sunrealtype lb, sunrealtype ub)
{
  return (ARKodeSetFixedStepBounds(arkode_mem, lb, ub));
}

int ARKStepSetAdaptivityMethod(void* arkode_mem, int imethod, int idefault,
                               int pq, sunrealtype adapt_params[3])
{
  return (arkSetAdaptivityMethod(arkode_mem, imethod, idefault, pq, adapt_params));
}

int ARKStepSetAdaptivityFn(void* arkode_mem, ARKAdaptFn hfun, void* h_data)
{
  return (arkSetAdaptivityFn(arkode_mem, hfun, h_data));
}

int ARKStepSetMaxFirstGrowth(void* arkode_mem, sunrealtype etamx1)
{
  return (ARKodeSetMaxFirstGrowth(arkode_mem, etamx1));
}

int ARKStepSetMaxEFailGrowth(void* arkode_mem, sunrealtype etamxf)
{
  return (ARKodeSetMaxEFailGrowth(arkode_mem, etamxf));
}

int ARKStepSetSmallNumEFails(void* arkode_mem, int small_nef)
{
  return (ARKodeSetSmallNumEFails(arkode_mem, small_nef));
}

int ARKStepSetMaxCFailGrowth(void* arkode_mem, sunrealtype etacf)
{
  return (ARKodeSetMaxCFailGrowth(arkode_mem, etacf));
}

int ARKStepSetNonlinCRDown(void* arkode_mem, sunrealtype crdown)
{
  return (ARKodeSetNonlinCRDown(arkode_mem, crdown));
}

int ARKStepSetNonlinRDiv(void* arkode_mem, sunrealtype rdiv)
{
  return (ARKodeSetNonlinRDiv(arkode_mem, rdiv));
}

int ARKStepSetDeltaGammaMax(void* arkode_mem, sunrealtype dgmax)
{
  return (ARKodeSetDeltaGammaMax(arkode_mem, dgmax));
}

int ARKStepSetLSetupFrequency(void* arkode_mem, int msbp)
{
  return (ARKodeSetLSetupFrequency(arkode_mem, msbp));
}

int ARKStepSetPredictorMethod(void* arkode_mem, int pred_method)
{
  return (ARKodeSetPredictorMethod(arkode_mem, pred_method));
}

int ARKStepSetStabilityFn(void* arkode_mem, ARKExpStabFn EStab, void* estab_data)
{
  return (ARKodeSetStabilityFn(arkode_mem, EStab, estab_data));
}

int ARKStepSetMaxErrTestFails(void* arkode_mem, int maxnef)
{
  return (ARKodeSetMaxErrTestFails(arkode_mem, maxnef));
}

int ARKStepSetMaxNonlinIters(void* arkode_mem, int maxcor)
{
  return (ARKodeSetMaxNonlinIters(arkode_mem, maxcor));
}

int ARKStepSetMaxConvFails(void* arkode_mem, int maxncf)
{
  return (ARKodeSetMaxConvFails(arkode_mem, maxncf));
}

int ARKStepSetNonlinConvCoef(void* arkode_mem, sunrealtype nlscoef)
{
  return (ARKodeSetNonlinConvCoef(arkode_mem, nlscoef));
}

int ARKStepSetConstraints(void* arkode_mem, N_Vector constraints)
{
  return (ARKodeSetConstraints(arkode_mem, constraints));
}

int ARKStepSetMaxNumSteps(void* arkode_mem, long int mxsteps)
{
  return (ARKodeSetMaxNumSteps(arkode_mem, mxsteps));
}

int ARKStepSetMaxHnilWarns(void* arkode_mem, int mxhnil)
{
  return (ARKodeSetMaxHnilWarns(arkode_mem, mxhnil));
}

int ARKStepSetInitStep(void* arkode_mem, sunrealtype hin)
{
  return (ARKodeSetInitStep(arkode_mem, hin));
}

int ARKStepSetMinStep(void* arkode_mem, sunrealtype hmin)
{
  return (ARKodeSetMinStep(arkode_mem, hmin));
}

int ARKStepSetMaxStep(void* arkode_mem, sunrealtype hmax)
{
  return (ARKodeSetMaxStep(arkode_mem, hmax));
}

int ARKStepSetInterpolateStopTime(void* arkode_mem, sunbooleantype interp)
{
  return (ARKodeSetInterpolateStopTime(arkode_mem, interp));
}

int ARKStepSetStopTime(void* arkode_mem, sunrealtype tstop)
{
  return (ARKodeSetStopTime(arkode_mem, tstop));
}

int ARKStepClearStopTime(void* arkode_mem)
{
  return (ARKodeClearStopTime(arkode_mem));
}

int ARKStepSetFixedStep(void* arkode_mem, sunrealtype hfixed)
{
  return (ARKodeSetFixedStep(arkode_mem, hfixed));
}

int ARKStepSetMaxNumConstrFails(void* arkode_mem, int maxfails)
{
  return (ARKodeSetMaxNumConstrFails(arkode_mem, maxfails));
}

int ARKStepSetRootDirection(void* arkode_mem, int* rootdir)
{
  return (ARKodeSetRootDirection(arkode_mem, rootdir));
}

int ARKStepSetNoInactiveRootWarn(void* arkode_mem)
{
  return (ARKodeSetNoInactiveRootWarn(arkode_mem));
}

int ARKStepSetUserData(void* arkode_mem, void* user_data)
{
  return (ARKodeSetUserData(arkode_mem, user_data));
}

int ARKStepSetPostprocessStepFn(void* arkode_mem, ARKPostProcessFn ProcessStep)
{
  return (ARKodeSetPostprocessStepFn(arkode_mem, ProcessStep));
}

int ARKStepSetPostprocessStageFn(void* arkode_mem, ARKPostProcessFn ProcessStage)
{
  return (ARKodeSetPostprocessStageFn(arkode_mem, ProcessStage));
}

int ARKStepSetStagePredictFn(void* arkode_mem, ARKStagePredictFn PredictStage)
{
  return (ARKodeSetStagePredictFn(arkode_mem, PredictStage));
}

int ARKStepSetJacFn(void* arkode_mem, ARKLsJacFn jac)
{
  return (ARKodeSetJacFn(arkode_mem, jac));
}

int ARKStepSetMassFn(void* arkode_mem, ARKLsMassFn mass)
{
  return (ARKodeSetMassFn(arkode_mem, mass));
}

int ARKStepSetJacEvalFrequency(void* arkode_mem, long int msbj)
{
  return (ARKodeSetJacEvalFrequency(arkode_mem, msbj));
}

int ARKStepSetLinearSolutionScaling(void* arkode_mem, sunbooleantype onoff)
{
  return (ARKodeSetLinearSolutionScaling(arkode_mem, onoff));
}

int ARKStepSetEpsLin(void* arkode_mem, sunrealtype eplifac)
{
  return (ARKodeSetEpsLin(arkode_mem, eplifac));
}

int ARKStepSetMassEpsLin(void* arkode_mem, sunrealtype eplifac)
{
  return (ARKodeSetMassEpsLin(arkode_mem, eplifac));
}

int ARKStepSetLSNormFactor(void* arkode_mem, sunrealtype nrmfac)
{
  return (ARKodeSetLSNormFactor(arkode_mem, nrmfac));
}

int ARKStepSetMassLSNormFactor(void* arkode_mem, sunrealtype nrmfac)
{
  return (ARKodeSetMassLSNormFactor(arkode_mem, nrmfac));
}

int ARKStepSetPreconditioner(void* arkode_mem, ARKLsPrecSetupFn psetup,
                             ARKLsPrecSolveFn psolve)
{
  return (ARKodeSetPreconditioner(arkode_mem, psetup, psolve));
}

int ARKStepSetMassPreconditioner(void* arkode_mem, ARKLsMassPrecSetupFn psetup,
                                 ARKLsMassPrecSolveFn psolve)
{
  return (ARKodeSetMassPreconditioner(arkode_mem, psetup, psolve));
}

int ARKStepSetJacTimes(void* arkode_mem, ARKLsJacTimesSetupFn jtsetup,
                       ARKLsJacTimesVecFn jtimes)
{
  return (ARKodeSetJacTimes(arkode_mem, jtsetup, jtimes));
}

int ARKStepSetJacTimesRhsFn(void* arkode_mem, ARKRhsFn jtimesRhsFn)
{
  return (ARKodeSetJacTimesRhsFn(arkode_mem, jtimesRhsFn));
}

int ARKStepSetMassTimes(void* arkode_mem, ARKLsMassTimesSetupFn msetup,
                        ARKLsMassTimesVecFn mtimes, void* mtimes_data)
{
  return (ARKodeSetMassTimes(arkode_mem, msetup, mtimes, mtimes_data));
}

int ARKStepSetLinSysFn(void* arkode_mem, ARKLsLinSysFn linsys)
{
  return (ARKodeSetLinSysFn(arkode_mem, linsys));
}

int ARKStepEvolve(void* arkode_mem, sunrealtype tout, N_Vector yout,
                  sunrealtype* tret, int itask)
{
  return (ARKodeEvolve(arkode_mem, tout, yout, tret, itask));
}

int ARKStepGetDky(void* arkode_mem, sunrealtype t, int k, N_Vector dky)
{
  return (ARKodeGetDky(arkode_mem, t, k, dky));
}

int ARKStepComputeState(void* arkode_mem, N_Vector zcor, N_Vector z)
{
  return (ARKodeComputeState(arkode_mem, zcor, z));
}

int ARKStepGetNumExpSteps(void* arkode_mem, long int* nsteps)
{
  return (ARKodeGetNumExpSteps(arkode_mem, nsteps));
}

int ARKStepGetNumAccSteps(void* arkode_mem, long int* nsteps)
{
  return (ARKodeGetNumAccSteps(arkode_mem, nsteps));
}

int ARKStepGetNumStepAttempts(void* arkode_mem, long int* nstep_attempts)
{
  return (ARKodeGetNumStepAttempts(arkode_mem, nstep_attempts));
}

int ARKStepGetNumLinSolvSetups(void* arkode_mem, long int* nlinsetups)
{
  return (ARKodeGetNumLinSolvSetups(arkode_mem, nlinsetups));
}

int ARKStepGetNumErrTestFails(void* arkode_mem, long int* netfails)
{
  return (ARKodeGetNumErrTestFails(arkode_mem, netfails));
}

int ARKStepGetEstLocalErrors(void* arkode_mem, N_Vector ele)
{
  return (ARKodeGetEstLocalErrors(arkode_mem, ele));
}

int ARKStepGetWorkSpace(void* arkode_mem, long int* lenrw, long int* leniw)
{
  return (ARKodeGetWorkSpace(arkode_mem, lenrw, leniw));
}

int ARKStepGetNumSteps(void* arkode_mem, long int* nsteps)
{
  return (ARKodeGetNumSteps(arkode_mem, nsteps));
}

int ARKStepGetActualInitStep(void* arkode_mem, sunrealtype* hinused)
{
  return (ARKodeGetActualInitStep(arkode_mem, hinused));
}

int ARKStepGetLastStep(void* arkode_mem, sunrealtype* hlast)
{
  return (ARKodeGetLastStep(arkode_mem, hlast));
}

int ARKStepGetCurrentStep(void* arkode_mem, sunrealtype* hcur)
{
  return (ARKodeGetCurrentStep(arkode_mem, hcur));
}

int ARKStepGetCurrentTime(void* arkode_mem, sunrealtype* tcur)
{
  return (ARKodeGetCurrentTime(arkode_mem, tcur));
}

int ARKStepGetCurrentState(void* arkode_mem, N_Vector* state)
{
  return (ARKodeGetCurrentState(arkode_mem, state));
}

int ARKStepGetCurrentGamma(void* arkode_mem, sunrealtype* gamma)
{
  return (ARKodeGetCurrentGamma(arkode_mem, gamma));
}

int ARKStepGetCurrentMassMatrix(void* arkode_mem, SUNMatrix* M)
{
  return (ARKodeGetCurrentMassMatrix(arkode_mem, M));
}

int ARKStepGetTolScaleFactor(void* arkode_mem, sunrealtype* tolsfact)
{
  return (ARKodeGetTolScaleFactor(arkode_mem, tolsfact));
}

int ARKStepGetErrWeights(void* arkode_mem, N_Vector eweight)
{
  return (ARKodeGetErrWeights(arkode_mem, eweight));
}

int ARKStepGetResWeights(void* arkode_mem, N_Vector rweight)
{
  return (ARKodeGetResWeights(arkode_mem, rweight));
}

int ARKStepGetNumGEvals(void* arkode_mem, long int* ngevals)
{
  return (ARKodeGetNumGEvals(arkode_mem, ngevals));
}

int ARKStepGetRootInfo(void* arkode_mem, int* rootsfound)
{
  return (ARKodeGetRootInfo(arkode_mem, rootsfound));
}

int ARKStepGetNumConstrFails(void* arkode_mem, long int* nconstrfails)
{
  return (ARKodeGetNumConstrFails(arkode_mem, nconstrfails));
}

int ARKStepGetUserData(void* arkode_mem, void** user_data)
{
  return (ARKodeGetUserData(arkode_mem, user_data));
}

int ARKStepPrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt)
{
  return (ARKodePrintAllStats(arkode_mem, outfile, fmt));
}

char* ARKStepGetReturnFlagName(long int flag)
{
  return (ARKodeGetReturnFlagName(flag));
}

int ARKStepWriteParameters(void* arkode_mem, FILE* fp)
{
  return (ARKodeWriteParameters(arkode_mem, fp));
}

int ARKStepWriteButcher(void* arkode_mem, FILE* fp)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeMem and ARKodeARKStepMem structures */
  retval = arkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* check that Butcher table is non-NULL (otherwise report error) */
  if ((step_mem->Be == NULL) && (step_mem->Bi == NULL))
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "Butcher table memory is NULL");
    return (ARK_MEM_NULL);
  }

  /* print Butcher tables to file */
  fprintf(fp, "\nARKStep Butcher tables (stages = %i):\n", step_mem->stages);
  if (step_mem->explicit && (step_mem->Be != NULL))
  {
    fprintf(fp, "  Explicit Butcher table:\n");
    ARKodeButcherTable_Write(step_mem->Be, fp);
  }
  fprintf(fp, "\n");
  if (step_mem->implicit && (step_mem->Bi != NULL))
  {
    fprintf(fp, "  Implicit Butcher table:\n");
    ARKodeButcherTable_Write(step_mem->Bi, fp);
  }
  fprintf(fp, "\n");

  return (ARK_SUCCESS);
}

int ARKStepGetStepStats(void* arkode_mem, long int* nsteps, sunrealtype* hinused,
                        sunrealtype* hlast, sunrealtype* hcur, sunrealtype* tcur)
{
  return (ARKodeGetStepStats(arkode_mem, nsteps, hinused, hlast, hcur, tcur));
}

int ARKStepGetNonlinearSystemData(void* arkode_mem, sunrealtype* tcur,
                                  N_Vector* zpred, N_Vector* z, N_Vector* Fi,
                                  sunrealtype* gamma, N_Vector* sdata,
                                  void** user_data)
{
  return (ARKodeGetNonlinearSystemData(arkode_mem, tcur, zpred, z, Fi, gamma,
                                       sdata, user_data));
}

int ARKStepGetNumNonlinSolvIters(void* arkode_mem, long int* nniters)
{
  return (ARKodeGetNumNonlinSolvIters(arkode_mem, nniters));
}

int ARKStepGetNumNonlinSolvConvFails(void* arkode_mem, long int* nnfails)
{
  return (ARKodeGetNumNonlinSolvConvFails(arkode_mem, nnfails));
}

int ARKStepGetNonlinSolvStats(void* arkode_mem, long int* nniters,
                              long int* nnfails)
{
  return (ARKodeGetNonlinSolvStats(arkode_mem, nniters, nnfails));
}

int ARKStepGetNumStepSolveFails(void* arkode_mem, long int* nncfails)
{
  return (ARKodeGetNumStepSolveFails(arkode_mem, nncfails));
}

int ARKStepGetJac(void* arkode_mem, SUNMatrix* J)
{
  return (ARKodeGetJac(arkode_mem, J));
}

int ARKStepGetJacTime(void* arkode_mem, sunrealtype* t_J)
{
  return (ARKodeGetJacTime(arkode_mem, t_J));
}

int ARKStepGetJacNumSteps(void* arkode_mem, long* nst_J)
{
  return (ARKodeGetJacNumSteps(arkode_mem, nst_J));
}

int ARKStepGetLinWorkSpace(void* arkode_mem, long int* lenrwLS, long int* leniwLS)
{
  return (ARKodeGetLinWorkSpace(arkode_mem, lenrwLS, leniwLS));
}

int ARKStepGetNumJacEvals(void* arkode_mem, long int* njevals)
{
  return (ARKodeGetNumJacEvals(arkode_mem, njevals));
}

int ARKStepGetNumPrecEvals(void* arkode_mem, long int* npevals)
{
  return (ARKodeGetNumPrecEvals(arkode_mem, npevals));
}

int ARKStepGetNumPrecSolves(void* arkode_mem, long int* npsolves)
{
  return (ARKodeGetNumPrecSolves(arkode_mem, npsolves));
}

int ARKStepGetNumLinIters(void* arkode_mem, long int* nliters)
{
  return (ARKodeGetNumLinIters(arkode_mem, nliters));
}

int ARKStepGetNumLinConvFails(void* arkode_mem, long int* nlcfails)
{
  return (ARKodeGetNumLinConvFails(arkode_mem, nlcfails));
}

int ARKStepGetNumJTSetupEvals(void* arkode_mem, long int* njtsetups)
{
  return (ARKodeGetNumJTSetupEvals(arkode_mem, njtsetups));
}

int ARKStepGetNumJtimesEvals(void* arkode_mem, long int* njvevals)
{
  return (ARKodeGetNumJtimesEvals(arkode_mem, njvevals));
}

int ARKStepGetNumLinRhsEvals(void* arkode_mem, long int* nfevalsLS)
{
  return (ARKodeGetNumLinRhsEvals(arkode_mem, nfevalsLS));
}

int ARKStepGetLastLinFlag(void* arkode_mem, long int* flag)
{
  return (ARKodeGetLastLinFlag(arkode_mem, flag));
}

int ARKStepGetMassWorkSpace(void* arkode_mem, long int* lenrwMLS,
                            long int* leniwMLS)
{
  return (ARKodeGetMassWorkSpace(arkode_mem, lenrwMLS, leniwMLS));
}

int ARKStepGetNumMassSetups(void* arkode_mem, long int* nmsetups)
{
  return (ARKodeGetNumMassSetups(arkode_mem, nmsetups));
}

int ARKStepGetNumMassMultSetups(void* arkode_mem, long int* nmvsetups)
{
  return (ARKodeGetNumMassMultSetups(arkode_mem, nmvsetups));
}

int ARKStepGetNumMassMult(void* arkode_mem, long int* nmvevals)
{
  return (ARKodeGetNumMassMult(arkode_mem, nmvevals));
}

int ARKStepGetNumMassSolves(void* arkode_mem, long int* nmsolves)
{
  return (ARKodeGetNumMassSolves(arkode_mem, nmsolves));
}

int ARKStepGetNumMassPrecEvals(void* arkode_mem, long int* nmpevals)
{
  return (ARKodeGetNumMassPrecEvals(arkode_mem, nmpevals));
}

int ARKStepGetNumMassPrecSolves(void* arkode_mem, long int* nmpsolves)
{
  return (ARKodeGetNumMassPrecSolves(arkode_mem, nmpsolves));
}

int ARKStepGetNumMassIters(void* arkode_mem, long int* nmiters)
{
  return (ARKodeGetNumMassIters(arkode_mem, nmiters));
}

int ARKStepGetNumMassConvFails(void* arkode_mem, long int* nmcfails)
{
  return (ARKodeGetNumMassConvFails(arkode_mem, nmcfails));
}

int ARKStepGetNumMTSetups(void* arkode_mem, long int* nmtsetups)
{
  return (ARKodeGetNumMTSetups(arkode_mem, nmtsetups));
}

int ARKStepGetLastMassFlag(void* arkode_mem, long int* flag)
{
  return (ARKodeGetLastMassFlag(arkode_mem, flag));
}

char* ARKStepGetLinReturnFlagName(long int flag)
{
  return (ARKodeGetLinReturnFlagName(flag));
}

void ARKStepFree(void** arkode_mem) { ARKodeFree(arkode_mem); }

void ARKStepPrintMem(void* arkode_mem, FILE* outfile)
{
  ARKodePrintMem(arkode_mem, outfile);
}

int ARKStepSetRelaxFn(void* arkode_mem, ARKRelaxFn rfn, ARKRelaxJacFn rjac)
{
  return (ARKodeSetRelaxFn(arkode_mem, rfn, rjac));
}

int ARKStepSetRelaxEtaFail(void* arkode_mem, sunrealtype eta_rf)
{
  return (ARKodeSetRelaxEtaFail(arkode_mem, eta_rf));
}

int ARKStepSetRelaxLowerBound(void* arkode_mem, sunrealtype lower)
{
  return (ARKodeSetRelaxLowerBound(arkode_mem, lower));
}

int ARKStepSetRelaxMaxFails(void* arkode_mem, int max_fails)
{
  return (ARKodeSetRelaxMaxFails(arkode_mem, max_fails));
}

int ARKStepSetRelaxMaxIters(void* arkode_mem, int max_iters)
{
  return (ARKodeSetRelaxMaxIters(arkode_mem, max_iters));
}

int ARKStepSetRelaxSolver(void* arkode_mem, ARKRelaxSolver solver)
{
  return (ARKodeSetRelaxSolver(arkode_mem, solver));
}

int ARKStepSetRelaxResTol(void* arkode_mem, sunrealtype res_tol)
{
  return (ARKodeSetRelaxResTol(arkode_mem, res_tol));
}

int ARKStepSetRelaxTol(void* arkode_mem, sunrealtype rel_tol, sunrealtype abs_tol)
{
  return (ARKodeSetRelaxTol(arkode_mem, rel_tol, abs_tol));
}

int ARKStepSetRelaxUpperBound(void* arkode_mem, sunrealtype upper)
{
  return (ARKodeSetRelaxUpperBound(arkode_mem, upper));
}

int ARKStepGetNumRelaxFnEvals(void* arkode_mem, long int* r_evals)
{
  return (ARKodeGetNumRelaxFnEvals(arkode_mem, r_evals));
}

int ARKStepGetNumRelaxJacEvals(void* arkode_mem, long int* J_evals)
{
  return (ARKodeGetNumRelaxJacEvals(arkode_mem, J_evals));
}

int ARKStepGetNumRelaxFails(void* arkode_mem, long int* relax_fails)
{
  return (ARKodeGetNumRelaxFails(arkode_mem, relax_fails));
}

int ARKStepGetNumRelaxBoundFails(void* arkode_mem, long int* fails)
{
  return (ARKodeGetNumRelaxBoundFails(arkode_mem, fails));
}

int ARKStepGetNumRelaxSolveFails(void* arkode_mem, long int* fails)
{
  return (ARKodeGetNumRelaxSolveFails(arkode_mem, fails));
}

int ARKStepGetNumRelaxSolveIters(void* arkode_mem, long int* iters)
{
  return (ARKodeGetNumRelaxSolveIters(arkode_mem, iters));
}

/*===============================================================
  EOF
  ===============================================================*/
