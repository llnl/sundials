/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban and Cosmin Petra @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * This is the implementation file for the optional input and output
 * functions for the adjoint module in the IDAS solver.
 * -----------------------------------------------------------------
 */

/*
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>

#include "idas_impl.h"

/*
 * =================================================================
 * IDAA PRIVATE CONSTANTS
 * =================================================================
 */

#define ONE SUN_RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * Optional input functions for ASA
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * IDAAdjSetNoSensi
 * -----------------------------------------------------------------
 * Disables the forward sensitivity analysis in IDASolveF.
 * -----------------------------------------------------------------
 */

int IDAAdjSetNoSensi(void* ida_mem)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;

  /* Is ida_mem valid? */
  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == SUNFALSE)
  {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, __LINE__, __func__, __FILE__,
                    MSGAM_NO_ADJ);
    return (IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  IDAADJ_mem->ia_storeSensi = SUNFALSE;

  return (IDA_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Optional input functions for backward integration
 * -----------------------------------------------------------------
 */

int IDASetNonlinearSolverB(void* ida_mem, int which, SUNNonlinearSolver NLS)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void* ida_memB;

  /* Check if ida_mem exists */
  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSGAM_NULL_IDAMEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Was ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == SUNFALSE)
  {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, __LINE__, __func__, __FILE__,
                    MSGAM_NO_ADJ);
    return (IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if (which >= IDAADJ_mem->ia_nbckpbs)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSGAM_BAD_WHICH);
    return (IDA_ILL_INPUT);
  }

  /* Find the IDABMem entry in the linked list corresponding to 'which' */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL)
  {
    if (which == IDAB_mem->ida_index) { break; }
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }

  ida_memB = (void*)(IDAB_mem->IDA_mem);

  return (IDASetNonlinearSolver(ida_memB, NLS));
}

int IDASetUserDataB(void* ida_mem, int which, void* user_dataB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;

  /* Is ida_mem valid? */
  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == SUNFALSE)
  {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, __LINE__, __func__, __FILE__,
                    MSGAM_NO_ADJ);
    return (IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if (which >= IDAADJ_mem->ia_nbckpbs)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSGAM_BAD_WHICH);
    return (IDA_ILL_INPUT);
  }

  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL)
  {
    if (which == IDAB_mem->ida_index) { break; }
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }

  /* Set user data for this backward problem. */
  IDAB_mem->ida_user_data = user_dataB;

  return (IDA_SUCCESS);
}

int IDASetMaxOrdB(void* ida_mem, int which, int maxordB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void* ida_memB;

  /* Is ida_mem valid? */
  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == SUNFALSE)
  {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, __LINE__, __func__, __FILE__,
                    MSGAM_NO_ADJ);
    return (IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if (which >= IDAADJ_mem->ia_nbckpbs)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSGAM_BAD_WHICH);
    return (IDA_ILL_INPUT);
  }

  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL)
  {
    if (which == IDAB_mem->ida_index) { break; }
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }
  ida_memB = (void*)IDAB_mem->IDA_mem;

  return IDASetMaxOrd(ida_memB, maxordB);
}

int IDASetMaxNumStepsB(void* ida_mem, int which, long int mxstepsB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void* ida_memB;

  /* Is ida_mem valid? */
  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == SUNFALSE)
  {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, __LINE__, __func__, __FILE__,
                    MSGAM_NO_ADJ);
    return (IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if (which >= IDAADJ_mem->ia_nbckpbs)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSGAM_BAD_WHICH);
    return (IDA_ILL_INPUT);
  }

  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL)
  {
    if (which == IDAB_mem->ida_index) { break; }
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }
  ida_memB = (void*)IDAB_mem->IDA_mem;

  return IDASetMaxNumSteps(ida_memB, mxstepsB);
}

int IDASetInitStepB(void* ida_mem, int which, sunrealtype hinB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void* ida_memB;

  /* Is ida_mem valid? */
  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == SUNFALSE)
  {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, __LINE__, __func__, __FILE__,
                    MSGAM_NO_ADJ);
    return (IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if (which >= IDAADJ_mem->ia_nbckpbs)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSGAM_BAD_WHICH);
    return (IDA_ILL_INPUT);
  }

  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL)
  {
    if (which == IDAB_mem->ida_index) { break; }
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }
  ida_memB = (void*)IDAB_mem->IDA_mem;

  return IDASetInitStep(ida_memB, hinB);
}

int IDASetMaxStepB(void* ida_mem, int which, sunrealtype hmaxB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void* ida_memB;

  /* Is ida_mem valid? */
  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == SUNFALSE)
  {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, __LINE__, __func__, __FILE__,
                    MSGAM_NO_ADJ);
    return (IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if (which >= IDAADJ_mem->ia_nbckpbs)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSGAM_BAD_WHICH);
    return (IDA_ILL_INPUT);
  }

  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL)
  {
    if (which == IDAB_mem->ida_index) { break; }
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }
  ida_memB = (void*)IDAB_mem->IDA_mem;

  return IDASetMaxStep(ida_memB, hmaxB);
}

int IDASetSuppressAlgB(void* ida_mem, int which, sunbooleantype suppressalgB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void* ida_memB;

  /* Is ida_mem valid? */
  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == SUNFALSE)
  {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, __LINE__, __func__, __FILE__,
                    MSGAM_NO_ADJ);
    return (IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if (which >= IDAADJ_mem->ia_nbckpbs)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSGAM_BAD_WHICH);
    return (IDA_ILL_INPUT);
  }

  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL)
  {
    if (which == IDAB_mem->ida_index) { break; }
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }
  ida_memB = (void*)IDAB_mem->IDA_mem;

  return IDASetSuppressAlg(ida_memB, suppressalgB);
}

int IDASetIdB(void* ida_mem, int which, N_Vector idB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void* ida_memB;

  /* Is ida_mem valid? */
  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == SUNFALSE)
  {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, __LINE__, __func__, __FILE__,
                    MSGAM_NO_ADJ);
    return (IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if (which >= IDAADJ_mem->ia_nbckpbs)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSGAM_BAD_WHICH);
    return (IDA_ILL_INPUT);
  }

  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL)
  {
    if (which == IDAB_mem->ida_index) { break; }
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }
  ida_memB = (void*)IDAB_mem->IDA_mem;

  return IDASetId(ida_memB, idB);
}

int IDASetConstraintsB(void* ida_mem, int which, N_Vector constraintsB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void* ida_memB;

  /* Is ida_mem valid? */
  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == SUNFALSE)
  {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, __LINE__, __func__, __FILE__,
                    MSGAM_NO_ADJ);
    return (IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if (which >= IDAADJ_mem->ia_nbckpbs)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSGAM_BAD_WHICH);
    return (IDA_ILL_INPUT);
  }

  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL)
  {
    if (which == IDAB_mem->ida_index) { break; }
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }
  ida_memB = (void*)IDAB_mem->IDA_mem;

  return IDASetConstraints(ida_memB, constraintsB);
}

/*
 * ----------------------------------------------------------------
 * Input quadrature functions for ASA
 * ----------------------------------------------------------------
 */

int IDASetQuadErrConB(void* ida_mem, int which, int errconQB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void* ida_memB;

  /* Is ida_mem valid? */
  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == SUNFALSE)
  {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, __LINE__, __func__, __FILE__,
                    MSGAM_NO_ADJ);
    return (IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if (which >= IDAADJ_mem->ia_nbckpbs)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSGAM_BAD_WHICH);
    return (IDA_ILL_INPUT);
  }

  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL)
  {
    if (which == IDAB_mem->ida_index) { break; }
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }
  ida_memB = (void*)IDAB_mem->IDA_mem;

  return IDASetQuadErrCon(ida_memB, errconQB);
}

/*
 * -----------------------------------------------------------------
 * Optional output functions for backward integration
 * -----------------------------------------------------------------
 */

/*
 * IDAGetAdjIDABmem
 *
 * This function returns a (void *) pointer to the IDAS
 * memory allocated for the backward problem. This pointer can
 * then be used to call any of the IDAGet* IDAS routines to
 * extract optional output for the backward integration phase.
 */

void* IDAGetAdjIDABmem(void* ida_mem, int which)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void* ida_memB;

  /* Is ida_mem valid? */
  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, 0, __LINE__, __func__, __FILE__, MSGAM_NULL_IDAMEM);
    return (NULL);
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == SUNFALSE)
  {
    IDAProcessError(IDA_mem, 0, __LINE__, __func__, __FILE__, MSGAM_NO_ADJ);
    return (NULL);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if (which >= IDAADJ_mem->ia_nbckpbs)
  {
    IDAProcessError(IDA_mem, 0, __LINE__, __func__, __FILE__, MSGAM_BAD_WHICH);
    return (NULL);
  }

  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL)
  {
    if (which == IDAB_mem->ida_index) { break; }
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }
  ida_memB = (void*)IDAB_mem->IDA_mem;

  return (ida_memB);
}

/*
 * IDAGetAdjCheckPointsInfo
 *
 * Loads an array of nckpnts structures of type IDAadjCheckPointRec
 * defined below.
 *
 * The user must allocate space for ckpnt (ncheck+1).
 */

int IDAGetAdjCheckPointsInfo(void* ida_mem, IDAadjCheckPointRec* ckpnt)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDAckpntMem ck_mem;
  int i;

  /* Is ida_mem valid? */
  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSGAM_NULL_IDAMEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == SUNFALSE)
  {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, __LINE__, __func__, __FILE__,
                    MSGAM_NO_ADJ);
    return (IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  i      = 0;
  ck_mem = IDAADJ_mem->ck_mem;
  while (ck_mem != NULL)
  {
    ckpnt[i].my_addr   = (void*)ck_mem;
    ckpnt[i].next_addr = (void*)ck_mem->ck_next;
    ckpnt[i].t0        = ck_mem->ck_t0;
    ckpnt[i].t1        = ck_mem->ck_t1;
    ckpnt[i].nstep     = ck_mem->ck_nst;
    ckpnt[i].order     = ck_mem->ck_kk;
    ckpnt[i].step      = ck_mem->ck_hh;

    ck_mem = ck_mem->ck_next;
    i++;
  }

  return (IDA_SUCCESS);
}

/* IDAGetConsistentICB
 *
 * Returns the consistent initial conditions computed by IDACalcICB or
 * IDACalcICBS
 *
 * It must be preceded by a successful call to IDACalcICB or IDACalcICBS
 * for 'which' backward problem.
 */

int IDAGetConsistentICB(void* ida_mem, int which, N_Vector yyB0_mod,
                        N_Vector ypB0_mod)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDABMem IDAB_mem;
  void* ida_memB;
  int flag;

  /* Is ida_mem valid? */
  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSGAM_NULL_IDAMEM);
    return IDA_MEM_NULL;
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == SUNFALSE)
  {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, __LINE__, __func__, __FILE__,
                    MSGAM_NO_ADJ);
    return (IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  /* Check the value of which */
  if (which >= IDAADJ_mem->ia_nbckpbs)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSGAM_BAD_WHICH);
    return (IDA_ILL_INPUT);
  }

  /* Find the IDABMem entry in the linked list corresponding to 'which'. */
  IDAB_mem = IDAADJ_mem->IDAB_mem;
  while (IDAB_mem != NULL)
  {
    if (which == IDAB_mem->ida_index) { break; }
    /* advance */
    IDAB_mem = IDAB_mem->ida_next;
  }
  ida_memB = (void*)IDAB_mem->IDA_mem;

  flag = IDAGetConsistentIC(ida_memB, yyB0_mod, ypB0_mod);

  return (flag);
}

/*-----------------------------------------------------------------*/

int IDAGetUserDataB(void* ida_mem, int which, void** user_dataB)
{
  IDAMem IDA_mem;
  IDAadjMem IDAa_mem;
  IDABMem IDAB_mem;

  /* Check if IDA_mem exists */
  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSGAM_NULL_IDAMEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Was ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == SUNFALSE)
  {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, __LINE__, __func__, __FILE__,
                    MSGAM_NO_ADJ);
    return (IDA_NO_ADJ);
  }
  IDAa_mem = IDA_mem->ida_adj_mem;

  /* Check which */
  if (which >= IDAa_mem->ia_nbckpbs)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSGAM_BAD_WHICH);
    return (IDA_ILL_INPUT);
  }

  /* Find the IDABMem entry in the linked list corresponding to which */
  IDAB_mem = IDAa_mem->IDAB_mem;
  while (IDAB_mem != NULL)
  {
    if (which == IDAB_mem->ida_index) { break; }
    IDAB_mem = IDAB_mem->ida_next;
  }

  *user_dataB = IDAB_mem->ida_user_data;

  return (IDA_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Undocumented development user-callable functions
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * IDAGetAdjDataPointHermite
 * -----------------------------------------------------------------
 * Returns the 2 vectors stored for cubic Hermite interpolation at
 * the data point 'which'. The user must allocate space for yy and
 * yd.
 *
 * Returns IDA_MEM_NULL if ida_mem is NULL, IDA_ILL_INPUT if the
 * interpolation type previously specified is not IDA_HERMITE or
 * IDA_SUCCESS otherwise.
 *
 */
int IDAGetAdjDataPointHermite(void* ida_mem, int which, sunrealtype* t,
                              N_Vector yy, N_Vector yd)

{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDAdtpntMem* dt_mem;
  IDAhermiteDataMem content;

  /* Is ida_mem valid? */
  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSGAM_NULL_IDAMEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == SUNFALSE)
  {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, __LINE__, __func__, __FILE__,
                    MSGAM_NO_ADJ);
    return (IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  dt_mem = IDAADJ_mem->dt_mem;

  if (IDAADJ_mem->ia_interpType != IDA_HERMITE)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSGAM_WRONG_INTERP);
    return (IDA_ILL_INPUT);
  }

  *t      = dt_mem[which]->t;
  content = (IDAhermiteDataMem)dt_mem[which]->content;

  if (yy != NULL) { N_VScale(ONE, content->y, yy); }
  if (yd != NULL) { N_VScale(ONE, content->yd, yd); }

  return (IDA_SUCCESS);
}

/*
 * IDAGetAdjDataPointPolynomial
 *
 * Returns the vector stored for polynomial interpolation at the
 * data point 'which'. The user must allocate space for y.
 *
 * Returns IDA_MEM_NULL if ida_mem is NULL, IDA_ILL_INPUT if the
 * interpolation type previously specified is not IDA_POLYNOMIAL or
 * IDA_SUCCESS otherwise.
 */

int IDAGetAdjDataPointPolynomial(void* ida_mem, int which, sunrealtype* t,
                                 int* order, N_Vector y)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;
  IDAdtpntMem* dt_mem;
  IDApolynomialDataMem content;
  /* Is ida_mem valid? */
  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSGAM_NULL_IDAMEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == SUNFALSE)
  {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, __LINE__, __func__, __FILE__,
                    MSGAM_NO_ADJ);
    return (IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  dt_mem = IDAADJ_mem->dt_mem;

  if (IDAADJ_mem->ia_interpType != IDA_POLYNOMIAL)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSGAM_WRONG_INTERP);
    return (IDA_ILL_INPUT);
  }

  *t      = dt_mem[which]->t;
  content = (IDApolynomialDataMem)dt_mem[which]->content;

  if (y != NULL) { N_VScale(ONE, content->y, y); }

  *order = content->order;

  return (IDA_SUCCESS);
}

/*
 * IDAGetAdjCurrentCheckPoint
 *
 * Returns the address of the 'active' check point.
 */

int IDAGetAdjCurrentCheckPoint(void* ida_mem, void** addr)
{
  IDAMem IDA_mem;
  IDAadjMem IDAADJ_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSGAM_NULL_IDAMEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem)ida_mem;

  /* Is ASA initialized? */
  if (IDA_mem->ida_adjMallocDone == SUNFALSE)
  {
    IDAProcessError(IDA_mem, IDA_NO_ADJ, __LINE__, __func__, __FILE__,
                    MSGAM_NO_ADJ);
    return (IDA_NO_ADJ);
  }
  IDAADJ_mem = IDA_mem->ida_adj_mem;

  *addr = (void*)IDAADJ_mem->ia_ckpntData;

  return (IDA_SUCCESS);
}
