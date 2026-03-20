/*
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban, Cody J. Balos @ LLNL
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
 * This is the implementation file for the CVDIAG linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvode_diag_impl.h"
#include "cvode_impl.h"

/* Other Constants */

#define FRACT SUN_RCONST(0.1)
#define ONE   SUN_RCONST(1.0)

/* CVDIAG linit, lsetup, lsolve, and lfree routines */

static int CVDiagInit(CVodeMem cv_mem);

static int CVDiagSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                       N_Vector fpred, sunbooleantype* jcurPtr, N_Vector vtemp1,
                       N_Vector vtemp2, N_Vector vtemp3);

static int CVDiagSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                       N_Vector ycur, N_Vector fcur);

static int CVDiagFree(CVodeMem cv_mem);

/*
 * -----------------------------------------------------------------
 * CVDiag
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the diagonal linear solver module.  CVDense first
 * calls the existing lfree routine if this is not NULL.  Then it sets
 * the cv_linit, cv_lsetup, cv_lsolve, cv_lfree fields in (*cvode_mem)
 * to be CVDiagInit, CVDiagSetup, CVDiagSolve, and CVDiagFree,
 * respectively.  It allocates memory for a structure of type
 * CVDiagMemRec and sets the cv_lmem field in (*cvode_mem) to the
 * address of this structure.  It sets setupNonNull in (*cvode_mem) to
 * SUNTRUE.  Finally, it allocates memory for M, bit, and bitcomp.
 * The CVDiag return value is SUCCESS = 0, LMEM_FAIL = -1, or
 * LIN_ILL_INPUT=-2.
 * -----------------------------------------------------------------
 */

int CVDiag(void* cvode_mem)
{
  CVodeMem cv_mem;
  CVDiagMem cvdiag_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CVDIAG_MEM_NULL, __LINE__, __func__, __FILE__,
                   MSGDG_CVMEM_NULL);
    return (CVDIAG_MEM_NULL);
  }
  cv_mem = (CVodeMem)cvode_mem;

  /* Check if N_VCompare and N_VInvTest are present */
  if (cv_mem->cv_tempv->ops->nvcompare == NULL ||
      cv_mem->cv_tempv->ops->nvinvtest == NULL)
  {
    cvProcessError(cv_mem, CVDIAG_ILL_INPUT, __LINE__, __func__, __FILE__,
                   MSGDG_BAD_NVECTOR);
    return (CVDIAG_ILL_INPUT);
  }

  if (cv_mem->cv_lfree != NULL) { cv_mem->cv_lfree(cv_mem); }

  /* Set four main function fields in cv_mem */
  cv_mem->cv_linit   = CVDiagInit;
  cv_mem->cv_lreinit = CVDiagInit;
  cv_mem->cv_lsetup  = CVDiagSetup;
  cv_mem->cv_lsolve  = CVDiagSolve;
  cv_mem->cv_lfree   = CVDiagFree;

  /* Get memory for CVDiagMemRec */
  cvdiag_mem = NULL;
  cvdiag_mem = (CVDiagMem)malloc(sizeof(CVDiagMemRec));
  if (cvdiag_mem == NULL)
  {
    cvProcessError(cv_mem, CVDIAG_MEM_FAIL, __LINE__, __func__, __FILE__,
                   MSGDG_MEM_FAIL);
    return (CVDIAG_MEM_FAIL);
  }

  cvdiag_mem->di_last_flag = CVDIAG_SUCCESS;

  /* Allocate memory for M, bit, and bitcomp */

  cvdiag_mem->di_M = N_VClone(cv_mem->cv_tempv);
  if (cvdiag_mem->di_M == NULL)
  {
    cvProcessError(cv_mem, CVDIAG_MEM_FAIL, __LINE__, __func__, __FILE__,
                   MSGDG_MEM_FAIL);
    free(cvdiag_mem);
    cvdiag_mem = NULL;
    return (CVDIAG_MEM_FAIL);
  }

  cvdiag_mem->di_bit = N_VClone(cv_mem->cv_tempv);
  if (cvdiag_mem->di_bit == NULL)
  {
    cvProcessError(cv_mem, CVDIAG_MEM_FAIL, __LINE__, __func__, __FILE__,
                   MSGDG_MEM_FAIL);
    N_VDestroy(cvdiag_mem->di_M);
    free(cvdiag_mem);
    cvdiag_mem = NULL;
    return (CVDIAG_MEM_FAIL);
  }

  cvdiag_mem->di_bitcomp = N_VClone(cv_mem->cv_tempv);
  if (cvdiag_mem->di_bitcomp == NULL)
  {
    cvProcessError(cv_mem, CVDIAG_MEM_FAIL, __LINE__, __func__, __FILE__,
                   MSGDG_MEM_FAIL);
    N_VDestroy(cvdiag_mem->di_M);
    N_VDestroy(cvdiag_mem->di_bit);
    free(cvdiag_mem);
    cvdiag_mem = NULL;
    return (CVDIAG_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  cv_mem->cv_lmem = cvdiag_mem;

  return (CVDIAG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVDiagGetWorkSpace
 * -----------------------------------------------------------------
 */

int CVDiagGetWorkSpace(void* cvode_mem, long int* lenrwLS, long int* leniwLS)
{
  CVodeMem cv_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CVDIAG_MEM_NULL, __LINE__, __func__, __FILE__,
                   MSGDG_CVMEM_NULL);
    return (CVDIAG_MEM_NULL);
  }
  cv_mem = (CVodeMem)cvode_mem;

  *lenrwLS = 3 * cv_mem->cv_lrw1;
  *leniwLS = 3 * cv_mem->cv_liw1;

  return (CVDIAG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVDiagGetNumRhsEvals
 * -----------------------------------------------------------------
 */

int CVDiagGetNumRhsEvals(void* cvode_mem, long int* nfevalsLS)
{
  CVodeMem cv_mem;
  CVDiagMem cvdiag_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CVDIAG_MEM_NULL, __LINE__, __func__, __FILE__,
                   MSGDG_CVMEM_NULL);
    return (CVDIAG_MEM_NULL);
  }
  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_lmem == NULL)
  {
    cvProcessError(cv_mem, CVDIAG_LMEM_NULL, __LINE__, __func__, __FILE__,
                   MSGDG_LMEM_NULL);
    return (CVDIAG_LMEM_NULL);
  }
  cvdiag_mem = (CVDiagMem)cv_mem->cv_lmem;

  *nfevalsLS = cvdiag_mem->di_nfeDI;

  return (CVDIAG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVDiagGetLastFlag
 * -----------------------------------------------------------------
 */

int CVDiagGetLastFlag(void* cvode_mem, long int* flag)
{
  CVodeMem cv_mem;
  CVDiagMem cvdiag_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CVDIAG_MEM_NULL, __LINE__, __func__, __FILE__,
                   MSGDG_CVMEM_NULL);
    return (CVDIAG_MEM_NULL);
  }
  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_lmem == NULL)
  {
    cvProcessError(cv_mem, CVDIAG_LMEM_NULL, __LINE__, __func__, __FILE__,
                   MSGDG_LMEM_NULL);
    return (CVDIAG_LMEM_NULL);
  }
  cvdiag_mem = (CVDiagMem)cv_mem->cv_lmem;

  *flag = cvdiag_mem->di_last_flag;

  return (CVDIAG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVDiagGetReturnFlagName
 * -----------------------------------------------------------------
 */

char* CVDiagGetReturnFlagName(long int flag)
{
  char* name;

  name = (char*)malloc(30 * sizeof(char));

  switch (flag)
  {
  case CVDIAG_SUCCESS: sprintf(name, "CVDIAG_SUCCESS"); break;
  case CVDIAG_MEM_NULL: sprintf(name, "CVDIAG_MEM_NULL"); break;
  case CVDIAG_LMEM_NULL: sprintf(name, "CVDIAG_LMEM_NULL"); break;
  case CVDIAG_ILL_INPUT: sprintf(name, "CVDIAG_ILL_INPUT"); break;
  case CVDIAG_MEM_FAIL: sprintf(name, "CVDIAG_MEM_FAIL"); break;
  case CVDIAG_INV_FAIL: sprintf(name, "CVDIAG_INV_FAIL"); break;
  case CVDIAG_RHSFUNC_UNRECVR: sprintf(name, "CVDIAG_RHSFUNC_UNRECVR"); break;
  case CVDIAG_RHSFUNC_RECVR: sprintf(name, "CVDIAG_RHSFUNC_RECVR"); break;
  default: sprintf(name, "NONE");
  }

  return (name);
}

/*
 * -----------------------------------------------------------------
 * CVDiagInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the diagonal
 * linear solver.
 * -----------------------------------------------------------------
 */

static int CVDiagInit(CVodeMem cv_mem)
{
  CVDiagMem cvdiag_mem;

  cvdiag_mem = (CVDiagMem)cv_mem->cv_lmem;

  cvdiag_mem->di_nfeDI = 0;

  cvdiag_mem->di_last_flag = CVDIAG_SUCCESS;
  return (0);
}

/*
 * -----------------------------------------------------------------
 * CVDiagSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the diagonal linear
 * solver.  It constructs a diagonal approximation to the Newton matrix
 * M = I - gamma*J, updates counters, and inverts M.
 * -----------------------------------------------------------------
 */

static int CVDiagSetup(CVodeMem cv_mem, SUNDIALS_MAYBE_UNUSED int convfail,
                       N_Vector ypred, N_Vector fpred, sunbooleantype* jcurPtr,
                       N_Vector vtemp1, N_Vector vtemp2,
                       SUNDIALS_MAYBE_UNUSED N_Vector vtemp3)
{
  sunrealtype r;
  N_Vector ftemp, y;
  sunbooleantype invOK;
  CVDiagMem cvdiag_mem;
  int retval;

  cvdiag_mem = (CVDiagMem)cv_mem->cv_lmem;

  /* Rename work vectors for use as temporary values of y and f */
  ftemp = vtemp1;
  y     = vtemp2;

  /* Form y with perturbation = FRACT*(func. iter. correction) */
  r = FRACT * cv_mem->cv_rl1;
#ifdef SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS
  if (cv_mem->cv_usefused)
  {
    cvDiagSetup_formY(cv_mem->cv_h, r, fpred, cv_mem->cv_zn[1], ypred, ftemp, y);
  }
  else
#endif
  {
    N_VLinearSum(cv_mem->cv_h, fpred, -ONE, cv_mem->cv_zn[1], ftemp);
    N_VLinearSum(r, ftemp, ONE, ypred, y);
  }

  /* Evaluate f at perturbed y */
  retval = cv_mem->cv_f(cv_mem->cv_tn, y, cvdiag_mem->di_M, cv_mem->cv_user_data);
  cvdiag_mem->di_nfeDI++;
  if (retval < 0)
  {
    cvProcessError(cv_mem, CVDIAG_RHSFUNC_UNRECVR, __LINE__, __func__, __FILE__,
                   MSGDG_RHSFUNC_FAILED);
    cvdiag_mem->di_last_flag = CVDIAG_RHSFUNC_UNRECVR;
    return (-1);
  }
  if (retval > 0)
  {
    cvdiag_mem->di_last_flag = CVDIAG_RHSFUNC_RECVR;
    return (1);
  }

  /* Construct M = I - gamma*J with J = diag(deltaf_i/deltay_i) */
#ifdef SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS
  if (cv_mem->cv_usefused)
  {
    cvDiagSetup_buildM(FRACT, cv_mem->cv_uround, cv_mem->cv_h, ftemp, fpred,
                       cv_mem->cv_ewt, cvdiag_mem->di_bit,
                       cvdiag_mem->di_bitcomp, y, cvdiag_mem->di_M);
  }
  else
#endif
  {
    N_VLinearSum(ONE, cvdiag_mem->di_M, -ONE, fpred, cvdiag_mem->di_M);
    N_VLinearSum(FRACT, ftemp, -(cv_mem->cv_h), cvdiag_mem->di_M,
                 cvdiag_mem->di_M);
    N_VProd(ftemp, cv_mem->cv_ewt, y);
    /* Protect against deltay_i being at roundoff level */
    N_VCompare(cv_mem->cv_uround, y, cvdiag_mem->di_bit);
    N_VAddConst(cvdiag_mem->di_bit, -ONE, cvdiag_mem->di_bitcomp);
    N_VProd(ftemp, cvdiag_mem->di_bit, y);
    N_VLinearSum(FRACT, y, -ONE, cvdiag_mem->di_bitcomp, y);
    N_VDiv(cvdiag_mem->di_M, y, cvdiag_mem->di_M);
    N_VProd(cvdiag_mem->di_M, cvdiag_mem->di_bit, cvdiag_mem->di_M);
    N_VLinearSum(ONE, cvdiag_mem->di_M, -ONE, cvdiag_mem->di_bitcomp,
                 cvdiag_mem->di_M);
  }

  /* Invert M with test for zero components */
  invOK = N_VInvTest(cvdiag_mem->di_M, cvdiag_mem->di_M);
  if (!invOK)
  {
    cvdiag_mem->di_last_flag = CVDIAG_INV_FAIL;
    return (1);
  }

  /* Set jcur = SUNTRUE, save gamma in gammasv, and return */
  *jcurPtr                 = SUNTRUE;
  cvdiag_mem->di_gammasv   = cv_mem->cv_gamma;
  cvdiag_mem->di_last_flag = CVDIAG_SUCCESS;
  return (0);
}

/*
 * -----------------------------------------------------------------
 * CVDiagSolve
 * -----------------------------------------------------------------
 * This routine performs the solve operation for the diagonal linear
 * solver.  If necessary it first updates gamma in M = I - gamma*J.
 * -----------------------------------------------------------------
 */

static int CVDiagSolve(CVodeMem cv_mem, N_Vector b,
                       SUNDIALS_MAYBE_UNUSED N_Vector weight,
                       SUNDIALS_MAYBE_UNUSED N_Vector ycur,
                       SUNDIALS_MAYBE_UNUSED N_Vector fcur)
{
  sunbooleantype invOK;
  sunrealtype r;
  CVDiagMem cvdiag_mem;

  cvdiag_mem = (CVDiagMem)cv_mem->cv_lmem;

  /* If gamma has changed, update factor in M, and save gamma value */

  if (cvdiag_mem->di_gammasv != cv_mem->cv_gamma)
  {
    r = cv_mem->cv_gamma / cvdiag_mem->di_gammasv;
#ifdef SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS
    if (cv_mem->cv_usefused) { cvDiagSolve_updateM(r, cvdiag_mem->di_M); }
    else
#endif
    {
      N_VInv(cvdiag_mem->di_M, cvdiag_mem->di_M);
      N_VAddConst(cvdiag_mem->di_M, -ONE, cvdiag_mem->di_M);
      N_VScale(r, cvdiag_mem->di_M, cvdiag_mem->di_M);
      N_VAddConst(cvdiag_mem->di_M, ONE, cvdiag_mem->di_M);
    }
    invOK = N_VInvTest(cvdiag_mem->di_M, cvdiag_mem->di_M);
    if (!invOK)
    {
      cvdiag_mem->di_last_flag = CVDIAG_INV_FAIL;
      return (1);
    }
    cvdiag_mem->di_gammasv = cv_mem->cv_gamma;
  }

  /* Apply M-inverse to b */
  N_VProd(b, cvdiag_mem->di_M, b);

  cvdiag_mem->di_last_flag = CVDIAG_SUCCESS;
  return (0);
}

/*
 * -----------------------------------------------------------------
 * CVDiagFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the diagonal linear solver.
 * -----------------------------------------------------------------
 */

static int CVDiagFree(CVodeMem cv_mem)
{
  CVDiagMem cvdiag_mem;

  cvdiag_mem = (CVDiagMem)cv_mem->cv_lmem;

  N_VDestroy(cvdiag_mem->di_M);
  N_VDestroy(cvdiag_mem->di_bit);
  N_VDestroy(cvdiag_mem->di_bitcomp);
  free(cvdiag_mem);
  cv_mem->cv_lmem = NULL;

  return (0);
}
