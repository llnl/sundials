/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ UMBC
 * Based on codes <solver>_klu.c, written by Carol Woodward @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * This is the implementation file for the KLU implementation of
 * the SUNLINSOL package.
 * -----------------------------------------------------------------*/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_math.h>
#include <sunlinsol/sunlinsol_klu.h>

#include "sundials_cli.h"
#include "sundials_macros.h"

#define ONE       SUN_RCONST(1.0)
#define TWOTHIRDS SUN_RCONST(0.666666666666666666666666666666667)

/*
 * -----------------------------------------------------------------
 * KLU solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define KLU_CONTENT(S)    ((SUNLinearSolverContent_KLU)(S->content))
#define LASTFLAG(S)       (KLU_CONTENT(S)->last_flag)
#define FIRSTFACTORIZE(S) (KLU_CONTENT(S)->first_factorize)
#define SYMBOLIC(S)       (KLU_CONTENT(S)->symbolic)
#define NUMERIC(S)        (KLU_CONTENT(S)->numeric)
#define COMMON(S)         (KLU_CONTENT(S)->common)
#define SOLVE(S)          (KLU_CONTENT(S)->klu_solver)

/*
 * ----------------------------------------------------------------------------
 * Un-exported implementation specific routines
 * ----------------------------------------------------------------------------
 */

static SUNErrCode setFromCommandLine_KLU(SUNLinearSolver S, const char* LSid,
                                         int argc, char* argv[]);

SUNErrCode SUNLinSolSetOptions_KLU(SUNLinearSolver S, const char* LSid,
                                   const char* file_name, int argc, char* argv[]);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new KLU linear solver
 */

SUNLinearSolver SUNLinSol_KLU(N_Vector y, SUNMatrix A, SUNContext sunctx)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_KLU content;
  int flag;

  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE) { return (NULL); }

  if (SUNSparseMatrix_Rows(A) != SUNSparseMatrix_Columns(A)) { return (NULL); }

  if ((N_VGetVectorID(y) != SUNDIALS_NVEC_SERIAL) &&
      (N_VGetVectorID(y) != SUNDIALS_NVEC_OPENMP) &&
      (N_VGetVectorID(y) != SUNDIALS_NVEC_PTHREADS))
  {
    return (NULL);
  }

  if (SUNSparseMatrix_Rows(A) != N_VGetLength(y)) { return (NULL); }

  /* Create an empty linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty(sunctx);
  if (S == NULL) { return (NULL); }

  /* Attach operations */
  S->ops->gettype    = SUNLinSolGetType_KLU;
  S->ops->getid      = SUNLinSolGetID_KLU;
  S->ops->setoptions = SUNLinSolSetOptions_KLU;
  S->ops->initialize = SUNLinSolInitialize_KLU;
  S->ops->setup      = SUNLinSolSetup_KLU;
  S->ops->solve      = SUNLinSolSolve_KLU;
  S->ops->lastflag   = SUNLinSolLastFlag_KLU;
  S->ops->space      = SUNLinSolSpace_KLU;
  S->ops->free       = SUNLinSolFree_KLU;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_KLU)malloc(sizeof *content);
  if (content == NULL)
  {
    SUNLinSolFree(S);
    return (NULL);
  }

  /* Attach content */
  S->content = content;

  /* Fill content */
  content->last_flag       = 0;
  content->first_factorize = 1;
  content->symbolic        = NULL;
  content->numeric         = NULL;

#if defined(SUNDIALS_INT64_T)
  if (SUNSparseMatrix_SparseType(A) == SUN_CSC_MAT)
  {
    content->klu_solver = (KLUSolveFn)&klu_l_solve;
  }
  else { content->klu_solver = (KLUSolveFn)&klu_l_tsolve; }
#elif defined(SUNDIALS_INT32_T)
  if (SUNSparseMatrix_SparseType(A) == SUN_CSC_MAT)
  {
    content->klu_solver = &klu_solve;
  }
  else { content->klu_solver = &klu_tsolve; }
#else /* incompatible sunindextype for KLU */
#error Incompatible sunindextype for KLU
#endif

  flag = sun_klu_defaults(&(content->common));
  if (flag == 0)
  {
    SUNLinSolFree(S);
    return (NULL);
  }
  (content->common).ordering = SUNKLU_ORDERING_DEFAULT;

  return (S);
}

/* ----------------------------------------------------------------------------
 * Function to reinitialize a KLU linear solver
 */

SUNErrCode SUNLinSol_KLUReInit(SUNLinearSolver S, SUNMatrix A, sunindextype nnz,
                               int reinit_type)
{
  /* Check for non-NULL SUNLinearSolver */
  if ((S == NULL) || (A == NULL)) { return SUN_ERR_ARG_CORRUPT; }

  /* Check for valid SUNMatrix */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE) { return SUN_ERR_ARG_INCOMPATIBLE; }

  /* Check for valid reinit_type */
  if ((reinit_type != SUNKLU_REINIT_FULL) &&
      (reinit_type != SUNKLU_REINIT_PARTIAL))
  {
    return SUN_ERR_ARG_INCOMPATIBLE;
  }

  /* Full re-initialization: reallocate matrix for updated storage */
  if (reinit_type == SUNKLU_REINIT_FULL)
  {
    if (SUNSparseMatrix_Reallocate(A, nnz) != 0) { return SUN_ERR_MEM_FAIL; }
  }

  /* Free the prior factorization and reset for first factorization */
  if (SYMBOLIC(S) != NULL) { sun_klu_free_symbolic(&SYMBOLIC(S), &COMMON(S)); }
  if (NUMERIC(S) != NULL) { sun_klu_free_numeric(&NUMERIC(S), &COMMON(S)); }
  FIRSTFACTORIZE(S) = 1;

  LASTFLAG(S) = SUN_SUCCESS;
  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Function to control set routines via the command line or file
 */

SUNErrCode SUNLinSolSetOptions_KLU(SUNLinearSolver S, const char* LSid,
                                   const char* file_name, int argc, char* argv[])
{
  SUNFunctionBegin(S->sunctx);

  /* File-based option control is currently unimplemented */
  SUNAssert((file_name == NULL || strlen(file_name) == 0),
            SUN_ERR_ARG_INCOMPATIBLE);

  if (argc > 0 && argv != NULL)
  {
    SUNCheckCall(setFromCommandLine_KLU(S, LSid, argc, argv));
  }

  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Function to control set routines via the command line
 */

static SUNErrCode setFromCommandLine_KLU(SUNLinearSolver S, const char* LSid,
                                         int argc, char* argv[])
{
  SUNFunctionBegin(S->sunctx);

  /* Prefix for options to set */
  const char* default_id = "sunlinearsolver";
  size_t offset          = strlen(default_id) + 1;
  if (LSid != NULL && strlen(LSid) > 0) { offset = strlen(LSid) + 1; }
  char* prefix = (char*)malloc(sizeof(char) * (offset + 1));
  if (LSid != NULL && strlen(LSid) > 0) { strcpy(prefix, LSid); }
  else { strcpy(prefix, default_id); }
  strcat(prefix, ".");

  for (int idx = 1; idx < argc; idx++)
  {
    int retval;

    /* skip command-line arguments that do not begin with correct prefix */
    if (strncmp(argv[idx], prefix, strlen(prefix)) != 0) { continue; }

    /* control over SetOrdering function */
    if (strcmp(argv[idx] + offset, "ordering") == 0)
    {
      idx += 1;
      int iarg = atoi(argv[idx]);
      retval   = SUNLinSol_KLUSetOrdering(S, iarg);
      if (retval != SUN_SUCCESS)
      {
        free(prefix);
        return retval;
      }
      continue;
    }
  }
  free(prefix);
  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Function to set the ordering type for a KLU linear solver
 */

SUNErrCode SUNLinSol_KLUSetOrdering(SUNLinearSolver S, int ordering_choice)
{
  /* Check for legal ordering_choice */
  if ((ordering_choice < 0) || (ordering_choice > 2))
  {
    return SUN_ERR_ARG_INCOMPATIBLE;
  }

  /* Check for non-NULL SUNLinearSolver */
  if (S == NULL) { return SUN_ERR_ARG_CORRUPT; }

  /* Set ordering_choice */
  COMMON(S).ordering = ordering_choice;

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * accessor functions
 * -----------------------------------------------------------------
 */

sun_klu_symbolic* SUNLinSol_KLUGetSymbolic(SUNLinearSolver S)
{
  return (SYMBOLIC(S));
}

sun_klu_numeric* SUNLinSol_KLUGetNumeric(SUNLinearSolver S)
{
  return (NUMERIC(S));
}

sun_klu_common* SUNLinSol_KLUGetCommon(SUNLinearSolver S)
{
  return (&(COMMON(S)));
}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_KLU(SUNDIALS_MAYBE_UNUSED SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_DIRECT);
}

SUNLinearSolver_ID SUNLinSolGetID_KLU(SUNDIALS_MAYBE_UNUSED SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_KLU);
}

SUNErrCode SUNLinSolInitialize_KLU(SUNLinearSolver S)
{
  /* Force factorization */
  FIRSTFACTORIZE(S) = 1;

  LASTFLAG(S) = SUN_SUCCESS;
  return (LASTFLAG(S));
}

int SUNLinSolSetup_KLU(SUNLinearSolver S, SUNMatrix A)
{
  int retval;
  sunrealtype uround_twothirds;

  uround_twothirds = SUNRpowerR(SUN_UNIT_ROUNDOFF, TWOTHIRDS);

  /* Ensure that A is a sparse matrix */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE)
  {
    LASTFLAG(S) = SUN_ERR_ARG_INCOMPATIBLE;
    return (LASTFLAG(S));
  }

  /* On first decomposition, get the symbolic factorization */
  if (FIRSTFACTORIZE(S))
  {
    /* Perform symbolic analysis of sparsity structure */
    if (SYMBOLIC(S)) { sun_klu_free_symbolic(&SYMBOLIC(S), &COMMON(S)); }
    SYMBOLIC(S) = sun_klu_analyze(SUNSparseMatrix_NP(A),
                                  SUNSparseMatrix_IndexPointers(A),
                                  SUNSparseMatrix_IndexValues(A), &COMMON(S));
    if (SYMBOLIC(S) == NULL)
    {
      LASTFLAG(S) = SUN_ERR_EXT_FAIL;
      return (LASTFLAG(S));
    }

    /* ------------------------------------------------------------
       Compute the LU factorization of the matrix
       ------------------------------------------------------------*/
    if (NUMERIC(S)) { sun_klu_free_numeric(&NUMERIC(S), &COMMON(S)); }
    NUMERIC(S) = sun_klu_factor(SUNSparseMatrix_IndexPointers(A),
                                SUNSparseMatrix_IndexValues(A),
                                SUNSparseMatrix_Data(A), SYMBOLIC(S), &COMMON(S));
    if (NUMERIC(S) == NULL)
    {
      LASTFLAG(S) = SUN_ERR_EXT_FAIL;
      return (LASTFLAG(S));
    }

    FIRSTFACTORIZE(S) = 0;
  }
  else
  { /* not the first decomposition, so just refactor */

    retval = sun_klu_refactor(SUNSparseMatrix_IndexPointers(A),
                              SUNSparseMatrix_IndexValues(A),
                              SUNSparseMatrix_Data(A), SYMBOLIC(S), NUMERIC(S),
                              &COMMON(S));
    if (retval == 0)
    {
      LASTFLAG(S) = SUNLS_PACKAGE_FAIL_REC;
      return (LASTFLAG(S));
    }

    /*-----------------------------------------------------------
      Check if a cheap estimate of the reciprocal of the condition
      number is getting too small.  If so, delete
      the prior numeric factorization and recompute it.
      -----------------------------------------------------------*/

    retval = sun_klu_rcond(SYMBOLIC(S), NUMERIC(S), &COMMON(S));
    if (retval == 0)
    {
      LASTFLAG(S) = SUNLS_PACKAGE_FAIL_REC;
      return (LASTFLAG(S));
    }

    if (COMMON(S).rcond < uround_twothirds)
    {
      /* Condition number may be getting large.
	 Compute more accurate estimate */
      retval = sun_klu_condest(SUNSparseMatrix_IndexPointers(A),
                               SUNSparseMatrix_Data(A), SYMBOLIC(S), NUMERIC(S),
                               &COMMON(S));
      if (retval == 0)
      {
        LASTFLAG(S) = SUNLS_PACKAGE_FAIL_REC;
        return (LASTFLAG(S));
      }

      if (COMMON(S).condest > (ONE / uround_twothirds))
      {
        /* More accurate estimate also says condition number is
           large, so recompute the numeric factorization */
        sun_klu_free_numeric(&NUMERIC(S), &COMMON(S));
        NUMERIC(S) = sun_klu_factor(SUNSparseMatrix_IndexPointers(A),
                                    SUNSparseMatrix_IndexValues(A),
                                    SUNSparseMatrix_Data(A), SYMBOLIC(S),
                                    &COMMON(S));
        if (NUMERIC(S) == NULL)
        {
          LASTFLAG(S) = SUN_ERR_EXT_FAIL;
          return (LASTFLAG(S));
        }
      }
    }
  }

  LASTFLAG(S) = SUN_SUCCESS;
  return (LASTFLAG(S));
}

int SUNLinSolSolve_KLU(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b,
                       SUNDIALS_MAYBE_UNUSED sunrealtype tol)
{
  int flag;
  sunrealtype* xdata;

  /* check for valid inputs */
  if ((A == NULL) || (S == NULL) || (x == NULL) || (b == NULL))
  {
    return SUN_ERR_ARG_CORRUPT;
  }

  /* copy b into x */
  N_VScale(ONE, b, x);

  /* access x data array */
  xdata = N_VGetArrayPointer(x);
  if (xdata == NULL)
  {
    LASTFLAG(S) = SUN_ERR_MEM_FAIL;
    return (LASTFLAG(S));
  }

  /* Call KLU to solve the linear system */
  flag = SOLVE(S)(SYMBOLIC(S), NUMERIC(S), SUNSparseMatrix_NP(A), 1, xdata,
                  &COMMON(S));
  if (flag == 0)
  {
    LASTFLAG(S) = SUNLS_PACKAGE_FAIL_REC;
    return (LASTFLAG(S));
  }

  LASTFLAG(S) = SUN_SUCCESS;
  return (LASTFLAG(S));
}

sunindextype SUNLinSolLastFlag_KLU(SUNLinearSolver S) { return (LASTFLAG(S)); }

SUNErrCode SUNLinSolSpace_KLU(SUNDIALS_MAYBE_UNUSED SUNLinearSolver S,
                              long int* lenrwLS, long int* leniwLS)
{
  /* since the klu structures are opaque objects, we
     omit those from these results */
  *leniwLS = 2;
  *lenrwLS = 0;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolFree_KLU(SUNLinearSolver S)
{
  /* return with success if already freed */
  if (S == NULL) { return SUN_SUCCESS; }

  /* delete items from the contents structure (if it exists) */
  if (S->content)
  {
    if (NUMERIC(S)) { sun_klu_free_numeric(&NUMERIC(S), &COMMON(S)); }
    if (SYMBOLIC(S)) { sun_klu_free_symbolic(&SYMBOLIC(S), &COMMON(S)); }
    free(S->content);
    S->content = NULL;
  }

  /* delete generic structures */
  if (S->ops)
  {
    free(S->ops);
    S->ops = NULL;
  }
  free(S);
  S = NULL;
  return SUN_SUCCESS;
}
