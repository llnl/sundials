/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ UMBC
 * Based on codes <solver>_lapack.c by: Radu Serban @ LLNL
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
 * This is the implementation file for the LAPACK dense
 * implementation of the SUNLINSOL package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_math.h>
#include <sunlinsol/sunlinsol_lapackdense.h>

#include "sundials_lapack_defs.h"
#include "sundials_macros.h"

/* Interfaces to match 'sunrealtype' with the correct LAPACK functions */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define xgetrf_f77 dgetrf_f77
#define xgetrs_f77 dgetrs_f77
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define xgetrf_f77 sgetrf_f77
#define xgetrs_f77 sgetrs_f77
#else
#error Incompatible sunrealtype for LAPACK; disable LAPACK and rebuild
#endif

#define ONE SUN_RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * LapackDense solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define LAPACKDENSE_CONTENT(S) \
  ((SUNLinearSolverContent_LapackDense)(S->content))
#define PIVOTS(S)   (LAPACKDENSE_CONTENT(S)->pivots)
#define LASTFLAG(S) (LAPACKDENSE_CONTENT(S)->last_flag)

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new LAPACK dense linear solver
 */

SUNLinearSolver SUNLinSol_LapackDense(N_Vector y, SUNMatrix A, SUNContext sunctx)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_LapackDense content;
  sunindextype MatrixRows;

  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_DENSE) { return (NULL); }

  if (SUNDenseMatrix_Rows(A) != SUNDenseMatrix_Columns(A)) { return (NULL); }

  if ((N_VGetVectorID(y) != SUNDIALS_NVEC_SERIAL) &&
      (N_VGetVectorID(y) != SUNDIALS_NVEC_OPENMP) &&
      (N_VGetVectorID(y) != SUNDIALS_NVEC_PTHREADS))
  {
    return (NULL);
  }

  MatrixRows = SUNDenseMatrix_Rows(A);
  if (MatrixRows != N_VGetLength(y)) { return (NULL); }

  /* Create linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty(sunctx);
  if (S == NULL) { return (NULL); }

  /* Attach operations */
  S->ops->gettype    = SUNLinSolGetType_LapackDense;
  S->ops->getid      = SUNLinSolGetID_LapackDense;
  S->ops->initialize = SUNLinSolInitialize_LapackDense;
  S->ops->setup      = SUNLinSolSetup_LapackDense;
  S->ops->solve      = SUNLinSolSolve_LapackDense;
  S->ops->lastflag   = SUNLinSolLastFlag_LapackDense;
  S->ops->space      = SUNLinSolSpace_LapackDense;
  S->ops->free       = SUNLinSolFree_LapackDense;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_LapackDense)malloc(sizeof *content);
  if (content == NULL)
  {
    SUNLinSolFree(S);
    return (NULL);
  }

  /* Attach content */
  S->content = content;

  /* Fill content */
  content->N         = MatrixRows;
  content->last_flag = 0;
  content->pivots    = NULL;

  /* Allocate content */
  content->pivots = (sunindextype*)malloc(MatrixRows * sizeof(sunindextype));
  if (content->pivots == NULL)
  {
    SUNLinSolFree(S);
    return (NULL);
  }

  return (S);
}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_LapackDense(
  SUNDIALS_MAYBE_UNUSED SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_DIRECT);
}

SUNLinearSolver_ID SUNLinSolGetID_LapackDense(SUNDIALS_MAYBE_UNUSED SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_LAPACKDENSE);
}

SUNErrCode SUNLinSolInitialize_LapackDense(SUNLinearSolver S)
{
  /* all solver-specific memory has already been allocated */
  LASTFLAG(S) = SUN_SUCCESS;
  return SUN_SUCCESS;
}

int SUNLinSolSetup_LapackDense(SUNLinearSolver S, SUNMatrix A)
{
  sunindextype n, ier;

  /* check for valid inputs */
  if ((A == NULL) || (S == NULL)) { return SUN_ERR_ARG_CORRUPT; }

  /* Ensure that A is a dense matrix */
  if (SUNMatGetID(A) != SUNMATRIX_DENSE)
  {
    LASTFLAG(S) = SUN_ERR_ARG_INCOMPATIBLE;
    return SUN_ERR_ARG_INCOMPATIBLE;
  }

  /* Call LAPACK to do LU factorization of A */
  n   = SUNDenseMatrix_Rows(A);
  ier = 0;
  xgetrf_f77(&n, &n, SUNDenseMatrix_Data(A), &n, PIVOTS(S), &ier);
  LASTFLAG(S) = ier;
  if (ier > 0) { return (SUNLS_LUFACT_FAIL); }
  if (ier < 0) { return SUN_ERR_EXT_FAIL; }
  return SUN_SUCCESS;
}

int SUNLinSolSolve_LapackDense(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                               N_Vector b, SUNDIALS_MAYBE_UNUSED sunrealtype tol)
{
  sunindextype n, one, ier;
  sunrealtype* xdata;

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
    return SUN_ERR_MEM_FAIL;
  }

  /* Call LAPACK to solve the linear system */
  n   = SUNDenseMatrix_Rows(A);
  one = 1;
  ier = 0;
  xgetrs_f77("N", &n, &one, SUNDenseMatrix_Data(A), &n, PIVOTS(S), xdata, &n,
             &ier);
  LASTFLAG(S) = ier;
  if (ier < 0) { return SUN_ERR_EXT_FAIL; }

  LASTFLAG(S) = SUN_SUCCESS;
  return SUN_SUCCESS;
}

sunindextype SUNLinSolLastFlag_LapackDense(SUNLinearSolver S)
{
  return (LASTFLAG(S));
}

SUNErrCode SUNLinSolSpace_LapackDense(SUNLinearSolver S, long int* lenrwLS,
                                      long int* leniwLS)
{
  *lenrwLS = 0;
  *leniwLS = 2 + LAPACKDENSE_CONTENT(S)->N;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolFree_LapackDense(SUNLinearSolver S)
{
  /* return if S is already free */
  if (S == NULL) { return SUN_SUCCESS; }

  /* delete items from contents, then delete generic structure */
  if (S->content)
  {
    if (PIVOTS(S))
    {
      free(PIVOTS(S));
      PIVOTS(S) = NULL;
    }
    free(S->content);
    S->content = NULL;
  }
  if (S->ops)
  {
    free(S->ops);
    S->ops = NULL;
  }
  free(S);
  S = NULL;
  return SUN_SUCCESS;
}
