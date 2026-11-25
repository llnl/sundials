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

#include <limits.h>
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

#if SUNDIALS_LAPACK_CONVERT_PIVOTS
typedef struct
{
  struct _SUNLinearSolverContent_LapackDense content;
  sunlapackinttype* pivots;
} SunLapackDenseStorage;
#define LAPACK_PIVOTS(S) (((SunLapackDenseStorage*)(S->content))->pivots)
#else
#define LAPACK_PIVOTS(S) PIVOTS(S)
#endif

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

  if (!sunLapackIndexInRange(MatrixRows) ||
      (uintmax_t)MatrixRows > SIZE_MAX / sizeof(sunindextype))
  {
    return NULL;
  }

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
#if SUNDIALS_LAPACK_CONVERT_PIVOTS
  SunLapackDenseStorage* storage = malloc(sizeof *storage);
  if (storage != NULL)
  {
    storage->pivots = NULL;
    content = &storage->content;
  }
#else
  content = (SUNLinearSolverContent_LapackDense)malloc(sizeof *content);
#endif
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
  content->pivots = calloc((size_t)MatrixRows, sizeof(*content->pivots));
  if (content->pivots == NULL)
  {
    SUNLinSolFree(S);
    return (NULL);
  }

#if SUNDIALS_LAPACK_CONVERT_PIVOTS
  if ((uintmax_t)MatrixRows <= SIZE_MAX / sizeof(sunlapackinttype))
  {
    LAPACK_PIVOTS(S) = malloc(MatrixRows * sizeof(sunlapackinttype));
  }
  if (LAPACK_PIVOTS(S) == NULL)
  {
    SUNLinSolFree(S);
    return NULL;
  }
#endif

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
  sunlapackinttype ier;

  /* check for valid inputs */
  if ((A == NULL) || (S == NULL)) { return SUN_ERR_ARG_CORRUPT; }

  /* Ensure that A is a dense matrix */
  if (SUNMatGetID(A) != SUNMATRIX_DENSE)
  {
    LASTFLAG(S) = SUN_ERR_ARG_INCOMPATIBLE;
    return SUN_ERR_ARG_INCOMPATIBLE;
  }

  /* Call LAPACK to do LU factorization of A */
  sunindextype rows = SUNDenseMatrix_Rows(A);
  if (rows != LAPACKDENSE_CONTENT(S)->N || SUNDenseMatrix_Columns(A) != rows ||
      !sunLapackIndexInRange(rows))
  {
    LASTFLAG(S) = SUN_ERR_ARG_INCOMPATIBLE;
    return SUN_ERR_ARG_INCOMPATIBLE;
  }
  sunlapackinttype n = (sunlapackinttype)rows;
  ier = 0;
  xgetrf_f77(&n, &n, SUNDenseMatrix_Data(A), &n, LAPACK_PIVOTS(S), &ier);
#if SUNDIALS_LAPACK_CONVERT_PIVOTS
  if (ier >= 0)
  {
    for (sunindextype i = 0; i < rows; i++)
    {
      PIVOTS(S)[i] = (sunindextype)LAPACK_PIVOTS(S)[i];
    }
  }
#endif
  LASTFLAG(S) = ier;
  if (ier > 0) { return (SUNLS_LUFACT_FAIL); }
  if (ier < 0) { return SUN_ERR_EXT_FAIL; }
  return SUN_SUCCESS;
}

int SUNLinSolSolve_LapackDense(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                               N_Vector b, SUNDIALS_MAYBE_UNUSED sunrealtype tol)
{
  sunlapackinttype one, ier;
  sunrealtype* xdata;

  if ((A == NULL) || (S == NULL) || (x == NULL) || (b == NULL))
  {
    return SUN_ERR_ARG_CORRUPT;
  }

  if (SUNMatGetID(A) != SUNMATRIX_DENSE)
  {
    LASTFLAG(S) = SUN_ERR_ARG_INCOMPATIBLE;
    return SUN_ERR_ARG_INCOMPATIBLE;
  }

  /* Call LAPACK to solve the linear system */
  sunindextype rows = SUNDenseMatrix_Rows(A);
  if (rows != LAPACKDENSE_CONTENT(S)->N || SUNDenseMatrix_Columns(A) != rows ||
      N_VGetLength(x) != rows || N_VGetLength(b) != rows || !sunLapackIndexInRange(rows))
  {
    LASTFLAG(S) = SUN_ERR_ARG_INCOMPATIBLE;
    return SUN_ERR_ARG_INCOMPATIBLE;
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

  sunlapackinttype n = (sunlapackinttype)rows;
  one = 1;
  ier = 0;
#if SUNDIALS_LAPACK_CONVERT_PIVOTS
  for (sunindextype i = 0; i < rows; i++)
  {
    sunindextype pivot = PIVOTS(S)[i];
    if (pivot < 1 || pivot > rows)
    {
      LASTFLAG(S) = SUN_ERR_ARG_OUTOFRANGE;
      return SUN_ERR_ARG_OUTOFRANGE;
    }
    LAPACK_PIVOTS(S)[i] = (sunlapackinttype)pivot;
  }
#endif
  xgetrs_f77("N", &n, &one, SUNDenseMatrix_Data(A), &n, LAPACK_PIVOTS(S),
             xdata, &n, &ier);
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
  const long int arrays = 1 + SUNDIALS_LAPACK_CONVERT_PIVOTS;
  const sunindextype n = LAPACKDENSE_CONTENT(S)->N;
  if ((uintmax_t)n > (uintmax_t)((LONG_MAX - 2) / arrays))
  {
    return SUN_ERR_OUTOFRANGE;
  }
  *lenrwLS = 0;
  *leniwLS = 2 + arrays * (long int)n;
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
#if SUNDIALS_LAPACK_CONVERT_PIVOTS
    free(LAPACK_PIVOTS(S));
#endif
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
