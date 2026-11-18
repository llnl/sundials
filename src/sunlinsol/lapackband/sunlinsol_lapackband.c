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
 * This is the implementation file for the LAPACK band
 * implementation of the SUNLINSOL package.
 * -----------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_math.h>
#include <sunlinsol/sunlinsol_lapackband.h>

#include "sundials_lapack_defs.h"
#include "sundials_macros.h"

/* Interfaces to match 'sunrealtype' with the correct LAPACK functions */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define xgbtrf_f77 dgbtrf_f77
#define xgbtrs_f77 dgbtrs_f77
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define xgbtrf_f77 sgbtrf_f77
#define xgbtrs_f77 sgbtrs_f77
#else
#error Incompatible sunrealtype for LAPACK; disable LAPACK and rebuild
#endif

#define ONE SUN_RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * Band solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define LAPACKBAND_CONTENT(S) ((SUNLinearSolverContent_LapackBand)(S->content))
#define PIVOTS(S)             (LAPACKBAND_CONTENT(S)->pivots)
#define LASTFLAG(S)           (LAPACKBAND_CONTENT(S)->last_flag)

#if SUNDIALS_LAPACK_CONVERT_PIVOTS
typedef struct
{
  struct _SUNLinearSolverContent_LapackBand content;
  sunlapackinttype* pivots;
} SunLapackBandStorage;
#define LAPACK_PIVOTS(S) (((SunLapackBandStorage*)(S->content))->pivots)
#else
#define LAPACK_PIVOTS(S) PIVOTS(S)
#endif

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new LAPACK band linear solver
 */

SUNLinearSolver SUNLinSol_LapackBand(N_Vector y, SUNMatrix A, SUNContext sunctx)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_LapackBand content;
  sunindextype MatrixRows;

  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_BAND) { return (NULL); }

  if (SUNBandMatrix_Rows(A) != SUNBandMatrix_Columns(A)) { return (NULL); }

  if ((N_VGetVectorID(y) != SUNDIALS_NVEC_SERIAL) &&
      (N_VGetVectorID(y) != SUNDIALS_NVEC_OPENMP) &&
      (N_VGetVectorID(y) != SUNDIALS_NVEC_PTHREADS))
  {
    return (NULL);
  }

  MatrixRows = SUNBandMatrix_Rows(A);
  if (MatrixRows != N_VGetLength(y)) { return (NULL); }

  if (!sunLapackIndexInRange(MatrixRows) ||
      (uintmax_t)MatrixRows > SIZE_MAX / sizeof(sunindextype))
  {
    return NULL;
  }
  if (!sunLapackIndexInRange(SUNBandMatrix_LDim(A))) { return NULL; }

  /* Create an empty linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty(sunctx);
  if (S == NULL) { return (NULL); }

  /* Attach operations */
  S->ops->gettype    = SUNLinSolGetType_LapackBand;
  S->ops->getid      = SUNLinSolGetID_LapackBand;
  S->ops->initialize = SUNLinSolInitialize_LapackBand;
  S->ops->setup      = SUNLinSolSetup_LapackBand;
  S->ops->solve      = SUNLinSolSolve_LapackBand;
  S->ops->lastflag   = SUNLinSolLastFlag_LapackBand;
  S->ops->space      = SUNLinSolSpace_LapackBand;
  S->ops->free       = SUNLinSolFree_LapackBand;

  /* Create content */
  content = NULL;
#if SUNDIALS_LAPACK_CONVERT_PIVOTS
  SunLapackBandStorage* storage = malloc(sizeof *storage);
  if (storage != NULL)
  {
    storage->pivots = NULL;
    content = &storage->content;
  }
#else
  content = (SUNLinearSolverContent_LapackBand)malloc(sizeof *content);
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

SUNLinearSolver_Type SUNLinSolGetType_LapackBand(
  SUNDIALS_MAYBE_UNUSED SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_DIRECT);
}

SUNLinearSolver_ID SUNLinSolGetID_LapackBand(SUNDIALS_MAYBE_UNUSED SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_LAPACKBAND);
}

SUNErrCode SUNLinSolInitialize_LapackBand(SUNLinearSolver S)
{
  /* all solver-specific memory has already been allocated */
  LASTFLAG(S) = SUN_SUCCESS;
  return SUN_SUCCESS;
}

int SUNLinSolSetup_LapackBand(SUNLinearSolver S, SUNMatrix A)
{
  sunlapackinttype ier;

  /* check for valid inputs */
  if ((A == NULL) || (S == NULL)) { return SUN_ERR_ARG_CORRUPT; }

  /* Ensure that A is a band matrix */
  if (SUNMatGetID(A) != SUNMATRIX_BAND)
  {
    LASTFLAG(S) = SUN_ERR_ARG_INCOMPATIBLE;
    return SUN_ERR_ARG_INCOMPATIBLE;
  }

  /* Call LAPACK to do LU factorization of A */
  ier  = 0;
  sunindextype rows = SUNBandMatrix_Rows(A);
  sunindextype lower = SUNBandMatrix_LowerBandwidth(A);
  sunindextype upper = SUNBandMatrix_UpperBandwidth(A);
  sunindextype leading_dim = SUNBandMatrix_LDim(A);
  if (rows != LAPACKBAND_CONTENT(S)->N || !sunLapackIndexInRange(rows) ||
      !sunLapackIndexInRange(lower) || !sunLapackIndexInRange(upper) ||
      !sunLapackIndexInRange(leading_dim))
  {
    LASTFLAG(S) = SUN_ERR_ARG_INCOMPATIBLE;
    return SUN_ERR_ARG_INCOMPATIBLE;
  }
  sunlapackinttype n = (sunlapackinttype)rows;
  sunlapackinttype ml = (sunlapackinttype)lower;
  sunlapackinttype mu = (sunlapackinttype)upper;
  sunlapackinttype ldim = (sunlapackinttype)leading_dim;
  xgbtrf_f77(&n, &n, &ml, &mu, SUNBandMatrix_Data(A), &ldim, LAPACK_PIVOTS(S), &ier);

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

int SUNLinSolSolve_LapackBand(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                              N_Vector b, SUNDIALS_MAYBE_UNUSED sunrealtype tol)
{
  sunlapackinttype one, ier;
  sunrealtype* xdata;

  /* check for valid inputs */
  if ((A == NULL) || (S == NULL) || (x == NULL) || (b == NULL))
  {
    return SUN_ERR_ARG_CORRUPT;
  }

  if (SUNMatGetID(A) != SUNMATRIX_BAND)
  {
    LASTFLAG(S) = SUN_ERR_ARG_INCOMPATIBLE;
    return SUN_ERR_ARG_INCOMPATIBLE;
  }

  /* Call LAPACK to solve the linear system */
  ier  = 0;
  sunindextype rows = SUNBandMatrix_Rows(A);
  sunindextype lower = SUNBandMatrix_LowerBandwidth(A);
  sunindextype upper = SUNBandMatrix_UpperBandwidth(A);
  sunindextype leading_dim = SUNBandMatrix_LDim(A);
  if (rows != LAPACKBAND_CONTENT(S)->N || N_VGetLength(x) != rows ||
      N_VGetLength(b) != rows || !sunLapackIndexInRange(rows) ||
      !sunLapackIndexInRange(lower) || !sunLapackIndexInRange(upper) ||
      !sunLapackIndexInRange(leading_dim))
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
  sunlapackinttype ml = (sunlapackinttype)lower;
  sunlapackinttype mu = (sunlapackinttype)upper;
  sunlapackinttype ldim = (sunlapackinttype)leading_dim;
  one  = 1;
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
  xgbtrs_f77("N", &n, &ml, &mu, &one, SUNBandMatrix_Data(A), &ldim, LAPACK_PIVOTS(S),
             xdata, &n, &ier);
  LASTFLAG(S) = ier;
  if (ier < 0) { return SUN_ERR_EXT_FAIL; }

  LASTFLAG(S) = SUN_SUCCESS;
  return SUN_SUCCESS;
}

sunindextype SUNLinSolLastFlag_LapackBand(SUNLinearSolver S)
{
  return (LASTFLAG(S));
}

SUNErrCode SUNLinSolSpace_LapackBand(SUNLinearSolver S, long int* lenrwLS,
                                     long int* leniwLS)
{
  const long int arrays = 1 + SUNDIALS_LAPACK_CONVERT_PIVOTS;
  const sunindextype n = LAPACKBAND_CONTENT(S)->N;
  if ((uintmax_t)n > (uintmax_t)((LONG_MAX - 2) / arrays))
  {
    return SUN_ERR_OUTOFRANGE;
  }
  *lenrwLS = 0;
  *leniwLS = 2 + arrays * (long int)n;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolFree_LapackBand(SUNLinearSolver S)
{
  /* return with success if already freed */
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
