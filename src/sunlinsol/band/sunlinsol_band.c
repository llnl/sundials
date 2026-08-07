/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ UMBC
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
 * This is the implementation file for the band implementation of
 * the SUNLINSOL package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_math.h>
#include <sunlinsol/sunlinsol_band.h>

#include "sundials_macros.h"

#define ONE SUN_RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * Band solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define BAND_CONTENT(S) ((SUNLinearSolverContent_Band)(S->content))
#define PIVOTS(S)       (BAND_CONTENT(S)->pivots)
#define LASTFLAG(S)     (BAND_CONTENT(S)->last_flag)

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new band linear solver
 */

SUNLinearSolver SUNLinSol_Band(SUNDIALS_MAYBE_UNUSED N_Vector y, SUNMatrix A,
                               SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNLinearSolver S;
  SUNLinearSolverContent_Band content;
  sunindextype MatrixRows;

  SUNAssertNull(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
  SUNAssertNull(SUNBandMatrix_Rows(A) == SUNBandMatrix_Columns(A),
                SUN_ERR_ARG_DIMSMISMATCH);
  SUNAssertNull(y->ops->nvgetarraypointer, SUN_ERR_ARG_INCOMPATIBLE);

  /* Check that A has appropriate storage upper bandwidth for factorization */
  MatrixRows = SUNBandMatrix_Rows(A);
  SUNAssertNull(SUNBandMatrix_StoredUpperBandwidth(A) >=
                  SUNMIN(MatrixRows - 1, SUNBandMatrix_LowerBandwidth(A) +
                                           SUNBandMatrix_UpperBandwidth(A)),
                SUN_ERR_ARG_INCOMPATIBLE);
  SUNAssertNull(MatrixRows == N_VGetLength(y), SUN_ERR_ARG_DIMSMISMATCH);

  /* Create an empty linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  S->ops->gettype    = SUNLinSolGetType_Band;
  S->ops->getid      = SUNLinSolGetID_Band;
  S->ops->initialize = SUNLinSolInitialize_Band;
  S->ops->setup      = SUNLinSolSetup_Band;
  S->ops->solve      = SUNLinSolSolve_Band;
  S->ops->lastflag   = SUNLinSolLastFlag_Band;
  S->ops->space      = SUNLinSolSpace_Band;
  S->ops->free       = SUNLinSolFree_Band;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_Band)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  S->content = content;

  /* Fill content */
  content->N         = MatrixRows;
  content->last_flag = 0;
  content->pivots    = NULL;

  /* Allocate content */
  content->pivots = (sunindextype*)malloc(MatrixRows * sizeof(sunindextype));
  SUNAssertNull(content->pivots, SUN_ERR_MALLOC_FAIL);

  return (S);
}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_Band(SUNDIALS_MAYBE_UNUSED SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_DIRECT);
}

SUNLinearSolver_ID SUNLinSolGetID_Band(SUNDIALS_MAYBE_UNUSED SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_BAND);
}

SUNErrCode SUNLinSolInitialize_Band(SUNLinearSolver S)
{
  /* all solver-specific memory has already been allocated */
  LASTFLAG(S) = SUN_SUCCESS;
  return SUN_SUCCESS;
}

int SUNLinSolSetup_Band(SUNLinearSolver S, SUNMatrix A)
{
  SUNFunctionBegin(S->sunctx);
  sunrealtype** A_cols;
  sunindextype* pivots;

  SUNAssert(A, SUN_ERR_ARG_CORRUPT);
  SUNAssert(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);

  /* access data pointers (return with failure on NULL) */
  A_cols = NULL;
  pivots = NULL;
  A_cols = SM_COLS_B(A);
  pivots = PIVOTS(S);
  SUNAssert(A_cols, SUN_ERR_ARG_CORRUPT);
  SUNAssert(pivots, SUN_ERR_ARG_CORRUPT);

  /* ensure that storage upper bandwidth is sufficient for fill-in */
  SUNAssert(SM_SUBAND_B(A) >=
              SUNMIN(SM_COLUMNS_B(A) - 1, SM_UBAND_B(A) + SM_LBAND_B(A)),
            SUN_ERR_ARG_INCOMPATIBLE);

  /* perform LU factorization of input matrix */
  LASTFLAG(S) = SUNDlsMat_bandGBTRF(A_cols, SM_COLUMNS_B(A), SM_UBAND_B(A),
                                    SM_LBAND_B(A), SM_SUBAND_B(A), pivots);

  /* store error flag (if nonzero, that row encountered zero-valued pivot) */
  if (LASTFLAG(S) > 0) { return (SUNLS_LUFACT_FAIL); }
  return SUN_SUCCESS;
}

int SUNLinSolSolve_Band(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b,
                        SUNDIALS_MAYBE_UNUSED sunrealtype tol)
{
  SUNFunctionBegin(S->sunctx);
  sunrealtype **A_cols, *xdata;
  sunindextype* pivots;

  /* copy b into x */
  N_VScale(ONE, b, x);
  SUNCheckLastErr();

  /* access data pointers (return with failure on NULL) */
  A_cols = NULL;
  xdata  = NULL;
  pivots = NULL;
  A_cols = SUNBandMatrix_Cols(A);
  SUNCheckLastErr();
  xdata = N_VGetArrayPointer(x);
  SUNCheckLastErr();
  pivots = PIVOTS(S);
  SUNAssert(pivots, SUN_ERR_ARG_CORRUPT);

  /* solve using LU factors */
  SUNDlsMat_bandGBTRS(A_cols, SM_COLUMNS_B(A), SM_SUBAND_B(A), SM_LBAND_B(A),
                      pivots, xdata);
  LASTFLAG(S) = SUN_SUCCESS;
  return SUN_SUCCESS;
}

sunindextype SUNLinSolLastFlag_Band(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  return LASTFLAG(S);
}

SUNErrCode SUNLinSolSpace_Band(SUNLinearSolver S, long int* lenrwLS,
                               long int* leniwLS)
{
  SUNFunctionBegin(S->sunctx);
  SUNAssert(SUNLinSolGetID(S) == SUNLINEARSOLVER_BAND, SUN_ERR_ARG_WRONGTYPE);
  *leniwLS = 2 + BAND_CONTENT(S)->N;
  *lenrwLS = 0;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolFree_Band(SUNLinearSolver S)
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
