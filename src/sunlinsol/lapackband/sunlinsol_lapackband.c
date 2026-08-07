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
  content = (SUNLinearSolverContent_LapackBand)malloc(sizeof *content);
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
  sunindextype n, ml, mu, ldim, ier;

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
  n    = SUNBandMatrix_Rows(A);
  ml   = SUNBandMatrix_LowerBandwidth(A);
  mu   = SUNBandMatrix_UpperBandwidth(A);
  ldim = SUNBandMatrix_LDim(A);
  xgbtrf_f77(&n, &n, &ml, &mu, SUNBandMatrix_Data(A), &ldim, PIVOTS(S), &ier);

  LASTFLAG(S) = ier;
  if (ier > 0) { return (SUNLS_LUFACT_FAIL); }
  if (ier < 0) { return SUN_ERR_EXT_FAIL; }
  return SUN_SUCCESS;
}

int SUNLinSolSolve_LapackBand(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                              N_Vector b, SUNDIALS_MAYBE_UNUSED sunrealtype tol)
{
  sunindextype n, ml, mu, ldim, one, ier;
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
    return SUN_ERR_MEM_FAIL;
  }

  /* Call LAPACK to solve the linear system */
  ier  = 0;
  n    = SUNBandMatrix_Rows(A);
  ml   = SUNBandMatrix_LowerBandwidth(A);
  mu   = SUNBandMatrix_UpperBandwidth(A);
  ldim = SUNBandMatrix_LDim(A);
  one  = 1;
  xgbtrs_f77("N", &n, &ml, &mu, &one, SUNBandMatrix_Data(A), &ldim, PIVOTS(S),
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
  *lenrwLS = 0;
  *leniwLS = 2 + LAPACKBAND_CONTENT(S)->N;
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
