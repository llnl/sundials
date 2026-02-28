/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ UMBC
 *                David Gardner @ LLNL
 * Based on code sundials_dense.c by: Scott D. Cohen,
 *     Alan C. Hindmarsh and Radu Serban @ LLNL
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
 * This is the implementation file for the dense implementation of
 * the SUNMATRIX package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_errors.h>
#include <sunmatrix/sunmatrix_dense.h>

#include "sundials_macros.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/* Private function prototypes */
static sunbooleantype compatibleMatrices(SUNMatrix A, SUNMatrix B);
static sunbooleantype compatibleMatrixAndVectors(SUNMatrix A, N_Vector x,
                                                 N_Vector y);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new dense matrix
 */

SUNMatrix SUNDenseMatrix(sunindextype M, sunindextype N, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNMatrix A;
  SUNMatrixContent_Dense content;
  sunindextype j;

  /* return with NULL matrix on illegal dimension input */
  SUNAssertNull(N > 0 && M > 0, SUN_ERR_ARG_OUTOFRANGE);

  /* Create an empty matrix object */
  A = NULL;
  A = SUNMatNewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  A->ops->getid                    = SUNMatGetID_Dense;
  A->ops->clone                    = SUNMatClone_Dense;
  A->ops->destroy                  = SUNMatDestroy_Dense;
  A->ops->zero                     = SUNMatZero_Dense;
  A->ops->copy                     = SUNMatCopy_Dense;
  A->ops->scaleadd                 = SUNMatScaleAdd_Dense;
  A->ops->scaleaddi                = SUNMatScaleAddI_Dense;
  A->ops->matvec                   = SUNMatMatvec_Dense;
  A->ops->mathermitiantransposevec = SUNMatHermitianTransposeVec_Dense;
  A->ops->space                    = SUNMatSpace_Dense;

  /* Create content */
  content = NULL;
  content = (SUNMatrixContent_Dense)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  A->content = content;

  /* Fill content */
  content->M     = M;
  content->N     = N;
  content->ldata = M * N;
  content->data  = NULL;
  content->cols  = NULL;

  /* Allocate content */
  content->data = (sunrealtype*)calloc(M * N, sizeof(sunrealtype));
  SUNAssertNull(content->data, SUN_ERR_MALLOC_FAIL);

  content->cols = (sunrealtype**)malloc(N * sizeof(sunrealtype*));
  SUNAssertNull(content->cols, SUN_ERR_MALLOC_FAIL);
  for (j = 0; j < N; j++) { content->cols[j] = content->data + j * M; }

  return (A);
}

/* ----------------------------------------------------------------------------
 * Function to print the dense matrix
 */

void SUNDenseMatrix_Print(SUNMatrix A, FILE* outfile)
{
  SUNFunctionBegin(A->sunctx);
  sunindextype i, j;

  SUNAssertVoid(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE);

  /* perform operation */
  fprintf(outfile, "\n");
  for (i = 0; i < SM_ROWS_D(A); i++)
  {
    for (j = 0; j < SM_COLUMNS_D(A); j++)
    {
      fprintf(outfile, SUN_FORMAT_E "  ", SM_ELEMENT_D(A, i, j));
    }
    fprintf(outfile, "\n");
  }
  return;
}

/* ----------------------------------------------------------------------------
 * Functions to access the contents of the dense matrix structure
 */

sunindextype SUNDenseMatrix_Rows(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNoRet(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE);
  return SM_ROWS_D(A);
}

sunindextype SUNDenseMatrix_Columns(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNoRet(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE);
  return SM_COLUMNS_D(A);
}

sunindextype SUNDenseMatrix_LData(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNoRet(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE);
  return SM_LDATA_D(A);
}

sunrealtype* SUNDenseMatrix_Data(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNull(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE);
  return SM_DATA_D(A);
}

sunrealtype** SUNDenseMatrix_Cols(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNull(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE);
  return SM_COLS_D(A);
}

sunrealtype* SUNDenseMatrix_Column(SUNMatrix A, sunindextype j)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNull(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE);
  return SM_COLUMN_D(A, j);
}

/*
 * -----------------------------------------------------------------
 * implementation of matrix operations
 * -----------------------------------------------------------------
 */

SUNMatrix_ID SUNMatGetID_Dense(SUNDIALS_MAYBE_UNUSED SUNMatrix A)
{
  return SUNMATRIX_DENSE;
}

SUNMatrix SUNMatClone_Dense(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNMatrix B = SUNDenseMatrix(SM_ROWS_D(A), SM_COLUMNS_D(A), A->sunctx);
  SUNCheckLastErrNull();
  return (B);
}

void SUNMatDestroy_Dense(SUNMatrix A)
{
  if (A == NULL) { return; }

  /* free content */
  if (A->content != NULL)
  {
    /* free data array */
    if (SM_DATA_D(A) != NULL)
    {
      free(SM_DATA_D(A));
      SM_DATA_D(A) = NULL;
    }
    /* free column pointers */
    if (SM_CONTENT_D(A)->cols != NULL)
    {
      free(SM_CONTENT_D(A)->cols);
      SM_CONTENT_D(A)->cols = NULL;
    }
    /* free content struct */
    free(A->content);
    A->content = NULL;
  }

  /* free ops and matrix */
  if (A->ops)
  {
    free(A->ops);
    A->ops = NULL;
  }
  free(A);
  A = NULL;

  return;
}

SUNErrCode SUNMatZero_Dense(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  sunindextype i;
  sunrealtype* Adata;

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE);

  /* Perform operation A_ij = 0 */
  Adata = SM_DATA_D(A);
  for (i = 0; i < SM_LDATA_D(A); i++) { Adata[i] = ZERO; }

  return SUN_SUCCESS;
}

SUNErrCode SUNMatCopy_Dense(SUNMatrix A, SUNMatrix B)
{
  SUNFunctionBegin(A->sunctx);
  sunindextype i, j;

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE);
  SUNAssert(SUNMatGetID(B) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE);
  SUNCheck(compatibleMatrices(A, B), SUN_ERR_ARG_DIMSMISMATCH);

  /* Perform operation B_ij = A_ij */
  for (j = 0; j < SM_COLUMNS_D(A); j++)
  {
    for (i = 0; i < SM_ROWS_D(A); i++)
    {
      SM_ELEMENT_D(B, i, j) = SM_ELEMENT_D(A, i, j);
    }
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNMatScaleAddI_Dense(sunrealtype c, SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  sunindextype i, j;

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE);

  /* Perform operation A = c*A + I */
  for (j = 0; j < SM_COLUMNS_D(A); j++)
  {
    for (i = 0; i < SM_ROWS_D(A); i++)
    {
      SM_ELEMENT_D(A, i, j) *= c;
      if (i == j) { SM_ELEMENT_D(A, i, j) += ONE; }
    }
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNMatScaleAdd_Dense(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  SUNFunctionBegin(A->sunctx);
  sunindextype i, j;

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE);
  SUNCheck(compatibleMatrices(A, B), SUN_ERR_ARG_DIMSMISMATCH);

  /* Perform operation A = c*A + B */
  for (j = 0; j < SM_COLUMNS_D(A); j++)
  {
    for (i = 0; i < SM_ROWS_D(A); i++)
    {
      SM_ELEMENT_D(A, i, j) = c * SM_ELEMENT_D(A, i, j) + SM_ELEMENT_D(B, i, j);
    }
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNMatMatvec_Dense(SUNMatrix A, N_Vector x, N_Vector y)
{
  SUNFunctionBegin(A->sunctx);
  sunrealtype *xd, *yd;

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE);
  SUNCheck(compatibleMatrixAndVectors(A, x, y), SUN_ERR_ARG_DIMSMISMATCH);

  /* access vector data (return if NULL data pointers) */
  xd = N_VGetArrayPointer(x);
  SUNCheckLastErr();
  yd = N_VGetArrayPointer(y);
  SUNCheckLastErr();

  SUNAssert(xd, SUN_ERR_MEM_FAIL);
  SUNAssert(yd, SUN_ERR_MEM_FAIL);
  SUNAssert(xd != yd, SUN_ERR_MEM_FAIL);

  /* Perform operation y = Ax */
  for (sunindextype i = 0; i < SM_ROWS_D(A); i++) { yd[i] = ZERO; }
  for (sunindextype j = 0; j < SM_COLUMNS_D(A); j++)
  {
    sunrealtype* col_j = SM_COLUMN_D(A, j);
    for (sunindextype i = 0; i < SM_ROWS_D(A); i++)
    {
      yd[i] += col_j[i] * xd[j];
    }
  }
  return SUN_SUCCESS;
}

SUNErrCode SUNMatHermitianTransposeVec_Dense(SUNMatrix A, N_Vector x, N_Vector y)
{
  SUNFunctionBegin(A->sunctx);

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE);
  SUNCheck(compatibleMatrixAndVectors(A, y, x), SUN_ERR_ARG_DIMSMISMATCH);

  /* access vector data (return if NULL data pointers) */
  sunrealtype *xd, *yd;
  xd = N_VGetArrayPointer(x);
  SUNCheckLastErr();
  yd = N_VGetArrayPointer(y);
  SUNCheckLastErr();

  SUNAssert(xd, SUN_ERR_MEM_FAIL);
  SUNAssert(yd, SUN_ERR_MEM_FAIL);
  SUNAssert(xd != yd, SUN_ERR_MEM_FAIL);

  /* Perform operation y = A^T x */
  for (sunindextype i = 0; i < SM_COLUMNS_D(A); i++) { yd[i] = ZERO; }
  for (sunindextype i = 0; i < SM_COLUMNS_D(A); i++)
  {
    sunrealtype* row_i = SM_COLUMN_D(A, i);
    for (sunindextype j = 0; j < SM_ROWS_D(A); j++)
    {
      yd[i] += row_i[j] * xd[j];
    }
  }
  return SUN_SUCCESS;
}

SUNErrCode SUNMatSpace_Dense(SUNMatrix A, long int* lenrw, long int* leniw)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssert(SUNMatGetID(A) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE);
  SUNAssert(lenrw, SUN_ERR_ARG_CORRUPT);
  SUNAssert(leniw, SUN_ERR_ARG_CORRUPT);
  *lenrw = SM_LDATA_D(A);
  *leniw = 3 + SM_COLUMNS_D(A);
  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

SUNDIALS_MAYBE_UNUSED
static sunbooleantype compatibleMatrices(SUNMatrix A, SUNMatrix B)
{
  /* both matrices must have the same shape */
  if ((SM_ROWS_D(A) != SM_ROWS_D(B)) || (SM_COLUMNS_D(A) != SM_COLUMNS_D(B)))
  {
    return SUNFALSE;
  }

  return SUNTRUE;
}

SUNDIALS_MAYBE_UNUSED
static sunbooleantype compatibleMatrixAndVectors(SUNMatrix A, N_Vector x,
                                                 N_Vector y)
{
  /* Vectors must provide nvgetarraypointer and cannot be a parallel vector */
  if (!x->ops->nvgetarraypointer || !y->ops->nvgetarraypointer)
  {
    return SUNFALSE;
  }

  /* Check that the dimensions agree */
  if ((N_VGetLength(x) != SM_COLUMNS_D(A)) || (N_VGetLength(y) != SM_ROWS_D(A)))
  {
    return SUNFALSE;
  }

  return SUNTRUE;
}
