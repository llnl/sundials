/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ UMBC
 *                David Gardner @ LLNL
 * Based on code sundials_band.c by: Alan C. Hindmarsh and
 *    Radu Serban @ LLNL
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
 * the SUNMATRIX package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_math.h>
#include <sunmatrix/sunmatrix_band.h>

#include "sundials_macros.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/* Private function prototypes */
static sunbooleantype compatibleMatrices(SUNMatrix A, SUNMatrix B);
static sunbooleantype compatibleMatrixAndVectors(SUNMatrix A, N_Vector x,
                                                 N_Vector y);
static SUNErrCode SMScaleAddNew_Band(sunrealtype c, SUNMatrix A, SUNMatrix B);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new band matrix with default storage upper bandwidth
 */

SUNMatrix SUNBandMatrix(sunindextype N, sunindextype mu, sunindextype ml,
                        SUNContext sunctx)
{
  return (SUNBandMatrixStorage(N, mu, ml, mu + ml, sunctx));
}

/* ----------------------------------------------------------------------------
 * Function to create a new band matrix with specified storage upper bandwidth
 */

SUNMatrix SUNBandMatrixStorage(sunindextype N, sunindextype mu, sunindextype ml,
                               sunindextype smu, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNMatrix A;
  SUNMatrixContent_Band content;
  sunindextype j, colSize;

  SUNAssertNull(N > 0, SUN_ERR_ARG_OUTOFRANGE);
  SUNAssertNull(smu >= 0, SUN_ERR_ARG_OUTOFRANGE);
  SUNAssertNull(ml >= 0, SUN_ERR_ARG_OUTOFRANGE);

  /* Create an empty matrix object */
  A = NULL;
  A = SUNMatNewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  A->ops->getid                    = SUNMatGetID_Band;
  A->ops->clone                    = SUNMatClone_Band;
  A->ops->destroy                  = SUNMatDestroy_Band;
  A->ops->zero                     = SUNMatZero_Band;
  A->ops->copy                     = SUNMatCopy_Band;
  A->ops->scaleadd                 = SUNMatScaleAdd_Band;
  A->ops->scaleaddi                = SUNMatScaleAddI_Band;
  A->ops->matvec                   = SUNMatMatvec_Band;
  A->ops->mathermitiantransposevec = SUNMatHermitianTransposeVec_Band;
  A->ops->space                    = SUNMatSpace_Band;

  /* Create content */
  content = NULL;
  content = (SUNMatrixContent_Band)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  A->content = content;

  /* Fill content */
  colSize        = smu + ml + 1;
  content->M     = N;
  content->N     = N;
  content->mu    = mu;
  content->ml    = ml;
  content->s_mu  = smu;
  content->ldim  = colSize;
  content->ldata = N * colSize;
  content->data  = NULL;
  content->cols  = NULL;

  /* Allocate content */
  content->data = (sunrealtype*)calloc(N * colSize, sizeof(sunrealtype));
  SUNAssertNull(content->data, SUN_ERR_MALLOC_FAIL);

  content->cols = (sunrealtype**)malloc(N * sizeof(sunrealtype*));
  SUNAssertNull(content->cols, SUN_ERR_MALLOC_FAIL);
  for (j = 0; j < N; j++) { content->cols[j] = content->data + j * colSize; }

  return (A);
}

/* ----------------------------------------------------------------------------
 * Function to print the band matrix
 */

void SUNBandMatrix_Print(SUNMatrix A, FILE* outfile)
{
  SUNFunctionBegin(A->sunctx);
  sunindextype i, j, start, finish;

  SUNAssertVoid(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);

  /* perform operation */
  fprintf(outfile, "\n");
  for (i = 0; i < SM_ROWS_B(A); i++)
  {
    start  = SUNMAX(0, i - SM_LBAND_B(A));
    finish = SUNMIN(SM_COLUMNS_B(A) - 1, i + SM_UBAND_B(A));
    for (j = 0; j < start; j++) { fprintf(outfile, "%12s  ", ""); }
    for (j = start; j <= finish; j++)
    {
      fprintf(outfile, SUN_FORMAT_E "  ", SM_ELEMENT_B(A, i, j));
    }
    fprintf(outfile, "\n");
  }
  return;
}

/* ----------------------------------------------------------------------------
 * Functions to access the contents of the band matrix structure
 */

sunindextype SUNBandMatrix_Rows(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNoRet(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
  return SM_ROWS_B(A);
}

sunindextype SUNBandMatrix_Columns(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNoRet(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
  return SM_COLUMNS_B(A);
}

sunindextype SUNBandMatrix_LowerBandwidth(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNoRet(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
  return SM_LBAND_B(A);
}

sunindextype SUNBandMatrix_UpperBandwidth(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNoRet(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
  return SM_UBAND_B(A);
}

sunindextype SUNBandMatrix_StoredUpperBandwidth(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNoRet(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
  return SM_SUBAND_B(A);
}

sunindextype SUNBandMatrix_LDim(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNoRet(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
  return SM_LDIM_B(A);
}

sunindextype SUNBandMatrix_LData(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNoRet(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
  return SM_LDATA_B(A);
}

sunrealtype* SUNBandMatrix_Data(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNull(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
  return SM_DATA_B(A);
}

sunrealtype** SUNBandMatrix_Cols(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNull(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
  return SM_COLS_B(A);
}

sunrealtype* SUNBandMatrix_Column(SUNMatrix A, sunindextype j)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNull(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
  return SM_COLUMN_B(A, j);
}

/*
 * -----------------------------------------------------------------
 * implementation of matrix operations
 * -----------------------------------------------------------------
 */

SUNMatrix_ID SUNMatGetID_Band(SUNDIALS_MAYBE_UNUSED SUNMatrix A)
{
  return SUNMATRIX_BAND;
}

SUNMatrix SUNMatClone_Band(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNMatrix B = SUNBandMatrixStorage(SM_COLUMNS_B(A), SM_UBAND_B(A),
                                     SM_LBAND_B(A), SM_SUBAND_B(A), A->sunctx);
  SUNCheckLastErrNull();
  return (B);
}

void SUNMatDestroy_Band(SUNMatrix A)
{
  if (A == NULL) { return; }

  /* free content */
  if (A->content != NULL)
  {
    /* free data array */
    if (SM_DATA_B(A))
    {
      free(SM_DATA_B(A));
      SM_DATA_B(A) = NULL;
    }
    /* free column pointers */
    if (SM_COLS_B(A))
    {
      free(SM_COLS_B(A));
      SM_COLS_B(A) = NULL;
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

SUNErrCode SUNMatZero_Band(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  sunindextype i;
  sunrealtype* Adata;

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);

  /* Perform operation */
  Adata = SM_DATA_B(A);
  for (i = 0; i < SM_LDATA_B(A); i++) { Adata[i] = ZERO; }

  return SUN_SUCCESS;
}

SUNErrCode SUNMatCopy_Band(SUNMatrix A, SUNMatrix B)
{
  SUNFunctionBegin(A->sunctx);
  sunindextype i, j, colSize, ml, mu, smu;
  sunrealtype *A_colj, *B_colj;

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
  SUNAssert(SUNMatGetID(B) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
  SUNCheck(compatibleMatrices(A, B), SUN_ERR_ARG_DIMSMISMATCH);

  /* Grow B if A's bandwidth is larger */
  if ((SM_UBAND_B(A) > SM_UBAND_B(B)) || (SM_LBAND_B(A) > SM_LBAND_B(B)))
  {
    ml                     = SUNMAX(SM_LBAND_B(B), SM_LBAND_B(A));
    mu                     = SUNMAX(SM_UBAND_B(B), SM_UBAND_B(A));
    smu                    = SUNMAX(SM_SUBAND_B(B), SM_SUBAND_B(A));
    colSize                = smu + ml + 1;
    SM_CONTENT_B(B)->mu    = mu;
    SM_CONTENT_B(B)->ml    = ml;
    SM_CONTENT_B(B)->s_mu  = smu;
    SM_CONTENT_B(B)->ldim  = colSize;
    SM_CONTENT_B(B)->ldata = SM_COLUMNS_B(B) * colSize;
    SM_CONTENT_B(B)->data =
      (sunrealtype*)realloc(SM_CONTENT_B(B)->data,
                            SM_COLUMNS_B(B) * colSize * sizeof(sunrealtype));
    for (j = 0; j < SM_COLUMNS_B(B); j++)
    {
      SM_CONTENT_B(B)->cols[j] = SM_CONTENT_B(B)->data + j * colSize;
    }
  }

  /* Perform operation */
  SUNCheckCall(SUNMatZero_Band(B));
  for (j = 0; j < SM_COLUMNS_B(B); j++)
  {
    B_colj = SM_COLUMN_B(B, j);
    A_colj = SM_COLUMN_B(A, j);
    for (i = -SM_UBAND_B(A); i <= SM_LBAND_B(A); i++) { B_colj[i] = A_colj[i]; }
  }
  return SUN_SUCCESS;
}

SUNErrCode SUNMatScaleAddI_Band(sunrealtype c, SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  sunindextype i, j;
  sunrealtype* A_colj;

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);

  /* Perform operation */
  for (j = 0; j < SM_COLUMNS_B(A); j++)
  {
    A_colj = SM_COLUMN_B(A, j);
    for (i = -SM_UBAND_B(A); i <= SM_LBAND_B(A); i++) { A_colj[i] *= c; }
    SM_ELEMENT_B(A, j, j) += ONE;
  }
  return SUN_SUCCESS;
}

SUNErrCode SUNMatScaleAdd_Band(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  SUNFunctionBegin(A->sunctx);
  sunindextype i, j;
  sunrealtype *A_colj, *B_colj;

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
  SUNAssert(SUNMatGetID(B) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
  SUNCheck(compatibleMatrices(A, B), SUN_ERR_ARG_DIMSMISMATCH);

  /* Call separate routine in B has larger bandwidth(s) than A */
  if ((SM_UBAND_B(B) > SM_UBAND_B(A)) || (SM_LBAND_B(B) > SM_LBAND_B(A)))
  {
    SUNCheckCall(SMScaleAddNew_Band(c, A, B));
    return SUN_SUCCESS;
  }

  /* Otherwise, perform operation in-place */
  for (j = 0; j < SM_COLUMNS_B(A); j++)
  {
    A_colj = SM_COLUMN_B(A, j);
    B_colj = SM_COLUMN_B(B, j);
    for (i = -SM_UBAND_B(B); i <= SM_LBAND_B(B); i++)
    {
      A_colj[i] = c * A_colj[i] + B_colj[i];
    }
  }
  return SUN_SUCCESS;
}

SUNErrCode SUNMatMatvec_Band(SUNMatrix A, N_Vector x, N_Vector y)
{
  SUNFunctionBegin(A->sunctx);
  sunindextype i, j, is, ie;
  sunrealtype *col_j, *xd, *yd;

  SUNCheck(compatibleMatrixAndVectors(A, x, y), SUN_ERR_ARG_DIMSMISMATCH);

  /* access vector data (return if failure) */
  xd = N_VGetArrayPointer(x);
  SUNCheckLastErr();
  yd = N_VGetArrayPointer(y);
  SUNCheckLastErr();

  SUNAssert(xd, SUN_ERR_MEM_FAIL);
  SUNAssert(yd, SUN_ERR_MEM_FAIL);

  /* Perform operation */
  for (i = 0; i < SM_ROWS_B(A); i++) { yd[i] = ZERO; }
  for (j = 0; j < SM_COLUMNS_B(A); j++)
  {
    col_j = SM_COLUMN_B(A, j);
    is    = SUNMAX(0, j - SM_UBAND_B(A));
    ie    = SUNMIN(SM_ROWS_B(A) - 1, j + SM_LBAND_B(A));
    for (i = is; i <= ie; i++) { yd[i] += col_j[i - j] * xd[j]; }
  }
  return SUN_SUCCESS;
}

SUNErrCode SUNMatHermitianTransposeVec_Band(SUNMatrix A, N_Vector x, N_Vector y)
{
  SUNFunctionBegin(A->sunctx);
  sunindextype i, j, is, ie;
  sunrealtype *col_j, *xd, *yd;

  SUNCheck(compatibleMatrixAndVectors(A, y, x), SUN_ERR_ARG_DIMSMISMATCH);

  /* access vector data (return if failure) */
  xd = N_VGetArrayPointer(x);
  SUNCheckLastErr();
  yd = N_VGetArrayPointer(y);
  SUNCheckLastErr();

  SUNAssert(xd, SUN_ERR_MEM_FAIL);
  SUNAssert(yd, SUN_ERR_MEM_FAIL);

  /* Perform operation */
  for (i = 0; i < SM_ROWS_B(A); i++) { yd[i] = ZERO; }
  for (j = 0; j < SM_COLUMNS_B(A); j++)
  {
    col_j = SM_COLUMN_B(A, j);
    is    = SUNMAX(0, j - SM_UBAND_B(A));
    ie    = SUNMIN(SM_ROWS_B(A) - 1, j + SM_LBAND_B(A));
    for (i = is; i <= ie; i++) { yd[j] += col_j[i - j] * xd[i]; }
  }
  return SUN_SUCCESS;
}

SUNErrCode SUNMatSpace_Band(SUNMatrix A, long int* lenrw, long int* leniw)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssert(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
  SUNAssert(lenrw, SUN_ERR_ARG_CORRUPT);
  SUNAssert(leniw, SUN_ERR_ARG_CORRUPT);
  *lenrw = SM_COLUMNS_B(A) * (SM_SUBAND_B(A) + SM_LBAND_B(A) + 1);
  *leniw = 7 + SM_COLUMNS_B(A);
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
  /* both matrices must have the same number of columns
     (note that we do not check for identical bandwidth) */
  if (SM_ROWS_B(A) != SM_ROWS_B(B)) { return SUNFALSE; }
  if (SM_COLUMNS_B(A) != SM_COLUMNS_B(B)) { return SUNFALSE; }
  return SUNTRUE;
}

SUNDIALS_MAYBE_UNUSED
static sunbooleantype compatibleMatrixAndVectors(SUNMatrix A, N_Vector x,
                                                 N_Vector y)
{
  /* Vectors must provide nvgetarraypointer */
  if (!x->ops->nvgetarraypointer || !y->ops->nvgetarraypointer)
  {
    return SUNFALSE;
  }

  /* Check that the dimensions agree */
  if ((N_VGetLength(x) != SM_COLUMNS_B(A)) || (N_VGetLength(y) != SM_ROWS_B(A)))
  {
    return SUNFALSE;
  }

  return SUNTRUE;
}

SUNErrCode SMScaleAddNew_Band(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  SUNFunctionBegin(A->sunctx);
  sunindextype i, j, ml, mu, smu;
  sunrealtype *A_colj, *B_colj, *C_colj;
  SUNMatrix C;

  /* create new matrix large enough to hold both A and B */
  ml  = SUNMAX(SM_LBAND_B(A), SM_LBAND_B(B));
  mu  = SUNMAX(SM_UBAND_B(A), SM_UBAND_B(B));
  smu = SUNMIN(SM_COLUMNS_B(A) - 1, mu + ml);
  C   = SUNBandMatrixStorage(SM_COLUMNS_B(A), mu, ml, smu, A->sunctx);
  SUNCheckLastErr();

  /* scale/add c*A into new matrix */
  for (j = 0; j < SM_COLUMNS_B(A); j++)
  {
    A_colj = SM_COLUMN_B(A, j);
    C_colj = SM_COLUMN_B(C, j);
    for (i = -SM_UBAND_B(A); i <= SM_LBAND_B(A); i++)
    {
      C_colj[i] = c * A_colj[i];
    }
  }

  /* add B into new matrix */
  for (j = 0; j < SM_COLUMNS_B(B); j++)
  {
    B_colj = SM_COLUMN_B(B, j);
    C_colj = SM_COLUMN_B(C, j);
    for (i = -SM_UBAND_B(B); i <= SM_LBAND_B(B); i++)
    {
      C_colj[i] += B_colj[i];
    }
  }

  /* replace A contents with C contents, nullify C content pointer, destroy C */
  free(SM_DATA_B(A));
  SM_DATA_B(A) = NULL;
  free(SM_COLS_B(A));
  SM_COLS_B(A) = NULL;
  free(A->content);
  A->content = NULL;
  A->content = C->content;
  C->content = NULL;
  SUNMatDestroy_Band(C);

  return SUN_SUCCESS;
}
