/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ UMBC
 *                David Gardner @ LLNL
 *                Steven B. Roberts @ LLNL
 * Based on code sundials_sparse.c by: Carol Woodward and
 *     Slaven Peles @ LLNL, and Daniel R. Reynolds @ UMBC
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
 * This is the implementation file for the sparse implementation of
 * the SUNMATRIX package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_math.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>

#include "sundials_macros.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/* Private function prototypes */
static sunbooleantype compatibleMatrices(SUNMatrix A, SUNMatrix B);
static sunbooleantype compatibleMatrixAndVectors(SUNMatrix A, N_Vector x,
                                                 N_Vector y);
static SUNErrCode Matvec_SparseCSC(SUNMatrix A, N_Vector x, N_Vector y);
static SUNErrCode Matvec_SparseCSR(SUNMatrix A, N_Vector x, N_Vector y);
static SUNErrCode MatTransposeVec_SparseCSC(SUNMatrix A, N_Vector x, N_Vector y);
static SUNErrCode MatTransposeVec_SparseCSR(SUNMatrix A, N_Vector x, N_Vector y);
static SUNErrCode format_convert(const SUNMatrix A, SUNMatrix B);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/*
 * ==================================================================
 * Private function prototypes (functions working on SlsMat)
 * ==================================================================
 */

/* ----------------------------------------------------------------------------
 * Function to create a new sparse matrix
 */

SUNMatrix SUNSparseMatrix(sunindextype M, sunindextype N, sunindextype NNZ,
                          int sparsetype, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNMatrix A;
  SUNMatrixContent_Sparse content;

  /* return with NULL matrix on illegal input */
  SUNAssertNull(M > 0 && N > 0, SUN_ERR_ARG_OUTOFRANGE);
  SUNAssertNull(NNZ >= 0, SUN_ERR_ARG_OUTOFRANGE);
  SUNAssertNull(sparsetype == SUN_CSC_MAT || sparsetype == SUN_CSR_MAT,
                SUN_ERR_ARG_OUTOFRANGE);

  /* Create an empty matrix object */
  A = NULL;
  A = SUNMatNewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  A->ops->getid                    = SUNMatGetID_Sparse;
  A->ops->clone                    = SUNMatClone_Sparse;
  A->ops->destroy                  = SUNMatDestroy_Sparse;
  A->ops->zero                     = SUNMatZero_Sparse;
  A->ops->copy                     = SUNMatCopy_Sparse;
  A->ops->scaleadd                 = SUNMatScaleAdd_Sparse;
  A->ops->scaleaddi                = SUNMatScaleAddI_Sparse;
  A->ops->matvec                   = SUNMatMatvec_Sparse;
  A->ops->mathermitiantransposevec = SUNMatHermitianTransposeVec_Sparse;
  A->ops->space                    = SUNMatSpace_Sparse;

  /* Create content */
  content = NULL;
  content = (SUNMatrixContent_Sparse)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  A->content = content;

  /* Fill content */
  content->sparsetype = sparsetype;
  content->M          = M;
  content->N          = N;
  content->NNZ        = NNZ;
  switch (sparsetype)
  {
  case SUN_CSC_MAT:
    content->NP      = N;
    content->rowvals = &(content->indexvals);
    content->colptrs = &(content->indexptrs);
    /* CSR indices */
    content->colvals = NULL;
    content->rowptrs = NULL;
    break;
  case SUN_CSR_MAT:
    content->NP      = M;
    content->colvals = &(content->indexvals);
    content->rowptrs = &(content->indexptrs);
    /* CSC indices */
    content->rowvals = NULL;
    content->colptrs = NULL;
  }
  content->data      = NULL;
  content->indexvals = NULL;
  content->indexptrs = NULL;

  /* Allocate content */
  content->data = (sunrealtype*)calloc(NNZ, sizeof(sunrealtype));
  SUNAssertNull(content->data, SUN_ERR_MALLOC_FAIL);

  content->indexvals = (sunindextype*)calloc(NNZ, sizeof(sunindextype));
  SUNAssertNull(content->indexvals, SUN_ERR_MALLOC_FAIL);

  content->indexptrs = (sunindextype*)calloc((content->NP + 1),
                                             sizeof(sunindextype));
  SUNAssertNull(content->indexptrs, SUN_ERR_MALLOC_FAIL);
  content->indexptrs[content->NP] = 0;

  return (A);
}

/* ----------------------------------------------------------------------------
 * Function to create a new sparse matrix from an existing dense matrix
 * by copying all nonzero values into the sparse matrix structure.  Returns NULL
 * if the request for matrix storage cannot be satisfied.
 */

SUNMatrix SUNSparseFromDenseMatrix(SUNMatrix Ad, sunrealtype droptol,
                                   int sparsetype)
{
  SUNFunctionBegin(Ad->sunctx);
  sunindextype i, j, nnz;
  sunindextype M, N;
  SUNMatrix As;

  SUNAssertNull(SUNMatGetID(Ad) == SUNMATRIX_DENSE, SUN_ERR_ARG_WRONGTYPE);
  SUNAssertNull(sparsetype == SUN_CSC_MAT || sparsetype == SUN_CSR_MAT,
                SUN_ERR_ARG_OUTOFRANGE);
  SUNAssertNull(droptol >= ZERO, SUN_ERR_ARG_OUTOFRANGE);

  /* set size of new matrix */
  M = SM_ROWS_D(Ad);
  N = SM_COLUMNS_D(Ad);

  /* determine total number of nonzeros */
  nnz = 0;
  for (j = 0; j < N; j++)
  {
    for (i = 0; i < M; i++)
    {
      nnz += (SUNRabs(SM_ELEMENT_D(Ad, i, j)) > droptol);
    }
  }

  /* allocate sparse matrix */
  As = NULL;
  As = SUNSparseMatrix(M, N, nnz, sparsetype, Ad->sunctx);
  SUNCheckLastErrNull();

  /* copy nonzeros from Ad into As, based on CSR/CSC type */
  nnz = 0;
  if (sparsetype == SUN_CSC_MAT)
  {
    for (j = 0; j < N; j++)
    {
      (SM_INDEXPTRS_S(As))[j] = nnz;
      for (i = 0; i < M; i++)
      {
        if (SUNRabs(SM_ELEMENT_D(Ad, i, j)) > droptol)
        {
          (SM_INDEXVALS_S(As))[nnz] = i;
          (SM_DATA_S(As))[nnz++]    = SM_ELEMENT_D(Ad, i, j);
        }
      }
    }
    (SM_INDEXPTRS_S(As))[N] = nnz;
  }
  else
  { /* SUN_CSR_MAT */
    for (i = 0; i < M; i++)
    {
      (SM_INDEXPTRS_S(As))[i] = nnz;
      for (j = 0; j < N; j++)
      {
        if (SUNRabs(SM_ELEMENT_D(Ad, i, j)) > droptol)
        {
          (SM_INDEXVALS_S(As))[nnz] = j;
          (SM_DATA_S(As))[nnz++]    = SM_ELEMENT_D(Ad, i, j);
        }
      }
    }
    (SM_INDEXPTRS_S(As))[M] = nnz;
  }

  return (As);
}

/* ----------------------------------------------------------------------------
 * Function to create a new sparse matrix from an existing band matrix
 * by copying all nonzero values into the sparse matrix structure.  Returns NULL
 * if the request for matrix storage cannot be satisfied.
 */

SUNMatrix SUNSparseFromBandMatrix(SUNMatrix Ab, sunrealtype droptol,
                                  int sparsetype)
{
  SUNFunctionBegin(Ab->sunctx);
  sunindextype i, j, nnz;
  sunindextype M, N;
  SUNMatrix As;

  SUNAssertNull(SUNMatGetID(Ab) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
  SUNAssertNull(sparsetype == SUN_CSC_MAT || sparsetype == SUN_CSR_MAT,
                SUN_ERR_ARG_OUTOFRANGE);
  SUNAssertNull(droptol >= ZERO, SUN_ERR_ARG_OUTOFRANGE);

  /* set size of new matrix */
  M = SM_ROWS_B(Ab);
  N = SM_COLUMNS_B(Ab);

  /* determine total number of nonzeros */
  nnz = 0;
  for (j = 0; j < N; j++)
  {
    for (i = SUNMAX(0, j - SM_UBAND_B(Ab));
         i <= SUNMIN(M - 1, j + SM_LBAND_B(Ab)); i++)
    {
      nnz += (SUNRabs(SM_ELEMENT_B(Ab, i, j)) > droptol);
    }
  }

  /* allocate sparse matrix */
  As = SUNSparseMatrix(M, N, nnz, sparsetype, Ab->sunctx);
  SUNCheckLastErrNull();

  /* copy nonzeros from Ab into As, based on CSR/CSC type */
  nnz = 0;
  if (sparsetype == SUN_CSC_MAT)
  {
    for (j = 0; j < N; j++)
    {
      (SM_INDEXPTRS_S(As))[j] = nnz;
      for (i = SUNMAX(0, j - SM_UBAND_B(Ab));
           i <= SUNMIN(M - 1, j + SM_LBAND_B(Ab)); i++)
      {
        if (SUNRabs(SM_ELEMENT_B(Ab, i, j)) > droptol)
        {
          (SM_INDEXVALS_S(As))[nnz] = i;
          (SM_DATA_S(As))[nnz++]    = SM_ELEMENT_B(Ab, i, j);
        }
      }
    }
    (SM_INDEXPTRS_S(As))[N] = nnz;
  }
  else
  { /* SUN_CSR_MAT */
    for (i = 0; i < M; i++)
    {
      (SM_INDEXPTRS_S(As))[i] = nnz;
      for (j = SUNMAX(0, i - SM_LBAND_B(Ab));
           j <= SUNMIN(N - 1, i + SM_UBAND_B(Ab)); j++)
      {
        if (SUNRabs(SM_ELEMENT_B(Ab, i, j)) > droptol)
        {
          (SM_INDEXVALS_S(As))[nnz] = j;
          (SM_DATA_S(As))[nnz++]    = SM_ELEMENT_B(Ab, i, j);
        }
      }
    }
    (SM_INDEXPTRS_S(As))[M] = nnz;
  }

  return (As);
}

/* ----------------------------------------------------------------------------
 * Function to create a new CSR matrix from a CSC matrix.
 */
SUNErrCode SUNSparseMatrix_ToCSR(const SUNMatrix A, SUNMatrix* Bout)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssert(SUNMatGetID(A) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);
  SUNAssert(SM_SPARSETYPE_S(A) == SUN_CSC_MAT, SUN_ERR_ARG_OUTOFRANGE);

  *Bout = SUNSparseMatrix(SM_ROWS_S(A), SM_COLUMNS_S(A), SM_NNZ_S(A),
                          SUN_CSR_MAT, A->sunctx);
  SUNCheckLastErr();
  SUNCheckCall(format_convert(A, *Bout));
  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Function to create a new CSC matrix from a CSR matrix.
 */
SUNErrCode SUNSparseMatrix_ToCSC(const SUNMatrix A, SUNMatrix* Bout)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssert(SUNMatGetID(A) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);
  SUNAssert(SM_SPARSETYPE_S(A) == SUN_CSR_MAT, SUN_ERR_ARG_OUTOFRANGE);

  *Bout = SUNSparseMatrix(SM_ROWS_S(A), SM_COLUMNS_S(A), SM_NNZ_S(A),
                          SUN_CSC_MAT, A->sunctx);
  SUNCheckLastErr();

  SUNCheckCall(format_convert(A, *Bout));
  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Function to reallocate internal sparse matrix storage arrays so that the
 * resulting sparse matrix holds indexptrs[NP] nonzeros.  Returns 0 on success
 * and 1 on failure (e.g. if A does not have sparse type, or if nnz is negative)
 */

SUNErrCode SUNSparseMatrix_Realloc(SUNMatrix A)
{
  sunindextype nzmax;
  SUNFunctionBegin(A->sunctx);

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);

  /* get total number of nonzeros (return with failure if illegal) */
  nzmax = (SM_INDEXPTRS_S(A))[SM_NP_S(A)];
  SUNAssert(nzmax >= 0, SUN_ERR_ARG_CORRUPT);

  /* perform reallocation */
  SM_INDEXVALS_S(A) = (sunindextype*)realloc(SM_INDEXVALS_S(A),
                                             nzmax * sizeof(sunindextype));
  SUNAssert(SM_INDEXVALS_S(A), SUN_ERR_MALLOC_FAIL);

  SM_DATA_S(A) = (sunrealtype*)realloc(SM_DATA_S(A), nzmax * sizeof(sunrealtype));
  SUNAssert(SM_DATA_S(A), SUN_ERR_MALLOC_FAIL);

  SM_NNZ_S(A) = nzmax;

  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Function to reallocate internal sparse matrix storage arrays so that the
 * resulting sparse matrix has storage for a specified number of nonzeros.
 * Returns 0 on success and 1 on failure (e.g. if A does not have sparse type,
 * or if nnz is negative)
 */

SUNErrCode SUNSparseMatrix_Reallocate(SUNMatrix A, sunindextype NNZ)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssert(SUNMatGetID(A) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);
  SUNAssert(NNZ >= 0, SUN_ERR_ARG_OUTOFRANGE);

  /* perform reallocation */
  SM_INDEXVALS_S(A) = (sunindextype*)realloc(SM_INDEXVALS_S(A),
                                             NNZ * sizeof(sunindextype));
  SUNAssert(SM_INDEXVALS_S(A), SUN_ERR_MALLOC_FAIL);

  SM_DATA_S(A) = (sunrealtype*)realloc(SM_DATA_S(A), NNZ * sizeof(sunrealtype));
  SUNAssert(SM_DATA_S(A), SUN_ERR_MALLOC_FAIL);

  SM_NNZ_S(A) = NNZ;

  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Function to print the sparse matrix
 */

void SUNSparseMatrix_Print(SUNMatrix A, FILE* outfile)
{
  SUNFunctionBegin(A->sunctx);
  sunindextype i, j;
  const char* matrixtype;
  const char* indexname;

  SUNAssertVoid(SUNMatGetID(A) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);

  /* perform operation */
  if (SM_SPARSETYPE_S(A) == SUN_CSC_MAT)
  {
    indexname  = "col";
    matrixtype = "CSC";
  }
  else
  {
    indexname  = "row";
    matrixtype = "CSR";
  }
  fprintf(outfile, "\n");
  fprintf(outfile, "%ld by %ld %s matrix, NNZ: %ld \n", (long int)SM_ROWS_S(A),
          (long int)SM_COLUMNS_S(A), matrixtype, (long int)SM_NNZ_S(A));
  for (j = 0; j < SM_NP_S(A); j++)
  {
    fprintf(outfile, "%s %ld : locations %ld to %ld\n", indexname, (long int)j,
            (long int)(SM_INDEXPTRS_S(A))[j],
            (long int)(SM_INDEXPTRS_S(A))[j + 1] - 1);
    fprintf(outfile, "  ");
    for (i = (SM_INDEXPTRS_S(A))[j]; i < (SM_INDEXPTRS_S(A))[j + 1]; i++)
    {
      fprintf(outfile, "%ld: " SUN_FORMAT_E "  ",
              (long int)(SM_INDEXVALS_S(A))[i], (SM_DATA_S(A))[i]);
    }
    fprintf(outfile, "\n");
  }
  return;
}

/* ----------------------------------------------------------------------------
 * Functions to access the contents of the sparse matrix structure
 */

sunindextype SUNSparseMatrix_Rows(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNoRet(SUNMatGetID(A) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);
  return SM_ROWS_S(A);
}

sunindextype SUNSparseMatrix_Columns(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNoRet(SUNMatGetID(A) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);
  return SM_COLUMNS_S(A);
}

sunindextype SUNSparseMatrix_NNZ(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNoRet(SUNMatGetID(A) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);
  return SM_NNZ_S(A);
}

sunindextype SUNSparseMatrix_NP(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNoRet(SUNMatGetID(A) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);
  return SM_NP_S(A);
}

int SUNSparseMatrix_SparseType(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNoRet(SUNMatGetID(A) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);
  return SM_SPARSETYPE_S(A);
}

sunrealtype* SUNSparseMatrix_Data(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNull(SUNMatGetID(A) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);
  return SM_DATA_S(A);
}

sunindextype* SUNSparseMatrix_IndexValues(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNull(SUNMatGetID(A) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);
  return SM_INDEXVALS_S(A);
}

sunindextype* SUNSparseMatrix_IndexPointers(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssertNull(SUNMatGetID(A) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);
  return SM_INDEXPTRS_S(A);
}

/*
 * -----------------------------------------------------------------
 * implementation of matrix operations
 * -----------------------------------------------------------------
 */

SUNMatrix_ID SUNMatGetID_Sparse(SUNDIALS_MAYBE_UNUSED SUNMatrix A)
{
  return SUNMATRIX_SPARSE;
}

SUNMatrix SUNMatClone_Sparse(SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  SUNMatrix B = SUNSparseMatrix(SM_ROWS_S(A), SM_COLUMNS_S(A), SM_NNZ_S(A),
                                SM_SPARSETYPE_S(A), A->sunctx);
  SUNCheckLastErrNull();
  return (B);
}

void SUNMatDestroy_Sparse(SUNMatrix A)
{
  if (A == NULL) { return; }

  /* free content */
  if (A->content != NULL)
  {
    /* free data array */
    if (SM_DATA_S(A))
    {
      free(SM_DATA_S(A));
      SM_DATA_S(A) = NULL;
    }
    /* free index values array */
    if (SM_INDEXVALS_S(A))
    {
      free(SM_INDEXVALS_S(A));
      SM_INDEXVALS_S(A)        = NULL;
      SM_CONTENT_S(A)->rowvals = NULL;
      SM_CONTENT_S(A)->colvals = NULL;
    }
    /* free index pointers array */
    if (SM_INDEXPTRS_S(A))
    {
      free(SM_INDEXPTRS_S(A));
      SM_INDEXPTRS_S(A)        = NULL;
      SM_CONTENT_S(A)->colptrs = NULL;
      SM_CONTENT_S(A)->rowptrs = NULL;
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

SUNErrCode SUNMatZero_Sparse(SUNMatrix A)
{
  sunindextype i;

  /* Perform operation */
  for (i = 0; i < SM_NNZ_S(A); i++)
  {
    (SM_DATA_S(A))[i]      = ZERO;
    (SM_INDEXVALS_S(A))[i] = 0;
  }
  for (i = 0; i < SM_NP_S(A); i++) { (SM_INDEXPTRS_S(A))[i] = 0; }
  (SM_INDEXPTRS_S(A))[SM_NP_S(A)] = 0;
  return SUN_SUCCESS;
}

SUNErrCode SUNMatCopy_Sparse(SUNMatrix A, SUNMatrix B)
{
  sunindextype i, A_nz;
  SUNFunctionBegin(A->sunctx);

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);
  SUNAssert(SUNMatGetID(B) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);
  SUNCheck(compatibleMatrices(A, B), SUN_ERR_ARG_DIMSMISMATCH);

  /* Perform operation */
  A_nz = (SM_INDEXPTRS_S(A))[SM_NP_S(A)];

  /* ensure that B is allocated with at least as
     much memory as we have nonzeros in A */
  if (SM_NNZ_S(B) < A_nz)
  {
    SM_INDEXVALS_S(B) = (sunindextype*)realloc(SM_INDEXVALS_S(B),
                                               A_nz * sizeof(sunindextype));
    SUNAssert(SM_INDEXVALS_S(B), SUN_ERR_MALLOC_FAIL);

    SM_DATA_S(B) = (sunrealtype*)realloc(SM_DATA_S(B),
                                         A_nz * sizeof(sunrealtype));
    SUNAssert(SM_DATA_S(B), SUN_ERR_MALLOC_FAIL);

    SM_NNZ_S(B) = A_nz;
  }

  /* zero out B so that copy works correctly */
  SUNCheckCall(SUNMatZero_Sparse(B));

  /* copy the data and row indices over */
  for (i = 0; i < A_nz; i++)
  {
    (SM_DATA_S(B))[i]      = (SM_DATA_S(A))[i];
    (SM_INDEXVALS_S(B))[i] = (SM_INDEXVALS_S(A))[i];
  }

  /* copy the column pointers over */
  for (i = 0; i < SM_NP_S(A); i++)
  {
    (SM_INDEXPTRS_S(B))[i] = (SM_INDEXPTRS_S(A))[i];
  }
  (SM_INDEXPTRS_S(B))[SM_NP_S(A)] = A_nz;

  return SUN_SUCCESS;
}

SUNErrCode SUNMatScaleAddI_Sparse(sunrealtype c, SUNMatrix A)
{
  SUNFunctionBegin(A->sunctx);
  const sunindextype N = SM_SPARSETYPE_S(A) == SUN_CSC_MAT ? SM_COLUMNS_S(A)
                                                           : SM_ROWS_S(A);
  const sunindextype M = SM_SPARSETYPE_S(A) == SUN_CSC_MAT ? SM_ROWS_S(A)
                                                           : SM_COLUMNS_S(A);

  sunindextype* Ap = SM_INDEXPTRS_S(A);
  SUNAssert(Ap, SUN_ERR_ARG_CORRUPT);
  sunindextype* Ai = SM_INDEXVALS_S(A);
  SUNAssert(Ai, SUN_ERR_ARG_CORRUPT);
  sunrealtype* Ax = SM_DATA_S(A);
  SUNAssert(Ax, SUN_ERR_ARG_CORRUPT);

  sunindextype newvals = 0;
  for (sunindextype j = 0; j < N; j++)
  {
    /* scan column (row if CSR) of A, searching for diagonal value */
    sunbooleantype found = SUNFALSE;
    for (sunindextype i = Ap[j]; i < Ap[j + 1]; i++)
    {
      if (Ai[i] == j)
      {
        found = SUNTRUE;
        Ax[i] = ONE + c * Ax[i];
      }
      else { Ax[i] *= c; }
    }
    /* If no diagonal element found and the current column (row) can actually
     * contain a diagonal element, increment the storage counter */
    if (!found && j < M) { newvals++; }
  }

  /* At this point, A has the correctly updated values except for any new
   * diagonal elements that need to be added (of which there are newvals). Now,
   * we allocate additional storage if needed */
  const sunindextype new_nnz = Ap[N] + newvals;
  if (new_nnz > SM_NNZ_S(A))
  {
    SUNCheckCall(SUNSparseMatrix_Reallocate(A, new_nnz));
    Ap = SM_INDEXPTRS_S(A);
    Ai = SM_INDEXVALS_S(A);
    Ax = SM_DATA_S(A);
  }

  for (sunindextype j = N - 1; newvals > 0; j--)
  {
    sunbooleantype found = SUNFALSE;
    for (sunindextype i = Ap[j + 1] - 1; i >= Ap[j]; i--)
    {
      if (Ai[i] == j) { found = SUNTRUE; }

      /* Shift elements to make room for diagonal elements */
      Ai[i + newvals] = Ai[i];
      Ax[i + newvals] = Ax[i];
    }

    Ap[j + 1] += newvals;
    if (!found && j < M)
    {
      /* This column (row) needs a diagonal element added */
      newvals--;
      Ai[Ap[j] + newvals] = j;
      Ax[Ap[j] + newvals] = ONE;
    }
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNMatScaleAdd_Sparse(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  sunindextype j, i, p, nz, newvals, M, N, cend;
  sunbooleantype newmat;
  sunindextype *w, *Ap, *Ai, *Bp, *Bi, *Cp, *Ci;
  sunrealtype *x, *Ax, *Bx, *Cx;
  SUNMatrix C;
  SUNFunctionBegin(A->sunctx);

  SUNAssert(SUNMatGetID(A) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);
  SUNAssert(SUNMatGetID(B) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);
  SUNCheck(compatibleMatrices(A, B), SUN_ERR_ARG_DIMSMISMATCH);

  /* store shortcuts to matrix dimensions (M is inner dimension, N is outer) */
  if (SM_SPARSETYPE_S(A) == SUN_CSC_MAT)
  {
    M = SM_ROWS_S(A);
    N = SM_COLUMNS_S(A);
  }
  else
  {
    M = SM_COLUMNS_S(A);
    N = SM_ROWS_S(A);
  }

  /* access data arrays from A and B (return if failure) */
  Ap = Ai = Bp = Bi = NULL;
  Ax = Bx = NULL;
  Ap      = SM_INDEXPTRS_S(A);
  SUNAssert(Ap, SUN_ERR_ARG_CORRUPT);
  Ai = SM_INDEXVALS_S(A);
  SUNAssert(Ai, SUN_ERR_ARG_CORRUPT);
  Ax = SM_DATA_S(A);
  SUNAssert(Ax, SUN_ERR_ARG_CORRUPT);
  Bp = SM_INDEXPTRS_S(B);
  SUNAssert(Bp, SUN_ERR_ARG_CORRUPT);
  Bi = SM_INDEXVALS_S(B);
  SUNAssert(Bi, SUN_ERR_ARG_CORRUPT);
  Bx = SM_DATA_S(B);
  SUNAssert(Bx, SUN_ERR_ARG_CORRUPT);

  /* create work arrays for row indices and nonzero column values */
  w = (sunindextype*)malloc(M * sizeof(sunindextype));
  SUNAssert(w, SUN_ERR_MALLOC_FAIL);
  x = (sunrealtype*)malloc(M * sizeof(sunrealtype));
  SUNAssert(x, SUN_ERR_MALLOC_FAIL);

  /* determine if A already contains the sparsity pattern of B */
  newvals = 0;
  for (j = 0; j < N; j++)
  {
    /* clear work array */
    for (i = 0; i < M; i++) { w[i] = 0; }

    /* scan column of A, incrementing w by one */
    for (i = Ap[j]; i < Ap[j + 1]; i++) { w[Ai[i]] += 1; }

    /* scan column of B, decrementing w by one */
    for (i = Bp[j]; i < Bp[j + 1]; i++) { w[Bi[i]] -= 1; }

    /* if any entry of w is negative, A doesn't contain B's sparsity,
       so increment necessary storage counter */
    for (i = 0; i < M; i++)
    {
      if (w[i] < 0) { newvals += 1; }
    }
  }

  /* If extra nonzeros required, check whether A has sufficient storage space
     for new nonzero entries (so B can be inserted into existing storage) */
  newmat = SUNFALSE; /* no reallocation needed */
  if (newvals > (SM_NNZ_S(A) - Ap[N])) { newmat = SUNTRUE; }

  /* perform operation based on existing/necessary structure */

  /*   case 1: A already contains sparsity pattern of B */
  if (newvals == 0)
  {
    /* iterate through columns, adding matrices */
    for (j = 0; j < N; j++)
    {
      /* clear work array */
      for (i = 0; i < M; i++) { x[i] = ZERO; }

      /* scan column of B, updating work array */
      for (i = Bp[j]; i < Bp[j + 1]; i++) { x[Bi[i]] = Bx[i]; }

      /* scan column of A, updating array entries appropriately */
      for (i = Ap[j]; i < Ap[j + 1]; i++) { Ax[i] = c * Ax[i] + x[Ai[i]]; }
    }

    /*   case 2: A has sufficient storage, but does not already contain B's
     * sparsity */
  }
  else if (!newmat)
  {
    /* determine storage location where last column (row) should end */
    nz = Ap[N] + newvals;

    /* store pointer past last column (row) from original A,
       and store updated value in revised A */
    cend  = Ap[N];
    Ap[N] = nz;

    /* iterate through columns (rows) backwards */
    for (j = N - 1; j >= 0; j--)
    {
      /* clear out temporary arrays for this column (row) */
      for (i = 0; i < M; i++)
      {
        w[i] = 0;
        x[i] = SUN_RCONST(0.0);
      }

      /* iterate down column (row) of A, collecting nonzeros */
      for (p = Ap[j]; p < cend; p++)
      {
        w[Ai[p]] += 1;        /* indicate that row (column) is filled */
        x[Ai[p]] = c * Ax[p]; /* collect/scale value */
      }

      /* iterate down column of B, collecting nonzeros */
      for (p = Bp[j]; p < Bp[j + 1]; p++)
      {
        w[Bi[p]] += 1;     /* indicate that row is filled */
        x[Bi[p]] += Bx[p]; /* collect value */
      }

      /* fill entries of A with this column's (row's) data */
      for (i = M - 1; i >= 0; i--)
      {
        if (w[i] > 0)
        {
          Ai[--nz] = i;
          Ax[nz]   = x[i];
        }
      }

      /* store ptr past this col (row) from orig A, update value for new A */
      cend  = Ap[j];
      Ap[j] = nz;
    }

    /*   case 3: A must be reallocated with sufficient storage */
  }
  else
  {
    /* create new matrix for sum */
    C = SUNSparseMatrix(SM_ROWS_S(A), SM_COLUMNS_S(A), Ap[N] + newvals,
                        SM_SPARSETYPE_S(A), A->sunctx);
    SUNCheckLastErr();

    /* access data from CSR structures (return if failure) */
    Cp = NULL;
    Ci = NULL;
    Cx = NULL;
    Cp = SM_INDEXPTRS_S(C);
    SUNAssert(Cp, SUN_ERR_ARG_CORRUPT);
    Ci = SM_INDEXVALS_S(C);
    SUNAssert(Ci, SUN_ERR_ARG_CORRUPT);
    Cx = SM_DATA_S(C);
    SUNAssert(Cx, SUN_ERR_ARG_CORRUPT);

    /* initialize total nonzero count */
    nz = 0;

    /* iterate through columns (rows) */
    for (j = 0; j < N; j++)
    {
      /* set current column (row) pointer to current # nonzeros */
      Cp[j] = nz;

      /* clear out temporary arrays for this column (row) */
      for (i = 0; i < M; i++)
      {
        w[i] = 0;
        x[i] = SUN_RCONST(0.0);
      }

      /* iterate down column of A, collecting nonzeros */
      for (p = Ap[j]; p < Ap[j + 1]; p++)
      {
        w[Ai[p]] += 1;        /* indicate that row is filled */
        x[Ai[p]] = c * Ax[p]; /* collect/scale value */
      }

      /* iterate down column of B, collecting nonzeros */
      for (p = Bp[j]; p < Bp[j + 1]; p++)
      {
        w[Bi[p]] += 1;     /* indicate that row is filled */
        x[Bi[p]] += Bx[p]; /* collect value */
      }

      /* fill entries of C with this column's data */
      for (i = 0; i < M; i++)
      {
        if (w[i] > 0)
        {
          Ci[nz]   = i;
          Cx[nz++] = x[i];
        }
      }
    }

    /* indicate end of data */
    Cp[N] = nz;

    /* update A's structure with C's values; nullify C's pointers */
    SM_NNZ_S(A) = SM_NNZ_S(C);

    free(SM_DATA_S(A));
    SM_DATA_S(A) = SM_DATA_S(C);
    SM_DATA_S(C) = NULL;

    free(SM_INDEXVALS_S(A));
    SM_INDEXVALS_S(A) = SM_INDEXVALS_S(C);
    SM_INDEXVALS_S(C) = NULL;

    free(SM_INDEXPTRS_S(A));
    SM_INDEXPTRS_S(A) = SM_INDEXPTRS_S(C);
    SM_INDEXPTRS_S(C) = NULL;

    /* clean up */
    SUNMatDestroy_Sparse(C);
  }

  /* clean up */
  free(w);
  free(x);

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode SUNMatMatvec_Sparse(SUNMatrix A, N_Vector x, N_Vector y)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssert(SUNMatGetID(A) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);
  SUNCheck(compatibleMatrixAndVectors(A, x, y), SUN_ERR_ARG_DIMSMISMATCH);

  /* Perform operation */
  if (SM_SPARSETYPE_S(A) == SUN_CSC_MAT)
  {
    SUNCheckCall(Matvec_SparseCSC(A, x, y));
  }
  else { SUNCheckCall(Matvec_SparseCSR(A, x, y)); }

  return SUN_SUCCESS;
}

SUNErrCode SUNMatHermitianTransposeVec_Sparse(SUNMatrix A, N_Vector x, N_Vector y)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssert(SUNMatGetID(A) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);
  SUNCheck(compatibleMatrixAndVectors(A, y, x), SUN_ERR_ARG_DIMSMISMATCH);

  /* Perform operation */
  if (SM_SPARSETYPE_S(A) == SUN_CSC_MAT)
  {
    SUNCheckCall(MatTransposeVec_SparseCSC(A, x, y));
  }
  else { SUNCheckCall(MatTransposeVec_SparseCSR(A, x, y)); }

  return SUN_SUCCESS;
}

SUNErrCode SUNMatSpace_Sparse(SUNMatrix A, long int* lenrw, long int* leniw)
{
  SUNFunctionBegin(A->sunctx);
  SUNAssert(SUNMatGetID(A) == SUNMATRIX_SPARSE, SUN_ERR_ARG_WRONGTYPE);
  SUNAssert(lenrw, SUN_ERR_ARG_CORRUPT);
  SUNAssert(leniw, SUN_ERR_ARG_CORRUPT);
  *lenrw = SM_NNZ_S(A);
  *leniw = 10 + SM_NP_S(A) + SM_NNZ_S(A);
  return SUN_SUCCESS;
}

/*
 * =================================================================
 * private functions
 * =================================================================
 */

/* -----------------------------------------------------------------
 * Function to check compatibility of two sparse SUNMatrix objects
 */

SUNDIALS_MAYBE_UNUSED
static sunbooleantype compatibleMatrices(SUNMatrix A, SUNMatrix B)
{
  /* both matrices must have the same shape and sparsity type */
  if (SUNSparseMatrix_Rows(A) != SUNSparseMatrix_Rows(B)) { return SUNFALSE; }
  if (SUNSparseMatrix_Columns(A) != SUNSparseMatrix_Columns(B))
  {
    return SUNFALSE;
  }
  if (SM_SPARSETYPE_S(A) != SM_SPARSETYPE_S(B)) { return SUNFALSE; }

  return SUNTRUE;
}

/* -----------------------------------------------------------------
 * Function to check compatibility of a SUNMatrix object with two
 * N_Vectors (A*x = b)
 */

SUNDIALS_MAYBE_UNUSED
static sunbooleantype compatibleMatrixAndVectors(SUNMatrix A, N_Vector x,
                                                 N_Vector y)
{
  /* vectors must implement N_VGetArrayPointer */
  if ((x->ops->nvgetarraypointer == NULL) || (y->ops->nvgetarraypointer == NULL))
  {
    return SUNFALSE;
  }

  /* Verify that the dimensions of A, x, and y agree */
  if ((SUNSparseMatrix_Columns(A) != N_VGetLength(x)) ||
      (SUNSparseMatrix_Rows(A) != N_VGetLength(y)))
  {
    return SUNFALSE;
  }

  return SUNTRUE;
}

/* -----------------------------------------------------------------
 * Computes y=A*x, where A is a CSC SUNMatrix_Sparse of dimension MxN, x is a
 * compatible N_Vector object of length N, and y is a compatible
 * N_Vector object of length M.
 *
 * Returns 0 if successful, 1 if unsuccessful (failed memory access, or both
 * x and y are the same vector).
 */
SUNErrCode Matvec_SparseCSC(SUNMatrix A, N_Vector x, N_Vector y)
{
  sunindextype i, j;
  sunindextype *Ap, *Ai;
  sunrealtype *Ax, *xd, *yd;
  SUNFunctionBegin(A->sunctx);

  /* access data from CSC structure (return if failure) */
  Ap = SM_INDEXPTRS_S(A);
  SUNAssert(Ap, SUN_ERR_ARG_CORRUPT);
  Ai = SM_INDEXVALS_S(A);
  SUNAssert(Ai, SUN_ERR_ARG_CORRUPT);
  Ax = SM_DATA_S(A);
  SUNAssert(Ax, SUN_ERR_ARG_CORRUPT);

  /* access vector data (return if failure) */
  xd = N_VGetArrayPointer(x);
  SUNCheckLastErr();
  yd = N_VGetArrayPointer(y);
  SUNCheckLastErr();

  /* initialize result */
  for (i = 0; i < SM_ROWS_S(A); i++) { yd[i] = ZERO; }

  /* iterate through matrix columns */
  for (j = 0; j < SM_COLUMNS_S(A); j++)
  {
    /* iterate down column of A, performing product */
    for (i = Ap[j]; i < Ap[j + 1]; i++) { yd[Ai[i]] += Ax[i] * xd[j]; }
  }

  return SUN_SUCCESS;
}

SUNErrCode MatTransposeVec_SparseCSC(SUNMatrix A, N_Vector x, N_Vector y)
{
  sunindextype i, j;
  sunindextype *Ap, *Ai;
  sunrealtype *Ax, *xd, *yd;
  SUNFunctionBegin(A->sunctx);

  /* access data from CSC structure (return if failure) */
  Ap = SM_INDEXPTRS_S(A);
  SUNAssert(Ap, SUN_ERR_ARG_CORRUPT);
  Ai = SM_INDEXVALS_S(A);
  SUNAssert(Ai, SUN_ERR_ARG_CORRUPT);
  Ax = SM_DATA_S(A);
  SUNAssert(Ax, SUN_ERR_ARG_CORRUPT);

  /* access vector data (return if failure) */
  xd = N_VGetArrayPointer(x);
  SUNCheckLastErr();
  yd = N_VGetArrayPointer(y);
  SUNCheckLastErr();

  /* initialize result vector */
  for (i = 0; i < SM_COLUMNS_S(A); i++) { yd[i] = ZERO; }

  /* iterate through matrix columns (rows of the transposed matrix) */
  for (j = 0; j < SM_COLUMNS_S(A); j++)
  {
    /* iterate through non-zero elements in the current column */
    for (i = Ap[j]; i < Ap[j + 1]; i++) { yd[j] += Ax[i] * xd[Ai[i]]; }
  }

  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * Computes y=A*x, where A is a CSR SUNMatrix_Sparse of dimension MxN, x is a
 * compatible N_Vector object of length N, and y is a compatible
 * N_Vector object of length M.
 *
 * Returns 0 if successful, -1 if unsuccessful (failed memory access).
 */
SUNErrCode Matvec_SparseCSR(SUNMatrix A, N_Vector x, N_Vector y)
{
  sunindextype i, j;
  sunindextype *Ap, *Aj;
  sunrealtype *Ax, *xd, *yd;
  SUNFunctionBegin(A->sunctx);

  /* access data from CSR structure (return if failure) */
  Ap = SM_INDEXPTRS_S(A);
  SUNAssert(Ap, SUN_ERR_ARG_CORRUPT);
  Aj = SM_INDEXVALS_S(A);
  SUNAssert(Aj, SUN_ERR_ARG_CORRUPT);
  Ax = SM_DATA_S(A);
  SUNAssert(Ax, SUN_ERR_ARG_CORRUPT);

  /* access vector data (return if failure) */
  xd = N_VGetArrayPointer(x);
  SUNCheckLastErr();
  yd = N_VGetArrayPointer(y);
  SUNCheckLastErr();
  SUNAssert(xd, SUN_ERR_ARG_CORRUPT);
  SUNAssert(yd, SUN_ERR_ARG_CORRUPT);
  SUNAssert(xd != yd, SUN_ERR_ARG_CORRUPT);

  /* initialize result */
  for (i = 0; i < SM_ROWS_S(A); i++) { yd[i] = ZERO; }

  /* iterate through matrix rows */
  for (i = 0; i < SM_ROWS_S(A); i++)
  {
    /* iterate along row of A, performing product */
    for (j = Ap[i]; j < Ap[i + 1]; j++) { yd[i] += Ax[j] * xd[Aj[j]]; }
  }

  return SUN_SUCCESS;
}

SUNErrCode MatTransposeVec_SparseCSR(SUNMatrix A, N_Vector x, N_Vector y)
{
  sunindextype i, j;
  sunindextype *Ap, *Aj;
  sunrealtype *Ax, *xd, *yd;
  SUNFunctionBegin(A->sunctx);

  /* access data from CSR structure (return if failure) */
  Ap = SM_INDEXPTRS_S(A);
  SUNAssert(Ap, SUN_ERR_ARG_CORRUPT);
  Aj = SM_INDEXVALS_S(A);
  SUNAssert(Aj, SUN_ERR_ARG_CORRUPT);
  Ax = SM_DATA_S(A);
  SUNAssert(Ax, SUN_ERR_ARG_CORRUPT);

  /* access vector data (return if failure) */
  xd = N_VGetArrayPointer(x);
  SUNCheckLastErr();
  yd = N_VGetArrayPointer(y);
  SUNCheckLastErr();
  SUNAssert(xd, SUN_ERR_ARG_CORRUPT);
  SUNAssert(yd, SUN_ERR_ARG_CORRUPT);
  SUNAssert(xd != yd, SUN_ERR_ARG_CORRUPT);

  /* initialize result vector */
  for (i = 0; i < SM_COLUMNS_S(A); i++) { yd[i] = ZERO; }

  /* iterate over rows of the original matrix (columns of the transposed matrix) */
  for (i = 0; i < SM_ROWS_S(A); i++)
  {
    for (j = Ap[i]; j < Ap[i + 1]; j++) { yd[Aj[j]] += Ax[j] * xd[i]; }
  }

  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * Copies A into a matrix B in the opposite format of A.
 * Returns 0 if successful, nonzero if unsuccessful.
 */
SUNErrCode format_convert(const SUNMatrix A, SUNMatrix B)
{
  SUNFunctionBegin(A->sunctx);

  sunrealtype *Ax, *Bx;
  sunindextype *Ap, *Aj;
  sunindextype *Bp, *Bi;
  sunindextype n_row, n_col, nnz;
  sunindextype n, col, csum, row, last;

  if (SM_SPARSETYPE_S(A) == SM_SPARSETYPE_S(B))
  {
    SUNCheckCall(SUNMatCopy_Sparse(A, B));
    return SUN_SUCCESS;
  }

  Ap = SM_INDEXPTRS_S(A);
  SUNAssert(Ap, SUN_ERR_ARG_CORRUPT);
  Aj = SM_INDEXVALS_S(A);
  SUNAssert(Aj, SUN_ERR_ARG_CORRUPT);
  Ax = SM_DATA_S(A);
  SUNAssert(Ax, SUN_ERR_ARG_CORRUPT);

  n_row = (SM_SPARSETYPE_S(A) == SUN_CSR_MAT) ? SM_ROWS_S(A) : SM_COLUMNS_S(A);
  n_col = (SM_SPARSETYPE_S(A) == SUN_CSR_MAT) ? SM_COLUMNS_S(A) : SM_ROWS_S(A);

  Bp = SM_INDEXPTRS_S(B);
  Bi = SM_INDEXVALS_S(B);
  Bx = SM_DATA_S(B);

  nnz = Ap[n_row];

  SUNCheckCall(SUNMatZero_Sparse(B));

  /* compute number of non-zero entries per column (if CSR) or per row (if CSC)
   * of A */
  for (n = 0; n < nnz; n++) { Bp[Aj[n]]++; }

  /* cumualtive sum the nnz per column to get Bp[] */
  for (col = 0, csum = 0; col < n_col; col++)
  {
    sunindextype temp = Bp[col];
    Bp[col]           = csum;
    csum += temp;
  }
  Bp[n_col] = nnz;

  for (row = 0; row < n_row; row++)
  {
    sunindextype jj;
    for (jj = Ap[row]; jj < Ap[row + 1]; jj++)
    {
      col               = Aj[jj];
      sunindextype dest = Bp[col];

      Bi[dest] = row;
      Bx[dest] = Ax[jj];

      Bp[col]++;
    }
  }

  for (col = 0, last = 0; col <= n_col; col++)
  {
    sunindextype temp = Bp[col];
    Bp[col]           = last;
    last              = temp;
  }

  return SUN_SUCCESS;
}
