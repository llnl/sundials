/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ UMBC
 *                David Gardner @ LLNL
 * Based on code sundials_sparse.h by: Carol Woodward and
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
 * This is the header file for the sparse implementation of the
 * SUNMATRIX module, SUNMATRIX_SPARSE.
 *
 * Notes:
 *   - The definition of the generic SUNMatrix structure can be found
 *     in the header file sundials_matrix.h.
 *   - The definition of the type 'sunrealtype' can be found in the
 *     header file sundials_types.h, and it may be changed (at the
 *     configuration stage) according to the user's needs.
 *     The sundials_types.h file also contains the definition
 *     for the type 'sunbooleantype' and 'indextype'.
 * -----------------------------------------------------------------
 */

#ifndef _SUNMATRIX_SPARSE_H
#define _SUNMATRIX_SPARSE_H

#include <stdio.h>
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_dense.h>

#include "sundials/sundials_types.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ------------------------
 * Matrix Type Definitions
 * ------------------------ */

#define SUN_CSC_MAT 0
#define SUN_CSR_MAT 1

/* DEPRECATION NOTICE: CSC_MAT and CSR_MAT will be removed in the next major release.
    use SUN_CSC_MAT and SUN_CSR_MAT instead. */
#define CSC_MAT 0
#define CSR_MAT 1

/* ------------------------------------------
 * Sparse Implementation of SUNMATRIX_SPARSE
 * ------------------------------------------ */

struct _SUNMatrixContent_Sparse
{
  sunindextype M;
  sunindextype N;
  sunindextype NNZ;
  sunindextype NP;
  sunrealtype* data;
  int sparsetype;
  sunindextype* indexvals;
  sunindextype* indexptrs;
  /* CSC indices */
  sunindextype** rowvals;
  sunindextype** colptrs;
  /* CSR indices */
  sunindextype** colvals;
  sunindextype** rowptrs;
};

typedef struct _SUNMatrixContent_Sparse* SUNMatrixContent_Sparse;

/* ---------------------------------------
 * Macros for access to SUNMATRIX_SPARSE
 * --------------------------------------- */

#define SM_CONTENT_S(A) ((SUNMatrixContent_Sparse)(A->content))

#define SM_ROWS_S(A) (SM_CONTENT_S(A)->M)

#define SM_COLUMNS_S(A) (SM_CONTENT_S(A)->N)

#define SM_NNZ_S(A) (SM_CONTENT_S(A)->NNZ)

#define SM_NP_S(A) (SM_CONTENT_S(A)->NP)

#define SM_SPARSETYPE_S(A) (SM_CONTENT_S(A)->sparsetype)

#define SM_DATA_S(A) (SM_CONTENT_S(A)->data)

#define SM_INDEXVALS_S(A) (SM_CONTENT_S(A)->indexvals)

#define SM_INDEXPTRS_S(A) (SM_CONTENT_S(A)->indexptrs)

/* ----------------------------------------
 * Exported Functions for SUNMATRIX_SPARSE
 * ---------------------------------------- */

SUNDIALS_EXPORT
SUNMatrix SUNSparseMatrix(sunindextype M, sunindextype N, sunindextype NNZ,
                          int sparsetype, SUNContext sunctx);

SUNDIALS_EXPORT
SUNMatrix SUNSparseFromDenseMatrix(SUNMatrix A, sunrealtype droptol,
                                   int sparsetype);

SUNDIALS_EXPORT
SUNMatrix SUNSparseFromBandMatrix(SUNMatrix A, sunrealtype droptol,
                                  int sparsetype);

SUNDIALS_EXPORT
SUNErrCode SUNSparseMatrix_ToCSR(const SUNMatrix A, SUNMatrix* Bout);

SUNDIALS_EXPORT
SUNErrCode SUNSparseMatrix_ToCSC(const SUNMatrix A, SUNMatrix* Bout);

SUNDIALS_EXPORT
SUNErrCode SUNSparseMatrix_Realloc(SUNMatrix A);

SUNDIALS_EXPORT
SUNErrCode SUNSparseMatrix_Reallocate(SUNMatrix A, sunindextype NNZ);

SUNDIALS_EXPORT
void SUNSparseMatrix_Print(SUNMatrix A, FILE* outfile);

SUNDIALS_EXPORT
sunindextype SUNSparseMatrix_Rows(SUNMatrix A);

SUNDIALS_EXPORT
sunindextype SUNSparseMatrix_Columns(SUNMatrix A);

SUNDIALS_EXPORT
sunindextype SUNSparseMatrix_NNZ(SUNMatrix A);

SUNDIALS_EXPORT
sunindextype SUNSparseMatrix_NP(SUNMatrix A);

SUNDIALS_EXPORT
int SUNSparseMatrix_SparseType(SUNMatrix A);

SUNDIALS_EXPORT
sunrealtype* SUNSparseMatrix_Data(SUNMatrix A);

SUNDIALS_EXPORT
sunindextype* SUNSparseMatrix_IndexValues(SUNMatrix A);

SUNDIALS_EXPORT
sunindextype* SUNSparseMatrix_IndexPointers(SUNMatrix A);

SUNDIALS_EXPORT
SUNMatrix_ID SUNMatGetID_Sparse(SUNMatrix A);

SUNDIALS_EXPORT
SUNMatrix SUNMatClone_Sparse(SUNMatrix A);

SUNDIALS_EXPORT
void SUNMatDestroy_Sparse(SUNMatrix A);

SUNDIALS_EXPORT
SUNErrCode SUNMatZero_Sparse(SUNMatrix A);

SUNDIALS_EXPORT
SUNErrCode SUNMatCopy_Sparse(SUNMatrix A, SUNMatrix B);

SUNDIALS_EXPORT
SUNErrCode SUNMatScaleAdd_Sparse(sunrealtype c, SUNMatrix A, SUNMatrix B);

SUNDIALS_EXPORT
SUNErrCode SUNMatScaleAddI_Sparse(sunrealtype c, SUNMatrix A);

SUNDIALS_EXPORT
SUNErrCode SUNMatMatvec_Sparse(SUNMatrix A, N_Vector x, N_Vector y);

SUNDIALS_EXPORT
SUNErrCode SUNMatHermitianTransposeVec_Sparse(SUNMatrix A, N_Vector x,
                                              N_Vector y);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
SUNErrCode SUNMatSpace_Sparse(SUNMatrix A, long int* lenrw, long int* leniw);

#ifdef __cplusplus
}
#endif

#endif
