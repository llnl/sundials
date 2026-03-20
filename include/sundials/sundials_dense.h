/* -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
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
 * This is the header file for a generic package of DENSE matrix
 * operations, based on the DlsMat type defined in sundials_direct.h.
 *
 * There are two sets of dense solver routines listed in
 * this file: one set uses type DlsMat defined below and the
 * other set uses the type sunrealtype ** for dense matrix arguments.
 * Routines that work with the type DlsMat begin with "Dense".
 * Routines that work with sunrealtype** begin with "dense".
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_DENSE_H
#define _SUNDIALS_DENSE_H

#include <sundials/sundials_direct.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * ----------------------------------------------------------------------------
 * Functions: SUNDlsMat_DenseGETRF and SUNDlsMat_DenseGETRS
 * ----------------------------------------------------------------------------
 * SUNDlsMat_DenseGETRF performs the LU factorization of the M by N dense matrix A.
 * This is done using standard Gaussian elimination with partial (row) pivoting.
 * Note that this applies only to matrices with M >= N and full column rank.
 *
 * A successful LU factorization leaves the matrix A and the pivot array p with
 * the following information:
 *
 * (1) p[k] contains the row number of the pivot element chosen at the beginning
 *     of elimination step k, k=0, 1, ..., N-1.
 *
 * (2) If the unique LU factorization of A is given by PA = LU, where P is a
 *     permutation matrix, L is a lower trapezoidal matrix with all 1's on the
 *     diagonal, and U is an upper triangular matrix, then the upper triangular
 *     part of A (including its diagonal) contains U and the strictly lower
 *     trapezoidal part of A contains the multipliers, I-L.
 *
 * For square matrices (M = N), L is unit lower triangular.
 *
 * SUNDlsMat_DenseGETRF returns 0 if successful. Otherwise it encountered a zero
 * diagonal element during the factorization. In this case it returns the column
 * index (numbered from one) at which it encountered the zero.
 *
 * SUNDlsMat_DenseGETRS solves the N-dimensional system A x = b using the LU
 * factorization in A and the pivot information in p computed in
 * SUNDlsMat_DenseGETRF. The solution x is returned in b. This routine cannot fail
 * if the corresponding call to SUNDlsMat_DenseGETRF did not fail.
 * SUNDlsMat_DenseGETRS does NOT check for a square matrix!
 *
 * ----------------------------------------------------------------------------
 * SUNDlsMat_DenseGETRF and SUNDlsMat_DenseGETRS are simply wrappers around
 * SUNDlsMat_denseGETRF and SUNDlsMat_denseGETRS, respectively, which perform all the
 * work by directly accessing the data in the SUNDlsMat A (i.e. in A->cols).
 * ----------------------------------------------------------------------------
 */

SUNDIALS_EXPORT
sunindextype SUNDlsMat_DenseGETRF(SUNDlsMat A, sunindextype* p);

SUNDIALS_EXPORT
void SUNDlsMat_DenseGETRS(SUNDlsMat A, sunindextype* p, sunrealtype* b);

SUNDIALS_EXPORT
sunindextype SUNDlsMat_denseGETRF(sunrealtype** a, sunindextype m,
                                  sunindextype n, sunindextype* p);

SUNDIALS_EXPORT
void SUNDlsMat_denseGETRS(sunrealtype** a, sunindextype n, sunindextype* p,
                          sunrealtype* b);

/*
 * ----------------------------------------------------------------------------
 * Functions : SUNDlsMat_DensePOTRF and SUNDlsMat_DensePOTRS
 * ----------------------------------------------------------------------------
 * SUNDlsMat_DensePOTRF computes the Cholesky factorization of a real symmetric
 * positive definite matrix A.
 * ----------------------------------------------------------------------------
 * SUNDlsMat_DensePOTRS solves a system of linear equations A*X = B with a
 * symmetric positive definite matrix A using the Cholesky factorization A =
 * L*L**T computed by SUNDlsMat_DensePOTRF.
 *
 * ----------------------------------------------------------------------------
 * SUNDlsMat_DensePOTRF and SUNDlsMat_DensePOTRS are simply wrappers around
 * SUNDlsMat_densePOTRF and SUNDlsMat_densePOTRS, respectively, which perform all the
 * work by directly accessing the data in the DlsMat A (i.e. the field cols)
 * ----------------------------------------------------------------------------
 */

SUNDIALS_EXPORT
sunindextype SUNDlsMat_DensePOTRF(SUNDlsMat A);

SUNDIALS_EXPORT
void SUNDlsMat_DensePOTRS(SUNDlsMat A, sunrealtype* b);

SUNDIALS_EXPORT
sunindextype SUNDlsMat_densePOTRF(sunrealtype** a, sunindextype m);

SUNDIALS_EXPORT
void SUNDlsMat_densePOTRS(sunrealtype** a, sunindextype m, sunrealtype* b);

/*
 * -----------------------------------------------------------------------------
 * Functions : SUNDlsMat_DenseGEQRF and SUNDlsMat_DenseORMQR
 * -----------------------------------------------------------------------------
 * SUNDlsMat_DenseGEQRF computes a QR factorization of a real M-by-N matrix A: A =
 * Q * R (with M>= N).
 *
 * SUNDlsMat_DenseGEQRF requires a temporary work vector wrk of length M.
 * -----------------------------------------------------------------------------
 * SUNDlsMat_DenseORMQR computes the product w = Q * v where Q is a real orthogonal
 * matrix defined as the product of k elementary reflectors
 *
 *        Q = H(1) H(2) . . . H(k)
 *
 * as returned by SUNDlsMat_DenseGEQRF. Q is an M-by-N matrix, v is a vector of
 * length N and w is a vector of length M (with M >= N).
 *
 * SUNDlsMat_DenseORMQR requires a temporary work vector wrk of length M.
 *
 * -----------------------------------------------------------------------------
 * SUNDlsMat_DenseGEQRF and SUNDlsMat_DenseORMQR are simply wrappers around
 * SUNDlsMat_denseGEQRF and SUNDlsMat_denseORMQR, respectively, which perform all the
 * work by directly accessing the data in the DlsMat A (i.e. the field cols)
 * -----------------------------------------------------------------------------
 */

SUNDIALS_EXPORT
int SUNDlsMat_DenseGEQRF(SUNDlsMat A, sunrealtype* beta, sunrealtype* wrk);

SUNDIALS_EXPORT
int SUNDlsMat_DenseORMQR(SUNDlsMat A, sunrealtype* beta, sunrealtype* vn,
                         sunrealtype* vm, sunrealtype* wrk);

SUNDIALS_EXPORT
int SUNDlsMat_denseGEQRF(sunrealtype** a, sunindextype m, sunindextype n,
                         sunrealtype* beta, sunrealtype* wrk);

SUNDIALS_EXPORT
int SUNDlsMat_denseORMQR(sunrealtype** a, sunindextype m, sunindextype n,
                         sunrealtype* beta, sunrealtype* v, sunrealtype* w,
                         sunrealtype* wrk);

/*
 * ----------------------------------------------------------------------------
 * Function : SUNDlsMat_DenseCopy
 * ----------------------------------------------------------------------------
 * SUNDlsMat_DenseCopy copies the contents of the M-by-N matrix A into the
 * M-by-N matrix B.
 *
 * SUNDlsMat_DenseCopy is a wrapper around SUNDlsMat_denseCopy which accesses
 * the data in the SUNDlsMat A and SUNDlsMat B (i.e. the fields cols)
 * -----------------------------------------------------------------------------
 */

SUNDIALS_EXPORT
void SUNDlsMat_DenseCopy(SUNDlsMat A, SUNDlsMat B);

SUNDIALS_EXPORT
void SUNDlsMat_denseCopy(sunrealtype** a, sunrealtype** b, sunindextype m,
                         sunindextype n);

/*
 * -----------------------------------------------------------------------------
 * Function: SUNDlsMat_DenseScale
 * -----------------------------------------------------------------------------
 * SUNDlsMat_DenseScale scales the elements of the M-by-N matrix A by the
 * constant c and stores the result back in A.
 *
 * SUNDlsMat_DenseScale is a wrapper around SUNDlsMat_denseScale which performs
 * the actual scaling by accessing the data in the SUNDlsMat A (i.e. in
 * A->cols).
 * -----------------------------------------------------------------------------
 */

SUNDIALS_EXPORT
void SUNDlsMat_DenseScale(sunrealtype c, SUNDlsMat A);

SUNDIALS_EXPORT
void SUNDlsMat_denseScale(sunrealtype c, sunrealtype** a, sunindextype m,
                          sunindextype n);

/*
 * -----------------------------------------------------------------------------
 * Function: SUNDlsMat_denseAddIdentity
 * -----------------------------------------------------------------------------
 * SUNDlsMat_denseAddIdentity adds the identity matrix to the n-by-n matrix
 * stored in a sunrealtype** array.
 * -----------------------------------------------------------------------------
 */

SUNDIALS_EXPORT
void SUNDlsMat_denseAddIdentity(sunrealtype** a, sunindextype n);

/*
 * -----------------------------------------------------------------------------
 * Function: SUNDlsMat_DenseMatvec
 * -----------------------------------------------------------------------------
 * SUNDlsMat_DenseMatvec computes the matrix-vector product y = A*x, where A is
 * an M-by-N matrix, x is a vector of length N, and y is a vector of length M.
 * No error checking is performed on the length of the arrays x and y.  Only y
 * is modified in this routine.
 *
 * SUNDlsMat_DenseMatvec is a wrapper around SUNDlsMat_denseMatvec which
 * performs the actual product by accessing the data in the SUNDlsMat A.
 * -----------------------------------------------------------------------------
 */

SUNDIALS_EXPORT
void SUNDlsMat_DenseMatvec(SUNDlsMat A, sunrealtype* x, sunrealtype* y);

SUNDIALS_EXPORT
void SUNDlsMat_denseMatvec(sunrealtype** a, sunrealtype* x, sunrealtype* y,
                           sunindextype m, sunindextype n);

#ifdef __cplusplus
}
#endif

#endif
