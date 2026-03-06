/*
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
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
 * This is the implementation file for a generic package of dense
 * matrix operations.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_math.h>

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)

/*
 * -----------------------------------------------------
 * Functions working on SUNDlsMat
 * -----------------------------------------------------
 */

sunindextype SUNDlsMat_DenseGETRF(SUNDlsMat A, sunindextype* p)
{
  return (SUNDlsMat_denseGETRF(A->cols, A->M, A->N, p));
}

void SUNDlsMat_DenseGETRS(SUNDlsMat A, sunindextype* p, sunrealtype* b)
{
  SUNDlsMat_denseGETRS(A->cols, A->N, p, b);
}

sunindextype SUNDlsMat_DensePOTRF(SUNDlsMat A)
{
  return (SUNDlsMat_densePOTRF(A->cols, A->M));
}

void SUNDlsMat_DensePOTRS(SUNDlsMat A, sunrealtype* b)
{
  SUNDlsMat_densePOTRS(A->cols, A->M, b);
}

int SUNDlsMat_DenseGEQRF(SUNDlsMat A, sunrealtype* beta, sunrealtype* wrk)
{
  return (SUNDlsMat_denseGEQRF(A->cols, A->M, A->N, beta, wrk));
}

int SUNDlsMat_DenseORMQR(SUNDlsMat A, sunrealtype* beta, sunrealtype* vn,
                         sunrealtype* vm, sunrealtype* wrk)
{
  return (SUNDlsMat_denseORMQR(A->cols, A->M, A->N, beta, vn, vm, wrk));
}

void SUNDlsMat_DenseCopy(SUNDlsMat A, SUNDlsMat B)
{
  SUNDlsMat_denseCopy(A->cols, B->cols, A->M, A->N);
}

void SUNDlsMat_DenseScale(sunrealtype c, SUNDlsMat A)
{
  SUNDlsMat_denseScale(c, A->cols, A->M, A->N);
}

void SUNDlsMat_DenseMatvec(SUNDlsMat A, sunrealtype* x, sunrealtype* y)
{
  SUNDlsMat_denseMatvec(A->cols, x, y, A->M, A->N);
}

sunindextype SUNDlsMat_denseGETRF(sunrealtype** a, sunindextype m,
                                  sunindextype n, sunindextype* p)
{
  sunindextype i, j, k, l;
  sunrealtype *col_j, *col_k;
  sunrealtype temp, mult, a_kj;

  /* k-th elimination step number */
  for (k = 0; k < n; k++)
  {
    col_k = a[k];

    /* find l = pivot row number */
    l = k;
    for (i = k + 1; i < m; i++)
    {
      if (SUNRabs(col_k[i]) > SUNRabs(col_k[l])) { l = i; }
    }
    p[k] = l;

    /* check for zero pivot element */
    if (col_k[l] == ZERO) { return (k + 1); }

    /* swap a(k,1:n) and a(l,1:n) if necessary */
    if (l != k)
    {
      for (i = 0; i < n; i++)
      {
        temp    = a[i][l];
        a[i][l] = a[i][k];
        a[i][k] = temp;
      }
    }

    /* Scale the elements below the diagonal in
     * column k by 1.0/a(k,k). After the above swap
     * a(k,k) holds the pivot element. This scaling
     * stores the pivot row multipliers a(i,k)/a(k,k)
     * in a(i,k), i=k+1, ..., m-1.
     */
    mult = ONE / col_k[k];
    for (i = k + 1; i < m; i++) { col_k[i] *= mult; }

    /* row_i = row_i - [a(i,k)/a(k,k)] row_k, i=k+1, ..., m-1 */
    /* row k is the pivot row after swapping with row l.      */
    /* The computation is done one column at a time,          */
    /* column j=k+1, ..., n-1.                                */

    for (j = k + 1; j < n; j++)
    {
      col_j = a[j];
      a_kj  = col_j[k];

      /* a(i,j) = a(i,j) - [a(i,k)/a(k,k)]*a(k,j)  */
      /* a_kj = a(k,j), col_k[i] = - a(i,k)/a(k,k) */

      if (a_kj != ZERO)
      {
        for (i = k + 1; i < m; i++) { col_j[i] -= a_kj * col_k[i]; }
      }
    }
  }

  /* return 0 to indicate success */

  return (0);
}

void SUNDlsMat_denseGETRS(sunrealtype** a, sunindextype n, sunindextype* p,
                          sunrealtype* b)
{
  sunindextype i, k, pk;
  sunrealtype *col_k, tmp;

  /* Permute b, based on pivot information in p */
  for (k = 0; k < n; k++)
  {
    pk = p[k];
    if (pk != k)
    {
      tmp   = b[k];
      b[k]  = b[pk];
      b[pk] = tmp;
    }
  }

  /* Solve Ly = b, store solution y in b */
  for (k = 0; k < n - 1; k++)
  {
    col_k = a[k];
    for (i = k + 1; i < n; i++) { b[i] -= col_k[i] * b[k]; }
  }

  /* Solve Ux = y, store solution x in b */
  for (k = n - 1; k > 0; k--)
  {
    col_k = a[k];
    b[k] /= col_k[k];
    for (i = 0; i < k; i++) { b[i] -= col_k[i] * b[k]; }
  }
  b[0] /= a[0][0];
}

/*
 * Cholesky decomposition of a symmetric positive-definite matrix
 * A = C^T*C: gaxpy version.
 * Only the lower triangle of A is accessed and it is overwritten with
 * the lower triangle of C.
 */

sunindextype SUNDlsMat_densePOTRF(sunrealtype** a, sunindextype m)
{
  sunrealtype *a_col_j, *a_col_k;
  sunrealtype a_diag;
  sunindextype i, j, k;

  for (j = 0; j < m; j++)
  {
    a_col_j = a[j];

    if (j > 0)
    {
      for (i = j; i < m; i++)
      {
        for (k = 0; k < j; k++)
        {
          a_col_k = a[k];
          a_col_j[i] -= a_col_k[i] * a_col_k[j];
        }
      }
    }

    a_diag = a_col_j[j];
    if (a_diag <= ZERO) { return (j + 1); }
    a_diag = SUNRsqrt(a_diag);

    for (i = j; i < m; i++) { a_col_j[i] /= a_diag; }
  }

  return (0);
}

/*
 * Solution of Ax=b, with A s.p.d., based on the Cholesky decomposition
 * obtained with denPOTRF.; A = C*C^T, C lower triangular
 *
 */

void SUNDlsMat_densePOTRS(sunrealtype** a, sunindextype m, sunrealtype* b)
{
  sunrealtype *col_j, *col_i;
  sunindextype i, j;

  /* Solve C y = b, forward substitution - column version.
     Store solution y in b */
  for (j = 0; j < m - 1; j++)
  {
    col_j = a[j];
    b[j] /= col_j[j];
    for (i = j + 1; i < m; i++) { b[i] -= b[j] * col_j[i]; }
  }
  col_j = a[m - 1];
  b[m - 1] /= col_j[m - 1];

  /* Solve C^T x = y, backward substitution - row version.
     Store solution x in b */
  col_j = a[m - 1];
  b[m - 1] /= col_j[m - 1];
  for (i = m - 2; i >= 0; i--)
  {
    col_i = a[i];
    for (j = i + 1; j < m; j++) { b[i] -= col_i[j] * b[j]; }
    b[i] /= col_i[i];
  }
}

/*
 * QR factorization of a rectangular matrix A of size m by n (m >= n)
 * using Householder reflections.
 *
 * On exit, the elements on and above the diagonal of A contain the n by n
 * upper triangular matrix R; the elements below the diagonal, with the array beta,
 * represent the orthogonal matrix Q as a product of elementary reflectors .
 *
 * v (of length m) must be provided as workspace.
 *
 */

int SUNDlsMat_denseGEQRF(sunrealtype** a, sunindextype m, sunindextype n,
                         sunrealtype* beta, sunrealtype* v)
{
  sunrealtype ajj, s, mu, v1, v1_2;
  sunrealtype *col_j, *col_k;
  sunindextype i, j, k;

  /* For each column...*/
  for (j = 0; j < n; j++)
  {
    col_j = a[j];

    ajj = col_j[j];

    /* Compute the j-th Householder vector (of length m-j) */
    v[0] = ONE;
    s    = ZERO;
    for (i = 1; i < m - j; i++)
    {
      v[i] = col_j[i + j];
      s += v[i] * v[i];
    }

    if (s != ZERO)
    {
      mu      = SUNRsqrt(ajj * ajj + s);
      v1      = (ajj <= ZERO) ? ajj - mu : -s / (ajj + mu);
      v1_2    = v1 * v1;
      beta[j] = TWO * v1_2 / (s + v1_2);
      for (i = 1; i < m - j; i++) { v[i] /= v1; }
    }
    else { beta[j] = ZERO; }

    /* Update upper triangle of A (load R) */
    for (k = j; k < n; k++)
    {
      col_k = a[k];
      s     = ZERO;
      for (i = 0; i < m - j; i++) { s += col_k[i + j] * v[i]; }
      s *= beta[j];
      for (i = 0; i < m - j; i++) { col_k[i + j] -= s * v[i]; }
    }

    /* Update A (load Householder vector) */
    if (j < m - 1)
    {
      for (i = 1; i < m - j; i++) { col_j[i + j] = v[i]; }
    }
  }

  return (0);
}

/*
 * Computes vm = Q * vn, where the orthogonal matrix Q is stored as
 * elementary reflectors in the m by n matrix A and in the vector beta.
 * (NOTE: It is assumed that an QR factorization has been previously
 * computed with denGEQRF).
 *
 * vn (IN) has length n, vm (OUT) has length m, and it's assumed that m >= n.
 *
 * v (of length m) must be provided as workspace.
 */

int SUNDlsMat_denseORMQR(sunrealtype** a, sunindextype m, sunindextype n,
                         sunrealtype* beta, sunrealtype* vn, sunrealtype* vm,
                         sunrealtype* v)
{
  sunrealtype *col_j, s;
  sunindextype i, j;

  /* Initialize vm */
  for (i = 0; i < n; i++) { vm[i] = vn[i]; }
  for (i = n; i < m; i++) { vm[i] = ZERO; }

  /* Accumulate (backwards) corrections into vm */
  for (j = n - 1; j >= 0; j--)
  {
    col_j = a[j];

    v[0] = ONE;
    s    = vm[j];
    for (i = 1; i < m - j; i++)
    {
      v[i] = col_j[i + j];
      s += v[i] * vm[i + j];
    }
    s *= beta[j];

    for (i = 0; i < m - j; i++) { vm[i + j] -= s * v[i]; }
  }

  return (0);
}

void SUNDlsMat_denseCopy(sunrealtype** a, sunrealtype** b, sunindextype m,
                         sunindextype n)
{
  sunindextype i, j;
  sunrealtype *a_col_j, *b_col_j;

  for (j = 0; j < n; j++)
  {
    a_col_j = a[j];
    b_col_j = b[j];
    for (i = 0; i < m; i++) { b_col_j[i] = a_col_j[i]; }
  }
}

void SUNDlsMat_denseScale(sunrealtype c, sunrealtype** a, sunindextype m,
                          sunindextype n)
{
  sunindextype i, j;
  sunrealtype* col_j;

  for (j = 0; j < n; j++)
  {
    col_j = a[j];
    for (i = 0; i < m; i++) { col_j[i] *= c; }
  }
}

void SUNDlsMat_denseAddIdentity(sunrealtype** a, sunindextype n)
{
  sunindextype i;

  for (i = 0; i < n; i++) { a[i][i] += ONE; }
}

void SUNDlsMat_denseMatvec(sunrealtype** a, sunrealtype* x, sunrealtype* y,
                           sunindextype m, sunindextype n)
{
  sunindextype i, j;
  sunrealtype* col_j;

  for (i = 0; i < m; i++) { y[i] = ZERO; }

  for (j = 0; j < n; j++)
  {
    col_j = a[j];
    for (i = 0; i < m; i++) { y[i] += col_j[i] * x[j]; }
  }
}
