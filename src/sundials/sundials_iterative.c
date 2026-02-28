/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 *                Shelby Lockhart @ LLNL
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
 * This is the implementation file for the iterative.h header
 * file. It contains the implementation of functions that may be
 * useful for many different iterative solvers of A x = b.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>

#include "sundials/sundials_errors.h"
#include "sundials_iterative_impl.h"

#define FACTOR SUN_RCONST(1000.0)
#define ZERO   SUN_RCONST(0.0)
#define ONE    SUN_RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * Function : SUNModifiedGS
 * -----------------------------------------------------------------
 * This implementation of SUNModifiedGS is a slight modification of
 * a previous modified Gram-Schmidt routine (called mgs) written by
 * Milo Dorr.
 * -----------------------------------------------------------------
 */

SUNErrCode SUNModifiedGS(N_Vector* v, sunrealtype** h, int k, int p,
                         sunrealtype* new_vk_norm)
{
  SUNFunctionBegin(v[0]->sunctx);

  int i, k_minus_1, i0;
  sunrealtype new_norm_2, new_product, vk_norm, temp;

  vk_norm = N_VDotProd(v[k], v[k]);
  SUNCheckLastErr();
  vk_norm   = SUNRsqrt(vk_norm);
  k_minus_1 = k - 1;
  i0        = SUNMAX(k - p, 0);

  /* Perform modified Gram-Schmidt */

  for (i = i0; i < k; i++)
  {
    h[i][k_minus_1] = N_VDotProd(v[i], v[k]);
    SUNCheckLastErr();
    N_VLinearSum(ONE, v[k], -h[i][k_minus_1], v[i], v[k]);
    SUNCheckLastErr();
  }

  /* Compute the norm of the new vector at v[k] */

  *new_vk_norm = N_VDotProd(v[k], v[k]);
  SUNCheckLastErr();
  *new_vk_norm = SUNRsqrt(*new_vk_norm);

  /* If the norm of the new vector at v[k] is less than
     FACTOR (== 1000) times unit roundoff times the norm of the
     input vector v[k], then the vector will be reorthogonalized
     in order to ensure that nonorthogonality is not being masked
     by a very small vector length. */

  temp = FACTOR * vk_norm;
  if ((temp + (*new_vk_norm)) != temp) { return SUN_SUCCESS; }

  new_norm_2 = ZERO;

  for (i = i0; i < k; i++)
  {
    new_product = N_VDotProd(v[i], v[k]);
    SUNCheckLastErr();
    temp = FACTOR * h[i][k_minus_1];
    if ((temp + new_product) == temp) { continue; }
    h[i][k_minus_1] += new_product;
    N_VLinearSum(ONE, v[k], -new_product, v[i], v[k]);
    SUNCheckLastErr();
    new_norm_2 += SUNSQR(new_product);
  }

  if (new_norm_2 != ZERO)
  {
    new_product  = SUNSQR(*new_vk_norm) - new_norm_2;
    *new_vk_norm = (new_product > ZERO) ? SUNRsqrt(new_product) : ZERO;
  }

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Function : SUNClassicalGS
 * -----------------------------------------------------------------
 * This implementation of SUNClassicalGS was contributed by Homer
 * Walker and Peter Brown.
 * -----------------------------------------------------------------
 */

SUNErrCode SUNClassicalGS(N_Vector* v, sunrealtype** h, int k, int p,
                          sunrealtype* new_vk_norm, sunrealtype* stemp,
                          N_Vector* vtemp)
{
  SUNFunctionBegin(v[0]->sunctx);
  int i, i0, k_minus_1;
  sunrealtype vk_norm;

  k_minus_1 = k - 1;
  i0        = SUNMAX(k - p, 0);

  /* Perform Classical Gram-Schmidt */

  SUNCheckCall(N_VDotProdMulti(k - i0 + 1, v[k], v + i0, stemp));

  vk_norm = SUNRsqrt(stemp[k - i0]);
  for (i = k - i0 - 1; i >= 0; i--)
  {
    h[i][k_minus_1] = stemp[i];
    stemp[i + 1]    = -stemp[i];
    vtemp[i + 1]    = v[i];
  }
  stemp[0] = ONE;
  vtemp[0] = v[k];

  SUNCheckCall(N_VLinearCombination(k - i0 + 1, stemp, vtemp, v[k]));

  /* Compute the norm of the new vector at v[k] */

  *new_vk_norm = SUNRsqrt(N_VDotProd(v[k], v[k]));
  SUNCheckLastErr();

  /* Reorthogonalize if necessary */

  if ((FACTOR * (*new_vk_norm)) < vk_norm)
  {
    SUNCheckCall(N_VDotProdMulti(k - i0, v[k], v + i0, stemp + 1));

    stemp[0] = ONE;
    vtemp[0] = v[k];
    for (i = i0; i < k; i++)
    {
      h[i][k_minus_1] += stemp[i - i0 + 1];
      stemp[i - i0 + 1] = -stemp[i - i0 + 1];
      vtemp[i - i0 + 1] = v[i - i0];
    }

    SUNCheckCall(N_VLinearCombination(k + 1, stemp, vtemp, v[k]));

    *new_vk_norm = SUNRsqrt(N_VDotProd(v[k], v[k]));
    SUNCheckLastErr();
  }

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRfact
 * -----------------------------------------------------------------
 * This implementation of SUNQRfact is a slight modification of a
 * previous routine (called qrfact) written by Milo Dorr.
 * -----------------------------------------------------------------
 */

int SUNQRfact(int n, sunrealtype** h, sunrealtype* q, int job)
{
  sunrealtype c, s, temp1, temp2, temp3;
  int i, j, k, q_ptr, n_minus_1, code = 0;

  switch (job)
  {
  case 0:

    /* Compute a new factorization of H */

    code = 0;
    for (k = 0; k < n; k++)
    {
      /* Multiply column k by the previous k-1 Givens rotations */

      for (j = 0; j < k - 1; j++)
      {
        i           = 2 * j;
        temp1       = h[j][k];
        temp2       = h[j + 1][k];
        c           = q[i];
        s           = q[i + 1];
        h[j][k]     = c * temp1 - s * temp2;
        h[j + 1][k] = s * temp1 + c * temp2;
      }

      /* Compute the Givens rotation components c and s */

      q_ptr = 2 * k;
      temp1 = h[k][k];
      temp2 = h[k + 1][k];
      if (temp2 == ZERO)
      {
        c = ONE;
        s = ZERO;
      }
      else if (SUNRabs(temp2) >= SUNRabs(temp1))
      {
        temp3 = temp1 / temp2;
        s     = -ONE / SUNRsqrt(ONE + SUNSQR(temp3));
        c     = -s * temp3;
      }
      else
      {
        temp3 = temp2 / temp1;
        c     = ONE / SUNRsqrt(ONE + SUNSQR(temp3));
        s     = -c * temp3;
      }
      q[q_ptr]     = c;
      q[q_ptr + 1] = s;
      if ((h[k][k] = c * temp1 - s * temp2) == ZERO) { code = k + 1; }
    }
    break;

  default:

    /* Update the factored H to which a new column has been added */

    n_minus_1 = n - 1;
    code      = 0;

    /* Multiply the new column by the previous n-1 Givens rotations */

    for (k = 0; k < n_minus_1; k++)
    {
      i                   = 2 * k;
      temp1               = h[k][n_minus_1];
      temp2               = h[k + 1][n_minus_1];
      c                   = q[i];
      s                   = q[i + 1];
      h[k][n_minus_1]     = c * temp1 - s * temp2;
      h[k + 1][n_minus_1] = s * temp1 + c * temp2;
    }

    /* Compute new Givens rotation and multiply it times the last two
       entries in the new column of H.  Note that the second entry of
       this product will be 0, so it is not necessary to compute it. */

    temp1 = h[n_minus_1][n_minus_1];
    temp2 = h[n][n_minus_1];
    if (temp2 == ZERO)
    {
      c = ONE;
      s = ZERO;
    }
    else if (SUNRabs(temp2) >= SUNRabs(temp1))
    {
      temp3 = temp1 / temp2;
      s     = -ONE / SUNRsqrt(ONE + SUNSQR(temp3));
      c     = -s * temp3;
    }
    else
    {
      temp3 = temp2 / temp1;
      c     = ONE / SUNRsqrt(ONE + SUNSQR(temp3));
      s     = -c * temp3;
    }
    q_ptr        = 2 * n_minus_1;
    q[q_ptr]     = c;
    q[q_ptr + 1] = s;
    if ((h[n_minus_1][n_minus_1] = c * temp1 - s * temp2) == ZERO) { code = n; }
  }

  return (code);
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRsol
 * -----------------------------------------------------------------
 * This implementation of SUNQRsol is a slight modification of a
 * previous routine (called qrsol) written by Milo Dorr.
 * -----------------------------------------------------------------
 */

int SUNQRsol(int n, sunrealtype** h, sunrealtype* q, sunrealtype* b)
{
  sunrealtype c, s, temp1, temp2;
  int i, k, q_ptr, code = 0;

  /* Compute Q*b */

  for (k = 0; k < n; k++)
  {
    q_ptr    = 2 * k;
    c        = q[q_ptr];
    s        = q[q_ptr + 1];
    temp1    = b[k];
    temp2    = b[k + 1];
    b[k]     = c * temp1 - s * temp2;
    b[k + 1] = s * temp1 + c * temp2;
  }

  /* Solve  R*x = Q*b */

  for (k = n - 1; k >= 0; k--)
  {
    if (h[k][k] == ZERO)
    {
      code = k + 1;
      break;
    }
    b[k] /= h[k][k];
    for (i = 0; i < k; i++) { b[i] -= b[k] * h[i][k]; }
  }

  return (code);
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRAdd_MGS
 * -----------------------------------------------------------------
 * Implementation of QRAdd to be called in Anderson Acceleration
 * -----------------------------------------------------------------
 */

SUNErrCode SUNQRAdd_MGS(N_Vector* Q, sunrealtype* R, N_Vector df, int m,
                        int mMax, void* QRdata)
{
  SUNFunctionBegin(Q[0]->sunctx);

  sunindextype j;
  SUNQRData qrdata = (SUNQRData)QRdata;

  N_VScale(ONE, df, qrdata->vtemp);
  SUNCheckLastErr();
  for (j = 0; j < m; j++)
  {
    R[m * mMax + j] = N_VDotProd(Q[j], qrdata->vtemp);
    SUNCheckLastErr();
    N_VLinearSum(ONE, qrdata->vtemp, -R[m * mMax + j], Q[j], qrdata->vtemp);
    SUNCheckLastErr();
  }
  R[m * mMax + m] = N_VDotProd(qrdata->vtemp, qrdata->vtemp);
  SUNCheckLastErr();
  R[m * mMax + m] = SUNRsqrt(R[m * mMax + m]);
  N_VScale((1 / R[m * mMax + m]), qrdata->vtemp, Q[m]);
  SUNCheckLastErr();

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRAdd_ICWY
 * -----------------------------------------------------------------
 * Low synchronous implementation of QRAdd to be called in
 * Anderson Acceleration.
 * -----------------------------------------------------------------
 */

SUNErrCode SUNQRAdd_ICWY(N_Vector* Q, sunrealtype* R, N_Vector df, int m,
                         int mMax, void* QRdata)
{
  SUNFunctionBegin(Q[0]->sunctx);
  sunindextype j, k;
  SUNQRData qrdata = (SUNQRData)QRdata;

  N_VScale(ONE, df, qrdata->vtemp);
  SUNCheckLastErr(); /* stores d_fi in temp */

  if (m > 0)
  {
    /* T(1:k-1,k-1)^T = Q(:,1:k-1)^T * Q(:,k-1) */
    SUNCheckCall(
      N_VDotProdMulti(m, Q[m - 1], Q, qrdata->temp_array + (m - 1) * mMax));

    /* T(k-1,k-1) = 1.0 */
    qrdata->temp_array[(m - 1) * mMax + (m - 1)] = ONE;

    /* R(1:k-1,k) = Q_k-1^T * df */
    SUNCheckCall(N_VDotProdMulti(m, qrdata->vtemp, Q, R + m * mMax));

    /* Solve T^T * R(1:k-1,k) = R(1:k-1,k) */
    for (k = 0; k < m; k++)
    {
      /* Skip setting the diagonal element because it doesn't change */
      for (j = k + 1; j < m; j++)
      {
        R[m * mMax + j] -= R[m * mMax + k] * qrdata->temp_array[j * mMax + k];
      }
    }
    /* end */

    /* Q(:,k-1) = df - Q_k-1 R(1:k-1,k) */
    SUNCheckCall(N_VLinearCombination(m, R + m * mMax, Q, qrdata->vtemp2));
    N_VLinearSum(ONE, qrdata->vtemp, -ONE, qrdata->vtemp2, qrdata->vtemp);
    SUNCheckLastErr();
  }

  /* R(k,k) = \| df \| */
  R[m * mMax + m] = N_VDotProd(qrdata->vtemp, qrdata->vtemp);
  SUNCheckLastErr();
  R[m * mMax + m] = SUNRsqrt(R[m * mMax + m]);
  /* Q(:,k) = df / \| df \| */
  N_VScale((1 / R[m * mMax + m]), qrdata->vtemp, Q[m]);
  SUNCheckLastErr();

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRAdd_ICWY_SB
 * -----------------------------------------------------------------
 * Low synchronous implementation of QRAdd to be called in
 * Anderson Acceleration which utilizes a single buffer reduction.
 * -----------------------------------------------------------------
 */

SUNErrCode SUNQRAdd_ICWY_SB(N_Vector* Q, sunrealtype* R, N_Vector df, int m,
                            int mMax, void* QRdata)
{
  SUNFunctionBegin(Q[0]->sunctx);
  sunindextype j, k;
  SUNQRData qrdata = (SUNQRData)QRdata;

  N_VScale(ONE, df, qrdata->vtemp);
  SUNCheckLastErr(); /* stores d_fi in temp */

  if (m > 0)
  {
    /* T(1:k-1,k-1)^T = Q(:,1:k-1)^T * Q(:,k-1) */
    SUNCheckCall(N_VDotProdMultiLocal(m, Q[m - 1], Q,
                                      qrdata->temp_array + (m - 1) * mMax));

    /* R(1:k-1,k) = Q_k-1^T * df */
    /* Put R values at end of temp_array */
    SUNCheckCall(N_VDotProdMultiLocal(m, qrdata->vtemp, Q,
                                      qrdata->temp_array + (m - 1) * mMax + m));

    SUNCheckCall(N_VDotProdMultiAllReduce(m + m, qrdata->vtemp,
                                          qrdata->temp_array + (m - 1) * mMax));

    /* Move the last values from temp array into R */
    for (k = 0; k < m; k++)
    {
      R[m * mMax + k] = qrdata->temp_array[(m - 1) * mMax + m + k];
    }

    /* T(k-1,k-1) = 1.0 */
    qrdata->temp_array[(m - 1) * mMax + (m - 1)] = ONE;

    /* Solve T^T * R(1:k-1,k) = R(1:k-1,k) */
    for (k = 0; k < m; k++)
    {
      /* Skip setting the diagonal element because it doesn't change */
      for (j = k + 1; j < m; j++)
      {
        R[m * mMax + j] -= R[m * mMax + k] * qrdata->temp_array[j * mMax + k];
      }
    }
    /* end */

    /* Q(:,k-1) = df - Q_k-1 R(1:k-1,k) */
    SUNCheckCall(N_VLinearCombination(m, R + m * mMax, Q, qrdata->vtemp2));
    N_VLinearSum(ONE, qrdata->vtemp, -ONE, qrdata->vtemp2, qrdata->vtemp);
    SUNCheckLastErr();
  }

  /* R(k,k) = \| df \| */
  R[m * mMax + m] = N_VDotProd(qrdata->vtemp, qrdata->vtemp);
  SUNCheckLastErr();
  R[m * mMax + m] = SUNRsqrt(R[m * mMax + m]);
  /* Q(:,k) = df / \| df \| */
  N_VScale((1 / R[m * mMax + m]), qrdata->vtemp, Q[m]);
  SUNCheckLastErr();

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRAdd_CGS2
 * -----------------------------------------------------------------
 * Low synchronous Implementation of QRAdd to be called in
 * Anderson Acceleration.
 * -----------------------------------------------------------------
 */

SUNErrCode SUNQRAdd_CGS2(N_Vector* Q, sunrealtype* R, N_Vector df, int m,
                         int mMax, void* QRdata)
{
  SUNFunctionBegin(Q[0]->sunctx);
  sunindextype j;
  SUNQRData qrdata = (SUNQRData)QRdata;

  N_VScale(ONE, df, qrdata->vtemp);
  SUNCheckLastErr(); /* temp = df */

  if (m > 0)
  {
    /* s_k = Q_k-1^T df_aa -- update with sdata as a sunrealtype* array */
    SUNCheckCall(N_VDotProdMulti(m, qrdata->vtemp, Q, R + m * mMax));

    /* y = df - Q_k-1 s_k */
    SUNCheckCall(N_VLinearCombination(m, R + m * mMax, Q, qrdata->vtemp2));
    N_VLinearSum(ONE, qrdata->vtemp, -ONE, qrdata->vtemp2, qrdata->vtemp2);
    SUNCheckLastErr();

    /* z_k = Q_k-1^T y */
    SUNCheckCall(N_VDotProdMulti(m, qrdata->vtemp2, Q, qrdata->temp_array));

    /* df = y - Q_k-1 z_k */
    SUNCheckCall(N_VLinearCombination(m, qrdata->temp_array, Q, Q[m]));
    N_VLinearSum(ONE, qrdata->vtemp2, -ONE, Q[m], qrdata->vtemp);
    SUNCheckLastErr();

    /* R(1:k-1,k) = s_k + z_k */
    for (j = 0; j < m; j++)
    {
      R[m * mMax + j] = R[m * mMax + j] + qrdata->temp_array[j];
    }
  }

  /* R(k,k) = \| df \| */
  R[m * mMax + m] = N_VDotProd(qrdata->vtemp, qrdata->vtemp);
  SUNCheckLastErr();
  R[m * mMax + m] = SUNRsqrt(R[m * mMax + m]);
  /* Q(:,k) = df / R(k,k) */
  N_VScale((1 / R[m * mMax + m]), qrdata->vtemp, Q[m]);
  SUNCheckLastErr();

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRAdd_DCGS2
 * -----------------------------------------------------------------
 * Low synchronous Implementation of QRAdd to be called in
 * Anderson Acceleration.
 * -----------------------------------------------------------------
 */

SUNErrCode SUNQRAdd_DCGS2(N_Vector* Q, sunrealtype* R, N_Vector df, int m,
                          int mMax, void* QRdata)
{
  SUNFunctionBegin(Q[0]->sunctx);
  sunindextype j;
  SUNQRData qrdata = (SUNQRData)QRdata;

  N_VScale(ONE, df, qrdata->vtemp);
  SUNCheckLastErr(); /* temp = df */

  if (m > 0)
  {
    /* R(1:k-1,k) = Q_k-1^T df_aa */
    SUNCheckCall(N_VDotProdMulti(m, qrdata->vtemp, Q, R + m * mMax));
    /* Delayed reorthogonalization */
    if (m > 1)
    {
      /* s = Q_k-2^T Q(:,k-1) */
      SUNCheckCall(N_VDotProdMulti(m - 1, Q[m - 1], Q, qrdata->temp_array));

      /* Q(:,k-1) = Q(:,k-1) - Q_k-2 s */
      SUNCheckCall(
        N_VLinearCombination(m - 1, qrdata->temp_array, Q, qrdata->vtemp2));
      N_VLinearSum(ONE, Q[m - 1], -ONE, qrdata->vtemp2, Q[m - 1]);
      SUNCheckLastErr();

      /* R(1:k-2,k-1) = R(1:k-2,k-1) + s */
      for (j = 0; j < m - 1; j++)
      {
        R[(m - 1) * mMax + j] = R[(m - 1) * mMax + j] + qrdata->temp_array[j];
      }
    }

    /* df = df - Q(:,k-1) R(1:k-1,k) */
    SUNCheckCall(N_VLinearCombination(m, R + m * mMax, Q, qrdata->vtemp2));
    N_VLinearSum(ONE, qrdata->vtemp, -ONE, qrdata->vtemp2, qrdata->vtemp);
    SUNCheckLastErr();
  }

  /* R(k,k) = \| df \| */
  R[m * mMax + m] = N_VDotProd(qrdata->vtemp, qrdata->vtemp);
  SUNCheckLastErr();
  R[m * mMax + m] = SUNRsqrt(R[m * mMax + m]);
  /* Q(:,k) = df / R(k,k) */
  N_VScale((1 / R[m * mMax + m]), qrdata->vtemp, Q[m]);
  SUNCheckLastErr();

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Function : SUNQRAdd_DCGS2_SB
 * -----------------------------------------------------------------
 * Low synchronous Implementation of QRAdd to be called in
 * Anderson Acceleration which utilizes a single buffer reduction.
 * -----------------------------------------------------------------
 */

SUNErrCode SUNQRAdd_DCGS2_SB(N_Vector* Q, sunrealtype* R, N_Vector df, int m,
                             int mMax, void* QRdata)
{
  SUNFunctionBegin(Q[0]->sunctx);
  sunindextype j;
  SUNQRData qrdata = (SUNQRData)QRdata;

  N_VScale(ONE, df, qrdata->vtemp);
  SUNCheckLastErr(); /* temp = df */

  if (m > 0)
  {
    if (m == 1)
    {
      /* R(1:k-1,k) = Q_k-1^T df_aa */
      SUNCheckCall(N_VDotProdMulti(m, qrdata->vtemp, Q, R + m * mMax));
    }
    /* Delayed reorthogonalization */
    else if (m > 1)
    {
      /* R(1:k-1,k) = Q_k-1^T df_aa */
      /* Put R values at beginning of temp array */
      SUNCheckCall(N_VDotProdMultiLocal(m, qrdata->vtemp, Q, qrdata->temp_array));

      /* s = Q_k-2^T Q(:,k-1) */
      SUNCheckCall(
        N_VDotProdMultiLocal(m - 1, Q[m - 1], Q, qrdata->temp_array + m));
      SUNCheckCall(
        N_VDotProdMultiAllReduce(m + m - 1, qrdata->vtemp, qrdata->temp_array));

      /* Move R values to R */
      for (j = 0; j < m; j++) { R[m * mMax + j] = qrdata->temp_array[j]; }

      /* Q(:,k-1) = Q(:,k-1) - Q_k-2 s */
      SUNCheckCall(
        N_VLinearCombination(m - 1, qrdata->temp_array + m, Q, qrdata->vtemp2));
      N_VLinearSum(ONE, Q[m - 1], -ONE, qrdata->vtemp2, Q[m - 1]);
      SUNCheckLastErr();

      /* R(1:k-2,k-1) = R(1:k-2,k-1) + s */
      for (j = 0; j < m - 1; j++)
      {
        R[(m - 1) * mMax + j] = R[(m - 1) * mMax + j] + qrdata->temp_array[m + j];
      }
    }

    /* df = df - Q(:,k-1) R(1:k-1,k) */
    SUNCheckCall(N_VLinearCombination(m, R + m * mMax, Q, qrdata->vtemp2));
    N_VLinearSum(ONE, qrdata->vtemp, -ONE, qrdata->vtemp2, qrdata->vtemp);
    SUNCheckLastErr();
  }

  /* R(k,k) = \| df \| */
  R[m * mMax + m] = N_VDotProd(qrdata->vtemp, qrdata->vtemp);
  SUNCheckLastErr();
  R[m * mMax + m] = SUNRsqrt(R[m * mMax + m]);
  /* Q(:,k) = df / R(k,k) */
  N_VScale((1 / R[m * mMax + m]), qrdata->vtemp, Q[m]);
  SUNCheckLastErr();

  return SUN_SUCCESS;
}
