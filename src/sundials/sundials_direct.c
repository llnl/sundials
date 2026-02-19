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
 * This is the implementation file for operations to be used by a
 * generic direct linear solver.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_math.h>

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

SUNDlsMat SUNDlsMat_NewDenseMat(sunindextype M, sunindextype N)
{
  SUNDlsMat A;
  sunindextype j;

  if ((M <= 0) || (N <= 0)) { return (NULL); }

  A = NULL;
  A = (SUNDlsMat)malloc(sizeof *A);
  if (A == NULL) { return (NULL); }

  A->data = (sunrealtype*)malloc(M * N * sizeof(sunrealtype));
  if (A->data == NULL)
  {
    free(A);
    A = NULL;
    return (NULL);
  }
  A->cols = (sunrealtype**)malloc(N * sizeof(sunrealtype*));
  if (A->cols == NULL)
  {
    free(A->data);
    A->data = NULL;
    free(A);
    A = NULL;
    return (NULL);
  }

  for (j = 0; j < N; j++) { A->cols[j] = A->data + j * M; }

  A->M     = M;
  A->N     = N;
  A->ldim  = M;
  A->ldata = M * N;

  A->type = SUNDIALS_DENSE;

  return (A);
}

sunrealtype** SUNDlsMat_newDenseMat(sunindextype m, sunindextype n)
{
  sunindextype j;
  sunrealtype** a;

  if ((n <= 0) || (m <= 0)) { return (NULL); }

  a = NULL;
  a = (sunrealtype**)malloc(n * sizeof(sunrealtype*));
  if (a == NULL) { return (NULL); }

  a[0] = NULL;
  a[0] = (sunrealtype*)malloc(m * n * sizeof(sunrealtype));
  if (a[0] == NULL)
  {
    free(a);
    a = NULL;
    return (NULL);
  }

  for (j = 1; j < n; j++) { a[j] = a[0] + j * m; }

  return (a);
}

SUNDlsMat SUNDlsMat_NewBandMat(sunindextype N, sunindextype mu, sunindextype ml,
                               sunindextype smu)
{
  SUNDlsMat A;
  sunindextype j, colSize;

  if (N <= 0) { return (NULL); }

  A = NULL;
  A = (SUNDlsMat)malloc(sizeof *A);
  if (A == NULL) { return (NULL); }

  colSize = smu + ml + 1;
  A->data = NULL;
  A->data = (sunrealtype*)malloc(N * colSize * sizeof(sunrealtype));
  if (A->data == NULL)
  {
    free(A);
    A = NULL;
    return (NULL);
  }

  A->cols = NULL;
  A->cols = (sunrealtype**)malloc(N * sizeof(sunrealtype*));
  if (A->cols == NULL)
  {
    free(A->data);
    free(A);
    A = NULL;
    return (NULL);
  }

  for (j = 0; j < N; j++) { A->cols[j] = A->data + j * colSize; }

  A->M     = N;
  A->N     = N;
  A->mu    = mu;
  A->ml    = ml;
  A->s_mu  = smu;
  A->ldim  = colSize;
  A->ldata = N * colSize;

  A->type = SUNDIALS_BAND;

  return (A);
}

sunrealtype** SUNDlsMat_newBandMat(sunindextype n, sunindextype smu,
                                   sunindextype ml)
{
  sunrealtype** a;
  sunindextype j, colSize;

  if (n <= 0) { return (NULL); }

  a = NULL;
  a = (sunrealtype**)malloc(n * sizeof(sunrealtype*));
  if (a == NULL) { return (NULL); }

  colSize = smu + ml + 1;
  a[0]    = NULL;
  a[0]    = (sunrealtype*)malloc(n * colSize * sizeof(sunrealtype));
  if (a[0] == NULL)
  {
    free(a);
    a = NULL;
    return (NULL);
  }

  for (j = 1; j < n; j++) { a[j] = a[0] + j * colSize; }

  return (a);
}

void SUNDlsMat_DestroyMat(SUNDlsMat A)
{
  free(A->data);
  A->data = NULL;
  free(A->cols);
  free(A);
  A = NULL;
}

void SUNDlsMat_destroyMat(sunrealtype** a)
{
  free(a[0]);
  a[0] = NULL;
  free(a);
  a = NULL;
}

int* SUNDlsMat_NewIntArray(int N)
{
  int* vec;

  if (N <= 0) { return (NULL); }

  vec = NULL;
  vec = (int*)malloc(N * sizeof(int));

  return (vec);
}

int* SUNDlsMat_newIntArray(int n)
{
  int* v;

  if (n <= 0) { return (NULL); }

  v = NULL;
  v = (int*)malloc(n * sizeof(int));

  return (v);
}

sunindextype* SUNDlsMat_NewIndexArray(sunindextype N)
{
  sunindextype* vec;

  if (N <= 0) { return (NULL); }

  vec = NULL;
  vec = (sunindextype*)malloc(N * sizeof(sunindextype));

  return (vec);
}

sunindextype* SUNDlsMat_newIndexArray(sunindextype n)
{
  sunindextype* v;

  if (n <= 0) { return (NULL); }

  v = NULL;
  v = (sunindextype*)malloc(n * sizeof(sunindextype));

  return (v);
}

sunrealtype* SUNDlsMat_NewRealArray(sunindextype N)
{
  sunrealtype* vec;

  if (N <= 0) { return (NULL); }

  vec = NULL;
  vec = (sunrealtype*)malloc(N * sizeof(sunrealtype));

  return (vec);
}

sunrealtype* SUNDlsMat_newRealArray(sunindextype m)
{
  sunrealtype* v;

  if (m <= 0) { return (NULL); }

  v = NULL;
  v = (sunrealtype*)malloc(m * sizeof(sunrealtype));

  return (v);
}

void SUNDlsMat_DestroyArray(void* V)
{
  free(V);
  V = NULL;
}

void SUNDlsMat_destroyArray(void* v)
{
  free(v);
  v = NULL;
}

void SUNDlsMat_AddIdentity(SUNDlsMat A)
{
  sunindextype i;

  switch (A->type)
  {
  case SUNDIALS_DENSE:
    for (i = 0; i < A->N; i++) { A->cols[i][i] += ONE; }
    break;

  case SUNDIALS_BAND:
    for (i = 0; i < A->M; i++) { A->cols[i][A->s_mu] += ONE; }
    break;
  }
}

void SUNDlsMat_SetToZero(SUNDlsMat A)
{
  sunindextype i, j, colSize;
  sunrealtype* col_j;

  switch (A->type)
  {
  case SUNDIALS_DENSE:

    for (j = 0; j < A->N; j++)
    {
      col_j = A->cols[j];
      for (i = 0; i < A->M; i++) { col_j[i] = ZERO; }
    }

    break;

  case SUNDIALS_BAND:

    colSize = A->mu + A->ml + 1;
    for (j = 0; j < A->M; j++)
    {
      col_j = A->cols[j] + A->s_mu - A->mu;
      for (i = 0; i < colSize; i++) { col_j[i] = ZERO; }
    }

    break;
  }
}

void SUNDlsMat_PrintMat(SUNDlsMat A, FILE* outfile)
{
  sunindextype i, j, start, finish;
  sunrealtype** a;

  switch (A->type)
  {
  case SUNDIALS_DENSE:

    fprintf(outfile, "\n");
    for (i = 0; i < A->M; i++)
    {
      for (j = 0; j < A->N; j++)
      {
        fprintf(outfile, SUN_FORMAT_E "  ", SUNDLS_DENSE_ELEM(A, i, j));
      }
      fprintf(outfile, "\n");
    }

    break;

  case SUNDIALS_BAND:

    a = A->cols;
    fprintf(outfile, "\n");
    for (i = 0; i < A->N; i++)
    {
      start  = SUNMAX(0, i - A->ml);
      finish = SUNMIN(A->N - 1, i + A->mu);
      for (j = 0; j < start; j++) { fprintf(outfile, "%12s  ", ""); }
      for (j = start; j <= finish; j++)
      {
        fprintf(outfile, SUN_FORMAT_E "  ", a[j][i - j + A->s_mu]);
      }
      fprintf(outfile, "\n");
    }

    break;
  }
}
