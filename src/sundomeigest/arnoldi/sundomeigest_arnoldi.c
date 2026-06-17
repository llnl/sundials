/* -----------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ SMU
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
 * This is the implementation file for the Arnoldi Iteration
 * implementation of the SUNDomEigEst package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sundomeigest/sundomeigest_arnoldi.h>

#include "sundials_lapack_defs.h"
#include "sundials_logger_impl.h"
#include "sundials_macros.h"

/* Interfaces to match 'sunscalartype' with the correct LAPACK functions 
   and to declare the appropriate sorting function for eigenvalues based
   on the scalar type. */
#if defined(SUNDIALS_SCALAR_TYPE_REAL)
static void sundomeigest_SortEigenvaluesByMagnitude(sunrealtype* wr, sunrealtype* wi, int n);
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define xgeev_f77 dgeev_f77
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define xgeev_f77 sgeev_f77
#else
#error Incompatible sunrealtype for LAPACK; disable LAPACK and rebuild
#endif
#elif defined(SUNDIALS_SCALAR_TYPE_COMPLEX)
static void sundomeigest_SortEigenvaluesByMagnitude(sunscalartype* wr, int n);
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define xgeev_f77 zgeev_f77
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define xgeev_f77 cgeev_f77
#else
#error Incompatible sunrealtype for LAPACK; disable LAPACK and rebuild
#endif
#endif

/* Default estimator parameters */
#define DEE_NUM_OF_WARMUPS_ARNOLDI_DEFAULT 100

/* Default Arnoldi Iteration parameters */
#define DEE_KRYLOV_DIM_DEFAULT 3

/*
 * -----------------------------------------------------------------
 * Arnoldi itetation structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define Arnoldi_CONTENT(DEE) ((SUNDomEigEstimatorContent_Arnoldi)(DEE->content))

/*
 * -----------------------------------------------------------------
 * Arnoldi module private function prototypes
 * -----------------------------------------------------------------
 */

int sundomeigest_Compare(const void* a, const void* b);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new Arnoldi estimator
 */

SUNDomEigEstimator SUNDomEigEstimator_Arnoldi(N_Vector q, int kry_dim,
                                              SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNDomEigEstimator DEE;
  SUNDomEigEstimatorContent_Arnoldi content;

  /* Check if kry_dim >= 2 */
  if (kry_dim < 3) { kry_dim = DEE_KRYLOV_DIM_DEFAULT; }

  /* Input vector must be non-NULL */
  SUNAssertNull(q, SUN_ERR_ARG_CORRUPT);
  SUNAssertNull(q->ops, SUN_ERR_ARG_CORRUPT);

  /* Check for required vector operations */
  SUNAssertNull(q->ops->nvclone, SUN_ERR_ARG_INCOMPATIBLE);
  SUNAssertNull(q->ops->nvdestroy, SUN_ERR_ARG_INCOMPATIBLE);
  SUNAssertNull(q->ops->nvdotprod, SUN_ERR_ARG_INCOMPATIBLE);
  SUNAssertNull(q->ops->nvscale, SUN_ERR_ARG_INCOMPATIBLE);

  /* Check if q != 0 vector */
  sunscalartype qdotq;
  SUNCheckCallNull(N_VDotProdComplex(q, q, &qdotq));
  SUNAssertNull(SUN_REAL(qdotq) > SUN_SMALL_REAL, SUN_ERR_ARG_INCOMPATIBLE);

  /* Create dominant eigenvalue estimator */
  DEE = NULL;
  DEE = SUNDomEigEstimator_NewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  DEE->ops->setatimes = SUNDomEigEstimator_SetATimes_Arnoldi;
  DEE->ops->setnumpreprocessiters =
    SUNDomEigEstimator_SetNumPreprocessIters_Arnoldi;
  DEE->ops->setinitialguess   = SUNDomEigEstimator_SetInitialGuess_Arnoldi;
  DEE->ops->initialize        = SUNDomEigEstimator_Initialize_Arnoldi;
  DEE->ops->estimate          = SUNDomEigEstimator_Estimate_Arnoldi;
  DEE->ops->getnumiters       = SUNDomEigEstimator_GetNumIters_Arnoldi;
  DEE->ops->getnumatimescalls = SUNDomEigEstimator_GetNumATimesCalls_Arnoldi;
  DEE->ops->write             = SUNDomEigEstimator_Write_Arnoldi;
  DEE->ops->destroy           = SUNDomEigEstimator_Destroy_Arnoldi;

  /* Create content */
  content = NULL;
  content = (SUNDomEigEstimatorContent_Arnoldi)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content  */
  DEE->content = content;

  /* Fill content */
  content->ATimes       = NULL;
  content->ATdata       = NULL;
  content->V            = NULL;
  content->q            = NULL;
  content->kry_dim      = kry_dim;
  content->num_warmups  = DEE_NUM_OF_WARMUPS_ARNOLDI_DEFAULT;
  content->num_iters    = 0;
  content->num_ATimes   = 0;
  content->LAPACK_A     = NULL;
  content->LAPACK_wr    = NULL;
  content->LAPACK_wi    = NULL;
  content->LAPACK_work  = NULL;
  content->LAPACK_lwork = 0;
  content->LAPACK_rwork = NULL;
  content->LAPACK_arr   = NULL;
  content->Hes          = NULL;

  /* Allocate content */
  content->q = N_VClone(q);
  SUNCheckLastErrNull();

  N_VScale(SUN_CCONST(1.0, 0.0), q, content->q);
  SUNCheckLastErrNull();

  content->V = N_VCloneVectorArray(kry_dim + 1, q);
  SUNCheckLastErrNull();

  return (DEE);
}

/*
 * -----------------------------------------------------------------
 * implementation of dominant eigenvalue estimator operations
 * -----------------------------------------------------------------
 */

SUNErrCode SUNDomEigEstimator_SetATimes_Arnoldi(SUNDomEigEstimator DEE,
                                                void* A_data, SUNATimesFn ATimes)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  /* set function pointers to integrator-supplied ATimes routine
     and data, and return with success */
  Arnoldi_CONTENT(DEE)->ATimes = ATimes;
  Arnoldi_CONTENT(DEE)->ATdata = A_data;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_Initialize_Arnoldi(SUNDomEigEstimator DEE)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  if (Arnoldi_CONTENT(DEE)->kry_dim < 2)
  {
    Arnoldi_CONTENT(DEE)->kry_dim = DEE_KRYLOV_DIM_DEFAULT;
  }
  if (Arnoldi_CONTENT(DEE)->num_warmups < 0)
  {
    Arnoldi_CONTENT(DEE)->num_warmups = DEE_NUM_OF_WARMUPS_ARNOLDI_DEFAULT;
  }

  SUNAssert(Arnoldi_CONTENT(DEE)->ATimes, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);

  if (Arnoldi_CONTENT(DEE)->LAPACK_A == NULL)
  {
    Arnoldi_CONTENT(DEE)->LAPACK_A = (sunscalartype*)malloc(
      (Arnoldi_CONTENT(DEE)->kry_dim * Arnoldi_CONTENT(DEE)->kry_dim) *
      sizeof(sunscalartype));
    SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_A, SUN_ERR_MALLOC_FAIL);
  }
  if (Arnoldi_CONTENT(DEE)->LAPACK_wr == NULL)
  {
    Arnoldi_CONTENT(DEE)->LAPACK_wr =
      malloc(Arnoldi_CONTENT(DEE)->kry_dim * sizeof(sunscalartype));
    SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_wr, SUN_ERR_MALLOC_FAIL);
  }
  if (Arnoldi_CONTENT(DEE)->LAPACK_wi == NULL)
  {
    Arnoldi_CONTENT(DEE)->LAPACK_wi =
      malloc(Arnoldi_CONTENT(DEE)->kry_dim * sizeof(sunscalartype));
    SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_wi, SUN_ERR_MALLOC_FAIL);
  }

  /* query the workspace size (call with lwork = -1) */
  char jobvl         = 'N';
  char jobvr         = 'N';
  Arnoldi_CONTENT(DEE)->LAPACK_rwork = NULL;
  sunindextype N     = Arnoldi_CONTENT(DEE)->kry_dim;
  sunindextype lda   = Arnoldi_CONTENT(DEE)->kry_dim;
  sunindextype ldvl  = 1;
  sunindextype ldvr  = 1;
  sunindextype info  = 0;
  sunindextype lwork = -1;
  sunscalartype work = SUN_CCONST(0.0, 0.0);

#if defined(SUNDIALS_SCALAR_TYPE_COMPLEX)
  Arnoldi_CONTENT(DEE)->LAPACK_rwork = (sunrealtype*)malloc(2 * N * sizeof(sunrealtype));
  SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_rwork, SUN_ERR_MALLOC_FAIL);

  xgeev_f77(&jobvl, &jobvr, &N, Arnoldi_CONTENT(DEE)->LAPACK_A, &lda,
            Arnoldi_CONTENT(DEE)->LAPACK_wr,
            NULL, &ldvl, NULL, &ldvr,
            &work, &lwork,
            Arnoldi_CONTENT(DEE)->LAPACK_rwork,
            &info);
#else
  xgeev_f77(&jobvl, &jobvr, &N, Arnoldi_CONTENT(DEE)->LAPACK_A, &lda,
            Arnoldi_CONTENT(DEE)->LAPACK_wr,
            Arnoldi_CONTENT(DEE)->LAPACK_wi,
            NULL, &ldvl, NULL, &ldvr,
            &work, &lwork,
            &info);
#endif

  if (info != 0) { return SUN_ERR_EXT_FAIL; }

  /* The workspace size is returned as the first entry of the work array */
  Arnoldi_CONTENT(DEE)->LAPACK_lwork = (sunindextype)work;

  Arnoldi_CONTENT(DEE)->LAPACK_work = (sunscalartype*)malloc(
    Arnoldi_CONTENT(DEE)->LAPACK_lwork * sizeof(sunscalartype));
  SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_work, SUN_ERR_MALLOC_FAIL);

  /* LAPACK array */
  Arnoldi_CONTENT(DEE)->LAPACK_arr =
    (sunscalartype**)malloc(Arnoldi_CONTENT(DEE)->kry_dim * sizeof(sunscalartype*));
  SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_arr, SUN_ERR_MALLOC_FAIL);

  for (int k = 0; k < Arnoldi_CONTENT(DEE)->kry_dim; k++)
  {
    Arnoldi_CONTENT(DEE)->LAPACK_arr[k] =
      (sunscalartype*)malloc(2 * sizeof(sunscalartype));
    SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_arr[k], SUN_ERR_MALLOC_FAIL);
  }

  /* Hessenberg matrix Hes */
  Arnoldi_CONTENT(DEE)->Hes = (sunscalartype**)malloc(
    (Arnoldi_CONTENT(DEE)->kry_dim + 1) * sizeof(sunscalartype*));
  SUNAssert(Arnoldi_CONTENT(DEE)->Hes, SUN_ERR_MALLOC_FAIL);

  for (int k = 0; k <= Arnoldi_CONTENT(DEE)->kry_dim; k++)
  {
    Arnoldi_CONTENT(DEE)->Hes[k] =
      (sunscalartype*)malloc(Arnoldi_CONTENT(DEE)->kry_dim * sizeof(sunscalartype));
    SUNAssert(Arnoldi_CONTENT(DEE)->Hes[k], SUN_ERR_MALLOC_FAIL);
  }

  /* Initialize the vector V */
  sunscalartype qdotq;
  SUNCheckCall(N_VDotProdComplex(Arnoldi_CONTENT(DEE)->q, Arnoldi_CONTENT(DEE)->q, &qdotq));
  sunrealtype normq = SUNRsqrt(SUN_REAL(qdotq));

  N_VScale(SUN_CCONST(1.0, 0.0) / normq, Arnoldi_CONTENT(DEE)->q, Arnoldi_CONTENT(DEE)->V[0]);
  SUNCheckLastErr();

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_SetNumPreprocessIters_Arnoldi(SUNDomEigEstimator DEE,
                                                            int num_iters)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  /* Check if num_iters >= 0 */
  if (num_iters < 0) { num_iters = DEE_NUM_OF_WARMUPS_ARNOLDI_DEFAULT; }

  /* set the number of warmups */
  Arnoldi_CONTENT(DEE)->num_warmups = num_iters;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_SetInitialGuess_Arnoldi(SUNDomEigEstimator DEE,
                                                      N_Vector q)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(q, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  sunscalartype qdotq;
  SUNCheckCall(N_VDotProdComplex(q, q, &qdotq));
  sunrealtype normq = SUNRsqrt(SUN_REAL(qdotq));

  /* set the initial guess */
  N_VScale(SUN_CCONST(1.0, 0.0) / normq, q, Arnoldi_CONTENT(DEE)->V[0]);
  SUNCheckLastErr();

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_Estimate_Arnoldi(SUNDomEigEstimator DEE,
                                               sunrealtype* lambdaR,
                                               sunrealtype* lambdaI)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(lambdaR, SUN_ERR_ARG_CORRUPT);
  SUNAssert(lambdaI, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE)->ATimes, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE)->Hes, SUN_ERR_ARG_CORRUPT);

  int retval;
  sunindextype n = Arnoldi_CONTENT(DEE)->kry_dim;
  sunrealtype normq;
  sunscalartype qdotq;
  Arnoldi_CONTENT(DEE)->num_ATimes = 0;
  Arnoldi_CONTENT(DEE)->num_iters  = 0;

  /* Set the initial q = A^{num_warmups}q/||A^{num_warmups}q|| */
  for (int i = 0; i < Arnoldi_CONTENT(DEE)->num_warmups; i++)
  {
    retval = Arnoldi_CONTENT(DEE)->ATimes(Arnoldi_CONTENT(DEE)->ATdata,
                                          Arnoldi_CONTENT(DEE)->V[0],
                                          Arnoldi_CONTENT(DEE)->q);
    Arnoldi_CONTENT(DEE)->num_ATimes++;
    Arnoldi_CONTENT(DEE)->num_iters++;
    if (retval != 0) { return SUN_ERR_USER_FCN_FAIL; }

    SUNCheckCall(N_VDotProdComplex(Arnoldi_CONTENT(DEE)->q, Arnoldi_CONTENT(DEE)->q, &qdotq));
    normq = SUNRsqrt(SUN_REAL(qdotq));

    N_VScale(SUN_CCONST(1.0, 0.0) / normq, Arnoldi_CONTENT(DEE)->q, Arnoldi_CONTENT(DEE)->V[0]);
    SUNCheckLastErr();
  }

  for (int i = 0; i < n; i++)
  {
    /* Compute the next Krylov vector */
    retval = Arnoldi_CONTENT(DEE)->ATimes(Arnoldi_CONTENT(DEE)->ATdata,
                                          Arnoldi_CONTENT(DEE)->V[i],
                                          Arnoldi_CONTENT(DEE)->V[i + 1]);
    Arnoldi_CONTENT(DEE)->num_ATimes++;
    Arnoldi_CONTENT(DEE)->num_iters++;
    if (retval != 0) { return SUN_ERR_USER_FCN_FAIL; }

    sunrealtype new_vk_norm;
    SUNCheckCall(SUNModifiedGS(Arnoldi_CONTENT(DEE)->V,
                               Arnoldi_CONTENT(DEE)->Hes, i + 1, (int)n,
                               &new_vk_norm));
    Arnoldi_CONTENT(DEE)->Hes[i + 1][i] = (sunscalartype)new_vk_norm;

    /* Unitize the computed orthogonal vector */
    N_VScale(SUN_CCONST(1.0, 0.0) / Arnoldi_CONTENT(DEE)->Hes[i + 1][i],
             Arnoldi_CONTENT(DEE)->V[i + 1], Arnoldi_CONTENT(DEE)->V[i + 1]);
    SUNCheckLastErr();
  }

  /* Pack the Hessenberg matrix in column-major order for LAPACK dgeev_ call */
  int k = 0;
  for (int j = 0; j < n; j++)
  {
    for (int i = 0; i < n; i++)
    {
      Arnoldi_CONTENT(DEE)->LAPACK_A[k] = Arnoldi_CONTENT(DEE)->Hes[i][j];
      k++;
    }
  }

  char jobvl = 'N'; // Do not compute left eigenvectors
  char jobvr = 'N'; // Do not compute right eigenvectors

  /* Call LAPACK's dgeev function
      return info values refer to
    = 0:  successful exit
    < 0:  if info = -i, the i-th argument had an illegal value.
    > 0:  if info = i, the QR algorithm failed to compute all the
          eigenvalues, and no eigenvectors have been computed;
          elements i+1:N of LAPACK_wr and LAPACK_wi contain
          eigenvalues which have converged.
  */
  sunindextype lda  = n;
  sunindextype ldvl = n;
  sunindextype ldvr = n;
  sunindextype info;
  sunindextype lwork = Arnoldi_CONTENT(DEE)->LAPACK_lwork;
  sunrealtype* rwork = Arnoldi_CONTENT(DEE)->LAPACK_rwork;
#if defined(SUNDIALS_SCALAR_TYPE_COMPLEX)
  xgeev_f77(&jobvl, &jobvr, &n, Arnoldi_CONTENT(DEE)->LAPACK_A, &lda,
            Arnoldi_CONTENT(DEE)->LAPACK_wr,
            NULL, &ldvl, NULL, &ldvr, Arnoldi_CONTENT(DEE)->LAPACK_work,
            &lwork, rwork, &info);
#else
  xgeev_f77(&jobvl, &jobvr, &n, Arnoldi_CONTENT(DEE)->LAPACK_A, &lda,
            Arnoldi_CONTENT(DEE)->LAPACK_wr, Arnoldi_CONTENT(DEE)->LAPACK_wi,
            NULL, &ldvl, NULL, &ldvr, Arnoldi_CONTENT(DEE)->LAPACK_work, &lwork,
            &info);
#endif

  if (info != 0) { return SUN_ERR_EXT_FAIL; }

  /* order the eigenvalues by their magnitude */
#if defined(SUNDIALS_SCALAR_TYPE_COMPLEX)
  sundomeigest_SortEigenvaluesByMagnitude(Arnoldi_CONTENT(DEE)->LAPACK_wr, n);
#else
  sundomeigest_SortEigenvaluesByMagnitude(Arnoldi_CONTENT(DEE)->LAPACK_wr,
                                          Arnoldi_CONTENT(DEE)->LAPACK_wi, n);
#endif

  /* Copy the dominant eigenvalue */
#if defined(SUNDIALS_SCALAR_TYPE_COMPLEX)
  *lambdaR = SUN_REAL(Arnoldi_CONTENT(DEE)->LAPACK_wr[0]);
  *lambdaI = SUN_IMAG(Arnoldi_CONTENT(DEE)->LAPACK_wr[0]);
#else
  *lambdaR = Arnoldi_CONTENT(DEE)->LAPACK_wr[0];
  *lambdaI = Arnoldi_CONTENT(DEE)->LAPACK_wi[0];
#endif

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_GetNumIters_Arnoldi(SUNDomEigEstimator DEE,
                                                  long int* num_iters)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(num_iters, SUN_ERR_ARG_CORRUPT);

  *num_iters = Arnoldi_CONTENT(DEE)->num_iters;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_GetNumATimesCalls_Arnoldi(SUNDomEigEstimator DEE,
                                                        long int* num_ATimes)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(num_ATimes, SUN_ERR_ARG_CORRUPT);

  *num_ATimes = Arnoldi_CONTENT(DEE)->num_ATimes;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_Write_Arnoldi(SUNDomEigEstimator DEE, FILE* outfile)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(outfile, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  fprintf(outfile, "\nArnoldi Iteration SUNDomEigEstimator:\n");
  fprintf(outfile, "Krylov dimension         = %d\n",
          Arnoldi_CONTENT(DEE)->kry_dim);
  fprintf(outfile, "Num. preprocessing iters = %d\n",
          Arnoldi_CONTENT(DEE)->num_warmups);
  fprintf(outfile, "Num. iters               = %ld\n",
          Arnoldi_CONTENT(DEE)->num_iters);
  fprintf(outfile, "Num. ATimes calls        = %ld\n\n",
          Arnoldi_CONTENT(DEE)->num_ATimes);

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_Destroy_Arnoldi(SUNDomEigEstimator* DEEptr)
{
  if ((*DEEptr) == NULL) { return SUN_SUCCESS; }

  SUNDomEigEstimator DEE = *DEEptr;

  if (DEE->content)
  {
    /* delete items from within the content structure */
    if (Arnoldi_CONTENT(DEE)->q)
    {
      N_VDestroy(Arnoldi_CONTENT(DEE)->q);
      Arnoldi_CONTENT(DEE)->q = NULL;
    }
    if (Arnoldi_CONTENT(DEE)->V)
    {
      N_VDestroyVectorArray(Arnoldi_CONTENT(DEE)->V,
                            Arnoldi_CONTENT(DEE)->kry_dim + 1);
      Arnoldi_CONTENT(DEE)->V = NULL;
    }
    if (Arnoldi_CONTENT(DEE)->LAPACK_A != NULL)
    {
      free(Arnoldi_CONTENT(DEE)->LAPACK_A);
      Arnoldi_CONTENT(DEE)->LAPACK_A = NULL;
    }
    if (Arnoldi_CONTENT(DEE)->LAPACK_wr != NULL)
    {
      free(Arnoldi_CONTENT(DEE)->LAPACK_wr);
      Arnoldi_CONTENT(DEE)->LAPACK_wr = NULL;
    }
    if (Arnoldi_CONTENT(DEE)->LAPACK_wi != NULL)
    {
      free(Arnoldi_CONTENT(DEE)->LAPACK_wi);
      Arnoldi_CONTENT(DEE)->LAPACK_wi = NULL;
    }
    if (Arnoldi_CONTENT(DEE)->LAPACK_work != NULL)
    {
      free(Arnoldi_CONTENT(DEE)->LAPACK_work);
      Arnoldi_CONTENT(DEE)->LAPACK_work = NULL;
    }

    /* free LAPACK_arr */
    if (Arnoldi_CONTENT(DEE)->LAPACK_arr != NULL)
    {
      for (int k = 0; k < Arnoldi_CONTENT(DEE)->kry_dim; k++)
      {
        free(Arnoldi_CONTENT(DEE)->LAPACK_arr[k]);
        Arnoldi_CONTENT(DEE)->LAPACK_arr[k] = NULL;
      }
      free(Arnoldi_CONTENT(DEE)->LAPACK_arr);
      Arnoldi_CONTENT(DEE)->LAPACK_arr = NULL;
    }
    /* free Hes */
    if (Arnoldi_CONTENT(DEE)->Hes != NULL)
    {
      for (int k = 0; k <= Arnoldi_CONTENT(DEE)->kry_dim; k++)
      {
        free(Arnoldi_CONTENT(DEE)->Hes[k]);
        Arnoldi_CONTENT(DEE)->Hes[k] = NULL;
      }
      free(Arnoldi_CONTENT(DEE)->Hes);
      Arnoldi_CONTENT(DEE)->Hes = NULL;
    }

    free(DEE->content);
    DEE->content = NULL;
  }
  if (DEE->ops)
  {
    free(DEE->ops);
    DEE->ops = NULL;
  }
  free(DEE);
  *DEEptr = NULL;
  return SUN_SUCCESS;
}

#if defined(SUNDIALS_SCALAR_TYPE_COMPLEX)
static void sundomeigest_SortEigenvaluesByMagnitude(sunscalartype* wr, int n)
#else
static void sundomeigest_SortEigenvaluesByMagnitude(sunrealtype* wr,
                                                    sunrealtype* wi, int n)
#endif
{
  for (int i = 0; i < n - 1; i++)
  {
    int imax = i;

#if defined(SUNDIALS_SCALAR_TYPE_COMPLEX)
    sunrealtype re_max = SUN_REAL(wr[i]);
    sunrealtype im_max = SUN_IMAG(wr[i]);
#else
    sunrealtype re_max = wr[i];
    sunrealtype im_max = wi[i];
#endif

    sunrealtype max_mag = re_max * re_max + im_max * im_max;

    for (int j = i + 1; j < n; j++)
    {
#if defined(SUNDIALS_SCALAR_TYPE_COMPLEX)
      sunrealtype re = SUN_REAL(wr[j]);
      sunrealtype im = SUN_IMAG(wr[j]);
#else
      sunrealtype re = wr[j];
      sunrealtype im = wi[j];
#endif

      sunrealtype mag = re * re + im * im;

      if (mag > max_mag)
      {
        max_mag = mag;
        imax    = j;
      }
    }

    if (imax != i)
    {
#if defined(SUNDIALS_SCALAR_TYPE_COMPLEX)
      sunscalartype tmp = wr[i];
      wr[i]             = wr[imax];
      wr[imax]          = tmp;
#else
      sunrealtype tmp_r = wr[i];
      sunrealtype tmp_i = wi[i];

      wr[i] = wr[imax];
      wi[i] = wi[imax];

      wr[imax] = tmp_r;
      wi[imax] = tmp_i;
#endif
    }
  }
}