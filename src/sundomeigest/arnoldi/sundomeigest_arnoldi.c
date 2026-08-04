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

/* Interfaces to match 'sunrealtype' with the correct LAPACK functions */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define xgeev_f77 dgeev_f77
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define xgeev_f77 sgeev_f77
#else
#error Incompatible sunrealtype for LAPACK; disable LAPACK and rebuild
#endif

#ifndef MAX_DQITERS
#define MAX_DQITERS 3
#endif

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/* Default estimator parameters */
#define DEE_NUM_OF_WARMUPS_ARNOLDI_DEFAULT 100
#define DEE_TOL_OF_WARMUPS_ARNOLDI_DEFAULT SUN_RCONST(0.005)

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

SUNErrCode dee_DQJtimes_Arnoldi(void* voidstarDEE, N_Vector v, N_Vector Jv);

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
  SUNAssertNull(N_VDotProd(q, q) > SUN_SMALL_REAL, SUN_ERR_ARG_INCOMPATIBLE);

  /* Create dominant eigenvalue estimator */
  DEE = NULL;
  DEE = SUNDomEigEstimator_NewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  DEE->ops->setatimes = SUNDomEigEstimator_SetATimes_Arnoldi;
  DEE->ops->setrhs    = SUNDomEigEstimator_SetRhs_Arnoldi;
  DEE->ops->setrhslinearizationpoint =
    SUNDomEigEstimator_SetRhsLinearizationPoint_Arnoldi;
  DEE->ops->setnumpreprocessiters =
    SUNDomEigEstimator_SetNumPreprocessIters_Arnoldi;
  DEE->ops->setreltol         = SUNDomEigEstimator_SetRelTol_Arnoldi;
  DEE->ops->setinitialguess   = SUNDomEigEstimator_SetInitialGuess_Arnoldi;
  DEE->ops->initialize        = SUNDomEigEstimator_Initialize_Arnoldi;
  DEE->ops->estimate          = SUNDomEigEstimator_Estimate_Arnoldi;
  DEE->ops->getnumiters       = SUNDomEigEstimator_GetNumIters_Arnoldi;
  DEE->ops->getnumrhsevals    = SUNDomEigEstimator_GetNumRhsEvals_Arnoldi;
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
  content->ATimes        = NULL;
  content->ATdata        = NULL;
  content->V             = NULL;
  content->q             = NULL;
  content->rhs_linY      = NULL;
  content->rhs_linT      = ZERO;
  content->Fy            = NULL;
  content->work          = NULL;
  content->kry_dim       = kry_dim;
  content->num_warmups   = DEE_NUM_OF_WARMUPS_ARNOLDI_DEFAULT;
  content->num_iters     = 0;
  content->num_ATimes    = 0;
  content->warmup_to_tol = SUNFALSE;
  content->tol_warmup    = DEE_TOL_OF_WARMUPS_ARNOLDI_DEFAULT;
  content->rhsfn         = NULL;
  content->rhs_data      = NULL;
  content->nfevals       = 0;
  content->LAPACK_A      = NULL;
  content->LAPACK_wr     = NULL;
  content->LAPACK_wi     = NULL;
  content->LAPACK_work   = NULL;
  content->LAPACK_lwork  = 0;
  content->LAPACK_arr    = NULL;
  content->Hes           = NULL;

  /* Allocate content */
  content->q = N_VClone(q);
  SUNCheckLastErrNull();

  content->V = N_VCloneVectorArray(kry_dim + 1, q);
  SUNCheckLastErrNull();

  /* Initialize the vector V[0] */
  sunrealtype normq = N_VDotProd(q, q);
  SUNCheckLastErrNull();

  normq = SUNRsqrt(normq);

  N_VScale(ONE / normq, q, content->V[0]);
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

SUNErrCode SUNDomEigEstimator_SetRhs_Arnoldi(SUNDomEigEstimator DEE,
                                             void* rhs_data, SUNRhsFn RHSfn)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(RHSfn, SUN_ERR_ARG_CORRUPT);

  /* set function pointers to integrator-supplied RHS routine
     and data, and return with success */
  Arnoldi_CONTENT(DEE)->rhsfn    = RHSfn;
  Arnoldi_CONTENT(DEE)->rhs_data = rhs_data;

  DEE->ops->setatimes(DEE, (void*)DEE, dee_DQJtimes_Arnoldi);
  SUNCheckLastErr();

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_SetRhsLinearizationPoint_Arnoldi(
  SUNDomEigEstimator DEE, sunrealtype t, N_Vector y)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(y, SUN_ERR_ARG_CORRUPT);

  if (Arnoldi_CONTENT(DEE)->rhs_linY == NULL)
  {
    Arnoldi_CONTENT(DEE)->rhs_linY = N_VClone(y);
    SUNCheckLastErr();
  }

  Arnoldi_CONTENT(DEE)->rhs_linT = t;

  N_VScale(ONE, y, Arnoldi_CONTENT(DEE)->rhs_linY);
  SUNCheckLastErr();

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
    Arnoldi_CONTENT(DEE)->LAPACK_A = (sunrealtype*)malloc(
      (Arnoldi_CONTENT(DEE)->kry_dim * Arnoldi_CONTENT(DEE)->kry_dim) *
      sizeof(sunrealtype));
    SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_A, SUN_ERR_MALLOC_FAIL);
  }
  if (Arnoldi_CONTENT(DEE)->LAPACK_wr == NULL)
  {
    Arnoldi_CONTENT(DEE)->LAPACK_wr =
      malloc(Arnoldi_CONTENT(DEE)->kry_dim * sizeof(sunrealtype));
    SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_wr, SUN_ERR_MALLOC_FAIL);
  }
  if (Arnoldi_CONTENT(DEE)->LAPACK_wi == NULL)
  {
    Arnoldi_CONTENT(DEE)->LAPACK_wi =
      malloc(Arnoldi_CONTENT(DEE)->kry_dim * sizeof(sunrealtype));
    SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_wi, SUN_ERR_MALLOC_FAIL);
  }

  /* query the workspace size (call with lwork = -1) */
  char jobvl         = 'N';
  char jobvr         = 'N';
  sunindextype N     = Arnoldi_CONTENT(DEE)->kry_dim;
  sunindextype lda   = Arnoldi_CONTENT(DEE)->kry_dim;
  sunindextype ldvl  = 1;
  sunindextype ldvr  = 1;
  sunindextype info  = 0;
  sunindextype lwork = -1;
  sunrealtype work   = SUN_RCONST(0.0);

  xgeev_f77(&jobvl, &jobvr, &N, Arnoldi_CONTENT(DEE)->LAPACK_A, &lda,
            Arnoldi_CONTENT(DEE)->LAPACK_wr, Arnoldi_CONTENT(DEE)->LAPACK_wi,
            NULL, &ldvl, NULL, &ldvr, &work, &lwork, &info);

  /* The workspace size is returned as the first entry of the work array */
  Arnoldi_CONTENT(DEE)->LAPACK_lwork = (sunindextype)work;

  Arnoldi_CONTENT(DEE)->LAPACK_work = (sunrealtype*)malloc(
    Arnoldi_CONTENT(DEE)->LAPACK_lwork * sizeof(sunrealtype));
  SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_work, SUN_ERR_MALLOC_FAIL);

  /* LAPACK array */
  Arnoldi_CONTENT(DEE)->LAPACK_arr =
    (sunrealtype**)malloc(Arnoldi_CONTENT(DEE)->kry_dim * sizeof(sunrealtype*));
  SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_arr, SUN_ERR_MALLOC_FAIL);

  for (int k = 0; k < Arnoldi_CONTENT(DEE)->kry_dim; k++)
  {
    Arnoldi_CONTENT(DEE)->LAPACK_arr[k] =
      (sunrealtype*)malloc(2 * sizeof(sunrealtype));
    SUNAssert(Arnoldi_CONTENT(DEE)->LAPACK_arr[k], SUN_ERR_MALLOC_FAIL);
  }

  /* Hessenberg matrix Hes */
  Arnoldi_CONTENT(DEE)->Hes = (sunrealtype**)malloc(
    (Arnoldi_CONTENT(DEE)->kry_dim + 1) * sizeof(sunrealtype*));
  SUNAssert(Arnoldi_CONTENT(DEE)->Hes, SUN_ERR_MALLOC_FAIL);

  for (int k = 0; k <= Arnoldi_CONTENT(DEE)->kry_dim; k++)
  {
    Arnoldi_CONTENT(DEE)->Hes[k] =
      (sunrealtype*)malloc(Arnoldi_CONTENT(DEE)->kry_dim * sizeof(sunrealtype));
    SUNAssert(Arnoldi_CONTENT(DEE)->Hes[k], SUN_ERR_MALLOC_FAIL);
  }

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

SUNErrCode SUNDomEigEstimator_SetRelTol_Arnoldi(SUNDomEigEstimator DEE,
                                                sunrealtype tol)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  /* set the tolerance for preprocessing iterations */
  if (tol < ZERO)
  {
    Arnoldi_CONTENT(DEE)->warmup_to_tol = SUNFALSE;
    return SUN_SUCCESS;
  }
  else if (tol == ZERO || tol > ONE - SUN_UNIT_ROUNDOFF)
  {
    tol = DEE_TOL_OF_WARMUPS_ARNOLDI_DEFAULT;
  }
  Arnoldi_CONTENT(DEE)->tol_warmup = tol;

  /* set the type of warmup iterations */
  Arnoldi_CONTENT(DEE)->warmup_to_tol = SUNTRUE;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_SetInitialGuess_Arnoldi(SUNDomEigEstimator DEE,
                                                      N_Vector q)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(q, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  sunrealtype normq = N_VDotProd(q, q);
  SUNCheckLastErr();

  normq = SUNRsqrt(normq);

  /* set the initial guess */
  N_VScale(ONE / normq, q, Arnoldi_CONTENT(DEE)->V[0]);
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
  Arnoldi_CONTENT(DEE)->num_ATimes = 0;
  Arnoldi_CONTENT(DEE)->num_iters  = 0;

  sunrealtype res;
  sunrealtype new_lambda = ZERO;
  sunrealtype old_lambda = ZERO;

  /* Set the initial q = A^{num_warmups}q/||A^{num_warmups}q|| */
  for (int i = 0; i < Arnoldi_CONTENT(DEE)->num_warmups; i++)
  {
    retval = Arnoldi_CONTENT(DEE)->ATimes(Arnoldi_CONTENT(DEE)->ATdata,
                                          Arnoldi_CONTENT(DEE)->V[0],
                                          Arnoldi_CONTENT(DEE)->q);
    Arnoldi_CONTENT(DEE)->num_ATimes++;
    Arnoldi_CONTENT(DEE)->num_iters++;
    if (retval != 0) { return SUN_ERR_USER_FCN_FAIL; }

    if (Arnoldi_CONTENT(DEE)->warmup_to_tol)
    {
      new_lambda = N_VDotProd(Arnoldi_CONTENT(DEE)->V[0],
                              Arnoldi_CONTENT(DEE)->q); //Rayleigh quotient
      SUNCheckLastErr();
    }

    normq = N_VDotProd(Arnoldi_CONTENT(DEE)->q, Arnoldi_CONTENT(DEE)->q);
    SUNCheckLastErr();

    normq = SUNRsqrt(normq);
    N_VScale(ONE / normq, Arnoldi_CONTENT(DEE)->q, Arnoldi_CONTENT(DEE)->V[0]);
    SUNCheckLastErr();

    if (Arnoldi_CONTENT(DEE)->warmup_to_tol)
    {
      res        = SUNRabs(new_lambda - old_lambda);
      old_lambda = new_lambda;
      if (res <= Arnoldi_CONTENT(DEE)->tol_warmup * SUNRabs(new_lambda))
      {
        break;
      }
    }
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

    SUNCheckCall(SUNModifiedGS(Arnoldi_CONTENT(DEE)->V,
                               Arnoldi_CONTENT(DEE)->Hes, i + 1, (int)n,
                               &(Arnoldi_CONTENT(DEE)->Hes[i + 1][i])));

    /* Unitize the computed orthogonal vector */
    N_VScale(SUN_RCONST(1.0) / Arnoldi_CONTENT(DEE)->Hes[i + 1][i],
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
  xgeev_f77(&jobvl, &jobvr, &n, Arnoldi_CONTENT(DEE)->LAPACK_A, &lda,
            Arnoldi_CONTENT(DEE)->LAPACK_wr, Arnoldi_CONTENT(DEE)->LAPACK_wi,
            NULL, &ldvl, NULL, &ldvr, Arnoldi_CONTENT(DEE)->LAPACK_work, &lwork,
            &info);

  if (info != 0) { return SUN_ERR_EXT_FAIL; }

  /* order the eigenvalues by their magnitude */
  for (int i = 0; i < n; i++)
  {
    Arnoldi_CONTENT(DEE)->LAPACK_arr[i][0] = Arnoldi_CONTENT(DEE)->LAPACK_wr[i];
    Arnoldi_CONTENT(DEE)->LAPACK_arr[i][1] = Arnoldi_CONTENT(DEE)->LAPACK_wi[i];
  }

  /* Sort the array using qsort */
  qsort(Arnoldi_CONTENT(DEE)->LAPACK_arr, n,
        sizeof(Arnoldi_CONTENT(DEE)->LAPACK_arr[0]), sundomeigest_Compare);

  /* Substitute the ordered eigenvalues back in LAPACK_w* */
  for (int i = 0; i < n; i++)
  {
    Arnoldi_CONTENT(DEE)->LAPACK_wr[i] = Arnoldi_CONTENT(DEE)->LAPACK_arr[i][0];
    Arnoldi_CONTENT(DEE)->LAPACK_wi[i] = Arnoldi_CONTENT(DEE)->LAPACK_arr[i][1];
  }

  /* Copy the dominant eigenvalue */
  *lambdaR = Arnoldi_CONTENT(DEE)->LAPACK_wr[0];
  *lambdaI = Arnoldi_CONTENT(DEE)->LAPACK_wi[0];

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

SUNErrCode SUNDomEigEstimator_GetNumRhsEvals_Arnoldi(SUNDomEigEstimator DEE,
                                                     long int* num_rhs_evals)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(num_rhs_evals, SUN_ERR_ARG_CORRUPT);

  *num_rhs_evals = Arnoldi_CONTENT(DEE)->nfevals;

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
    if (Arnoldi_CONTENT(DEE)->rhs_linY)
    {
      N_VDestroy(Arnoldi_CONTENT(DEE)->rhs_linY);
      Arnoldi_CONTENT(DEE)->rhs_linY = NULL;
    }
    if (Arnoldi_CONTENT(DEE)->Fy)
    {
      N_VDestroy(Arnoldi_CONTENT(DEE)->Fy);
      Arnoldi_CONTENT(DEE)->Fy = NULL;
    }
    if (Arnoldi_CONTENT(DEE)->work)
    {
      N_VDestroy(Arnoldi_CONTENT(DEE)->work);
      Arnoldi_CONTENT(DEE)->work = NULL;
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

// Comparison function for qsort
int sundomeigest_Compare(const void* a, const void* b)
{
  const sunrealtype* cplx_a = *(const sunrealtype* const*)a;
  const sunrealtype* cplx_b = *(const sunrealtype* const*)b;

  sunrealtype mag_a = SUNRsqrt(cplx_a[0] * cplx_a[0] + cplx_a[1] * cplx_a[1]);
  sunrealtype mag_b = SUNRsqrt(cplx_b[0] * cplx_b[0] + cplx_b[1] * cplx_b[1]);
  return (mag_b > mag_a) - (mag_b < mag_a); // Descending order
}

/*---------------------------------------------------------------
  dee_DQJtimes_Arnoldi:

  This routine generates a difference quotient approximation to
  the Jacobian-vector product f_y(t,y) * v. The approximation is
  Jv = [f(y + v*sig) - f(y)]/sig, where
      sig = sign(y^T v) * sqrt(unit roundoff)
            * max(|y^T v|, ||v||_1) / (v^T v).
  ---------------------------------------------------------------*/
SUNErrCode dee_DQJtimes_Arnoldi(void* voidstarDEE, N_Vector v, N_Vector Jv)
{
  SUNDomEigEstimator DEE = (SUNDomEigEstimator)voidstarDEE;
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(v, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Jv, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE)->rhsfn, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Arnoldi_CONTENT(DEE)->rhs_linY, SUN_ERR_ARG_CORRUPT);

  sunrealtype vdotv = N_VDotProd(v, v);
  if (vdotv <= SUN_SMALL_REAL)
  {
    N_VScale(ZERO, v, Jv);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  if (Arnoldi_CONTENT(DEE)->work == NULL)
  {
    Arnoldi_CONTENT(DEE)->work = N_VClone(v);
    SUNCheckLastErr();
  }
  if (Arnoldi_CONTENT(DEE)->Fy == NULL)
  {
    Arnoldi_CONTENT(DEE)->Fy = N_VClone(v);
    SUNCheckLastErr();
  }

  sunrealtype sig, siginv;
  int iter, retval;

  N_Vector y    = Arnoldi_CONTENT(DEE)->rhs_linY;
  N_Vector work = Arnoldi_CONTENT(DEE)->work;
  N_Vector Fy   = Arnoldi_CONTENT(DEE)->Fy;

  retval = Arnoldi_CONTENT(DEE)->rhsfn(Arnoldi_CONTENT(DEE)->rhs_linT, y, Fy,
                                       Arnoldi_CONTENT(DEE)->rhs_data);
  Arnoldi_CONTENT(DEE)->nfevals++;
  if (retval != 0) { return SUN_ERR_USER_FCN_FAIL; }

  /* Initialize perturbation */
  sunrealtype ydotv   = N_VDotProd(y, v);
  sunrealtype sq1norm = N_VL1Norm(v);
  sunrealtype sign    = (ydotv >= ZERO) ? ONE : -ONE;
  sunrealtype sqrteps = SUNRsqrt(SUN_UNIT_ROUNDOFF);
  sig = sign * sqrteps * SUNMAX(SUNRabs(ydotv), sq1norm) / vdotv;

  for (iter = 0; iter < MAX_DQITERS; iter++)
  {
    /* Set work = y + sig*v */
    N_VLinearSum(sig, v, ONE, y, work);

    /* Set Jv = f(tn, y+sig*v) */
    retval = Arnoldi_CONTENT(DEE)->rhsfn(Arnoldi_CONTENT(DEE)->rhs_linT, work,
                                         Jv, Arnoldi_CONTENT(DEE)->rhs_data);
    Arnoldi_CONTENT(DEE)->nfevals++;
    if (retval == 0) { break; }
    if (retval < 0) { return SUN_ERR_USER_FCN_FAIL; }

    /* If f failed recoverably, shrink sig and retry */
    sig *= SUN_RCONST(0.25);
  }

  /* If retval still isn't 0, return with a recoverable failure */
  if (retval > 0) { return (+1); }

  /* Replace Jv by (Jv - fn)/sig */
  siginv = ONE / sig;
  N_VLinearSum(siginv, Jv, -siginv, Fy, Jv);

  return SUN_SUCCESS;
}