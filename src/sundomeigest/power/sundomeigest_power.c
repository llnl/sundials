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
 * This is the implementation file for the Power Iteration (PI)
 * implementation of the SUNDomEigEst package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundomeigest/sundomeigest_power.h>

#include "sundials_logger_impl.h"
#include "sundials_macros.h"

#ifndef MAX_DQITERS
#define MAX_DQITERS 3
#endif

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/* Default estimator parameters */
#define DEE_NUM_OF_WARMUPS_PI_DEFAULT 100

/* Default Power Iteration parameters */
#define DEE_TOL_DEFAULT      SUN_RCONST(0.005)
#define DEE_MAX_ITER_DEFAULT 100

/*
 * -----------------------------------------------------------------
 * Power itetation structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define PI_CONTENT(DEE) ((SUNDomEigEstimatorContent_Power)(DEE->content))

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/*
* --------------------------------------------------------------------------
* private functions
* --------------------------------------------------------------------------
*/

SUNErrCode sundomeigestimator_complex_dom_eigs_from_PI(
  SUNDomEigEstimator DEE, sunrealtype lambdaR, sunrealtype h21, N_Vector v_prev,
  N_Vector v, sunrealtype* lambdaR_out, sunrealtype* lambdaI_out);

SUNErrCode dee_DQJtimes_Power(void* voidstarDEE, N_Vector v, N_Vector Jv);

/* ----------------------------------------------------------------------------
 * Function to create a new PI estimator
 */

SUNDomEigEstimator SUNDomEigEstimator_Power(N_Vector q, long int max_iters,
                                            sunrealtype rel_tol,
                                            SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNDomEigEstimator DEE;
  SUNDomEigEstimatorContent_Power content;

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

  /* check for max_iters values; if illegal use defaults */
  if (max_iters <= 0) { max_iters = DEE_MAX_ITER_DEFAULT; }

  /* Check if rel_tol > 0 and < 1 */
  if (rel_tol < SUN_SMALL_REAL || rel_tol > ONE - SUN_UNIT_ROUNDOFF)
  {
    rel_tol = DEE_TOL_DEFAULT;
  }

  /* Create dominant eigenvalue estimator */
  DEE = NULL;
  DEE = SUNDomEigEstimator_NewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  DEE->ops->setatimes = SUNDomEigEstimator_SetATimes_Power;
  DEE->ops->setrhs    = SUNDomEigEstimator_SetRhs_Power;
  DEE->ops->setrhslinearizationpoint =
    SUNDomEigEstimator_SetRhsLinearizationPoint_Power;
  DEE->ops->setmaxiters = SUNDomEigEstimator_SetMaxIters_Power;
  DEE->ops->setnumpreprocessiters = SUNDomEigEstimator_SetNumPreprocessIters_Power;
  DEE->ops->setreltol         = SUNDomEigEstimator_SetRelTol_Power;
  DEE->ops->setinitialguess   = SUNDomEigEstimator_SetInitialGuess_Power;
  DEE->ops->initialize        = SUNDomEigEstimator_Initialize_Power;
  DEE->ops->estimate          = SUNDomEigEstimator_Estimate_Power;
  DEE->ops->getres            = SUNDomEigEstimator_GetRes_Power;
  DEE->ops->getnumiters       = SUNDomEigEstimator_GetNumIters_Power;
  DEE->ops->getnumrhsevals    = SUNDomEigEstimator_GetNumRhsEvals_Power;
  DEE->ops->getnumatimescalls = SUNDomEigEstimator_GetNumATimesCalls_Power;
  DEE->ops->write             = SUNDomEigEstimator_Write_Power;
  DEE->ops->destroy           = SUNDomEigEstimator_Destroy_Power;

  /* Create content */
  content = NULL;
  content = (SUNDomEigEstimatorContent_Power)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content  */
  DEE->content = content;

  /* Fill content */
  content->ATimes      = NULL;
  content->ATdata      = NULL;
  content->V           = NULL;
  content->q           = NULL;
  content->q_prev      = NULL;
  content->rhs_linY    = NULL;
  content->rhs_linT    = ZERO;
  content->Fy          = NULL;
  content->work        = NULL;
  content->is_complex  = SUNTRUE;
  content->max_iters   = max_iters;
  content->num_warmups = DEE_NUM_OF_WARMUPS_PI_DEFAULT;
  content->rel_tol     = rel_tol;
  content->res         = ZERO;
  content->rhsfn       = NULL;
  content->rhs_data    = NULL;
  content->nfevals     = 0;
  content->num_iters   = 0;
  content->num_ATimes  = 0;

  /* Allocate content */
  content->q = N_VClone(q);
  SUNCheckLastErrNull();

  content->V = N_VClone(q);
  SUNCheckLastErrNull();

  /* Initialize the vector V */
  sunrealtype normq = N_VDotProd(q, q);
  SUNCheckLastErrNull();

  normq = SUNRsqrt(normq);

  N_VScale(ONE / normq, q, PI_CONTENT(DEE)->V);
  SUNCheckLastErrNull();

  return (DEE);
}

/*
 * -----------------------------------------------------------------
 * implementation of dominant eigenvalue estimator operations
 * -----------------------------------------------------------------
 */

SUNErrCode SUNDomEigEstimator_SetATimes_Power(SUNDomEigEstimator DEE,
                                              void* A_data, SUNATimesFn ATimes)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(ATimes, SUN_ERR_ARG_CORRUPT);

  /* set function pointers to integrator-supplied ATimes routine
     and data, and return with success */
  PI_CONTENT(DEE)->ATimes = ATimes;
  PI_CONTENT(DEE)->ATdata = A_data;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_SetRhs_Power(SUNDomEigEstimator DEE,
                                           void* rhs_data, SUNRhsFn RHSfn)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(RHSfn, SUN_ERR_ARG_CORRUPT);

  /* set function pointers to integrator-supplied RHS routine
     and data, and return with success */
  PI_CONTENT(DEE)->rhsfn    = RHSfn;
  PI_CONTENT(DEE)->rhs_data = rhs_data;

  DEE->ops->setatimes(DEE, (void*)DEE, dee_DQJtimes_Power);
  SUNCheckLastErr();

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_SetRhsLinearizationPoint_Power(
  SUNDomEigEstimator DEE, sunrealtype t, N_Vector y)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(y, SUN_ERR_ARG_CORRUPT);

  if (PI_CONTENT(DEE)->rhs_linY == NULL)
  {
    PI_CONTENT(DEE)->rhs_linY = N_VClone(y);
    SUNCheckLastErr();
  }

  PI_CONTENT(DEE)->rhs_linT = t;

  N_VScale(ONE, y, PI_CONTENT(DEE)->rhs_linY);
  SUNCheckLastErr();

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_SetIsReal_Power(SUNDomEigEstimator DEE,
                                              sunbooleantype real)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  /* set the complex flag to the opposite of the real flag */
  PI_CONTENT(DEE)->is_complex = !real;

  /* q_prev is allocated in SUNDomEigEstimator_Initialize_Power, which is expected to be 
  called after this routine. If the user calls this routine after initialization, we need 
  to free q_prev here. */
  if (!(PI_CONTENT(DEE)->is_complex) && PI_CONTENT(DEE)->q_prev)
  {
    N_VDestroy(PI_CONTENT(DEE)->q_prev);
    PI_CONTENT(DEE)->q_prev = NULL;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_Initialize_Power(SUNDomEigEstimator DEE)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  if (PI_CONTENT(DEE)->rel_tol < SUN_SMALL_REAL ||
      PI_CONTENT(DEE)->rel_tol > ONE - SUN_UNIT_ROUNDOFF)
  {
    PI_CONTENT(DEE)->rel_tol = DEE_TOL_DEFAULT;
  }
  if (PI_CONTENT(DEE)->num_warmups < 0)
  {
    PI_CONTENT(DEE)->num_warmups = DEE_NUM_OF_WARMUPS_PI_DEFAULT;
  }
  if (PI_CONTENT(DEE)->max_iters <= 0)
  {
    PI_CONTENT(DEE)->max_iters = DEE_MAX_ITER_DEFAULT;
  }

  SUNAssert(PI_CONTENT(DEE)->ATimes, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);

  if (PI_CONTENT(DEE)->is_complex)
  {
    SUNAssert(PI_CONTENT(DEE)->q_prev == NULL, SUN_ERR_ARG_CORRUPT);
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_SetNumPreprocessIters_Power(SUNDomEigEstimator DEE,
                                                          int num_iters)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  /* Check if num_iters >= 0 */
  if (num_iters < 0) { num_iters = DEE_NUM_OF_WARMUPS_PI_DEFAULT; }

  /* set the number of warmups */
  PI_CONTENT(DEE)->num_warmups = num_iters;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_SetRelTol_Power(SUNDomEigEstimator DEE,
                                              sunrealtype rel_tol)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  /* Check if rel_tol > 0 and < 1 */
  if (rel_tol < SUN_SMALL_REAL || rel_tol > ONE - SUN_UNIT_ROUNDOFF)
  {
    rel_tol = DEE_TOL_DEFAULT;
  }

  /* set the tolerance */
  PI_CONTENT(DEE)->rel_tol = rel_tol;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_SetMaxIters_Power(SUNDomEigEstimator DEE,
                                                long int max_iters)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  /* Check for legal number of iters */
  if (max_iters <= 0) { max_iters = DEE_MAX_ITER_DEFAULT; }

  /* Set max iters */
  PI_CONTENT(DEE)->max_iters = max_iters;
  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_SetInitialGuess_Power(SUNDomEigEstimator DEE,
                                                    N_Vector q)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(q, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  sunrealtype normq = N_VDotProd(q, q);
  SUNCheckLastErr();

  normq = SUNRsqrt(normq);

  /* set the initial guess */
  N_VScale(ONE / normq, q, PI_CONTENT(DEE)->V);
  SUNCheckLastErr();

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_Estimate_Power(SUNDomEigEstimator DEE,
                                             sunrealtype* lambdaR,
                                             sunrealtype* lambdaI)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(lambdaR, SUN_ERR_ARG_CORRUPT);
  SUNAssert(lambdaI, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->ATimes, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->V, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->q, SUN_ERR_ARG_CORRUPT);
  SUNAssert((PI_CONTENT(DEE)->max_iters >= 0), SUN_ERR_ARG_CORRUPT);
  if (PI_CONTENT(DEE)->is_complex && (PI_CONTENT(DEE)->q_prev == NULL))
  {
    /* allocate q_prev vector */
    PI_CONTENT(DEE)->q_prev = N_VClone(PI_CONTENT(DEE)->q);
    SUNCheckLastErr();
  }

  sunrealtype newlambdaR = ZERO;
  sunrealtype oldlambdaR = ZERO;

  int retval;
  sunbooleantype converged;
  sunrealtype normq;
  PI_CONTENT(DEE)->num_ATimes = 0;
  PI_CONTENT(DEE)->num_iters  = 0;

  /* Set the initial q = A^{num_warmups}q/||A^{num_warmups}q|| */
  for (int i = 0; i < PI_CONTENT(DEE)->num_warmups; i++)
  {
    retval = PI_CONTENT(DEE)->ATimes(PI_CONTENT(DEE)->ATdata,
                                     PI_CONTENT(DEE)->V, PI_CONTENT(DEE)->q);
    PI_CONTENT(DEE)->num_ATimes++;
    PI_CONTENT(DEE)->num_iters++;
    if (retval != 0) { return SUN_ERR_USER_FCN_FAIL; }

    normq = N_VDotProd(PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->q);
    SUNCheckLastErr();

    normq = SUNRsqrt(normq);
    N_VScale(ONE / normq, PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->V);
    SUNCheckLastErr();
  }

  for (int k = 0; k < PI_CONTENT(DEE)->max_iters; k++)
  {
    if (PI_CONTENT(DEE)->is_complex)
    {
      N_VScale(ONE, PI_CONTENT(DEE)->V, PI_CONTENT(DEE)->q_prev);
      SUNCheckLastErr();
    }

    retval = PI_CONTENT(DEE)->ATimes(PI_CONTENT(DEE)->ATdata,
                                     PI_CONTENT(DEE)->V, PI_CONTENT(DEE)->q);
    PI_CONTENT(DEE)->num_ATimes++;
    PI_CONTENT(DEE)->num_iters++;
    if (retval != 0) { return SUN_ERR_USER_FCN_FAIL; }

    newlambdaR = N_VDotProd(PI_CONTENT(DEE)->V,
                            PI_CONTENT(DEE)->q); //Rayleigh quotient
    SUNCheckLastErr();

    PI_CONTENT(DEE)->res = SUNRabs(newlambdaR - oldlambdaR);
    converged =
      (PI_CONTENT(DEE)->res <= PI_CONTENT(DEE)->rel_tol * SUNRabs(newlambdaR));

    if (converged && !PI_CONTENT(DEE)->is_complex) { break; }

    normq = N_VDotProd(PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->q);
    SUNCheckLastErr();

    normq = SUNRsqrt(normq);
    N_VScale(ONE / normq, PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->V);
    SUNCheckLastErr();

    if (converged) { break; }

    oldlambdaR = newlambdaR;
  }

  if (PI_CONTENT(DEE)->is_complex)
  {
    retval = sundomeigestimator_complex_dom_eigs_from_PI(DEE, newlambdaR, normq,
                                                         PI_CONTENT(DEE)->q_prev,
                                                         PI_CONTENT(DEE)->V,
                                                         lambdaR, lambdaI);
    if (retval != 0) { return SUN_ERR_USER_FCN_FAIL; }
  }
  else
  {
    *lambdaR = newlambdaR;
    *lambdaI = ZERO;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_GetRes_Power(SUNDomEigEstimator DEE,
                                           sunrealtype* res)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(res, SUN_ERR_ARG_CORRUPT);

  *res = PI_CONTENT(DEE)->res;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_GetNumIters_Power(SUNDomEigEstimator DEE,
                                                long int* num_iters)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(num_iters, SUN_ERR_ARG_CORRUPT);

  *num_iters = PI_CONTENT(DEE)->num_iters;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_GetNumRhsEvals_Power(SUNDomEigEstimator DEE,
                                                   long int* num_rhs_evals)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(num_rhs_evals, SUN_ERR_ARG_CORRUPT);

  *num_rhs_evals = PI_CONTENT(DEE)->nfevals;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_GetNumATimesCalls_Power(SUNDomEigEstimator DEE,
                                                      long int* num_ATimes)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(num_ATimes, SUN_ERR_ARG_CORRUPT);

  *num_ATimes = PI_CONTENT(DEE)->num_ATimes;

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_Write_Power(SUNDomEigEstimator DEE, FILE* outfile)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(outfile, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  if (DEE == NULL || outfile == NULL) { return SUN_ERR_ARG_CORRUPT; }

  fprintf(outfile, "\nPower Iteration SUNDomEigEstimator:\n");
  fprintf(outfile, "Max. iters               = %ld\n",
          PI_CONTENT(DEE)->max_iters);
  fprintf(outfile, "Num. preprocessing iters = %d\n",
          PI_CONTENT(DEE)->num_warmups);
  fprintf(outfile, "Relative tolerance       = " SUN_FORMAT_G "\n",
          PI_CONTENT(DEE)->rel_tol);
  fprintf(outfile, "Residual                 = " SUN_FORMAT_G "\n",
          PI_CONTENT(DEE)->res);
  fprintf(outfile, "Num. iters               = %ld\n",
          PI_CONTENT(DEE)->num_iters);
  fprintf(outfile, "Num. ATimes calls        = %ld\n\n",
          PI_CONTENT(DEE)->num_ATimes);

  return SUN_SUCCESS;
}

SUNErrCode sundomeigestimator_complex_dom_eigs_from_PI(
  SUNDomEigEstimator DEE, sunrealtype lambdaR, sunrealtype h21, N_Vector v_prev,
  N_Vector v, sunrealtype* lambdaR_out, sunrealtype* lambdaI_out)
{
  SUNFunctionBegin(DEE->sunctx);

  int retval;
  sunrealtype cos_qs, gram_det, det_G_inv, h11, h12, h22, p11, p12, p21, p22;
  /* The threshold for identifying real or complex DEE is experimentally 
  determined based on the relative tolerance PI_CONTENT(DEE)->rel_tol */
  sunrealtype gram_det_tol = SUNMAX(SUN_RCONST(10.0) * SUN_UNIT_ROUNDOFF,
                                    SUN_RCONST(10.0) * PI_CONTENT(DEE)->rel_tol);
  cos_qs                   = N_VDotProd(v_prev, v);
  SUNCheckLastErr();

  /* Safety against roundoff in dot product */
  if (cos_qs > ONE) { cos_qs = ONE; }
  if (cos_qs < -ONE) { cos_qs = -ONE; }

  /* Use Gram determinant as the near-dependence measure:
     G = [ [1, cos_qs], [cos_qs, 1] ], det(G) = 1 - cos_qs^2
     This assumes v_prev and v are normalized. */
  gram_det = ONE - cos_qs * cos_qs;

  if (gram_det <= gram_det_tol)
  {
    /* Dominant eigenvalue is real */
    *lambdaR_out = lambdaR;
    *lambdaI_out = ZERO;
    return SUN_SUCCESS;
  }
  else
  {
    det_G_inv = ONE / gram_det;

    /* Solve for G = [v_prev v]' * [v_prev v] and compute
       projected matrix P = G^{-1} * [v_prev v]' * A * [v_prev v] */

    h11 = lambdaR;

    retval = PI_CONTENT(DEE)->ATimes(PI_CONTENT(DEE)->ATdata, v,
                                     PI_CONTENT(DEE)->q);
    PI_CONTENT(DEE)->num_ATimes++;
    if (retval != 0) { return SUN_ERR_USER_FCN_FAIL; }

    h12 = N_VDotProd(v_prev, PI_CONTENT(DEE)->q);
    SUNCheckLastErr();
    h22 = N_VDotProd(v, PI_CONTENT(DEE)->q);
    SUNCheckLastErr();

    p11 = det_G_inv * (h11 - cos_qs * h21);
    p12 = det_G_inv * (h12 - cos_qs * h22);
    p21 = det_G_inv * (h21 - cos_qs * h11);
    p22 = det_G_inv * (h22 - cos_qs * h12);

    /* Compute eigenvalues of P */
    sunrealtype traceP  = p11 + p22;
    sunrealtype detP    = p11 * p22 - p12 * p21;
    sunrealtype discrim = traceP * traceP - SUN_RCONST(4.0) * detP;
    if (discrim >= ZERO)
    {
      /* Dominant eigenvalue is real */
      sunrealtype sqrt_discrim = SUNRsqrt(discrim);
      sunrealtype lam_plus     = (traceP + sqrt_discrim) / SUN_RCONST(2.0);
      sunrealtype lam_minus    = (traceP - sqrt_discrim) / SUN_RCONST(2.0);
      if (SUNRabs(lam_plus) >= SUNRabs(lam_minus)) { *lambdaR_out = lam_plus; }
      else { *lambdaR_out = lam_minus; }
      *lambdaI_out = ZERO;
      *lambdaI_out = ZERO;
    }
    else
    {
      /* Dominant eigenvalue is complex */
      *lambdaR_out = traceP / SUN_RCONST(2.0);
      *lambdaI_out = SUNRsqrt(-discrim) / SUN_RCONST(2.0);
    }
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_Destroy_Power(SUNDomEigEstimator* DEEptr)
{
  SUNFunctionBegin((*DEEptr)->sunctx);

  if ((*DEEptr) == NULL) { return SUN_SUCCESS; }

  SUNDomEigEstimator DEE = *DEEptr;

  if (DEE->content)
  {
    /* delete items from within the content structure */
    if (PI_CONTENT(DEE)->q)
    {
      N_VDestroy(PI_CONTENT(DEE)->q);
      PI_CONTENT(DEE)->q = NULL;
    }
    if (PI_CONTENT(DEE)->q_prev)
    {
      N_VDestroy(PI_CONTENT(DEE)->q_prev);
      PI_CONTENT(DEE)->q_prev = NULL;
    }
    if (PI_CONTENT(DEE)->V)
    {
      N_VDestroy(PI_CONTENT(DEE)->V);
      PI_CONTENT(DEE)->V = NULL;
    }
    if (PI_CONTENT(DEE)->rhs_linY)
    {
      N_VDestroy(PI_CONTENT(DEE)->rhs_linY);
      PI_CONTENT(DEE)->rhs_linY = NULL;
    }
    if (PI_CONTENT(DEE)->Fy)
    {
      N_VDestroy(PI_CONTENT(DEE)->Fy);
      PI_CONTENT(DEE)->Fy = NULL;
    }
    if (PI_CONTENT(DEE)->work)
    {
      N_VDestroy(PI_CONTENT(DEE)->work);
      PI_CONTENT(DEE)->work = NULL;
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

/*---------------------------------------------------------------
  dee_DQJtimes_Power:

  This routine generates a difference quotient approximation to
  the Jacobian-vector product f_y(t,y) * v. The approximation is
  Jv = [f(y + v*sig) - f(y)]/sig, where
      sig = sign(y^T v) * sqrt(unit roundoff)
            * max(|y^T v|, ||v||_1) / (v^T v).
  ---------------------------------------------------------------*/
SUNErrCode dee_DQJtimes_Power(void* voidstarDEE, N_Vector v, N_Vector Jv)
{
  SUNDomEigEstimator DEE = (SUNDomEigEstimator)voidstarDEE;
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(v, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Jv, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->rhsfn, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->rhs_linY, SUN_ERR_ARG_CORRUPT);

  sunrealtype vdotv = N_VDotProd(v, v);
  if (vdotv <= SUN_SMALL_REAL)
  {
    N_VScale(ZERO, v, Jv);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  if (PI_CONTENT(DEE)->work == NULL)
  {
    PI_CONTENT(DEE)->work = N_VClone(v);
    SUNCheckLastErr();
  }
  if (PI_CONTENT(DEE)->Fy == NULL)
  {
    PI_CONTENT(DEE)->Fy = N_VClone(v);
    SUNCheckLastErr();
  }

  sunrealtype sig, siginv;
  int iter, retval;

  N_Vector y    = PI_CONTENT(DEE)->rhs_linY;
  N_Vector work = PI_CONTENT(DEE)->work;
  N_Vector Fy   = PI_CONTENT(DEE)->Fy;

  retval = PI_CONTENT(DEE)->rhsfn(PI_CONTENT(DEE)->rhs_linT, y, Fy,
                                  PI_CONTENT(DEE)->rhs_data);
  PI_CONTENT(DEE)->nfevals++;
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
    retval = PI_CONTENT(DEE)->rhsfn(PI_CONTENT(DEE)->rhs_linT, work, Jv,
                                    PI_CONTENT(DEE)->rhs_data);
    PI_CONTENT(DEE)->nfevals++;
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
