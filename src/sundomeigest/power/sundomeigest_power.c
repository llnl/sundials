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
  SUNDomEigEstimator DEE, sunrealtype lambdaR, N_Vector v_prev, N_Vector v,
  sunrealtype* lambdaR_out, sunrealtype* lambdaI_out);

SUNErrCode dee_DQJtimes(void* voidstarDEE, N_Vector v, N_Vector Jv);

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

  /* Check if rel_tol > 0 */
  if (rel_tol < SUN_SMALL_REAL) { rel_tol = DEE_TOL_DEFAULT; }

  /* Create dominant eigenvalue estimator */
  DEE = NULL;
  DEE = SUNDomEigEstimator_NewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  DEE->ops->setatimes = SUNDomEigEstimator_SetATimes_Power;
  DEE->ops->setrhs    = SUNDomEigEstimator_SetRHS_Power;
  DEE->ops->setrhslinearizationvector =
    SUNDomEigEstimator_SetRHSLinearizationVector_Power;
  DEE->ops->setmaxiters = SUNDomEigEstimator_SetMaxIters_Power;
  DEE->ops->setnumpreprocessiters = SUNDomEigEstimator_SetNumPreprocessIters_Power;
  DEE->ops->setreltol         = SUNDomEigEstimator_SetRelTol_Power;
  DEE->ops->setinitialguess   = SUNDomEigEstimator_SetInitialGuess_Power;
  DEE->ops->initialize        = SUNDomEigEstimator_Initialize_Power;
  DEE->ops->estimate          = SUNDomEigEstimator_Estimate_Power;
  DEE->ops->getres            = SUNDomEigEstimator_GetRes_Power;
  DEE->ops->getnumiters       = SUNDomEigEstimator_GetNumIters_Power;
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
  content->rhs_linV    = NULL;
  content->Fv          = NULL;
  content->work        = NULL;
  content->complex     = SUNTRUE;
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

  N_VScale(ONE, q, content->q);
  SUNCheckLastErrNull();

  content->V = N_VClone(q);
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

SUNErrCode SUNDomEigEstimator_SetRHS_Power(SUNDomEigEstimator DEE,
                                           void* rhs_data, DEERhsFn RHSfn)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(RHSfn, SUN_ERR_ARG_CORRUPT);

  /* set function pointers to integrator-supplied RHS routine
     and data, and return with success */
  PI_CONTENT(DEE)->rhsfn    = RHSfn;
  PI_CONTENT(DEE)->rhs_data = rhs_data;

  DEE->ops->setatimes(DEE, (void*)DEE, dee_DQJtimes);
  SUNCheckLastErr();

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_SetRHSLinearizationVector_Power(SUNDomEigEstimator DEE,
                                                              N_Vector v)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(v, SUN_ERR_ARG_CORRUPT);

  if (PI_CONTENT(DEE)->rhs_linV == NULL)
  {
    PI_CONTENT(DEE)->rhs_linV = N_VClone(v);
    SUNCheckLastErr();
  }

  N_VScale(ONE, v, PI_CONTENT(DEE)->rhs_linV);
  SUNCheckLastErr();

  return SUN_SUCCESS;
}

SUNErrCode SUNDomEigEstimator_SetDEEisReal_Power(SUNDomEigEstimator DEE,
                                                 sunbooleantype real)
{
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);

  /* set the complex flag to the opposite of the real flag */
  PI_CONTENT(DEE)->complex = !real;

  /* q_prev is allocated in SUNDomEigEstimator_Initialize_Power, which is expected to be 
  called after this routine. If the user calls this routine after initialization, we need 
  to free q_prev here. */
  if (!(PI_CONTENT(DEE)->complex) && PI_CONTENT(DEE)->q_prev)
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

  if (PI_CONTENT(DEE)->rel_tol < SUN_SMALL_REAL)
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

  if (PI_CONTENT(DEE)->complex)
  {
    SUNAssert(PI_CONTENT(DEE)->q_prev == NULL, SUN_ERR_ARG_CORRUPT);
  }

  /* Initialize the vector V */
  sunrealtype normq = N_VDotProd(PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->q);
  SUNCheckLastErr();

  normq = SUNRsqrt(normq);

  N_VScale(ONE / normq, PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->V);
  SUNCheckLastErr();

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

  /* Check if rel_tol > 0 */
  if (rel_tol < SUN_SMALL_REAL) { rel_tol = DEE_TOL_DEFAULT; }

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
  if (PI_CONTENT(DEE)->complex && (PI_CONTENT(DEE)->q_prev == NULL))
  {
    /* allocate q_prev vector */
    PI_CONTENT(DEE)->q_prev = N_VClone(PI_CONTENT(DEE)->q);
    SUNCheckLastErr();
  }

  sunrealtype newlambdaR = ZERO;
  sunrealtype oldlambdaR = ZERO;

  int retval;
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
    if (PI_CONTENT(DEE)->complex)
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

    normq = N_VDotProd(PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->q);
    SUNCheckLastErr();

    normq = SUNRsqrt(normq);
    N_VScale(ONE / normq, PI_CONTENT(DEE)->q, PI_CONTENT(DEE)->V);
    SUNCheckLastErr();

    PI_CONTENT(DEE)->res = SUNRabs(newlambdaR - oldlambdaR);
    if (PI_CONTENT(DEE)->res <= PI_CONTENT(DEE)->rel_tol * SUNRabs(newlambdaR))
    {
      break;
    }

    oldlambdaR = newlambdaR;
  }

  if (PI_CONTENT(DEE)->complex)
  {
    retval = sundomeigestimator_complex_dom_eigs_from_PI(DEE, newlambdaR,
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
  SUNDomEigEstimator DEE, sunrealtype lambdaR, N_Vector v_prev, N_Vector v,
  sunrealtype* lambdaR_out, sunrealtype* lambdaI_out)
{
  SUNFunctionBegin(DEE->sunctx);

  int retval;
  sunrealtype cos_qs, det_G_inv, h11, h12, h21, h22, p11, p12, p21, p22;
  sunrealtype proj_cond_inv = SUN_RCONST(1e-10);
  cos_qs                    = N_VDotProd(v_prev, v);
  SUNCheckLastErr();

  if (SUNRabs(SUNRabs(cos_qs) - ONE) < proj_cond_inv)
  {
    /* Dominant eigenvalue is real */
    *lambdaR_out = lambdaR;
    *lambdaI_out = ZERO;
    return SUN_SUCCESS;
  }
  else
  {
    det_G_inv = ONE / (ONE - cos_qs * cos_qs);

    /* Solve for G = [v_prev v]' * [v_prev v] and compute 
    projected matrix P = G^{-1} * [v_prev v]' * A * [v_prev v] */

    retval = PI_CONTENT(DEE)->ATimes(PI_CONTENT(DEE)->ATdata, v_prev,
                                     PI_CONTENT(DEE)->q);
    PI_CONTENT(DEE)->num_ATimes++;
    PI_CONTENT(DEE)->num_iters++;
    if (retval != 0) { return SUN_ERR_USER_FCN_FAIL; }

    h11 = N_VDotProd(v_prev, PI_CONTENT(DEE)->q);
    h21 = N_VDotProd(v, PI_CONTENT(DEE)->q);

    retval = PI_CONTENT(DEE)->ATimes(PI_CONTENT(DEE)->ATdata, v,
                                     PI_CONTENT(DEE)->q);
    PI_CONTENT(DEE)->num_ATimes++;
    PI_CONTENT(DEE)->num_iters++;
    if (retval != 0) { return SUN_ERR_USER_FCN_FAIL; }

    h12 = N_VDotProd(v_prev, PI_CONTENT(DEE)->q);
    h22 = N_VDotProd(v, PI_CONTENT(DEE)->q);

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
      *lambdaR_out = (traceP + SUNRsqrt(discrim)) / SUN_RCONST(2.0);
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
    if (PI_CONTENT(DEE)->rhs_linV)
    {
      N_VDestroy(PI_CONTENT(DEE)->rhs_linV);
      PI_CONTENT(DEE)->rhs_linV = NULL;
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
  dee_DQJtimes:

  This routine generates a difference quotient approximation to
  the Jacobian-vector product f_y(t,y) * v. The approximation is
  Jv = [f(y + v*sig) - f(y)]/sig, where sig = 1 / ||v||_WRMS,
  i.e. the WRMS norm of v*sig is 1.
  ---------------------------------------------------------------*/
SUNErrCode dee_DQJtimes(void* voidstarDEE, N_Vector v, N_Vector Jv)
{
  SUNDomEigEstimator DEE = (SUNDomEigEstimator)voidstarDEE;
  SUNFunctionBegin(DEE->sunctx);

  SUNAssert(DEE, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE), SUN_ERR_ARG_CORRUPT);
  SUNAssert(v, SUN_ERR_ARG_CORRUPT);
  SUNAssert(Jv, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->rhsfn, SUN_ERR_ARG_CORRUPT);
  SUNAssert(PI_CONTENT(DEE)->rhs_linV, SUN_ERR_ARG_CORRUPT);

  if (PI_CONTENT(DEE)->work == NULL)
  {
    PI_CONTENT(DEE)->work = N_VClone(v);
    SUNCheckLastErr();
  }
  if (PI_CONTENT(DEE)->Fv == NULL)
  {
    PI_CONTENT(DEE)->Fv = N_VClone(v);
    SUNCheckLastErr();
  }

  sunrealtype sig, siginv;
  int iter, retval;

  N_Vector y    = PI_CONTENT(DEE)->rhs_linV;
  N_Vector work = PI_CONTENT(DEE)->work;
  N_Vector Fv   = PI_CONTENT(DEE)->Fv;

  retval = PI_CONTENT(DEE)->rhsfn(y, PI_CONTENT(DEE)->Fv,
                                  PI_CONTENT(DEE)->rhs_data);
  PI_CONTENT(DEE)->nfevals++;
  if (retval != 0) { return SUN_ERR_USER_FCN_FAIL; }

  /* Initialize perturbation to 1/||v|| */
  // sunrealtype sig = ONE / N_VWrmsNorm(v, ark_mem->ewt);
  sig = ONE / N_VWrmsNorm(v, v);

  for (iter = 0; iter < MAX_DQITERS; iter++)
  {
    /* Set work = y + sig*v */
    N_VLinearSum(sig, v, ONE, y, work);

    /* Set Jv = f(tn, y+sig*v) */
    retval = PI_CONTENT(DEE)->rhsfn(work, Jv, PI_CONTENT(DEE)->rhs_data);
    PI_CONTENT(DEE)->nfevals++;
    if (retval == 0) { break; }
    if (retval < 0) { return (-1); }

    /* If f failed recoverably, shrink sig and retry */
    sig *= SUN_RCONST(0.25);
  }

  /* Replace Jv by (Jv - fn)/sig */
  siginv = ONE / sig;
  N_VLinearSum(siginv, Jv, -siginv, Fv, Jv);

  return SUN_SUCCESS;
}