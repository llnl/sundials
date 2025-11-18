/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ UMBC
 * Based on sundials_sptfqmr.c code, written by Aaron Collier @ LLNL
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
 * This is the implementation file for the SPTFQMR implementation of
 * the SUNLINSOL package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_math.h>
#include <sunlinsol/sunlinsol_sptfqmr.h>
#include "sundials_logger_impl.h"

#include "sundials_cli.h"
#include "sundials_macros.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * SPTFQMR solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define SPTFQMR_CONTENT(S) ((SUNLinearSolverContent_SPTFQMR)(S->content))
#define LASTFLAG(S)        (SPTFQMR_CONTENT(S)->last_flag)

/*
 * ----------------------------------------------------------------------------
 * Un-exported implementation specific routines
 * ----------------------------------------------------------------------------
 */

static SUNErrCode setFromCommandLine_SPTFQMR(SUNLinearSolver S, const char* LSid,
                                             int argc, char* argv[]);

SUNErrCode SUNLinSolSetOptions_SPTFQMR(SUNLinearSolver S, const char* LSid,
                                       const char* file_name, int argc,
                                       char* argv[]);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new SPTFQMR linear solver
 */

SUNLinearSolver SUNLinSol_SPTFQMR(N_Vector y, int pretype, int maxl,
                                  SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNLinearSolver S;
  SUNLinearSolverContent_SPTFQMR content;

  /* check for legal pretype and maxl values; if illegal use defaults */
  if ((pretype != SUN_PREC_NONE) && (pretype != SUN_PREC_LEFT) &&
      (pretype != SUN_PREC_RIGHT) && (pretype != SUN_PREC_BOTH))
  {
    pretype = SUN_PREC_NONE;
  }
  if (maxl <= 0) { maxl = SUNSPTFQMR_MAXL_DEFAULT; }

  /* check that the supplied N_Vector supports all requisite operations */
  SUNAssertNull((y->ops->nvclone) && (y->ops->nvdestroy) &&
                  (y->ops->nvlinearsum) && (y->ops->nvconst) && (y->ops->nvprod) &&
                  (y->ops->nvdiv) && (y->ops->nvscale) && (y->ops->nvdotprod),
                SUN_ERR_ARG_INCOMPATIBLE);

  /* Create linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  S->ops->gettype           = SUNLinSolGetType_SPTFQMR;
  S->ops->getid             = SUNLinSolGetID_SPTFQMR;
  S->ops->setatimes         = SUNLinSolSetATimes_SPTFQMR;
  S->ops->setoptions        = SUNLinSolSetOptions_SPTFQMR;
  S->ops->setpreconditioner = SUNLinSolSetPreconditioner_SPTFQMR;
  S->ops->setscalingvectors = SUNLinSolSetScalingVectors_SPTFQMR;
  S->ops->setzeroguess      = SUNLinSolSetZeroGuess_SPTFQMR;
  S->ops->initialize        = SUNLinSolInitialize_SPTFQMR;
  S->ops->setup             = SUNLinSolSetup_SPTFQMR;
  S->ops->solve             = SUNLinSolSolve_SPTFQMR;
  S->ops->numiters          = SUNLinSolNumIters_SPTFQMR;
  S->ops->resnorm           = SUNLinSolResNorm_SPTFQMR;
  S->ops->resid             = SUNLinSolResid_SPTFQMR;
  S->ops->lastflag          = SUNLinSolLastFlag_SPTFQMR;
  S->ops->space             = SUNLinSolSpace_SPTFQMR;
  S->ops->free              = SUNLinSolFree_SPTFQMR;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_SPTFQMR)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  S->content = content;

  /* Fill content */
  content->last_flag = 0;
  content->maxl      = maxl;
  content->pretype   = pretype;
  content->zeroguess = SUNFALSE;
  content->numiters  = 0;
  content->resnorm   = ZERO;
  content->r_star    = NULL;
  content->q         = NULL;
  content->d         = NULL;
  content->v         = NULL;
  content->p         = NULL;
  content->r         = NULL;
  content->u         = NULL;
  content->vtemp1    = NULL;
  content->vtemp2    = NULL;
  content->vtemp3    = NULL;
  content->s1        = NULL;
  content->s2        = NULL;
  content->ATimes    = NULL;
  content->ATData    = NULL;
  content->Psetup    = NULL;
  content->Psolve    = NULL;
  content->PData     = NULL;

  /* Allocate content */
  content->r_star = N_VClone(y);
  SUNCheckLastErrNull();
  content->q = N_VClone(y);
  SUNCheckLastErrNull();
  content->d = N_VClone(y);
  SUNCheckLastErrNull();
  content->v = N_VClone(y);
  SUNCheckLastErrNull();
  content->p = N_VClone(y);
  SUNCheckLastErrNull();
  content->r = N_VCloneVectorArray(2, y);
  SUNCheckLastErrNull();
  content->u = N_VClone(y);
  SUNCheckLastErrNull();
  content->vtemp1 = N_VClone(y);
  SUNCheckLastErrNull();
  content->vtemp2 = N_VClone(y);
  SUNCheckLastErrNull();
  content->vtemp3 = N_VClone(y);
  SUNCheckLastErrNull();

  return (S);
}

/* ----------------------------------------------------------------------------
 * Function to control set routines via the command line or file
 */

SUNErrCode SUNLinSolSetOptions_SPTFQMR(SUNLinearSolver S, const char* LSid,
                                       SUNDIALS_MAYBE_UNUSED const char* file_name,
                                       int argc, char* argv[])
{
  SUNFunctionBegin(S->sunctx);

  /* File-based option control is currently unimplemented */
  SUNAssert((file_name == NULL || strlen(file_name) == 0),
            SUN_ERR_ARG_INCOMPATIBLE);

  if (argc > 0 && argv != NULL)
  {
    SUNCheckCall(setFromCommandLine_SPTFQMR(S, LSid, argc, argv));
  }

  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Function to control set routines via the command line
 */

static SUNErrCode setFromCommandLine_SPTFQMR(SUNLinearSolver S, const char* LSid,
                                             int argc, char* argv[])
{
  SUNFunctionBegin(S->sunctx);

  /* Prefix for options to set */
  const char* default_id = "sunlinearsolver";
  size_t offset          = strlen(default_id) + 1;
  if (LSid != NULL && strlen(LSid) > 0) { offset = strlen(LSid) + 1; }
  char* prefix = (char*)malloc(sizeof(char) * (offset + 1));
  if (LSid != NULL && strlen(LSid) > 0) { strcpy(prefix, LSid); }
  else { strcpy(prefix, default_id); }
  strcat(prefix, ".");

  for (int idx = 1; idx < argc; idx++)
  {
    int retval;

    /* skip command-line arguments that do not begin with correct prefix */
    if (strncmp(argv[idx], prefix, strlen(prefix)) != 0) { continue; }

    /* control over PrecType function */
    if (strcmp(argv[idx] + offset, "prec_type") == 0)
    {
      idx += 1;
      int iarg = atoi(argv[idx]);
      retval   = SUNLinSol_SPTFQMRSetPrecType(S, iarg);
      if (retval != SUN_SUCCESS)
      {
        free(prefix);
        return retval;
      }
      continue;
    }

    /* control over Maxl function */
    if (strcmp(argv[idx] + offset, "maxl") == 0)
    {
      idx += 1;
      int iarg = atoi(argv[idx]);
      retval   = SUNLinSol_SPTFQMRSetMaxl(S, iarg);
      if (retval != SUN_SUCCESS)
      {
        free(prefix);
        return retval;
      }
      continue;
    }
  }
  free(prefix);
  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Function to set the type of preconditioning for SPTFQMR to use
 */

SUNErrCode SUNLinSol_SPTFQMRSetPrecType(SUNLinearSolver S, int pretype)
{
  SUNFunctionBegin(S->sunctx);
  /* Check for legal pretype */
  SUNAssert((pretype == SUN_PREC_NONE) || (pretype == SUN_PREC_LEFT) ||
              (pretype == SUN_PREC_RIGHT) || (pretype == SUN_PREC_BOTH),
            SUN_ERR_ARG_OUTOFRANGE);

  /* Set pretype */
  SPTFQMR_CONTENT(S)->pretype = pretype;
  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Function to set the maximum number of iterations for SPTFQMR to use
 */

SUNErrCode SUNLinSol_SPTFQMRSetMaxl(SUNLinearSolver S, int maxl)
{
  /* Check for legal pretype */
  if (maxl <= 0) { maxl = SUNSPTFQMR_MAXL_DEFAULT; }

  /* Set pretype */
  SPTFQMR_CONTENT(S)->maxl = maxl;
  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_SPTFQMR(SUNDIALS_MAYBE_UNUSED SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_ITERATIVE);
}

SUNLinearSolver_ID SUNLinSolGetID_SPTFQMR(SUNDIALS_MAYBE_UNUSED SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_SPTFQMR);
}

SUNErrCode SUNLinSolInitialize_SPTFQMR(SUNLinearSolver S)
{
  SUNFunctionBegin(S->sunctx);
  SUNLinearSolverContent_SPTFQMR content;

  /* set shortcut to SPTFQMR memory structure */
  content = SPTFQMR_CONTENT(S);

  /* ensure valid options */
  if (content->maxl <= 0) { content->maxl = SUNSPTFQMR_MAXL_DEFAULT; }

  SUNAssert(content->ATimes, SUN_ERR_ARG_CORRUPT);

  if ((content->pretype != SUN_PREC_LEFT) &&
      (content->pretype != SUN_PREC_RIGHT) && (content->pretype != SUN_PREC_BOTH))
  {
    content->pretype = SUN_PREC_NONE;
  }

  SUNAssert((content->pretype == SUN_PREC_NONE) || (content->Psolve != NULL),
            SUN_ERR_ARG_CORRUPT);

  /* no additional memory to allocate */

  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetATimes_SPTFQMR(SUNLinearSolver S, void* ATData,
                                      SUNATimesFn ATimes)
{
  /* set function pointers to integrator-supplied ATimes routine
     and data, and return with success */
  SPTFQMR_CONTENT(S)->ATimes = ATimes;
  SPTFQMR_CONTENT(S)->ATData = ATData;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetPreconditioner_SPTFQMR(SUNLinearSolver S, void* PData,
                                              SUNPSetupFn Psetup,
                                              SUNPSolveFn Psolve)
{
  /* set function pointers to integrator-supplied Psetup and PSolve
     routines and data, and return with success */
  SPTFQMR_CONTENT(S)->Psetup = Psetup;
  SPTFQMR_CONTENT(S)->Psolve = Psolve;
  SPTFQMR_CONTENT(S)->PData  = PData;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetScalingVectors_SPTFQMR(SUNLinearSolver S, N_Vector s1,
                                              N_Vector s2)
{
  /* set N_Vector pointers to integrator-supplied scaling vectors,
     and return with success */
  SPTFQMR_CONTENT(S)->s1 = s1;
  SPTFQMR_CONTENT(S)->s2 = s2;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetZeroGuess_SPTFQMR(SUNLinearSolver S, sunbooleantype onoff)
{
  /* set flag indicating a zero initial guess */
  SPTFQMR_CONTENT(S)->zeroguess = onoff;
  return SUN_SUCCESS;
}

int SUNLinSolSetup_SPTFQMR(SUNLinearSolver S, SUNDIALS_MAYBE_UNUSED SUNMatrix A)
{
  SUNFunctionBegin(S->sunctx);

  int status = SUN_SUCCESS;
  SUNPSetupFn Psetup;
  void* PData;

  /* Set shortcuts to SPTFQMR memory structures */
  Psetup = SPTFQMR_CONTENT(S)->Psetup;
  PData  = SPTFQMR_CONTENT(S)->PData;

  /* no solver-specific setup is required, but if user-supplied
     Psetup routine exists, call that here */
  if (Psetup != NULL)
  {
    status = Psetup(PData);
    if (status != 0)
    {
      LASTFLAG(S) = (status < 0) ? SUNLS_PSET_FAIL_UNREC : SUNLS_PSET_FAIL_REC;
      return (LASTFLAG(S));
    }
  }

  /* return with success */
  LASTFLAG(S) = SUN_SUCCESS;
  return SUN_SUCCESS;
}

int SUNLinSolSolve_SPTFQMR(SUNLinearSolver S, SUNDIALS_MAYBE_UNUSED SUNMatrix A,
                           N_Vector x, N_Vector b, sunrealtype delta)
{
  SUNFunctionBegin(S->sunctx);

  /* local data and shortcut variables */
  sunrealtype alpha, tau, eta, beta, c, sigma, v_bar, omega;
  sunrealtype rho[2];
  sunrealtype r_init_norm, r_curr_norm;
  sunrealtype temp_val;
  sunbooleantype preOnLeft, preOnRight, scale_x, scale_b, converged, b_ok;
  sunbooleantype* zeroguess;
  int n, m, l_max;
  void *A_data, *P_data;
  SUNATimesFn atimes;
  SUNPSolveFn psolve;
  sunrealtype* res_norm;
  int* nli;
  N_Vector sx, sb, r_star, q, d, v, p, *r, u, vtemp1, vtemp2, vtemp3;
  sunrealtype cv[3];
  N_Vector Xv[3];
  int status = SUN_SUCCESS;

  /* Make local shortcuts to solver variables. */
  l_max     = SPTFQMR_CONTENT(S)->maxl;
  r_star    = SPTFQMR_CONTENT(S)->r_star;
  q         = SPTFQMR_CONTENT(S)->q;
  d         = SPTFQMR_CONTENT(S)->d;
  v         = SPTFQMR_CONTENT(S)->v;
  p         = SPTFQMR_CONTENT(S)->p;
  r         = SPTFQMR_CONTENT(S)->r;
  u         = SPTFQMR_CONTENT(S)->u;
  vtemp1    = SPTFQMR_CONTENT(S)->vtemp1;
  vtemp2    = SPTFQMR_CONTENT(S)->vtemp2;
  vtemp3    = SPTFQMR_CONTENT(S)->vtemp3;
  sb        = SPTFQMR_CONTENT(S)->s1;
  sx        = SPTFQMR_CONTENT(S)->s2;
  A_data    = SPTFQMR_CONTENT(S)->ATData;
  P_data    = SPTFQMR_CONTENT(S)->PData;
  atimes    = SPTFQMR_CONTENT(S)->ATimes;
  psolve    = SPTFQMR_CONTENT(S)->Psolve;
  zeroguess = &(SPTFQMR_CONTENT(S)->zeroguess);
  nli       = &(SPTFQMR_CONTENT(S)->numiters);
  res_norm  = &(SPTFQMR_CONTENT(S)->resnorm);

  /* Initialize counters and convergence flag */
  temp_val = r_curr_norm = -ONE;
  *nli                   = 0;
  converged              = SUNFALSE;
  b_ok                   = SUNFALSE;

  /* set sunbooleantype flags for internal solver options */
  preOnLeft  = ((SPTFQMR_CONTENT(S)->pretype == SUN_PREC_LEFT) ||
               (SPTFQMR_CONTENT(S)->pretype == SUN_PREC_BOTH));
  preOnRight = ((SPTFQMR_CONTENT(S)->pretype == SUN_PREC_RIGHT) ||
                (SPTFQMR_CONTENT(S)->pretype == SUN_PREC_BOTH));
  scale_x    = (sx != NULL);
  scale_b    = (sb != NULL);

  SUNLogInfo(S->sunctx->logger, "linear-solver", "solver = sptfqmr");

  SUNLogInfo(S->sunctx->logger, "begin-iterations-list", "");

  /* Check for unsupported use case */
  if (preOnRight && !(*zeroguess))
  {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUN_ERR_ARG_INCOMPATIBLE;

    SUNLogInfo(S->sunctx->logger, "end-iterations-list",
               "status = failed unsupported configuration");

    return SUN_ERR_ARG_INCOMPATIBLE;
  }

  /* Check if Atimes function has been set */
  SUNAssert(atimes, SUN_ERR_ARG_CORRUPT);

  /* If preconditioning, check if psolve has been set */
  SUNAssert(!(preOnLeft || preOnRight) || psolve, SUN_ERR_ARG_CORRUPT);

  /* Set r_star to initial (unscaled) residual r_star = r_0 = b - A*x_0 */
  /* NOTE: if x == 0 then just set residual to b and continue */
  if (*zeroguess)
  {
    N_VScale(ONE, b, r_star);
    SUNCheckLastErr();
  }
  else
  {
    status = atimes(A_data, x, r_star);
    if (status != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (status < 0) ? SUNLS_ATIMES_FAIL_UNREC
                                 : SUNLS_ATIMES_FAIL_REC;

      SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                 "status = failed matvec, retval = %d", status);

      return (LASTFLAG(S));
    }
    N_VLinearSum(ONE, b, -ONE, r_star, r_star);
    SUNCheckLastErr();
  }

  /* Apply left preconditioner and b-scaling to r_star (or really just r_0) */
  if (preOnLeft)
  {
    status = psolve(P_data, r_star, vtemp1, delta, SUN_PREC_LEFT);
    if (status != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                 : SUNLS_PSOLVE_FAIL_REC;

      SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                 "status = failed preconditioner solve, retval = %d", status);

      return (LASTFLAG(S));
    }
  }
  else
  {
    N_VScale(ONE, r_star, vtemp1);
    SUNCheckLastErr();
  }

  if (scale_b)
  {
    N_VProd(sb, vtemp1, r_star);
    SUNCheckLastErr();
  }
  else
  {
    N_VScale(ONE, vtemp1, r_star);
    SUNCheckLastErr();
  }

  /* Initialize rho[0] */
  /* NOTE: initialized here to reduce number of computations - avoid need
           to compute r_star^T*r_star twice, and avoid needlessly squaring
           values */
  rho[0] = N_VDotProd(r_star, r_star);
  SUNCheckLastErr();

  /* Compute norm of initial residual (r_0) to see if we really need
     to do anything */
  *res_norm = r_init_norm = SUNRsqrt(rho[0]);

  if (r_init_norm <= delta)
  {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUN_SUCCESS;

    SUNLogInfo(S->sunctx->logger, "end-iterations-list",
               "cur-iter = 0, res-norm = " SUN_FORMAT_G ", status = success",
               *res_norm);

    return (LASTFLAG(S));
  }

  SUNLogInfo(S->sunctx->logger, "linear-iterate",
             "cur-iter = 0, res-norm = " SUN_FORMAT_G, *res_norm);

  /* Set v = A*r_0 (preconditioned and scaled) */
  if (scale_x)
  {
    N_VDiv(r_star, sx, vtemp1);
    SUNCheckLastErr();
  }
  else
  {
    N_VScale(ONE, r_star, vtemp1);
    SUNCheckLastErr();
  }

  if (preOnRight)
  {
    N_VScale(ONE, vtemp1, v);
    SUNCheckLastErr();
    status = psolve(P_data, v, vtemp1, delta, SUN_PREC_RIGHT);
    if (status != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                 : SUNLS_PSOLVE_FAIL_REC;

      SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                 "status = failed preconditioner solve, retval = %d", status);

      return (LASTFLAG(S));
    }
  }

  status = atimes(A_data, vtemp1, v);
  if (status != 0)
  {
    *zeroguess = SUNFALSE;
    LASTFLAG(S) = (status < 0) ? SUNLS_ATIMES_FAIL_UNREC : SUNLS_ATIMES_FAIL_REC;

    SUNLogInfo(S->sunctx->logger, "end-iterations-list",
               "status = failed matvec, retval = %d", status);

    return (LASTFLAG(S));
  }

  if (preOnLeft)
  {
    status = psolve(P_data, v, vtemp1, delta, SUN_PREC_LEFT);
    if (status != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                 : SUNLS_PSOLVE_FAIL_REC;

      SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                 "status = failed preconditioner solve, retval = %d", status);

      return (LASTFLAG(S));
    }
  }
  else
  {
    N_VScale(ONE, v, vtemp1);
    SUNCheckLastErr();
  }

  if (scale_b)
  {
    N_VProd(sb, vtemp1, v);
    SUNCheckLastErr();
  }
  else
  {
    N_VScale(ONE, vtemp1, v);
    SUNCheckLastErr();
  }

  /* Initialize remaining variables */
  N_VScale(ONE, r_star, r[0]);
  SUNCheckLastErr();
  N_VScale(ONE, r_star, u);
  SUNCheckLastErr();
  N_VScale(ONE, r_star, p);
  SUNCheckLastErr();
  N_VConst(ZERO, d);
  SUNCheckLastErr();

  /* Set x = sx x if non-zero guess */
  if (scale_x && !(*zeroguess))
  {
    N_VProd(sx, x, x);
    SUNCheckLastErr();
  }

  tau   = r_init_norm;
  v_bar = eta = ZERO;

  SUNLogInfo(S->sunctx->logger, "end-iterations-list", "status = continue");

  /* START outer loop */
  for (n = 0; n < l_max; ++n)
  {
    SUNLogInfo(S->sunctx->logger, "begin-iterations-list", "");

    /* Increment linear iteration counter */
    (*nli)++;

    /* sigma = r_star^T*v */
    sigma = N_VDotProd(r_star, v);
    SUNCheckLastErr();

    /* alpha = rho[0]/sigma */
    alpha = rho[0] / sigma;

    /* q = u-alpha*v */
    N_VLinearSum(ONE, u, -alpha, v, q);
    SUNCheckLastErr();

    /* r[1] = r[0]-alpha*A*(u+q) */
    N_VLinearSum(ONE, u, ONE, q, r[1]);
    SUNCheckLastErr();
    if (scale_x)
    {
      N_VDiv(r[1], sx, r[1]);
      SUNCheckLastErr();
    }

    if (preOnRight)
    {
      N_VScale(ONE, r[1], vtemp1);
      SUNCheckLastErr();
      status = psolve(P_data, vtemp1, r[1], delta, SUN_PREC_RIGHT);
      if (status != 0)
      {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                   : SUNLS_PSOLVE_FAIL_REC;

        SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                   "status = failed preconditioner solve, retval = %d", status);

        return (LASTFLAG(S));
      }
    }

    status = atimes(A_data, r[1], vtemp1);
    if (status != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (status < 0) ? SUNLS_ATIMES_FAIL_UNREC
                                 : SUNLS_ATIMES_FAIL_REC;

      SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                 "status = failed matvec, retval = %d", status);

      return (LASTFLAG(S));
    }

    if (preOnLeft)
    {
      status = psolve(P_data, vtemp1, r[1], delta, SUN_PREC_LEFT);
      if (status != 0)
      {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                   : SUNLS_PSOLVE_FAIL_REC;

        SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                   "status = failed preconditioner solve, retval = %d", status);

        return (LASTFLAG(S));
      }
    }
    else
    {
      N_VScale(ONE, vtemp1, r[1]);
      SUNCheckLastErr();
    }

    if (scale_b)
    {
      N_VProd(sb, r[1], vtemp1);
      SUNCheckLastErr();
    }
    else
    {
      N_VScale(ONE, r[1], vtemp1);
      SUNCheckLastErr();
    }
    N_VLinearSum(ONE, r[0], -alpha, vtemp1, r[1]);
    SUNCheckLastErr();

    /* START inner loop */
    for (m = 0; m < 2; ++m)
    {
      /* d = [*]+(v_bar^2*eta/alpha)*d */
      /* NOTES:
       *   (1) [*] = u if m == 0, and q if m == 1
       *   (2) using temp_val reduces the number of required computations
       *       if the inner loop is executed twice
       */
      if (m == 0)
      {
        temp_val = N_VDotProd(r[1], r[1]);
        SUNCheckLastErr();
        temp_val = SUNRsqrt(temp_val);
        omega    = N_VDotProd(r[0], r[0]);
        SUNCheckLastErr();
        omega = SUNRsqrt(SUNRsqrt(omega) * temp_val);
        N_VLinearSum(ONE, u, SUNSQR(v_bar) * eta / alpha, d, d);
        SUNCheckLastErr();
      }
      else
      {
        omega = temp_val;
        N_VLinearSum(ONE, q, SUNSQR(v_bar) * eta / alpha, d, d);
        SUNCheckLastErr();
      }

      /* v_bar = omega/tau */
      v_bar = omega / tau;

      /* c = (1+v_bar^2)^(-1/2) */
      c = ONE / SUNRsqrt(ONE + SUNSQR(v_bar));

      /* tau = tau*v_bar*c */
      tau = tau * v_bar * c;

      /* eta = c^2*alpha */
      eta = SUNSQR(c) * alpha;

      /* x = x+eta*d */
      if (n == 0 && m == 0 && *zeroguess)
      {
        N_VScale(eta, d, x);
        SUNCheckLastErr();
      }
      else
      {
        N_VLinearSum(ONE, x, eta, d, x);
        SUNCheckLastErr();
      }

      /* Check for convergence... */
      /* NOTE: just use approximation to norm of residual, if possible */
      *res_norm = r_curr_norm = tau * SUNRsqrt(m + 1);

      /* Exit inner loop if iteration has converged based upon approximation
         to norm of current residual */
      if (r_curr_norm <= delta)
      {
        SUNLogInfo(S->sunctx->logger, "linear-iterate",
                   "cur-iter = %i, res-norm = " SUN_FORMAT_G, n + 1, *nli,
                   *res_norm);

        converged = SUNTRUE;
        break;
      }

      /* Decide if actual norm of residual vector should be computed */
      /* NOTES:
       *   (1) if r_curr_norm > delta, then check if actual residual norm
       *       is OK (recall we first compute an approximation)
       *   (2) if r_curr_norm >= r_init_norm and m == 1 and n == l_max, then
       *       compute actual residual norm to see if the iteration can be
       *       saved
       *   (3) the scaled and preconditioned right-hand side of the given
       *       linear system (denoted by b) is only computed once, and the
       *       result is stored in vtemp3 so it can be reused - reduces the
       *       number of psovles if using left preconditioning
       */
      if ((r_curr_norm > delta) ||
          (r_curr_norm >= r_init_norm && m == 1 && n == l_max))
      {
        /* Compute norm of residual ||b-A*x||_2 (preconditioned and scaled) */
        if (scale_x)
        {
          N_VDiv(x, sx, vtemp1);
          SUNCheckLastErr();
        }
        else
        {
          N_VScale(ONE, x, vtemp1);
          SUNCheckLastErr();
        }

        if (preOnRight)
        {
          status = psolve(P_data, vtemp1, vtemp2, delta, SUN_PREC_RIGHT);
          if (status != 0)
          {
            *zeroguess  = SUNFALSE;
            LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                       : SUNLS_PSOLVE_FAIL_REC;

            SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                       "status = failed preconditioner solve, retval = %d",
                       status);

            return (LASTFLAG(S));
          }
          N_VScale(ONE, vtemp2, vtemp1);
          SUNCheckLastErr();
        }

        status = atimes(A_data, vtemp1, vtemp2);
        if (status != 0)
        {
          *zeroguess  = SUNFALSE;
          LASTFLAG(S) = (status < 0) ? SUNLS_ATIMES_FAIL_UNREC
                                     : SUNLS_ATIMES_FAIL_REC;

          SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                     "status = failed matvec, retval = %d", status);

          return (LASTFLAG(S));
        }

        if (preOnLeft)
        {
          status = psolve(P_data, vtemp2, vtemp1, delta, SUN_PREC_LEFT);
          if (status != 0)
          {
            *zeroguess  = SUNFALSE;
            LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                       : SUNLS_PSOLVE_FAIL_REC;

            SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                       "status = failed preconditioner solve, retval = %d",
                       status);

            return (LASTFLAG(S));
          }
        }
        else
        {
          N_VScale(ONE, vtemp2, vtemp1);
          SUNCheckLastErr();
        }

        if (scale_b)
        {
          N_VProd(sb, vtemp1, vtemp2);
          SUNCheckLastErr();
        }
        else
        {
          N_VScale(ONE, vtemp1, vtemp2);
          SUNCheckLastErr();
        }

        /* Only precondition and scale b once (result saved for reuse) */
        if (!b_ok)
        {
          b_ok = SUNTRUE;
          if (preOnLeft)
          {
            status = psolve(P_data, b, vtemp3, delta, SUN_PREC_LEFT);
            if (status != 0)
            {
              *zeroguess  = SUNFALSE;
              LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                         : SUNLS_PSOLVE_FAIL_REC;

              SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                         "status = failed preconditioner solve, retval = %d",
                         status);

              return (LASTFLAG(S));
            }
          }
          else
          {
            N_VScale(ONE, b, vtemp3);
            SUNCheckLastErr();
          }

          if (scale_b)
          {
            N_VProd(sb, vtemp3, vtemp3);
            SUNCheckLastErr();
          }
        }
        N_VLinearSum(ONE, vtemp3, -ONE, vtemp2, vtemp1);
        SUNCheckLastErr();
        r_curr_norm = N_VDotProd(vtemp1, vtemp1);
        SUNCheckLastErr();
        *res_norm = r_curr_norm = SUNRsqrt(r_curr_norm);

        /* Exit inner loop if inequality condition is satisfied
           (meaning exit if we have converged) */
        if (r_curr_norm <= delta)
        {
          SUNLogInfo(S->sunctx->logger, "linear-iterate",
                     "cur-iter = %i, res-norm = " SUN_FORMAT_G, n + 1, *nli,
                     *res_norm);

          converged = SUNTRUE;
          break;
        }
      }

      SUNLogInfo(S->sunctx->logger, "linear-iterate",
                 "cur-iter = %i, res-norm = " SUN_FORMAT_G, n + 1, *nli,
                 *res_norm);

    } /* END inner loop */

    /* If converged, then exit outer loop as well */
    if (converged == SUNTRUE) { break; }

    /* rho[1] = r_star^T*r_[1] */
    rho[1] = N_VDotProd(r_star, r[1]);
    SUNCheckLastErr();

    /* beta = rho[1]/rho[0] */
    beta = rho[1] / rho[0];

    /* u = r[1]+beta*q */
    N_VLinearSum(ONE, r[1], beta, q, u);
    SUNCheckLastErr();

    /* p = u+beta*(q+beta*p) = beta*beta*p + beta*q + u */
    cv[0] = SUNSQR(beta);
    Xv[0] = p;

    cv[1] = beta;
    Xv[1] = q;

    cv[2] = ONE;
    Xv[2] = u;

    SUNCheckCall(N_VLinearCombination(3, cv, Xv, p));

    /* v = A*p */
    if (scale_x)
    {
      N_VDiv(p, sx, vtemp1);
      SUNCheckLastErr();
    }
    else
    {
      N_VScale(ONE, p, vtemp1);
      SUNCheckLastErr();
    }

    if (preOnRight)
    {
      N_VScale(ONE, vtemp1, v);
      SUNCheckLastErr();
      status = psolve(P_data, v, vtemp1, delta, SUN_PREC_RIGHT);
      if (status != 0)
      {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                   : SUNLS_PSOLVE_FAIL_REC;

        SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                   "status = failed preconditioner solve, retval = %d", status);

        return (LASTFLAG(S));
      }
    }

    status = atimes(A_data, vtemp1, v);
    if (status != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (status < 0) ? SUNLS_ATIMES_FAIL_UNREC
                                 : SUNLS_ATIMES_FAIL_REC;

      SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                 "status = failed matvec, retval = %d", status);

      return (LASTFLAG(S));
    }

    if (preOnLeft)
    {
      status = psolve(P_data, v, vtemp1, delta, SUN_PREC_LEFT);
      if (status != 0)
      {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                   : SUNLS_PSOLVE_FAIL_REC;

        SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                   "status = failed preconditioner solve, retval = %d", status);

        return (LASTFLAG(S));
      }
    }
    else
    {
      N_VScale(ONE, v, vtemp1);
      SUNCheckLastErr();
    }

    if (scale_b)
    {
      N_VProd(sb, vtemp1, v);
      SUNCheckLastErr();
    }
    else
    {
      N_VScale(ONE, vtemp1, v);
      SUNCheckLastErr();
    }

    /* Shift variable values */
    /* NOTE: reduces storage requirements */
    N_VScale(ONE, r[1], r[0]);
    SUNCheckLastErr();
    rho[0] = rho[1];

    SUNLogInfoIf(n < l_max - 1, S->sunctx->logger, "end-iterations-list",
                 "status = continue");

  } /* END outer loop */

  /* Determine return value */
  /* If iteration converged or residual was reduced, then return current iterate
   * (x) */
  if ((converged == SUNTRUE) || (r_curr_norm < r_init_norm))
  {
    if (scale_x)
    {
      N_VDiv(x, sx, x);
      SUNCheckLastErr();
    }

    if (preOnRight)
    {
      status = psolve(P_data, x, vtemp1, delta, SUN_PREC_RIGHT);
      if (status != 0)
      {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                   : SUNLS_PSOLVE_FAIL_REC;

        SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                   "status = failed preconditioner solve, retval = %d", status);

        return (LASTFLAG(S));
      }
      N_VScale(ONE, vtemp1, x);
      SUNCheckLastErr();
    }

    *zeroguess = SUNFALSE;
    if (converged == SUNTRUE)
    {
      SUNLogInfo(S->sunctx->logger, "end-iterations-list", "status = success");

      LASTFLAG(S) = SUN_SUCCESS;
    }
    else
    {
      SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                 "status = failed residual reduced");

      LASTFLAG(S) = SUNLS_RES_REDUCED;
    }
    return (LASTFLAG(S));
  }
  else
  {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUNLS_CONV_FAIL;

    SUNLogInfo(S->sunctx->logger, "end-iterations-list",
               "status = failed max iterations");

    return (LASTFLAG(S));
  }
}

int SUNLinSolNumIters_SPTFQMR(SUNLinearSolver S)
{
  return (SPTFQMR_CONTENT(S)->numiters);
}

sunrealtype SUNLinSolResNorm_SPTFQMR(SUNLinearSolver S)
{
  return (SPTFQMR_CONTENT(S)->resnorm);
}

N_Vector SUNLinSolResid_SPTFQMR(SUNLinearSolver S)
{
  return (SPTFQMR_CONTENT(S)->vtemp1);
}

sunindextype SUNLinSolLastFlag_SPTFQMR(SUNLinearSolver S)
{
  return (LASTFLAG(S));
}

SUNErrCode SUNLinSolSpace_SPTFQMR(SUNLinearSolver S, long int* lenrwLS,
                                  long int* leniwLS)
{
  SUNFunctionBegin(S->sunctx);
  sunindextype liw1, lrw1;
  if (SPTFQMR_CONTENT(S)->vtemp1->ops->nvspace)
  {
    N_VSpace(SPTFQMR_CONTENT(S)->vtemp1, &lrw1, &liw1);
    SUNCheckLastErr();
  }
  else { lrw1 = liw1 = 0; }
  *lenrwLS = lrw1 * 11;
  *leniwLS = liw1 * 11;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolFree_SPTFQMR(SUNLinearSolver S)
{
  if (S == NULL) { return SUN_SUCCESS; }

  if (S->content)
  {
    /* delete items from within the content structure */
    if (SPTFQMR_CONTENT(S)->r_star)
    {
      N_VDestroy(SPTFQMR_CONTENT(S)->r_star);
      SPTFQMR_CONTENT(S)->r_star = NULL;
    }
    if (SPTFQMR_CONTENT(S)->q)
    {
      N_VDestroy(SPTFQMR_CONTENT(S)->q);
      SPTFQMR_CONTENT(S)->q = NULL;
    }
    if (SPTFQMR_CONTENT(S)->d)
    {
      N_VDestroy(SPTFQMR_CONTENT(S)->d);
      SPTFQMR_CONTENT(S)->d = NULL;
    }
    if (SPTFQMR_CONTENT(S)->v)
    {
      N_VDestroy(SPTFQMR_CONTENT(S)->v);
      SPTFQMR_CONTENT(S)->v = NULL;
    }
    if (SPTFQMR_CONTENT(S)->p)
    {
      N_VDestroy(SPTFQMR_CONTENT(S)->p);
      SPTFQMR_CONTENT(S)->p = NULL;
    }
    if (SPTFQMR_CONTENT(S)->r)
    {
      N_VDestroyVectorArray(SPTFQMR_CONTENT(S)->r, 2);
      SPTFQMR_CONTENT(S)->r = NULL;
    }
    if (SPTFQMR_CONTENT(S)->u)
    {
      N_VDestroy(SPTFQMR_CONTENT(S)->u);
      SPTFQMR_CONTENT(S)->u = NULL;
    }
    if (SPTFQMR_CONTENT(S)->vtemp1)
    {
      N_VDestroy(SPTFQMR_CONTENT(S)->vtemp1);
      SPTFQMR_CONTENT(S)->vtemp1 = NULL;
    }
    if (SPTFQMR_CONTENT(S)->vtemp2)
    {
      N_VDestroy(SPTFQMR_CONTENT(S)->vtemp2);
      SPTFQMR_CONTENT(S)->vtemp2 = NULL;
    }
    if (SPTFQMR_CONTENT(S)->vtemp3)
    {
      N_VDestroy(SPTFQMR_CONTENT(S)->vtemp3);
      SPTFQMR_CONTENT(S)->vtemp3 = NULL;
    }
    free(S->content);
    S->content = NULL;
  }
  if (S->ops)
  {
    free(S->ops);
    S->ops = NULL;
  }
  free(S);
  S = NULL;
  return SUN_SUCCESS;
}
