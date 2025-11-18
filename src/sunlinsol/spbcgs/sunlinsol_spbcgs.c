/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ UMBC
 * Based on sundials_spbcgs.c code, written by Peter Brown and
 *                Aaron Collier @ LLNL
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
 * This is the implementation file for the SPBCGS implementation of
 * the SUNLINSOL package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_math.h>
#include <sunlinsol/sunlinsol_spbcgs.h>
#include "sundials_logger_impl.h"

#include "sundials_cli.h"
#include "sundials_macros.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * SPBCGS solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define SPBCGS_CONTENT(S) ((SUNLinearSolverContent_SPBCGS)(S->content))
#define PRETYPE(S)        (SPBCGS_CONTENT(S)->pretype)
#define LASTFLAG(S)       (SPBCGS_CONTENT(S)->last_flag)

/*
 * ----------------------------------------------------------------------------
 * Un-exported implementation specific routines
 * ----------------------------------------------------------------------------
 */

static SUNErrCode setFromCommandLine_SPBCGS(SUNLinearSolver S, const char* LSid,
                                            int argc, char* argv[]);

SUNErrCode SUNLinSolSetOptions_SPBCGS(SUNLinearSolver S, const char* LSid,
                                      const char* file_name, int argc,
                                      char* argv[]);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new SPBCGS linear solver
 */

SUNLinearSolver SUNLinSol_SPBCGS(N_Vector y, int pretype, int maxl,
                                 SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNLinearSolver S;
  SUNLinearSolverContent_SPBCGS content;

  /* check for legal pretype and maxl values; if illegal use defaults */
  if ((pretype != SUN_PREC_NONE) && (pretype != SUN_PREC_LEFT) &&
      (pretype != SUN_PREC_RIGHT) && (pretype != SUN_PREC_BOTH))
  {
    pretype = SUN_PREC_NONE;
  }
  if (maxl <= 0) { maxl = SUNSPBCGS_MAXL_DEFAULT; }

  /* check that the supplied N_Vector supports all requisite operations */
  SUNAssertNull((y->ops->nvclone) && (y->ops->nvdestroy) &&
                  (y->ops->nvlinearsum) && (y->ops->nvprod) &&
                  (y->ops->nvdiv) && (y->ops->nvscale) && (y->ops->nvdotprod),
                SUN_ERR_ARG_INCOMPATIBLE);

  /* Create linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  S->ops->gettype           = SUNLinSolGetType_SPBCGS;
  S->ops->getid             = SUNLinSolGetID_SPBCGS;
  S->ops->setatimes         = SUNLinSolSetATimes_SPBCGS;
  S->ops->setoptions        = SUNLinSolSetOptions_SPBCGS;
  S->ops->setpreconditioner = SUNLinSolSetPreconditioner_SPBCGS;
  S->ops->setscalingvectors = SUNLinSolSetScalingVectors_SPBCGS;
  S->ops->setzeroguess      = SUNLinSolSetZeroGuess_SPBCGS;
  S->ops->initialize        = SUNLinSolInitialize_SPBCGS;
  S->ops->setup             = SUNLinSolSetup_SPBCGS;
  S->ops->solve             = SUNLinSolSolve_SPBCGS;
  S->ops->numiters          = SUNLinSolNumIters_SPBCGS;
  S->ops->resnorm           = SUNLinSolResNorm_SPBCGS;
  S->ops->resid             = SUNLinSolResid_SPBCGS;
  S->ops->lastflag          = SUNLinSolLastFlag_SPBCGS;
  S->ops->space             = SUNLinSolSpace_SPBCGS;
  S->ops->free              = SUNLinSolFree_SPBCGS;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_SPBCGS)malloc(sizeof *content);
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
  content->r         = NULL;
  content->p         = NULL;
  content->q         = NULL;
  content->u         = NULL;
  content->Ap        = NULL;
  content->vtemp     = NULL;
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

  content->r = N_VClone(y);
  SUNCheckLastErrNull();

  content->p = N_VClone(y);
  SUNCheckLastErrNull();

  content->q = N_VClone(y);
  SUNCheckLastErrNull();

  content->u = N_VClone(y);
  SUNCheckLastErrNull();

  content->Ap = N_VClone(y);
  SUNCheckLastErrNull();

  content->vtemp = N_VClone(y);
  SUNCheckLastErrNull();

  return (S);
}

/* ----------------------------------------------------------------------------
 * Function to control set routines via the command line or file
 */

SUNErrCode SUNLinSolSetOptions_SPBCGS(SUNLinearSolver S, const char* LSid,
                                      SUNDIALS_MAYBE_UNUSED const char* file_name,
                                      int argc, char* argv[])
{
  SUNFunctionBegin(S->sunctx);

  /* File-based option control is currently unimplemented */
  SUNAssert((file_name == NULL || strlen(file_name) == 0),
            SUN_ERR_ARG_INCOMPATIBLE);

  if (argc > 0 && argv != NULL)
  {
    SUNCheckCall(setFromCommandLine_SPBCGS(S, LSid, argc, argv));
  }

  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Function to control set routines via the command line
 */

static SUNErrCode setFromCommandLine_SPBCGS(SUNLinearSolver S, const char* LSid,
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
      retval   = SUNLinSol_SPBCGSSetPrecType(S, iarg);
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
      retval   = SUNLinSol_SPBCGSSetMaxl(S, iarg);
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
 * Function to set the type of preconditioning for SPBCGS to use
 */

SUNErrCode SUNLinSol_SPBCGSSetPrecType(SUNLinearSolver S, int pretype)
{
  SUNFunctionBegin(S->sunctx);
  /* Check for legal pretype */
  SUNAssert((pretype == SUN_PREC_NONE) || (pretype == SUN_PREC_LEFT) ||
              (pretype == SUN_PREC_RIGHT) || (pretype == SUN_PREC_BOTH),
            SUN_ERR_ARG_CORRUPT);

  /* Set pretype */
  PRETYPE(S) = pretype;
  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Function to set the maximum number of iterations for SPBCGS to use
 */

SUNErrCode SUNLinSol_SPBCGSSetMaxl(SUNLinearSolver S, int maxl)
{
  SUNFunctionBegin(S->sunctx);

  /* Check for legal pretype */
  if (maxl <= 0) { maxl = SUNSPBCGS_MAXL_DEFAULT; }

  /* Set pretype */
  SPBCGS_CONTENT(S)->maxl = maxl;
  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_SPBCGS(SUNDIALS_MAYBE_UNUSED SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_ITERATIVE);
}

SUNLinearSolver_ID SUNLinSolGetID_SPBCGS(SUNDIALS_MAYBE_UNUSED SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_SPBCGS);
}

SUNErrCode SUNLinSolInitialize_SPBCGS(SUNLinearSolver S)
{
  SUNFunctionBegin(S->sunctx);

  if (SPBCGS_CONTENT(S)->maxl <= 0)
  {
    SPBCGS_CONTENT(S)->maxl = SUNSPBCGS_MAXL_DEFAULT;
  }

  SUNAssert(SPBCGS_CONTENT(S)->ATimes, SUN_ERR_ARG_CORRUPT);

  if ((PRETYPE(S) != SUN_PREC_LEFT) && (PRETYPE(S) != SUN_PREC_RIGHT) &&
      (PRETYPE(S) != SUN_PREC_BOTH))
  {
    PRETYPE(S) = SUN_PREC_NONE;
  }

  SUNAssert((PRETYPE(S) == SUN_PREC_NONE) || SPBCGS_CONTENT(S)->Psolve,
            SUN_ERR_ARG_CORRUPT);

  /* no additional memory to allocate */

  /* return with success */
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetATimes_SPBCGS(SUNLinearSolver S, void* ATData,
                                     SUNATimesFn ATimes)
{
  SUNFunctionBegin(S->sunctx);
  /* set function pointers to integrator-supplied ATimes routine
     and data, and return with success */
  SPBCGS_CONTENT(S)->ATimes = ATimes;
  SPBCGS_CONTENT(S)->ATData = ATData;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetPreconditioner_SPBCGS(SUNLinearSolver S, void* PData,
                                             SUNPSetupFn Psetup,
                                             SUNPSolveFn Psolve)
{
  SUNFunctionBegin(S->sunctx);
  /* set function pointers to integrator-supplied Psetup and PSolve
     routines and data, and return with success */
  SPBCGS_CONTENT(S)->Psetup = Psetup;
  SPBCGS_CONTENT(S)->Psolve = Psolve;
  SPBCGS_CONTENT(S)->PData  = PData;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetScalingVectors_SPBCGS(SUNLinearSolver S, N_Vector s1,
                                             N_Vector s2)
{
  SUNFunctionBegin(S->sunctx);
  /* set N_Vector pointers to integrator-supplied scaling vectors,
     and return with success */
  SPBCGS_CONTENT(S)->s1 = s1;
  SPBCGS_CONTENT(S)->s2 = s2;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetZeroGuess_SPBCGS(SUNLinearSolver S, sunbooleantype onoff)
{
  SUNFunctionBegin(S->sunctx);
  /* set flag indicating a zero initial guess */
  SPBCGS_CONTENT(S)->zeroguess = onoff;
  return SUN_SUCCESS;
}

int SUNLinSolSetup_SPBCGS(SUNLinearSolver S, SUNDIALS_MAYBE_UNUSED SUNMatrix A)
{
  SUNFunctionBegin(S->sunctx);

  int status;
  SUNPSetupFn Psetup;
  void* PData;

  /* Set shortcuts to SPBCGS memory structures */
  Psetup = SPBCGS_CONTENT(S)->Psetup;
  PData  = SPBCGS_CONTENT(S)->PData;

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
  return (LASTFLAG(S));
}

int SUNLinSolSolve_SPBCGS(SUNLinearSolver S, SUNDIALS_MAYBE_UNUSED SUNMatrix A,
                          N_Vector x, N_Vector b, sunrealtype delta)
{
  SUNFunctionBegin(S->sunctx);

  /* local data and shortcut variables */
  sunrealtype alpha, beta, omega, omega_denom, beta_num, beta_denom, r_norm, rho;
  N_Vector r_star, r, p, q, u, Ap, vtemp;
  sunbooleantype preOnLeft, preOnRight, scale_x, scale_b, converged;
  sunbooleantype* zeroguess;
  int l, l_max;
  void *A_data, *P_data;
  N_Vector sx, sb;
  SUNATimesFn atimes;
  SUNPSolveFn psolve;
  sunrealtype* res_norm;
  int* nli;
  int status;

  /* local variables for fused vector operations */
  sunrealtype cv[3];
  N_Vector Xv[3];

  /* Make local shortcuts to solver variables. */
  l_max     = SPBCGS_CONTENT(S)->maxl;
  r_star    = SPBCGS_CONTENT(S)->r_star;
  r         = SPBCGS_CONTENT(S)->r;
  p         = SPBCGS_CONTENT(S)->p;
  q         = SPBCGS_CONTENT(S)->q;
  u         = SPBCGS_CONTENT(S)->u;
  Ap        = SPBCGS_CONTENT(S)->Ap;
  vtemp     = SPBCGS_CONTENT(S)->vtemp;
  sb        = SPBCGS_CONTENT(S)->s1;
  sx        = SPBCGS_CONTENT(S)->s2;
  A_data    = SPBCGS_CONTENT(S)->ATData;
  P_data    = SPBCGS_CONTENT(S)->PData;
  atimes    = SPBCGS_CONTENT(S)->ATimes;
  psolve    = SPBCGS_CONTENT(S)->Psolve;
  zeroguess = &(SPBCGS_CONTENT(S)->zeroguess);
  nli       = &(SPBCGS_CONTENT(S)->numiters);
  res_norm  = &(SPBCGS_CONTENT(S)->resnorm);

  /* Initialize counters and convergence flag */
  *nli      = 0;
  converged = SUNFALSE;

  /* set sunbooleantype flags for internal solver options */
  preOnLeft = ((PRETYPE(S) == SUN_PREC_LEFT) || (PRETYPE(S) == SUN_PREC_BOTH));
  preOnRight = ((PRETYPE(S) == SUN_PREC_RIGHT) || (PRETYPE(S) == SUN_PREC_BOTH));
  scale_x = (sx != NULL);
  scale_b = (sb != NULL);

  SUNLogInfo(S->sunctx->logger, "linear-solver", "solver = spbcgs");

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

  /* Set r_star to initial (unscaled) residual r_0 = b - A*x_0 */

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

  /* Apply left preconditioner and b-scaling to r_star = r_0 */

  if (preOnLeft)
  {
    status = psolve(P_data, r_star, r, delta, SUN_PREC_LEFT);
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
    N_VScale(ONE, r_star, r);
    SUNCheckLastErr();
  }

  if (scale_b)
  {
    N_VProd(sb, r, r_star);
    SUNCheckLastErr();
  }
  else
  {
    N_VScale(ONE, r, r_star);
    SUNCheckLastErr();
  }

  /* Initialize beta_denom to the dot product of r0 with r0 */

  beta_denom = N_VDotProd(r_star, r_star);
  SUNCheckLastErr();

  /* Set r_norm to L2 norm of r_star = sb P1_inv r_0, and
     return if small */

  *res_norm = r_norm = rho = SUNRsqrt(beta_denom);

  if (r_norm <= delta)
  {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUN_SUCCESS;

    SUNLogInfo(S->sunctx->logger, "end-iterations-list",
               "cur-iter = 0, res-norm = " SUN_FORMAT_G ", status = success",
               *res_norm);

    return (LASTFLAG(S));
  }

  SUNLogInfo(S->sunctx->logger, "end-iterations-list",
             "cur-iter = 0, res-norm = " SUN_FORMAT_G ", status = continue",
             *res_norm);

  /* Copy r_star to r and p */

  N_VScale(ONE, r_star, r);
  SUNCheckLastErr();
  N_VScale(ONE, r_star, p);
  SUNCheckLastErr();

  /* Set x = sx x if non-zero guess */
  if (scale_x && !(*zeroguess))
  {
    N_VProd(sx, x, x);
    SUNCheckLastErr();
  }

  /* Begin main iteration loop */

  for (l = 0; l < l_max; l++)
  {
    SUNLogInfo(S->sunctx->logger, "begin-iterations-list", "");

    (*nli)++;

    /* Generate Ap = A-tilde p, where A-tilde = sb P1_inv A P2_inv sx_inv */

    /*   Apply x-scaling: vtemp = sx_inv p */

    if (scale_x)
    {
      N_VDiv(p, sx, vtemp);
      SUNCheckLastErr();
    }
    else
    {
      N_VScale(ONE, p, vtemp);
      SUNCheckLastErr();
    }

    /*   Apply right preconditioner: vtemp = P2_inv sx_inv p */

    if (preOnRight)
    {
      N_VScale(ONE, vtemp, Ap);
      SUNCheckLastErr();
      status = psolve(P_data, Ap, vtemp, delta, SUN_PREC_RIGHT);
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

    /*   Apply A: Ap = A P2_inv sx_inv p */

    status = atimes(A_data, vtemp, Ap);
    if (status != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (status < 0) ? SUNLS_ATIMES_FAIL_UNREC
                                 : SUNLS_ATIMES_FAIL_REC;

      SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                 "status = failed matvec, retval = ", status);

      return (LASTFLAG(S));
    }

    /*   Apply left preconditioner: vtemp = P1_inv A P2_inv sx_inv p */

    if (preOnLeft)
    {
      status = psolve(P_data, Ap, vtemp, delta, SUN_PREC_LEFT);
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
      N_VScale(ONE, Ap, vtemp);
      SUNCheckLastErr();
    }

    /*   Apply b-scaling: Ap = sb P1_inv A P2_inv sx_inv p */

    if (scale_b)
    {
      N_VProd(sb, vtemp, Ap);
      SUNCheckLastErr();
    }
    else
    {
      N_VScale(ONE, vtemp, Ap);
      SUNCheckLastErr();
    }

    /* Calculate alpha = <r,r_star>/<Ap,r_star> */

    alpha = N_VDotProd(Ap, r_star);
    SUNCheckLastErr();
    alpha = beta_denom / alpha;

    /* Update q = r - alpha*Ap = r - alpha*(sb P1_inv A P2_inv sx_inv p) */

    N_VLinearSum(ONE, r, -alpha, Ap, q);
    SUNCheckLastErr();

    /* Generate u = A-tilde q */

    /*   Apply x-scaling: vtemp = sx_inv q */

    if (scale_x)
    {
      N_VDiv(q, sx, vtemp);
      SUNCheckLastErr();
    }
    else
    {
      N_VScale(ONE, q, vtemp);
      SUNCheckLastErr();
    }

    /*   Apply right preconditioner: vtemp = P2_inv sx_inv q */

    if (preOnRight)
    {
      N_VScale(ONE, vtemp, u);
      SUNCheckLastErr();
      status = psolve(P_data, u, vtemp, delta, SUN_PREC_RIGHT);
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

    /*   Apply A: u = A P2_inv sx_inv u */

    status = atimes(A_data, vtemp, u);
    if (status != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (status < 0) ? SUNLS_ATIMES_FAIL_UNREC
                                 : SUNLS_ATIMES_FAIL_REC;

      SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                 "status = failed matvec, retval = %d", status);

      return (LASTFLAG(S));
    }

    /*   Apply left preconditioner: vtemp = P1_inv A P2_inv sx_inv p */

    if (preOnLeft)
    {
      status = psolve(P_data, u, vtemp, delta, SUN_PREC_LEFT);
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
      N_VScale(ONE, u, vtemp);
      SUNCheckLastErr();
    }

    /*   Apply b-scaling: u = sb P1_inv A P2_inv sx_inv u */

    if (scale_b)
    {
      N_VProd(sb, vtemp, u);
      SUNCheckLastErr();
    }
    else
    {
      N_VScale(ONE, vtemp, u);
      SUNCheckLastErr();
    }

    /* Calculate omega = <u,q>/<u,u> */

    omega_denom = N_VDotProd(u, u);
    SUNCheckLastErr();
    if (omega_denom == ZERO) { omega_denom = ONE; }
    omega = N_VDotProd(u, q);
    SUNCheckLastErr();
    omega /= omega_denom;

    /* Update x = x + alpha*p + omega*q */
    if (l == 0 && *zeroguess)
    {
      N_VLinearSum(alpha, p, omega, q, x);
      SUNCheckLastErr();
    }
    else
    {
      cv[0] = ONE;
      Xv[0] = x;

      cv[1] = alpha;
      Xv[1] = p;

      cv[2] = omega;
      Xv[2] = q;

      SUNCheckCall(N_VLinearCombination(3, cv, Xv, x));
    }

    /* Update the residual r = q - omega*u */

    N_VLinearSum(ONE, q, -omega, u, r);
    SUNCheckLastErr();

    /* Set rho = norm(r) and check convergence */

    *res_norm = rho = SUNRsqrt(N_VDotProd(r, r));
    SUNCheckLastErr();

    SUNLogInfo(S->sunctx->logger, "linear-iterate",
               "cur-iter = %i, res-norm = " SUN_FORMAT_G, *nli, *res_norm);

    if (rho <= delta)
    {
      converged = SUNTRUE;
      break;
    }

    /* Not yet converged, continue iteration */
    /* Update beta = <rnew,r_star> / <rold,r_start> * alpha / omega */

    beta_num = N_VDotProd(r, r_star);
    SUNCheckLastErr();
    beta = ((beta_num / beta_denom) * (alpha / omega));

    /* Update p = r + beta*(p - omega*Ap) = beta*p - beta*omega*Ap + r */
    cv[0] = beta;
    Xv[0] = p;

    cv[1] = -alpha * (beta_num / beta_denom);
    Xv[1] = Ap;

    cv[2] = ONE;
    Xv[2] = r;

    SUNCheckCall(N_VLinearCombination(3, cv, Xv, p));

    /* update beta_denom for next iteration */
    beta_denom = beta_num;

    SUNLogInfoIf(l < l_max - 1, S->sunctx->logger, "end-iterations-list",
                 "status = continue");
  }

  /* Main loop finished */

  if ((converged == SUNTRUE) || (rho < r_norm))
  {
    /* Apply the x-scaling and right preconditioner: x = P2_inv sx_inv x */

    if (scale_x)
    {
      N_VDiv(x, sx, x);
      SUNCheckLastErr();
    }
    if (preOnRight)
    {
      status = psolve(P_data, x, vtemp, delta, SUN_PREC_RIGHT);
      if (status != 0)
      {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                   : SUNLS_PSOLVE_FAIL_REC;

        SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                   "status = failed preconditioner solve, retval = %d", status);

        return (LASTFLAG(S));
      }
      N_VScale(ONE, vtemp, x);
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

int SUNLinSolNumIters_SPBCGS(SUNLinearSolver S)
{
  /* return the stored 'numiters' value */
  return (SPBCGS_CONTENT(S)->numiters);
}

sunrealtype SUNLinSolResNorm_SPBCGS(SUNLinearSolver S)
{
  /* return the stored 'resnorm' value */
  return (SPBCGS_CONTENT(S)->resnorm);
}

N_Vector SUNLinSolResid_SPBCGS(SUNLinearSolver S)
{
  /* return the stored 'r' vector */
  return (SPBCGS_CONTENT(S)->r);
}

sunindextype SUNLinSolLastFlag_SPBCGS(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  return (LASTFLAG(S));
}

SUNErrCode SUNLinSolSpace_SPBCGS(SUNLinearSolver S, long int* lenrwLS,
                                 long int* leniwLS)
{
  SUNFunctionBegin(S->sunctx);
  sunindextype liw1, lrw1;
  if (SPBCGS_CONTENT(S)->vtemp->ops->nvspace)
  {
    N_VSpace(SPBCGS_CONTENT(S)->vtemp, &lrw1, &liw1);
    SUNCheckLastErr();
  }
  else { lrw1 = liw1 = 0; }
  *lenrwLS = lrw1 * 9;
  *leniwLS = liw1 * 9;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolFree_SPBCGS(SUNLinearSolver S)
{
  if (S->content)
  {
    /* delete items from within the content structure */
    if (SPBCGS_CONTENT(S)->r_star)
    {
      N_VDestroy(SPBCGS_CONTENT(S)->r_star);
      SPBCGS_CONTENT(S)->r_star = NULL;
    }
    if (SPBCGS_CONTENT(S)->r)
    {
      N_VDestroy(SPBCGS_CONTENT(S)->r);
      SPBCGS_CONTENT(S)->r = NULL;
    }
    if (SPBCGS_CONTENT(S)->p)
    {
      N_VDestroy(SPBCGS_CONTENT(S)->p);
      SPBCGS_CONTENT(S)->p = NULL;
    }
    if (SPBCGS_CONTENT(S)->q)
    {
      N_VDestroy(SPBCGS_CONTENT(S)->q);
      SPBCGS_CONTENT(S)->q = NULL;
    }
    if (SPBCGS_CONTENT(S)->u)
    {
      N_VDestroy(SPBCGS_CONTENT(S)->u);
      SPBCGS_CONTENT(S)->u = NULL;
    }
    if (SPBCGS_CONTENT(S)->Ap)
    {
      N_VDestroy(SPBCGS_CONTENT(S)->Ap);
      SPBCGS_CONTENT(S)->Ap = NULL;
    }
    if (SPBCGS_CONTENT(S)->vtemp)
    {
      N_VDestroy(SPBCGS_CONTENT(S)->vtemp);
      SPBCGS_CONTENT(S)->vtemp = NULL;
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
