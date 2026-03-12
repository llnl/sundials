/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ UMBC
 * Based on sundials_spfgmr.c code, written by Daniel R. Reynolds
 *                and Hilari C. Tiedeman @ UMBC
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
 * This is the implementation file for the SPFGMR implementation of
 * the SUNLINSOL package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_math.h>
#include <sunlinsol/sunlinsol_spfgmr.h>
#include "sundials_logger_impl.h"

#include "sundials_cli.h"
#include "sundials_macros.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * SPFGMR solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define SPFGMR_CONTENT(S) ((SUNLinearSolverContent_SPFGMR)(S->content))
#define LASTFLAG(S)       (SPFGMR_CONTENT(S)->last_flag)

/*
 * ----------------------------------------------------------------------------
 * Un-exported implementation specific routines
 * ----------------------------------------------------------------------------
 */

static SUNErrCode setFromCommandLine_SPFGMR(SUNLinearSolver S, const char* LSid,
                                            int argc, char* argv[]);

SUNErrCode SUNLinSolSetOptions_SPFGMR(SUNLinearSolver S, const char* LSid,
                                      const char* file_name, int argc,
                                      char* argv[]);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new SPFGMR linear solver
 */

SUNLinearSolver SUNLinSol_SPFGMR(N_Vector y, int pretype, int maxl,
                                 SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNLinearSolver S;
  SUNLinearSolverContent_SPFGMR content;

  /* set preconditioning flag (enabling any preconditioner implies right
     preconditioning, since SPFGMR does not support left preconditioning) */
  pretype = ((pretype == SUN_PREC_LEFT) || (pretype == SUN_PREC_RIGHT) ||
             (pretype == SUN_PREC_BOTH))
              ? SUN_PREC_RIGHT
              : SUN_PREC_NONE;

  /* if maxl input is illegal, set to default */
  if (maxl <= 0) { maxl = SUNSPFGMR_MAXL_DEFAULT; }

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
  S->ops->gettype           = SUNLinSolGetType_SPFGMR;
  S->ops->getid             = SUNLinSolGetID_SPFGMR;
  S->ops->setatimes         = SUNLinSolSetATimes_SPFGMR;
  S->ops->setoptions        = SUNLinSolSetOptions_SPFGMR;
  S->ops->setpreconditioner = SUNLinSolSetPreconditioner_SPFGMR;
  S->ops->setscalingvectors = SUNLinSolSetScalingVectors_SPFGMR;
  S->ops->setzeroguess      = SUNLinSolSetZeroGuess_SPFGMR;
  S->ops->initialize        = SUNLinSolInitialize_SPFGMR;
  S->ops->setup             = SUNLinSolSetup_SPFGMR;
  S->ops->solve             = SUNLinSolSolve_SPFGMR;
  S->ops->numiters          = SUNLinSolNumIters_SPFGMR;
  S->ops->resnorm           = SUNLinSolResNorm_SPFGMR;
  S->ops->resid             = SUNLinSolResid_SPFGMR;
  S->ops->lastflag          = SUNLinSolLastFlag_SPFGMR;
  S->ops->space             = SUNLinSolSpace_SPFGMR;
  S->ops->free              = SUNLinSolFree_SPFGMR;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_SPFGMR)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  S->content = content;

  /* Fill content */
  content->last_flag    = 0;
  content->maxl         = maxl;
  content->pretype      = pretype;
  content->gstype       = SUNSPFGMR_GSTYPE_DEFAULT;
  content->max_restarts = SUNSPFGMR_MAXRS_DEFAULT;
  content->zeroguess    = SUNFALSE;
  content->numiters     = 0;
  content->resnorm      = ZERO;
  content->xcor         = NULL;
  content->vtemp        = NULL;
  content->s1           = NULL;
  content->s2           = NULL;
  content->ATimes       = NULL;
  content->ATData       = NULL;
  content->Psetup       = NULL;
  content->Psolve       = NULL;
  content->PData        = NULL;
  content->V            = NULL;
  content->Z            = NULL;
  content->Hes          = NULL;
  content->givens       = NULL;
  content->yg           = NULL;
  content->cv           = NULL;
  content->Xv           = NULL;

  /* Allocate content */
  content->xcor = N_VClone(y);
  SUNCheckLastErrNull();
  content->vtemp = N_VClone(y);
  SUNCheckLastErrNull();

  return (S);
}

/* ----------------------------------------------------------------------------
 * Function to control set routines via the command line or file
 */

SUNErrCode SUNLinSolSetOptions_SPFGMR(SUNLinearSolver S, const char* LSid,
                                      SUNDIALS_MAYBE_UNUSED const char* file_name,
                                      int argc, char* argv[])
{
  SUNFunctionBegin(S->sunctx);

  /* File-based option control is currently unimplemented */
  SUNAssert((file_name == NULL || strlen(file_name) == 0),
            SUN_ERR_ARG_INCOMPATIBLE);

  if (argc > 0 && argv != NULL)
  {
    SUNCheckCall(setFromCommandLine_SPFGMR(S, LSid, argc, argv));
  }

  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Function to control set routines via the command line
 */

static SUNErrCode setFromCommandLine_SPFGMR(SUNLinearSolver S, const char* LSid,
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
      retval   = SUNLinSol_SPFGMRSetPrecType(S, iarg);
      if (retval != SUN_SUCCESS)
      {
        free(prefix);
        return retval;
      }
      continue;
    }

    /* control over GSType function */
    if (strcmp(argv[idx] + offset, "gs_type") == 0)
    {
      idx += 1;
      int iarg = atoi(argv[idx]);
      retval   = SUNLinSol_SPFGMRSetGSType(S, iarg);
      if (retval != SUN_SUCCESS)
      {
        free(prefix);
        return retval;
      }
      continue;
    }

    /* control over MaxRestarts function */
    if (strcmp(argv[idx] + offset, "max_restarts") == 0)
    {
      idx += 1;
      int iarg = atoi(argv[idx]);
      retval   = SUNLinSol_SPFGMRSetMaxRestarts(S, iarg);
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
 * Function to toggle preconditioning on/off -- turns on if pretype is any
 * one of SUN_PREC_LEFT, SUN_PREC_RIGHT or SUN_PREC_BOTH; otherwise turns off
 */

SUNErrCode SUNLinSol_SPFGMRSetPrecType(SUNLinearSolver S, int pretype)
{
  /* Check for legal pretype */
  pretype = ((pretype == SUN_PREC_LEFT) || (pretype == SUN_PREC_RIGHT) ||
             (pretype == SUN_PREC_BOTH))
              ? SUN_PREC_RIGHT
              : SUN_PREC_NONE;

  /* Set pretype */
  SPFGMR_CONTENT(S)->pretype = pretype;
  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Function to set the type of Gram-Schmidt orthogonalization for SPFGMR to use
 */

SUNErrCode SUNLinSol_SPFGMRSetGSType(SUNLinearSolver S, int gstype)
{
  SUNFunctionBegin(S->sunctx);

  /* Check for legal gstype */
  SUNAssert(gstype == SUN_MODIFIED_GS || gstype == SUN_CLASSICAL_GS,
            SUN_ERR_ARG_OUTOFRANGE);

  /* Set pretype */
  SPFGMR_CONTENT(S)->gstype = gstype;
  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Function to set the maximum number of FGMRES restarts to allow
 */

SUNErrCode SUNLinSol_SPFGMRSetMaxRestarts(SUNLinearSolver S, int maxrs)
{
  /* Illegal maxrs implies use of default value */
  if (maxrs < 0) { maxrs = SUNSPFGMR_MAXRS_DEFAULT; }

  /* Set max_restarts */
  SPFGMR_CONTENT(S)->max_restarts = maxrs;
  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_SPFGMR(SUNDIALS_MAYBE_UNUSED SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_ITERATIVE);
}

SUNLinearSolver_ID SUNLinSolGetID_SPFGMR(SUNDIALS_MAYBE_UNUSED SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_SPFGMR);
}

SUNErrCode SUNLinSolInitialize_SPFGMR(SUNLinearSolver S)
{
  SUNFunctionBegin(S->sunctx);
  int k;
  SUNLinearSolverContent_SPFGMR content;

  /* set shortcut to SPFGMR memory structure */
  content = SPFGMR_CONTENT(S);

  /* ensure valid options */
  if (content->max_restarts < 0)
  {
    content->max_restarts = SUNSPFGMR_MAXRS_DEFAULT;
  }

  SUNAssert(content->ATimes, SUN_ERR_ARG_CORRUPT);

  if ((content->pretype != SUN_PREC_LEFT) &&
      (content->pretype != SUN_PREC_RIGHT) && (content->pretype != SUN_PREC_BOTH))
  {
    content->pretype = SUN_PREC_NONE;
  }

  SUNAssert((content->pretype == SUN_PREC_NONE) || (content->Psolve != NULL),
            SUN_ERR_ARG_CORRUPT);

  /* allocate solver-specific memory (where the size depends on the
     choice of maxl) here */

  /*   Krylov subspace vectors */
  if (content->V == NULL)
  {
    content->V = N_VCloneVectorArray(content->maxl + 1, content->vtemp);
    SUNCheckLastErr();
  }

  /*   Preconditioned basis vectors */
  if (content->Z == NULL)
  {
    content->Z = N_VCloneVectorArray(content->maxl + 1, content->vtemp);
    SUNCheckLastErr();
  }

  /*   Hessenberg matrix Hes */
  if (content->Hes == NULL)
  {
    content->Hes =
      (sunrealtype**)malloc((content->maxl + 1) * sizeof(sunrealtype*));
    SUNAssert(content->Hes, SUN_ERR_MALLOC_FAIL);

    for (k = 0; k <= content->maxl; k++)
    {
      content->Hes[k] = NULL;
      content->Hes[k] = (sunrealtype*)malloc(content->maxl * sizeof(sunrealtype));
      SUNAssert(content->Hes[k], SUN_ERR_MALLOC_FAIL);
    }
  }

  /*   Givens rotation components */
  if (content->givens == NULL)
  {
    content->givens =
      (sunrealtype*)malloc(2 * content->maxl * sizeof(sunrealtype));
    SUNAssert(content->givens, SUN_ERR_MALLOC_FAIL);
  }

  /*    y and g vectors */
  if (content->yg == NULL)
  {
    content->yg = (sunrealtype*)malloc((content->maxl + 1) * sizeof(sunrealtype));
    SUNAssert(content->yg, SUN_ERR_MALLOC_FAIL);
  }

  /*    cv vector for fused vector ops */
  if (content->cv == NULL)
  {
    content->cv = (sunrealtype*)malloc((content->maxl + 1) * sizeof(sunrealtype));
    SUNAssert(content->cv, SUN_ERR_MALLOC_FAIL);
  }

  /*    Xv vector for fused vector ops */
  if (content->Xv == NULL)
  {
    content->Xv = (N_Vector*)malloc((content->maxl + 1) * sizeof(N_Vector));
    SUNAssert(content->Xv, SUN_ERR_MALLOC_FAIL);
  }

  /* return with success */
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetATimes_SPFGMR(SUNLinearSolver S, void* ATData,
                                     SUNATimesFn ATimes)
{
  /* set function pointers to integrator-supplied ATimes routine
     and data, and return with success */
  SPFGMR_CONTENT(S)->ATimes = ATimes;
  SPFGMR_CONTENT(S)->ATData = ATData;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetPreconditioner_SPFGMR(SUNLinearSolver S, void* PData,
                                             SUNPSetupFn Psetup,
                                             SUNPSolveFn Psolve)
{
  /* set function pointers to integrator-supplied Psetup and PSolve
     routines and data, and return with success */
  SPFGMR_CONTENT(S)->Psetup = Psetup;
  SPFGMR_CONTENT(S)->Psolve = Psolve;
  SPFGMR_CONTENT(S)->PData  = PData;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetScalingVectors_SPFGMR(SUNLinearSolver S, N_Vector s1,
                                             N_Vector s2)
{
  /* set N_Vector pointers to integrator-supplied scaling vectors,
     and return with success */
  SPFGMR_CONTENT(S)->s1 = s1;
  SPFGMR_CONTENT(S)->s2 = s2;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetZeroGuess_SPFGMR(SUNLinearSolver S, sunbooleantype onoff)
{
  /* set flag indicating a zero initial guess */
  SPFGMR_CONTENT(S)->zeroguess = onoff;
  return SUN_SUCCESS;
}

int SUNLinSolSetup_SPFGMR(SUNLinearSolver S, SUNDIALS_MAYBE_UNUSED SUNMatrix A)
{
  SUNFunctionBegin(S->sunctx);

  int status;
  SUNPSetupFn Psetup;
  void* PData;

  /* Set shortcuts to SPFGMR memory structures */
  Psetup = SPFGMR_CONTENT(S)->Psetup;
  PData  = SPFGMR_CONTENT(S)->PData;

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

int SUNLinSolSolve_SPFGMR(SUNLinearSolver S, SUNDIALS_MAYBE_UNUSED SUNMatrix A,
                          N_Vector x, N_Vector b, sunrealtype delta)
{
  SUNFunctionBegin(S->sunctx);

  /* local data and shortcut variables */
  N_Vector *V, *Z, xcor, vtemp, s1, s2;
  sunrealtype **Hes, *givens, *yg, *res_norm;
  sunrealtype beta, rotation_product, r_norm, s_product, rho;
  sunbooleantype preOnRight, scale1, scale2, converged;
  sunbooleantype* zeroguess;
  int i, j, k, l, l_max, krydim, ntries, max_restarts, gstype;
  int* nli;
  void *A_data, *P_data;
  SUNATimesFn atimes;
  SUNPSolveFn psolve;
  int status;

  /* local shortcuts for fused vector operations */
  sunrealtype* cv;
  N_Vector* Xv;

  /* Initialize some variables */
  krydim = 0;

  /* Make local shortcuts to solver variables. */
  l_max        = SPFGMR_CONTENT(S)->maxl;
  max_restarts = SPFGMR_CONTENT(S)->max_restarts;
  gstype       = SPFGMR_CONTENT(S)->gstype;
  V            = SPFGMR_CONTENT(S)->V;
  Z            = SPFGMR_CONTENT(S)->Z;
  Hes          = SPFGMR_CONTENT(S)->Hes;
  givens       = SPFGMR_CONTENT(S)->givens;
  xcor         = SPFGMR_CONTENT(S)->xcor;
  yg           = SPFGMR_CONTENT(S)->yg;
  vtemp        = SPFGMR_CONTENT(S)->vtemp;
  s1           = SPFGMR_CONTENT(S)->s1;
  s2           = SPFGMR_CONTENT(S)->s2;
  A_data       = SPFGMR_CONTENT(S)->ATData;
  P_data       = SPFGMR_CONTENT(S)->PData;
  atimes       = SPFGMR_CONTENT(S)->ATimes;
  psolve       = SPFGMR_CONTENT(S)->Psolve;
  zeroguess    = &(SPFGMR_CONTENT(S)->zeroguess);
  nli          = &(SPFGMR_CONTENT(S)->numiters);
  res_norm     = &(SPFGMR_CONTENT(S)->resnorm);
  cv           = SPFGMR_CONTENT(S)->cv;
  Xv           = SPFGMR_CONTENT(S)->Xv;

  /* Initialize counters and convergence flag */
  *nli      = 0;
  converged = SUNFALSE;

  /* Set sunbooleantype flags for internal solver options */
  preOnRight = ((SPFGMR_CONTENT(S)->pretype == SUN_PREC_LEFT) ||
                (SPFGMR_CONTENT(S)->pretype == SUN_PREC_RIGHT) ||
                (SPFGMR_CONTENT(S)->pretype == SUN_PREC_BOTH));
  scale1     = (s1 != NULL);
  scale2     = (s2 != NULL);

  /* Check if Atimes function has been set */
  SUNAssert(atimes, SUN_ERR_ARG_CORRUPT);

  /* If preconditioning, check if psolve has been set */
  SUNAssert(!preOnRight || psolve, SUN_ERR_ARG_CORRUPT);

  SUNLogInfo(S->sunctx->logger, "linear-solver", "solver = spfgmr");

  SUNLogInfo(S->sunctx->logger, "begin-iterations-list", "");

  /* Set vtemp and V[0] to initial (unscaled) residual r_0 = b - A*x_0 */
  if (*zeroguess)
  {
    N_VScale(ONE, b, vtemp);
    SUNCheckLastErr();
  }
  else
  {
    status = atimes(A_data, x, vtemp);
    if (status != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (status < 0) ? SUNLS_ATIMES_FAIL_UNREC
                                 : SUNLS_ATIMES_FAIL_REC;

      SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                 "status = failed matvec, retval = %d", status);

      return (LASTFLAG(S));
    }
    N_VLinearSum(ONE, b, -ONE, vtemp, vtemp);
    SUNCheckLastErr();
  }

  /* Apply left scaling to vtemp = r_0 to fill V[0]. */
  if (scale1)
  {
    N_VProd(s1, vtemp, V[0]);
    SUNCheckLastErr();
  }
  else
  {
    N_VScale(ONE, vtemp, V[0]);
    SUNCheckLastErr();
  }

  /* Set r_norm = beta to L2 norm of V[0] = s1 r_0, and return if small */
  r_norm = N_VDotProd(V[0], V[0]);
  SUNCheckLastErr();
  *res_norm = r_norm = beta = SUNRsqrt(r_norm);

  if (r_norm <= delta)
  {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUN_SUCCESS;

    SUNLogInfo(S->sunctx->logger, "end-iterations-list",
               "cur-iter = 0, total-iters = 0, res-norm = " SUN_FORMAT_G
               ", status = success",
               *res_norm);

    return (LASTFLAG(S));
  }

  SUNLogInfo(S->sunctx->logger, "end-iterations-list",
             "cur-iter = 0, total-iters = 0, res-norm = " SUN_FORMAT_G
             ", status = continue",
             *res_norm);

  /* Initialize rho to avoid compiler warning message */
  rho = beta;

  /* Set xcor = 0. */
  N_VConst(ZERO, xcor);
  SUNCheckLastErr();

  /* Begin outer iterations: up to (max_restarts + 1) attempts. */
  for (ntries = 0; ntries <= max_restarts; ntries++)
  {
    /* Initialize the Hessenberg matrix Hes and Givens rotation
       product.  Normalize the initial vector V[0].             */
    for (i = 0; i <= l_max; i++)
    {
      for (j = 0; j < l_max; j++) { Hes[i][j] = ZERO; }
    }
    rotation_product = ONE;
    N_VScale(ONE / r_norm, V[0], V[0]);
    SUNCheckLastErr();

    /* Inner loop: generate Krylov sequence and Arnoldi basis. */
    for (l = 0; l < l_max; l++)
    {
      SUNLogInfo(S->sunctx->logger, "begin-iterations-list", "");

      (*nli)++;

      krydim = l + 1;

      /* Generate A-tilde V[l], where A-tilde = s1 A P_inv s2_inv. */

      /*   Apply right scaling: vtemp = s2_inv V[l]. */
      if (scale2)
      {
        N_VDiv(V[l], s2, vtemp);
        SUNCheckLastErr();
      }
      else
      {
        N_VScale(ONE, V[l], vtemp);
        SUNCheckLastErr();
      }

      /*   Apply right preconditioner: vtemp = Z[l] = P_inv s2_inv V[l]. */
      if (preOnRight)
      {
        N_VScale(ONE, vtemp, V[l + 1]);
        SUNCheckLastErr();
        status = psolve(P_data, V[l + 1], vtemp, delta, SUN_PREC_RIGHT);
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
      N_VScale(ONE, vtemp, Z[l]);
      SUNCheckLastErr();

      /*   Apply A: V[l+1] = A P_inv s2_inv V[l]. */
      status = atimes(A_data, vtemp, V[l + 1]);
      if (status != 0)
      {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (status < 0) ? SUNLS_ATIMES_FAIL_UNREC
                                   : SUNLS_ATIMES_FAIL_REC;

        SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                   "status = failed matvec, retval = %d", status);

        return (LASTFLAG(S));
      }

      /*   Apply left scaling: V[l+1] = s1 A P_inv s2_inv V[l]. */
      if (scale1)
      {
        N_VProd(s1, V[l + 1], V[l + 1]);
        SUNCheckLastErr();
      }

      /* Orthogonalize V[l+1] against previous V[i]: V[l+1] = w_tilde. */
      if (gstype == SUN_CLASSICAL_GS)
      {
        SUNCheckCall(
          SUNClassicalGS(V, Hes, l + 1, l_max, &(Hes[l + 1][l]), cv, Xv));
      }
      else
      {
        SUNCheckCall(SUNModifiedGS(V, Hes, l + 1, l_max, &(Hes[l + 1][l])));
      }

      /* Update the QR factorization of Hes. */
      if (SUNQRfact(krydim, Hes, givens, l) != 0)
      {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = SUNLS_QRFACT_FAIL;

        SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                   "status = failed QR factorization");

        return (LASTFLAG(S));
      }

      /* Update residual norm estimate; break if convergence test passes. */
      rotation_product *= givens[2 * l + 1];
      *res_norm = rho = SUNRabs(rotation_product * r_norm);

      SUNLogInfo(S->sunctx->logger, "linear-iterate",
                 "cur-iter = %i, total-iters = %i, res-norm = " SUN_FORMAT_G,
                 l + 1, *nli, *res_norm);

      if (rho <= delta)
      {
        converged = SUNTRUE;
        break;
      }

      /* Normalize V[l+1] with norm value from the Gram-Schmidt routine. */
      N_VScale(ONE / Hes[l + 1][l], V[l + 1], V[l + 1]);
      SUNCheckLastErr();

      SUNLogInfoIf(l < l_max - 1, S->sunctx->logger, "end-iterations-list",
                   "status = continue");
    }

    /* Inner loop is done.  Compute the new correction vector xcor. */

    /*   Construct g, then solve for y. */
    yg[0] = r_norm;
    for (i = 1; i <= krydim; i++) { yg[i] = ZERO; }
    if (SUNQRsol(krydim, Hes, givens, yg) != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = SUNLS_QRSOL_FAIL;

      SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                 "status = failed QR solve");

      return (LASTFLAG(S));
    }

    /*   Add correction vector Z_l y to xcor. */
    cv[0] = ONE;
    Xv[0] = xcor;

    for (k = 0; k < krydim; k++)
    {
      cv[k + 1] = yg[k];
      Xv[k + 1] = Z[k];
    }
    SUNCheckCall(N_VLinearCombination(krydim + 1, cv, Xv, xcor));

    /* If converged, construct the final solution vector x and return. */
    if (converged)
    {
      if (*zeroguess)
      {
        N_VScale(ONE, xcor, x);
        SUNCheckLastErr();
      }
      else
      {
        N_VLinearSum(ONE, x, ONE, xcor, x);
        SUNCheckLastErr();
      }
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = SUN_SUCCESS;

      SUNLogInfo(S->sunctx->logger, "end-iterations-list", "status = success");

      return (LASTFLAG(S));
    }

    /* Not yet converged; if allowed, prepare for restart. */
    if (ntries == max_restarts) { break; }

    /* Construct last column of Q in yg. */
    s_product = ONE;
    for (i = krydim; i > 0; i--)
    {
      yg[i] = s_product * givens[2 * i - 2];
      s_product *= givens[2 * i - 1];
    }
    yg[0] = s_product;

    /* Scale r_norm and yg. */
    r_norm *= s_product;
    for (i = 0; i <= krydim; i++) { yg[i] *= r_norm; }
    r_norm = SUNRabs(r_norm);

    /* Multiply yg by V_(krydim+1) to get last residual vector; restart. */
    for (k = 0; k <= krydim; k++)
    {
      cv[k] = yg[k];
      Xv[k] = V[k];
    }
    SUNCheckCall(N_VLinearCombination(krydim + 1, cv, Xv, V[0]));

    SUNLogInfo(S->sunctx->logger, "end-iterations-list", "status = continue");
  }

  /* Failed to converge, even after allowed restarts.
     If the residual norm was reduced below its initial value, compute
     and return x anyway.  Otherwise return failure flag. */
  if (rho < beta)
  {
    if (*zeroguess)
    {
      N_VScale(ONE, xcor, x);
      SUNCheckLastErr();
    }
    else
    {
      N_VLinearSum(ONE, x, ONE, xcor, x);
      SUNCheckLastErr();
    }
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUNLS_RES_REDUCED;

    SUNLogInfo(S->sunctx->logger, "end-iterations-list",
               "status = failed residual reduced");

    return (LASTFLAG(S));
  }

  *zeroguess  = SUNFALSE;
  LASTFLAG(S) = SUNLS_CONV_FAIL;

  SUNLogInfo(S->sunctx->logger, "end-iterations-list",
             "status = failed max iterations");

  return (LASTFLAG(S));
}

int SUNLinSolNumIters_SPFGMR(SUNLinearSolver S)
{
  return (SPFGMR_CONTENT(S)->numiters);
}

sunrealtype SUNLinSolResNorm_SPFGMR(SUNLinearSolver S)
{
  return (SPFGMR_CONTENT(S)->resnorm);
}

N_Vector SUNLinSolResid_SPFGMR(SUNLinearSolver S)
{
  return (SPFGMR_CONTENT(S)->vtemp);
}

sunindextype SUNLinSolLastFlag_SPFGMR(SUNLinearSolver S)
{
  return (LASTFLAG(S));
}

SUNErrCode SUNLinSolSpace_SPFGMR(SUNLinearSolver S, long int* lenrwLS,
                                 long int* leniwLS)
{
  SUNFunctionBegin(S->sunctx);
  int maxl;
  sunindextype liw1, lrw1;
  maxl = SPFGMR_CONTENT(S)->maxl;
  if (SPFGMR_CONTENT(S)->vtemp->ops->nvspace)
  {
    N_VSpace(SPFGMR_CONTENT(S)->vtemp, &lrw1, &liw1);
    SUNCheckLastErr();
  }
  else { lrw1 = liw1 = 0; }
  *lenrwLS = lrw1 * (2 * maxl + 4) + maxl * (maxl + 5) + 2;
  *leniwLS = liw1 * (2 * maxl + 4);
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolFree_SPFGMR(SUNLinearSolver S)
{
  int k;

  if (S == NULL) { return SUN_SUCCESS; }

  if (S->content)
  {
    /* delete items from within the content structure */
    if (SPFGMR_CONTENT(S)->xcor)
    {
      N_VDestroy(SPFGMR_CONTENT(S)->xcor);
      SPFGMR_CONTENT(S)->xcor = NULL;
    }
    if (SPFGMR_CONTENT(S)->vtemp)
    {
      N_VDestroy(SPFGMR_CONTENT(S)->vtemp);
      SPFGMR_CONTENT(S)->vtemp = NULL;
    }
    if (SPFGMR_CONTENT(S)->V)
    {
      N_VDestroyVectorArray(SPFGMR_CONTENT(S)->V, SPFGMR_CONTENT(S)->maxl + 1);
      SPFGMR_CONTENT(S)->V = NULL;
    }
    if (SPFGMR_CONTENT(S)->Z)
    {
      N_VDestroyVectorArray(SPFGMR_CONTENT(S)->Z, SPFGMR_CONTENT(S)->maxl + 1);
      SPFGMR_CONTENT(S)->Z = NULL;
    }
    if (SPFGMR_CONTENT(S)->Hes)
    {
      for (k = 0; k <= SPFGMR_CONTENT(S)->maxl; k++)
      {
        if (SPFGMR_CONTENT(S)->Hes[k])
        {
          free(SPFGMR_CONTENT(S)->Hes[k]);
          SPFGMR_CONTENT(S)->Hes[k] = NULL;
        }
      }
      free(SPFGMR_CONTENT(S)->Hes);
      SPFGMR_CONTENT(S)->Hes = NULL;
    }
    if (SPFGMR_CONTENT(S)->givens)
    {
      free(SPFGMR_CONTENT(S)->givens);
      SPFGMR_CONTENT(S)->givens = NULL;
    }
    if (SPFGMR_CONTENT(S)->yg)
    {
      free(SPFGMR_CONTENT(S)->yg);
      SPFGMR_CONTENT(S)->yg = NULL;
    }
    if (SPFGMR_CONTENT(S)->cv)
    {
      free(SPFGMR_CONTENT(S)->cv);
      SPFGMR_CONTENT(S)->cv = NULL;
    }
    if (SPFGMR_CONTENT(S)->Xv)
    {
      free(SPFGMR_CONTENT(S)->Xv);
      SPFGMR_CONTENT(S)->Xv = NULL;
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
