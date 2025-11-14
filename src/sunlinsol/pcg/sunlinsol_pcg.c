/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds, Ashley Crawford @ UMBC
 * Based on sundials_pcg.c code, written by Daniel Reynolds @ UMBC
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
 * This is the implementation file for the PCG implementation of
 * the SUNLINSOL package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_math.h>
#include <sunlinsol/sunlinsol_pcg.h>
#include "sundials_logger_impl.h"

#include "sundials_cli.h"
#include "sundials_macros.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * PCG solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define PCG_CONTENT(S) ((SUNLinearSolverContent_PCG)(S->content))
#define PRETYPE(S)     (PCG_CONTENT(S)->pretype)
#define LASTFLAG(S)    (PCG_CONTENT(S)->last_flag)

/*
 * -----------------------------------------------------------------
 * unexported functions
 * -----------------------------------------------------------------
 */

static SUNErrCode setFromCommandLine_PCG(SUNLinearSolver S, const char* LSid,
                                         int argc, char* argv[]);

SUNErrCode SUNLinSolSetOptions_PCG(SUNLinearSolver S, const char* LSid,
                                   const char* file_name, int argc, char* argv[]);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new PCG linear solver
 */

SUNLinearSolver SUNLinSol_PCG(N_Vector y, int pretype, int maxl, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNLinearSolver S;
  SUNLinearSolverContent_PCG content;

  /* check for legal pretype and maxl values; if illegal use defaults */
  if ((pretype != SUN_PREC_NONE) && (pretype != SUN_PREC_LEFT) &&
      (pretype != SUN_PREC_RIGHT) && (pretype != SUN_PREC_BOTH))
  {
    pretype = SUN_PREC_NONE;
  }
  if (maxl <= 0) { maxl = SUNPCG_MAXL_DEFAULT; }

  /* Create linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  S->ops->gettype           = SUNLinSolGetType_PCG;
  S->ops->getid             = SUNLinSolGetID_PCG;
  S->ops->setatimes         = SUNLinSolSetATimes_PCG;
  S->ops->setoptions        = SUNLinSolSetOptions_PCG;
  S->ops->setpreconditioner = SUNLinSolSetPreconditioner_PCG;
  S->ops->setscalingvectors = SUNLinSolSetScalingVectors_PCG;
  S->ops->setzeroguess      = SUNLinSolSetZeroGuess_PCG;
  S->ops->initialize        = SUNLinSolInitialize_PCG;
  S->ops->setup             = SUNLinSolSetup_PCG;
  S->ops->solve             = SUNLinSolSolve_PCG;
  S->ops->numiters          = SUNLinSolNumIters_PCG;
  S->ops->resnorm           = SUNLinSolResNorm_PCG;
  S->ops->resid             = SUNLinSolResid_PCG;
  S->ops->lastflag          = SUNLinSolLastFlag_PCG;
  S->ops->space             = SUNLinSolSpace_PCG;
  S->ops->free              = SUNLinSolFree_PCG;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_PCG)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content  */
  S->content = content;

  /* Fill content */
  content->last_flag = 0;
  content->maxl      = maxl;
  content->pretype   = pretype;
  content->zeroguess = SUNFALSE;
  content->numiters  = 0;
  content->resnorm   = ZERO;
  content->r         = NULL;
  content->p         = NULL;
  content->z         = NULL;
  content->Ap        = NULL;
  content->s         = NULL;
  content->ATimes    = NULL;
  content->ATData    = NULL;
  content->Psetup    = NULL;
  content->Psolve    = NULL;
  content->PData     = NULL;

  /* Allocate content */
  content->r = N_VClone(y);
  SUNCheckLastErrNull();

  content->p = N_VClone(y);
  SUNCheckLastErrNull();

  content->z = N_VClone(y);
  SUNCheckLastErrNull();

  content->Ap = N_VClone(y);
  SUNCheckLastErrNull();

  return (S);
}

/* ----------------------------------------------------------------------------
 * Function to control set routines via the command line or file
 */

SUNErrCode SUNLinSolSetOptions_PCG(SUNLinearSolver S, const char* LSid,
                                   SUNDIALS_MAYBE_UNUSED const char* file_name,
                                   int argc, char* argv[])
{
  SUNFunctionBegin(S->sunctx);

  /* File-based option control is currently unimplemented */
  SUNAssert((file_name == NULL || strlen(file_name) == 0),
            SUN_ERR_ARG_INCOMPATIBLE);

  if (argc > 0 && argv != NULL)
  {
    SUNCheckCall(setFromCommandLine_PCG(S, LSid, argc, argv));
  }

  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Function to control set routines via the command line
 */

static SUNErrCode setFromCommandLine_PCG(SUNLinearSolver S, const char* LSid,
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
      retval   = SUNLinSol_PCGSetPrecType(S, iarg);
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
      retval   = SUNLinSol_PCGSetMaxl(S, iarg);
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
 * Function to set the type of preconditioning for PCG to use
 */

SUNErrCode SUNLinSol_PCGSetPrecType(SUNLinearSolver S, int pretype)
{
  SUNFunctionBegin(S->sunctx);
  /* Check for legal pretype */
  SUNAssert((pretype == SUN_PREC_NONE) || (pretype == SUN_PREC_LEFT) ||
              (pretype == SUN_PREC_RIGHT) || (pretype == SUN_PREC_BOTH),
            SUN_ERR_ARG_OUTOFRANGE);

  /* Set pretype */
  PRETYPE(S) = pretype;
  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Function to set the maximum number of iterations for PCG to use
 */

SUNErrCode SUNLinSol_PCGSetMaxl(SUNLinearSolver S, int maxl)
{
  /* Check for legal number of iters */
  if (maxl <= 0) { maxl = SUNPCG_MAXL_DEFAULT; }

  /* Set max iters */
  PCG_CONTENT(S)->maxl = maxl;
  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_PCG(SUNDIALS_MAYBE_UNUSED SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_ITERATIVE);
}

SUNLinearSolver_ID SUNLinSolGetID_PCG(SUNDIALS_MAYBE_UNUSED SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_PCG);
}

SUNErrCode SUNLinSolInitialize_PCG(SUNLinearSolver S)
{
  SUNFunctionBegin(S->sunctx);
  if (PCG_CONTENT(S)->maxl <= 0) { PCG_CONTENT(S)->maxl = SUNPCG_MAXL_DEFAULT; }

  SUNAssert(PCG_CONTENT(S)->ATimes, SUN_ERR_ARG_CORRUPT);

  if ((PRETYPE(S) != SUN_PREC_LEFT) && (PRETYPE(S) != SUN_PREC_RIGHT) &&
      (PRETYPE(S) != SUN_PREC_BOTH))
  {
    PRETYPE(S) = SUN_PREC_NONE;
  }

  SUNAssert((PCG_CONTENT(S)->pretype == SUN_PREC_NONE) ||
              (PCG_CONTENT(S)->Psolve != NULL),
            SUN_ERR_ARG_CORRUPT);

  /* no additional memory to allocate */

  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetATimes_PCG(SUNLinearSolver S, void* ATData,
                                  SUNATimesFn ATimes)
{
  /* set function pointers to integrator-supplied ATimes routine
     and data, and return with success */
  PCG_CONTENT(S)->ATimes = ATimes;
  PCG_CONTENT(S)->ATData = ATData;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetPreconditioner_PCG(SUNLinearSolver S, void* PData,
                                          SUNPSetupFn Psetup, SUNPSolveFn Psolve)
{
  /* set function pointers to integrator-supplied Psetup and PSolve
     routines and data, and return with success */
  PCG_CONTENT(S)->Psetup = Psetup;
  PCG_CONTENT(S)->Psolve = Psolve;
  PCG_CONTENT(S)->PData  = PData;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetScalingVectors_PCG(SUNLinearSolver S, N_Vector s,
                                          SUNDIALS_MAYBE_UNUSED N_Vector nul)
{
  /* set N_Vector pointer to integrator-supplied scaling vector
     (only use the first one), and return with success */
  PCG_CONTENT(S)->s = s;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetZeroGuess_PCG(SUNLinearSolver S, sunbooleantype onoff)
{
  /* set flag indicating a zero initial guess */
  PCG_CONTENT(S)->zeroguess = onoff;
  return SUN_SUCCESS;
}

int SUNLinSolSetup_PCG(SUNLinearSolver S, SUNDIALS_MAYBE_UNUSED SUNMatrix nul)
{
  SUNFunctionBegin(S->sunctx);

  int status;
  SUNPSetupFn Psetup;
  void* PData;

  /* Set shortcuts to PCG memory structures */
  Psetup = PCG_CONTENT(S)->Psetup;
  PData  = PCG_CONTENT(S)->PData;

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

int SUNLinSolSolve_PCG(SUNLinearSolver S, SUNDIALS_MAYBE_UNUSED SUNMatrix nul,
                       N_Vector x, N_Vector b, sunrealtype delta)
{
  SUNFunctionBegin(S->sunctx);

  /* local data and shortcut variables */
  sunrealtype alpha, beta, r0_norm, rho, rz, rz_old;
  N_Vector r, p, z, Ap, w;
  sunbooleantype UsePrec, UseScaling, converged;
  sunbooleantype* zeroguess;
  int l, l_max, pretype;
  void *A_data, *P_data;
  SUNATimesFn atimes;
  SUNPSolveFn psolve;
  sunrealtype* res_norm;
  int* nli;
  int status;

  /* Make local shortcuts to solver variables. */
  l_max     = PCG_CONTENT(S)->maxl;
  r         = PCG_CONTENT(S)->r;
  p         = PCG_CONTENT(S)->p;
  z         = PCG_CONTENT(S)->z;
  Ap        = PCG_CONTENT(S)->Ap;
  w         = PCG_CONTENT(S)->s;
  A_data    = PCG_CONTENT(S)->ATData;
  P_data    = PCG_CONTENT(S)->PData;
  atimes    = PCG_CONTENT(S)->ATimes;
  psolve    = PCG_CONTENT(S)->Psolve;
  pretype   = PCG_CONTENT(S)->pretype;
  zeroguess = &(PCG_CONTENT(S)->zeroguess);
  nli       = &(PCG_CONTENT(S)->numiters);
  res_norm  = &(PCG_CONTENT(S)->resnorm);

  /* Initialize counters and convergence flag */
  *nli      = 0;
  converged = SUNFALSE;

  /* set sunbooleantype flags for internal solver options */
  UsePrec    = ((pretype == SUN_PREC_BOTH) || (pretype == SUN_PREC_LEFT) ||
             (pretype == SUN_PREC_RIGHT));
  UseScaling = (w != NULL);

  /* Check if Atimes function has been set */
  SUNAssert(atimes, SUN_ERR_ARG_CORRUPT);

  /* If preconditioning, check if psolve has been set */
  SUNAssert(!UsePrec || psolve, SUN_ERR_ARG_CORRUPT);

  SUNLogInfo(S->sunctx->logger, "linear-solver", "solver = pcg");

  SUNLogInfo(S->sunctx->logger, "begin-iterations-list", "");

  /* Set r to initial residual r_0 = b - A*x_0 */
  if (*zeroguess)
  {
    N_VScale(ONE, b, r);
    SUNCheckLastErr();
  }
  else
  {
    status = atimes(A_data, x, r);
    if (status != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (status < 0) ? SUNLS_ATIMES_FAIL_UNREC
                                 : SUNLS_ATIMES_FAIL_REC;

      SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                 "status = failed matvec");

      return (LASTFLAG(S));
    }
    N_VLinearSum(ONE, b, -ONE, r, r);
    SUNCheckLastErr();
  }

  /* Set rho to scaled L2 norm of r, and return if small */
  if (UseScaling)
  {
    N_VProd(r, w, Ap);
    SUNCheckLastErr();
  }
  else
  {
    N_VScale(ONE, r, Ap);
    SUNCheckLastErr();
  }
  rho = N_VDotProd(Ap, Ap);
  SUNCheckLastErr();
  *res_norm = r0_norm = rho = SUNRsqrt(rho);

  if (rho <= delta)
  {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUN_SUCCESS;

    SUNLogInfo(S->sunctx->logger, "end-iterations-list",
               "cur-iter = 0, total-iters = 0, res-norm = " SUN_FORMAT_G
               ", status = success",
               *res_norm);

    return (LASTFLAG(S));
  }

  /* Apply preconditioner and b-scaling to r = r_0 */
  if (UsePrec)
  {
    status = psolve(P_data, r, z, delta, SUN_PREC_LEFT); /* z = P^{-1}r */
    if (status != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                 : SUNLS_PSOLVE_FAIL_REC;

      SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                 "cur-iter = 0, total-iters = 0, res-norm = " SUN_FORMAT_G
                 ", status = failed preconditioner solve",
                 *res_norm);

      return (LASTFLAG(S));
    }
  }
  else
  {
    N_VScale(ONE, r, z);
    SUNCheckLastErr();
  }

  /* Initialize rz to <r,z> */
  rz = N_VDotProd(r, z);
  SUNCheckLastErr();

  /* Copy z to p */
  N_VScale(ONE, z, p);
  SUNCheckLastErr();

  SUNLogInfo(S->sunctx->logger, "end-iterations-list",
             "cur-iter = 0, total-iters = 0, res-norm = " SUN_FORMAT_G
             ", status = continue",
             *res_norm);

  /* Begin main iteration loop */
  for (l = 0; l < l_max; l++)
  {
    SUNLogInfo(S->sunctx->logger, "begin-iterations-list", "");

    /* increment counter */
    (*nli)++;

    /* Generate Ap = A*p */
    status = atimes(A_data, p, Ap);
    if (status != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (status < 0) ? SUNLS_ATIMES_FAIL_UNREC
                                 : SUNLS_ATIMES_FAIL_REC;

      SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                 "status = failed matvec");

      return (LASTFLAG(S));
    }

    /* Calculate alpha = <r,z> / <Ap,p> */
    alpha = N_VDotProd(Ap, p);
    SUNCheckLastErr();
    alpha = rz / alpha;

    /* Update x = x + alpha*p */
    if (l == 0 && *zeroguess)
    {
      N_VScale(alpha, p, x);
      SUNCheckLastErr();
    }
    else
    {
      N_VLinearSum(ONE, x, alpha, p, x);
      SUNCheckLastErr();
    }

    /* Update r = r - alpha*Ap */
    N_VLinearSum(ONE, r, -alpha, Ap, r);
    SUNCheckLastErr();

    /* Set rho and check convergence */
    if (UseScaling)
    {
      N_VProd(r, w, Ap);
      SUNCheckLastErr();
    }
    else
    {
      N_VScale(ONE, r, Ap);
      SUNCheckLastErr();
    }
    rho = N_VDotProd(Ap, Ap);
    SUNCheckLastErr();
    *res_norm = rho = SUNRsqrt(rho);

    SUNLogInfo(S->sunctx->logger, "linear-iterate",
               "cur-iter = %i, res-norm = " SUN_FORMAT_G, *nli, *res_norm);

    if (rho <= delta)
    {
      converged = SUNTRUE;
      break;
    }

    /* Exit early on last iteration */
    if (l == l_max - 1) { break; }

    /* Apply preconditioner:  z = P^{-1}*r */
    if (UsePrec)
    {
      status = psolve(P_data, r, z, delta, SUN_PREC_LEFT);
      if (status != 0)
      {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                   : SUNLS_PSOLVE_FAIL_REC;

        SUNLogInfo(S->sunctx->logger, "end-iterations-list",
                   "status = failed preconditioner solve");

        return (LASTFLAG(S));
      }
    }
    else
    {
      N_VScale(ONE, r, z);
      SUNCheckLastErr();
    }

    /* update rz */
    rz_old = rz;
    rz     = N_VDotProd(r, z);
    SUNCheckLastErr();

    /* Calculate beta = <r,z> / <r_old,z_old> */
    beta = rz / rz_old;

    /* Update p = z + beta*p */
    N_VLinearSum(ONE, z, beta, p, p);
    SUNCheckLastErr();

    SUNLogInfo(S->sunctx->logger, "end-iterations-list", "status = continue");
  }

  /* Main loop finished, return with result */
  *zeroguess = SUNFALSE;
  if (converged == SUNTRUE)
  {
    SUNLogInfo(S->sunctx->logger, "end-iterations-list", "status = success");

    LASTFLAG(S) = SUN_SUCCESS;
  }
  else if (rho < r0_norm)
  {
    SUNLogInfo(S->sunctx->logger, "end-iterations-list",
               "status = failed residual reduced");

    LASTFLAG(S) = SUNLS_RES_REDUCED;
  }
  else
  {
    SUNLogInfo(S->sunctx->logger, "end-iterations-list",
               "status = failed max iterations");

    LASTFLAG(S) = SUNLS_CONV_FAIL;
  }
  return (LASTFLAG(S));
}

int SUNLinSolNumIters_PCG(SUNLinearSolver S)
{
  /* return the stored 'numiters' value */
  return (PCG_CONTENT(S)->numiters);
}

sunrealtype SUNLinSolResNorm_PCG(SUNLinearSolver S)
{
  /* return the stored 'resnorm' value */
  return (PCG_CONTENT(S)->resnorm);
}

N_Vector SUNLinSolResid_PCG(SUNLinearSolver S)
{
  /* return the stored 'r' vector */
  return (PCG_CONTENT(S)->r);
}

sunindextype SUNLinSolLastFlag_PCG(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  return (LASTFLAG(S));
}

SUNErrCode SUNLinSolSpace_PCG(SUNLinearSolver S, long int* lenrwLS,
                              long int* leniwLS)
{
  SUNFunctionBegin(S->sunctx);
  sunindextype liw1, lrw1;
  N_VSpace(PCG_CONTENT(S)->r, &lrw1, &liw1);
  SUNCheckLastErr();
  *lenrwLS = 1 + lrw1 * 4;
  *leniwLS = 4 + liw1 * 4;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolFree_PCG(SUNLinearSolver S)
{
  if (S == NULL) { return SUN_SUCCESS; }

  if (S->content)
  {
    /* delete items from within the content structure */
    if (PCG_CONTENT(S)->r)
    {
      N_VDestroy(PCG_CONTENT(S)->r);
      PCG_CONTENT(S)->r = NULL;
    }
    if (PCG_CONTENT(S)->p)
    {
      N_VDestroy(PCG_CONTENT(S)->p);
      PCG_CONTENT(S)->p = NULL;
    }
    if (PCG_CONTENT(S)->z)
    {
      N_VDestroy(PCG_CONTENT(S)->z);
      PCG_CONTENT(S)->z = NULL;
    }
    if (PCG_CONTENT(S)->Ap)
    {
      N_VDestroy(PCG_CONTENT(S)->Ap);
      PCG_CONTENT(S)->Ap = NULL;
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
