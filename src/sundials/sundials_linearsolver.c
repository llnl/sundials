/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ UMBC
 *                David J. Gardner, Carol S. Woodward, and
 *                Slaven Peles @ LLNL
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
 * This is the implementation file for a generic SUNLINEARSOLVER
 * package.  It contains the implementation of the SUNLinearSolver
 * operations listed in sundials_linearsolver.h
 * -----------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_errors.h>
#include "sundials_logger_impl.h"

#if defined(SUNDIALS_BUILD_WITH_PROFILING)
static SUNProfiler getSUNProfiler(SUNLinearSolver S)
{
  return (S->sunctx->profiler);
}
#endif

/* Forward declaration of function used to destroy any data allocated for Python */
#if defined(SUNDIALS_ENABLE_PYTHON)
void SUNLinearSolverFunctionTable_Destroy(void* ptr);
#endif

/* internal function prototypes */
SUNErrCode sunlsSetFromCommandLine(SUNLinearSolver S, const char* LSid,
                                   int argc, char* argv[]);

/* -----------------------------------------------------------------
 * Create a new empty SUNLinearSolver object
 * ----------------------------------------------------------------- */

SUNLinearSolver SUNLinSolNewEmpty(SUNContext sunctx)
{
  SUNLinearSolver LS;
  SUNLinearSolver_Ops ops;

  if (sunctx == NULL) { return NULL; }

  SUNFunctionBegin(sunctx);

  /* create linear solver object */
  LS = NULL;
  LS = (SUNLinearSolver)malloc(sizeof *LS);
  SUNAssertNull(LS, SUN_ERR_MALLOC_FAIL);

  /* create linear solver ops structure */
  ops = NULL;
  ops = (SUNLinearSolver_Ops)malloc(sizeof *ops);
  SUNAssertNull(ops, SUN_ERR_MALLOC_FAIL);

  /* initialize operations to NULL */
  ops->gettype           = NULL;
  ops->getid             = NULL;
  ops->setatimes         = NULL;
  ops->setpreconditioner = NULL;
  ops->setscalingvectors = NULL;
  ops->setoptions        = NULL;
  ops->setzeroguess      = NULL;
  ops->initialize        = NULL;
  ops->setup             = NULL;
  ops->solve             = NULL;
  ops->numiters          = NULL;
  ops->resnorm           = NULL;
  ops->resid             = NULL;
  ops->lastflag          = NULL;
  ops->space             = NULL;
  ops->free              = NULL;

  /* attach ops and initialize content and context to NULL */
  LS->ops     = ops;
  LS->content = NULL;
  LS->python  = NULL;
  LS->sunctx  = sunctx;

  return (LS);
}

/* -----------------------------------------------------------------
 * Free a generic SUNLinearSolver (assumes content is already empty)
 * ----------------------------------------------------------------- */

void SUNLinSolFreeEmpty(SUNLinearSolver S)
{
  if (S == NULL) { return; }

  free(S->ops);
  S->ops = NULL;

#if defined(SUNDIALS_ENABLE_PYTHON)
  SUNLinearSolverFunctionTable_Destroy(S->python);
#endif
  S->python = NULL;

  free(S);
  return;
}

/* -----------------------------------------------------------------
 * internal utility routines
 * ----------------------------------------------------------------- */

SUNErrCode sunlsSetFromCommandLine(SUNLinearSolver S, const char* LSid,
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

    /* control over ZeroGuess function */
    if (strcmp(argv[idx] + offset, "zero_guess") == 0)
    {
      idx += 1;
      int iarg = atoi(argv[idx]);
      retval   = SUNLinSolSetZeroGuess(S, iarg);
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

/* -----------------------------------------------------------------
 * Functions in the 'ops' structure
 * -----------------------------------------------------------------*/

SUNLinearSolver_Type SUNLinSolGetType(SUNLinearSolver S)
{
  return (S->ops->gettype(S));
}

SUNLinearSolver_ID SUNLinSolGetID(SUNLinearSolver S)
{
  if (S->ops->getid) { return (S->ops->getid(S)); }
  else { return (SUNLINEARSOLVER_CUSTOM); }
}

SUNErrCode SUNLinSolSetATimes(SUNLinearSolver S, void* A_data, SUNATimesFn ATimes)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(S));
  if (S->ops->setatimes) { ier = S->ops->setatimes(S, A_data, ATimes); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(S));
  return (ier);
}

SUNErrCode SUNLinSolSetPreconditioner(SUNLinearSolver S, void* P_data,
                                      SUNPSetupFn Pset, SUNPSolveFn Psol)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(S));
  if (S->ops->setpreconditioner)
  {
    ier = S->ops->setpreconditioner(S, P_data, Pset, Psol);
  }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(S));
  return (ier);
}

SUNErrCode SUNLinSolSetScalingVectors(SUNLinearSolver S, N_Vector s1, N_Vector s2)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(S));
  if (S->ops->setscalingvectors) { ier = S->ops->setscalingvectors(S, s1, s2); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(S));
  return (ier);
}

SUNErrCode SUNLinSolSetOptions(SUNLinearSolver S, const char* LSid,
                               const char* file_name, int argc, char* argv[])
{
  SUNErrCode ier = SUN_SUCCESS;
  if (S == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(S->sunctx);
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(S));

  /* File-based option control is currently unimplemented */
  SUNAssert((file_name == NULL || strlen(file_name) == 0),
            SUN_ERR_ARG_INCOMPATIBLE);

  /* First, process all base-class options */
  if (argc > 0 && argv != NULL)
  {
    SUNCheckCall(sunlsSetFromCommandLine(S, LSid, argc, argv));
  }

  /* Second, ask the implementation to process any remaining options */
  if (S->ops->setoptions)
  {
    ier = S->ops->setoptions(S, LSid, file_name, argc, argv);
  }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(S));
  return (ier);
}

SUNErrCode SUNLinSolSetZeroGuess(SUNLinearSolver S, sunbooleantype onoff)
{
  if (S->ops->setzeroguess) { return (S->ops->setzeroguess(S, onoff)); }
  else { return SUN_SUCCESS; }
}

SUNErrCode SUNLinSolInitialize(SUNLinearSolver S)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(S));
  if (S->ops->initialize) { ier = S->ops->initialize(S); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(S));
  return (ier);
}

int SUNLinSolSetup(SUNLinearSolver S, SUNMatrix A)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(S));
  if (S->ops->setup) { ier = S->ops->setup(S, A); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(S));
  return (ier);
}

int SUNLinSolSolve(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b,
                   sunrealtype tol)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(S));
  ier = S->ops->solve(S, A, x, b, tol);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(S));
  return (ier);
}

int SUNLinSolNumIters(SUNLinearSolver S)
{
  int result;
  if (S->ops->numiters) { result = S->ops->numiters(S); }
  else { result = 0; }
  return (result);
}

sunrealtype SUNLinSolResNorm(SUNLinearSolver S)
{
  sunrealtype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(S));
  if (S->ops->resnorm) { result = S->ops->resnorm(S); }
  else { result = SUN_RCONST(0.0); }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(S));
  return (result);
}

N_Vector SUNLinSolResid(SUNLinearSolver S)
{
  N_Vector resid;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(S));
  if (S->ops->resid) { resid = S->ops->resid(S); }
  else { resid = NULL; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(S));
  return (resid);
}

sunindextype SUNLinSolLastFlag(SUNLinearSolver S)
{
  if (S->ops->lastflag) { return ((sunindextype)S->ops->lastflag(S)); }
  else { return 0; }
}

SUNErrCode SUNLinSolSpace(SUNLinearSolver S, long int* lenrwLS, long int* leniwLS)
{
  if (S->ops->space) { return (S->ops->space(S, lenrwLS, leniwLS)); }
  else
  {
    *lenrwLS = 0;
    *leniwLS = 0;
    return SUN_SUCCESS;
  }
}

SUNErrCode SUNLinSolFree(SUNLinearSolver S)
{
  if (S == NULL) { return SUN_SUCCESS; }

  /* if the free operation exists use it */
  if (S->ops)
  {
    if (S->ops->free) { return (S->ops->free(S)); }
  }

  /* if we reach this point, either ops == NULL or free == NULL,
     try to cleanup by freeing the content, ops, and solver */
  free(S->content);
  S->content = NULL;

  free(S->ops);
  S->ops = NULL;

#if defined(SUNDIALS_ENABLE_PYTHON)
  SUNLinearSolverFunctionTable_Destroy(S->python);
#endif
  S->python = NULL;

  free(S);
  S = NULL;

  return SUN_SUCCESS;
}
