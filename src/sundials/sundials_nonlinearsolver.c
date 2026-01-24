/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 * This is the implementation file for a generic SUNNonlinerSolver package. It
 * contains the implementation of the SUNNonlinearSolver operations listed in
 * the 'ops' structure in sundials_nonlinearsolver.h
 * ---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>

#include <sundials/priv/sundials_context_impl.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>
#include "sundials_logger_impl.h"

#if defined(SUNDIALS_BUILD_WITH_PROFILING)
static SUNProfiler getSUNProfiler(SUNNonlinearSolver NLS)
{
  return (NLS->sunctx->profiler);
}
#endif

/* Forward declaration of function used to destroy any data allocated for Python */
#if defined(SUNDIALS_ENABLE_PYTHON)
void SUNNonlinearSolverFunctionTable_Destroy(void* ptr);
#endif

/* internal function prototypes */
SUNErrCode sunnlsSetFromCommandLine(SUNNonlinearSolver NLS, const char* NLSid,
                                    int argc, char* argv[]);

/* -----------------------------------------------------------------------------
 * Create a new empty SUNLinearSolver object
 * ---------------------------------------------------------------------------*/

SUNNonlinearSolver SUNNonlinSolNewEmpty(SUNContext sunctx)
{
  if (sunctx == NULL) { return NULL; }

  SUNFunctionBegin(sunctx);
  SUNNonlinearSolver NLS;
  SUNNonlinearSolver_Ops ops;

  /* create nonlinear solver object */
  NLS = NULL;
  NLS = (SUNNonlinearSolver)malloc(sizeof *NLS);
  SUNAssertNull(NLS, SUN_ERR_MALLOC_FAIL);

  /* create nonlinear solver ops structure */
  ops = NULL;
  ops = (SUNNonlinearSolver_Ops)malloc(sizeof *ops);
  SUNAssertNull(ops, SUN_ERR_MALLOC_FAIL);

  /* initialize operations to NULL */
  ops->gettype         = NULL;
  ops->initialize      = NULL;
  ops->setup           = NULL;
  ops->solve           = NULL;
  ops->free            = NULL;
  ops->setsysfn        = NULL;
  ops->setlsetupfn     = NULL;
  ops->setlsolvefn     = NULL;
  ops->setctestfn      = NULL;
  ops->setoptions      = NULL;
  ops->setmaxiters     = NULL;
  ops->getnumiters     = NULL;
  ops->getcuriter      = NULL;
  ops->getnumconvfails = NULL;

  /* attach context and ops, initialize content to NULL */
  NLS->sunctx  = sunctx;
  NLS->ops     = ops;
  NLS->content = NULL;
  NLS->python  = NULL;

  return (NLS);
}

/* -----------------------------------------------------------------------------
 * Free a generic SUNNonlinearSolver (assumes content is already empty)
 * ---------------------------------------------------------------------------*/

void SUNNonlinSolFreeEmpty(SUNNonlinearSolver NLS)
{
  if (NLS == NULL) { return; }

  /* free non-NULL ops structure */
  free(NLS->ops);
  NLS->ops = NULL;

#if defined(SUNDIALS_ENABLE_PYTHON)
  SUNNonlinearSolverFunctionTable_Destroy(NLS->python);
#endif
  NLS->python = NULL;

  /* free overall N_Vector object and return */
  free(NLS);
  return;
}

/* -----------------------------------------------------------------------------
 * core functions
 * ---------------------------------------------------------------------------*/

SUNNonlinearSolver_Type SUNNonlinSolGetType(SUNNonlinearSolver NLS)
{
  return (NLS->ops->gettype(NLS));
}

SUNErrCode SUNNonlinSolInitialize(SUNNonlinearSolver NLS)
{
  SUNErrCode ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(NLS));
  if (NLS->ops->initialize) { ier = NLS->ops->initialize(NLS); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(NLS));
  return (ier);
}

int SUNNonlinSolSetup(SUNNonlinearSolver NLS, N_Vector y, void* mem)
{
  int ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(NLS));
  if (NLS->ops->setup) { ier = NLS->ops->setup(NLS, y, mem); }
  else { ier = SUN_SUCCESS; }
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(NLS));
  return (ier);
}

int SUNNonlinSolSolve(SUNNonlinearSolver NLS, N_Vector y0, N_Vector y, N_Vector w,
                      sunrealtype tol, sunbooleantype callLSetup, void* mem)
{
  int status;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(NLS));
  status = NLS->ops->solve(NLS, y0, y, w, tol, callLSetup, mem);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(NLS));
  return (status);
}

SUNErrCode SUNNonlinSolFree(SUNNonlinearSolver NLS)
{
  if (NLS == NULL) { return (SUN_SUCCESS); }

  /* if the free operation exists use it */
  if (NLS->ops)
  {
    if (NLS->ops->free) { return (NLS->ops->free(NLS)); }
  }

  /* if we reach this point, either ops == NULL or free == NULL,
     try to cleanup by freeing the content, ops, and solver */
  free(NLS->content);
  NLS->content = NULL;
  free(NLS->ops);
  NLS->ops = NULL;
#if defined(SUNDIALS_ENABLE_PYTHON)
  SUNNonlinearSolverFunctionTable_Destroy(NLS->python);
#endif
  NLS->python = NULL;
  free(NLS);
  NLS = NULL;

  return (SUN_SUCCESS);
}

/* -----------------------------------------------------------------------------
 * internal utility routines
 * ---------------------------------------------------------------------------*/

SUNErrCode sunnlsSetFromCommandLine(SUNNonlinearSolver NLS, const char* NLSid,
                                    int argc, char* argv[])
{
  SUNFunctionBegin(NLS->sunctx);

  /* Prefix for options to set */
  const char* default_id = "sunnonlinearsolver";
  size_t offset          = strlen(default_id) + 1;
  if (NLSid != NULL && strlen(NLSid) > 0) { offset = strlen(NLSid) + 1; }
  char* prefix = (char*)malloc(sizeof(char) * (offset + 1));
  if (NLSid != NULL && strlen(NLSid) > 0) { strcpy(prefix, NLSid); }
  else { strcpy(prefix, default_id); }
  strcat(prefix, ".");

  for (int idx = 1; idx < argc; idx++)
  {
    int retval;

    /* skip command-line arguments that do not begin with correct prefix */
    if (strncmp(argv[idx], prefix, strlen(prefix)) != 0) { continue; }

    /* control over MaxIters function */
    if (strcmp(argv[idx] + offset, "max_iters") == 0)
    {
      idx += 1;
      int iarg = atoi(argv[idx]);
      retval   = SUNNonlinSolSetMaxIters(NLS, iarg);
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

/* -----------------------------------------------------------------------------
 * set functions
 * ---------------------------------------------------------------------------*/

/* set the nonlinear system function (required) */
SUNErrCode SUNNonlinSolSetSysFn(SUNNonlinearSolver NLS, SUNNonlinSolSysFn SysFn)
{
  return (NLS->ops->setsysfn(NLS, SysFn));
}

/* set the linear solver setup function (optional) */
SUNErrCode SUNNonlinSolSetLSetupFn(SUNNonlinearSolver NLS,
                                   SUNNonlinSolLSetupFn LSetupFn)
{
  if (NLS->ops->setlsetupfn) { return (NLS->ops->setlsetupfn(NLS, LSetupFn)); }
  else { return (SUN_SUCCESS); }
}

/* set the linear solver solve function (optional) */
SUNErrCode SUNNonlinSolSetLSolveFn(SUNNonlinearSolver NLS,
                                   SUNNonlinSolLSolveFn LSolveFn)
{
  if (NLS->ops->setlsolvefn) { return (NLS->ops->setlsolvefn(NLS, LSolveFn)); }
  else { return (SUN_SUCCESS); }
}

/* set the convergence test function (optional) */
SUNErrCode SUNNonlinSolSetConvTestFn(SUNNonlinearSolver NLS,
                                     SUNNonlinSolConvTestFn CTestFn,
                                     void* ctest_data)
{
  if (NLS->ops->setctestfn)
  {
    return (NLS->ops->setctestfn(NLS, CTestFn, ctest_data));
  }
  else { return (SUN_SUCCESS); }
}

SUNErrCode SUNNonlinSolSetOptions(SUNNonlinearSolver NLS, const char* NLSid,
                                  const char* file_name, int argc, char* argv[])
{
  if (NLS == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(NLS->sunctx);

  /* File-based option control is currently unimplemented */
  SUNAssert((file_name == NULL || strlen(file_name) == 0),
            SUN_ERR_ARG_INCOMPATIBLE);

  /* First, process all base-class options */
  if (argc > 0 && argv != NULL)
  {
    SUNCheckCall(sunnlsSetFromCommandLine(NLS, NLSid, argc, argv));
  }

  /* Second, ask the implementation to process any remaining options */
  if (NLS->ops->setoptions)
  {
    return (NLS->ops->setoptions(NLS, NLSid, file_name, argc, argv));
  }
  else { return (SUN_SUCCESS); }
}

SUNErrCode SUNNonlinSolSetMaxIters(SUNNonlinearSolver NLS, int maxiters)
{
  if (NLS->ops->setmaxiters) { return (NLS->ops->setmaxiters(NLS, maxiters)); }
  else { return (SUN_SUCCESS); }
}

/* -----------------------------------------------------------------------------
 * get functions
 * ---------------------------------------------------------------------------*/

/* get the total number on nonlinear iterations (optional) */
SUNErrCode SUNNonlinSolGetNumIters(SUNNonlinearSolver NLS, long int* niters)
{
  if (NLS->ops->getnumiters) { return (NLS->ops->getnumiters(NLS, niters)); }
  else
  {
    *niters = 0;
    return (SUN_SUCCESS);
  }
}

/* get the iteration count for the current nonlinear solve */
SUNErrCode SUNNonlinSolGetCurIter(SUNNonlinearSolver NLS, int* iter)
{
  if (NLS->ops->getcuriter) { return (NLS->ops->getcuriter(NLS, iter)); }
  else
  {
    *iter = -1;
    return (SUN_SUCCESS);
  }
}

/* get the total number on nonlinear solve convergence failures (optional) */
SUNErrCode SUNNonlinSolGetNumConvFails(SUNNonlinearSolver NLS,
                                       long int* nconvfails)
{
  if (NLS->ops->getnumconvfails)
  {
    return (NLS->ops->getnumconvfails(NLS, nconvfails));
  }
  else
  {
    *nconvfails = 0;
    return (SUN_SUCCESS);
  }
}
