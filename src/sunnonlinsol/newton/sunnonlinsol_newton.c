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
 * This is the implementation file for the SUNNonlinearSolver module
 * implementation of Newton's method.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_nvector_senswrapper.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#include "sundials_logger_impl.h"
#include "sundials_macros.h"

/* Content structure accessibility macros  */
#define NEWTON_CONTENT(S) ((SUNNonlinearSolverContent_Newton)(S->content))

/* Constant macros */
#define ZERO SUN_RCONST(0.0) /* real 0.0 */
#define ONE  SUN_RCONST(1.0) /* real 1.0 */

/*==============================================================================
  Constructor to create a new Newton solver
  ============================================================================*/

SUNNonlinearSolver SUNNonlinSol_Newton(N_Vector y, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNNonlinearSolver NLS;
  SUNNonlinearSolverContent_Newton content;

  /* Check that the supplied N_Vector supports all required operations */
  SUNAssertNull(y->ops->nvclone && y->ops->nvdestroy && y->ops->nvscale &&
                  y->ops->nvlinearsum,
                SUN_ERR_ARG_INCOMPATIBLE);

  /* Create an empty nonlinear linear solver object */
  NLS = SUNNonlinSolNewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  NLS->ops->gettype         = SUNNonlinSolGetType_Newton;
  NLS->ops->initialize      = SUNNonlinSolInitialize_Newton;
  NLS->ops->solve           = SUNNonlinSolSolve_Newton;
  NLS->ops->free            = SUNNonlinSolFree_Newton;
  NLS->ops->setsysfn        = SUNNonlinSolSetSysFn_Newton;
  NLS->ops->setlsetupfn     = SUNNonlinSolSetLSetupFn_Newton;
  NLS->ops->setlsolvefn     = SUNNonlinSolSetLSolveFn_Newton;
  NLS->ops->setctestfn      = SUNNonlinSolSetConvTestFn_Newton;
  NLS->ops->setmaxiters     = SUNNonlinSolSetMaxIters_Newton;
  NLS->ops->getnumiters     = SUNNonlinSolGetNumIters_Newton;
  NLS->ops->getcuriter      = SUNNonlinSolGetCurIter_Newton;
  NLS->ops->getnumconvfails = SUNNonlinSolGetNumConvFails_Newton;

  /* Create content */
  content = NULL;
  content = (SUNNonlinearSolverContent_Newton)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Initialize all components of content to 0/NULL */
  memset(content, 0, sizeof(struct _SUNNonlinearSolverContent_Newton));

  /* Attach content */
  NLS->content = content;

  /* Fill general content */
  content->Sys        = NULL;
  content->LSetup     = NULL;
  content->LSolve     = NULL;
  content->CTest      = NULL;
  content->jcur       = SUNFALSE;
  content->curiter    = 0;
  content->maxiters   = 3;
  content->niters     = 0;
  content->nconvfails = 0;
  content->ctest_data = NULL;

  /* Fill allocatable content */
  content->delta = N_VClone(y);
  SUNCheckLastErrNull();

  return (NLS);
}

/*==============================================================================
  Constructor wrapper to create a new Newton solver for sensitivity solvers
  ============================================================================*/

SUNNonlinearSolver SUNNonlinSol_NewtonSens(int count, N_Vector y,
                                           SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNNonlinearSolver NLS;
  N_Vector w;

  /* create sensitivity vector wrapper */
  w = N_VNew_SensWrapper(count, y);
  SUNCheckLastErrNull();

  /* create nonlinear solver using sensitivity vector wrapper */
  NLS = SUNNonlinSol_Newton(w, sunctx);
  SUNCheckLastErrNull();

  /* free sensitivity vector wrapper */
  N_VDestroy(w);
  SUNCheckLastErrNull();

  /* return NLS object */
  return (NLS);
}

/*==============================================================================
  GetType, Initialize, Setup, Solve, and Free operations
  ============================================================================*/

SUNNonlinearSolver_Type SUNNonlinSolGetType_Newton(
  SUNDIALS_MAYBE_UNUSED SUNNonlinearSolver NLS)
{
  return (SUNNONLINEARSOLVER_ROOTFIND);
}

SUNErrCode SUNNonlinSolInitialize_Newton(SUNNonlinearSolver NLS)
{
  SUNFunctionBegin(NLS->sunctx);
  /* check that all required function pointers have been set */
  SUNAssert(NEWTON_CONTENT(NLS)->Sys && NEWTON_CONTENT(NLS)->CTest &&
              NEWTON_CONTENT(NLS)->LSolve,
            SUN_ERR_ARG_CORRUPT);

  /* reset the total number of iterations and convergence failures */
  NEWTON_CONTENT(NLS)->niters     = 0;
  NEWTON_CONTENT(NLS)->nconvfails = 0;

  /* reset the Jacobian status */
  NEWTON_CONTENT(NLS)->jcur = SUNFALSE;

  return SUN_SUCCESS;
}

/*------------------------------------------------------------------------------
  SUNNonlinSolSolve_Newton: Performs the nonlinear solve F(y) = 0

  Successful solve return code:
    SUN_SUCCESS = 0

  Recoverable failure return codes (positive):
    SUN_NLS_CONV_RECVR
    *_RHSFUNC_RECVR (ODEs) or *_RES_RECVR (DAEs)
    *_LSETUP_RECVR
    *_LSOLVE_RECVR

  Unrecoverable failure return codes (negative):
    SUN_ERR_*
    *_RHSFUNC_FAIL (ODEs) or *_RES_FAIL (DAEs)
    *_LSETUP_FAIL
    *_LSOLVE_FAIL

  Note return values beginning with * are package specific values returned by
  the Sys, LSetup, and LSolve functions provided to the nonlinear solver.
  ----------------------------------------------------------------------------*/
int SUNNonlinSolSolve_Newton(SUNNonlinearSolver NLS,
                             SUNDIALS_MAYBE_UNUSED N_Vector y0, N_Vector ycor,
                             N_Vector w, sunrealtype tol,
                             sunbooleantype callLSetup, void* mem)
{
  SUNFunctionBegin(NLS->sunctx);
  /* local variables */
  int retval;
  sunbooleantype jbad;
  N_Vector delta;

  /* check that all required function pointers have been set */
  SUNAssert(NEWTON_CONTENT(NLS)->Sys && NEWTON_CONTENT(NLS)->CTest &&
              NEWTON_CONTENT(NLS)->LSolve,
            SUN_ERR_ARG_CORRUPT);
  SUNAssert(!callLSetup || (callLSetup && NEWTON_CONTENT(NLS)->LSetup),
            SUN_ERR_ARG_CORRUPT);

  SUNLogInfo(NLS->sunctx->logger, "nonlinear-solver", "solver = Newton");

  /* set local shortcut variables */
  delta = NEWTON_CONTENT(NLS)->delta;

  /* assume the Jacobian is good */
  jbad = SUNFALSE;

  /* initialize iteration and convergence fail counters for this solve */
  NEWTON_CONTENT(NLS)->niters     = 0;
  NEWTON_CONTENT(NLS)->nconvfails = 0;

  /* looping point for attempts at solution of the nonlinear system:
       Evaluate the nonlinear residual function (store in delta)
       Setup the linear solver if necessary
       Perform Newton iteration */
  for (;;)
  {
    SUNLogInfo(NLS->sunctx->logger, "begin-iterations-list", "");

    /* initialize current iteration counter for this solve attempt */
    NEWTON_CONTENT(NLS)->curiter = 0;

    /* compute the nonlinear residual, store in delta */
    retval = NEWTON_CONTENT(NLS)->Sys(ycor, delta, mem);
    if (retval != SUN_SUCCESS) { break; }

    /* if indicated, setup the linear system */
    if (callLSetup)
    {
      retval = NEWTON_CONTENT(NLS)->LSetup(jbad, &(NEWTON_CONTENT(NLS)->jcur),
                                           mem);
      if (retval != SUN_SUCCESS) { break; }
    }

    /* looping point for Newton iteration. Break out on any error. */
    for (;;)
    {
      /* increment nonlinear solver iteration counter */
      NEWTON_CONTENT(NLS)->niters++;

      /* compute the negative of the residual for the linear system rhs */
      N_VScale(-ONE, delta, delta);
      SUNCheckLastErr();

      /* solve the linear system to get Newton update delta */
      retval = NEWTON_CONTENT(NLS)->LSolve(delta, mem);
      if (retval != SUN_SUCCESS) { break; }

      /* update the Newton iterate */
      N_VLinearSum(ONE, ycor, ONE, delta, ycor);
      SUNCheckLastErr();

      /* test for convergence */
      retval = NEWTON_CONTENT(NLS)->CTest(NLS, ycor, delta, tol, w,
                                          NEWTON_CONTENT(NLS)->ctest_data);

      NEWTON_CONTENT(NLS)->curiter++;

      SUNLogInfo(NLS->sunctx->logger, "nonlinear-iterate",
                 "cur-iter = %i, total-iters = %li, update-norm = %.16g",
                 NEWTON_CONTENT(NLS)->curiter, NEWTON_CONTENT(NLS)->niters,
                 N_VWrmsNorm(delta, w));

      /* if successful update Jacobian status and return */
      if (retval == SUN_SUCCESS)
      {
        SUNLogInfo(NLS->sunctx->logger, "end-iterations-list",
                   "status = success");
        NEWTON_CONTENT(NLS)->jcur = SUNFALSE;
        return SUN_SUCCESS;
      }

      /* check if the iteration should continue; otherwise exit Newton loop */
      if (retval != SUN_NLS_CONTINUE) { break; }

      /* not yet converged, test for max allowed iterations. */
      if (NEWTON_CONTENT(NLS)->curiter >= NEWTON_CONTENT(NLS)->maxiters)
      {
        retval = SUN_NLS_CONV_RECVR;
        break;
      }

      SUNLogInfo(NLS->sunctx->logger, "end-iterations-list", "status = continue");
      SUNLogInfo(NLS->sunctx->logger, "begin-iterations-list", "");

      /* compute the nonlinear residual, store in delta */
      retval = NEWTON_CONTENT(NLS)->Sys(ycor, delta, mem);
      if (retval != SUN_SUCCESS) { break; }

    } /* end of Newton iteration loop */

    /* all errors go here */

    /* If there is a recoverable convergence failure and the Jacobian-related
       data appears not to be current, increment the convergence failure count,
       reset the initial correction to zero, and loop again with a call to
       lsetup in which jbad is TRUE. Otherwise break out and return. */
    if ((retval > 0) && !(NEWTON_CONTENT(NLS)->jcur) &&
        (NEWTON_CONTENT(NLS)->LSetup))
    {
      NEWTON_CONTENT(NLS)->nconvfails++;
      callLSetup = SUNTRUE;
      jbad       = SUNTRUE;
      N_VConst(ZERO, ycor);
      SUNCheckLastErr();

      SUNLogInfo(NLS->sunctx->logger, "end-iterations-list", "status = continue");

      continue;
    }
    else
    {
      SUNLogInfo(NLS->sunctx->logger, "end-iterations-list",
                 "status = failure, retval = %i", retval);
      break;
    }

  } /* end of setup loop */

  /* increment number of convergence failures */
  NEWTON_CONTENT(NLS)->nconvfails++;

  /* all error returns exit here */
  return (retval);
}

SUNErrCode SUNNonlinSolFree_Newton(SUNNonlinearSolver NLS)
{
  /* return if NLS is already free */
  if (NLS == NULL) { return SUN_SUCCESS; }

  SUNFunctionBegin(NLS->sunctx);
  /* free items from contents, then the generic structure */
  if (NLS->content)
  {
    if (NEWTON_CONTENT(NLS)->delta)
    {
      N_VDestroy(NEWTON_CONTENT(NLS)->delta);
      SUNCheckLastErr();
    }
    NEWTON_CONTENT(NLS)->delta = NULL;

    free(NLS->content);
    NLS->content = NULL;
  }

  /* free the ops structure */
  if (NLS->ops)
  {
    free(NLS->ops);
    NLS->ops = NULL;
  }

  /* free the nonlinear solver */
  free(NLS);

  return SUN_SUCCESS;
}

/*==============================================================================
  Set functions
  ============================================================================*/

SUNErrCode SUNNonlinSolSetSysFn_Newton(SUNNonlinearSolver NLS,
                                       SUNNonlinSolSysFn SysFn)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNAssert(SysFn, SUN_ERR_ARG_CORRUPT);
  NEWTON_CONTENT(NLS)->Sys = SysFn;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetLSetupFn_Newton(SUNNonlinearSolver NLS,
                                          SUNNonlinSolLSetupFn LSetupFn)
{
  NEWTON_CONTENT(NLS)->LSetup = LSetupFn;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetLSolveFn_Newton(SUNNonlinearSolver NLS,
                                          SUNNonlinSolLSolveFn LSolveFn)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNAssert(LSolveFn, SUN_ERR_ARG_CORRUPT);
  NEWTON_CONTENT(NLS)->LSolve = LSolveFn;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetConvTestFn_Newton(SUNNonlinearSolver NLS,
                                            SUNNonlinSolConvTestFn CTestFn,
                                            void* ctest_data)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNAssert(CTestFn, SUN_ERR_ARG_CORRUPT);

  NEWTON_CONTENT(NLS)->CTest = CTestFn;

  /* attach convergence test data */
  NEWTON_CONTENT(NLS)->ctest_data = ctest_data;

  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetMaxIters_Newton(SUNNonlinearSolver NLS, int maxiters)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNAssert(maxiters >= 1, SUN_ERR_ARG_OUTOFRANGE);
  NEWTON_CONTENT(NLS)->maxiters = maxiters;
  return SUN_SUCCESS;
}

/*==============================================================================
  Get functions
  ============================================================================*/

SUNErrCode SUNNonlinSolGetNumIters_Newton(SUNNonlinearSolver NLS, long int* niters)
{
  /* return the number of nonlinear iterations in the last solve */
  *niters = NEWTON_CONTENT(NLS)->niters;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolGetCurIter_Newton(SUNNonlinearSolver NLS, int* iter)
{
  /* return the current nonlinear solver iteration count */
  *iter = NEWTON_CONTENT(NLS)->curiter;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolGetNumConvFails_Newton(SUNNonlinearSolver NLS,
                                              long int* nconvfails)
{
  /* return the total number of nonlinear convergence failures */
  *nconvfails = NEWTON_CONTENT(NLS)->nconvfails;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolGetSysFn_Newton(SUNNonlinearSolver NLS,
                                       SUNNonlinSolSysFn* SysFn)
{
  /* return the nonlinear system defining function */
  *SysFn = NEWTON_CONTENT(NLS)->Sys;
  return SUN_SUCCESS;
}
