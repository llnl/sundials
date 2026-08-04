/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 * This is the implementation file for the SUNNonlinearSolver module
 * implementation that automatically switches between Newton and fixed-point
 * iterations.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>
#include <sunnonlinsol/sunnonlinsol_auto.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#include "sundials_logger_impl.h"
#include "sundials_macros.h"

/* Content structure accessibility macros */
#define AUTO_CONTENT(S) ((SUNNonlinearSolverContent_Auto)(S->content))

/* this is effectively the maximum number of times we would allow 
   switching to happen within a single solve. It is unlikely to 
   ever be reached, but it prevents an infinite switchign loop from
   being possible. */
#define SUNNLS_AUTO_MAX_SOLVE_ATTEMPTS 3

/* Default switching parameters 
   2.0 and 0.8 come from the numerical experiments in Norsett & Thomsen 1986,
   as does the Newton to fixed-point switch delay of 10 solves. 
   The idea behind setting the fixed-point to Newton switch delay to 1 is to
   allow the integrator to see at least one failed convergence test, cut the time step,
   and try fixed-point again before switching to Newton. 
*/
#define SUNNLS_AUTO_DEFAULT_NEWT_TO_FP_THRESHOLD SUN_RCONST(2.0)
#define SUNNLS_AUTO_DEFAULT_FP_TO_NEWT_THRESHOLD SUN_RCONST(0.8)
#define SUNNLS_AUTO_DEFAULT_NEWT_TO_FP_DELAY     10
#define SUNNLS_AUTO_DEFAULT_FP_TO_NEWT_DELAY     1

static SUNErrCode setFromCommandLine_Auto(SUNNonlinearSolver NLS,
                                          const char* NLSid, int argc,
                                          char* argv[]);

typedef struct
{
  SUNNonlinearSolver auto_nls;
  SUNNonlinSolConvTestFn user_ctest_fn;
  void* user_ctest_data;
} SUNNonlinSolAutoConvTestData;

static int SUNNonlinSolConvTest_Auto(SUNNonlinearSolver sub_nls, N_Vector y,
                                     N_Vector del, sunrealtype tol,
                                     N_Vector ewt, void* mem);

SUNDIALS_MAYBE_UNUSED
static const char* SUNNonlinSolAutoType_ToString(SUNNonlinSolAutoType type)
{
  switch (type)
  {
  case SUNNONLINSOL_AUTO_FIXEDPOINT: return "Fixed-Point";
  case SUNNONLINSOL_AUTO_NEWTON: return "Newton";
  default: return "Unknown";
  }
}

SUNNonlinearSolver SUNNonlinSol_Auto(N_Vector y, int m,
                                     SUNNonlinSolAutoType initial_solver_type,
                                     SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNNonlinearSolver NLS                 = NULL;
  SUNNonlinearSolverContent_Auto content = NULL;

  NLS = SUNNonlinSolNewEmpty(sunctx);
  SUNCheckLastErrNull();

  NLS->ops->gettype            = SUNNonlinSolGetType_Auto;
  NLS->ops->initialize         = SUNNonlinSolInitialize_Auto;
  NLS->ops->solve              = SUNNonlinSolSolve_Auto;
  NLS->ops->free               = SUNNonlinSolFree_Auto;
  NLS->ops->setsysfns          = SUNNonlinSolSetSysFns_Auto;
  NLS->ops->setctestfn         = SUNNonlinSolSetConvTestFn_Auto;
  NLS->ops->setnormfn          = SUNNonlinSolSetNormFn_Auto;
  NLS->ops->setgetupdatenormfn = SUNNonlinSolSetGetUpdateNormFn_Auto;
  NLS->ops->setgetconvratefn   = SUNNonlinSolSetGetConvRateFn_Auto;
  NLS->ops->setlsetupfn        = SUNNonlinSolSetLSetupFn_Auto;
  NLS->ops->setlsolvefn        = SUNNonlinSolSetLSolveFn_Auto;
  NLS->ops->setoptions         = SUNNonlinSolSetOptions_Auto;
  NLS->ops->setmaxiters        = SUNNonlinSolSetMaxIters_Auto;
  NLS->ops->getnumiters        = SUNNonlinSolGetNumIters_Auto;
  NLS->ops->getcuriter         = SUNNonlinSolGetCurIter_Auto;
  NLS->ops->getnumconvfails    = SUNNonlinSolGetNumConvFails_Auto;

  content = (SUNNonlinearSolverContent_Auto)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  NLS->content = content;

  content->active_solver_type      = initial_solver_type;
  content->getconvrate_fn          = NULL;
  content->getconvrate_data        = NULL;
  content->num_iters               = 0;
  content->num_conv_fails          = 0;
  content->fp_to_newt_delay        = SUNNLS_AUTO_DEFAULT_FP_TO_NEWT_DELAY;
  content->newt_to_fp_delay        = SUNNLS_AUTO_DEFAULT_NEWT_TO_FP_DELAY;
  content->fp_to_newt_threshold    = SUNNLS_AUTO_DEFAULT_FP_TO_NEWT_THRESHOLD;
  content->newt_to_fp_threshold    = SUNNLS_AUTO_DEFAULT_NEWT_TO_FP_THRESHOLD;
  content->num_solves_since_switch = 0;
  content->switch_count            = 0;
  content->fp_num_iters_total      = 0;
  content->newton_num_iters_total  = 0;
  content->fp_num_conv_fails_total = 0;
  content->newton_num_conv_fails_total = 0;
  content->auto_ctest_data             = NULL;
  content->fp_solver                   = SUNNonlinSol_FixedPoint(y, m, sunctx);
  content->newton_solver               = SUNNonlinSol_Newton(y, sunctx);
  if (content->fp_solver == NULL || content->newton_solver == NULL)
  {
    if (content->fp_solver) { SUNNonlinSolFree(content->fp_solver); }
    if (content->newton_solver) { SUNNonlinSolFree(content->newton_solver); }
    free(content);
    NLS->content = NULL;
    SUNNonlinSolFreeEmpty(NLS);
    return NULL;
  }

  content->auto_ctest_data = malloc(sizeof(SUNNonlinSolAutoConvTestData));
  if (content->auto_ctest_data == NULL)
  {
    SUNNonlinSolFree(content->fp_solver);
    SUNNonlinSolFree(content->newton_solver);
    free(content);
    NLS->content = NULL;
    SUNNonlinSolFreeEmpty(NLS);
    return NULL;
  }

  ((SUNNonlinSolAutoConvTestData*)content->auto_ctest_data)->auto_nls = NLS;
  ((SUNNonlinSolAutoConvTestData*)content->auto_ctest_data)->user_ctest_fn = NULL;
  ((SUNNonlinSolAutoConvTestData*)content->auto_ctest_data)->user_ctest_data =
    NULL;

  if (SUNNonlinSolSetConvTestFn(content->fp_solver, SUNNonlinSolConvTest_Auto,
                                content->auto_ctest_data) != SUN_SUCCESS ||
      SUNNonlinSolSetConvTestFn(content->newton_solver, SUNNonlinSolConvTest_Auto,
                                content->auto_ctest_data) != SUN_SUCCESS)
  {
    free(content->auto_ctest_data);
    SUNNonlinSolFree(content->fp_solver);
    SUNNonlinSolFree(content->newton_solver);
    free(content);
    NLS->content = NULL;
    SUNNonlinSolFreeEmpty(NLS);
    return NULL;
  }

  if (SUNNonlinSolSetComputeStiffnessRatio_Newton(content->newton_solver,
                                                  SUNTRUE) != SUN_SUCCESS)
  {
    free(content->auto_ctest_data);
    SUNNonlinSolFree(content->fp_solver);
    SUNNonlinSolFree(content->newton_solver);
    free(content);
    NLS->content = NULL;
    SUNNonlinSolFreeEmpty(NLS);
    return NULL;
  }

  return NLS;
}

SUNNonlinearSolver_Type SUNNonlinSolGetType_Auto(
  SUNDIALS_MAYBE_UNUSED SUNNonlinearSolver NLS)
{
  return SUNNONLINEARSOLVER_HYBRID;
}

SUNErrCode SUNNonlinSolInitialize_Auto(SUNNonlinearSolver NLS)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNCheckCall(SUNNonlinSolInitialize(AUTO_CONTENT(NLS)->fp_solver));
  SUNCheckCall(SUNNonlinSolInitialize(AUTO_CONTENT(NLS)->newton_solver));
  return SUN_SUCCESS;
}

int SUNNonlinSolSolve_Auto(SUNNonlinearSolver NLS, N_Vector y0, N_Vector ycor,
                           N_Vector w, sunrealtype tol,
                           sunbooleantype callSetup, void* mem)
{
  SUNFunctionBegin(NLS->sunctx);
  int retval;
  long int attempts, iters, nconvfails;
  SUNNonlinearSolver subsolver;
  SUNNonlinearSolverContent_Auto C = AUTO_CONTENT(NLS);
  SUNNonlinSolAutoType solve_solver_type;

  C->num_iters      = 0;
  C->num_conv_fails = 0;

  for (attempts = 0; attempts < SUNNLS_AUTO_MAX_SOLVE_ATTEMPTS; attempts++)
  {
    solve_solver_type = C->active_solver_type;
    subsolver         = (solve_solver_type == SUNNONLINSOL_AUTO_FIXEDPOINT)
                          ? C->fp_solver
                          : C->newton_solver;

    SUNLogInfo(NLS->sunctx->logger, "begin-subsolver-solves-list",
               "requested = %s",
               SUNNonlinSolAutoType_ToString(C->active_solver_type));

    retval = SUNNonlinSolSolve(subsolver, y0, ycor, w, tol, callSetup, mem);

    iters = 0;
    if (SUNNonlinSolGetNumIters(subsolver, &iters) == SUN_SUCCESS)
    {
      C->num_iters += iters;
      if (solve_solver_type == SUNNONLINSOL_AUTO_FIXEDPOINT)
      {
        C->fp_num_iters_total += iters;
      }
      else { C->newton_num_iters_total += iters; }
    }

    nconvfails = 0;
    if (SUNNonlinSolGetNumConvFails(subsolver, &nconvfails) == SUN_SUCCESS)
    {
      C->num_conv_fails += nconvfails;
      if (solve_solver_type == SUNNONLINSOL_AUTO_FIXEDPOINT)
      {
        C->fp_num_conv_fails_total += nconvfails;
      }
      else { C->newton_num_conv_fails_total += nconvfails; }
    }

    SUNLogInfo(NLS->sunctx->logger, "end-subsolver-solves-list",
               "status = complete, retval = %i, iters = %li, conv-fails = %li",
               retval, iters, nconvfails);

    if (retval == SUN_NLS_SWITCH) { continue; }

    C->num_solves_since_switch++;

    return retval;
  }

  return SUN_SUCCESS;
}

int SUNNonlinSolConvTest_Auto(SUNNonlinearSolver sub_nls, N_Vector y,
                              N_Vector del, sunrealtype tol, N_Vector ewt,
                              void* mem)
{
  SUNNonlinSolAutoConvTestData* data = (SUNNonlinSolAutoConvTestData*)mem;
  if (data == NULL || data->user_ctest_fn == NULL)
  {
    return SUN_ERR_ARG_CORRUPT;
  }
  SUNNonlinearSolver auto_nls      = data->auto_nls;
  SUNNonlinearSolverContent_Auto C = AUTO_CONTENT(auto_nls);

  int retval = data->user_ctest_fn(sub_nls, y, del, tol, ewt,
                                   data->user_ctest_data);
  /* Return early if convergence test error is unrecoverable */
  if (retval < 0) { return retval; }

  /* We follow the switching strategy outlined in
     Norsett, Syvert P., and Per G. Thomsen. "Switching between modified Newton
     and fix-point iteration for implicit ODE-solvers." BIT Numerical Mathematics
     26, no. 3 (1986): 339-348. https://doi.org/10.1007/BF01933714. */

  if (C->active_solver_type == SUNNONLINSOL_AUTO_FIXEDPOINT)
  {
    sunrealtype crate;
    SUNErrCode crate_retval;

    /* If the integrator-provided convergence test passed, exit with success and
       don't consider switching, since fixed-point is still converging fine. */
    if (retval == SUN_SUCCESS) { return SUN_SUCCESS; }

    /* Get the convergence rate from the user-provided function */
    if (C->getconvrate_fn == NULL) { return SUN_ERR_NOT_IMPLEMENTED; }
    crate_retval = C->getconvrate_fn(&crate, C->getconvrate_data);
    if (crate_retval != SUN_SUCCESS) { return crate_retval; }

    sunbooleantype diverging = crate >= C->fp_to_newt_threshold;
    sunbooleantype dont_delay = C->num_solves_since_switch >= C->fp_to_newt_delay;

    if (diverging && dont_delay)
    {
      C->num_solves_since_switch = 0;
      C->active_solver_type      = SUNNONLINSOL_AUTO_NEWTON;
      C->switch_count++;
      SUNLogInfo(auto_nls->sunctx->logger, "auto-nonlinear-solver-switch",
                 "from = Fixed-Point, to = Newton, crate = " SUN_FORMAT_G
                 ", threshold = " SUN_FORMAT_G
                 ", delay = %li, user_ctest_retval = %i",
                 crate, C->fp_to_newt_threshold, C->fp_to_newt_delay, retval);
      /* Return SUN_NLS_SWITCH so that the solver loop continues but with Newton */
      return SUN_NLS_SWITCH;
    }
  }
  else
  {
    /* Since Newton is active, check if we should switch to fixed-point regardless
       of if the convergence test passed. */

    /* Get the stiffness ratio from the Newton solver */
    sunrealtype stiffr = SUN_RCONST(0.0);
    SUNErrCode stiffr_retval =
      SUNNonlinSolGetStiffnessRatio_Newton(C->newton_solver, &stiffr);
    if (stiffr_retval != SUN_SUCCESS) { return stiffr_retval; }

    sunbooleantype contraction = stiffr < C->newt_to_fp_threshold;
    sunbooleantype dont_delay = C->num_solves_since_switch >= C->newt_to_fp_delay;

    if (contraction && dont_delay)
    {
      C->num_solves_since_switch = 0;
      C->active_solver_type      = SUNNONLINSOL_AUTO_FIXEDPOINT;
      C->switch_count++;
      SUNLogInfo(auto_nls->sunctx->logger, "auto-nonlinear-solver-switch",
                 "from = Newton, to = Fixed-Point, stiffr = " SUN_FORMAT_G
                 ", threshold = " SUN_FORMAT_G
                 ", delay = %li, user_ctest_retval = %i",
                 stiffr, C->newt_to_fp_threshold, C->newt_to_fp_delay, retval);

      /* If the integrator-provided convergence test passed, then we return with success
         so that the solve loop is stopped. The switch to fixed-point will happen on
         the next time step. If the convergence test failed, we return SUN_NLS_SWITCH
         so that the solver loop continues but with fixed-point iteration. */
      if (retval == SUN_SUCCESS) { return SUN_SUCCESS; }
      else { return SUN_NLS_SWITCH; }
    }
  }

  return retval;
}

SUNErrCode SUNNonlinSolFree_Auto(SUNNonlinearSolver NLS)
{
  if (NLS == NULL) { return SUN_SUCCESS; }

  if (NLS->content)
  {
    if (AUTO_CONTENT(NLS)->auto_ctest_data)
    {
      free(AUTO_CONTENT(NLS)->auto_ctest_data);
      AUTO_CONTENT(NLS)->auto_ctest_data = NULL;
    }
    SUNNonlinSolFree(AUTO_CONTENT(NLS)->fp_solver);
    SUNNonlinSolFree(AUTO_CONTENT(NLS)->newton_solver);
    free(NLS->content);
    NLS->content = NULL;
  }
  if (NLS->ops)
  {
    free(NLS->ops);
    NLS->ops = NULL;
  }
  free(NLS);
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetSysFns_Auto(SUNNonlinearSolver NLS,
                                      SUNNonlinSolSysFn root_sys_fn,
                                      SUNNonlinSolSysFn fixed_point_fn)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNCheckCall(
    SUNNonlinSolSetSysFn(AUTO_CONTENT(NLS)->newton_solver, root_sys_fn));
  SUNCheckCall(SUNNonlinSolSetSysFn(AUTO_CONTENT(NLS)->fp_solver, fixed_point_fn));
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetConvTestFn_Auto(SUNNonlinearSolver NLS,
                                          SUNNonlinSolConvTestFn CTestFn,
                                          void* ctest_data)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNAssert(CTestFn, SUN_ERR_ARG_CORRUPT);

  SUNNonlinSolAutoConvTestData* data =
    (SUNNonlinSolAutoConvTestData*)AUTO_CONTENT(NLS)->auto_ctest_data;
  data->auto_nls        = NLS;
  data->user_ctest_fn   = CTestFn;
  data->user_ctest_data = ctest_data;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetLSetupFn_Auto(SUNNonlinearSolver NLS,
                                        SUNNonlinSolLSetupFn LSetupFn)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNCheckCall(
    SUNNonlinSolSetLSetupFn(AUTO_CONTENT(NLS)->newton_solver, LSetupFn));
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetLSolveFn_Auto(SUNNonlinearSolver NLS,
                                        SUNNonlinSolLSolveFn LSolveFn)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNCheckCall(
    SUNNonlinSolSetLSolveFn(AUTO_CONTENT(NLS)->newton_solver, LSolveFn));
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetNormFn_Auto(SUNNonlinearSolver NLS,
                                      SUNNonlinSolNormFn NormFn,
                                      void* norm_fn_data)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNCheckCall(
    SUNNonlinSolSetNormFn(AUTO_CONTENT(NLS)->fp_solver, NormFn, norm_fn_data));
  SUNCheckCall(SUNNonlinSolSetNormFn(AUTO_CONTENT(NLS)->newton_solver, NormFn,
                                     norm_fn_data));
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetGetUpdateNormFn_Auto(
  SUNNonlinearSolver NLS, SUNNonlinSolGetUpdateNormFn GetUpdateNormFn,
  void* getupdatenorm_data)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNCheckCall(SUNNonlinSolSetGetUpdateNormFn(AUTO_CONTENT(NLS)->fp_solver,
                                              GetUpdateNormFn,
                                              getupdatenorm_data));
  SUNCheckCall(SUNNonlinSolSetGetUpdateNormFn(AUTO_CONTENT(NLS)->newton_solver,
                                              GetUpdateNormFn,
                                              getupdatenorm_data));
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetGetConvRateFn_Auto(SUNNonlinearSolver NLS,
                                             SUNNonlinSolGetConvRateFn GetConvRateFn,
                                             void* getconvrate_data)
{
  SUNFunctionBegin(NLS->sunctx);
  AUTO_CONTENT(NLS)->getconvrate_fn   = GetConvRateFn;
  AUTO_CONTENT(NLS)->getconvrate_data = getconvrate_data;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetMaxIters_Auto(SUNNonlinearSolver NLS, int maxiters)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNCheckCall(SUNNonlinSolSetMaxIters(AUTO_CONTENT(NLS)->fp_solver, maxiters));
  SUNCheckCall(
    SUNNonlinSolSetMaxIters(AUTO_CONTENT(NLS)->newton_solver, maxiters));
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetSwitchingParameters_Auto(
  SUNNonlinearSolver NLS, sunrealtype newt_to_fp_threshold,
  long int newt_to_fp_delay, sunrealtype fp_to_newt_threshold,
  long int fp_to_newt_delay)
{
  SUNFunctionBegin(NLS->sunctx);

  AUTO_CONTENT(NLS)->newt_to_fp_threshold =
    (newt_to_fp_threshold < SUN_RCONST(0.0))
      ? SUNNLS_AUTO_DEFAULT_NEWT_TO_FP_THRESHOLD
      : newt_to_fp_threshold;
  AUTO_CONTENT(NLS)->newt_to_fp_delay = (newt_to_fp_delay < 0)
                                          ? SUNNLS_AUTO_DEFAULT_NEWT_TO_FP_DELAY
                                          : newt_to_fp_delay;
  AUTO_CONTENT(NLS)->fp_to_newt_threshold =
    (fp_to_newt_threshold < SUN_RCONST(0.0))
      ? SUNNLS_AUTO_DEFAULT_FP_TO_NEWT_THRESHOLD
      : fp_to_newt_threshold;
  AUTO_CONTENT(NLS)->fp_to_newt_delay = (fp_to_newt_delay < 0)
                                          ? SUNNLS_AUTO_DEFAULT_FP_TO_NEWT_DELAY
                                          : fp_to_newt_delay;

  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolGetNumIters_Auto(SUNNonlinearSolver NLS, long int* niters)
{
  SUNFunctionBegin(NLS->sunctx);
  *niters = AUTO_CONTENT(NLS)->num_iters;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolGetFixedPointSolver_Auto(SUNNonlinearSolver NLS,
                                                SUNNonlinearSolver* fp_nls)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNAssert(fp_nls, SUN_ERR_ARG_CORRUPT);
  *fp_nls = AUTO_CONTENT(NLS)->fp_solver;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolGetNewtonSolver_Auto(SUNNonlinearSolver NLS,
                                            SUNNonlinearSolver* newton_nls)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNAssert(newton_nls, SUN_ERR_ARG_CORRUPT);
  *newton_nls = AUTO_CONTENT(NLS)->newton_solver;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolGetActiveSolverType_Auto(
  SUNNonlinearSolver NLS, SUNNonlinSolAutoType* active_solver_type)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNAssert(active_solver_type, SUN_ERR_ARG_CORRUPT);
  *active_solver_type = AUTO_CONTENT(NLS)->active_solver_type;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolGetSwitchCount_Auto(SUNNonlinearSolver NLS,
                                           long int* switch_count)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNAssert(switch_count, SUN_ERR_ARG_CORRUPT);
  *switch_count = AUTO_CONTENT(NLS)->switch_count;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolGetTotalNumItersByType_Auto(SUNNonlinearSolver NLS,
                                                   long int* fp_iters,
                                                   long int* newt_iters)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNAssert(fp_iters, SUN_ERR_ARG_CORRUPT);
  SUNAssert(newt_iters, SUN_ERR_ARG_CORRUPT);

  *fp_iters   = AUTO_CONTENT(NLS)->fp_num_iters_total;
  *newt_iters = AUTO_CONTENT(NLS)->newton_num_iters_total;

  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolGetCurIter_Auto(SUNNonlinearSolver NLS, int* iter)
{
  if (AUTO_CONTENT(NLS)->active_solver_type == SUNNONLINSOL_AUTO_FIXEDPOINT)
  {
    return SUNNonlinSolGetCurIter(AUTO_CONTENT(NLS)->fp_solver, iter);
  }
  else
  {
    return SUNNonlinSolGetCurIter(AUTO_CONTENT(NLS)->newton_solver, iter);
  }
}

SUNErrCode SUNNonlinSolGetNumConvFails_Auto(SUNNonlinearSolver NLS,
                                            long int* nconvfails)
{
  SUNFunctionBegin(NLS->sunctx);
  *nconvfails = AUTO_CONTENT(NLS)->num_conv_fails;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolGetTotalNumConvFailsByType_Auto(SUNNonlinearSolver NLS,
                                                       long int* fp_nconvfails,
                                                       long int* newt_nconvfails)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNAssert(fp_nconvfails, SUN_ERR_ARG_CORRUPT);
  SUNAssert(newt_nconvfails, SUN_ERR_ARG_CORRUPT);
  *fp_nconvfails   = AUTO_CONTENT(NLS)->fp_num_conv_fails_total;
  *newt_nconvfails = AUTO_CONTENT(NLS)->newton_num_conv_fails_total;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetOptions_Auto(SUNNonlinearSolver NLS, const char* NLSid,
                                       SUNDIALS_MAYBE_UNUSED const char* file_name,
                                       int argc, char* argv[])
{
  SUNFunctionBegin(NLS->sunctx);

  SUNAssert((file_name == NULL || strlen(file_name) == 0),
            SUN_ERR_ARG_INCOMPATIBLE);

  if (argc > 0 && argv != NULL)
  {
    SUNCheckCall(setFromCommandLine_Auto(NLS, NLSid, argc, argv));
  }

  return SUN_SUCCESS;
}

static SUNErrCode setFromCommandLine_Auto(SUNNonlinearSolver NLS,
                                          const char* NLSid, int argc,
                                          char* argv[])
{
  SUNFunctionBegin(NLS->sunctx);

  const char* default_id = "sunnonlinearsolver";
  size_t offset          = strlen(default_id) + 1;
  if (NLSid != NULL && strlen(NLSid) > 0) { offset = strlen(NLSid) + 1; }

  char* prefix = (char*)malloc(sizeof(char) * (offset + 1));
  SUNAssert(prefix, SUN_ERR_MALLOC_FAIL);
  if (NLSid != NULL && strlen(NLSid) > 0) { strcpy(prefix, NLSid); }
  else { strcpy(prefix, default_id); }
  strcat(prefix, ".");

  for (int idx = 1; idx < argc; idx++)
  {
    int retval;

    if (strncmp(argv[idx], prefix, strlen(prefix)) != 0) { continue; }

    if (strcmp(argv[idx] + offset, "switching_parameters") == 0)
    {
      if (idx + 4 >= argc)
      {
        free(prefix);
        return SUN_ERR_ARG_INCOMPATIBLE;
      }
      retval =
        SUNNonlinSolSetSwitchingParameters_Auto(NLS,
                                                (sunrealtype)atof(argv[idx + 1]),
                                                atol(argv[idx + 2]),
                                                (sunrealtype)atof(argv[idx + 3]),
                                                atol(argv[idx + 4]));
      if (retval != SUN_SUCCESS)
      {
        free(prefix);
        return retval;
      }
      idx += 4;
      continue;
    }
  }

  free(prefix);
  return SUN_SUCCESS;
}
