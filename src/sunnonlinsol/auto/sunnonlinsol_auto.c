#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#include <sunnonlinsol/sunnonlinsol_auto.h>

#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#include "sundials/priv/sundials_errors_impl.h"
#include "sundials/sundials_errors.h"
#include "sundials_logger_impl.h"
#include "sundials/sundials_nonlinearsolver.h"

/* Content structure accessibility macros */
#define AUTO_CONTENT(S) ((SUNNonlinearSolverContent_Auto)(S->content))

typedef struct
{
  SUNNonlinearSolver auto_nls;
  SUNNonlinSolConvTestFn user_ctest_fn;
  void* user_ctest_data;
} SUNNonlinSolAutoConvTestData;

static int SUNNonlinSolConvTest_Auto(SUNNonlinearSolver sub_nls, N_Vector y,
                                     N_Vector del, sunrealtype tol,
                                     N_Vector ewt, void* mem);

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
                                     SUNNonlinSolAutoType active_solver_type,
                                     SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  SUNNonlinearSolver NLS                 = NULL;
  SUNNonlinearSolverContent_Auto content = NULL;

  NLS = SUNNonlinSolNewEmpty(sunctx);
  SUNCheckLastErrNull();

  NLS->ops->gettype         = SUNNonlinSolGetType_Auto;
  NLS->ops->initialize      = SUNNonlinSolInitialize_Auto;
  NLS->ops->solve           = SUNNonlinSolSolve_Auto;
  NLS->ops->free            = SUNNonlinSolFree_Auto;
  NLS->ops->setsysfn        = SUNNonlinSolSetSysFn_Auto;
  NLS->ops->setctestfn      = SUNNonlinSolSetConvTestFn_Auto;
  NLS->ops->setlsetupfn     = SUNNonlinSolSetLSetupFn_Auto;
  NLS->ops->setlsolvefn     = SUNNonlinSolSetLSolveFn_Auto;
  NLS->ops->setmaxiters     = SUNNonlinSolSetMaxIters_Auto;
  NLS->ops->getnumiters     = SUNNonlinSolGetNumIters_Auto;
  NLS->ops->getcuriter      = SUNNonlinSolGetCurIter_Auto;
  NLS->ops->getnumconvfails = SUNNonlinSolGetNumConvFails_Auto;
  NLS->ops->getdelnrm       = SUNNonlinSolGetDelNrm_Auto;

  content = (SUNNonlinearSolverContent_Auto)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  NLS->content = content;

  content->active_solver_type   = active_solver_type;
  content->user_ctest_fn        = NULL;
  content->user_ctest_data      = NULL;
  content->maxiters             = 3;
  content->curiter              = 0;
  content->niters               = 0;
  content->nconvfails           = 0;
  content->fp_to_newt_delay     = 0;
  content->newt_to_fp_delay     = 10;
  content->fp_to_newt_threshold = SUN_RCONST(0.8);
  content->newt_to_fp_threshold = SUN_RCONST(2.0);
  content->nsolves_since_switch = 0;
  content->switch_count         = 0;
  content->fp_niters_total      = 0;
  content->newt_niters_total    = 0;
  content->auto_ctest_data      = NULL;
  content->fp_solver            = SUNNonlinSol_FixedPoint(y, m, sunctx);
  content->newton_solver        = SUNNonlinSol_Newton(y, sunctx);

  return NLS;
}

SUNNonlinearSolver_Type SUNNonlinSolGetType_Auto(SUNNonlinearSolver NLS)
{
  return SUNNONLINEARSOLVER_HYBRID;
}

SUNErrCode SUNNonlinSolInitialize_Auto(SUNNonlinearSolver NLS)
{
  SUNFunctionBegin(NLS->sunctx);
  if (AUTO_CONTENT(NLS)->active_solver_type == SUNNONLINSOL_AUTO_FIXEDPOINT)
  {
    SUNCheckCall(SUNNonlinSolInitialize(AUTO_CONTENT(NLS)->fp_solver));
  }
  else
  {
    SUNCheckCall(SUNNonlinSolInitialize(AUTO_CONTENT(NLS)->newton_solver));
  }
  return SUN_SUCCESS;
}

int SUNNonlinSolSolve_Auto(SUNNonlinearSolver NLS, N_Vector y0, N_Vector ycor,
                           N_Vector w, sunrealtype tol,
                           sunbooleantype callSetup, void* mem)
{
  SUNFunctionBegin(NLS->sunctx);
  int retval;
  long int iters;

  SUNLogInfo(NLS->sunctx->logger, "nonlinear-solver",
             "solver = Auto, active = %s",
             SUNNonlinSolAutoType_ToString(AUTO_CONTENT(NLS)->active_solver_type));

  if (AUTO_CONTENT(NLS)->active_solver_type == SUNNONLINSOL_AUTO_FIXEDPOINT)
  {
    retval = SUNNonlinSolSolve(AUTO_CONTENT(NLS)->fp_solver, y0, ycor, w, tol,
                               callSetup, mem);
    iters = 0;
    if (SUNNonlinSolGetNumIters(AUTO_CONTENT(NLS)->fp_solver, &iters) ==
        SUN_SUCCESS)
    {
      AUTO_CONTENT(NLS)->fp_niters_total += iters;
    }
  }
  else
  {
    retval = SUNNonlinSolSolve(AUTO_CONTENT(NLS)->newton_solver, y0, ycor, w,
                               tol, callSetup, mem);
    iters = 0;
    if (SUNNonlinSolGetNumIters(AUTO_CONTENT(NLS)->newton_solver, &iters) ==
        SUN_SUCCESS)
    {
      AUTO_CONTENT(NLS)->newt_niters_total += iters;
    }
  }

  /* increment solve counter used for switch-delay gating */
  AUTO_CONTENT(NLS)->nsolves_since_switch++;

  return retval;
}

static int SUNNonlinSolConvTest_Auto(SUNNonlinearSolver sub_nls, N_Vector y,
                                     N_Vector del, sunrealtype tol,
                                     N_Vector ewt, void* mem)
{
  SUNNonlinSolAutoConvTestData* data = (SUNNonlinSolAutoConvTestData*)mem;
  SUNNonlinearSolver auto_nls        = data->auto_nls;
  SUNNonlinearSolverContent_Auto C   = AUTO_CONTENT(auto_nls);

  int retval = data->user_ctest_fn(sub_nls, y, del, tol, ewt,
                                   data->user_ctest_data);

  /* We follow the switching strategy outlined in 
     Nørsett, Syvert P., and Per G. Thomsen. "Switching between modified Newton 
     and fix-point iteration for implicit ODE-solvers." BIT Numerical Mathematics
     26, no. 3 (1986): 339-348. https://doi.org/10.1007/BF01933714. */

  if (C->active_solver_type == SUNNONLINSOL_AUTO_FIXEDPOINT)
  {
    /* In the unlikely event that we are 'diverging', but the integrator-provided
       convergence test passed, we just exit with success and don't consider any switching. */
    if (retval == SUN_SUCCESS) { return retval; }

    SUNNonlinearSolverContent_FixedPoint fp_content =
      (SUNNonlinearSolverContent_FixedPoint)C->fp_solver->content;

    /* Check if we are diverging */
    sunbooleantype diverging  = fp_content->crate >= C->fp_to_newt_threshold;
    sunbooleantype dont_delay = C->nsolves_since_switch >= C->fp_to_newt_delay;

    if (diverging && dont_delay)
    {
      C->nsolves_since_switch = 0;
      C->active_solver_type   = SUNNONLINSOL_AUTO_NEWTON;
      C->switch_count++;
      SUNLogInfo(auto_nls->sunctx->logger, "auto-nonlinear-solver-switch",
                 "from = Fixed-Point, to = Newton, crate = " SUN_FORMAT_G
                 ", threshold = " SUN_FORMAT_G ", delay = %li",
                 fp_content->crate, C->fp_to_newt_threshold, C->fp_to_newt_delay);
      return SUN_NLS_SWITCH;
    }
  }
  else
  {
    SUNNonlinearSolverContent_Newton newton_content =
      (SUNNonlinearSolverContent_Newton)C->newton_solver->content;

    sunbooleantype contraction = newton_content->stiffr < C->newt_to_fp_threshold;
    sunbooleantype dont_delay = C->nsolves_since_switch >= C->newt_to_fp_delay;

    if (contraction && dont_delay)
    {
      C->nsolves_since_switch = 0;
      C->active_solver_type   = SUNNONLINSOL_AUTO_FIXEDPOINT;
      C->switch_count++;
      SUNLogInfo(auto_nls->sunctx->logger, "auto-nonlinear-solver-switch",
                 "from = Newton, to = Fixed-Point, stiffr = " SUN_FORMAT_G
                 ", threshold = " SUN_FORMAT_G ", delay = %li",
                 newton_content->stiffr, C->newt_to_fp_threshold,
                 C->newt_to_fp_delay);
      return SUN_NLS_SWITCH;
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
    if (AUTO_CONTENT(NLS)->fp_solver)
    {
      AUTO_CONTENT(NLS)->fp_solver->ops->free(AUTO_CONTENT(NLS)->fp_solver);
    }
    if (AUTO_CONTENT(NLS)->newton_solver)
    {
      AUTO_CONTENT(NLS)->newton_solver->ops->free(
        AUTO_CONTENT(NLS)->newton_solver);
    }
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

SUNErrCode SUNNonlinSolSetSysFn_Auto(SUNNonlinearSolver NLS,
                                     SUNNonlinSolSysFn SysFn)
{
  SUNFunctionBegin(NLS->sunctx);
  if (AUTO_CONTENT(NLS)->active_solver_type == SUNNONLINSOL_AUTO_FIXEDPOINT)
  {
    SUNCheckCall(SUNNonlinSolSetSysFn(AUTO_CONTENT(NLS)->fp_solver, SysFn));
  }
  else
  {
    SUNCheckCall(SUNNonlinSolSetSysFn(AUTO_CONTENT(NLS)->newton_solver, SysFn));
  }
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

  AUTO_CONTENT(NLS)->user_ctest_fn   = CTestFn;
  AUTO_CONTENT(NLS)->user_ctest_data = ctest_data;

  if (!AUTO_CONTENT(NLS)->auto_ctest_data)
  {
    AUTO_CONTENT(NLS)->auto_ctest_data =
      malloc(sizeof(SUNNonlinSolAutoConvTestData));
    SUNAssert(AUTO_CONTENT(NLS)->auto_ctest_data, SUN_ERR_MALLOC_FAIL);
  }

  SUNNonlinSolAutoConvTestData* data =
    (SUNNonlinSolAutoConvTestData*)AUTO_CONTENT(NLS)->auto_ctest_data;
  data->auto_nls        = NLS;
  data->user_ctest_fn   = CTestFn;
  data->user_ctest_data = ctest_data;

  SUNCheckCall(SUNNonlinSolSetConvTestFn(AUTO_CONTENT(NLS)->fp_solver,
                                         SUNNonlinSolConvTest_Auto, data));
  SUNCheckCall(SUNNonlinSolSetConvTestFn(AUTO_CONTENT(NLS)->newton_solver,
                                         SUNNonlinSolConvTest_Auto, data));
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetLSetupFn_Auto(SUNNonlinearSolver NLS,
                                        SUNNonlinSolLSetupFn LSetupFn)
{
  SUNFunctionBegin(NLS->sunctx);
  if (AUTO_CONTENT(NLS)->active_solver_type == SUNNONLINSOL_AUTO_NEWTON)
  {
    SUNCheckCall(
      SUNNonlinSolSetLSetupFn(AUTO_CONTENT(NLS)->newton_solver, LSetupFn));
  }
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetLSolveFn_Auto(SUNNonlinearSolver NLS,
                                        SUNNonlinSolLSolveFn LSolveFn)
{
  SUNFunctionBegin(NLS->sunctx);
  if (AUTO_CONTENT(NLS)->active_solver_type == SUNNONLINSOL_AUTO_NEWTON)
  {
    SUNCheckCall(
      SUNNonlinSolSetLSolveFn(AUTO_CONTENT(NLS)->newton_solver, LSolveFn));
  }
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolSetMaxIters_Auto(SUNNonlinearSolver NLS, int maxiters)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNErrCode retval = SUN_SUCCESS;
  if (AUTO_CONTENT(NLS)->active_solver_type == SUNNONLINSOL_AUTO_FIXEDPOINT)
  {
    retval = SUNNonlinSolSetMaxIters(AUTO_CONTENT(NLS)->fp_solver, maxiters);
  }
  else
  {
    retval = SUNNonlinSolSetMaxIters(AUTO_CONTENT(NLS)->newton_solver, maxiters);
  }
  return retval;
}

SUNErrCode SUNNonlinSolSetSwitchingParameters_Auto(
  SUNNonlinearSolver NLS, sunrealtype newt_to_fp_threshold,
  long int newt_to_fp_delay, sunrealtype fp_to_newt_threshold,
  long int fp_to_newt_delay)
{
  SUNFunctionBegin(NLS->sunctx);

  SUNAssert(newt_to_fp_threshold <= SUN_RCONST(2.0), SUN_ERR_ARG_OUTOFRANGE);
  SUNAssert(fp_to_newt_threshold <= SUN_RCONST(1.0), SUN_ERR_ARG_OUTOFRANGE);

  AUTO_CONTENT(NLS)->newt_to_fp_threshold =
    (newt_to_fp_threshold < SUN_RCONST(0.0)) ? SUN_RCONST(2.0)
                                             : newt_to_fp_threshold;
  AUTO_CONTENT(NLS)->newt_to_fp_delay = (newt_to_fp_delay < 0) ? 10
                                                               : newt_to_fp_delay;
  AUTO_CONTENT(NLS)->fp_to_newt_threshold =
    (fp_to_newt_threshold < SUN_RCONST(0.0)) ? SUN_RCONST(0.8)
                                             : fp_to_newt_threshold;
  AUTO_CONTENT(NLS)->fp_to_newt_delay = (fp_to_newt_delay < 0) ? 0
                                                               : fp_to_newt_delay;

  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolGetNumIters_Auto(SUNNonlinearSolver NLS, long int* niters)
{
  SUNFunctionBegin(NLS->sunctx);
  long int fp_iters   = 0;
  long int newt_iters = 0;
  SUNCheckCall(SUNNonlinSolGetNumIters(AUTO_CONTENT(NLS)->fp_solver, &fp_iters));
  SUNCheckCall(
    SUNNonlinSolGetNumIters(AUTO_CONTENT(NLS)->newton_solver, &newt_iters));
  *niters = fp_iters + newt_iters;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolGetNumItersByType_Auto(SUNNonlinearSolver NLS,
                                              long int* fp_iters,
                                              long int* newt_iters)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNAssert(fp_iters, SUN_ERR_ARG_CORRUPT);
  SUNAssert(newt_iters, SUN_ERR_ARG_CORRUPT);

  *fp_iters   = AUTO_CONTENT(NLS)->fp_niters_total;
  *newt_iters = AUTO_CONTENT(NLS)->newt_niters_total;

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
  long int fp_nvconvfails  = 0;
  long int newt_nconvfails = 0;
  SUNCheckCall(
    SUNNonlinSolGetNumConvFails(AUTO_CONTENT(NLS)->fp_solver, &fp_nvconvfails));
  SUNCheckCall(SUNNonlinSolGetNumConvFails(AUTO_CONTENT(NLS)->newton_solver,
                                           &newt_nconvfails));
  *nconvfails = fp_nvconvfails + newt_nconvfails;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolGetDelNrm_Auto(SUNNonlinearSolver NLS, sunrealtype* delnrm)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNAssert(delnrm, SUN_ERR_ARG_CORRUPT);
  if (AUTO_CONTENT(NLS)->active_solver_type == SUNNONLINSOL_AUTO_FIXEDPOINT)
  {
    return SUNNonlinSolGetDelNrm(AUTO_CONTENT(NLS)->fp_solver, delnrm);
  }
  else
  {
    return SUNNonlinSolGetDelNrm(AUTO_CONTENT(NLS)->newton_solver, delnrm);
  }
}
