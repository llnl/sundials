#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#include <sunnonlinsol/sunnonlinsol_auto.h>

#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#include "sundials_logger_impl.h"
#include "sundials_macros.h"

/* Content structure accessibility macros */
#define AUTO_CONTENT(S) ((SUNNonlinearSolverContent_Auto)(S->content))

SUNNonlinearSolver SUNNonlinSol_Auto(N_Vector y, int m,
                                     SUNNonlinSolAutoType type, SUNContext sunctx)
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
  NLS->ops->setmaxiters     = SUNNonlinSolSetMaxIters_Auto;
  NLS->ops->getnumiters     = SUNNonlinSolGetNumIters_Auto;
  NLS->ops->getcuriter      = SUNNonlinSolGetCurIter_Auto;
  NLS->ops->getnumconvfails = SUNNonlinSolGetNumConvFails_Auto;

  content = (SUNNonlinearSolverContent_Auto)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  NLS->content = content;

  content->type          = type;
  content->maxiters      = 3;
  content->curiter       = 0;
  content->niters        = 0;
  content->nconvfails    = 0;
  content->fp_solver     = SUNNonlinSol_FixedPoint(y, m, sunctx);
  content->newton_solver = SUNNonlinSol_Newton(y, sunctx);

  return NLS;
}

SUNNonlinearSolver_Type SUNNonlinSolGetType_Auto(SUNNonlinearSolver NLS)
{
  if (AUTO_CONTENT(NLS)->type == SUNNONLINSOL_AUTO_FIXEDPOINT)
  {
    return SUNNONLINEARSOLVER_FIXEDPOINT;
  }
  else { return SUNNONLINEARSOLVER_ROOTFIND; }
}

SUNErrCode SUNNonlinSolInitialize_Auto(SUNNonlinearSolver NLS)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNErrCode retval = SUN_SUCCESS;
  if (AUTO_CONTENT(NLS)->type == SUNNONLINSOL_AUTO_FIXEDPOINT)
  {
    retval = SUNNonlinSolInitialize(AUTO_CONTENT(NLS)->fp_solver);
  }
  else { retval = SUNNonlinSolInitialize(AUTO_CONTENT(NLS)->newton_solver); }
  return retval;
}

int SUNNonlinSolSolve_Auto(SUNNonlinearSolver NLS, N_Vector y0, N_Vector ycor,
                           N_Vector w, sunrealtype tol,
                           sunbooleantype callSetup, void* mem)
{
  SUNFunctionBegin(NLS->sunctx);
  int retval;
  if (AUTO_CONTENT(NLS)->type == SUNNONLINSOL_AUTO_FIXEDPOINT)
  {
    retval = SUNNonlinSolSolve(AUTO_CONTENT(NLS)->fp_solver, y0, ycor, w, tol,
                               callSetup, mem);

    SUNNonlinearSolverContent_FixedPoint fp_content =
      (SUNNonlinearSolverContent_FixedPoint)AUTO_CONTENT(NLS)->fp_solver->content;
    const sunrealtype alpha = 0.8;
    if (fp_content->crate >= alpha)
    {
      printf(">>>> crate=%g\n", fp_content->crate);
    }
  }
  else
  {
    retval = SUNNonlinSolSolve(AUTO_CONTENT(NLS)->newton_solver, y0, ycor, w,
                               tol, callSetup, mem);
  }
  return retval;
}

SUNErrCode SUNNonlinSolFree_Auto(SUNNonlinearSolver NLS)
{
  if (NLS == NULL) { return SUN_SUCCESS; }
  if (NLS->content)
  {
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
  SUNErrCode retval = SUN_SUCCESS;
  if (AUTO_CONTENT(NLS)->type == SUNNONLINSOL_AUTO_FIXEDPOINT)
  {
    retval = SUNNonlinSolSetSysFn(AUTO_CONTENT(NLS)->fp_solver, SysFn);
  }
  else
  {
    retval = SUNNonlinSolSetSysFn(AUTO_CONTENT(NLS)->newton_solver, SysFn);
  }
  return retval;
}

SUNErrCode SUNNonlinSolSetConvTestFn_Auto(SUNNonlinearSolver NLS,
                                          SUNNonlinSolConvTestFn CTestFn,
                                          void* ctest_data)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNErrCode retval = SUN_SUCCESS;
  if (AUTO_CONTENT(NLS)->type == SUNNONLINSOL_AUTO_FIXEDPOINT)
  {
    retval = SUNNonlinSolSetConvTestFn(AUTO_CONTENT(NLS)->fp_solver, CTestFn,
                                       ctest_data);
  }
  else
  {
    retval = SUNNonlinSolSetConvTestFn(AUTO_CONTENT(NLS)->newton_solver,
                                       CTestFn, ctest_data);
  }
  return retval;
}

SUNErrCode SUNNonlinSolSetMaxIters_Auto(SUNNonlinearSolver NLS, int maxiters)
{
  SUNFunctionBegin(NLS->sunctx);
  SUNErrCode retval = SUN_SUCCESS;
  if (AUTO_CONTENT(NLS)->type == SUNNONLINSOL_AUTO_FIXEDPOINT)
  {
    retval = SUNNonlinSolSetMaxIters(AUTO_CONTENT(NLS)->fp_solver, maxiters);
  }
  else
  {
    retval = SUNNonlinSolSetMaxIters(AUTO_CONTENT(NLS)->newton_solver, maxiters);
  }
  return retval;
}

SUNErrCode SUNNonlinSolGetNumIters_Auto(SUNNonlinearSolver NLS, long int* niters)
{
  long int fp_iters   = 0;
  long int newt_iters = 0;
  SUNCheckCall(SUNNonlinSolGetNumIters(AUTO_CONTENT(NLS)->fp_solver, &fp_iters));
  SUNCheckCall(
    SUNNonlinSolGetNumIters(AUTO_CONTENT(NLS)->newton_solver, &newt_iters));
  *niters = fp_iters + newt_iters;
  return SUN_SUCCESS;
}

SUNErrCode SUNNonlinSolGetCurIter_Auto(SUNNonlinearSolver NLS, int* iter)
{
  if (AUTO_CONTENT(NLS)->type == SUNNONLINSOL_AUTO_FIXEDPOINT)
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
  long int fp_nvconvfails  = 0;
  long int newt_nconvfails = 0;
  SUNCheckCall(
    SUNNonlinSolGetNumIters(AUTO_CONTENT(NLS)->fp_solver, &fp_nvconvfails));
  SUNCheckCall(SUNNonlinSolGetNumIters(AUTO_CONTENT(NLS)->newton_solver,
                                       &newt_nconvfails));
  *nconvfails = fp_nvconvfails + newt_nconvfails;
  return SUN_SUCCESS;
}

/* Optionally, add a setter to switch type at runtime */
SUNErrCode SUNNonlinSolSetType_Auto(SUNNonlinearSolver NLS,
                                    SUNNonlinSolAutoType type)
{
  AUTO_CONTENT(NLS)->type = type;
  return SUN_SUCCESS;
}
