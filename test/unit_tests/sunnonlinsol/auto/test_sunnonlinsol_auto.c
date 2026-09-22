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
 * Focused unit test for the SUNNonlinearSolver Auto module.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "nvector/nvector_serial.h"
#include "sundials/sundials_types.h"
#include "sunnonlinsol/sunnonlinsol_auto.h"
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"
#include "sunnonlinsol/sunnonlinsol_newton.h"

#define ZERO SUN_RCONST(0.0)
#define HALF SUN_RCONST(0.5)
#define ONE  SUN_RCONST(1.0)

typedef struct
{
  int convtest_calls;
  int convrate_calls;
  int continue_calls;
}* TestMem;

static int check_retval(void* flagvalue, const char* funcname, int opt);
static int SetupProblem(SUNNonlinearSolver NLS, TestMem mem, int continue_calls);
static int FPFunction(N_Vector ycor, N_Vector gy, void* mem);
static int Res(N_Vector ycor, N_Vector f, void* mem);
static int LSolve(N_Vector b, void* mem);
static int ConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del,
                    sunrealtype tol, N_Vector ewt, void* mem);
static SUNErrCode GetConvRate(sunrealtype* crate, void* mem);

int main(void)
{
  int retval               = 0;
  SUNContext sunctx        = NULL;
  SUNNonlinearSolver NLS   = NULL;
  N_Vector y0              = NULL;
  N_Vector ycor            = NULL;
  N_Vector w               = NULL;
  TestMem mem              = NULL;
  long int switch_count    = 0;
  long int fp_iters        = 0;
  long int newt_iters      = 0;
  long int nconvfails      = 0;
  long int fp_nconvfails   = 0;
  long int newt_nconvfails = 0;
  SUNNonlinSolAutoType active_solver_type;
  SUNNonlinearSolverContent_Auto content;
  SUNNonlinearSolverContent_FixedPoint fp_content;
  SUNNonlinearSolverContent_Newton newton_content;

  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return 1; }

  mem = (TestMem)malloc(sizeof(*mem));
  if (check_retval(mem, "malloc", 0)) { return 1; }
  mem->convtest_calls = 0;
  mem->convrate_calls = 0;
  mem->continue_calls = 0;

  y0 = N_VNew_Serial(1, sunctx);
  if (check_retval(y0, "N_VNew_Serial", 0)) { return 1; }

  ycor = N_VClone(y0);
  if (check_retval(ycor, "N_VClone", 0)) { return 1; }

  w = N_VClone(y0);
  if (check_retval(w, "N_VClone", 0)) { return 1; }

  NV_Ith_S(y0, 0)   = ZERO;
  NV_Ith_S(ycor, 0) = ZERO;
  NV_Ith_S(w, 0)    = ONE;

  NLS = SUNNonlinSol_Auto(y0, 0, SUNNONLINSOL_AUTO_FIXEDPOINT, sunctx);
  if (check_retval(NLS, "SUNNonlinSol_Auto", 0)) { return 1; }

  content = (SUNNonlinearSolverContent_Auto)NLS->content;
  fp_content = (SUNNonlinearSolverContent_FixedPoint)content->fp_solver->content;
  newton_content =
    (SUNNonlinearSolverContent_Newton)content->newton_solver->content;

  if (fp_content->CTest == NULL || fp_content->ctest_data == NULL)
  {
    printf("ERROR: fixed-point subsolver wrapper callback was not installed\n");
    return 1;
  }

  if (newton_content->CTest == NULL || newton_content->ctest_data == NULL)
  {
    printf("ERROR: Newton subsolver wrapper callback was not installed\n");
    return 1;
  }

  retval = SetupProblem(NLS, mem, 1);
  if (check_retval(&retval, "SetupProblem", 1)) { return 1; }

  retval = SUNNonlinSolSetSwitchingParameters_Auto(NLS, HALF, 1, HALF, 0);
  if (check_retval(&retval, "SUNNonlinSolSetSwitchingParameters_Auto", 1))
  {
    return 1;
  }

  retval = SUNNonlinSolSetMaxIters(NLS, 1);
  if (check_retval(&retval, "SUNNonlinSolSetMaxIters", 1)) { return 1; }

  if (fp_content->maxiters != 1 || newton_content->maxiters != 1)
  {
    printf("ERROR: Auto max-iter setter did not update both subsolvers\n");
    return 1;
  }

  retval = SUNNonlinSolSolve(NLS, y0, ycor, w, ONE, SUNFALSE, mem);
  if (retval != SUN_NLS_CONV_RECVR)
  {
    printf("ERROR: expected fixed-point convergence failure without a switch, "
           "got %d\n",
           retval);
    return 1;
  }

  retval = SUNNonlinSolGetActiveSolverType_Auto(NLS, &active_solver_type);
  if (check_retval(&retval, "SUNNonlinSolGetActiveSolverType_Auto", 1))
  {
    return 1;
  }

  if (active_solver_type != SUNNONLINSOL_AUTO_FIXEDPOINT)
  {
    printf(
      "ERROR: expected fixed-point to remain active after first iteration\n");
    return 1;
  }

  retval = SUNNonlinSolGetSwitchCount_Auto(NLS, &switch_count);
  if (check_retval(&retval, "SUNNonlinSolGetSwitchCount_Auto", 1)) { return 1; }

  if (switch_count != 0)
  {
    printf("ERROR: expected no automatic solver switch, got %ld\n", switch_count);
    return 1;
  }

  if (mem->convrate_calls != 0)
  {
    printf("ERROR: expected no first-iteration convergence rate query, got "
           "%d\n",
           mem->convrate_calls);
    return 1;
  }

  SUNNonlinSolFree(NLS);

  NLS = SUNNonlinSol_Auto(y0, 0, SUNNONLINSOL_AUTO_FIXEDPOINT, sunctx);
  if (check_retval(NLS, "SUNNonlinSol_Auto", 0)) { return 1; }

  retval = SetupProblem(NLS, mem, 2);
  if (check_retval(&retval, "SetupProblem", 1)) { return 1; }

  retval = SUNNonlinSolSetSwitchingParameters_Auto(NLS, HALF, 1, HALF, 0);
  if (check_retval(&retval, "SUNNonlinSolSetSwitchingParameters_Auto", 1))
  {
    return 1;
  }

  retval = SUNNonlinSolSetMaxIters(NLS, 2);
  if (check_retval(&retval, "SUNNonlinSolSetMaxIters", 1)) { return 1; }

  retval = SUNNonlinSolSolve(NLS, y0, ycor, w, ONE, SUNFALSE, mem);
  if (check_retval(&retval, "SUNNonlinSolSolve", 1)) { return 1; }

  retval = SUNNonlinSolGetActiveSolverType_Auto(NLS, &active_solver_type);
  if (check_retval(&retval, "SUNNonlinSolGetActiveSolverType_Auto", 1))
  {
    return 1;
  }

  if (active_solver_type != SUNNONLINSOL_AUTO_NEWTON)
  {
    printf("ERROR: expected Newton to be active after the switch\n");
    return 1;
  }

  retval = SUNNonlinSolGetSwitchCount_Auto(NLS, &switch_count);
  if (check_retval(&retval, "SUNNonlinSolGetSwitchCount_Auto", 1)) { return 1; }

  if (switch_count != 1)
  {
    printf("ERROR: expected one automatic solver switch, got %ld\n",
           switch_count);
    return 1;
  }

  retval = SUNNonlinSolGetTotalNumItersByType_Auto(NLS, &fp_iters, &newt_iters);
  if (check_retval(&retval, "SUNNonlinSolGetTotalNumItersByType_Auto", 1))
  {
    return 1;
  }

  if (fp_iters != 2 || newt_iters != 1)
  {
    printf("ERROR: expected two fixed-point iterations and one Newton "
           "iteration, got fp=%ld newt=%ld\n",
           fp_iters, newt_iters);
    return 1;
  }

  if (mem->convrate_calls != 1)
  {
    printf("ERROR: expected one convergence rate query, got %d\n",
           mem->convrate_calls);
    return 1;
  }

  retval = SUNNonlinSolGetNumConvFails(NLS, &nconvfails);
  if (check_retval(&retval, "SUNNonlinSolGetNumConvFails", 1)) { return 1; }

  retval = SUNNonlinSolGetTotalNumConvFailsByType_Auto(NLS, &fp_nconvfails,
                                                       &newt_nconvfails);
  if (check_retval(&retval, "SUNNonlinSolGetTotalNumConvFailsByType_Auto", 1))
  {
    return 1;
  }

  if (nconvfails != 0 || fp_nconvfails != 0 || newt_nconvfails != 0)
  {
    printf("ERROR: expected zero convergence failures, got total=%ld fp=%ld "
           "newt=%ld\n",
           nconvfails, fp_nconvfails, newt_nconvfails);
    return 1;
  }

  SUNNonlinSolFree(NLS);
  N_VDestroy(y0);
  N_VDestroy(ycor);
  N_VDestroy(w);
  free(mem);
  SUNContext_Free(&sunctx);

  printf("SUCCESS\n");
  return 0;
}

static int SetupProblem(SUNNonlinearSolver NLS, TestMem mem, int continue_calls)
{
  int retval;

  mem->convtest_calls = 0;
  mem->convrate_calls = 0;
  mem->continue_calls = continue_calls;

  retval = SUNNonlinSolSetSysFns(NLS, Res, FPFunction);
  if (retval != SUN_SUCCESS) { return retval; }

  retval = SUNNonlinSolSetLSolveFn(NLS, LSolve);
  if (retval != SUN_SUCCESS) { return retval; }

  retval = SUNNonlinSolSetConvTestFn(NLS, ConvTest, mem);
  if (retval != SUN_SUCCESS) { return retval; }

  return SUNNonlinSolSetGetConvRateFn(NLS, GetConvRate, mem);
}

static int FPFunction(N_Vector ycor, N_Vector gy, void* mem)
{
  (void)ycor;
  (void)mem;
  NV_Ith_S(gy, 0) = ZERO;
  return 0;
}

static int Res(N_Vector ycor, N_Vector f, void* mem)
{
  (void)mem;
  NV_Ith_S(f, 0) = NV_Ith_S(ycor, 0) - ONE;
  return 0;
}

static int LSolve(N_Vector b, void* mem)
{
  (void)b;
  (void)mem;
  return 0;
}

static int ConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del,
                    sunrealtype tol, N_Vector ewt, void* mem)
{
  TestMem test_mem = (TestMem)mem;

  (void)NLS;
  (void)y;
  (void)del;
  (void)tol;
  (void)ewt;

  test_mem->convtest_calls++;
  return (test_mem->convtest_calls <= test_mem->continue_calls) ? SUN_NLS_CONTINUE
                                                                : SUN_SUCCESS;
}

static SUNErrCode GetConvRate(sunrealtype* crate, void* mem)
{
  TestMem test_mem = (TestMem)mem;
  test_mem->convrate_calls++;
  *crate = ONE;
  return SUN_SUCCESS;
}

static int check_retval(void* flagvalue, const char* funcname, int opt)
{
  if (opt == 0 && flagvalue == NULL)
  {
    printf("ERROR: %s() returned NULL\n", funcname);
    return 1;
  }

  if (opt == 1)
  {
    int err = *((int*)flagvalue);
    if (err != SUN_SUCCESS)
    {
      printf("ERROR: %s() returned %d\n", funcname, err);
      return 1;
    }
  }

  return 0;
}
