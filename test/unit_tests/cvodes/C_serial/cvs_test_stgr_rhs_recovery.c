/* -----------------------------------------------------------------------------
 * Programmer(s): Bill Hlavacek
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
 * Unit test for a recoverable RHS failure at the converged state in the
 * staggered sensitivity corrector.
 *
 * With CV_STAGGERED or CV_STAGGERED1, cvStep evaluates f at the converged
 * state before solving for the sensitivities. A recoverable failure of that
 * call must be handled like any other recoverable failure: the prediction is
 * restored, the step is retried with a smaller step size, and the failure is
 * counted as a step solve failure (CVodeGetNumStepSolveFails). Retrying from
 * the advanced time instead moves tn by two steps, and the first of them is
 * never corrected.
 *
 * The test solves y' = -p y with one sensitivity in two passes, for each of
 * the two staggered methods. The first pass records which f calls are made at
 * a converged state; that is the last f call before the sensitivity RHS is
 * called. The second pass fails one of those calls once with a recoverable
 * error. Every step must then start where the previous one ended, the solution
 * and sensitivity must still match the exact solution, and the second pass
 * must count exactly one more step solve failure than the first.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "cvodes/cvodes.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sundials/sundials_matrix.h"
#include "sundials/sundials_nvector.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#else
#define GSYM "g"
#endif

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

#define MAX_RECORDED 1000

/* Which converged-state f call (in the order of the first pass) to fail */
#define FAIL_NTH 10

typedef struct
{
  sunrealtype p;    /* decay rate, the sensitivity parameter          */
  long int nf;      /* f calls so far                                 */
  long int last_f;  /* index of the most recent f call                */
  long int fail_at; /* index of the f call to fail, or -1             */
  int failed;       /* whether that call was failed                   */
  int record;       /* whether to record converged-state f calls      */
  long int recorded[MAX_RECORDED];
  int nrecorded;
}* UserData;

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  UserData udata = (UserData)user_data;

  udata->nf++;
  udata->last_f = udata->nf;
  if (udata->nf == udata->fail_at)
  {
    udata->failed = 1;
    return 1;
  }

  NV_Ith_S(ydot, 0) = -udata->p * NV_Ith_S(y, 0);
  return 0;
}

static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData udata = (UserData)user_data;

  SM_ELEMENT_D(J, 0, 0) = -udata->p;
  return 0;
}

static int fS(int Ns, sunrealtype t, N_Vector y, N_Vector ydot, int iS,
              N_Vector yS, N_Vector ySdot, void* user_data, N_Vector tmp1,
              N_Vector tmp2)
{
  UserData udata = (UserData)user_data;

  /* The most recent f call was made at the converged state, unless it was
     already recorded (more than one sensitivity RHS call per f call) */
  if (udata->record && udata->nrecorded < MAX_RECORDED &&
      (udata->nrecorded == 0 ||
       udata->recorded[udata->nrecorded - 1] != udata->last_f))
  {
    udata->recorded[udata->nrecorded++] = udata->last_f;
  }

  NV_Ith_S(ySdot, 0) = -udata->p * NV_Ith_S(yS, 0) - NV_Ith_S(y, 0);
  return 0;
}

/* Check the return flag of a SUNDIALS routine, and report the routine and the
   flag if it failed. Returns 0 on success and 1 on failure. */
static int check_flag(int flag, const char* routine)
{
  if (flag == CV_SUCCESS) { return 0; }
  fprintf(stderr, "ERROR: %s returned %i\n", routine, flag);
  return 1;
}

/* Check a pointer returned by a SUNDIALS routine, and report the routine if it
   returned NULL. Returns 0 on success and 1 on failure. */
static int check_ptr(const void* ptr, const char* routine)
{
  if (ptr) { return 0; }
  fprintf(stderr, "ERROR: %s returned NULL\n", routine);
  return 1;
}

/* Objects created for one integration, freed together by free_resources */
typedef struct
{
  N_Vector y;
  N_Vector* yS;
  SUNMatrix A;
  SUNLinearSolver LS;
  void* cvode_mem;
} Resources;

static void free_resources(Resources* r)
{
  if (r->cvode_mem) { CVodeFree(&r->cvode_mem); }
  if (r->LS) { SUNLinSolFree(r->LS); }
  if (r->A) { SUNMatDestroy(r->A); }
  if (r->yS) { N_VDestroyVectorArray(r->yS, 1); }
  if (r->y) { N_VDestroy(r->y); }
}

/* Integrate to tout one step at a time with the sensitivity method ism, and
   return the number of step solve failures in ncfn. Returns 0 if every step
   started where the previous one ended and the final state and sensitivity are
   accurate. Objects are created in r; the caller frees them. */
static int run_integration(UserData udata, SUNContext sunctx, int ism,
                           Resources* r, long int* ncfn)
{
  const sunrealtype rtol = SUN_RCONST(1.0e-6);
  const sunrealtype atol = SUN_RCONST(1.0e-10);
  const sunrealtype tout = SUN_RCONST(2.0);

  int flag            = 0;
  long int nstep      = 0;
  sunrealtype tret    = ZERO;
  sunrealtype tprev   = ZERO;
  sunrealtype tcur    = ZERO;
  sunrealtype hlast   = ZERO;
  sunrealtype tstart  = ZERO;
  sunrealtype fuzz    = ZERO;
  sunrealtype y_exact = ZERO;
  sunrealtype s_exact = ZERO;
  sunrealtype y_err   = ZERO;
  sunrealtype s_err   = ZERO;

  r->y = N_VNew_Serial(1, sunctx);
  if (check_ptr(r->y, "N_VNew_Serial")) { return 1; }
  NV_Ith_S(r->y, 0) = ONE;

  r->yS = N_VCloneVectorArray(1, r->y);
  if (check_ptr(r->yS, "N_VCloneVectorArray")) { return 1; }
  N_VConst(ZERO, r->yS[0]);

  r->cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_ptr(r->cvode_mem, "CVodeCreate")) { return 1; }

  flag = CVodeInit(r->cvode_mem, f, ZERO, r->y);
  if (check_flag(flag, "CVodeInit")) { return 1; }

  flag = CVodeSetUserData(r->cvode_mem, udata);
  if (check_flag(flag, "CVodeSetUserData")) { return 1; }

  flag = CVodeSStolerances(r->cvode_mem, rtol, atol);
  if (check_flag(flag, "CVodeSStolerances")) { return 1; }

  r->A = SUNDenseMatrix(1, 1, sunctx);
  if (check_ptr(r->A, "SUNDenseMatrix")) { return 1; }

  r->LS = SUNLinSol_Dense(r->y, r->A, sunctx);
  if (check_ptr(r->LS, "SUNLinSol_Dense")) { return 1; }

  flag = CVodeSetLinearSolver(r->cvode_mem, r->LS, r->A);
  if (check_flag(flag, "CVodeSetLinearSolver")) { return 1; }

  flag = CVodeSetJacFn(r->cvode_mem, jac);
  if (check_flag(flag, "CVodeSetJacFn")) { return 1; }

  flag = CVodeSensInit1(r->cvode_mem, 1, ism, fS, r->yS);
  if (check_flag(flag, "CVodeSensInit1")) { return 1; }

  flag = CVodeSensEEtolerances(r->cvode_mem);
  if (check_flag(flag, "CVodeSensEEtolerances")) { return 1; }

  flag = CVodeSetSensErrCon(r->cvode_mem, SUNTRUE);
  if (check_flag(flag, "CVodeSetSensErrCon")) { return 1; }

  while (tret < tout)
  {
    flag = CVode(r->cvode_mem, tout, r->y, &tret, CV_ONE_STEP);
    if (check_flag(flag, "CVode")) { return 1; }
    nstep++;

    flag = CVodeGetCurrentTime(r->cvode_mem, &tcur);
    if (check_flag(flag, "CVodeGetCurrentTime")) { return 1; }

    flag = CVodeGetLastStep(r->cvode_mem, &hlast);
    if (check_flag(flag, "CVodeGetLastStep")) { return 1; }

    /* The step just taken must start where the previous one ended */
    tstart = tcur - hlast;
    fuzz   = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF *
           (SUNRabs(tcur) + SUNRabs(hlast));
    if (SUNRabs(tstart - tprev) > fuzz)
    {
      fprintf(stderr,
              "ERROR: step %li ends at t = %" GSYM " with h = %" GSYM
              ", so it starts at t = %" GSYM
              ", but the previous step ended at t = %" GSYM "\n",
              nstep, tcur, hlast, tstart, tprev);
      return 1;
    }
    tprev = tcur;
  }

  flag = CVodeGetSens(r->cvode_mem, &tret, r->yS);
  if (check_flag(flag, "CVodeGetSens")) { return 1; }

  flag = CVodeGetNumStepSolveFails(r->cvode_mem, ncfn);
  if (check_flag(flag, "CVodeGetNumStepSolveFails")) { return 1; }

  y_exact = SUNRexp(-udata->p * tret);
  s_exact = -tret * y_exact;
  y_err   = SUNRabs(NV_Ith_S(r->y, 0) - y_exact) / SUNRabs(y_exact);
  s_err   = SUNRabs(NV_Ith_S(r->yS[0], 0) - s_exact) / SUNRabs(s_exact);
  if (y_err > SUN_RCONST(1.0e-4) || s_err > SUN_RCONST(1.0e-4))
  {
    fprintf(stderr,
            "ERROR: at t = %" GSYM ", relative error %" GSYM " in y and %" GSYM
            " in dy/dp\n",
            tret, y_err, s_err);
    return 1;
  }

  return 0;
}

/* Run one integration and free everything it created. Returns 0 on success. */
static int integrate(UserData udata, SUNContext sunctx, int ism, long int* ncfn)
{
  Resources r = {NULL, NULL, NULL, NULL, NULL};
  int status  = run_integration(udata, sunctx, ism, &r, ncfn);
  free_resources(&r);
  return status;
}

/* Run both passes for the sensitivity method ism. Returns 0 on success. */
static int test_method(UserData udata, SUNContext sunctx, int ism)
{
  long int ncfn1 = 0;
  long int ncfn2 = 0;

  /* First pass: record the converged-state f calls */
  udata->p         = ONE;
  udata->nf        = 0;
  udata->last_f    = 0;
  udata->fail_at   = -1;
  udata->failed    = 0;
  udata->record    = 1;
  udata->nrecorded = 0;
  if (integrate(udata, sunctx, ism, &ncfn1))
  {
    fprintf(stderr, "ERROR: first pass, with no failure, did not pass\n");
    return 1;
  }
  if (udata->nrecorded <= FAIL_NTH)
  {
    fprintf(stderr, "ERROR: only %i converged-state f calls recorded\n",
            udata->nrecorded);
    return 1;
  }

  /* Second pass: the same run, failing one converged-state f call once */
  udata->fail_at   = udata->recorded[FAIL_NTH];
  udata->record    = 0;
  udata->nf        = 0;
  udata->last_f    = 0;
  udata->failed    = 0;
  udata->nrecorded = 0;
  if (integrate(udata, sunctx, ism, &ncfn2))
  {
    fprintf(stderr, "ERROR: second pass, failing f call %li, did not pass\n",
            udata->fail_at);
    return 1;
  }
  if (!udata->failed)
  {
    fprintf(stderr, "ERROR: f call %li was never made, so nothing was tested\n",
            udata->fail_at);
    return 1;
  }
  if (ncfn2 != ncfn1 + 1)
  {
    fprintf(stderr,
            "ERROR: %li step solve failures with the failed call, %li "
            "without; expected one more\n",
            ncfn2, ncfn1);
    return 1;
  }

  return 0;
}

int main(int argc, char* argv[])
{
  const int methods[2]        = {CV_STAGGERED, CV_STAGGERED1};
  const char* method_names[2] = {"CV_STAGGERED", "CV_STAGGERED1"};

  SUNContext sunctx = NULL;
  UserData udata    = NULL;
  int flag          = 0;
  int failed        = 0;
  int nfailed       = 0;
  int m             = 0;

  flag = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_flag(flag, "SUNContext_Create")) { return 1; }

  udata = (UserData)calloc(1, sizeof(*udata));
  if (check_ptr(udata, "calloc"))
  {
    SUNContext_Free(&sunctx);
    return 1;
  }

  /* Test each method, and report every result */
  for (m = 0; m < 2; m++)
  {
    failed = test_method(udata, sunctx, methods[m]);
    printf("%s: %s\n", method_names[m], failed ? "FAIL" : "PASS");
    nfailed += failed;
  }

  free(udata);
  SUNContext_Free(&sunctx);
  return nfailed ? 1 : 0;
}
