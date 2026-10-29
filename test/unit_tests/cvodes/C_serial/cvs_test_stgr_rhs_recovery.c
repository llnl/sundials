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
 * With CV_STAGGERED, cvStep evaluates f at the converged state before solving
 * for the sensitivities. A recoverable failure of that call must be handled
 * like any other recoverable failure: the prediction is restored and the step
 * is retried with a smaller step size. Retrying from the advanced time instead
 * moves tn by two steps, and the first of them is never corrected.
 *
 * The test solves y' = -p y with one sensitivity in two passes. The first pass
 * records which f calls are made at a converged state; that is the last f call
 * before the sensitivity RHS is called. The second pass fails one of those
 * calls once with a recoverable error. Every step must then start where the
 * previous one ended, and the solution and sensitivity must still match the
 * exact solution.
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

/* Integrate to tout one step at a time. Returns 0 if every step started where
   the previous one ended and the final state and sensitivity are accurate. */
static int integrate(UserData udata, SUNContext sunctx)
{
  const sunrealtype rtol = SUN_RCONST(1.0e-6);
  const sunrealtype atol = SUN_RCONST(1.0e-10);
  const sunrealtype tout = SUN_RCONST(2.0);

  N_Vector y         = NULL;
  N_Vector* yS       = NULL;
  SUNMatrix A        = NULL;
  SUNLinearSolver LS = NULL;
  void* cvode_mem    = NULL;

  int flag            = 0;
  int status          = 1;
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

  y = N_VNew_Serial(1, sunctx);
  if (!y) { goto cleanup; }
  NV_Ith_S(y, 0) = ONE;

  yS = N_VCloneVectorArray(1, y);
  if (!yS) { goto cleanup; }
  N_VConst(ZERO, yS[0]);

  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (!cvode_mem) { goto cleanup; }

  if (CVodeInit(cvode_mem, f, ZERO, y)) { goto cleanup; }
  if (CVodeSetUserData(cvode_mem, udata)) { goto cleanup; }
  if (CVodeSStolerances(cvode_mem, rtol, atol)) { goto cleanup; }

  A = SUNDenseMatrix(1, 1, sunctx);
  if (!A) { goto cleanup; }
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (!LS) { goto cleanup; }
  if (CVodeSetLinearSolver(cvode_mem, LS, A)) { goto cleanup; }
  if (CVodeSetJacFn(cvode_mem, jac)) { goto cleanup; }

  if (CVodeSensInit1(cvode_mem, 1, CV_STAGGERED, fS, yS)) { goto cleanup; }
  if (CVodeSensEEtolerances(cvode_mem)) { goto cleanup; }
  if (CVodeSetSensErrCon(cvode_mem, SUNTRUE)) { goto cleanup; }

  tprev = ZERO;
  while (tret < tout)
  {
    flag = CVode(cvode_mem, tout, y, &tret, CV_ONE_STEP);
    if (flag < 0)
    {
      printf("ERROR: CVode returned %i\n", flag);
      goto cleanup;
    }
    nstep++;

    if (CVodeGetCurrentTime(cvode_mem, &tcur)) { goto cleanup; }
    if (CVodeGetLastStep(cvode_mem, &hlast)) { goto cleanup; }

    /* The step just taken must start where the previous one ended */
    tstart = tcur - hlast;
    fuzz   = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF *
           (SUNRabs(tcur) + SUNRabs(hlast));
    if (SUNRabs(tstart - tprev) > fuzz)
    {
      printf("ERROR: step %li ends at t = %" GSYM " with h = %" GSYM
             ", so it starts at t = %" GSYM
             ", but the previous step ended at t = %" GSYM "\n",
             nstep, tcur, hlast, tstart, tprev);
      goto cleanup;
    }
    tprev = tcur;
  }

  if (CVodeGetSens(cvode_mem, &tret, yS)) { goto cleanup; }

  y_exact = SUNRexp(-udata->p * tret);
  s_exact = -tret * y_exact;
  y_err   = SUNRabs(NV_Ith_S(y, 0) - y_exact) / SUNRabs(y_exact);
  s_err   = SUNRabs(NV_Ith_S(yS[0], 0) - s_exact) / SUNRabs(s_exact);
  if (y_err > SUN_RCONST(1.0e-4) || s_err > SUN_RCONST(1.0e-4))
  {
    printf("ERROR: at t = %" GSYM ", relative error %" GSYM " in y and %" GSYM
           " in dy/dp\n",
           tret, y_err, s_err);
    goto cleanup;
  }

  status = 0;

cleanup:
  N_VDestroy(y);
  if (yS) { N_VDestroyVectorArray(yS, 1); }
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  CVodeFree(&cvode_mem);
  return status;
}

int main(int argc, char* argv[])
{
  SUNContext sunctx = NULL;
  UserData udata    = NULL;
  int status        = 1;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx))
  {
    fprintf(stderr, "SUNContext_Create failed\n");
    return 1;
  }

  udata = (UserData)calloc(1, sizeof(*udata));
  if (!udata) { goto cleanup; }
  udata->p = ONE;

  /* First pass: record the converged-state f calls */
  udata->fail_at = -1;
  udata->record  = 1;
  if (integrate(udata, sunctx))
  {
    printf("ERROR: the first pass, with no failure, did not pass\n");
    goto cleanup;
  }
  if (udata->nrecorded <= FAIL_NTH)
  {
    printf("ERROR: only %i converged-state f calls recorded\n", udata->nrecorded);
    goto cleanup;
  }

  /* Second pass: the same run, failing one converged-state f call once */
  udata->fail_at   = udata->recorded[FAIL_NTH];
  udata->record    = 0;
  udata->nf        = 0;
  udata->last_f    = 0;
  udata->failed    = 0;
  udata->nrecorded = 0;
  if (integrate(udata, sunctx)) { goto cleanup; }
  if (!udata->failed)
  {
    printf("ERROR: f call %li was never made, so nothing was tested\n",
           udata->fail_at);
    goto cleanup;
  }

  printf("SUCCESS\n");
  status = 0;

cleanup:
  free(udata);
  SUNContext_Free(&sunctx);
  return status;
}
