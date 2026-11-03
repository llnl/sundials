/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
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
 * Unit test for the STS temporal error estimate in ExtSTS methods.  The test
 * problem is 1D advection-diffusion on (0,1),
 *
 *   u_t = D u_xx + c u_x,   u(0,t) = u(1,t) = 0,   u(x,0) = sin(pi x),
 *
 * discretized with second-order centered differences on N interior points,
 * with diffusion as f^D (STS) and advection as f^E (explicit, slow).  For c = 0
 * the semi-discrete solution is u(t) = exp(lambda t) sin(pi x), with
 * lambda = -(4 D / dx^2) sin^2(pi dx / 2).
 *
 * Tests 1-4 enable the STS error estimate with
 * MRIStepSetUseExtSTSErrorEstimate.  The tests are:
 *   1. c = 0 (no slow dynamics): the adaptive ExtSTS solution must be accurate
 *      and use a number of steps comparable to LSRKStep run directly.
 *   2. c = 1e-2: the adaptive ExtSTS solution must be accurate.
 *   3. D = 0 (no STS dynamics): ExtSTS must take the same steps as plain
 *      MRIStep with an exact inner integrator, so the STS estimate does not
 *      perturb the MRI estimate.
 *   4. c = 0, one step with the 2-stage default ExtSTS table: the ExtSTS error
 *      estimate (from fixed-step error accumulation) must approximately equal
 *      the LSRKStep embedding estimate for the same step.
 *   5. MRI H-Tol controllers must be rejected for ExtSTS methods.
 *   6. Option handling: the STS error estimate is disabled by default; the
 *      setter rejects non-ExtSTS objects; the command-line key matches the
 *      setter; and ARKodeSetDefaults disables the estimate again.
 *
 * The STS method is chosen by the first command-line argument (0 = RKC,
 * 1 = RKL; default 1).
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <arkode/arkode_lsrkstep.h>
#include <arkode/arkode_mristep.h>
#include <nvector/nvector_serial.h>
#include <sunadaptcontroller/sunadaptcontroller_mrihtol.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>

#include "arkode/arkode_mristep_impl.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)
#define FOUR SUN_RCONST(4.0)

#define NX   50
#define DIFF SUN_RCONST(0.1)

/* Precision-dependent parameters.  TOL is the integration tolerance and REFTOL
   the tolerance for the Test 2 reference solution.  HSTEP is the Test 4 step
   size: the RKC/RKL embedding estimates in Test 4 must be dominated by
   truncation error (which scales as HSTEP^3), not by roundoff in the STS
   stages (which grows with the stage count, so with HSTEP and NX). */
#if defined(SUNDIALS_SINGLE_PRECISION)
#define TOL    SUN_RCONST(1.0e-3)
#define REFTOL SUN_RCONST(1.0e-6)
#define HSTEP  SUN_RCONST(0.2)
#else
#define TOL    SUN_RCONST(1.0e-6)
#define REFTOL SUN_RCONST(1.0e-10)
#define HSTEP  SUN_RCONST(0.02)
#endif

/* Pass thresholds.  ERRTOL bounds the solution error in Tests 1 and 2.
   SAMETOL bounds the relative difference between ExtSTS and MRIStep in Test 3:
   with f^D = 0 the STS error estimate is pure roundoff, of size
   SUN_UNIT_ROUNDOFF / TOL in the weighted norm, so it perturbs the step sizes
   (and hence the solution) by an amount that scales with that ratio. */
#define ERRTOL  (SUN_RCONST(100.0) * TOL)
#define SAMETOL (SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF / TOL)

/* problem data */
typedef struct
{
  sunrealtype dx; /* mesh spacing */
  sunrealtype d;  /* diffusion coefficient */
  sunrealtype c;  /* advection coefficient */
  N_Vector work;  /* workspace for the full RHS */
} UserData;

/* f^D: diffusion */
static int f_diffusion(SUNDIALS_MAYBE_UNUSED sunrealtype t, N_Vector y,
                       N_Vector f, void* user_data)
{
  UserData* udata      = (UserData*)user_data;
  const sunrealtype* u = N_VGetArrayPointer(y);
  sunrealtype* r       = N_VGetArrayPointer(f);
  const sunrealtype k  = udata->d / (udata->dx * udata->dx);
  for (int i = 0; i < NX; i++)
  {
    const sunrealtype ul = (i > 0) ? u[i - 1] : ZERO;
    const sunrealtype ur = (i < NX - 1) ? u[i + 1] : ZERO;
    r[i]                 = k * (ul - TWO * u[i] + ur);
  }
  return 0;
}

/* f^E: advection */
static int f_advection(SUNDIALS_MAYBE_UNUSED sunrealtype t, N_Vector y,
                       N_Vector f, void* user_data)
{
  UserData* udata      = (UserData*)user_data;
  const sunrealtype* u = N_VGetArrayPointer(y);
  sunrealtype* r       = N_VGetArrayPointer(f);
  const sunrealtype k  = udata->c / (TWO * udata->dx);
  for (int i = 0; i < NX; i++)
  {
    const sunrealtype ul = (i > 0) ? u[i - 1] : ZERO;
    const sunrealtype ur = (i < NX - 1) ? u[i + 1] : ZERO;
    r[i]                 = k * (ur - ul);
  }
  return 0;
}

/* f^D + f^E, for LSRKStep direct run */
static int f_full(sunrealtype t, N_Vector y, N_Vector f, void* user_data)
{
  UserData* udata = (UserData*)user_data;
  int retval      = f_diffusion(t, y, f, user_data);
  if (retval != 0) { return retval; }
  retval = f_advection(t, y, udata->work, user_data);
  if (retval != 0) { return retval; }
  N_VLinearSum(ONE, f, ONE, udata->work, f);
  return 0;
}

/* dominant eigenvalue of f^D (advection only adds a small imaginary part) */
static int dom_eig(SUNDIALS_MAYBE_UNUSED sunrealtype t,
                   SUNDIALS_MAYBE_UNUSED N_Vector y,
                   SUNDIALS_MAYBE_UNUSED N_Vector fn, sunrealtype* lambdaR,
                   sunrealtype* lambdaI, void* user_data,
                   SUNDIALS_MAYBE_UNUSED N_Vector temp1,
                   SUNDIALS_MAYBE_UNUSED N_Vector temp2,
                   SUNDIALS_MAYBE_UNUSED N_Vector temp3)
{
  UserData* udata = (UserData*)user_data;
  *lambdaR        = -FOUR * udata->d / (udata->dx * udata->dx);
  *lambdaI        = ZERO;
  return 0;
}

/* initial condition and exact semi-discrete solution for c = 0 */
static void exact_solution(sunrealtype t, N_Vector y, const UserData* udata)
{
  const sunrealtype pi     = SUN_RCONST(4.0) * atan(ONE);
  const sunrealtype s      = sin(pi * udata->dx / TWO);
  const sunrealtype lambda = -FOUR * udata->d / (udata->dx * udata->dx) * s * s;
  sunrealtype* u           = N_VGetArrayPointer(y);
  for (int i = 0; i < NX; i++)
  {
    u[i] = exp(lambda * t) * sin(pi * (i + 1) * udata->dx);
  }
}

/* Exact inner stepper for f^D = 0: y' = sum_k F_k ((t - tshift)/tscale)^k */
static int exact_inner_evolve(MRIStepInnerStepper stepper, sunrealtype t0,
                              sunrealtype tout, N_Vector y)
{
  sunrealtype tshift, tscale;
  N_Vector* forcing;
  int nforcing;
  int retval = MRIStepInnerStepper_GetForcingData(stepper, &tshift, &tscale,
                                                  &forcing, &nforcing);
  if (retval != 0) { return retval; }

  /* int_t0^tout tau^k dt = tscale/(k+1) (tau(tout)^(k+1) - tau(t0)^(k+1)) */
  const sunrealtype tau0 = (t0 - tshift) / tscale;
  const sunrealtype tauf = (tout - tshift) / tscale;
  sunrealtype pow0 = tau0, powf = tauf;
  for (int k = 0; k < nforcing; k++)
  {
    const sunrealtype w = tscale / (k + 1) * (powf - pow0);
    N_VLinearSum(ONE, y, w, forcing[k], y);
    pow0 *= tau0;
    powf *= tauf;
  }
  return 0;
}

static int exact_inner_fullrhs(SUNDIALS_MAYBE_UNUSED MRIStepInnerStepper stepper,
                               SUNDIALS_MAYBE_UNUSED sunrealtype t,
                               SUNDIALS_MAYBE_UNUSED N_Vector y, N_Vector f,
                               SUNDIALS_MAYBE_UNUSED int mode)
{
  N_VConst(ZERO, f);
  return 0;
}

static int exact_inner_reset(SUNDIALS_MAYBE_UNUSED MRIStepInnerStepper stepper,
                             SUNDIALS_MAYBE_UNUSED sunrealtype tR,
                             SUNDIALS_MAYBE_UNUSED N_Vector yR)
{
  return 0;
}

/* Check a return flag; print an error message and return 1 on failure */
static int check_flag(int flag, const char* funcname)
{
  if (flag < 0)
  {
    fprintf(stderr, "ERROR: %s returned %i\n", funcname, flag);
    return 1;
  }
  return 0;
}

/* Integrators compared by run_adaptive: LSRKStep run directly, or ExtSTS with
   the STS error estimate left at its default, disabled or enabled with
   MRIStepSetUseExtSTSErrorEstimate, or enabled from the command line */
typedef enum
{
  RUN_LSRK,
  RUN_EXTSTS_DEFAULT,
  RUN_EXTSTS_OFF,
  RUN_EXTSTS_ON,
  RUN_EXTSTS_CLI
} RunType;

/* Integrate to tf with the integrator selected by run; return the step count
   and the final solution in y */
static int run_adaptive(RunType run, ARKODE_LSRKMethodType sts, UserData* udata,
                        sunrealtype tf, N_Vector y, long int* nst, SUNContext ctx)
{
  void* arkode_mem = NULL;
  void* sts_mem    = NULL;
  sunrealtype t    = ZERO;
  int fails        = 0;

  exact_solution(ZERO, y, udata);

  /* create the integrator; for ExtSTS, set the STS error option as requested */
  if (run != RUN_LSRK)
  {
    arkode_mem = MRIStepCreateExtSTS(f_diffusion, f_advection, NULL, ZERO, y,
                                     ctx);
    if (arkode_mem == NULL) { return 1; }
    fails += check_flag(MRIStepGetSTS(arkode_mem, &sts_mem), "MRIStepGetSTS");
    if (run == RUN_EXTSTS_OFF || run == RUN_EXTSTS_ON)
    {
      fails += check_flag(MRIStepSetUseExtSTSErrorEstimate(arkode_mem,
                                                           run == RUN_EXTSTS_ON),
                          "MRIStepSetUseExtSTSErrorEstimate");
    }
    else if (run == RUN_EXTSTS_CLI)
    {
      /* mimic "<program> arkode.use_extsts_error_estimate 1" */
      char arg0[] = "ark_test_extsts_stserror";
      char arg1[] = "arkode.use_extsts_error_estimate";
      char arg2[] = "1";
      char* cli[] = {arg0, arg1, arg2};
      fails += check_flag(ARKodeSetOptions(arkode_mem, NULL, NULL, 3, cli),
                          "ARKodeSetOptions");
    }
  }
  else
  {
    arkode_mem = LSRKStepCreateSTS(f_full, ZERO, y, ctx);
    if (arkode_mem == NULL) { return 1; }
    sts_mem = arkode_mem;
  }
  /* common settings: STS method, analytic dominant eigenvalue (computed once),
     tolerances and step limit on the outer integrator */
  fails += check_flag(ARKodeSetUserData(arkode_mem, udata),
                      "ARKodeSetUserData");
  fails += check_flag(ARKodeSetUserData(sts_mem, udata), "ARKodeSetUserData");
  fails += check_flag(LSRKStepSetSTSMethod(sts_mem, sts),
                      "LSRKStepSetSTSMethod");
  fails += check_flag(LSRKStepSetDomEigFn(sts_mem, dom_eig),
                      "LSRKStepSetDomEigFn");
  fails += check_flag(LSRKStepSetDomEigFrequency(sts_mem, 0),
                      "LSRKStepSetDomEigFrequency");
  fails += check_flag(ARKodeSStolerances(arkode_mem, TOL, TOL),
                      "ARKodeSStolerances");
  fails += check_flag(ARKodeSetMaxNumSteps(arkode_mem, 100000),
                      "ARKodeSetMaxNumSteps");
  /* integrate to tf (skipped if any setup call failed) */
  if (fails == 0)
  {
    fails += check_flag(ARKodeEvolve(arkode_mem, tf, y, &t, ARK_NORMAL),
                        "ARKodeEvolve");
    fails += check_flag(ARKodeGetNumSteps(arkode_mem, nst),
                        "ARKodeGetNumSteps");
  }
  ARKodeFree(&arkode_mem);
  return fails;
}

/* Test 1: c = 0, accuracy and step count relative to LSRKStep */
static int test_no_slow(ARKODE_LSRKMethodType sts, UserData* udata, SUNContext ctx)
{
  const sunrealtype tf = SUN_RCONST(5.0);
  long int nst_extsts = 0, nst_lsrk = 0;
  int fails = 0;

  N_Vector y      = N_VNew_Serial(NX, ctx);
  N_Vector yexact = N_VClone(y);

  udata->d = DIFF;
  udata->c = ZERO;
  exact_solution(tf, yexact, udata);

  /* max-norm errors of LSRKStep and ExtSTS against the exact solution */
  fails += run_adaptive(RUN_LSRK, sts, udata, tf, y, &nst_lsrk, ctx);
  N_VLinearSum(ONE, y, -ONE, yexact, y);
  const sunrealtype err_lsrk = N_VMaxNorm(y);

  fails += run_adaptive(RUN_EXTSTS_ON, sts, udata, tf, y, &nst_extsts, ctx);
  N_VLinearSum(ONE, y, -ONE, yexact, y);
  const sunrealtype err_extsts = N_VMaxNorm(y);

  printf("Test 1 (c = 0): LSRKStep steps = %li, error = %.3e; "
         "ExtSTS steps = %li, error = %.3e\n",
         nst_lsrk, (double)err_lsrk, nst_extsts, (double)err_extsts);

  /* ExtSTS must meet the tolerance with an LSRKStep-like step count */
  if (err_extsts > ERRTOL)
  {
    fprintf(stderr, "  FAIL: ExtSTS error %.3e exceeds %.1e\n",
            (double)err_extsts, (double)ERRTOL);
    fails++;
  }
  if ((2 * nst_extsts < nst_lsrk) || (nst_extsts > 2 * nst_lsrk))
  {
    fprintf(stderr, "  FAIL: ExtSTS steps %li not within a factor of 2 of %li\n",
            nst_extsts, nst_lsrk);
    fails++;
  }

  N_VDestroy(y);
  N_VDestroy(yexact);
  return fails;
}

/* Test 2: c = 1e-2, accuracy against a tight-tolerance LSRKStep reference */
static int test_small_slow(ARKODE_LSRKMethodType sts, UserData* udata,
                           SUNContext ctx)
{
  const sunrealtype tf = SUN_RCONST(5.0);
  long int nst_extsts  = 0;
  int fails            = 0;

  N_Vector y    = N_VNew_Serial(NX, ctx);
  N_Vector yref = N_VClone(y);

  udata->d = DIFF;
  udata->c = SUN_RCONST(1.0e-2);

  /* reference: LSRKStep on f^D + f^E with tolerances REFTOL */
  {
    sunrealtype t = ZERO;
    exact_solution(ZERO, yref, udata);
    void* arkode_mem = LSRKStepCreateSTS(f_full, ZERO, yref, ctx);
    fails += check_flag(ARKodeSetUserData(arkode_mem, udata),
                        "ARKodeSetUserData");
    fails += check_flag(LSRKStepSetSTSMethod(arkode_mem, sts),
                        "LSRKStepSetSTSMethod");
    fails += check_flag(LSRKStepSetDomEigFn(arkode_mem, dom_eig),
                        "LSRKStepSetDomEigFn");
    fails += check_flag(LSRKStepSetDomEigFrequency(arkode_mem, 0),
                        "LSRKStepSetDomEigFrequency");
    fails += check_flag(ARKodeSStolerances(arkode_mem, REFTOL, REFTOL),
                        "ARKodeSStolerances");
    fails += check_flag(ARKodeSetMaxNumSteps(arkode_mem, 1000000),
                        "ARKodeSetMaxNumSteps");
    fails += check_flag(ARKodeEvolve(arkode_mem, tf, yref, &t, ARK_NORMAL),
                        "ARKodeEvolve");
    ARKodeFree(&arkode_mem);
  }

  /* ExtSTS solution error against the reference */
  fails += run_adaptive(RUN_EXTSTS_ON, sts, udata, tf, y, &nst_extsts, ctx);
  N_VLinearSum(ONE, y, -ONE, yref, y);
  const sunrealtype err = N_VMaxNorm(y);

  printf("Test 2 (c = 1e-2): ExtSTS steps = %li, error = %.3e\n", nst_extsts,
         (double)err);
  if (err > ERRTOL)
  {
    fprintf(stderr, "  FAIL: ExtSTS error %.3e exceeds %.1e\n", (double)err,
            (double)ERRTOL);
    fails++;
  }

  N_VDestroy(y);
  N_VDestroy(yref);
  return fails;
}

/* Test 3: D = 0, ExtSTS vs. plain MRIStep with an exact inner integrator */
static int test_no_sts(ARKODE_LSRKMethodType sts, UserData* udata, SUNContext ctx)
{
  const sunrealtype tf = SUN_RCONST(0.5);
  long int nst_extsts = 0, nst_mri = 0;
  sunrealtype t = ZERO;
  int fails     = 0;

  N_Vector y_extsts         = N_VNew_Serial(NX, ctx);
  N_Vector y_mri            = N_VClone(y_extsts);
  MRIStepInnerStepper inner = NULL;
  MRIStepCoupling MRIC      = NULL;
  void* sts_mem             = NULL;

  udata->d = ZERO;
  udata->c = ONE;
  exact_solution(ZERO, y_extsts, udata);
  N_VScale(ONE, y_extsts, y_mri);

  /* ExtSTS with f^D = 0 (the dominant eigenvalue is then zero) */
  void* extsts_mem = MRIStepCreateExtSTS(f_diffusion, f_advection, NULL, ZERO,
                                         y_extsts, ctx);
  if (extsts_mem == NULL) { return 1; }
  fails += check_flag(MRIStepGetSTS(extsts_mem, &sts_mem), "MRIStepGetSTS");
  fails += check_flag(MRIStepSetUseExtSTSErrorEstimate(extsts_mem, SUNTRUE),
                      "MRIStepSetUseExtSTSErrorEstimate");
  fails += check_flag(ARKodeSetUserData(extsts_mem, udata),
                      "ARKodeSetUserData");
  fails += check_flag(ARKodeSetUserData(sts_mem, udata), "ARKodeSetUserData");
  fails += check_flag(LSRKStepSetSTSMethod(sts_mem, sts),
                      "LSRKStepSetSTSMethod");
  fails += check_flag(LSRKStepSetDomEigFn(sts_mem, dom_eig),
                      "LSRKStepSetDomEigFn");
  fails += check_flag(ARKodeSStolerances(extsts_mem, TOL, TOL),
                      "ARKodeSStolerances");
  fails += check_flag(ARKodeSetMaxNumSteps(extsts_mem, 100000),
                      "ARKodeSetMaxNumSteps");

  /* plain MRIStep with the same coupling table and an exact inner stepper */
  fails += check_flag(MRIStepInnerStepper_Create(ctx, &inner),
                      "MRIStepInnerStepper_Create");
  fails += check_flag(MRIStepInnerStepper_SetEvolveFn(inner, exact_inner_evolve),
                      "MRIStepInnerStepper_SetEvolveFn");
  fails += check_flag(MRIStepInnerStepper_SetFullRhsFn(inner, exact_inner_fullrhs),
                      "MRIStepInnerStepper_SetFullRhsFn");
  fails += check_flag(MRIStepInnerStepper_SetResetFn(inner, exact_inner_reset),
                      "MRIStepInnerStepper_SetResetFn");
  void* mri_mem = MRIStepCreate(f_advection, NULL, ZERO, y_mri, inner, ctx);
  if (mri_mem == NULL) { return fails + 1; }
  MRIC = MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_EXP_GKC21);
  fails += check_flag(MRIStepSetCoupling(mri_mem, MRIC), "MRIStepSetCoupling");
  MRIStepCoupling_Free(MRIC);
  fails += check_flag(ARKodeSetUserData(mri_mem, udata), "ARKodeSetUserData");
  fails += check_flag(ARKodeSStolerances(mri_mem, TOL, TOL),
                      "ARKodeSStolerances");
  fails += check_flag(ARKodeSetMaxNumSteps(mri_mem, 100000),
                      "ARKodeSetMaxNumSteps");

  /* integrate both to tf and compare step counts and final solutions */
  if (fails == 0)
  {
    fails += check_flag(ARKodeEvolve(extsts_mem, tf, y_extsts, &t, ARK_NORMAL),
                        "ARKodeEvolve");
    fails += check_flag(ARKodeEvolve(mri_mem, tf, y_mri, &t, ARK_NORMAL),
                        "ARKodeEvolve");
    fails += check_flag(ARKodeGetNumSteps(extsts_mem, &nst_extsts),
                        "ARKodeGetNumSteps");
    fails += check_flag(ARKodeGetNumSteps(mri_mem, &nst_mri),
                        "ARKodeGetNumSteps");
  }

  N_VLinearSum(ONE, y_extsts, -ONE, y_mri, y_mri);
  const sunrealtype diff = N_VMaxNorm(y_mri) / N_VMaxNorm(y_extsts);
  printf("Test 3 (D = 0): MRIStep steps = %li, ExtSTS steps = %li, "
         "relative difference = %.3e\n",
         nst_mri, nst_extsts, (double)diff);
  if (nst_extsts != nst_mri)
  {
    fprintf(stderr, "  FAIL: step counts differ\n");
    fails++;
  }
  if (diff > SAMETOL)
  {
    fprintf(stderr, "  FAIL: solutions differ by %.3e (limit %.3e)\n",
            (double)diff, (double)SAMETOL);
    fails++;
  }

  ARKodeFree(&extsts_mem);
  ARKodeFree(&mri_mem);
  MRIStepInnerStepper_Free(&inner);
  N_VDestroy(y_extsts);
  N_VDestroy(y_mri);
  return fails;
}

/* Test 4: c = 0, one step of size h: ExtSTS (fixed step, accumulated error)
   vs. LSRKStep (adaptive, accumulated error over one accepted attempt) */
static int test_one_step(ARKODE_LSRKMethodType sts, UserData* udata,
                         SUNContext ctx)
{
  const sunrealtype h = HSTEP;
  sunrealtype t = ZERO, dsm_extsts = ZERO, dsm_lsrk = ZERO;
  long int nattempts = 0;
  void* sts_mem      = NULL;
  int fails          = 0;

  N_Vector y = N_VNew_Serial(NX, ctx);
  udata->d   = DIFF;
  udata->c   = ZERO;

  /* ExtSTS (default GKC21 explicit table), fixed step with error estimation */
  exact_solution(ZERO, y, udata);
  void* extsts_mem = MRIStepCreateExtSTS(f_diffusion, f_advection, NULL, ZERO,
                                         y, ctx);
  if (extsts_mem == NULL) { return 1; }
  fails += check_flag(MRIStepGetSTS(extsts_mem, &sts_mem), "MRIStepGetSTS");
  fails += check_flag(ARKodeSetUserData(extsts_mem, udata),
                      "ARKodeSetUserData");
  fails += check_flag(ARKodeSetUserData(sts_mem, udata), "ARKodeSetUserData");
  fails += check_flag(LSRKStepSetSTSMethod(sts_mem, sts),
                      "LSRKStepSetSTSMethod");
  fails += check_flag(LSRKStepSetDomEigFn(sts_mem, dom_eig),
                      "LSRKStepSetDomEigFn");
  fails += check_flag(ARKodeSStolerances(extsts_mem, TOL, TOL),
                      "ARKodeSStolerances");
  fails += check_flag(MRIStepSetUseExtSTSErrorEstimate(extsts_mem, SUNTRUE),
                      "MRIStepSetUseExtSTSErrorEstimate");
  fails += check_flag(ARKodeSetFixedStep(extsts_mem, h), "ARKodeSetFixedStep");
  fails +=
    check_flag(ARKodeSetAccumulatedErrorType(extsts_mem, ARK_ACCUMERROR_MAX),
               "ARKodeSetAccumulatedErrorType");
  if (fails == 0)
  {
    fails += check_flag(ARKodeEvolve(extsts_mem, h, y, &t, ARK_ONE_STEP),
                        "ARKodeEvolve");
    fails += check_flag(ARKodeGetAccumulatedError(extsts_mem, &dsm_extsts),
                        "ARKodeGetAccumulatedError");
    /* convert the accumulated error to the WRMS dsm scale */
    dsm_extsts /= TOL;
  }
  ARKodeFree(&extsts_mem);

  /* LSRKStep run directly, one adaptive step with initial step h */
  exact_solution(ZERO, y, udata);
  void* lsrk_mem = LSRKStepCreateSTS(f_full, ZERO, y, ctx);
  if (lsrk_mem == NULL) { return fails + 1; }
  fails += check_flag(ARKodeSetUserData(lsrk_mem, udata), "ARKodeSetUserData");
  fails += check_flag(LSRKStepSetSTSMethod(lsrk_mem, sts),
                      "LSRKStepSetSTSMethod");
  fails += check_flag(LSRKStepSetDomEigFn(lsrk_mem, dom_eig),
                      "LSRKStepSetDomEigFn");
  fails += check_flag(ARKodeSStolerances(lsrk_mem, TOL, TOL),
                      "ARKodeSStolerances");
  fails += check_flag(ARKodeSetInitStep(lsrk_mem, h), "ARKodeSetInitStep");
  fails += check_flag(ARKodeSetAccumulatedErrorType(lsrk_mem, ARK_ACCUMERROR_MAX),
                      "ARKodeSetAccumulatedErrorType");
  if (fails == 0)
  {
    fails += check_flag(ARKodeEvolve(lsrk_mem, h, y, &t, ARK_ONE_STEP),
                        "ARKodeEvolve");
    fails += check_flag(ARKodeGetNumStepAttempts(lsrk_mem, &nattempts),
                        "ARKodeGetNumStepAttempts");
    fails += check_flag(ARKodeGetAccumulatedError(lsrk_mem, &dsm_lsrk),
                        "ARKodeGetAccumulatedError");
    dsm_lsrk /= TOL;
  }
  ARKodeFree(&lsrk_mem);

  /* The RKC/RKL embedded estimate for a step of size dt scales as dt^3, and
     the GKC21 fast solves have lengths theta_k h with theta_1 = 2 - sqrt(2)
     and theta_2 = sqrt(2) - 1. The time-weighted RMS of the two estimates is
     therefore sqrt(theta_1^7 + theta_2^7) times the LSRKStep estimate for a
     single step of size h. This must be updated if the LSRKStep embeddings
     change order. */
  const sunrealtype theta1 = SUN_RCONST(2.0) - SUNRsqrt(SUN_RCONST(2.0));
  const sunrealtype theta2 = ONE - theta1;
  const sunrealtype expected =
    SUNRsqrt(SUNRpowerI(theta1, 7) + SUNRpowerI(theta2, 7));
  const sunrealtype ratio = dsm_extsts / dsm_lsrk;
  printf("Test 4 (one step, h = %g): LSRKStep dsm = %.3e (%li attempt), "
         "ExtSTS dsm = %.3e, ratio = %.4f (expected %.4f)\n",
         (double)h, (double)dsm_lsrk, nattempts, (double)dsm_extsts,
         (double)ratio, (double)expected);
  /* the LSRKStep estimate is only comparable if its first step was accepted */
  if (nattempts != 1)
  {
    fprintf(stderr, "  FAIL: LSRKStep step was rejected; reduce h\n");
    fails++;
  }
  if (SUNRabs(ratio - expected) > SUN_RCONST(0.1) * expected)
  {
    fprintf(stderr,
            "  FAIL: estimate ratio %.4f differs from expected %.4f by more "
            "than 10%%\n",
            (double)ratio, (double)expected);
    fails++;
  }

  N_VDestroy(y);
  return fails;
}

/* Test 5: MRI H-Tol controllers are rejected for ExtSTS methods */
static int test_htol_rejected(UserData* udata, SUNContext ctx)
{
  sunrealtype t = ZERO;
  void* sts_mem = NULL;
  int fails     = 0;

  N_Vector y = N_VNew_Serial(NX, ctx);
  udata->d   = DIFF;
  udata->c   = ZERO;
  exact_solution(ZERO, y, udata);

  void* extsts_mem = MRIStepCreateExtSTS(f_diffusion, f_advection, NULL, ZERO,
                                         y, ctx);
  if (extsts_mem == NULL) { return 1; }
  fails += check_flag(MRIStepGetSTS(extsts_mem, &sts_mem), "MRIStepGetSTS");
  fails += check_flag(ARKodeSetUserData(extsts_mem, udata),
                      "ARKodeSetUserData");
  fails += check_flag(ARKodeSetUserData(sts_mem, udata), "ARKodeSetUserData");
  fails += check_flag(LSRKStepSetDomEigFn(sts_mem, dom_eig),
                      "LSRKStepSetDomEigFn");
  fails += check_flag(ARKodeSStolerances(extsts_mem, TOL, TOL),
                      "ARKodeSStolerances");

  /* attach an MRI H-Tol controller; the failure is reported at the first step */
  SUNAdaptController HControl   = SUNAdaptController_I(ctx);
  SUNAdaptController TolControl = SUNAdaptController_I(ctx);
  SUNAdaptController htol = SUNAdaptController_MRIHTol(HControl, TolControl, ctx);
  if (HControl == NULL || TolControl == NULL || htol == NULL) { fails++; }
  else
  {
    fails += check_flag(ARKodeSetAdaptController(extsts_mem, htol),
                        "ARKodeSetAdaptController");
  }

  if (fails == 0)
  {
    const int flag = ARKodeEvolve(extsts_mem, ONE, y, &t, ARK_ONE_STEP);
    printf("Test 5 (H-Tol controller): ARKodeEvolve returned %i\n", flag);
    if (flag != ARK_ILL_INPUT)
    {
      fprintf(stderr, "  FAIL: expected ARK_ILL_INPUT (%i)\n", ARK_ILL_INPUT);
      fails++;
    }
  }

  ARKodeFree(&extsts_mem);
  SUNAdaptController_Destroy(htol);
  SUNAdaptController_Destroy(HControl);
  SUNAdaptController_Destroy(TolControl);
  N_VDestroy(y);
  return fails;
}

/* Test 6: option handling for MRIStepSetUseExtSTSErrorEstimate */
static int test_option(ARKODE_LSRKMethodType sts, UserData* udata, SUNContext ctx)
{
  const sunrealtype tf = SUN_RCONST(5.0);
  long int nst_default = 0, nst_off = 0, nst_on = 0, nst_cli = 0;
  MRIStepInnerStepper inner = NULL;
  void* sts_mem             = NULL;
  int fails                 = 0;
  int flag;

  N_Vector y_default = N_VNew_Serial(NX, ctx);
  N_Vector y_off     = N_VClone(y_default);
  N_Vector y_on      = N_VClone(y_default);
  N_Vector y_cli     = N_VClone(y_default);

  udata->d = DIFF;
  udata->c = ZERO;

  /* (a) default off: never calling the setter matches disabling it, and
     (c) the command-line key matches enabling it with the setter.  Identical
     code paths are expected here, so the solutions are compared exactly. */
  fails += run_adaptive(RUN_EXTSTS_DEFAULT, sts, udata, tf, y_default,
                        &nst_default, ctx);
  fails += run_adaptive(RUN_EXTSTS_OFF, sts, udata, tf, y_off, &nst_off, ctx);
  fails += run_adaptive(RUN_EXTSTS_ON, sts, udata, tf, y_on, &nst_on, ctx);
  fails += run_adaptive(RUN_EXTSTS_CLI, sts, udata, tf, y_cli, &nst_cli, ctx);

  /* differences between runs that should be identical */
  N_VLinearSum(ONE, y_default, -ONE, y_off, y_off);
  N_VLinearSum(ONE, y_on, -ONE, y_cli, y_cli);
  const sunrealtype diff_off = N_VMaxNorm(y_off);
  const sunrealtype diff_cli = N_VMaxNorm(y_cli);
  printf("Test 6 (option): steps default = %li, off = %li, on = %li, "
         "command line = %li; max |default - off| = %.3e, "
         "max |on - command line| = %.3e\n",
         nst_default, nst_off, nst_on, nst_cli, (double)diff_off,
         (double)diff_cli);
  if (nst_default != nst_off || diff_off != ZERO)
  {
    fprintf(stderr, "  FAIL: default run differs from disabled run\n");
    fails++;
  }
  if (nst_on != nst_cli || diff_cli != ZERO)
  {
    fprintf(stderr, "  FAIL: command-line run differs from enabled run\n");
    fails++;
  }
  if (nst_on == nst_default)
  {
    fprintf(stderr, "  FAIL: enabled run has the same step count as default\n");
    fails++;
  }

  /* (b) setter return values: NULL memory, ExtSTS object, plain MRIStep */
  flag = MRIStepSetUseExtSTSErrorEstimate(NULL, SUNTRUE);
  if (flag != ARK_MEM_NULL)
  {
    fprintf(stderr, "  FAIL: setter with NULL memory returned %i\n", flag);
    fails++;
  }

  /* the setter succeeds on an ExtSTS object */
  exact_solution(ZERO, y_on, udata);
  void* extsts_mem = MRIStepCreateExtSTS(f_diffusion, f_advection, NULL, ZERO,
                                         y_on, ctx);
  if (extsts_mem == NULL) { return fails + 1; }
  fails += check_flag(MRIStepGetSTS(extsts_mem, &sts_mem), "MRIStepGetSTS");
  flag = MRIStepSetUseExtSTSErrorEstimate(extsts_mem, SUNTRUE);
  if (flag != ARK_SUCCESS)
  {
    fprintf(stderr, "  FAIL: setter on ExtSTS object returned %i\n", flag);
    fails++;
  }

  /* (d) ARKodeSetDefaults disables the estimate again */
  fails += check_flag(ARKodeSetDefaults(extsts_mem), "ARKodeSetDefaults");
  ARKodeMRIStepMem step_mem = (ARKodeMRIStepMem)((ARKodeMem)extsts_mem)->step_mem;
  if (step_mem->extsts_use_err_est)
  {
    fprintf(stderr, "  FAIL: ARKodeSetDefaults did not disable the option\n");
    fails++;
  }
  ARKodeFree(&extsts_mem);

  /* both the setter and the command-line key fail on a plain MRIStep object */
  fails += check_flag(MRIStepInnerStepper_Create(ctx, &inner),
                      "MRIStepInnerStepper_Create");
  fails += check_flag(MRIStepInnerStepper_SetEvolveFn(inner, exact_inner_evolve),
                      "MRIStepInnerStepper_SetEvolveFn");
  fails += check_flag(MRIStepInnerStepper_SetFullRhsFn(inner, exact_inner_fullrhs),
                      "MRIStepInnerStepper_SetFullRhsFn");
  void* mri_mem = MRIStepCreate(f_advection, NULL, ZERO, y_on, inner, ctx);
  if (mri_mem == NULL) { return fails + 1; }
  flag = MRIStepSetUseExtSTSErrorEstimate(mri_mem, SUNTRUE);
  if (flag != ARK_ILL_INPUT)
  {
    fprintf(stderr, "  FAIL: setter on plain MRIStep object returned %i\n", flag);
    fails++;
  }
  {
    char arg0[] = "ark_test_extsts_stserror";
    char arg1[] = "arkode.use_extsts_error_estimate";
    char arg2[] = "1";
    char* cli[] = {arg0, arg1, arg2};
    flag        = ARKodeSetOptions(mri_mem, NULL, NULL, 3, cli);
    if (flag != ARK_ILL_INPUT)
    {
      fprintf(stderr,
              "  FAIL: command-line key on plain MRIStep object returned %i\n",
              flag);
      fails++;
    }
  }
  ARKodeFree(&mri_mem);
  MRIStepInnerStepper_Free(&inner);

  N_VDestroy(y_default);
  N_VDestroy(y_off);
  N_VDestroy(y_on);
  N_VDestroy(y_cli);
  return fails;
}

int main(int argc, char* argv[])
{
  SUNContext ctx = NULL;
  UserData udata;
  int fails = 0;

  const ARKODE_LSRKMethodType sts =
    (argc > 1 && atoi(argv[1]) == 0) ? ARKODE_LSRK_RKC_2 : ARKODE_LSRK_RKL_2;
  printf("ExtSTS STS-error estimate test, %s\n",
         (sts == ARKODE_LSRK_RKC_2) ? "RKC" : "RKL");

  if (SUNContext_Create(SUN_COMM_NULL, &ctx)) { return 1; }

  udata.dx   = ONE / (NX + 1);
  udata.d    = DIFF;
  udata.c    = ZERO;
  udata.work = N_VNew_Serial(NX, ctx);

  fails += test_no_slow(sts, &udata, ctx);
  fails += test_small_slow(sts, &udata, ctx);
  fails += test_no_sts(sts, &udata, ctx);
  fails += test_one_step(sts, &udata, ctx);
  fails += test_htol_rejected(&udata, ctx);
  fails += test_option(sts, &udata, ctx);

  N_VDestroy(udata.work);
  SUNContext_Free(&ctx);

  if (fails) { printf("FAIL: %i test failure(s)\n", fails); }
  else { printf("SUCCESS\n"); }
  return (fails) ? 1 : 0;
}
