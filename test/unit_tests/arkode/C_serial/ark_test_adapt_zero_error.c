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
 * Unit test for the handling of exactly-zero local error estimates by the
 * history-based SUNAdaptController modules, and for the ARKODE check that a
 * controller does not return NaN.
 *
 * A zero error stored in the controller history must restart the history, so
 * that the I controller is used until enough positive errors have been stored
 * again. A zero current error must give an infinite step (bounded by ARKODE).
 *
 * The test has three parts:
 *   1. Controller-level checks that call the controllers directly.
 *   2. ARKODE-level checks that wrap a controller in a pass-through controller
 *      which replaces the error passed to UpdateH by zero on one accepted step,
 *      and then verify the controller output on the following steps.
 *   3. An ARKODE-level check that a controller returning NaN causes
 *      ARKodeEvolve to return ARK_CONTROLLER_ERR.
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>
#include <sunadaptcontroller/sunadaptcontroller_imexgus.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>
#include <sundials/sundials_math.h>

/* relative tolerance for comparing the controller output against the formulas */
#define RTOL (SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF)

/* error bias used for all controllers (not the default, so that it is tested) */
#define BIAS SUN_RCONST(1.5)

/* method order passed to the controller in the controller-level checks */
#define ORDER 3

/* -----------------------------------------------------------------------------
 * Description of each controller: the constructor and the gains, from which
 * the expected step sizes are computed independently of the implementation.
 * ---------------------------------------------------------------------------*/

typedef SUNAdaptController (*CreateFn)(SUNContext sunctx);

typedef struct
{
  const char* name;
  CreateFn create;
  sunbooleantype imex; /* SUNFALSE: Soderlind family, SUNTRUE: ImExGus     */
  sunrealtype k[5];    /* Soderlind: k1..k5, ImExGus: k1e, k2e, k1i, k2i  */
  int historysize;     /* number of stored errors used by the controller */
} ControllerCase;

/* Soderlind controller with all five gains nonzero and of mixed sign */
static SUNAdaptController create_soderlind_custom(SUNContext sunctx)
{
  SUNAdaptController C = SUNAdaptController_Soderlind(sunctx);
  if (C == NULL) { return NULL; }
  if (SUNAdaptController_SetParams_Soderlind(C, SUN_RCONST(0.123),
                                             SUN_RCONST(-0.456), SUN_RCONST(0.789),
                                             SUN_RCONST(-1.0), SUN_RCONST(-2.0)))
  {
    SUNAdaptController_Destroy(C);
    return NULL;
  }
  return C;
}

/* ImExGus controller whose k2 gains have the opposite sign of its k1 gains, so
   that a zero current error would give 0 * Inf = NaN in the history formula */
static SUNAdaptController create_imexgus_mixed(SUNContext sunctx)
{
  SUNAdaptController C = SUNAdaptController_ImExGus(sunctx);
  if (C == NULL) { return NULL; }
  if (SUNAdaptController_SetParams_ImExGus(C, SUN_RCONST(0.4), SUN_RCONST(-0.2),
                                           SUN_RCONST(0.95), SUN_RCONST(-0.5)))
  {
    SUNAdaptController_Destroy(C);
    return NULL;
  }
  return C;
}

/* The gains below repeat the defaults of each constructor, mapped onto the
   Soderlind form h' = h e1^(-k1/q) ep^(-k2/q) epp^(-k3/q) (h/hp)^k4
   (hp/hpp)^k5 with q = p + 1 (see sunadaptcontroller_soderlind.c). */
static const ControllerCase controller_cases[] = {
  {"PID",
   SUNAdaptController_PID,
   SUNFALSE,
   {SUN_RCONST(0.58), SUN_RCONST(-0.21), SUN_RCONST(0.1), SUN_RCONST(0.0),
    SUN_RCONST(0.0)},
   2},
  {"PI",
   SUNAdaptController_PI,
   SUNFALSE,
   {SUN_RCONST(0.8), SUN_RCONST(-0.31), SUN_RCONST(0.0), SUN_RCONST(0.0),
    SUN_RCONST(0.0)},
   1},
  {"I",
   SUNAdaptController_I,
   SUNFALSE,
   {SUN_RCONST(1.0), SUN_RCONST(0.0), SUN_RCONST(0.0), SUN_RCONST(0.0),
    SUN_RCONST(0.0)},
   0},
  /* ExpGus(k1, k2) = Soderlind(k1 + k2, -k2, 0, 0, 0) */
  {"ExpGus",
   SUNAdaptController_ExpGus,
   SUNFALSE,
   {SUN_RCONST(0.367) + SUN_RCONST(0.268), -SUN_RCONST(0.268), SUN_RCONST(0.0),
    SUN_RCONST(0.0), SUN_RCONST(0.0)},
   1},
  /* ImpGus(k1, k2) = Soderlind(k1 + k2, -k2, 0, 1, 0) */
  {"ImpGus",
   SUNAdaptController_ImpGus,
   SUNFALSE,
   {SUN_RCONST(0.98) + SUN_RCONST(0.95), -SUN_RCONST(0.95), SUN_RCONST(0.0),
    SUN_RCONST(1.0), SUN_RCONST(0.0)},
   1},
  {"H0211",
   SUNAdaptController_H0211,
   SUNFALSE,
   {SUN_RCONST(0.5), SUN_RCONST(0.5), SUN_RCONST(0.0), SUN_RCONST(-0.5),
    SUN_RCONST(0.0)},
   1},
  {"H0321",
   SUNAdaptController_H0321,
   SUNFALSE,
   {SUN_RCONST(1.25), SUN_RCONST(0.5), SUN_RCONST(-0.75), SUN_RCONST(0.25),
    SUN_RCONST(0.75)},
   2},
  {"H211",
   SUNAdaptController_H211,
   SUNFALSE,
   {SUN_RCONST(0.25), SUN_RCONST(0.25), SUN_RCONST(0.0), SUN_RCONST(-0.25),
    SUN_RCONST(0.0)},
   1},
  {"H312",
   SUNAdaptController_H312,
   SUNFALSE,
   {SUN_RCONST(1.0) / SUN_RCONST(8.0), SUN_RCONST(0.25),
    SUN_RCONST(1.0) / SUN_RCONST(8.0), -SUN_RCONST(3.0) / SUN_RCONST(8.0),
    -SUN_RCONST(1.0) / SUN_RCONST(8.0)},
   2},
  {"Soderlind (custom)",
   create_soderlind_custom,
   SUNFALSE,
   {SUN_RCONST(0.123), SUN_RCONST(-0.456), SUN_RCONST(0.789), SUN_RCONST(-1.0),
    SUN_RCONST(-2.0)},
   2},
  {"ImExGus",
   SUNAdaptController_ImExGus,
   SUNTRUE,
   {SUN_RCONST(0.367), SUN_RCONST(0.268), SUN_RCONST(0.95), SUN_RCONST(0.95),
    SUN_RCONST(0.0)},
   1},
  {"ImExGus (mixed-sign gains)",
   create_imexgus_mixed,
   SUNTRUE,
   {SUN_RCONST(0.4), SUN_RCONST(-0.2), SUN_RCONST(0.95), SUN_RCONST(-0.5),
    SUN_RCONST(0.0)},
   1},
};

static const int num_controller_cases =
  (int)(sizeof(controller_cases) / sizeof(controller_cases[0]));

/* The controllers run inside ARKODE in part 2 (a subset of the above that is
   reasonable for a nonstiff ERK integration) */
static const char* const arkode_case_names[] = {"PID",  "PI",   "ExpGus",
                                                "H211", "H312", "ImExGus"};

static const int num_arkode_cases =
  (int)(sizeof(arkode_case_names) / sizeof(arkode_case_names[0]));

/* -----------------------------------------------------------------------------
 * Expected step sizes
 * ---------------------------------------------------------------------------*/

/* I controller, used when there is no usable history:
   h' = h (bias dsm)^(-1/(p+1)) */
static sunrealtype expected_i_law(sunrealtype h, int p, sunrealtype dsm)
{
  return h * SUNRpowerR(BIAS * dsm, -SUN_RCONST(1.0) / (p + 1));
}

/* Controller formula with full history: (h, dsm) are the current step and
   error, and (hp, dsmp), (hpp, dsmpp) the previous two stored ones */
static sunrealtype expected_full(const ControllerCase* cc, sunrealtype h, int p,
                                 sunrealtype dsm, sunrealtype hp,
                                 sunrealtype dsmp, sunrealtype hpp,
                                 sunrealtype dsmpp)
{
  const sunrealtype ord = (sunrealtype)(p + 1);
  const sunrealtype e1  = BIAS * dsm;
  const sunrealtype ep  = BIAS * dsmp;

  if (cc->imex)
  {
    const sunrealtype e2   = e1 / ep;
    const sunrealtype hexp = SUNRpowerR(e1, -cc->k[0] / ord) *
                             SUNRpowerR(e2, -cc->k[1] / ord);
    const sunrealtype himp = (h / hp) * SUNRpowerR(e1, -cc->k[2] / ord) *
                             SUNRpowerR(e2, -cc->k[3] / ord);
    return h * SUNMIN(himp, hexp);
  }

  sunrealtype hnew = h * SUNRpowerR(e1, -cc->k[0] / ord);
  if (cc->historysize > 0)
  {
    hnew *= SUNRpowerR(ep, -cc->k[1] / ord) * SUNRpowerR(h / hp, cc->k[3]);
  }
  if (cc->historysize > 1)
  {
    const sunrealtype epp = BIAS * dsmpp;
    hnew *= SUNRpowerR(epp, -cc->k[2] / ord) * SUNRpowerR(hp / hpp, cc->k[4]);
  }
  return hnew;
}

/* Compare a step size against its expected value; returns 1 on a mismatch */
static int check_value(const char* name, const char* what, sunrealtype actual,
                       sunrealtype expected)
{
  if (!isfinite(actual) || actual <= SUN_RCONST(0.0) ||
      SUNRCompareTol(actual, expected, RTOL))
  {
    fprintf(stderr,
            "FAIL [%s] %s: hnew = " SUN_FORMAT_G ", expected " SUN_FORMAT_G "\n",
            name, what, actual, expected);
    return 1;
  }
  return 0;
}

/* Check that a step size is +Inf; returns 1 otherwise */
static int check_inf(const char* name, const char* what, sunrealtype actual)
{
  if (!(isinf(actual) && actual > SUN_RCONST(0.0)))
  {
    fprintf(stderr, "FAIL [%s] %s: hnew = " SUN_FORMAT_G ", expected +Inf\n",
            name, what, actual);
    return 1;
  }
  return 0;
}

/* -----------------------------------------------------------------------------
 * Part 1: controller-level checks
 * ---------------------------------------------------------------------------*/

/* Fixed sequence of positive step sizes and errors */
#define NSEQ 8
static const sunrealtype hseq[NSEQ]   = {SUN_RCONST(0.1),  SUN_RCONST(0.13),
                                         SUN_RCONST(0.11), SUN_RCONST(0.12),
                                         SUN_RCONST(0.09), SUN_RCONST(0.105),
                                         SUN_RCONST(0.14), SUN_RCONST(0.1)};
static const sunrealtype dsmseq[NSEQ] = {SUN_RCONST(0.31), SUN_RCONST(0.72),
                                         SUN_RCONST(0.45), SUN_RCONST(0.88),
                                         SUN_RCONST(0.27), SUN_RCONST(0.63),
                                         SUN_RCONST(0.5),  SUN_RCONST(0.95)};

/* Reset the controller and store the first three positive errors, which fills
   the history of every controller */
static int fill_history(SUNAdaptController C)
{
  if (SUNAdaptController_Reset(C)) { return 1; }
  for (int i = 0; i < 3; i++)
  {
    if (SUNAdaptController_UpdateH(C, hseq[i], dsmseq[i])) { return 1; }
  }
  return 0;
}

static int test_controller(const ControllerCase* cc, SUNContext sunctx)
{
  int fails = 0;
  sunrealtype hnew;

  SUNAdaptController C = cc->create(sunctx);
  if (C == NULL)
  {
    fprintf(stderr, "FAIL [%s] could not create controller\n", cc->name);
    return 1;
  }
  if (SUNAdaptController_SetErrorBias(C, BIAS))
  {
    fprintf(stderr, "FAIL [%s] could not set error bias\n", cc->name);
    SUNAdaptController_Destroy(C);
    return 1;
  }

  /* (a) No change for positive errors: with a full history of positive
     errors, the step matches the controller formula */
  fails += fill_history(C);
  SUNAdaptController_EstimateStep(C, hseq[3], ORDER, dsmseq[3], &hnew);
  fails += check_value(cc->name, "positive errors", hnew,
                       expected_full(cc, hseq[3], ORDER, dsmseq[3], hseq[2],
                                     dsmseq[2], hseq[1], dsmseq[1]));

  /* (b) Zero current error with a positive history gives +Inf */
  SUNAdaptController_EstimateStep(C, hseq[3], ORDER, SUN_RCONST(0.0), &hnew);
  fails += check_inf(cc->name, "zero current error", hnew);

  /* (c) Zero stored error, then positive errors: the I controller is used
     until historysize positive errors have been stored, after which the
     controller formula applies to the new history */
  fails += fill_history(C);
  SUNAdaptController_UpdateH(C, hseq[3], SUN_RCONST(0.0));
  for (int j = 0; j <= cc->historysize; j++)
  {
    const int i = 4 + j;
    char what[64];
    SUNAdaptController_EstimateStep(C, hseq[i], ORDER, dsmseq[i], &hnew);
    if (j < cc->historysize)
    {
      snprintf(what, sizeof(what), "step %i after zero error (I law)", j + 1);
      fails += check_value(cc->name, what, hnew,
                           expected_i_law(hseq[i], ORDER, dsmseq[i]));
    }
    else
    {
      snprintf(what, sizeof(what), "step %i after zero error (full history)",
               j + 1);
      fails +=
        check_value(cc->name, what, hnew,
                    expected_full(cc, hseq[i], ORDER, dsmseq[i], hseq[i - 1],
                                  dsmseq[i - 1], hseq[i - 2], dsmseq[i - 2]));
    }
    SUNAdaptController_UpdateH(C, hseq[i], dsmseq[i]);
  }

  /* (d) Repeated zero errors: a zero current error after a stored zero gives
     +Inf, and a second stored zero still restarts the history */
  fails += fill_history(C);
  SUNAdaptController_UpdateH(C, hseq[3], SUN_RCONST(0.0));
  SUNAdaptController_EstimateStep(C, hseq[4], ORDER, SUN_RCONST(0.0), &hnew);
  fails += check_inf(cc->name, "repeated zero error", hnew);
  SUNAdaptController_UpdateH(C, hseq[4], SUN_RCONST(0.0));
  SUNAdaptController_EstimateStep(C, hseq[5], ORDER, dsmseq[5], &hnew);
  fails += check_value(cc->name, "step after two zero errors (I law)", hnew,
                       expected_i_law(hseq[5], ORDER, dsmseq[5]));

  SUNAdaptController_Destroy(C);
  return fails;
}

/* -----------------------------------------------------------------------------
 * Pass-through controller for the ARKODE-level checks: forwards all calls to
 * an inner controller, records the inputs and output of the EstimateStep call
 * for each accepted step, and can either replace the error passed to UpdateH
 * by zero on one accepted step, or replace the step returned by EstimateStep
 * by NaN on one step.
 * ---------------------------------------------------------------------------*/

#define MAX_RECORDS 64

typedef struct
{
  sunrealtype h, dsm, hnew;
  int p;
  int attempts; /* number of EstimateStep calls for this step (1 = no failure) */
} StepRecord;

typedef struct
{
  SUNAdaptController inner;
  int zero_step;    /* accepted step whose stored error is replaced by 0 */
  int nan_step;     /* step whose EstimateStep output is replaced by NaN */
  int num_accepted; /* accepted steps so far */
  StepRecord last;  /* EstimateStep record for the current step attempt */
  StepRecord records[MAX_RECORDS + 1]; /* records[k] for accepted step k */
} ShimContent;

#define SHIM(C) ((ShimContent*)((C)->content))

static SUNAdaptController_Type shim_gettype(SUNAdaptController C)
{
  (void)C;
  return SUN_ADAPTCONTROLLER_H;
}

static SUNErrCode shim_estimatestep(SUNAdaptController C, sunrealtype h, int p,
                                    sunrealtype dsm, sunrealtype* hnew)
{
  ShimContent* content = SHIM(C);
  SUNErrCode ier = SUNAdaptController_EstimateStep(content->inner, h, p, dsm,
                                                   hnew);
  if (ier) { return ier; }
  if (content->num_accepted + 1 == content->nan_step) { *hnew = NAN; }

  content->last.h    = h;
  content->last.p    = p;
  content->last.dsm  = dsm;
  content->last.hnew = *hnew;
  content->last.attempts++;
  return SUN_SUCCESS;
}

static SUNErrCode shim_updateh(SUNAdaptController C, sunrealtype h,
                               sunrealtype dsm)
{
  ShimContent* content = SHIM(C);
  content->num_accepted++;
  if (content->num_accepted <= MAX_RECORDS)
  {
    content->records[content->num_accepted] = content->last;
  }
  content->last.attempts = 0;

  if (content->num_accepted == content->zero_step) { dsm = SUN_RCONST(0.0); }
  return SUNAdaptController_UpdateH(content->inner, h, dsm);
}

static SUNErrCode shim_reset(SUNAdaptController C)
{
  ShimContent* content   = SHIM(C);
  content->num_accepted  = 0;
  content->last.attempts = 0;
  return SUNAdaptController_Reset(content->inner);
}

static SUNErrCode shim_destroy(SUNAdaptController C)
{
  free(C->content);
  C->content = NULL;
  SUNAdaptController_DestroyEmpty(C);
  return SUN_SUCCESS;
}

/* The shim does not own the inner controller */
static SUNAdaptController create_shim(SUNContext sunctx, SUNAdaptController inner,
                                      int zero_step, int nan_step)
{
  SUNAdaptController C = SUNAdaptController_NewEmpty(sunctx);
  if (C == NULL) { return NULL; }

  ShimContent* content = (ShimContent*)calloc(1, sizeof(ShimContent));
  if (content == NULL)
  {
    SUNAdaptController_DestroyEmpty(C);
    return NULL;
  }
  content->inner     = inner;
  content->zero_step = zero_step;
  content->nan_step  = nan_step;

  C->content           = content;
  C->ops->gettype      = shim_gettype;
  C->ops->estimatestep = shim_estimatestep;
  C->ops->updateh      = shim_updateh;
  C->ops->reset        = shim_reset;
  C->ops->destroy      = shim_destroy;
  return C;
}

/* -----------------------------------------------------------------------------
 * Test problem and error handler for the ARKODE-level checks
 * ---------------------------------------------------------------------------*/

/* RHS for y' = -y: smooth, with nonzero local error estimates on every step */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  (void)t;
  (void)user_data;
  N_VScale(SUN_RCONST(-1.0), y, ydot);
  return 0;
}

/* Print errors unless one is expected; record the last error code */
typedef struct
{
  sunbooleantype expect_error;
  SUNErrCode last_code;
} ErrData;

static void err_fn(int line, const char* func, const char* file, const char* msg,
                   SUNErrCode err_code, void* err_user_data, SUNContext sunctx)
{
  (void)sunctx;
  ErrData* err_data   = (ErrData*)err_user_data;
  err_data->last_code = err_code;
  if (!err_data->expect_error)
  {
    fprintf(stderr, "Error at line %i of %s in %s: %s\n", line, func, file, msg);
  }
}

/* Create an ERKStep integrator for y' = -y, y(0) = 1, using controller C */
static void* create_integrator(SUNContext sunctx, N_Vector y, SUNAdaptController C)
{
  N_VConst(SUN_RCONST(1.0), y);
  void* arkode_mem = ERKStepCreate(f, SUN_RCONST(0.0), y, sunctx);
  if (arkode_mem == NULL) { return NULL; }
  if (ARKodeSStolerances(arkode_mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-10)) ||
      ARKodeSetAdaptController(arkode_mem, C))
  {
    ARKodeFree(&arkode_mem);
    return NULL;
  }
  return arkode_mem;
}

/* -----------------------------------------------------------------------------
 * Part 2: ARKODE-level check of a zero error inserted into the history
 * ---------------------------------------------------------------------------*/

/* accepted step whose stored error is replaced by zero */
#define ZERO_STEP 6

static int test_arkode_zero_error(const ControllerCase* cc, SUNContext sunctx,
                                  N_Vector y)
{
  int fails       = 0;
  const int H     = cc->historysize;
  const int nlast = ZERO_STEP + H + 1; /* last step to check */

  SUNAdaptController inner = cc->create(sunctx);
  if (inner == NULL || SUNAdaptController_SetErrorBias(inner, BIAS))
  {
    fprintf(stderr, "FAIL [%s, ARKODE] could not create controller\n", cc->name);
    SUNAdaptController_Destroy(inner);
    return 1;
  }
  SUNAdaptController C = create_shim(sunctx, inner, ZERO_STEP, -1);
  void* arkode_mem     = (C != NULL) ? create_integrator(sunctx, y, C) : NULL;
  if (arkode_mem == NULL)
  {
    fprintf(stderr, "FAIL [%s, ARKODE] could not create integrator\n", cc->name);
    SUNAdaptController_Destroy(C);
    SUNAdaptController_Destroy(inner);
    return 1;
  }

  /* take single steps until the steps to check have been accepted */
  const ShimContent* content = SHIM(C);
  sunrealtype tret;
  while (content->num_accepted < nlast)
  {
    if (ARKodeEvolve(arkode_mem, SUN_RCONST(10.0), y, &tret, ARK_ONE_STEP) < 0)
    {
      fprintf(stderr, "FAIL [%s, ARKODE] ARKodeEvolve failed at step %i\n",
              cc->name, content->num_accepted + 1);
      fails++;
      break;
    }
  }

  /* Steps up to and including ZERO_STEP see a full positive history; the H
     steps after it use the I law; step ZERO_STEP + H + 1 again has a full
     positive history. */
  for (int k = H + 1; fails == 0 && k <= nlast; k++)
  {
    const StepRecord* r = &content->records[k];
    char what[64];
    if (k > ZERO_STEP && r->attempts != 1)
    {
      fprintf(stderr, "FAIL [%s, ARKODE] step %i took %i attempts, expected 1\n",
              cc->name, k, r->attempts);
      fails++;
    }
    if (k > ZERO_STEP && k <= ZERO_STEP + H)
    {
      snprintf(what, sizeof(what), "ARKODE step %i (I law)", k);
      fails += check_value(cc->name, what, r->hnew,
                           expected_i_law(r->h, r->p, r->dsm));
    }
    else
    {
      const StepRecord* rp  = &content->records[k - 1];
      const StepRecord* rpp = (k >= 2) ? &content->records[k - 2] : rp;
      snprintf(what, sizeof(what), "ARKODE step %i (full history)", k);
      fails += check_value(cc->name, what, r->hnew,
                           expected_full(cc, r->h, r->p, r->dsm, rp->h, rp->dsm,
                                         rpp->h, rpp->dsm));
    }
  }

  ARKodeFree(&arkode_mem);
  SUNAdaptController_Destroy(C);
  SUNAdaptController_Destroy(inner);
  return fails;
}

/* -----------------------------------------------------------------------------
 * Part 3: ARKODE-level check that a NaN step gives ARK_CONTROLLER_ERR
 * ---------------------------------------------------------------------------*/

#define NAN_STEP 4

static int test_arkode_nan(SUNContext sunctx, N_Vector y, ErrData* err_data)
{
  int fails = 0;

  SUNAdaptController inner = SUNAdaptController_PID(sunctx);
  SUNAdaptController C =
    (inner != NULL) ? create_shim(sunctx, inner, -1, NAN_STEP) : NULL;
  void* arkode_mem = (C != NULL) ? create_integrator(sunctx, y, C) : NULL;
  if (arkode_mem == NULL)
  {
    fprintf(stderr, "FAIL [NaN check] could not create integrator\n");
    SUNAdaptController_Destroy(C);
    SUNAdaptController_Destroy(inner);
    return 1;
  }

  const ShimContent* content = SHIM(C);
  sunrealtype tret;
  int flag               = ARK_SUCCESS;
  err_data->expect_error = SUNTRUE;
  while (flag >= 0 && content->num_accepted < NAN_STEP)
  {
    flag = ARKodeEvolve(arkode_mem, SUN_RCONST(10.0), y, &tret, ARK_ONE_STEP);
  }
  err_data->expect_error = SUNFALSE;

  if (flag != ARK_CONTROLLER_ERR)
  {
    fprintf(stderr,
            "FAIL [NaN check] ARKodeEvolve returned %i (%s), expected "
            "ARK_CONTROLLER_ERR\n",
            flag, ARKodeGetReturnFlagName(flag));
    fails++;
  }
  else if (content->num_accepted != NAN_STEP - 1)
  {
    fprintf(stderr,
            "FAIL [NaN check] failure after %i accepted steps, expected %i\n",
            content->num_accepted, NAN_STEP - 1);
    fails++;
  }

  ARKodeFree(&arkode_mem);
  SUNAdaptController_Destroy(C);
  SUNAdaptController_Destroy(inner);
  return fails;
}

/* -----------------------------------------------------------------------------
 * Main program
 * ---------------------------------------------------------------------------*/

static const ControllerCase* find_case(const char* name)
{
  for (int i = 0; i < num_controller_cases; i++)
  {
    if (strcmp(controller_cases[i].name, name) == 0)
    {
      return &controller_cases[i];
    }
  }
  return NULL;
}

int main(void)
{
  int fails = 0;
  ErrData err_data;
  err_data.expect_error = SUNFALSE;
  err_data.last_code    = SUN_SUCCESS;

  SUNContext sunctx;
  if (SUNContext_Create(SUN_COMM_NULL, &sunctx))
  {
    fprintf(stderr, "SUNContext_Create failed\n");
    return 1;
  }
  if (SUNContext_ClearErrHandlers(sunctx) ||
      SUNContext_PushErrHandler(sunctx, err_fn, &err_data))
  {
    fprintf(stderr, "Setting the error handler failed\n");
    SUNContext_Free(&sunctx);
    return 1;
  }

  N_Vector y = N_VNew_Serial(1, sunctx);
  if (y == NULL)
  {
    fprintf(stderr, "N_VNew_Serial failed\n");
    SUNContext_Free(&sunctx);
    return 1;
  }

  /* Part 1: controller-level checks */
  for (int i = 0; i < num_controller_cases; i++)
  {
    fails += test_controller(&controller_cases[i], sunctx);
  }

  /* Part 2: ARKODE-level checks of a zero error inserted into the history */
  for (int i = 0; i < num_arkode_cases; i++)
  {
    const ControllerCase* cc = find_case(arkode_case_names[i]);
    if (cc == NULL)
    {
      fprintf(stderr, "FAIL unknown controller %s\n", arkode_case_names[i]);
      fails++;
      continue;
    }
    fails += test_arkode_zero_error(cc, sunctx, y);
  }

  /* Part 3: ARKODE-level check of a NaN step */
  fails += test_arkode_nan(sunctx, y, &err_data);

  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  if (fails)
  {
    printf("FAILED: %i check(s) failed\n", fails);
    return 1;
  }
  printf("SUCCESS\n");
  return 0;
}
