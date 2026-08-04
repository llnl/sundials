/* -----------------------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
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
 * Unit tests for PDAEStep nonlinear solver APIs.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>

#include "arkode/arkode.h"
#include "arkode/arkode_pdaestep.h"
#include "arkode/arkode_pdaestep_impl.h"
#include "nvector/nvector_serial.h"
#include "sunnonlinsol/sunnonlinsol_auto.h"
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

#define CHECK_RET(flag, expected, name)                               \
  do {                                                                \
    int _flag = (flag);                                               \
    if (_flag != (expected))                                          \
    {                                                                 \
      fprintf(stderr, "%s returned %i, expected %i\n", (name), _flag, \
              (expected));                                            \
      return 1;                                                       \
    }                                                                 \
  }                                                                   \
  while (0)

#define CHECK_TRUE(cond, name)                \
  do {                                        \
    if (!(cond))                              \
    {                                         \
      fprintf(stderr, "%s failed\n", (name)); \
      return 1;                               \
    }                                         \
  }                                           \
  while (0)

static int component_res(sunrealtype t, N_Vector y, N_Vector w, N_Vector yp,
                         N_Vector res, void* user_data)
{
  N_VConst(ZERO, res);
  return 0;
}

static int algebraic_res(sunrealtype t, N_Vector y, N_Vector w, N_Vector res,
                         void* user_data)
{
  N_VConst(ZERO, res);
  return 0;
}

static int create_problem(SUNContext sunctx, void** arkode_mem, N_Vector* x,
                          N_Vector* z, N_Vector* w, N_Vector* xp, N_Vector* zp,
                          N_Vector* wp, N_Vector* y, N_Vector* yp)
{
  PDAEStepComponentResFn component_res_fns[1] = {component_res};

  *x  = N_VNew_Serial(1, sunctx);
  *z  = N_VNew_Serial(1, sunctx);
  *w  = N_VNew_Serial(1, sunctx);
  *xp = N_VNew_Serial(1, sunctx);
  *zp = N_VNew_Serial(1, sunctx);
  *wp = N_VNew_Serial(1, sunctx);
  if (!*x || !*z || !*w || !*xp || !*zp || !*wp) { return 1; }

  N_VConst(ONE, *x);
  N_VConst(ONE, *z);
  N_VConst(ZERO, *w);
  N_VConst(ZERO, *xp);
  N_VConst(ZERO, *zp);
  N_VConst(ZERO, *wp);

  *y  = PDAEStepManyVector(x, z, *w, 1);
  *yp = PDAEStepManyVector(xp, zp, *wp, 1);
  if (!*y || !*yp) { return 1; }

  *arkode_mem = PDAEStepCreate(component_res_fns, algebraic_res, ZERO, *y, *yp,
                               1, sunctx);
  if (!*arkode_mem) { return 1; }

  return 0;
}

static void free_problem(void** arkode_mem, N_Vector x, N_Vector z, N_Vector w,
                         N_Vector xp, N_Vector zp, N_Vector wp, N_Vector y,
                         N_Vector yp)
{
  ARKodeFree(arkode_mem);
  N_VDestroy(x);
  N_VDestroy(z);
  N_VDestroy(w);
  N_VDestroy(xp);
  N_VDestroy(zp);
  N_VDestroy(wp);
  N_VDestroy(y);
  N_VDestroy(yp);
}

int main(void)
{
  SUNContext sunctx = NULL;
  void* arkode_mem  = NULL;
  N_Vector x = NULL, z = NULL, w = NULL, xp = NULL, zp = NULL, wp = NULL;
  N_Vector y = NULL, yp = NULL, alg_template = NULL;
  SUNNonlinearSolver fp_nls   = NULL;
  SUNNonlinearSolver auto_nls = NULL;
  sunrealtype tret            = ZERO;
  long int nsetups = -1, nniters = -1, nnfails = -1;

  CHECK_RET(SUNContext_Create(SUN_COMM_NULL, &sunctx), 0, "SUNContext_Create");

  CHECK_RET(create_problem(sunctx, &arkode_mem, &x, &z, &w, &xp, &zp, &wp, &y,
                           &yp),
            0, "create_problem");

  CHECK_RET(PDAEStepGetAlgebraicVectorTemplate(arkode_mem, &alg_template),
            ARK_SUCCESS, "PDAEStepGetAlgebraicVectorTemplate");
  CHECK_TRUE(alg_template != NULL, "alg_template != NULL");

  CHECK_RET(PDAEStepSetNonlinearSolver(arkode_mem, NULL), ARK_ILL_INPUT,
            "PDAEStepSetNonlinearSolver(NULL)");

  fp_nls = SUNNonlinSol_FixedPoint(alg_template, 1, sunctx);
  CHECK_TRUE(fp_nls != NULL, "SUNNonlinSol_FixedPoint");
  CHECK_RET(PDAEStepSetNonlinearSolver(arkode_mem, fp_nls), ARK_ILL_INPUT,
            "PDAEStepSetNonlinearSolver(fixed-point)");
  SUNNonlinSolFree(fp_nls);
  fp_nls = NULL;

  auto_nls = SUNNonlinSol_Auto(alg_template, 1, SUNNONLINSOL_AUTO_FIXEDPOINT,
                               sunctx);
  CHECK_TRUE(auto_nls != NULL, "SUNNonlinSol_Auto");
  CHECK_RET(PDAEStepSetNonlinearSolver(arkode_mem, auto_nls), ARK_ILL_INPUT,
            "PDAEStepSetNonlinearSolver(auto)");
  SUNNonlinSolFree(auto_nls);
  auto_nls = NULL;

  CHECK_RET(PDAEStepSetMaxNonlinIters(arkode_mem, 7), ARK_SUCCESS,
            "PDAEStepSetMaxNonlinIters");
  CHECK_TRUE(((ARKodePDAEStepMem)((ARKodeMem)arkode_mem)->step_mem)->maxcor == 7,
             "maxcor update");
  CHECK_RET(PDAEStepSetMaxNonlinIters(arkode_mem, 0), ARK_SUCCESS,
            "PDAEStepSetMaxNonlinIters reset");
  CHECK_TRUE(((ARKodePDAEStepMem)((ARKodeMem)arkode_mem)->step_mem)->maxcor ==
               DEFAULT_MAX_COR,
             "maxcor reset");

  CHECK_RET(PDAEStepSetNonlinConvCoef(arkode_mem, SUN_RCONST(0.25)),
            ARK_SUCCESS, "PDAEStepSetNonlinConvCoef");
  CHECK_TRUE(((ARKodePDAEStepMem)((ARKodeMem)arkode_mem)->step_mem)->nlscoef ==
               SUN_RCONST(0.25),
             "nlscoef update");
  CHECK_RET(PDAEStepSetNonlinConvCoef(arkode_mem, ZERO), ARK_SUCCESS,
            "PDAEStepSetNonlinConvCoef reset");
  CHECK_TRUE(((ARKodePDAEStepMem)((ARKodeMem)arkode_mem)->step_mem)->nlscoef ==
               DEFAULT_NLSCOEF,
             "nlscoef reset");

  CHECK_RET(PDAEStepSetNonlinCRDown(arkode_mem, SUN_RCONST(0.45)), ARK_SUCCESS,
            "PDAEStepSetNonlinCRDown");
  CHECK_TRUE(((ARKodePDAEStepMem)((ARKodeMem)arkode_mem)->step_mem)->crdown ==
               SUN_RCONST(0.45),
             "crdown update");
  CHECK_RET(PDAEStepSetNonlinCRDown(arkode_mem, ZERO), ARK_SUCCESS,
            "PDAEStepSetNonlinCRDown reset");
  CHECK_TRUE(((ARKodePDAEStepMem)((ARKodeMem)arkode_mem)->step_mem)->crdown ==
               DEFAULT_CRDOWN,
             "crdown reset");

  CHECK_RET(PDAEStepSetNonlinRDiv(arkode_mem, SUN_RCONST(3.0)), ARK_SUCCESS,
            "PDAEStepSetNonlinRDiv");
  CHECK_TRUE(((ARKodePDAEStepMem)((ARKodeMem)arkode_mem)->step_mem)->rdiv ==
               SUN_RCONST(3.0),
             "rdiv update");
  CHECK_RET(PDAEStepSetNonlinRDiv(arkode_mem, ZERO), ARK_SUCCESS,
            "PDAEStepSetNonlinRDiv reset");
  CHECK_TRUE(((ARKodePDAEStepMem)((ARKodeMem)arkode_mem)->step_mem)->rdiv ==
               DEFAULT_RDIV,
             "rdiv reset");

  CHECK_RET(PDAEStepGetNumLinSolvSetups(arkode_mem, &nsetups), ARK_SUCCESS,
            "PDAEStepGetNumLinSolvSetups");
  CHECK_RET(PDAEStepGetNonlinSolvStats(arkode_mem, &nniters, &nnfails),
            ARK_SUCCESS, "PDAEStepGetNonlinSolvStats");
  CHECK_TRUE(nsetups == 0 && nniters == 0 && nnfails == 0, "initial stats");

  CHECK_RET(ARKodeSetFixedStep(arkode_mem, SUN_RCONST(0.01)), ARK_SUCCESS,
            "ARKodeSetFixedStep");
  CHECK_RET(ARKodeEvolve(arkode_mem, SUN_RCONST(0.01), y, &tret, ARK_NORMAL),
            ARK_NLS_INIT_FAIL, "ARKodeEvolve missing PDAEStep LS");

  free_problem(&arkode_mem, x, z, w, xp, zp, wp, y, yp);
  SUNContext_Free(&sunctx);

  printf("SUCCESS\n");
  return 0;
}
