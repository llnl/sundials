/* -----------------------------------------------------------------
 * Originating author: Ajrat Makhmutov, 2026-11-18
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------- */

#include <kinsol/kinsol.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sundials/sundials_math.h>

typedef struct
{
  int strategy;
  int action;
  long int calls;
  sunrealtype observed;
  int invalid;
} TestData;

static int System(N_Vector u, N_Vector value, void* user_data)
{
  TestData* data = user_data;
  sunrealtype x  = NV_Ith_S(u, 0);

  NV_Ith_S(value, 0) = data->strategy == KIN_FP
                        ? (x + SUN_RCONST(2.0)) / SUN_RCONST(2.0)
                        : x * x - SUN_RCONST(2.0);
  return 0;
}

static int PicardLinearJacobian(N_Vector u, N_Vector value, SUNMatrix matrix,
                                void* user_data, N_Vector tmp1, N_Vector tmp2)
{
  (void)u;
  (void)value;
  (void)user_data;
  (void)tmp1;
  (void)tmp2;
  SM_ELEMENT_D(matrix, 0, 0) = SUN_RCONST(3.0);
  return 0;
}

static int Observe(long int iteration, N_Vector u, N_Vector value, void* user_data)
{
  TestData* data = user_data;
  sunrealtype x  = NV_Ith_S(u, 0);

  data->calls++;
  data->observed = x;
  if (iteration != data->calls) { data->invalid = 1; }

  if (value)
  {
    sunrealtype expected = data->strategy == KIN_FP
                             ? (x + SUN_RCONST(2.0)) / SUN_RCONST(2.0)
                             : x * x - SUN_RCONST(2.0);
    if (SUNRabs(NV_Ith_S(value, 0) - expected) > SUN_RCONST(20.0) * SUN_UNIT_ROUNDOFF)
    {
      data->invalid = 1;
    }
  }
  else if (data->strategy != KIN_FP) { data->invalid = 1; }

  return data->action;
}

static int Run(SUNContext ctx, int strategy)
{
  TestData data      = {strategy, 1, 0, SUN_RCONST(0.0), 0};
  void* mem          = KINCreate(ctx);
  N_Vector u         = N_VNew_Serial(1, ctx);
  N_Vector scale     = N_VClone(u);
  SUNMatrix matrix   = SUNDenseMatrix(1, 1, ctx);
  SUNLinearSolver ls = SUNLinSol_Dense(u, matrix, ctx);
  int status         = 1;

  N_VConst(SUN_RCONST(1.0), scale);
  if (KINInit(mem, System, u) || KINSetUserData(mem, &data) ||
      KINSetIterationFn(mem, Observe) ||
      (strategy != KIN_FP && KINSetLinearSolver(mem, ls, matrix)) ||
      (strategy == KIN_PICARD && KINSetJacFn(mem, PicardLinearJacobian)))
  {
    goto cleanup;
  }

  NV_Ith_S(u, 0) = strategy == KIN_LINESEARCH ? SUN_RCONST(0.1) : SUN_RCONST(1.0);
  if (KINSol(mem, u, strategy, scale, scale) != KIN_USER_STOP ||
      data.calls != 1 || data.invalid || data.observed != NV_Ith_S(u, 0))
  {
    goto cleanup;
  }

  data.action   = -1;
  data.calls    = 0;
  NV_Ith_S(u, 0) = SUN_RCONST(1.0);
  if (KINSol(mem, u, strategy, scale, scale) != KIN_ITERATION_FN_FAIL ||
      data.calls != 1 || data.invalid || data.observed != NV_Ith_S(u, 0))
  {
    goto cleanup;
  }

  if (KINSetIterationFn(mem, NULL)) { goto cleanup; }
  NV_Ith_S(u, 0) = SUN_RCONST(1.0);
  if (KINSol(mem, u, strategy, scale, scale) < 0 || data.calls != 1)
  {
    goto cleanup;
  }

  status = 0;

cleanup:
  KINFree(&mem);
  SUNLinSolFree(ls);
  SUNMatDestroy(matrix);
  N_VDestroy(scale);
  N_VDestroy(u);
  if (status) { fprintf(stderr, "Iteration callback failed for strategy %d\n", strategy); }
  return status;
}

int main(void)
{
  SUNContext ctx;
  int status = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (status) { return 1; }

  const int strategies[] = {KIN_NONE, KIN_LINESEARCH, KIN_PICARD, KIN_FP};
  for (size_t i = 0; i < sizeof(strategies) / sizeof(*strategies); ++i)
  {
    status |= Run(ctx, strategies[i]);
  }

  SUNContext_Free(&ctx);
  return status;
}
