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

#include <idas/idas.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sundials/sundials_math.h>

static int Residual(sunrealtype t, N_Vector y, N_Vector yp, N_Vector residual,
                    void* user_data)
{
  (void)t;
  (void)user_data;
  sunrealtype x          = NV_Ith_S(y, 0);
  NV_Ith_S(residual, 0) = x * x + NV_Ith_S(yp, 0) - SUN_RCONST(2.0);
  return 0;
}

int main(void)
{
  SUNContext ctx;
  if (SUNContext_Create(SUN_COMM_NULL, &ctx)) { return 1; }

  N_Vector y         = N_VNew_Serial(1, ctx);
  N_Vector yp        = N_VClone(y);
  N_Vector guess     = N_VClone(y);
  SUNMatrix matrix   = SUNDenseMatrix(1, 1, ctx);
  SUNLinearSolver ls = SUNLinSol_Dense(y, matrix, ctx);
  void* mem          = IDACreate(ctx);
  int result         = 1;

  N_VConst(SUN_RCONST(1.0), y);
  N_VConst(SUN_RCONST(0.0), yp);
  if (IDAInit(mem, Residual, SUN_RCONST(0.0), y, yp) ||
      IDASStolerances(mem, SUN_RCONST(1e-4), SUN_RCONST(1e-6)) ||
      IDASetLinearSolver(mem, ls, matrix) || IDASetMaxNumItersIC(mem, 1) ||
      IDASetMaxNumJacsIC(mem, 1))
  {
    goto cleanup;
  }

  int flag = IDACalcICWithOutput(mem, IDA_Y_INIT, SUN_RCONST(1.0), y, yp);
  if (flag >= 0 ||
      NV_Ith_S(y, 0) <= SUN_RCONST(1.0) || NV_Ith_S(yp, 0) != SUN_RCONST(0.0) ||
      IDAGetConsistentIC(mem, guess, NULL) ||
      NV_Ith_S(guess, 0) != SUN_RCONST(1.0))
  {
    fprintf(stderr, "Unexpected failed IC results (flag %d, candidate %g, derivative %g, stored %g)\n",
            flag, (double)NV_Ith_S(y, 0), (double)NV_Ith_S(yp, 0),
            (double)NV_Ith_S(guess, 0));
    goto cleanup;
  }

  if (IDAReInit(mem, SUN_RCONST(0.0), y, yp) ||
      IDASetMaxNumItersIC(mem, 10) || IDASetMaxNumJacsIC(mem, 4) ||
      IDACalcICWithOutput(mem, IDA_Y_INIT, SUN_RCONST(1.0), y, yp) ||
      SUNRabs(NV_Ith_S(y, 0) * NV_Ith_S(y, 0) - SUN_RCONST(2.0)) > SUN_RCONST(1e-4) ||
      IDAGetConsistentIC(mem, guess, NULL) ||
      NV_Ith_S(guess, 0) != NV_Ith_S(y, 0))
  {
    goto cleanup;
  }

  result = 0;

cleanup:
  IDAFree(&mem);
  SUNLinSolFree(ls);
  SUNMatDestroy(matrix);
  N_VDestroy(guess);
  N_VDestroy(yp);
  N_VDestroy(y);
  SUNContext_Free(&ctx);
  return result;
}
