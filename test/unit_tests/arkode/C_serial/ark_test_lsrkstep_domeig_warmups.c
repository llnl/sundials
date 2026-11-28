/* -----------------------------------------------------------------------------
 * Programmer(s): Vladimir Ivannikov @ Hereon
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
 * Unit test for the number of preprocessing iterations (warmups) LSRKStep
 * applies to an attached SUNDomEigEstimator: the initial count should be used
 * only for the first estimate, whereas the value sent to
 * LSRKStepSetNumDomEigEstPreprocessIters should be used for all subsequent
 * estimates.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode/arkode_lsrkstep.h"
#include "nvector/nvector_serial.h"
#include "sundomeigest/sundomeigest_power.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

#define INIT_WARMUPS    7
#define STEADY_WARMUPS  3
#define UPDATED_WARMUPS 5

/* y' = -diag(1, 2) y */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* yd  = N_VGetArrayPointer(y);
  sunrealtype* ydd = N_VGetArrayPointer(ydot);
  ydd[0]           = -yd[0];
  ydd[1]           = -SUN_RCONST(2.0) * yd[1];
  return 0;
}

/* Take one step and check the work done by the dominant eigenvalue estimate in
   that step. With a real estimate and max_iters = 1, each estimate performs
   exactly warmups + 1 ATimes calls. */
static int step_and_check(void* arkode_mem, SUNDomEigEstimator DEE, N_Vector y,
                          sunrealtype* t, int expected_warmups, const char* label)
{
  long int num_evals_before = 0, num_evals_after = 0, num_ATimes = 0;

  if (LSRKStepGetNumDomEigUpdates(arkode_mem, &num_evals_before)) { return 1; }

  if (ARKodeEvolve(arkode_mem, ONE, y, t, ARK_ONE_STEP) < 0)
  {
    fprintf(stderr, "%s: ARKodeEvolve failed\n", label);
    return 1;
  }

  if (LSRKStepGetNumDomEigUpdates(arkode_mem, &num_evals_after)) { return 1; }
  if (num_evals_after != num_evals_before + 1)
  {
    fprintf(stderr, "%s: expected one estimate in the step, got %li\n", label,
            num_evals_after - num_evals_before);
    return 1;
  }

  if (SUNDomEigEstimator_GetNumATimesCalls(DEE, &num_ATimes)) { return 1; }
  if (num_ATimes != expected_warmups + 1)
  {
    fprintf(stderr, "%s: expected %i warmups, estimate used %li\n", label,
            expected_warmups, num_ATimes - 1);
    return 1;
  }

  printf("%s: %i warmups (PASSED)\n", label, expected_warmups);
  return 0;
}

int main(int argc, char* argv[])
{
  int fails              = 0;
  SUNContext sunctx      = NULL;
  N_Vector y             = NULL;
  void* arkode_mem       = NULL;
  SUNDomEigEstimator DEE = NULL;
  sunrealtype t          = ZERO;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx))
  {
    fprintf(stderr, "SUNContext_Create failed\n");
    return 1;
  }

  y = N_VNew_Serial(2, sunctx);
  if (!y)
  {
    fprintf(stderr, "N_VNew_Serial returned NULL\n");
    return 1;
  }
  N_VConst(ONE, y);

  arkode_mem = LSRKStepCreateSTS(f, ZERO, y, sunctx);
  if (!arkode_mem)
  {
    fprintf(stderr, "LSRKStepCreateSTS returned NULL\n");
    return 1;
  }

  DEE = SUNDomEigEstimator_Power(y, 1, ZERO, sunctx);
  if (!DEE)
  {
    fprintf(stderr, "SUNDomEigEstimator_Power returned NULL\n");
    return 1;
  }

  if (LSRKStepSetDomEigEstimator(arkode_mem, DEE))
  {
    fprintf(stderr, "LSRKStepSetDomEigEstimator failed\n");
    return 1;
  }

  /* The complex-eigenvalue path may add an ATimes call, so request a real one */
  if (SUNDomEigEstimator_SetIsReal_Power(DEE, SUNTRUE))
  {
    fprintf(stderr, "SUNDomEigEstimator_SetIsReal_Power failed\n");
    return 1;
  }

  if (LSRKStepSetNumDomEigEstInitPreprocessIters(arkode_mem, INIT_WARMUPS))
  {
    fprintf(stderr, "LSRKStepSetNumDomEigEstInitPreprocessIters failed\n");
    return 1;
  }

  if (LSRKStepSetNumDomEigEstPreprocessIters(arkode_mem, STEADY_WARMUPS))
  {
    fprintf(stderr, "LSRKStepSetNumDomEigEstPreprocessIters failed\n");
    return 1;
  }

  if (LSRKStepSetDomEigFrequency(arkode_mem, 1))
  {
    fprintf(stderr, "LSRKStepSetDomEigFrequency failed\n");
    return 1;
  }

  if (ARKodeSetFixedStep(arkode_mem, SUN_RCONST(0.01)))
  {
    fprintf(stderr, "ARKodeSetFixedStep failed\n");
    return 1;
  }

  if (ARKodeInit(arkode_mem))
  {
    fprintf(stderr, "ARKodeInit failed\n");
    return 1;
  }

  fails += step_and_check(arkode_mem, DEE, y, &t, INIT_WARMUPS, "First estimate");
  fails += step_and_check(arkode_mem, DEE, y, &t, STEADY_WARMUPS,
                          "Second estimate");

  /* Change the steady-state count after the first estimate */
  if (LSRKStepSetNumDomEigEstPreprocessIters(arkode_mem, UPDATED_WARMUPS))
  {
    fprintf(stderr, "LSRKStepSetNumDomEigEstPreprocessIters failed\n");
    return 1;
  }
  fails += step_and_check(arkode_mem, DEE, y, &t, UPDATED_WARMUPS,
                          "Estimate after setter");

  ARKodeFree(&arkode_mem);
  SUNDomEigEstimator_Destroy(&DEE);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  if (fails)
  {
    printf("FAIL: %i check(s) failed\n", fails);
    return 1;
  }

  printf("SUCCESS\n");
  return 0;
}
