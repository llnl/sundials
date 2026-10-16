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
 * Unit test checking that LSRKStep recomputes the dominant eigenvalue on the
 * first step after ARKodeResize, independent of the update frequency set with
 * LSRKStepSetDomEigFrequency (including a constant dominant eigenvalue).
 *
 * The problem is y' = lambda(N) y with lambda(N) = -100 N, where N is the
 * problem size, so growing the problem increases the stiffness and the number
 * of stages required for stability.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode/arkode_lsrkstep.h"
#include "nvector/nvector_serial.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

#define N_INITIAL 2
#define N_RESIZED 20

static sunrealtype lambda(N_Vector y)
{
  return SUN_RCONST(-100.0) * (sunrealtype)N_VGetLength(y);
}

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VScale(lambda(y), y, ydot);
  return 0;
}

static int dom_eig(sunrealtype t, N_Vector y, N_Vector fn, sunrealtype* lambdaR,
                   sunrealtype* lambdaI, void* user_data, N_Vector temp1,
                   N_Vector temp2, N_Vector temp3)
{
  *lambdaR = lambda(y);
  *lambdaI = ZERO;
  return 0;
}

/* Returns the number of failed checks for the given update frequency */
static int run_test(long int dom_eig_freq, SUNContext sunctx)
{
  int fails                   = 0;
  void* arkode_mem            = NULL;
  N_Vector y                  = NULL;
  N_Vector y_big              = NULL;
  sunrealtype t               = ZERO;
  long int num_updates_before = 0, num_updates_after = 0;
  int stages_before = 0, stages_after = 0;

  y     = N_VNew_Serial(N_INITIAL, sunctx);
  y_big = N_VNew_Serial(N_RESIZED, sunctx);
  if (!y || !y_big)
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

  if (LSRKStepSetDomEigFn(arkode_mem, dom_eig) ||
      LSRKStepSetDomEigFrequency(arkode_mem, dom_eig_freq) ||
      ARKodeSetFixedStep(arkode_mem, SUN_RCONST(0.01)))
  {
    fprintf(stderr, "Setting LSRKStep options failed\n");
    return 1;
  }

  /* With this frequency, only the first step computes an estimate */
  for (int i = 0; i < 3; i++)
  {
    if (ARKodeEvolve(arkode_mem, ONE, y, &t, ARK_ONE_STEP) < 0)
    {
      fprintf(stderr, "ARKodeEvolve failed before resize\n");
      return 1;
    }
  }

  if (LSRKStepGetNumDomEigUpdates(arkode_mem, &num_updates_before) ||
      LSRKStepGetMaxNumStages(arkode_mem, &stages_before))
  {
    return 1;
  }
  if (num_updates_before != 1)
  {
    fprintf(stderr, "freq %li: expected 1 update before resize, got %li\n",
            dom_eig_freq, num_updates_before);
    fails++;
  }

  N_VConst(N_VGetArrayPointer(y)[0], y_big);
  if (ARKodeResize(arkode_mem, y_big, ONE, t, NULL, NULL))
  {
    fprintf(stderr, "ARKodeResize failed\n");
    return 1;
  }

  if (ARKodeEvolve(arkode_mem, ONE, y_big, &t, ARK_ONE_STEP) < 0)
  {
    fprintf(stderr, "ARKodeEvolve failed after resize\n");
    return 1;
  }

  if (LSRKStepGetNumDomEigUpdates(arkode_mem, &num_updates_after) ||
      LSRKStepGetMaxNumStages(arkode_mem, &stages_after))
  {
    return 1;
  }

  if (num_updates_after != num_updates_before + 1)
  {
    fprintf(stderr, "freq %li: expected a new estimate after resize, got %li update(s)\n",
            dom_eig_freq, num_updates_after - num_updates_before);
    fails++;
  }

  /* The stiffer resized problem must be stepped with more stages */
  if (stages_after <= stages_before)
  {
    fprintf(stderr,
            "freq %li: expected more than %i stages after resize, got %i\n",
            dom_eig_freq, stages_before, stages_after);
    fails++;
  }

  if (!fails)
  {
    printf("freq %li: estimate refreshed after resize, stages %i -> %i "
           "(PASSED)\n",
           dom_eig_freq, stages_before, stages_after);
  }

  ARKodeFree(&arkode_mem);
  N_VDestroy(y);
  N_VDestroy(y_big);

  return fails;
}

int main(int argc, char* argv[])
{
  int fails         = 0;
  SUNContext sunctx = NULL;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx))
  {
    fprintf(stderr, "SUNContext_Create failed\n");
    return 1;
  }

  /* Infrequent updates and a constant dominant eigenvalue */
  fails += run_test(1000, sunctx);
  fails += run_test(0, sunctx);

  SUNContext_Free(&sunctx);

  if (fails)
  {
    printf("FAIL: %i check(s) failed\n", fails);
    return 1;
  }

  printf("SUCCESS\n");
  return 0;
}
