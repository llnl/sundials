/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * Unit test for creating a minimal MRIStepInnerStepper using the multirate
 * Dahlquist problem problem y' = lambda_s y + lambda_f y and a custom explicit
 * Euler inner stepper.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode/arkode.h"
#include "arkode/arkode_mristep.h"
#include "nvector/nvector_serial.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

static int ode_slow_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* y_data    = N_VGetArrayPointer(ydot);
  sunrealtype* ydot_data = N_VGetArrayPointer(ydot);
  ydot_data[0]           = -ONE * y_data[0];
  return 0;
}

static int fast_evolve(MRIStepInnerStepper fast_mem, sunrealtype t0,
                       sunrealtype tf, N_Vector y)
{
  int i               = 0;
  sunrealtype h_fast  = (t0 - tf) / SUN_RCONST(10.0);
  sunrealtype* y_data = N_VGetArrayPointer(y);

  for (i = 0; i < 10; i++) { y_data[0] += (h_fast * -ONE * y_data[0]); }

  return 0;
}

int main(int argc, char* argv[])
{
  SUNContext sunctx            = NULL;
  N_Vector y                   = NULL;
  void* arkode_mem             = NULL;
  MRIStepInnerStepper fast_mem = NULL;

  int flag         = 0;
  int arkode_flag  = 0;
  sunrealtype tout = SUN_RCONST(0.10);
  sunrealtype tret = ZERO;

  /* --------------
   * Create context
   * -------------- */

  flag = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (flag) { return 1; }

  /* -----------------------
   * Setup initial condition
   * ----------------------- */

  y = N_VNew_Serial(1, sunctx);
  if (!y) { return 1; }
  N_VConst(ONE, y);

  /* ---------------------
   * Setup fast integrator
   * --------------------- */

  flag = MRIStepInnerStepper_Create(sunctx, &fast_mem);
  if (flag) { return 1; }

  flag = MRIStepInnerStepper_SetEvolveFn(fast_mem, fast_evolve);
  if (flag) { return 1; }

  /* ---------------------
   * Setup slow integrator
   * --------------------- */

  arkode_mem = MRIStepCreate(ode_slow_rhs, NULL, ZERO, y, fast_mem, sunctx);
  if (!arkode_mem) { return 1; }

  flag = ARKodeSetFixedStep(arkode_mem, SUN_RCONST(0.01));
  if (flag) { return 1; }

  flag = ARKodeSetInterpolantType(arkode_mem, ARK_INTERP_HERMITE);
  if (flag) { return 1; }

  /* ---------------
   * Advance in time
   * --------------- */

  /* Evolve should return a failure when using Hermite interpolation */
  arkode_flag = ARKodeEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);
  printf("ARKodeEvolve returned %i\n", arkode_flag);
  if (arkode_flag != ARK_RHSFUNC_FAIL) { return 1; }

  /* -----------------------
   * Reinitialize integrator
   * ----------------------- */

  N_VConst(ONE, y);

  flag = MRIStepReInit(arkode_mem, ode_slow_rhs, NULL, ZERO, y);
  if (flag) { return 1; }

  flag = ARKodeSetInterpolantType(arkode_mem, ARK_INTERP_LAGRANGE);
  if (flag) { return 1; }

  /* ---------------
   * Advance in time
   * --------------- */

  /* Evolve should succeed when using Lagrange interpolation */
  arkode_flag = ARKodeEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);
  printf("ARKodeEvolve returned %i\n", arkode_flag);
  if (arkode_flag != ARK_SUCCESS) { return 1; }

  /* --------
   * Clean up
   * -------- */

  MRIStepInnerStepper_Free(&fast_mem);
  ARKodeFree(&arkode_mem);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  if (!flag) { printf("SUCCESS\n"); }

  return flag;
}

/*---- end of file ----*/
