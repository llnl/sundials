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
 * Unit test for creating a minimal SUNStepper using the multirate
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

static sunrealtype fast_tcur        = ZERO;
static int fast_evolve_status       = 0;
static sunrealtype fast_accum_error = ZERO;
static sunrealtype fast_rtol        = ZERO;

static int fast_evolve(SUNStepper fast_mem, sunrealtype tf, N_Vector y,
                       sunrealtype* tret)
{
  int i               = 0;
  sunrealtype h_fast  = (tf - fast_tcur) / SUN_RCONST(10.0);
  sunrealtype* y_data = N_VGetArrayPointer(y);

  for (i = 0; i < 10; i++) { y_data[0] += (h_fast * -ONE * y_data[0]); }

  fast_tcur = tf;
  *tret     = tf;

  return fast_evolve_status;
}

static SUNErrCode fast_reset(SUNStepper fast_mem, sunrealtype tR, N_Vector yR)
{
  fast_tcur = tR;
  return SUN_SUCCESS;
}

static SUNErrCode fast_set_forcing(SUNStepper fast_mem, sunrealtype tshift,
                                   sunrealtype tscale, N_Vector* forcing,
                                   int nforcing)
{
  return SUN_SUCCESS;
}

static SUNErrCode fast_get_accumulated_error(SUNStepper fast_mem,
                                             sunrealtype* accum_error)
{
  *accum_error = fast_accum_error;
  return SUN_SUCCESS;
}

static SUNErrCode fast_reset_accumulated_error(SUNStepper fast_mem)
{
  fast_accum_error = ZERO;
  return SUN_SUCCESS;
}

static SUNErrCode fast_set_rtol(SUNStepper fast_mem, sunrealtype rtol)
{
  fast_rtol = rtol;
  return SUN_SUCCESS;
}

int main(int argc, char* argv[])
{
  SUNContext sunctx   = NULL;
  N_Vector y          = NULL;
  void* arkode_mem    = NULL;
  SUNStepper fast_mem = NULL;

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

  flag = SUNStepper_Create(sunctx, &fast_mem);
  if (flag) { return 1; }

  flag = SUNStepper_SetEvolveFn(fast_mem, fast_evolve);
  if (flag) { return 1; }

  flag = SUNStepper_SetResetFn(fast_mem, fast_reset);
  if (flag) { return 1; }

  flag = SUNStepper_SetForcingFn(fast_mem, fast_set_forcing);
  if (flag) { return 1; }

  flag = SUNStepper_SetGetAccumulatedErrorFn(fast_mem,
                                             fast_get_accumulated_error);
  if (flag) { return 1; }

  flag = SUNStepper_SetResetAccumulatedErrorFn(fast_mem,
                                               fast_reset_accumulated_error);
  if (flag) { return 1; }

  flag = SUNStepper_SetRTolFn(fast_mem, fast_set_rtol);
  if (flag) { return 1; }

  /* Verify evolve statuses are propagated without translation. */
  fast_evolve_status = 1;
  if (SUNStepper_Evolve(fast_mem, ZERO, y, &tret) != 1) { return 1; }
  fast_evolve_status = -1;
  if (SUNStepper_Evolve(fast_mem, ZERO, y, &tret) != -1) { return 1; }
  fast_evolve_status = 0;

  /* Verify the accumulated-error and tolerance operations. */
  fast_accum_error = SUN_RCONST(0.25);
  sunrealtype accum_error;
  flag = SUNStepper_GetAccumulatedError(fast_mem, &accum_error);
  if (flag || accum_error != fast_accum_error) { return 1; }
  flag = SUNStepper_ResetAccumulatedError(fast_mem);
  if (flag || fast_accum_error != ZERO) { return 1; }
  flag = SUNStepper_SetRTol(fast_mem, SUN_RCONST(1.0e-4));
  if (flag || fast_rtol != SUN_RCONST(1.0e-4)) { return 1; }

  /* Verify the stateless forcing helper. */
  N_Vector forcing[2] = {N_VClone(y), N_VClone(y)};
  if (!forcing[0] || !forcing[1]) { return 1; }
  N_VConst(SUN_RCONST(2.0), forcing[0]);
  N_VConst(SUN_RCONST(3.0), forcing[1]);

  N_VConst(ONE, y);
  flag = SUNStepper_AddForcing(SUN_RCONST(3.0), ONE, SUN_RCONST(2.0),
                               forcing, 2, y);
  if (flag || N_VGetArrayPointer(y)[0] != SUN_RCONST(6.0)) { return 1; }

  N_VConst(ONE, y);
  flag = SUNStepper_AddForcing(ZERO, ZERO, ZERO, forcing, 1, y);
  if (flag || N_VGetArrayPointer(y)[0] != SUN_RCONST(3.0)) { return 1; }

  flag = SUNStepper_AddForcing(ZERO, ZERO, ZERO, NULL, 0, y);
  if (flag || N_VGetArrayPointer(y)[0] != SUN_RCONST(3.0)) { return 1; }

  N_VDestroy(forcing[0]);
  N_VDestroy(forcing[1]);
  N_VConst(ONE, y);

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

  SUNStepper_Destroy(&fast_mem);
  ARKodeFree(&arkode_mem);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  if (!flag) { printf("SUCCESS\n"); }

  return flag;
}

/*---- end of file ----*/
