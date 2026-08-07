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
 * Unit test for step size growth and handling of inf/nan within a controller
 * due to 0 local error.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode/arkode_erkstep.h"
#include "nvector/nvector_serial.h"

/* If an error occurs, print to stderr and exit */
static void err_fn(const int line, const char* const func, const char* const file,
                   const char* const msg, const SUNErrCode err_code,
                   void* const err_user_data, const SUNContext sunctx)
{
  fprintf(stderr, "Error at line %i of %s in %s: %s\n", line, func, file, msg);
  exit(err_code);
}

/* RHS for the simple ODE y'=0 */
static int f(const sunrealtype t, const N_Vector y, const N_Vector ydot,
             void* const user_data)
{
  N_VConst(SUN_RCONST(0.0), ydot);
  return 0;
}

/* Take a single time step solving y'=0 and check the max growth was attained */
static int check_step(void* const arkode_mem, const N_Vector y,
                      const sunrealtype h_expected, const int step)
{
  /* Integration much farther than expected step */
  const sunrealtype tout = SUN_RCONST(100.0) * h_expected;
  sunrealtype tret;
  /* The ERK method should be able to take the maximum possible timestep without
     any rejected steps or local error */
  ARKodeEvolve(arkode_mem, tout, y, &tret, ARK_ONE_STEP);

  long int local_err_fails;
  ARKodeGetNumErrTestFails(arkode_mem, &local_err_fails);
  if (local_err_fails != 0)
  {
    fprintf(stderr, "Expected 0 local error failures at step %i but is %li\n",
            step, local_err_fails);
    return 1;
  }

  const N_Vector err = N_VClone(y);
  ARKodeGetEstLocalErrors(arkode_mem, err);
  const sunrealtype err_norm = N_VMaxNorm(err);
  N_VDestroy(err);

  if (err_norm != SUN_RCONST(0.0))
  {
    fprintf(stderr,
            "Expected local error at step %i to be 0 but is " SUN_FORMAT_G "\n",
            step, err_norm);
    return 1;
  }

  sunrealtype h_actual;
  ARKodeGetCurrentStep(arkode_mem, &h_actual);
  if (SUNRCompare(h_expected, h_actual))
  {
    fprintf(stderr,
            "Expected h at step %i to be " SUN_FORMAT_G " but is " SUN_FORMAT_G
            "\n",
            step, h_expected, h_actual);
    return 1;
  }

  return 0;
}

/* Take several steps solving y'=0 and check the max growth was attained */
static int check_steps(void* const arkode_mem, const N_Vector y,
                       const sunrealtype h0, const sunrealtype first_growth,
                       const sunrealtype growth)
{
  ARKodeSetInitStep(arkode_mem, h0);
  ARKodeSetMaxFirstGrowth(arkode_mem, first_growth);
  ARKodeSetMaxGrowth(arkode_mem, growth);

  sunrealtype h_expect = first_growth * h0;
  /* Take enough steps to allow the controller history to fill up */
  const int num_steps = 5;
  for (int step = 1; step <= num_steps; step++)
  {
    const int retval = check_step(arkode_mem, y, h_expect, step);
    if (retval != 0) { return retval; }
    h_expect *= growth;
  }

  return 0;
}

int main(void)
{
  SUNContext sunctx;
  int retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (retval != 0)
  {
    fprintf(stderr, "SUNContext_Create returned %i\n", retval);
    return 1;
  }

  retval = SUNContext_PushErrHandler(sunctx, err_fn, NULL);
  if (retval != 0)
  {
    fprintf(stderr, "SUNContext_PushErrHandler returned %i\n", retval);
    return 1;
  }

  const N_Vector y = N_VNew_Serial(1, sunctx);
  N_VConst(SUN_RCONST(1.0), y);

  /* Forward integration from 0 */
  void* arkode_mem = ERKStepCreate(f, SUN_RCONST(0.0), y, sunctx);
  retval = check_steps(arkode_mem, y, SUN_RCONST(1.0e-4), SUN_RCONST(1234.0),
                       SUN_RCONST(3.0));
  if (retval != 0) { return retval; }

  /* Backward integration from positive time */
  ERKStepReInit(arkode_mem, f, SUN_RCONST(999.0), y);
  retval = check_steps(arkode_mem, y, SUN_RCONST(-1.0e-2), SUN_RCONST(1.6),
                       SUN_RCONST(2.3));
  if (retval != 0) { return retval; }

  /* Forward integration from a negative time */
  ERKStepReInit(arkode_mem, f, SUN_RCONST(-999.0), y);
  retval = check_steps(arkode_mem, y, SUN_RCONST(20.0), SUN_RCONST(1.0e5),
                       SUN_RCONST(1.1e3));
  if (retval != 0) { return retval; }

  /* Try a non-default controller */
  ERKStepReInit(arkode_mem, f, SUN_RCONST(0.0), y);
  SUNAdaptController controller1 = SUNAdaptController_ExpGus(sunctx);
  ARKodeSetAdaptController(arkode_mem, controller1);
  retval = check_steps(arkode_mem, y, SUN_RCONST(0.1), SUN_RCONST(4.0),
                       SUN_RCONST(10.0));
  if (retval != 0) { return retval; }

  /* Try another non-default controller */
  ERKStepReInit(arkode_mem, f, SUN_RCONST(0.0), y);
  SUNAdaptController controller2 = SUNAdaptController_Soderlind(sunctx);
  SUNAdaptController_SetParams_Soderlind(controller2, SUN_RCONST(0.123),
                                         SUN_RCONST(-0.456), SUN_RCONST(0.789),
                                         SUN_RCONST(-1.0), SUN_RCONST(-2.0));
  ARKodeSetAdaptController(arkode_mem, controller2);
  retval = check_steps(arkode_mem, y, SUN_RCONST(-0.1), SUN_RCONST(4.0),
                       SUN_RCONST(10.0));
  if (retval != 0) { return retval; }

  SUNAdaptController_Destroy(controller1);
  SUNAdaptController_Destroy(controller2);
  ARKodeFree(&arkode_mem);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  printf("SUCCESS\n");

  return 0;
}
