/* -----------------------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ SMU
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
 * Unit test for the adaptive step size controller memory and its behavior with
 * fixed and adaptive stepping switches.
 *
 * Solve y' = -y with ERKStep and a PID controller. After taking two adaptive
 * steps, enable fixed stepping, take two fixed steps, and re-enable adaptivity.
 * Check that enabling fixed stepping resets the controller history, fixed
 * steps leave it unchanged, and adaptive steps populate it again. Turning on
 * adaptivity after fixed stepping uses the step size from the last fixed step.
 * Expected errors are computed from ARKODE's local error estimates and weights.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode/arkode_erkstep.h"
#include "nvector/nvector_serial.h"
#include "sunadaptcontroller/sunadaptcontroller_soderlind.h"
#include "sundials/sundials_math.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/* User-supplied RHS function */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/* Private functions to check return values and controller history */
static int check_retval(void* flagvalue, const char* funcname, int opt);

/* Check the controller history against expected values */
static int check_history(SUNAdaptController controller, sunrealtype ep,
                         sunrealtype epp, sunrealtype hp, sunrealtype hpp,
                         int firststeps);

/* Main program */
int main(void)
{
  SUNContext sunctx                           = NULL;
  N_Vector y                                  = NULL;
  N_Vector ele                                = NULL;
  N_Vector ewt                                = NULL;
  void* arkode_mem                            = NULL;
  SUNAdaptController controller               = NULL;
  SUNAdaptControllerContent_Soderlind content = NULL;
  sunrealtype t                               = ZERO;
  sunrealtype hlast[3]                        = {ONE, ONE, ONE};
  sunrealtype elast[3]                        = {ONE, ONE, ONE};
  sunrealtype hfixed                          = SUN_RCONST(0.001);
  int retval;
  int fails = 0;
  int i;

  printf("\nARKODE step size and error controller test:\n");

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return 1; }

  /* Create solution vector and set initial condition */
  y = N_VNew_Serial(1, sunctx);
  if (check_retval((void*)y, "N_VNew_Serial", 0)) { return 1; }
  N_VConst(ONE, y);

  /* Allocate work vectors for the local error estimate and error weights */
  ele = N_VClone(y);
  if (check_retval((void*)ele, "N_VClone", 0)) { return 1; }
  ewt = N_VClone(y);
  if (check_retval((void*)ewt, "N_VClone", 0)) { return 1; }

  /* Create ERKStep memory and set tolerances */
  arkode_mem = ERKStepCreate(f, ZERO, y, sunctx);
  if (check_retval(arkode_mem, "ERKStepCreate", 0)) { return 1; }

  retval = ARKodeSStolerances(arkode_mem, SUN_RCONST(1.0e-5), SUN_RCONST(1.0e-9));
  if (check_retval(&retval, "ARKodeSStolerances", 1)) { return 1; }

  /* Attach a PID controller, which uses two previous error entries */
  controller = SUNAdaptController_PID(sunctx);
  if (check_retval((void*)controller, "SUNAdaptController_PID", 0))
  {
    return 1;
  }

  retval = ARKodeSetAdaptController(arkode_mem, controller);
  if (check_retval(&retval, "ARKodeSetAdaptController", 1)) { return 1; }

  retval = ARKodeSetAccumulatedErrorType(arkode_mem, ARK_ACCUMERROR_NONE);
  if (check_retval(&retval, "ARKodeSetAccumulatedErrorType", 1)) { return 1; }

  printf("\nSetting initial step to " SUN_FORMAT_G "\n", hfixed);
  retval = ARKodeSetInitStep(arkode_mem, hfixed);
  if (check_retval(&retval, "ARKodeSetInitStep", 1)) { return 1; }

  content = (SUNAdaptControllerContent_Soderlind)controller->content;
  if (content->historysize != 2)
  {
    fprintf(stderr, "Expected PID history size 2 but is %i\n",
            content->historysize);
    return 1;
  }

  printf("\nInitial controller history check\n");
  fails += check_history(controller, ONE, ONE, ONE, ONE, 0);

  /* Take adaptive steps and record their step sizes and normalized errors */
  for (i = 0; i < 2; i++)
  {
    retval = ARKodeEvolve(arkode_mem, t + ONE, y, &t, ARK_ONE_STEP);
    if (check_retval(&retval, "ARKodeEvolve", 1)) { return 1; }

    retval = ARKodeGetLastStep(arkode_mem, &hlast[i + 1]);
    if (check_retval(&retval, "ARKodeGetLastStep", 1)) { return 1; }

    retval = ARKodeGetEstLocalErrors(arkode_mem, ele);
    if (check_retval(&retval, "ARKodeGetEstLocalErrors", 1)) { return 1; }
    retval = ARKodeGetErrWeights(arkode_mem, ewt);
    if (check_retval(&retval, "ARKodeGetErrWeights", 1)) { return 1; }
    elast[i + 1] = content->bias * N_VWrmsNorm(ele, ewt);

    printf("\nChecking controller history after adaptive step %i\n", i + 1);
    fails += check_history(controller, elast[i + 1], elast[i], hlast[i + 1],
                           hlast[i], i + 1);
  }

  printf("\nFixed-step set to " SUN_FORMAT_G "\n", hfixed);

  /* Enable fixed stepping and check that the controller is reset */
  retval = ARKodeSetFixedStep(arkode_mem, hfixed);
  if (check_retval(&retval, "ARKodeSetFixedStep", 1)) { return 1; }

  printf("\nChecking controller reset after enabling fixed stepping\n");
  fails += check_history(controller, ONE, ONE, ONE, ONE, 0);

  /* Fixed steps must leave the controller history unchanged */
  for (i = 0; i < 2; i++)
  {
    retval = ARKodeEvolve(arkode_mem, t + ONE, y, &t, ARK_ONE_STEP);
    if (check_retval(&retval, "ARKodeEvolve", 1)) { return 1; }

    printf("\nChecking controller history after fixed step %i\n", i + 1);
    fails += check_history(controller, ONE, ONE, ONE, ONE, 0);
  }

  /* Re-enable adaptivity without any intervening reset calls */
  retval = ARKodeSetFixedStep(arkode_mem, ZERO);
  if (check_retval(&retval, "ARKodeSetFixedStep", 1)) { return 1; }

  printf("\nChecking controller history after re-enabling adaptivity\n");
  fails += check_history(controller, ONE, ONE, ONE, ONE, 0);

  /* Take adaptive steps and record their step sizes and normalized errors */
  for (i = 0; i < 2; i++)
  {
    retval = ARKodeEvolve(arkode_mem, t + ONE, y, &t, ARK_ONE_STEP);
    if (check_retval(&retval, "ARKodeEvolve", 1)) { return 1; }

    retval = ARKodeGetLastStep(arkode_mem, &hlast[i + 1]);
    if (check_retval(&retval, "ARKodeGetLastStep", 1)) { return 1; }

    retval = ARKodeGetEstLocalErrors(arkode_mem, ele);
    if (check_retval(&retval, "ARKodeGetEstLocalErrors", 1)) { return 1; }
    retval = ARKodeGetErrWeights(arkode_mem, ewt);
    if (check_retval(&retval, "ARKodeGetErrWeights", 1)) { return 1; }
    elast[i + 1] = content->bias * N_VWrmsNorm(ele, ewt);

    printf("\nChecking controller history after adaptive step %i\n", i + 1);
    fails += check_history(controller, elast[i + 1], elast[i], hlast[i + 1],
                           hlast[i], i + 1);
  }

  /* Free integrator, controller, vectors, and context */
  ARKodeFree(&arkode_mem);
  SUNAdaptController_Destroy(controller);
  N_VDestroy(ele);
  N_VDestroy(ewt);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  if (fails)
  {
    fprintf(stderr, "FAIL: %i checks failed\n", fails);
    return 1;
  }

  printf("SUCCESS\n");
  return 0;
}

/* RHS for the scalar ODE y' = -y */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  (void)t;
  (void)user_data;
  N_VScale(-ONE, y, ydot);
  return 0;
}

/* Check function return values */
static int check_retval(void* flagvalue, const char* funcname, int opt)
{
  int* errflag;

  /* Check if SUNDIALS function returned NULL pointer */
  if (opt == 0 && flagvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  /* Check if flag < 0 */
  else if (opt == 1)
  {
    errflag = (int*)flagvalue;
    if (*errflag < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *errflag);
      return 1;
    }
  }

  return 0;
}

/* Check the PID implementation's stored errors, step sizes, and history count */
static int check_history(SUNAdaptController controller, sunrealtype ep,
                         sunrealtype epp, sunrealtype hp, sunrealtype hpp,
                         int firststeps)
{
  SUNAdaptControllerContent_Soderlind content =
    (SUNAdaptControllerContent_Soderlind)controller->content;

  printf("  Expected: ep=" SUN_FORMAT_G ", epp=" SUN_FORMAT_G
         ", hp=" SUN_FORMAT_G ", hpp=" SUN_FORMAT_G ", firststeps=%i\n",
         ep, epp, hp, hpp, firststeps);
  printf("  Actual:   ep=" SUN_FORMAT_G ", epp=" SUN_FORMAT_G
         ", hp=" SUN_FORMAT_G ", hpp=" SUN_FORMAT_G ", firststeps=%i\n",
         content->ep, content->epp, content->hp, content->hpp,
         content->firststeps);

  /* Reconstructed error norms may differ by roundoff */
  if (SUNRCompare(content->ep, ep) || SUNRCompare(content->epp, epp) ||
      content->hp != hp || content->hpp != hpp ||
      content->firststeps != firststeps)
  {
    fprintf(stderr, "     FAIL: Unexpected controller history:\n");
    return 1;
  }

  return 0;
}
