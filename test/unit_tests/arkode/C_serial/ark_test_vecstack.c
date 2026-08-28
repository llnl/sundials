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
 * Unit test for attaching and retrieving SUNVecStack objects
 * ---------------------------------------------------------------------------*/

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>

#include "arkode/arkode_erkstep.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_vecstack.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VConst(ZERO, ydot);
  return 0;
}

static int check_flag(const char* name, int flag)
{
  if (flag)
  {
    fprintf(stderr, "%s returned %i\n", name, flag);
    return 1;
  }

  return 0;
}

int main(int argc, char* argv[])
{
  int flag                = 0;
  int64_t num_active_vecs = 0;
  int64_t num_idle_vecs   = 0;
  int64_t num_vecs        = 0;
  sunrealtype tret        = ZERO;
  SUNContext sunctx       = NULL;
  SUNVecStack stack       = NULL;
  SUNVecStack stack_check = NULL;
  SUNVecStack stack_other = NULL;
  N_Vector ele            = NULL;
  N_Vector tmp            = NULL;
  N_Vector y              = NULL;
  N_Vector y_resize       = NULL;
  void* arkode_mem        = NULL;

  /* Create the shared SUNDIALS context used by ARKODE, the N_Vector, and the
     user-owned vector stack created later in this test. */
  flag = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_flag("SUNContext_Create", flag)) { return 1; }

  /* Create a minimal serial vector to use as the initial condition and as the
     template for any vectors allocated by SUNVecStack. */
  y = N_VNew_Serial(1, sunctx);
  if (!y)
  {
    fprintf(stderr, "N_VNew_Serial returned NULL\n");
    return 1;
  }
  N_VConst(ONE, y);

  /* Check the default path through a public ARKODE stepper constructor. ARKODE
     should create an internal vector stack when the stepper is initialized with
     a template vector. */
  arkode_mem = ERKStepCreate(f, ZERO, y, sunctx);
  if (!arkode_mem)
  {
    fprintf(stderr, "ERKStepCreate returned NULL\n");
    return 1;
  }

  flag = ARKodeGetVecStack(arkode_mem, &stack_check);
  if (check_flag("ARKodeGetVecStack", flag)) { return 1; }
  if (stack_check == NULL)
  {
    fprintf(stderr, "Expected non-NULL stack after ERKStepCreate\n");
    return 1;
  }

  /* Advance the problem once so that the stack is exercised by the ERKStep
     stage computations and interpolation infrastructure. */
  flag = ARKodeSStolerances(arkode_mem, SUN_RCONST(1.0e-4), SUN_RCONST(1.0e-8));
  if (check_flag("ARKodeSStolerances", flag)) { return 1; }

  flag = ARKodeSetMaxNumSteps(arkode_mem, 1000);
  if (check_flag("ARKodeSetMaxNumSteps", flag)) { return 1; }

  flag = ARKodeEvolve(arkode_mem, SUN_RCONST(0.1), y, &tret, ARK_NORMAL);
  if (check_flag("ARKodeEvolve", flag)) { return 1; }

  flag = SUNVecStack_GetNumActiveVecs(stack_check, &num_active_vecs);
  if (check_flag("SUNVecStack_GetNumActiveVecs", flag)) { return 1; }
  if (num_active_vecs <= 0)
  {
    fprintf(stderr, "Expected active workspace vectors after ARKodeEvolve\n");
    return 1;
  }

  /* Query the local error estimate. This checks that ARKODE can keep lte alive
     after a step when it has been checked out from the vector stack. */
  ele = N_VClone(y);
  if (!ele)
  {
    fprintf(stderr, "N_VClone returned NULL\n");
    return 1;
  }

  flag = ARKodeGetEstLocalErrors(arkode_mem, ele);
  if (check_flag("ARKodeGetEstLocalErrors", flag)) { return 1; }

  /* Resize the problem while using ARKODE's internal stack. This should refresh
     the stack template so any future temporary vectors match the new problem
     size. */
  y_resize = N_VNew_Serial(2, sunctx);
  if (!y_resize)
  {
    fprintf(stderr, "N_VNew_Serial returned NULL\n");
    return 1;
  }
  N_VConst(ONE, y_resize);

  flag = ARKodeResize(arkode_mem, y_resize, ONE, tret, NULL, NULL);
  if (check_flag("ARKodeResize", flag)) { return 1; }

  flag = ARKodeGetVecStack(arkode_mem, &stack_check);
  if (check_flag("ARKodeGetVecStack", flag)) { return 1; }

  flag = SUNVecStack_Pop(stack_check, &tmp);
  if (check_flag("SUNVecStack_Pop", flag)) { return 1; }
  if (N_VGetLength(tmp) != 2)
  {
    fprintf(stderr, "Expected resized stack vector length 2\n");
    return 1;
  }
  flag = SUNVecStack_Push(stack_check, &tmp);
  if (check_flag("SUNVecStack_Push", flag)) { return 1; }

  flag = ARKodeEvolve(arkode_mem, SUN_RCONST(0.2), y_resize, &tret, ARK_NORMAL);
  if (check_flag("ARKodeEvolve", flag)) { return 1; }

  /* Free the ARKODE memory block. Since ARKODE created the stack above,
     ARKodeFree is responsible for destroying it. */
  ARKodeFree(&arkode_mem);

  /* Next check the user-supplied path. ARKODE public stepper constructors
     create an internal stack, but no vectors are checked out from it until
     ARKodeInit/ARKodeEvolve performs the stepper setup. */
  arkode_mem = ERKStepCreate(f, ZERO, y, sunctx);
  if (!arkode_mem)
  {
    fprintf(stderr, "ERKStepCreate returned NULL\n");
    return 1;
  }

  /* Seed the user-owned stack with one vector so we can verify after ARKodeFree
     that the stack still exists and has not been destroyed by ARKODE. */
  flag = SUNVecStack_Create(y, 1, sunctx, &stack);
  if (check_flag("SUNVecStack_Create", flag)) { return 1; }

  /* Attach the user-owned stack. ARKODE should store this pointer but leave
     ownership with the caller. */
  flag = ARKodeSetVecStack(arkode_mem, stack);
  if (check_flag("ARKodeSetVecStack", flag)) { return 1; }

  flag = ARKodeGetVecStack(arkode_mem, &stack_check);
  if (check_flag("ARKodeGetVecStack", flag)) { return 1; }
  if (stack_check != stack)
  {
    fprintf(stderr, "ARKodeGetVecStack did not return the user stack\n");
    return 1;
  }

  /* Initializing with a user-provided stack should reuse that stack instead of
     the initial internal one. */
  flag = ARKodeSStolerances(arkode_mem, SUN_RCONST(1.0e-4), SUN_RCONST(1.0e-8));
  if (check_flag("ARKodeSStolerances", flag)) { return 1; }

  flag = ARKodeSetMaxNumSteps(arkode_mem, 1000);
  if (check_flag("ARKodeSetMaxNumSteps", flag)) { return 1; }

  flag = ARKodeEvolve(arkode_mem, SUN_RCONST(0.1), y, &tret, ARK_NORMAL);
  if (check_flag("ARKodeEvolve", flag)) { return 1; }

  flag = SUNVecStack_GetNumActiveVecs(stack, &num_active_vecs);
  if (check_flag("SUNVecStack_GetNumActiveVecs", flag)) { return 1; }
  if (num_active_vecs <= 0)
  {
    fprintf(stderr, "Expected active workspace vectors after ARKodeEvolve\n");
    return 1;
  }

  /* Resizing with a user-owned stack is rejected until SUNVecStack supports
     caller-managed resizing of idle cached vectors. */
  flag = ARKodeResize(arkode_mem, y_resize, ONE, tret, NULL, NULL);
  if (flag != ARK_ILL_INPUT)
  {
    fprintf(stderr, "Expected ARKodeResize to reject a user-owned stack\n");
    return 1;
  }

  /* Attempting to replace the stack while vectors are checked out should fail.
     This keeps ARKODE from losing track of vectors it must return on free. */
  flag = SUNVecStack_Create(y, 0, sunctx, &stack_other);
  if (check_flag("SUNVecStack_Create", flag)) { return 1; }

  flag = ARKodeSetVecStack(arkode_mem, stack_other);
  if (flag != ARK_ILL_INPUT)
  {
    fprintf(stderr, "Expected ARKodeSetVecStack to reject an active stack\n");
    return 1;
  }

  flag = SUNVecStack_Destroy(&stack_other);
  if (check_flag("SUNVecStack_Destroy", flag)) { return 1; }

  /* Free the ARKODE memory block. Since the caller owns the stack, it should
     remain usable after ARKodeFree. */
  ARKodeFree(&arkode_mem);

  flag = SUNVecStack_GetNumVecs(stack, &num_vecs);
  if (check_flag("SUNVecStack_GetNumVecs", flag)) { return 1; }
  if (num_vecs <= 1)
  {
    fprintf(stderr,
            "Expected workspace vectors in the user stack after ARKodeFree\n");
    return 1;
  }

  flag = SUNVecStack_GetNumActiveVecs(stack, &num_active_vecs);
  if (check_flag("SUNVecStack_GetNumActiveVecs", flag)) { return 1; }
  if (num_active_vecs != 0)
  {
    fprintf(stderr,
            "Expected no active vectors after ARKodeFree, got %" PRId64 "\n",
            num_active_vecs);
    return 1;
  }

  flag = SUNVecStack_GetNumIdleVecs(stack, &num_idle_vecs);
  if (check_flag("SUNVecStack_GetNumIdleVecs", flag)) { return 1; }
  if (num_idle_vecs != num_vecs)
  {
    fprintf(stderr,
            "Expected all user stack vectors to be idle after ARKodeFree\n");
    return 1;
  }

  flag = SUNVecStack_Destroy(&stack);
  if (check_flag("SUNVecStack_Destroy", flag)) { return 1; }

  N_VDestroy(ele);
  N_VDestroy(y);
  N_VDestroy(y_resize);
  SUNContext_Free(&sunctx);

  return 0;
}
