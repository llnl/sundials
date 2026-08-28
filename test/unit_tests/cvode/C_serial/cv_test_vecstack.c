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
 * Unit test for attaching and retrieving SUNVecStack objects
 * ---------------------------------------------------------------------------*/

#include <stdint.h>
#include <stdio.h>

#include "cvode/cvode.h"
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
  int64_t num_vecs        = 0;
  SUNContext sunctx       = NULL;
  SUNVecStack stack       = NULL;
  SUNVecStack stack_check = NULL;
  N_Vector y              = NULL;
  void* cvode_mem         = NULL;

  /* Create the shared SUNDIALS context used by CVODE, the N_Vector, and the
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

  /* First check the default path. CVODE should not create an internal vector
     stack until CVodeInit has a template vector available. */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (!cvode_mem)
  {
    fprintf(stderr, "CVodeCreate returned NULL\n");
    return 1;
  }

  flag = CVodeGetVecStack(cvode_mem, &stack_check);
  if (check_flag("CVodeGetVecStack", flag)) { return 1; }
  if (stack_check != NULL)
  {
    fprintf(stderr, "Expected NULL stack before CVodeInit\n");
    return 1;
  }

  /* CVodeInit should create an internal vector stack when the user has not
     already attached one. */
  flag = CVodeInit(cvode_mem, f, ZERO, y);
  if (check_flag("CVodeInit", flag)) { return 1; }

  flag = CVodeGetVecStack(cvode_mem, &stack_check);
  if (check_flag("CVodeGetVecStack", flag)) { return 1; }
  if (stack_check == NULL)
  {
    fprintf(stderr, "Expected non-NULL stack after CVodeInit\n");
    return 1;
  }

  /* Free the first CVODE memory block. Since CVODE created the stack above,
     CVodeFree is responsible for destroying it. */
  CVodeFree(&cvode_mem);

  /* Next check the user-supplied path. This new CVODE memory block will receive
     a stack owned by the caller before CVodeInit is called. */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (!cvode_mem)
  {
    fprintf(stderr, "CVodeCreate returned NULL\n");
    return 1;
  }

  /* Seed the user-owned stack with one vector so we can verify after CVodeFree
     that the stack still exists and has not been destroyed by CVODE. */
  flag = SUNVecStack_Create(y, 1, sunctx, &stack);
  if (check_flag("SUNVecStack_Create", flag)) { return 1; }

  /* Attach the user-owned stack. CVODE should store this pointer but leave
     ownership with the caller. */
  flag = CVodeSetVecStack(cvode_mem, stack);
  if (check_flag("CVodeSetVecStack", flag)) { return 1; }

  /* Initializing with a user-provided stack should reuse that stack instead of
     allocating a new internal one. */
  flag = CVodeInit(cvode_mem, f, ZERO, y);
  if (check_flag("CVodeInit", flag)) { return 1; }

  flag = CVodeGetVecStack(cvode_mem, &stack_check);
  if (check_flag("CVodeGetVecStack", flag)) { return 1; }
  if (stack_check != stack)
  {
    fprintf(stderr, "CVodeGetVecStack did not return the user stack\n");
    return 1;
  }

  /* CVodeFree must not destroy a user-owned stack. The caller should still be
     able to query and destroy it after the CVODE memory is gone. */
  CVodeFree(&cvode_mem);

  flag = SUNVecStack_GetNumVecs(stack, &num_vecs);
  if (check_flag("SUNVecStack_GetNumVecs", flag)) { return 1; }
  if (num_vecs != 1)
  {
    fprintf(stderr, "Expected user stack to remain valid after CVodeFree\n");
    return 1;
  }

  flag = SUNVecStack_Destroy(&stack);
  if (check_flag("SUNVecStack_Destroy", flag)) { return 1; }

  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  printf("SUCCESS\n");

  return 0;
}

/*---- end of file ----*/
