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
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

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
  SUNLinearSolver LS      = NULL;
  SUNMatrix A             = NULL;
  SUNVecStack stack       = NULL;
  SUNVecStack stack_check = NULL;
  N_Vector ele            = NULL;
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

  /* The checked-out CVODE workspace vectors should now be active in the
     user-owned stack. This includes acor so local-error queries remain valid
     while CVODE is alive. */
  flag = SUNVecStack_GetNumActiveVecs(stack, &num_active_vecs);
  if (check_flag("SUNVecStack_GetNumActiveVecs", flag)) { return 1; }
  if (num_active_vecs != 6)
  {
    fprintf(stderr, "Expected six active workspace vectors after CVodeInit\n");
    return 1;
  }

  /* Advance the problem once, then query the estimated local error. This checks
     that cv_acor still has the required post-CVode lifetime even though its
     storage came from SUNVecStack. */
  flag = CVodeSStolerances(cvode_mem, SUN_RCONST(1.0e-4), SUN_RCONST(1.0e-8));
  if (check_flag("CVodeSStolerances", flag)) { return 1; }

  A = SUNDenseMatrix(1, 1, sunctx);
  if (!A)
  {
    fprintf(stderr, "SUNDenseMatrix returned NULL\n");
    return 1;
  }

  LS = SUNLinSol_Dense(y, A, sunctx);
  if (!LS)
  {
    fprintf(stderr, "SUNLinSol_Dense returned NULL\n");
    return 1;
  }

  flag = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_flag("CVodeSetLinearSolver", flag)) { return 1; }

  ele = N_VClone(y);
  if (!ele)
  {
    fprintf(stderr, "N_VClone returned NULL\n");
    return 1;
  }

  flag = CVode(cvode_mem, SUN_RCONST(0.1), y, &tret, CV_NORMAL);
  if (check_flag("CVode", flag)) { return 1; }

  flag = CVodeGetEstLocalErrors(cvode_mem, ele);
  if (check_flag("CVodeGetEstLocalErrors", flag)) { return 1; }

  /* CVodeFree must not destroy a user-owned stack. The caller should still be
     able to query and destroy it after the CVODE memory is gone. */
  CVodeFree(&cvode_mem);

  flag = SUNVecStack_GetNumVecs(stack, &num_vecs);
  if (check_flag("SUNVecStack_GetNumVecs", flag)) { return 1; }
  if (num_vecs != 6)
  {
    fprintf(stderr, "Expected user stack to remain valid after CVodeFree\n");
    return 1;
  }

  flag = SUNVecStack_GetNumActiveVecs(stack, &num_active_vecs);
  if (check_flag("SUNVecStack_GetNumActiveVecs", flag)) { return 1; }
  if (num_active_vecs != 0)
  {
    fprintf(stderr, "Expected all user stack vectors to be returned\n");
    return 1;
  }

  flag = SUNVecStack_GetNumIdleVecs(stack, &num_idle_vecs);
  if (check_flag("SUNVecStack_GetNumIdleVecs", flag)) { return 1; }
  if (num_idle_vecs != num_vecs)
  {
    fprintf(stderr, "Expected all user stack vectors to be idle\n");
    return 1;
  }

  flag = SUNVecStack_Destroy(&stack);
  if (check_flag("SUNVecStack_Destroy", flag)) { return 1; }

  N_VDestroy(ele);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  printf("SUCCESS\n");

  return 0;
}

/*---- end of file ----*/
