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

#include <stdint.h>
#include <stdio.h>

#include "kinsol/kinsol.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_vecstack.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

static int sysfn(N_Vector u, N_Vector fval, void* user_data)
{
  sunrealtype* u_data    = N_VGetArrayPointer(u);
  sunrealtype* fval_data = N_VGetArrayPointer(fval);

  fval_data[0] = u_data[0] - ONE;

  return 0;
}

static int jacfn(N_Vector u, N_Vector fu, SUNMatrix J, void* user_data,
                 N_Vector tmp1, N_Vector tmp2)
{
  sunrealtype* J_data = SUNDenseMatrix_Data(J);

  J_data[0] = ONE;

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

static int setup_linear_solver(void* kin_mem, N_Vector u, SUNContext sunctx,
                               SUNMatrix* A, SUNLinearSolver* LS)
{
  *A = SUNDenseMatrix(1, 1, sunctx);
  if (!(*A))
  {
    fprintf(stderr, "SUNDenseMatrix returned NULL\n");
    return 1;
  }

  *LS = SUNLinSol_Dense(u, *A, sunctx);
  if (!(*LS))
  {
    fprintf(stderr, "SUNLinSol_Dense returned NULL\n");
    return 1;
  }

  if (check_flag("KINSetLinearSolver", KINSetLinearSolver(kin_mem, *LS, *A)))
  {
    return 1;
  }

  if (check_flag("KINSetJacFn", KINSetJacFn(kin_mem, jacfn))) { return 1; }

  return 0;
}

int main(int argc, char* argv[])
{
  int flag                = 0;
  int64_t num_active_vecs = 0;
  int64_t num_idle_vecs   = 0;
  int64_t num_vecs        = 0;
  SUNContext sunctx       = NULL;
  SUNLinearSolver LS      = NULL;
  SUNMatrix A             = NULL;
  SUNVecStack stack       = NULL;
  SUNVecStack stack_check = NULL;
  N_Vector fscale         = NULL;
  N_Vector u              = NULL;
  N_Vector uscale         = NULL;
  void* kin_mem           = NULL;

  /* Create the shared SUNDIALS context used by KINSOL, the N_Vectors, and the
     user-owned vector stack created later in this test. */
  flag = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_flag("SUNContext_Create", flag)) { return 1; }

  /* Create minimal serial vectors for the nonlinear solve.  The u vector also
     acts as the template for vectors allocated by SUNVecStack. */
  u = N_VNew_Serial(1, sunctx);
  if (!u)
  {
    fprintf(stderr, "N_VNew_Serial returned NULL\n");
    return 1;
  }
  N_VConst(ZERO, u);

  uscale = N_VClone(u);
  if (!uscale)
  {
    fprintf(stderr, "N_VClone returned NULL\n");
    return 1;
  }
  N_VConst(ONE, uscale);

  fscale = N_VClone(u);
  if (!fscale)
  {
    fprintf(stderr, "N_VClone returned NULL\n");
    return 1;
  }
  N_VConst(ONE, fscale);

  /* First check the default path. KINSOL should not create an internal vector
     stack until KINInit has a template vector available. */
  kin_mem = KINCreate(sunctx);
  if (!kin_mem)
  {
    fprintf(stderr, "KINCreate returned NULL\n");
    return 1;
  }

  flag = KINGetVecStack(kin_mem, &stack_check);
  if (check_flag("KINGetVecStack", flag)) { return 1; }
  if (stack_check != NULL)
  {
    fprintf(stderr, "Expected NULL stack before KINInit\n");
    return 1;
  }

  /* KINInit should create an internal vector stack when the user has not already
     attached one. */
  flag = KINInit(kin_mem, sysfn, u);
  if (check_flag("KINInit", flag)) { return 1; }

  flag = KINGetVecStack(kin_mem, &stack_check);
  if (check_flag("KINGetVecStack", flag)) { return 1; }
  if (stack_check == NULL)
  {
    fprintf(stderr, "Expected non-NULL stack after KINInit\n");
    return 1;
  }

  /* Free the first KINSOL memory block. Since KINSOL created the stack above,
     KINFree is responsible for destroying it. */
  KINFree(&kin_mem);

  /* Next check the user-supplied path. This new KINSOL memory block will receive
     a stack owned by the caller before KINInit is called. */
  kin_mem = KINCreate(sunctx);
  if (!kin_mem)
  {
    fprintf(stderr, "KINCreate returned NULL\n");
    return 1;
  }

  /* Seed the user-owned stack with one vector so we can verify after KINFree
     that the stack still exists and has not been destroyed by KINSOL. */
  flag = SUNVecStack_Create(u, 1, sunctx, &stack);
  if (check_flag("SUNVecStack_Create", flag)) { return 1; }

  /* Attach the user-owned stack. KINSOL should store this pointer but leave
     ownership with the caller. */
  flag = KINSetVecStack(kin_mem, stack);
  if (check_flag("KINSetVecStack", flag)) { return 1; }

  /* Initializing with a user-provided stack should reuse that stack instead of
     allocating a new internal one. */
  flag = KINInit(kin_mem, sysfn, u);
  if (check_flag("KINInit", flag)) { return 1; }

  flag = KINGetVecStack(kin_mem, &stack_check);
  if (check_flag("KINGetVecStack", flag)) { return 1; }
  if (stack_check != stack)
  {
    fprintf(stderr, "KINGetVecStack did not return the user stack\n");
    return 1;
  }

  /* KINInit checks out five core workspace vectors from the stack: unew, fval,
     pp, vtemp1, and vtemp2. */
  flag = SUNVecStack_GetNumActiveVecs(stack, &num_active_vecs);
  if (check_flag("SUNVecStack_GetNumActiveVecs", flag)) { return 1; }
  if (num_active_vecs != 5)
  {
    fprintf(stderr, "Expected five active workspace vectors after KINInit\n");
    return 1;
  }

  /* Attempting to replace the stack while KINSOL has checked-out vectors should
     fail, leaving the original user-owned stack attached. */
  flag = KINSetVecStack(kin_mem, stack);
  if (flag != KIN_ILL_INPUT)
  {
    fprintf(stderr, "Expected KINSetVecStack to reject an active stack\n");
    return 1;
  }

  flag = KINGetVecStack(kin_mem, &stack_check);
  if (check_flag("KINGetVecStack", flag)) { return 1; }
  if (stack_check != stack)
  {
    fprintf(stderr, "KINSetVecStack changed the active stack\n");
    return 1;
  }

  /* Run a tiny Picard solve. This exercises the strategy-specific gval checkout
     in addition to the five core vectors from KINInit. */
  if (setup_linear_solver(kin_mem, u, sunctx, &A, &LS)) { return 1; }

  flag = KINSol(kin_mem, u, KIN_PICARD, uscale, fscale);
  if (check_flag("KINSol", flag)) { return 1; }

  flag = SUNVecStack_GetNumActiveVecs(stack, &num_active_vecs);
  if (check_flag("SUNVecStack_GetNumActiveVecs", flag)) { return 1; }
  if (num_active_vecs != 6)
  {
    fprintf(stderr, "Expected six active workspace vectors after KINSol\n");
    return 1;
  }

  /* KINFree must not destroy a user-owned stack. The caller should still be able
     to query and destroy it after the KINSOL memory is gone. */
  KINFree(&kin_mem);

  flag = SUNVecStack_GetNumVecs(stack, &num_vecs);
  if (check_flag("SUNVecStack_GetNumVecs", flag)) { return 1; }
  if (num_vecs != 6)
  {
    fprintf(stderr, "Expected user stack to remain valid after KINFree\n");
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

  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(fscale);
  N_VDestroy(uscale);
  N_VDestroy(u);
  SUNContext_Free(&sunctx);

  printf("SUCCESS\n");

  return 0;
}

/*---- end of file ----*/
