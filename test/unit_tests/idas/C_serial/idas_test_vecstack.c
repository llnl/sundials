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

#include "idas/idas.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_vecstack.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

static int dae_res(sunrealtype t, N_Vector y, N_Vector ydot, N_Vector res,
                   void* user_data)
{
  sunrealtype* ydot_data = N_VGetArrayPointer(ydot);
  sunrealtype* res_data  = N_VGetArrayPointer(res);
  res_data[0]            = ydot_data[0] - ONE;
  return 0;
}

static int dae_jac(sunrealtype t, sunrealtype cj, N_Vector y, N_Vector yp,
                   N_Vector rr, SUNMatrix J, void* user_data, N_Vector tempv1,
                   N_Vector tempv2, N_Vector tempv3)
{
  sunrealtype* J_data = SUNDenseMatrix_Data(J);
  J_data[0]           = cj;
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

static int setup_linear_solver(void* ida_mem, N_Vector y, SUNContext sunctx,
                               SUNMatrix* A, SUNLinearSolver* LS)
{
  *A = SUNDenseMatrix(1, 1, sunctx);
  if (!(*A))
  {
    fprintf(stderr, "SUNDenseMatrix returned NULL\n");
    return 1;
  }

  *LS = SUNLinSol_Dense(y, *A, sunctx);
  if (!(*LS))
  {
    fprintf(stderr, "SUNLinSol_Dense returned NULL\n");
    return 1;
  }

  if (check_flag("IDASetLinearSolver", IDASetLinearSolver(ida_mem, *LS, *A)))
  {
    return 1;
  }

  if (check_flag("IDASetJacFn", IDASetJacFn(ida_mem, dae_jac))) { return 1; }

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
  N_Vector yp             = NULL;
  void* ida_mem           = NULL;

  /* Create the shared SUNDIALS context used by IDAS, the N_Vectors, and the
     user-owned vector stack created later in this test. */
  flag = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_flag("SUNContext_Create", flag)) { return 1; }

  /* Create minimal serial vectors for the initial condition.  The y vector
     also acts as the template for vectors allocated by SUNVecStack. */
  y = N_VNew_Serial(1, sunctx);
  if (!y)
  {
    fprintf(stderr, "N_VNew_Serial returned NULL\n");
    return 1;
  }
  N_VConst(ZERO, y);

  yp = N_VClone(y);
  if (!yp)
  {
    fprintf(stderr, "N_VClone returned NULL\n");
    return 1;
  }
  N_VConst(ONE, yp);

  /* First check the default path. IDAS should not create an internal vector
     stack until IDAInit has a template vector available. */
  ida_mem = IDACreate(sunctx);
  if (!ida_mem)
  {
    fprintf(stderr, "IDACreate returned NULL\n");
    return 1;
  }

  flag = IDAGetVecStack(ida_mem, &stack_check);
  if (check_flag("IDAGetVecStack", flag)) { return 1; }
  if (stack_check != NULL)
  {
    fprintf(stderr, "Expected NULL stack before IDAInit\n");
    return 1;
  }

  /* IDAInit should create an internal vector stack when the user has not
     already attached one. */
  flag = IDAInit(ida_mem, dae_res, ZERO, y, yp);
  if (check_flag("IDAInit", flag)) { return 1; }

  flag = IDAGetVecStack(ida_mem, &stack_check);
  if (check_flag("IDAGetVecStack", flag)) { return 1; }
  if (stack_check == NULL)
  {
    fprintf(stderr, "Expected non-NULL stack after IDAInit\n");
    return 1;
  }

  /* Free the first IDAS memory block. Since IDAS created the stack above,
     IDAFree is responsible for destroying it. */
  IDAFree(&ida_mem);

  /* Next check the user-supplied path. This new IDAS memory block will receive
     a stack owned by the caller before IDAInit is called. */
  ida_mem = IDACreate(sunctx);
  if (!ida_mem)
  {
    fprintf(stderr, "IDACreate returned NULL\n");
    return 1;
  }

  /* Seed the user-owned stack with one vector so we can verify after IDAFree
     that the stack still exists and has not been destroyed by IDAS. */
  flag = SUNVecStack_Create(y, 1, sunctx, &stack);
  if (check_flag("SUNVecStack_Create", flag)) { return 1; }

  /* Attach the user-owned stack. IDAS should store this pointer but leave
     ownership with the caller. */
  flag = IDASetVecStack(ida_mem, stack);
  if (check_flag("IDASetVecStack", flag)) { return 1; }

  /* Initializing with a user-provided stack should reuse that stack instead of
     allocating a new internal one. */
  flag = IDAInit(ida_mem, dae_res, ZERO, y, yp);
  if (check_flag("IDAInit", flag)) { return 1; }

  flag = IDAGetVecStack(ida_mem, &stack_check);
  if (check_flag("IDAGetVecStack", flag)) { return 1; }
  if (stack_check != stack)
  {
    fprintf(stderr, "IDAGetVecStack did not return the user stack\n");
    return 1;
  }

  /* The checked-out IDAS workspace vectors should now be active in the
     user-owned stack. This includes ee so local-error queries remain valid
     while IDAS is alive. */
  flag = SUNVecStack_GetNumActiveVecs(stack, &num_active_vecs);
  if (check_flag("SUNVecStack_GetNumActiveVecs", flag)) { return 1; }
  if (num_active_vecs != 8)
  {
    fprintf(stderr, "Expected eight active workspace vectors after IDAInit\n");
    return 1;
  }

  /* Attempting to replace the stack while IDAS has checked-out vectors should
     fail, leaving the original user-owned stack attached. */
  flag = IDASetVecStack(ida_mem, stack);
  if (flag != IDA_ILL_INPUT)
  {
    fprintf(stderr, "Expected IDASetVecStack to reject an active stack\n");
    return 1;
  }

  flag = IDAGetVecStack(ida_mem, &stack_check);
  if (check_flag("IDAGetVecStack", flag)) { return 1; }
  if (stack_check != stack)
  {
    fprintf(stderr, "IDASetVecStack changed the active stack\n");
    return 1;
  }

  /* Advance the problem once, then query the estimated local error. This checks
     that ida_ee still has the required post-IDASolve lifetime even though its
     storage came from SUNVecStack. */
  flag = IDASStolerances(ida_mem, SUN_RCONST(1.0e-4), SUN_RCONST(1.0e-8));
  if (check_flag("IDASStolerances", flag)) { return 1; }

  if (setup_linear_solver(ida_mem, y, sunctx, &A, &LS)) { return 1; }

  ele = N_VClone(y);
  if (!ele)
  {
    fprintf(stderr, "N_VClone returned NULL\n");
    return 1;
  }

  flag = IDASolve(ida_mem, SUN_RCONST(0.1), &tret, y, yp, IDA_NORMAL);
  if (check_flag("IDASolve", flag)) { return 1; }

  flag = IDAGetEstLocalErrors(ida_mem, ele);
  if (check_flag("IDAGetEstLocalErrors", flag)) { return 1; }

  /* IDAFree must not destroy a user-owned stack. The caller should still be
     able to query and destroy it after the IDAS memory is gone. */
  IDAFree(&ida_mem);

  flag = SUNVecStack_GetNumVecs(stack, &num_vecs);
  if (check_flag("SUNVecStack_GetNumVecs", flag)) { return 1; }
  if (num_vecs != 8)
  {
    fprintf(stderr, "Expected user stack to remain valid after IDAFree\n");
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
  N_VDestroy(yp);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  printf("SUCCESS\n");

  return 0;
}

/*---- end of file ----*/
