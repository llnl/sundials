/* -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Unit test for missing linear solver validation in ARKStep
 * ---------------------------------------------------------------------------*/

#include <stdio.h>

#include "arkode/arkode_arkstep.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_context.h"
#include "sundials/sundials_matrix.h"
#include "sundials/sundials_nvector.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"
#include "sunnonlinsol/sunnonlinsol_newton.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

enum
{
  DEFAULT_NLS,
  CUSTOM_FIXEDPOINT_NLS,
  CUSTOM_NEWTON_NLS
};

static int zero_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  (void)t;
  (void)y;
  (void)user_data;
  N_VConst(ZERO, ydot);
  return 0;
}

static int setup_stepper(void** arkode_mem, N_Vector y, SUNContext sunctx,
                         int use_explicit_method)
{
  int flag;

  *arkode_mem = ARKStepCreate(zero_rhs, zero_rhs, ZERO, y, sunctx);
  if (*arkode_mem == NULL) { return 1; }

  flag = ARKodeSetFixedStep(*arkode_mem, SUN_RCONST(0.1));
  if (flag != ARK_SUCCESS) { return 1; }

  if (use_explicit_method)
  {
    flag = ARKStepSetTableNum(*arkode_mem, ARKODE_DIRK_NONE,
                              ARKODE_HEUN_EULER_2_1_2);
  }
  else
  {
    flag = ARKStepSetTableNum(*arkode_mem, ARKODE_BACKWARD_EULER_1_1,
                              ARKODE_ERK_NONE);
  }
  if (flag != ARK_SUCCESS) { return 1; }

  return 0;
}

static int check_case(N_Vector y, SUNContext sunctx, int use_explicit_method,
                      int custom_nls_type, int attach_linear_solver,
                      int expected_init_flag, const char* description)
{
  SUNMatrix A               = NULL;
  SUNLinearSolver LS        = NULL;
  SUNNonlinearSolver NLS    = NULL;
  void* arkode_mem          = NULL;
  sunrealtype tret          = ZERO;
  int flag                  = 0;
  int test_result           = 0;

  N_VConst(ZERO, y);

  if (setup_stepper(&arkode_mem, y, sunctx, use_explicit_method))
  {
    fprintf(stderr, "%s: failed to set up ARKStep\n", description);
    test_result = 1;
    goto cleanup;
  }

  if (custom_nls_type == CUSTOM_FIXEDPOINT_NLS ||
      custom_nls_type == CUSTOM_NEWTON_NLS)
  {
    NLS = custom_nls_type == CUSTOM_FIXEDPOINT_NLS
            ? SUNNonlinSol_FixedPoint(y, 0, sunctx)
            : SUNNonlinSol_Newton(y, sunctx);
    if (NLS == NULL)
    {
      fprintf(stderr, "%s: failed to create custom nonlinear solver\n",
              description);
      test_result = 1;
      goto cleanup;
    }

    flag = ARKodeSetNonlinearSolver(arkode_mem, NLS);
    if (flag != ARK_SUCCESS)
    {
      fprintf(stderr, "%s: failed to attach fixed-point solver\n", description);
      test_result = 1;
      goto cleanup;
    }
  }

  if (attach_linear_solver)
  {
    A = SUNDenseMatrix(1, 1, sunctx);
    if (A == NULL)
    {
      fprintf(stderr, "%s: failed to create dense matrix\n", description);
      test_result = 1;
      goto cleanup;
    }

    LS = SUNLinSol_Dense(y, A, sunctx);
    if (LS == NULL)
    {
      fprintf(stderr, "%s: failed to create dense linear solver\n", description);
      test_result = 1;
      goto cleanup;
    }

    flag = ARKodeSetLinearSolver(arkode_mem, LS, A);
    if (flag != ARK_SUCCESS)
    {
      fprintf(stderr, "%s: failed to attach dense linear solver\n", description);
      test_result = 1;
      goto cleanup;
    }
  }

  flag = ARKodeInit(arkode_mem);
  if (flag != expected_init_flag)
  {
    fprintf(stderr, "%s: ARKodeInit returned %d, expected %d\n", description,
            flag, expected_init_flag);
    test_result = 1;
    goto cleanup;
  }

  if (expected_init_flag == ARK_SUCCESS)
  {
    flag = ARKodeEvolve(arkode_mem, ONE, y, &tret, ARK_NORMAL);
    if (flag < ARK_SUCCESS)
    {
      fprintf(stderr, "%s: ARKodeEvolve returned %d\n", description, flag);
      test_result = 1;
      goto cleanup;
    }
  }

cleanup:
  if (arkode_mem != NULL) { ARKodeFree(&arkode_mem); }
  if (NLS != NULL) { SUNNonlinSolFree(NLS); }
  if (LS != NULL) { SUNLinSolFree(LS); }
  if (A != NULL) { SUNMatDestroy(A); }

  return test_result;
}

int main(void)
{
  SUNContext sunctx = NULL;
  N_Vector y        = NULL;
  int flag          = 0;

  flag = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (flag != 0)
  {
    fprintf(stderr, "SUNContext_Create returned %d\n", flag);
    return 1;
  }

  y = N_VNew_Serial(1, sunctx);
  if (y == NULL)
  {
    fprintf(stderr, "N_VNew_Serial failed\n");
    SUNContext_Free(&sunctx);
    return 1;
  }

  if (check_case(y, sunctx, 0, DEFAULT_NLS, 0, ARK_NLS_INIT_FAIL,
                 "default Newton without a linear solver") ||
      check_case(y, sunctx, 0, DEFAULT_NLS, 1, ARK_SUCCESS,
                 "default Newton with a dense linear solver") ||
      check_case(y, sunctx, 0, CUSTOM_FIXEDPOINT_NLS, 0, ARK_SUCCESS,
                 "custom fixed-point solver without a linear solver") ||
      check_case(y, sunctx, 0, CUSTOM_NEWTON_NLS, 0, ARK_NLS_INIT_FAIL,
                 "custom Newton without a linear solver") ||
      check_case(y, sunctx, 1, DEFAULT_NLS, 0, ARK_SUCCESS,
                 "explicit method without a linear solver"))
  {
    flag = 1;
  }

  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  return flag;
}
