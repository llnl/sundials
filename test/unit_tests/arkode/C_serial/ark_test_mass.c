/* -----------------------------------------------------------------------------
 * Programmer(s): Steven Roberts @ LLNL
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
 * Unit test that checks the number of mass matrix solves taken in 1 step of an
 * ARKODE integrator.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode/arkode_arkstep.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunnonlinsol/sunnonlinsol_newton.h"

/* A simple nonlinear RHS function */
static int f_explicit(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VProd(y, y, ydot);
  return 0;
}

/* A simple nonlinear RHS function */
static int f_implicit(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VInv(y, ydot);
  return 0;
}

/* A simple mass matrix which can be fixed or time-dependent based on user
 * data */
static int mass(sunrealtype t, SUNMatrix mass, void* user_data, N_Vector tmp1,
                N_Vector tmp2, N_Vector tmp3)
{
  sunbooleantype time_dep  = *(sunbooleantype*)user_data;
  SM_ELEMENT_D(mass, 0, 0) = time_dep ? SUNRexp(t) : 2.0;
  return 0;
}

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a value so check if
             retval >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_retval(void* flagvalue, const char* funcname, int opt)
{
  int* errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
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

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}

static int solve(const char* im, const char* ex, int steps,
                 sunbooleantype time_dep, sunbooleantype deduce_implicit_rhs,
                 long int expected_mass_solves)
{
  int retval = 0;
  int s;
  sunrealtype t               = SUN_RCONST(1.0);
  SUNContext sunctx           = NULL;
  N_Vector y                  = NULL;
  SUNMatrix jacobian_mat      = NULL;
  SUNMatrix mass_mat          = NULL;
  SUNLinearSolver jacobian_ls = NULL;
  SUNLinearSolver mass_ls     = NULL;
  SUNNonlinearSolver nls      = NULL;
  void* arkode_mem            = NULL;
  long int actual_mass_solves;

  /* Create the SUNDIALS context object for this simulation. */
  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return 1;

  /* Create solution vector */
  y = N_VNew_Serial(1, sunctx);
  if (check_retval(y, "N_VNew_Serial", 0)) return 1;
  N_VConst(1.0, y);

  /* Create the matrices */
  jacobian_mat = SUNDenseMatrix(1, 1, sunctx);
  if (check_retval(jacobian_mat, "N_VNew_Serial", 0)) return 1;

  mass_mat = SUNDenseMatrix(1, 1, sunctx);
  if (check_retval(mass_mat, "N_VNew_Serial", 0)) return 1;

  /* Create ARKODE mem structure */
  arkode_mem = ARKStepCreate(f_explicit, f_implicit, SUN_RCONST(0.0), y, sunctx);
  if (check_retval(arkode_mem, "ARKStepCreate", 0)) return 1;

  retval = ARKStepSetTableName(arkode_mem, im, ex);
  if (check_retval(&retval, "ARKStepSetTableNum", 1)) return 1;

  retval = ARKodeSetUserData(arkode_mem, &time_dep);
  if (check_retval(&retval, "ARKodeSetUserData", 1)) return 1;

  /* Specify a time step so no extra mass evaluations occur from initial step
   * size procedure */
  retval = ARKodeSetInitStep(arkode_mem, SUN_RCONST(1.0e-5));
  if (check_retval(&retval, "ARKodeSetInitStep", 1)) return 1;

  /* Use Lagrange interpolation because Hermite may need addition mass matrix
   * solves to compute the interpolant */
  retval = ARKodeSetInterpolantType(arkode_mem, ARK_INTERP_LAGRANGE);
  if (check_retval(&retval, "ARKodeSetInterpolantType", 1)) return 1;

  /* Configure the solvers */
  jacobian_ls = SUNLinSol_Dense(y, jacobian_mat, sunctx);
  if (check_retval(jacobian_ls, "SUNLinSol_Dense", 0)) return 1;

  retval = ARKodeSetLinearSolver(arkode_mem, jacobian_ls, jacobian_mat);
  if (check_retval(&retval, "ARKodeSetLinearSolver", 1)) return 1;

  mass_ls = SUNLinSol_Dense(y, mass_mat, sunctx);
  if (check_retval(mass_ls, "SUNLinSol_Dense", 0)) return 1;

  retval = ARKodeSetMassLinearSolver(arkode_mem, mass_ls, mass_mat, time_dep);
  if (check_retval(&retval, "ARKodeSetMassLinearSolver", 0)) return 1;

  nls = SUNNonlinSol_Newton(y, sunctx);
  if (check_retval(nls, "SUNNonlinSol_Newton", 0)) return 1;

  retval = ARKodeSetDeduceImplicitRhs(arkode_mem, deduce_implicit_rhs);
  if (check_retval(&retval, "ARKodeSetDeduceImplicitRhs", 1)) return 1;

  /* Let ARKODE estimate the Jacobian with finite differences */
  retval = ARKodeSetJacFn(arkode_mem, NULL);
  if (check_retval(&retval, "ARKodeSetJacFn", 1)) return 1;

  retval = ARKodeSetMassFn(arkode_mem, mass);
  if (check_retval(&retval, "ARKodeSetMassFn", 1)) return 1;

  /* Take time step(s) */
  for (s = 0; s < steps; s++)
  {
    retval = ARKodeEvolve(arkode_mem, t, y, &t, ARK_ONE_STEP);
    if (check_retval(&retval, "ARKodeEvolve", 1)) return 1;
  }

  retval = ARKodeGetNumMassSolves(arkode_mem, &actual_mass_solves);
  if (check_retval(&retval, "ARKodeGetNumMassSolves", 1)) return 1;

  /* Free integrator memory */
  ARKodeFree(&arkode_mem);

  /* Clean up */
  N_VDestroy(y);
  SUNMatDestroy(jacobian_mat);
  SUNMatDestroy(mass_mat);
  SUNLinSolFree(jacobian_ls);
  SUNLinSolFree(mass_ls);
  SUNNonlinSolFree(nls);
  SUNContext_Free(&sunctx);

  retval = actual_mass_solves != expected_mass_solves;
  printf("%6s | %-29s| %-30s| %-6d| %-9s| %-14s| %-7ld| %-9ld\n",
         retval ? "Fail" : "Pass", im, ex, steps, time_dep ? "Yes" : "No ",
         deduce_implicit_rhs ? "Yes" : "No ", actual_mass_solves,
         expected_mass_solves);
  return retval;
}

/* Main program */
int main(int argc, char* argv[])
{
  int retval = 0;

  printf("Mass Matrix Solve Count Test\n\n");
  printf("Result | Implicit Method              | Explicit Method              "
         " | Steps | Time-Dep | Deduce Im RHS | Actual | Expected\n");
  printf("-------+------------------------------+------------------------------"
         "-+-------+----------+---------------+--------+---------\n");

  /* Tests taking 1 timestep */

  /* SDIRK */
  /* 1 for primary and 1 for embedded solution. Optimally, could be 0 solves */
  retval += solve("ARKODE_SDIRK_2_1_2", "ARKODE_ERK_NONE", 1, SUNFALSE,
                  SUNFALSE, 2);
  /* 1 for primary and 1 for embedded solution. Optimally, could be 0 solves */
  retval += solve("ARKODE_SDIRK_2_1_2", "ARKODE_ERK_NONE", 1, SUNFALSE, SUNTRUE,
                  2);
  /* 1 per stage */
  retval += solve("ARKODE_SDIRK_2_1_2", "ARKODE_ERK_NONE", 1, SUNTRUE, SUNFALSE,
                  2);
  /* Optimal */
  retval += solve("ARKODE_SDIRK_2_1_2", "ARKODE_ERK_NONE", 1, SUNTRUE, SUNTRUE,
                  0);

  /* Stiffly-accurate SDIRK */
  /* 1 for embedded solution. 0 needed for primary due to stiff accuracy
   * property. Optimally, could be 0 solves */
  retval += solve("ARKODE_SDIRK_5_3_4", "ARKODE_ERK_NONE", 1, SUNFALSE,
                  SUNFALSE, 1);
  /* 1 for embedded solution. 0 needed for primary due to stiff accuracy
   * property. Optimally, could be 0 solves */
  retval += solve("ARKODE_SDIRK_5_3_4", "ARKODE_ERK_NONE", 1, SUNFALSE, SUNTRUE,
                  1);
  /* 1 per stage */
  retval += solve("ARKODE_SDIRK_5_3_4", "ARKODE_ERK_NONE", 1, SUNTRUE, SUNFALSE,
                  5);
  /* Optimal */
  retval += solve("ARKODE_SDIRK_5_3_4", "ARKODE_ERK_NONE", 1, SUNTRUE, SUNTRUE,
                  0);

  /* FSAL ESDIRK */
  /* 1 for embedded solution. 0 needed for primary due to FSAL property. The
   * first step is computing y'_0 = M^{-1} f(y_0), which is needed for Hermite
   * interpolation, but not Lagrange. Practically, could be 1 solve. Technically
   * the optimal is 0 solves if we express embedded solution as linear
   * combination of Y_i. This requires d to be in the rowspace of A. Since this
   * method has the FSAL property this small optimization would only save 1
   * solve on the first step. */
  retval += solve("ARKODE_ESDIRK324L2SA_4_2_3", "ARKODE_ERK_NONE", 1, SUNFALSE,
                  SUNFALSE, 2);
  /* Same comment as previous */
  retval += solve("ARKODE_ESDIRK324L2SA_4_2_3", "ARKODE_ERK_NONE", 1, SUNFALSE,
                  SUNTRUE, 2);
  /* 1 per stage */
  retval += solve("ARKODE_ESDIRK324L2SA_4_2_3", "ARKODE_ERK_NONE", 1, SUNTRUE,
                  SUNFALSE, 4);
  /* Optimal */
  retval += solve("ARKODE_ESDIRK324L2SA_4_2_3", "ARKODE_ERK_NONE", 1, SUNTRUE,
                  SUNTRUE, 1);

  /* ERK */
  /* 1 per stage, 1 for primary solution, and 1 for embedded solution.
   * Optimally, could be 3 solves */
  retval += solve("ARKODE_DIRK_NONE", "ARKODE_SHU_OSHER_3_2_3", 1, SUNFALSE,
                  SUNFALSE, 5);
  /* 1 per stage, 1 for primary solution, and 1 for embedded solution.
   * Optimally, could be 3 solves */
  retval += solve("ARKODE_DIRK_NONE", "ARKODE_SHU_OSHER_3_2_3", 1, SUNFALSE,
                  SUNTRUE, 5);
  /* 1 per stage. Optimal */
  retval += solve("ARKODE_DIRK_NONE", "ARKODE_SHU_OSHER_3_2_3", 1, SUNTRUE,
                  SUNFALSE, 3);
  /* 1 per stage. Optimal */
  retval += solve("ARKODE_DIRK_NONE", "ARKODE_SHU_OSHER_3_2_3", 1, SUNTRUE,
                  SUNTRUE, 3);

  /* FSAL ERK */
  /* 1 per stage and 1 for embedded solution. 0 needed for primary due to FSAL
   * property. Optimally, could be 4 solves */
  retval += solve("ARKODE_DIRK_NONE", "ARKODE_BOGACKI_SHAMPINE_4_2_3", 1,
                  SUNFALSE, SUNFALSE, 5);
  /* 1 per stage and 1 for embedded solution. 0 needed for primary due to FSAL
   * property. Optimally, could be 4 solves */
  retval += solve("ARKODE_DIRK_NONE", "ARKODE_BOGACKI_SHAMPINE_4_2_3", 1,
                  SUNFALSE, SUNTRUE, 5);
  /* 1 per stage. Optimal */
  retval += solve("ARKODE_DIRK_NONE", "ARKODE_BOGACKI_SHAMPINE_4_2_3", 1,
                  SUNTRUE, SUNFALSE, 4);
  /* 1 per stage. Optimal */
  retval += solve("ARKODE_DIRK_NONE", "ARKODE_BOGACKI_SHAMPINE_4_2_3", 1,
                  SUNTRUE, SUNTRUE, 4);

  /* IMEX 1 */
  /* 1 for primary solution and 1 for embedded solution. The first step is
   * computing y'_0 = M^{-1} f(y_0), which is needed for Hermite interpolation,
   * but not Lagrange. Optimally, could be 2 solves */
  retval += solve("ARKODE_ARK2_DIRK_3_1_2", "ARKODE_ARK2_ERK_3_1_2", 1,
                  SUNFALSE, SUNFALSE, 3);
  /* Same comment as previous */
  retval += solve("ARKODE_ARK2_DIRK_3_1_2", "ARKODE_ARK2_ERK_3_1_2", 1,
                  SUNFALSE, SUNTRUE, 3);
  /* 1 per implicit and explicit stage */
  retval += solve("ARKODE_ARK2_DIRK_3_1_2", "ARKODE_ARK2_ERK_3_1_2", 1, SUNTRUE,
                  SUNFALSE, 6);
  /* 1 per explicit stage and 1 for the first stage of implicit method.
   * Optimal */
  retval += solve("ARKODE_ARK2_DIRK_3_1_2", "ARKODE_ARK2_ERK_3_1_2", 1, SUNTRUE,
                  SUNTRUE, 4);

  /* IMEX 2 */
  /* 1 for primary solution and 1 for embedded solution. The first step is
   * computing y'_0 = M^{-1} f(y_0), which is needed for Hermite interpolation,
   * but not Lagrange. Optimally, could be 2 solves */
  retval += solve("ARKODE_ARK548L2SA_DIRK_8_4_5", "ARKODE_ARK548L2SA_ERK_8_4_5",
                  1, SUNFALSE, SUNFALSE, 3);
  /* Same comment as previous */
  retval += solve("ARKODE_ARK548L2SA_DIRK_8_4_5", "ARKODE_ARK548L2SA_ERK_8_4_5",
                  1, SUNFALSE, SUNTRUE, 3);
  /* 1 per implicit and explicit stage */
  retval += solve("ARKODE_ARK548L2SA_DIRK_8_4_5", "ARKODE_ARK548L2SA_ERK_8_4_5",
                  1, SUNTRUE, SUNFALSE, 16);
  /* 1 per explicit stage and 1 for the first stage of implicit method.
   * Optimal */
  retval += solve("ARKODE_ARK548L2SA_DIRK_8_4_5", "ARKODE_ARK548L2SA_ERK_8_4_5",
                  1, SUNTRUE, SUNTRUE, 9);

  /* 2 timestep tests to check FSAL optimizations. This assumes there are no
   * rejected steps */

  /* FSAL ESDIRK */
  /* 2 per step. The fixed mass matrix implementation cannot benefit from the
   * FSAL property */
  retval += solve("ARKODE_ESDIRK324L2SA_4_2_3", "ARKODE_ERK_NONE", 2, SUNFALSE,
                  SUNFALSE, 4);
  /* 2 per step. The fixed mass matrix implementation cannot benefit from the
   * FSAL property */
  retval += solve("ARKODE_ESDIRK324L2SA_4_2_3", "ARKODE_ERK_NONE", 2, SUNFALSE,
                  SUNTRUE, 4);
  /* 1 per stage except the first stage of step 2 due to FSAL optimization */
  retval += solve("ARKODE_ESDIRK324L2SA_4_2_3", "ARKODE_ERK_NONE", 2, SUNTRUE,
                  SUNFALSE, 7);
  /* Optimal */
  retval += solve("ARKODE_ESDIRK324L2SA_4_2_3", "ARKODE_ERK_NONE", 2, SUNTRUE,
                  SUNTRUE, 1);

  /* FSAL ERK */
  /* 1 per stage and 1 for embedded solution. 0 needed for primary due to FSAL
   * property. Optimally, could be 7 solves */
  retval += solve("ARKODE_DIRK_NONE", "ARKODE_BOGACKI_SHAMPINE_4_2_3", 2,
                  SUNFALSE, SUNFALSE, 10);
  /* 1 per stage and 1 for embedded solution. 0 needed for primary due to FSAL
   * property. Optimally, could be 7 solves */
  retval += solve("ARKODE_DIRK_NONE", "ARKODE_BOGACKI_SHAMPINE_4_2_3", 2,
                  SUNFALSE, SUNTRUE, 10);
  /* Optimal */
  retval += solve("ARKODE_DIRK_NONE", "ARKODE_BOGACKI_SHAMPINE_4_2_3", 2,
                  SUNTRUE, SUNFALSE, 7);
  /* Optimal */
  retval += solve("ARKODE_DIRK_NONE", "ARKODE_BOGACKI_SHAMPINE_4_2_3", 2,
                  SUNTRUE, SUNTRUE, 7);

  /* IMEX */
  /* While the implicit part has the FSAL property, the overall IMEX method does
   * not. The mass solves are all double the 1 step results */
  retval += solve("ARKODE_ARK548L2SA_DIRK_8_4_5", "ARKODE_ARK548L2SA_ERK_8_4_5",
                  2, SUNFALSE, SUNFALSE, 6);
  retval += solve("ARKODE_ARK548L2SA_DIRK_8_4_5", "ARKODE_ARK548L2SA_ERK_8_4_5",
                  2, SUNFALSE, SUNTRUE, 6);
  retval += solve("ARKODE_ARK548L2SA_DIRK_8_4_5", "ARKODE_ARK548L2SA_ERK_8_4_5",
                  2, SUNTRUE, SUNFALSE, 32);
  retval += solve("ARKODE_ARK548L2SA_DIRK_8_4_5", "ARKODE_ARK548L2SA_ERK_8_4_5",
                  2, SUNTRUE, SUNTRUE, 18);

  return retval;
}
