/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
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
 * Test inequality constraints in ARKStep
 * ---------------------------------------------------------------------------*/

#include <cstdio>
#include <iostream>

// Include desired integrators, vectors, linear solvers, and nonlinear solvers
#include "arkode/arkode.h"
#include "arkode/arkode_arkstep.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_context.hpp"
#include "sundials/sundials_math.h"
#include "sundials/sundials_matrix.h"
#include "sundials/sundials_nvector.h"
#include "sundials/sundials_types.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#include "problems/constant.hpp"
#include "utilities/check_return.hpp"

using namespace std;
using namespace problems::constant;

int main(int argc, char* argv[])
{
  cout << "Start ARKStep inequality constraints test" << endl;

  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Create right hand side vector
  N_Vector rhs = N_VNew_Serial(2, sunctx);
  if (check_ptr(rhs, "N_VNew_Serial")) { return 1; }

  sunrealtype* rhs_data = N_VGetArrayPointer(rhs);
  rhs_data[0]           = SUN_RCONST(-2.0);
  rhs_data[1]           = SUN_RCONST(1.0);

  // Create initial condition
  sunrealtype t0 = SUN_RCONST(0.0);
  N_Vector y     = N_VClone(rhs);
  if (check_ptr(y, "N_VClone")) { return 1; }

  sunrealtype* y_data = N_VGetArrayPointer(y);
  y_data[0]           = SUN_RCONST(1.0);
  y_data[1]           = SUN_RCONST(-1.0);

  // Create integrator
  void* arkode_mem = ARKStepCreate(nullptr, ode_rhs, t0, y, sunctx);
  if (check_ptr(arkode_mem, "ARKStepCreate")) { return 1; }

  // Create user data structure
  UserData user_data{t0, y, rhs};

  int flag = ARKodeSetUserData(arkode_mem, &user_data);
  if (check_flag(flag, "ARKodeSetUserData")) { return 1; }

  // Set relative and absolute tolerances
  const sunrealtype rtol = SUN_RCONST(1.0e-6);
  const sunrealtype atol = SUN_RCONST(1.0e-10);

  flag = ARKodeSStolerances(arkode_mem, rtol, atol);
  if (check_flag(flag, "ARKodeSStolerances")) { return 1; }

  // Create and attach matrix and linear solver
  SUNMatrix A = SUNDenseMatrix(2, 2, sunctx);
  if (check_ptr(A, "SUNDenseMatrix")) { return 1; }

  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_ptr(LS, "SUNLinSol_Dense")) { return 1; }

  flag = ARKodeSetLinearSolver(arkode_mem, LS, A);
  if (check_flag(flag, "ARKodeSetLinearSolver")) { return 1; }

  flag = ARKodeSetJacFn(arkode_mem, ode_rhs_jac);
  if (check_flag(flag, "ARKodeSetJacFn")) { return 1; }

  // Set method order
  flag = ARKodeSetOrder(arkode_mem, 2);
  if (check_flag(flag, "ARKodeSetOrder")) { return 1; }

  // Create constraint vector
  N_Vector constraints = N_VClone(rhs);
  if (check_ptr(y, "N_VClone")) { return 1; }

  sunrealtype* c_data = N_VGetArrayPointer(constraints);
  c_data[0]           = SUN_RCONST(1.0);  // >= 0.0
  c_data[1]           = SUN_RCONST(-1.0); // <= 0.0

  flag = ARKodeSetConstraints(arkode_mem, constraints);
  if (check_flag(flag, "ARKodeSetConstraints")) { return 1; }

  // Set the initial time step size such that y = -y0 at t = t0 + h
  N_Vector temp = N_VClone(rhs);
  if (check_ptr(y, "N_VClone")) { return 1; }

  N_VDiv(y, rhs, temp);
  N_VAbs(temp, temp);
  sunrealtype exact_step = N_VMin(temp);                // step to reach y = 0
  sunrealtype init_step = SUN_RCONST(2.0) * exact_step; // step to reach y = -y0

  flag = ARKodeSetInitStep(arkode_mem, init_step);
  if (check_flag(flag, "ARKodeSetInitStep")) { return 1; }

  // Initial time and fist output time
  sunrealtype tret = SUN_RCONST(0.0);
  sunrealtype tout = t0 + init_step;

  // Advance one step
  flag = ARKodeEvolve(arkode_mem, tout, y, &tret, ARK_ONE_STEP);
  if (check_flag(flag, "ARKode")) { return 1; }

  sunrealtype h_init;
  flag = ARKodeGetActualInitStep(arkode_mem, &h_init);
  if (check_flag(flag, "ARKodeGetActualInitStep")) { return 1; }

  sunrealtype h_last;
  flag = ARKodeGetLastStep(arkode_mem, &h_last);
  if (check_flag(flag, "ARKodeGetLastStep")) { return 1; }

  long int num_constraint_fails;
  flag = ARKodeGetNumConstrFails(arkode_mem, &num_constraint_fails);

  cout << endl;
  cout << "y[0]          = " << y_data[0] << endl;
  cout << "y[1]          = " << y_data[1] << endl;
  cout << "Initial step  = " << h_init << endl;
  cout << "Exact step    = " << exact_step << endl;
  cout << "Last step     = " << h_last << endl;
  cout << "Expected step = " << SUN_RCONST(0.9) * exact_step << endl;
  cout << "Failed steps  = " << num_constraint_fails << endl;
  cout << endl;

  // Check solution
  if (y_data[0] < SUN_RCONST(0.0))
  {
    cout << "FAILED: First solution component is negative!" << endl;
    return 1;
  }

  if (y_data[1] > SUN_RCONST(0.0))
  {
    cout << "FAILED: Second solution component is positive!" << endl;
    return 1;
  }

  // Expected step = safety factor * exact step
  if (SUNRCompare(SUN_RCONST(0.9) * exact_step, h_last))
  {
    cout << "FAILED: Expected step and actual step do not match!" << endl;
    return 1;
  }

  // Should only have one constraint failure
  if (num_constraint_fails != 1)
  {
    cout << "FAILED: Expected one constraint failure!" << endl;
    return 1;
  }

  // Print final statistics
  flag = ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(flag, "ARKodePrintAllStats")) { return 1; }

  // Clean up and return
  N_VDestroy(rhs);
  N_VDestroy(y);
  N_VDestroy(temp);
  N_VDestroy(constraints);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  ARKodeFree(&arkode_mem);

  cout << endl << "End ARKStep inequality constraints test" << endl;

  return 0;
}

/*---- end of file ----*/
