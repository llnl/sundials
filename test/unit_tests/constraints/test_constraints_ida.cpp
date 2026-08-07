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
 * Test inequality constraints in IDA
 * ---------------------------------------------------------------------------*/

#include <cstdio>
#include <iostream>

// Include desired integrators, vectors, linear solvers, and nonlinear solvers
#include "ida/ida.h"
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
  cout << "Start IDA inequality constraints test" << endl;

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

  N_Vector yp = N_VClone(rhs);
  if (check_ptr(y, "N_VClone")) { return 1; }

  N_VScale(SUN_RCONST(1.0), rhs, yp);

  // Create integrator
  void* ida_mem = IDACreate(sunctx);
  if (check_ptr(ida_mem, "IDACreate")) { return 1; }

  int flag = IDAInit(ida_mem, dae_res, SUN_RCONST(0.0), y, yp);
  if (check_flag(flag, "IDAInit")) { return 1; }

  // Create user data structure
  UserData user_data{t0, y, rhs};

  flag = IDASetUserData(ida_mem, &user_data);
  if (check_flag(flag, "IDASetUserData")) { return 1; }

  // Set relative and absolute tolerances
  const sunrealtype rtol = SUN_RCONST(1.0e-6);
  const sunrealtype atol = SUN_RCONST(1.0e-10);

  flag = IDASStolerances(ida_mem, rtol, atol);
  if (check_flag(flag, "IDASStolerances")) { return 1; }

  // Create and attach matrix and linear solver
  SUNMatrix A = SUNDenseMatrix(2, 2, sunctx);
  if (check_ptr(A, "SUNDenseMatrix")) { return 1; }

  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_ptr(LS, "SUNLinSol_Dense")) { return 1; }

  flag = IDASetLinearSolver(ida_mem, LS, A);
  if (check_flag(flag, "IDASetLinearSolver")) { return 1; }

  flag = IDASetJacFn(ida_mem, dae_res_jac);
  if (check_flag(flag, "IDASetJacFn")) { return 1; }

  // Create constraint vector
  N_Vector constraints = N_VClone(rhs);
  if (check_ptr(y, "N_VClone")) { return 1; }

  sunrealtype* c_data = N_VGetArrayPointer(constraints);
  c_data[0]           = SUN_RCONST(1.0);  // >= 0.0
  c_data[1]           = SUN_RCONST(-1.0); // <= 0.0

  flag = IDASetConstraints(ida_mem, constraints);
  if (check_flag(flag, "IDASetConstraints")) { return 1; }

  // Set the initial time step size such that y = -y0 at t = t0 + h
  N_Vector temp = N_VClone(rhs);
  if (check_ptr(y, "N_VClone")) { return 1; }

  N_VDiv(y, rhs, temp);
  N_VAbs(temp, temp);
  sunrealtype exact_step = N_VMin(temp);                // step to reach y = 0
  sunrealtype init_step = SUN_RCONST(2.0) * exact_step; // step to reach y = -y0

  flag = IDASetInitStep(ida_mem, init_step);
  if (check_flag(flag, "IDASetInitStep")) { return 1; }

  // Initial time and fist output time
  sunrealtype tret = SUN_RCONST(0.0);
  sunrealtype tout = t0 + init_step;

  // Advance one step
  flag = IDASolve(ida_mem, tout, &tret, y, yp, IDA_ONE_STEP);
  if (check_flag(flag, "IDA")) { return 1; }

  sunrealtype h_init;
  flag = IDAGetActualInitStep(ida_mem, &h_init);
  if (check_flag(flag, "IDAGetActualInitStep")) { return 1; }

  sunrealtype h_last;
  flag = IDAGetLastStep(ida_mem, &h_last);
  if (check_flag(flag, "IDAGetLastStep")) { return 1; }

  long int num_constraint_fails = 0;
  flag = IDAGetNumConstraintFails(ida_mem, &num_constraint_fails);

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
  flag = IDAPrintAllStats(ida_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(flag, "IDAPrintAllStats")) { return 1; }

  // Clean up and return
  N_VDestroy(rhs);
  N_VDestroy(y);
  N_VDestroy(yp);
  N_VDestroy(temp);
  N_VDestroy(constraints);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  IDAFree(&ida_mem);

  cout << endl << "End IDA inequality constraints test" << endl;

  return 0;
}

/*---- end of file ----*/
