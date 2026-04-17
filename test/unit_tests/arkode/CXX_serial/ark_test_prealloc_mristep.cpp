/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 * based on test_logging_arkode_mristep.cpp by David J. Gardner @ LLNL
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
 * Test data preallocation in MRIStep
 * ---------------------------------------------------------------------------*/

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>

// Include desired integrators, vectors, linear solvers, and nonlinear solvers
#include "arkode/arkode.h"
#include "arkode/arkode_arkstep.h"
#include "arkode/arkode_mristep.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_context.hpp"
#include "sundials/sundials_iterative.h"
#include "sundials/sundials_matrix.h"
#include "sundials/sundials_nonlinearsolver.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#include "problems/kpr.hpp"
#include "utilities/check_return.hpp"

using namespace std;
using namespace problems::kpr;

int main(int argc, char* argv[])
{
  cout << "Start MRIStep preallocation test" << endl;

  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Preallocate data
  bool preallocate_data = false;
  if (argc > 1) { preallocate_data = stoi(argv[1]); }

  // Create initial condition
  N_Vector y = N_VNew_Serial(2, sunctx);
  if (check_ptr(y, "N_VNew_Serial")) { return 1; }

  sunrealtype utrue, vtrue;
  int flag = true_sol(zero, &utrue, &vtrue);
  if (check_flag(flag, "true_sol")) { return 1; }

  sunrealtype* ydata = N_VGetArrayPointer(y);
  ydata[0]           = utrue;
  ydata[1]           = vtrue;

  // Create fast stepper
  void* inner_arkode_mem = ARKStepCreate(ode_rhs_ff, nullptr, zero, y, sunctx);
  if (check_ptr(inner_arkode_mem, "ARKStepCreate")) { return 1; }

  flag = ARKodeSetUserData(inner_arkode_mem, &problem_data);
  if (check_flag(flag, "ARKodeSetUserData")) { return 1; }

  // Relative and absolute tolerances
  const sunrealtype inner_rtol = SUN_RCONST(1.0e-6);
  const sunrealtype inner_atol = SUN_RCONST(1.0e-10);
  flag = ARKodeSStolerances(inner_arkode_mem, inner_rtol, inner_atol);
  if (check_flag(flag, "ARKodeSStolerances")) { return 1; }

  MRIStepInnerStepper stepper = nullptr;
  flag = ARKStepCreateMRIStepInnerStepper(inner_arkode_mem, &stepper);
  if (check_flag(flag, "ARKStepCreateMRIStepInnerStepper")) { return 1; }

  // Create MRIStep memory structure
  void* arkode_mem = nullptr;
  cout << "Using Im-MRI-GARK method" << endl;
  arkode_mem = MRIStepCreate(nullptr, ode_rhs_s, zero, y, stepper, sunctx);
  if (check_ptr(arkode_mem, "MRIStepCreate")) { return 1; }
  flag = ARKodeSetUserData(arkode_mem, &problem_data);
  if (check_flag(flag, "ARKodeSetUserData")) { return 1; }

  // Relative and absolute tolerances
  const sunrealtype rtol = SUN_RCONST(1.0e-6);
  const sunrealtype atol = SUN_RCONST(1.0e-10);
  flag                   = ARKodeSStolerances(arkode_mem, rtol, atol);
  if (check_flag(flag, "ARKodeSStolerances")) { return 1; }

  // Implicit algebraic solvers
  cout << "Using Newton nonlinear solver" << endl;
  cout << "Using dense direct linear solver" << endl;
  SUNMatrix A            = nullptr;
  SUNLinearSolver LS     = nullptr;
  SUNNonlinearSolver NLS = nullptr;
  A                      = SUNDenseMatrix(2, 2, sunctx);
  if (check_ptr(A, "SUNDenseMatrix")) { return 1; }
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_ptr(LS, "SUNLinSol_Dense")) { return 1; }
  flag = ARKodeSetLinearSolver(arkode_mem, LS, A);
  if (check_flag(flag, "ARKodeSetLinearSolver")) { return 1; }
  flag = ARKodeSetJacFn(arkode_mem, ode_rhs_jac);
  if (check_flag(flag, "ARKodeSetJacFn")) { return 1; }

  // Data preallocation (all steppers)
  if (preallocate_data)
  {
    flag = ARKodeInit(inner_arkode_mem);
    if (check_flag(flag, "ARKodeInit")) { return 1; }
    flag = ARKodeInit(arkode_mem);
    if (check_flag(flag, "ARKodeInit")) { return 1; }
  }

  // Initial time and fist output time
  const sunrealtype dtout = one; // output interval
  const int nout          = 3;   // number of outputs
  sunrealtype tret        = zero;
  sunrealtype tout        = tret + dtout;

  // Output initial condition
  cout << scientific;
  cout << setprecision(numeric_limits<sunrealtype>::digits10);
  cout << "           t              ";
  cout << "          u              ";
  cout << "          v              ";
  cout << "        u err            ";
  cout << "        v err            " << endl;
  for (int i = 0; i < 9; i++) { cout << "--------------"; }
  cout << endl;

  cout << setw(22) << tret << setw(25) << ydata[0] << setw(25) << ydata[1]
       << setw(25) << abs(ydata[0] - utrue) << setw(25) << abs(ydata[1] - vtrue)
       << endl;

  // Advance in time
  for (int i = 0; i < nout; i++)
  {
    flag = ARKodeEvolve(arkode_mem, tout, y, &tret, ARK_ONE_STEP);
    if (check_flag(flag, "ARKode")) { return 1; }

    flag = true_sol(tret, &utrue, &vtrue);
    if (check_flag(flag, "true_sol")) { return 1; }

    cout << setw(22) << tret << setw(25) << ydata[0] << setw(25) << ydata[1]
         << setw(25) << abs(ydata[0] - utrue) << setw(25)
         << abs(ydata[1] - vtrue) << endl;

    // update output time
    tout += dtout;
  }
  for (int i = 0; i < 9; i++) { cout << "--------------"; }
  cout << endl;

  // Print some final statistics
  cout << endl << "Outer integrator statistics:" << endl;
  flag = ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(flag, "ARKodePrintAllStats")) { return 1; }
  cout << endl << "Inner integrator statistics:" << endl;
  flag = ARKodePrintAllStats(inner_arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(flag, "ARKodePrintAllStats")) { return 1; }

  // Clean up and return with successful completion
  N_VDestroy(y);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  SUNNonlinSolFree(NLS);
  MRIStepInnerStepper_Free(&stepper);
  ARKodeFree(&inner_arkode_mem);
  ARKodeFree(&arkode_mem);

  cout << "End MRIStep preallocation test" << endl;

  return 0;
}

/*---- end of file ----*/
