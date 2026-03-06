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
 * Test logging output in MRIStep
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
#include "sundials/sundials_logger.h"
#include "sundials/sundials_matrix.h"
#include "sundials/sundials_nonlinearsolver.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"

#include "problems/kpr.hpp"
#include "utilities/check_return.hpp"

using namespace std;
using namespace problems::kpr;

int main(int argc, char* argv[])
{
  cout << "Start MRIStep Logging test" << endl;

  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Method to use:
  //   0 = Ex-MRI-GARK
  //   1 = Im-MRI-GARK
  //   2 = ImEx-MRI-GARK
  //   3 = Ex-MRI-SR
  //   4 = Im-MRI-SR
  //   5 = ImEx-MRI-SR
  //   6 = MERK
  int method_type = 0;
  if (argc > 1) { method_type = stoi(argv[1]); }

  // Use Newton (1) otherwise use fixed-point
  bool newton = true;
  if (argc > 2) { newton = stoi(argv[2]); }

  // Use direct dense solver (1) otherwise use GMRES
  bool direct = true;
  if (argc > 3) { direct = stoi(argv[3]); }

  // Ensure logging output goes to stdout
  SUNLogger logger;
  int flag = SUNContext_GetLogger(sunctx, &logger);
  if (check_flag(flag, "SUNContext_GetLogger")) { return 1; }

  SUNLogger_SetErrorFilename(logger, "stdout");
  SUNLogger_SetWarningFilename(logger, "stdout");
  SUNLogger_SetInfoFilename(logger, "stdout");
  SUNLogger_SetDebugFilename(logger, "stdout");

  // Create initial condition
  N_Vector y = N_VNew_Serial(2, sunctx);
  if (check_ptr(y, "N_VNew_Serial")) { return 1; }

  sunrealtype utrue, vtrue;
  flag = true_sol(zero, &utrue, &vtrue);
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
  if (method_type == 0)
  {
    cout << "Using Ex-MRI-GARK method" << endl;
    arkode_mem = MRIStepCreate(ode_rhs_s, nullptr, zero, y, stepper, sunctx);
    if (check_ptr(arkode_mem, "MRIStepCreate")) { return 1; }
  }
  else if (method_type == 1)
  {
    cout << "Using Im-MRI-GARK method" << endl;
    arkode_mem = MRIStepCreate(nullptr, ode_rhs_s, zero, y, stepper, sunctx);
    if (check_ptr(arkode_mem, "MRIStepCreate")) { return 1; }
  }
  else if (method_type == 2)
  {
    cout << "Using ImEx-MRI-GARK method" << endl;
    arkode_mem = MRIStepCreate(ode_rhs_se, ode_rhs_si, zero, y, stepper, sunctx);
    if (check_ptr(arkode_mem, "MRIStepCreate")) { return 1; }
  }
  else if (method_type == 3)
  {
    cout << "Using Ex-MRI-SR method" << endl;
    arkode_mem = MRIStepCreate(ode_rhs_s, nullptr, zero, y, stepper, sunctx);
    if (check_ptr(arkode_mem, "MRIStepCreate")) { return 1; }
  }
  else if (method_type == 4)
  {
    cout << "Using Im-MRI-SR method" << endl;
    arkode_mem = MRIStepCreate(nullptr, ode_rhs_s, zero, y, stepper, sunctx);
    if (check_ptr(arkode_mem, "MRIStepCreate")) { return 1; }
  }
  else if (method_type == 5)
  {
    cout << "Using ImEx-MRI-SR method" << endl;
    arkode_mem = MRIStepCreate(ode_rhs_se, ode_rhs_si, zero, y, stepper, sunctx);
    if (check_ptr(arkode_mem, "MRIStepCreate")) { return 1; }
  }
  else if (method_type == 6)
  {
    cout << "Using MERK method" << endl;
    arkode_mem = MRIStepCreate(ode_rhs_s, nullptr, zero, y, stepper, sunctx);
    if (check_ptr(arkode_mem, "MRIStepCreate")) { return 1; }
  }
  else
  {
    cerr << "ERROR: Invalid method type option " << method_type;
    return 1;
  }

  // Set MRI-SR or MERK method
  MRIStepCoupling C = nullptr;
  if (method_type == 3 || method_type == 4 || method_type == 5)
  {
    C = MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_SR21);
    if (check_ptr(C, "MRIStepCoupling_LoadTable")) { return 1; }
  }
  else if (method_type == 6)
  {
    C = MRIStepCoupling_LoadTable(ARKODE_MERK21);
    if (check_ptr(C, "MRIStepCoupling_LoadTable")) { return 1; }
  }

  if (C)
  {
    flag = MRIStepSetCoupling(arkode_mem, C);
    if (check_flag(flag, "MRIStepSetCoupling")) { return 1; }
    MRIStepCoupling_Free(C);
  }

  flag = ARKodeSetUserData(arkode_mem, &problem_data);
  if (check_flag(flag, "ARKodeSetUserData")) { return 1; }

  // Relative and absolute tolerances
  const sunrealtype rtol = SUN_RCONST(1.0e-6);
  const sunrealtype atol = SUN_RCONST(1.0e-10);

  flag = ARKodeSStolerances(arkode_mem, rtol, atol);
  if (check_flag(flag, "ARKodeSStolerances")) { return 1; }

  // Set the slow step size
  const sunrealtype hfixed = SUN_RCONST(0.001);

  flag = ARKodeSetFixedStep(arkode_mem, hfixed);
  if (check_flag(flag, "ARKodeSetFixedStep")) { return 1; }

  SUNMatrix A            = nullptr;
  SUNLinearSolver LS     = nullptr;
  SUNNonlinearSolver NLS = nullptr;

  if (method_type == 1 || method_type == 2 || method_type == 4 || method_type == 5)
  {
    if (!newton)
    {
      cout << "Using fixed-point nonlinear solver" << endl;

      NLS = SUNNonlinSol_FixedPoint(y, 0, sunctx);
      if (check_ptr(NLS, "SUNNonlinSol_FixedPoint")) { return 1; }

      flag = ARKodeSetNonlinearSolver(arkode_mem, NLS);
      if (check_flag(flag, "ARKodeSetLinearSolver")) { return 1; }
    }
    else
    {
      cout << "Using Newton nonlinear solver" << endl;
      if (direct)
      {
        cout << "Using dense direct linear solver" << endl;

        A = SUNDenseMatrix(2, 2, sunctx);
        if (check_ptr(A, "SUNDenseMatrix")) { return 1; }

        LS = SUNLinSol_Dense(y, A, sunctx);
        if (check_ptr(LS, "SUNLinSol_Dense")) { return 1; }

        flag = ARKodeSetLinearSolver(arkode_mem, LS, A);
        if (check_flag(flag, "ARKodeSetLinearSolver")) { return 1; }

        if (method_type == 1 || method_type == 4)
        {
          // implicit slow
          flag = ARKodeSetJacFn(arkode_mem, ode_rhs_jac);
          if (check_flag(flag, "ARKodeSetJacFn")) { return 1; }
        }
        else if (method_type == 2 || method_type == 5)
        {
          // implicit-explicit slow
          flag = ARKodeSetJacFn(arkode_mem, ode_rhs_jac_im);
          if (check_flag(flag, "ARKodeSetJacFn")) { return 1; }
        }
        else
        {
          cerr << "ERROR: Unknown splitting option";
          return 1;
        }
      }
      else
      {
        cout << "Using GMRES iterative linear solver" << endl;

        LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, sunctx);
        if (check_ptr(LS, "SUNLinSol_SPGMR")) { return 1; }

        flag = ARKodeSetLinearSolver(arkode_mem, LS, A);
        if (check_flag(flag, "ARKodeSetLinearSolver")) { return 1; }
      }
    }
  }

  // Initial time and fist output time
  const sunrealtype dtout = one; // output interval
  const int nout          = 3;   // number of outputs
  sunrealtype tret        = zero;
  sunrealtype tout        = tret + dtout;

  // Output initial contion
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
  flag = ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(flag, "ARKodePrintAllStats")) { return 1; }

  // Clean up and return with successful completion
  N_VDestroy(y);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  SUNNonlinSolFree(NLS);
  MRIStepInnerStepper_Free(&stepper);
  ARKodeFree(&inner_arkode_mem);
  ARKodeFree(&arkode_mem);

  cout << "End MRIStep Logging test" << endl;

  return 0;
}

/*---- end of file ----*/
