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
 * Test logging output in IDA
 * ---------------------------------------------------------------------------*/

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>

// Include desired integrators, vectors, linear solvers, and nonlinear solvers
#include "idas/idas.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_context.hpp"
#include "sundials/sundials_iterative.h"
#include "sundials/sundials_logger.h"
#include "sundials/sundials_matrix.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunmatrix/sunmatrix_dense.h"

#include "problems/kpr.hpp"
#include "utilities/check_return.hpp"

using namespace std;
using namespace problems::kpr;

int main(int argc, char* argv[])
{
  cout << "Start IDAS Logging test" << endl;

  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Use direct dense solver (1) otherwise use GMRES
  bool direct = true;
  if (argc > 1) { direct = stoi(argv[1]); }

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

  N_Vector yp = N_VNew_Serial(2, sunctx);
  if (check_ptr(y, "N_VNew_Serial")) { return 1; }

  sunrealtype uptrue, vptrue;
  flag = true_sol_p(zero, &uptrue, &vptrue);
  if (check_flag(flag, "true_sol")) { return 1; }

  sunrealtype* ypdata = N_VGetArrayPointer(yp);
  ypdata[0]           = uptrue;
  ypdata[1]           = vptrue;

  // Create IDA memory structure
  void* ida_mem = IDACreate(sunctx);
  if (check_ptr(ida_mem, "IDACreate")) { return 1; }

  flag = IDAInit(ida_mem, dae_res, zero, y, yp);
  if (check_flag(flag, "IDAInit")) { return 1; }

  flag = IDASetUserData(ida_mem, &problem_data);
  if (check_flag(flag, "IDASetUserData")) { return 1; }

  // Relative and absolute tolerances
  const sunrealtype rtol = SUN_RCONST(1.0e-6);
  const sunrealtype atol = SUN_RCONST(1.0e-10);

  flag = IDASStolerances(ida_mem, rtol, atol);
  if (check_flag(flag, "IDASStolerances")) { return 1; }

  cout << "Using Newton nonlinear solver" << endl;

  SUNMatrix A        = nullptr;
  SUNLinearSolver LS = nullptr;

  if (direct)
  {
    cout << "Using dense direct linear solver" << endl;

    A = SUNDenseMatrix(2, 2, sunctx);
    if (check_ptr(A, "SUNDenseMatrix")) { return 1; }

    LS = SUNLinSol_Dense(y, A, sunctx);
    if (check_ptr(LS, "SUNLinSol_Dense")) { return 1; }

    flag = IDASetLinearSolver(ida_mem, LS, A);
    if (check_flag(flag, "IDASetLinearSolver")) { return 1; }

    flag = IDASetJacFn(ida_mem, dae_res_jac);
    if (check_flag(flag, "IDASetJacFn")) { return 1; }
  }
  else
  {
    cout << "Using GMRES iterative linear solver" << endl;

    LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, sunctx);
    if (check_ptr(LS, "SUNLinSol_SPGMR")) { return 1; }

    flag = IDASetLinearSolver(ida_mem, LS, A);
    if (check_flag(flag, "IDASetLinearSolver")) { return 1; }
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
    flag = IDASolve(ida_mem, tout, &tret, y, yp, IDA_ONE_STEP);
    if (check_flag(flag, "IDA")) { return 1; }

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
  flag = IDAPrintAllStats(ida_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(flag, "IDAPrintAllStats")) { return 1; }

  // Clean up and return with successful completion
  N_VDestroy(y);
  N_VDestroy(yp);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  IDAFree(&ida_mem);

  cout << "End IDAS Logging test" << endl;

  return 0;
}

/*---- end of file ----*/
