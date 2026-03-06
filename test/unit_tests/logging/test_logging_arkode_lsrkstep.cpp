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
 * Test logging output in LSRKStep
 * ---------------------------------------------------------------------------*/

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>

// Include desired integrators, vectors, linear solvers, and nonlinear solvers
#include "arkode/arkode_lsrkstep.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_logger.h"

#include "problems/prv.hpp"
#include "sundials/sundials_nvector.h"
#include "utilities/check_return.hpp"

using namespace std;
using namespace problems::prv;

int main(int argc, char* argv[])
{
  cout << "Start LSRKStep Logging test" << endl;

  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Use RKC (0), RKL (1), SSPs2 (2), SSPs3 (3), SSP43 (4), SSP104 (5)
  int method = 0;
  if (argc > 1) { method = stoi(argv[1]); }

  // Ensure logging output goes to stdout
  SUNLogger logger;
  int flag = SUNContext_GetLogger(sunctx, &logger);
  if (check_flag(flag, "SUNContext_GetLogger")) { return 1; }

  SUNLogger_SetErrorFilename(logger, "stdout");
  SUNLogger_SetWarningFilename(logger, "stdout");
  SUNLogger_SetInfoFilename(logger, "stdout");
  SUNLogger_SetDebugFilename(logger, "stdout");

  // Create initial condition
  N_Vector y = N_VNew_Serial(1, sunctx);
  if (check_ptr(y, "N_VNew_Serial")) { return 1; }
  N_VConst(true_solution(zero), y);

  // Create LSRKStep memory structure
  void* arkode_mem = nullptr;
  if (method < 2) { arkode_mem = LSRKStepCreateSTS(ode_rhs, zero, y, sunctx); }
  else { arkode_mem = LSRKStepCreateSSP(ode_rhs, zero, y, sunctx); }
  if (check_ptr(arkode_mem, "LSRKStepCreate")) { return 1; }

  // Select method
  if (method == 0)
  {
    flag = LSRKStepSetSTSMethodByName(arkode_mem, "ARKODE_LSRK_RKC_2");
  }
  else if (method == 1)
  {
    flag = LSRKStepSetSTSMethodByName(arkode_mem, "ARKODE_LSRK_RKL_2");
  }
  else if (method == 2)
  {
    flag = LSRKStepSetSSPMethodByName(arkode_mem, "ARKODE_LSRK_SSP_S_2");
  }
  else if (method == 3)
  {
    flag = LSRKStepSetSSPMethodByName(arkode_mem, "ARKODE_LSRK_SSP_S_3");
  }
  else if (method == 4)
  {
    flag = LSRKStepSetSSPMethodByName(arkode_mem, "ARKODE_LSRK_SSP_S_3");
    if (flag == 0) { flag = LSRKStepSetNumSSPStages(arkode_mem, 4); }
  }
  else if (method == 5)
  {
    flag = LSRKStepSetSSPMethodByName(arkode_mem, "ARKODE_LSRK_SSP_10_4");
  }
  else
  {
    cerr << "Invalid method option\n";
    return 1;
  }
  if (check_flag(flag, "LSRKStepSetMethodByName")) { return 1; }

  flag = ARKodeSetUserData(arkode_mem, &problem_data);
  if (check_flag(flag, "ARKodeSetUserData")) { return 1; }

  // Relative and absolute tolerances
  const sunrealtype rtol = SUN_RCONST(1.0e-6);
  const sunrealtype atol = SUN_RCONST(1.0e-10);

  flag = ARKodeSStolerances(arkode_mem, rtol, atol);
  if (check_flag(flag, "ARKodeSStolerances")) { return 1; }

  // Specify dominant eigenvalue function
  flag = LSRKStepSetDomEigFn(arkode_mem, ode_dom_eig);
  if (check_flag(flag, "LSRKStepSetDomEigFn")) { return 1; }

  // Initial time and fist output time
  const sunrealtype dtout = one; // output interval
  const int nout          = 3;   // number of outputs
  sunrealtype tret        = zero;
  sunrealtype tout        = tret + dtout;

  const int width = numeric_limits<sunrealtype>::digits10 + 8;

  // Output initial contion
  cout << scientific;
  cout << setprecision(numeric_limits<sunrealtype>::digits10);
  cout << setw(width) << " t";
  cout << setw(width) << " y";
  cout << setw(width) << " y err" << endl;
  for (int i = 0; i < 3 * width; i++) { cout << "-"; }
  cout << endl;

  sunrealtype* y_data = N_VGetArrayPointer(y);

  cout << setw(width) << tret << setw(width) << y_data[0] << setw(width)
       << abs(y_data[0] - true_solution(tret)) << endl;

  // Advance in time
  for (int i = 0; i < nout; i++)
  {
    flag = ARKodeEvolve(arkode_mem, tout, y, &tret, ARK_ONE_STEP);
    if (check_flag(flag, "ARKodeEvolve")) { return 1; }

    cout << setw(width) << tret << setw(width) << y_data[0] << setw(width)
         << abs(y_data[0] - true_solution(tret)) << endl;

    // update output time
    tout += dtout;
  }
  for (int i = 0; i < 3 * width; i++) { cout << "-"; }
  cout << endl;

  // Print some final statistics
  flag = ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(flag, "ARKodePrintAllStats")) { return 1; }

  // Clean up and return with successful completion
  N_VDestroy(y);
  ARKodeFree(&arkode_mem);

  cout << "End LSRKStep Logging test" << endl;

  return 0;
}

/*---- end of file ----*/
