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
 * Test logging output in SPRKStep
 * ---------------------------------------------------------------------------*/

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>

// Include desired integrators, vectors, linear solvers, and nonlinear solvers
#include "arkode/arkode_sprkstep.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_context.hpp"
#include "sundials/sundials_logger.h"

#include "problems/kepler.hpp"
#include "utilities/check_return.hpp"

using namespace std;
using namespace problems::kepler;

int main(int argc, char* argv[])
{
  cout << "Start SPRKStep Logging test" << endl;

  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Use compensated summation
  bool use_compensated_sum = false;
  if (argc > 1) { use_compensated_sum = stoi(argv[1]); }

  // Ensure logging output goes to stdout
  SUNLogger logger;
  int flag = SUNContext_GetLogger(sunctx, &logger);
  if (check_flag(flag, "SUNContext_GetLogger")) { return 1; }

  SUNLogger_SetErrorFilename(logger, "stdout");
  SUNLogger_SetWarningFilename(logger, "stdout");
  SUNLogger_SetInfoFilename(logger, "stdout");
  SUNLogger_SetDebugFilename(logger, "stdout");

  // Create initial condition
  N_Vector y = N_VNew_Serial(4, sunctx);
  if (check_ptr(y, "N_VNew_Serial")) { return 1; }

  flag = initial_condition(y, eccentricity);
  if (check_flag(flag, "initial_condition")) { return 1; }

  // Create SPRKStep memory structure
  void* arkode_mem = SPRKStepCreate(ode_rhs_force, ode_rhs_velocity, zero, y,
                                    sunctx);
  if (check_ptr(arkode_mem, "SPKStepCreate")) { return 1; }

  // Step size
  const sunrealtype dt = SUN_RCONST(0.001);
  flag                 = ARKodeSetFixedStep(arkode_mem, dt);
  if (check_flag(flag, "ARKodeSetFixedStep")) { return 1; }

  // Compensated summation
  flag = SPRKStepSetUseCompensatedSums(arkode_mem, use_compensated_sum);
  if (check_flag(flag, "SPRKStepSetUseCompensatedSums")) { return 1; }

  // Initial time and fist output time
  const sunrealtype dtout = dt; // output interval
  const int nout          = 3;  // number of outputs
  sunrealtype tret        = zero;
  sunrealtype tout        = tret + dtout;

  // Output initial contion
  sunrealtype* ydata = N_VGetArrayPointer(y);
  if (check_ptr(y, "N_VGetArrayPointer")) { return 1; }

  cout << scientific;
  cout << setprecision(numeric_limits<sunrealtype>::digits10);
  cout << "           t              ";
  cout << "          q1             ";
  cout << "          q2             ";
  cout << "        q3               ";
  cout << "        q4               " << endl;
  for (int i = 0; i < 9; i++) { cout << "--------------"; }
  cout << endl;

  cout << setw(22) << tret << setw(25) << ydata[0] << setw(25) << ydata[1]
       << setw(25) << ydata[2] << setw(25) << ydata[3] << endl;

  // Advance in time
  for (int i = 0; i < nout; i++)
  {
    flag = ARKodeEvolve(arkode_mem, tout, y, &tret, ARK_ONE_STEP);
    if (check_flag(flag, "ARKodeEvolve")) { return 1; }

    cout << setw(22) << tret << setw(25) << ydata[0] << setw(25) << ydata[1]
         << setw(25) << ydata[2] << setw(25) << ydata[3] << endl;

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
  ARKodeFree(&arkode_mem);

  cout << "End SPRKStep Logging test" << endl;

  return 0;
}

/*---- end of file ----*/
