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
 * Test logging output in ForcingStep
 * ---------------------------------------------------------------------------*/

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>

// Include desired integrators and vectors
#include "arkode/arkode_erkstep.h"
#include "arkode/arkode_forcingstep.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_context.hpp"
#include "sundials/sundials_logger.h"

#include "problems/estep.hpp"
#include "utilities/check_return.hpp"

using namespace std;
using namespace problems::estep;

int main(int argc, char* argv[])
{
  cout << "Start ForcingStep Logging test" << endl;

  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Ensure logging output goes to stdout
  SUNLogger logger;
  int flag = SUNContext_GetLogger(sunctx, &logger);
  if (check_flag(flag, "SUNContext_GetLogger")) { return 1; }

  SUNLogger_SetErrorFilename(logger, "stdout");
  SUNLogger_SetWarningFilename(logger, "stdout");
  SUNLogger_SetInfoFilename(logger, "stdout");
  SUNLogger_SetDebugFilename(logger, "stdout");

  // Step sizes: overall, partition 1, partition 2
  sunrealtype dt   = SUN_RCONST(0.001);
  sunrealtype dt_1 = dt / 2;
  sunrealtype dt_2 = dt / 4;

  // Create initial condition
  N_Vector y = N_VNew_Serial(1, sunctx);
  if (check_ptr(y, "N_VNew_Serial")) { return 1; }

  flag = initial_condition(y);
  if (check_flag(flag, "initial_condition")) { return 1; }

  // Create partition 1 integrator
  void* stepper_1 = ERKStepCreate(ode_rhs_1, zero, y, sunctx);
  if (check_ptr(stepper_1, "ERKStepCreate")) { return 1; }

  flag = ARKodeSetUserData(stepper_1, &problem_data);
  if (check_flag(flag, "ARKodeSetUserData")) { return 1; }

  flag = ARKodeSetFixedStep(stepper_1, dt_1);
  if (check_flag(flag, "ARKodeSetFixedStep")) { return 1; }

  // Create partition 1 integrator
  void* stepper_2 = ERKStepCreate(ode_rhs_2, zero, y, sunctx);
  if (check_ptr(stepper_2, "ERKStepCreate")) { return 1; }

  flag = ARKodeSetFixedStep(stepper_2, dt_2);
  if (check_flag(flag, "ARKodeSetFixedStep")) { return 1; }

  // Create the overall integrator
  SUNStepper steppers[2];

  flag = ARKodeCreateSUNStepper(stepper_1, &steppers[0]);
  if (check_flag(flag, "ARKodeSetFixedStep")) { return 1; }

  flag = ARKodeCreateSUNStepper(stepper_2, &steppers[1]);
  if (check_flag(flag, "ARKodeSetFixedStep")) { return 1; }

  void* arkode_mem = ForcingStepCreate(steppers[0], steppers[1], zero, y, sunctx);
  if (check_ptr(arkode_mem, "ForcingStepCreate")) { return 1; }

  flag = ARKodeSetFixedStep(arkode_mem, dt);
  if (check_flag(flag, "ARKodeSetFixedStep")) { return 1; }

  // True solution vector
  N_Vector yt = N_VClone(y);

  flag = true_solution(zero, problem_data, yt);
  if (check_flag(flag, "true_solution")) { return 1; }

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

  sunrealtype* y_data  = N_VGetArrayPointer(y);
  sunrealtype* yt_data = N_VGetArrayPointer(yt);

  cout << setw(width) << tret << setw(width) << y_data[0] << setw(width)
       << abs(y_data[0] - yt_data[0]) << endl;

  // Advance in time
  for (int i = 0; i < nout; i++)
  {
    flag = ARKodeEvolve(arkode_mem, tout, y, &tret, ARK_ONE_STEP);
    if (check_flag(flag, "ARKodeEvolve")) { return 1; }

    flag = true_solution(tret, problem_data, yt);
    if (check_flag(flag, "true_solution")) { return 1; }

    cout << setw(width) << tret << setw(width) << y_data[0] << setw(width)
         << abs(y_data[0] - yt_data[0]) << endl;

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
  N_VDestroy(yt);
  ARKodeFree(&arkode_mem);
  ARKodeFree(&stepper_1);
  ARKodeFree(&stepper_2);
  SUNStepper_Destroy(&steppers[0]);
  SUNStepper_Destroy(&steppers[1]);

  cout << "End ForcingStep Logging test" << endl;

  return 0;
}

/*---- end of file ----*/
