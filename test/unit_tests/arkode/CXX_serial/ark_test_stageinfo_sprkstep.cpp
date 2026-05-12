/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
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
 * Test stage information and pre/postprocessing routines in SPRKStep
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

#include "problems/kepler.hpp"
#include "utilities/check_return.hpp"

using namespace std;
using namespace problems::kepler;

// Store the main integrator global memory for these tests only, so that we
// can call ARKodeGetLastTime etc in the step and stage pre/postprocessing
// callback functions below.  This would normally be stored in user_data, but
// here we reuse problem definitions from other tests.
void* arkode_mem = nullptr;
static int pre_step(sunrealtype t, N_Vector y, long int step, int attempt,
                    void* user_data);
static int post_step(sunrealtype t, N_Vector y, long int step, void* user_data);
static int pre_rhs(sunrealtype t, N_Vector y, void* user_data);
static int postprocess_step(sunrealtype t, N_Vector y, void* user_data);
static int postprocess_stage(sunrealtype t, N_Vector y, void* user_data);

int main(int argc, char* argv[])
{
  cout << "Start SPRKStep StageInfo test" << endl;

  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Create initial condition
  N_Vector y = N_VNew_Serial(4, sunctx);
  if (check_ptr(y, "N_VNew_Serial")) { return 1; }

  int flag = initial_condition(y, eccentricity);
  if (check_flag(flag, "initial_condition")) { return 1; }

  // Create SPRKStep memory structure
  arkode_mem = SPRKStepCreate(ode_rhs_force, ode_rhs_velocity, zero, y, sunctx);
  if (check_ptr(arkode_mem, "SPKStepCreate")) { return 1; }

  // Step size
  const sunrealtype dt = SUN_RCONST(0.001);
  flag                 = ARKodeSetFixedStep(arkode_mem, dt);
  if (check_flag(flag, "ARKodeSetFixedStep")) { return 1; }

  // Set pre/post step and stage routines
  flag = ARKodeSetPreStepFn(arkode_mem, pre_step);
  if (check_flag(flag, "ARKodeSetPreStepFn")) { return 1; }
  flag = ARKodeSetPostStepFn(arkode_mem, post_step);
  if (check_flag(flag, "ARKodeSetPostStepFn")) { return 1; }
  flag = ARKodeSetPostprocessStepFn(arkode_mem, postprocess_step);
  if (check_flag(flag, "ARKodeSetPostprocessStepFn")) { return 1; }
  flag = ARKodeSetPreRhsFn(arkode_mem, pre_rhs);
  if (check_flag(flag, "ARKodeSetPreRhsFn")) { return 1; }
  flag = ARKodeSetPostprocessStageFn(arkode_mem, postprocess_stage);
  if (check_flag(flag, "ARKodeSetPostprocessStageFn")) { return 1; }

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

  cout << "End SPRKStep StageInfo test" << endl;

  return 0;
}

// Callback functions
static int pre_step(sunrealtype t, N_Vector y, long int step, int attempt,
                    void* user_data)
{
  sunrealtype tn, tcur;
  if (ARKodeGetLastTime(arkode_mem, &tn) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetLastTime" << std::endl;
    return -1;
  }
  if (ARKodeGetCurrentTime(arkode_mem, &tcur) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetCurrentTime" << std::endl;
    return -1;
  }
  std::cout << "    [Pre-step at t = " << std::setprecision(2) << t
            << " (tn = " << tn << " , tcur = " << tcur << "),"
            << " step = " << step << ", attempt = " << attempt << ", "
            << std::setprecision(10)
            << "||y||_2 = " << SUNRsqrt(N_VDotProd(y, y)) << "]" << std::endl
            << std::flush;
  return 0;
}

static int post_step(sunrealtype t, N_Vector y, long int step, void* user_data)
{
  sunrealtype tn, tcur;
  if (ARKodeGetLastTime(arkode_mem, &tn) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetLastTime" << std::endl;
    return -1;
  }
  if (ARKodeGetCurrentTime(arkode_mem, &tcur) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetCurrentTime" << std::endl;
    return -1;
  }
  std::cout << "    [Post-step at t = " << std::setprecision(2) << t
            << " (tn = " << tn << " , tcur = " << tcur << "),"
            << " step = " << step << ", " << std::setprecision(10)
            << "||y||_2 = " << SUNRsqrt(N_VDotProd(y, y)) << "]" << std::endl
            << std::flush;
  return 0;
}

static int postprocess_step(sunrealtype t, N_Vector y, void* user_data)
{
  sunrealtype tn, tcur;
  if (ARKodeGetLastTime(arkode_mem, &tn) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetLastTime" << std::endl;
    return -1;
  }
  if (ARKodeGetCurrentTime(arkode_mem, &tcur) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetCurrentTime" << std::endl;
    return -1;
  }
  std::cout << "    [Postprocess-step at t = " << std::setprecision(2) << t
            << " (tn = " << tn << " , tcur = " << tcur << "),"
            << std::setprecision(10)
            << "||y||_2 = " << SUNRsqrt(N_VDotProd(y, y)) << "]" << std::endl
            << std::flush;
  return 0;
}

static int pre_rhs(sunrealtype t, N_Vector y, void* user_data)
{
  int stage, max_stages;
  sunrealtype tn, tcur;
  if (ARKodeGetLastTime(arkode_mem, &tn) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetLastTime" << std::endl;
    return -1;
  }
  if (ARKodeGetCurrentTime(arkode_mem, &tcur) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetCurrentTime" << std::endl;
    return -1;
  }
  if (ARKodeGetStageIndex(arkode_mem, &stage, &max_stages) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetStageIndex" << std::endl;
    return -1;
  }
  std::cout << "    [Pre-RHS (stage " << stage << " of " << max_stages
            << ") at t = " << std::setprecision(2) << t << " (tn = " << tn
            << " , tcur = " << tcur << "), " << std::setprecision(10)
            << "||y||_2 = " << SUNRsqrt(N_VDotProd(y, y)) << "]" << std::endl
            << std::flush;
  return 0;
}

static int postprocess_stage(sunrealtype t, N_Vector y, void* user_data)
{
  int stage, max_stages;
  sunrealtype tn, tcur;
  if (ARKodeGetLastTime(arkode_mem, &tn) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetLastTime" << std::endl;
    return -1;
  }
  if (ARKodeGetCurrentTime(arkode_mem, &tcur) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetCurrentTime" << std::endl;
    return -1;
  }
  if (ARKodeGetStageIndex(arkode_mem, &stage, &max_stages) != ARK_SUCCESS)
  {
    std::cerr << "Error in ARKodeGetStageIndex" << std::endl;
    return -1;
  }
  std::cout << "    [Post-stage processing (stage " << stage << " of "
            << max_stages << ") at t = " << std::setprecision(2) << t
            << " (tn = " << tn << " , tcur = " << tcur << "), "
            << std::setprecision(10)
            << "||y||_2 = " << SUNRsqrt(N_VDotProd(y, y)) << "]" << std::endl
            << std::flush;
  return 0;
}

/*---- end of file ----*/
