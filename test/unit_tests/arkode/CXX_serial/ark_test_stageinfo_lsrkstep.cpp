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
 * Test stage information and pre/postprocessing routines in LSRKStep
 * ---------------------------------------------------------------------------*/

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>

// Include desired integrators, vectors, linear solvers, and nonlinear solvers
#include "arkode/arkode_lsrkstep.h"
#include "nvector/nvector_serial.h"

#include "problems/prv.hpp"
#include "sundials/sundials_nvector.h"
#include "utilities/check_return.hpp"

using namespace std;
using namespace problems::prv;

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
  cout << "Start LSRKStep StageInfo test" << endl;

  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Use RKC (0), RKL (1), SSPs2 (2), SSPs3 (3), SSP43 (4), SSP104 (5)
  int method = 0;
  if (argc > 1) { method = stoi(argv[1]); }

  // Create initial condition
  N_Vector y = N_VNew_Serial(1, sunctx);
  if (check_ptr(y, "N_VNew_Serial")) { return 1; }
  N_VConst(true_solution(zero), y);

  // Create LSRKStep memory structure
  if (method < 2) { arkode_mem = LSRKStepCreateSTS(ode_rhs, zero, y, sunctx); }
  else { arkode_mem = LSRKStepCreateSSP(ode_rhs, zero, y, sunctx); }
  if (check_ptr(arkode_mem, "LSRKStepCreate")) { return 1; }

  // Select method
  int flag;
  if (method == 0)
  {
    cout << "Using RKC method" << endl;
    flag = LSRKStepSetSTSMethodByName(arkode_mem, "ARKODE_LSRK_RKC_2");
  }
  else if (method == 1)
  {
    cout << "Using RKL method" << endl;
    flag = LSRKStepSetSTSMethodByName(arkode_mem, "ARKODE_LSRK_RKL_2");
  }
  else if (method == 2)
  {
    cout << "Using SSP(s,2) method" << endl;
    flag = LSRKStepSetSSPMethodByName(arkode_mem, "ARKODE_LSRK_SSP_S_2");
  }
  else if (method == 3)
  {
    cout << "Using SSP(9,3) method" << endl;
    flag = LSRKStepSetSSPMethodByName(arkode_mem, "ARKODE_LSRK_SSP_S_3");
    if (flag == 0) { flag = LSRKStepSetNumSSPStages(arkode_mem, 4); }
  }
  else if (method == 4)
  {
    cout << "Using SSP(4,3) method" << endl;
    flag = LSRKStepSetSSPMethodByName(arkode_mem, "ARKODE_LSRK_SSP_S_3");
    if (flag == 0) { flag = LSRKStepSetNumSSPStages(arkode_mem, 4); }
  }
  else if (method == 5)
  {
    cout << "Using SSP(10,4) method" << endl;
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

  cout << "End LSRKStep StageInfo test" << endl;

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
