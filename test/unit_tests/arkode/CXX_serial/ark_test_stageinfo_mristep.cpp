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
 * Test stage information and pre/postprocessing routines in MRIStep
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
#include "sunlinsol/sunlinsol_spgmr.h"

#include "problems/kpr.hpp"
#include "utilities/check_return.hpp"

using namespace std;
using namespace problems::kpr;

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
  cout << "Start MRIStep StageInfo test" << endl;

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
    C = MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_SR32);
    if (check_ptr(C, "MRIStepCoupling_LoadTable")) { return 1; }
  }
  else if (method_type == 6)
  {
    C = MRIStepCoupling_LoadTable(ARKODE_MERK32);
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

  SUNMatrix A        = nullptr;
  SUNLinearSolver LS = nullptr;

  if (method_type == 1 || method_type == 2 || method_type == 4 || method_type == 5)
  {
    cout << "Using Newton nonlinear solver" << endl;
    cout << "Using GMRES iterative linear solver" << endl;

    LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, sunctx);
    if (check_ptr(LS, "SUNLinSol_SPGMR")) { return 1; }

    flag = ARKodeSetLinearSolver(arkode_mem, LS, A);
    if (check_flag(flag, "ARKodeSetLinearSolver")) { return 1; }
  }

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
  MRIStepInnerStepper_Free(&stepper);
  ARKodeFree(&inner_arkode_mem);
  ARKodeFree(&arkode_mem);

  cout << "End MRIStep StageInfo test" << endl;

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
