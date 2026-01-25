/* ---------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ---------------------------------------------------------------------------
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
 * ---------------------------------------------------------------------------
 * IMEX multirate Dahlquist problem:
 *
 * y' = lambda_e * y + lambda_i * y + lambda_f * y
 * ---------------------------------------------------------------------------*/

// Header files
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_mristep.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <map>
#include <nvector/nvector_serial.h>
#include <string>
#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>

#include "arkode/arkode_mri_tables_impl.h"

enum class interp_type
{
  none = -1,
  hermite,
  lagrange,
};

// Problem parameters
struct ProblemData
{
  sunrealtype lambda_e = SUN_RCONST(-1.0);
  sunrealtype lambda_i = SUN_RCONST(-1.0);
  sunrealtype lambda_f = SUN_RCONST(-1.0);
};

// Problem options
struct ProblemOptions
{
  // Initial time
  sunrealtype t0 = SUN_RCONST(0.0);

  // Number of time steps
  int nsteps = 1;

  // Relative and absolute tolerances
  sunrealtype reltol = SUN_RCONST(1.0e-4);
  sunrealtype abstol = SUN_RCONST(1.0e-6);

  // Slow and fast step sizes
  sunrealtype hs = SUN_RCONST(0.01);
  sunrealtype hf = SUN_RCONST(0.01);

  // Interpolant type
  // -1 = None
  // 0  = Hermite
  // 1  = Lagrange
  interp_type i_type = interp_type::hermite;
};

// User-supplied Functions called by the solver
static int fe(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int ff(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Ji(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Private function to check function return values
static int check_flag(void* flagvalue, const std::string funcname, int opt);

// Test drivers
static int run_tests(MRISTEP_METHOD_TYPE type, ProblemOptions& prob_opts,
                     ProblemData& prob_data, SUNContext ctx);

// -----------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  // Problem data and options structures
  ProblemData prob_data;
  ProblemOptions prob_opts;

  // Check for inputs
  if (argc > 1)
  {
    if (std::stoi(argv[1]) == -1) { prob_opts.i_type = interp_type::none; }
    else if (std::stoi(argv[1]) == 0)
    {
      prob_opts.i_type = interp_type::hermite;
    }
    else if (std::stoi(argv[1]) == 1)
    {
      prob_opts.i_type = interp_type::lagrange;
    }
    else
    {
      std::cerr << "ERROR: Invalid interpolation type option" << std::endl;
      return 1;
    }
  }

  // Output problem setup
  std::cout << "\nDahlquist ODE test problem:\n";
  std::cout << "   lambda expl  = " << prob_data.lambda_e << "\n";
  std::cout << "   lambda impl  = " << prob_data.lambda_i << "\n";
  std::cout << "   lambda fast  = " << prob_data.lambda_f << "\n";
  std::cout << "   h slow       = " << prob_opts.hs << "\n";
  std::cout << "   h fast       = " << prob_opts.hf << "\n";
  std::cout << "   relative tol = " << prob_opts.reltol << "\n";
  std::cout << "   absolute tol = " << prob_opts.abstol << "\n";
  if (prob_opts.i_type == interp_type::hermite)
  {
    std::cout << "   interp type  = Hermite\n";
  }
  else if (prob_opts.i_type == interp_type::lagrange)
  {
    std::cout << "   interp type  = Lagrange\n";
  }
  else { std::cout << "   interp type  = None\n"; }
  // Create SUNDIALS context
  sundials::Context sunctx;

  // Test methods
  int numfails = 0;

  numfails += run_tests(MRISTEP_EXPLICIT, prob_opts, prob_data, sunctx);

  numfails += run_tests(MRISTEP_IMPLICIT, prob_opts, prob_data, sunctx);

  numfails += run_tests(MRISTEP_IMEX, prob_opts, prob_data, sunctx);

  numfails += run_tests(MRISTEP_MERK, prob_opts, prob_data, sunctx);

  numfails += run_tests(MRISTEP_SR, prob_opts, prob_data, sunctx);

  if (numfails) { std::cout << "\n\nFailed " << numfails << " tests!\n"; }
  else { std::cout << "\n\nAll tests passed!\n"; }

  // Return test status
  return numfails;
}

// -----------------------------------------------------------------------------
// Test drivers
// -----------------------------------------------------------------------------

int run_tests(MRISTEP_METHOD_TYPE type, ProblemOptions& prob_opts,
              ProblemData& prob_data, SUNContext sunctx)
{
  // Reusable error-checking flag
  int flag;

  // Test failure counter
  int numfails = 0;

  // Create initial condition vector
  N_Vector y = N_VNew_Serial(1, sunctx);
  if (check_flag((void*)y, "N_VNew_Serial", 0)) { return 1; }

  N_VConst(SUN_RCONST(1.0), y);

  // Create matrix and linear solver (if necessary)
  SUNMatrix A        = nullptr;
  SUNLinearSolver LS = nullptr;

  if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX || type == MRISTEP_SR)
  {
    // Initialize dense matrix data structures and solvers
    A = SUNDenseMatrix(1, 1, sunctx);
    if (check_flag((void*)A, "SUNDenseMatrix", 0)) { return 1; }

    LS = SUNLinSol_Dense(y, A, sunctx);
    if (check_flag((void*)LS, "SUNLinSol_Dense", 0)) { return 1; }
  }

  // ----------------------
  // Create fast integrator
  // ----------------------

  // Create explicit fast integrator
  void* arkstep_mem = ARKStepCreate(ff, nullptr, prob_opts.t0, y, sunctx);
  if (check_flag((void*)arkstep_mem, "ARKStepCreate", 0)) { return 1; }

  // Set user data
  flag = ARKodeSetUserData(arkstep_mem, &prob_data);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }

  // Specify tolerances
  flag = ARKodeSStolerances(arkstep_mem, prob_opts.reltol, prob_opts.abstol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) { return 1; }

  // Specify fixed time step size
  flag = ARKodeSetFixedStep(arkstep_mem, prob_opts.hf);
  if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }

  // Lagrange interpolant (removes additional RHS evaluation with DIRK methods)
  if (prob_opts.i_type == interp_type::lagrange)
  {
    flag = ARKodeSetInterpolantType(arkstep_mem, ARK_INTERP_LAGRANGE);
    if (check_flag(&flag, "ARKStepSetInterpolantType", 1)) { return 1; }
  }
  else if (prob_opts.i_type == interp_type::none)
  {
    flag = ARKodeSetInterpolantType(arkstep_mem, ARK_INTERP_NONE);
    if (check_flag(&flag, "ARKStepSetInterpolantType", 1)) { return 1; }
  }

  // Wrap ARKStep integrator as fast integrator object
  MRIStepInnerStepper inner_stepper = nullptr;
  flag = ARKodeCreateMRIStepInnerStepper(arkstep_mem, &inner_stepper);
  if (check_flag(&flag, "ARKodeCreateMRIStepInnerStepper", 1)) { return 1; }

  // ----------------------
  // Create slow integrator
  // ----------------------

  // Create slow integrator based on MRI type
  void* mristep_mem = nullptr;

  if ((type == MRISTEP_EXPLICIT) || (type == MRISTEP_MERK))
  {
    mristep_mem = MRIStepCreate(fe, nullptr, prob_opts.t0, y, inner_stepper,
                                sunctx);
  }
  else if (type == MRISTEP_IMPLICIT)
  {
    mristep_mem = MRIStepCreate(nullptr, fi, prob_opts.t0, y, inner_stepper,
                                sunctx);
  }
  else if ((type == MRISTEP_IMEX) || (type == MRISTEP_SR))
  {
    mristep_mem = MRIStepCreate(fe, fi, prob_opts.t0, y, inner_stepper, sunctx);
  }
  else { return 1; }
  if (check_flag((void*)mristep_mem, "MRIStepCreate", 0)) { return 1; }

  // Set user data
  flag = ARKodeSetUserData(mristep_mem, &prob_data);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }

  // Specify tolerances
  flag = ARKodeSStolerances(mristep_mem, prob_opts.reltol, prob_opts.abstol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) { return 1; }

  // Specify fixed time step sizes
  flag = ARKodeSetFixedStep(mristep_mem, prob_opts.hs);
  if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }

  if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX || type == MRISTEP_SR)
  {
    // Attach linear solver
    flag = ARKodeSetLinearSolver(mristep_mem, LS, A);
    if (check_flag(&flag, "ARKodeSetLinearSolver", 1)) { return 1; }

    // Set Jacobian function
    flag = ARKodeSetJacFn(mristep_mem, Ji);
    if (check_flag(&flag, "ARKodeSetJacFn", 1)) { return 1; }

    // Specify linearly implicit RHS, with non-time-dependent Jacobian
    flag = ARKodeSetLinear(mristep_mem, 0);
    if (check_flag(&flag, "ARKodeSetLinear", 1)) { return 1; }
  }

  // Lagrange interpolant (removes additional RHS evaluation with DIRK methods)
  if (prob_opts.i_type == interp_type::lagrange)
  {
    flag = ARKodeSetInterpolantType(mristep_mem, ARK_INTERP_LAGRANGE);
    if (check_flag(&flag, "ARKStepSetInterpolantType", 1)) { return 1; }
  }
  else if (prob_opts.i_type == interp_type::none)
  {
    flag = ARKodeSetInterpolantType(mristep_mem, ARK_INTERP_NONE);
    if (check_flag(&flag, "ARKStepSetInterpolantType", 1)) { return 1; }
  }

  // -------------------------------
  // Evolve with various MRI methods
  // -------------------------------

  // Methods to test paired with whether they are stiffly accurate
  std::map<std::string, bool> methods;

  if (type == MRISTEP_EXPLICIT)
  {
    std::cout << "\n=========================\n";
    std::cout << "Test explicit MRI methods\n";
    std::cout << "=========================\n";

    methods.insert({{"ARKODE_MRI_GARK_FORWARD_EULER", false},
                    {"ARKODE_MRI_GARK_ERK22a", false},
                    {"ARKODE_MRI_GARK_ERK22b", false},
                    {"ARKODE_MRI_GARK_RALSTON2", false},
                    {"ARKODE_MIS_KW3", false},
                    {"ARKODE_MRI_GARK_ERK33a", false},
                    {"ARKODE_MRI_GARK_RALSTON3", false},
                    {"ARKODE_MRI_GARK_ERK45a", false}});
  }
  else if (type == MRISTEP_MERK)
  {
    std::cout << "\n=================\n";
    std::cout << "Test MERK methods\n";
    std::cout << "=================\n";

    methods.insert({{"ARKODE_MERK21", true},
                    {"ARKODE_MERK32", true},
                    {"ARKODE_MERK43", true},
                    {"ARKODE_MERK54", true}});
  }
  else if (type == MRISTEP_IMPLICIT)
  {
    std::cout << "\n=========================\n";
    std::cout << "Test implicit MRI methods\n";
    std::cout << "=========================\n";

    methods.insert({{"ARKODE_MRI_GARK_BACKWARD_EULER", true},
                    {"ARKODE_MRI_GARK_IRK21a", true},
                    {"ARKODE_MRI_GARK_IMPLICIT_MIDPOINT", false},
                    {"ARKODE_MRI_GARK_ESDIRK34a", true},
                    {"ARKODE_MRI_GARK_ESDIRK46a", true}});
  }
  else if (type == MRISTEP_IMEX)
  {
    std::cout << "\n=====================\n";
    std::cout << "Test IMEX MRI methods\n";
    std::cout << "=====================\n";

    methods.insert({{"ARKODE_IMEX_MRI_GARK_EULER", true},
                    {"ARKODE_IMEX_MRI_GARK_TRAPEZOIDAL", false},
                    {"ARKODE_IMEX_MRI_GARK_MIDPOINT", false},
                    {"ARKODE_IMEX_MRI_GARK3a", false},
                    {"ARKODE_IMEX_MRI_GARK3b", false},
                    {"ARKODE_IMEX_MRI_GARK4", false}});
  }
  else if (type == MRISTEP_SR)
  {
    std::cout << "\n========================\n";
    std::cout << "Test IMEX MRI SR methods\n";
    std::cout << "========================\n";

    methods.insert({{"ARKODE_IMEX_MRI_SR21", true},
                    {"ARKODE_IMEX_MRI_SR32", true},
                    {"ARKODE_IMEX_MRI_SR43", true}});
  }
  else { return 1; }

  for (const auto& pair : methods)
  {
    int methodfails       = 0;
    std::string id        = pair.first;
    bool stiffly_accurate = pair.second;
    std::cout << "\nTesting method " << id << "\n";

    // -------------
    // Select method
    // -------------

    // Load method table
    MRIStepCoupling C = MRIStepCoupling_LoadTableByName(id.c_str());
    if (check_flag((void*)C, "MRIStepCoupling_LoadTable", 0)) { return 1; }

    MRIStepCoupling_Write(C, stdout);

    // Get the number of stored stages
    int* stage_map = new int[C->stages];
    int nstages_stored;

    flag = mriStepCoupling_GetStageMap(C, stage_map, &nstages_stored);
    if (check_flag(&flag, "mriStepCoupling_GetStageMap", 1)) { return 1; }

    std::cout << "  Stored stages = " << nstages_stored << "\n";
    delete[] stage_map;

    // Set coupling table
    flag = MRIStepSetCoupling(mristep_mem, C);
    if (check_flag(&flag, "MRIStepSetCoupling", 1)) { return 1; }

    // -----------------
    // Output statistics
    // -----------------

    sunrealtype t  = prob_opts.t0;
    sunrealtype tf = prob_opts.nsteps * prob_opts.hs;

    for (int i = 0; i < prob_opts.nsteps; i++)
    {
      // Advance in time
      flag = ARKodeEvolve(mristep_mem, tf, y, &t, ARK_ONE_STEP);
      if (check_flag(&flag, "ARKodeEvolve", 1)) { return 1; }

      // Update output time
      tf += prob_opts.hs;
    }

    // -----------------
    // Output statistics
    // -----------------

    long int mri_nst, mri_nfse, mri_nfsi;     // integrator
    long int mri_nni, mri_ncfn;               // nonlinear solver
    long int mri_nsetups, mri_nje, mri_nfeLS; // linear solver

    flag = ARKodeGetNumSteps(mristep_mem, &mri_nst);
    if (check_flag(&flag, "ARKodeGetNumSteps", 1)) { return 1; }

    flag = ARKodeGetNumRhsEvals(mristep_mem, 0, &mri_nfse);
    if (check_flag(&flag, "ARKodeGetNumRhsEvals", 1)) { return 1; }

    flag = ARKodeGetNumRhsEvals(mristep_mem, 1, &mri_nfsi);
    if (check_flag(&flag, "ARKodeGetNumRhsEvals", 1)) { return 1; }

    if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX || type == MRISTEP_SR)
    {
      flag = ARKodeGetNumNonlinSolvIters(mristep_mem, &mri_nni);
      if (check_flag(&flag, "ARKodeGetNumNonlinSolvIters", 1)) { return 1; }

      flag = ARKodeGetNumNonlinSolvConvFails(mristep_mem, &mri_ncfn);
      if (check_flag(&flag, "ARKodeGetNumNonlinSolvConvFails", 1)) { return 1; }

      flag = ARKodeGetNumLinSolvSetups(mristep_mem, &mri_nsetups);
      if (check_flag(&flag, "ARKodeGetNumLinSolvSetups", 1)) { return 1; }

      flag = ARKodeGetNumJacEvals(mristep_mem, &mri_nje);
      if (check_flag(&flag, "ARKodeGetNumJacEvals", 1)) { return 1; }

      flag = ARKodeGetNumLinRhsEvals(mristep_mem, &mri_nfeLS);
      check_flag(&flag, "ARKodeGetNumLinRhsEvals", 1);
    }

    sunrealtype pow = prob_data.lambda_f;
    if (type != MRISTEP_IMPLICIT) { pow += prob_data.lambda_e; }
    if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX || type == MRISTEP_SR)
    {
      pow += prob_data.lambda_i;
    }
    sunrealtype ytrue = exp(pow * t);

    sunrealtype* ydata = N_VGetArrayPointer(y);
    sunrealtype error  = ytrue - ydata[0];

    std::cout << "\nMRIStep Statistics:\n";
    std::cout << "   Time        = " << t << "\n";
    std::cout << "   y(t)        = " << ytrue << "\n";
    std::cout << "   y_n         = " << ydata[0] << "\n";
    std::cout << "   Error       = " << error << "\n";
    std::cout << "   Steps       = " << mri_nst << "\n";
    std::cout << "   Fe evals    = " << mri_nfse << "\n";
    std::cout << "   Fi evals    = " << mri_nfsi << "\n";

    if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX || type == MRISTEP_SR)
    {
      std::cout << "   NLS iters   = " << mri_nni << "\n";
      std::cout << "   NLS fails   = " << mri_ncfn << "\n";
      std::cout << "   LS setups   = " << mri_nsetups << "\n";
      std::cout << "   LS Fi evals = " << mri_nfeLS << "\n";
      std::cout << "   Ji evals    = " << mri_nje << "\n";
    }

    // ----------------
    // Check statistics
    // ----------------

    std::cout << "\nComparing Solver Statistics:\n";

    int nstages_evaluated = nstages_stored;
    if (stiffly_accurate) nstages_evaluated--;
    long int fe_evals = 0;
    if (type != MRISTEP_IMPLICIT) { fe_evals = mri_nst * nstages_evaluated; }

    if (mri_nfse != fe_evals)
    {
      methodfails++;
      std::cout << "Fe RHS evals: " << mri_nfse << " vs " << fe_evals << "\n";
    }

    long int fi_evals = 0;
    if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX || type == MRISTEP_SR)
    {
      fi_evals = mri_nst * nstages_evaluated + mri_nni;
    }

    if (mri_nfsi != fi_evals)
    {
      methodfails++;
      std::cout << "Fi RHS evals: " << mri_nfsi << " vs " << fi_evals << "\n";
    }

    if (methodfails) { std::cout << "Failed " << methodfails << " tests\n"; }
    else { std::cout << "All checks passed\n"; }
    numfails += methodfails;

    // -------------------
    // Setup for next test
    // -------------------

    // Free coupling table
    MRIStepCoupling_Free(C);

    // Reset state vector to the initial condition
    N_VConst(SUN_RCONST(1.0), y);

    // Re-initialize fast integrator
    flag = ARKStepReInit(arkstep_mem, ff, nullptr, prob_opts.t0, y);
    if (check_flag(&flag, "ARKStepReInit", 1)) { return 1; }

    // Re-initialize slow integrator based on MRI type
    if ((type == MRISTEP_EXPLICIT) || (type == MRISTEP_MERK))
    {
      flag = MRIStepReInit(mristep_mem, fe, nullptr, prob_opts.t0, y);
    }
    else if (type == MRISTEP_IMPLICIT)
    {
      flag = MRIStepReInit(mristep_mem, nullptr, fi, prob_opts.t0, y);
    }
    else if (type == MRISTEP_IMEX || type == MRISTEP_SR)
    {
      flag = MRIStepReInit(mristep_mem, fe, fi, prob_opts.t0, y);
    }
    else { return 1; }
    if (check_flag(&flag, "MRIStepReInit", 1)) { return 1; }
  }

  // Clean up
  MRIStepInnerStepper_Free(&inner_stepper);
  ARKodeFree(&mristep_mem);
  ARKodeFree(&arkstep_mem);
  if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX || type == MRISTEP_SR)
  {
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
  }
  N_VDestroy(y);

  return numfails;
}

// -----------------------------------------------------------------------------
// Functions called by the solver
// -----------------------------------------------------------------------------

// Explicit ODE RHS function fe(t,y)
static int fe(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* y_data    = N_VGetArrayPointer(y);
  sunrealtype* yd_data   = N_VGetArrayPointer(ydot);
  ProblemData* prob_data = static_cast<ProblemData*>(user_data);

  yd_data[0] = prob_data->lambda_e * y_data[0];

  return 0;
}

// Implicit ODE RHS function fi(t,y)
static int fi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* y_data    = N_VGetArrayPointer(y);
  sunrealtype* yd_data   = N_VGetArrayPointer(ydot);
  ProblemData* prob_data = static_cast<ProblemData*>(user_data);

  yd_data[0] = prob_data->lambda_i * y_data[0];

  return 0;
}

// Fast ODE RHS function ff(t,y)
static int ff(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* y_data    = N_VGetArrayPointer(y);
  sunrealtype* yd_data   = N_VGetArrayPointer(ydot);
  ProblemData* prob_data = static_cast<ProblemData*>(user_data);

  yd_data[0] = prob_data->lambda_f * y_data[0];

  return 0;
}

// Jacobian routine to compute J(t,y) = dfi/dy.
static int Ji(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* J_data    = SUNDenseMatrix_Data(J);
  ProblemData* prob_data = static_cast<ProblemData*>(user_data);

  J_data[0] = prob_data->lambda_i;

  return 0;
}

// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Check function return value
static int check_flag(void* flagvalue, const std::string funcname, int opt)
{
  int* errflag;

  // Check if function returned NULL pointer - no memory allocated
  if (opt == 0 && flagvalue == nullptr)
  {
    std::cerr << "\nMEMORY_ERROR: " << funcname
              << " failed - returned NULL pointer\n\n";
    return 1;
  }
  // Check if flag < 0
  else if (opt == 1)
  {
    errflag = (int*)flagvalue;
    if (*errflag < 0)
    {
      std::cerr << "\nSUNDIALS_ERROR: " << funcname
                << " failed with flag = " << *errflag << "\n\n";
      return 1;
    }
  }

  return 0;
}
