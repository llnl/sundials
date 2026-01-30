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
 * Dahlquist problem:
 *
 * y' = lambda_e * y
 * ---------------------------------------------------------------------------*/

// Header files
#include <arkode/arkode_erkstep.h>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <nvector/nvector_serial.h>
#include <string>
#include <sundials/sundials_core.hpp>

#include "arkode/arkode_butcher.h"

// Constants
#define NEG_ONE SUN_RCONST(-1.0)
#define ZERO    SUN_RCONST(0.0)
#define TWO     SUN_RCONST(2.0)

enum class interp_type
{
  none = -1,
  hermite,
  lagrange
};

// Problem parameters
struct ProblemData
{
  sunrealtype lambda_e = NEG_ONE;
  sunrealtype lambda_i = NEG_ONE;
};

// Problem options
struct ProblemOptions
{
  // Initial time
  sunrealtype t0 = ZERO;

  // Relative and absolute tolerances
  sunrealtype reltol = SUN_RCONST(1.0e-4);
  sunrealtype abstol = SUN_RCONST(1.0e-6);

  // Step size
  sunrealtype h = SUN_RCONST(0.01);

  // Interpolant type
  // -1 = None
  // 0  = Hermite
  // 1  = Lagrange
  interp_type i_type = interp_type::hermite;
};

// User-supplied Functions called by the solver
int fe(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

// Private function to check function return values
int check_flag(void* flagvalue, const std::string funcname, int opt);

// Test drivers
int run_tests(ARKodeButcherTable Be, ProblemData& prob_data,
              ProblemOptions& prob_opts, sundials::Context& sunctx);

int get_method_properties(ARKodeButcherTable Be, int& stages, int& order,
                          bool& explicit_first_stage, bool& stiffly_accurate,
                          bool& fsal);

int expected_rhs_evals(interp_type i_type, int stages,
                       bool explicit_first_stage, bool stiffly_accurate,
                       bool fsal, void* erkstep_mem, long int& nfe_expected);

int check_rhs_evals(void* erkstep_mem, long int nfe_expected);

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
  std::cout << "\nDahlquist ODE test problem:\n"
            << "  lambda expl  = " << prob_data.lambda_e << "\n"
            << "  step size    = " << prob_opts.h << "\n"
            << "  relative tol = " << prob_opts.reltol << "\n"
            << "  absolute tol = " << prob_opts.abstol << "\n";
  if (prob_opts.i_type == interp_type::hermite)
  {
    std::cout << "  interp type  = Hermite\n";
  }
  else if (prob_opts.i_type == interp_type::lagrange)
  {
    std::cout << "  interp type  = Lagrange\n";
  }
  else { std::cout << "  interp type  = None\n"; }

  // Create SUNDIALS context
  sundials::Context sunctx;

  // Test methods
  int flag, numfails = 0;

  ARKodeButcherTable Be = nullptr;

  int stages, order;
  bool explicit_first_stage, stiffly_accurate, fsal;

  // --------
  // Explicit
  // --------

  std::cout << "\n========================\n"
            << "Test explicit RK methods\n"
            << "========================\n";

  for (int i = ARKODE_MIN_ERK_NUM; i <= ARKODE_MAX_ERK_NUM; i++)
  {
    Be   = ARKodeButcherTable_LoadERK(static_cast<ARKODE_ERKTableID>(i));
    flag = get_method_properties(Be, stages, order, explicit_first_stage,
                                 stiffly_accurate, fsal);
    if (check_flag(&flag, "get_method_properties", 1)) { return 1; }

    std::cout << "\n========================" << std::endl;
    std::cout << ARKodeButcherTable_ERKIDToName(static_cast<ARKODE_ERKTableID>(i))
              << std::endl;
    std::cout << "  stages:             " << stages << std::endl;
    std::cout << "  order:              " << order << std::endl;
    std::cout << "  explicit 1st stage: " << explicit_first_stage << std::endl;
    std::cout << "  stiffly accurate:   " << stiffly_accurate << std::endl;
    std::cout << "  first same as last: " << fsal << std::endl;
    std::cout << "========================" << std::endl;

    numfails += run_tests(Be, prob_data, prob_opts, sunctx);

    ARKodeButcherTable_Free(Be);
    Be = nullptr;
  }

  if (numfails) { std::cout << "\n\nFailed " << numfails << " tests!\n"; }
  else { std::cout << "\n\nAll tests passed!\n"; }

  std::cout << std::endl;
  SUNContext_PrintAllocatorStats(sunctx, stdout, SUN_OUTPUTFORMAT_TABLE);

  // Return test status
  return numfails;
}

// -----------------------------------------------------------------------------
// Test drivers
// -----------------------------------------------------------------------------

int run_tests(ARKodeButcherTable Be, ProblemData& prob_data,
              ProblemOptions& prob_opts, sundials::Context& sunctx)
{
  // Reusable error-checking flag
  int flag;

  // Test failure counter
  int numfails = 0;

  // Get method properties
  int stages, order;
  bool explicit_first_stage, stiffly_accurate, fsal;
  flag = get_method_properties(Be, stages, order, explicit_first_stage,
                               stiffly_accurate, fsal);
  if (check_flag(&flag, "get_method_properties", 1)) { return 1; }

  // Create initial condition vector
  N_Vector y = N_VNew_Serial(1, sunctx);
  if (check_flag((void*)y, "N_VNew_Serial", 0)) { return 1; }

  N_VConst(SUN_RCONST(1.0), y);

  // -----------------
  // Create integrator
  // -----------------

  // Create integrator based on type
  void* erkstep_mem = nullptr;

  erkstep_mem = ERKStepCreate(fe, prob_opts.t0, y, sunctx);
  if (check_flag((void*)erkstep_mem, "ERKStepCreate", 0)) { return 1; }

  // Set user data
  flag = ARKodeSetUserData(erkstep_mem, &prob_data);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }

  // Specify tolerances
  flag = ARKodeSStolerances(erkstep_mem, prob_opts.reltol, prob_opts.abstol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) { return 1; }

  // Specify fixed time step size
  flag = ARKodeSetFixedStep(erkstep_mem, prob_opts.h);
  if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }

  // Lagrange interpolant (removes additional RHS evaluation with DIRK methods)
  if (prob_opts.i_type == interp_type::lagrange)
  {
    flag = ARKodeSetInterpolantType(erkstep_mem, ARK_INTERP_LAGRANGE);
    if (check_flag(&flag, "ARKodeSetInterpolantType", 1)) { return 1; }
  }
  else if (prob_opts.i_type == interp_type::none)
  {
    flag = ARKodeSetInterpolantType(erkstep_mem, ARK_INTERP_NONE);
    if (check_flag(&flag, "ERKodeSetInterpolantType", 1)) { return 1; }
  }

  // Attach Butcher tables
  flag = ERKStepSetTable(erkstep_mem, Be);
  if (check_flag(&flag, "ERKStepSetTables", 1)) { return 1; }

  // --------------
  // Evolve in time
  // --------------

  sunrealtype t_ret = prob_opts.t0;
  sunrealtype t_out = 3 * prob_opts.h;

  long int nfe_expected;

  for (int i = 0; i < 3; i++)
  {
    std::cout << "--------------------" << std::endl;

    // Advance in time
    flag = ARKodeEvolve(erkstep_mem, t_out, y, &t_ret, ARK_ONE_STEP);
    if (check_flag(&flag, "ARKodeEvolve", 1)) { return 1; }

    // Update output time
    t_out += prob_opts.h;

    // Check statistics
    flag = expected_rhs_evals(prob_opts.i_type, stages, explicit_first_stage,
                              stiffly_accurate, fsal, erkstep_mem, nfe_expected);
    if (check_flag(&flag, "expected_rhs_evals", 1)) { return 1; }

    numfails += check_rhs_evals(erkstep_mem, nfe_expected);

    if (numfails)
    {
      std::cout << "Failed " << numfails << " checks\n";
      break;
    }
  }

  // ----------------
  // Get dense output
  // ----------------

  long int extra_fe_evals = 0;

  if (numfails == 0)
  {
    std::cout << "--------------------" << std::endl;
    std::cout << "Dense Output" << std::endl;

    sunrealtype h_last;
    flag = ARKodeGetLastStep(erkstep_mem, &h_last);
    if (check_flag(&flag, "ARKodeGetLastStep", 1)) { return 1; }

    if (prob_opts.i_type != interp_type::none)
    {
      flag = ARKodeGetDky(erkstep_mem, t_ret - h_last / TWO, 0, y);
      if (check_flag(&flag, "ARKodeGetDky", 1)) { return 1; }
    }

    // Stiffly accurate (and FSAL) methods do not require an additional RHS
    // evaluation to get the new RHS value at the end of a step for dense
    // output. However, for methods with an explicit first stage this evaluation
    // can be used at the start of the next step. For methods with an implicit
    // first stage that are not stiffly accurate this evaluation replaces one
    // that would happen at the end of the next step (this is accounted for in
    // expected_rhs_evals after the next step is taken below).
    if (prob_opts.i_type == interp_type::hermite && !stiffly_accurate)
    {
      nfe_expected++;
    }

    // Higher order methods require additional RHS evaluations for dense output
    // with the Hermite interpolant (note default degree is order - 1, except
    // for first order where the degree is 1. These are not accounted for in
    // explicit_rhs_evals and must be carried forward.
    int degree = (order == 1) ? 1 : order - 1;
    if (prob_opts.i_type == interp_type::hermite && degree > 3)
    {
      extra_fe_evals += (degree == 4) ? 1 : 4;
    }

    numfails += check_rhs_evals(erkstep_mem, nfe_expected + extra_fe_evals);

    std::cout << "--------------------" << std::endl;
  }

  // --------------------
  // Additional time step
  // --------------------

  if (numfails == 0)
  {
    // Advance in time
    flag = ARKodeEvolve(erkstep_mem, t_out, y, &t_ret, ARK_ONE_STEP);
    if (check_flag(&flag, "ARKodeEvolve", 1)) { return 1; }

    // Update output time
    t_out += prob_opts.h;

    // Check statistics
    flag = expected_rhs_evals(prob_opts.i_type, stages, explicit_first_stage,
                              stiffly_accurate, fsal, erkstep_mem, nfe_expected);
    if (check_flag(&flag, "expected_rhs_evals", 1)) { return 1; }

    numfails += check_rhs_evals(erkstep_mem, nfe_expected + extra_fe_evals);

    std::cout << "--------------------" << std::endl;
  }

  std::cout << std::endl;
  SUNContext_PrintAllocatorStats(sunctx, stdout, SUN_OUTPUTFORMAT_TABLE);

  // --------
  // Clean up
  // --------

  ARKodeFree(&erkstep_mem);
  N_VDestroy(y);

  return numfails;
}

int get_method_properties(ARKodeButcherTable Be, int& stages, int& order,
                          bool& explicit_first_stage, bool& stiffly_accurate,
                          bool& fsal)
{
  stages = 0;
  if (Be) { stages = Be->stages; }
  else
  {
    std::cerr << "ERROR: Butcher table is NULL!" << std::endl;
    return 1;
  }

  order = 0;
  if (Be) { order = Be->q; }
  else
  {
    std::cerr << "ERROR: Both Butcher tables are NULL!" << std::endl;
    return 1;
  }

  // Check for explicit first stage
  explicit_first_stage = true;
  if (std::abs(Be->A[0][0]) > ZERO) { explicit_first_stage = false; }

  // Check for stiffly accurate method
  stiffly_accurate = ARKodeButcherTable_IsStifflyAccurate(Be);

  // Check for first same as last (FSAL) property
  fsal = explicit_first_stage && stiffly_accurate;

  return 0;
}

int expected_rhs_evals(interp_type i_type, int stages,
                       bool explicit_first_stage, bool stiffly_accurate,
                       bool fsal, void* erkstep_mem, long int& nfe_expected)
{
  int flag = 0;

  // Get number of steps and nonlinear solver iterations
  long int nst = 0;
  flag         = ARKodeGetNumSteps(erkstep_mem, &nst);
  if (check_flag(&flag, "ARKodeGetNumSteps", 1)) { return 1; }

  // Expected number of explicit functions evaluations
  nfe_expected = 0;
  if (fsal)
  {
    // Save one function evaluation after first step
    nfe_expected = stages + (stages - 1) * (nst - 1);
  }
  else { nfe_expected = stages * nst; }

  if (i_type == interp_type::hermite && !explicit_first_stage)
  {
    if (stiffly_accurate)
    {
      // One extra evaluation in the first step only
      nfe_expected++;
    }
    else
    {
      // One extra evaluation in each step
      nfe_expected += nst;
    }
  }

  std::cout << "Steps: " << nst << std::endl;

  return 0;
}

int check_rhs_evals(void* erkstep_mem, long int nfe_expected)
{
  int flag = 0;

  long int nst = 0;
  flag         = ARKodeGetNumSteps(erkstep_mem, &nst);
  if (check_flag(&flag, "ARKodeGetNumSteps", 1)) { return 1; }

  long int nfe;
  flag = ARKodeGetNumRhsEvals(erkstep_mem, 0, &nfe);
  if (check_flag(&flag, "ARKodeGetNumRhsEvals", 1)) { return 1; }

  std::cout << "Fe RHS evals:\n"
            << "  actual:   " << nfe << "\n"
            << "  expected: " << nfe_expected << "\n";

  if (nfe != nfe_expected)
  {
    std::cout << ">>> Check failed <<<" << std::endl;
    return 1;
  }

  return 0;
}

// -----------------------------------------------------------------------------
// Functions called by the solver
// -----------------------------------------------------------------------------

// Explicit ODE RHS function fe(t,y)
int fe(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* y_data    = N_VGetArrayPointer(y);
  sunrealtype* yd_data   = N_VGetArrayPointer(ydot);
  ProblemData* prob_data = static_cast<ProblemData*>(user_data);

  yd_data[0] = prob_data->lambda_e * y_data[0];

  return 0;
}

// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Check function return value
int check_flag(void* flagvalue, const std::string funcname, int opt)
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
