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
 * IMEX Dahlquist problem:
 *
 * 1) y' = lambda_e * y + lambda_i * y
 * 2) M y' = M (lambda_e * y + lambda_i * y)
 * 3) M(t) y' = M(t) (lambda_e * y + lambda_i * y)
 * ---------------------------------------------------------------------------*/

// Header files
#include <arkode/arkode_arkstep.h>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <nvector/nvector_serial.h>
#include <string>
#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>

#include "arkode/arkode.h"
#include "arkode/arkode_butcher.h"

// Constants
#define NEG_ONE SUN_RCONST(-1.0)
#define ZERO    SUN_RCONST(0.0)
#define ONE     SUN_RCONST(1.0)
#define TWO     SUN_RCONST(2.0)

// Method types
enum class mass_matrix_type
{
  identity,
  fixed,
  time_dependent
};

enum class rk_type
{
  expl,
  impl,
  imex
};

enum class interp_type
{
  none = -1,
  hermite,
  lagrange
};

// Problem parameters
struct ProblemData
{
  sunrealtype lambda_e    = NEG_ONE;
  sunrealtype lambda_i    = NEG_ONE;
  mass_matrix_type m_type = mass_matrix_type::identity;
};

// Problem options
struct ProblemOptions
{
  mass_matrix_type m_type = mass_matrix_type::identity;

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

  // Predictor type
  // 0 = trivial
  // 1 = maximum order (dense output)
  int p_type = 0;

  // Runge-Kutta type
  rk_type r_type;
};

// User-supplied Functions called by the solver
int fe(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
int fi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
int Ji(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int MassMatrix(sunrealtype t, SUNMatrix M, void* user_data, N_Vector tmp1,
               N_Vector tmp2, N_Vector tmp3);

// Private function to check function return values
int check_flag(void* flagvalue, const std::string funcname, int opt);

// Test drivers
int run_tests(ARKodeButcherTable Be, ARKodeButcherTable Bi,
              ProblemData& prob_data, ProblemOptions& prob_opts,
              sundials::Context& sunctx);

int get_method_properties(ARKodeButcherTable Be, ARKodeButcherTable Bi,
                          int& stages, int& order, bool& explicit_first_stage,
                          bool& stiffly_accurate, bool& fsal);

int expected_rhs_evals(ProblemOptions& prob_opts, int stages, int order,
                       bool explicit_first_stage, bool stiffly_accurate,
                       bool fsal, void* arkstep_mem, long int& nfe_expected,
                       long int& nfi_expected);

int check_rhs_evals(rk_type r_type, void* arkstep_mem, long int nfe_expected,
                    long int nfi_expected);

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
    if (std::stoi(argv[1]) == 1)
    {
      prob_data.m_type = mass_matrix_type::fixed;
      prob_opts.m_type = mass_matrix_type::fixed;
    }
    else if (std::stoi(argv[1]) == 2)
    {
      prob_data.m_type = mass_matrix_type::time_dependent;
      prob_opts.m_type = mass_matrix_type::time_dependent;
    }
    else
    {
      prob_data.m_type = mass_matrix_type::identity;
      prob_opts.m_type = mass_matrix_type::identity;
    }
  }

  if (argc > 2)
  {
    if (std::stoi(argv[2]) == -1) { prob_opts.i_type = interp_type::none; }
    else if (std::stoi(argv[2]) == 0)
    {
      prob_opts.i_type = interp_type::hermite;
    }
    else if (std::stoi(argv[2]) == 1)
    {
      prob_opts.i_type = interp_type::lagrange;
    }
    else
    {
      std::cerr << "ERROR: Invalid interpolation type option" << std::endl;
      return 1;
    }
  }

  if (argc > 3)
  {
    prob_opts.p_type = std::stoi(argv[3]);
    if (prob_opts.p_type > 1)
    {
      std::cerr << "ERROR: Only the trivial (0) and max order (1) "
                << "predictors are supported" << std::endl;
      return 1;
    }
  }

  // Output problem setup
  std::cout << "\nDahlquist ODE test problem:\n";
  if (prob_data.m_type == mass_matrix_type::identity)
  {
    std::cout << "  problem type = Identity\n";
  }
  else if (prob_data.m_type == mass_matrix_type::fixed)
  {
    std::cout << "  problem type = Fixed mass matrix\n";
  }
  else { std::cout << "  problem type = Time-dependent mass matrix\n"; }
  std::cout << "  lambda expl  = " << prob_data.lambda_e << "\n"
            << "  lambda impl  = " << prob_data.lambda_i << "\n"
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
  if (prob_opts.p_type == 0) { std::cout << "  pred type    = Trivial (0)\n"; }
  else { std::cout << "  pred type    = Max order (1)\n"; }

  // Create SUNDIALS context
  sundials::Context sunctx;

  // Test methods
  int flag, numfails = 0;

  ARKodeButcherTable Be = nullptr;
  ARKodeButcherTable Bi = nullptr;

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
    flag = get_method_properties(Be, Bi, stages, order, explicit_first_stage,
                                 stiffly_accurate, fsal);
    if (check_flag(&flag, "get_method_properties", 1)) { return 1; }

    std::cout << "\n========================" << std::endl;
    std::cout << "ERK: "
              << ARKodeButcherTable_ERKIDToName(static_cast<ARKODE_ERKTableID>(i))
              << std::endl;
    std::cout << "  stages:             " << stages << std::endl;
    std::cout << "  order:              " << order << std::endl;
    std::cout << "  explicit 1st stage: " << explicit_first_stage << std::endl;
    std::cout << "  stiffly accurate:   " << stiffly_accurate << std::endl;
    std::cout << "  first same as last: " << fsal << std::endl;
    std::cout << "========================" << std::endl;

    numfails += run_tests(Be, Bi, prob_data, prob_opts, sunctx);

    ARKodeButcherTable_Free(Be);
    Be = nullptr;
  }

  // --------
  // Implicit
  // --------

  std::cout << "\n========================\n"
            << "Test implicit RK methods\n"
            << "========================\n";

  for (int i = ARKODE_MIN_DIRK_NUM; i <= ARKODE_MAX_DIRK_NUM; i++)
  {
    Bi   = ARKodeButcherTable_LoadDIRK(static_cast<ARKODE_DIRKTableID>(i));
    flag = get_method_properties(Be, Bi, stages, order, explicit_first_stage,
                                 stiffly_accurate, fsal);
    if (check_flag(&flag, "get_method_properties", 1)) { return 1; }

    std::cout << "\n========================" << std::endl;
    std::cout << "DIRK: "
              << ARKodeButcherTable_DIRKIDToName(
                   static_cast<ARKODE_DIRKTableID>(i))
              << std::endl;
    std::cout << "  stages:             " << stages << std::endl;
    std::cout << "  order:              " << order << std::endl;
    std::cout << "  explicit 1st stage: " << explicit_first_stage << std::endl;
    std::cout << "  stiffly accurate:   " << stiffly_accurate << std::endl;
    std::cout << "  first same as last: " << fsal << std::endl;
    std::cout << "========================" << std::endl;

    numfails += run_tests(Be, Bi, prob_data, prob_opts, sunctx);

    ARKodeButcherTable_Free(Bi);
    Bi = nullptr;
  }

  // ----
  // IMEX
  // ----

  std::cout << "\n=====================\n"
            << "Test IMEX ARK methods\n"
            << "=====================\n";

  Be          = ARKodeButcherTable_Alloc(2, SUNFALSE);
  Be->A[1][0] = ONE;
  Be->b[0]    = ONE;
  Be->c[1]    = ONE;
  Be->q       = 1;

  Bi          = ARKodeButcherTable_Alloc(2, SUNFALSE);
  Bi->A[1][1] = ONE;
  Bi->b[1]    = ONE;
  Bi->c[1]    = ONE;
  Bi->q       = 1;

  flag = get_method_properties(Be, Bi, stages, order, explicit_first_stage,
                               stiffly_accurate, fsal);
  if (check_flag(&flag, "get_method_properties", 1)) { return 1; }

  std::cout << "\n========================" << std::endl;
  std::cout << "IMEX Euler" << std::endl;
  std::cout << "  stages:             " << stages << std::endl;
  std::cout << "  order:              " << order << std::endl;
  std::cout << "  explicit 1st stage: " << explicit_first_stage << std::endl;
  std::cout << "  stiffly accurate:   " << stiffly_accurate << std::endl;
  std::cout << "  first same as last: " << fsal << std::endl;
  std::cout << "========================" << std::endl;

  numfails += run_tests(Be, Bi, prob_data, prob_opts, sunctx);

  ARKodeButcherTable_Free(Be);
  ARKodeButcherTable_Free(Bi);
  Be = nullptr;
  Bi = nullptr;

  const char* ark_methods_erk[6] = {"ARKODE_ARK2_ERK_3_1_2",
                                    "ARKODE_ARK324L2SA_ERK_4_2_3",
                                    "ARKODE_ARK436L2SA_ERK_6_3_4",
                                    "ARKODE_ARK437L2SA_ERK_7_3_4",
                                    "ARKODE_ARK548L2SA_ERK_8_4_5",
                                    "ARKODE_ARK548L2SAb_ERK_8_4_5"};

  const char* ark_methods_dirk[6] = {"ARKODE_ARK2_DIRK_3_1_2",
                                     "ARKODE_ARK324L2SA_DIRK_4_2_3",
                                     "ARKODE_ARK436L2SA_DIRK_6_3_4",
                                     "ARKODE_ARK437L2SA_DIRK_7_3_4",
                                     "ARKODE_ARK548L2SA_DIRK_8_4_5",
                                     "ARKODE_ARK548L2SAb_DIRK_8_4_5"};

  for (int i = 0; i < 6; i++)
  {
    Be = ARKodeButcherTable_LoadERKByName(ark_methods_erk[i]);
    Bi = ARKodeButcherTable_LoadDIRKByName(ark_methods_dirk[i]);

    flag = get_method_properties(Be, Bi, stages, order, explicit_first_stage,
                                 stiffly_accurate, fsal);
    if (check_flag(&flag, "get_method_properties", 1)) { return 1; }

    std::cout << "\n========================" << std::endl;
    std::cout << "DIRK: " << ark_methods_dirk[i] << std::endl;
    std::cout << "ERK: " << ark_methods_erk[i] << std::endl;
    std::cout << "  stages:             " << stages << std::endl;
    std::cout << "  order:              " << order << std::endl;
    std::cout << "  explicit 1st stage: " << explicit_first_stage << std::endl;
    std::cout << "  stiffly accurate:   " << stiffly_accurate << std::endl;
    std::cout << "  first same as last: " << fsal << std::endl;
    std::cout << "========================" << std::endl;

    numfails += run_tests(Be, Bi, prob_data, prob_opts, sunctx);

    ARKodeButcherTable_Free(Be);
    ARKodeButcherTable_Free(Bi);
    Be = nullptr;
    Bi = nullptr;
  }

  if (numfails) { std::cout << "\n\nFailed " << numfails << " tests!\n"; }
  else { std::cout << "\n\nAll tests passed!\n"; }

  // Return test status
  return numfails;
}

// -----------------------------------------------------------------------------
// Test drivers
// -----------------------------------------------------------------------------

int run_tests(ARKodeButcherTable Be, ARKodeButcherTable Bi,
              ProblemData& prob_data, ProblemOptions& prob_opts,
              sundials::Context& sunctx)
{
  // Reusable error-checking flag
  int flag;

  // Test failure counter
  int numfails = 0;

  // Determine method type
  if (Be && !Bi) { prob_opts.r_type = rk_type::expl; }
  else if (!Be && Bi) { prob_opts.r_type = rk_type::impl; }
  else if (Be && Bi) { prob_opts.r_type = rk_type::imex; }
  else
  {
    std::cerr << "ERROR: Both Butcher tables are NULL" << std::endl;
    return 1;
  }

  // Get method properties
  int stages, order;
  bool explicit_first_stage, stiffly_accurate, fsal;
  flag = get_method_properties(Be, Bi, stages, order, explicit_first_stage,
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
  void* arkstep_mem = nullptr;

  if (prob_opts.r_type == rk_type::expl)
  {
    arkstep_mem = ARKStepCreate(fe, nullptr, prob_opts.t0, y, sunctx);
  }
  else if (prob_opts.r_type == rk_type::impl)
  {
    arkstep_mem = ARKStepCreate(nullptr, fi, prob_opts.t0, y, sunctx);
  }
  else { arkstep_mem = ARKStepCreate(fe, fi, prob_opts.t0, y, sunctx); }
  if (check_flag((void*)arkstep_mem, "ARKStepCreate", 0)) { return 1; }

  // Set user data
  flag = ARKodeSetUserData(arkstep_mem, &prob_data);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }

  // Specify tolerances
  flag = ARKodeSStolerances(arkstep_mem, prob_opts.reltol, prob_opts.abstol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) { return 1; }

  // Specify fixed time step size
  flag = ARKodeSetFixedStep(arkstep_mem, prob_opts.h);
  if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }

  // Attach Butcher tables (ignore actual method order)
  flag = ARKStepSetTables(arkstep_mem, 1, 0, Bi, Be);
  if (check_flag(&flag, "ARKStepSetTables", 1)) { return 1; }

  // Lagrange interpolant (removes additional RHS evaluation with DIRK methods)
  if (prob_opts.i_type == interp_type::lagrange)
  {
    flag = ARKodeSetInterpolantType(arkstep_mem, ARK_INTERP_LAGRANGE);
    if (check_flag(&flag, "ARKodeSetInterpolantType", 1)) { return 1; }
  }
  else if (prob_opts.i_type == interp_type::none)
  {
    flag = ARKodeSetInterpolantType(arkstep_mem, ARK_INTERP_NONE);
    if (check_flag(&flag, "ARKodeSetInterpolantType", 1)) { return 1; }
  }

  // Create matrix and linear solver (if necessary)
  SUNMatrix A        = nullptr;
  SUNLinearSolver LS = nullptr;

  if (prob_opts.r_type == rk_type::impl || prob_opts.r_type == rk_type::imex)
  {
    // Initialize dense matrix data structures and solvers
    A = SUNDenseMatrix(1, 1, sunctx);
    if (check_flag((void*)A, "SUNDenseMatrix", 0)) { return 1; }

    LS = SUNLinSol_Dense(y, A, sunctx);
    if (check_flag((void*)LS, "SUNLinSol_Dense", 0)) { return 1; }

    // Attach linear solver
    flag = ARKodeSetLinearSolver(arkstep_mem, LS, A);
    if (check_flag(&flag, "ARKodeSetLinearSolver", 1)) { return 1; }

    // Set Jacobian function
    flag = ARKodeSetJacFn(arkstep_mem, Ji);
    if (check_flag(&flag, "ARKodeSetJacFn", 1)) { return 1; }

    // Specify linearly implicit RHS, with non-time-dependent Jacobian
    flag = ARKodeSetLinear(arkstep_mem, 0);
    if (check_flag(&flag, "ARKodeSetLinear", 1)) { return 1; }

    // Specify an autonomous RHS
    flag = ARKodeSetAutonomous(arkstep_mem, SUNTRUE);
    if (check_flag(&flag, "ARKodeSetAutonomous", 1)) { return 1; }

    // Specify implicit predictor method
    flag = ARKodeSetPredictorMethod(arkstep_mem, prob_opts.p_type);
    if (check_flag(&flag, "ARKodeSetPredictorMethod", 1)) { return 1; }
  }

  // Create mass matrix and linear solver (if necessary)
  SUNMatrix M         = nullptr;
  SUNLinearSolver MLS = nullptr;

  if (prob_data.m_type == mass_matrix_type::fixed ||
      prob_data.m_type == mass_matrix_type::time_dependent)
  {
    M = SUNDenseMatrix(1, 1, sunctx);
    if (check_flag((void*)M, "SUNDenseMatrix", 0)) { return 1; }

    MLS = SUNLinSol_Dense(y, M, sunctx);
    if (check_flag((void*)MLS, "SUNLinSol_Dense", 0)) { return 1; }

    int time_dep = 0;
    if (prob_data.m_type == mass_matrix_type::time_dependent) { time_dep = 1; }

    flag = ARKodeSetMassLinearSolver(arkstep_mem, MLS, M, time_dep);
    if (check_flag(&flag, "ARKodeSetMassLinearSolver", 1)) { return 1; }

    flag = ARKodeSetMassFn(arkstep_mem, MassMatrix);
    if (check_flag(&flag, "ARKodeSetMassFn", 1)) { return 1; }
  }

  // --------------
  // Evolve in time
  // --------------

  sunrealtype t_ret = prob_opts.t0;
  sunrealtype t_out = 3 * prob_opts.h;

  long int nfe_expected, nfi_expected;

  for (int i = 0; i < 3; i++)
  {
    std::cout << "--------------------" << std::endl;

    // Advance in time
    flag = ARKodeEvolve(arkstep_mem, t_out, y, &t_ret, ARK_ONE_STEP);
    if (check_flag(&flag, "ARKodeEvolve", 1)) { return 1; }

    // Update output time
    t_out += prob_opts.h;

    // Check statistics
    flag = expected_rhs_evals(prob_opts, stages, order, explicit_first_stage,
                              stiffly_accurate, fsal, arkstep_mem, nfe_expected,
                              nfi_expected);
    if (check_flag(&flag, "expected_rhs_evals", 1)) { return 1; }

    numfails += check_rhs_evals(prob_opts.r_type, arkstep_mem, nfe_expected,
                                nfi_expected);

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
  long int extra_fi_evals = 0;

  if (numfails == 0)
  {
    std::cout << "--------------------" << std::endl;
    std::cout << "Dense Output" << std::endl;

    sunrealtype h_last;
    flag = ARKodeGetLastStep(arkstep_mem, &h_last);
    if (check_flag(&flag, "ARKodeGetLastStep", 1)) { return 1; }

    if (prob_opts.i_type != interp_type::none)
    {
      flag = ARKodeGetDky(arkstep_mem, t_ret - h_last / TWO, 0, y);
      if (check_flag(&flag, "ARKodeGetDky", 1)) { return 1; }
    }

    // Stiffly accurate (and FSAL) methods do not require an additional RHS
    // evaluation to get the new RHS value at the end of a step for dense
    // output. Methods with an explicit first stage can use this RHS value at
    // the start of the next step.
    //
    // If a method has an implicit first stage but is not stiffly accurate, the
    // evaluation needed for interpolation replaces one that would happen at the
    // end of the next step (this is accounted for in expected_rhs_evals after
    // the next step is taken below). However, if the trivial predictor is used
    // then the RHS evaluation can be reused across the stage solves in the
    // initial evaluation of the residual.
    if (prob_opts.i_type == interp_type::hermite && !stiffly_accurate)
    {
      if (prob_opts.r_type == rk_type::expl || prob_opts.r_type == rk_type::imex)
      {
        nfe_expected++;
      }
      if (prob_opts.r_type == rk_type::impl || prob_opts.r_type == rk_type::imex)
      {
        nfi_expected++;
      }
    }

    // Higher order methods require additional RHS evaluations for dense output
    // with the Hermite interpolant (note default degree is order - 1, except
    // for first order where the degree is 1. These are not accounted for in
    // explicit_rhs_evals and must be carried forward.
    int degree = (order == 1) ? 1 : order - 1;
    if (prob_opts.i_type == interp_type::hermite && degree > 3)
    {
      if (prob_opts.r_type == rk_type::expl || prob_opts.r_type == rk_type::imex)
      {
        extra_fe_evals += (degree == 4) ? 1 : 4;
      }
      if (prob_opts.r_type == rk_type::impl || prob_opts.r_type == rk_type::imex)
      {
        extra_fi_evals += (degree == 4) ? 1 : 4;
      }
    }

    numfails += check_rhs_evals(prob_opts.r_type, arkstep_mem,
                                nfe_expected + extra_fe_evals,
                                nfi_expected + extra_fi_evals);

    std::cout << "--------------------" << std::endl;
  }

  // --------------------
  // Additional time step
  // --------------------

  if (numfails == 0)
  {
    // Advance in time
    flag = ARKodeEvolve(arkstep_mem, t_out, y, &t_ret, ARK_ONE_STEP);
    if (check_flag(&flag, "ARKodeEvolve", 1)) { return 1; }

    // Update output time
    t_out += prob_opts.h;

    // Check statistics
    flag = expected_rhs_evals(prob_opts, stages, order, explicit_first_stage,
                              stiffly_accurate, fsal, arkstep_mem, nfe_expected,
                              nfi_expected);
    if (check_flag(&flag, "expected_rhs_evals", 1)) { return 1; }

    numfails += check_rhs_evals(prob_opts.r_type, arkstep_mem,
                                nfe_expected + extra_fe_evals,
                                nfi_expected + extra_fi_evals);

    std::cout << "--------------------" << std::endl;
  }

  // --------
  // Clean up
  // --------

  ARKodeFree(&arkstep_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  SUNLinSolFree(MLS);
  SUNMatDestroy(M);
  N_VDestroy(y);

  return numfails;
}

int get_method_properties(ARKodeButcherTable Be, ARKodeButcherTable Bi,
                          int& stages, int& order, bool& explicit_first_stage,
                          bool& stiffly_accurate, bool& fsal)
{
  stages = 0;
  if (Bi) { stages = Bi->stages; }
  else if (Be) { stages = Be->stages; }
  else
  {
    std::cerr << "ERROR: Both Butcher tables are NULL!" << std::endl;
    return 1;
  }

  // Built-in ARK methods have the same order for Bi and Be
  order = 0;
  if (Bi) { order = Bi->q; }
  else if (Be) { order = Be->q; }
  else
  {
    std::cerr << "ERROR: Both Butcher tables are NULL!" << std::endl;
    return 1;
  }

  // Check for explicit first stage
  explicit_first_stage = true;
  if (Bi)
  {
    if (std::abs(Bi->A[0][0]) > ZERO) { explicit_first_stage = false; }
  }
  if (Be)
  {
    if (std::abs(Be->A[0][0]) > ZERO) { explicit_first_stage = false; }
  }

  // Check for stiffly accurate method
  stiffly_accurate = true;
  if (Bi)
  {
    if (!ARKodeButcherTable_IsStifflyAccurate(Bi)) { stiffly_accurate = false; }
  }
  if (Be)
  {
    if (!ARKodeButcherTable_IsStifflyAccurate(Be)) { stiffly_accurate = false; }
  }

  // Check for first same as last (FSAL) property
  fsal = explicit_first_stage && stiffly_accurate;

  return 0;
}

int expected_rhs_evals(ProblemOptions& prob_opts, int stages, int order,
                       bool explicit_first_stage, bool stiffly_accurate,
                       bool fsal, void* arkstep_mem, long int& nfe_expected,
                       long int& nfi_expected)
{
  int flag = 0;

  // Get number of steps and nonlinear solver iterations
  long int nst = 0;

  flag = ARKodeGetNumSteps(arkstep_mem, &nst);
  if (check_flag(&flag, "ARKodeGetNumSteps", 1)) { return 1; }

  long int nni            = 0;
  long int extra_fe_evals = 0;
  long int extra_fi_evals = 0;

  bool save_fn_for_residual = prob_opts.p_type == 0 &&
                              prob_opts.m_type !=
                                mass_matrix_type::time_dependent &&
                              prob_opts.r_type != rk_type::expl;

  if (prob_opts.r_type == rk_type::impl || prob_opts.r_type == rk_type::imex)
  {
    flag = ARKodeGetNumNonlinSolvIters(arkstep_mem, &nni);
    if (check_flag(&flag, "ARKodeGetNumNonlinSolvIters", 1)) { return 1; }

    // Remove one evaluation per implicit stage
    int implicit_stages = (explicit_first_stage) ? stages - 1 : stages;

    if (save_fn_for_residual) { extra_fi_evals -= implicit_stages * nst; }

    // With higher order methods some predictors require additional RHS when
    // using Hermite interpolation (note default degree is order - 1, except
    // for first order where the degree is 1.
    int degree = (order == 1) ? 1 : order - 1;

    if (prob_opts.p_type != 0 && prob_opts.i_type == interp_type::hermite &&
        degree > 3)
    {
      if (prob_opts.r_type == rk_type::expl || prob_opts.r_type == rk_type::imex)
      {
        extra_fe_evals = (degree == 4) ? 1 : 4;
      }
      if (prob_opts.r_type == rk_type::impl || prob_opts.r_type == rk_type::imex)
      {
        extra_fi_evals = (degree == 4) ? 1 : 4;
      }
      extra_fe_evals *= implicit_stages * (nst - 1);
      extra_fi_evals *= implicit_stages * (nst - 1);
    }
  }

  // Expected number of explicit functions evaluations
  nfe_expected = extra_fe_evals;

  if (prob_opts.r_type == rk_type::expl || prob_opts.r_type == rk_type::imex)
  {
    if (fsal)
    {
      // Save one function evaluation after first step
      nfe_expected += stages + (stages - 1) * (nst - 1);
    }
    else { nfe_expected += stages * nst; }

    if (prob_opts.i_type == interp_type::hermite && !explicit_first_stage)
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
  }

  // Expected number of implicit functions evaluations
  nfi_expected = extra_fi_evals;

  if (prob_opts.r_type == rk_type::impl || prob_opts.r_type == rk_type::imex)
  {
    if (fsal)
    {
      // Save one function evaluation after first step
      nfi_expected += stages + (stages - 1) * (nst - 1) + nni;
    }
    else { nfi_expected += stages * nst + nni; }

    if (prob_opts.i_type == interp_type::hermite && !explicit_first_stage)
    {
      if (stiffly_accurate)
      {
        // One extra evaluation in the first step only
        nfi_expected++;
      }
      else
      {
        // One extra evaluation in each step
        nfi_expected += nst;
      }
    }

    if (prob_opts.i_type != interp_type::hermite && save_fn_for_residual &&
        !explicit_first_stage)
    {
      if (stiffly_accurate) { nfi_expected++; }
      else { nfi_expected += nst; }
    }
  }

  std::cout << "Steps: " << nst << std::endl;

  if (prob_opts.r_type == rk_type::impl || prob_opts.r_type == rk_type::imex)
  {
    std::cout << "NLS iters: " << nni << std::endl;
  }

  return 0;
}

int check_rhs_evals(rk_type r_type, void* arkstep_mem, long int nfe_expected,
                    long int nfi_expected)
{
  int flag = 0;

  long int nst;
  flag = ARKodeGetNumSteps(arkstep_mem, &nst);
  if (check_flag(&flag, "ARKodeGetNumSteps", 1)) { return 1; }

  long int nfe, nfi;
  flag = ARKodeGetNumRhsEvals(arkstep_mem, 0, &nfe);
  if (check_flag(&flag, "ARKodeGetNumRhsEvals", 1)) { return 1; }

  flag = ARKodeGetNumRhsEvals(arkstep_mem, 1, &nfi);
  if (check_flag(&flag, "ARKodeGetNumRhsEvals", 1)) { return 1; }

  if (r_type == rk_type::expl || r_type == rk_type::imex)
  {
    std::cout << "Fe RHS evals:\n"
              << "  actual:   " << nfe << "\n"
              << "  expected: " << nfe_expected << "\n";
  }
  if (r_type == rk_type::impl || r_type == rk_type::imex)
  {
    std::cout << "Fi RHS evals:\n"
              << "  actual:   " << nfi << "\n"
              << "  expected: " << nfi_expected << "\n";
  }

  if (nfe != nfe_expected || nfi != nfi_expected)
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

  if (prob_data->m_type == mass_matrix_type::fixed) { yd_data[0] *= TWO; }
  else if (prob_data->m_type == mass_matrix_type::time_dependent)
  {
    yd_data[0] *= TWO + std::cos(t);
  }

  return 0;
}

// Implicit ODE RHS function fi(t,y)
int fi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* y_data    = N_VGetArrayPointer(y);
  sunrealtype* yd_data   = N_VGetArrayPointer(ydot);
  ProblemData* prob_data = static_cast<ProblemData*>(user_data);

  yd_data[0] = prob_data->lambda_i * y_data[0];

  if (prob_data->m_type == mass_matrix_type::fixed) { yd_data[0] *= TWO; }
  else if (prob_data->m_type == mass_matrix_type::time_dependent)
  {
    yd_data[0] *= TWO + std::cos(t);
  }

  return 0;
}

// Jacobian routine to compute J(t,y) = dfi/dy.
int Ji(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* J_data    = SUNDenseMatrix_Data(J);
  ProblemData* prob_data = static_cast<ProblemData*>(user_data);

  J_data[0] = prob_data->lambda_i;

  return 0;
}

int MassMatrix(sunrealtype t, SUNMatrix M, void* user_data, N_Vector tmp1,
               N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* M_data    = SUNDenseMatrix_Data(M);
  ProblemData* prob_data = static_cast<ProblemData*>(user_data);

  if (prob_data->m_type == mass_matrix_type::fixed) { M_data[0] = TWO; }
  else { M_data[0] = TWO + std::cos(t); }

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
    if (*errflag != 0)
    {
      std::cerr << "\nSUNDIALS_ERROR: " << funcname
                << " failed with flag = " << *errflag << "\n\n";
      return 1;
    }
  }

  return 0;
}
