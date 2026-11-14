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
 * This example simulates the 1D advection-diffusion-reaction equation,
 *
 *   u_t = -c u_x + d u_xx + A - (w + 1) * u + v * u^2
 *   v_t = -c v_x + d u_xx + w * u - v * u^2
 *   w_t = -c w_x + d w_xx + (B - w) / eps - w * u
 *
 * where u, v, and w represent the concentrations of chemical species, c = 0.001
 * is the advection speed, d = 0.01 is the diffusion rate, and the species with
 * constant concentration over time are A = 0.6 and B = 2.0.
 *
 * The problem is evolved for t in [0, 3] and x in [0, 1], with initial
 * conditions given by
 *
 *   u(0,x) =  A  + 0.1 * sin(pi * x)
 *   v(0,x) = B/A + 0.1 * sin(pi * x)
 *   w(0,x) =  B  + 0.1 * sin(pi * x)
 *
 * and stationary boundary conditions i.e.,
 *
 *   u_t(t,0) = u_t(t,1) = 0,
 *   v_t(t,0) = v_t(t,1) = 0,
 *   w_t(t,0) = w_t(t,1) = 0.
 *
 * The system is advanced in time using one of the following approaches based on
 * the --integrator <int> flag value. The following options are available:
 *
 *   0. An explicit Runge-Kutta method with ERKStep.
 *
 *   1. An explicit, diagonally implicit, or IMEX Runge-Kutta method with
 *      ARKStep. The method used depends on the value of the --splitting <int>
 *      flag (denoted by IDX in the tables below).
 *
 *      Advection-Diffusion-Reaction splittings
 *
 *      | IDX | A | D | R | Description                                     |
 *      +-----+---+---+---+-------------------------------------------------+
 *      |  0  | E | E | E | fully explicit                                  |
 *      |  1  | E | E | I | implicit reaction, explicit advection-diffusion |
 *      |  2  | E | I | E | implicit diffusion, explicit advection-reaction |
 *      |  3  | E | I | I | implicit diffusion-reaction, explicit advection |
 *      |  4  | I | E | E | implicit advection, explicit diffusion-reaction |
 *      |  5  | I | E | I | implicit advection-reaction, explicit diffusion |
 *      |  6  | I | I | E | implicit advection-diffusion, explicit reaction |
 *      |  7  | I | I | I | fully implicit                                  |
 *      +-----+---+---+---+-------------------------------------------------+
 *
 *      Advection-Reaction splittings (use the --no-diffusion flag)
 *
 *      | IDX | A | R | Description                           |
 *      +-----+---+---+---------------------------------------+
 *      |  0  | E | E | fully explicit                        |
 *      |  1  | E | I | implicit reaction, explicit advection |
 *      |  2  | I | E | implicit advection, explicit reaction |
 *      |  3  | I | I | fully implicit                        |
 *      +-----+---+---+---------------------------------------+
 *
 *      Diffusion-Reaction splittings (use the --no-advection flag)
 *
 *      | IDX | D | R | Description                           |
 *      +-----+---+---+---------------------------------------+
 *      |  0  | E | E | fully explicit                        |
 *      |  1  | E | I | explicit diffusion, implicit reaction |
 *      |  2  | I | E | explicit reaction, implicit diffusion |
 *      |  3  | I | I | fully implicit                        |
 *      +-----+---+---+---------------------------------------+
 *
 *   2. An explicit, implicit-solve-decoupled, or IMEX-solve-decoupled MRI-GARK
 *      method with MRIStep. Advection is treated explicitly at the slow time
 *      scale, diffusion implicitly at the slow time scale, and the reaction
 *      terms are integrated using an explicit or implicit method from ARKStep
 *      at the fast time scale depending on the value of the --splitting <int>
 *      flag (denoted by IDX in the tables below).
 *
 *      Advection-Diffusion-Reaction splittings
 *
 *      | IDX |  A  |  D  |  R  | Description                                  |
 *      +-----+-----+-----+-----+----------------------------------------------+
 *      |  0  | S-E | S-I | F-E | Slow-explicit advection,                     |
 *      |     |     |     |     | Slow-Implicit diffusion,                     |
 *      |     |     |     |     | Fast-explicit reaction                       |
 *      |  1  | S-E | S-I | F-I | Slow-explicit advection,                     |
 *      |     |     |     |     | Slow-Implicit diffusion,                     |
 *      |     |     |     |     | Fast-Implicit reaction                       |
 *      +-----+-----+-----+-----+----------------------------------------------+
 *
 *      Advection-Reaction splittings (use the --no-diffusion flag)
 *
 *      | IDX |  A  |  R  | Description                                     |
 *      +-----+-----+-----+-------------------------------------------------+
 *      |  0  | S-E | F-E | Slow-explicit advection, Fast-explicit reaction |
 *      |  1  | S-E | F-I | Slow-explicit advection, Fast-Implicit reaction |
 *      +-----+-----+-----+-------------------------------------------------+
 *
 *      Diffusion-Reaction splittings (use the --no-advection flag)
 *
 *      | IDX |  D  |  R  | Description                                     |
 *      +-----+-----+-----+-------------------------------------------------+
 *      |  0  | S-I | F-E | Slow-implicit diffusion, Fast-explicit reaction |
 *      |  1  | S-I | F-I | Slow-implicit diffusion, Fast-Implicit reaction |
 *      +-----+-----+-----+-------------------------------------------------+
 *
 *   3. An explicit, implicit-solve-decoupled, or IMEX-solve-decoupled MRI-GARK
 *      method with MRIStep. Advection is treated explicitly at the slow time
 *      scale, diffusion implicitly at the slow time scale, and the reaction
 *      terms are integrated implicitly using a custom inner stepper wrapping
 *      CVODE,
 *
 *      Advection-Diffusion-Reaction splitting
 *
 *      |  A  |  D  |  R  | Description                                  |
 *      +-----+-----+-----+----------------------------------------------+
 *      | S-E | S-I | F-I | Slow-explicit advection,                     |
 *      |     |     |     | Slow-Implicit diffusion,                     |
 *      |     |     |     | Fast-Implicit reaction                       |
 *      +-----+-----+-----+----------------------------------------------+
 *
 *      Advection-Reaction splitting (use the --no-diffusion flag)
 *
 *      |  A  |  R  | Description                                     |
 *      +-----+-----+-------------------------------------------------+
 *      | S-E | F-I | Slow-explicit advection, Fast-Implicit reaction |
 *      +-----+-----+-------------------------------------------------+
 *
 *      Diffusion-Reaction splitting (use the --no-advection flag)
 *
 *      |  D  |  R  | Description                                     |
 *      +-----+-----+-------------------------------------------------+
 *      | S-I | F-I | Slow-implicit diffusion, Fast-Implicit reaction |
 *      +-----+-----+-------------------------------------------------+
 *
 * Several command line options are available to change the problem parameters
 * and integrator settings. Use the flag --help for more information.
 * ---------------------------------------------------------------------------*/

#include "ark_advection_diffusion_reaction.hpp"

int main(int argc, char* argv[])
{
  // SUNDIALS context object for this simulation
  sundials::Context ctx;

  // -----------------
  // Setup the problem
  // -----------------

  UserData udata;
  UserOptions uopts;

  vector<string> args(argv + 1, argv + argc);

  int flag = ReadInputs(args, udata, uopts, ctx);
  if (check_flag(flag, "ReadInputs")) { return 1; }

  flag = PrintSetup(udata, uopts);
  if (check_flag(flag, "PrintSetup")) { return 1; }

  // Create state vector and set initial condition
  N_Vector y = N_VNew_Serial(udata.neq, ctx);
  if (check_ptr(y, "N_VNew_Serial")) { return 1; }

  flag = SetIC(y, udata);
  if (check_flag(flag, "SetIC")) { return 1; }

  // --------------------
  // Setup the integrator
  // --------------------

  // ARKODE memory structure
  void* arkode_mem = nullptr;

  // Matrix and linear solver for DIRK, IMEX, or MRI slow integrators
  SUNMatrix A        = nullptr;
  SUNLinearSolver LS = nullptr;

  // Matrix and linear solver for MRI fast integrator
  SUNMatrix A_fast        = nullptr;
  SUNLinearSolver LS_fast = nullptr;

  // Fast integrator for MRIStep
  MRIStepInnerStepper fast_mem = nullptr;

  // Create integrator
  switch (uopts.integrator)
  {
  case (0): flag = SetupERK(ctx, udata, uopts, y, &arkode_mem); break;
  case (1): flag = SetupARK(ctx, udata, uopts, y, &A, &LS, &arkode_mem); break;
  case (2):
    flag = SetupMRIARK(ctx, udata, uopts, y, &A, &LS, &A_fast, &LS_fast,
                       &fast_mem, &arkode_mem);
    break;
  case (3):
    flag = SetupMRICVODE(ctx, udata, uopts, y, &A, &LS, &A_fast, &LS_fast,
                         &fast_mem, &arkode_mem);
    break;
  default: flag = -1;
  }
  if (check_flag(flag, "Integrator setup")) { return 1; }

  // ----------------------
  // Evolve problem in time
  // ----------------------

  // Initial time, time between outputs, output time
  sunrealtype t     = ZERO;
  sunrealtype dTout = udata.tf / uopts.nout;
  sunrealtype tout  = dTout;

  // Initial output
  flag = OpenOutput(udata, uopts);
  if (check_flag(flag, "OpenOutput")) { return 1; }

  flag = WriteOutput(t, y, udata, uopts);
  if (check_flag(flag, "WriteOutput")) { return 1; }

  // Loop over output times
  for (int iout = 0; iout < uopts.nout; iout++)
  {
    // Evolve
    if (uopts.output == 3)
    {
      // Stop at output time (do not interpolate output)
      flag = ARKodeSetStopTime(arkode_mem, tout);
      if (check_flag(flag, "ARKodeSetStopTime")) { return 1; }
    }

    //   Advance in time
    flag = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_flag(flag, "ARKodeEvolve")) { break; }

    // Output solution
    flag = WriteOutput(t, y, udata, uopts);
    if (check_flag(flag, "WriteOutput")) { return 1; }

    // Update output time
    tout += dTout;
    tout = (tout > udata.tf) ? udata.tf : tout;
  }

  // Close output
  flag = CloseOutput(uopts);
  if (check_flag(flag, "CloseOutput")) { return 1; }

  // ------------
  // Output stats
  // ------------

  if (uopts.output)
  {
    cout << "Final integrator statistics:" << endl;
    switch (uopts.integrator)
    {
    case (0): flag = OutputStatsERK(arkode_mem, udata); break;
    case (1): flag = OutputStatsARK(arkode_mem, udata); break;
    case (2): flag = OutputStatsMRIARK(arkode_mem, fast_mem, udata); break;
    case (3): flag = OutputStatsMRICVODE(arkode_mem, fast_mem, udata); break;
    default: flag = -1;
    }
    if (check_flag(flag, "OutputStats")) { return 1; }
  }

  // --------
  // Clean up
  // --------

  switch (uopts.integrator)
  {
  case (0): ARKodeFree(&arkode_mem); break;
  case (1): ARKodeFree(&arkode_mem); break;
  case (2):
  {
    void* inner_arkode_mem = nullptr;
    MRIStepInnerStepper_GetContent(fast_mem, &inner_arkode_mem);
    ARKodeFree(&inner_arkode_mem);
    MRIStepInnerStepper_Free(&fast_mem);
    ARKodeFree(&arkode_mem);
    break;
  }
  case (3):
  {
    void* inner_content = nullptr;
    MRIStepInnerStepper_GetContent(fast_mem, &inner_content);
    CVodeInnerStepperContent* content = (CVodeInnerStepperContent*)inner_content;
    CVodeFree(&(content->cvode_mem));
    delete content;
    MRIStepInnerStepper_Free(&fast_mem);
    ARKodeFree(&arkode_mem);
    break;
  }
  }

  N_VDestroy(y);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  SUNMatDestroy(A_fast);
  SUNLinSolFree(LS_fast);

  return 0;
}

// -----------------------------------------------------------------------------
// Setup the integrator
// -----------------------------------------------------------------------------

int SetupERK(SUNContext ctx, UserData& udata, UserOptions& uopts, N_Vector y,
             void** arkode_mem)
{
  // Problem configuration
  ARKRhsFn f_RHS; // explicit RHS function

  if (udata.diffusion && udata.advection)
  {
    // Explicit -- advection-diffusion-reaction
    f_RHS = f_adv_diff_react;
  }
  else if (!udata.diffusion && udata.advection)
  {
    // Explicit -- advection-reaction
    f_RHS = f_adv_react;
  }
  else if (udata.diffusion && !udata.advection)
  {
    // Explicit -- diffusion-reaction
    f_RHS = f_diff_react;
  }
  else
  {
    cerr << "ERROR: Invalid problem configuration" << endl;
    return -1;
  }

  // Create ERKStep memory
  *arkode_mem = ERKStepCreate(f_RHS, ZERO, y, ctx);
  if (check_ptr(arkode_mem, "ERKStepCreate")) { return 1; }

  // Specify tolerances
  int flag = ARKodeSStolerances(*arkode_mem, uopts.rtol, uopts.atol);
  if (check_flag(flag, "ARKodeSStolerances")) { return 1; }

  // Attach user data
  flag = ARKodeSetUserData(*arkode_mem, &udata);
  if (check_flag(flag, "ARKodeSetUserData")) { return 1; }

  // Select method order
  flag = ARKodeSetOrder(*arkode_mem, uopts.order);
  if (check_flag(flag, "ARKodeSetOrder")) { return 1; }

  // Set fixed step size or adaptivity method
  if (uopts.fixed_h > ZERO)
  {
    flag = ARKodeSetFixedStep(*arkode_mem, uopts.fixed_h);
    if (check_flag(flag, "ARKodeSetFixedStep")) { return 1; }
  }
  else
  {
    flag = ARKodeSetAdaptControllerByName(*arkode_mem, uopts.controller.c_str());
    if (check_flag(flag, "ARKodeSetAdaptControllerByName")) { return 1; }
  }

  // Set max steps between outputs
  flag = ARKodeSetMaxNumSteps(*arkode_mem, uopts.maxsteps);
  if (check_flag(flag, "ARKodeSetMaxNumSteps")) { return 1; }

  // Set stopping time
  flag = ARKodeSetStopTime(*arkode_mem, udata.tf);
  if (check_flag(flag, "ARKodeSetStopTime")) { return 1; }

  return 0;
}

int SetupARK(SUNContext ctx, UserData& udata, UserOptions& uopts, N_Vector y,
             SUNMatrix* A, SUNLinearSolver* LS, void** arkode_mem)
{
  // Problem configuration
  ARKRhsFn fe_RHS;   // explicit RHS function
  ARKRhsFn fi_RHS;   // implicit RHS function
  ARKLsJacFn Ji_RHS; // Jacobian of implicit RHS function

  // advection-diffusion-reaction
  if (udata.diffusion && udata.advection)
  {
    switch (udata.splitting)
    {
    case (0):
      // ERK -- fully explicit
      fe_RHS = f_adv_diff_react;
      fi_RHS = nullptr;
      Ji_RHS = nullptr;
      break;
    case (1):
      // IMEX -- explicit advection-diffusion, implicit reaction
      fe_RHS = f_adv_diff;
      fi_RHS = f_reaction;
      Ji_RHS = J_reaction;
      break;
    case (2):
      // IMEX -- explicit advection-reaction, implicit diffusion
      fe_RHS = f_adv_react;
      fi_RHS = f_diffusion;
      Ji_RHS = J_diffusion;
      break;
    case (3):
      // IMEX -- explicit advection, implicit diffusion-reaction
      fe_RHS = f_advection;
      fi_RHS = f_diff_react;
      Ji_RHS = J_diff_react;
      break;
    case (4):
      // IMEX -- explicit diffusion-reaction, implicit advection
      fe_RHS = f_diff_react;
      fi_RHS = f_advection;
      Ji_RHS = J_advection;
      break;
    case (5):
      // IMEX -- explicit diffusion, implicit advection-reaction
      fe_RHS = f_diffusion;
      fi_RHS = f_adv_react;
      Ji_RHS = J_adv_react;
      break;
    case (6):
      // IMEX -- explicit reaction, implicit advection-diffusion
      fe_RHS = f_reaction;
      fi_RHS = f_adv_diff;
      Ji_RHS = J_adv_diff;
      break;
    case (7):
      // DIRK -- fully implicit
      fe_RHS = nullptr;
      fi_RHS = f_adv_diff_react;
      Ji_RHS = J_adv_diff_react;
      break;
    default:
      cerr << "ERROR: Invalid splitting option" << endl;
      return -1;
      break;
    }
  }
  // advection-reaction
  else if (!udata.diffusion && udata.advection)
  {
    switch (udata.splitting)
    {
    case (0):
      // ERK -- fully explicit
      fe_RHS = f_adv_react;
      fi_RHS = nullptr;
      Ji_RHS = nullptr;
      break;
    case (1):
      // IMEX -- explicit advection, implicit reaction
      fe_RHS = f_advection;
      fi_RHS = f_reaction;
      Ji_RHS = J_reaction;
      break;
    case (2):
      // IMEX -- explicit reaction, implicit advection
      fe_RHS = f_reaction;
      fi_RHS = f_advection;
      Ji_RHS = J_advection;
      break;
    case (3):
      // DIRK -- fully implicit
      fe_RHS = nullptr;
      fi_RHS = f_adv_react;
      Ji_RHS = J_adv_react;
      break;
    default:
      cerr << "ERROR: Invalid splitting option" << endl;
      return -1;
      break;
    }
  }
  // diffusion-reaction
  else if (udata.diffusion && !udata.advection)
  {
    switch (udata.splitting)
    {
    case (0):
      // ERK -- fully explicit
      fe_RHS = f_diff_react;
      fi_RHS = nullptr;
      Ji_RHS = nullptr;
      break;
    case (1):
      // IMEX -- explicit diffusion, implicit reaction
      fe_RHS = f_diffusion;
      fi_RHS = f_reaction;
      Ji_RHS = J_reaction;
      break;
    case (2):
      // IMEX -- explicit reaction, implicit diffusion
      fe_RHS = f_reaction;
      fi_RHS = f_diffusion;
      Ji_RHS = J_diffusion;
      break;
    case (4):
      // DIRK -- fully implicit
      fe_RHS = nullptr;
      fi_RHS = f_diff_react;
      Ji_RHS = J_diff_react;
      break;
    default:
      cerr << "ERROR: Invalid splitting option" << endl;
      return -1;
      break;
    }
  }
  else
  {
    cerr << "ERROR: Invalid problem configuration" << endl;
    return -1;
  }

  // Create ARKStep memory
  *arkode_mem = ARKStepCreate(fe_RHS, fi_RHS, ZERO, y, ctx);
  if (check_ptr(arkode_mem, "ARKStepCreate")) { return 1; }

  // Specify tolerances
  int flag = ARKodeSStolerances(*arkode_mem, uopts.rtol, uopts.atol);
  if (check_flag(flag, "ARKodeSStolerances")) { return 1; }

  // Attach user data
  flag = ARKodeSetUserData(*arkode_mem, &udata);
  if (check_flag(flag, "ARKodeSetUserData")) { return 1; }

  // If implicit, setup solvers
  if (fi_RHS)
  {
    // Create banded matrix
    *A = SUNBandMatrix(udata.neq, 3, 3, ctx);
    if (check_ptr(*A, "SUNBandMatrix")) { return 1; }

    // Create linear solver
    *LS = SUNLinSol_Band(y, *A, ctx);
    if (check_ptr(*LS, "SUNLinSol_Band")) { return 1; }

    // Attach linear solver
    flag = ARKodeSetLinearSolver(*arkode_mem, *LS, *A);
    if (check_flag(flag, "ARKodeSetLinearSolver")) { return 1; }

    // Attach Jacobian function
    flag = ARKodeSetJacFn(*arkode_mem, Ji_RHS);
    if (check_flag(flag, "ARKodeSetJacFn")) { return 1; }

    // Set the predictor method
    flag = ARKodeSetPredictorMethod(*arkode_mem, uopts.predictor);
    if (check_flag(flag, "ARKodeSetPredictorMethod")) { return 1; }

    // Set linear solver setup frequency
    flag = ARKodeSetLSetupFrequency(*arkode_mem, uopts.ls_setup_freq);
    if (check_flag(flag, "ARKodeSetLSetupFrequency")) { return 1; }

    if (uopts.linear)
    {
      // Specify linearly implicit non-time-dependent RHS
      flag = ARKodeSetLinear(*arkode_mem, SUNFALSE);
      if (check_flag(flag, "ARKodeSetLinear")) { return 1; }
    }
  }

  // Select method
  if (!fe_RHS && uopts.ark_dirk)
  {
    // Use the DIRK method from the default ARK method
    switch (uopts.order)
    {
    case (3):
      flag = ARKStepSetTableName(*arkode_mem, "ARKODE_ARK324L2SA_DIRK_4_2_3",
                                 "ARKODE_ERK_NONE");
      break;
    case (4):
      flag = ARKStepSetTableName(*arkode_mem, "ARKODE_ARK436L2SA_DIRK_6_3_4",
                                 "ARKODE_ERK_NONE");
      break;
    case (5):
      flag = ARKStepSetTableName(*arkode_mem, "ARKODE_ARK548L2SA_DIRK_8_4_5",
                                 "ARKODE_ERK_NONE");
      break;
    default:
      cerr << "ERROR: Invalid order to use ARK DIRK method" << endl;
      return -1;
      break;
    }
    if (check_flag(flag, "ARKStepSetTableNum")) { return 1; }
  }
  else
  {
    // Select default method of a given order
    flag = ARKodeSetOrder(*arkode_mem, uopts.order);
    if (check_flag(flag, "ARKodeSetOrder")) { return 1; }
  }

  // Set fixed step size or adaptivity method
  if (uopts.fixed_h > ZERO)
  {
    flag = ARKodeSetFixedStep(*arkode_mem, uopts.fixed_h);
    if (check_flag(flag, "ARKodeSetFixedStep")) { return 1; }
  }
  else
  {
    flag = ARKodeSetAdaptControllerByName(*arkode_mem, uopts.controller.c_str());
    if (check_flag(flag, "ARKodeSetAdaptControllerByName")) { return 1; }
  }

  // Set max steps between outputs
  flag = ARKodeSetMaxNumSteps(*arkode_mem, uopts.maxsteps);
  if (check_flag(flag, "ARKodeSetMaxNumSteps")) { return 1; }

  // Set stopping time
  flag = ARKodeSetStopTime(*arkode_mem, udata.tf);
  if (check_flag(flag, "ARKodeSetStopTime")) { return 1; }

  return 0;
}

int SetupMRIARK(SUNContext ctx, UserData& udata, UserOptions& uopts, N_Vector y,
                SUNMatrix* A, SUNLinearSolver* LS, SUNMatrix* A_fast,
                SUNLinearSolver* LS_fast, MRIStepInnerStepper* fast_mem,
                void** arkode_mem)
{
  // Problem configuration
  ARKRhsFn fse_RHS;   // slow explicit RHS function
  ARKRhsFn fsi_RHS;   // slow implicit RHS function
  ARKLsJacFn Jsi_RHS; // Jacobian of slow implicit RHS function
  ARKRhsFn ffe_RHS;   // fast explicit RHS function
  ARKRhsFn ffi_RHS;   // fast implicit RHS function
  ARKLsJacFn Jfi_RHS; // Jacobian of fast implicit RHS function

  // Slow time scale
  if (udata.diffusion && udata.advection)
  {
    // IMEX slow -- advection-diffusion
    fse_RHS = f_advection;
    fsi_RHS = f_diffusion;
    Jsi_RHS = J_diffusion;
  }
  else if (!udata.diffusion && udata.advection)
  {
    // Explicit slow -- advection
    fse_RHS = f_advection;
    fsi_RHS = nullptr;
    Jsi_RHS = nullptr;
  }
  else if (udata.diffusion && !udata.advection)
  {
    // Implicit slow -- diffusion
    fse_RHS = nullptr;
    fsi_RHS = f_diffusion;
    Jsi_RHS = J_diffusion;
  }
  else
  {
    // No slow time scale
    cerr << "ERROR: Invalid problem configuration" << endl;
    return -1;
  }

  // Fast time scale
  if (udata.splitting)
  {
    // Implicit fast -- reaction
    ffe_RHS = nullptr;
    ffi_RHS = f_reaction;
    Jfi_RHS = J_reaction;
  }
  else
  {
    // Explicit fast -- reaction
    ffe_RHS = f_reaction;
    ffi_RHS = nullptr;
    Jfi_RHS = nullptr;
  }

  // -------------------------
  // Setup the fast integrator
  // -------------------------

  // Create ARKStep memory
  void* fast_arkode_mem = ARKStepCreate(ffe_RHS, ffi_RHS, ZERO, y, ctx);
  if (check_ptr(arkode_mem, "ARKStepCreate")) { return 1; }

  // Specify tolerances
  int flag = ARKodeSStolerances(fast_arkode_mem, uopts.rtol_fast,
                                uopts.atol_fast);
  if (check_flag(flag, "ARKodeSStolerances")) { return 1; }

  // Attach user data
  flag = ARKodeSetUserData(fast_arkode_mem, &udata);
  if (check_flag(flag, "ARKodeSetUserData")) { return 1; }

  // If implicit, setup solvers
  if (ffi_RHS)
  {
    // Create banded matrix
    *A_fast = SUNBandMatrix(udata.neq, 2, 2, ctx);
    if (check_ptr(*A_fast, "SUNBandMatrix")) { return 1; }

    // Create linear solver
    *LS_fast = SUNLinSol_Band(y, *A_fast, ctx);
    if (check_ptr(*LS_fast, "SUNLinSol_Band")) { return 1; }

    // Attach linear solver
    flag = ARKodeSetLinearSolver(fast_arkode_mem, *LS_fast, *A_fast);
    if (check_flag(flag, "ARKodeSetLinearSolver")) { return 1; }

    // Attach Jacobian function
    flag = ARKodeSetJacFn(fast_arkode_mem, Jfi_RHS);
    if (check_flag(flag, "ARKodeSetJacFn")) { return 1; }

    // Set the predictor method
    flag = ARKodeSetPredictorMethod(fast_arkode_mem, uopts.predictor_fast);
    if (check_flag(flag, "ARKodeSetPredictorMethod")) { return 1; }

    // Set linear solver setup frequency
    flag = ARKodeSetLSetupFrequency(fast_arkode_mem, uopts.ls_setup_freq_fast);
    if (check_flag(flag, "ARKodeSetLSetupFrequency")) { return 1; }
  }

  // Select method order
  flag = ARKodeSetOrder(fast_arkode_mem, uopts.order_fast);
  if (check_flag(flag, "ARKodeSetOrder")) { return 1; }

  // Set fixed step size or adaptivity method
  if (uopts.fixed_h_fast > ZERO)
  {
    flag = ARKodeSetFixedStep(fast_arkode_mem, uopts.fixed_h_fast);
    if (check_flag(flag, "ARKodeSetFixedStep")) { return 1; }
  }
  else
  {
    flag = ARKodeSetAdaptControllerByName(fast_arkode_mem,
                                          uopts.controller_fast.c_str());
    if (check_flag(flag, "ARKodeSetAdaptControllerByName")) { return 1; }
  }

  // Set max steps between outputs
  flag = ARKodeSetMaxNumSteps(fast_arkode_mem, uopts.maxsteps);
  if (check_flag(flag, "ARKodeSetMaxNumSteps")) { return 1; }

  // Wrap ARKODE as an MRIStepInnerStepper
  flag = ARKodeCreateMRIStepInnerStepper(fast_arkode_mem, fast_mem);
  if (check_flag(flag, "ARKodeCreateMRIStepInnerStepper")) { return 1; }

  // -------------------------
  // Setup the slow integrator
  // -------------------------

  // Create slow integrator for diffusion and attach fast integrator
  *arkode_mem = MRIStepCreate(fse_RHS, fsi_RHS, ZERO, y, *fast_mem, ctx);
  if (check_ptr(*arkode_mem, "MRIStepCreate")) { return 1; }

  // Set the slow step size
  flag = ARKodeSetFixedStep(*arkode_mem, uopts.fixed_h);
  if (check_flag(flag, "ARKodeSetFixedStep")) { return 1; }

  // Specify tolerances
  flag = ARKodeSStolerances(*arkode_mem, uopts.rtol, uopts.atol);
  if (check_flag(flag, "ARKodeSStolerances")) { return 1; }

  // Attach user data
  flag = ARKodeSetUserData(*arkode_mem, &udata);
  if (check_flag(flag, "ARKodeSetUserData")) { return 1; }

  // If implicit, setup solvers
  if (fsi_RHS)
  {
    // Create banded matrix
    *A = SUNBandMatrix(udata.neq, 3, 3, ctx);
    if (check_ptr(*A, "SUNBandMatrix")) { return 1; }

    // Create linear solver
    *LS = SUNLinSol_Band(y, *A, ctx);
    if (check_ptr(*LS, "SUNLinSol_Band")) { return 1; }

    // Attach linear solver
    flag = ARKodeSetLinearSolver(*arkode_mem, *LS, *A);
    if (check_flag(flag, "ARKodeSetLinearSolver")) { return 1; }

    // Attach Jacobian function
    flag = ARKodeSetJacFn(*arkode_mem, Jsi_RHS);
    if (check_flag(flag, "ARKodeSetJacFn")) { return 1; }

    // Set linear solver setup frequency
    flag = ARKodeSetLSetupFrequency(*arkode_mem, uopts.ls_setup_freq);
    if (check_flag(flag, "ARKodeSetLSetupFrequency")) { return 1; }

    // Set the predictor method
    flag = ARKodeSetPredictorMethod(*arkode_mem, uopts.predictor);
    if (check_flag(flag, "ARKodeSetPredictorMethod")) { return 1; }

    if (uopts.linear)
    {
      // Specify linearly implicit non-time-dependent RHS
      flag = ARKodeSetLinear(*arkode_mem, SUNFALSE);
      if (check_flag(flag, "ARKodeSetLinear")) { return 1; }
    }
  }

  // Select method order
  flag = ARKodeSetOrder(*arkode_mem, uopts.order);
  if (check_flag(flag, "ARKodeSetOrder")) { return 1; }

  // Set max steps between outputs
  flag = ARKodeSetMaxNumSteps(*arkode_mem, uopts.maxsteps);
  if (check_flag(flag, "ARKodeSetMaxNumSteps")) { return 1; }

  // Set stopping time
  flag = ARKodeSetStopTime(*arkode_mem, udata.tf);
  if (check_flag(flag, "ARKodeSetStopTime")) { return 1; }

  return 0;
}

int SetupMRICVODE(SUNContext ctx, UserData& udata, UserOptions& uopts,
                  N_Vector y, SUNMatrix* A, SUNLinearSolver* LS,
                  SUNMatrix* A_fast, SUNLinearSolver* LS_fast,
                  MRIStepInnerStepper* fast_mem, void** arkode_mem)
{
  // Problem configuration
  ARKRhsFn fse_RHS;   // slow explicit RHS function
  ARKRhsFn fsi_RHS;   // slow implicit RHS function
  ARKLsJacFn Jsi_RHS; // Jacobian of slow implicit RHS function
  ARKRhsFn ff_RHS;    // fast RHS function
  ARKLsJacFn Jf_RHS;  // Jacobian of fast RHS function

  // Slow time scale
  if (udata.diffusion && udata.advection)
  {
    // IMEX slow -- advection-diffusion
    fse_RHS = f_advection;
    fsi_RHS = f_diffusion;
    Jsi_RHS = J_diffusion;
  }
  else if (!udata.diffusion && udata.advection)
  {
    // Explicit slow -- advection
    fse_RHS = f_advection;
    fsi_RHS = nullptr;
    Jsi_RHS = nullptr;
  }
  else if (udata.diffusion && !udata.advection)
  {
    // Implicit slow -- diffusion
    fse_RHS = nullptr;
    fsi_RHS = f_diffusion;
    Jsi_RHS = J_diffusion;
  }
  else
  {
    // No slow time scale
    cerr << "ERROR: Invalid problem configuration" << endl;
    return -1;
  }

  // Fast time scale -- Implicit fast reaction
  ff_RHS = f_react_forcing;
  Jf_RHS = J_reaction;

  // -------------------------
  // Setup the fast integrator
  // -------------------------

  // Create the solver memory and specify the Adams methods
  void* cvode_mem = CVodeCreate(CV_BDF, ctx);
  if (check_ptr(cvode_mem, "CVodeCreate")) { return 1; }

  // Initialize the integrator memory
  int flag = CVodeInit(cvode_mem, ff_RHS, ZERO, y);
  if (check_flag(flag, "CVodeInit")) { return 1; }

  // Specify tolerances
  flag = CVodeSStolerances(cvode_mem, uopts.rtol_fast, uopts.atol_fast);
  if (check_flag(flag, "CVodeSVtolerances")) { return 1; }

  // Attach user data
  flag = CVodeSetUserData(cvode_mem, &udata);
  if (check_flag(flag, "CVodeSetUserData")) { return 1; }

  // Create banded matrix
  *A_fast = SUNBandMatrix(udata.neq, 2, 2, ctx);
  if (check_ptr(*A_fast, "SUNBandMatrix")) { return 1; }

  // Create linear solver
  *LS_fast = SUNLinSol_Band(y, *A_fast, ctx);
  if (check_ptr(*LS_fast, "SUNLinSol_Band")) { return 1; }

  // Attach linear solver
  flag = CVodeSetLinearSolver(cvode_mem, *LS_fast, *A_fast);
  if (check_flag(flag, "CVodeSetLinearSolver")) { return 1; }

  // Attach Jacobian function
  flag = CVodeSetJacFn(cvode_mem, Jf_RHS);
  if (check_flag(flag, "CVodeSetJacFn")) { return 1; }

  // Set linear solver setup frequency
  flag = CVodeSetLSetupFrequency(cvode_mem, uopts.ls_setup_freq_fast);
  if (check_flag(flag, "CVodeSetLSetupFrequency")) { return 1; }

  // Set max step size change in first step
  flag = CVodeSetEtaMaxFirstStep(cvode_mem, uopts.etamx1_fast);
  if (check_flag(flag, "CVodeSetEtaMaxFirstStep")) { return 1; }

  // Set max steps between outputs
  flag = CVodeSetMaxNumSteps(cvode_mem, uopts.maxsteps);
  if (check_flag(flag, "CVodeSetMaxNumSteps")) { return 1; }

  // Create the inner stepper wrapper
  flag = MRIStepInnerStepper_Create(ctx, fast_mem);
  if (check_flag(flag, "MRIStepInnerStepper_Create")) { return 1; }

  // Attach memory and operations
  CVodeInnerStepperContent* inner_content = new CVodeInnerStepperContent;

  inner_content->cvode_mem   = cvode_mem;
  inner_content->user_data   = &udata;
  inner_content->save_hinit  = uopts.save_hinit;
  inner_content->save_hcur   = uopts.save_hcur;
  inner_content->hcur_factor = uopts.hcur_factor;

  flag = MRIStepInnerStepper_SetContent(*fast_mem, inner_content);
  if (check_flag(flag, "MRIStepInnerStepper_SetContent")) { return 1; }

  flag = MRIStepInnerStepper_SetEvolveFn(*fast_mem, CVodeInnerStepper_Evolve);
  if (check_flag(flag, "MRIStepInnerStepper_SetEvolve")) { return 1; }

  flag = MRIStepInnerStepper_SetFullRhsFn(*fast_mem, CVodeInnerStepper_FullRhs);
  if (check_flag(flag, "MRIStepInnerStepper_SetFullRhsFn")) { return 1; }

  flag = MRIStepInnerStepper_SetResetFn(*fast_mem, CVodeInnerStepper_Reset);
  if (check_flag(flag, "MRIStepInnerStepper_SetResetFn")) { return 1; }

  // Attach inner stepper memory to user data
  udata.fast_mem = *fast_mem;

  // -------------------------
  // Setup the slow integrator
  // -------------------------

  // Create slow integrator for diffusion and attach fast integrator
  *arkode_mem = MRIStepCreate(fse_RHS, fsi_RHS, ZERO, y, *fast_mem, ctx);
  if (check_ptr(*arkode_mem, "MRIStepCreate")) { return 1; }

  // Set the slow step size
  flag = ARKodeSetFixedStep(*arkode_mem, uopts.fixed_h);
  if (check_flag(flag, "ARKodeSetFixedStep")) { return 1; }

  // Specify tolerances
  flag = ARKodeSStolerances(*arkode_mem, uopts.rtol, uopts.atol);
  if (check_flag(flag, "ARKodeSStolerances")) { return 1; }

  // Attach user data
  flag = ARKodeSetUserData(*arkode_mem, &udata);
  if (check_flag(flag, "ARKodeSetUserData")) { return 1; }

  // If implicit, setup solvers
  if (fsi_RHS)
  {
    // Create banded matrix
    *A = SUNBandMatrix(udata.neq, 3, 3, ctx);
    if (check_ptr(*A, "SUNBandMatrix")) { return 1; }

    // Create linear solver
    *LS = SUNLinSol_Band(y, *A, ctx);
    if (check_ptr(*LS, "SUNLinSol_Band")) { return 1; }

    // Attach linear solver
    flag = ARKodeSetLinearSolver(*arkode_mem, *LS, *A);
    if (check_flag(flag, "ARKodeSetLinearSolver")) { return 1; }

    // Attach Jacobian function
    flag = ARKodeSetJacFn(*arkode_mem, Jsi_RHS);
    if (check_flag(flag, "ARKodeSetJacFn")) { return 1; }

    // Set linear solver setup frequency
    flag = ARKodeSetLSetupFrequency(*arkode_mem, uopts.ls_setup_freq);
    if (check_flag(flag, "ARKodeSetLSetupFrequency")) { return 1; }

    // Set the predictor method
    flag = ARKodeSetPredictorMethod(*arkode_mem, uopts.predictor);
    if (check_flag(flag, "ARKodeSetPredictorMethod")) { return 1; }

    if (uopts.linear)
    {
      // Specify linearly implicit non-time-dependent RHS
      flag = ARKodeSetLinear(*arkode_mem, SUNFALSE);
      if (check_flag(flag, "ARKodeSetLinear")) { return 1; }
    }
  }

  // Select method order
  flag = ARKodeSetOrder(*arkode_mem, uopts.order);
  if (check_flag(flag, "ARKodeSetOrder")) { return 1; }

  // Set max steps between outputs
  flag = ARKodeSetMaxNumSteps(*arkode_mem, uopts.maxsteps);
  if (check_flag(flag, "ARKodeSetMaxNumSteps")) { return 1; }

  // Set stopping time
  flag = ARKodeSetStopTime(*arkode_mem, udata.tf);
  if (check_flag(flag, "ARKodeSetStopTime")) { return 1; }

  return 0;
}

// -----------------------------------------------------------------------------
// Custom inner stepper functions
// -----------------------------------------------------------------------------

// Advance the fast ODE in time
int CVodeInnerStepper_Evolve(MRIStepInnerStepper fast_mem, sunrealtype t0,
                             sunrealtype tout, N_Vector y)
{
  void* inner_content = nullptr;
  int flag = MRIStepInnerStepper_GetContent(fast_mem, &inner_content);
  if (check_flag(flag, "MRIStepInnerStepper_GetContent")) { return -1; }

  CVodeInnerStepperContent* content = (CVodeInnerStepperContent*)inner_content;

  // Set initial step size (if saved)
  if (content->save_hinit && content->hinit > ZERO)
  {
    flag = CVodeSetInitStep(content->cvode_mem, content->hinit);
    if (flag) { return -1; }
  }
  else if (content->save_hcur && content->hcur > ZERO)
  {
    flag = CVodeSetInitStep(content->cvode_mem,
                            content->hcur_factor * content->hcur);
    if (flag) { return -1; }
  }

  // Evolve in time
  flag = CVodeSetStopTime(content->cvode_mem, tout);
  if (check_flag(flag, "CVodeSetStopTime")) { return -1; }

  sunrealtype tret;
  flag = CVode(content->cvode_mem, tout, y, &tret, CV_NORMAL);
  if (flag < 0) { return -1; }

  // Save the initial step size
  if (content->save_hinit)
  {
    flag = CVodeGetActualInitStep(content->cvode_mem, &(content->hinit));
    if (flag) { return -1; }
  }

  // Save the current step size
  if (content->save_hcur)
  {
    flag = CVodeGetCurrentStep(content->cvode_mem, &(content->hcur));
    if (flag) { return -1; }
  }

  return 0;
}

// Compute the RHS of the fast ODE
int CVodeInnerStepper_FullRhs(MRIStepInnerStepper fast_mem, sunrealtype t,
                              N_Vector y, N_Vector f, int mode)
{
  void* inner_content = nullptr;
  int flag = MRIStepInnerStepper_GetContent(fast_mem, &inner_content);
  if (check_flag(flag, "MRIStepInnerStepper_GetContent")) { return -1; }

  CVodeInnerStepperContent* content = (CVodeInnerStepperContent*)inner_content;

  flag = f_reaction(t, y, f, content->user_data);
  if (flag) { return -1; }

  return 0;
}

// Reset the fast integrator to the given time and state
int CVodeInnerStepper_Reset(MRIStepInnerStepper fast_mem, sunrealtype tR,
                            N_Vector yR)
{
  void* inner_content = nullptr;
  int flag = MRIStepInnerStepper_GetContent(fast_mem, &inner_content);
  if (check_flag(flag, "MRIStepInnerStepper_GetContent")) { return -1; }

  CVodeInnerStepperContent* content = (CVodeInnerStepperContent*)inner_content;

  // Save current stats before reinit
  flag = UpdateCVodeStats(content);
  if (check_flag(flag, "UpdateCVodeStats")) { return -1; }

  // Reinitialize CVODE with new state
  flag = CVodeReInit(content->cvode_mem, tR, yR);
  if (check_flag(flag, "CVodeReInit")) { return -1; }

  return 0;
}

// -----------------------------------------------------------------------------
// Functions called by the integrator
// -----------------------------------------------------------------------------

// Advection RHS function
int f_advection(sunrealtype t, N_Vector y, N_Vector f, void* user_data)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Access data arrays
  sunrealtype* ydata = N_VGetArrayPointer(y);
  if (check_ptr(ydata, "N_VGetArrayPointer")) { return -1; }

  sunrealtype* fdata = N_VGetArrayPointer(f);
  if (check_ptr(fdata, "N_VGetArrayPointer")) { return -1; }

  // Compute advection RHS
  sunrealtype ul, ur;
  sunrealtype vl, vr;
  sunrealtype wl, wr;

  sunrealtype c = -ONE * udata->c / (TWO * udata->dx);

  fdata[0] = fdata[1] = fdata[2] = ZERO;

  for (sunindextype i = 1; i < udata->nx - 1; i++)
  {
    ul = ydata[UIDX(i - 1)];
    ur = ydata[UIDX(i + 1)];

    vl = ydata[VIDX(i - 1)];
    vr = ydata[VIDX(i + 1)];

    wl = ydata[WIDX(i - 1)];
    wr = ydata[WIDX(i + 1)];

    fdata[UIDX(i)] = c * (ur - ul);
    fdata[VIDX(i)] = c * (vr - vl);
    fdata[WIDX(i)] = c * (wr - wl);
  }

  fdata[udata->neq - 3] = fdata[udata->neq - 2] = fdata[udata->neq - 1] = ZERO;

  return 0;
}

// Advection Jacobian function
int J_advection(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  sunrealtype c = -ONE * udata->c / (TWO * udata->dx);

  for (sunindextype i = 1; i < udata->nx - 1; i++)
  {
    SM_ELEMENT_B(J, UIDX(i), UIDX(i - 1)) = -c;
    SM_ELEMENT_B(J, UIDX(i), UIDX(i + 1)) = c;

    SM_ELEMENT_B(J, VIDX(i), VIDX(i - 1)) = -c;
    SM_ELEMENT_B(J, VIDX(i), VIDX(i + 1)) = c;

    SM_ELEMENT_B(J, WIDX(i), WIDX(i - 1)) = -c;
    SM_ELEMENT_B(J, WIDX(i), WIDX(i + 1)) = c;
  }

  return 0;
}

// Diffusion RHS function
int f_diffusion(sunrealtype t, N_Vector y, N_Vector f, void* user_data)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Access data arrays
  sunrealtype* ydata = N_VGetArrayPointer(y);
  if (check_ptr(ydata, "N_VGetArrayPointer")) { return -1; }

  sunrealtype* fdata = N_VGetArrayPointer(f);
  if (check_ptr(fdata, "N_VGetArrayPointer")) { return -1; }

  // Compute diffusion RHS
  sunrealtype ul, uc, ur;
  sunrealtype vl, vc, vr;
  sunrealtype wl, wc, wr;

  sunrealtype d = udata->d / (udata->dx * udata->dx);

  fdata[0] = fdata[1] = fdata[2] = ZERO;

  for (sunindextype i = 1; i < udata->nx - 1; i++)
  {
    ul = ydata[UIDX(i - 1)];
    uc = ydata[UIDX(i)];
    ur = ydata[UIDX(i + 1)];

    vl = ydata[VIDX(i - 1)];
    vc = ydata[VIDX(i)];
    vr = ydata[VIDX(i + 1)];

    wl = ydata[WIDX(i - 1)];
    wc = ydata[WIDX(i)];
    wr = ydata[WIDX(i + 1)];

    fdata[UIDX(i)] = d * (ul - TWO * uc + ur);
    fdata[VIDX(i)] = d * (vl - TWO * vc + vr);
    fdata[WIDX(i)] = d * (wl - TWO * wc + wr);
  }

  fdata[udata->neq - 3] = fdata[udata->neq - 2] = fdata[udata->neq - 1] = ZERO;

  return 0;
}

// Diffusion Jacobian function
int J_diffusion(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  sunrealtype d = udata->d / (udata->dx * udata->dx);

  for (sunindextype i = 1; i < udata->nx - 1; i++)
  {
    SM_ELEMENT_B(J, UIDX(i), UIDX(i - 1)) = d;
    SM_ELEMENT_B(J, UIDX(i), UIDX(i))     = -d * TWO;
    SM_ELEMENT_B(J, UIDX(i), UIDX(i + 1)) = d;

    SM_ELEMENT_B(J, VIDX(i), VIDX(i - 1)) = d;
    SM_ELEMENT_B(J, VIDX(i), VIDX(i))     = -d * TWO;
    SM_ELEMENT_B(J, VIDX(i), VIDX(i + 1)) = d;

    SM_ELEMENT_B(J, WIDX(i), WIDX(i - 1)) = d;
    SM_ELEMENT_B(J, WIDX(i), WIDX(i))     = -d * TWO;
    SM_ELEMENT_B(J, WIDX(i), WIDX(i + 1)) = d;
  }

  return 0;
}

// Reaction RHS function
int f_reaction(sunrealtype t, N_Vector y, N_Vector f, void* user_data)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Access data arrays
  sunrealtype* ydata = N_VGetArrayPointer(y);
  if (check_ptr(ydata, "N_VGetArrayPointer")) { return -1; }

  sunrealtype* fdata = N_VGetArrayPointer(f);
  if (check_ptr(fdata, "N_VGetArrayPointer")) { return -1; }

  // Compute reaction RHS
  sunrealtype u, v, w;

  fdata[0] = fdata[1] = fdata[2] = ZERO;

  for (sunindextype i = 1; i < udata->nx - 1; i++)
  {
    u = ydata[UIDX(i)];
    v = ydata[VIDX(i)];
    w = ydata[WIDX(i)];

    fdata[UIDX(i)] = udata->A - (w + ONE) * u + v * u * u;
    fdata[VIDX(i)] = w * u - v * u * u;
    fdata[WIDX(i)] = ((udata->B - w) / udata->eps) - w * u;
  }

  fdata[udata->neq - 3] = fdata[udata->neq - 2] = fdata[udata->neq - 1] = ZERO;

  return 0;
}

// Diffusion Jacobian function
int J_reaction(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Access data array
  sunrealtype* ydata = N_VGetArrayPointer(y);
  if (check_ptr(ydata, "N_VGetArrayPointer")) { return 1; }

  sunrealtype u, v, w;

  for (sunindextype i = 1; i < udata->nx - 1; i++)
  {
    u = ydata[UIDX(i)];
    v = ydata[VIDX(i)];
    w = ydata[WIDX(i)];

    // all vars wrt u
    SM_ELEMENT_B(J, UIDX(i), UIDX(i)) = -(w + ONE) + TWO * u * v;
    SM_ELEMENT_B(J, VIDX(i), UIDX(i)) = w - TWO * u * v;
    SM_ELEMENT_B(J, WIDX(i), UIDX(i)) = -w;

    // all vars wrt v
    SM_ELEMENT_B(J, UIDX(i), VIDX(i)) = u * u;
    SM_ELEMENT_B(J, VIDX(i), VIDX(i)) = -u * u;

    // all vars wrt w
    SM_ELEMENT_B(J, UIDX(i), WIDX(i)) = -u;
    SM_ELEMENT_B(J, VIDX(i), WIDX(i)) = u;
    SM_ELEMENT_B(J, WIDX(i), WIDX(i)) = (-ONE / udata->eps) - u;
  }

  return 0;
}

// Advection-diffusion RHS function
int f_adv_diff(sunrealtype t, N_Vector y, N_Vector f, void* user_data)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Compute advection
  int flag = f_advection(t, y, f, user_data);
  if (flag) { return flag; }

  // Compute diffusion
  flag = f_diffusion(t, y, udata->temp_v, user_data);
  if (flag) { return flag; }

  // Combine advection and reaction
  N_VLinearSum(ONE, f, ONE, udata->temp_v, f);

  return 0;
}

// Advection-diffusion Jacobian function
int J_adv_diff(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Compute diffusion Jacobian
  int flag = J_advection(t, y, fy, J, user_data, tmp1, tmp2, tmp3);
  if (flag) { return flag; }

  // Compute diffusion Jacobian
  flag = SUNMatZero(udata->temp_J);
  if (flag) { return flag; }

  flag = J_diffusion(t, y, fy, udata->temp_J, user_data, tmp1, tmp2, tmp3);
  if (flag) { return flag; }

  // Combine Jacobians
  flag = SUNMatScaleAdd(ONE, J, udata->temp_J);
  if (flag) { return -1; }

  return 0;
}

// Advection-reaction RHS function
int f_adv_react(sunrealtype t, N_Vector y, N_Vector f, void* user_data)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Compute advection
  int flag = f_advection(t, y, f, user_data);
  if (flag) { return flag; }

  // Compute reactions
  flag = f_reaction(t, y, udata->temp_v, user_data);
  if (flag) { return flag; }

  // Combine advection and reaction
  N_VLinearSum(ONE, f, ONE, udata->temp_v, f);

  return 0;
}

// Diffusion-reaction Jacobian function
int J_adv_react(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Compute advection Jacobian
  int flag = J_advection(t, y, fy, J, user_data, tmp1, tmp2, tmp3);
  if (flag) { return flag; }

  // Compute reaction Jacobian
  flag = SUNMatZero(udata->temp_J);
  if (flag) { return flag; }

  flag = J_reaction(t, y, fy, udata->temp_J, user_data, tmp1, tmp2, tmp3);
  if (flag) { return flag; }

  // Combine Jacobians
  flag = SUNMatScaleAdd(ONE, J, udata->temp_J);
  if (flag) { return -1; }

  return 0;
}

// Diffusion-reaction RHS function
int f_diff_react(sunrealtype t, N_Vector y, N_Vector f, void* user_data)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Compute diffusion
  int flag = f_diffusion(t, y, f, user_data);
  if (flag) { return flag; }

  // Compute reactions
  flag = f_reaction(t, y, udata->temp_v, user_data);
  if (flag) { return flag; }

  // Combine advection and reaction
  N_VLinearSum(ONE, f, ONE, udata->temp_v, f);

  return 0;
}

// Diffusion-reaction Jacobian function
int J_diff_react(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                 void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Compute diffusion Jacobian
  int flag = J_diffusion(t, y, fy, J, user_data, tmp1, tmp2, tmp3);
  if (flag) { return flag; }

  // Compute reaction Jacobian
  flag = SUNMatZero(udata->temp_J);
  if (flag) { return flag; }

  flag = J_reaction(t, y, fy, udata->temp_J, user_data, tmp1, tmp2, tmp3);
  if (flag) { return flag; }

  // Combine Jacobians
  flag = SUNMatScaleAdd(ONE, J, udata->temp_J);
  if (flag) { return -1; }

  return 0;
}

// Advection-diffusion-reaction RHS function
int f_adv_diff_react(sunrealtype t, N_Vector y, N_Vector f, void* user_data)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Compute advection
  int flag = f_advection(t, y, f, user_data);
  if (flag) { return flag; }

  // Compute diffusion
  flag = f_diffusion(t, y, udata->temp_v, user_data);
  if (flag) { return flag; }

  // Combine advection and reaction
  N_VLinearSum(ONE, f, ONE, udata->temp_v, f);

  // Compute reactions
  flag = f_reaction(t, y, udata->temp_v, user_data);
  if (flag) { return flag; }

  // Combine advection and reaction
  N_VLinearSum(ONE, f, ONE, udata->temp_v, f);

  return 0;
}

// Diffusion-reaction Jacobian function
int J_adv_diff_react(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                     void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Compute diffusion Jacobian
  int flag = J_advection(t, y, fy, J, user_data, tmp1, tmp2, tmp3);
  if (flag) { return flag; }

  // Compute reaction Jacobian
  flag = SUNMatZero(udata->temp_J);
  if (flag) { return flag; }

  flag = J_diffusion(t, y, fy, udata->temp_J, user_data, tmp1, tmp2, tmp3);
  if (flag) { return flag; }

  // Combine Jacobians
  flag = SUNMatScaleAdd(ONE, J, udata->temp_J);
  if (flag) { return -1; }

  // Compute reaction Jacobian
  flag = SUNMatZero(udata->temp_J);
  if (flag) { return flag; }

  flag = J_reaction(t, y, fy, udata->temp_J, user_data, tmp1, tmp2, tmp3);
  if (flag) { return flag; }

  // Combine Jacobians
  flag = SUNMatScaleAdd(ONE, J, udata->temp_J);
  if (flag) { return -1; }

  return 0;
}

// Reaction RHS function with MRI forcing
int f_react_forcing(sunrealtype t, N_Vector y, N_Vector f, void* user_data)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Compute reaction RHS
  int flag = f_reaction(t, y, f, user_data);
  if (flag) { return flag; }

  // Apply inner forcing for MRI + CVODE
  flag = MRIStepInnerStepper_AddForcing(udata->fast_mem, t, f);
  if (check_flag(flag, "MRIStepInnerStepper_AddForcing")) { return -1; }

  return 0;
}

// Compute the initial condition
int SetIC(N_Vector y, UserData& udata)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);
  if (check_ptr(ydata, "N_VGetArrayPointer")) { return -1; }

  sunrealtype x, p;

  for (sunindextype i = 0; i < udata.nx; i++)
  {
    x              = udata.xl + i * udata.dx;
    p              = SUN_RCONST(0.1) * sin(PI * x);
    ydata[UIDX(i)] = udata.A + p;
    ydata[VIDX(i)] = udata.B / udata.A + p;
    ydata[WIDX(i)] = udata.B + p;
  }

  return 0;
}

//---- end of file ----
