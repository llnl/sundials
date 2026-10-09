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
 * This example problem is adapted from:
 *
 * H. Ranocha, M. Sayyari, L. Dalcin, M. Parsani, and D.I. Ketcheson,
 * "Relaxation Runge-Kutta Methods: Fully-Discrete Explicit Entropy-Stable
 * Schemes for the Compressible Euler and Navier-Stokes Equations," SIAM Journal
 * on Scientific Computing, 42(2), 2020, https://doi.org/10.1137/19M1263480.
 * -----------------------------------------------------------------------------
 * This example evolves system y = [u, v]^T
 *
 *   du/dt = -sin(v)
 *   dv/dt = u
 *
 * for t in the interval [0, 10] with the initial condition y0 = [1.5 1.0]^T.
 * The conserved energy and its Jacobian for the system are
 *
 *   e(y) = 0.5 u^2 - cos(v) and e'(y) = [u, sin(v)]^T
 *
 * The problem is advanced in time with an explicit or implicit relaxed
 * Runge-Kutta method to ensure conservation of the energy.
 * ---------------------------------------------------------------------------*/

// Standard headers
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

// Common utility functions
#include <example_utilities.hpp>

// SUNDIALS headers
#include <arkode/arkode_arkstep.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>

#include "arkode/arkode_butcher.h"
#include "sundials/sundials_types.h"

/* ----------------------- *
 * User-supplied functions *
 * ----------------------- */

// ODE RHS function
int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

// ODE RHS Jacobian function
int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Energy function
int Eng(N_Vector y, sunrealtype* e, void* user_data);

// Energy Jacobian function
int JacEng(N_Vector y, N_Vector J, void* user_data);

/* ------------ *
 * Main Program *
 * ------------ */

int main(int argc, char* argv[])
{
  // Create the SUNDIALS context object for this simulation
  sundials::Context ctx;

  // Initial and final times
  sunrealtype t0 = SUN_RCONST(0.0);
  sunrealtype tf = SUN_RCONST(10.0);

  // Relative and absolute tolerances
  sunrealtype reltol = SUN_RCONST(1.0e-6);
  sunrealtype abstol = SUN_RCONST(1.0e-10);

  // Command line options
  int relax           = 1;               // enable relaxation
  int implicit        = 1;               // implicit
  sunrealtype fixed_h = SUN_RCONST(0.0); // adaptive

  /* -------------------- *
   * Output Problem Setup *
   * -------------------- */

  if (argc > 1) { relax = atoi(argv[1]); }
  if (argc > 2) { implicit = atoi(argv[2]); }
  if (argc > 3) { fixed_h = atof(argv[3]); }

  std::cout << "Nonlinear Pendulum problem:\n";
  if (implicit) { std::cout << "   method     = DIRK\n"; }
  else { std::cout << "   method     = ERK\n"; }
  if (relax) { std::cout << "   relaxation = ON\n"; }
  else { std::cout << "   relaxation = OFF\n"; }
  std::cout << "   reltol     = " << reltol << "\n";
  std::cout << "   abstol     = " << abstol << "\n";
  if (fixed_h > SUN_RCONST(0.0))
  {
    std::cout << "   fixed h    = " << fixed_h << "\n";
  }
  std::cout << std::endl;

  /* ------------ *
   * Setup ARKODE *
   * ------------ */

  // Create serial vector and set the initial condition values
  N_Vector y = N_VNew_Serial(2, ctx);
  if (check_ptr(y, "N_VNew_Serial")) { return 1; }

  sunrealtype* ydata = N_VGetArrayPointer(y);
  if (check_ptr(ydata, "N_VGetArrayPointer")) { return 1; }

  ydata[0] = SUN_RCONST(1.5);
  ydata[1] = SUN_RCONST(1.0);

  // Initial energy
  sunrealtype eng0;
  int flag = Eng(y, &eng0, nullptr);
  if (check_flag(flag, "Eng")) { return 1; }

  // Initialize ARKStep
  void* arkode_mem = nullptr;
  if (implicit) { arkode_mem = ARKStepCreate(nullptr, f, t0, y, ctx); }
  else { arkode_mem = ARKStepCreate(f, nullptr, t0, y, ctx); }
  if (check_ptr(arkode_mem, "ARKStepCreate")) { return 1; }

  // Specify tolerances
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(flag, "ARKodeSStolerances")) { return 1; }

  if (relax)
  {
    // Enable relaxation methods
    flag = ARKodeSetRelaxFn(arkode_mem, Eng, JacEng);
    if (check_flag(flag, "ARKodeSetRelaxFn")) { return 1; }
  }

  SUNMatrix A        = nullptr;
  SUNLinearSolver LS = nullptr;

  if (implicit)
  {
    // Create dense matrix and linear solver
    A = SUNDenseMatrix(2, 2, ctx);
    if (check_ptr(A, "SUNDenseMatrix")) { return 1; }

    LS = SUNLinSol_Dense(y, A, ctx);
    if (check_ptr(LS, "SUNLinSol_Dense")) { return 1; }

    // Attach the matrix and linear solver
    flag = ARKodeSetLinearSolver(arkode_mem, LS, A);
    if (check_flag(flag, "ARKodeSetLinearSolver")) { return 1; }

    // Set Jacobian routine
    flag = ARKodeSetJacFn(arkode_mem, Jac);
    if (check_flag(flag, "ARKodeSetJacFn")) { return 1; }

    if (fixed_h > SUN_RCONST(0.0))
    {
      // 3rd-order SDIRK method of Norsett
      ARKodeButcherTable B = ARKodeButcherTable_Alloc(2, SUNFALSE);

      const sunrealtype gamma =
        ((SUN_RCONST(3.0) + std::sqrt(SUN_RCONST(3.0))) / SUN_RCONST(6.0));

      B->A[0][0] = gamma;
      B->A[1][0] = SUN_RCONST(1.0) - SUN_RCONST(2.0) * gamma;
      B->A[1][1] = gamma;

      B->c[0] = gamma;
      B->c[1] = SUN_RCONST(1.0) - gamma;

      B->b[0] = SUN_RCONST(0.5);
      B->b[1] = SUN_RCONST(0.5);

      B->q = 3;
      B->p = 0;

      flag = ARKStepSetTables(arkode_mem, 3, 0, B, nullptr);
      if (check_flag(flag, "ARKStepSetTables")) { return 1; }

      ARKodeButcherTable_Free(B);
    }
    else
    {
      // Select a Butcher table with non-negative b values
      flag = ARKStepSetTableName(arkode_mem, "ARKODE_SDIRK_2_1_2",
                                 "ARKODE_ERK_NONE");
      if (check_flag(flag, "ARKStepSetTableName")) { return 1; }
    }
  }

  if (fixed_h > SUN_RCONST(0.0))
  {
    flag = ARKodeSetFixedStep(arkode_mem, fixed_h);
    if (check_flag(flag, "ARKodeSetFixedStep")) { return 1; }
  }

  flag = ARKodeSetNonlinConvCoef(arkode_mem, SUN_RCONST(0.01));
  if (check_flag(flag, "ARKodeSetNonlinConvCoef")) { return 1; }

  flag = ARKodeSetAutonomous(arkode_mem, SUNTRUE);
  if (check_flag(flag, "ARKodeSetAutonomous")) { return 1; }

  /* --------------- *
   * Advance in Time *
   * --------------- */

  // Initial time
  sunrealtype t = t0;

  // Output the initial condition and energy
  int swidth = 8;
  int rwidth = std::numeric_limits<sunrealtype>::digits10 + 8;

  std::ofstream outfile("ark_pendulum.txt");
  outfile << "# vars: t u v energy energy_err\n";
  outfile << std::scientific;
  outfile << std::setprecision(std::numeric_limits<sunrealtype>::digits10);
  outfile << t << " " << ydata[0] << " " << ydata[1] << " " << eng0 << " "
          << SUN_RCONST(0.0) << std::endl;

  std::cout << std::setw(swidth) << "step" << std::setw(rwidth) << "t"
            << std::setw(rwidth) << "u" << std::setw(rwidth) << "v"
            << std::setw(rwidth) << "e" << std::setw(rwidth) << "e err"
            << std::endl;
  for (int i = 0; i < swidth + 5 * rwidth; i++) { std::cout << "-"; }

  std::cout << std::endl;
  std::cout << std::scientific;
  std::cout << std::setprecision(std::numeric_limits<sunrealtype>::digits10);
  std::cout << std::setw(swidth) << 0 << std::setw(rwidth) << t
            << std::setw(rwidth) << ydata[0] << std::setw(rwidth) << ydata[1]
            << std::setw(rwidth) << eng0 << std::setw(rwidth) << SUN_RCONST(0.0);
  std::cout << std::endl;

  while (t < tf)
  {
    // Evolve in time
    flag = ARKodeEvolve(arkode_mem, tf, y, &t, ARK_ONE_STEP);
    if (check_flag(flag, "ARKodeEvolve")) { break; }

    // Output solution and errors
    sunrealtype eng;
    flag = Eng(y, &eng, nullptr);
    if (check_flag(flag, "Eng")) { return 1; }

    sunrealtype eng_chng = eng - eng0;

    /* Output to the screen periodically */
    long int nst;
    flag = ARKodeGetNumSteps(arkode_mem, &nst);
    check_flag(flag, "ARKodeGetNumSteps");

    if (nst % 1000 == 0)
    {
      std::cout << std::setw(swidth) << nst << std::setw(rwidth) << t
                << std::setw(rwidth) << ydata[0] << std::setw(rwidth)
                << ydata[1] << std::setw(rwidth) << eng << std::setw(rwidth)
                << eng_chng << std::endl;
    }

    /* Write all steps to file */
    outfile << t << " " << ydata[0] << " " << ydata[1] << " " << eng << " "
            << eng_chng << std::endl;
  }

  for (int i = 0; i < swidth + 5 * rwidth; i++) { std::cout << "-"; }
  std::cout << std::endl;
  outfile.close();

  /* ------------ *
   * Output Stats *
   * ------------ */

  // ARKODE statistics
  long int nst, nst_a, netf, nfe, nfi;

  // Get final statistics on how the solve progressed
  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  check_flag(flag, "ARKodeGetNumSteps");

  flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  check_flag(flag, "ARKodeGetNumStepAttempts");

  flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
  check_flag(flag, "ARKodeGetNumErrTestFails");

  flag = ARKodeGetNumRhsEvals(arkode_mem, 0, &nfe);
  check_flag(flag, "ARKodeGetNumRhsEvals");

  flag = ARKodeGetNumRhsEvals(arkode_mem, 1, &nfi);
  check_flag(flag, "ARKodeGetNumRhsEvals");

  std::cout << std::endl;
  std::cout << "Final Solver Statistics:\n";
  std::cout << "  Internal solver steps = " << nst << " (attempted = " << nst_a
            << ")\n";
  std::cout << "  Total number of error test failures = " << netf << "\n";
  std::cout << "  Total RHS evals:  Fe = " << nfe << ",  Fi = " << nfi << "\n";

  if (implicit)
  {
    long int nsetups, nje, nfeLS, nni, ncfn;

    flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
    check_flag(flag, "ARKodeGetNumNonlinSolvIters");

    flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
    check_flag(flag, "ARKodeGetNumNonlinSolvConvFails");

    flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
    check_flag(flag, "ARKodeGetNumLinSolvSetups");

    flag = ARKodeGetNumJacEvals(arkode_mem, &nje);
    check_flag(flag, "ARKodeGetNumJacEvals");

    flag = ARKodeGetNumLinRhsEvals(arkode_mem, &nfeLS);
    check_flag(flag, "ARKodeGetNumLinRhsEvals");

    std::cout << "  Total number of Newton iterations = " << nni << "\n";
    std::cout << "  Total number of linear solver convergence failures = " << ncfn
              << "\n";
    std::cout << "  Total linear solver setups = " << nsetups << "\n";
    std::cout << "  Total number of Jacobian evaluations = " << nje << "\n";
    std::cout << "  Total RHS evals for setting up the linear system = " << nfeLS
              << "\n";
  }

  if (relax)
  {
    long int nre, nrje, nrf, nrbf, nrnlsi, nrnlsf;

    flag = ARKodeGetNumRelaxFnEvals(arkode_mem, &nre);
    check_flag(flag, "ARKodeGetNumRelaxFnEvals");

    flag = ARKodeGetNumRelaxJacEvals(arkode_mem, &nrje);
    check_flag(flag, "ARKodeGetNumRelaxJacEvals");

    flag = ARKodeGetNumRelaxFails(arkode_mem, &nrf);
    check_flag(flag, "ARKodeGetNumRelaxFails");

    flag = ARKodeGetNumRelaxBoundFails(arkode_mem, &nrbf);
    check_flag(flag, "ARKodeGetNumRelaxBoundFails");

    flag = ARKodeGetNumRelaxSolveFails(arkode_mem, &nrnlsf);
    check_flag(flag, "ARKodeGetNumRelaxSolveFails");

    flag = ARKodeGetNumRelaxSolveIters(arkode_mem, &nrnlsi);
    check_flag(flag, "ARKodeGetNumRelaxSolveIters");

    std::cout << "  Total Relaxation Fn evals    = " << nre << "\n";
    std::cout << "  Total Relaxation Jac evals   = " << nrje << "\n";
    std::cout << "  Total Relaxation fails       = " << nrf << "\n";
    std::cout << "  Total Relaxation bound fails = " << nrbf << "\n";
    std::cout << "  Total Relaxation NLS fails   = " << nrnlsf << "\n";
    std::cout << "  Total Relaxation NLS iters   = " << nrnlsi << "\n";
  }
  std::cout << "\n";

  /* -------- *
   * Clean up *
   * -------- */

  // Free ARKODE integrator and SUNDIALS objects
  ARKodeFree(&arkode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);

  return flag;
}

/* ----------------------- *
 * User-supplied functions *
 * ----------------------- */

// ODE RHS function f(t,y).
int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* fdata = N_VGetArrayPointer(ydot);

  fdata[0] = -std::sin(ydata[1]);
  fdata[1] = ydata[0];

  return 0;
}

// ODE RHS Jacobian function J(t,y) = df/dy.
int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* Jdata = SUNDenseMatrix_Data(J);

  // column 0
  Jdata[0] = SUN_RCONST(0.0);
  Jdata[1] = SUN_RCONST(1.0);

  // column 1
  Jdata[2] = -std::cos(ydata[1]);
  Jdata[3] = SUN_RCONST(0.0);

  return 0;
}

// Energy function e(y)
int Eng(N_Vector y, sunrealtype* e, void* user_data)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);

  *e = SUN_RCONST(0.5) * ydata[0] * ydata[0] - std::cos(ydata[1]);

  return 0;
}

// Energy function Jacobian Je(y) = de/dy
int JacEng(N_Vector y, N_Vector J, void* user_data)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* jdata = N_VGetArrayPointer(J);

  jdata[0] = ydata[0];
  jdata[1] = std::sin(ydata[1]);

  return 0;
}
