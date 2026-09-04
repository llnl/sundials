/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ UMBC
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
 * Example problem:
 *
 * The following test simulates a simple anisotropic 2D heat equation,
 *
 *   u_t = kx u_xx + ky u_yy + b,
 *
 * for t in [0, 1] and (x,y) in [0, 1]^2, with initial conditions
 *
 *   u(0,x,y) = sin^2(pi x) sin^2(pi y) + 1,
 *
 * stationary boundary conditions
 *
 *   u_t(t,0,y) = u_t(t,1,y) = u_t(t,x,0) = u_t(t,x,1) = 0,
 *
 * and the heat source
 *
 *   b(t,x,y) = -2 pi sin^2(pi x) sin^2(pi y) sin(pi t) cos(pi t)
 *              - kx 2 pi^2 (cos^2(pi x) - sin^2(pi x)) sin^2(pi y) cos^2(pi t)
 *              - ky 2 pi^2 (cos^2(pi y) - sin^2(pi y)) sin^2(pi x) cos^2(pi t).
 *
 * Under this setup, the problem has the analytical solution
 *
 *    u(t,x,y) = sin^2(pi x) sin^2(pi y) cos^2(pi t) + 1.
 *
 * The spatial derivatives are computed using second-order centered differences,
 * with the data distributed over nx * ny points on a uniform spatial grid. The
 * problem is advanced in time with BDF methods using an inexact Newton method
 * paired with the PCG or SPGMR linear solver. Several command line options are
 * available to change the problem parameters and CVODE settings. Use the flag
 * --help for more information.
 * ---------------------------------------------------------------------------*/

// Include user data structure and utility functions for this problem
#include "cv_heat2D.hpp"

// Include integrator, vector, matrix, and linear solver headers
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_pcg.h>
#include <sunlinsol/sunlinsol_spgmr.h>

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrator
// -----------------------------------------------------------------------------

// ODE right hand side function
int f(sunrealtype t, N_Vector u, N_Vector f, void* user_data);

// Preconditioner Setup and Solve functions
int PSetup(sunrealtype t, N_Vector u, N_Vector f, sunbooleantype jok,
           sunbooleantype* jcurPtr, sunrealtype gamma, void* user_data);

int PSolve(sunrealtype t, N_Vector u, N_Vector f, N_Vector r, N_Vector z,
           sunrealtype gamma, sunrealtype delta, int lr, void* user_data);

// -----------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  // SUNDIALS context
  sundials::Context sunctx;

  // ---------------
  // Setup UserData
  // ---------------

  // Read input options
  UserData udata;
  std::vector<std::string> args(argv + 1, argv + argc);
  if (ReadInputs(args, udata)) { return 1; }
  PrintUserData(udata);

  // ---------------
  // Create vectors
  // ---------------

  // Create solution vector
  N_Vector u = N_VNew_Serial(udata.nodes, sunctx);
  if (check_ptr(u, "N_VNew_Serial")) { return 1; }

  // Set initial condition
  int flag = Solution(ZERO, u, udata);
  if (check_flag(flag, "Solution")) { return 1; }

  // Create error vector
  N_Vector e = N_VClone(u);
  if (check_ptr(e, "N_VClone")) { return 1; }

  // ---------------------
  // Create linear solver
  // ---------------------

  // Create linear solver
  int prectype = (udata.prec) ? SUN_PREC_RIGHT : SUN_PREC_NONE;

  SUNLinearSolver LS = nullptr;
  if (udata.pcg)
  {
    LS = SUNLinSol_PCG(u, prectype, udata.liniters, sunctx);
    if (check_ptr(LS, "SUNLinSol_PCG")) { return 1; }
  }
  else
  {
    LS = SUNLinSol_SPGMR(u, prectype, udata.liniters, sunctx);
    if (check_ptr(LS, "SUNLinSol_SPGMR")) { return 1; }
  }

  // Allocate preconditioner workspace
  if (udata.prec)
  {
    udata.d = N_VClone(u);
    if (check_ptr((udata.d), "N_VClone")) { return 1; }
  }

  // --------------
  // Setup CVODE
  // --------------

  // Create integrator
  void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_ptr(cvode_mem, "CVodeCreate")) { return 1; }

  // Initialize integrator
  flag = CVodeInit(cvode_mem, f, ZERO, u);
  if (check_flag(flag, "CVodeInit")) { return 1; }

  // Specify tolerances
  flag = CVodeSStolerances(cvode_mem, udata.rtol, udata.atol);
  if (check_flag(flag, "CVodeSStolerances")) { return 1; }

  // Attach user data
  flag = CVodeSetUserData(cvode_mem, (void*)&udata);
  if (check_flag(flag, "CVodeSetUserData")) { return 1; }

  // Attach linear solver
  flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
  if (check_flag(flag, "CVodeSetLinearSolver")) { return 1; }

  if (udata.prec)
  {
    // Attach preconditioner
    flag = CVodeSetPreconditioner(cvode_mem, PSetup, PSolve);
    if (check_flag(flag, "CVodeSetPreconditioner")) { return 1; }

    // Set linear solver setup frequency (update preconditioner)
    flag = CVodeSetLSetupFrequency(cvode_mem, udata.msbp);
    if (check_flag(flag, "CVodeSetLSetupFrequency")) { return 1; }
  }

  // Set linear solver tolerance factor
  flag = CVodeSetEpsLin(cvode_mem, udata.epslin);
  if (check_flag(flag, "CVodeSetEpsLin")) { return 1; }

  // Set max steps between outputs
  flag = CVodeSetMaxNumSteps(cvode_mem, udata.maxsteps);
  if (check_flag(flag, "CVodeSetMaxNumSteps")) { return 1; }

  // Set stopping time
  flag = CVodeSetStopTime(cvode_mem, udata.tf);
  if (check_flag(flag, "CVodeSetStopTime")) { return 1; }

  // -----------------------
  // Loop over output times
  // -----------------------

  auto t     = static_cast<sunrealtype>(ZERO);
  auto dTout = static_cast<sunrealtype>(udata.tf / udata.nout);
  auto tout  = dTout;

  // Initial output
  flag = OpenOutput(udata);
  if (check_flag(flag, "OpenOutput")) { return 1; }

  flag = WriteOutput(t, u, e, udata);
  if (check_flag(flag, "WriteOutput")) { return 1; }

  for (int iout = 0; iout < udata.nout; iout++)
  {
    // Evolve in time
    flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if (check_flag(flag, "CVode")) { break; }

    // Output solution and error
    flag = WriteOutput(t, u, e, udata);
    if (check_flag(flag, "WriteOutput")) { return 1; }

    // Update output time
    tout += dTout;
    tout = (tout > udata.tf) ? udata.tf : tout;
  }

  // Close output
  flag = CloseOutput(udata);
  if (check_flag(flag, "CloseOutput")) { return 1; }

  // --------------
  // Final outputs
  // --------------

  // Print integrator and solver stats
  long int nst, netf, nf, nni, ncfn, nli, nlcf, nsetups, nf_ls, nJv;
  flag = CVodeGetNumSteps(cvode_mem, &nst);
  if (check_flag(flag, "CVodeGetNumSteps")) { return -1; }
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  if (check_flag(flag, "CVodeGetNumErrTestFails")) { return -1; }
  flag = CVodeGetNumRhsEvals(cvode_mem, &nf);
  if (check_flag(flag, "CVodeGetNumRhsEvals")) { return -1; }
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  if (check_flag(flag, "CVodeGetNumNonlinSolvIters")) { return -1; }
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  if (check_flag(flag, "CVodeGetNumNonlinSolvConvFails")) { return -1; }
  flag = CVodeGetNumLinIters(cvode_mem, &nli);
  if (check_flag(flag, "CVodeGetNumLinIters")) { return -1; }
  flag = CVodeGetNumLinConvFails(cvode_mem, &nlcf);
  if (check_flag(flag, "CVodeGetNumLinConvFails")) { return -1; }
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  if (check_flag(flag, "CVodeGetNumLinSolvSetups")) { return -1; }
  flag = CVodeGetNumLinRhsEvals(cvode_mem, &nf_ls);
  if (check_flag(flag, "CVodeGetNumLinRhsEvals")) { return -1; }
  flag = CVodeGetNumJtimesEvals(cvode_mem, &nJv);
  if (check_flag(flag, "CVodeGetNumJtimesEvals")) { return -1; }

  std::cout << std::fixed << std::setprecision(6)
            << "Final integrator statistics:\n"
            << "  Steps            = " << nst << "\n"
            << "  Error test fails = " << netf << "\n"
            << "  RHS evals        = " << nf << "\n"
            << "  NLS iters        = " << nni << "\n"
            << "  NLS fails        = " << ncfn << "\n"
            << "  LS iters         = " << nli << "\n"
            << "  LS fails         = " << nlcf << "\n"
            << "  LS setups        = " << nsetups << "\n"
            << "  LS RHS evals     = " << nf_ls << "\n"
            << "  Jv products      = " << nJv << "\n"
            << std::endl;

  // Compute average nls iters per step attempt and ls iters per nls iter
  auto avgnli = static_cast<sunrealtype>(nni) / static_cast<sunrealtype>(nst);
  auto avgli  = static_cast<sunrealtype>(nli) / static_cast<sunrealtype>(nni);
  std::cout << "  Avg NLS iters per step    = " << avgnli << "\n";
  std::cout << "  Avg LS iters per NLS iter = " << avgli << "\n";
  std::cout << std::endl;

  // Get preconditioner stats
  if (udata.prec)
  {
    long int npe, nps;
    flag = CVodeGetNumPrecEvals(cvode_mem, &npe);
    if (check_flag(flag, "CVodeGetNumPrecEvals")) { return -1; }
    flag = CVodeGetNumPrecSolves(cvode_mem, &nps);
    if (check_flag(flag, "CVodeGetNumPrecSolves")) { return -1; }

    std::cout << "  Preconditioner setups = " << npe << "\n";
    std::cout << "  Preconditioner solves = " << nps << "\n";
    std::cout << std::endl;
  }

  // Output final error
  flag = SolutionError(t, u, e, udata);
  if (check_flag(flag, "SolutionError")) { return 1; }

  sunrealtype maxerr = N_VMaxNorm(e);

  std::cout << std::scientific
            << std::setprecision(std::numeric_limits<sunrealtype>::digits10)
            << "  Max error = " << maxerr << std::endl;

  // --------------------
  // Clean up and return
  // --------------------

  CVodeFree(&cvode_mem); // Free integrator memory
  SUNLinSolFree(LS);     // Free linear solver
  N_VDestroy(u);         // Free vectors
  N_VDestroy(e);

  return 0;
}

// -----------------------------------------------------------------------------
// Functions called by the integrator
// -----------------------------------------------------------------------------

// f routine to compute the ODE RHS function f(t,y).
int f(sunrealtype t, N_Vector u, N_Vector f, void* user_data)
{
  // Access problem data and set shortcuts
  auto udata    = static_cast<UserData*>(user_data);
  const auto nx = udata->nx;
  const auto ny = udata->ny;
  const auto dx = udata->dx;
  const auto dy = udata->dy;
  const auto kx = udata->kx;
  const auto ky = udata->ky;

  // Access data arrays
  auto uarray = N_VGetArrayPointer(u);
  if (check_ptr(uarray, "N_VGetArrayPointer")) { return -1; }

  auto farray = N_VGetArrayPointer(f);
  if (check_ptr(farray, "N_VGetArrayPointer")) { return -1; }

  // Constants for computing f(t,y)
  const auto cx = kx / (dx * dx);
  const auto cy = ky / (dy * dy);
  const auto cc = -TWO * (cx + cy);

  const auto bx = kx * TWO * PI * PI;
  const auto by = ky * TWO * PI * PI;

  const auto sin_t_cos_t = sin(PI * t) * cos(PI * t);
  const auto cos_sqr_t   = cos(PI * t) * cos(PI * t);

  // Initialize RHS vector to zero (handles boundary conditions)
  N_VConst(ZERO, f);

  // Iterate over domain interior and fill the RHS vector
  for (sunindextype j = 1; j < ny - 1; j++)
  {
    for (sunindextype i = 1; i < nx - 1; i++)
    {
      auto x = i * dx;
      auto y = j * dy;

      auto sin_sqr_x = sin(PI * x) * sin(PI * x);
      auto sin_sqr_y = sin(PI * y) * sin(PI * y);

      auto cos_sqr_x = cos(PI * x) * cos(PI * x);
      auto cos_sqr_y = cos(PI * y) * cos(PI * y);

      // center, north, south, east, and west indices
      auto idx_c = i + j * nx;
      auto idx_n = i + (j + 1) * nx;
      auto idx_s = i + (j - 1) * nx;
      auto idx_e = (i + 1) + j * nx;
      auto idx_w = (i - 1) + j * nx;

      farray[idx_c] = cc * uarray[idx_c] + cx * (uarray[idx_w] + uarray[idx_e]) +
                      cy * (uarray[idx_s] + uarray[idx_n]) -
                      TWO * PI * sin_sqr_x * sin_sqr_y * sin_t_cos_t -
                      bx * (cos_sqr_x - sin_sqr_x) * sin_sqr_y * cos_sqr_t -
                      by * (cos_sqr_y - sin_sqr_y) * sin_sqr_x * cos_sqr_t;
    }
  }

  // Return success
  return 0;
}

// Preconditioner setup routine
int PSetup(sunrealtype t, N_Vector u, N_Vector f, sunbooleantype jok,
           sunbooleantype* jcurPtr, sunrealtype gamma, void* user_data)
{
  // Access problem data
  auto udata = static_cast<UserData*>(user_data);

  // Access data array
  sunrealtype* diag = N_VGetArrayPointer(udata->d);
  if (check_ptr(diag, "N_VGetArrayPointer")) { return -1; }

  // Constants for computing diffusion
  auto cx = udata->kx / (udata->dx * udata->dx);
  auto cy = udata->ky / (udata->dy * udata->dy);
  auto cc = -TWO * (cx + cy);

  // Set all entries of d to the inverse diagonal values of interior
  // (since boundary RHS is 0, set boundary diagonals to the same)
  sunrealtype c = ONE / (ONE - gamma * cc);
  N_VConst(c, udata->d);

  // Return success
  return 0;
}

// Preconditioner solve routine for Pz = r
int PSolve(sunrealtype t, N_Vector u, N_Vector f, N_Vector r, N_Vector z,
           sunrealtype gamma, sunrealtype delta, int lr, void* user_data)
{
  // Access user_data structure
  UserData* udata = static_cast<UserData*>(user_data);

  // Perform Jacobi iteration
  N_VProd(udata->d, r, z);

  // Return success
  return 0;
}

//---- end of file ----
