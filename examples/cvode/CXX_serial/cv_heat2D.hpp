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
 * Header file for 2D heat equation example problem.
 * See cv_heat2D.cpp for more information.
 * ---------------------------------------------------------------------------*/

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

// SUNDIALS types
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

// Common utility functions
#include <example_utilities.hpp>

// Macros for problem constants
#define PI   SUN_RCONST(3.141592653589793238462643383279502884197169)
#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)

// -----------------------------------------------------------------------------
// User data structure
// -----------------------------------------------------------------------------

struct UserData
{
  // Diffusion coefficients in the x and y directions
  sunrealtype kx = ONE;
  sunrealtype ky = ONE;

  // Final time
  sunrealtype tf = ONE;

  // Upper bounds in x and y directions
  sunrealtype xu = ONE;
  sunrealtype yu = ONE;

  // Number of nodes in the x and y directions
  sunindextype nx = 32;
  sunindextype ny = 32;

  // Total number of nodes
  sunindextype nodes = nx * ny;

  // Mesh spacing in the x and y directions
  sunrealtype dx = xu / (nx - 1);
  sunrealtype dy = yu / (ny - 1);

  // Integrator settings
  sunrealtype rtol = SUN_RCONST(1.0e-4); // relative tolerance
  sunrealtype atol = SUN_RCONST(1.0e-8); // absolute tolerance
  int maxsteps     = 0;                  // max number of steps between outputs

  // Linear solver and preconditioner settings
  bool pcg           = true; // use PCG (true) or GMRES (false)
  bool prec          = true; // preconditioner on/off
  int liniters       = 20;   // number of linear iterations
  int msbp           = 0; // max number of steps between preconditioner setups
  sunrealtype epslin = ZERO; // linear solver tolerance factor

  // Inverse of Jacobian diagonal for preconditioner
  N_Vector d = nullptr;

  // Output variables
  bool output = false; // write solution to disk
  int nout    = 20;    // number of output times
  std::ofstream uout;  // output file stream
  std::ofstream eout;  // error file stream

  // Destructor
  ~UserData();
};

// Free memory allocated within Userdata
UserData::~UserData()
{
  // Free preconditioner data
  if (d)
  {
    N_VDestroy(d);
    d = nullptr;
  }
}

// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Compute the exact solution
static int Solution(sunrealtype t, N_Vector u, UserData& udata)
{
  auto uarray = N_VGetArrayPointer(u);
  if (check_ptr(uarray, "N_VGetArrayPointer")) { return -1; }

  // Initialize u to one (handles boundary conditions)
  N_VConst(ONE, u);

  // Compute the true solution
  auto cos_sqr_t = cos(PI * t) * cos(PI * t);

  for (sunindextype j = 1; j < udata.ny - 1; j++)
  {
    for (sunindextype i = 1; i < udata.nx - 1; i++)
    {
      auto x = i * udata.dx;
      auto y = j * udata.dy;

      auto sin_sqr_x = sin(PI * x) * sin(PI * x);
      auto sin_sqr_y = sin(PI * y) * sin(PI * y);

      auto idx    = i + j * udata.nx;
      uarray[idx] = sin_sqr_x * sin_sqr_y * cos_sqr_t + ONE;
    }
  }

  return 0;
}

// Compute the solution error
static int SolutionError(sunrealtype t, N_Vector u, N_Vector e, UserData& udata)
{
  // Compute true solution
  auto flag = Solution(t, e, udata);
  if (flag != 0) { return -1; }

  // Compute absolute error
  N_VLinearSum(ONE, u, -ONE, e, e);
  N_VAbs(e, e);

  return 0;
}

// Print command line options
static void InputHelp()
{
  std::cout << std::endl
            << "Command line options:\n"
            << "  --nx <nx>          : number of x mesh points\n"
            << "  --nx <nx>          : number of y mesh points\n"
            << "  --xu <xu>          : x upper bound\n"
            << "  --yu <yu>          : y upper bound\n"
            << "  --kx <kx>          : x diffusion coefficient\n"
            << "  --kx <ky>          : y diffusion coefficient\n"
            << "  --tf <time>        : final time\n"
            << "  --rtol <rtol>      : relative tolerance\n"
            << "  --atol <atol>      : absolute tolerance\n"
            << "  --gmres            : use GMRES linear solver\n"
            << "  --noprec           : disable preconditioner\n"
            << "  --liniters <iters> : max number of iterations\n"
            << "  --epslin <factor>  : linear tolerance factor\n"
            << "  --msbp <steps>     : max steps between prec setups\n"
            << "  --output           : write solution to disk\n"
            << "  --nout <nout>      : number of outputs\n"
            << "  --maxsteps <steps> : max steps between outputs\n"
            << "  --help             : print this message and exit\n";
}

// Read command line inputs
static int ReadInputs(std::vector<std::string>& args, UserData& udata)
{
  if (find(args.begin(), args.end(), "--help") != args.end())
  {
    InputHelp();
    return 1;
  }

  find_arg(args, "--nx", udata.nx);
  find_arg(args, "--ny", udata.ny);
  find_arg(args, "--xu", udata.xu);
  find_arg(args, "--yu", udata.yu);
  find_arg(args, "--kx", udata.kx);
  find_arg(args, "--ky", udata.ky);
  find_arg(args, "--tf", udata.tf);
  find_arg(args, "--rtol", udata.rtol);
  find_arg(args, "--atol", udata.atol);
  find_arg(args, "--gmres", udata.pcg, false);
  find_arg(args, "--noprec", udata.prec, false);
  find_arg(args, "--liniters", udata.liniters);
  find_arg(args, "--msbp", udata.msbp);
  find_arg(args, "--epslin", udata.epslin);
  find_arg(args, "--output", udata.output);
  find_arg(args, "--nout", udata.nout);
  find_arg(args, "--maxsteps", udata.maxsteps);

  // Recompute total number of nodes and mesh spacing
  udata.nodes = udata.nx * udata.ny;
  udata.dx    = udata.xu / (udata.nx - 1);
  udata.dy    = udata.yu / (udata.ny - 1);

  // Return success
  return 0;
}

// Print user data
static void PrintUserData(UserData& udata)
{
  std::cout << std::endl
            << "2D Heat problem:\n"
            << " ----------------------------\n"
            << "  kx        = " << udata.kx << "\n"
            << "  ky        = " << udata.ky << "\n"
            << "  tf        = " << udata.tf << "\n"
            << "  xu        = " << udata.xu << "\n"
            << "  yu        = " << udata.yu << "\n"
            << "  nx        = " << udata.nx << "\n"
            << "  ny        = " << udata.ny << "\n"
            << "  dx        = " << udata.dx << "\n"
            << "  dy        = " << udata.dy << "\n"
            << " ----------------------------\n";
  if (udata.pcg) { std::cout << "  linear solver  = PCG\n"; }
  else { std::cout << "  linear solver  = GMRES\n"; }
  std::cout << "  rtol      = " << udata.rtol << "\n"
            << "  atol      = " << udata.atol << "\n"
            << " ----------------------------\n"
            << "  lin iters = " << udata.liniters << "\n"
            << "  eps lin   = " << udata.epslin << "\n"
            << "  prec      = " << udata.prec << "\n"
            << "  msbp      = " << udata.msbp << "\n"
            << " ----------------------------\n"
            << "  output    = " << udata.output << "\n"
            << " ----------------------------\n"
            << std::endl;
}

// Initialize output
static int OpenOutput(UserData& udata)
{
  // Header for status output
  std::cout << std::scientific
            << std::setprecision(std::numeric_limits<sunrealtype>::digits10)
            << "          t                     ||u||_rms      "
            << "          max error\n"
            << " ----------------------------------------------"
            << "-------------------------\n";

  // Output problem information and open output streams
  if (udata.output)
  {
    // Each processor outputs subdomain information
    std::ofstream dout;
    dout.open("heat2d_info.txt");
    dout << "xu  " << udata.xu << std::endl;
    dout << "yu  " << udata.yu << std::endl;
    dout << "nx  " << udata.nx << std::endl;
    dout << "ny  " << udata.ny << std::endl;
    dout << "nt  " << udata.nout + 1 << std::endl;
    dout.close();

    // Open output streams for solution and error
    udata.uout.open("heat2d_solution.txt");
    udata.uout << std::scientific
               << std::setprecision(std::numeric_limits<sunrealtype>::digits10);

    udata.eout.open("heat2d_error.txt");
    udata.eout << std::scientific
               << std::setprecision(std::numeric_limits<sunrealtype>::digits10);
  }

  return 0;
}

// Write output
static int WriteOutput(sunrealtype t, N_Vector u, N_Vector e, UserData& udata)
{
  // Compute the error
  auto flag = SolutionError(t, u, e, udata);
  if (check_flag(flag, "SolutionError")) { return 1; }

  // Compute max error
  sunrealtype max = N_VMaxNorm(e);

  // Compute rms norm of the state
  sunrealtype urms = sqrt(N_VDotProd(u, u) / udata.nx / udata.ny);

  // Output current status
  std::cout << std::setw(22) << t << std::setw(25) << urms << std::setw(25)
            << max << std::endl;

  // Write solution and error to disk
  if (udata.output)
  {
    sunrealtype* uarray = N_VGetArrayPointer(u);
    if (check_ptr(uarray, "N_VGetArrayPointer")) { return -1; }

    udata.uout << t << " ";
    for (sunindextype i = 0; i < udata.nodes; i++)
    {
      udata.uout << uarray[i] << " ";
    }
    udata.uout << std::endl;

    // Output error to disk
    sunrealtype* earray = N_VGetArrayPointer(e);
    if (check_ptr(earray, "N_VGetArrayPointer")) { return -1; }

    udata.eout << t << " ";
    for (sunindextype i = 0; i < udata.nodes; i++)
    {
      udata.eout << earray[i] << " ";
    }
    udata.eout << std::endl;
  }

  return 0;
}

// Finalize output
static int CloseOutput(UserData& udata)
{
  // Footer for status output
  std::cout << " ----------------------------------------------"
            << "-------------------------\n\n";

  if (udata.output)
  {
    // Close output streams
    udata.uout.close();
    udata.eout.close();
  }

  return 0;
}
