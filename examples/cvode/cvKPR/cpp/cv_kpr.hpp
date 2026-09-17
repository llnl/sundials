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
 * Kvaerno-Prothero-Robinson ODE test problem, see .cpp file for details
 * ---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

// SUNDIALS types
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

// Common utility functions
#include <example_utilities.hpp>

// Macros for problem constants
#define ZERO   SUN_RCONST(0.0)
#define HALF   SUN_RCONST(0.5)
#define ONE    SUN_RCONST(1.0)
#define TWO    SUN_RCONST(2.0)
#define TWENTY SUN_RCONST(20.0)

// -----------------------------------------------------------------------------
// Problem options
// -----------------------------------------------------------------------------

struct Options
{
  // Relative and absolute tolerances
  sunrealtype rtol = SUN_RCONST(1.0e-6);
  sunrealtype atol = SUN_RCONST(1.0e-10);

  // Finite difference Jacobian
  bool fd_jac = false;

  // Output options
  sunrealtype dtout = ONE; // output interval
  int nout          = 10;  // number of outputs
};

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// Compute r(t)
inline sunrealtype r(sunrealtype t) { return HALF * cos(t); }

// Compute the derivative of r(t)
inline sunrealtype rdot(sunrealtype t) { return -HALF * sin(t); }

// Compute s(t)
inline sunrealtype s(sunrealtype t) { return cos(TWENTY * t); }

// Compute the derivative of s(t)
inline sunrealtype sdot(sunrealtype t) { return -TWENTY * sin(TWENTY * t); }

// Compute the true solution
inline int true_sol(sunrealtype t, sunrealtype* u, sunrealtype* v)
{
  *u = sqrt(ONE + r(t));
  *v = sqrt(TWO + s(t));

  return 0;
}

// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Print command line options
static void InputHelp()
{
  std::cout << std::endl;
  std::cout << "Command line options:" << std::endl;
  std::cout << "  --help         : print options and exit\n";
  std::cout << "  --rtol         : relative tolerance\n";
  std::cout << "  --atol         : absolute tolerance\n";
  std::cout << "  --fdjac        : finite-difference Jacobian\n";
  std::cout << "  --dtout        : output interval\n";
  std::cout << "  --nout         : number of outputs\n";
}

static int ReadInputs(std::vector<std::string>& args, Options& opts,
                      SUNContext ctx)
{
  if (find(args.begin(), args.end(), "--help") != args.end())
  {
    InputHelp();
    return 1;
  }

  // Problem options
  find_arg(args, "--rtol", opts.rtol);
  find_arg(args, "--atol", opts.atol);
  find_arg(args, "--fdjac", opts.fd_jac);
  find_arg(args, "--dtout", opts.dtout);
  find_arg(args, "--nout", opts.nout);

  return 0;
}

/*---- end of file ----*/
