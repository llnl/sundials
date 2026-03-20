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

// Include desired integrators, vectors, linear solvers, and nonlinear solvers
#include "ida/ida.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_core.hpp"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

// Macros for problem constants
#define ZERO   SUN_RCONST(0.0)
#define HALF   SUN_RCONST(0.5)
#define ONE    SUN_RCONST(1.0)
#define TWO    SUN_RCONST(2.0)
#define TWENTY SUN_RCONST(20.0)

using namespace std;

// -----------------------------------------------------------------------------
// Test options
// -----------------------------------------------------------------------------

struct TestOptions
{
  // Relative and absolute tolerances
  sunrealtype rtol = SUN_RCONST(1.0e-6);
  sunrealtype atol = SUN_RCONST(1.0e-10);

  // Fixed step size eta bounds (use defaults = 0.0 and 1.5)
  sunrealtype eta_min_fx = -ONE;
  sunrealtype eta_max_fx = -ONE;

  // Max and min eta bounds on a general step (use default = 10)
  sunrealtype eta_max = -ONE;
  sunrealtype eta_min = -ONE;

  // Min eta bound after an error test fail (use default = 0.25)
  sunrealtype eta_min_ef = -ONE;

  // Eta value on a nonlinear solver convergence failure (use default = 0.25)
  sunrealtype eta_cf = -ONE;

  // Parameter for if a change in c_j needs a call lsetup (use defaults 0.25)
  sunrealtype dcj = -ONE;

  // Output options
  sunrealtype dtout = ONE; // output interval
  int nout          = 10;  // number of outputs
};

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrators
// -----------------------------------------------------------------------------

// DAE residual function
int res(sunrealtype t, N_Vector y, N_Vector yp, N_Vector rr, void* user_data);

// Jacobian of RHS function
int J(sunrealtype t, sunrealtype cj, N_Vector y, N_Vector yp, N_Vector rr,
      SUNMatrix J, void* user_data, N_Vector tempv1, N_Vector tempv2,
      N_Vector tempv3);

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// Compute r(t)
static sunrealtype r(sunrealtype t) { return HALF * cos(t); }

// Compute the derivative of r(t)
static sunrealtype rdot(sunrealtype t) { return -HALF * sin(t); }

// Compute s(t)
static sunrealtype s(sunrealtype t) { return cos(TWENTY * t); }

// Compute the derivative of s(t)
static sunrealtype sdot(sunrealtype t) { return -TWENTY * sin(TWENTY * t); }

// Compute the true solution
static int true_sol(sunrealtype t, sunrealtype* u, sunrealtype* v)
{
  *u = sqrt(ONE + r(t));
  *v = sqrt(TWO + s(t));

  return 0;
}

// Compute the true solution derivative
static int true_sol_p(sunrealtype t, sunrealtype* up, sunrealtype* vp)
{
  *up = rdot(t) / (TWO * sqrt(ONE + r(t)));
  *vp = sdot(t) / (TWO * sqrt(TWO + s(t)));

  return 0;
}

// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Check function return flag
static int check_flag(int flag, const string funcname)
{
  if (!flag) { return 0; }
  if (flag < 0) { cerr << "ERROR: "; }
  cerr << funcname << " returned " << flag << endl;
  return 1;
}

// Check if a function returned a NULL pointer
static int check_ptr(void* ptr, const string funcname)
{
  if (ptr) { return 0; }
  cerr << "ERROR: " << funcname << " returned NULL" << endl;
  return 1;
}

inline void find_arg(vector<string>& args, const string key, sunrealtype& dest)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
#if defined(SUNDIALS_SINGLE_PRECISION)
    dest = stof(*(it + 1));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    dest = stod(*(it + 1));
#elif defined(SUNDIALS_EXTENDED_PRECISION)
    dest = stold(*(it + 1));
#endif
    args.erase(it, it + 2);
  }
}

inline void find_arg(vector<string>& args, const string key, long int& dest)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stoll(*(it + 1));
    args.erase(it, it + 2);
  }
}

inline void find_arg(vector<string>& args, const string key, int& dest)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stoi(*(it + 1));
    args.erase(it, it + 2);
  }
}

inline void find_arg(vector<string>& args, const string key, bool& dest,
                     bool store = true)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = store;
    args.erase(it);
  }
}

// Print command line options
static void InputHelp()
{
  cout << endl;
  cout << "Command line options:" << endl;
  cout << "  --help        : print options and exit\n";
  cout << "  --rtol        : relative tolerance\n";
  cout << "  --atol        : absolute tolerance\n";
  cout << "  --eta_min_fx  : fixed step lower bound\n";
  cout << "  --eta_max_fx  : fixed step upper bound\n";
  cout << "  --eta_max     : general step max eta\n";
  cout << "  --eta_min     : general step min eta\n";
  cout << "  --eta_min_ef  : min eta on err test fail\n";
  cout << "  --eta_cf      : eta on an error test fail\n";
  cout << "  --dcj         : change in cj triggers lsetup\n";
  cout << "  --dtout       : output interval\n";
  cout << "  --nout        : number of outputs\n";
}

static int ReadInputs(vector<string>& args, TestOptions& opts, SUNContext ctx)
{
  if (find(args.begin(), args.end(), "--help") != args.end())
  {
    InputHelp();
    return 1;
  }

  // Test options
  find_arg(args, "--rtol", opts.rtol);
  find_arg(args, "--atol", opts.atol);
  find_arg(args, "--eta_min_fx", opts.eta_min_fx);
  find_arg(args, "--eta_max_fx", opts.eta_max_fx);
  find_arg(args, "--eta_max", opts.eta_max);
  find_arg(args, "--eta_min", opts.eta_min);
  find_arg(args, "--eta_min_ef", opts.eta_min_ef);
  find_arg(args, "--eta_cf", opts.eta_cf);
  find_arg(args, "--dcj", opts.dcj);
  find_arg(args, "--dtout", opts.dtout);
  find_arg(args, "--nout", opts.nout);

  return 0;
}

/*---- end of file ----*/
