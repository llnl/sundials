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
 * Header file for extsts advection-diffusion-reaction equation example, see
 * extsts_advection_diffusion_reaction.cpp for more details.
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
#include "arkode/arkode_lsrkstep.h"
#include "arkode/arkode_mristep.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_core.hpp"
#include "sunlinsol/sunlinsol_band.h"
#include "sunmatrix/sunmatrix_band.h"

// Macros for problem constants
#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)

#define NSPECIES 3

#define WIDTH (10 + numeric_limits<sunrealtype>::digits10)

// Macro to access each species at an x location
#define UIDX(i) (NSPECIES * (i))
#define VIDX(i) (NSPECIES * (i) + 1)
#define WIDX(i) (NSPECIES * (i) + 2)

using namespace std;

// -----------------------------------------------------------------------------
// Problem parameters
// -----------------------------------------------------------------------------

struct UserData
{
  // RHS options
  bool reaction = true;
  bool advection = true;

  // Advection and diffusion coefficients
  sunrealtype c = SUN_RCONST(1.0e-2);
  sunrealtype d = SUN_RCONST(1.0e-1);

  // Feed and reaction rates
  sunrealtype A = SUN_RCONST(0.6);
  sunrealtype B = SUN_RCONST(2.0);

  // Stiffness parameter
  sunrealtype eps = SUN_RCONST(1.0e-2);

  // Final simulation time
  sunrealtype tf = SUN_RCONST(3.0);

  // Domain boundaries
  sunrealtype xl = ZERO;
  sunrealtype xu = ONE;

  // Number of nodes
  sunindextype nx = 512;

  // Mesh spacing
  sunrealtype dx = (xu - xl) / (nx - 1);

  // Number of equations
  sunindextype neq = NSPECIES * nx;

  // Temporary workspace vector and matrix
  N_Vector temp_v  = nullptr;
  SUNMatrix temp_J = nullptr;

  // Inner stepper memory
  MRIStepInnerStepper sts_mem = nullptr;

  ~UserData();
};

UserData::~UserData()
{
  if (temp_v)
  {
    N_VDestroy(temp_v);
    temp_v = nullptr;
  }

  if (temp_J)
  {
    SUNMatDestroy(temp_J);
    temp_J = nullptr;
  }
}

// -----------------------------------------------------------------------------
// Problem options
// -----------------------------------------------------------------------------

struct UserOptions
{
  // STS method type: 0 = RKC, 1 = RKL
  int sts_method = 0;

  // MRI method name, or Butcher table name plus flag to build MRI table
  // as an MIS method
  string mri_method = "ARKODE_IMEX_MRI_GARK_ARS222";
  bool build_mri_table = false;

  // Relative and absolute tolerances
  sunrealtype rtol = SUN_RCONST(1.e-4);
  sunrealtype atol = SUN_RCONST(1.e-9);

  // Step size selection (ZERO = adaptive steps)
  sunrealtype fixed_h = ZERO;

  int maxsteps      = 10000; // max steps between outputs
  int ls_setup_freq = 0;     // linear solver setup frequency

  bool calc_error = false;
  bool write_solution = false;

  int output = 1;  // 0 = none, 1 = stats, 2 = disk, 3 = disk with tstop
  int nout   = 10; // number of output times
  ofstream uout;   // output file stream
};

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrators
// -----------------------------------------------------------------------------

// ODE right hand side (RHS) functions
int f_advection(sunrealtype t, N_Vector y, N_Vector f, void* user_data);
int f_diffusion(sunrealtype t, N_Vector y, N_Vector f, void* user_data);
int f_reaction(sunrealtype t, N_Vector y, N_Vector f, void* user_data);

// Jacobian of RHS functions
int J_reaction(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Dominant eigenvalue function (for diffusion operator in LSRKStep)
int diffusion_domeig(sunrealtype t, N_Vector y, N_Vector fn,
                     sunrealtype* lambdaR, sunrealtype* lambdaI,
                     void* user_data, N_Vector temp1, N_Vector temp2,
                     N_Vector temp3);

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// Compute the initial condition
int SetIC(N_Vector y, UserData& udata);

// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------

// Check function return flag
static int check_flag(int flag, const string funcname)
{
  if (flag < 0)
  {
    cerr << "ERROR: " << funcname << " returned " << flag << endl;
    return 1;
  }
  return 0;
}

// Check if a function returned a NULL pointer
static int check_ptr(void* ptr, const string funcname)
{
  if (ptr) { return 0; }
  cerr << "ERROR: " << funcname << " returned NULL" << endl;
  return 1;
}

// Print command line options
static void InputHelp()
{
  cout << endl;
  cout << "Command line options:" << endl;
  cout << "  --no-advection           : disable advection\n";
  cout << "  --no-reaction            : disable reaction\n";
  cout << "  --c <real>               : advection coefficient\n";
  cout << "  --d <real>               : diffusion coefficient\n";
  cout << "  --A <real>               : species A concentration\n";
  cout << "  --B <real>               : species B concentration\n";
  cout << "  --eps <real>             : stiffness parameter\n";
  cout << "  --tf <real>              : final time\n";
  cout << "  --xl <real>              : domain lower boundary\n";
  cout << "  --xu <real>              : domain upper boundary\n";
  cout << "  --nx <int>               : number of mesh points\n";
  cout << "  --sts_method <int>       : STS method type (0=RKC, 1=RKL)\n";
  cout << "  --mri_method <string>    : MRI method or Butcher table name\n";
  cout << "  --build_mri_table        : build MRI table from Butcher table\n";
  cout << "  --rtol <real>            : relative tolerance\n";
  cout << "  --atol <real>            : absolute tolerance\n";
  cout << "  --fixed_h <real>         : fixed step size\n";
  cout << "  --lssetupfreq <int>      : LS setup frequency\n";
  cout << "  --maxsteps <int>         : max steps between outputs\n";
  cout << "  --output <int>           : output level\n";
  cout << "  --nout <int>             : number of outputs\n";
  cout << "  --help                   : print options and exit\n";
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

#if defined(SUNDIALS_INT64_T)
inline void find_arg(vector<string>& args, const string key, sunindextype& dest)
{
  auto it = find(args.begin(), args.end(), key);
  if (it != args.end())
  {
    dest = stoll(*(it + 1));
    args.erase(it, it + 2);
  }
}
#endif

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

inline void find_arg(vector<string>& args, const string key, string& dest)
{
  auto it = find(args.cbegin(), args.cend(), key);
  if (it != args.end())
  {
    dest = std::move(*(it + 1));
    args.erase(it, it + 2);
  }
}

static int ReadInputs(vector<string>& args, UserData& udata, UserOptions& uopts,
                      SUNContext ctx)
{
  if (find(args.begin(), args.end(), "--help") != args.end())
  {
    InputHelp();
    return 1;
  }

  // Problem parameters
  find_arg(args, "--no-advection", udata.advection, false);
  find_arg(args, "--no-reaction", udata.reaction, false);
  find_arg(args, "--c", udata.c);
  find_arg(args, "--d", udata.d);
  find_arg(args, "--A", udata.A);
  find_arg(args, "--B", udata.B);
  find_arg(args, "--eps", udata.eps);
  find_arg(args, "--tf", udata.tf);
  find_arg(args, "--xl", udata.xl);
  find_arg(args, "--xu", udata.xu);
  find_arg(args, "--nx", udata.nx);

  // Integrator options
  find_arg(args, "--sts_method", uopts.sts_method);
  find_arg(args, "--mri_method", uopts.mri_method);
  find_arg(args, "--build_mri_table", uopts.build_mri_table);
  find_arg(args, "--rtol", uopts.rtol);
  find_arg(args, "--atol", uopts.atol);
  find_arg(args, "--fixed_h", uopts.fixed_h);
  find_arg(args, "--lssetupfreq", uopts.ls_setup_freq);
  find_arg(args, "--maxsteps", uopts.maxsteps);
  find_arg(args, "--write_solution", uopts.write_solution);
  find_arg(args, "--output", uopts.output);
  find_arg(args, "--nout", uopts.nout);

  // Recompute mesh spacing and total number of nodes
  udata.dx  = (udata.xu - udata.xl) / (udata.nx - 1);
  udata.neq = NSPECIES * udata.nx;

  // Input checks
  if (!udata.reaction && !udata.advection)
  {
    cerr << "ERROR: Invalid problem configuration" << endl;
    return -1;
  }

  // Ensure that build_mri_table is only set for non-ImEx problem
  if (uopts.build_mri_table && !udata.advection && !udata.reaction)
  {
    cerr << "ERROR: Cannot build MRI table for non-ImEx problem" << endl;
    return -1;
  }

  return 0;
}

// Print user data
static int PrintSetup(UserData& udata, UserOptions& uopts)
{
  cout << endl;
  cout << "Problem parameters and options:" << endl;
  cout << " --------------------------------- " << endl;
  cout << "  c                = " << udata.c << endl;
  cout << "  d                = " << udata.d << endl;
  cout << "  A                = " << udata.A << endl;
  cout << "  B                = " << udata.B << endl;
  cout << "  eps              = " << udata.eps << endl;
  cout << " --------------------------------- " << endl;
  cout << "  tf               = " << udata.tf << endl;
  cout << "  xl               = " << udata.xl << endl;
  cout << "  xu               = " << udata.xu << endl;
  cout << "  nx               = " << udata.nx << endl;
  cout << "  dx               = " << udata.dx << endl;
  cout << " --------------------------------- " << endl;

  cout << "  integrator       = ExtSTS" << endl;
  if (udata.advection) { cout << "  advection        = Explicit" << endl; }
  else { cout << "  advection        = OFF" << endl; }
  if (udata.reaction) {  cout << "  reaction         = Implicit" << endl; }
  else { cout << "  reaction         = OFF" << endl; }
  cout << "  diffusion        = Explicit" << endl;

  cout << "  rtol             = " << uopts.rtol << endl;
  cout << "  atol             = " << uopts.atol << endl;
  cout << "  fixed h          = " << uopts.fixed_h << endl;
  cout << "  ls setup freq    = " << uopts.ls_setup_freq << endl;
  cout << " --------------------------------- " << endl;
  if (uopts.build_mri_table)
  {
    cout << "  MRI method constructed from Butcher table " << uopts.mri_method << endl;
  }
  else
  {
    cout << "  MRI method       = " << uopts.mri_method << endl;
  }
  if (uopts.sts_method == 0)
  { cout << "  STS method       = RKC" << endl; }
  else
  { cout << "  STS method       = RKL" << endl; }
  cout << " --------------------------------- " << endl;
  cout << "  output           = " << uopts.output << endl;
  cout << " --------------------------------- " << endl;
  cout << endl;

  return 0;
}

// Initialize output
static int OpenOutput(UserData& udata, UserOptions& uopts)
{
  // Header for status output
  if (uopts.output)
  {
    cout << scientific;
    cout << setprecision(numeric_limits<sunrealtype>::digits10);
    cout << "          t           ";
    cout << "          ||y||_rms      ";
    if (uopts.calc_error)
    {
      cout << "   ||yerr||_rms";
    }
    cout << endl;
    cout << " ---------------------";
    if (uopts.calc_error)
    {
      cout << "---------------";
    }
    cout << "-------------------------" << endl;
  }

  // Open output stream and output problem information
  if (uopts.output >= 2)
  {
    // Open output stream
    stringstream fname;
    fname << "advection_diffusion_reaction.out";
    uopts.uout.open(fname.str());

    uopts.uout << scientific;
    uopts.uout << setprecision(numeric_limits<sunrealtype>::digits10);
    uopts.uout << "# title Advection-Diffusion-Reaction (Brusselator)" << endl;
    uopts.uout << "# nvar 3" << endl;
    uopts.uout << "# vars u v w" << endl;
    uopts.uout << "# nt " << uopts.nout + 1 << endl;
    uopts.uout << "# nx " << udata.nx << endl;
    uopts.uout << "# xl " << udata.xl << endl;
    uopts.uout << "# xu " << udata.xu << endl;
  }

  return 0;
}

// Write output
static int WriteOutput(sunrealtype t, N_Vector y, UserData& udata,
                       UserOptions& uopts)
{
  if (uopts.output)
  {
    // Compute rms norm of the state
    sunrealtype urms = sqrt(N_VDotProd(y, y) / udata.nx);
    cout << setw(22) << t << setw(25) << urms << endl;

    // Write solution to disk
    if (uopts.output >= 2)
    {
      sunrealtype* ydata = N_VGetArrayPointer(y);
      if (check_ptr(ydata, "N_VGetArrayPointer")) { return -1; }

      uopts.uout << t;
      for (sunindextype i = 0; i < udata.nx; i++)
      {
        uopts.uout << setw(WIDTH) << ydata[UIDX(i)];
        uopts.uout << setw(WIDTH) << ydata[VIDX(i)];
        uopts.uout << setw(WIDTH) << ydata[WIDX(i)];
      }
      uopts.uout << endl;
    }
  }

  return 0;
}

// Finalize output
static int CloseOutput(UserOptions& uopts)
{
  // Footer for status output
  if (uopts.output)
  {
    cout << " ---------------------";
    if (uopts.calc_error)
    {
      cout << "---------------";
    }
    cout << "-------------------------" << endl;
    cout << endl;
  }

  // Close output streams
  if (uopts.output >= 2) { uopts.uout.close(); }

  return 0;
}

//---- end of file ----
