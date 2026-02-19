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
 * Header file for ARKODE advection-diffusion-reaction equation example, see
 * ark_advection_diffusion_reaction.cpp for more details.
 * ---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

// Include desired integrators, vectors, linear solvers, and nonlinear solvers
#include "arkode/arkode_arkstep.h"
#include "arkode/arkode_erkstep.h"
#include "arkode/arkode_mristep.h"
#include "cvode/cvode.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_core.hpp"
#include "sunlinsol/sunlinsol_band.h"
#include "sunmatrix/sunmatrix_band.h"

// Macros for problem constants
#define PI   SUN_RCONST(3.141592653589793238462643383279502884197169)
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
  bool diffusion = true;
  bool advection = true;

  // Splitting option
  int splitting = 3;

  // Advection and diffusion coefficients
  sunrealtype c = SUN_RCONST(1.0e-3);
  sunrealtype d = SUN_RCONST(1.0e-2);

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
  MRIStepInnerStepper fast_mem = nullptr;

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
  // Integration method (0 = ERK, 1 = ARK, 2 = MRIARK, 3 = MRICVODE)
  int integrator = 1;

  // Method order
  int order      = 3;
  int order_fast = 3;
  bool ark_dirk  = false;

  // Relative and absolute tolerances
  sunrealtype rtol      = SUN_RCONST(1.e-4);
  sunrealtype atol      = SUN_RCONST(1.e-9);
  sunrealtype rtol_fast = SUN_RCONST(1.e-4);
  sunrealtype atol_fast = SUN_RCONST(1.e-9);

  // Step size selection (ZERO = adaptive steps)
  sunrealtype fixed_h      = ZERO;
  sunrealtype fixed_h_fast = ZERO;

  // First step growth factor
  sunrealtype etamx1_fast = ZERO;

  int maxsteps      = 10000; // max steps between outputs
  string controller = "I";   // step size adaptivity method
  int predictor     = 0;     // predictor for nonlinear systems
  int ls_setup_freq = 0;     // linear solver setup frequency

  string controller_fast = "I"; // fast step size adaptivity method
  int predictor_fast     = 0;   // predictor for fast nonlinear systems
  int ls_setup_freq_fast = 0;   // fast linear solver setup frequency

  bool linear = false; // signal that the problem is linearly implicit

  // save and reuse prior fast time step sizes
  bool save_hinit = false;
  bool save_hcur  = false;

  sunrealtype hcur_factor = SUN_RCONST(0.7);

  int output = 1;  // 0 = none, 1 = stats, 2 = disk, 3 = disk with tstop
  int nout   = 10; // number of output times
  ofstream uout;   // output file stream
};

// -----------------------------------------------------------------------------
// Custom inner stepper content and functions
// -----------------------------------------------------------------------------

struct CVodeInnerStepperContent
{
  void* cvode_mem = nullptr; // CVODE memory structure
  void* user_data = nullptr; // user data pointer

  // saved step sizes
  bool save_hinit = false;
  bool save_hcur  = false;

  sunrealtype hcur_factor = ONE;

  sunrealtype hinit = ZERO; // initial step size
  sunrealtype hcur  = ZERO; // current step size

  // saved integrator stats
  long int nst     = 0; // time steps
  long int netf    = 0; // error test fails
  long int nfe     = 0; // rhs evals
  long int nni     = 0; // nonlinear iterations
  long int nncf    = 0; // nonlinear convergence failures
  long int nsetups = 0; // linear solver setups
  long int nje     = 0; // Jacobian evals
};

int CVodeInnerStepper_Evolve(MRIStepInnerStepper fast_mem, sunrealtype t0,
                             sunrealtype tout, N_Vector y);

int CVodeInnerStepper_FullRhs(MRIStepInnerStepper fast_mem, sunrealtype t,
                              N_Vector y, N_Vector f, int mode);

int CVodeInnerStepper_Reset(MRIStepInnerStepper fast_mem, sunrealtype tR,
                            N_Vector yR);

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrators
// -----------------------------------------------------------------------------

// ODE right hand side (RHS) functions
int f_advection(sunrealtype t, N_Vector y, N_Vector f, void* user_data);
int f_diffusion(sunrealtype t, N_Vector y, N_Vector f, void* user_data);
int f_reaction(sunrealtype t, N_Vector y, N_Vector f, void* user_data);

int f_adv_diff(sunrealtype t, N_Vector y, N_Vector f, void* user_data);
int f_adv_react(sunrealtype t, N_Vector y, N_Vector f, void* user_data);
int f_diff_react(sunrealtype t, N_Vector y, N_Vector f, void* user_data);
int f_adv_diff_react(sunrealtype t, N_Vector y, N_Vector f, void* user_data);

int f_react_forcing(sunrealtype t, N_Vector y, N_Vector f, void* user_data);

// Jacobian of RHS functions
int J_advection(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int J_diffusion(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int J_reaction(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int J_adv_diff(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int J_adv_react(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int J_diff_react(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                 void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int J_adv_diff_react(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                     void* user_data, N_Vector tmp1, N_Vector tmp2,
                     N_Vector tmp3);

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// Integrator setup functions
int SetupERK(SUNContext ctx, UserData& udata, UserOptions& uopts, N_Vector y,
             void** arkode_mem);

int SetupARK(SUNContext ctx, UserData& udata, UserOptions& uopts, N_Vector y,
             SUNMatrix* A, SUNLinearSolver* LS, void** arkode_mem);

int SetupMRIARK(SUNContext ctx, UserData& udata, UserOptions& uopts, N_Vector y,
                SUNMatrix* A, SUNLinearSolver* LS, SUNMatrix* A_fast,
                SUNLinearSolver* LS_fast, MRIStepInnerStepper* fast_mem,
                void** arkode_mem);

int SetupMRICVODE(SUNContext ctx, UserData& udata, UserOptions& uopts,
                  N_Vector y, SUNMatrix* A, SUNLinearSolver* LS,
                  SUNMatrix* A_fast, SUNLinearSolver* LS_fast,
                  MRIStepInnerStepper* fast_mem, void** arkode_mem);

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

// Print ERK integrator statistics
static int OutputStatsERK(void* arkode_mem, UserData& udata)
{
  int flag;

  // Get integrator and solver stats
  long int nst, nst_a, netf, nfe;
  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  if (check_flag(flag, "ARKodeGetNumSteps")) { return -1; }
  flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  if (check_flag(flag, "ARKodeGetNumStepAttempts")) { return -1; }
  flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
  if (check_flag(flag, "ARKodeGetNumErrTestFails")) { return -1; }
  flag = ARKodeGetNumRhsEvals(arkode_mem, 0, &nfe);
  if (check_flag(flag, "ARKodeGetNumRhsEvals")) { return -1; }

  cout << "  Steps            = " << nst << endl;
  cout << "  Step attempts    = " << nst_a << endl;
  cout << "  Error test fails = " << netf << endl;
  cout << "  RHS evals        = " << nfe << endl;

  return 0;
}

// Print ARK integrator statistics
static int OutputStatsARK(void* arkode_mem, UserData& udata)
{
  int flag;

  // Get integrator and solver stats
  long int nst, nst_a, netf, nfe, nfi;
  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  if (check_flag(flag, "ARKodeGetNumSteps")) { return -1; }
  flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  if (check_flag(flag, "ARKodeGetNumStepAttempts")) { return -1; }
  flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
  if (check_flag(flag, "ARKodeGetNumErrTestFails")) { return -1; }
  flag = ARKodeGetNumRhsEvals(arkode_mem, 0, &nfe);
  if (check_flag(flag, "ARKodeGetNumRhsEvals")) { return -1; }
  flag = ARKodeGetNumRhsEvals(arkode_mem, 1, &nfi);
  if (check_flag(flag, "ARKodeGetNumRhsEvals")) { return -1; }

  cout << fixed << setprecision(6);
  cout << "  Steps              = " << nst << endl;
  cout << "  Step attempts      = " << nst_a << endl;
  cout << "  Error test fails   = " << netf << endl;
  cout << "  Explicit RHS evals = " << nfe << endl;
  cout << "  Implicit RHS evals = " << nfi << endl;

  if (udata.splitting)
  {
    long int nni, ncfn;
    flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
    if (check_flag(flag, "ARKodeGetNumNonlinSolvIters")) { return -1; }
    flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
    if (check_flag(flag, "ARKodeGetNumNonlinSolvConvFails")) { return -1; }

    long int nsetups, nje;
    flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
    if (check_flag(flag, "ARKodeGetNumLinSolvSetups")) { return -1; }
    flag = ARKodeGetNumJacEvals(arkode_mem, &nje);
    if (check_flag(flag, "ARKodeGetNumJacEvals")) { return -1; }

    cout << "  NLS iters          = " << nni << endl;
    cout << "  NLS fails          = " << ncfn << endl;
    cout << "  LS setups          = " << nsetups << endl;
    cout << "  J evals            = " << nje << endl;
    cout << endl;

    sunrealtype avgnli = (sunrealtype)nni / (sunrealtype)nst_a;
    sunrealtype avgls  = (sunrealtype)nsetups / (sunrealtype)nni;
    cout << "  Avg NLS iters per step attempt = " << avgnli << endl;
    cout << "  Avg LS setups per NLS iter     = " << avgls << endl;
  }
  cout << endl;

  return 0;
}

// Print MRI integrator statistics
static int OutputStatsMRIARK(void* arkode_mem, MRIStepInnerStepper fast_mem,
                             UserData& udata)
{
  int flag;

  // Get slow integrator and solver stats
  long int nst, nst_a, netf, nfe, nfi;
  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  if (check_flag(flag, "ARKodeGetNumSteps")) { return -1; }
  flag = ARKodeGetNumRhsEvals(arkode_mem, 0, &nfe);
  if (check_flag(flag, "ARKodeGetNumRhsEvals")) { return -1; }
  flag = ARKodeGetNumRhsEvals(arkode_mem, 1, &nfi);
  if (check_flag(flag, "ARKodeGetNumRhsEvals")) { return -1; }

  cout << fixed << setprecision(6);
  cout << endl << "Slow Integrator:" << endl;
  cout << "  Steps                   = " << nst << endl;
  cout << "  Slow explicit RHS evals = " << nfe << endl;
  cout << "  Slow implicit RHS evals = " << nfi << endl;

  if (udata.diffusion)
  {
    long int nni, ncfn;
    flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
    if (check_flag(flag, "ARKodeGetNumNonlinSolvIters")) { return -1; }
    flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
    if (check_flag(flag, "ARKodeGetNumNonlinSolvConvFails")) { return -1; }

    long int nsetups, nje;
    flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
    if (check_flag(flag, "ARKodeGetNumLinSolvSetups")) { return -1; }
    flag = ARKodeGetNumJacEvals(arkode_mem, &nje);
    if (check_flag(flag, "ARKodeGetNumJacEvals")) { return -1; }

    cout << "  NLS iters               = " << nni << endl;
    cout << "  NLS fails               = " << ncfn << endl;
    cout << "  LS setups               = " << nsetups << endl;
    cout << "  J evals                 = " << nje << endl;
    cout << endl;

    // Compute average nls iters per step and ls iters per nls iter
    sunrealtype avgnli = (sunrealtype)nni / (sunrealtype)nst;
    sunrealtype avgls  = (sunrealtype)nsetups / (sunrealtype)nni;
    cout << "  Avg NLS iters per slow step = " << avgnli << endl;
    cout << "  Avg LS setups per NLS iter  = " << avgls << endl;
  }
  cout << endl;

  // Get fast integrator stats and solver stats
  void* fast_arkode_mem;
  MRIStepInnerStepper_GetContent(fast_mem, &fast_arkode_mem);

  // Get fast integrator and solver stats
  flag = ARKodeGetNumSteps(fast_arkode_mem, &nst);
  if (check_flag(flag, "ARKodeGetNumSteps")) { return -1; }
  flag = ARKodeGetNumStepAttempts(fast_arkode_mem, &nst_a);
  if (check_flag(flag, "ARKodeGetNumStepAttempts")) { return -1; }
  flag = ARKodeGetNumErrTestFails(fast_arkode_mem, &netf);
  if (check_flag(flag, "ARKodeGetNumErrTestFails")) { return -1; }
  flag = ARKodeGetNumRhsEvals(fast_arkode_mem, 0, &nfe);
  if (check_flag(flag, "ARKodeGetNumRhsEvals")) { return -1; }
  flag = ARKodeGetNumRhsEvals(fast_arkode_mem, 1, &nfi);
  if (check_flag(flag, "ARKodeGetNumRhsEvals")) { return -1; }

  cout << fixed << setprecision(6);
  cout << endl << "Fast Integrator:" << endl;
  cout << "  Steps                   = " << nst << endl;
  cout << "  Step attempts           = " << nst_a << endl;
  cout << "  Error test fails        = " << netf << endl;
  cout << "  Slow explicit RHS evals = " << nfe << endl;
  cout << "  Slow implicit RHS evals = " << nfi << endl;

  if (udata.splitting)
  {
    long int nni, ncfn;
    flag = ARKodeGetNumNonlinSolvIters(fast_arkode_mem, &nni);
    if (check_flag(flag, "ARKodeGetNumNonlinSolvIters")) { return -1; }
    flag = ARKodeGetNumNonlinSolvConvFails(fast_arkode_mem, &ncfn);
    if (check_flag(flag, "ARKodeGetNumNonlinSolvConvFails")) { return -1; }

    long int nsetups, nje;
    flag = ARKodeGetNumLinSolvSetups(fast_arkode_mem, &nsetups);
    if (check_flag(flag, "ARKodeGetNumLinSolvSetups")) { return -1; }
    flag = ARKodeGetNumJacEvals(fast_arkode_mem, &nje);
    if (check_flag(flag, "ARKodeGetNumJacEvals")) { return -1; }

    cout << "  NLS iters               = " << nni << endl;
    cout << "  NLS fails               = " << ncfn << endl;
    cout << "  LS setups               = " << nsetups << endl;
    cout << "  J evals                 = " << nje << endl;
    cout << endl;

    sunrealtype avgnli = (sunrealtype)nni / (sunrealtype)nst;
    sunrealtype avgls  = (sunrealtype)nsetups / (sunrealtype)nni;
    cout << "  Avg NLS iters per fast step = " << avgnli << endl;
    cout << "  Avg LS setups per NLS iter  = " << avgls << endl;
  }
  cout << endl;

  return 0;
}

// Save current stats
static int UpdateCVodeStats(CVodeInnerStepperContent* content)
{
  int flag;
  long int nst, netf, nfe, nni, nncf, nsetups, nje;

  flag = CVodeGetNumSteps(content->cvode_mem, &nst);
  if (check_flag(flag, "CVodeGetNumSteps")) { return -1; }
  content->nst += nst;

  flag = CVodeGetNumErrTestFails(content->cvode_mem, &netf);
  if (check_flag(flag, "CVodeGetNumErrTestFails")) { return -1; }
  content->netf += netf;

  flag = CVodeGetNumRhsEvals(content->cvode_mem, &nfe);
  if (check_flag(flag, "CVodeGetNumRhsEvals")) { return -1; }
  content->nfe += nfe;

  flag = CVodeGetNumNonlinSolvIters(content->cvode_mem, &nni);
  if (check_flag(flag, "CVodeGetNumNonlinSolvIters")) { return -1; }
  content->nni += nni;

  flag = CVodeGetNumNonlinSolvConvFails(content->cvode_mem, &nncf);
  if (check_flag(flag, "CVodeGetNumNonlinSolvConvFails")) { return -1; }
  content->nncf += nncf;

  flag = CVodeGetNumLinSolvSetups(content->cvode_mem, &nsetups);
  if (check_flag(flag, "CVodeGetNumLinSolveSetups")) { return -1; }
  content->nsetups += nsetups;

  flag = CVodeGetNumJacEvals(content->cvode_mem, &nje);
  if (check_flag(flag, "CVodeGetNumJacEvals")) { return -1; }
  content->nje += nje;

  return 0;
}

// Print MRI integrator statistics
static int OutputStatsMRICVODE(void* arkode_mem, MRIStepInnerStepper fast_mem,
                               UserData& udata)
{
  int flag;

  // Get slow integrator and solver stats
  long int nsts, nfse, nfsi;
  flag = ARKodeGetNumSteps(arkode_mem, &nsts);
  if (check_flag(flag, "ARKodeGetNumSteps")) { return -1; }
  flag = ARKodeGetNumRhsEvals(arkode_mem, 0, &nfse);
  if (check_flag(flag, "ARKodeGetNumRhsEvals")) { return -1; }
  flag = ARKodeGetNumRhsEvals(arkode_mem, 1, &nfsi);
  if (check_flag(flag, "ARKodeGetNumRhsEvals")) { return -1; }

  long int nni, ncfn;
  flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
  if (check_flag(flag, "ARKodeGetNumNonlinSolvIters")) { return -1; }
  flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  if (check_flag(flag, "ARKodeGetNumNonlinSolvConvFails")) { return -1; }

  long int nsetups, nje;
  flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
  if (check_flag(flag, "ARKodeGetNumLinSolvSetups")) { return -1; }
  flag = ARKodeGetNumJacEvals(arkode_mem, &nje);
  if (check_flag(flag, "ARKodeGetNumJacEvals")) { return -1; }

  cout << fixed << setprecision(6);
  cout << endl << "Slow Integrator:" << endl;
  cout << "  Steps                   = " << nsts << endl;
  cout << "  Slow explicit RHS evals = " << nfse << endl;
  cout << "  Slow implicit RHS evals = " << nfsi << endl;
  cout << "  NLS iters               = " << nni << endl;
  cout << "  NLS fails               = " << ncfn << endl;
  cout << "  LS setups               = " << nsetups << endl;
  cout << "  J evals                 = " << nje << endl;
  cout << endl;

  // Compute average nls iters per step and ls iters per nls iter
  sunrealtype avgnli = (sunrealtype)nni / (sunrealtype)nsts;
  sunrealtype avgls  = (sunrealtype)nsetups / (sunrealtype)nni;
  cout << "  Avg NLS iters per slow step = " << avgnli << endl;
  cout << "  Avg LS setups per NLS iter  = " << avgls << endl;
  cout << endl;

  // Get fast integrator stats and solver stats
  void* inner_content;
  MRIStepInnerStepper_GetContent(fast_mem, &inner_content);
  CVodeInnerStepperContent* content = (CVodeInnerStepperContent*)inner_content;

  // Update CVODE stats
  flag = UpdateCVodeStats(content);
  if (check_flag(flag, "UpdateCVodeStats")) { return -1; }

  cout << fixed << setprecision(6);
  cout << "Fast Integrator:" << endl;
  cout << "  Steps            = " << content->nst << endl;
  cout << "  Error test fails = " << content->netf << endl;
  cout << "  Fast RHS evals   = " << content->nfe << endl;
  cout << "  NLS iters        = " << content->nni << endl;
  cout << "  NLS fails        = " << content->nncf << endl;
  cout << "  LS setups        = " << content->nsetups << endl;
  cout << "  J evals          = " << content->nje << endl;
  cout << endl;

  avgnli = (sunrealtype)content->nni / (sunrealtype)content->nst;
  avgls  = (sunrealtype)content->nsetups / (sunrealtype)content->nni;
  cout << "  Avg NLS iters per fast step = " << avgnli << endl;
  cout << "  Avg LS setups per NLS iter  = " << avgls << endl;
  cout << endl;

  return 0;
}

// Print command line options
static void InputHelp()
{
  cout << endl;
  cout << "Command line options:" << endl;
  cout << "  --no-advection           : disable advection\n";
  cout << "  --no-diffusion           : disable diffusion\n";
  cout << "  --splitting <int>        : RHS splitting\n";
  cout << "  --c <real>               : advection coefficient\n";
  cout << "  --d <real>               : diffusion coefficient\n";
  cout << "  --A <real>               : species A concentration\n";
  cout << "  --B <real>               : species B concentration\n";
  cout << "  --eps <real>             : stiffness parameter\n";
  cout << "  --tf <real>              : final time\n";
  cout << "  --xl <real>              : domain lower boundary\n";
  cout << "  --xu <real>              : domain upper boundary\n";
  cout << "  --nx <int>               : number of mesh points\n";
  cout << "  --integrator <int>       : integrator option\n";
  cout << "  --order <int>            : method order\n";
  cout << "  --order_fast <int>       : MRI fast method order\n";
  cout << "  --ark_dirk               : Use DIRK method from ARK method\n";
  cout << "  --rtol <real>            : relative tolerance\n";
  cout << "  --atol <real>            : absolute tolerance\n";
  cout << "  --rtol_fast <real>       : MRI fast relative tolerance\n";
  cout << "  --atol_fast <real>       : MRI fast absolute tolerance\n";
  cout << "  --fixed_h <real>         : fixed step size\n";
  cout << "  --fixed_h_fast <real>    : MRI fast fixed step size\n";
  cout << "  --controller <name>      : time step adaptivity\n";
  cout << "  --predictor <int>        : nonlinear solver predictor\n";
  cout << "  --lssetupfreq <int>      : LS setup frequency\n";
  cout << "  --controller_fast <name> : MRI fast time step adaptivity\n";
  cout << "  --predictor_fast <int>   : MRI fast nonlinear solver predictor\n";
  cout << "  --lssetupfreq_fast <int> : MRI fast LS setup frequency\n";
  cout << "  --maxsteps <int>         : max steps between outputs\n";
  cout << "  --etamx1_fast <real>     : max step size growth in first step\n";
  cout << "  --linear                 : linearly implicit\n";
  cout << "  --save_hinit             : reuse initial fast step\n";
  cout << "  --save_hcur              : reuse current fast step\n";
  cout << "  --hcur_factor            : current fast step safety factor\n";
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
  find_arg(args, "--no-diffusion", udata.diffusion, false);
  find_arg(args, "--splitting", udata.splitting);
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
  find_arg(args, "--integrator", uopts.integrator);
  find_arg(args, "--order", uopts.order);
  find_arg(args, "--order_fast", uopts.order_fast);
  find_arg(args, "--ark_dirk", uopts.ark_dirk);
  find_arg(args, "--rtol", uopts.rtol);
  find_arg(args, "--atol", uopts.atol);
  find_arg(args, "--rtol_fast", uopts.rtol_fast);
  find_arg(args, "--atol_fast", uopts.atol_fast);
  find_arg(args, "--fixed_h", uopts.fixed_h);
  find_arg(args, "--fixed_h_fast", uopts.fixed_h_fast);
  find_arg(args, "--predictor", uopts.predictor);
  find_arg(args, "--controller", uopts.controller);
  find_arg(args, "--lssetupfreq", uopts.ls_setup_freq);
  find_arg(args, "--predictor_fast", uopts.predictor_fast);
  find_arg(args, "--controller_fast", uopts.controller_fast);
  find_arg(args, "--lssetupfreq_fast", uopts.ls_setup_freq_fast);
  find_arg(args, "--maxsteps", uopts.maxsteps);
  find_arg(args, "--etamx1_fast", uopts.etamx1_fast);
  find_arg(args, "--linear", uopts.linear);
  find_arg(args, "--save_hinit", uopts.save_hinit);
  find_arg(args, "--save_hcur", uopts.save_hcur);
  find_arg(args, "--hcur_factor", uopts.hcur_factor);
  find_arg(args, "--output", uopts.output);
  find_arg(args, "--nout", uopts.nout);

  // Recompute mesh spacing and total number of nodes
  udata.dx  = (udata.xu - udata.xl) / (udata.nx - 1);
  udata.neq = NSPECIES * udata.nx;

  // Create workspace
  switch (uopts.integrator)
  {
  case (0):
    // Create workspace vector
    udata.temp_v = N_VNew_Serial(udata.neq, ctx);
    if (check_ptr(udata.temp_v, "N_VNew_Serial")) { return 1; }
    N_VConst(ZERO, udata.temp_v);
    break;
  case (1):
    // Create workspace vector
    udata.temp_v = N_VNew_Serial(udata.neq, ctx);
    if (check_ptr(udata.temp_v, "N_VNew_Serial")) { return 1; }
    N_VConst(ZERO, udata.temp_v);
    // Create workspace matrix
    udata.temp_J = SUNBandMatrix(udata.neq, 3, 3, ctx);
    if (check_ptr(udata.temp_J, "SUNBandMatrix")) { return 1; }
    SUNMatZero(udata.temp_J);
    break;
  case (2): break;
  case (3): break;
  default: cerr << "Invalid integrator option" << endl; return 1;
  }

  // Input checks
  if (!udata.diffusion && !udata.advection)
  {
    cerr << "ERROR: Invalid problem configuration" << endl;
    return -1;
  }

  if (udata.diffusion && udata.advection)
  {
    if (udata.splitting < 0 || udata.splitting > 7)
    {
      cerr << "ERROR: Invalid splitting option" << endl;
      return -1;
    }
  }
  else
  {
    if (udata.splitting < 0 || udata.splitting > 4)
    {
      cerr << "ERROR: Invalid splitting option" << endl;
      return -1;
    }
  }

  if (uopts.integrator < 0 || uopts.integrator > 3)
  {
    cerr << "ERROR: Invalid integrator option" << endl;
    return -1;
  }

  // ERKStep setup requires splitting 0 (fully explicit)
  if (uopts.integrator == 0) { udata.splitting = 0; }

  // MRIStep + ARKStep requires splitting 0 or 1 (explicit v implicit reactions)
  if (uopts.integrator == 2 && udata.splitting) { udata.splitting = 1; }

  // MRIStep + CVODE requires splitting 1 (implicit reactions)
  if (uopts.integrator == 3) { udata.splitting = 1; }

  // MRIStep + CVODE requires adaptive fast step sizes
  if (uopts.integrator == 3) { uopts.fixed_h_fast = ZERO; }

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
  cout << " --------------------------------- " << endl;
  cout << "  tf               = " << udata.tf << endl;
  cout << "  xl               = " << udata.xl << endl;
  cout << "  xu               = " << udata.xu << endl;
  cout << "  nx               = " << udata.nx << endl;
  cout << "  dx               = " << udata.dx << endl;
  cout << " --------------------------------- " << endl;

  if (uopts.integrator == 0)
  {
    cout << "  integrator       = ERK" << endl;
    if (udata.advection) { cout << "  advection        = Explicit" << endl; }
    else { cout << "  advection        = OFF" << endl; }
    if (udata.diffusion) { cout << "  diffusion        = Explicit" << endl; }
    else { cout << "  diffusion        = OFF" << endl; }
    cout << "  reaction         = Explicit" << endl;
  }
  else if (uopts.integrator == 1)
  {
    cout << "  integrator       = ARK" << endl;
    // advection-diffusion-reaction
    if (udata.diffusion && udata.advection)
    {
      switch (udata.splitting)
      {
      case (0):
        // ERK -- fully explicit
        cout << "  advection        = Explicit" << endl;
        cout << "  diffusion        = Explicit" << endl;
        cout << "  reaction         = Explicit" << endl;
        break;
      case (1):
        // IMEX -- explicit advection-diffusion, implicit reaction
        cout << "  advection        = Explicit" << endl;
        cout << "  diffusion        = Explicit" << endl;
        cout << "  reaction         = Implicit" << endl;
        break;
      case (2):
        // IMEX -- explicit advection-reaction, implicit diffusion
        cout << "  advection        = Explicit" << endl;
        cout << "  diffusion        = Implicit" << endl;
        cout << "  reaction         = Explicit" << endl;
        break;
      case (3):
        // IMEX -- explicit advection, implicit diffusion-reaction
        cout << "  advection        = Explicit" << endl;
        cout << "  diffusion        = Implicit" << endl;
        cout << "  reaction         = Implicit" << endl;
        break;
      case (4):
        // IMEX -- explicit diffusion-reaction, implicit advection
        cout << "  advection        = Implicit" << endl;
        cout << "  diffusion        = Explicit" << endl;
        cout << "  reaction         = Explicit" << endl;
        break;
      case (5):
        // IMEX -- explicit diffusion, implicit advection-reaction
        cout << "  advection        = Implicit" << endl;
        cout << "  diffusion        = Explicit" << endl;
        cout << "  reaction         = Implicit" << endl;
        break;
      case (6):
        // IMEX -- explicit reaction, implicit advection-diffusion
        cout << "  advection        = Implicit" << endl;
        cout << "  diffusion        = Implicit" << endl;
        cout << "  reaction         = Explicit" << endl;
        break;
      case (7):
        // DIRK -- fully implicit
        cout << "  advection        = Implicit" << endl;
        cout << "  diffusion        = Implicit" << endl;
        cout << "  reaction         = Implicit" << endl;
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
        cout << "  advection        = Explicit" << endl;
        cout << "  diffusion        = OFF" << endl;
        cout << "  reaction         = Explicit" << endl;
        break;
      case (1):
        // IMEX -- explicit advection, implicit reaction
        cout << "  advection        = Explicit" << endl;
        cout << "  diffusion        = OFF" << endl;
        cout << "  reaction         = Implicit" << endl;
        break;
      case (2):
        // IMEX -- explicit reaction, implicit advection
        cout << "  advection        = Implicit" << endl;
        cout << "  diffusion        = OFF" << endl;
        cout << "  reaction         = Explicit" << endl;
        break;
      case (3):
        // DIRK -- fully implicit
        cout << "  advection        = Implicit" << endl;
        cout << "  diffusion        = OFF" << endl;
        cout << "  reaction         = Implicit" << endl;
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
        cout << "  advection        = OFF" << endl;
        cout << "  diffusion        = Explicit" << endl;
        cout << "  reaction         = Explicit" << endl;
        break;
      case (1):
        // IMEX -- explicit diffusion, implicit reaction
        cout << "  advection        = OFF" << endl;
        cout << "  diffusion        = Explicit" << endl;
        cout << "  reaction         = Implicit" << endl;
        break;
      case (2):
        // IMEX -- explicit reaction, implicit diffusion
        cout << "  advection        = OFF" << endl;
        cout << "  diffusion        = Implicit" << endl;
        cout << "  reaction         = Explicit" << endl;
        break;
      case (4):
        // DIRK -- fully implicit
        cout << "  advection        = OFF" << endl;
        cout << "  diffusion        = Implicit" << endl;
        cout << "  reaction         = Implicit" << endl;
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
  }
  else if (uopts.integrator == 2)
  {
    cout << "  integrator       = MRI + ARK" << endl;
    if (udata.diffusion && udata.advection)
    {
      // IMEX slow -- advection-diffusion
      cout << "  advection        = Slow-Explicit" << endl;
      cout << "  diffusion        = Slow-Implicit" << endl;
      if (udata.splitting)
      {
        cout << "  reaction         = Fast-Implicit" << endl;
      }
      else { cout << "  reaction         = Fast-Explicit" << endl; }
    }
    else if (!udata.diffusion && udata.advection)
    {
      // Explicit slow -- advection
      cout << "  advection        = Slow-Explicit" << endl;
      cout << "  diffusion        = OFF" << endl;
      if (udata.splitting)
      {
        cout << "  reaction         = Fast-Implicit" << endl;
      }
      else { cout << "  reaction         = Fast-Explicit" << endl; }
    }
    else if (udata.diffusion && !udata.advection)
    {
      cout << "  advection        = OFF" << endl;
      cout << "  diffusion        = Implicit" << endl;
      if (udata.splitting)
      {
        cout << "  reaction         = Fast-Implicit" << endl;
      }
      else { cout << "  reaction         = Fast-Explicit" << endl; }
    }
    else
    {
      // No slow time scale
      cerr << "ERROR: Invalid problem configuration" << endl;
      return -1;
    }
  }
  else if (uopts.integrator == 3)
  {
    cout << "  integrator       = MRI + CVODE" << endl;
    // Slow time scale
    if (udata.diffusion && udata.advection)
    {
      // IMEX slow -- advection-diffusion
      cout << "  advection        = Slow-Explicit" << endl;
      cout << "  diffusion        = Slow-Implicit" << endl;
      cout << "  reaction         = Fast-Implicit" << endl;
    }
    else if (!udata.diffusion && udata.advection)
    {
      // Explicit slow -- advection
      cout << "  advection        = Slow-Explicit" << endl;
      cout << "  diffusion        = OFF" << endl;
      cout << "  reaction         = Fast-Implicit" << endl;
    }
    else if (udata.diffusion && !udata.advection)
    {
      // Implicit slow -- diffusion
      cout << "  advection        = OFF" << endl;
      cout << "  diffusion        = Slow-Implicit" << endl;
      cout << "  reaction         = Fast-Implicit" << endl;
    }
    else
    {
      // No slow time scale
      cerr << "ERROR: Invalid problem configuration" << endl;
      return -1;
    }
  }
  else
  {
    cerr << "ERROR: Invalid integrator option" << endl;
    return -1;
  }
  if (uopts.ark_dirk)
  {
    cout << "  order (ark_dirk) = " << uopts.order << endl;
  }
  else { cout << "  order            = " << uopts.order << endl; }
  cout << "  rtol             = " << uopts.rtol << endl;
  cout << "  atol             = " << uopts.atol << endl;
  cout << "  fixed h          = " << uopts.fixed_h << endl;
  cout << "  controller       = " << uopts.controller << endl;
  if (uopts.integrator > 0)
  {
    if (uopts.predictor == 0)
    {
      cout << "  predictor        = trivial" << endl;
    }
    else if (uopts.predictor == 1)
    {
      cout << "  predictor        = max order" << endl;
    }
    else if (uopts.predictor == 2)
    {
      cout << "  predictor        = variable order" << endl;
    }
    else if (uopts.predictor == 3)
    {
      cout << "  predictor        = cutoff order" << endl;
    }
    else { cout << "  predictor        = " << uopts.predictor << endl; }
    cout << "  ls setup freq    = " << uopts.ls_setup_freq << endl;
    cout << "  linear           = " << uopts.linear << endl;
  }
  if (uopts.integrator == 2)
  {
    cout << " --------------------------------- " << endl;
    cout << "  fast integrator  = ARK" << endl;
    cout << "  rtol             = " << uopts.rtol_fast << endl;
    cout << "  atol             = " << uopts.atol_fast << endl;
    cout << "  order            = " << uopts.order_fast << endl;
    cout << "  fixed h          = " << uopts.fixed_h_fast << endl;
    cout << "  fast controller  = PID" << uopts.controller_fast << endl;
    if (uopts.predictor_fast == 0)
    {
      cout << "  predictor        = trivial" << endl;
    }
    else if (uopts.predictor_fast == 1)
    {
      cout << "  predictor        = max order" << endl;
    }
    else if (uopts.predictor_fast == 2)
    {
      cout << "  predictor        = variable order" << endl;
    }
    else if (uopts.predictor_fast == 3)
    {
      cout << "  predictor        = cutoff order" << endl;
    }
    else { cout << "  predictor        = " << uopts.predictor_fast << endl; }
    cout << "  ls setup freq    = " << uopts.ls_setup_freq_fast << endl;
  }
  else if (uopts.integrator == 3)
  {
    cout << " --------------------------------- " << endl;
    cout << "  fast integrator  = CVODE" << endl;
    cout << "  rtol             = " << uopts.rtol_fast << endl;
    cout << "  atol             = " << uopts.atol_fast << endl;
    cout << "  ls setup freq    = " << uopts.ls_setup_freq << endl;
    cout << "  etamx first step = " << uopts.etamx1_fast << endl;
    cout << "  reuse initial h  = " << uopts.save_hinit << endl;
    cout << "  reuse current h  = " << uopts.save_hcur << endl;
    cout << "  current h factor = " << uopts.hcur_factor << endl;
  }
  cout << " --------------------------------- " << endl;
  cout << "  output          = " << uopts.output << endl;
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
    cout << "          ||y||_rms      " << endl;
    cout << " ---------------------";
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
    cout << "-------------------------" << endl;
    cout << endl;
  }

  // Close output streams
  if (uopts.output >= 2) { uopts.uout.close(); }

  return 0;
}

//---- end of file ----
