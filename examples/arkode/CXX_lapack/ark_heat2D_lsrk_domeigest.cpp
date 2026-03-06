/* -----------------------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ SMU
 *
 * (adapted from ark_heat2D.cpp and adapted from ark_heat2D_lsrk.cpp,
 * co-authored by Daniel Reynolds and David Gardner (LLNL))
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
 * for t in [0, 1] and (x,y) in [0, 1]^2, with initial condition
 *
 *   u(0,x,y) = sin^2(pi x) sin^2(pi y),
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
 *    u(t,x,y) = sin^2(pi x) sin^2(pi y) cos^2(pi t).
 *
 * The spatial derivatives are computed using second-order centered differences,
 * with the data distributed over nx * ny points on a uniform spatial grid. The
 * problem is advanced in time with the LSRKStep module in ARKODE.
 * Several command line options are available to change the problem parameters
 * and ARKODE settings. Use the flag --help for more information.
 * ---------------------------------------------------------------------------*/

#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "arkode/arkode_lsrkstep.h" // access to LSRKStep
#include "nvector/nvector_serial.h" // access to the serial N_Vector
#include "sunadaptcontroller/sunadaptcontroller_imexgus.h"
#include "sunadaptcontroller/sunadaptcontroller_soderlind.h"

#include <sundomeigest/sundomeigest_arnoldi.h> // access to Arnoldi Iteration module
#include <sundomeigest/sundomeigest_power.h> // access to Power Iteration module

// Macros for problem constants
#define PI   SUN_RCONST(3.141592653589793238462643383279502884197169)
#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)

// Macro to access (x,y) location in 1D NVector array
#define IDX(x, y, n) ((n) * (y) + (x))

using namespace std;

// -----------------------------------------------------------------------------
// User data structure
// -----------------------------------------------------------------------------

struct UserData
{
  // Diffusion coefficients in the x and y directions
  sunrealtype kx;
  sunrealtype ky;

  // Enable/disable forcing
  bool forcing;

  // Final time
  sunrealtype tf;

  // Upper bounds in x and y directions
  sunrealtype xu;
  sunrealtype yu;

  // Number of nodes in the x and y directions
  sunindextype nx;
  sunindextype ny;

  // Total number of nodes
  sunindextype nodes;

  // Mesh spacing in the x and y directions
  sunrealtype dx;
  sunrealtype dy;

  // Integrator settings
  sunrealtype rtol;       // relative tolerance
  sunrealtype atol;       // absolute tolerance
  sunrealtype hfixed;     // fixed step size
  std::string controller; // step size adaptivity method
  int maxsteps;           // max number of steps between outputs
  bool diagnostics;       // output diagnostics

  // LSRKStep options
  ARKODE_LSRKMethodType method; // LSRK method choice
  long int eigfrequency;        // dominant eigenvalue update frequency
  int stage_max_limit;          // maximum number of stages per step
  sunrealtype eigsafety;        // dominant eigenvalue safety factor

  // Output variables
  int output;    // output level
  int nout;      // number of output times
  ofstream uout; // output file stream
  ofstream eout; // error file stream
  N_Vector e;    // error vector

  // DEE options
  int dee_id;            // DEE ID
  int dee_num_init_wups; // number of initial warmups before the first estimate
  int dee_num_succ_wups; // number of succeeding warmups before each estimate
  int dee_max_iters;     // max number of iterations
  int dee_krylov_dim;    // Krylov dimension for DEE
  double dee_reltol;     // tolerance

  // Timing variables
  bool timing; // print timings
  double evolvetime;
  double rhstime;
};

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrator
// -----------------------------------------------------------------------------

// ODE right hand side function
static int f(sunrealtype t, N_Vector u, N_Vector f, void* user_data);

// -----------------------------------------------------------------------------
// UserData and input functions
// -----------------------------------------------------------------------------

// Set the default values in the UserData structure
static int InitUserData(UserData* udata);

// Free memory allocated within UserData
static int FreeUserData(UserData* udata);

// Read the command line inputs and set UserData values
static int ReadInputs(int* argc, char*** argv, UserData* udata);

// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------

// Compute the true solution
static int Solution(sunrealtype t, N_Vector u, UserData* udata);

// Compute the solution error solution
static int SolutionError(sunrealtype t, N_Vector u, N_Vector e, UserData* udata);

// Print the command line options
static void InputHelp();

// Print some UserData information
static int PrintUserData(UserData* udata);

// Output solution and error
static int OpenOutput(UserData* udata);
static int WriteOutput(sunrealtype t, N_Vector u, UserData* udata);
static int CloseOutput(UserData* udata);

// Print integration timing
static int OutputTiming(UserData* udata);

// Check function return values
static int check_flag(void* flagvalue, const string funcname, int opt);

// -----------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int flag;                // reusable error-checking flag
  UserData* udata  = NULL; // user data structure
  N_Vector u       = NULL; // vector for storing solution
  void* arkode_mem = NULL; // ARKODE memory structure

  // Dominant Eigenvalue Estimator (DEE) pointers and variables
  SUNDomEigEstimator DEE = NULL; /* domeig estimator object */
  N_Vector q             = NULL; /* random initial eigenvector */

  // Timing variables
  chrono::time_point<chrono::steady_clock> t1;
  chrono::time_point<chrono::steady_clock> t2;

  // Create the SUNDIALS context object for this simulation
  SUNContext ctx;
  flag = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_flag(&flag, "SUNContext_Create", 1)) { return 1; }

  // ---------------
  // Setup UserData
  // ---------------

  // Allocate and initialize user data structure with default values. The
  // defaults may be overwritten by command line inputs in ReadInputs below.
  udata = new UserData;
  flag  = InitUserData(udata);
  if (check_flag(&flag, "InitUserData", 1)) { return 1; }

  // Parse command line inputs
  flag = ReadInputs(&argc, &argv, udata);
  if (flag != 0) { return 1; }

  // Output problem setup/options
  flag = PrintUserData(udata);
  if (check_flag(&flag, "PrintUserData", 1)) { return 1; }

  if (udata->diagnostics)
  {
    SUNLogger logger = NULL;

    flag = SUNContext_GetLogger(ctx, &logger);
    if (check_flag(&flag, "SUNContext_GetLogger", 1)) { return 1; }

    flag = SUNLogger_SetInfoFilename(logger, "diagnostics.txt");
    if (check_flag(&flag, "SUNLogger_SetInfoFilename", 1)) { return 1; }

    flag = SUNLogger_SetDebugFilename(logger, "diagnostics.txt");
    if (check_flag(&flag, "SUNLogger_SetDebugFilename", 1)) { return 1; }
  }

  // ----------------------
  // Create serial vectors
  // ----------------------

  // Create vector for solution
  u = N_VNew_Serial(udata->nodes, ctx);
  if (check_flag((void*)u, "N_VNew_Serial", 0)) { return 1; }

  // Set initial condition
  flag = Solution(ZERO, u, udata);
  if (check_flag(&flag, "Solution", 1)) { return 1; }

  // Create vector for error
  udata->e = N_VClone(u);
  if (check_flag((void*)(udata->e), "N_VClone", 0)) { return 1; }

  // --------------
  // Setup ARKODE
  // --------------

  // Create integrator
  arkode_mem = LSRKStepCreateSTS(f, ZERO, u, ctx);
  if (check_flag((void*)arkode_mem, "LSRKStepCreateSTS", 0)) { return 1; }

  // Specify tolerances
  flag = ARKodeSStolerances(arkode_mem, udata->rtol, udata->atol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) { return 1; }

  // Attach user data
  flag = ARKodeSetUserData(arkode_mem, (void*)udata);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }

  // Select LSRK method
  flag = LSRKStepSetSTSMethod(arkode_mem, udata->method);
  if (check_flag(&flag, "LSRKStepSetSTSMethod", 1)) { return 1; }

  /* Set the initial random eigenvector for the DEE */
  q = N_VClone(u);
  if (check_flag(q, "N_VClone", 0)) { return 1; }

  sunrealtype* qd = N_VGetArrayPointer(q);
  for (int i = 0; i < udata->nodes; i++)
  {
    qd[i] = (sunrealtype)rand() / (sunrealtype)RAND_MAX;
  }

  if (udata->dee_id == 0)
  {
    /* Create power iteration dominant eigenvalue estimator */
    DEE = SUNDomEigEstimator_Power(q, udata->dee_max_iters, udata->dee_reltol,
                                   ctx);
    if (check_flag(DEE, "SUNDomEigEstimator_Power", 0)) { return 1; }
  }
  else if (udata->dee_id == 1)
  {
    /* Create Arnoldi dominant eigenvalue estimator */
    DEE = SUNDomEigEstimator_Arnoldi(q, udata->dee_krylov_dim, ctx);
    if (check_flag(DEE, "SUNDomEigEstimator_Arnoldi", 0)) { return 1; }
  }
  else
  {
    cerr << "Invalid DEE ID: " << udata->dee_id << endl;
    return 1;
  }

  /* After the DEE creation, random q vector is no longer needed.
     It is used only to initialize the DEE */
  N_VDestroy(q);

  /* Attach the DEE to the LSRKStep module.
  There is no need to set Atimes or initialize since these are all
  performed after attaching the DEE by LSRKStep. */
  flag = LSRKStepSetDomEigEstimator(arkode_mem, DEE);
  if (check_flag(&flag, "LSRKStepSetDomEigEstimator", 2)) { return 1; }

  /* Specify the number of preprocessing warmups before the first estimate call. */
  flag = LSRKStepSetNumDomEigEstInitPreprocessIters(arkode_mem,
                                                    udata->dee_num_init_wups);
  if (check_flag(&flag, "LSRKStepSetNumDomEigEstInitPreprocessIters", 2))
  {
    return 1;
  }

  /* Specify the succeeding warmups before each estimate call.
     This is the number of warmups that will be performed after the first
     estimate call and before each subsequent estimate calls. */
  flag = LSRKStepSetNumDomEigEstPreprocessIters(arkode_mem,
                                                udata->dee_num_succ_wups);
  if (check_flag(&flag, "LSRKStepSetNumDomEigEstPreprocessIters", 1))
  {
    return 1;
  }

  flag = LSRKStepSetDomEigFrequency(arkode_mem, udata->eigfrequency);
  if (check_flag(&flag, "LSRKStepSetDomEigFrequency", 1)) { return 1; }

  // Set maximum number of stages per step
  flag = LSRKStepSetMaxNumStages(arkode_mem, udata->stage_max_limit);
  if (check_flag(&flag, "LSRKStepSetMaxNumStages", 1)) { return 1; }

  // Set spectral radius safety factor
  flag = LSRKStepSetDomEigSafetyFactor(arkode_mem, udata->eigsafety);
  if (check_flag(&flag, "LSRKStepSetDomEigSafetyFactor", 1)) { return 1; }

  // Set fixed step size or adaptivity method
  if (udata->hfixed > ZERO)
  {
    flag = ARKodeSetFixedStep(arkode_mem, udata->hfixed);
    if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }
  }
  else
  {
    flag = ARKodeSetAdaptControllerByName(arkode_mem, udata->controller.c_str());
    if (check_flag(&flag, "ARKodeSetAdaptControllerByName", 1)) { return 1; }
  }

  // Set max steps between outputs
  flag = ARKodeSetMaxNumSteps(arkode_mem, udata->maxsteps);
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) { return 1; }

  // Set stopping time
  flag = ARKodeSetStopTime(arkode_mem, udata->tf);
  if (check_flag(&flag, "ARKodeSetStopTime", 1)) { return 1; }

  // -----------------------
  // Loop over output times
  // -----------------------

  sunrealtype t     = ZERO;
  sunrealtype dTout = udata->tf / udata->nout;
  sunrealtype tout  = dTout;

  // initial output
  flag = OpenOutput(udata);
  if (check_flag(&flag, "OpenOutput", 1)) { return 1; }

  flag = WriteOutput(t, u, udata);
  if (check_flag(&flag, "WriteOutput", 1)) { return 1; }

  for (int iout = 0; iout < udata->nout; iout++)
  {
    // Start timer
    t1 = chrono::steady_clock::now();

    // Evolve in time
    flag = ARKodeEvolve(arkode_mem, tout, u, &t, ARK_NORMAL);
    if (check_flag(&flag, "ARKodeEvolve", 1)) { break; }

    // Stop timer
    t2 = chrono::steady_clock::now();

    // Update timer
    udata->evolvetime += chrono::duration<double>(t2 - t1).count();

    // Output solution and error
    flag = WriteOutput(t, u, udata);
    if (check_flag(&flag, "WriteOutput", 1)) { return 1; }

    // Update output time
    tout += dTout;
    tout = (tout > udata->tf) ? udata->tf : tout;
  }

  // Close output
  flag = CloseOutput(udata);
  if (check_flag(&flag, "CloseOutput", 1)) { return 1; }

  // --------------
  // Final outputs
  // --------------

  // Print final integrator stats
  if (udata->output > 0)
  {
    cout << "Final integrator statistics:" << endl;
    flag = ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
    if (check_flag(&flag, "ARKodePrintAllStats", 1)) { return 1; }
  }

  if (udata->forcing)
  {
    // Output final error
    flag = SolutionError(t, u, udata->e, udata);
    if (check_flag(&flag, "SolutionError", 1)) { return 1; }

    sunrealtype maxerr = N_VMaxNorm(udata->e);

    cout << scientific;
    cout << setprecision(numeric_limits<sunrealtype>::digits10);
    cout << "  Max error = " << maxerr << endl;
  }

  // Print timing
  if (udata->timing)
  {
    flag = OutputTiming(udata);
    if (check_flag(&flag, "OutputTiming", 1)) { return 1; }
  }

  // --------------------
  // Clean up and return
  // --------------------

  ARKodeFree(&arkode_mem);          // Free integrator memory
  N_VDestroy(u);                    // Free vectors
  SUNDomEigEstimator_Destroy(&DEE); /* Free DEE object */
  FreeUserData(udata);              // Free user data
  delete udata;
  SUNContext_Free(&ctx); // Free context

  return 0;
}

// -----------------------------------------------------------------------------
// Functions called by the integrator
// -----------------------------------------------------------------------------

// f routine to compute the ODE RHS function f(t,y).
static int f(sunrealtype t, N_Vector u, N_Vector f, void* user_data)
{
  // Timing variables
  chrono::time_point<chrono::steady_clock> t1;
  chrono::time_point<chrono::steady_clock> t2;

  // Start timer
  t1 = chrono::steady_clock::now();

  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Shortcuts to number of nodes
  sunindextype nx = udata->nx;
  sunindextype ny = udata->ny;

  // Constants for computing diffusion term
  sunrealtype cx = udata->kx / (udata->dx * udata->dx);
  sunrealtype cy = udata->ky / (udata->dy * udata->dy);
  sunrealtype cc = -TWO * (cx + cy);

  // Access data arrays
  sunrealtype* uarray = N_VGetArrayPointer(u);
  if (check_flag((void*)uarray, "N_VGetArrayPointer", 0)) { return -1; }

  sunrealtype* farray = N_VGetArrayPointer(f);
  if (check_flag((void*)farray, "N_VGetArrayPointer", 0)) { return -1; }

  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, f);

  // Iterate over domain interior and compute rhs forcing term
  if (udata->forcing)
  {
    sunrealtype x, y;
    sunrealtype sin_sqr_x, sin_sqr_y;
    sunrealtype cos_sqr_x, cos_sqr_y;

    sunrealtype bx = (udata->kx) * TWO * PI * PI;
    sunrealtype by = (udata->ky) * TWO * PI * PI;

    sunrealtype sin_t_cos_t = sin(PI * t) * cos(PI * t);
    sunrealtype cos_sqr_t   = cos(PI * t) * cos(PI * t);

    for (sunindextype j = 1; j < ny - 1; j++)
    {
      for (sunindextype i = 1; i < nx - 1; i++)
      {
        x = i * udata->dx;
        y = j * udata->dy;

        sin_sqr_x = sin(PI * x) * sin(PI * x);
        sin_sqr_y = sin(PI * y) * sin(PI * y);

        cos_sqr_x = cos(PI * x) * cos(PI * x);
        cos_sqr_y = cos(PI * y) * cos(PI * y);

        farray[IDX(i, j, nx)] =
          -TWO * PI * sin_sqr_x * sin_sqr_y * sin_t_cos_t -
          bx * (cos_sqr_x - sin_sqr_x) * sin_sqr_y * cos_sqr_t -
          by * (cos_sqr_y - sin_sqr_y) * sin_sqr_x * cos_sqr_t;
      }
    }
  }

  // Iterate over domain interior and add rhs diffusion term
  for (sunindextype j = 1; j < ny - 1; j++)
  {
    for (sunindextype i = 1; i < nx - 1; i++)
    {
      farray[IDX(i, j, nx)] +=
        cc * uarray[IDX(i, j, nx)] +
        cx * (uarray[IDX(i - 1, j, nx)] + uarray[IDX(i + 1, j, nx)]) +
        cy * (uarray[IDX(i, j - 1, nx)] + uarray[IDX(i, j + 1, nx)]);
    }
  }

  // Stop timer
  t2 = chrono::steady_clock::now();

  // Update timer
  udata->rhstime += chrono::duration<double>(t2 - t1).count();

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// UserData and input functions
// -----------------------------------------------------------------------------

// Initialize memory allocated within Userdata
static int InitUserData(UserData* udata)
{
  // Diffusion coefficient
  udata->kx = SUN_RCONST(10.0);
  udata->ky = SUN_RCONST(10.0);

  // Enable forcing
  udata->forcing = true;

  // Final time
  udata->tf = ONE;

  // Upper bounds in x and y directions
  udata->xu = ONE;
  udata->yu = ONE;

  // Number of nodes in the x and y directions
  udata->nx    = 64;
  udata->ny    = 64;
  udata->nodes = udata->nx * udata->ny;

  // Mesh spacing in the x and y directions
  udata->dx = udata->xu / (udata->nx - 1);
  udata->dy = udata->yu / (udata->ny - 1);

  // Integrator settings
  udata->rtol        = SUN_RCONST(1.e-5);  // relative tolerance
  udata->atol        = SUN_RCONST(1.e-10); // absolute tolerance
  udata->hfixed      = ZERO;               // using adaptive step sizes
  udata->controller  = "I";                // I controller
  udata->maxsteps    = 0;                  // use default
  udata->diagnostics = false;              // output diagnostics

  // LSRKStep options
  udata->method          = ARKODE_LSRK_RKC_2; // RKC
  udata->eigfrequency    = 25;   // update eigenvalue at least every 25 steps
  udata->stage_max_limit = 1000; // allow up to 1000 stages/step
  udata->eigsafety       = SUN_RCONST(1.01); // 1% safety factor

  // Output variables
  udata->output = 1;  // 0 = no output, 1 = stats output, 2 = output to disk
  udata->nout   = 20; // Number of output times
  udata->e      = NULL;

  // DEE options
  udata->dee_id            = 1; // DEE ID (0 for PI and 1 for Arnoldi)
  udata->dee_num_init_wups = 20;
  udata->dee_num_succ_wups = 5;
  udata->dee_max_iters     = 100;
  udata->dee_krylov_dim    = 3;
  udata->dee_reltol        = 0.005;

  // Timing variables
  udata->timing     = false;
  udata->evolvetime = 0.0;
  udata->rhstime    = 0.0;

  // Return success
  return 0;
}

// Free memory allocated within Userdata
static int FreeUserData(UserData* udata)
{
  // Free error vector
  if (udata->e)
  {
    N_VDestroy(udata->e);
    udata->e = NULL;
  }

  // Return success
  return 0;
}

// Read command line inputs
static int ReadInputs(int* argc, char*** argv, UserData* udata)
{
  // Check for input args
  int arg_idx = 1;

  while (arg_idx < (*argc))
  {
    string arg = (*argv)[arg_idx++];

    // Mesh points
    if (arg == "--mesh")
    {
      udata->nx = stoi((*argv)[arg_idx++]);
      udata->ny = stoi((*argv)[arg_idx++]);
    }
    // Domain upper bounds
    else if (arg == "--domain")
    {
      udata->xu = stoi((*argv)[arg_idx++]);
      udata->yu = stoi((*argv)[arg_idx++]);
    }
    // Diffusion parameters
    else if (arg == "--k")
    {
      udata->kx = stod((*argv)[arg_idx++]);
      udata->ky = stod((*argv)[arg_idx++]);
    }
    // Disable forcing
    else if (arg == "--noforcing") { udata->forcing = false; }
    // Temporal domain settings
    else if (arg == "--tf") { udata->tf = stod((*argv)[arg_idx++]); }
    // Integrator settings
    else if (arg == "--rtol") { udata->rtol = stod((*argv)[arg_idx++]); }
    else if (arg == "--atol") { udata->atol = stod((*argv)[arg_idx++]); }
    else if (arg == "--fixedstep") { udata->hfixed = stod((*argv)[arg_idx++]); }
    else if (arg == "--controller") { udata->controller = (*argv)[arg_idx++]; }
    else if (arg == "--diagnostics") { udata->diagnostics = true; }
    else if (arg == "--method")
    {
      udata->method = (ARKODE_LSRKMethodType)stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--eigfrequency")
    {
      udata->eigfrequency = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--stage_max_limit")
    {
      udata->stage_max_limit = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--eigsafety")
    {
      udata->eigsafety = stod((*argv)[arg_idx++]);
    }
    // Output settings
    else if (arg == "--output") { udata->output = stoi((*argv)[arg_idx++]); }
    else if (arg == "--nout") { udata->nout = stoi((*argv)[arg_idx++]); }
    else if (arg == "--maxsteps")
    {
      udata->maxsteps = stoi((*argv)[arg_idx++]);
    }
    // DEE options
    else if (arg == "--dee_id") { udata->dee_id = stoi((*argv)[arg_idx++]); }
    else if (arg == "--dee_num_init_wups")
    {
      udata->dee_num_init_wups = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--dee_num_succ_wups")
    {
      udata->dee_num_succ_wups = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--dee_max_iters")
    {
      udata->dee_max_iters = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--dee_krylov_dim")
    {
      udata->dee_krylov_dim = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--dee_reltol")
    {
      udata->dee_reltol = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--timing") { udata->timing = true; }
    // Help
    else if (arg == "--help")
    {
      InputHelp();
      return -1;
    }
    // Unknown input
    else
    {
      cerr << "ERROR: Invalid input " << arg << endl;
      InputHelp();
      return -1;
    }
  }

  // Recompute total number of nodes
  udata->nodes = udata->nx * udata->ny;

  // Recompute x and y mesh spacing
  udata->dx = (udata->xu) / (udata->nx - 1);
  udata->dy = (udata->yu) / (udata->ny - 1);

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------

// Compute the exact solution
static int Solution(sunrealtype t, N_Vector u, UserData* udata)
{
  sunrealtype x, y;
  sunrealtype cos_sqr_t;
  sunrealtype sin_sqr_x, sin_sqr_y;

  // Constants for computing solution
  cos_sqr_t = cos(PI * t) * cos(PI * t);

  // Initialize u to zero (handles boundary conditions)
  N_VConst(ZERO, u);

  sunrealtype* uarray = N_VGetArrayPointer(u);
  if (check_flag((void*)uarray, "N_VGetArrayPointer", 0)) { return -1; }

  for (sunindextype j = 1; j < udata->ny - 1; j++)
  {
    for (sunindextype i = 1; i < udata->nx - 1; i++)
    {
      x = i * udata->dx;
      y = j * udata->dy;

      sin_sqr_x = sin(PI * x) * sin(PI * x);
      sin_sqr_y = sin(PI * y) * sin(PI * y);

      uarray[IDX(i, j, udata->nx)] = sin_sqr_x * sin_sqr_y * cos_sqr_t;
    }
  }

  return 0;
}

// Compute the solution error
static int SolutionError(sunrealtype t, N_Vector u, N_Vector e, UserData* udata)
{
  // Compute true solution
  int flag = Solution(t, e, udata);
  if (flag != 0) { return -1; }

  // Compute absolute error
  N_VLinearSum(ONE, u, -ONE, e, e);
  N_VAbs(e, e);

  return 0;
}

// Print command line options
static void InputHelp()
{
  cout << endl;
  cout << "Command line options:" << endl;
  cout
    << "  --mesh <nx> <ny>            : mesh points in the x and y directions"
    << endl;
  cout << "  --domain <xu> <yu>             : domain upper bound in the x and "
          "y direction"
       << endl;
  cout << "  --k <kx> <ky>               : diffusion coefficients" << endl;
  cout << "  --noforcing                 : disable forcing term" << endl;
  cout << "  --tf <time>                 : final time" << endl;
  cout << "  --rtol <rtol>               : relative tolerance" << endl;
  cout << "  --atol <atol>               : absolute tolerance" << endl;
  cout << "  --fixedstep <step>          : used fixed step size" << endl;
  cout << "  --controller <ctr>          : time step adaptivity controller"
       << endl;
  cout << "  --method <mth>              : LSRK method choice" << endl;
  cout << "  --eigfrequency <nst>        : dominant eigenvalue update frequency"
       << endl;
  cout << "  --stage_max_limit <smax>    : maximum number of stages per step"
       << endl;
  cout << "  --eigsafety <safety>        : dominant eigenvalue safety factor"
       << endl;
  cout << "  --diagnostics               : output diagnostics" << endl;
  cout << "  --output <level>            : output level" << endl;
  cout << "  --nout <nout>               : number of outputs" << endl;
  cout << "  --maxsteps <steps>          : max steps between outputs" << endl;
  cout << "  --dee_id <id>               : DomEig Estimator (DEE) id (PI: 0, "
          "Arnoldi: 1)"
       << endl;
  cout << "  --dee_num_init_wups <num>   : number of DEE initial warmups" << endl;
  cout << "  --dee_num_succ_wups <num>   : number of DEE succeeding warmups"
       << endl;
  cout << "  --dee_max_iters <num>       : max iterations in DEE" << endl;
  cout << "  --dee_krylov_dim <dim>      : Krylov dimension for DEE" << endl;
  cout << "  --dee_reltol <tol>          : DEE tolerance" << endl;
  cout << "  --timing                    : print timing data" << endl;
  cout << "  --help                      : print this message and exit" << endl;
}

// Print user data
static int PrintUserData(UserData* udata)
{
  cout << endl;
  cout << "2D Heat PDE test problem:" << endl;
  cout << " ------------------------------------ " << endl;
  cout << " kx                = " << udata->kx << endl;
  cout << " ky                = " << udata->ky << endl;
  cout << " forcing           = " << udata->forcing << endl;
  cout << " tf                = " << udata->tf << endl;
  cout << " xu                = " << udata->xu << endl;
  cout << " yu                = " << udata->yu << endl;
  cout << " nx                = " << udata->nx << endl;
  cout << " ny                = " << udata->ny << endl;
  cout << " dx                = " << udata->dx << endl;
  cout << " dy                = " << udata->dy << endl;
  cout << " ------------------------------------ " << endl;
  cout << " rtol              = " << udata->rtol << endl;
  cout << " atol              = " << udata->atol << endl;
  cout << " fixed h           = " << udata->hfixed << endl;
  cout << " controller        = " << udata->controller << endl;
  cout << " method            = " << udata->method << endl;
  cout << " eigfrequency      = " << udata->eigfrequency << endl;
  cout << " stage_max_limit   = " << udata->stage_max_limit << endl;
  cout << " eigsafety         = " << udata->eigsafety << endl;
  cout << " ------------------------------------ " << endl;
  cout << " output            = " << udata->output << endl;
  cout << " max steps         = " << udata->maxsteps << endl;
  cout << " ------------------------------------ " << endl;
  cout << " dee ID            = " << udata->dee_id << endl;
  cout << " dee num_init_wups = " << udata->dee_num_init_wups << endl;
  cout << " dee num_succ_wups = " << udata->dee_num_succ_wups << endl;
  cout << " dee_max_iters     = " << udata->dee_max_iters << endl;
  cout << " dee_krylov_dim    = " << udata->dee_krylov_dim << endl;
  cout << " dee_reltol        = " << udata->dee_reltol << endl;
  cout << " ------------------------------------ " << endl;
  cout << endl;

  return 0;
}

// Initialize output
static int OpenOutput(UserData* udata)
{
  // Header for status output
  if (udata->output > 0)
  {
    cout << scientific;
    cout << setprecision(numeric_limits<sunrealtype>::digits10);
    if (udata->forcing)
    {
      cout << "          t           ";
      cout << "          ||u||_rms      ";
      cout << "          max error      " << endl;
      cout << " ---------------------";
      cout << "-------------------------";
      cout << "-------------------------" << endl;
    }
    else
    {
      cout << "          t           ";
      cout << "          ||u||_rms      " << endl;
      cout << " ---------------------";
      cout << "-------------------------" << endl;
    }
  }

  // Output problem information and open output streams
  if (udata->output == 2)
  {
    // Each processor outputs subdomain information
    ofstream dout;
    dout.open("heat2d_info.txt");
    dout << "xu  " << udata->xu << endl;
    dout << "yu  " << udata->yu << endl;
    dout << "nx  " << udata->nx << endl;
    dout << "ny  " << udata->ny << endl;
    dout << "nt  " << udata->nout + 1 << endl;
    dout.close();

    // Open output streams for solution and error
    udata->uout.open("heat2d_solution.txt");
    udata->uout << scientific;
    udata->uout << setprecision(numeric_limits<sunrealtype>::digits10);

    if (udata->forcing)
    {
      udata->eout.open("heat2d_error.txt");
      udata->eout << scientific;
      udata->eout << setprecision(numeric_limits<sunrealtype>::digits10);
    }
  }

  return 0;
}

// Write output
static int WriteOutput(sunrealtype t, N_Vector u, UserData* udata)
{
  int flag;

  if (udata->output > 0)
  {
    // Compute rms norm of the state
    sunrealtype urms = sqrt(N_VDotProd(u, u) / udata->nx / udata->ny);

    // Output current status
    if (udata->forcing)
    {
      // Compute the error
      flag = SolutionError(t, u, udata->e, udata);
      if (check_flag(&flag, "SolutionError", 1)) { return 1; }

      // Compute max error
      sunrealtype max = N_VMaxNorm(udata->e);

      cout << setw(22) << t << setw(25) << urms << setw(25) << max << endl;
    }
    else { cout << setw(22) << t << setw(25) << urms << endl; }

    // Write solution and error to disk
    if (udata->output == 2)
    {
      sunrealtype* uarray = N_VGetArrayPointer(u);
      if (check_flag((void*)uarray, "N_VGetArrayPointer", 0)) { return -1; }

      udata->uout << t << " ";
      for (sunindextype i = 0; i < udata->nodes; i++)
      {
        udata->uout << uarray[i] << " ";
      }
      udata->uout << endl;

      if (udata->forcing)
      {
        // Output error to disk
        sunrealtype* earray = N_VGetArrayPointer(udata->e);
        if (check_flag((void*)earray, "N_VGetArrayPointer", 0)) { return -1; }

        udata->eout << t << " ";
        for (sunindextype i = 0; i < udata->nodes; i++)
        {
          udata->eout << earray[i] << " ";
        }
        udata->eout << endl;
      }
    }
  }

  return 0;
}

// Finalize output
static int CloseOutput(UserData* udata)
{
  // Footer for status output
  if (udata->output > 0)
  {
    if (udata->forcing)
    {
      cout << " ---------------------";
      cout << "-------------------------";
      cout << "-------------------------" << endl;
      cout << endl;
    }
    else
    {
      cout << " ---------------------";
      cout << "-------------------------" << endl;
      cout << endl;
    }
  }

  if (udata->output == 2)
  {
    // Close output streams
    udata->uout.close();
    if (udata->forcing) { udata->eout.close(); }
  }

  return 0;
}

static int OutputTiming(UserData* udata)
{
  cout << scientific;
  cout << setprecision(6);

  cout << "  Evolve time = " << udata->evolvetime << " sec" << endl;
  cout << "  RHS time    = " << udata->rhstime << " sec" << endl;
  cout << endl;

  return 0;
}

// Check function return value
static int check_flag(void* flagvalue, const string funcname, int opt)
{
  // Check if the function returned a NULL pointer
  if (opt == 0)
  {
    if (flagvalue == NULL)
    {
      cerr << endl
           << "ERROR: " << funcname << " returned NULL pointer" << endl
           << endl;
      return 1;
    }
  }
  // Check the function return flag value
  else if (opt == 1 || opt == 2)
  {
    int errflag = *((int*)flagvalue);
    if ((opt == 1 && errflag < 0) || (opt == 2 && errflag != 0))
    {
      cerr << endl
           << "ERROR: " << funcname << " returned with flag = " << errflag << endl
           << endl;
      return 1;
    }
  }
  else
  {
    cerr << endl
         << "ERROR: check_flag called with an invalid option value" << endl;
    return 1;
  }

  return 0;
}

//---- end of file ----
