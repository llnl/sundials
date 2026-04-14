/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *
 * (adapted from ark_heat2D_p.cpp, co-authored by Daniel Reynolds and David
 * Gardner (LLNL))
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

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

#include "arkode/arkode_lsrkstep.h"   // access to LSRKStep
#include "mpi.h"                      // MPI header file
#include "nvector/nvector_parallel.h" // access to the MPI N_Vector
#include "sunadaptcontroller/sunadaptcontroller_imexgus.h"
#include "sunadaptcontroller/sunadaptcontroller_soderlind.h"

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

  // Global number of nodes in the x and y directions
  sunindextype nx;
  sunindextype ny;

  // Global total number of nodes
  sunindextype nodes;

  // Mesh spacing in the x and y directions
  sunrealtype dx;
  sunrealtype dy;

  // Local number of nodes in the x and y directions
  sunindextype nx_loc;
  sunindextype ny_loc;

  // Overall number of local nodes
  sunindextype nodes_loc;

  // Global x and y indices of this subdomain
  sunindextype is; // x starting index
  sunindextype ie; // x ending index
  sunindextype js; // y starting index
  sunindextype je; // y ending index

  // MPI variables
  MPI_Comm comm_c; // Cartesian communicator in space

  int nprocs_w; // total number of MPI processes in Comm world
  int npx;      // number of MPI processes in the x-direction
  int npy;      // number of MPI processes in the y-direction

  int myid_c; // process ID in Cartesian communicator

  // Flags denoting if this process has a neighbor
  bool HaveNbrW;
  bool HaveNbrE;
  bool HaveNbrS;
  bool HaveNbrN;

  // Neighbor IDs for exchange
  int ipW;
  int ipE;
  int ipS;
  int ipN;

  // Receive buffers for neighbor exchange
  sunrealtype* Wrecv;
  sunrealtype* Erecv;
  sunrealtype* Srecv;
  sunrealtype* Nrecv;

  // Receive requests for neighbor exchange
  MPI_Request reqRW;
  MPI_Request reqRE;
  MPI_Request reqRS;
  MPI_Request reqRN;

  // Send buffers for neighbor exchange
  sunrealtype* Wsend;
  sunrealtype* Esend;
  sunrealtype* Ssend;
  sunrealtype* Nsend;

  // Send requests for neighbor exchange
  MPI_Request reqSW;
  MPI_Request reqSE;
  MPI_Request reqSS;
  MPI_Request reqSN;

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

  // Timing variables
  bool timing; // print timings
  double evolvetime;
  double rhstime;
  double exchangetime;
};

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrator
// -----------------------------------------------------------------------------

// ODE right hand side function
static int f(sunrealtype t, N_Vector u, N_Vector f, void* user_data);

// Spectral radius estimation routine
static int eig(sunrealtype t, N_Vector y, N_Vector fn, sunrealtype* lambdaR,
               sunrealtype* lambdaI, void* user_data, N_Vector temp1,
               N_Vector temp2, N_Vector temp3);

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// Setup the parallel decomposition
static int SetupDecomp(MPI_Comm comm_w, UserData* udata);

// Perform neighbor exchange
static int PostRecv(UserData* udata);
static int SendData(N_Vector y, UserData* udata);
static int WaitRecv(UserData* udata);

// -----------------------------------------------------------------------------
// UserData and input functions
// -----------------------------------------------------------------------------

// Set the default values in the UserData structure
static int InitUserData(UserData* udata);

// Free memory allocated within UserData
static int FreeUserData(UserData* udata);

// Read the command line inputs and set UserData values
static int ReadInputs(int* argc, char*** argv, UserData* udata, bool outproc);

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

  // Timing variables
  double t1 = 0.0;
  double t2 = 0.0;

  // MPI variables
  MPI_Comm comm_w = MPI_COMM_WORLD; // MPI communicator
  int myid;                         // MPI process ID

  // Initialize MPI
  flag = MPI_Init(&argc, &argv);
  if (check_flag(&flag, "MPI_Init", 1)) { return 1; }

  flag = MPI_Comm_rank(comm_w, &myid);
  if (check_flag(&flag, "MPI_Comm_rank", 1)) { return 1; }

  // Create the SUNDIALS context object for this simulation
  SUNContext ctx;
  flag = SUNContext_Create(comm_w, &ctx);
  if (check_flag(&flag, "SUNContext_Create", 1)) { return 1; }

  // Set output process flag
  bool outproc = (myid == 0);

  // ------------------------------------------
  // Setup UserData and parallel decomposition
  // ------------------------------------------

  // Allocate and initialize user data structure with default values. The
  // defaults may be overwritten by command line inputs in ReadInputs below.
  udata = new UserData;
  flag  = InitUserData(udata);
  if (check_flag(&flag, "InitUserData", 1)) { return 1; }

  // Parse command line inputs
  flag = ReadInputs(&argc, &argv, udata, outproc);
  if (flag != 0) { return 1; }

  // Setup parallel decomposition
  flag = SetupDecomp(comm_w, udata);
  if (check_flag(&flag, "SetupDecomp", 1)) { return 1; }

  // Output problem setup/options
  if (outproc)
  {
    flag = PrintUserData(udata);
    if (check_flag(&flag, "PrintUserData", 1)) { return 1; }
  }

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

  // ------------------------
  // Create parallel vectors
  // ------------------------

  // Create vector for solution
  u = N_VNew_Parallel(udata->comm_c, udata->nodes_loc, udata->nodes, ctx);
  if (check_flag((void*)u, "N_VNew_Parallel", 0)) { return 1; }

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

  // Select LSRK spectral radius function and options
  flag = LSRKStepSetDomEigFn(arkode_mem, eig);
  if (check_flag(&flag, "LSRKStepSetDomEigFn", 1)) { return 1; }
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

  // Initial output
  flag = OpenOutput(udata);
  if (check_flag(&flag, "OpenOutput", 1)) { return 1; }

  flag = WriteOutput(t, u, udata);
  if (check_flag(&flag, "WriteOutput", 1)) { return 1; }

  for (int iout = 0; iout < udata->nout; iout++)
  {
    // Start timer
    t1 = MPI_Wtime();

    // Evolve in time
    flag = ARKodeEvolve(arkode_mem, tout, u, &t, ARK_NORMAL);
    if (check_flag(&flag, "ARKodeEvolve", 1)) { break; }

    // Stop timer
    t2 = MPI_Wtime();

    // Update timer
    udata->evolvetime += t2 - t1;

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
  if (udata->output > 0 && outproc)
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

    if (outproc)
    {
      cout << scientific;
      cout << setprecision(numeric_limits<sunrealtype>::digits10);
      cout << "  Max error = " << maxerr << endl;
    }
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

  ARKodeFree(&arkode_mem); // Free integrator memory
  N_VDestroy(u);           // Free vectors
  FreeUserData(udata);     // Free user data
  delete udata;
  SUNContext_Free(&ctx); // Free context
  flag = MPI_Finalize(); // Finalize MPI
  return 0;
}

// -----------------------------------------------------------------------------
// Setup the parallel decomposition
// -----------------------------------------------------------------------------

static int SetupDecomp(MPI_Comm comm_w, UserData* udata)
{
  int flag;

  // Check that this has not been called before
  if (udata->Erecv != NULL || udata->Wrecv != NULL || udata->Srecv != NULL ||
      udata->Nrecv != NULL)
  {
    cerr << "SetupDecomp error: parallel decomposition already set up" << endl;
    return -1;
  }

  // Get the number of processes
  flag = MPI_Comm_size(comm_w, &(udata->nprocs_w));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_size = " << flag << endl;
    return -1;
  }

  // Check the processor grid
  if ((udata->npx * udata->npy) != udata->nprocs_w)
  {
    cerr << "Error: npx * npy != nproc" << endl;
    return -1;
  }

  // Set up 2D Cartesian communicator
  int dims[2];
  dims[0] = udata->npx;
  dims[1] = udata->npy;

  int periods[2];
  periods[0] = 0;
  periods[1] = 0;

  flag = MPI_Cart_create(comm_w, 2, dims, periods, 0, &(udata->comm_c));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_create = " << flag << endl;
    return -1;
  }

  // Get my rank in the new Cartesian communicator
  flag = MPI_Comm_rank(udata->comm_c, &(udata->myid_c));
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_rank = " << flag << endl;
    return -1;
  }

  // Get dimension of the Cartesian communicator and my coordinates
  int coords[2];
  flag = MPI_Cart_get(udata->comm_c, 2, dims, periods, coords);
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_get = " << flag << endl;
    return -1;
  }

  // Determine local extents in x-direction
  int idx         = coords[0];
  sunindextype qx = udata->nx / dims[0];
  sunindextype rx = udata->nx % dims[0];

  udata->is = qx * idx + (idx < rx ? idx : rx);
  udata->ie = udata->is + qx - 1 + (idx < rx ? 1 : 0);

  // Sanity check
  if (udata->ie > (udata->nx - 1))
  {
    cerr << "Error ie > nx - 1" << endl;
    return -1;
  }

  // Determine local extents in y-direction
  int idy         = coords[1];
  sunindextype qy = udata->ny / dims[1];
  sunindextype ry = udata->ny % dims[1];

  udata->js = qy * idy + (idy < ry ? idy : ry);
  udata->je = udata->js + qy - 1 + (idy < ry ? 1 : 0);

  // Sanity check
  if (udata->je > (udata->ny - 1))
  {
    cerr << "Error je > ny - 1" << endl;
    return -1;
  }

  // Number of local nodes
  udata->nx_loc = (udata->ie) - (udata->is) + 1;
  udata->ny_loc = (udata->je) - (udata->js) + 1;

  // Initialize global and local vector lengths
  udata->nodes     = udata->nx * udata->ny;
  udata->nodes_loc = udata->nx_loc * udata->ny_loc;

  // Determine if this proc has neighbors
  udata->HaveNbrW = (udata->is != 0);
  udata->HaveNbrE = (udata->ie != udata->nx - 1);
  udata->HaveNbrS = (udata->js != 0);
  udata->HaveNbrN = (udata->je != udata->ny - 1);

  // Allocate exchange buffers if necessary
  if (udata->HaveNbrW)
  {
    udata->Wrecv = new sunrealtype[udata->ny_loc];
    udata->Wsend = new sunrealtype[udata->ny_loc];
  }
  if (udata->HaveNbrE)
  {
    udata->Erecv = new sunrealtype[udata->ny_loc];
    udata->Esend = new sunrealtype[udata->ny_loc];
  }
  if (udata->HaveNbrS)
  {
    udata->Srecv = new sunrealtype[udata->nx_loc];
    udata->Ssend = new sunrealtype[udata->nx_loc];
  }
  if (udata->HaveNbrN)
  {
    udata->Nrecv = new sunrealtype[udata->nx_loc];
    udata->Nsend = new sunrealtype[udata->nx_loc];
  }

  // MPI neighborhood information
  int nbcoords[2];

  // West neighbor
  if (udata->HaveNbrW)
  {
    nbcoords[0] = coords[0] - 1;
    nbcoords[1] = coords[1];
    flag        = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipW));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // East neighbor
  if (udata->HaveNbrE)
  {
    nbcoords[0] = coords[0] + 1;
    nbcoords[1] = coords[1];
    flag        = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipE));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // South neighbor
  if (udata->HaveNbrS)
  {
    nbcoords[0] = coords[0];
    nbcoords[1] = coords[1] - 1;
    flag        = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipS));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // North neighbor
  if (udata->HaveNbrN)
  {
    nbcoords[0] = coords[0];
    nbcoords[1] = coords[1] + 1;
    flag        = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipN));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// Functions called by the integrator
// -----------------------------------------------------------------------------

// f routine to compute the ODE RHS function f(t,y).
static int f(sunrealtype t, N_Vector u, N_Vector f, void* user_data)
{
  int flag;
  sunindextype i, j;

  // Start timer
  double t1 = MPI_Wtime();

  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Open exchange receives
  flag = PostRecv(udata);
  if (check_flag(&flag, "PostRecv", 1)) { return -1; }

  // Send exchange data
  flag = SendData(u, udata);
  if (check_flag(&flag, "SendData", 1)) { return -1; }

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0 : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0 : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

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

  // Iterate over subdomain and compute rhs forcing term
  if (udata->forcing)
  {
    sunrealtype x, y;
    sunrealtype sin_sqr_x, sin_sqr_y;
    sunrealtype cos_sqr_x, cos_sqr_y;

    sunrealtype bx = (udata->kx) * TWO * PI * PI;
    sunrealtype by = (udata->ky) * TWO * PI * PI;

    sunrealtype sin_t_cos_t = sin(PI * t) * cos(PI * t);
    sunrealtype cos_sqr_t   = cos(PI * t) * cos(PI * t);

    for (j = jstart; j < jend; j++)
    {
      for (i = istart; i < iend; i++)
      {
        x = (udata->is + i) * udata->dx;
        y = (udata->js + j) * udata->dy;

        sin_sqr_x = sin(PI * x) * sin(PI * x);
        sin_sqr_y = sin(PI * y) * sin(PI * y);

        cos_sqr_x = cos(PI * x) * cos(PI * x);
        cos_sqr_y = cos(PI * y) * cos(PI * y);

        farray[IDX(i, j, nx_loc)] =
          -TWO * PI * sin_sqr_x * sin_sqr_y * sin_t_cos_t -
          bx * (cos_sqr_x - sin_sqr_x) * sin_sqr_y * cos_sqr_t -
          by * (cos_sqr_y - sin_sqr_y) * sin_sqr_x * cos_sqr_t;
      }
    }
  }

  // Iterate over subdomain interior and add rhs diffusion term
  for (j = 1; j < ny_loc - 1; j++)
  {
    for (i = 1; i < nx_loc - 1; i++)
    {
      farray[IDX(i, j, nx_loc)] +=
        cc * uarray[IDX(i, j, nx_loc)] +
        cx * (uarray[IDX(i - 1, j, nx_loc)] + uarray[IDX(i + 1, j, nx_loc)]) +
        cy * (uarray[IDX(i, j - 1, nx_loc)] + uarray[IDX(i, j + 1, nx_loc)]);
    }
  }

  // Wait for exchange receives
  flag = WaitRecv(udata);
  if (check_flag(&flag, "WaitRecv", 1)) { return -1; }

  // Iterate over subdomain boundaries and add rhs diffusion term
  sunrealtype* Warray = udata->Wrecv;
  sunrealtype* Earray = udata->Erecv;
  sunrealtype* Sarray = udata->Srecv;
  sunrealtype* Narray = udata->Nrecv;

  // West face (updates south-west and north-west corners if necessary)
  if (udata->HaveNbrW)
  {
    i = 0;
    if (udata->HaveNbrS) // South-West corner
    {
      j = 0;
      farray[IDX(i, j, nx_loc)] +=
        cc * uarray[IDX(i, j, nx_loc)] +
        cx * (Warray[j] + uarray[IDX(i + 1, j, nx_loc)]) +
        cy * (Sarray[i] + uarray[IDX(i, j + 1, nx_loc)]);
    }

    for (j = 1; j < ny_loc - 1; j++)
    {
      farray[IDX(i, j, nx_loc)] +=
        cc * uarray[IDX(i, j, nx_loc)] +
        cx * (Warray[j] + uarray[IDX(i + 1, j, nx_loc)]) +
        cy * (uarray[IDX(i, j - 1, nx_loc)] + uarray[IDX(i, j + 1, nx_loc)]);
    }

    if (udata->HaveNbrN) // North-West corner
    {
      j = ny_loc - 1;
      farray[IDX(i, j, nx_loc)] +=
        cc * uarray[IDX(i, j, nx_loc)] +
        cx * (Warray[j] + uarray[IDX(i + 1, j, nx_loc)]) +
        cy * (uarray[IDX(i, j - 1, nx_loc)] + Narray[i]);
    }
  }

  // East face (updates south-east and north-east corners if necessary)
  if (udata->HaveNbrE)
  {
    i = nx_loc - 1;
    if (udata->HaveNbrS) // South-East corner
    {
      j = 0;
      farray[IDX(i, j, nx_loc)] +=
        cc * uarray[IDX(i, j, nx_loc)] +
        cx * (uarray[IDX(i - 1, j, nx_loc)] + Earray[j]) +
        cy * (Sarray[i] + uarray[IDX(i, j + 1, nx_loc)]);
    }

    for (j = 1; j < ny_loc - 1; j++)
    {
      farray[IDX(i, j, nx_loc)] +=
        cc * uarray[IDX(i, j, nx_loc)] +
        cx * (uarray[IDX(i - 1, j, nx_loc)] + Earray[j]) +
        cy * (uarray[IDX(i, j - 1, nx_loc)] + uarray[IDX(i, j + 1, nx_loc)]);
    }

    if (udata->HaveNbrN) // North-East corner
    {
      j = ny_loc - 1;
      farray[IDX(i, j, nx_loc)] +=
        cc * uarray[IDX(i, j, nx_loc)] +
        cx * (uarray[IDX(i - 1, j, nx_loc)] + Earray[j]) +
        cy * (uarray[IDX(i, j - 1, nx_loc)] + Narray[i]);
    }
  }

  // South face (excludes corners)
  if (udata->HaveNbrS)
  {
    j = 0;
    for (i = 1; i < nx_loc - 1; i++)
    {
      farray[IDX(i, j, nx_loc)] +=
        cc * uarray[IDX(i, j, nx_loc)] +
        cx * (uarray[IDX(i - 1, j, nx_loc)] + uarray[IDX(i + 1, j, nx_loc)]) +
        cy * (Sarray[i] + uarray[IDX(i, j + 1, nx_loc)]);
    }
  }

  // North face (excludes corners)
  if (udata->HaveNbrN)
  {
    j = udata->ny_loc - 1;
    for (i = 1; i < nx_loc - 1; i++)
    {
      farray[IDX(i, j, nx_loc)] +=
        cc * uarray[IDX(i, j, nx_loc)] +
        cx * (uarray[IDX(i - 1, j, nx_loc)] + uarray[IDX(i + 1, j, nx_loc)]) +
        cy * (uarray[IDX(i, j - 1, nx_loc)] + Narray[i]);
    }
  }

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->rhstime += t2 - t1;

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// RHS helper functions
// -----------------------------------------------------------------------------

// Post exchange receives
static int PostRecv(UserData* udata)
{
  int flag;

  // Start timer
  double t1 = MPI_Wtime();

  // Open Irecv buffers
  if (udata->HaveNbrW)
  {
    flag = MPI_Irecv(udata->Wrecv, (int)udata->ny_loc, MPI_SUNREALTYPE,
                     udata->ipW, MPI_ANY_TAG, udata->comm_c, &(udata->reqRW));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrE)
  {
    flag = MPI_Irecv(udata->Erecv, (int)udata->ny_loc, MPI_SUNREALTYPE,
                     udata->ipE, MPI_ANY_TAG, udata->comm_c, &(udata->reqRE));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrS)
  {
    flag = MPI_Irecv(udata->Srecv, (int)udata->nx_loc, MPI_SUNREALTYPE,
                     udata->ipS, MPI_ANY_TAG, udata->comm_c, &(udata->reqRS));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrN)
  {
    flag = MPI_Irecv(udata->Nrecv, (int)udata->nx_loc, MPI_SUNREALTYPE,
                     udata->ipN, MPI_ANY_TAG, udata->comm_c, &(udata->reqRN));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->exchangetime += t2 - t1;

  // Return success
  return 0;
}

// Send exchange data
static int SendData(N_Vector y, UserData* udata)
{
  int flag, i;
  sunindextype ny_loc = udata->ny_loc;
  sunindextype nx_loc = udata->nx_loc;

  // Start timer
  double t1 = MPI_Wtime();

  // Access data array
  sunrealtype* Y = N_VGetArrayPointer(y);
  if (check_flag((void*)Y, "N_VGetArrayPointer", 0)) { return -1; }

  // Send data
  if (udata->HaveNbrW)
  {
    for (i = 0; i < ny_loc; i++) { udata->Wsend[i] = Y[IDX(0, i, nx_loc)]; }
    flag = MPI_Isend(udata->Wsend, (int)udata->ny_loc, MPI_SUNREALTYPE,
                     udata->ipW, 0, udata->comm_c, &(udata->reqSW));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrE)
  {
    for (i = 0; i < ny_loc; i++)
    {
      udata->Esend[i] = Y[IDX(nx_loc - 1, i, nx_loc)];
    }
    flag = MPI_Isend(udata->Esend, (int)udata->ny_loc, MPI_SUNREALTYPE,
                     udata->ipE, 1, udata->comm_c, &(udata->reqSE));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrS)
  {
    for (i = 0; i < nx_loc; i++) { udata->Ssend[i] = Y[IDX(i, 0, nx_loc)]; }
    flag = MPI_Isend(udata->Ssend, (int)udata->nx_loc, MPI_SUNREALTYPE,
                     udata->ipS, 2, udata->comm_c, &(udata->reqSS));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrN)
  {
    for (i = 0; i < nx_loc; i++)
    {
      udata->Nsend[i] = Y[IDX(i, ny_loc - 1, nx_loc)];
    }
    flag = MPI_Isend(udata->Nsend, (int)udata->nx_loc, MPI_SUNREALTYPE,
                     udata->ipN, 3, udata->comm_c, &(udata->reqSN));
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->exchangetime += t2 - t1;

  // Return success
  return 0;
}

// Wait for exchange data
static int WaitRecv(UserData* udata)
{
  // Local variables
  int flag;
  MPI_Status stat;

  // Start timer
  double t1 = MPI_Wtime();

  // Wait for messages to finish
  if (udata->HaveNbrW)
  {
    flag = MPI_Wait(&(udata->reqRW), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
    flag = MPI_Wait(&(udata->reqSW), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrE)
  {
    flag = MPI_Wait(&(udata->reqRE), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
    flag = MPI_Wait(&(udata->reqSE), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrS)
  {
    flag = MPI_Wait(&(udata->reqRS), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
    flag = MPI_Wait(&(udata->reqSS), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  if (udata->HaveNbrN)
  {
    flag = MPI_Wait(&(udata->reqRN), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
    flag = MPI_Wait(&(udata->reqSN), &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->exchangetime += t2 - t1;

  // Return success
  return 0;
}

// Spectral radius estimation routine
static int eig(sunrealtype t, N_Vector y, N_Vector fn, sunrealtype* lambdaR,
               sunrealtype* lambdaI, void* user_data, N_Vector temp1,
               N_Vector temp2, N_Vector temp3)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Fill in spectral radius value
  *lambdaR = -SUN_RCONST(8.0) * max(udata->kx / udata->dx / udata->dx,
                                    udata->ky / udata->dy / udata->dy);
  *lambdaI = SUN_RCONST(0.0);

  // return with success
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

  // Global number of nodes in the x and y directions
  udata->nx    = 64;
  udata->ny    = 64;
  udata->nodes = udata->nx * udata->ny;

  // Mesh spacing in the x and y directions
  udata->dx = udata->xu / (udata->nx - 1);
  udata->dy = udata->yu / (udata->ny - 1);

  // Locals number of nodes in the x and y directions (set in SetupDecomp)
  udata->nx_loc    = 0;
  udata->ny_loc    = 0;
  udata->nodes_loc = 0;

  // Global indices of this subdomain (set in SetupDecomp)
  udata->is = 0;
  udata->ie = 0;
  udata->js = 0;
  udata->je = 0;

  // MPI variables (set in SetupDecomp)
  udata->comm_c = MPI_COMM_NULL;

  udata->nprocs_w = 1;
  udata->npx      = 1;
  udata->npy      = 1;

  udata->myid_c = 0;

  // Flags denoting neighbors (set in SetupDecomp)
  udata->HaveNbrW = true;
  udata->HaveNbrE = true;
  udata->HaveNbrS = true;
  udata->HaveNbrN = true;

  // Exchange receive buffers (allocated in SetupDecomp)
  udata->Erecv = NULL;
  udata->Wrecv = NULL;
  udata->Nrecv = NULL;
  udata->Srecv = NULL;

  // Exchange send buffers (allocated in SetupDecomp)
  udata->Esend = NULL;
  udata->Wsend = NULL;
  udata->Nsend = NULL;
  udata->Ssend = NULL;

  // Neighbors IDs (set in SetupDecomp)
  udata->ipW = -1;
  udata->ipE = -1;
  udata->ipS = -1;
  udata->ipN = -1;

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

  // Timing variables
  udata->timing       = false;
  udata->evolvetime   = 0.0;
  udata->rhstime      = 0.0;
  udata->exchangetime = 0.0;

  // Return success
  return 0;
}

// Free memory allocated within Userdata
static int FreeUserData(UserData* udata)
{
  // Free exchange buffers
  if (udata->Wrecv != NULL) { delete[] udata->Wrecv; }
  if (udata->Wsend != NULL) { delete[] udata->Wsend; }
  if (udata->Erecv != NULL) { delete[] udata->Erecv; }
  if (udata->Esend != NULL) { delete[] udata->Esend; }
  if (udata->Srecv != NULL) { delete[] udata->Srecv; }
  if (udata->Ssend != NULL) { delete[] udata->Ssend; }
  if (udata->Nrecv != NULL) { delete[] udata->Nrecv; }
  if (udata->Nsend != NULL) { delete[] udata->Nsend; }

  // Free MPI Cartesian communicator
  if (udata->comm_c != MPI_COMM_NULL) { MPI_Comm_free(&(udata->comm_c)); }

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
static int ReadInputs(int* argc, char*** argv, UserData* udata, bool outproc)
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
    // MPI processes
    else if (arg == "--np")
    {
      udata->npx = stoi((*argv)[arg_idx++]);
      udata->npy = stoi((*argv)[arg_idx++]);
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
    else if (arg == "--timing") { udata->timing = true; }
    // Help
    else if (arg == "--help")
    {
      if (outproc) { InputHelp(); }
      return -1;
    }
    // Unknown input
    else
    {
      if (outproc)
      {
        cerr << "ERROR: Invalid input " << arg << endl;
        InputHelp();
      }
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

  // Iterative over domain interior
  sunindextype istart = (udata->HaveNbrW) ? 0 : 1;
  sunindextype iend   = (udata->HaveNbrE) ? udata->nx_loc : udata->nx_loc - 1;

  sunindextype jstart = (udata->HaveNbrS) ? 0 : 1;
  sunindextype jend   = (udata->HaveNbrN) ? udata->ny_loc : udata->ny_loc - 1;

  sunrealtype* uarray = N_VGetArrayPointer(u);
  if (check_flag((void*)uarray, "N_VGetArrayPointer", 0)) { return -1; }

  for (sunindextype j = jstart; j < jend; j++)
  {
    for (sunindextype i = istart; i < iend; i++)
    {
      x = (udata->is + i) * udata->dx;
      y = (udata->js + j) * udata->dy;

      sin_sqr_x = sin(PI * x) * sin(PI * x);
      sin_sqr_y = sin(PI * y) * sin(PI * y);

      uarray[IDX(i, j, udata->nx_loc)] = sin_sqr_x * sin_sqr_y * cos_sqr_t;
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
  cout << "  --mesh <nx> <ny>        : mesh points in the x and y directions"
       << endl;
  cout << "  --np <npx> <npy>        : number of MPI processes in the x and y "
          "directions"
       << endl;
  cout
    << "  --domain <xu> <yu>      : domain upper bound in the x and y direction"
    << endl;
  cout << "  --k <kx> <ky>           : diffusion coefficients" << endl;
  cout << "  --noforcing             : disable forcing term" << endl;
  cout << "  --tf <time>             : final time" << endl;
  cout << "  --rtol <rtol>           : relative tolerance" << endl;
  cout << "  --atol <atol>           : absolute tolerance" << endl;
  cout << "  --fixedstep <step>      : used fixed step size" << endl;
  cout << "  --controller <ctr>      : time step adaptivity controller" << endl;
  cout << "  --method <mth>          : LSRK method choice" << endl;
  cout << "  --eigfrequency <nst>    : dominant eigenvalue update frequency"
       << endl;
  cout << "  --stage_max_limit <smax>  : maximum number of stages per step"
       << endl;
  cout << "  --eigsafety <safety>    : dominant eigenvalue safety factor" << endl;
  cout << "  --diagnostics           : output diagnostics" << endl;
  cout << "  --output <level>        : output level" << endl;
  cout << "  --nout <nout>           : number of outputs" << endl;
  cout << "  --maxsteps <steps>      : max steps between outputs" << endl;
  cout << "  --timing                : print timing data" << endl;
  cout << "  --help                  : print this message and exit" << endl;
}

// Print user data
static int PrintUserData(UserData* udata)
{
  cout << endl;
  cout << "2D Heat PDE test problem:" << endl;
  cout << " --------------------------------- " << endl;
  cout << "  nprocs         = " << udata->nprocs_w << endl;
  cout << "  npx            = " << udata->npx << endl;
  cout << "  npy            = " << udata->npy << endl;
  cout << " --------------------------------- " << endl;
  cout << "  kx             = " << udata->kx << endl;
  cout << "  ky             = " << udata->ky << endl;
  cout << "  forcing        = " << udata->forcing << endl;
  cout << "  tf             = " << udata->tf << endl;
  cout << "  xu             = " << udata->xu << endl;
  cout << "  yu             = " << udata->yu << endl;
  cout << "  nx             = " << udata->nx << endl;
  cout << "  ny             = " << udata->ny << endl;
  cout << "  nxl (proc 0)   = " << udata->nx_loc << endl;
  cout << "  nyl (proc 0)   = " << udata->ny_loc << endl;
  cout << "  dx             = " << udata->dx << endl;
  cout << "  dy             = " << udata->dy << endl;
  cout << " --------------------------------- " << endl;
  cout << "  rtol           = " << udata->rtol << endl;
  cout << "  atol           = " << udata->atol << endl;
  cout << "  fixed h        = " << udata->hfixed << endl;
  cout << "  controller     = " << udata->controller << endl;
  cout << "  method         = " << udata->method << endl;
  cout << "  eigfrequency   = " << udata->eigfrequency << endl;
  cout << "  stage_max_limit  = " << udata->stage_max_limit << endl;
  cout << "  eigsafety      = " << udata->eigsafety << endl;
  cout << " --------------------------------- " << endl;
  cout << "  output         = " << udata->output << endl;
  cout << "  max steps      = " << udata->maxsteps << endl;
  cout << " --------------------------------- " << endl;
  cout << endl;

  return 0;
}

// Initialize output
static int OpenOutput(UserData* udata)
{
  bool outproc = (udata->myid_c == 0);

  // Header for status output
  if (udata->output > 0 && outproc)
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
    stringstream fname;
    fname << "heat2d_info." << setfill('0') << setw(5) << udata->myid_c << ".txt";

    ofstream dout;
    dout.open(fname.str());
    dout << "xu  " << udata->xu << endl;
    dout << "yu  " << udata->yu << endl;
    dout << "nx  " << udata->nx << endl;
    dout << "ny  " << udata->ny << endl;
    dout << "px  " << udata->npx << endl;
    dout << "py  " << udata->npy << endl;
    dout << "np  " << udata->nprocs_w << endl;
    dout << "is  " << udata->is << endl;
    dout << "ie  " << udata->ie << endl;
    dout << "js  " << udata->js << endl;
    dout << "je  " << udata->je << endl;
    dout << "nt  " << udata->nout + 1 << endl;
    dout.close();

    // Open output streams for solution and error
    fname.str("");
    fname.clear();
    fname << "heat2d_solution." << setfill('0') << setw(5) << udata->myid_c
          << ".txt";
    udata->uout.open(fname.str());

    udata->uout << scientific;
    udata->uout << setprecision(numeric_limits<sunrealtype>::digits10);

    if (udata->forcing)
    {
      fname.str("");
      fname.clear();
      fname << "heat2d_error." << setfill('0') << setw(5) << udata->myid_c
            << ".txt";
      udata->eout.open(fname.str());

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
  sunrealtype max = ZERO;
  bool outproc    = (udata->myid_c == 0);

  if (udata->output > 0)
  {
    if (udata->forcing)
    {
      // Compute the error
      flag = SolutionError(t, u, udata->e, udata);
      if (check_flag(&flag, "SolutionError", 1)) { return 1; }

      // Compute max error
      max = N_VMaxNorm(udata->e);
    }

    // Compute rms norm of the state
    sunrealtype urms = sqrt(N_VDotProd(u, u) / udata->nx / udata->ny);

    // Output current status
    if (outproc)
    {
      if (udata->forcing)
      {
        cout << setw(22) << t << setw(25) << urms << setw(25) << max << endl;
      }
      else { cout << setw(22) << t << setw(25) << urms << endl; }
    }

    // Write solution and error to disk
    if (udata->output == 2)
    {
      sunrealtype* uarray = N_VGetArrayPointer(u);
      if (check_flag((void*)uarray, "N_VGetArrayPointer", 0)) { return -1; }

      udata->uout << t << " ";
      for (sunindextype i = 0; i < udata->nodes_loc; i++)
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
        for (sunindextype i = 0; i < udata->nodes_loc; i++)
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
  bool outproc = (udata->myid_c == 0);

  // Footer for status output
  if (outproc && (udata->output > 0))
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
  bool outproc = (udata->myid_c == 0);

  if (outproc)
  {
    cout << scientific;
    cout << setprecision(6);
  }

  double maxtime = 0.0;

  MPI_Reduce(&(udata->evolvetime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc) { cout << "  Evolve time   = " << maxtime << " sec" << endl; }

  MPI_Reduce(&(udata->rhstime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc) { cout << "  RHS time      = " << maxtime << " sec" << endl; }

  MPI_Reduce(&(udata->exchangetime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  Exchange time = " << maxtime << " sec" << endl;
    cout << endl;
  }

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
