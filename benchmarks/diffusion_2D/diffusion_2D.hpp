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
 * Shared header file for 2D diffusion benchmark problem
 * ---------------------------------------------------------------------------*/

#ifndef DIFFUSION_2D_HPP
#define DIFFUSION_2D_HPP

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <sundials/sundials_core.hpp>
#include <vector>

#include "mpi.h"

#if defined(USE_HIP)
#include "nvector/nvector_hip.h"
#include "nvector/nvector_mpiplusx.h"
#elif defined(USE_CUDA)
#include "nvector/nvector_cuda.h"
#include "nvector/nvector_mpiplusx.h"
#else
#include "nvector/nvector_parallel.h"
#endif

#include "sunlinsol/sunlinsol_pcg.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#if defined(USE_SUPERLU_DIST)
#include "sunlinsol/sunlinsol_superludist.h"
#include "sunmatrix/sunmatrix_slunrloc.h"
#endif

// Macros for problem constants
#define PI   SUN_RCONST(3.141592653589793238462643383279502884197169)
#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)

// Macro to access (x,y) location in 1D NVector array
#define IDX(x, y, n) ((n) * (y) + (x))

// Ceiling for integers ceil(a/b) = ceil((a + b - 1) / b)
#define ICEIL(a, b) (((a) + (b)-1) / (b))

using namespace std;

// -----------------------------------------------------------------------------
// UserData structure
// -----------------------------------------------------------------------------

struct UserData
{
  SUNProfiler prof = NULL;

  // Diffusion coefficients in the x and y directions
  sunrealtype kx = ONE;
  sunrealtype ky = ONE;

  // Enable/disable forcing
  bool forcing = true;

  // Final time
  sunrealtype tf = ONE;

  // Lower and Upper bounds in x and y directions
  sunrealtype xl = ZERO;
  sunrealtype yl = ZERO;
  sunrealtype xu = ONE;
  sunrealtype yu = ONE;

  // Global number of nodes in the x and y directions
  sunindextype nx = 32;
  sunindextype ny = 32;

  // Global total number of nodes
  sunindextype nodes = nx * ny;

  // Mesh spacing in the x and y directions
  sunrealtype dx = xu / (nx - 1);
  sunrealtype dy = yu / (ny - 1);

  // Minimum number of local nodes in the x and y directions
  sunindextype qx = 0;
  sunindextype qy = 0;

  // Leftover nodes in the x and y directions
  sunindextype rx = 0;
  sunindextype ry = 0;

  // Local number of nodes in the x and y directions
  sunindextype nx_loc = 0;
  sunindextype ny_loc = 0;

  // Overall number of local nodes
  sunindextype nodes_loc = 0;

  // Global x and y indices of this subdomain
  sunindextype is = 0; // x starting index
  sunindextype ie = 0; // x ending index
  sunindextype js = 0; // y starting index
  sunindextype je = 0; // y ending index

  // MPI variables
  MPI_Comm comm_c = MPI_COMM_NULL; // Cartesian communicator in space

  int np  = 1; // total number of MPI processes in Comm world
  int npx = 0; // number of MPI processes in the x-direction
  int npy = 0; // number of MPI processes in the y-direction

  int myid_c = 0; // process ID in Cartesian communicator
  int idx    = 0; // process x-coordinate
  int idy    = 0; // process y-coordinate

  // Flags denoting if this process has a neighbor
  bool HaveNbrW = true;
  bool HaveNbrE = true;
  bool HaveNbrS = true;
  bool HaveNbrN = true;

  // Neighbor IDs for exchange
  int ipW = -1;
  int ipE = -1;
  int ipS = -1;
  int ipN = -1;

  // Receive buffers for neighbor exchange
  sunrealtype* Wrecv = NULL;
  sunrealtype* Erecv = NULL;
  sunrealtype* Srecv = NULL;
  sunrealtype* Nrecv = NULL;

  // Receive requests for neighbor exchange
  MPI_Request reqRW;
  MPI_Request reqRE;
  MPI_Request reqRS;
  MPI_Request reqRN;

  // Send buffers for neighbor exchange
  sunrealtype* Wsend = NULL;
  sunrealtype* Esend = NULL;
  sunrealtype* Ssend = NULL;
  sunrealtype* Nsend = NULL;

  // Send requests for neighbor exchange
  MPI_Request reqSW;
  MPI_Request reqSE;
  MPI_Request reqSS;
  MPI_Request reqSN;

  // Inverse of Jacobian diagonal for preconditioner
  N_Vector diag = NULL;

  UserData(SUNProfiler prof_) : prof(prof_) {}

  ~UserData();

  // Helper functions
  int parse_args(vector<string>& args, bool outproc);
  void help();
  void print();
  int setup();
  int start_exchange(const N_Vector u);
  int end_exchange();

private:
  int allocate_buffers();
  int pack_buffers(const N_Vector u);
  int free_buffers();
};

// -----------------------------------------------------------------------------
// UserOutput structure
// -----------------------------------------------------------------------------

struct UserOutput
{
  // Output variables
  int output     = 1;    // 0 = no output, 1 = stats output, 2 = output to disk
  int nout       = 20;   // number of output times
  N_Vector error = NULL; // error vector
  ofstream uoutstream;   // output file stream
  ofstream eoutstream;   // error file stream

  // Helper functions
  int parse_args(vector<string>& args, bool outproc);
  void help();
  void print();

  // Output functions
  int open(UserData* udata);
  int write(sunrealtype t, N_Vector u, UserData* udata);
  int close(UserData* udata);
};

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrator
// -----------------------------------------------------------------------------

// Common function for computing the Laplacian
int laplacian(sunrealtype t, N_Vector u, N_Vector f, UserData* udata);

#if defined(BENCHMARK_ODE)

#if defined(USE_SUPERLU_DIST)
int laplacian_matrix_sludist(N_Vector u, SUNMatrix L, UserData* udata);
#endif

// ODE right hand side function
int diffusion(sunrealtype t, N_Vector u, N_Vector f, void* user_data);

int diffusion_jac(sunrealtype t, N_Vector u, N_Vector f, SUNMatrix Jac,
                  void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Preconditioner setup and solve functions
int PSetup(sunrealtype t, N_Vector u, N_Vector f, sunbooleantype jok,
           sunbooleantype* jcurPtr, sunrealtype gamma, void* user_data);

int PSolve(sunrealtype t, N_Vector u, N_Vector f, N_Vector r, N_Vector z,
           sunrealtype gamma, sunrealtype delta, int lr, void* user_data);

#elif defined(BENCHMARK_DAE)

#if defined(USE_SUPERLU_DIST)
int laplacian_matrix_sludist(N_Vector u, sunrealtype cj, SUNMatrix L,
                             UserData* udata);
#endif

// DAE residual function
int diffusion(sunrealtype t, N_Vector u, N_Vector up, N_Vector res,
              void* user_data);

int diffusion_jac(sunrealtype t, sunrealtype cj, N_Vector u, N_Vector up,
                  N_Vector res, SUNMatrix Jac, void* user_data, N_Vector tmp1,
                  N_Vector tmp2, N_Vector tmp3);

// Preconditioner setup and solve functions
int PSetup(sunrealtype t, N_Vector u, N_Vector up, N_Vector res, sunrealtype cj,
           void* user_data);

int PSolve(sunrealtype t, N_Vector u, N_Vector up, N_Vector res, N_Vector r,
           N_Vector z, sunrealtype cj, sunrealtype delta, void* user_data);

#else
#error "Missing ODE/DAE preprocessor directive"
#endif

// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Synchronize devices
int DeviceSynchronize();

// Copy data from devices
int CopyDataFromDevice(N_Vector y);

// Compute the true solution and its derivative
int Solution(sunrealtype t, N_Vector u, UserData* udata);
int SolutionDerivative(sunrealtype t, N_Vector up, UserData* udata);

// Compute the solution error
int SolutionError(sunrealtype t, N_Vector u, N_Vector e, UserData* udata);

// Check function return values
int check_flag(const void* flagvalue, const string funcname, int opt);

#endif
