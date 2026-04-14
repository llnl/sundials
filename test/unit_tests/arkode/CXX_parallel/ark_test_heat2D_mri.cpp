/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *---------------------------------------------------------------
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
 *---------------------------------------------------------------
 * Routine to check that MRIStep and ARKStep exhibit the same
 * solver statistics when both run with fixed-steps, the same
 * solver parameters, and MRIStep runs using a solve-decoupled
 * DIRK method at the slow time scale.
 *
 * This routine will switch between the default Newton nonlinear
 * solver and the 'linear' version based on a 0/1 command-line
 * argument (1 => linear).
 *---------------------------------------------------------------*/

// Header files
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode/arkode_arkstep.h"    // prototypes for ARKStep fcts., consts
#include "arkode/arkode_mristep.h"    // prototypes for MRIStep fcts., consts
#include "mpi.h"                      // MPI header file
#include "nvector/nvector_parallel.h" // parallel N_Vector types, fcts., macros
#include "sundials/sundials_types.h"  // def. of type 'sunrealtype'
#include "sunlinsol/sunlinsol_pcg.h"  // access to PCG SUNLinearSolver

using namespace std;

// accessor macros between (x,y) location and 1D NVector array
#define IDX(x, y, n) ((n) * (y) + (x))
#define PI           SUN_RCONST(3.141592653589793238462643383279502884197169)
#define ZERO         SUN_RCONST(0.0)
#define ONE          SUN_RCONST(1.0)
#define TWO          SUN_RCONST(2.0)

// user data structure
typedef struct
{
  sunindextype nx; // global number of x grid points
  sunindextype ny; // global number of y grid points
  sunindextype is; // global x indices of this subdomain
  sunindextype ie;
  sunindextype js; // global y indices of this subdomain
  sunindextype je;
  sunindextype nxl;    // local number of x grid points
  sunindextype nyl;    // local number of y grid points
  sunrealtype dx;      // x-directional mesh spacing
  sunrealtype dy;      // y-directional mesh spacing
  sunrealtype kx;      // x-directional diffusion coefficient
  sunrealtype ky;      // y-directional diffusion coefficient
  N_Vector h;          // heat source vector
  N_Vector d;          // inverse of Jacobian diagonal
  MPI_Comm comm;       // communicator object
  int myid;            // MPI process ID
  int nprocs;          // total number of MPI processes
  bool HaveBdry[2][2]; // flags denoting if on physical boundary
  sunrealtype* Erecv;  // receive buffers for neighbor exchange
  sunrealtype* Wrecv;
  sunrealtype* Nrecv;
  sunrealtype* Srecv;
  sunrealtype* Esend; // send buffers for neighbor exchange
  sunrealtype* Wsend;
  sunrealtype* Nsend;
  sunrealtype* Ssend;
} UserData;

// User-supplied Functions Called by the Solver
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int PSet(sunrealtype t, N_Vector y, N_Vector fy, sunbooleantype jok,
                sunbooleantype* jcurPtr, sunrealtype gamma, void* user_data);
static int PSol(sunrealtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z,
                sunrealtype gamma, sunrealtype delta, int lr, void* user_data);

// Private functions
//    checks function return values
static int check_flag(void* flagvalue, const string funcname, int opt);
//    sets default values into UserData structure
static int InitUserData(UserData* udata);
//    sets up parallel decomposition
static int SetupDecomp(UserData* udata);
//    performs neighbor exchange
static int Exchange(N_Vector y, UserData* udata);
//    frees memory allocated within UserData
static int FreeUserData(UserData* udata);
//    check if relative difference is within tolerance
static bool Compare(long int a, long int b, sunrealtype tol);

// Main Program
int main(int argc, char* argv[])
{
  /* Create the SUNDIALS context object for this simulation. */
  SUNContext ctx = NULL;
  SUNContext_Create(SUN_COMM_NULL, &ctx);

  // general problem parameters
  sunrealtype T0   = SUN_RCONST(0.0); // initial time
  sunrealtype Tf   = SUN_RCONST(0.3); // final time
  int Nt           = 1000;            // total number of internal steps
  sunindextype nx  = 60;              // spatial mesh size
  sunindextype ny  = 120;
  sunrealtype kx   = SUN_RCONST(0.5); // heat conductivity coefficients
  sunrealtype ky   = SUN_RCONST(0.75);
  sunrealtype rtol = SUN_RCONST(1.e-5); // relative and absolute tolerances
  sunrealtype atol = SUN_RCONST(1.e-10);
  UserData* udata  = NULL;
  sunrealtype* data;
  sunindextype N, Ntot, i, j;
  int numfails;
  sunbooleantype linear;
  sunrealtype t;
  long int ark_nst, ark_nfe, ark_nfi, ark_nsetups, ark_nli, ark_nJv, ark_nlcf,
    ark_nni, ark_ncfn, ark_npe, ark_nps;
  long int mri_nst, mri_nfse, mri_nfsi, mri_nsetups, mri_nli, mri_nJv, mri_nlcf,
    mri_nni, mri_ncfn, mri_npe, mri_nps;

  // general problem variables
  int flag;                   // reusable error-checking flag
  int myid;                   // MPI process ID
  N_Vector y          = NULL; // empty vector for storing solution
  SUNLinearSolver LSa = NULL; // empty linear solver memory structures
  SUNLinearSolver LSm = NULL;
  void* arkstep_mem   = NULL; // empty ARKStep memory structure
  void* mristep_mem   = NULL; // empty MRIStep memory structure
  void* inner_mem     = NULL; // empty inner ARKStep memory structure

  // initialize MPI
  flag = MPI_Init(&argc, &argv);
  if (check_flag(&flag, "MPI_Init", 1)) { return 1; }
  flag = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (check_flag(&flag, "MPI_Comm_rank", 1)) { return 1; }

  // if an argument supplied, set linear (otherwise use SUNFALSE)
  linear = SUNFALSE;
  if (argc > 1) { linear = stoi(argv[1], NULL); }

  // allocate and fill udata structure
  udata = new UserData;
  flag  = InitUserData(udata);
  if (check_flag(&flag, "InitUserData", 1)) { return 1; }
  udata->nx = nx;
  udata->ny = ny;
  udata->kx = kx;
  udata->ky = ky;
  udata->dx = ONE / (ONE * nx - ONE); // x mesh spacing
  udata->dy = ONE / (ONE * ny - ONE); // y mesh spacing

  // Set up parallel decomposition
  flag = SetupDecomp(udata);
  if (check_flag(&flag, "SetupDecomp", 1)) { return 1; }

  // Initial problem output
  bool outproc = (udata->myid == 0);
  if (outproc)
  {
    cout << "\n2D Heat PDE test problem:\n";
    cout << "   nprocs = " << udata->nprocs << "\n";
    cout << "   nx = " << udata->nx << "\n";
    cout << "   ny = " << udata->ny << "\n";
    cout << "   kx = " << udata->kx << "\n";
    cout << "   ky = " << udata->ky << "\n";
    cout << "   rtol = " << rtol << "\n";
    cout << "   atol = " << atol << "\n";
    cout << "   nxl (proc 0) = " << udata->nxl << "\n";
    cout << "   nyl (proc 0) = " << udata->nyl << "\n";
    if (linear) { cout << "   Linearly implicit solver\n\n"; }
    else { cout << "   Nonlinear implicit solver\n\n"; }
  }

  // Initialize vector data structures
  N    = (udata->nxl) * (udata->nyl);
  Ntot = nx * ny;
  y    = N_VNew_Parallel(udata->comm, N, Ntot,
                         ctx); // Create parallel vector for solution
  if (check_flag((void*)y, "N_VNew_Parallel", 0)) { return 1; }
  N_VConst(ZERO, y);      // Set initial conditions
  udata->h = N_VClone(y); // Create vector for heat source
  if (check_flag((void*)udata->h, "N_VNew_Parallel", 0)) { return 1; }
  udata->d = N_VClone(y); // Create vector for Jacobian diagonal
  if (check_flag((void*)udata->d, "N_VNew_Parallel", 0)) { return 1; }

  // Initialize linear solver data structures
  LSa = SUNLinSol_PCG(y, 1, 20, ctx);
  if (check_flag((void*)LSa, "SUNLinSol_PCG", 0)) { return 1; }
  LSm = SUNLinSol_PCG(y, 1, 20, ctx);
  if (check_flag((void*)LSm, "SUNLinSol_PCG", 0)) { return 1; }

  // fill in the heat source array
  data = N_VGetArrayPointer(udata->h);
  for (j = 0; j < udata->nyl; j++)
  {
    for (i = 0; i < udata->nxl; i++)
    {
      data[IDX(i, j, udata->nxl)] = sin(PI * (udata->is + i) * udata->dx) *
                                    sin(TWO * PI * (udata->js + j) * udata->dy);
    }
  }

  /* Call ARKStepCreate and MRIStepCreate to initialize the timesteppers */
  arkstep_mem = ARKStepCreate(NULL, f, T0, y, ctx);
  if (check_flag((void*)arkstep_mem, "ARKStepCreate", 0)) { return 1; }

  inner_mem = ARKStepCreate(f0, NULL, T0, y, ctx);
  if (check_flag((void*)inner_mem, "ARKStepCreate", 0)) { return 1; }

  MRIStepInnerStepper inner_stepper = NULL;
  flag = ARKodeCreateMRIStepInnerStepper(inner_mem, &inner_stepper);
  if (check_flag(&flag, "ARKodeCreateMRIStepInnerStepper", 1)) { return 1; }

  mristep_mem = MRIStepCreate(NULL, f, T0, y, inner_stepper, ctx);
  if (check_flag((void*)mristep_mem, "MRIStepCreate", 0)) { return 1; }

  // Create solve-decoupled DIRK2 (trapezoidal) Butcher table
  ARKodeButcherTable B = ARKodeButcherTable_Alloc(2, SUNFALSE);
  if (check_flag((void*)B, "ARKodeButcherTable_Alloc", 0)) { return 1; }
  B->A[1][0] = SUN_RCONST(0.5);
  B->A[1][1] = SUN_RCONST(0.5);
  B->b[0]    = SUN_RCONST(0.5);
  B->b[1]    = SUN_RCONST(0.5);
  B->c[1]    = ONE;
  B->q       = 2;

  // Create solve-decoupled DIRK2 (trapezoidal) Butcher table
  ARKodeButcherTable Bc = ARKodeButcherTable_Alloc(3, SUNFALSE);
  if (check_flag((void*)Bc, "ARKodeButcherTable_Alloc", 0)) { return 1; }
  Bc->A[1][0] = ONE;
  Bc->A[2][0] = SUN_RCONST(0.5);
  Bc->A[2][2] = SUN_RCONST(0.5);
  Bc->b[0]    = SUN_RCONST(0.5);
  Bc->b[2]    = SUN_RCONST(0.5);
  Bc->c[1]    = ONE;
  Bc->c[2]    = ONE;
  Bc->q       = 2;

  // Create the MIS coupling table
  MRIStepCoupling C = MRIStepCoupling_MIStoMRI(Bc, 2, 0);
  if (check_flag((void*)C, "MRIStepCoupling_MIStoMRI", 0)) { return 1; }

  // Set routines
  flag = ARKodeSetUserData(arkstep_mem,
                           (void*)udata); // Pass udata to user functions
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }
  flag = ARKodeSetNonlinConvCoef(arkstep_mem,
                                 SUN_RCONST(1.e-7)); // Update solver convergence coeff.
  if (check_flag(&flag, "ARKodeSetNonlinConvCoef", 1)) { return 1; }
  flag = ARKodeSStolerances(arkstep_mem, rtol, atol); // Specify tolerances
  if (check_flag(&flag, "ARKodeSStolerances", 1)) { return 1; }
  flag = ARKodeSetFixedStep(arkstep_mem, Tf / Nt); // Specify fixed time step size
  if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }
  flag = ARKStepSetTables(arkstep_mem, 2, 0, B, NULL); // Specify Butcher table
  if (check_flag(&flag, "ARKStepSetTables", 1)) { return 1; }
  flag = ARKodeSetMaxNumSteps(arkstep_mem, 2 * Nt); // Increase num internal steps
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) { return 1; }

  flag = ARKodeSetUserData(mristep_mem,
                           (void*)udata); // Pass udata to user functions
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }
  flag = ARKodeSetNonlinConvCoef(mristep_mem,
                                 SUN_RCONST(1.e-7)); // Update solver convergence coeff.
  if (check_flag(&flag, "ARKodeSetNonlinConvCoef", 1)) { return 1; }
  flag = ARKodeSStolerances(mristep_mem, rtol, atol); // Specify tolerances
  if (check_flag(&flag, "ARKodeSStolerances", 1)) { return 1; }
  flag = ARKodeSetFixedStep(mristep_mem, Tf / Nt); // Specify fixed time step sizes
  if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }
  flag = ARKodeSetFixedStep(inner_mem, Tf / Nt / 10);
  if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }
  flag = MRIStepSetCoupling(mristep_mem, C); // Specify coupling table
  if (check_flag(&flag, "MRIStepSetCoupling", 1)) { return 1; }
  flag = ARKodeSetMaxNumSteps(mristep_mem, 2 * Nt); // Increase num internal steps
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) { return 1; }

  // Linear solver interface
  flag = ARKodeSetLinearSolver(arkstep_mem, LSa, NULL); // Attach linear solver
  if (check_flag(&flag, "ARKodeSetLinearSolver", 1)) { return 1; }
  flag = ARKodeSetPreconditioner(arkstep_mem, PSet,
                                 PSol); // Specify the Preconditioner
  if (check_flag(&flag, "ARKodeSetPreconditioner", 1)) { return 1; }

  flag = ARKodeSetLinearSolver(mristep_mem, LSm, NULL); // Attach linear solver
  if (check_flag(&flag, "ARKodeSetLinearSolver", 1)) { return 1; }
  flag = ARKodeSetPreconditioner(mristep_mem, PSet,
                                 PSol); // Specify the Preconditioner
  if (check_flag(&flag, "ARKodeSetPreconditioner", 1)) { return 1; }

  // Optionally specify linearly implicit RHS, with non-time-dependent preconditioner
  if (linear)
  {
    flag = ARKodeSetLinear(arkstep_mem, 0);
    if (check_flag(&flag, "ARKodeSetLinear", 1)) { return 1; }

    flag = ARKodeSetLinear(mristep_mem, 0);
    if (check_flag(&flag, "ARKodeSetLinear", 1)) { return 1; }
  }

  // First call ARKodeEvolve to evolve the full problem, and print results
  t = T0;
  N_VConst(ZERO, y);
  flag = ARKodeEvolve(arkstep_mem, Tf, y, &t, ARK_NORMAL);
  if (check_flag(&flag, "ARKodeEvolve", 1)) { return 1; }
  flag = ARKodeGetNumSteps(arkstep_mem, &ark_nst);
  if (check_flag(&flag, "ARKodeGetNumSteps", 1)) { return 1; }
  flag = ARKodeGetNumRhsEvals(arkstep_mem, 0, &ark_nfe);
  if (check_flag(&flag, "ARKodeGetNumRhsEvals", 1)) { return 1; }
  flag = ARKodeGetNumRhsEvals(arkstep_mem, 1, &ark_nfi);
  if (check_flag(&flag, "ARKodeGetNumRhsEvals", 1)) { return 1; }
  flag = ARKodeGetNumLinSolvSetups(arkstep_mem, &ark_nsetups);
  if (check_flag(&flag, "ARKodeGetNumLinSolvSetups", 1)) { return 1; }
  flag = ARKodeGetNumNonlinSolvIters(arkstep_mem, &ark_nni);
  if (check_flag(&flag, "ARKodeGetNumNonlinSolvIters", 1)) { return 1; }
  flag = ARKodeGetNumNonlinSolvConvFails(arkstep_mem, &ark_ncfn);
  if (check_flag(&flag, "ARKodeGetNumNonlinSolvConvFails", 1)) { return 1; }
  flag = ARKodeGetNumLinIters(arkstep_mem, &ark_nli);
  if (check_flag(&flag, "ARKodeGetNumLinIters", 1)) { return 1; }
  flag = ARKodeGetNumJtimesEvals(arkstep_mem, &ark_nJv);
  if (check_flag(&flag, "ARKodeGetNumJtimesEvals", 1)) { return 1; }
  flag = ARKodeGetNumLinConvFails(arkstep_mem, &ark_nlcf);
  if (check_flag(&flag, "ARKodeGetNumLinConvFails", 1)) { return 1; }
  flag = ARKodeGetNumPrecEvals(arkstep_mem, &ark_npe);
  if (check_flag(&flag, "ARKodeGetNumPrecEvals", 1)) { return 1; }
  flag = ARKodeGetNumPrecSolves(arkstep_mem, &ark_nps);
  if (check_flag(&flag, "ARKodeGetNumPrecSolves", 1)) { return 1; }
  if (outproc)
  {
    cout << "\nARKStep Solver Statistics:\n";
    cout << "   Internal solver steps = " << ark_nst << "\n";
    cout << "   Total RHS evals:  Fe = " << ark_nfe << ",  Fi = " << ark_nfi
         << "\n";
    cout << "   Total linear solver setups = " << ark_nsetups << "\n";
    cout << "   Total linear iterations = " << ark_nli << "\n";
    cout << "   Total number of Jacobian-vector products = " << ark_nJv << "\n";
    cout << "   Total number of Preconditioner setups = " << ark_npe << "\n";
    cout << "   Total number of Preconditioner solves = " << ark_nps << "\n";
    cout << "   Total number of linear solver convergence failures = " << ark_nlcf
         << "\n";
    cout << "   Total number of Newton iterations = " << ark_nni << "\n";
    cout << "   Total number of nonlinear solver convergence failures = "
         << ark_ncfn << "\n";
  }

  // Second call ARKodeEvolve to evolve the full problem, and print results
  t = T0;
  N_VConst(ZERO, y);
  flag = ARKodeEvolve(mristep_mem, Tf, y, &t, ARK_NORMAL);
  if (check_flag(&flag, "ARKodeEvolve", 1)) { return 1; }
  flag = ARKodeGetNumSteps(mristep_mem, &mri_nst);
  if (check_flag(&flag, "ARKodeGetNumSteps", 1)) { return 1; }
  flag = ARKodeGetNumRhsEvals(mristep_mem, 0, &mri_nfse);
  if (check_flag(&flag, "ARKodeGetNumRhsEvals", 1)) { return 1; }
  flag = ARKodeGetNumRhsEvals(mristep_mem, 1, &mri_nfsi);
  if (check_flag(&flag, "ARKodeGetNumRhsEvals", 1)) { return 1; }
  flag = ARKodeGetNumLinSolvSetups(mristep_mem, &mri_nsetups);
  if (check_flag(&flag, "ARKodeGetNumLinSolvSetups", 1)) { return 1; }
  flag = ARKodeGetNumNonlinSolvIters(mristep_mem, &mri_nni);
  if (check_flag(&flag, "ARKodeGetNumNonlinSolvIters", 1)) { return 1; }
  flag = ARKodeGetNumNonlinSolvConvFails(mristep_mem, &mri_ncfn);
  if (check_flag(&flag, "ARKodeGetNumNonlinSolvConvFails", 1)) { return 1; }
  flag = ARKodeGetNumLinIters(mristep_mem, &mri_nli);
  if (check_flag(&flag, "ARKodeGetNumLinIters", 1)) { return 1; }
  flag = ARKodeGetNumJtimesEvals(mristep_mem, &mri_nJv);
  if (check_flag(&flag, "ARKodeGetNumJtimesEvals", 1)) { return 1; }
  flag = ARKodeGetNumLinConvFails(mristep_mem, &mri_nlcf);
  if (check_flag(&flag, "ARKodeGetNumLinConvFails", 1)) { return 1; }
  flag = ARKodeGetNumPrecEvals(mristep_mem, &mri_npe);
  if (check_flag(&flag, "ARKodeGetNumPrecEvals", 1)) { return 1; }
  flag = ARKodeGetNumPrecSolves(mristep_mem, &mri_nps);
  if (check_flag(&flag, "ARKodeGetNumPrecSolves", 1)) { return 1; }
  if (outproc)
  {
    cout << "\nMRIStep Solver Statistics:\n";
    cout << "   Internal solver steps = " << mri_nst << "\n";
    cout << "   Total RHS evals:  Fe = " << mri_nfsi << "\n";
    cout << "   Total linear solver setups = " << mri_nsetups << "\n";
    cout << "   Total linear iterations = " << mri_nli << "\n";
    cout << "   Total number of Jacobian-vector products = " << mri_nJv << "\n";
    cout << "   Total number of Preconditioner setups = " << mri_npe << "\n";
    cout << "   Total number of Preconditioner solves = " << mri_nps << "\n";
    cout << "   Total number of linear solver convergence failures = " << mri_nlcf
         << "\n";
    cout << "   Total number of Newton iterations = " << mri_nni << "\n";
    cout << "   Total number of nonlinear solver convergence failures = "
         << mri_ncfn << "\n";
  }

  // Compare solver statistics
  numfails = 0;
  if (outproc) { cout << "\nComparing Solver Statistics:\n"; }
  if (ark_nst != mri_nst)
  {
    numfails += 1;
    if (outproc)
    {
      cout << "  Internal solver steps error: " << ark_nst << " vs " << mri_nst
           << "\n";
    }
  }
  if (ark_nfi != mri_nfsi)
  {
    numfails += 1;
    if (outproc)
    {
      cout << "  RHS evals error: " << ark_nfi << " vs " << mri_nfsi << "\n";
    }
  }
  if (ark_nsetups != mri_nsetups)
  {
    numfails += 1;
    if (outproc)
    {
      cout << "  Linear solver setups error: " << ark_nsetups << " vs "
           << mri_nsetups << "\n";
    }
  }
  if (!Compare(ark_nli, mri_nli, ONE))
  {
    numfails += 1;
    if (outproc)
    {
      cout << "  Linear iterations error: " << ark_nli << " vs " << mri_nli
           << "\n";
    }
  }
  if (!Compare(ark_nJv, mri_nJv, ONE))
  {
    numfails += 1;
    if (outproc)
    {
      cout << "  Jacobian-vector products error: " << ark_nJv << " vs "
           << mri_nJv << "\n";
    }
  }
  if (!Compare(ark_nps, mri_nps, ONE))
  {
    numfails += 1;
    if (outproc)
    {
      cout << "  Preconditioner solves error: " << ark_nps << " vs " << mri_nps
           << "\n";
    }
  }
  if (ark_nlcf != mri_nlcf)
  {
    numfails += 1;
    if (outproc)
    {
      cout << "  Linear convergence failures error: " << ark_nlcf << " vs "
           << mri_nlcf << "\n";
    }
  }
  if (ark_nni != mri_nni)
  {
    numfails += 1;
    if (outproc)
    {
      cout << "  Newton iterations error: " << ark_nni << " vs " << mri_nni
           << "\n";
    }
  }
  if (ark_ncfn != mri_ncfn)
  {
    numfails += 1;
    if (outproc)
    {
      cout << "  Nonlinear convergence failures error: " << ark_ncfn << " vs "
           << mri_ncfn << "\n";
    }
  }
  if (outproc)
  {
    if (numfails) { cout << "Failed " << numfails << " tests\n"; }
    else { cout << "All tests pass!\n"; }
  }

  // Clean up and return with successful completion
  ARKodeButcherTable_Free(B);  // Free Butcher table
  ARKodeButcherTable_Free(Bc); // Free Butcher table
  MRIStepCoupling_Free(C);     // Free MRI coupling table
  ARKodeFree(&arkstep_mem);    // Free integrator memory
  ARKodeFree(&mristep_mem);
  ARKodeFree(&inner_mem);
  MRIStepInnerStepper_Free(&inner_stepper);
  SUNLinSolFree(LSa); // Free linear solver
  SUNLinSolFree(LSm);
  N_VDestroy(y); // Free vectors
  N_VDestroy(udata->h);
  N_VDestroy(udata->d);
  FreeUserData(udata); // Free user data
  delete udata;

  SUNContext_Free(&ctx);
  flag = MPI_Finalize(); // Finalize MPI
  return (numfails);
}

/*--------------------------------
 * Functions called by the solver
 *--------------------------------*/

// f routine to compute the ODE RHS function f(t,y).
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VConst(ZERO, ydot);                    // Initialize ydot to zero
  UserData* udata  = (UserData*)user_data; // access problem data
  sunindextype nxl = udata->nxl;           // set variable shortcuts
  sunindextype nyl = udata->nyl;
  sunrealtype kx   = udata->kx;
  sunrealtype ky   = udata->ky;
  sunrealtype dx   = udata->dx;
  sunrealtype dy   = udata->dy;
  sunrealtype* Y   = N_VGetArrayPointer(y); // access data arrays
  if (check_flag((void*)Y, "N_VGetArrayPointer", 0)) { return -1; }
  sunrealtype* Ydot = N_VGetArrayPointer(ydot);
  if (check_flag((void*)Ydot, "N_VGetArrayPointer", 0)) { return -1; }

  // Exchange boundary data with neighbors
  int ierr = Exchange(y, udata);
  if (check_flag(&ierr, "Exchange", 1)) { return -1; }

  // iterate over subdomain interior, computing approximation to RHS
  sunrealtype c1 = kx / dx / dx;
  sunrealtype c2 = ky / dy / dy;
  sunrealtype c3 = -TWO * (c1 + c2);
  sunindextype i, j;
  for (j = 1; j < nyl - 1; j++)
  { // diffusive terms
    for (i = 1; i < nxl - 1; i++)
    {
      Ydot[IDX(i, j, nxl)] =
        c1 * (Y[IDX(i - 1, j, nxl)] + Y[IDX(i + 1, j, nxl)]) +
        c2 * (Y[IDX(i, j - 1, nxl)] + Y[IDX(i, j + 1, nxl)]) +
        c3 * Y[IDX(i, j, nxl)];
    }
  }

  // iterate over subdomain boundaries (if not at overall domain boundary)
  if (!udata->HaveBdry[0][0])
  { // West face
    i = 0;
    for (j = 1; j < nyl - 1; j++)
    {
      Ydot[IDX(i, j, nxl)] = c1 * (udata->Wrecv[j] + Y[IDX(i + 1, j, nxl)]) +
                             c2 * (Y[IDX(i, j - 1, nxl)] + Y[IDX(i, j + 1, nxl)]) +
                             c3 * Y[IDX(i, j, nxl)];
    }
  }
  if (!udata->HaveBdry[0][1])
  { // East face
    i = nxl - 1;
    for (j = 1; j < nyl - 1; j++)
    {
      Ydot[IDX(i, j, nxl)] = c1 * (Y[IDX(i - 1, j, nxl)] + udata->Erecv[j]) +
                             c2 * (Y[IDX(i, j - 1, nxl)] + Y[IDX(i, j + 1, nxl)]) +
                             c3 * Y[IDX(i, j, nxl)];
    }
  }
  if (!udata->HaveBdry[1][0])
  { // South face
    j = 0;
    for (i = 1; i < nxl - 1; i++)
    {
      Ydot[IDX(i, j, nxl)] =
        c1 * (Y[IDX(i - 1, j, nxl)] + Y[IDX(i + 1, j, nxl)]) +
        c2 * (udata->Srecv[i] + Y[IDX(i, j + 1, nxl)]) + c3 * Y[IDX(i, j, nxl)];
    }
  }
  if (!udata->HaveBdry[1][1])
  { // West face
    j = nyl - 1;
    for (i = 1; i < nxl - 1; i++)
    {
      Ydot[IDX(i, j, nxl)] =
        c1 * (Y[IDX(i - 1, j, nxl)] + Y[IDX(i + 1, j, nxl)]) +
        c2 * (Y[IDX(i, j - 1, nxl)] + udata->Nrecv[i]) + c3 * Y[IDX(i, j, nxl)];
    }
  }
  if (!udata->HaveBdry[0][0] && !udata->HaveBdry[1][0])
  { // South-West corner
    i                    = 0;
    j                    = 0;
    Ydot[IDX(i, j, nxl)] = c1 * (udata->Wrecv[j] + Y[IDX(i + 1, j, nxl)]) +
                           c2 * (udata->Srecv[i] + Y[IDX(i, j + 1, nxl)]) +
                           c3 * Y[IDX(i, j, nxl)];
  }
  if (!udata->HaveBdry[0][0] && !udata->HaveBdry[1][1])
  { // North-West corner
    i                    = 0;
    j                    = nyl - 1;
    Ydot[IDX(i, j, nxl)] = c1 * (udata->Wrecv[j] + Y[IDX(i + 1, j, nxl)]) +
                           c2 * (Y[IDX(i, j - 1, nxl)] + udata->Nrecv[i]) +
                           c3 * Y[IDX(i, j, nxl)];
  }
  if (!udata->HaveBdry[0][1] && !udata->HaveBdry[1][0])
  { // South-East corner
    i                    = nxl - 1;
    j                    = 0;
    Ydot[IDX(i, j, nxl)] = c1 * (Y[IDX(i - 1, j, nxl)] + udata->Erecv[j]) +
                           c2 * (udata->Srecv[i] + Y[IDX(i, j + 1, nxl)]) +
                           c3 * Y[IDX(i, j, nxl)];
  }
  if (!udata->HaveBdry[0][1] && !udata->HaveBdry[1][1])
  { // North-East corner
    i                    = nxl - 1;
    j                    = nyl - 1;
    Ydot[IDX(i, j, nxl)] = c1 * (Y[IDX(i - 1, j, nxl)] + udata->Erecv[j]) +
                           c2 * (Y[IDX(i, j - 1, nxl)] + udata->Nrecv[i]) +
                           c3 * Y[IDX(i, j, nxl)];
  }

  // add in heat source
  N_VLinearSum(ONE, ydot, ONE, udata->h, ydot);
  return 0; // Return with success
}

// f0 routine to compute a zero-valued ODE RHS function f(t,y).
static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  // Initialize ydot to zero and return
  N_VConst(ZERO, ydot);
  return 0;
}

// Preconditioner setup routine (fills inverse of Jacobian diagonal)
static int PSet(sunrealtype t, N_Vector y, N_Vector fy, sunbooleantype jok,
                sunbooleantype* jcurPtr, sunrealtype gamma, void* user_data)
{
  UserData* udata   = (UserData*)user_data; // variable shortcuts
  sunrealtype kx    = udata->kx;
  sunrealtype ky    = udata->ky;
  sunrealtype dx    = udata->dx;
  sunrealtype dy    = udata->dy;
  sunrealtype* diag = N_VGetArrayPointer(udata->d); // access data arrays
  if (check_flag((void*)diag, "N_VGetArrayPointer", 0)) { return -1; }

  // set all entries of d to the diagonal values of interior
  // (since boundary RHS is 0, set boundary diagonals to the same)
  sunrealtype c = ONE + gamma * TWO * (kx / dx / dx + ky / dy / dy);
  N_VConst(c, udata->d);
  N_VInv(udata->d, udata->d); // invert diagonal
  return 0;                   // Return with success
}

// Preconditioner solve routine
static int PSol(sunrealtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z,
                sunrealtype gamma, sunrealtype delta, int lr, void* user_data)
{
  UserData* udata = (UserData*)user_data; // access user_data structure
  N_VProd(r, udata->d, z);                // perform Jacobi iteration
  return 0;                               // Return with success
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_flag(void* flagvalue, const string funcname, int opt)
{
  int* errflag;

  // Check if SUNDIALS function returned NULL pointer - no memory allocated
  if (opt == 0 && flagvalue == NULL)
  {
    cerr << "\nSUNDIALS_ERROR: " << funcname
         << " failed - returned NULL pointer\n\n";
    return 1;
  }

  // Check if flag < 0
  else if (opt == 1)
  {
    errflag = (int*)flagvalue;
    if (*errflag < 0)
    {
      cerr << "\nSUNDIALS_ERROR: " << funcname
           << " failed with flag = " << *errflag << "\n\n";
      return 1;
    }
  }

  // Check if function returned NULL pointer - no memory allocated
  else if (opt == 2 && flagvalue == NULL)
  {
    cerr << "\nMEMORY_ERROR: " << funcname
         << " failed - returned NULL pointer\n\n";
    return 1;
  }

  return 0;
}

// Set up parallel decomposition
static int SetupDecomp(UserData* udata)
{
  // check that this has not been called before
  if (udata->Erecv != NULL || udata->Wrecv != NULL || udata->Srecv != NULL ||
      udata->Nrecv != NULL)
  {
    cerr << "SetupDecomp warning: parallel decomposition already set up\n";
    return 1;
  }

  // get suggested parallel decomposition
  int ierr, dims[] = {0, 0};
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &(udata->nprocs));
  if (ierr != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_size = " << ierr << "\n";
    return -1;
  }
  ierr = MPI_Dims_create(udata->nprocs, 2, dims);
  if (ierr != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Dims_create = " << ierr << "\n";
    return -1;
  }

  // set up 2D Cartesian communicator
  int periods[] = {0, 0};
  ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &(udata->comm));
  if (ierr != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_create = " << ierr << "\n";
    return -1;
  }
  ierr = MPI_Comm_rank(udata->comm, &(udata->myid));
  if (ierr != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_rank = " << ierr << "\n";
    return -1;
  }

  // determine local extents
  int coords[2];
  ierr = MPI_Cart_get(udata->comm, 2, dims, periods, coords);
  if (ierr != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_get = " << ierr << "\n";
    return -1;
  }
  udata->is  = (udata->nx) * (coords[0]) / (dims[0]);
  udata->ie  = (udata->nx) * (coords[0] + 1) / (dims[0]) - 1;
  udata->js  = (udata->ny) * (coords[1]) / (dims[1]);
  udata->je  = (udata->ny) * (coords[1] + 1) / (dims[1]) - 1;
  udata->nxl = (udata->ie) - (udata->is) + 1;
  udata->nyl = (udata->je) - (udata->js) + 1;

  // determine if I have neighbors, and allocate exchange buffers
  udata->HaveBdry[0][0] = (udata->is == 0);
  udata->HaveBdry[0][1] = (udata->ie == udata->nx - 1);
  udata->HaveBdry[1][0] = (udata->js == 0);
  udata->HaveBdry[1][1] = (udata->je == udata->ny - 1);
  if (!udata->HaveBdry[0][0])
  {
    udata->Wrecv = new sunrealtype[udata->nyl];
    udata->Wsend = new sunrealtype[udata->nyl];
  }
  if (!udata->HaveBdry[0][1])
  {
    udata->Erecv = new sunrealtype[udata->nyl];
    udata->Esend = new sunrealtype[udata->nyl];
  }
  if (!udata->HaveBdry[1][0])
  {
    udata->Srecv = new sunrealtype[udata->nxl];
    udata->Ssend = new sunrealtype[udata->nxl];
  }
  if (!udata->HaveBdry[1][1])
  {
    udata->Nrecv = new sunrealtype[udata->nxl];
    udata->Nsend = new sunrealtype[udata->nxl];
  }

  return 0; // return with success flag
}

// Perform neighbor exchange
static int Exchange(N_Vector y, UserData* udata)
{
  // local variables
  MPI_Request reqSW, reqSE, reqSS, reqSN, reqRW, reqRE, reqRS, reqRN;
  MPI_Status stat;
  int ierr, i, ipW = -1, ipE = -1, ipS = -1, ipN = -1;
  int coords[2], dims[2], periods[2], nbcoords[2];
  sunindextype nyl = udata->nyl;
  sunindextype nxl = udata->nxl;

  // access data array
  sunrealtype* Y = N_VGetArrayPointer(y);
  if (check_flag((void*)Y, "N_VGetArrayPointer", 0)) { return -1; }

  // MPI neighborhood information
  ierr = MPI_Cart_get(udata->comm, 2, dims, periods, coords);
  if (ierr != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_get = " << ierr << "\n";
    return -1;
  }
  if (!udata->HaveBdry[0][0])
  {
    nbcoords[0] = coords[0] - 1;
    nbcoords[1] = coords[1];
    ierr        = MPI_Cart_rank(udata->comm, nbcoords, &ipW);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[0][1])
  {
    nbcoords[0] = coords[0] + 1;
    nbcoords[1] = coords[1];
    ierr        = MPI_Cart_rank(udata->comm, nbcoords, &ipE);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][0])
  {
    nbcoords[0] = coords[0];
    nbcoords[1] = coords[1] - 1;
    ierr        = MPI_Cart_rank(udata->comm, nbcoords, &ipS);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][1])
  {
    nbcoords[0] = coords[0];
    nbcoords[1] = coords[1] + 1;
    ierr        = MPI_Cart_rank(udata->comm, nbcoords, &ipN);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << ierr << "\n";
      return -1;
    }
  }

  // open Irecv buffers
  if (!udata->HaveBdry[0][0])
  {
    ierr = MPI_Irecv(udata->Wrecv, (int)udata->nyl, MPI_SUNREALTYPE, ipW,
                     MPI_ANY_TAG, udata->comm, &reqRW);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[0][1])
  {
    ierr = MPI_Irecv(udata->Erecv, (int)udata->nyl, MPI_SUNREALTYPE, ipE,
                     MPI_ANY_TAG, udata->comm, &reqRE);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][0])
  {
    ierr = MPI_Irecv(udata->Srecv, (int)udata->nxl, MPI_SUNREALTYPE, ipS,
                     MPI_ANY_TAG, udata->comm, &reqRS);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][1])
  {
    ierr = MPI_Irecv(udata->Nrecv, (int)udata->nxl, MPI_SUNREALTYPE, ipN,
                     MPI_ANY_TAG, udata->comm, &reqRN);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << ierr << "\n";
      return -1;
    }
  }

  // send data
  if (!udata->HaveBdry[0][0])
  {
    for (i = 0; i < nyl; i++) { udata->Wsend[i] = Y[IDX(0, i, nxl)]; }
    ierr = MPI_Isend(udata->Wsend, (int)udata->nyl, MPI_SUNREALTYPE, ipW, 0,
                     udata->comm, &reqSW);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[0][1])
  {
    for (i = 0; i < nyl; i++) { udata->Esend[i] = Y[IDX(nxl - 1, i, nxl)]; }
    ierr = MPI_Isend(udata->Esend, (int)udata->nyl, MPI_SUNREALTYPE, ipE, 1,
                     udata->comm, &reqSE);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][0])
  {
    for (i = 0; i < nxl; i++) { udata->Ssend[i] = Y[IDX(i, 0, nxl)]; }
    ierr = MPI_Isend(udata->Ssend, (int)udata->nxl, MPI_SUNREALTYPE, ipS, 2,
                     udata->comm, &reqSS);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][1])
  {
    for (i = 0; i < nxl; i++) { udata->Nsend[i] = Y[IDX(i, nyl - 1, nxl)]; }
    ierr = MPI_Isend(udata->Nsend, (int)udata->nxl, MPI_SUNREALTYPE, ipN, 3,
                     udata->comm, &reqSN);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << ierr << "\n";
      return -1;
    }
  }

  // wait for messages to finish
  if (!udata->HaveBdry[0][0])
  {
    ierr = MPI_Wait(&reqRW, &stat);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
    ierr = MPI_Wait(&reqSW, &stat);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[0][1])
  {
    ierr = MPI_Wait(&reqRE, &stat);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
    ierr = MPI_Wait(&reqSE, &stat);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][0])
  {
    ierr = MPI_Wait(&reqRS, &stat);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
    ierr = MPI_Wait(&reqSS, &stat);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
  }
  if (!udata->HaveBdry[1][1])
  {
    ierr = MPI_Wait(&reqRN, &stat);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
    ierr = MPI_Wait(&reqSN, &stat);
    if (ierr != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << ierr << "\n";
      return -1;
    }
  }

  return 0; // return with success flag
}

// Initialize memory allocated within Userdata
static int InitUserData(UserData* udata)
{
  udata->nx             = 0;
  udata->ny             = 0;
  udata->is             = 0;
  udata->ie             = 0;
  udata->js             = 0;
  udata->je             = 0;
  udata->nxl            = 0;
  udata->nyl            = 0;
  udata->dx             = ZERO;
  udata->dy             = ZERO;
  udata->kx             = ZERO;
  udata->ky             = ZERO;
  udata->h              = NULL;
  udata->d              = NULL;
  udata->comm           = MPI_COMM_WORLD;
  udata->myid           = 0;
  udata->nprocs         = 0;
  udata->HaveBdry[0][0] = 1;
  udata->HaveBdry[0][1] = 1;
  udata->HaveBdry[1][0] = 1;
  udata->HaveBdry[1][1] = 1;
  udata->Erecv          = NULL;
  udata->Wrecv          = NULL;
  udata->Nrecv          = NULL;
  udata->Srecv          = NULL;
  udata->Esend          = NULL;
  udata->Wsend          = NULL;
  udata->Nsend          = NULL;
  udata->Ssend          = NULL;

  return 0; // return with success flag
}

// Free memory allocated within Userdata
static int FreeUserData(UserData* udata)
{
  // free exchange buffers
  if (udata->Wrecv != NULL) { delete[] udata->Wrecv; }
  if (udata->Wsend != NULL) { delete[] udata->Wsend; }
  if (udata->Erecv != NULL) { delete[] udata->Erecv; }
  if (udata->Esend != NULL) { delete[] udata->Esend; }
  if (udata->Srecv != NULL) { delete[] udata->Srecv; }
  if (udata->Ssend != NULL) { delete[] udata->Ssend; }
  if (udata->Nrecv != NULL) { delete[] udata->Nrecv; }
  if (udata->Nsend != NULL) { delete[] udata->Nsend; }

  if (udata->comm != MPI_COMM_WORLD) { MPI_Comm_free(&(udata->comm)); }

  return 0; // return with success flag
}

// Check if relative difference of a and b is less than tolerance
static bool Compare(long int a, long int b, sunrealtype tol)
{
  sunrealtype rel_diff = SUN_RCONST(100.0) * abs(static_cast<sunrealtype>(a - b) /
                                                 static_cast<sunrealtype>(a));

  return (rel_diff > tol) ? false : true;
}

//---- end of file ----
