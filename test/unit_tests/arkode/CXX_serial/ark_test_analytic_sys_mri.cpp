/*-----------------------------------------------------------------
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
 * solver and the fixed-point solver based on a 0/1 command-line
 * argument (1 => fixed-point).
 *-----------------------------------------------------------------*/

// Header files
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_mristep.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <nvector/nvector_serial.h>
#include <string>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

using namespace std;

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

// User-supplied Functions Called by the Solver
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Private function to perform matrix-matrix product
static int dense_MM(SUNMatrix A, SUNMatrix B, SUNMatrix C);

// Private function to check function return values
static int check_flag(void* flagvalue, const string funcname, int opt);

//    check if relative difference is within tolerance
static bool Compare(long int a, long int b, sunrealtype tol);

// SUNContext for the simulation
static SUNContext sunctx = NULL;

// Main Program
int main(int argc, char* argv[])
{
  // Create the SUNDIALS context object for this simulation.
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  // general problem parameters
  sunrealtype T0     = SUN_RCONST(0.0);    // initial time
  sunrealtype Tf     = SUN_RCONST(0.05);   // final time
  int Nt             = 1000;               // total number of internal steps
  sunindextype NEQ   = 3;                  // number of dependent vars.
  sunrealtype reltol = SUN_RCONST(1.0e-6); // tolerances
  sunrealtype abstol = SUN_RCONST(1.0e-10);
  sunrealtype lambda = SUN_RCONST(-100.0); // stiffness parameter

  // general problem variables
  int flag;                       // reusable error-checking flag
  N_Vector y              = NULL; // empty vector for storing solution
  SUNMatrix Aa            = NULL; // empty dense matrices for solvers
  SUNMatrix Am            = NULL;
  SUNNonlinearSolver NLSa = NULL; // empty nonlinear solvers
  SUNNonlinearSolver NLSm = NULL;
  SUNLinearSolver LSa     = NULL; // empty dense linear solvers
  SUNLinearSolver LSm     = NULL;
  void* arkstep_mem       = NULL; // empty ARKStep memory structure
  void* mristep_mem       = NULL; // empty MRIStep memory structure
  void* inner_mem         = NULL; // empty inner ARKStep memory structure
  int numfails;
  sunbooleantype fixedpoint;
  sunrealtype t, tcur;
  long int ark_nst, ark_nfe, ark_nfi, ark_nsetups, ark_nje, ark_nfeLS, ark_nni,
    ark_ncfn;
  long int mri_nst, mri_nfse, mri_nfsi, mri_nsetups, mri_nje, mri_nfeLS,
    mri_nni, mri_ncfn;

  // if an argument supplied, set fixedpoint (otherwise use SUNFALSE)
  fixedpoint = SUNFALSE;
  if (argc > 1) { fixedpoint = std::stoi(argv[1], NULL); }
  if (argc > 2) { Nt = std::stoi(argv[2], NULL); }

  // Initial problem output
  cout << "\nAnalytical ODE test problem:\n";
  cout << "   lambda = " << lambda << "\n";
  cout << "   reltol = " << reltol << "\n";
  cout << "   abstol = " << abstol << "\n\n";
  if (fixedpoint) { cout << "   Fixed-point nonlinear solver\n"; }
  else { cout << "   Newton nonlinear solver\n"; }

  // Initialize vector data structure and specify initial condition
  y = N_VNew_Serial(NEQ, sunctx);
  if (check_flag((void*)y, "N_VNew_Serial", 0)) { return 1; }
  N_VConst(ONE, y);

  /* Call ARKStepCreate and MRIStepCreate to initialize the timesteppers */
  arkstep_mem = ARKStepCreate(NULL, f, T0, y, sunctx);
  if (check_flag((void*)arkstep_mem, "ARKStepCreate", 0)) { return 1; }

  inner_mem = ARKStepCreate(f0, NULL, T0, y, sunctx);
  if (check_flag((void*)inner_mem, "ARKStepCreate", 0)) { return 1; }

  MRIStepInnerStepper inner_stepper = NULL;
  flag = ARKodeCreateMRIStepInnerStepper(inner_mem, &inner_stepper);
  if (check_flag(&flag, "ARKodeCreateMRIStepInnerStepper", 1)) { return 1; }

  mristep_mem = MRIStepCreate(NULL, f, T0, y, inner_stepper, sunctx);
  if (check_flag((void*)mristep_mem, "MRIStepCreate", 0)) { return 1; }

  // Create DIRK2 (trapezoidal) Butcher table
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
                           (void*)&lambda); // Pass lambda to user functions
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }
  flag = ARKodeSStolerances(arkstep_mem, reltol, abstol); // Specify tolerances
  if (check_flag(&flag, "ARKodeSStolerances", 1)) { return 1; }
  flag = ARKodeSetFixedStep(arkstep_mem, Tf / Nt); // Specify fixed time step size
  if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }
  flag = ARKStepSetTables(arkstep_mem, 2, 0, B, NULL); // Specify Butcher table
  if (check_flag(&flag, "ARKStepSetTables", 1)) { return 1; }
  flag = ARKodeSetMaxNumSteps(arkstep_mem, 2 * Nt); // Increase num internal steps
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) { return 1; }

  flag = ARKodeSetUserData(mristep_mem,
                           (void*)&lambda); // Pass lambda to user functions
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }
  flag = ARKodeSStolerances(mristep_mem, reltol, abstol); // Specify tolerances
  if (check_flag(&flag, "ARKodeSStolerances", 1)) { return 1; }
  flag = ARKodeSetFixedStep(mristep_mem, Tf / Nt); // Specify fixed time step sizes
  if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }
  flag = ARKodeSetFixedStep(inner_mem, Tf / Nt / 10);
  if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }
  flag = MRIStepSetCoupling(mristep_mem, C); // Specify coupling table
  if (check_flag(&flag, "MRIStepSetCoupling", 1)) { return 1; }
  flag = ARKodeSetMaxNumSteps(mristep_mem, 2 * Nt); // Increase num internal steps
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) { return 1; }

  // Initialize implicit solver data structures
  if (fixedpoint)
  {
    // Initialize fixed-point solvers and attach to integrators
    NLSa = SUNNonlinSol_FixedPoint(y, 50, sunctx);
    if (check_flag((void*)NLSa, "SUNNonlinSol_FixedPoint", 0)) { return 1; }
    flag = ARKodeSetNonlinearSolver(arkstep_mem, NLSa);
    if (check_flag(&flag, "ARKodeSetNonlinearSolver", 1)) { return 1; }

    NLSm = SUNNonlinSol_FixedPoint(y, 50, sunctx);
    if (check_flag((void*)NLSm, "SUNNonlinSol_FixedPoint", 0)) { return 1; }
    flag = ARKodeSetNonlinearSolver(mristep_mem, NLSm);
    if (check_flag(&flag, "ARKodeSetNonlinearSolver", 1)) { return 1; }
  }
  else
  {
    // Initialize dense matrix data structures and solvers
    Aa = SUNDenseMatrix(NEQ, NEQ, sunctx);
    if (check_flag((void*)Aa, "SUNDenseMatrix", 0)) { return 1; }
    LSa = SUNLinSol_Dense(y, Aa, sunctx);
    if (check_flag((void*)LSa, "SUNLinSol_Dense", 0)) { return 1; }

    Am = SUNDenseMatrix(NEQ, NEQ, sunctx);
    if (check_flag((void*)Am, "SUNDenseMatrix", 0)) { return 1; }
    LSm = SUNLinSol_Dense(y, Am, sunctx);
    if (check_flag((void*)LSm, "SUNLinSol_Dense", 0)) { return 1; }

    // Linear solver interface
    flag = ARKodeSetLinearSolver(arkstep_mem, LSa, Aa);
    if (check_flag(&flag, "ARKodeSetLinearSolver", 1)) { return 1; }
    flag = ARKodeSetJacFn(arkstep_mem, Jac);
    if (check_flag(&flag, "ARKodeSetJacFn", 1)) { return 1; }

    flag = ARKodeSetLinearSolver(mristep_mem, LSm, Am);
    if (check_flag(&flag, "ARKodeSetLinearSolver", 1)) { return 1; }
    flag = ARKodeSetJacFn(mristep_mem, Jac);
    if (check_flag(&flag, "ARKodeSetJacFn", 1)) { return 1; }

    // Specify linearly implicit RHS, with non-time-dependent Jacobian
    flag = ARKodeSetLinear(arkstep_mem, 0);
    if (check_flag(&flag, "ARKodeSetLinear", 1)) { return 1; }

    flag = ARKodeSetLinear(mristep_mem, 0);
    if (check_flag(&flag, "ARKodeSetLinear", 1)) { return 1; }
  }

  // First call ARKodeEvolve to evolve the full problem, and print results
  t = T0;
  N_VConst(ONE, y);
  flag = ARKodeEvolve(arkstep_mem, Tf, y, &t, ARK_NORMAL);
  if (check_flag(&flag, "ARKodeEvolve", 1)) { return 1; }
  flag = ARKodeGetCurrentTime(arkstep_mem, &tcur);
  if (check_flag(&flag, "ARKodeGetCurrentTime", 1)) { return 1; }
  flag = ARKodeGetNumSteps(arkstep_mem, &ark_nst);
  if (check_flag(&flag, "ARKodeGetNumSteps", 1)) { return 1; }
  flag = ARKodeGetNumRhsEvals(arkstep_mem, 0, &ark_nfe);
  if (check_flag(&flag, "ARKodeGetNumRhsEvals", 1)) { return 1; }
  flag = ARKodeGetNumRhsEvals(arkstep_mem, 1, &ark_nfi);
  if (check_flag(&flag, "ARKodeGetNumRhsEvals", 1)) { return 1; }
  flag = ARKodeGetNumNonlinSolvIters(arkstep_mem, &ark_nni);
  if (check_flag(&flag, "ARKodeGetNumNonlinSolvIters", 1)) { return 1; }
  flag = ARKodeGetNumNonlinSolvConvFails(arkstep_mem, &ark_ncfn);
  if (check_flag(&flag, "ARKodeGetNumNonlinSolvConvFails", 1)) { return 1; }
  if (!fixedpoint)
  {
    flag = ARKodeGetNumLinSolvSetups(arkstep_mem, &ark_nsetups);
    if (check_flag(&flag, "ARKodeGetNumLinSolvSetups", 1)) { return 1; }
    flag = ARKodeGetNumJacEvals(arkstep_mem, &ark_nje);
    if (check_flag(&flag, "ARKodeGetNumJacEvals", 1)) { return 1; }
    flag = ARKodeGetNumLinRhsEvals(arkstep_mem, &ark_nfeLS);
    check_flag(&flag, "ARKodeGetNumLinRhsEvals", 1);
  }
  cout << "\nARKStep Solver Statistics:\n";
  cout << "   Return time = " << t << "\n";
  cout << "   Internal final time = " << tcur << "\n";
  cout << "   Internal solver steps = " << ark_nst << "\n";
  cout << "   Total RHS evals:  Fe = " << ark_nfe << ",  Fi = " << ark_nfi
       << "\n";
  cout << "   Total number of nonlinear iterations = " << ark_nni << "\n";
  cout << "   Total number of nonlinear solver convergence failures = "
       << ark_ncfn << "\n";
  if (!fixedpoint)
  {
    cout << "   Total linear solver setups = " << ark_nsetups << "\n";
    cout << "   Total RHS evals for setting up the linear system = " << ark_nfeLS
         << "\n";
    cout << "   Total number of Jacobian evaluations = " << ark_nje << "\n";
  }

  // Second call ARKodeEvolve to evolve the full problem, and print results
  t = T0;
  N_VConst(ZERO, y);
  flag = ARKodeEvolve(mristep_mem, Tf, y, &t, ARK_NORMAL);
  if (check_flag(&flag, "ARKodeEvolve", 1)) { return 1; }
  flag = ARKodeGetCurrentTime(arkstep_mem, &tcur);
  if (check_flag(&flag, "ARKodeGetCurrentTime", 1)) { return 1; }
  flag = ARKodeGetNumSteps(mristep_mem, &mri_nst);
  if (check_flag(&flag, "ARKodeGetNumSteps", 1)) { return 1; }
  flag = ARKodeGetNumRhsEvals(mristep_mem, 0, &mri_nfse);
  if (check_flag(&flag, "ARKodeGetNumRhsEvals", 1)) { return 1; }
  flag = ARKodeGetNumRhsEvals(mristep_mem, 1, &mri_nfsi);
  if (check_flag(&flag, "ARKodeGetNumRhsEvals", 1)) { return 1; }
  flag = ARKodeGetNumNonlinSolvIters(mristep_mem, &mri_nni);
  if (check_flag(&flag, "ARKodeGetNumNonlinSolvIters", 1)) { return 1; }
  flag = ARKodeGetNumNonlinSolvConvFails(mristep_mem, &mri_ncfn);
  if (check_flag(&flag, "ARKodeGetNumNonlinSolvConvFails", 1)) { return 1; }
  if (!fixedpoint)
  {
    flag = ARKodeGetNumLinSolvSetups(mristep_mem, &mri_nsetups);
    if (check_flag(&flag, "ARKodeGetNumLinSolvSetups", 1)) { return 1; }
    flag = ARKodeGetNumJacEvals(mristep_mem, &mri_nje);
    if (check_flag(&flag, "ARKodeGetNumJacEvals", 1)) { return 1; }
    flag = ARKodeGetNumLinRhsEvals(mristep_mem, &mri_nfeLS);
    check_flag(&flag, "ARKodeGetNumLinRhsEvals", 1);
  }
  cout << "\nMRIStep Solver Statistics:\n";
  cout << "   Return time = " << t << "\n";
  cout << "   Internal final time = " << tcur << "\n";
  cout << "   Internal solver steps = " << mri_nst << "\n";
  cout << "   Total RHS evals:  Fs = " << mri_nfsi << "\n";
  cout << "   Total number of nonlinear iterations = " << mri_nni << "\n";
  cout << "   Total number of nonlinear solver convergence failures = "
       << mri_ncfn << "\n";
  if (!fixedpoint)
  {
    cout << "   Total linear solver setups = " << mri_nsetups << "\n";
    cout << "   Total RHS evals for setting up the linear system = " << mri_nfeLS
         << "\n";
    cout << "   Total number of Jacobian evaluations = " << mri_nje << "\n";
  }

  // Compare solver statistics
  numfails = 0;
  cout << "\nComparing Solver Statistics:\n";
  if (ark_nst != mri_nst)
  {
    numfails += 1;
    cout << "  Internal solver steps error: " << ark_nst << " vs " << mri_nst
         << "\n";
  }
  if (!Compare(ark_nfi, mri_nfsi, ONE))
  {
    numfails += 1;
    cout << "  RHS evals error: " << ark_nfi << " vs " << mri_nfsi << "\n";
  }
  if (!Compare(ark_nni, mri_nni, ONE))
  {
    numfails += 1;
    cout << "  Nonlinear iterations error: " << ark_nni << " vs " << mri_nni
         << "\n";
  }
  if (ark_ncfn != mri_ncfn)
  {
    numfails += 1;
    cout << "  Nonlinear convergence failures error: " << ark_ncfn << " vs "
         << mri_ncfn << "\n";
  }
  if (!fixedpoint)
  {
    if (ark_nsetups != mri_nsetups)
    {
      numfails += 1;
      cout << "  Linear solver setups error: " << ark_nsetups << " vs "
           << mri_nsetups << "\n";
    }
    if (ark_nfeLS != mri_nfeLS)
    {
      numfails += 1;
      cout << "  RHS evals for LS error: " << ark_nfeLS << " vs " << mri_nfeLS
           << "\n";
    }
    if (ark_nje != mri_nje)
    {
      numfails += 1;
      cout << "  Jacobian evaluations error: " << ark_nje << " vs " << mri_nje
           << "\n";
    }
  }
  if (numfails) { cout << "Failed " << numfails << " tests\n"; }
  else { cout << "All tests pass!\n"; }

  // Clean up and return with successful completion
  ARKodeButcherTable_Free(B);  // Free Butcher table
  ARKodeButcherTable_Free(Bc); // Free Butcher table
  MRIStepCoupling_Free(C);     // Free MRI coupling table
  ARKodeFree(&arkstep_mem);    // Free integrator memory
  ARKodeFree(&mristep_mem);
  ARKodeFree(&inner_mem);
  MRIStepInnerStepper_Free(&inner_stepper);
  if (fixedpoint)
  {
    SUNNonlinSolFree(NLSa); // Free nonlinear solvers
    SUNNonlinSolFree(NLSm);
  }
  else
  {
    SUNLinSolFree(LSa); // Free linear solvers
    SUNLinSolFree(LSm);
    SUNMatDestroy(Aa); // Free A matrices
    SUNMatDestroy(Am);
  }
  N_VDestroy(y); // Free y vector

  SUNContext_Free(&sunctx);
  return (numfails);
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

// f routine to compute the ODE RHS function f(t,y).
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* rdata = (sunrealtype*)user_data; // cast user_data to sunrealtype
  sunrealtype lam    = rdata[0];       // set shortcut for stiffness parameter
  sunrealtype y0     = NV_Ith_S(y, 0); // access current solution values
  sunrealtype y1     = NV_Ith_S(y, 1);
  sunrealtype y2     = NV_Ith_S(y, 2);
  sunrealtype yd0, yd1, yd2;

  // fill in the RHS function: f(t,y) = V*D*Vi*y
  yd0 = SUN_RCONST(0.25) * (SUN_RCONST(5.0) * y0 + SUN_RCONST(1.0) * y1 -
                            SUN_RCONST(3.0) * y2); // yd = Vi*y
  yd1 = SUN_RCONST(0.25) *
        (SUN_RCONST(2.0) * y0 + SUN_RCONST(2.0) * y1 - SUN_RCONST(2.0) * y2);
  yd2 = SUN_RCONST(0.25) *
        (SUN_RCONST(1.0) * y0 + SUN_RCONST(1.0) * y1 + SUN_RCONST(1.0) * y2);
  y0  = -SUN_RCONST(0.5) * yd0; //  y = D*yd
  y1  = -SUN_RCONST(0.1) * yd1;
  y2  = lam * yd2;
  yd0 = SUN_RCONST(1.0) * y0 - SUN_RCONST(1.0) * y1 +
        SUN_RCONST(1.0) * y2; // yd = V*y
  yd1 = -SUN_RCONST(1.0) * y0 + SUN_RCONST(2.0) * y1 + SUN_RCONST(1.0) * y2;
  yd2 = SUN_RCONST(0.0) * y0 - SUN_RCONST(1.0) * y1 + SUN_RCONST(2.0) * y2;
  NV_Ith_S(ydot, 0) = yd0;
  NV_Ith_S(ydot, 1) = yd1;
  NV_Ith_S(ydot, 2) = yd2;

  return 0; // Return with success
}

// f0 routine to compute a zero-valued ODE RHS function f(t,y).
static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  // Initialize ydot to zero and return
  N_VConst(ZERO, ydot);
  return 0;
}

// Jacobian routine to compute J(t,y) = df/dy.
static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* rdata = (sunrealtype*)user_data; // cast user_data to sunrealtype
  sunrealtype lam    = rdata[0]; // set shortcut for stiffness parameter
  SUNMatrix V = SUNDenseMatrix(3, 3, sunctx); // create temporary SUNMatrix objects
  SUNMatrix D = SUNDenseMatrix(3, 3, sunctx); // create temporary SUNMatrix objects
  SUNMatrix Vi = SUNDenseMatrix(3, 3, sunctx); // create temporary SUNMatrix objects

  SUNMatZero(V); // initialize temporary matrices to zero
  SUNMatZero(D); // (not technically required)
  SUNMatZero(Vi);

  // Fill in temporary matrices:
  //    V = [1 -1 1; -1 2 1; 0 -1 2]
  SM_ELEMENT_D(V, 0, 0) = SUN_RCONST(1.0);
  SM_ELEMENT_D(V, 0, 1) = -SUN_RCONST(1.0);
  SM_ELEMENT_D(V, 0, 2) = SUN_RCONST(1.0);
  SM_ELEMENT_D(V, 1, 0) = -SUN_RCONST(1.0);
  SM_ELEMENT_D(V, 1, 1) = SUN_RCONST(2.0);
  SM_ELEMENT_D(V, 1, 2) = SUN_RCONST(1.0);
  SM_ELEMENT_D(V, 2, 0) = SUN_RCONST(0.0);
  SM_ELEMENT_D(V, 2, 1) = -SUN_RCONST(1.0);
  SM_ELEMENT_D(V, 2, 2) = SUN_RCONST(2.0);

  //    Vi = 0.25*[5 1 -3; 2 2 -2; 1 1 1]
  SM_ELEMENT_D(Vi, 0, 0) = SUN_RCONST(0.25) * SUN_RCONST(5.0);
  SM_ELEMENT_D(Vi, 0, 1) = SUN_RCONST(0.25) * SUN_RCONST(1.0);
  SM_ELEMENT_D(Vi, 0, 2) = -SUN_RCONST(0.25) * SUN_RCONST(3.0);
  SM_ELEMENT_D(Vi, 1, 0) = SUN_RCONST(0.25) * SUN_RCONST(2.0);
  SM_ELEMENT_D(Vi, 1, 1) = SUN_RCONST(0.25) * SUN_RCONST(2.0);
  SM_ELEMENT_D(Vi, 1, 2) = -SUN_RCONST(0.25) * SUN_RCONST(2.0);
  SM_ELEMENT_D(Vi, 2, 0) = SUN_RCONST(0.25) * SUN_RCONST(1.0);
  SM_ELEMENT_D(Vi, 2, 1) = SUN_RCONST(0.25) * SUN_RCONST(1.0);
  SM_ELEMENT_D(Vi, 2, 2) = SUN_RCONST(0.25) * SUN_RCONST(1.0);

  //    D = [-0.5 0 0; 0 -0.1 0; 0 0 lam]
  SM_ELEMENT_D(D, 0, 0) = -SUN_RCONST(0.5);
  SM_ELEMENT_D(D, 1, 1) = -SUN_RCONST(0.1);
  SM_ELEMENT_D(D, 2, 2) = lam;

  // Compute J = V*D*Vi
  if (dense_MM(D, Vi, J) != 0)
  { // J = D*Vi
    cerr << "matmul error\n";
    return 1;
  }
  if (dense_MM(V, J, D) != 0)
  { // D = V*J [= V*D*Vi]
    cerr << "matmul error\n";
    return 1;
  }
  SUNMatCopy(D, J);

  SUNMatDestroy(V);  // Free V matrix
  SUNMatDestroy(D);  // Free D matrix
  SUNMatDestroy(Vi); // Free Vi matrix

  return 0; // Return with success
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

// SUNDenseMatrix matrix-multiply utility routine: C = A*B
static int dense_MM(SUNMatrix A, SUNMatrix B, SUNMatrix C)
{
  // check for legal dimensions
  if ((SUNDenseMatrix_Columns(A) != SUNDenseMatrix_Rows(B)) ||
      (SUNDenseMatrix_Rows(C) != SUNDenseMatrix_Rows(A)) ||
      (SUNDenseMatrix_Columns(C) != SUNDenseMatrix_Columns(B)))
  {
    cerr << "\n matmul error: dimension mismatch\n\n";
    return 1;
  }

  sunrealtype** adata = SUNDenseMatrix_Cols(A); // access data and extents
  sunrealtype** bdata = SUNDenseMatrix_Cols(B);
  sunrealtype** cdata = SUNDenseMatrix_Cols(C);
  sunindextype m      = SUNDenseMatrix_Rows(C);
  sunindextype n      = SUNDenseMatrix_Columns(C);
  sunindextype l      = SUNDenseMatrix_Columns(A);
  sunindextype i, j, k;
  SUNMatZero(C); // initialize output

  // perform multiply (not optimal, but fine for 3x3 matrices)
  for (i = 0; i < m; i++)
  {
    for (j = 0; j < n; j++)
    {
      for (k = 0; k < l; k++) { cdata[i][j] += adata[i][k] * bdata[k][j]; }
    }
  }

  return 0; // Return with success
}

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

// Check if relative difference of a and b is less than tolerance
static bool Compare(long int a, long int b, sunrealtype tol)
{
  sunrealtype rel_diff = SUN_RCONST(100.0) * abs(static_cast<sunrealtype>(a - b) /
                                                 static_cast<sunrealtype>(a));

  return (rel_diff > tol) ? false : true;
}

//---- end of file ----
