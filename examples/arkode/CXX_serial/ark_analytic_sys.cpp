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
 * Example problem:
 *
 * The following is a simple example problem with analytical
 * solution,
 *    dy/dt = A*y
 * where A = V*D*Vi,
 *      V = [1 -1 1; -1 2 1; 0 -1 2];
 *      Vi = 0.25*[5 1 -3; 2 2 -2; 1 1 1];
 *      D = [-0.5 0 0; 0 -0.1 0; 0 0 lam];
 * where lam is a large negative number. The analytical solution to
 * this problem is
 *   Y(t) = V*exp(D*t)*Vi*Y0
 * for t in the interval [0.0, 0.05], with initial condition:
 * y(0) = [1,1,1]'.
 *
 * The stiffness of the problem is directly proportional to the
 * value of "lambda".  The value of lambda should be negative to
 * result in a well-posed ODE; for values with magnitude larger than
 * 100 the problem becomes quite stiff.
 *
 * In this example, we choose lambda = -100.
 *
 * This program solves the problem with the DIRK method,
 * Newton iteration with the dense linear solver, and a
 * user-supplied Jacobian routine.
 * Output is printed every 1.0 units of time (10 total).
 * Run statistics (optional outputs) are printed at the end.
 *-----------------------------------------------------------------*/

// Header files
#include <arkode/arkode_arkstep.h> // prototypes for ARKStep fcts., consts.
#include <cmath>
#include <iostream>
#include <nvector/nvector_serial.h> // access to serial N_Vector
#include <stdio.h>
#include <string.h>
#include <sundials/sundials_core.hpp>
#include <sundials/sundials_logger.h>
#include <sundials/sundials_types.h>   // def. of type 'sunrealtype'
#include <sunlinsol/sunlinsol_dense.h> // access to dense SUNLinearSolver
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define ESYM "Le"
#define FSYM "Lf"
#else
#define ESYM "e"
#define FSYM "f"
#endif

using namespace std;

// User-supplied Functions Called by the Solver
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Private function to perform matrix-matrix product
static int SUNDlsMat_dense_MM(SUNMatrix A, SUNMatrix B, SUNMatrix C);

// Private function to check function return values
static int check_flag(void* flagvalue, const string funcname, int opt);

// Main Program
int main()
{
  // general problem parameters
  sunrealtype T0     = SUN_RCONST(0.0);    // initial time
  sunrealtype Tf     = SUN_RCONST(0.05);   // final time
  sunrealtype dTout  = SUN_RCONST(0.005);  // time between outputs
  sunindextype NEQ   = 3;                  // number of dependent vars.
  sunrealtype reltol = SUN_RCONST(1.0e-6); // tolerances
  sunrealtype abstol = SUN_RCONST(1.0e-10);
  sunrealtype lambda = SUN_RCONST(-100.0); // stiffness parameter

  // general problem variables
  int flag;                  // reusable error-checking flag
  N_Vector y         = NULL; // empty vector for storing solution
  SUNMatrix A        = NULL; // empty dense matrix for solver
  SUNLinearSolver LS = NULL; // empty dense linear solver
  void* arkode_mem   = NULL; // empty ARKode memory structure

  // Initial problem output
  cout << "\nAnalytical ODE test problem:\n";
  cout << "   lambda = " << lambda << "\n";
  cout << "   reltol = " << reltol << "\n";
  cout << "   abstol = " << abstol << "\n\n";

  // Create the SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Configure SUNDIALS logging via the runtime API.
  // This requires that SUNDIALS was configured with the CMake options
  //   SUNDIALS_LOGGING_LEVEL=n
  // where n is one of:
  //   1 --> log only errors,
  //   2 --> log errors + warnings,
  //   3 --> log errors + warnings + info output
  //   4 --> all of the above plus debugging output like internal integrator values
  //   5 --> all of the above and even more
  // SUNDIALS will only log up to the max level n, but a lesser level can
  // be configured at runtime by only providing output files for the
  // desired levels. We will enable all logging here on the condition
  // that SUNDIALS was built with the correct logging level enabled.
  SUNLogger logger = NULL;
  if (SUNDIALS_LOGGING_LEVEL >= SUN_LOGLEVEL_ERROR)
  {
    flag = SUNLogger_Create(SUN_COMM_NULL, // no MPI communicator
                            0,             // output on process 0 (the only one)
                            &logger);
    if (check_flag(&flag, "SUNLogger_Create", 1)) { return 1; }

    // Attach the logger
    flag = SUNContext_SetLogger(sunctx, logger);
    if (check_flag(&flag, "SUNContext_SetLogger", 1)) { return 1; }

    // Setup log files
    flag = SUNLogger_SetErrorFilename(logger, "stderr");
    if (check_flag(&flag, "SUNLogger_SetErrorFilename", 1)) { return 1; }
  }

  if (SUNDIALS_LOGGING_LEVEL >= SUN_LOGLEVEL_WARNING)
  {
    flag = SUNLogger_SetWarningFilename(logger, "stderr");
    if (check_flag(&flag, "SUNLogger_SetWarningFilename", 1)) { return 1; }
  }

  if (SUNDIALS_LOGGING_LEVEL >= SUN_LOGLEVEL_INFO)
  {
    flag = SUNLogger_SetInfoFilename(logger, "ark_analytic_sys.info.log");
    if (check_flag(&flag, "SUNLogger_SetInfoFilename", 1)) { return 1; }
  }

  if (SUNDIALS_LOGGING_LEVEL >= SUN_LOGLEVEL_DEBUG)
  {
    flag = SUNLogger_SetDebugFilename(logger, "stderr");
    if (check_flag(&flag, "SUNLogger_SetDebugFilename", 1)) { return 1; }
  }

  // Initialize vector data structure and specify initial condition
  y = N_VNew_Serial(NEQ, sunctx);
  if (check_flag((void*)y, "N_VNew_Serial", 0)) { return 1; }
  NV_Ith_S(y, 0) = 1.0;
  NV_Ith_S(y, 1) = 1.0;
  NV_Ith_S(y, 2) = 1.0;

  // Initialize dense matrix data structure and solver
  A = SUNDenseMatrix(NEQ, NEQ, sunctx);
  if (check_flag((void*)A, "SUNDenseMatrix", 0)) { return 1; }

  LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_flag((void*)LS, "SUNLinSol_Dense", 0)) { return 1; }

  /* Call ARKStepCreate to initialize the ARK timestepper memory and
     specify the right-hand side function in y'=f(t,y), the initial time
     T0, and the initial dependent variable vector y.  Note: since
     this problem is fully implicit, we set f_E to NULL and f_I to f. */
  arkode_mem = ARKStepCreate(NULL, f, T0, y, sunctx);
  if (check_flag((void*)arkode_mem, "ARKStepCreate", 0)) { return 1; }

  // Set routines
  flag = ARKodeSetUserData(arkode_mem,
                           (void*)&lambda); // Pass lambda to user functions
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }

  flag = ARKodeSStolerances(arkode_mem, reltol, abstol); // Specify tolerances
  if (check_flag(&flag, "ARKodeSStolerances", 1)) { return 1; }

  // Linear solver interface
  flag = ARKodeSetLinearSolver(arkode_mem, LS,
                               A); // Attach matrix and linear solver
  if (check_flag(&flag, "ARKodeSetLinearSolver", 1)) { return 1; }

  flag = ARKodeSetJacFn(arkode_mem, Jac); // Set Jacobian routine
  if (check_flag(&flag, "ARKodeSetJacFn", 1)) { return 1; }

  // Specify linearly implicit RHS, with non-time-dependent Jacobian
  flag = ARKodeSetLinear(arkode_mem, 0);
  if (check_flag(&flag, "ARKodeSetLinear", 1)) { return 1; }

  flag = ARKodeSetAutonomous(arkode_mem, SUNTRUE);
  if (check_flag(&flag, "ARKodeSetAutonomous", 1)) { return 1; }

  // Open output stream for results, output comment line
  FILE* UFID = fopen("solution.txt", "w");
  fprintf(UFID, "# t y1 y2 y3\n");

  // output initial condition to disk
  fprintf(UFID, " %.16" ESYM " %.16" ESYM " %.16" ESYM " %.16" ESYM "\n", T0,
          NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));

  /* Main time-stepping loop: calls ARKodeEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached */
  sunrealtype t    = T0;
  sunrealtype tout = T0 + dTout;
  cout << "      t        y0        y1        y2\n";
  cout << "   --------------------------------------\n";
  while (Tf - t > 1.0e-15)
  {
    flag = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_NORMAL); // call integrator
    if (check_flag(&flag, "ARKodeEvolve", 1)) { break; }
    printf("  %8.4" FSYM "  %8.5" FSYM "  %8.5" FSYM "  %8.5" FSYM
           "\n", // access/print solution
           t, NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));
    fprintf(UFID, " %.16" ESYM " %.16" ESYM " %.16" ESYM " %.16" ESYM "\n", t,
            NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));
    if (flag >= 0)
    { // successful solve: update time
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
    else
    { // unsuccessful solve: break
      fprintf(stderr, "Solver failure, stopping integration\n");
      break;
    }
  }
  cout << "   --------------------------------------\n";
  fclose(UFID);

  // Print some final statistics
  long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;
  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  check_flag(&flag, "ARKodeGetNumSteps", 1);
  flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  check_flag(&flag, "ARKodeGetNumStepAttempts", 1);
  flag = ARKodeGetNumRhsEvals(arkode_mem, 0, &nfe);
  check_flag(&flag, "ARKodeGetNumRhsEvals", 1);
  flag = ARKodeGetNumRhsEvals(arkode_mem, 1, &nfi);
  check_flag(&flag, "ARKodeGetNumRhsEvals", 1);
  flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
  check_flag(&flag, "ARKodeGetNumLinSolvSetups", 1);
  flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
  check_flag(&flag, "ARKodeGetNumErrTestFails", 1);
  flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
  check_flag(&flag, "ARKodeGetNumNonlinSolvIters", 1);
  flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  check_flag(&flag, "ARKodeGetNumNonlinSolvConvFails", 1);
  flag = ARKodeGetNumJacEvals(arkode_mem, &nje);
  check_flag(&flag, "ARKodeGetNumJacEvals", 1);
  flag = ARKodeGetNumLinRhsEvals(arkode_mem, &nfeLS);
  check_flag(&flag, "ARKodeGetNumLinRhsEvals", 1);

  cout << "\nFinal Solver Statistics:\n";
  cout << "   Internal solver steps = " << nst << " (attempted = " << nst_a
       << ")\n";
  cout << "   Total RHS evals:  Fe = " << nfe << ",  Fi = " << nfi << "\n";
  cout << "   Total linear solver setups = " << nsetups << "\n";
  cout << "   Total RHS evals for setting up the linear system = " << nfeLS
       << "\n";
  cout << "   Total number of Jacobian evaluations = " << nje << "\n";
  cout << "   Total number of Newton iterations = " << nni << "\n";
  cout << "   Total number of linear solver convergence failures = " << ncfn
       << "\n";
  cout << "   Total number of error test failures = " << netf << "\n\n";

  // Clean up and return with successful completion
  ARKodeFree(&arkode_mem); // Free integrator memory
  SUNLinSolFree(LS);       // Free linear solver
  SUNMatDestroy(A);        // Free A matrix
  N_VDestroy(y);           // Free y vector
  if (logger)
  {
    SUNLogger_Destroy(&logger); // Free logger
  }
  return 0;
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
  yd0               = 0.25 * (5.0 * y0 + 1.0 * y1 - 3.0 * y2); // yd = Vi*y
  yd1               = 0.25 * (2.0 * y0 + 2.0 * y1 - 2.0 * y2);
  yd2               = 0.25 * (1.0 * y0 + 1.0 * y1 + 1.0 * y2);
  y0                = -0.5 * yd0; //  y = D*yd
  y1                = -0.1 * yd1;
  y2                = lam * yd2;
  yd0               = 1.0 * y0 - 1.0 * y1 + 1.0 * y2; // yd = V*y
  yd1               = -1.0 * y0 + 2.0 * y1 + 1.0 * y2;
  yd2               = 0.0 * y0 - 1.0 * y1 + 2.0 * y2;
  NV_Ith_S(ydot, 0) = yd0;
  NV_Ith_S(ydot, 1) = yd1;
  NV_Ith_S(ydot, 2) = yd2;

  return 0; // Return with success
}

// Jacobian routine to compute J(t,y) = df/dy.
static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* rdata = (sunrealtype*)user_data; // cast user_data to sunrealtype
  sunrealtype lam    = rdata[0]; // set shortcut for stiffness parameter

  // Get Jacobian context
  SUNContext sunctx = J->sunctx;

  // create temporary SUNMatrix objects
  SUNMatrix V  = SUNDenseMatrix(3, 3, sunctx);
  SUNMatrix D  = SUNDenseMatrix(3, 3, sunctx);
  SUNMatrix Vi = SUNDenseMatrix(3, 3, sunctx);

  // initialize temporary matrices to zero (not technically required)
  SUNMatZero(V);
  SUNMatZero(D);
  SUNMatZero(Vi);

  // Fill in temporary matrices:
  //    V = [1 -1 1; -1 2 1; 0 -1 2]
  SM_ELEMENT_D(V, 0, 0) = 1.0;
  SM_ELEMENT_D(V, 0, 1) = -1.0;
  SM_ELEMENT_D(V, 0, 2) = 1.0;
  SM_ELEMENT_D(V, 1, 0) = -1.0;
  SM_ELEMENT_D(V, 1, 1) = 2.0;
  SM_ELEMENT_D(V, 1, 2) = 1.0;
  SM_ELEMENT_D(V, 2, 0) = 0.0;
  SM_ELEMENT_D(V, 2, 1) = -1.0;
  SM_ELEMENT_D(V, 2, 2) = 2.0;

  //    Vi = 0.25*[5 1 -3; 2 2 -2; 1 1 1]
  SM_ELEMENT_D(Vi, 0, 0) = 0.25 * 5.0;
  SM_ELEMENT_D(Vi, 0, 1) = 0.25 * 1.0;
  SM_ELEMENT_D(Vi, 0, 2) = -0.25 * 3.0;
  SM_ELEMENT_D(Vi, 1, 0) = 0.25 * 2.0;
  SM_ELEMENT_D(Vi, 1, 1) = 0.25 * 2.0;
  SM_ELEMENT_D(Vi, 1, 2) = -0.25 * 2.0;
  SM_ELEMENT_D(Vi, 2, 0) = 0.25 * 1.0;
  SM_ELEMENT_D(Vi, 2, 1) = 0.25 * 1.0;
  SM_ELEMENT_D(Vi, 2, 2) = 0.25 * 1.0;

  //    D = [-0.5 0 0; 0 -0.1 0; 0 0 lam]
  SM_ELEMENT_D(D, 0, 0) = -0.5;
  SM_ELEMENT_D(D, 1, 1) = -0.1;
  SM_ELEMENT_D(D, 2, 2) = lam;

  // Compute J = V*D*Vi
  if (SUNDlsMat_dense_MM(D, Vi, J) != 0)
  { // J = D*Vi
    cerr << "matmul error\n";
    return 1;
  }
  if (SUNDlsMat_dense_MM(V, J, D) != 0)
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

// SUNDenseMatrix matrix-multiply utility routine: C = A*B.
static int SUNDlsMat_dense_MM(SUNMatrix A, SUNMatrix B, SUNMatrix C)
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
      for (k = 0; k < l; k++) { cdata[j][i] += adata[k][i] * bdata[j][k]; }
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

//---- end of file ----
