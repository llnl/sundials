/* -----------------------------------------------------------------
 * Programmer(s): Carol Woodward.
 *      Based on cvRoberts_dns.c and modified to use KLU.
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
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
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODE. The problem is from
 * chemical kinetics, and consists of the following three rate
 * equations:
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2
 *    dy3/dt = 3.e7*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 * While integrating the system, we also use the rootfinding
 * feature to find the points at which y1 = 1e-4 or at which
 * y3 = 0.01. This program solves the problem with the BDF method,
 * Newton iteration with the KLU sparse direct linear solver, and a
 * user-supplied Jacobian routine.
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance. Output is printed in decades from t = .4 to t = 4.e10.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------*/

#include <cvodes/cvodes.h>          /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h> /* access to serial N_Vector            */
#include <stdio.h>
#include <sunlinsol/sunlinsol_klu.h> /* access to KLU sparse direct solver   */
#include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix           */

/* User-defined vector accessor macro: Ith */

/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0. */

#define Ith(v, i) NV_Ith_S(v, i - 1) /* i-th vector component, i=1..NEQ */

/* Problem Constants */

#define NEQ   3               /* number of equations  */
#define NNZ   7               /* number of non-zero entries in the Jacobian */
#define Y1    SUN_RCONST(1.0) /* initial y components */
#define Y2    SUN_RCONST(0.0)
#define Y3    SUN_RCONST(0.0)
#define RTOL  SUN_RCONST(1.0e-4) /* scalar relative tolerance            */
#define ATOL1 SUN_RCONST(1.0e-8) /* vector absolute tolerance components */
#define ATOL2 SUN_RCONST(1.0e-14)
#define ATOL3 SUN_RCONST(1.0e-6)
#define T0    SUN_RCONST(0.0)  /* initial time           */
#define T1    SUN_RCONST(0.4)  /* first output time      */
#define TMULT SUN_RCONST(10.0) /* output time factor     */
#define NOUT  12               /* number of output times */

/* Functions Called by the Solver */

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

static int g(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data);

static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private functions to output results */

static void PrintOutput(sunrealtype t, sunrealtype y1, sunrealtype y2,
                        sunrealtype y3);
static void PrintRootInfo(int root_f1, int root_f2);

/* Private function to print final statistics */

static void PrintFinalStats(void* cvode_mem);

/* Private function to check function return values */

static int check_retval(void* returnvalue, const char* funcname, int opt);

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main(void)
{
  SUNContext sunctx;
  sunrealtype t, tout;
  N_Vector y;
  N_Vector abstol;
  SUNMatrix A;
  SUNLinearSolver LS;
  void* cvode_mem;
  int retval, iout;
  int retvalr;
  int rootsfound[2];

  y         = NULL;
  abstol    = NULL;
  A         = NULL;
  LS        = NULL;
  cvode_mem = NULL;

  /* Create the SUNDIALS context */
  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return (1); }

  /* Initial conditions */
  y = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)y, "N_VNew_Serial", 0)) { return (1); }

  /* Initialize y */
  Ith(y, 1) = Y1;
  Ith(y, 2) = Y2;
  Ith(y, 3) = Y3;

  /* Set the vector absolute tolerance */
  abstol = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)abstol, "N_VNew_Serial", 0)) { return (1); }

  Ith(abstol, 1) = ATOL1;
  Ith(abstol, 2) = ATOL2;
  Ith(abstol, 3) = ATOL3;

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) { return (1); }

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the initial time T0, and
   * the initial dependent variable vector y. */
  retval = CVodeInit(cvode_mem, f, T0, y);
  if (check_retval(&retval, "CVodeInit", 1)) { return (1); }

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  retval = CVodeSVtolerances(cvode_mem, RTOL, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) { return (1); }

  /* Call CVodeRootInit to specify the root function g with 2 components */
  retval = CVodeRootInit(cvode_mem, 2, g);
  if (check_retval(&retval, "CVodeRootInit", 1)) { return (1); }

  /* Create sparse SUNMatrix for use in linear solves */
  A = SUNSparseMatrix(NEQ, NEQ, NNZ, SUN_CSC_MAT, sunctx);
  if (check_retval((void*)A, "SUNSparseMatrix", 0)) { return (1); }

  /* Create KLU solver object for use by CVode */
  LS = SUNLinSol_KLU(y, A, sunctx);
  if (check_retval((void*)LS, "SUNLinSol_KLU", 0)) { return (1); }

  /* Attach the matrix and linear solver */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) { return (1); }

  /* Set the user-supplied Jacobian routine Jac */
  retval = CVodeSetJacFn(cvode_mem, Jac);
  if (check_retval(&retval, "CVodeSetJacFn", 1)) { return (1); }

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
  printf(" \n3-species kinetics problem\n\n");

  iout = 0;
  tout = T1;
  while (1)
  {
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    PrintOutput(t, Ith(y, 1), Ith(y, 2), Ith(y, 3));

    if (retval == CV_ROOT_RETURN)
    {
      retvalr = CVodeGetRootInfo(cvode_mem, rootsfound);
      if (check_retval(&retvalr, "CVodeGetRootInfo", 1)) { return (1); }
      PrintRootInfo(rootsfound[0], rootsfound[1]);
    }

    if (check_retval(&retval, "CVode", 1)) { break; }
    if (retval == CV_SUCCESS)
    {
      iout++;
      tout *= TMULT;
    }

    if (iout == NOUT) { break; }
  }

  /* Print some final statistics */
  PrintFinalStats(cvode_mem);

  /* Free memory */
  N_VDestroy(y);            /* Free y vector */
  N_VDestroy(abstol);       /* Free abstol vector */
  CVodeFree(&cvode_mem);    /* Free CVODES memory */
  SUNLinSolFree(LS);        /* Free the linear solver memory */
  SUNMatDestroy(A);         /* Free the matrix memory */
  SUNContext_Free(&sunctx); /* Free the SUNDIALS context */

  return (0);
}

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/*
 * f routine. Compute function f(t,y).
 */

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype y1, y2, y3, yd1, yd3;

  y1 = Ith(y, 1);
  y2 = Ith(y, 2);
  y3 = Ith(y, 3);

  yd1 = Ith(ydot, 1) = SUN_RCONST(-0.04) * y1 + SUN_RCONST(1.0e4) * y2 * y3;
  yd3 = Ith(ydot, 3) = SUN_RCONST(3.0e7) * y2 * y2;
  Ith(ydot, 2)       = -yd1 - yd3;

  return (0);
}

/*
 * g routine. Compute functions g_i(t,y) for i = 0,1.
 */

static int g(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data)
{
  sunrealtype y1, y3;

  y1      = Ith(y, 1);
  y3      = Ith(y, 3);
  gout[0] = y1 - SUN_RCONST(0.0001);
  gout[1] = y3 - SUN_RCONST(0.01);

  return (0);
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  /* State at which to evaluate the Jacobian */
  sunrealtype* yval = N_VGetArrayPointer(y);

  /* J is stored in CSC format:
     data    = non-zero matrix entries stored column-wise (length NNZ)
     rowvals = row index for each non-zero matrix entry (length NNZ)
     colptrs = i-th entry is the index in data where the first non-zero matrix
               entry of the i-th column is stored (length NEQ + 1) */
  sunrealtype* data     = SUNSparseMatrix_Data(J);
  sunindextype* rowvals = SUNSparseMatrix_IndexValues(J);
  sunindextype* colptrs = SUNSparseMatrix_IndexPointers(J);

  /* first column entries start at data[0], two entries (rows 0 and 1) */
  colptrs[0] = 0;

  rowvals[0] = 0;
  data[0]    = SUN_RCONST(-0.04);

  rowvals[1] = 1;
  data[1]    = SUN_RCONST(0.04);

  /* second column entries start at data[2], three entries (rows 0, 1, and 2) */
  colptrs[1] = 2;

  rowvals[2] = 0;
  data[2]    = SUN_RCONST(1.0e4) * yval[2];

  rowvals[3] = 1;
  data[3]    = (SUN_RCONST(-1.0e4) * yval[2]) - (SUN_RCONST(6.0e7) * yval[1]);

  rowvals[4] = 2;
  data[4]    = SUN_RCONST(6.0e7) * yval[1];

  /* third column entries start at data[5], two entries (rows 0 and 1) */
  colptrs[2] = 5;

  rowvals[5] = 0;
  data[5]    = SUN_RCONST(1.0e4) * yval[1];

  rowvals[6] = 1;
  data[6]    = SUN_RCONST(-1.0e4) * yval[1];

  /* number of non-zeros */
  colptrs[3] = 7;

  return (0);
}

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

static void PrintOutput(sunrealtype t, sunrealtype y1, sunrealtype y2,
                        sunrealtype y3)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %0.4Le      y =%14.6Le  %14.6Le  %14.6Le\n", t, y1, y2, y3);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#else
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#endif

  return;
}

static void PrintRootInfo(int root_f1, int root_f2)
{
  printf("    rootsfound[] = %3d %3d\n", root_f1, root_f2);

  return;
}

/*
 * Get and print some final statistics
 */

static void PrintFinalStats(void* cvode_mem)
{
  long int nst, nfe, nsetups, nje, nni, nnf, ncfn, netf, nge;
  int retval;

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &nnf);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);
  retval = CVodeGetNumStepSolveFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumStepSolveFails", 1);

  retval = CVodeGetNumJacEvals(cvode_mem, &nje);
  check_retval(&retval, "CVodeGetNumJacEvals", 1);

  retval = CVodeGetNumGEvals(cvode_mem, &nge);
  check_retval(&retval, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe = %-6ld nsetups = %-6ld nje = %ld\n", nst, nfe,
         nsetups, nje);
  printf("nni = %-6ld nnf = %-6ld netf = %-6ld    ncfn = %-6ld nge = %ld\n\n",
         nni, nnf, netf, ncfn, nge);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return (1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  return (0);
}
