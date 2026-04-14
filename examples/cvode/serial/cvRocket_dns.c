/* -----------------------------------------------------------------
 * Programmer(s): Alan C. Hindmarsh @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem, with the coding needed
 * for its solution by CVODE. The problem is a simpliflied model of a
 * rocket, ascending vertically, with mass decreasing over time. The
 * system (of size 2) is given by
 *    y_1 = rocket height H, y_1(0) = 0,
 *    y_2 = rocket velocity v, y_2(0) = 0,
 *    dH/dt = v,
 *    dv/dt = a(t,v).
 * The upward acceleration a(t,v) is given by
 *    a(t,v) = F/(M_r + M_f) - Dv - g,
 * where F = engine thrust force (constant) M_r = rocket mass without
 * fuel, M_f = fuel mass = M_f0 - r*t, r = fuel burn rate,
 * D = drag coefficient, g = gravitational acceleration.
 * The engine force is reset to 0 when the fuel mass reaches 0, or
 * when H reaches a preset height H_c, whichever happens first.
 * Rootfinding is used to locate the time at which M_f = 0 or H =
 * H_c, and also the time at which the rocket reaches its maximum
 * height, given by the condition v = 0, t > 0.
 *
 * The problem is solved with the BDF method and Dense linear solver.
 *
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------*/

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>

#include "sundials/sundials_nvector.h"

/* Problem Constants */

#define NEQ 2 /* number of equations  */

#define Force  SUN_RCONST(2200.0) /* engine force */
#define massr  SUN_RCONST(10.0)   /* rocket mass (empty) */
#define massf0 SUN_RCONST(1.0)    /* initial fuel mass */
#define brate  SUN_RCONST(0.1)    /* fuel burn rate */
#define Drag   SUN_RCONST(0.3)    /* Drag coefficient */
#define grav   SUN_RCONST(32.0)   /* acceleration due to gravity */
#define Hcut   SUN_RCONST(4000.0) /* height of engine cutoff */

#define Y1    SUN_RCONST(0.0) /* initial y components */
#define Y2    SUN_RCONST(0.0)
#define RTOL  SUN_RCONST(1.0e-5) /* scalar relative tolerance            */
#define ATOL1 SUN_RCONST(1.0e-2) /* vector absolute tolerance components */
#define ATOL2 SUN_RCONST(1.0e-1)
#define T0    SUN_RCONST(0.0) /* initial time           */
#define T1    SUN_RCONST(1.0) /* first output time      */
#define TINC  SUN_RCONST(1.0) /* output time increment  */
#define NOUT  70              /* number of output times */

#define ZERO SUN_RCONST(0.0)

/* Functions Called by the Solver */

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

static int g(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data);

static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private functions to output results */
static void PrintOutput(sunrealtype t, sunrealtype y1, sunrealtype y2);
static void PrintRootInfo(int root_f1, int root_f2, int numroot);

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
  sunrealtype reltol, t, tout;
  N_Vector y, abstol;
  SUNMatrix A;
  SUNLinearSolver LS;
  void* cvode_mem;
  int retval, retvalr, iout;
  int numroot;
  int rootsfound[2];
  sunbooleantype engine_on;
  sunrealtype *ydata, *adata;

  y = abstol = NULL;
  A          = NULL;
  LS         = NULL;
  cvode_mem  = NULL;
  ydata      = NULL;

  /* Create the SUNDIALS context */
  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return (1); }

  /* Create serial vector of length NEQ for I.C. and abstol */
  y = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)y, "N_VNew_Serial", 0)) { return (1); }
  abstol = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)abstol, "N_VNew_Serial", 0)) { return (1); }

  /* Initialize y */
  ydata    = N_VGetArrayPointer(y);
  ydata[0] = Y1;
  ydata[1] = Y2;

  /* Set the scalar relative tolerance */
  reltol = RTOL;
  /* Set the vector absolute tolerance */
  adata    = N_VGetArrayPointer(abstol);
  adata[0] = ATOL1;
  adata[1] = ATOL2;

  /* Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) { return (1); }

  /* Call CVodeInit to initialize the integrator memory and specify the right-hand side function in
   * y'=f(t,y), the initial time T0, and the initial dependent variable vector y. */
  retval = CVodeInit(cvode_mem, f, T0, y);
  if (check_retval(&retval, "CVodeInit", 1)) { return (1); }

  /* Call CVodeSVtolerances to specify the scalar relative tolerance and vector absolute tolerances */
  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) { return (1); }

  /* Provide sunbooleantype engine_on as user data for use in f and g routines */
  retval = CVodeSetUserData(cvode_mem, &engine_on);
  if (check_retval((void*)&retval, "CVodeSetUserData", 1)) { return (1); }

  /* Call CVodeRootInit to specify the root function g with 2 components */
  retval = CVodeRootInit(cvode_mem, 2, g);
  if (check_retval(&retval, "CVodeRootInit", 1)) { return (1); }

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ, sunctx);
  if (check_retval((void*)A, "SUNDenseMatrix", 0)) { return (1); }

  /* Create dense SUNLinearSolver object for use by CVode */
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) { return (1); }

  /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) { return (1); }

  /* Set the user-supplied Jacobian routine Jac */
  retval = CVodeSetJacFn(cvode_mem, Jac);
  if (check_retval(&retval, "CVodeSetJacFn", 1)) { return (1); }

  /* In loop, call CVode, print results, check for root stops, and test for error.  On the first
     root return, restart with engine turned off. Break out of loop when NOUT preset output times
     have been reached, or when the returned value of H is negative.  */
  printf(" \nAccelerating rocket problem\n\n");

  iout      = 0;
  tout      = T1;
  engine_on = SUNTRUE;
  numroot   = 2;
  while (1)
  {
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    if (check_retval(&retval, "CVode", 1)) { return (1); }

    PrintOutput(t, ydata[0], ydata[1]);

    if (engine_on && (retval == CV_ROOT_RETURN))
    { /* engine cutoff */
      retvalr = CVodeGetRootInfo(cvode_mem, rootsfound);
      if (check_retval(&retvalr, "CVodeGetRootInfo", 1)) { return (1); }
      PrintRootInfo(rootsfound[0], rootsfound[1], numroot);
      engine_on = SUNFALSE;
      numroot   = 1;
      /* Call CVodeRootInit to specify the root function g with 1 component */
      retval = CVodeRootInit(cvode_mem, 1, g);
      if (check_retval(&retval, "CVodeRootInit", 1)) { return (1); }
      /* Reinitialize the solver with current t and y values. */
      retval = CVodeReInit(cvode_mem, t, y);
      if (check_retval((void*)&retval, "CVodeReInit", 1)) { return (1); }
    }
    else if ((!engine_on) && (retval == CV_ROOT_RETURN))
    { /* max.  height */
      retvalr = CVodeGetRootInfo(cvode_mem, rootsfound);
      if (check_retval(&retvalr, "CVodeGetRootInfo", 1)) { return (1); }
      PrintRootInfo(rootsfound[0], rootsfound[1], numroot);
    }

    if (retval == CV_SUCCESS)
    {
      iout++;
      tout += TINC;
    }

    if (iout == NOUT) { break; }
    if (ydata[0] < ZERO) { break; }
  }

  /* Print some final statistics */
  PrintFinalStats(cvode_mem);

  /* Free y and abstol vectors */
  N_VDestroy(y);
  N_VDestroy(abstol);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  /* Free the linear solver memory */
  SUNLinSolFree(LS);

  /* Free the matrix memory */
  SUNMatDestroy(A);

  /* Free the SUNDIALS context */
  SUNContext_Free(&sunctx);

  return (retval);
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
  sunrealtype v, acc;
  sunbooleantype engine_on;
  sunrealtype *ydata, *yddata;

  engine_on = *((sunbooleantype*)user_data);
  ydata     = N_VGetArrayPointer(y);
  yddata    = N_VGetArrayPointer(ydot);

  v         = ydata[1];
  yddata[0] = v;

  acc = engine_on ? Force / (massr + massf0 - brate * t) : ZERO;

  yddata[1] = acc - Drag * v - grav;

  return (0);
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  /* Jdata is column-major */
  sunrealtype* Jdata = SUNDenseMatrix_Data(J);

  Jdata[1] = SUN_RCONST(1.0);
  Jdata[3] = -Drag;

  return (0);
}

/*
 * g routine. Compute functions g_i(t,y).
 */

static int g(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data)
{
  sunrealtype H, v;
  sunbooleantype engine_on;
  sunrealtype* ydata;

  engine_on = *((sunbooleantype*)user_data);

  ydata = N_VGetArrayPointer(y);

  if (engine_on)
  {
    gout[0] = massf0 - brate * t;
    H       = ydata[0];
    gout[1] = H - Hcut;
  }
  else
  {
    v       = ydata[1];
    gout[0] = v;
  }

  return (0);
}

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

static void PrintOutput(sunrealtype t, sunrealtype y1, sunrealtype y2)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %0.4Le      y =%14.6Le  %14.6Le\n", t, y1, y2);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %0.4e      y =%14.6e  %14.6e\n", t, y1, y2);
#else
  printf("At t = %0.4e      y =%14.6e  %14.6e\n", t, y1, y2);
#endif

  return;
}

static void PrintRootInfo(int root_f1, int root_f2, int numroot)
{
  if (numroot == 2)
  {
    printf("    rootsfound[] = %3d %3d\n", root_f1, root_f2);
  }
  if (numroot == 1) { printf("    rootsfound[] = %3d\n", root_f1); }

  return;
}

/*
 * Get and print some final statistics
 */

static void PrintFinalStats(void* cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
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
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

  retval = CVodeGetNumJacEvals(cvode_mem, &nje);
  check_retval(&retval, "CVodeGetNumJacEvals", 1);
  retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
  check_retval(&retval, "CVodeGetNumLinRhsEvals", 1);

  retval = CVodeGetNumGEvals(cvode_mem, &nge);
  check_retval(&retval, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = % ld\n",
         nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n", nni, ncfn,
         netf, nge);
}

/*
 * Check function return value... opt == 0 means SUNDIALS function allocates memory so check if
 *   returned NULL pointer opt == 1 means SUNDIALS function returns an integer value so check if
 *   retval < 0 opt == 2 means function allocates memory so check if returned NULL pointer
 */

static int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n ",
            funcname);
    return (1);
  }
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
  else if (opt == 2 && returnvalue == NULL)
  {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n ",
            funcname);
    return (1);
  }

  return (0);
}
