/* -----------------------------------------------------------------
 * Programmers: Radu Serban and Alan Hindmarsh @ LLNL
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
 * Modification of the CVODE example cvRoberts_dns to illustrate
 * the treatment of unphysical solution components through the RHS
 * function return retval.
 *
 * Note that, to make possible negative solution components, the
 * absolute tolerances had to be loosened a bit from their values
 * in cvRoberts_dns.
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
 * -----------------------------------------------------------------*/

#include <cvode/cvode.h>            /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h> /* access to serial N_Vector            */
#include <stdio.h>
#include <sundials/sundials_types.h> /* defs. of sunrealtype, sunindextype      */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */

/* Problem Constants */

#define NEQ   3               /* number of equations  */
#define Y1    SUN_RCONST(1.0) /* initial y components */
#define Y2    SUN_RCONST(0.0)
#define Y3    SUN_RCONST(0.0)
#define RTOL  SUN_RCONST(1.0e-4) /* scalar relative tolerance            */
#define ATOL1 SUN_RCONST(1.0e-7) /* vector absolute tolerance components */
#define ATOL2 SUN_RCONST(1.0e-13)
#define ATOL3 SUN_RCONST(1.0e-5)
#define T0    SUN_RCONST(0.0)  /* initial time           */
#define T1    SUN_RCONST(0.4)  /* first output time      */
#define TMULT SUN_RCONST(10.0) /* output time factor     */
#define NOUT  14               /* number of output times */

/* Functions Called by the Solver */

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/* Private functions to output results */

static void PrintOutput(sunrealtype t, sunrealtype y1, sunrealtype y2,
                        sunrealtype y3);

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
  sunbooleantype check_negative;

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
  NV_Ith_S(y, 0) = Y1;
  NV_Ith_S(y, 1) = Y2;
  NV_Ith_S(y, 2) = Y3;

  /* Set the vector absolute tolerance */
  abstol = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)abstol, "N_VNew_Serial", 0)) { return (1); }

  NV_Ith_S(abstol, 0) = ATOL1;
  NV_Ith_S(abstol, 1) = ATOL2;
  NV_Ith_S(abstol, 2) = ATOL3;

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

  /* Call CVodeSetUserData to pass the check negative retval as user data */
  retval = CVodeSetUserData(cvode_mem, &check_negative);
  if (check_retval(&retval, "CVodeSetUserData", 1)) { return (1); }

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ, sunctx);
  if (check_retval((void*)A, "SUNDenseMatrix", 0)) { return (1); }

  /* Create dense SUNLinearSolver object for use by CVode */
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) { return (1); }

  /* Attach the matrix and linear solver */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) { return (1); }

  /* Case 1: ignore negative solution components */
  printf("Ignore negative solution components\n\n");
  check_negative = SUNFALSE;
  /* In loop, call CVode in CV_NORMAL mode */
  iout = 0;
  tout = T1;
  while (1)
  {
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    PrintOutput(t, NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));
    iout++;
    tout *= TMULT;
    if (iout == NOUT) { break; }
  }
  /* Print some final statistics */
  PrintFinalStats(cvode_mem);

  /* Case 2: intercept negative solution components */
  printf("Intercept negative solution components\n\n");
  check_negative = SUNTRUE;
  /* Reinitialize solver */
  NV_Ith_S(y, 0) = Y1;
  NV_Ith_S(y, 1) = Y2;
  NV_Ith_S(y, 2) = Y3;
  retval         = CVodeReInit(cvode_mem, T0, y);
  /* In loop, call CVode in CV_NORMAL mode */
  iout = 0;
  tout = T1;
  while (1)
  {
    CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    PrintOutput(t, NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));
    iout++;
    tout *= TMULT;
    if (iout == NOUT) { break; }
  }
  /* Print some final statistics */
  PrintFinalStats(cvode_mem);

  /* Free memory */
  N_VDestroy(y);            /* Free y vector */
  N_VDestroy(abstol);       /* Free abstol vector */
  CVodeFree(&cvode_mem);    /* Free CVODE memory */
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
  sunbooleantype* check_negative;

  check_negative = (sunbooleantype*)user_data;

  y1 = NV_Ith_S(y, 0);
  y2 = NV_Ith_S(y, 1);
  y3 = NV_Ith_S(y, 2);

  if (*check_negative && (y1 < 0 || y2 < 0 || y3 < 0)) { return (1); }

  yd1 = NV_Ith_S(ydot, 0) = SUN_RCONST(-0.04) * y1 + SUN_RCONST(1.0e4) * y2 * y3;
  yd3 = NV_Ith_S(ydot, 2) = SUN_RCONST(3.0e7) * y2 * y2;
  NV_Ith_S(ydot, 1)       = -yd1 - yd3;

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

static void PrintFinalStats(void* cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, nnf, ncfn, netf;
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
  retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
  check_retval(&retval, "CVodeGetNumLinRhsEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
         nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld nnf = %-6ld netf = %-6ld    ncfn = %-6ld\n\n", nni, nnf,
         netf, ncfn);
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
