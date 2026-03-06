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
 * The following test simulates a brusselator problem from chemical
 * kinetics.  This is an ODE system with 3 components, Y = [u,v,w],
 * satisfying the equations,
 *    du/dt = a - (w+1)*u + v*u^2
 *    dv/dt = w*u - v*u^2
 *    dw/dt = (b-w)/ep - w*u
 * for t in the interval [0.0, 10.0], with initial conditions
 * Y0 = [u0,v0,w0].
 *
 * We have 3 different testing scenarios:
 *
 * Test 1:  u0=3.9,  v0=1.1,  w0=2.8,  a=1.2,  b=2.5,  ep=1.0e-5
 *    Here, all three components exhibit a rapid transient change
 *    during the first 0.2 time units, followed by a slow and
 *    smooth evolution.
 *
 * Test 2:  u0=1.2,  v0=3.1,  w0=3,  a=1,  b=3.5,  ep=5.0e-6
 *    Here, w experiences a fast initial transient, jumping 0.5
 *    within a few steps.  All values proceed smoothly until
 *    around t=6.5, when both u and v undergo a sharp transition,
 *    with u increaseing from around 0.5 to 5 and v decreasing
 *    from around 6 to 1 in less than 0.5 time units.  After this
 *    transition, both u and v continue to evolve somewhat
 *    rapidly for another 1.4 time units, and finish off smoothly.
 *
 * Test 3:  u0=3,  v0=3,  w0=3.5,  a=0.5,  b=3,  ep=5.0e-4
 *    Here, all components undergo very rapid initial transients
 *    during the first 0.3 time units, and all then proceed very
 *    smoothly for the remainder of the simulation.
 *
 * This file is hard-coded to use test 3.
 *
 * This program solves the problem with the ARK method, using an
 * accelerated fixed-point iteration for the nonlinear solver.
 *
 * 100 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 *-----------------------------------------------------------------*/

/* Header files */
#include <arkode/arkode_arkstep.h> /* prototypes for ARKStep fcts., consts */
#include <math.h>
#include <nvector/nvector_serial.h> /* serial N_Vector types, fcts., macros */
#include <stdio.h>
#include <sundials/sundials_logger.h>
#include <sundials/sundials_types.h> /* def. of type 'sunrealtype'              */
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h> /* access to FP nonlinear solver        */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* User-supplied Functions Called by the Solver */
static int fe(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/* Private function to check function return values */
static int check_flag(void* flagvalue, const char* funcname, int opt);

/* Main Program */
int main(int argc, char* argv[])
{
  /* general problem parameters */
  sunrealtype T0     = SUN_RCONST(0.0);       /* initial time */
  sunrealtype Tf     = SUN_RCONST(10.0);      /* final time */
  sunrealtype dTout  = SUN_RCONST(1.0);       /* time between outputs */
  sunindextype NEQ   = 3;                     /* number of dependent vars. */
  int Nt             = (int)ceil(Tf / dTout); /* number of output times */
  int test           = 3;                     /* test problem to run */
  sunrealtype reltol = 1.0e-6;                /* tolerances */
  sunrealtype abstol = 1.0e-10;
  int fp_m           = 3;  /* dimension of acceleration subspace */
  int maxcor         = 10; /* maximum # of nonlinear iterations/step */
  sunrealtype a, b, ep, u0, v0, w0;
  sunrealtype rdata[3];
  int monitor = 0; /* turn on/off monitoring */

  /* general problem variables */
  SUNContext ctx         = NULL;
  SUNLogger logger       = NULL;
  const char* info_fname = "ark_brusselator_fp-info.txt";
  int flag;                      /* reusable error-checking flag */
  N_Vector y             = NULL; /* empty vector for storing solution */
  SUNNonlinearSolver NLS = NULL; /* empty nonlinear solver object */
  void* arkode_mem       = NULL; /* empty ARKode memory structure */
  FILE* UFID;
  sunrealtype t, tout;
  int iout;
  long int nst, nst_a, nfe, nfi, nni, ncfn, netf;

  /* read inputs */
  if (argc == 2) { monitor = atoi(argv[1]); }

  /* create SUNDIALS context and a logger which will record
     nonlinear solver info (e.g., residual) amongst other things. */

  flag = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_flag(&flag, "SUNContext_Create", 1)) { return 1; }

  flag = SUNLogger_Create(SUN_COMM_NULL, 0, &logger);
  if (check_flag(&flag, "SUNLogger_Create", 1)) { return 1; }

  if (monitor)
  {
    flag = SUNLogger_SetInfoFilename(logger, info_fname);
    if (check_flag(&flag, "SUNLogger_SetInfoFilename", 1)) { return 1; }
  }

  flag = SUNContext_SetLogger(ctx, logger);
  if (check_flag(&flag, "SUNContext_SetLogger", 1)) { return 1; }

  /* set up the test problem according to the desired test */
  if (test == 1)
  {
    u0 = SUN_RCONST(3.9);
    v0 = SUN_RCONST(1.1);
    w0 = SUN_RCONST(2.8);
    a  = SUN_RCONST(1.2);
    b  = SUN_RCONST(2.5);
    ep = SUN_RCONST(1.0e-5);
  }
  else if (test == 3)
  {
    u0 = SUN_RCONST(3.0);
    v0 = SUN_RCONST(3.0);
    w0 = SUN_RCONST(3.5);
    a  = SUN_RCONST(0.5);
    b  = SUN_RCONST(3.0);
    ep = SUN_RCONST(5.0e-4);
  }
  else
  {
    u0 = SUN_RCONST(1.2);
    v0 = SUN_RCONST(3.1);
    w0 = SUN_RCONST(3.0);
    a  = SUN_RCONST(1.0);
    b  = SUN_RCONST(3.5);
    ep = SUN_RCONST(5.0e-6);
  }

  /* Initial problem output */
  printf("\nBrusselator ODE test problem, fixed-point solver:\n");
  printf("    initial conditions:  u0 = %" GSYM ",  v0 = %" GSYM
         ",  w0 = %" GSYM "\n",
         u0, v0, w0);
  printf("    problem parameters:  a = %" GSYM ",  b = %" GSYM ",  ep = %" GSYM
         "\n",
         a, b, ep);
  printf("    reltol = %.1" ESYM ",  abstol = %.1" ESYM "\n\n", reltol, abstol);

  /* Initialize data structures */
  rdata[0] = a; /* set user data  */
  rdata[1] = b;
  rdata[2] = ep;
  y        = N_VNew_Serial(NEQ, ctx); /* Create serial vector for solution */
  if (check_flag((void*)y, "N_VNew_Serial", 0)) { return 1; }
  NV_Ith_S(y, 0) = u0; /* Set initial conditions */
  NV_Ith_S(y, 1) = v0;
  NV_Ith_S(y, 2) = w0;

  /* Call ARKStepCreate to initialize the ARK timestepper module and
     specify the right-hand side functions in y'=fe(t,y)+fi(t,y),
     the initial time T0, and the initial dependent variable vector y. */
  arkode_mem = ARKStepCreate(fe, fi, T0, y, ctx);
  if (check_flag((void*)arkode_mem, "ARKStepCreate", 0)) { return 1; }

  /* Initialize fixed-point nonlinear solver and attach to ARKODE */
  NLS = SUNNonlinSol_FixedPoint(y, fp_m, ctx);
  if (check_flag((void*)NLS, "SUNNonlinSol_FixedPoint", 0)) { return 1; }

  flag = ARKodeSetNonlinearSolver(arkode_mem, NLS);
  if (check_flag(&flag, "ARKodeSetNonlinearSolver", 1)) { return 1; }

  /* Set routines */
  flag = ARKodeSetUserData(arkode_mem,
                           (void*)rdata); /* Pass rdata to user functions */
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }

  flag = ARKodeSStolerances(arkode_mem, reltol, abstol); /* Specify tolerances */
  if (check_flag(&flag, "ARKodeSStolerances", 1)) { return 1; }

  flag = ARKodeSetMaxNonlinIters(arkode_mem,
                                 maxcor); /* Increase default iterations */
  if (check_flag(&flag, "ARKodeSetMaxNonlinIters", 1)) { return 1; }

  flag = ARKodeSetAutonomous(arkode_mem, SUNTRUE);
  if (check_flag(&flag, "ARKodeSetAutonomous", 1)) { return 1; }

  /* Open output stream for results, output comment line */
  UFID = fopen("solution.txt", "w");
  fprintf(UFID, "# t u v w\n");

  /* output initial condition to disk */
  fprintf(UFID, " %.16" ESYM " %.16" ESYM " %.16" ESYM " %.16" ESYM "\n", T0,
          NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));

  /* Main time-stepping loop: calls ARKodeEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t    = T0;
  tout = T0 + dTout;
  printf("        t           u           v           w\n");
  printf("   ----------------------------------------------\n");
  for (iout = 0; iout < Nt; iout++)
  {
    flag = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_NORMAL); /* call integrator */
    if (check_flag(&flag, "ARKodeEvolve", 1)) { break; }
    printf("  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM
           "\n", /* access/print solution */
           t, NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));
    fprintf(UFID, " %.16" ESYM " %.16" ESYM " %.16" ESYM " %.16" ESYM "\n", t,
            NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));
    if (flag >= 0)
    { /* successful solve: update time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
    else
    { /* unsuccessful solve: break */
      fprintf(stderr, "Solver failure, stopping integration\n");
      break;
    }
  }
  printf("   ----------------------------------------------\n");
  fclose(UFID);

  /* Print some final statistics */
  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  check_flag(&flag, "ARKodeGetNumSteps", 1);
  flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  check_flag(&flag, "ARKodeGetNumStepAttempts", 1);
  flag = ARKodeGetNumRhsEvals(arkode_mem, 0, &nfe);
  check_flag(&flag, "ARKodeGetNumRhsEvals", 1);
  flag = ARKodeGetNumRhsEvals(arkode_mem, 1, &nfi);
  check_flag(&flag, "ARKodeGetNumRhsEvals", 1);
  flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
  check_flag(&flag, "ARKodeGetNumErrTestFails", 1);
  flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
  check_flag(&flag, "ARKodeGetNumNonlinSolvIters", 1);
  flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  check_flag(&flag, "ARKodeGetNumNonlinSolvConvFails", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total number of fixed-point iterations = %li\n", nni);
  printf("   Total number of nonlinear solver convergence failures = %li\n",
         ncfn);
  printf("   Total number of error test failures = %li\n\n", netf);

  /* Clean up and return with successful completion */
  N_VDestroy(y);
  ARKodeFree(&arkode_mem);
  SUNNonlinSolFree(NLS);
  SUNLogger_Destroy(&logger);
  SUNContext_Free(&ctx);

  return 0;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* fi routine to compute the implicit portion of the ODE RHS. */
static int fi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* rdata = (sunrealtype*)user_data; /* cast user_data to sunrealtype */
  sunrealtype b  = rdata[1];                    /* access data entries */
  sunrealtype ep = rdata[2];
  sunrealtype w  = NV_Ith_S(y, 2); /* access solution values */

  /* fill in the RHS function */
  NV_Ith_S(ydot, 0) = 0.0;
  NV_Ith_S(ydot, 1) = 0.0;
  NV_Ith_S(ydot, 2) = (b - w) / ep;

  return 0; /* Return with success */
}

/* fe routine to compute the explicit portion of the ODE RHS. */
static int fe(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* rdata = (sunrealtype*)user_data; /* cast user_data to sunrealtype */
  sunrealtype a = rdata[0];                     /* access data entries */
  sunrealtype u = NV_Ith_S(y, 0);               /* access solution values */
  sunrealtype v = NV_Ith_S(y, 1);
  sunrealtype w = NV_Ith_S(y, 2);

  /* fill in the RHS function */
  NV_Ith_S(ydot, 0) = a - (w + 1.0) * u + v * u * u;
  NV_Ith_S(ydot, 1) = w * u - v * u * u;
  NV_Ith_S(ydot, 2) = -w * u;

  return 0; /* Return with success */
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
static int check_flag(void* flagvalue, const char* funcname, int opt)
{
  int* errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  /* Check if flag < 0 */
  else if (opt == 1)
  {
    errflag = (int*)flagvalue;
    if (*errflag < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return 1;
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}

/*---- end of file ----*/
