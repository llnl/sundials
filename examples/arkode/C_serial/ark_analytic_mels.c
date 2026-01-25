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
 *    dy/dt = lambda*y + 1/(1+t^2) - lambda*atan(t)
 * for t in the interval [0.0, 10.0], with initial condition: y=0.
 *
 * The stiffness of the problem is directly proportional to the
 * value of "lambda".  The value of lambda should be negative to
 * result in a well-posed ODE; for values with magnitude larger
 * than 100 the problem becomes quite stiff.
 *
 * This program solves the problem with the DIRK method and a
 * custom 'matrix-embedded' SUNLinearSolver. Output is printed
 * every 1.0 units of time (10 total).  Run statistics (optional
 * outputs) are printed at the end.
 *-----------------------------------------------------------------*/

/* Header files */
#include <arkode/arkode_arkstep.h> /* prototypes for ARKStep fcts., consts */
#include <math.h>
#include <nvector/nvector_serial.h> /* serial N_Vector types, fcts., macros */
#include <stdio.h>
#include <sundials/sundials_types.h> /* definition of type sunrealtype          */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* User-supplied functions called by ARKODE */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/* Custom linear solver data structure, accessor macros, and routines */
static SUNLinearSolver MatrixEmbeddedLS(void* arkode_mem, SUNContext ctx);
static SUNLinearSolver_Type MatrixEmbeddedLSType(SUNLinearSolver S);
static int MatrixEmbeddedLSSolve(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                                 N_Vector b, sunrealtype tol);
static SUNErrCode MatrixEmbeddedLSFree(SUNLinearSolver S);

/* Private function to check function return values */
static int check_retval(void* returnvalue, const char* funcname, int opt);

/* Private function to check computed solution */
static int check_ans(N_Vector y, sunrealtype t, sunrealtype rtol,
                     sunrealtype atol);

/* Main Program */
int main(void)
{
  /* general problem parameters */
  sunrealtype T0     = SUN_RCONST(0.0);    /* initial time */
  sunrealtype Tf     = SUN_RCONST(10.0);   /* final time */
  sunrealtype dTout  = SUN_RCONST(1.0);    /* time between outputs */
  sunindextype NEQ   = 1;                  /* number of dependent vars. */
  sunrealtype reltol = SUN_RCONST(1.0e-5); /* tolerances */
  sunrealtype abstol = SUN_RCONST(1.0e-10);
  sunrealtype lambda = SUN_RCONST(-100.0); /* stiffness parameter */

  /* general problem variables */
  int retval;                /* reusable error-checking flag */
  N_Vector y         = NULL; /* empty vector for storing solution */
  SUNLinearSolver LS = NULL; /* empty linear solver object */
  void* arkode_mem   = NULL; /* empty ARKODE memory structure */
  sunrealtype t, tout;
  long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;

  /* Create the SUNDIALS context object for this simulation */
  SUNContext ctx;
  retval = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return 1; }

  /* Initial diagnostics output */
  printf("\nAnalytical ODE test problem:\n");
  printf("   lambda = %" GSYM "\n", lambda);
  printf("   reltol = %.1" ESYM "\n", reltol);
  printf("   abstol = %.1" ESYM "\n\n", abstol);

  /* Initialize data structures */
  y = N_VNew_Serial(NEQ, ctx); /* Create serial vector for solution */
  if (check_retval((void*)y, "N_VNew_Serial", 0)) { return 1; }
  N_VConst(SUN_RCONST(0.0), y); /* Specify initial condition */

  /* Call ARKStepCreate to initialize the ARK timestepper module and
     specify the right-hand side function in y'=f(t,y), the initial time
     T0, and the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. */
  arkode_mem = ARKStepCreate(NULL, f, T0, y, ctx);
  if (check_retval((void*)arkode_mem, "ARKStepCreate", 0)) { return 1; }

  /* Set routines */
  retval = ARKodeSetUserData(arkode_mem,
                             (void*)&lambda); /* Pass lambda to user functions */
  if (check_retval(&retval, "ARKodeSetUserData", 1)) { return 1; }
  retval = ARKodeSStolerances(arkode_mem, reltol, abstol); /* Specify tolerances */
  if (check_retval(&retval, "ARKodeSStolerances", 1)) { return 1; }

  /* Initialize custom matrix-embedded linear solver */
  LS = MatrixEmbeddedLS(arkode_mem, ctx);
  if (check_retval((void*)LS, "MatrixEmbeddedLS", 0)) { return 1; }
  retval = ARKodeSetLinearSolver(arkode_mem, LS, NULL); /* Attach linear solver */
  if (check_retval(&retval, "ARKodeSetLinearSolver", 1)) { return 1; }

  /* Specify linearly implicit RHS, with non-time-dependent Jacobian */
  retval = ARKodeSetLinear(arkode_mem, 0);
  if (check_retval(&retval, "ARKodeSetLinear", 1)) { return 1; }

  /* Main time-stepping loop: calls ARKodeEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached. */
  t    = T0;
  tout = T0 + dTout;
  printf("        t           u\n");
  printf("   ---------------------\n");
  while (Tf - t > 1.0e-15)
  {
    retval = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_NORMAL); /* call integrator */
    if (check_retval(&retval, "ARKodeEvolve", 1)) { break; }
    printf("  %10.6" FSYM "  %10.6" FSYM "\n", t,
           NV_Ith_S(y, 0)); /* access/print solution */
    if (retval >= 0)
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
  printf("   ---------------------\n");

  /* Get/print some final statistics on how the solve progressed */
  retval = ARKodeGetNumSteps(arkode_mem, &nst);
  check_retval(&retval, "ARKodeGetNumSteps", 1);
  retval = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  check_retval(&retval, "ARKodeGetNumStepAttempts", 1);
  retval = ARKodeGetNumRhsEvals(arkode_mem, 0, &nfe);
  check_retval(&retval, "ARKodeGetNumRhsEvals", 1);
  retval = ARKodeGetNumRhsEvals(arkode_mem, 1, &nfi);
  check_retval(&retval, "ARKodeGetNumRhsEvals", 1);
  retval = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
  check_retval(&retval, "ARKodeGetNumLinSolvSetups", 1);
  retval = ARKodeGetNumErrTestFails(arkode_mem, &netf);
  check_retval(&retval, "ARKodeGetNumErrTestFails", 1);
  retval = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
  check_retval(&retval, "ARKodeGetNumNonlinSolvIters", 1);
  retval = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  check_retval(&retval, "ARKodeGetNumNonlinSolvConvFails", 1);
  retval = ARKodeGetNumJacEvals(arkode_mem, &nje);
  check_retval(&retval, "ARKodeGetNumJacEvals", 1);
  retval = ARKodeGetNumLinRhsEvals(arkode_mem, &nfeLS);
  check_retval(&retval, "ARKodeGetNumLinRhsEvals", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total number of linear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n\n", netf);

  /* check the solution error */
  retval = check_ans(y, t, reltol, abstol);

  /* Clean up and return */
  N_VDestroy(y);           /* Free y vector */
  ARKodeFree(&arkode_mem); /* Free integrator memory */
  SUNLinSolFree(LS);       /* Free linear solver */
  SUNContext_Free(&ctx);   /* Free the SUNContext */

  return retval;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* rdata = (sunrealtype*)user_data; /* cast user_data to sunrealtype */
  sunrealtype lambda = rdata[0]; /* set shortcut for stiffness parameter */
  sunrealtype u      = NV_Ith_S(y, 0); /* access current solution value */

  /* fill in the RHS function: "NV_Ith_S" accesses the 0th entry of ydot */
  NV_Ith_S(ydot, 0) = lambda * u + SUN_RCONST(1.0) / (SUN_RCONST(1.0) + t * t) -
                      lambda * atan(t);

  return 0; /* return with success */
}

/*-------------------------------------
 * Custom matrix-embedded linear solver
 *-------------------------------------*/

/* constructor */
static SUNLinearSolver MatrixEmbeddedLS(void* arkode_mem, SUNContext ctx)
{
  /* Create an empty linear solver */
  SUNLinearSolver LS = SUNLinSolNewEmpty(ctx);
  if (LS == NULL) { return NULL; }

  /* Attach operations */
  LS->ops->gettype = MatrixEmbeddedLSType;
  LS->ops->solve   = MatrixEmbeddedLSSolve;
  LS->ops->free    = MatrixEmbeddedLSFree;

  /* Set content pointer to ARKODE memory */
  LS->content = arkode_mem;

  /* Return solver */
  return (LS);
}

/* type descriptor */
static SUNLinearSolver_Type MatrixEmbeddedLSType(SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_MATRIX_EMBEDDED);
}

/* linear solve routine */
static int MatrixEmbeddedLSSolve(SUNLinearSolver LS, SUNMatrix A, N_Vector x,
                                 N_Vector b, sunrealtype tol)
{
  /* temporary variables */
  int retval;
  N_Vector z, zpred, Fi, sdata;
  sunrealtype tcur, gamma;
  void* user_data;
  sunrealtype* rdata;
  sunrealtype lambda;

  /* retrieve implicit system data from ARKODE */
  retval = ARKodeGetNonlinearSystemData(LS->content, &tcur, &zpred, &z, &Fi,
                                        &gamma, &sdata, &user_data);
  if (check_retval((void*)&retval, "ARKodeGetNonlinearSystemData", 1))
  {
    return (-1);
  }

  /* extract stiffness parameter from user_data */
  rdata  = (sunrealtype*)user_data;
  lambda = rdata[0];

  /* perform linear solve: (1-gamma*lambda)*x = b */
  NV_Ith_S(x, 0) = NV_Ith_S(b, 0) / (1 - gamma * lambda);

  /* return with success */
  return (SUN_SUCCESS);
}

/* destructor */
static SUNErrCode MatrixEmbeddedLSFree(SUNLinearSolver LS)
{
  if (LS == NULL) { return (SUN_SUCCESS); }
  LS->content = NULL;
  SUNLinSolFreeEmpty(LS);
  return (SUN_SUCCESS);
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
static int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  /* Check if flag < 0 */
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *retval);
      return 1;
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}

/* check the computed solution */
static int check_ans(N_Vector y, sunrealtype t, sunrealtype rtol, sunrealtype atol)
{
  int passfail = 0;          /* answer pass (0) or fail (1) flag     */
  sunrealtype ans, err, ewt; /* answer data, error, and error weight */

  /* compute solution error */
  ans = atan(t);
  ewt = SUN_RCONST(1.0) / (rtol * SUNRabs(ans) + atol);
  err = ewt * SUNRabs(NV_Ith_S(y, 0) - ans);

  /* The local errors accumulate from step to step so that the global error is
   * not quite within the local error tolerances. This factor accounts for
   * this. */
  sunrealtype global_bound = SUN_RCONST(1.5);

  /* is the solution within the tolerances? */
  passfail = (err < global_bound) ? 0 : 1;

  if (passfail)
  {
    fprintf(stdout, "\nSUNDIALS_WARNING: check_ans error=%" GSYM "\n\n", err);
  }

  return (passfail);
}

/*---- end of file ----*/
