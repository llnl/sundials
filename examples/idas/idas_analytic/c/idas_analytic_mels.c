/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
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
 * The following is a simple example problem with analytical
 * solution adapted from example 10.2 of Ascher & Petzold, "Computer
 * Methods for Ordinary Differential Equations and
 * Differential-Algebraic Equations," SIAM, 1998, page 267:
 *    x1'(t) = (1-alpha)/(t-2)*x1 - x1 + (alpha-1)*x2 + 2*exp(t)
 *         0 = (t+2)*x1 - (t+2)*exp(t)
 * for t in the interval [0.0, 1.0], with initial condition:
 *    x1(0) = 1   and   x2(0) = -1/2.
 * The problem has true solution
 *    x1(t) = exp(t)  and  x2(t) = exp(t)/(t-2)
 *
 * This program solves the problem with IDAS using a custom
 * 'matrix-embedded' SUNLinearSolver. Output is printed
 * every 0.1 units of time (10 total).  Run statistics (optional
 * outputs) are printed at the end.
 * -----------------------------------------------------------------*/

#include <idas/idas.h>              /* prototypes for IDAS fcts., consts.   */
#include <nvector/nvector_serial.h> /* access to serial N_Vector            */
#include <stdio.h>
#include <sundials/sundials_math.h> /* defs. of SUNRabs, SUNRexp, etc.      */
#include <sundials/sundials_types.h> /* defs. of sunrealtype, sunindextype      */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* User-supplied functions called by IDA */
int fres(sunrealtype tres, N_Vector yy, N_Vector yp, N_Vector resval,
         void* user_data);

/* Custom linear solver data structure, accessor macros, and routines */
static SUNLinearSolver MatrixEmbeddedLS(void* ida_mem, SUNContext ctx);
static SUNLinearSolver_Type MatrixEmbeddedLSType(SUNLinearSolver S);
static int MatrixEmbeddedLSSolve(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                                 N_Vector b, sunrealtype tol);
static SUNErrCode MatrixEmbeddedLSFree(SUNLinearSolver S);

/* Private function to check function return values */
static int check_retval(void* returnvalue, const char* funcname, int opt);

/* Private function to check computed solution */
static void analytical_solution(sunrealtype t, N_Vector y, N_Vector yp);
static int check_ans(N_Vector y, sunrealtype t, sunrealtype rtol,
                     sunrealtype atol);

/* Main Program */
int main(int argc, char* argv[])
{
  /* SUNDIALS context object */
  SUNContext ctx;

  /* general problem parameters */
  sunrealtype T0     = SUN_RCONST(0.0);    /* initial time */
  sunrealtype Tf     = SUN_RCONST(1.0);    /* final time */
  sunrealtype dTout  = SUN_RCONST(0.1);    /* time between outputs */
  sunindextype NEQ   = 2;                  /* number of dependent vars. */
  sunrealtype reltol = SUN_RCONST(1.0e-4); /* tolerances */
  sunrealtype abstol = SUN_RCONST(1.0e-9);
  sunrealtype alpha  = SUN_RCONST(10.0); /* stiffness parameter */

  /* general problem variables */
  int retval;                /* reusable error-checking flag */
  N_Vector yy        = NULL; /* empty vector for storing solution */
  N_Vector yp        = NULL; /* empty vector for storing solution derivative */
  SUNLinearSolver LS = NULL; /* empty linear solver object */
  void* ida_mem      = NULL; /* empty IDA memory structure */
  sunrealtype t, tout, h0;
  long int nst, nre, nni, netf, ncfn, nreLS;

  /* Initial diagnostics output */
  printf("\nAnalytical DAE test problem:\n");
  printf("    alpha = %" GSYM "\n", alpha);
  printf("   reltol = %.1" ESYM "\n", reltol);
  printf("   abstol = %.1" ESYM "\n\n", abstol);

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return 1; }

  /* Initialize data structures */
  yy = N_VNew_Serial(NEQ, ctx); /* Create serial vector for solution */
  if (check_retval((void*)yy, "N_VNew_Serial", 0)) { return 1; }
  yp = N_VClone(yy); /* Create serial vector for solution derivative */
  if (check_retval((void*)yp, "N_VClone", 0)) { return 1; }
  analytical_solution(T0, yy, yp); /* Specify initial conditions */

  /* Call IDACreate and IDAInit to initialize IDA memory */
  ida_mem = IDACreate(ctx);
  if (check_retval((void*)ida_mem, "IDACreate", 0)) { return (1); }
  retval = IDAInit(ida_mem, fres, T0, yy, yp);
  if (check_retval(&retval, "IDAInit", 1)) { return (1); }

  /* Set routines */
  retval = IDASetUserData(ida_mem, (void*)&alpha);
  if (check_retval(&retval, "IDASetUserData", 1)) { return (1); }
  retval = IDASStolerances(ida_mem, reltol, abstol);
  if (check_retval(&retval, "IDASStolerances", 1)) { return (1); }

  /* Create custom matrix-embedded linear solver */
  LS = MatrixEmbeddedLS(ida_mem, ctx);
  if (check_retval((void*)LS, "MatrixEmbeddedLS", 0)) { return 1; }

  /* Attach the linear solver */
  retval = IDASetLinearSolver(ida_mem, LS, NULL);
  if (check_retval(&retval, "IDASetLinearSolver", 1)) { return (1); }

  /* Override any current settings with command-line options */
  retval = IDASetOptions(ida_mem, NULL, NULL, argc, argv);
  if (check_retval(&retval, "IDASetOptions", 1)) { return 1; }

  /* In loop, call IDASolve, print results, and test for error.
     Stops when the final time has been reached. */
  t    = T0;
  tout = T0 + dTout;
  printf("        t          x1         x2\n");
  printf("   ----------------------------------\n");
  while (Tf - t > 1.0e-15)
  {
    retval = IDASolve(ida_mem, tout, &t, yy, yp, IDA_NORMAL); /* call integrator */
    if (check_retval(&retval, "IDASolve", 1)) { return (1); }
    printf("  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM "\n", t,
           NV_Ith_S(yy, 0), NV_Ith_S(yy, 1)); /* access/print solution */
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
  printf("   ----------------------------------\n");

  /* Get/print some final statistics on how the solve progressed */
  retval = IDAGetActualInitStep(ida_mem, &h0);
  check_retval(&retval, "IDAGetActualInitStep", 1);
  retval = IDAGetNumSteps(ida_mem, &nst);
  check_retval(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetNumResEvals(ida_mem, &nre);
  check_retval(&retval, "IDAGetNumResEvals", 1);
  retval = IDAGetNumNonlinSolvIters(ida_mem, &nni);
  check_retval(&retval, "IDAGetNumNonlinSolvIters", 1);
  retval = IDAGetNumErrTestFails(ida_mem, &netf);
  check_retval(&retval, "IDAGetNumErrTestFails", 1);
  retval = IDAGetNumNonlinSolvConvFails(ida_mem, &ncfn);
  check_retval(&retval, "IDAGetNumNonlinSolvConvFails", 1);
  retval = IDAGetNumLinResEvals(ida_mem, &nreLS);
  check_retval(&retval, "IDAGetNumLinResEvals", 1);

  printf("\nFinal Solver Statistics: \n\n");
  printf("Initial time step                  = %8.6" FSYM "\n", h0);
  printf("Number of steps                    = %ld\n", nst);
  printf("Number of residual evaluations     = %ld\n", nre + nreLS);
  printf("Number of nonlinear iterations     = %ld\n", nni);
  printf("Number of error test failures      = %ld\n", netf);
  printf("Number of nonlinear conv. failures = %ld\n", ncfn);

  /* check the solution error */
  retval = check_ans(yy, t, reltol, abstol);

  /* Clean up and return */
  IDAFree(&ida_mem);
  SUNLinSolFree(LS);
  N_VDestroy(yy);
  N_VDestroy(yp);
  SUNContext_Free(&ctx);

  return (retval);
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* System residual function:
      0 = (1-alpha)/(t-2)*x1 - x1 + (alpha-1)*x2 + 2*exp(t) - x1'(t)
      0 = (t+2)*x1 - (t+2)*exp(t)
*/
int fres(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data)
{
  sunrealtype* rdata = (sunrealtype*)user_data; /* cast user_data to sunrealtype */
  sunrealtype alpha = rdata[0]; /* set shortcut for stiffness parameter */
  sunrealtype x1    = NV_Ith_S(yy, 0); /* access current solution values */
  sunrealtype x2    = NV_Ith_S(yy, 1);
  sunrealtype x1p   = NV_Ith_S(yp, 0); /* access current derivative values */
  sunrealtype ONE   = SUN_RCONST(1.0);
  sunrealtype TWO   = SUN_RCONST(2.0);

  NV_Ith_S(rr, 0) = (ONE - alpha) / (t - TWO) * x1 - x1 + (alpha - ONE) * x2 +
                    TWO * exp(t) - x1p;
  NV_Ith_S(rr, 1) = (t + TWO) * x1 - (t + TWO) * SUNRexp(t);

  return (0);
}

/*-------------------------------------
 * Custom matrix-embedded linear solver
 *-------------------------------------*/

/* constructor */
static SUNLinearSolver MatrixEmbeddedLS(void* ida_mem, SUNContext ctx)
{
  /* Create an empty linear solver */
  SUNLinearSolver LS = SUNLinSolNewEmpty(ctx);
  if (LS == NULL) { return NULL; }

  /* Attach operations */
  LS->ops->gettype = MatrixEmbeddedLSType;
  LS->ops->solve   = MatrixEmbeddedLSSolve;
  LS->ops->free    = MatrixEmbeddedLSFree;

  /* Set content pointer to IDA memory */
  LS->content = ida_mem;

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
  N_Vector yypred, yppred, yyn, ypn, res;
  sunrealtype tcur, cj;
  void* user_data;
  sunrealtype* rdata;
  sunrealtype alpha;
  sunrealtype a11, a12, a21, b1, b2;
  sunrealtype ONE = SUN_RCONST(1.0);
  sunrealtype TWO = SUN_RCONST(2.0);

  /* retrieve implicit system data from IDA */
  retval = IDAGetNonlinearSystemData(LS->content, &tcur, &yypred, &yppred, &yyn,
                                     &ypn, &res, &cj, &user_data);
  if (check_retval((void*)&retval, "IDAGetNonlinearSystemData", 1))
  {
    return (-1);
  }

  /* extract stiffness parameter from user_data */
  rdata = (sunrealtype*)user_data;
  alpha = rdata[0];

  /* perform linear solve: A*x=b
         A = df/dy + cj*df/dyp
      =>
         A = [ - cj - (alpha - 1)/(t - 2) - 1, alpha - 1]
             [                          t + 2,         0]

   */
  a11            = -cj - (alpha - ONE) / (tcur - TWO) - ONE;
  a12            = alpha - ONE;
  a21            = tcur + TWO;
  b1             = NV_Ith_S(b, 0);
  b2             = NV_Ith_S(b, 1);
  NV_Ith_S(x, 0) = b2 / a21;
  NV_Ith_S(x, 1) = -(a11 * b2 - a21 * b1) / (a12 * a21);

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

/* routine to fill analytical solution and its derivative */
static void analytical_solution(sunrealtype t, N_Vector y, N_Vector yp)
{
  NV_Ith_S(y, 0)  = SUNRexp(t);
  NV_Ith_S(y, 1)  = SUNRexp(t) / (t - SUN_RCONST(2.0));
  NV_Ith_S(yp, 0) = SUNRexp(t);
  NV_Ith_S(yp, 1) = SUNRexp(t) / (t - SUN_RCONST(2.0)) -
                    SUNRexp(t) / (t - SUN_RCONST(2.0)) / (t - SUN_RCONST(2.0));
}

/* check the computed solution */
static int check_ans(N_Vector y, sunrealtype t, sunrealtype rtol, sunrealtype atol)
{
  int passfail = 0; /* answer pass (0) or fail (1) retval */
  N_Vector ytrue;   /* true solution vector               */
  N_Vector ewt;     /* error weight vector                */
  N_Vector abstol;  /* absolute tolerance vector          */
  sunrealtype err;  /* wrms error                         */
  sunrealtype ONE = SUN_RCONST(1.0);

  /* create solution and error weight vectors */
  ytrue  = N_VClone(y);
  ewt    = N_VClone(y);
  abstol = N_VClone(y);

  /* set the solution data */
  analytical_solution(t, ytrue, abstol);

  /* compute the error weight vector, loosen atol */
  N_VConst(atol, abstol);
  N_VAbs(ytrue, ewt);
  N_VLinearSum(rtol, ewt, SUN_RCONST(10.0), abstol, ewt);
  if (N_VMin(ewt) <= SUN_RCONST(0.0))
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: check_ans failed - ewt <= 0\n\n");
    return (-1);
  }
  N_VInv(ewt, ewt);

  /* compute the solution error */
  N_VLinearSum(ONE, y, -ONE, ytrue, ytrue);
  err = N_VWrmsNorm(ytrue, ewt);

  /* is the solution within the tolerances? */
  passfail = (err < ONE) ? 0 : 1;

  if (passfail)
  {
    fprintf(stdout, "\nSUNDIALS_WARNING: check_ans error=%" GSYM "\n\n", err);
  }

  /* Free vectors */
  N_VDestroy(ytrue);
  N_VDestroy(abstol);
  N_VDestroy(ewt);

  return (passfail);
}

/*---- end of file ----*/
