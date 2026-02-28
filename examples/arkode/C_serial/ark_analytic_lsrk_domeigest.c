/*-----------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ SMU
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
 * The following is a simple example problem that solves the ODE
 *
 *    dy/dt = (lambda - alpha*cos((10 - t)/10*pi)*y + 1/(1+t^2)
 *          - (lambda - alpha*cos((10 - t)/10*pi)*atan(t),
 *
 * for t in the interval [0.0, 10.0], with an initial condition: y=0.
 * This initial value problem has the analytical solution
 *
 *    y(t) = arctan(t).
 *
 * The stiffness of the problem depends on both lambda and alpha together.
 * While lambda determines the center of the stiffness parameter,
 * the value of alpha determines the radius of the interval in which
 * the stiffness parameter lies.
 *
 * The value of lambda - alpha*cos((10 - t)/10*pi) should
 * be negative to result in a well-posed ODE; for values with magnitude
 * larger than 100 the problem becomes quite stiff.
 *
 * This program solves the problem with the LSRK method using internal
 * SUNDIALS dominant eigenvalue estimation (DEE) module.
 * Output is printed every 1.0 units of time (10 total).
 * Run statistics (optional outputs) are printed at the end.
 *-----------------------------------------------------------------*/

/* Header files */
#include <arkode/arkode_lsrkstep.h> /* prototypes for LSRKStep fcts., consts */
#include <math.h>
#include <nvector/nvector_serial.h> /* serial N_Vector types, fcts., macros */
#include <stdio.h>
#include <sundials/sundials_math.h> /* def. of SUNRsqrt, etc. */
#include <sundials/sundials_types.h> /* definition of type sunrealtype          */
#include <sundomeigest/sundomeigest_power.h> /* access to Power Iteration module */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

#if defined(SUNDIALS_DOUBLE_PRECISION)
#define ATAN(x) (atan((x)))
#define ACOS(x) (acos((x)))
#define COS(x)  (cos((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define ATAN(x) (atanf((x)))
#define ACOS(x) (acosf((x)))
#define COS(x)  (cosf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define ATAN(x) (atanl((x)))
#define ACOS(x) (acosl((x)))
#define COS(x)  (cosl((x)))
#endif

/* User-supplied Functions Called by the Solver */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/* Private function to check function return values */
static int check_flag(void* flagvalue, const char* funcname, int opt);

/* Private function to check computed solution */
static int check_ans(N_Vector y, sunrealtype t, sunrealtype rtol,
                     sunrealtype atol);

/* Private function to compute error */
static int compute_error(N_Vector y, sunrealtype t);

/* Main Program */
int main(int argc, char* argv[])
{
  /* general problem parameters */
  sunrealtype T0    = SUN_RCONST(0.0);  /* initial time */
  sunrealtype Tf    = SUN_RCONST(10.0); /* final time */
  sunrealtype dTout = SUN_RCONST(1.0);  /* time between outputs */
  sunindextype NEQ  = 1;                /* number of dependent vars. */

#if defined(SUNDIALS_DOUBLE_PRECISION)
  sunrealtype reltol = SUN_RCONST(1.0e-8); /* tolerances */
  sunrealtype abstol = SUN_RCONST(1.0e-8);
  sunrealtype lambda = SUN_RCONST(-1.0e+6); /* stiffness parameter 1 */
  sunrealtype alpha  = SUN_RCONST(1.0e+2);  /* stiffness parameter 2 */
#elif defined(SUNDIALS_SINGLE_PRECISION)
  sunrealtype reltol = SUN_RCONST(1.0e-4); /* tolerances */
  sunrealtype abstol = SUN_RCONST(1.0e-8);
  sunrealtype lambda = SUN_RCONST(-1.0e+3); /* stiffness parameter 1 */
  sunrealtype alpha  = SUN_RCONST(1.0e+1);  /* stiffness parameter 2 */
#elif defined(SUNDIALS_EXTENDED_PRECISION)
  sunrealtype reltol = SUN_RCONST(1.0e-8); /* tolerances */
  sunrealtype abstol = SUN_RCONST(1.0e-8);
  sunrealtype lambda = SUN_RCONST(-1.0e+6); /* stiffness parameter 1 */
  sunrealtype alpha  = SUN_RCONST(1.0e+2);  /* stiffness parameter 2 */
#endif

  sunrealtype UserData[2];
  UserData[0] = lambda;
  UserData[1] = alpha;

  /* general problem variables */
  int flag;                /* reusable error-checking flag */
  N_Vector y       = NULL; /* empty vector for storing solution */
  void* arkode_mem = NULL; /* empty ARKode memory structure */
  FILE *UFID, *FID;
  sunrealtype t, tout;

  /* Dominant Eigenvalue Estimator (DEE) pointers and variables */
  SUNDomEigEstimator DEE = NULL; /* domeig estimator object */
  sunindextype max_iters = 100;  /* max number of power iterations (PI)*/
  sunindextype numwarmup = 10;   /* number of preprocessing warmups */
  sunrealtype rel_tol    = SUN_RCONST(5.0e-3); /* relative error for PI*/
  N_Vector q             = NULL;               /* random initial eigenvector */

  /* Create the SUNDIALS context object for this simulation */
  SUNContext ctx;
  flag = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_flag(&flag, "SUNContext_Create", 1)) { return 1; }

  /* Initial diagnostics output */
  printf("\nAnalytical ODE test problem with a variable Jacobian:");
  printf("\nThe stiffness of the problem is directly proportional to");
  printf("\n\"lambda - alpha*cos((10 - t)/10*pi)\"\n\n");
  printf("    lambda = %" GSYM "\n", lambda);
  printf("     alpha = %" GSYM "\n", alpha);
  printf("    reltol = %.1" ESYM "\n", reltol);
  printf("    abstol = %.1" ESYM "\n\n", abstol);

  /* Initialize data structures */
  y = N_VNew_Serial(NEQ, ctx); /* Create serial vector for solution */
  if (check_flag((void*)y, "N_VNew_Serial", 0)) { return 1; }
  N_VConst(SUN_RCONST(0.0), y); /* Specify initial condition */

  /* Call LSRKStepCreateSTS to initialize the ARK timestepper module and
     specify the right-hand side function in y'=f(t,y), the initial time
     T0, and the initial dependent variable vector y. */
  arkode_mem = LSRKStepCreateSTS(f, T0, y, ctx);
  if (check_flag((void*)arkode_mem, "LSRKStepCreateSTS", 0)) { return 1; }

  /* Set routines */
  flag = ARKodeSetUserData(arkode_mem,
                           (void*)&UserData); /* Pass lambda to user functions */
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }

  /* Specify tolerances */
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) { return 1; }

  /* Set the initial random eigenvector for the DEE */
  q = N_VClone(y);
  if (check_flag(q, "N_VClone", 0)) { return 1; }

  sunrealtype* qd = N_VGetArrayPointer(q);
  for (int i = 0; i < NEQ; i++)
  {
    qd[i] = (sunrealtype)rand() / (sunrealtype)RAND_MAX;
  }

  /* Create power iteration dominant eigenvalue estimator (DEE) */
  DEE = SUNDomEigEstimator_Power(q, max_iters, rel_tol, ctx);
  if (check_flag(DEE, "SUNDomEigEstimator_Power", 0)) { return 1; }

  /* After the DEE creation, random q vector is no longer needed.
     It is used only to initialize the DEE */
  N_VDestroy(q);

  /* Attach the DEE to the LSRKStep module.
  There is no need to set Atimes or initialize since these are all
  performed after attaching the DEE by LSRKStep. */
  flag = LSRKStepSetDomEigEstimator(arkode_mem, DEE);
  if (check_flag(&flag, "LSRKStepSetDomEigEstimator", 1)) { return 1; }

  /* Set the number of preprocessing warmups. The warmup
  is used to compute a "better" initial eigenvector and so an initial
  eigenvalue. The warmup is performed only once by the LSRKStep module
  internally unless LSRKStepSetNumDomEigEstPreprocessIters is called to set
  a new number of succeeding warmups that would be executed before
  every dominant eigenvalue estimate calls */
  flag = LSRKStepSetNumDomEigEstInitPreprocessIters(arkode_mem, numwarmup);
  if (check_flag(&flag, "LSRKStepSetNumDomEigEstInitPreprocessIters", 1))
  {
    return 1;
  }

  /* Specify after how many successful steps dom_eig is recomputed */
  flag = LSRKStepSetDomEigFrequency(arkode_mem, 25);
  if (check_flag(&flag, "LSRKStepSetDomEigFrequency", 1)) { return 1; }

  /* Specify max number of stages allowed */
  flag = LSRKStepSetMaxNumStages(arkode_mem, 200);
  if (check_flag(&flag, "LSRKStepSetMaxNumStages", 1)) { return 1; }

  /* Specify max number of steps allowed */
  flag = ARKodeSetMaxNumSteps(arkode_mem, 2000);
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) { return 1; }

  /* Specify safety factor for user provided dom_eig */
  flag = LSRKStepSetDomEigSafetyFactor(arkode_mem, SUN_RCONST(1.01));
  if (check_flag(&flag, "LSRKStepSetDomEigSafetyFactor", 1)) { return 1; }

  /* Specify the number of preprocessing warmups before each estimate call
     succeeding the very first estimate call. */
  flag = LSRKStepSetNumDomEigEstPreprocessIters(arkode_mem, 0);
  if (check_flag(&flag, "LSRKStepSetNumDomEigEstPreprocessIters", 1))
  {
    return 1;
  }

  /* Specify the Runge--Kutta--Chebyshev LSRK method by name */
  flag = LSRKStepSetSTSMethodByName(arkode_mem, "ARKODE_LSRK_RKC_2");
  if (check_flag(&flag, "LSRKStepSetSTSMethodByName", 1)) { return 1; }

  /* Override any current settings with command-line options */
  flag = SUNDomEigEstimator_SetOptions(DEE, NULL, NULL, argc, argv);
  if (check_flag(&flag, "SUNDomEigEstimator_SetOptions", 1)) { return 1; }

  /* Override any current settings with command-line options */
  flag = ARKodeSetOptions(arkode_mem, NULL, NULL, argc, argv);
  if (check_flag(&flag, "ARKodeSetOptions", 1)) { return 1; }

  /* Open output stream for results, output comment line */
  UFID = fopen("solution.txt", "w");
  fprintf(UFID, "# t u\n");

  /* output initial condition to disk */
  fprintf(UFID, " %.16" ESYM " %.16" ESYM "\n", T0, N_VGetArrayPointer(y)[0]);

  /* Main time-stepping loop: calls ARKodeEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t    = T0;
  tout = T0 + dTout;
  printf("        t           u\n");
  printf("   ---------------------\n");
  while (Tf - t > SUN_RCONST(1.0e-15))
  {
    flag = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_NORMAL); /* call integrator */
    if (check_flag(&flag, "ARKodeEvolve", 1)) { break; }
    printf("  %10.6" FSYM "  %10.6" FSYM "\n", t,
           N_VGetArrayPointer(y)[0]); /* access/print solution */
    fprintf(UFID, " %.16" ESYM " %.16" ESYM "\n", t, N_VGetArrayPointer(y)[0]);
    if (flag < 0)
    { /* unsuccessful solve: break */
      fprintf(stderr, "Solver failure, stopping integration\n");
      break;
    }
    else
    { /* successful solve: update time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
  }
  printf("   ---------------------\n");
  fclose(UFID);

  /* Print final statistics */
  printf("\nFinal Statistics:\n");
  flag = ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  /* Print final statistics to a file in CSV format */
  FID  = fopen("ark_analytic_nonlin_stats.csv", "w");
  flag = ARKodePrintAllStats(arkode_mem, FID, SUN_OUTPUTFORMAT_CSV);
  fclose(FID);

  /* check the solution error */
  flag = check_ans(y, t, reltol, abstol);
  flag = compute_error(y, t);

  /* Clean up and return */
  N_VDestroy(y);                    /* Free y vector */
  ARKodeFree(&arkode_mem);          /* Free integrator memory */
  SUNDomEigEstimator_Destroy(&DEE); /* Free DEE object */
  SUNContext_Free(&ctx);            /* Free context */

  return flag;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* rdata = (sunrealtype*)user_data; /* cast user_data to sunrealtype */
  sunrealtype lambda = rdata[0]; /* set shortcut for stiffness parameter 1 */
  sunrealtype alpha  = rdata[1]; /* set shortcut for stiffness parameter 2 */
  sunrealtype u = N_VGetArrayPointer(y)[0]; /* access current solution value */

  /* fill in the RHS function: "N_VGetArrayPointer" accesses the 0th entry of ydot */
  N_VGetArrayPointer(ydot)[0] =
    (lambda - alpha * COS((SUN_RCONST(10.0) - t) / SUN_RCONST(10.0) *
                          ACOS(SUN_RCONST(-1.0)))) *
      u +
    SUN_RCONST(1.0) / (SUN_RCONST(1.0) + t * t) -
    (lambda - alpha * COS((SUN_RCONST(10.0) - t) / SUN_RCONST(10.0) *
                          ACOS(SUN_RCONST(-1.0)))) *
      ATAN(t);

  return 0; /* return with success */
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

/* check the computed solution */
static int check_ans(N_Vector y, sunrealtype t, sunrealtype rtol, sunrealtype atol)
{
  int passfail = 0;          /* answer pass (0) or fail (1) flag     */
  sunrealtype ans, err, ewt; /* answer data, error, and error weight */

  /* compute solution error */
  ans = ATAN(t);
  ewt = SUN_RCONST(1.0) / (rtol * SUNRabs(ans) + atol);
  err = ewt * SUNRabs(N_VGetArrayPointer(y)[0] - ans);

  /* is the solution within the tolerances? */
  passfail = (err < SUN_RCONST(1.0)) ? 0 : 1;

  if (passfail)
  {
    fprintf(stdout, "\nSUNDIALS_WARNING: check_ans error=%" GSYM "\n\n", err);
  }

  return (passfail);
}

/* check the error */
static int compute_error(N_Vector y, sunrealtype t)
{
  sunrealtype ans, err; /* answer data, error */

  /* compute solution error */
  ans = ATAN(t);
  err = SUNRabs(N_VGetArrayPointer(y)[0] - ans);

  fprintf(stdout, "\nACCURACY at the final time   = %" GSYM "\n", err);
  return 0;
}

/*---- end of file ----*/
