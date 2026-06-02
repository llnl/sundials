/*-----------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ SMU
 * Based on
 * ark_analytic_lsrk_domeigest.c by Mustafa Aggul @ SMU and
 * ark_brusselator.c by Daniel R. Reynolds @ UMBC
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
 * We use u0=1.2,  v0=3.1,  w0=3,  a=1,  b=3.5,  ep=5.0e-6.
 *  
 * In this case, w experiences a fast initial transient, jumping 0.5
 * within a few steps. All values proceed smoothly until
 * around t=6.5, when both u and v undergo a sharp transition,
 * with u increasing from around 0.5 to 5 and v decreasing
 * from around 6 to 1 in less than 0.5 time units. After this
 * transition, both u and v continue to evolve somewhat
 * rapidly for another 1.4 time units, and finish off smoothly.
 *
 * This program solves the problem with an STS method from LSRKStep using a
 * SUNDIALS dominant eigenvalue estimation (DEE) module.
 *
 * 100 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
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

/* User-supplied Functions Called by the Solver */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/* Private function to check function return values */
static int check_flag(void* flagvalue, const char* funcname, int opt);

/* Main Program */
int main(int argc, char* argv[])
{
  /* general problem parameters */
  sunrealtype T0    = SUN_RCONST(0.0);       /* initial time */
  sunrealtype Tf    = SUN_RCONST(10.0);      /* final time */
  sunrealtype dTout = SUN_RCONST(1.0);       /* time between outputs */
  sunindextype NEQ  = 3;                     /* number of dependent vars. */
  int Nt            = (int)ceil(Tf / dTout); /* number of output times */
  sunrealtype a, b, ep, u0, v0, w0;

#if defined(SUNDIALS_DOUBLE_PRECISION)
  sunrealtype reltol = SUN_RCONST(1.0e-6); /* tolerances */
  sunrealtype abstol = SUN_RCONST(1.0e-10);
#elif defined(SUNDIALS_SINGLE_PRECISION)
  sunrealtype reltol = SUN_RCONST(1.0e-4); /* tolerances */
  sunrealtype abstol = SUN_RCONST(1.0e-8);
#elif defined(SUNDIALS_EXTENDED_PRECISION)
  sunrealtype reltol = SUN_RCONST(1.0e-6); /* tolerances */
  sunrealtype abstol = SUN_RCONST(1.0e-10);
#endif

  /* general problem variables */
  int flag;                /* reusable error-checking flag */
  N_Vector y       = NULL; /* empty vector for storing solution */
  void* arkode_mem = NULL; /* empty ARKode memory structure */
  sunrealtype rdata[3];
  sunrealtype t, tout;
  int iout;

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

  /* set up the test problem */
  u0 = SUN_RCONST(1.2);
  v0 = SUN_RCONST(3.1);
  w0 = SUN_RCONST(3.0);
  a  = SUN_RCONST(1.0);
  b  = SUN_RCONST(3.5);
  ep = SUN_RCONST(5.0e-6);

  /* Initial problem output */
  printf("\nBrusselator ODE test problem:\n");
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

  sunrealtype* ydata = N_VGetArrayPointer(y);
  ydata[0]           = u0; /* Set initial conditions */
  ydata[1]           = v0;
  ydata[2]           = w0;

  /* Call LSRKStepCreateSTS to initialize the STS timestepper module and
     specify the right-hand side function in y'=f(t,y), the initial time
     T0, and the initial dependent variable vector y. */
  arkode_mem = LSRKStepCreateSTS(f, T0, y, ctx);
  if (check_flag((void*)arkode_mem, "LSRKStepCreateSTS", 0)) { return 1; }

  /* Set routines */
  flag = ARKodeSetUserData(arkode_mem,
                           (void*)rdata); /* Pass rdata to user functions */
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }

  flag = ARKodeSStolerances(arkode_mem, reltol, abstol); /* Specify tolerances */
  if (check_flag(&flag, "ARKodeSStolerances", 1)) { return 1; }

  flag = ARKodeSetInterpolantType(arkode_mem,
                                  ARK_INTERP_LAGRANGE); /* Specify stiff interpolant */
  if (check_flag(&flag, "ARKodeSetInterpolantType", 1)) { return 1; }

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
  There is no need to set Atimes or initialize since LSRKStep provides
  a default Atimes, and initialized the DEE, after it is attached. */
  flag = LSRKStepSetDomEigEstimator(arkode_mem, DEE);
  if (check_flag(&flag, "LSRKStepSetDomEigEstimator", 1)) { return 1; }

  /* Set the number of preprocessing warmups. The warmup
  is used to compute a "better" initial eigenvector and so an initial
  eigenvalue. The warmup is performed only once by the LSRKStep module
  internally unless LSRKStepSetNumDomEigEstPreprocessIters is called to set
  a new number of succeeding warmups that would be executed before
  every dominant eigenvalue estimation call */
  flag = LSRKStepSetNumDomEigEstInitPreprocessIters(arkode_mem, numwarmup);
  if (check_flag(&flag, "LSRKStepSetNumDomEigEstInitPreprocessIters", 1))
  {
    return 1;
  }

  /* Specify max number of stages allowed */
  flag = LSRKStepSetMaxNumStages(arkode_mem, 200);
  if (check_flag(&flag, "LSRKStepSetMaxNumStages", 1)) { return 1; }

  /* Specify max number of steps allowed */
  flag = ARKodeSetMaxNumSteps(arkode_mem, 2000);
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) { return 1; }

  /* Specify safety factor for user provided dom_eig */
  flag = LSRKStepSetDomEigSafetyFactor(arkode_mem, SUN_RCONST(1.01));
  if (check_flag(&flag, "LSRKStepSetDomEigSafetyFactor", 1)) { return 1; }

  /* Specify the Runge--Kutta--Legendre LSRK method by name */
  flag = LSRKStepSetSTSMethodByName(arkode_mem, "ARKODE_LSRK_RKL_2");
  if (check_flag(&flag, "LSRKStepSetSTSMethodByName", 1)) { return 1; }

  /* Override any current settings with command-line options */
  flag = SUNDomEigEstimator_SetOptions(DEE, NULL, NULL, argc, argv);
  if (check_flag(&flag, "SUNDomEigEstimator_SetOptions", 1)) { return 1; }

  /* Override any current settings with command-line options */
  flag = ARKodeSetOptions(arkode_mem, NULL, NULL, argc, argv);
  if (check_flag(&flag, "ARKodeSetOptions", 1)) { return 1; }

  /* Main time-stepping loop: calls ARKodeEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t    = T0;
  tout = T0 + dTout;
  printf("        t           u           v           w\n");
  printf("   -------------------------------------------\n");
  printf("  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM "\n", t,
         NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));

  for (iout = 0; iout < Nt; iout++)
  {
    flag = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_NORMAL); /* call integrator */
    if (check_flag(&flag, "ARKodeEvolve", 1)) { break; }
    printf("  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM
           "\n", /* access/print solution */
           t, NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));
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
  printf("   -------------------------------------------\n");

  /* Print final statistics */
  printf("\nFinal Statistics:\n");
  flag = ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(&flag, "ARKodePrintAllStats", 1)) { return 1; }

  /* Clean up and return with successful completion */
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
  sunrealtype a  = rdata[0];                    /* access data entries */
  sunrealtype b  = rdata[1];
  sunrealtype ep = rdata[2];
  sunrealtype u  = NV_Ith_S(y, 0); /* access solution values */
  sunrealtype v  = NV_Ith_S(y, 1);
  sunrealtype w  = NV_Ith_S(y, 2);

  /* fill in the RHS function */
  NV_Ith_S(ydot, 0) = a - (w + 1.0) * u + v * u * u;
  NV_Ith_S(ydot, 1) = w * u - v * u * u;
  NV_Ith_S(ydot, 2) = (b - w) / ep - w * u;

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
