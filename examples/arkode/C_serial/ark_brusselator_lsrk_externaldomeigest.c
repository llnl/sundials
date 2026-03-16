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
 * The following test simulates the same problem as 
 * ark_brusselator_lsrk_domeigest.c.
 * 
 * The only difference with that example is that this one
 * attaches a "user provided" dominant eigenvalue estimate,
 * and where that estimate is provided by a user-constructed
 * SUNDomEigEstimator object (see the dom_eig function below).
 * 
 * While this example and ark_brusselator_lsrk_domeigest.c should provide
 * the same results, the purpose of this example is to show 
 * users how they may manually construct SUNDomEigEstimator
 * objects to query the dominant eigenvalues of a desired 
 * function.  In particular, we note that there is no 
 * requirement for the SUNDomEigEstimator to be used purely
 * for super-time-stepping methods in LSRKStep, so users 
 * interested in other purposes may focus on the contents 
 * of the dom_eig function below as a template.
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

/* User-supplied Dominated Eigenvalue Called by the Solver */
static int dom_eig(sunrealtype t, N_Vector y, N_Vector fn, sunrealtype* lambdaR,
                   sunrealtype* lambdaI, void* user_data, N_Vector temp1,
                   N_Vector temp2, N_Vector temp3);

/* Private function to check function return values */
static int check_flag(void* flagvalue, const char* funcname, int opt);

/* user data structure */
typedef struct
{
  SUNContext ctx;
  sunrealtype rdata[3];
  SUNDomEigEstimator DEE;
  sunrealtype rel_tol;
  sunindextype max_iters;
} UserData;

/* Main Program */
int main(int argc, char* argv[])
{
  /* general problem parameters */
  sunrealtype T0     = SUN_RCONST(0.0);       /* initial time */
  sunrealtype Tf     = SUN_RCONST(10.0);      /* final time */
  sunrealtype dTout  = SUN_RCONST(1.0);       /* time between outputs */
  sunindextype NEQ   = 3;                     /* number of dependent vars. */
  int Nt             = (int)ceil(Tf / dTout); /* number of output times */
  int test           = 2;                     /* test problem to run */
  sunrealtype reltol = SUN_RCONST(1.0e-6);    /* tolerances */
  sunrealtype abstol = SUN_RCONST(1.0e-10);
  sunrealtype a, b, ep, u0, v0, w0;

  /* general problem variables */
  int flag;                /* reusable error-checking flag */
  N_Vector y       = NULL; /* empty vector for storing solution */
  void* arkode_mem = NULL; /* empty ARKode memory structure */
  UserData ProbData;       /* problem data structure     */
  sunrealtype t, tout;
  int iout;

  /* Create the SUNDIALS context object for this simulation */
  SUNContext ctx;
  flag = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_flag(&flag, "SUNContext_Create", 1)) { return 1; }

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
  printf("\nBrusselator ODE test problem:\n");
  printf("    initial conditions:  u0 = %" GSYM ",  v0 = %" GSYM
         ",  w0 = %" GSYM "\n",
         u0, v0, w0);
  printf("    problem parameters:  a = %" GSYM ",  b = %" GSYM ",  ep = %" GSYM
         "\n",
         a, b, ep);
  printf("    reltol = %.1" ESYM ",  abstol = %.1" ESYM "\n\n", reltol, abstol);

  /* Initialize data structures */
  ProbData.ctx       = ctx;
  ProbData.rdata[0]  = a; /* set user data  */
  ProbData.rdata[1]  = b;
  ProbData.rdata[2]  = ep;
  ProbData.DEE       = NULL;
  ProbData.rel_tol   = SUN_RCONST(5.0e-3);
  ProbData.max_iters = 100;
  y = N_VNew_Serial(NEQ, ctx); /* Create serial vector for solution */
  if (check_flag((void*)y, "N_VNew_Serial", 0)) { return 1; }
  NV_Ith_S(y, 0) = u0; /* Set initial conditions */
  NV_Ith_S(y, 1) = v0;
  NV_Ith_S(y, 2) = w0;

  /* Call LSRKStepCreateSTS to initialize the ARK timestepper module and
     specify the right-hand side function in y'=f(t,y), the initial time
     T0, and the initial dependent variable vector y. */
  arkode_mem = LSRKStepCreateSTS(f, T0, y, ctx);
  if (check_flag((void*)arkode_mem, "LSRKStepCreateSTS", 0)) { return 1; }

  /* Set routines */
  flag = ARKodeSetUserData(arkode_mem,
                           (void*)&ProbData); /* Pass rdata to user functions */
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }

  flag = ARKodeSStolerances(arkode_mem, reltol, abstol); /* Specify tolerances */
  if (check_flag(&flag, "ARKodeSStolerances", 1)) { return 1; }

  flag = ARKodeSetInterpolantType(arkode_mem,
                                  ARK_INTERP_LAGRANGE); /* Specify stiff interpolant */
  if (check_flag(&flag, "ARKodeSetInterpolantType", 1)) { return 1; }

  /* Specify user provided spectral radius */
  flag = LSRKStepSetDomEigFn(arkode_mem, dom_eig);
  if (check_flag(&flag, "LSRKStepSetDomEigFn", 1)) { return 1; }

  /* Specify after how many successful steps dom_eig is recomputed */
  flag = LSRKStepSetDomEigFrequency(arkode_mem, 0);
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

  /* Specify the Runge--Kutta--Chebyshev LSRK method by name */
  flag = LSRKStepSetSTSMethodByName(arkode_mem, "ARKODE_LSRK_RKL_2");
  if (check_flag(&flag, "LSRKStepSetSTSMethodByName", 1)) { return 1; }

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
  N_VDestroy(y);                             /* Free y vector */
  ARKodeFree(&arkode_mem);                   /* Free integrator memory */
  SUNDomEigEstimator_Destroy(&ProbData.DEE); /* Free DEE object */
  SUNContext_Free(&ctx);                     /* Free context */

  return flag;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  UserData* data     = (UserData*)user_data; /* cast user_data to UserData */
  sunrealtype* rdata = data->rdata;          /* access rdata from UserData */
  sunrealtype a      = rdata[0];             /* access data entries */
  sunrealtype b      = rdata[1];
  sunrealtype ep     = rdata[2];
  sunrealtype u      = NV_Ith_S(y, 0); /* access solution values */
  sunrealtype v      = NV_Ith_S(y, 1);
  sunrealtype w      = NV_Ith_S(y, 2);

  /* fill in the RHS function */
  NV_Ith_S(ydot, 0) = a - (w + 1.0) * u + v * u * u;
  NV_Ith_S(ydot, 1) = w * u - v * u * u;
  NV_Ith_S(ydot, 2) = (b - w) / ep - w * u;

  return 0; /* Return with success */
}

/* dom_eig routine to estimate the dominated eigenvalue */
static int dom_eig(sunrealtype t, N_Vector y, N_Vector fn, sunrealtype* lambdaR,
                   sunrealtype* lambdaI, void* user_data, N_Vector temp1,
                   N_Vector temp2, N_Vector temp3)
{
  int flag;
  UserData* data = (UserData*)user_data; /* cast user_data to UserData */

  SUNContext ctx         = data->ctx; /* access context from UserData */
  SUNDomEigEstimator DEE = data->DEE; /* access DEE from UserData */

  /* DEE is initialized to NULL, so on the first dom_eig call we need
     to create and initialize this object */
  if (DEE == NULL)
  {
    /* Create random initial vector for power iteration */
    sunrealtype* qd  = N_VGetArrayPointer(temp1);
    sunindextype NEQ = N_VGetLength(temp1);
    for (int i = 0; i < NEQ; i++)
    {
      qd[i] = (sunrealtype)rand() / (sunrealtype)RAND_MAX;
    }

    /* Create power iteration dominant eigenvalue estimator (DEE) */
    DEE = SUNDomEigEstimator_Power(temp1, data->max_iters, data->rel_tol, ctx);
    if (check_flag(DEE, "SUNDomEigEstimator_Power", 0)) { return 1; }

    /* Set the ODE right-hand side function at t for the Jacobian-vector products */
    flag = SUNDomEigEstimator_SetRHS(DEE, user_data, f);
    if (check_flag(&flag, "SUNDomEigEstimator_SetRHS", 1)) { return 1; }

    /* Set the linearization vector for the Jacobian-vector products */
    flag = SUNDomEigEstimator_SetRHSLinearizationPoint(DEE, t,
                                                       y); // set the time t as well
    if (check_flag(&flag, "SUNDomEigEstimator_SetRHSLinearizationPoint", 1))
    {
      return 1;
    }

    flag = SUNDomEigEstimator_Initialize(DEE);
    if (check_flag(&flag, "SUNDomEigEstimator_Initialize", 1)) { return 1; }
  }

  /* Update the linearization vector and time for the Jacobian-vector products */
  flag = SUNDomEigEstimator_SetRHSLinearizationPoint(DEE, t, y);
  if (check_flag(&flag, "SUNDomEigEstimator_SetRHSLinearizationPoint", 1))
  {
    return 1;
  }

  /* Estimate the dominant eigenvalue with power iteration */
  flag = SUNDomEigEstimator_Estimate(DEE, lambdaR, lambdaI);
  if (check_flag(&flag, "SUNDomEigEstimator_Estimate", 1)) { return 1; }

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

/*---- end of file ----*/
