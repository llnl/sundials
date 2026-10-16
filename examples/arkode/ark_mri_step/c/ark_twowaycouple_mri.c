/* ----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ----------------------------------------------------------------
 * Based a linear example program by Rujeko Chinomona @ UMBC.
 * ----------------------------------------------------------------
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
 * ----------------------------------------------------------------
 * Example problem:
 *
 * This example simulates an ODE system with 3 components,
 * Y = [u,v,w], given by the equations,
 *
 *   du/dt =  100v+w
 *   dv/dt = -100u
 *   dw/dt = -w+u
 *
 * for t in the interval [0.0, 2.0] with initial conditions
 * u(0)=9001/10001, v(0)=-1e-5/10001, and w(0)=1000. In this problem
 * the slow (w) and fast (u and v) components depend on one another.
 *
 * This program solves the problem with the MRI stepper. Outputs are
 * printed at equal intervals of 0.1 and run statistics are printed
 * at the end.
 * ----------------------------------------------------------------*/

/* Header files */
#include <arkode/arkode_arkstep.h> /* prototypes for ARKStep fcts., consts */
#include <arkode/arkode_mristep.h> /* prototypes for MRIStep fcts., consts */
#include <math.h>
#include <nvector/nvector_serial.h> /* serial N_Vector types, fcts., macros */
#include <stdio.h>
#include <sundials/sundials_types.h> /* def. of type 'sunrealtype'              */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* User-supplied functions called by the solver */
static int fs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int ff(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/* Private function to check function return values */
static int check_retval(void* returnvalue, const char* funcname, int opt);

/* Main Program */
int main(void)
{
  /* general problem parameters */
  sunrealtype T0    = SUN_RCONST(0.0);       /* initial time */
  sunrealtype Tf    = SUN_RCONST(2.0);       /* final time */
  sunrealtype dTout = SUN_RCONST(0.1);       /* time between outputs */
  sunindextype NEQ  = 3;                     /* number of dependent vars. */
  int Nt            = (int)ceil(Tf / dTout); /* number of output times */
  sunrealtype hs    = SUN_RCONST(0.001);     /* slow step size */
  sunrealtype hf    = SUN_RCONST(0.00002);   /* fast step size */
  sunrealtype u0, v0, w0;                    /* initial conditions */

  /* general problem variables */
  int retval;                               /* reusable error-checking flag */
  N_Vector y                        = NULL; /* vector for the solution      */
  void* arkode_mem                  = NULL; /* ARKode memory structure      */
  void* inner_arkode_mem            = NULL; /* ARKode memory structure      */
  MRIStepInnerStepper inner_stepper = NULL; /* inner stepper                */
  FILE* UFID;
  sunrealtype t, tout;
  int iout;
  long int nsts, nstf, nfse, nff;

  /* Create the SUNDIALS context object for this simulation */
  SUNContext ctx;
  retval = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return 1; }

  /*
   * Initialization
   */

  /* Set the initial contions */
  u0 = SUN_RCONST(9001.0) / SUN_RCONST(10001.0);
  v0 = SUN_RCONST(-1.0e5) / SUN_RCONST(10001.0);
  w0 = SUN_RCONST(1000.0);

  /* Initial problem output */
  printf("\nTwo way coupling ODE test problem:\n");
  printf("    initial conditions:  u0 = %" GSYM ",  v0 = %" GSYM
         ",  w0 = %" GSYM "\n",
         u0, v0, w0);
  printf("    hs = %" GSYM ",  hf = %" GSYM "\n\n", hs, hf);

  /* Create and initialize serial vector for the solution */
  y = N_VNew_Serial(NEQ, ctx);
  if (check_retval((void*)y, "N_VNew_Serial", 0)) { return 1; }
  NV_Ith_S(y, 0) = u0;
  NV_Ith_S(y, 1) = v0;
  NV_Ith_S(y, 2) = w0;

  /*
   * Create the fast integrator and set options
   */

  /* Initialize the fast integrator. Specify the explicit fast right-hand side
     function in y'=fe(t,y)+fi(t,y)+ff(t,y), the initial time T0, and the
     initial dependent variable vector y. */
  inner_arkode_mem = ARKStepCreate(ff, NULL, T0, y, ctx);
  if (check_retval((void*)inner_arkode_mem, "ARKStepCreate", 0)) { return 1; }

  /* Set the fast method */
  retval = ARKStepSetTableNum(inner_arkode_mem, -1, ARKODE_KNOTH_WOLKE_3_3);
  if (check_retval(&retval, "ARKStepSetTableNum", 1)) { return 1; }

  /* Set the fast step size */
  retval = ARKodeSetFixedStep(inner_arkode_mem, hf);
  if (check_retval(&retval, "ARKodeSetFixedStep", 1)) { return 1; }

  /* Create inner stepper */
  retval = ARKodeCreateMRIStepInnerStepper(inner_arkode_mem, &inner_stepper);
  if (check_retval(&retval, "ARKodeCreateMRIStepInnerStepper", 1)) { return 1; }

  /*
   * Create the slow integrator and set options
   */

  /* Initialize the slow integrator. Specify the explicit slow right-hand side
     function in y'=fe(t,y)+fi(t,y)+ff(t,y), the initial time T0, the
     initial dependent variable vector y, and the fast integrator. */
  arkode_mem = MRIStepCreate(fs, NULL, T0, y, inner_stepper, ctx);
  if (check_retval((void*)arkode_mem, "MRIStepCreate", 0)) { return 1; }

  /* Set the slow step size */
  retval = ARKodeSetFixedStep(arkode_mem, hs);
  if (check_retval(&retval, "ARKodeSetFixedStep", 1)) { return 1; }

  /*
   * Integrate ODE
   */

  /* Open output stream for results, output comment line */
  UFID = fopen("ark_twowaycouple_mri_solution.txt", "w");
  fprintf(UFID, "# t u v w\n");

  /* output initial condition to disk */
  fprintf(UFID, " %.16" ESYM " %.16" ESYM " %.16" ESYM " %.16" ESYM "\n", T0,
          NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));

  /* Main time-stepping loop: calls ARKodeEvolve to perform the
     integration, then prints results. Stops when the final time
     has been reached */
  t    = T0;
  tout = T0 + dTout;
  printf("        t           u           v           w\n");
  printf("   -----------------------------------------------\n");
  printf("  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM "\n", t,
         NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));

  for (iout = 0; iout < Nt; iout++)
  {
    /* call integrator */
    retval = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_retval(&retval, "ARKodeEvolve", 1)) { break; }

    /* access/print solution and error */
    printf("  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM "\n", t,
           NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));
    fprintf(UFID, " %.16" ESYM " %.16" ESYM " %.16" ESYM " %.16" ESYM "\n", t,
            NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));

    /* successful solve: update time */
    tout += dTout;
    tout = (tout > Tf) ? Tf : tout;
  }
  printf("   -----------------------------------------------\n");
  fclose(UFID);

  /*
   * Finalize
   */

  /* Get some slow integrator statistics */
  retval = ARKodeGetNumSteps(arkode_mem, &nsts);
  check_retval(&retval, "ARKodeGetNumSteps", 1);
  retval = ARKodeGetNumRhsEvals(arkode_mem, 0, &nfse);
  check_retval(&retval, "ARKodeGetNumRhsEvals", 1);

  /* Get some fast integrator statistics */
  retval = ARKodeGetNumSteps(inner_arkode_mem, &nstf);
  check_retval(&retval, "ARKodeGetNumSteps", 1);
  retval = ARKodeGetNumRhsEvals(inner_arkode_mem, 0, &nff);
  check_retval(&retval, "ARKodeGetNumRhsEvals", 1);

  /* Print some final statistics */
  printf("\nFinal Solver Statistics:\n");
  printf("   Steps: nsts = %li, nstf = %li\n", nsts, nstf);
  printf("   Total RHS evals:  Fs = %li,  Ff = %li\n", nfse, nff);

  /* Clean up and return */
  N_VDestroy(y);                            /* Free y vector */
  ARKodeFree(&inner_arkode_mem);            /* Free integrator memory */
  MRIStepInnerStepper_Free(&inner_stepper); /* Free inner stepper */
  ARKodeFree(&arkode_mem);                  /* Free integrator memory */
  SUNContext_Free(&ctx);                    /* Free context */

  return 0;
}

/* ------------------------------
 * Functions called by the solver
 * ------------------------------*/

/* ff routine to compute the fast portion of the ODE RHS. */
static int ff(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype c1 = SUN_RCONST(100.0); /* problem constant */
  sunrealtype u  = NV_Ith_S(y, 0);    /* access solution values */
  sunrealtype v  = NV_Ith_S(y, 1);

  /* fill in the RHS function */
  NV_Ith_S(ydot, 0) = c1 * v;
  NV_Ith_S(ydot, 1) = -c1 * u;
  NV_Ith_S(ydot, 2) = u;

  /* Return with success */
  return 0;
}

/* fs routine to compute the slow portion of the ODE RHS. */
static int fs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype w = NV_Ith_S(y, 2); /* access solution values */

  /* fill in the RHS function */
  NV_Ith_S(ydot, 0) = w;
  NV_Ith_S(ydot, 1) = SUN_RCONST(0.0);
  NV_Ith_S(ydot, 2) = -w;

  /* Return with success */
  return 0;
}

/* ------------------------------
 * Private helper functions
 * ------------------------------*/

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a retval so check if
             retval < 0
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

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
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

/*---- end of file ----*/
