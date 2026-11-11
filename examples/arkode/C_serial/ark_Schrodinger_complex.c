/*---------------------------------------------------------------
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
 * The following test simulates a simple 1D Schrodinger equation,
 *    i*hbar*psi_t = -(hbar^2/(2*m))*psi_xx + V0*psi
 * for t in [0, 1], x in [0, L), with initial conditions
 *    psi(0,x) = [exp(i*k*x) + i*exp(2*i*k*x)] / sqrt(2*L),
 * where k = 2*pi/L, periodic boundary conditions, i.e.
 *    psi(t,0) = psi(t,L),
 * and a constant, real potential V0 (V0 = 0 is a free particle).
 *
 * The spatial derivatives are computed using second-order
 * centered differences, with the data distributed over N points
 * on a uniform spatial grid.  The solution is complex-valued, so
 * this example requires a SUNDIALS build with a complex scalar
 * type (sunscalartype = suncomplextype).
 *
 * This program solves the problem with either an ERK or DIRK
 * method, selected by the optional command-line argument
 *    0 - ERK (default)
 *    1 - DIRK
 * For the DIRK method, we use a Newton iteration with the
 * SUNLinSol_SPGMR linear solver, and a user-supplied Jacobian-vector
 * product routine.  (The Newton system I - gamma*J is not Hermitian
 * since J = i*H for the Hermitian discrete Hamiltonian H, so a
 * GMRES-type solver is used in place of PCG.)
 *
 * Since each Fourier mode of the semidiscrete system evolves with
 * the known frequency
 *    omega_h(k) = 2*hbar*sin(k*dx/2)^2/(m*dx^2) + V0/hbar,
 * the time integration error is measured against this exact
 * semidiscrete solution.  The total probability
 *    P = dx*sum_j |psi_j|^2
 * is also monitored.
 *
 * 10 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 *---------------------------------------------------------------*/

/* Header files */
#include <arkode/arkode_arkstep.h> /* prototypes for ARKStep fcts., consts */
#include <math.h>
#include <nvector/nvector_serial.h> /* serial N_Vector types, fcts., macros */
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>  /* SUN_REAL, SUN_IMAG, SUNCabs, SUNCexp */
#include <sundials/sundials_types.h> /* defs. of sunscalartype, SUN_I, etc   */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver    */

#if !defined(SUNDIALS_SCALAR_TYPE_COMPLEX)
#error "This example requires SUNDIALS_SCALAR_TYPE=COMPLEX"
#endif

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM       "Lg"
#define ESYM       "Le"
#define FSYM       "Lf"
#define SUNRsin(x) (sinl((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define GSYM       "g"
#define ESYM       "e"
#define FSYM       "f"
#define SUNRsin(x) (sinf((x)))
#else
#define GSYM       "g"
#define ESYM       "e"
#define FSYM       "f"
#define SUNRsin(x) (sin((x)))
#endif

#define PI SUN_RCONST(3.14159265358979323846)

/* user data structure */
typedef struct
{
  sunindextype N;   /* number of intervals   */
  sunrealtype L;    /* domain length         */
  sunrealtype dx;   /* mesh spacing          */
  sunrealtype hbar; /* reduced Planck const. */
  sunrealtype m;    /* particle mass         */
  sunrealtype V0;   /* constant potential    */
}* UserData;

/* User-supplied Functions Called by the Solver */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Jac(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector y, N_Vector fy,
               void* user_data, N_Vector tmp);

/* Private functions to compute the exact semidiscrete solution and
   the probability (discrete L2 norm squared) of the solution */
static int Solution(sunrealtype t, N_Vector y, UserData udata);
static sunrealtype Probability(N_Vector y, UserData udata);

/* Private function to check function return values */
static int check_flag(void* flagvalue, const char* funcname, int opt);

/* Main Program */
int main(int argc, char* argv[])
{
  /* general problem parameters */
  sunrealtype T0   = SUN_RCONST(0.0); /* initial time */
  sunrealtype Tf   = SUN_RCONST(1.0); /* final time */
  int Nt           = 10;              /* total number of output times */
  sunrealtype rtol = 1.e-6;           /* relative tolerance */
  sunrealtype atol = 1.e-10;          /* absolute tolerance */
  UserData udata   = NULL;
  sunscalartype* data;
  sunindextype N   = 128;                  /* spatial mesh size */
  sunrealtype L    = SUN_RCONST(2.0) * PI; /* domain length */
  sunrealtype hbar = SUN_RCONST(1.0);      /* reduced Planck constant */
  sunrealtype m    = SUN_RCONST(1.0);      /* particle mass */
  sunrealtype V0   = SUN_RCONST(0.0);      /* constant potential */
  int method       = 0;                    /* 0 = ERK, 1 = DIRK */
  sunindextype i;

  /* general problem variables */
  int flag;                  /* reusable error-checking flag */
  N_Vector y         = NULL; /* empty vector for storing solution */
  N_Vector ytrue     = NULL; /* empty vector for storing exact solution */
  SUNLinearSolver LS = NULL; /* empty linear solver object */
  void* arkode_mem   = NULL; /* empty ARKode memory structure */
  FILE *FID, *UFID;
  sunrealtype t, dTout, tout, P0;
  int iout;
  long int nst, nst_a, nfe, nfi, nsetups, nli, nJv, nlcf, nni, ncfn, netf;

  /* Retrieve the time integration method from the command line */
  if (argc > 1) { method = atoi(argv[1]); }
  if (method != 0 && method != 1)
  {
    fprintf(stderr, "ERROR: method must be 0 (ERK) or 1 (DIRK)\n");
    return 1;
  }

  /* Create the SUNDIALS context object for this simulation */
  SUNContext ctx;
  flag = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_flag(&flag, "SUNContext_Create", 1)) { return 1; }

  /* allocate and fill udata structure */
  udata       = (UserData)malloc(sizeof(*udata));
  udata->N    = N;
  udata->L    = L;
  udata->hbar = hbar;
  udata->m    = m;
  udata->V0   = V0;
  udata->dx   = L / N; /* mesh spacing (periodic) */

  /* Initial problem output */
  printf("\n1D Schrodinger PDE test problem:\n");
  printf("  N = %li\n", (long int)udata->N);
  printf("  domain length:  L = %" GSYM "\n", udata->L);
  printf("  hbar = %" GSYM ",  m = %" GSYM ",  V0 = %" GSYM "\n", udata->hbar,
         udata->m, udata->V0);
  printf("  time integration method:  %s\n", (method == 0) ? "ERK" : "DIRK");

  /* Initialize data structures */
  y = N_VNew_Serial(N, ctx); /* Create serial vector for solution */
  if (check_flag((void*)y, "N_VNew_Serial", 0)) { return 1; }
  ytrue = N_VClone(y); /* Create serial vector for exact solution */
  if (check_flag((void*)ytrue, "N_VClone", 0)) { return 1; }
  flag = Solution(T0, y, udata); /* Set initial conditions */
  if (check_flag(&flag, "Solution", 1)) { return 1; }
  P0 = Probability(y, udata); /* Initial probability */

  /* Call ARKStepCreate to initialize the ARK timestepper module and
     specify the right-hand side function in y'=f(t,y), the initial time
     T0, and the initial dependent variable vector y.  For the ERK method
     we set f_E to f and f_I to NULL; for the DIRK method the problem is
     fully implicit, so we set f_E to NULL and f_I to f. */
  if (method == 0)
  {
    arkode_mem = ARKStepCreate(f, NULL, T0, y, ctx);
  }
  else
  {
    arkode_mem = ARKStepCreate(NULL, f, T0, y, ctx);
  }
  if (check_flag((void*)arkode_mem, "ARKStepCreate", 0)) { return 1; }

  /* Set routines */
  flag = ARKodeSetUserData(arkode_mem,
                           (void*)udata); /* Pass udata to user functions */
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }
  flag = ARKodeSetMaxNumSteps(arkode_mem, 100000); /* Increase max num steps */
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) { return 1; }
  flag = ARKodeSStolerances(arkode_mem, rtol, atol); /* Specify tolerances */
  if (check_flag(&flag, "ARKodeSStolerances", 1)) { return 1; }

  if (method == 1)
  {
    flag = ARKodeSetPredictorMethod(arkode_mem,
                                    1); /* Specify maximum-order predictor */
    if (check_flag(&flag, "ARKodeSetPredictorMethod", 1)) { return 1; }

    /* Initialize SPGMR solver -- no preconditioning, with up to N iterations  */
    LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, (int)N, ctx);
    if (check_flag((void*)LS, "SUNLinSol_SPGMR", 0)) { return 1; }

    /* Linear solver interface -- set user-supplied J*v routine (no 'jtsetup' required) */
    flag = ARKodeSetLinearSolver(arkode_mem, LS,
                                 NULL); /* Attach linear solver to ARKODE */
    if (check_flag(&flag, "ARKodeSetLinearSolver", 1)) { return 1; }
    flag = ARKodeSetJacTimes(arkode_mem, NULL, Jac); /* Set the Jacobian routine */
    if (check_flag(&flag, "ARKodeSetJacTimes", 1)) { return 1; }

    /* Specify linearly implicit RHS, with non-time-dependent Jacobian */
    flag = ARKodeSetLinear(arkode_mem, 0);
    if (check_flag(&flag, "ARKodeSetLinear", 1)) { return 1; }
  }

  /* Override any current settings with command-line options */
  flag = ARKodeSetOptions(arkode_mem, NULL, NULL, argc, argv);
  if (check_flag(&flag, "ARKodeSetOptions", 1)) { return 1; }
  if (method == 1)
  {
    flag = SUNLinSolSetOptions(LS, NULL, NULL, argc, argv);
    if (check_flag(&flag, "SUNLinSolSetOptions", 1)) { return 1; }
  }

  /* output mesh to disk */
  FID = fopen("schrodinger_mesh.txt", "w");
  for (i = 0; i < N; i++) { fprintf(FID, "  %.16" ESYM "\n", udata->dx * i); }
  fclose(FID);

  /* Open output stream for results, access data array */
  UFID = fopen("schrodinger1D.txt", "w");
  data = N_VGetArrayPointer(y);

  /* output initial condition to disk (real and imaginary parts) */
  for (i = 0; i < N; i++)
  {
    fprintf(UFID, " %.16" ESYM " %.16" ESYM "", SUN_REAL(data[i]),
            SUN_IMAG(data[i]));
  }
  fprintf(UFID, "\n");

  /* Main time-stepping loop: calls ARKodeEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t     = T0;
  dTout = (Tf - T0) / Nt;
  tout  = T0 + dTout;
  printf("        t      ||psi||_rms   probability    ||error||_rms\n");
  printf("   ---------------------------------------------------------\n");
  printf("  %10.6" FSYM "  %10.6" FSYM "  %14.10" FSYM "  %12.4" ESYM "\n", t,
         SUNRsqrt(Probability(y, udata) / L), Probability(y, udata),
         SUN_RCONST(0.0));
  for (iout = 0; iout < Nt; iout++)
  {
    flag = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_NORMAL); /* call integrator */
    if (check_flag(&flag, "ARKodeEvolve", 1)) { break; }
    flag = Solution(t, ytrue, udata); /* compute exact semidiscrete solution */
    if (check_flag(&flag, "Solution", 1)) { break; }
    N_VLinearSum(SUN_RCONST(1.0), y, -SUN_RCONST(1.0), ytrue, ytrue); /* error */
    printf("  %10.6" FSYM "  %10.6" FSYM "  %14.10" FSYM "  %12.4" ESYM "\n", t,
           SUNRsqrt(Probability(y, udata) / L), Probability(y, udata),
           SUNRsqrt(Probability(ytrue, udata) / L)); /* print solution stats */
    if (flag >= 0)
    { /* successful solve: update output time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
    else
    { /* unsuccessful solve: break */
      fprintf(stderr, "Solver failure, stopping integration\n");
      break;
    }

    /* output results to disk (real and imaginary parts) */
    for (i = 0; i < N; i++)
    {
      fprintf(UFID, " %.16" ESYM " %.16" ESYM "", SUN_REAL(data[i]),
              SUN_IMAG(data[i]));
    }
    fprintf(UFID, "\n");
  }
  printf("   ---------------------------------------------------------\n");
  printf("  probability drift |P - P0| = %.4" ESYM "\n",
         SUNRabs(Probability(y, udata) - P0));
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

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  if (method == 1)
  {
    flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
    check_flag(&flag, "ARKodeGetNumLinSolvSetups", 1);
    flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
    check_flag(&flag, "ARKodeGetNumNonlinSolvIters", 1);
    flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
    check_flag(&flag, "ARKodeGetNumNonlinSolvConvFails", 1);
    flag = ARKodeGetNumLinIters(arkode_mem, &nli);
    check_flag(&flag, "ARKodeGetNumLinIters", 1);
    flag = ARKodeGetNumJtimesEvals(arkode_mem, &nJv);
    check_flag(&flag, "ARKodeGetNumJtimesEvals", 1);
    flag = ARKodeGetNumLinConvFails(arkode_mem, &nlcf);
    check_flag(&flag, "ARKodeGetNumLinConvFails", 1);

    printf("   Total linear solver setups = %li\n", nsetups);
    printf("   Total linear iterations = %li\n", nli);
    printf("   Total number of Jacobian-vector products = %li\n", nJv);
    printf("   Total number of linear solver convergence failures = %li\n", nlcf);
    printf("   Total number of Newton iterations = %li\n", nni);
    printf("   Total number of nonlinear solver convergence failures = %li\n",
           ncfn);
  }
  printf("   Total number of error test failures = %li\n", netf);

  /* Clean up and return with successful completion */
  N_VDestroy(y); /* Free vectors */
  N_VDestroy(ytrue);
  free(udata);             /* Free user data */
  ARKodeFree(&arkode_mem); /* Free integrator memory */
  SUNLinSolFree(LS);       /* Free linear solver */
  SUNContext_Free(&ctx);   /* Free context */

  return 0;
}

/*--------------------------------
 * Functions called by the solver
 *--------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  UserData udata   = (UserData)user_data; /* access problem data */
  sunindextype N   = udata->N;            /* set variable shortcuts */
  sunrealtype hbar = udata->hbar;
  sunrealtype m    = udata->m;
  sunrealtype V0   = udata->V0;
  sunrealtype dx   = udata->dx;
  sunscalartype *Y = NULL, *Ydot = NULL;
  sunscalartype c1, c2;
  sunindextype i;

  Y = N_VGetArrayPointer(y); /* access data arrays */
  if (check_flag((void*)Y, "N_VGetArrayPointer", 0)) { return 1; }
  Ydot = N_VGetArrayPointer(ydot);
  if (check_flag((void*)Ydot, "N_VGetArrayPointer", 0)) { return 1; }
  N_VConst(0.0, ydot); /* Initialize ydot to zero */

  /* iterate over domain, computing all equations; dividing the PDE by
     i*hbar gives a positive imaginary coefficient for the Laplacian and
     a negative imaginary coefficient for the potential term */
  c1      = SUN_I * hbar / (SUN_RCONST(2.0) * m) / dx / dx;
  c2      = -SUN_RCONST(2.0) * c1 - SUN_I * V0 / hbar;
  Ydot[0] = c1 * Y[N - 1] + c2 * Y[0] + c1 * Y[1]; /* left (periodic) */
  for (i = 1; i < N - 1; i++)
  {
    Ydot[i] = c1 * Y[i - 1] + c2 * Y[i] + c1 * Y[i + 1];
  }
  Ydot[N - 1] = c1 * Y[N - 2] + c2 * Y[N - 1] + c1 * Y[0]; /* right (periodic) */

  return 0; /* Return with success */
}

/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(N_Vector v, N_Vector Jv, sunrealtype t, N_Vector y, N_Vector fy,
               void* user_data, N_Vector tmp)
{
  UserData udata   = (UserData)user_data; /* variable shortcuts */
  sunindextype N   = udata->N;
  sunrealtype hbar = udata->hbar;
  sunrealtype m    = udata->m;
  sunrealtype V0   = udata->V0;
  sunrealtype dx   = udata->dx;
  sunscalartype *V = NULL, *JV = NULL;
  sunscalartype c1, c2;
  sunindextype i;

  V = N_VGetArrayPointer(v); /* access data arrays */
  if (check_flag((void*)V, "N_VGetArrayPointer", 0)) { return 1; }
  JV = N_VGetArrayPointer(Jv);
  if (check_flag((void*)JV, "N_VGetArrayPointer", 0)) { return 1; }
  N_VConst(0.0, Jv); /* initialize Jv product to zero */

  /* iterate over domain, computing all Jacobian-vector products */
  c1    = SUN_I * hbar / (SUN_RCONST(2.0) * m) / dx / dx;
  c2    = -SUN_RCONST(2.0) * c1 - SUN_I * V0 / hbar;
  JV[0] = c1 * V[N - 1] + c2 * V[0] + c1 * V[1]; /* left (periodic) */
  for (i = 1; i < N - 1; i++)
  {
    JV[i] = c1 * V[i - 1] + c2 * V[i] + c1 * V[i + 1];
  }
  JV[N - 1] = c1 * V[N - 2] + c2 * V[N - 1] + c1 * V[0]; /* right (periodic) */

  return 0; /* Return with success */
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Compute the exact solution of the semidiscrete system at time t,
   psi_j(t) = [exp(i*(k*x_j - omega_h(k)*t)) + i*exp(i*(2k*x_j - omega_h(2k)*t))]
              / sqrt(2*L),
   with omega_h(k) = 2*hbar*sin(k*dx/2)^2/(m*dx^2) + V0/hbar. */
static int Solution(sunrealtype t, N_Vector y, UserData udata)
{
  sunindextype N   = udata->N; /* set variable shortcuts */
  sunrealtype L    = udata->L;
  sunrealtype hbar = udata->hbar;
  sunrealtype m    = udata->m;
  sunrealtype V0   = udata->V0;
  sunrealtype dx   = udata->dx;
  sunscalartype* Y = NULL;
  sunrealtype k, s, omega, x;
  sunscalartype coef;
  sunindextype i;
  int mode;

  Y = N_VGetArrayPointer(y); /* access data array */
  if (check_flag((void*)Y, "N_VGetArrayPointer", 0)) { return 1; }
  N_VConst(0.0, y); /* Initialize y to zero */

  /* accumulate the two Fourier modes */
  for (mode = 1; mode <= 2; mode++)
  {
    k     = mode * SUN_RCONST(2.0) * PI / L;
    s     = SUNRsin(SUN_RCONST(0.5) * k * dx);
    omega = SUN_RCONST(2.0) * hbar * s * s / (m * dx * dx) + V0 / hbar;
    coef  = (mode == 1) ? SUN_RCONST(1.0) : SUN_I;
    for (i = 0; i < N; i++)
    {
      x = dx * i;
      Y[i] += coef * SUNCexp(SUN_I * (k * x - omega * t));
    }
  }
  N_VScale(SUN_RCONST(1.0) / SUNRsqrt(SUN_RCONST(2.0) * L), y, y); /* normalize */

  return 0; /* Return with success */
}

/* Compute the total probability P = dx*sum_j |psi_j|^2 */
static sunrealtype Probability(N_Vector y, UserData udata)
{
  sunindextype N   = udata->N;
  sunrealtype dx   = udata->dx;
  sunscalartype* Y = N_VGetArrayPointer(y);
  sunrealtype P    = SUN_RCONST(0.0);
  sunrealtype a;
  sunindextype i;

  for (i = 0; i < N; i++)
  {
    a = SUNCabs(Y[i]);
    P += a * a;
  }

  return dx * P;
}

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
