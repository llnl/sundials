/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 * This example problem is adapted from:
 *
 * H. Ranocha, M. Sayyari, L. Dalcin, M. Parsani, and D.I. Ketcheson,
 * "Relaxation Runge-Kutta Methods: Fully-Discrete Explicit Entropy-Stable
 * Schemes for the Compressible Euler and Navier-Stokes Equations," SIAM Journal
 * on Scientific Computing, 42(2), 2020, https://doi.org/10.1137/19M1263480.
 * -----------------------------------------------------------------------------
 * This example evolves the equation du/dt = -exp(u) for t in the interval
 * [0, 5] with the initial condition u(0) = 0.5. The equation has the analytic
 * solution u(t) = -log(e^{-0.5} + t) and the dissipated exponential entropy is
 * given by ent(u) = exp(u) with Jacobian ent'(u) = de/du = exp(u).
 *
 * The problem is advanced in time with an explicit or implicit relaxed
 * Runge-Kutta method to ensure dissipation of the entropy.
 * ---------------------------------------------------------------------------*/

/* Header files */
#include <math.h>
#include <stdio.h>

/* SUNDIALS headers */
#include <arkode/arkode_arkstep.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>

/* Convince macros for calling precision-specific math functions */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define EXP(x) (exp((x)))
#define LOG(x) (log((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define EXP(x) (expf((x)))
#define LOG(x) (logf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define EXP(x) (expl((x)))
#define LOG(x) (logl((x)))
#endif

/* Convince macros for using precision-specific format specifiers */
#if defined(SUNDIALS_EXTENDED_PRECISION)
#define ESYM "Le"
#else
#define ESYM "e"
#endif

/* ----------------------- *
 * User-supplied functions *
 * ----------------------- */

/* ODE RHS function */
int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/* ODE RHS Jacobian function */
int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Entropy function */
int Ent(N_Vector y, sunrealtype* e, void* user_data);

/* Entropy Jacobian function */
int JacEnt(N_Vector y, N_Vector J, void* user_data);

/* ----------------- *
 * Utility functions *
 * ----------------- */

/* Analytic solution */
int ans(sunrealtype t, N_Vector y);

/* Check for an unrecoverable (negative) return flag from a SUNDIALS function */
int check_flag(int flag, const char* funcname);

/* Check if a function returned a NULL pointer */
int check_ptr(void* ptr, const char* funcname);

/* ------------ *
 * Main Program *
 * ------------ */

int main(int argc, char* argv[])
{
  /* Error-checking flag */
  int flag;

  /* Initial and final times */
  sunrealtype t0 = SUN_RCONST(0.0);
  sunrealtype tf = SUN_RCONST(5.0);

  /* Relative and absolute tolerances */
  sunrealtype reltol = SUN_RCONST(1.0e-6);
  sunrealtype abstol = SUN_RCONST(1.0e-10);

  /* SUNDIALS context, vector, matrix, and linear solver objects */
  SUNContext ctx     = NULL;
  N_Vector y         = NULL;
  N_Vector ytrue     = NULL;
  SUNMatrix A        = NULL;
  SUNLinearSolver LS = NULL;

  /* Pointer to vector data array */
  sunrealtype* ydata;
  sunrealtype* ytdata;

  /* Initial, current, and change in entropy value */
  sunrealtype ent0, ent, delta_ent;

  /* Solution errors */
  sunrealtype u_err;

  /* ARKODE memory structure */
  void* arkode_mem = NULL;

  /* ARKODE statistics */
  long int nst, nst_a, nfe, nfi;
  long int nrf, nrbf, nre, nrje, nrnlsi, nrnlsf;
  long int nsetups, nje, nfeLS, nni, ncfn, netf;

  /* Output time and file */
  sunrealtype t;
  FILE* UFID;

  /* Command line options */
  int relax           = 1;               /* enable relaxation */
  int implicit        = 1;               /* implicit          */
  sunrealtype fixed_h = SUN_RCONST(0.0); /* adaptive stepping */

  /* -------------------- *
   * Output Problem Setup *
   * -------------------- */

  if (argc > 1) { relax = atoi(argv[1]); }
  if (argc > 2) { implicit = atoi(argv[2]); }
  if (argc > 3) { fixed_h = atof(argv[3]); }

  printf("\nDissipated Exponential Entropy problem:\n");
  if (implicit) { printf("   method     = DIRK\n"); }
  else { printf("   method     = ERK\n"); }
  printf("   reltol     = %.1" ESYM "\n", reltol);
  printf("   abstol     = %.1" ESYM "\n", abstol);
  if (fixed_h > SUN_RCONST(0.0))
  {
    printf("   fixed h    = %.1" ESYM "\n", fixed_h);
  }
  if (relax) { printf("   relaxation = ON\n"); }
  else { printf("   relaxation = OFF\n"); }
  printf("\n");

  /* ------------ *
   * Setup ARKODE *
   * ------------ */

  /* Create the SUNDIALS context object for this simulation */
  flag = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_flag(flag, "SUNContext_Create")) { return 1; }

  /* Create serial vector and set the initial condition values */
  y = N_VNew_Serial(1, ctx);
  if (check_ptr(y, "N_VNew_Serial")) { return 1; }

  ydata = N_VGetArrayPointer(y);
  if (check_ptr(ydata, "N_VGetArrayPointer")) { return 1; }

  ydata[0] = SUN_RCONST(0.5);

  ytrue = N_VClone(y);
  if (check_ptr(ytrue, "N_VClone")) { return 1; }

  ytdata = N_VGetArrayPointer(ytrue);
  if (check_ptr(ytdata, "N_VGetArrayPointer")) { return 1; }

  /* Initialize ARKStep */
  if (implicit) { arkode_mem = ARKStepCreate(NULL, f, t0, y, ctx); }
  else { arkode_mem = ARKStepCreate(f, NULL, t0, y, ctx); }
  if (check_ptr(arkode_mem, "ARKStepCreate")) { return 1; }

  /* Set order */
  flag = ARKodeSetOrder(arkode_mem, 2);
  if (check_flag(flag, "ARKodeSetOrder")) { return 1; }

  /* Specify tolerances */
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(flag, "ARKodeSStolerances")) { return 1; }

  if (relax)
  {
    /* Enable relaxation methods */
    flag = ARKodeSetRelaxFn(arkode_mem, Ent, JacEnt);
    if (check_flag(flag, "ARKodeSetRelaxFn")) { return 1; }
  }

  if (implicit)
  {
    /* Create dense matrix and linear solver */
    A = SUNDenseMatrix(1, 1, ctx);
    if (check_ptr(A, "SUNDenseMatrix")) { return 1; }

    LS = SUNLinSol_Dense(y, A, ctx);
    if (check_ptr(LS, "SUNLinSol_Dense")) { return 1; }

    /* Attach the matrix and linear solver */
    flag = ARKodeSetLinearSolver(arkode_mem, LS, A);
    if (check_flag(flag, "ARKodeSetLinearSolver")) { return 1; }

    /* Set Jacobian routine */
    flag = ARKodeSetJacFn(arkode_mem, Jac);
    if (check_flag(flag, "ARKodeSetJacFn")) { return 1; }

    /* Tighten nonlinear solver tolerance */
    flag = ARKodeSetNonlinConvCoef(arkode_mem, SUN_RCONST(0.01));
    if (check_flag(flag, "ARKodeSetNonlinConvCoef")) { return 1; }
  }

  if (fixed_h > SUN_RCONST(0.0))
  {
    flag = ARKodeSetFixedStep(arkode_mem, fixed_h);
    if (check_flag(flag, "ARKodeSetFixedStep")) { return 1; }
  }

  /* Open output stream for results, output comment line */
  UFID = fopen("ark_dissipated_exp_entropy.txt", "w");
  fprintf(UFID, "# vars: t u entropy u_err delta_entropy\n");

  /* --------------- *
   * Advance in Time *
   * --------------- */

  /* Initial time */
  t = t0;

  /* Output the initial condition and entropy */
  flag = Ent(y, &ent0, NULL);
  if (check_flag(flag, "Ent")) { return 1; }

  fprintf(UFID,
          "%23.16" ESYM " %23.16" ESYM " %23.16" ESYM " %23.16" ESYM
          " %23.16" ESYM "\n",
          t0, ydata[0], ent0, SUN_RCONST(0.0), SUN_RCONST(0.0));

  printf(" step   t              u              e              u_err          "
         "delta e\n");
  printf(" --------------------------------------------------------------------"
         "-----------\n");
  printf("%5d %14.6" ESYM " %14.6" ESYM " %14.6" ESYM " %14.6" ESYM
         " %14.6" ESYM "\n",
         0, t, ydata[0], ent0, SUN_RCONST(0.0), SUN_RCONST(0.0));

  while (t < tf)
  {
    /* Evolve in time */
    flag = ARKodeEvolve(arkode_mem, tf, y, &t, ARK_ONE_STEP);
    if (check_flag(flag, "ARKodeEvolve")) { break; }

    /* Output solution and errors */
    flag = Ent(y, &ent, NULL);
    if (check_flag(flag, "Ent")) { return 1; }

    flag = ans(t, ytrue);
    if (check_flag(flag, "ans")) { return 1; }

    delta_ent = ent - ent0;
    u_err     = ydata[0] - ytdata[0];

    /* Output to the screen periodically */
    flag = ARKodeGetNumSteps(arkode_mem, &nst);
    check_flag(flag, "ARKodeGetNumSteps");

    if (nst % 40 == 0)
    {
      printf("%5ld %14.6" ESYM " %14.6" ESYM " %14.6" ESYM " %14.6" ESYM
             " %14.6" ESYM "\n",
             nst, t, ydata[0], ent, u_err, delta_ent);
    }

    /* Write all steps to file */
    fprintf(UFID,
            "%23.16" ESYM " %23.16" ESYM " %23.16" ESYM " %23.16" ESYM
            " %23.16" ESYM "\n",
            t, ydata[0], ent, u_err, delta_ent);
  }

  printf(" --------------------------------------------------------------------"
         "-----------\n");
  fclose(UFID);

  /* ------------ *
   * Output Stats *
   * ------------ */

  /* Get final statistics on how the solve progressed */
  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  check_flag(flag, "ARKodeGetNumSteps");

  flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  check_flag(flag, "ARKodeGetNumStepAttempts");

  flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
  check_flag(flag, "ARKodeGetNumErrTestFails");

  flag = ARKodeGetNumRhsEvals(arkode_mem, 0, &nfe);
  check_flag(flag, "ARKodeGetNumRhsEvals");

  flag = ARKodeGetNumRhsEvals(arkode_mem, 1, &nfi);
  check_flag(flag, "ARKodeGetNumRhsEvals");

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total number of error test failures = %li\n", netf);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);

  if (implicit)
  {
    flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
    check_flag(flag, "ARKodeGetNumNonlinSolvIters");

    flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
    check_flag(flag, "ARKodeGetNumNonlinSolvConvFails");

    flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
    check_flag(flag, "ARKodeGetNumLinSolvSetups");

    flag = ARKodeGetNumJacEvals(arkode_mem, &nje);
    check_flag(flag, "ARKodeGetNumJacEvals");

    flag = ARKodeGetNumLinRhsEvals(arkode_mem, &nfeLS);
    check_flag(flag, "ARKodeGetNumLinRhsEvals");

    printf("   Total number of Newton iterations = %li\n", nni);
    printf("   Total number of linear solver convergence failures = %li\n", ncfn);
    printf("   Total linear solver setups = %li\n", nsetups);
    printf("   Total number of Jacobian evaluations = %li\n", nje);
    printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
  }

  if (relax)
  {
    flag = ARKodeGetNumRelaxFnEvals(arkode_mem, &nre);
    check_flag(flag, "ARKodeGetNumRelaxFnEvals");

    flag = ARKodeGetNumRelaxJacEvals(arkode_mem, &nrje);
    check_flag(flag, "ARKodeGetNumRelaxJacEvals");

    flag = ARKodeGetNumRelaxFails(arkode_mem, &nrf);
    check_flag(flag, "ARKodeGetNumRelaxFails");

    flag = ARKodeGetNumRelaxBoundFails(arkode_mem, &nrbf);
    check_flag(flag, "ARKodeGetNumRelaxBoundFails");

    flag = ARKodeGetNumRelaxSolveFails(arkode_mem, &nrnlsf);
    check_flag(flag, "ARKodeGetNumRelaxSolveFails");

    flag = ARKodeGetNumRelaxSolveIters(arkode_mem, &nrnlsi);
    check_flag(flag, "ARKodeGetNumRelaxSolveIters");

    printf("   Total Relaxation Fn evals    = %li\n", nre);
    printf("   Total Relaxation Jac evals   = %li\n", nrje);
    printf("   Total Relaxation fails       = %li\n", nrf);
    printf("   Total Relaxation bound fails = %li\n", nrbf);
    printf("   Total Relaxation NLS fails   = %li\n", nrnlsf);
    printf("   Total Relaxation NLS iters   = %li\n", nrnlsi);
  }
  printf("\n");

  /* -------- *
   * Clean up *
   * -------- */

  /* Free ARKODE integrator and SUNDIALS objects */
  ARKodeFree(&arkode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  N_VDestroy(ytrue);
  SUNContext_Free(&ctx);

  return flag;
}

/* ----------------------- *
 * User-supplied functions *
 * ----------------------- */

/* ODE RHS function f(t,y). */
int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* fdata = N_VGetArrayPointer(ydot);
  fdata[0]           = -EXP(ydata[0]);
  return 0;
}

/* ODE RHS Jacobian function J(t,y) = df/dy. */
int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data,
        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* Jdata = SUNDenseMatrix_Data(J);
  Jdata[0]           = -EXP(ydata[0]);
  return 0;
}

/* Entropy function e(y) */
int Ent(N_Vector y, sunrealtype* e, void* user_data)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);
  *e                 = EXP(ydata[0]);
  return 0;
}

/* Entropy function Jacobian Je(y) = de/dy */
int JacEnt(N_Vector y, N_Vector J, void* user_data)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* jdata = N_VGetArrayPointer(J);
  jdata[0]           = EXP(ydata[0]);
  return 0;
}

/* ----------------- *
 * Utility functions *
 * ----------------- */

/* Analytic solution */
int ans(sunrealtype t, N_Vector y)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);
  ydata[0]           = -LOG(EXP(SUN_RCONST(-0.5)) + t);
  return 0;
}

/* Check for an unrecoverable (negative) return flag from a SUNDIALS function */
int check_flag(int flag, const char* funcname)
{
  if (flag < 0)
  {
    fprintf(stderr, "ERROR: %s() returned %d\n", funcname, flag);
    return 1;
  }
  return 0;
}

/* Check if a function returned a NULL pointer */
int check_ptr(void* ptr, const char* funcname)
{
  if (!ptr)
  {
    fprintf(stderr, "ERROR: %s() returned NULL\n", funcname);
    return 1;
  }
  return 0;
}
