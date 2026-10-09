/* clang-format off */
/* ----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * ----------------------------------------------------------------------------
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
 * ----------------------------------------------------------------------------
 * We consider the Kepler problem. We choose one body to be the center of our
 * coordinate system and then we use the coordinates q = (q1, q2) to represent
 * the position of the second body relative to the first (center). This yields
 * the ODE:
 *    dq/dt = [ p1 ]
 *            [ p2 ]
 *    dp/dt = [ -q1 / (q1^2 + q2^2)^(3/2) ]
 *          = [ -q2 / (q1^2 + q2^2)^(3/2) ]
 * with the initial conditions
 *    q(0) = [ 1 - e ],  p(0) = [        0          ]
 *           [   0   ]          [ sqrt((1+e)/(1-e)) ]
 * where e = 0.6 is the eccentricity.
 *
 * The Hamiltonian for the system,
 *    H(p,q) = 1/2 * (p1^2 + p2^2) - 1/sqrt(q1^2 + q2^2)
 * is conserved as well as the angular momentum,
 *    L(p,q) = q1*p2 - q2*p1.
 *
 * By default we solve the problem by letting y = [ q, p ]^T then using a 4th
 * order symplectic integrator via the SPRKStep time-stepper of ARKODE with a
 * fixed time-step size.
 *
 * The rootfinding feature of SPRKStep is used to count the number of complete orbits.
 * This is done by defining the function,
 *    g(q) = q2
 * and providing it to SPRKStep as the function to find the roots for g(q).
 *
 * The program also accepts command line arguments to change the method
 * used and time-stepping strategy. The program has the following CLI arguments:
 *
 *   --step-mode <fixed, adapt>  should we use a fixed time-step or adaptive time-step (default fixed)
 *   --stepper <SPRK, ERK>       should we use SPRKStep or ARKStep with an ERK method (default SPRK)
 *   --method <string>           which method to use (default ARKODE_SPRK_MCLACHLAN_4_4)
 *   --use-compensated-sums      turns on compensated summation in ARKODE where applicable
 *   --disable-tstop             turns off tstop mode
 *   --dt <Real>                 the fixed-time step size to use if fixed time stepping is turned on (default 0.01)
 *   --tf <Real>                 the final time for the simulation (default 100)
 *   --nout                      number of output times
 *   --count-orbits              use rootfinding to count the number of completed orbits
 *   --check-order               compute the order of the method used and check if it is within the expected range
 *
 * References:
 *    Ernst Hairer, Christain Lubich, Gerhard Wanner
 *    Geometric Numerical Integration: Structure-Preserving
 *    Algorithms for Ordinary Differential Equations
 *    Springer, 2006,
 *    ISSN 0179-3632
 * --------------------------------------------------------------------------*/
/* clang-format on */

#include "ark_kepler.h"

#include <arkode/arkode.h>
#include <arkode/arkode_arkstep.h> /* prototypes for ARKStep fcts., consts */
#include <arkode/arkode_sprk.h>
#include <arkode/arkode_sprkstep.h> /* prototypes for SPRKStep fcts., consts */
#include <math.h>
#include <nvector/nvector_serial.h> /* serial N_Vector type, fcts., macros  */
#include <stdio.h>
#include <string.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_math.h> /* def. math fcns, 'sunrealtype'           */
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#define NUM_DT 8

typedef struct
{
  sunrealtype ecc;
}* UserData;

typedef struct
{
  N_Vector sol;
  sunrealtype energy_error;
  int method_order;
} ProblemResult;

/* Helper functions */
static int SolveProblem(ProgramArgs* args, ProblemResult* result,
                        SUNContext sunctx);
static void InitialConditions(N_Vector y0, sunrealtype ecc);
static sunrealtype Hamiltonian(N_Vector yvec);
static sunrealtype AngularMomentum(N_Vector y);

/* RHS callback functions */
static int dydt(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int velocity(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int force(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/* g(q) callback function for rootfinding */
static int rootfn(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data);

int SolveProblem(ProgramArgs* args, ProblemResult* result, SUNContext sunctx)
{
  void* arkode_mem       = NULL;
  N_Vector y             = NULL;
  SUNNonlinearSolver NLS = NULL;
  UserData udata         = NULL;
  sunrealtype* ydata     = NULL;
  sunrealtype num_orbits = 0;
  FILE* conserved_fp     = NULL;
  FILE* solution_fp      = NULL;
  FILE* times_fp         = NULL;
  int rootsfound         = 0;
  int iout               = 0;
  int retval             = 0;

  const int count_orbits     = args->count_orbits;
  const int step_mode        = args->step_mode;
  const int stepper          = args->stepper;
  const int use_compsums     = args->use_compsums;
  const int num_output_times = args->num_output_times;
  const char* method_name    = args->method_name;
  const sunrealtype dt       = args->dt;
  sunrealtype Tf             = args->tf;

  /* Default problem parameters */
  const sunrealtype T0    = SUN_RCONST(0.0);
  const sunrealtype dTout = (Tf - T0) / ((sunrealtype)num_output_times);
  const sunrealtype ecc   = SUN_RCONST(0.6);

  printf("\n   Begin Kepler Problem\n\n");
  PrintArgs(args);

  /* Allocate and fill udata structure */
  udata      = (UserData)malloc(sizeof(*udata));
  udata->ecc = ecc;

  /* Allocate our state vector */
  y     = N_VNew_Serial(4, sunctx);
  ydata = N_VGetArrayPointer(y);

  /* Fill the initial conditions */
  InitialConditions(y, ecc);

  /* Create SPRKStep integrator */
  if (stepper == 0)
  {
    arkode_mem = SPRKStepCreate(force, velocity, T0, y, sunctx);

    /* Optional: enable temporal root-finding */
    if (count_orbits)
    {
      ARKodeRootInit(arkode_mem, 1, rootfn);
      if (check_retval(&retval, "ARKodeRootInit", 1)) { return 1; }
    }

    retval = SPRKStepSetMethodName(arkode_mem, method_name);
    if (check_retval(&retval, "SPRKStepSetMethodName", 1)) { return 1; }

    retval = ARKodeSetUseCompensatedSums(arkode_mem, use_compsums);
    if (check_retval(&retval, "ARKodeSetUseCompensatedSums", 1)) { return 1; }

    if (step_mode == 0)
    {
      retval = ARKodeSetFixedStep(arkode_mem, dt);
      if (check_retval(&retval, "ARKodeSetFixedStep", 1)) { return 1; }

      retval = ARKodeSetMaxNumSteps(arkode_mem, ((long int)ceil(Tf / dt)) + 1);
      if (check_retval(&retval, "ARKodeSetMaxNumSteps", 1)) { return 1; }
    }
    else
    {
      fprintf(stderr,
              "ERROR: adaptive time-steps are not supported with SPRKStep\n");
      return 1;
    }

    retval = ARKodeSetUserData(arkode_mem, (void*)udata);
    if (check_retval(&retval, "ARKodeSetUserData", 1)) { return 1; }
  }
  else if (stepper == 1)
  {
    arkode_mem = ARKStepCreate(dydt, NULL, T0, y, sunctx);

    retval = ARKStepSetTableName(arkode_mem, "ARKODE_DIRK_NONE", method_name);
    if (check_retval(&retval, "ARKStepSetTableName", 1)) { return 1; }

    if (count_orbits)
    {
      ARKodeRootInit(arkode_mem, 1, rootfn);
      if (check_retval(&retval, "ARKodeRootInit", 1)) { return 1; }
    }

    retval = ARKodeSetUserData(arkode_mem, (void*)udata);
    if (check_retval(&retval, "ARKodeSetUserData", 1)) { return 1; }

    retval = ARKodeSetMaxNumSteps(arkode_mem, ((long int)ceil(Tf / dt)) + 1);
    if (check_retval(&retval, "ARKodeSetMaxNumSteps", 1)) { return 1; }

    if (step_mode == 0) { retval = ARKodeSetFixedStep(arkode_mem, dt); }
    else
    {
      retval = ARKodeSStolerances(arkode_mem, dt, dt);
      if (check_retval(&retval, "ARKodeSStolerances", 1)) { return 1; }
    }
  }

  /* Open output files */
  if (stepper == 0)
  {
    const char* fmt1 = "ark_kepler_conserved_%s-dt-%.2e.txt";
    const char* fmt2 = "ark_kepler_solution_%s-dt-%.2e.txt";
    const char* fmt3 = "ark_kepler_times_%s-dt-%.2e.txt";
    char fname[256];
    sprintf(fname, fmt1, method_name, dt);
    conserved_fp = fopen(fname, "w+");
    sprintf(fname, fmt2, method_name, dt);
    solution_fp = fopen(fname, "w+");
    sprintf(fname, fmt3, method_name, dt);
    times_fp = fopen(fname, "w+");
  }
  else
  {
    const char* fmt1 = "ark_kepler_conserved_%s-dt-%.2e.txt";
    const char* fmt2 = "ark_kepler_solution_%s-dt-%.2e.txt";
    const char* fmt3 = "ark_kepler_times_%s-dt-%.2e.txt";
    char fname[256];
    sprintf(fname, fmt1, method_name, dt);
    conserved_fp = fopen(fname, "w+");
    sprintf(fname, fmt2, method_name, dt);
    solution_fp = fopen(fname, "w+");
    sprintf(fname, fmt3, method_name, dt);
    times_fp = fopen(fname, "w+");
  }

  /* Print out starting energy, momentum before integrating */
  sunrealtype tret = T0;
  sunrealtype tout = T0 + dTout;
  sunrealtype H0   = Hamiltonian(y);
  sunrealtype L0   = AngularMomentum(y);
  fprintf(stdout, "t = %.4Lf, H(p,q) = %.16Lf, L(p,q) = %.16Lf\n",
          (long double)tret, (long double)H0, (long double)L0);
  fprintf(times_fp, "%.16Lf\n", (long double)tret);
  fprintf(conserved_fp, "%.16Lf, %.16Lf\n", (long double)H0, (long double)L0);
  N_VPrintFile(y, solution_fp);

  /* Do integration */
  if (stepper == 0)
  {
    while (iout < num_output_times)
    {
      /* Optional: if the stop time is not set, then its possible that the
         exact requested output time will not be hit (even with a fixed
         time-step due to roundoff error accumulation) and interpolation will be
         used to get the solution at the output time. */
      if (args->use_tstop) { ARKodeSetStopTime(arkode_mem, tout); }
      retval = ARKodeEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);

      if (retval == ARK_ROOT_RETURN)
      {
        num_orbits += SUN_RCONST(0.5);

        fprintf(stdout, "ROOT RETURN:\t");
        ARKodeGetRootInfo(arkode_mem, &rootsfound);
        fprintf(stdout, "  g[0] = %3d, y[0] = %3Lg, y[1] = %3Lg, num. orbits is now %.2Lf\n",
                rootsfound, (long double)ydata[0], (long double)ydata[1],
                (long double)num_orbits);
        fprintf(stdout, "t = %.4Lf, H(p,q)-H0 = %.16Le, L(p,q)-L0 = %.16Le\n",
                (long double)tret, (long double)(Hamiltonian(y) - H0),
                (long double)(AngularMomentum(y) - L0));
      }
      else if (retval >= 0)
      {
        /* Output current integration status */
        fprintf(stdout, "t = %.4Lf, H(p,q)-H0 = %.16Le, L(p,q)-L0 = %.16Le\n",
                (long double)tret, (long double)(Hamiltonian(y) - H0),
                (long double)(AngularMomentum(y) - L0));
        fprintf(times_fp, "%.16Lf\n", (long double)tret);
        fprintf(conserved_fp, "%.16Lf, %.16Lf\n", (long double)Hamiltonian(y),
                (long double)AngularMomentum(y));

        N_VPrintFile(y, solution_fp);

        tout += dTout;
        tout = (tout > Tf) ? Tf : tout;
        iout++;
      }
      else
      {
        fprintf(stderr, "Solver failure, stopping integration\n");
        break;
      }
    }
  }
  else
  {
    while (iout < num_output_times)
    {
      /* Optional: if the stop time is not set, then its possible that the the
         exact requested output time will not be hit (even with a fixed
         time-step due to roundoff error accumulation) and interpolation will be
         used to get the solution at the output time. */
      if (args->use_tstop) { ARKodeSetStopTime(arkode_mem, tout); }
      retval = ARKodeEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);

      if (retval == ARK_ROOT_RETURN)
      {
        num_orbits += SUN_RCONST(0.5);

        fprintf(stdout, "ROOT RETURN:\t");
        ARKodeGetRootInfo(arkode_mem, &rootsfound);
        fprintf(stdout, "  g[0] = %3d, y[0] = %3Lg, y[1] = %3Lg, num. orbits is now %.2Lf\n",
                rootsfound, (long double)ydata[0], (long double)ydata[1],
                (long double)num_orbits);
        fprintf(stdout, "t = %.4Lf, H(p,q)-H0 = %.16Le, L(p,q)-L0 = %.16Le\n",
                (long double)tret, (long double)(Hamiltonian(y) - H0),
                (long double)(AngularMomentum(y) - L0));
      }
      else if (retval >= 0)
      {
        /* Output current integration status */
        fprintf(stdout, "t = %.4Lf, H(p,q)-H0 = %.16Le, L(p,q)-L0 = %.16Le\n",
                (long double)tret, (long double)(Hamiltonian(y) - H0),
                (long double)(AngularMomentum(y) - L0));
        fprintf(times_fp, "%.16Lf\n", (long double)tret);
        fprintf(conserved_fp, "%.16Lf, %.16Lf\n", (long double)Hamiltonian(y),
                (long double)AngularMomentum(y));

        N_VPrintFile(y, solution_fp);

        tout += dTout;
        tout = (tout > Tf) ? Tf : tout;
        iout++;
      }
      else
      {
        fprintf(stderr, "Solver failure, stopping integration\n");
        break;
      }
    }
  }

  /* Copy results */
  N_VScale(SUN_RCONST(1.0), y, result->sol);
  result->energy_error = Hamiltonian(y) - H0;

  free(udata);
  fclose(times_fp);
  fclose(conserved_fp);
  fclose(solution_fp);
  if (NLS) { SUNNonlinSolFree(NLS); }
  N_VDestroy(y);
  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  ARKodeFree(&arkode_mem);
  return 0;
}

void InitialConditions(N_Vector y0vec, sunrealtype ecc)
{
  const sunrealtype zero = SUN_RCONST(0.0);
  const sunrealtype one  = SUN_RCONST(1.0);
  sunrealtype* y0        = N_VGetArrayPointer(y0vec);

  y0[0] = one - ecc;
  y0[1] = zero;
  y0[2] = zero;
  y0[3] = SUNRsqrt((one + ecc) / (one - ecc));
}

sunrealtype Hamiltonian(N_Vector yvec)
{
  sunrealtype H              = 0.0;
  sunrealtype* y             = N_VGetArrayPointer(yvec);
  const sunrealtype sqrt_qTq = SUNRsqrt(y[0] * y[0] + y[1] * y[1]);
  const sunrealtype pTp      = y[2] * y[2] + y[3] * y[3];

  H = SUN_RCONST(0.5) * pTp - SUN_RCONST(1.0) / sqrt_qTq;

  return H;
}

sunrealtype AngularMomentum(N_Vector yvec)
{
  sunrealtype L        = 0.0;
  sunrealtype* y       = N_VGetArrayPointer(yvec);
  const sunrealtype q1 = y[0];
  const sunrealtype q2 = y[1];
  const sunrealtype p1 = y[2];
  const sunrealtype p2 = y[3];

  L = q1 * p2 - q2 * p1;

  return L;
}

int dydt(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  int retval = 0;

  retval += force(t, yvec, ydotvec, user_data);
  retval += velocity(t, yvec, ydotvec, user_data);

  return retval;
}

int velocity(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  sunrealtype* y       = N_VGetArrayPointer(yvec);
  sunrealtype* ydot    = N_VGetArrayPointer(ydotvec);
  const sunrealtype p1 = y[2];
  const sunrealtype p2 = y[3];

  ydot[0] = p1;
  ydot[1] = p2;

  return 0;
}

int force(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  sunrealtype* y             = N_VGetArrayPointer(yvec);
  sunrealtype* ydot          = N_VGetArrayPointer(ydotvec);
  const sunrealtype q1       = y[0];
  const sunrealtype q2       = y[1];
  const sunrealtype sqrt_qTq = SUNRsqrt(q1 * q1 + q2 * q2);

  ydot[2] = -q1 / SUNRpowerR(sqrt_qTq, SUN_RCONST(3.0));
  ydot[3] = -q2 / SUNRpowerR(sqrt_qTq, SUN_RCONST(3.0));

  return 0;
}

int rootfn(sunrealtype t, N_Vector yvec, sunrealtype* gout, void* user_data)
{
  sunrealtype* y       = N_VGetArrayPointer(yvec);
  const sunrealtype q2 = y[1];

  gout[0] = q2;

  return 0;
}

int main(int argc, char* argv[])
{
  ProgramArgs args;
  ProblemResult result;
  SUNContext sunctx = NULL;
  int retval        = 0;

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return 1; }

  /* Parse the command line arguments */
  if (ParseArgs(argc, argv, &args)) { return 1; };

  /* Allocate space for result variables */
  result.sol = N_VNew_Serial(4, sunctx);

  if (!args.check_order)
  {
    /* SolveProblem calls a stepper to evolve the problem to Tf */
    retval = SolveProblem(&args, &result, sunctx);
    if (check_retval(&retval, "SolveProblem", 1)) { return 1; }
  }
  else
  {
    int i = 0;
    /* Compute the order of accuracy of the method by testing
       it with different step sizes. */
    sunrealtype acc_orders[NUM_DT];
    sunrealtype con_orders[NUM_DT];
    sunrealtype acc_errors[NUM_DT];
    sunrealtype con_errors[NUM_DT];
    ARKodeSPRKTable method = ARKodeSPRKTable_LoadByName(args.method_name);
    int expected_order     = method->q;
    N_Vector ref_sol       = N_VClone(result.sol);
    N_Vector error         = N_VClone(result.sol);
    sunrealtype a11 = 0, a12 = 0, a21 = 0, a22 = 0;
    sunrealtype b1 = 0, b2 = 0, b1e = 0, b2e = 0;
    sunrealtype ord_max_acc = 0, ord_max_conv = 0, ord_avg = 0, ord_est = 0;
    sunrealtype refine = SUN_RCONST(.5);
    sunrealtype dt = (expected_order >= 3) ? SUN_RCONST(1e-1) : SUN_RCONST(1e-3);
    sunrealtype dts[NUM_DT];

    /* Create a reference solution using 8th order ERK with a small time step */
    const int old_step_mode     = args.step_mode;
    const int old_stepper       = args.stepper;
    const char* old_method_name = args.method_name;
    args.dt                     = SUN_RCONST(1e-3);
    args.step_mode              = 0;
    args.stepper                = 1;
    args.method_name            = "ARKODE_ARK548L2SAb_ERK_8_4_5";

    /* Free method, we just needed it to get its order */
    ARKodeSPRKTable_Free(method);

    /* SolveProblem calls a stepper to evolve the problem to Tf */
    retval = SolveProblem(&args, &result, sunctx);
    if (check_retval(&retval, "SolveProblem", 1)) { return 1; }

    /* Store the reference solution */
    N_VScale(SUN_RCONST(1.0), result.sol, ref_sol);

    /* Restore the program args */
    args.step_mode   = old_step_mode;
    args.stepper     = old_stepper;
    args.method_name = old_method_name;

    for (i = 0; i < NUM_DT; i++) { dts[i] = dt * pow(refine, i); }

    /* Compute the error with various step sizes */
    for (i = 0; i < NUM_DT; i++)
    {
      /* Set the dt to use for this solve */
      args.dt = dts[i];

      /* SolveProblem calls a stepper to evolve the problem to Tf */
      retval = SolveProblem(&args, &result, sunctx);
      if (check_retval(&retval, "SolveProblem", 1)) { return 1; }

      printf("\n");

      /* Compute the error */
      N_VLinearSum(SUN_RCONST(1.0), result.sol, -SUN_RCONST(1.0), ref_sol, error);
      acc_errors[i] = SUNRsqrt(N_VDotProd(error, error)) /
                      ((sunrealtype)N_VGetLength(error));
      con_errors[i] = SUNRabs(result.energy_error);

      a11 += 1;
      a12 += log(dts[i]);
      a21 += log(dts[i]);
      a22 += (log(dts[i]) * log(dts[i]));
      b1 += log(acc_errors[i]);
      b2 += (log(acc_errors[i]) * log(dts[i]));
      b1e += log(con_errors[i]);
      b2e += (log(con_errors[i]) * log(dts[i]));

      if (i >= 1)
      {
        acc_orders[i - 1] = log(acc_errors[i] / acc_errors[i - 1]) /
                            log(dts[i] / dts[i - 1]);
        con_orders[i - 1] = log(con_errors[i] / con_errors[i - 1]) /
                            log(dts[i] / dts[i - 1]);
      }
    }

    /* Compute the order of accuracy */
    retval = ComputeConvergence(NUM_DT, acc_orders, expected_order, a11, a12, a21,
                                a22, b1, b2, &ord_avg, &ord_max_acc, &ord_est);
    printf("Order of accuracy wrt solution:    expected = %d, max = %.4Lf,  "
           "avg "
           "= %.4Lf,  "
           "overall = %.4Lf\n",
           expected_order, (long double)ord_max_acc, (long double)ord_avg,
           (long double)ord_est);

    /* Compute the order of accuracy with respect to conservation */
    retval = ComputeConvergence(NUM_DT, con_orders, expected_order, a11, a12,
                                a21, a22, b1e, b2e, &ord_avg, &ord_max_conv,
                                &ord_est);

    printf("Order of accuracy wrt Hamiltonian: expected = %d, max = %.4Lf,  "
           "avg = %.4Lf,  overall = %.4Lf\n",
           expected_order, (long double)ord_max_conv, (long double)ord_avg,
           (long double)ord_est);

    if (ord_max_acc < (expected_order - SUN_RCONST(0.5)))
    {
      printf(">>> FAILURE: computed order of accuracy wrt solution is below "
             "expected (%d)\n",
             expected_order);
      return 1;
    }

    if (ord_max_conv < (expected_order - SUN_RCONST(0.5)))
    {
      printf(">>> FAILURE: computed order of accuracy wrt Hamiltonian is below "
             "expected (%d)\n",
             expected_order);
      return 1;
    }

    N_VDestroy(ref_sol);
    N_VDestroy(error);
  }

  N_VDestroy(result.sol);
  SUNContext_Free(&sunctx);

  return 0;
}
