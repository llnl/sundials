/* -----------------------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
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
 * Unit tests on several ODEs with analytical solutions to verify the
 * ForcingStep module.
 * ---------------------------------------------------------------------------*/

#include <arkode/arkode_erkstep.h>
#include <arkode/arkode_forcingstep.h>
#include <nvector/nvector_serial.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#else
#define GSYM "g"
#endif

/* RHS function for f_1(t, y) = t / y */
static int f_forward_1(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VInv(y, ydot);
  N_VScale(t, ydot, ydot);
  return 0;
}

/* RHS function for f_2(t, y) = 1 / y */
static int f_forward_2(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VInv(y, ydot);
  return 0;
}

/* Integrates the ODE
 *
 * y' = [t / y] + [1 / y]
 *
 * with initial condition y(0) = 1 and partitioning specified by the square
 * brackets. We integrate to t = 1 and check the error against the exact
 * solution y(t) = |t + 1|.
 */
static int test_forward(SUNContext ctx)
{
  sunrealtype t0         = SUN_RCONST(0.0);
  sunrealtype tf         = SUN_RCONST(1.0);
  sunrealtype dt         = SUN_RCONST(1.0e-4);
  sunrealtype local_tol  = SUN_RCONST(1.0e-6);
  sunrealtype global_tol = SUN_RCONST(10.0) * local_tol;

  N_Vector y = N_VNew_Serial(1, ctx);
  /* Use the wrong initial condition for the partitions to ensure it gets reset
   * to the correct value of 1 by ForcingStep */
  N_VConst(SUN_RCONST(0.0), y);

  void* parititon_mem[] = {ERKStepCreate(f_forward_1, t0, y, ctx),
                           ERKStepCreate(f_forward_2, t0, y, ctx)};
  ARKodeSStolerances(parititon_mem[0], local_tol, local_tol);
  ARKodeSStolerances(parititon_mem[1], local_tol, local_tol);

  SUNStepper steppers[] = {NULL, NULL};
  ARKodeCreateSUNStepper(parititon_mem[0], &steppers[0]);
  ARKodeCreateSUNStepper(parititon_mem[1], &steppers[1]);

  N_VConst(SUN_RCONST(1.0), y);
  void* arkode_mem = ForcingStepCreate(steppers[0], steppers[1], t0, y, ctx);
  ARKodeSetFixedStep(arkode_mem, dt);
  ARKodeSetMaxNumSteps(arkode_mem, -1);

  sunrealtype tret = t0;
  ARKodeEvolve(arkode_mem, tf, y, &tret, ARK_NORMAL);

  sunrealtype exact_solution     = SUN_RCONST(2.0);
  sunrealtype numerical_solution = N_VGetArrayPointer(y)[0];
  sunrealtype err                = numerical_solution - exact_solution;

  printf("Forward direction solution completed with an error of %" GSYM "\n",
         err);
  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  sunbooleantype fail = SUNRCompareTol(exact_solution, numerical_solution,
                                       global_tol);
  if (fail)
  {
    fprintf(stderr, "Error exceeded tolerance of %" GSYM "\n", global_tol);
  }
  printf("\n");

  N_VDestroy(y);
  ARKodeFree(&parititon_mem[0]);
  SUNStepper_Destroy(&steppers[0]);
  ARKodeFree(&parititon_mem[1]);
  SUNStepper_Destroy(&steppers[1]);
  ARKodeFree(&arkode_mem);

  return fail;
}

static int f_mixed_direction_1(sunrealtype t, N_Vector z, N_Vector zdot,
                               void* data)
{
  N_VGetArrayPointer(zdot)[0] = N_VGetArrayPointer(z)[1] - t;
  N_VGetArrayPointer(zdot)[1] = SUN_RCONST(0.0);
  return 0;
}

static int f_mixed_direction_2(sunrealtype t, N_Vector z, N_Vector zdot,
                               void* data)
{
  N_VGetArrayPointer(zdot)[0] = SUN_RCONST(0.0);
  N_VGetArrayPointer(zdot)[1] = t - N_VGetArrayPointer(z)[0];
  return 0;
}

/* Integrates the ODE
 *
 * y_1' = y_2 - t
 * y_2' = t - y_1
 *
 * with initial condition y(0) = [1, 1]^T backward in time to t = -1, then
 * forward to t = 0.4, and backward to t = 0. We apply a splitting method using
 * a component partitioning and check that the numerical solution is close to
 * the original initial condition.
 */
static int test_mixed_directions(SUNContext ctx)
{
  sunrealtype t0         = SUN_RCONST(0.0);
  sunrealtype t1         = SUN_RCONST(-1.0);
  sunrealtype t2         = SUN_RCONST(0.4);
  sunrealtype t3         = t0;
  sunrealtype dt         = -SUN_RCONST(1.23e-4);
  sunrealtype local_tol  = SUN_RCONST(1.0e-4);
  sunrealtype global_tol = SUN_RCONST(10.0) * local_tol;
  N_Vector y             = N_VNew_Serial(2, ctx);
  N_VConst(SUN_RCONST(1.0), y);
  N_Vector err = N_VClone(y);
  N_VConst(SUN_RCONST(1.0), err);

  void* parititon_mem[] = {ERKStepCreate(f_mixed_direction_1, t0, y, ctx),
                           ERKStepCreate(f_mixed_direction_2, t0, y, ctx)};
  ARKodeSStolerances(parititon_mem[0], local_tol, local_tol);
  ARKodeSStolerances(parititon_mem[1], local_tol, local_tol);

  SUNStepper steppers[] = {NULL, NULL};
  ARKodeCreateSUNStepper(parititon_mem[0], &steppers[0]);
  ARKodeCreateSUNStepper(parititon_mem[1], &steppers[1]);

  void* arkode_mem = ForcingStepCreate(steppers[0], steppers[1], t0, y, ctx);
  ARKodeSetFixedStep(arkode_mem, dt);
  ARKodeSetInterpolantType(arkode_mem, ARK_INTERP_HERMITE);
  ARKodeSetMaxNumSteps(arkode_mem, -1);

  /* Integrate from 0 to -1 */
  sunrealtype tret = t0;
  ARKodeEvolve(arkode_mem, t1, y, &tret, ARK_NORMAL);

  /* Integrate from -1 to 0.4 */
  ARKodeReset(arkode_mem, t1, y);
  ARKodeSetStepDirection(arkode_mem, t2 - t1);
  ARKodeEvolve(arkode_mem, t2, y, &tret, ARK_NORMAL);

  /* Integrate from 0.4 to 0 */
  ARKodeReset(arkode_mem, t2, y);
  ARKodeSetStepDirection(arkode_mem, t3 - t2);
  ARKodeEvolve(arkode_mem, t3, y, &tret, ARK_NORMAL);

  N_VLinearSum(SUN_RCONST(1.0), err, -SUN_RCONST(1.0), y, err);
  sunrealtype max_err = N_VMaxNorm(err);

  printf("Mixed direction solution completed with an error of %" GSYM "\n",
         max_err);
  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  sunbooleantype fail = max_err > global_tol;
  if (fail)
  {
    fprintf(stderr, "Error exceeded tolerance of %" GSYM "\n", global_tol);
  }
  printf("\n");

  N_VDestroy(y);
  N_VDestroy(err);
  ARKodeFree(&parititon_mem[0]);
  SUNStepper_Destroy(&steppers[0]);
  ARKodeFree(&parititon_mem[1]);
  SUNStepper_Destroy(&steppers[1]);
  ARKodeFree(&arkode_mem);

  return fail;
}

/* Integrates the ODE
 *
 * y' = [t / y] + [1 / y]
 *
 * with initial condition y(0) = 1 and partitioning specified by the square
 * brackets. We integrate to t = 1 then reinitialize the forcing method by
 * swapping the SUNSteppers and updating the initial condition to y(1) = -1.
 * Next we integrate to t = 2 and check the error against the exact solution
 * y(2) = -sqrt(6).
 */
static int test_reinit(SUNContext ctx)
{
  sunrealtype t0         = SUN_RCONST(0.0);
  sunrealtype t1         = SUN_RCONST(1.0);
  sunrealtype t2         = SUN_RCONST(2.0);
  sunrealtype dt         = SUN_RCONST(1.0e-4);
  sunrealtype local_tol  = SUN_RCONST(1.0e-6);
  sunrealtype global_tol = SUN_RCONST(10.0) * local_tol;

  N_Vector y = N_VNew_Serial(1, ctx);
  N_VConst(SUN_RCONST(1.0), y);

  void* parititon_mem[] = {ERKStepCreate(f_forward_1, t0, y, ctx),
                           ERKStepCreate(f_forward_2, t0, y, ctx)};
  ARKodeSStolerances(parititon_mem[0], local_tol, local_tol);
  ARKodeSStolerances(parititon_mem[1], local_tol, local_tol);

  SUNStepper steppers[] = {NULL, NULL};
  ARKodeCreateSUNStepper(parititon_mem[0], &steppers[0]);
  ARKodeCreateSUNStepper(parititon_mem[1], &steppers[1]);

  void* arkode_mem = ForcingStepCreate(steppers[0], steppers[1], t0, y, ctx);
  ARKodeSetFixedStep(arkode_mem, dt);
  ARKodeSetMaxNumSteps(arkode_mem, -1);

  sunrealtype tret = t0;
  ARKodeEvolve(arkode_mem, t1, y, &tret, ARK_NORMAL);

  /* Change the state and swap the steppers so now forcing is applied to
   * steppers[0] rather than steppers[1] */
  N_VConst(SUN_RCONST(-1.0), y);
  ERKStepReInit(parititon_mem[0], f_forward_1, t1, y);
  ForcingStepReInit(arkode_mem, steppers[1], steppers[0], t1, y);
  ERKStepReInit(parititon_mem[1], f_forward_2, t1, y);

  ARKodeEvolve(arkode_mem, t2, y, &tret, ARK_NORMAL);

  sunrealtype exact_solution     = -SUNRsqrt(SUN_RCONST(6.0));
  sunrealtype numerical_solution = N_VGetArrayPointer(y)[0];
  sunrealtype err                = numerical_solution - exact_solution;

  printf("Reinitialized solution completed with an error of %" GSYM "\n", err);
  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  sunbooleantype fail = SUNRCompareTol(exact_solution, numerical_solution,
                                       global_tol);
  if (fail)
  {
    fprintf(stderr, "Error exceeded tolerance of %" GSYM "\n", global_tol);
  }
  printf("\n");

  N_VDestroy(y);
  ARKodeFree(&parititon_mem[0]);
  SUNStepper_Destroy(&steppers[0]);
  ARKodeFree(&parititon_mem[1]);
  SUNStepper_Destroy(&steppers[1]);
  ARKodeFree(&arkode_mem);

  return fail;
}

int main(void)
{
  SUNContext ctx = NULL;
  SUNErrCode err = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (err != SUN_SUCCESS)
  {
    fprintf(stderr, "Failed to create the SUNContext\n");
    return 1;
  }

  err = SUNContext_PushErrHandler(ctx, SUNAbortErrHandlerFn, NULL);
  if (err != SUN_SUCCESS)
  {
    fprintf(stderr, "Failed to add error handler\n");
    return 1;
  }

  int errors = 0;
  errors += test_forward(ctx);
  errors += test_mixed_directions(ctx);
  errors += test_reinit(ctx);

  SUNContext_Free(&ctx);

  if (errors == 0) { printf("Success\n"); }
  else { printf("%d Test Failures\n", errors); }

  return 0;
}
