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
 * SplittingStep module.
 * ---------------------------------------------------------------------------*/

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <arkode/arkode_splittingstep.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.hpp>
#include <sundials/sundials_stepper.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#include <cmath>
#include <iostream>
#include <vector>

/* Integrates the ODE
 *
 * y' = \sum_{i=0}^{P-1} 2^i / (1 - 2^P) * y,    y(0) = 1
 *
 * with a 3rd order operator splitting method. Each partition is solved by an
 * ERK method. Using the exact solution y(t) = exp(-t), we confirm the numerical
 * solution is sufficiently accurate.
 */
static int test_forward(sundials::Context& ctx, int order, int partitions)
{
  constexpr auto t0         = SUN_RCONST(0.0);
  constexpr auto tf         = SUN_RCONST(1.0);
  constexpr auto dt         = SUN_RCONST(8.0e-3);
  constexpr auto local_tol  = SUN_RCONST(1.0e-6);
  constexpr auto global_tol = SUN_RCONST(10.0) * local_tol;
  auto y                    = N_VNew_Serial(1, ctx);
  N_VConst(SUN_RCONST(1.0), y);

  ARKRhsFn f = [](sunrealtype, N_Vector z, N_Vector zdot, void* user_data)
  {
    auto lambda = *static_cast<sunrealtype*>(user_data);
    N_VScale(lambda, z, zdot);
    return 0;
  };

  std::vector<void*> partition_mem(partitions);
  std::vector<sunrealtype> lambda(partitions);
  std::vector<SUNStepper> steppers(partitions);
  for (int i = 0; i < partitions; i++)
  {
    partition_mem[i] = ERKStepCreate(f, t0, y, ctx);
    /* The lambdas sum up to 1 */
    lambda[i] = std::pow(SUN_RCONST(2.0), i) /
                (1 - std::pow(SUN_RCONST(2.0), partitions));
    ARKodeSetUserData(partition_mem[i], &lambda[i]);
    ARKodeSStolerances(partition_mem[i], local_tol, local_tol);
    ARKodeCreateSUNStepper(partition_mem[i], &steppers[i]);
  }

  auto arkode_mem = SplittingStepCreate(steppers.data(), partitions, t0, y, ctx);
  ARKodeSetFixedStep(arkode_mem, dt);
  ARKodeSetOrder(arkode_mem, order);
  auto tret = t0;
  ARKodeEvolve(arkode_mem, tf, y, &tret, ARK_NORMAL);

  /* Check that the solution matches the exact solution */
  auto exact_solution     = std::exp(t0 - tf);
  auto numerical_solution = N_VGetArrayPointer(y)[0];
  auto err                = numerical_solution - exact_solution;

  std::cout << "Forward solution of order " << order << " with " << partitions
            << " partitions completed with an error of " << err << "\n";
  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  sunbooleantype fail = SUNRCompareTol(exact_solution, numerical_solution,
                                       global_tol);
  if (fail)
  {
    std::cerr << "Error exceeded tolerance of " << global_tol << "\n";
  }
  std::cout << "\n";

  N_VDestroy(y);
  for (int i = 0; i < partitions; i++)
  {
    ARKodeFree(&partition_mem[i]);
    SUNStepper_Destroy(&steppers[i]);
  }
  ARKodeFree(&arkode_mem);

  return fail;
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
static int test_mixed_directions(const sundials::Context& ctx, const char* name)
{
  constexpr auto t0         = SUN_RCONST(0.0);
  constexpr auto t1         = SUN_RCONST(-1.0);
  constexpr auto t2         = SUN_RCONST(0.4);
  constexpr auto t3         = t0;
  constexpr auto dt         = -SUN_RCONST(1.23e-3);
  constexpr auto local_tol  = SUN_RCONST(1.0e-6);
  constexpr auto global_tol = SUN_RCONST(25.0) * local_tol;
  auto y                    = N_VNew_Serial(2, ctx);
  N_VConst(SUN_RCONST(1.0), y);
  auto err = N_VClone(y);
  N_VConst(SUN_RCONST(1.0), err);

  ARKRhsFn f1 = [](sunrealtype t, N_Vector z, N_Vector zdot, void*)
  {
    N_VGetArrayPointer(zdot)[0] = N_VGetArrayPointer(z)[1] - t;
    N_VGetArrayPointer(zdot)[1] = SUN_RCONST(0.0);
    return 0;
  };

  ARKRhsFn f2 = [](sunrealtype t, N_Vector z, N_Vector zdot, void*)
  {
    N_VGetArrayPointer(zdot)[0] = SUN_RCONST(0.0);
    N_VGetArrayPointer(zdot)[1] = t - N_VGetArrayPointer(z)[0];
    return 0;
  };

  void* parititon_mem[] = {ARKStepCreate(nullptr, f1, t0, y, ctx),
                           ERKStepCreate(f2, t0, y, ctx)};
  ARKodeSStolerances(parititon_mem[0], local_tol, local_tol);
  ARKodeSStolerances(parititon_mem[1], local_tol, local_tol);

  SUNStepper steppers[] = {nullptr, nullptr};
  ARKodeCreateSUNStepper(parititon_mem[0], &steppers[0]);
  ARKodeCreateSUNStepper(parititon_mem[1], &steppers[1]);

  SUNNonlinearSolver nls = SUNNonlinSol_FixedPoint(y, 1, ctx);
  ARKodeSetNonlinearSolver(parititon_mem[0], nls);

  auto arkode_mem = SplittingStepCreate(steppers, 2, t0, y, ctx);
  ARKodeSetFixedStep(arkode_mem, dt);
  ARKodeSetInterpolantType(arkode_mem, ARK_INTERP_HERMITE);
  ARKodeSetMaxNumSteps(arkode_mem, -1);
  auto coefficients = SplittingStepCoefficients_LoadCoefficientsByName(name);
  SplittingStepSetCoefficients(arkode_mem, coefficients);
  SplittingStepCoefficients_Destroy(&coefficients);

  /* Integrate from 0 to -1 */
  auto tret = t0;
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
  auto max_err = N_VMaxNorm(err);

  std::cout << "Mixed direction solution using " << name
            << " completed with an error of " << max_err << "\n";
  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  sunbooleantype fail = max_err > global_tol;
  if (fail)
  {
    std::cerr << "Error exceeded tolerance of " << global_tol << "\n";
  }
  std::cout << "\n";

  N_VDestroy(y);
  N_VDestroy(err);
  ARKodeFree(&parititon_mem[0]);
  SUNStepper_Destroy(&steppers[0]);
  ARKodeFree(&parititon_mem[1]);
  SUNStepper_Destroy(&steppers[1]);
  ARKodeFree(&arkode_mem);
  SUNNonlinSolFree(nls);

  return fail;
}

/* Integrates the ODE
 *
 * y' = [y^2] - [y*(y + 1)] = -y
 *
 * with initial condition y(0) = 1 and partitioning specified by the square
 * brackets. Powers and products are done componentwise. At t = 0.5, we resize
 * from 1D to 2D:
 *
 * y_new = [1, y_old]^T
 *
 * Then we integrate to t = 1 and check the error against the exact solution
 * y(1) = [exp(-0.5), exp(-1)].
 */
static int test_resize(const sundials::Context& ctx)
{
  constexpr auto t0         = SUN_RCONST(0.0);
  constexpr auto t1         = SUN_RCONST(0.5);
  constexpr auto t2         = SUN_RCONST(1.0);
  constexpr auto dt         = SUN_RCONST(8.0e-3);
  constexpr auto local_tol  = SUN_RCONST(1.0e-5);
  constexpr auto global_tol = local_tol;
  auto y                    = N_VNew_Serial(1, ctx);
  N_VConst(SUN_RCONST(1.0), y);

  ARKRhsFn f1 = [](sunrealtype t, N_Vector z, N_Vector zdot, void*)
  {
    N_VProd(z, z, zdot);
    return 0;
  };

  ARKRhsFn f2 = [](sunrealtype t, N_Vector z, N_Vector zdot, void*)
  {
    N_VConst(-SUN_RCONST(1.0), zdot);
    N_VLinearSum(SUN_RCONST(1.0), zdot, -SUN_RCONST(1.0), z, zdot);
    N_VProd(zdot, z, zdot);
    return 0;
  };

  void* parititon_mem[] = {ERKStepCreate(f1, t0, y, ctx),
                           ERKStepCreate(f2, t0, y, ctx)};
  ARKodeSStolerances(parititon_mem[0], local_tol, local_tol);
  ARKodeSStolerances(parititon_mem[1], local_tol, local_tol);

  SUNStepper steppers[] = {nullptr, nullptr};
  ARKodeCreateSUNStepper(parititon_mem[0], &steppers[0]);
  ARKodeCreateSUNStepper(parititon_mem[1], &steppers[1]);

  auto arkode_mem = SplittingStepCreate(steppers, 2, t0, y, ctx);
  ARKodeSetFixedStep(arkode_mem, dt);
  auto coefficients = SplittingStepCoefficients_SymmetricParallel(2);
  SplittingStepSetCoefficients(arkode_mem, coefficients);
  SplittingStepCoefficients_Destroy(&coefficients);

  /* Integrate from 0 to 0.5 */
  auto tret = t0;
  ARKodeEvolve(arkode_mem, t1, y, &tret, ARK_NORMAL);

  /* Resize */
  auto y_new                   = N_VNew_Serial(2, ctx);
  N_VGetArrayPointer(y_new)[0] = SUN_RCONST(1.0);
  N_VGetArrayPointer(y_new)[1] = N_VGetArrayPointer(y)[0];
  N_VDestroy(y);
  ARKodeResize(arkode_mem, y_new, SUN_RCONST(1.0), t1, nullptr, nullptr);
  ARKodeResize(parititon_mem[0], y_new, SUN_RCONST(1.0), t1, nullptr, nullptr);
  ARKodeResize(parititon_mem[1], y_new, SUN_RCONST(1.0), t1, nullptr, nullptr);

  /* Integrate from 0.5 to 1 */
  ARKodeEvolve(arkode_mem, t2, y_new, &tret, ARK_NORMAL);

  auto err                   = N_VClone(y_new);
  N_VGetArrayPointer(err)[0] = std::exp(t1 - t2);
  N_VGetArrayPointer(err)[1] = std::exp(t0 - t2);
  N_VLinearSum(SUN_RCONST(1.0), err, -SUN_RCONST(1.0), y_new, err);
  auto max_err = N_VMaxNorm(err);

  std::cout << "Resized solution completed with an error of " << max_err << "\n";
  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  sunbooleantype fail = max_err > global_tol;
  if (fail)
  {
    std::cerr << "Error exceeded tolerance of " << global_tol << "\n";
  }
  std::cout << "\n";

  N_VDestroy(y_new);
  N_VDestroy(err);
  ARKodeFree(&parititon_mem[0]);
  SUNStepper_Destroy(&steppers[0]);
  ARKodeFree(&parititon_mem[1]);
  SUNStepper_Destroy(&steppers[1]);
  ARKodeFree(&arkode_mem);

  return fail;
}

/* Creates a custom SUNStepper for the linear, scalar ODE y' = lambda*y */
static SUNStepper create_exp_stepper(const sundials::Context& ctx,
                                     const sunrealtype& lam, N_Vector y)
{
  SUNStepper stepper = nullptr;
  SUNStepper_Create(ctx, &stepper);

  struct Content
  {
    sunrealtype lambda;
    sunrealtype t{};
    N_Vector v;

    Content(sunrealtype l, N_Vector tmpl) : lambda(l), v(N_VClone(tmpl)) {}

    ~Content() { N_VDestroy(v); }

    static Content& from_stepper(SUNStepper s)
    {
      void* content = nullptr;
      SUNStepper_GetContent(s, &content);
      return *static_cast<Content*>(content);
    }
  };

  SUNStepper_SetContent(stepper, new Content(lam, y));

  auto reset = [](SUNStepper s, sunrealtype tR, N_Vector vR)
  {
    auto& content = Content::from_stepper(s);
    content.t     = tR;
    N_VScale(SUN_RCONST(1.0), vR, content.v);
    return 0;
  };
  SUNStepper_SetResetFn(stepper, reset);

  auto empty_func = [](auto...) { return 0; };
  SUNStepper_SetStopTimeFn(stepper, empty_func);
  SUNStepper_SetStepDirectionFn(stepper, empty_func);
  SUNStepper_SetFullRhsFn(stepper, empty_func);

  auto evolve =
    [](SUNStepper s, sunrealtype tout, N_Vector vret, sunrealtype* tret)
  {
    auto& content = Content::from_stepper(s);
    N_VScale(std::exp(content.lambda * (tout - content.t)), content.v, vret);
    *tret = tout;
    return 0;
  };
  SUNStepper_SetEvolveFn(stepper, evolve);

  auto destroy = [](SUNStepper s)
  {
    delete &Content::from_stepper(s);
    return 0;
  };
  SUNStepper_SetDestroyFn(stepper, destroy);
  return stepper;
}

/* Integrates the ODE
 *
 * y' = -0.6*y - 0.4*y
 *
 * with initial condition y(0) = 1 with a sixth order splitting method and
 * exact, custom SUNSteppers for the two partitions. We integrate to t = 1 and
 * compare the numerical solution to the exact solution y(t) = exp(-t).
 */
static int test_custom_stepper(const sundials::Context& ctx, int order)
{
  constexpr auto lambda1 = SUN_RCONST(-0.6);
  constexpr auto lambda2 = SUN_RCONST(-0.4);
  constexpr auto t0      = SUN_RCONST(0.0);
  constexpr auto tf      = SUN_RCONST(1.0);
  constexpr auto dt      = SUN_RCONST(0.1);
  auto y                 = N_VNew_Serial(1, ctx);
  N_VConst(SUN_RCONST(1.0), y);

  SUNStepper steppers[] = {create_exp_stepper(ctx, lambda1, y),
                           create_exp_stepper(ctx, lambda2, y)};

  auto arkode_mem = SplittingStepCreate(steppers, 2, t0, y, ctx);
  ARKodeSetFixedStep(arkode_mem, dt);
  auto coefficients = SplittingStepCoefficients_SuzukiFractal(2, order);
  SplittingStepSetCoefficients(arkode_mem, coefficients);
  SplittingStepCoefficients_Destroy(&coefficients);

  auto tret = t0;
  ARKodeEvolve(arkode_mem, tf, y, &tret, ARK_NORMAL);

  /* Check that the solution matches the exact solution */
  auto exact_solution     = std::exp(t0 - tf);
  auto numerical_solution = N_VGetArrayPointer(y)[0];
  auto err                = numerical_solution - exact_solution;

  std::cout << "Custom SUNStepper solution of order " << order
            << " completed with an error of " << err << "\n";
  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  sunbooleantype fail = SUNRCompare(exact_solution, numerical_solution);
  if (fail) { std::cerr << "Error exceeded tolerance\n"; }
  std::cout << "\n";

  N_VDestroy(y);
  SUNStepper_Destroy(&steppers[0]);
  SUNStepper_Destroy(&steppers[1]);
  ARKodeFree(&arkode_mem);
  return fail;
}

/* Integrates the ODE
 *
 * y' = [y / (1 + t)] + [-y / (2 + t)] + [y / (3 + t)]
 *
 * with initial condition y(0) = 1 and partitioning specified by the square
 * brackets. After integration to t = 1, the third partition and its
 * corresponding SUNStepper are dropped from the ODE so there are only two
 * partitions. Then we integrate to t = 2 and check the error against the exact
 * solution y(2) = 2.
 */
static int test_reinit(const sundials::Context& ctx)
{
  constexpr auto t0         = SUN_RCONST(0.0);
  constexpr auto t1         = SUN_RCONST(1.0);
  constexpr auto t2         = SUN_RCONST(2.0);
  constexpr auto dt         = SUN_RCONST(8.0e-3);
  constexpr auto local_tol  = SUN_RCONST(1.0e-6);
  constexpr auto global_tol = SUN_RCONST(10.0) * local_tol;
  auto y                    = N_VNew_Serial(1, ctx);
  N_VConst(SUN_RCONST(1.0), y);

  ARKRhsFn fns[] = {[](sunrealtype t, N_Vector z, N_Vector zdot, void*)
                    {
                      N_VScale(SUN_RCONST(1.0) / (SUN_RCONST(1.0) + t), z, zdot);
                      return 0;
                    },
                    [](sunrealtype t, N_Vector z, N_Vector zdot, void*)
                    {
                      N_VScale(-SUN_RCONST(1.0) / (SUN_RCONST(2.0) + t), z, zdot);
                      return 0;
                    },
                    [](sunrealtype t, N_Vector z, N_Vector zdot, void*)
                    {
                      N_VScale(SUN_RCONST(1.0) / (SUN_RCONST(3.0) + t), z, zdot);
                      return 0;
                    }};

  void* partition_mem[] = {nullptr, nullptr, nullptr};
  SUNStepper steppers[] = {nullptr, nullptr, nullptr};

  for (std::size_t i = 0; i < 3; i++)
  {
    partition_mem[i] = ERKStepCreate(fns[i], t0, y, ctx);
    ARKodeSStolerances(partition_mem[i], local_tol, local_tol);
    ARKodeCreateSUNStepper(partition_mem[i], &steppers[i]);
  }

  auto arkode_mem = SplittingStepCreate(steppers, 3, t0, y, ctx);
  ARKodeSetFixedStep(arkode_mem, dt);
  ARKodeSetOrder(arkode_mem, 2);
  auto tret = t0;
  ARKodeEvolve(arkode_mem, t1, y, &tret, ARK_NORMAL);

  SplittingStepReInit(arkode_mem, steppers, 2, t1, y);
  for (std::size_t i = 0; i < 2; i++)
  {
    ERKStepReInit(partition_mem[i], fns[i], t1, y);
  }

  ARKodeEvolve(arkode_mem, t2, y, &tret, ARK_NORMAL);

  /* Check that the solution matches the exact solution */
  auto exact_solution     = SUN_RCONST(2.0);
  auto numerical_solution = N_VGetArrayPointer(y)[0];
  auto err                = numerical_solution - exact_solution;

  std::cout << "Reinitialized solution completed with an error of " << err
            << "\n";
  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  sunbooleantype fail = SUNRCompareTol(exact_solution, numerical_solution,
                                       global_tol);
  if (fail)
  {
    std::cerr << "Error exceeded tolerance of " << global_tol << "\n";
  }
  std::cout << "\n";

  N_VDestroy(y);
  for (std::size_t i = 0; i < 3; i++)
  {
    ARKodeFree(&partition_mem[i]);
    SUNStepper_Destroy(&steppers[i]);
  }
  ARKodeFree(&arkode_mem);

  return fail;
}

int main()
{
  sundials::Context ctx;
  SUNContext_PushErrHandler(ctx, SUNAbortErrHandlerFn, nullptr);

  int errors = 0;

  /* Run the tests */
  constexpr auto min_partitions = 2;
  constexpr auto max_partitions = 5;
  for (auto p = min_partitions; p <= max_partitions; p++)
  {
    constexpr auto min_order = 1;
    constexpr auto max_order = 4;
    for (auto order = min_order; order <= max_order; order++)
    {
      errors += test_forward(ctx, order, p);
    }
  }

  auto names = {"ARKODE_SPLITTING_STRANG_2_2_2", "ARKODE_SPLITTING_BEST_2_2_2",
                "ARKODE_SPLITTING_SUZUKI_3_3_2", "ARKODE_SPLITTING_RUTH_3_3_2"};
  for (auto name : names) { errors += test_mixed_directions(ctx, name); }

  errors += test_resize(ctx);
  errors += test_custom_stepper(ctx, 4);
  errors += test_custom_stepper(ctx, 6);
  errors += test_reinit(ctx);

  if (errors == 0) { std::cout << "Success\n"; }
  else { std::cout << errors << " Test Failures\n"; }

  return errors;
}
