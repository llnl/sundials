/* -----------------------------------------------------------------------------
 * Programmer(s): Demetrios Agourakis
 * Created: 2026-10-31
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2026, the SUNDIALS contributors.
 * All rights reserved.
 * See the top-level LICENSE and NOTICE files for details.
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Exercise residual failures after state convergence with staggered FSA.
 * F = y' + p*y, p = 1, y(0) = 1; s = dy/dp = -t*exp(-t).
 * Optional q' = y and qS' = s exercise auxiliary-history restoration.
 * ---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

#include "idas/idas.h"
#include "idas/idas_ls.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunnonlinsol/sunnonlinsol_newton.h"

namespace {
enum class Failure
{
  None,
  Once,
  Persistent,
  Fatal
};

using SolveFn = int (*)(SUNNonlinearSolver, N_Vector, N_Vector, N_Vector,
                        sunrealtype, sunbooleantype, void*);

struct Probe
{
  void* memory              = nullptr;
  SolveFn solve             = nullptr;
  Failure failure           = Failure::None;
  int inject_at             = 8;
  int converged_calls       = 0;
  int rejections            = 0;
  long int injected_nst     = -1;
  bool after_state_solve    = false;
  bool watchdog             = false;
  sunrealtype rejected_time = SUN_RCONST(0.0);
};

/* Observe the public Newton solver without changing its computations. The next
   state residual after a successful state solve is the callback under test. */
int observe_solve(SUNNonlinearSolver nls, N_Vector y0, N_Vector ycor,
                  N_Vector weights, sunrealtype tolerance,
                  sunbooleantype call_setup, void* memory)
{
  void* user_data = nullptr;
  if (IDAGetUserData(memory, &user_data) != IDA_SUCCESS) { return -1; }
  auto* probe              = static_cast<Probe*>(user_data);
  probe->after_state_solve = false;
  const int flag = probe->solve(nls, y0, ycor, weights, tolerance, call_setup,
                                memory);
  probe->after_state_solve = flag == SUN_SUCCESS;
  return flag;
}

// Inject the selected failure only after the state nonlinear solve converges.
int residual(sunrealtype t, N_Vector y, N_Vector yp, N_Vector r, void* user_data)
{
  auto* probe = static_cast<Probe*>(user_data);
  if (probe->after_state_solve)
  {
    probe->after_state_solve = false;
    ++probe->converged_calls;
    const bool selected = probe->converged_calls >= probe->inject_at;
    const bool reject   = selected && (probe->failure == Failure::Persistent ||
                                     ((probe->failure == Failure::Once ||
                                       probe->failure == Failure::Fatal) &&
                                      probe->rejections == 0));
    if (reject)
    {
      if (probe->rejections == 0)
      {
        probe->rejected_time = t;
        if (IDAGetNumSteps(probe->memory, &probe->injected_nst) != IDA_SUCCESS)
        {
          return -1;
        }
      }
      ++probe->rejections;
      /* Bound the broken implementation without using a wall-clock timeout. */
      if (probe->rejections > 32)
      {
        probe->watchdog = true;
        return -1;
      }
      return probe->failure == Failure::Fatal ? -7 : 7;
    }
  }
  NV_Ith_S(r, 0) = NV_Ith_S(yp, 0) + NV_Ith_S(y, 0);
  return 0;
}

// Evaluate the parameter sensitivity residual for y' + p*y = 0 at p = 1.
int sensitivity(int, sunrealtype, N_Vector y, N_Vector, N_Vector, N_Vector* s,
                N_Vector* sp, N_Vector* r, void*, N_Vector, N_Vector, N_Vector)
{
  NV_Ith_S(r[0], 0) = NV_Ith_S(sp[0], 0) + NV_Ith_S(s[0], 0) + NV_Ith_S(y, 0);
  return 0;
}

// Provide the exact scalar Jacobian of the state residual.
int jacobian(sunrealtype, sunrealtype cj, N_Vector, N_Vector, N_Vector,
             SUNMatrix matrix, void*, N_Vector, N_Vector, N_Vector)
{
  SM_ELEMENT_D(matrix, 0, 0) = cj + SUN_RCONST(1.0);
  return 0;
}

// Integrate the state as an auxiliary variable.
int quadrature(sunrealtype, N_Vector y, N_Vector, N_Vector qdot, void*)
{
  NV_Ith_S(qdot, 0) = NV_Ith_S(y, 0);
  return 0;
}

// Integrate the sensitivity of the auxiliary variable.
int quad_sensitivity(int, sunrealtype, N_Vector, N_Vector, N_Vector* s,
                     N_Vector*, N_Vector, N_Vector* qsdot, void*, N_Vector,
                     N_Vector, N_Vector)
{
  NV_Ith_S(qsdot[0], 0) = NV_Ith_S(s[0], 0);
  return 0;
}

struct Resources
{
  SUNContext context = nullptr;
  void* memory       = nullptr;
  N_Vector y = nullptr, yp = nullptr, s = nullptr, sp = nullptr;
  N_Vector q = nullptr, qs = nullptr;
  SUNMatrix matrix             = nullptr;
  SUNLinearSolver linear       = nullptr;
  SUNNonlinearSolver nonlinear = nullptr;

  ~Resources()
  {
    if (memory) { IDAFree(&memory); }
    if (nonlinear) { SUNNonlinSolFree(nonlinear); }
    if (linear) { SUNLinSolFree(linear); }
    if (matrix) { SUNMatDestroy(matrix); }
    for (N_Vector v : {qs, q, sp, s, yp, y})
    {
      if (v) { N_VDestroy(v); }
    }
    if (context) { SUNContext_Free(&context); }
  }
};

struct Case
{
  const char* name;
  int method;
  Failure failure;
  int inject_at;
  bool auxiliary;
};

bool run_case(const Case& test)
{
  Probe probe;
  Resources r;
  probe.failure   = test.failure;
  probe.inject_at = test.inject_at;
  auto check_flag = [](int flag, const char* routine)
  {
    if (flag == IDA_SUCCESS) { return true; }
    std::fprintf(stderr, "%s returned %d\n", routine, flag);
    return false;
  };
  auto check_ptr = [](const void* pointer, const char* routine)
  {
    if (pointer) { return true; }
    std::fprintf(stderr, "%s returned NULL\n", routine);
    return false;
  };
  // Create the context and suppress expected solver-failure diagnostics.
  if (!check_flag(SUNContext_Create(SUN_COMM_NULL, &r.context),
                  "SUNContext_Create") ||
      !check_flag(SUNContext_ClearErrHandlers(r.context),
                  "SUNContext_ClearErrHandlers"))
  {
    return false;
  }
  // Allocate state, sensitivity, quadrature, and solver resources.
  r.y = N_VNew_Serial(1, r.context);
  if (!check_ptr(r.y, "N_VNew_Serial")) { return false; }
  r.yp     = N_VClone(r.y);
  r.s      = N_VClone(r.y);
  r.sp     = N_VClone(r.y);
  r.q      = N_VClone(r.y);
  r.qs     = N_VClone(r.y);
  r.matrix = SUNDenseMatrix(1, 1, r.context);
  r.linear = r.matrix ? SUNLinSol_Dense(r.y, r.matrix, r.context) : nullptr;
  r.memory = IDACreate(r.context);
  if (!check_ptr(r.yp, "N_VClone(yp)") || !check_ptr(r.s, "N_VClone(s)") ||
      !check_ptr(r.sp, "N_VClone(sp)") || !check_ptr(r.q, "N_VClone(q)") ||
      !check_ptr(r.qs, "N_VClone(qs)") || !check_ptr(r.matrix, "SUNDenseMatrix") ||
      !check_ptr(r.linear, "SUNLinSol_Dense") ||
      !check_ptr(r.memory, "IDACreate"))
  {
    return false;
  }
  // Set consistent initial conditions for state and sensitivities.
  probe.memory      = r.memory;
  NV_Ith_S(r.y, 0)  = SUN_RCONST(1.0);
  NV_Ith_S(r.yp, 0) = SUN_RCONST(-1.0);
  NV_Ith_S(r.s, 0)  = SUN_RCONST(0.0);
  NV_Ith_S(r.sp, 0) = SUN_RCONST(-1.0);
  NV_Ith_S(r.q, 0) = NV_Ith_S(r.qs, 0) = SUN_RCONST(0.0);
  const sunrealtype eps = std::numeric_limits<sunrealtype>::epsilon();
  const sunrealtype rtol = std::max(SUN_RCONST(1.0e-7), SUN_RCONST(100.0) * eps);
  const sunrealtype atol = rtol * SUN_RCONST(0.01);
  auto check             = check_flag;
  if (!check(IDAInit(r.memory, residual, SUN_RCONST(0.0), r.y, r.yp),
             "IDAInit") ||
      !check(IDASetUserData(r.memory, &probe), "IDASetUserData") ||
      !check(IDASStolerances(r.memory, rtol, atol), "IDASStolerances") ||
      !check(IDASetLinearSolver(r.memory, r.linear, r.matrix),
             "IDASetLinearSolver") ||
      !check(IDASetJacFn(r.memory, jacobian), "IDASetJacFn") ||
      !check(IDASetInitStep(r.memory, SUN_RCONST(0.001)), "IDASetInitStep") ||
      !check(IDASetMaxStep(r.memory, SUN_RCONST(0.05)), "IDASetMaxStep") ||
      !check(IDASetMaxConvFails(r.memory, 4), "IDASetMaxConvFails"))
  {
    return false;
  }
  // Observe only the staggered state solve to target the extra residual call.
  if (test.method == IDA_STAGGERED)
  {
    r.nonlinear = SUNNonlinSol_Newton(r.y, r.context);
    if (!check_ptr(r.nonlinear, "SUNNonlinSol_Newton")) { return false; }
    probe.solve             = r.nonlinear->ops->solve;
    r.nonlinear->ops->solve = observe_solve;
    if (!check(IDASetNonlinearSolver(r.memory, r.nonlinear),
               "IDASetNonlinearSolver"))
    {
      return false;
    }
  }
  // Configure sensitivity integration and optional auxiliary error control.
  if (!check(IDASensInit(r.memory, 1, test.method, sensitivity, &r.s, &r.sp),
             "IDASensInit") ||
      !check(IDASensEEtolerances(r.memory), "IDASensEEtolerances") ||
      !check(IDASetSensErrCon(r.memory, test.auxiliary ? SUNTRUE : SUNFALSE),
             "IDASetSensErrCon"))
  {
    return false;
  }
  if (test.auxiliary)
  {
    if (!check(IDAQuadInit(r.memory, quadrature, r.q), "IDAQuadInit") ||
        !check(IDAQuadSStolerances(r.memory, rtol, atol),
               "IDAQuadSStolerances") ||
        !check(IDASetQuadErrCon(r.memory, SUNTRUE), "IDASetQuadErrCon") ||
        !check(IDAQuadSensInit(r.memory, quad_sensitivity, &r.qs),
               "IDAQuadSensInit") ||
        !check(IDAQuadSensEEtolerances(r.memory), "IDAQuadSensEEtolerances") ||
        !check(IDASetQuadSensErrCon(r.memory, SUNTRUE), "IDASetQuadSensErrCon"))
    {
      return false;
    }
  }
  sunrealtype previous = SUN_RCONST(0.0), time = SUN_RCONST(0.0);
  bool passed = true, reduced = false;
  int flag  = IDA_SUCCESS;
  int steps = 0;
  // Advance one accepted step at a time and detect gaps after rejected attempts.
  while (time < SUN_RCONST(1.0) && steps++ < 10000)
  {
    flag = IDASolve(r.memory, SUN_RCONST(1.0), &time, r.y, r.yp, IDA_ONE_STEP);
    if (flag != IDA_SUCCESS) { break; }
    sunrealtype h = 0;
    if (!check(IDAGetLastStep(r.memory, &h), "IDAGetLastStep"))
    {
      return false;
    }
    const sunrealtype continuity_tol = SUN_RCONST(128.0) * eps *
                                       std::max(SUN_RCONST(1.0), std::abs(time));
    if (std::abs((time - h) - previous) > continuity_tol)
    {
      std::fprintf(stderr, "step gap: previous=%.17Lg time=%.17Lg h=%.17Lg\n",
                   static_cast<long double>(previous),
                   static_cast<long double>(time), static_cast<long double>(h));
      passed = false;
    }
    if (probe.rejections && previous < probe.rejected_time &&
        time < probe.rejected_time)
    {
      reduced = true;
    }
    previous = time;
  }
  if (steps >= 10000 || probe.watchdog) { passed = false; }
  // Verify failure propagation, retry limits, and step reduction.
  long int failures = -1;
  if (!check(IDAGetNumStepSolveFails(r.memory, &failures),
             "IDAGetNumStepSolveFails"))
  {
    return false;
  }
  const int expected = test.failure == Failure::Persistent ? IDA_REP_RES_ERR
                       : test.failure == Failure::Fatal    ? IDA_RES_FAIL
                                                           : IDA_SUCCESS;
  const int expected_rejections    = test.failure == Failure::Persistent ? 4
                                     : test.failure == Failure::None     ? 0
                                                                         : 1;
  const long int expected_failures = test.failure == Failure::Persistent ? 4
                                     : test.failure == Failure::Once     ? 1
                                                                         : 0;
  passed                           = passed && flag == expected &&
           probe.rejections == expected_rejections &&
           failures == expected_failures;
  if (test.failure == Failure::Once) { passed = passed && reduced; }
  if (test.failure != Failure::None && test.inject_at == 1)
  {
    passed = passed && probe.injected_nst == 0;
  }
  // Check rollback to the last accepted time, including terminal failures.
  sunrealtype internal_time = 0;
  if (!check(IDAGetCurrentTime(r.memory, &internal_time), "IDAGetCurrentTime"))
  {
    return false;
  }
  passed = passed &&
           std::abs(internal_time - previous) <= SUN_RCONST(128.0) * eps;
  if (!check(IDAGetSens(r.memory, &time, &r.s), "IDAGetSens")) { return false; }
  long int accepted_steps = 0;
  if (!check(IDAGetNumSteps(r.memory, &accepted_steps), "IDAGetNumSteps"))
  {
    return false;
  }
  // First-step failure leaves interpolation at order zero.
  if (accepted_steps > 0 &&
      !check(IDAGetSensDky(r.memory, internal_time, 1, &r.sp), "IDAGetSensDky"))
  {
    return false;
  }
  // Compare state and sensitivity outputs against the analytic solution.
  const sunrealtype exact     = std::exp(-internal_time);
  const sunrealtype tolerance = SUN_RCONST(500.0) * rtol;
  passed = passed && std::abs(NV_Ith_S(r.y, 0) - exact) <= tolerance &&
           std::abs(NV_Ith_S(r.yp, 0) + exact) <= tolerance &&
           std::abs(NV_Ith_S(r.s, 0) + internal_time * exact) <= tolerance;
  if (accepted_steps > 0)
  {
    passed = passed &&
             std::abs(NV_Ith_S(r.sp, 0) -
                      (internal_time - SUN_RCONST(1.0)) * exact) <= tolerance;
  }
  if (test.auxiliary)
  {
    if (!check(IDAGetQuad(r.memory, &time, r.q), "IDAGetQuad") ||
        !check(IDAGetQuadSens(r.memory, &time, &r.qs), "IDAGetQuadSens"))
    {
      return false;
    }
    passed =
      passed &&
      std::abs(NV_Ith_S(r.q, 0) - (SUN_RCONST(1.0) - exact)) <= tolerance &&
      std::abs(NV_Ith_S(r.qs, 0) - ((internal_time + SUN_RCONST(1.0)) * exact -
                                    SUN_RCONST(1.0))) <= tolerance;
  }
  if (!passed)
  {
    std::fprintf(stderr, "flag=%d expected=%d rejected=%d failures=%ld reduced=%d watchdog=%d\n",
                 flag, expected, probe.rejections, failures, reduced,
                 probe.watchdog);
  }
  return passed;
}
} // namespace

int main()
{
  // Cover both integration methods and recoverable, persistent, and fatal errors.
  const Case cases[] =
    {{"staggered control", IDA_STAGGERED, Failure::None, 8, false},
     {"simultaneous control", IDA_SIMULTANEOUS, Failure::None, 8, false},
     {"recoverable residual", IDA_STAGGERED, Failure::Once, 8, false},
     {"first-step residual", IDA_STAGGERED, Failure::Once, 1, false},
     {"recoverable with quadratures", IDA_STAGGERED, Failure::Once, 8, true},
     {"persistent residual", IDA_STAGGERED, Failure::Persistent, 8, true},
     {"first-step persistent", IDA_STAGGERED, Failure::Persistent, 1, true},
     {"unrecoverable residual", IDA_STAGGERED, Failure::Fatal, 8, true}};
  int failures = 0;
  for (const auto& test : cases)
  {
    const bool passed = run_case(test);
    std::printf("%s: %s\n", test.name, passed ? "PASS" : "FAIL");
    if (!passed) { ++failures; }
  }
  return failures ? 1 : 0;
}
