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

int sensitivity(int, sunrealtype, N_Vector y, N_Vector, N_Vector, N_Vector* s,
                N_Vector* sp, N_Vector* r, void*, N_Vector, N_Vector, N_Vector)
{
  NV_Ith_S(r[0], 0) = NV_Ith_S(sp[0], 0) + NV_Ith_S(s[0], 0) + NV_Ith_S(y, 0);
  return 0;
}

int jacobian(sunrealtype, sunrealtype cj, N_Vector, N_Vector, N_Vector,
             SUNMatrix matrix, void*, N_Vector, N_Vector, N_Vector)
{
  SM_ELEMENT_D(matrix, 0, 0) = cj + SUN_RCONST(1.0);
  return 0;
}

int quadrature(sunrealtype, N_Vector y, N_Vector, N_Vector qdot, void*)
{
  NV_Ith_S(qdot, 0) = NV_Ith_S(y, 0);
  return 0;
}

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
  if (SUNContext_Create(SUN_COMM_NULL, &r.context) != SUN_SUCCESS ||
      SUNContext_ClearErrHandlers(r.context) != SUN_SUCCESS)
  {
    return false;
  }
  r.y = N_VNew_Serial(1, r.context);
  if (!r.y) { return false; }
  r.yp     = N_VClone(r.y);
  r.s      = N_VClone(r.y);
  r.sp     = N_VClone(r.y);
  r.q      = N_VClone(r.y);
  r.qs     = N_VClone(r.y);
  r.matrix = SUNDenseMatrix(1, 1, r.context);
  r.linear = r.matrix ? SUNLinSol_Dense(r.y, r.matrix, r.context) : nullptr;
  r.memory = IDACreate(r.context);
  if (!r.yp || !r.s || !r.sp || !r.q || !r.qs || !r.matrix || !r.linear ||
      !r.memory)
  {
    return false;
  }
  probe.memory      = r.memory;
  NV_Ith_S(r.y, 0)  = SUN_RCONST(1.0);
  NV_Ith_S(r.yp, 0) = SUN_RCONST(-1.0);
  NV_Ith_S(r.s, 0)  = SUN_RCONST(0.0);
  NV_Ith_S(r.sp, 0) = SUN_RCONST(-1.0);
  NV_Ith_S(r.q, 0) = NV_Ith_S(r.qs, 0) = SUN_RCONST(0.0);
  const sunrealtype eps = std::numeric_limits<sunrealtype>::epsilon();
  const sunrealtype rtol = std::max(SUN_RCONST(1.0e-7), SUN_RCONST(100.0) * eps);
  const sunrealtype atol = rtol * SUN_RCONST(0.01);
  auto check             = [&](int flag)
  {
    if (flag != IDA_SUCCESS) { std::fprintf(stderr, "setup flag=%d\n", flag); }
    return flag == IDA_SUCCESS;
  };
  if (!check(IDAInit(r.memory, residual, SUN_RCONST(0.0), r.y, r.yp)) ||
      !check(IDASetUserData(r.memory, &probe)) ||
      !check(IDASStolerances(r.memory, rtol, atol)) ||
      !check(IDASetLinearSolver(r.memory, r.linear, r.matrix)) ||
      !check(IDASetJacFn(r.memory, jacobian)) ||
      !check(IDASetInitStep(r.memory, SUN_RCONST(0.001))) ||
      !check(IDASetMaxStep(r.memory, SUN_RCONST(0.05))) ||
      !check(IDASetMaxConvFails(r.memory, 4)))
  {
    return false;
  }
  if (test.method == IDA_STAGGERED)
  {
    r.nonlinear = SUNNonlinSol_Newton(r.y, r.context);
    if (!r.nonlinear) { return false; }
    probe.solve             = r.nonlinear->ops->solve;
    r.nonlinear->ops->solve = observe_solve;
    if (!check(IDASetNonlinearSolver(r.memory, r.nonlinear))) { return false; }
  }
  if (!check(IDASensInit(r.memory, 1, test.method, sensitivity, &r.s, &r.sp)) ||
      !check(IDASensEEtolerances(r.memory)) ||
      !check(IDASetSensErrCon(r.memory, test.auxiliary ? SUNTRUE : SUNFALSE)))
  {
    return false;
  }
  if (test.auxiliary)
  {
    if (!check(IDAQuadInit(r.memory, quadrature, r.q)) ||
        !check(IDAQuadSStolerances(r.memory, rtol, atol)) ||
        !check(IDASetQuadErrCon(r.memory, SUNTRUE)) ||
        !check(IDAQuadSensInit(r.memory, quad_sensitivity, &r.qs)) ||
        !check(IDAQuadSensEEtolerances(r.memory)) ||
        !check(IDASetQuadSensErrCon(r.memory, SUNTRUE)))
    {
      return false;
    }
  }
  sunrealtype previous = SUN_RCONST(0.0), time = SUN_RCONST(0.0);
  bool passed = true, reduced = false;
  int flag  = IDA_SUCCESS;
  int steps = 0;
  while (time < SUN_RCONST(1.0) && steps++ < 10000)
  {
    flag = IDASolve(r.memory, SUN_RCONST(1.0), &time, r.y, r.yp, IDA_ONE_STEP);
    if (flag != IDA_SUCCESS) { break; }
    sunrealtype h = 0;
    if (!check(IDAGetLastStep(r.memory, &h))) { return false; }
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
  long int failures = -1;
  if (!check(IDAGetNumStepSolveFails(r.memory, &failures))) { return false; }
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
  sunrealtype internal_time = 0;
  if (!check(IDAGetCurrentTime(r.memory, &internal_time))) { return false; }
  passed = passed &&
           std::abs(internal_time - previous) <= SUN_RCONST(128.0) * eps;
  if (!check(IDAGetSens(r.memory, &time, &r.s))) { return false; }
  long int accepted_steps = 0;
  if (!check(IDAGetNumSteps(r.memory, &accepted_steps))) { return false; }
  // First-step failure leaves interpolation at order zero.
  if (accepted_steps > 0 &&
      !check(IDAGetSensDky(r.memory, internal_time, 1, &r.sp)))
  {
    return false;
  }
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
    if (!check(IDAGetQuad(r.memory, &time, r.q)) ||
        !check(IDAGetQuadSens(r.memory, &time, &r.qs)))
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
