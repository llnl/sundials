/* -----------------------------------------------------------------------------
 * Programmer(s): Demetrios Agourakis
 * Origin date: 2026-10-26
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2026, SUNDIALS contributors.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * At zero parameters, y' = 1 + sum_i p_i a_i d_i t^(d_i-1) has y = t,
 * while dy/dp_i = a_i t^d_i. State errors alone cannot control the
 * higher-degree sensitivities. A sine sensitivity also checks a
 * non-polynomial case. No quadratures are used in this test.
 * ---------------------------------------------------------------------------*/

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

#include "idas/idas.h"
#include "nvector/nvector_serial.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

struct Problem
{
  std::vector<int> degrees;
  std::vector<sunrealtype> factors;
};

static constexpr int Sine = -1;
#if defined(SUNDIALS_SINGLE_PRECISION)
static constexpr sunrealtype RelativeTolerance = SUN_RCONST(1e-4);
static constexpr sunrealtype AbsoluteTolerance = SUN_RCONST(1e-6);
static constexpr sunrealtype ErrorBound        = SUN_RCONST(1e-2);
#else
static constexpr sunrealtype RelativeTolerance = SUN_RCONST(1e-6);
static constexpr sunrealtype AbsoluteTolerance = SUN_RCONST(1e-10);
static constexpr sunrealtype ErrorBound        = SUN_RCONST(1e-4);
#endif

static sunrealtype SensitivityValue(sunrealtype t, int n)
{
  if (n == Sine) { return std::sin(t); }
  sunrealtype value = SUN_RCONST(1.0);
  for (int i = 0; i < n; ++i) { value *= t; }
  return value;
}

static int Residual(sunrealtype, N_Vector, N_Vector yp, N_Vector r, void*)
{
  NV_Ith_S(r, 0) = NV_Ith_S(yp, 0) - SUN_RCONST(1.0);
  return 0;
}

static int Sensitivity(int ns, sunrealtype t, N_Vector, N_Vector, N_Vector,
                       N_Vector*, N_Vector* sp, N_Vector* r, void* user_data,
                       N_Vector, N_Vector, N_Vector)
{
  const auto& problem = *static_cast<Problem*>(user_data);
  for (int i = 0; i < ns; ++i)
  {
    const int d                  = problem.degrees[i];
    const sunrealtype derivative = d == Sine ? std::cos(t)
                                             : d * SensitivityValue(t, d - 1);
    NV_Ith_S(r[i], 0) = NV_Ith_S(sp[i], 0) - problem.factors[i] * derivative;
  }
  return 0;
}

static void Check(int flag)
{
  if (flag < 0) { throw std::runtime_error("SUNDIALS call failed"); }
}

struct Resources
{
  SUNContext context = nullptr;
  N_Vector y = nullptr, yp = nullptr;
  std::vector<N_Vector> s, sp;
  void* mem              = nullptr;
  SUNMatrix matrix       = nullptr;
  SUNLinearSolver solver = nullptr;

  ~Resources()
  {
    if (mem) { IDAFree(&mem); }
    if (solver) { SUNLinSolFree(solver); }
    if (matrix) { SUNMatDestroy(matrix); }
    for (auto vector : sp)
    {
      if (vector) { N_VDestroy(vector); }
    }
    for (auto vector : s)
    {
      if (vector) { N_VDestroy(vector); }
    }
    if (yp) { N_VDestroy(yp); }
    if (y) { N_VDestroy(y); }
    if (context) { SUNContext_Free(&context); }
  }
};

static void Run(Problem problem)
{
  Resources r;
  Check(SUNContext_Create(SUN_COMM_NULL, &r.context));
  r.y = N_VNew_Serial(1, r.context);
  if (!r.y) { throw std::runtime_error("vector allocation failed"); }
  r.yp = N_VClone(r.y);
  if (!r.yp) { throw std::runtime_error("vector allocation failed"); }
  N_VConst(SUN_RCONST(0.0), r.y);
  N_VConst(SUN_RCONST(1.0), r.yp);
  const int ns = static_cast<int>(problem.degrees.size());
  r.s.resize(ns, nullptr);
  r.sp.resize(ns, nullptr);
  for (int i = 0; i < ns; ++i)
  {
    r.s[i]  = N_VClone(r.y);
    r.sp[i] = N_VClone(r.y);
    if (!r.s[i] || !r.sp[i])
    {
      throw std::runtime_error("vector allocation failed");
    }
    N_VConst(SUN_RCONST(0.0), r.s[i]);
    N_VConst(problem.degrees[i] == 1 || problem.degrees[i] == Sine
               ? problem.factors[i]
               : SUN_RCONST(0.0),
             r.sp[i]);
  }
  r.mem    = IDACreate(r.context);
  r.matrix = SUNDenseMatrix(1, 1, r.context);
  r.solver = SUNLinSol_Dense(r.y, r.matrix, r.context);
  if (!r.mem || !r.matrix || !r.solver)
  {
    throw std::runtime_error("allocation failed");
  }
  Check(IDAInit(r.mem, Residual, SUN_RCONST(0.0), r.y, r.yp));
  Check(IDASetUserData(r.mem, &problem));
  Check(IDASStolerances(r.mem, RelativeTolerance, AbsoluteTolerance));
  Check(IDASetLinearSolver(r.mem, r.solver, r.matrix));
  Check(IDASensInit(r.mem, ns, IDA_SIMULTANEOUS, Sensitivity, r.s.data(),
                    r.sp.data()));
  std::vector<sunrealtype> atol(ns, AbsoluteTolerance);
  Check(IDASensSStolerances(r.mem, RelativeTolerance, atol.data()));
  Check(IDASetSensErrCon(r.mem, SUNTRUE));
  Check(IDASetStopTime(r.mem, SUN_RCONST(1.0)));
  sunrealtype t = SUN_RCONST(0.0);
  Check(IDASolve(r.mem, SUN_RCONST(1.0), &t, r.y, r.yp, IDA_NORMAL));
  Check(IDAGetSens(r.mem, &t, r.s.data()));
  const sunrealtype state_error = std::abs(NV_Ith_S(r.y, 0) - t);
  if (!std::isfinite(state_error) || state_error > ErrorBound)
  {
    throw std::runtime_error("incorrect state");
  }
  for (int i = 0; i < ns; ++i)
  {
    const sunrealtype expected = problem.factors[i] *
                                 SensitivityValue(t, problem.degrees[i]);
    const sunrealtype error = std::abs(NV_Ith_S(r.s[i], 0) - expected);
    // A loose global-error bound, not an assertion that local rtol bounds global error.
    if (!std::isfinite(error) || error > ErrorBound)
    {
      std::fprintf(stderr, "sensitivity %d (degree %d) failed\n", i,
                   problem.degrees[i]);
      throw std::runtime_error("incorrect sensitivity");
    }
  }
}

int main()
{
  const std::vector<Problem> problems =
    {{{3}, {SUN_RCONST(1.0)}},
     {{5}, {SUN_RCONST(1.0)}},
     {{1, 3}, {SUN_RCONST(1.0), SUN_RCONST(-2.0)}},
     {{3, 5}, {SUN_RCONST(1.0), SUN_RCONST(-2.0)}},
     {{Sine}, {SUN_RCONST(1.0)}},
     {{1}, {SUN_RCONST(1.0)}},
     {{1, 3, 5, 2},
      {SUN_RCONST(1.0), SUN_RCONST(-2.0), SUN_RCONST(0.5), SUN_RCONST(-0.25)}}};
  int failures = 0;
  for (const auto& problem : problems)
  {
    try
    {
      Run(problem);
    }
    catch (const std::exception& e)
    {
      std::fprintf(stderr, "%s\n", e.what());
      ++failures;
    }
  }
  if (failures) { return 1; }
  std::puts("Simultaneous sensitivity error control passed");
  return 0;
}
