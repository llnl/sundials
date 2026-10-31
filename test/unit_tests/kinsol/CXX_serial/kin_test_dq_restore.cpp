/* -----------------------------------------------------------------------------
 * Programmer(s): Demetrios Agourakis
 * Created: 2026-10-31
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2026, the SUNDIALS contributors.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Check that a failed dense DQ Jacobian probe restores the solution vector.
 * The linear-solver setup failure must remain a failure.
 * ---------------------------------------------------------------------------*/

#include <array>
#include <cmath>
#include <cstdio>
#include <limits>

#include "kinsol/kinsol.h"
#include "kinsol/kinsol_ls.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

namespace {
constexpr sunindextype length = 3;

struct ProbeData
{
  std::array<sunrealtype, length> initial = {SUN_RCONST(0.995), SUN_RCONST(0.5),
                                             SUN_RCONST(-0.25)};
  std::array<sunrealtype, length> root    = {SUN_RCONST(0.75), SUN_RCONST(0.25),
                                             SUN_RCONST(-0.5)};
  sunindextype column                     = 0;
  int failure                             = 0;
  int rejections                          = 0;
};

int residual(N_Vector u, N_Vector f, void* user_data)
{
  auto* data   = static_cast<ProbeData*>(user_data);
  auto* values = N_VGetArrayPointer(u);
  if (data->failure != 0 && values[data->column] != data->initial[data->column])
  {
    ++data->rejections;
    return data->failure;
  }
  auto* output = N_VGetArrayPointer(f);
  for (sunindextype i = 0; i < length; ++i)
  {
    output[i] = values[i] - data->root[i];
  }
  return 0;
}

int jacobian(N_Vector, N_Vector, SUNMatrix matrix, void*, N_Vector, N_Vector)
{
  for (sunindextype j = 0; j < length; ++j)
  {
    for (sunindextype i = 0; i < length; ++i)
    {
      SM_ELEMENT_D(matrix, i, j) = i == j ? SUN_RCONST(1.0) : SUN_RCONST(0.0);
    }
  }
  return 0;
}

bool run_case(sunindextype column, int failure, bool exact)
{
  SUNContext context = nullptr;
  if (SUNContext_Create(SUN_COMM_NULL, &context)) { return false; }
  /* Expected setup failures are checked below; omit path-dependent diagnostics. */
  if (SUNContext_ClearErrHandlers(context))
  {
    SUNContext_Free(&context);
    return false;
  }
  N_Vector u             = N_VNew_Serial(length, context);
  N_Vector scale         = N_VNew_Serial(length, context);
  SUNMatrix matrix       = SUNDenseMatrix(length, length, context);
  SUNLinearSolver solver = u && matrix ? SUNLinSol_Dense(u, matrix, context)
                                       : nullptr;
  void* memory           = KINCreate(context);
  bool passed            = u && scale && matrix && solver && memory;
  ProbeData data;
  data.column  = column;
  data.failure = failure;
  if (passed)
  {
    auto* values = N_VGetArrayPointer(u);
    for (sunindextype i = 0; i < length; ++i) { values[i] = data.initial[i]; }
    N_VConst(SUN_RCONST(1.0), scale);
    int flag = KINInit(memory, residual, u);
    if (flag == KIN_SUCCESS) { flag = KINSetUserData(memory, &data); }
    if (flag == KIN_SUCCESS)
    {
      flag = KINSetRelErrFunc(memory, SUN_RCONST(1.0e-4));
    }
    if (flag == KIN_SUCCESS)
    {
      flag = KINSetLinearSolver(memory, solver, matrix);
    }
    if (flag == KIN_SUCCESS && exact) { flag = KINSetJacFn(memory, jacobian); }
    if (flag == KIN_SUCCESS)
    {
      flag = KINSol(memory, u, KIN_NONE, scale, scale);
    }
    long int evaluations = -1;
    passed = KINGetNumLinFuncEvals(memory, &evaluations) == KINLS_SUCCESS;
    if (failure != 0)
    {
      passed = passed && flag == KIN_LSETUP_FAIL && data.rejections == 1 &&
               evaluations == column + 1;
      for (sunindextype i = 0; i < length; ++i)
      {
        if (values[i] != data.initial[i])
        {
          std::fprintf(stderr, "Unrestored solution component %lld\n",
                       static_cast<long long>(i));
          passed = false;
        }
      }
    }
    else
    {
      const sunrealtype tolerance =
        SUN_RCONST(100.0) *
        std::sqrt(std::numeric_limits<sunrealtype>::epsilon());
      passed = passed && flag == KIN_SUCCESS && data.rejections == 0 &&
               (exact ? evaluations == 0 : evaluations >= length);
      for (sunindextype i = 0; i < length; ++i)
      {
        passed = passed && std::abs(values[i] - data.root[i]) <= tolerance;
      }
    }
    if (!passed)
    {
      std::fprintf(stderr, "flag=%d, rejected=%d, DQ evaluations=%ld\n", flag,
                   data.rejections, evaluations);
    }
  }
  if (memory) { KINFree(&memory); }
  if (solver) { SUNLinSolFree(solver); }
  if (matrix) { SUNMatDestroy(matrix); }
  if (scale) { N_VDestroy(scale); }
  if (u) { N_VDestroy(u); }
  SUNContext_Free(&context);
  return passed;
}
} // namespace

int main()
{
  int failures = 0;
  for (int failure : {1, -1})
  {
    for (sunindextype column = 0; column < length; ++column)
    {
      const bool passed = run_case(column, failure, false);
      std::printf("DQ rejection column %lld, callback %d: %s\n",
                  static_cast<long long>(column), failure,
                  passed ? "PASS" : "FAIL");
      if (!passed) { ++failures; }
    }
  }
  for (bool exact : {false, true})
  {
    const bool passed = run_case(0, 0, exact);
    std::printf("%s Jacobian success control: %s\n", exact ? "Explicit" : "DQ",
                passed ? "PASS" : "FAIL");
    if (!passed) { ++failures; }
  }
  return failures == 0 ? 0 : 1;
}
