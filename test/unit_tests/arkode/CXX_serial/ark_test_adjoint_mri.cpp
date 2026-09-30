/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2026, Lawrence Livermore National Security,
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
 * -----------------------------------------------------------------*/

#include <algorithm>
#include <arkode/arkode_erkstep.h>
#include <arkode/arkode_mristep.h>
#include <cmath>
#include <cstdio>
#include <nvector/nvector_manyvector.h>
#include <nvector/nvector_serial.h>
#include <sunadjointcheckpointscheme/sunadjointcheckpointscheme_fixed.h>
#include <sundials/sundials_adjointstepper.h>
#include <sundials/sundials_context.h>
#include <sunmemory/sunmemory_system.h>

#include "problems/lotka_volterra.hpp"

using namespace problems::lotka_volterra;

struct UserData
{
  sunrealtype* parameters;
  sunbooleantype has_mu;
};

static int slow_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  auto* data = static_cast<UserData*>(user_data);
  return ode_rhs_1(t, y, ydot, data->parameters);
}

static int fast_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  auto* data = static_cast<UserData*>(user_data);
  return ode_rhs_2(t, y, ydot, data->parameters);
}

static int slow_adj_rhs(sunrealtype t, N_Vector y, N_Vector sens,
                        N_Vector sens_dot, void* user_data)
{
  auto* data = static_cast<UserData*>(user_data);
  if (N_VGetNumSubvectors_ManyVector(sens) != 1 + data->has_mu ||
      N_VGetNumSubvectors_ManyVector(sens_dot) != 1 + data->has_mu)
  {
    return -1;
  }
  N_Vector lambda     = N_VGetSubvector_ManyVector(sens, 0);
  N_Vector lambda_dot = N_VGetSubvector_ManyVector(sens_dot, 0);
  int retval = ode_vjp_1(lambda, lambda_dot, t, y, nullptr, data->parameters,
                         nullptr);
  if (retval != 0 || !data->has_mu) { return retval; }
  return parameter_vjp_1(lambda, N_VGetSubvector_ManyVector(sens_dot, 1), t, y,
                         nullptr, data->parameters, nullptr);
}

static int fast_adj_rhs(sunrealtype t, N_Vector y, N_Vector sens,
                        N_Vector sens_dot, void* user_data)
{
  auto* data = static_cast<UserData*>(user_data);
  if (N_VGetNumSubvectors_ManyVector(sens) != 1 + data->has_mu ||
      N_VGetNumSubvectors_ManyVector(sens_dot) != 1 + data->has_mu)
  {
    return -1;
  }
  N_Vector lambda     = N_VGetSubvector_ManyVector(sens, 0);
  N_Vector lambda_dot = N_VGetSubvector_ManyVector(sens_dot, 0);
  int retval = ode_vjp_2(lambda, lambda_dot, t, y, nullptr, data->parameters,
                         nullptr);
  if (retval != 0 || !data->has_mu) { return retval; }
  return parameter_vjp_2(lambda, N_VGetSubvector_ManyVector(sens_dot, 1), t, y,
                         nullptr, data->parameters, nullptr);
}

static int run_test(sunbooleantype keep, sunbooleantype has_mu)
{
  constexpr sunrealtype t0  = SUN_RCONST(0.0);
  constexpr sunrealtype tf  = SUN_RCONST(0.3);
  constexpr sunrealtype hs  = SUN_RCONST(0.1);
  constexpr sunrealtype hf  = SUN_RCONST(1.0) / SUN_RCONST(300.0);
  sunrealtype parameters[4] = {SUN_RCONST(1.5), SUN_RCONST(1.0),
                               SUN_RCONST(3.0), SUN_RCONST(1.0)};
  UserData data{parameters, has_mu};
  SUNContext sunctx                        = nullptr;
  SUNMemoryHelper helper                   = nullptr;
  N_Vector y                               = nullptr;
  N_Vector lambda                          = nullptr;
  N_Vector mu                              = nullptr;
  N_Vector sf                              = nullptr;
  SUNStepper inner_sunstepper              = nullptr;
  SUNAdjointStepper inner_adjoint          = nullptr;
  SUNAdjointStepper outer_adjoint          = nullptr;
  SUNAdjointCheckpointScheme outer_scheme  = nullptr;
  SUNAdjointCheckpointScheme inner_scheme  = nullptr;
  MRIStepInnerAdjointProblem inner_problem = nullptr;
  void* inner_mem                          = nullptr;
  void* outer_mem                          = nullptr;
  int failures                             = 0;

  SUNContext_Create(SUN_COMM_NULL, &sunctx);
  helper = SUNMemoryHelper_Sys(sunctx);
  y      = N_VNew_Serial(2, sunctx);
  N_VConst(SUN_RCONST(1.0), y);

  inner_mem = ERKStepCreate(fast_rhs, t0, y, sunctx);
  ARKodeSetUserData(inner_mem, &data);
  ARKodeSetFixedStep(inner_mem, hf);
  ARKodeCreateSUNStepper(inner_mem, &inner_sunstepper);

  outer_mem = MRIStepCreate(slow_rhs, nullptr, t0, y, inner_sunstepper, sunctx);
  ARKodeSetUserData(outer_mem, &data);
  ARKodeSetFixedStep(outer_mem, hs);
  MRIStepCoupling coupling = MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_ERK33a);
  MRIStepSetCoupling(outer_mem, coupling);
  MRIStepCoupling_Free(coupling);

  SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM, helper, 2, 4,
                                          keep, sunctx, &outer_scheme);
  ARKodeSetAdjointCheckpointScheme(outer_mem, outer_scheme);
  sunrealtype tret = t0;
  if (ARKodeEvolve(outer_mem, tf, y, &tret, ARK_NORMAL) < 0) { failures++; }

  lambda = N_VNew_Serial(2, sunctx);
  if (has_mu) { mu = N_VNew_Serial(4, sunctx); }
  sunrealtype* y_data      = N_VGetArrayPointer(y);
  sunrealtype* lambda_data = N_VGetArrayPointer(lambda);
  lambda_data[0]           = y_data[0] - SUN_RCONST(1.0);
  lambda_data[1]           = y_data[1] - SUN_RCONST(1.0);
  if (has_mu) { N_VConst(SUN_RCONST(0.0), mu); }
  N_Vector subvectors[2] = {lambda, mu};
  sf = N_VNew_ManyVector(has_mu ? 2 : 1, subvectors, sunctx);

  if (MRIStepInnerAdjointProblem_Create(outer_mem, fast_adj_rhs, sf, &data,
                                        &inner_problem) != ARK_SUCCESS)
  {
    failures++;
  }
  N_Vector inner_sf         = nullptr;
  SUNAdjRhsFn inner_adj_rhs = nullptr;
  MRIStepInnerAdjointProblem_GetTerminalState(inner_problem, &inner_sf);
  MRIStepInnerAdjointProblem_GetAdjRhsFn(inner_problem, &inner_adj_rhs);
  N_Vector inner_parameters = N_VGetSubvector_ManyVector(inner_sf, 1);
  if (!inner_parameters ||
      N_VGetNumSubvectors_ManyVector(inner_parameters) != 2 + has_mu)
  {
    failures++;
  }

  SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM, helper, 2, 64,
                                          keep, sunctx, &inner_scheme);
  ARKodeSetAdjointCheckpointScheme(inner_mem, inner_scheme);
  if (ERKStepCreateAdjointStepper(inner_mem, inner_adj_rhs, tf, inner_sf,
                                  sunctx, &inner_adjoint) != ARK_SUCCESS)
  {
    failures++;
  }

  if (MRIStepCreateAdjointStepper(outer_mem, inner_adjoint, slow_adj_rhs,
                                  nullptr, inner_problem, tf, sf, sunctx,
                                  &outer_adjoint) != ARK_SUCCESS)
  {
    failures++;
  }
  if (!failures &&
      SUNAdjointStepper_Evolve(outer_adjoint, t0, sf, &tret) != SUN_SUCCESS)
  {
    failures++;
  }
  suncountertype adj_nst = 0;
  SUNAdjointStepper_GetNumSteps(outer_adjoint, &adj_nst);

  // Reference values computed with ERKStep using a 1e-5 fixed step.
  constexpr sunrealtype reference_lambda[2] = {SUN_RCONST(2.204924168753578e-01),
                                               SUN_RCONST(-3.050319456561253e-01)};
  constexpr sunrealtype reference_mu[4] = {SUN_RCONST(7.974709697811572e-02),
                                           SUN_RCONST(-5.950658873178583e-02),
                                           SUN_RCONST(8.224828446314870e-02),
                                           SUN_RCONST(-9.076588019985569e-02)};
  sunrealtype* mu_data     = has_mu ? N_VGetArrayPointer(mu) : nullptr;
  sunrealtype lambda_error = SUN_RCONST(0.0);
  sunrealtype mu_error     = SUN_RCONST(0.0);
  for (int i = 0; i < 2; i++)
  {
    lambda_error = std::max(lambda_error,
                            std::abs(lambda_data[i] - reference_lambda[i]));
  }
  for (int i = 0; has_mu && i < 4; i++)
  {
    mu_error = std::max(mu_error, std::abs(mu_data[i] - reference_mu[i]));
  }
  std::printf("mu = %d, keep = %d, steps = %ld, lambda error = %.3e, "
              "parameter error = %.3e\n",
              static_cast<int>(has_mu), static_cast<int>(keep),
              static_cast<long>(adj_nst), static_cast<double>(lambda_error),
              static_cast<double>(mu_error));
  if (lambda_error > SUN_RCONST(2.0e-4) || mu_error > SUN_RCONST(2.0e-4))
  {
    failures++;
  }

  SUNAdjointStepper_Destroy(&outer_adjoint);
  SUNAdjointStepper_Destroy(&inner_adjoint);
  MRIStepInnerAdjointProblem_Free(inner_problem);
  SUNAdjointCheckpointScheme_Destroy(&inner_scheme);
  SUNAdjointCheckpointScheme_Destroy(&outer_scheme);
  ARKodeFree(&outer_mem);
  SUNStepper_Destroy(&inner_sunstepper);
  ARKodeFree(&inner_mem);
  N_VDestroy(sf);
  N_VDestroy(lambda);
  if (mu) { N_VDestroy(mu); }
  N_VDestroy(y);
  SUNMemoryHelper_Destroy(helper);
  SUNContext_Free(&sunctx);
  return failures;
}

int main()
{
  int failures = run_test(SUNTRUE, SUNTRUE) + run_test(SUNFALSE, SUNTRUE) +
                 run_test(SUNTRUE, SUNFALSE) + run_test(SUNFALSE, SUNFALSE);
  std::printf("MRIStep adjoint failures = %d\n", failures);
  return failures;
}
