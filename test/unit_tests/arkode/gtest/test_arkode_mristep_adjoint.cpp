/* -----------------------------------------------------------------
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
 * -----------------------------------------------------------------*/

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>

#include <gtest/gtest.h>

#include <arkode/arkode.h>
#include <arkode/arkode_butcher_erk.h>
#include <arkode/arkode_erkstep.h>
#include <arkode/arkode_mristep.h>
#include <nvector/nvector_manyvector.h>
#include <nvector/nvector_serial.h>
#include <sunadjointcheckpointscheme/sunadjointcheckpointscheme_fixed.h>
#include <sundials/sundials_adjointstepper.h>
#include <sundials/sundials_context.hpp>
#include <sundials/sundials_math.h>
#include <sunmemory/sunmemory_system.h>

#include "arkode/arkode_mristep_impl.h"
#include "problems/lotka_volterra.hpp"

using namespace problems::lotka_volterra;

namespace {

constexpr sunrealtype t0 = SUN_RCONST(0.0);
constexpr sunrealtype tf = SUN_RCONST(0.3);

struct TestConfig
{
  const char* name;
  sunbooleantype keep;
  sunbooleantype has_mu;
  suncountertype outer_frequency;
  suncountertype inner_frequency;
  sunrealtype slow_step;
  sunrealtype fast_step;
};

struct UserData // TODO(SBR): this is inconsistent with problem user data
{
  std::array<sunrealtype, 4> parameters = {SUN_RCONST(1.5), SUN_RCONST(1.0),
                                           SUN_RCONST(3.0), SUN_RCONST(1.0)};
  sunbooleantype has_mu                 = SUNFALSE;
  sunbooleantype zero_fast              = SUNFALSE;
  bool slow_rhs_enabled                 = true;
  long int slow_rhs_calls               = 0;
  long int fast_rhs_calls               = 0;
  long int slow_adj_rhs_calls           = 0;
  long int fast_adj_rhs_calls           = 0;
};

struct Result
{
  std::array<sunrealtype, 2> state{};
  std::array<sunrealtype, 2> lambda{};
  std::array<sunrealtype, 2> interpolated_lambda{};
  std::array<sunrealtype, 4> mu{};
  sunrealtype forward_time     = t0;
  sunrealtype adjoint_time     = tf;
  long int forward_steps       = 0;
  suncountertype adjoint_steps = 0;
  bool has_interpolated_state  = false;
  bool valid                   = false;
  sunrealtype cost             = SUN_RCONST(0.0);
  std::array<long int, 2> forward_rhs_calls{};
  std::array<long int, 4> reverse_rhs_calls{};
};

struct RunOptions
{
  MRIStepCoupling coupling = nullptr;
  int perturb              = -1;
  sunrealtype delta        = SUN_RCONST(0.0);
  bool forward_only        = false;
};

struct Resources
{
  sundials::Context sunctx;
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

  ~Resources()
  {
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
    N_VDestroy(mu);
    N_VDestroy(y);
    SUNMemoryHelper_Destroy(helper);
  }
};

static int slow_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  auto* data = static_cast<UserData*>(user_data);
  data->slow_rhs_calls++;
  if (!data->slow_rhs_enabled) { return -1; }
  return data->zero_fast ? ode_rhs(t, y, ydot, data->parameters.data())
                         : ode_rhs_1(t, y, ydot, data->parameters.data());
}

static int fast_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  auto* data = static_cast<UserData*>(user_data);
  data->fast_rhs_calls++;
  if (data->zero_fast)
  {
    N_VConst(SUN_RCONST(0.0), ydot);
    return 0;
  }
  return ode_rhs_2(t, y, ydot, data->parameters.data());
}

static int full_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  auto* data = static_cast<UserData*>(user_data);
  return ode_rhs(t, y, ydot, data->parameters.data());
}

static int adj_rhs_impl(sunrealtype t, N_Vector y, N_Vector sens,
                        N_Vector sens_dot, UserData* data, bool fast)
{
  if (N_VGetNumSubvectors_ManyVector(sens) != 1 + data->has_mu ||
      N_VGetNumSubvectors_ManyVector(sens_dot) != 1 + data->has_mu)
  {
    return -1;
  }

  N_Vector lambda     = N_VGetSubvector_ManyVector(sens, 0);
  N_Vector lambda_dot = N_VGetSubvector_ManyVector(sens_dot, 0);
  int retval;

  if (fast && data->zero_fast)
  {
    N_VConst(SUN_RCONST(0.0), lambda_dot);
    retval = 0;
  }
  else if (fast)
  {
    retval = ode_vjp_2(lambda, lambda_dot, t, y, nullptr,
                       data->parameters.data(), nullptr);
  }
  else if (data->zero_fast)
  {
    retval = ode_vjp(lambda, lambda_dot, t, y, nullptr, data->parameters.data(),
                     nullptr);
  }
  else
  {
    retval = ode_vjp_1(lambda, lambda_dot, t, y, nullptr,
                       data->parameters.data(), nullptr);
  }
  if (retval != 0 || !data->has_mu) { return retval; }

  N_Vector mu_dot = N_VGetSubvector_ManyVector(sens_dot, 1);
  if (fast && data->zero_fast)
  {
    N_VConst(SUN_RCONST(0.0), mu_dot);
    return 0;
  }
  if (fast)
  {
    return parameter_vjp_2(lambda, mu_dot, t, y, nullptr,
                           data->parameters.data(), nullptr);
  }
  if (data->zero_fast)
  {
    return parameter_vjp(lambda, mu_dot, t, y, nullptr, data->parameters.data(),
                         nullptr);
  }
  return parameter_vjp_1(lambda, mu_dot, t, y, nullptr, data->parameters.data(),
                         nullptr);
}

static int slow_adj_rhs(sunrealtype t, N_Vector y, N_Vector sens,
                        N_Vector sens_dot, void* user_data)
{
  static_cast<UserData*>(user_data)->slow_adj_rhs_calls++;
  return adj_rhs_impl(t, y, sens, sens_dot, static_cast<UserData*>(user_data),
                      false);
}

static int fast_adj_rhs(sunrealtype t, N_Vector y, N_Vector sens,
                        N_Vector sens_dot, void* user_data)
{
  static_cast<UserData*>(user_data)->fast_adj_rhs_calls++;
  return adj_rhs_impl(t, y, sens, sens_dot, static_cast<UserData*>(user_data),
                      true);
}

static int full_adj_rhs(sunrealtype t, N_Vector y, N_Vector sens,
                        N_Vector sens_dot, void* user_data)
{
  auto* data = static_cast<UserData*>(user_data);
  if (N_VGetNumSubvectors_ManyVector(sens) != 1 + data->has_mu ||
      N_VGetNumSubvectors_ManyVector(sens_dot) != 1 + data->has_mu)
  {
    return -1;
  }

  N_Vector lambda = N_VGetSubvector_ManyVector(sens, 0);
  int retval = ode_vjp(lambda, N_VGetSubvector_ManyVector(sens_dot, 0), t, y,
                       nullptr, data->parameters.data(), nullptr);
  if (retval != 0 || !data->has_mu) { return retval; }
  return parameter_vjp(lambda, N_VGetSubvector_ManyVector(sens_dot, 1), t, y,
                       nullptr, data->parameters.data(), nullptr);
}

static void create_terminal_state(Resources& resources, UserData& data)
{
  resources.lambda = N_VNew_Serial(2, resources.sunctx);
  ASSERT_NE(resources.lambda, nullptr);
  if (data.has_mu)
  {
    resources.mu = N_VNew_Serial(4, resources.sunctx);
    ASSERT_NE(resources.mu, nullptr);
    N_VConst(SUN_RCONST(0.0), resources.mu);
  }

  sunrealtype* y_data      = N_VGetArrayPointer(resources.y);
  sunrealtype* lambda_data = N_VGetArrayPointer(resources.lambda);
  lambda_data[0]           = y_data[0] - SUN_RCONST(1.0);
  lambda_data[1]           = y_data[1] - SUN_RCONST(1.0);

  N_Vector subvectors[2] = {resources.lambda, resources.mu};
  resources.sf           = N_VNew_ManyVector(data.has_mu ? 2 : 1, subvectors,
                                   resources.sunctx);
  ASSERT_NE(resources.sf, nullptr);
}

static void save_result(Resources& resources, Result& result,
                        sunrealtype forward_time, sunrealtype adjoint_time)
{
  const sunrealtype* y_data      = N_VGetArrayPointer(resources.y);
  const sunrealtype* lambda_data = N_VGetArrayPointer(resources.lambda);
  std::copy_n(y_data, result.state.size(), result.state.begin());
  std::copy_n(lambda_data, result.lambda.size(), result.lambda.begin());
  if (resources.mu)
  {
    const sunrealtype* mu_data = N_VGetArrayPointer(resources.mu);
    std::copy_n(mu_data, result.mu.size(), result.mu.begin());
  }
  result.forward_time = forward_time;
  result.adjoint_time = adjoint_time;
  result.valid        = true;
}

template<typename Actual, typename Expected>
static sunrealtype max_error(const Actual& actual, const Expected& expected);

static void run_mristep(const TestConfig& config, bool zero_fast, Result& result,
                        bool reinit = false, bool interpolate = false,
                        bool forbid_slow_rhs_during_adjoint = false,
                        const RunOptions& options           = {})
{
  Resources resources;
  UserData data;
  data.has_mu    = config.has_mu;
  data.zero_fast = zero_fast;
  if (options.perturb >= 2)
  {
    data.parameters[options.perturb - 2] += options.delta;
  }

  resources.helper = SUNMemoryHelper_Sys(resources.sunctx);
  ASSERT_NE(resources.helper, nullptr);
  resources.y = N_VNew_Serial(2, resources.sunctx);
  ASSERT_NE(resources.y, nullptr);
  N_VConst(SUN_RCONST(1.0), resources.y);
  if (options.perturb >= 0 && options.perturb < 2)
  {
    N_VGetArrayPointer(resources.y)[options.perturb] += options.delta;
  }

  resources.inner_mem = ERKStepCreate(fast_rhs, t0, resources.y,
                                      resources.sunctx);
  ASSERT_NE(resources.inner_mem, nullptr);
  ASSERT_EQ(ARKodeSetUserData(resources.inner_mem, &data), ARK_SUCCESS);
  ASSERT_EQ(ARKodeSetFixedStep(resources.inner_mem, config.fast_step),
            ARK_SUCCESS);
  ASSERT_EQ(ARKodeCreateSUNStepper(resources.inner_mem,
                                   &resources.inner_sunstepper),
            ARK_SUCCESS);

  resources.outer_mem = MRIStepCreate(slow_rhs, nullptr, t0, resources.y,
                                      resources.inner_sunstepper,
                                      resources.sunctx);
  ASSERT_NE(resources.outer_mem, nullptr);
  ASSERT_EQ(ARKodeSetUserData(resources.outer_mem, &data), ARK_SUCCESS);
  ASSERT_EQ(ARKodeSetFixedStep(resources.outer_mem, config.slow_step),
            ARK_SUCCESS);

  MRIStepCoupling coupling = nullptr;
  if (options.coupling) { coupling = MRIStepCoupling_Copy(options.coupling); }
  else if (zero_fast)
  {
    ARKodeButcherTable table =
      ARKodeButcherTable_LoadERK(ARKODE_FORWARD_EULER_1_1);
    ASSERT_NE(table, nullptr);
    coupling = MRIStepCoupling_MIStoMRI(table, table->q, 0);
    ARKodeButcherTable_Free(table);
  }
  else { coupling = MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_ERK33a); }
  ASSERT_NE(coupling, nullptr);
  ASSERT_EQ(MRIStepSetCoupling(resources.outer_mem, coupling), ARK_SUCCESS);
  MRIStepCoupling_Free(coupling);

  ASSERT_EQ(SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM,
                                                    resources.helper,
                                                    config.outer_frequency, 16,
                                                    config.keep, resources.sunctx,
                                                    &resources.outer_scheme),
            SUN_SUCCESS);
  ASSERT_EQ(ARKodeSetAdjointCheckpointScheme(resources.outer_mem,
                                             resources.outer_scheme),
            ARK_SUCCESS);

  sunrealtype forward_time = t0;
  ASSERT_GE(ARKodeEvolve(resources.outer_mem, tf, resources.y, &forward_time,
                         ARK_NORMAL),
            ARK_SUCCESS);
  ASSERT_EQ(ARKodeGetNumSteps(resources.outer_mem, &result.forward_steps),
            ARK_SUCCESS);
  EXPECT_GT(result.forward_steps, 0);
  const long int forward_slow_rhs_calls = data.slow_rhs_calls;
  result.forward_rhs_calls     = {data.slow_rhs_calls, data.fast_rhs_calls};
  const sunrealtype* forward_y = N_VGetArrayPointer(resources.y);
  for (int i = 0; i < 2; i++)
  {
    sunrealtype difference = forward_y[i] - SUN_RCONST(1.0);
    result.cost += SUN_RCONST(0.5) * difference * difference;
  }
  if (options.forward_only)
  {
    result.valid = true;
    return;
  }
  if (forbid_slow_rhs_during_adjoint) { data.slow_rhs_enabled = false; }

  create_terminal_state(resources, data);
  ASSERT_EQ(MRIStepInnerAdjointProblem_Create(resources.outer_mem, fast_adj_rhs,
                                              resources.sf, &data,
                                              &resources.inner_problem),
            ARK_SUCCESS);

  N_Vector inner_sf         = nullptr;
  SUNAdjRhsFn inner_adj_rhs = nullptr;
  ASSERT_EQ(MRIStepInnerAdjointProblem_GetTerminalState(resources.inner_problem,
                                                        &inner_sf),
            ARK_SUCCESS);
  N_Vector inner_parameters = N_VGetSubvector_ManyVector(inner_sf, 1);
  ASSERT_NE(inner_parameters, nullptr);
  const sunindextype nparameters =
    N_VGetNumSubvectors_ManyVector(inner_parameters);
  ASSERT_GE(nparameters, data.has_mu ? 1 : 0);
  const sunindextype nomega = nparameters - (data.has_mu ? 1 : 0);
  for (sunindextype k = 0; k < nomega; k++)
  {
    EXPECT_EQ(N_VGetLength(N_VGetSubvector_ManyVector(inner_parameters, k)),
              N_VGetLength(resources.lambda));
  }
  if (data.has_mu)
  {
    EXPECT_EQ(N_VGetLength(N_VGetSubvector_ManyVector(inner_parameters, nomega)),
              N_VGetLength(resources.mu));
  }
  ASSERT_EQ(MRIStepInnerAdjointProblem_GetAdjRhsFn(resources.inner_problem,
                                                   &inner_adj_rhs),
            ARK_SUCCESS);

  ASSERT_EQ(SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM,
                                                    resources.helper,
                                                    config.inner_frequency, 128,
                                                    config.keep, resources.sunctx,
                                                    &resources.inner_scheme),
            SUN_SUCCESS);
  ASSERT_EQ(ARKodeSetAdjointCheckpointScheme(resources.inner_mem,
                                             resources.inner_scheme),
            ARK_SUCCESS);
  ASSERT_EQ(ERKStepCreateAdjointStepper(resources.inner_mem, inner_adj_rhs, tf,
                                        inner_sf, resources.sunctx,
                                        &resources.inner_adjoint),
            ARK_SUCCESS);
  ASSERT_EQ(MRIStepCreateAdjointStepper(resources.outer_mem,
                                        resources.inner_adjoint, slow_adj_rhs,
                                        nullptr, resources.inner_problem, tf,
                                        resources.sf, resources.sunctx,
                                        &resources.outer_adjoint),
            ARK_SUCCESS);

  sunrealtype adjoint_time = tf;
  if (interpolate)
  {
    const sunrealtype tout = tf - config.slow_step / SUN_RCONST(2.0);
    ASSERT_EQ(SUNAdjointStepper_Evolve(resources.outer_adjoint, tout,
                                       resources.sf, &adjoint_time),
              SUN_SUCCESS);
    EXPECT_NEAR(adjoint_time, tout, SUN_RCONST(10.0) * SUN_UNIT_ROUNDOFF);
    std::copy_n(N_VGetArrayPointer(resources.lambda),
                result.interpolated_lambda.size(),
                result.interpolated_lambda.begin());
    result.has_interpolated_state = true;
    ASSERT_EQ(SUNAdjointStepper_Evolve(resources.outer_adjoint, t0,
                                       resources.sf, &adjoint_time),
              SUN_SUCCESS);
  }
  else
  {
    for (long int i = 0; i < result.forward_steps; i++)
    {
      ASSERT_EQ(SUNAdjointStepper_OneStep(resources.outer_adjoint, t0,
                                          resources.sf, &adjoint_time),
                SUN_SUCCESS);
    }
  }

  if (reinit)
  {
    std::array<sunrealtype, 2> first_lambda;
    std::array<sunrealtype, 4> first_mu{};
    std::copy_n(N_VGetArrayPointer(resources.lambda), first_lambda.size(),
                first_lambda.begin());
    if (resources.mu)
    {
      std::copy_n(N_VGetArrayPointer(resources.mu), first_mu.size(),
                  first_mu.begin());
      N_VConst(SUN_RCONST(0.0), resources.mu);
    }
    const sunrealtype* y_data = N_VGetArrayPointer(resources.y);
    sunrealtype* lambda_data  = N_VGetArrayPointer(resources.lambda);
    lambda_data[0]            = y_data[0] - SUN_RCONST(1.0);
    lambda_data[1]            = y_data[1] - SUN_RCONST(1.0);

    ASSERT_EQ(SUNAdjointStepper_ReInit(resources.outer_adjoint, tf, resources.sf,
                                       static_cast<suncountertype>(
                                         result.forward_steps) -
                                         1),
              SUN_SUCCESS);
    adjoint_time = tf;
    for (long int i = 0; i < result.forward_steps; i++)
    {
      ASSERT_EQ(SUNAdjointStepper_OneStep(resources.outer_adjoint, t0,
                                          resources.sf, &adjoint_time),
                SUN_SUCCESS);
    }
    EXPECT_LE(max_error(first_lambda, std::array<sunrealtype, 2>{lambda_data[0],
                                                                 lambda_data[1]}),
              SUN_RCONST(10.0) * SUN_UNIT_ROUNDOFF);
    if (resources.mu)
    {
      const sunrealtype* mu_data = N_VGetArrayPointer(resources.mu);
      EXPECT_LE(max_error(first_mu,
                          std::array<sunrealtype, 4>{mu_data[0], mu_data[1],
                                                     mu_data[2], mu_data[3]}),
                SUN_RCONST(10.0) * SUN_UNIT_ROUNDOFF);
    }
  }
  ASSERT_EQ(SUNAdjointStepper_GetNumSteps(resources.outer_adjoint,
                                          &result.adjoint_steps),
            SUN_SUCCESS);
  EXPECT_EQ(result.adjoint_steps,
            static_cast<suncountertype>(result.forward_steps));
  if (forbid_slow_rhs_during_adjoint)
  {
    EXPECT_EQ(data.slow_rhs_calls, forward_slow_rhs_calls);
  }
  result.reverse_rhs_calls = {data.slow_rhs_calls - result.forward_rhs_calls[0],
                              data.fast_rhs_calls - result.forward_rhs_calls[1],
                              data.slow_adj_rhs_calls, data.fast_adj_rhs_calls};
  save_result(resources, result, forward_time, adjoint_time);
}

static void run_erkstep(const TestConfig& config, Result& result,
                        bool interpolate = false)
{
  Resources resources;
  UserData data;
  data.has_mu    = config.has_mu;
  data.zero_fast = SUNTRUE;

  resources.helper = SUNMemoryHelper_Sys(resources.sunctx);
  ASSERT_NE(resources.helper, nullptr);
  resources.y = N_VNew_Serial(2, resources.sunctx);
  ASSERT_NE(resources.y, nullptr);
  N_VConst(SUN_RCONST(1.0), resources.y);

  resources.outer_mem = ERKStepCreate(full_rhs, t0, resources.y,
                                      resources.sunctx);
  ASSERT_NE(resources.outer_mem, nullptr);
  ASSERT_EQ(ARKodeSetUserData(resources.outer_mem, &data), ARK_SUCCESS);
  ASSERT_EQ(ARKodeSetFixedStep(resources.outer_mem, config.slow_step),
            ARK_SUCCESS);
  ARKodeButcherTable table = ARKodeButcherTable_LoadERK(ARKODE_FORWARD_EULER_1_1);
  ASSERT_NE(table, nullptr);
  ASSERT_EQ(ERKStepSetTable(resources.outer_mem, table), ARK_SUCCESS);
  ARKodeButcherTable_Free(table);

  ASSERT_EQ(SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM,
                                                    resources.helper,
                                                    config.outer_frequency, 16,
                                                    SUNTRUE, resources.sunctx,
                                                    &resources.outer_scheme),
            SUN_SUCCESS);
  ASSERT_EQ(ARKodeSetAdjointCheckpointScheme(resources.outer_mem,
                                             resources.outer_scheme),
            ARK_SUCCESS);

  sunrealtype forward_time = t0;
  ASSERT_GE(ARKodeEvolve(resources.outer_mem, tf, resources.y, &forward_time,
                         ARK_NORMAL),
            ARK_SUCCESS);
  ASSERT_EQ(ARKodeGetNumSteps(resources.outer_mem, &result.forward_steps),
            ARK_SUCCESS);
  EXPECT_GT(result.forward_steps, 0);

  create_terminal_state(resources, data);
  ASSERT_EQ(ERKStepCreateAdjointStepper(resources.outer_mem, full_adj_rhs, tf,
                                        resources.sf, resources.sunctx,
                                        &resources.outer_adjoint),
            ARK_SUCCESS);

  sunrealtype adjoint_time = tf;
  if (interpolate)
  {
    const sunrealtype tout = tf - config.slow_step / SUN_RCONST(2.0);
    ASSERT_EQ(SUNAdjointStepper_Evolve(resources.outer_adjoint, tout,
                                       resources.sf, &adjoint_time),
              SUN_SUCCESS);
    EXPECT_NEAR(adjoint_time, tout, SUN_RCONST(10.0) * SUN_UNIT_ROUNDOFF);
    std::copy_n(N_VGetArrayPointer(resources.lambda),
                result.interpolated_lambda.size(),
                result.interpolated_lambda.begin());
    result.has_interpolated_state = true;
    ASSERT_EQ(SUNAdjointStepper_Evolve(resources.outer_adjoint, t0,
                                       resources.sf, &adjoint_time),
              SUN_SUCCESS);
  }
  else
  {
    for (long int i = 0; i < result.forward_steps; i++)
    {
      ASSERT_EQ(SUNAdjointStepper_OneStep(resources.outer_adjoint, t0,
                                          resources.sf, &adjoint_time),
                SUN_SUCCESS);
    }
  }
  ASSERT_EQ(SUNAdjointStepper_GetNumSteps(resources.outer_adjoint,
                                          &result.adjoint_steps),
            SUN_SUCCESS);
  EXPECT_EQ(result.adjoint_steps,
            static_cast<suncountertype>(result.forward_steps));
  save_result(resources, result, forward_time, adjoint_time);
}

template<typename Actual, typename Expected>
static sunrealtype max_error(const Actual& actual, const Expected& expected)
{
  sunrealtype error = SUN_RCONST(0.0);
  for (std::size_t i = 0; i < actual.size(); i++)
  {
    error = std::max(error, SUNRabs(actual[i] - expected[i]));
  }
  return error;
}

class MRIStepAdjointTest : public testing::TestWithParam<TestConfig>
{};

TEST(MRIStepAdjointValidationTest, RejectsNullOutputPointer)
{
  EXPECT_EQ(MRIStepCreateAdjointStepper(nullptr, nullptr, nullptr, nullptr,
                                        nullptr, tf, nullptr, nullptr, nullptr),
            ARK_ILL_INPUT);
}

TEST(MRIStepAdjointValidationTest, RejectsAnEmptyForwardHistory)
{
  Resources resources;
  UserData data;
  resources.helper = SUNMemoryHelper_Sys(resources.sunctx);
  ASSERT_NE(resources.helper, nullptr);
  resources.y = N_VNew_Serial(2, resources.sunctx);
  ASSERT_NE(resources.y, nullptr);
  N_VConst(SUN_RCONST(1.0), resources.y);

  resources.inner_mem = ERKStepCreate(fast_rhs, t0, resources.y,
                                      resources.sunctx);
  ASSERT_NE(resources.inner_mem, nullptr);
  ASSERT_EQ(ARKodeSetFixedStep(resources.inner_mem, SUN_RCONST(0.01)),
            ARK_SUCCESS);
  ASSERT_EQ(ARKodeCreateSUNStepper(resources.inner_mem,
                                   &resources.inner_sunstepper),
            ARK_SUCCESS);
  resources.outer_mem = MRIStepCreate(slow_rhs, nullptr, t0, resources.y,
                                      resources.inner_sunstepper,
                                      resources.sunctx);
  ASSERT_NE(resources.outer_mem, nullptr);
  ASSERT_EQ(ARKodeSetFixedStep(resources.outer_mem, SUN_RCONST(0.1)),
            ARK_SUCCESS);
  MRIStepCoupling coupling = MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_ERK33a);
  ASSERT_NE(coupling, nullptr);
  ASSERT_EQ(MRIStepSetCoupling(resources.outer_mem, coupling), ARK_SUCCESS);
  MRIStepCoupling_Free(coupling);

  ASSERT_EQ(SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM,
                                                    resources.helper, 1, 16,
                                                    SUNTRUE, resources.sunctx,
                                                    &resources.outer_scheme),
            SUN_SUCCESS);
  ASSERT_EQ(ARKodeSetAdjointCheckpointScheme(resources.outer_mem,
                                             resources.outer_scheme),
            ARK_SUCCESS);

  EXPECT_EQ(MRIStepCreateAdjointStepper(resources.outer_mem, nullptr,
                                        slow_adj_rhs, nullptr, nullptr, tf,
                                        nullptr, resources.sunctx,
                                        &resources.outer_adjoint),
            ARK_ILL_INPUT);
  EXPECT_EQ(resources.outer_adjoint, nullptr);
}

/* Fail one clone at a time to exercise ownership before and after the
   parameter vectors are transferred into the augmented ManyVector. */
static int clones_before_failure;

static N_Vector clone_with_failure(N_Vector w)
{
  if (clones_before_failure-- == 0) { return nullptr; }
  return N_VClone_Serial(w);
}

TEST(MRIStepAdjointLifecycleTest, InnerCreationCleansUpAfterCloneFailure)
{
  for (sunbooleantype has_mu : {SUNFALSE, SUNTRUE})
  {
    SCOPED_TRACE(has_mu);
    Resources resources;
    UserData data;
    data.has_mu      = has_mu;
    resources.helper = SUNMemoryHelper_Sys(resources.sunctx);
    ASSERT_NE(resources.helper, nullptr);
    resources.y = N_VNew_Serial(2, resources.sunctx);
    ASSERT_NE(resources.y, nullptr);
    N_VConst(SUN_RCONST(1.0), resources.y);
    resources.inner_mem = ERKStepCreate(fast_rhs, t0, resources.y,
                                        resources.sunctx);
    ASSERT_NE(resources.inner_mem, nullptr);
    ASSERT_EQ(ARKodeSetUserData(resources.inner_mem, &data), ARK_SUCCESS);
    ASSERT_EQ(ARKodeSetFixedStep(resources.inner_mem, SUN_RCONST(0.005)),
              ARK_SUCCESS);
    ASSERT_EQ(ARKodeCreateSUNStepper(resources.inner_mem,
                                     &resources.inner_sunstepper),
              ARK_SUCCESS);
    resources.outer_mem = MRIStepCreate(slow_rhs, nullptr, t0, resources.y,
                                        resources.inner_sunstepper,
                                        resources.sunctx);
    ASSERT_NE(resources.outer_mem, nullptr);
    ASSERT_EQ(ARKodeSetUserData(resources.outer_mem, &data), ARK_SUCCESS);
    ASSERT_EQ(ARKodeSetFixedStep(resources.outer_mem, SUN_RCONST(0.1)),
              ARK_SUCCESS);
    MRIStepCoupling coupling = MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_ERK45a);
    ASSERT_NE(coupling, nullptr);
    const int nforcing = coupling->nmat;
    ASSERT_EQ(MRIStepSetCoupling(resources.outer_mem, coupling), ARK_SUCCESS);
    MRIStepCoupling_Free(coupling);
    sunrealtype tret;
    ASSERT_GE(ARKodeEvolve(resources.outer_mem, tf, resources.y, &tret,
                           ARK_NORMAL),
              ARK_SUCCESS);
    create_terminal_state(resources, data);
    ASSERT_NE(resources.sf, nullptr);

    for (int fail_at = 0; fail_at < nforcing + has_mu + 1; fail_at++)
    {
      SCOPED_TRACE(fail_at);
      clones_before_failure          = fail_at;
      resources.lambda->ops->nvclone = clone_with_failure;
      if (has_mu) { resources.mu->ops->nvclone = clone_with_failure; }
      const int retval =
        MRIStepInnerAdjointProblem_Create(resources.outer_mem, fast_adj_rhs,
                                          resources.sf, &data,
                                          &resources.inner_problem);
      resources.lambda->ops->nvclone = N_VClone_Serial;
      if (has_mu) { resources.mu->ops->nvclone = N_VClone_Serial; }
      EXPECT_EQ(retval, ARK_MEM_FAIL);
      ASSERT_EQ(resources.inner_problem, nullptr);
      EXPECT_EQ(N_VGetSubvector_ManyVector(resources.sf, 0), resources.lambda);
      EXPECT_EQ(N_VGetArrayPointer(resources.lambda)[0],
                N_VGetArrayPointer(resources.y)[0] - SUN_RCONST(1.0));
      if (has_mu) { EXPECT_EQ(N_VMaxNorm(resources.mu), SUN_RCONST(0.0)); }
    }
    EXPECT_EQ(MRIStepInnerAdjointProblem_Create(resources.outer_mem,
                                                fast_adj_rhs, resources.sf,
                                                &data, &resources.inner_problem),
              ARK_SUCCESS);
    EXPECT_NE(resources.inner_problem, nullptr);
  }
}

TEST(MRIStepAdjointLifecycleTest, ReInitReusesWorkspaceAndCheckpointSequence)
{
  const TestConfig config{"ReInit",
                          SUNTRUE,
                          SUNTRUE,
                          1,
                          2,
                          SUN_RCONST(0.1),
                          SUN_RCONST(1.0) / SUN_RCONST(300.0)};
  Result result;
  run_mristep(config, false, result, true);
  ASSERT_TRUE(result.valid);
}

TEST(MRIStepAdjointInterpolationTest, EvolvesThroughNonStepAlignedOutput)
{
  const TestConfig config{"Interpolate",
                          SUNTRUE,
                          SUNFALSE,
                          1,
                          2,
                          SUN_RCONST(0.1),
                          SUN_RCONST(1.0) / SUN_RCONST(300.0)};
  Result interpolated;
  Result reference;
  run_mristep(config, true, interpolated, false, true);
  run_erkstep(config, reference, true);
  ASSERT_TRUE(interpolated.valid);
  ASSERT_TRUE(reference.valid);
  ASSERT_TRUE(interpolated.has_interpolated_state);
  ASSERT_TRUE(reference.has_interpolated_state);
  const sunrealtype tolerance = SUN_RCONST(10.0) * SUNRsqrt(SUN_UNIT_ROUNDOFF);
  EXPECT_LE(max_error(interpolated.interpolated_lambda,
                      reference.interpolated_lambda),
            tolerance);
  EXPECT_LE(max_error(interpolated.lambda, reference.lambda), tolerance);
}

TEST(MRIStepAdjointCheckpointTest, DoesNotReevaluateStoredSlowRhs)
{
  const TestConfig config{"StoredSlowRhs",
                          SUNTRUE,
                          SUNFALSE,
                          1,
                          2,
                          SUN_RCONST(0.1),
                          SUN_RCONST(1.0) / SUN_RCONST(300.0)};
  Result result;
  run_mristep(config, false, result, false, false, true);
  ASSERT_TRUE(result.valid);
}

TEST_P(MRIStepAdjointTest, SplitProblemIntegratesForwardAndBackward)
{
  const TestConfig& config = GetParam();
  Result result;
  run_mristep(config, false, result);
  ASSERT_TRUE(result.valid);

  constexpr std::array<sunrealtype, 2> reference_state =
    {SUN_RCONST(1.249171740707807), SUN_RCONST(0.566892772969817)};
  constexpr std::array<sunrealtype, 2> reference_lambda =
    {SUN_RCONST(2.204924168753578e-01), SUN_RCONST(-3.050319456561253e-01)};
  constexpr std::array<sunrealtype, 4> reference_mu =
    {SUN_RCONST(7.974709697811572e-02), SUN_RCONST(-5.950658873178583e-02),
     SUN_RCONST(8.224828446314870e-02), SUN_RCONST(-9.076588019985569e-02)};
  const sunrealtype tolerance =
    std::max(SUN_RCONST(5.0e-4), SUN_RCONST(10.0) * SUNRsqrt(SUN_UNIT_ROUNDOFF));

  EXPECT_NEAR(result.forward_time, tf, tolerance);
  EXPECT_NEAR(result.adjoint_time, t0, tolerance);
  EXPECT_LE(max_error(result.state, reference_state), tolerance);
  EXPECT_LE(max_error(result.lambda, reference_lambda), tolerance);
  if (config.has_mu)
  {
    EXPECT_LE(max_error(result.mu, reference_mu), tolerance);
  }
}

TEST_P(MRIStepAdjointTest, ZeroFastRhsMatchesErkAdjoint)
{
  const TestConfig& config = GetParam();
  Result mri_result;
  Result erk_result;
  run_mristep(config, true, mri_result);
  ASSERT_TRUE(mri_result.valid);
  run_erkstep(config, erk_result);
  ASSERT_TRUE(erk_result.valid);

  const sunrealtype tolerance = SUN_RCONST(20.0) * SUNRsqrt(SUN_UNIT_ROUNDOFF);
  EXPECT_NEAR(mri_result.forward_time, erk_result.forward_time, tolerance);
  EXPECT_NEAR(mri_result.adjoint_time, erk_result.adjoint_time, tolerance);
  EXPECT_EQ(mri_result.forward_steps, erk_result.forward_steps);
  EXPECT_EQ(mri_result.adjoint_steps, erk_result.adjoint_steps);
  EXPECT_LE(max_error(mri_result.state, erk_result.state), tolerance);
  EXPECT_LE(max_error(mri_result.lambda, erk_result.lambda), tolerance);
  if (config.has_mu)
  {
    EXPECT_LE(max_error(mri_result.mu, erk_result.mu), tolerance);
  }
}

/* These tables exercise repeated abscissae in the middle of a step. The
 * identity variant retains a dependency on stage 1 across stage 2. */
static MRIStepCoupling repeated_stage_table(bool identity)
{
  MRIStepCoupling table = MRIStepCoupling_Alloc(2, 4, MRISTEP_EXPLICIT);
  if (!table) { return nullptr; }
  table->q    = 1;
  table->p    = 0;
  table->c[1] = table->c[2] = SUN_RCONST(0.5);
  table->c[3]               = SUN_RCONST(1.0);
  table->W[0][1][0]         = SUN_RCONST(0.5);
  if (identity)
  {
    table->W[0][3][0] = -SUN_RCONST(0.5);
    table->W[0][3][1] = SUN_RCONST(1.0);
  }
  else
  {
    table->W[0][2][0] = -SUN_RCONST(0.5);
    table->W[1][2][0] = SUN_RCONST(0.5);
    table->W[0][2][1] = SUN_RCONST(0.5);
    table->W[1][2][1] = -SUN_RCONST(0.5);
    table->W[0][3][0] = -SUN_RCONST(0.25);
    table->W[0][3][2] = SUN_RCONST(0.75);
  }
  return table;
}

using CouplingPtr =
  std::unique_ptr<MRIStepCouplingMem, decltype(&MRIStepCoupling_Free)>;

TEST(MRIStepCouplingTest, StifflyAccurateStageIsFinalSolutionStage)
{
  for (auto type : {MRISTEP_EXPLICIT, MRISTEP_IMPLICIT, MRISTEP_IMEX})
  {
    SCOPED_TRACE(type);
    CouplingPtr table(MRIStepCoupling_Alloc(2, 4, type), &MRIStepCoupling_Free);
    ASSERT_NE(table, nullptr);
    table->q    = 1;
    table->c[1] = table->c[2] = table->c[3] = SUN_RCONST(1.0);
    if (table->W) { table->W[0][1][0] = SUN_RCONST(1.0); }
    if (table->G) { table->G[0][1][0] = SUN_RCONST(1.0); }

    EXPECT_EQ(mriStepCoupling_GetStageType(table.get(), 0), MRISTAGE_FIRST);
    EXPECT_EQ(mriStepCoupling_GetStageType(table.get(), 1), MRISTAGE_ERK_FAST);
    EXPECT_EQ(mriStepCoupling_GetStageType(table.get(), 2), MRISTAGE_ERK_NOFAST);
    EXPECT_EQ(mriStepCoupling_GetStageType(table.get(), 3), MRISTAGE_STIFF_ACC);
    EXPECT_EQ(mriStepCoupling_GetStageType(table.get(), 4), MRISTAGE_ERK_NOFAST);

    // A nonzero algebraic correction at the final time is not an identity.
    if (table->W) { table->W[1][3][0] = SUN_RCONST(1.0); }
    if (table->G) { table->G[1][3][0] = SUN_RCONST(1.0); }
    EXPECT_EQ(mriStepCoupling_GetStageType(table.get(), 3), MRISTAGE_ERK_NOFAST);
  }
}

TEST_P(MRIStepAdjointTest, AlgebraicStagesMatchDiscreteFiniteDifferences)
{
  const TestConfig& config    = GetParam();
  const sunrealtype delta     = std::cbrt(SUN_UNIT_ROUNDOFF);
  const sunrealtype tolerance = SUN_RCONST(100.0) * delta * delta;
  for (int variant = 0; variant < 3; variant++)
  {
    SCOPED_TRACE(variant);
    CouplingPtr table(variant == 0
                        ? MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_ERK22b)
                        : repeated_stage_table(variant == 2),
                      MRIStepCoupling_Free);
    ASSERT_NE(table, nullptr);
    RunOptions options;
    options.coupling = table.get();
    Result adjoint;
    run_mristep(config, false, adjoint, false, false, false, options);
    ASSERT_TRUE(adjoint.valid);
    for (int component = 0; component < (config.has_mu ? 6 : 2); component++)
    {
      SCOPED_TRACE(component);
      options.forward_only = true;
      options.perturb      = component;
      options.delta        = delta;
      Result plus;
      run_mristep(config, false, plus, false, false, false, options);
      ASSERT_TRUE(plus.valid);
      options.delta = -delta;
      Result minus;
      run_mristep(config, false, minus, false, false, false, options);
      ASSERT_TRUE(minus.valid);
      sunrealtype difference = (plus.cost - minus.cost) /
                               (SUN_RCONST(2.0) * delta);
      sunrealtype gradient = component < 2 ? adjoint.lambda[component]
                                           : adjoint.mu[component - 2];
      EXPECT_NEAR(gradient, difference, tolerance);
    }
  }
}

TEST_P(MRIStepAdjointTest, RedundantFinalStageAddsNoRhsEvaluations)
{
  const TestConfig& config = GetParam();
  CouplingPtr base(MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_ERK22b),
                   MRIStepCoupling_Free);
  ASSERT_NE(base, nullptr);
  CouplingPtr padded(MRIStepCoupling_Alloc(base->nmat, base->stages + 1,
                                           MRISTEP_EXPLICIT),
                     MRIStepCoupling_Free);
  ASSERT_NE(padded, nullptr);
  padded->q = base->q;
  padded->p = base->p;
  for (int i = 0; i < base->stages; i++)
  {
    padded->c[i] = base->c[i];
    for (int k = 0; k < base->nmat; k++)
    {
      for (int j = 0; j < base->stages; j++)
      {
        padded->W[k][i][j] = base->W[k][i][j];
      }
    }
  }
  padded->c[base->stages] = SUN_RCONST(1.0);
  // An embedding-only dependency must not trigger an adjoint evaluation.
  padded->W[0][padded->stages][base->stages - 1] = SUN_RCONST(0.25);
  padded->W[0][padded->stages][0]                = -SUN_RCONST(0.25);

  RunOptions options;
  options.coupling = base.get();
  Result original;
  run_mristep(config, false, original, false, false, false, options);
  ASSERT_TRUE(original.valid);
  options.coupling = padded.get();
  Result redundant;
  run_mristep(config, false, redundant, false, false, false, options);
  ASSERT_TRUE(redundant.valid);
  const sunrealtype tolerance = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF;
  EXPECT_LE(max_error(original.state, redundant.state), tolerance);
  EXPECT_LE(max_error(original.lambda, redundant.lambda), tolerance);
  EXPECT_LE(max_error(original.mu, redundant.mu), tolerance);
  EXPECT_EQ(original.forward_rhs_calls, redundant.forward_rhs_calls);
  EXPECT_EQ(original.reverse_rhs_calls, redundant.reverse_rhs_calls);
  EXPECT_EQ(original.reverse_rhs_calls[2], 2 * original.forward_steps);
}

static std::string config_name(const testing::TestParamInfo<TestConfig>& info)
{
  return info.param.name;
}

INSTANTIATE_TEST_SUITE_P(
  CheckpointAndStepConfigurations, MRIStepAdjointTest,
  testing::Values(
    TestConfig{"KeepLambdaCoarse", SUNTRUE, SUNFALSE, 1, 2, SUN_RCONST(0.1),
               SUN_RCONST(1.0) / SUN_RCONST(300.0)},
    TestConfig{"DropLambdaCoarse", SUNFALSE, SUNFALSE, 1, 2, SUN_RCONST(0.1),
               SUN_RCONST(1.0) / SUN_RCONST(300.0)},
    TestConfig{"KeepMuCoarse", SUNTRUE, SUNTRUE, 1, 2, SUN_RCONST(0.1),
               SUN_RCONST(1.0) / SUN_RCONST(300.0)},
    TestConfig{"DropMuCoarse", SUNFALSE, SUNTRUE, 1, 2, SUN_RCONST(0.1),
               SUN_RCONST(1.0) / SUN_RCONST(300.0)},
    TestConfig{"KeepLambdaFine", SUNTRUE, SUNFALSE, 2, 5, SUN_RCONST(0.05),
               SUN_RCONST(1.0) / SUN_RCONST(600.0)},
    TestConfig{"DropLambdaFine", SUNFALSE, SUNFALSE, 2, 5, SUN_RCONST(0.05),
               SUN_RCONST(1.0) / SUN_RCONST(600.0)},
    TestConfig{"KeepMuFine", SUNTRUE, SUNTRUE, 2, 5, SUN_RCONST(0.05),
               SUN_RCONST(1.0) / SUN_RCONST(600.0)},
    TestConfig{"DropMuFine", SUNFALSE, SUNTRUE, 2, 5, SUN_RCONST(0.05),
               SUN_RCONST(1.0) / SUN_RCONST(600.0)}),
  config_name);

} // namespace
