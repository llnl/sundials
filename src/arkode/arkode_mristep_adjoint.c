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
 * -----------------------------------------------------------------
 * MRIStep inner adjoint problem implementation.
 * -----------------------------------------------------------------*/

#include <nvector/nvector_manyvector.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>

#include "arkode_impl.h"
#include "arkode_mristep_impl.h"

static int mriStepInnerAdjointProblem_SetView(N_Vector view, N_Vector lambda,
                                              N_Vector mu)
{
  N_VectorContent_ManyVector content;

  if (!view || N_VGetVectorID(view) != SUNDIALS_NVEC_MANYVECTOR)
  {
    return ARK_ILL_INPUT;
  }

  content = (N_VectorContent_ManyVector)view->content;
  if (!content || content->num_subvectors != 2) { return ARK_ILL_INPUT; }

  content->subvec_array[0] = lambda;
  content->subvec_array[1] = mu;
  return ARK_SUCCESS;
}

static int mriStepInnerAdjointProblem_Rhs(sunrealtype t, N_Vector y,
                                          N_Vector sens, N_Vector sens_dot,
                                          void* user_data)
{
  MRIStepInnerAdjointProblem problem = (MRIStepInnerAdjointProblem)user_data;
  N_Vector lambda, lambda_dot, parameters, parameters_dot, mu, mu_dot;
  sunrealtype tau       = SUN_RCONST(0.0);
  sunrealtype tau_power = SUN_RCONST(1.0);
  int retval;

  if (!problem || !sens || !sens_dot) { return ARK_ILL_INPUT; }
  if (N_VGetVectorID(sens) != SUNDIALS_NVEC_MANYVECTOR ||
      N_VGetVectorID(sens_dot) != SUNDIALS_NVEC_MANYVECTOR ||
      N_VGetNumSubvectors_ManyVector(sens) != 2 ||
      N_VGetNumSubvectors_ManyVector(sens_dot) != 2)
  {
    return ARK_ILL_INPUT;
  }

  lambda         = N_VGetSubvector_ManyVector(sens, 0);
  lambda_dot     = N_VGetSubvector_ManyVector(sens_dot, 0);
  parameters     = N_VGetSubvector_ManyVector(sens, 1);
  parameters_dot = N_VGetSubvector_ManyVector(sens_dot, 1);

  if (!lambda || !lambda_dot || !parameters || !parameters_dot ||
      N_VGetVectorID(parameters) != SUNDIALS_NVEC_MANYVECTOR ||
      N_VGetVectorID(parameters_dot) != SUNDIALS_NVEC_MANYVECTOR ||
      N_VGetNumSubvectors_ManyVector(parameters) != problem->nforcing + 1 ||
      N_VGetNumSubvectors_ManyVector(parameters_dot) != problem->nforcing + 1)
  {
    return ARK_ILL_INPUT;
  }

  if (problem->nforcing > 1)
  {
    if (problem->tscale == SUN_RCONST(0.0)) { return ARK_ILL_INPUT; }
    tau = (t - problem->tshift) / problem->tscale;
  }

  mu     = N_VGetSubvector_ManyVector(parameters, 0);
  mu_dot = N_VGetSubvector_ManyVector(parameters_dot, 0);
  retval = mriStepInnerAdjointProblem_SetView(problem->sens_view, lambda, mu);
  if (retval != ARK_SUCCESS) { return retval; }
  retval = mriStepInnerAdjointProblem_SetView(problem->sens_dot_view,
                                              lambda_dot, mu_dot);
  if (retval != ARK_SUCCESS) { return retval; }

  retval = problem->adj_f(t, y, problem->sens_view, problem->sens_dot_view,
                          problem->user_data);
  if (retval != ARK_SUCCESS) { return retval; }

  for (int k = 0; k < problem->nforcing; k++)
  {
    N_Vector omega_dot = N_VGetSubvector_ManyVector(parameters_dot, k + 1);
    if (!omega_dot) { return ARK_ILL_INPUT; }
    N_VScale(tau_power, lambda, omega_dot);
    tau_power *= tau;
  }

  return ARK_SUCCESS;
}

int MRIStepInnerAdjointProblem_Create(void* arkode_mem, SUNAdjRhsFn adj_f,
                                      N_Vector sf, void* user_data,
                                      MRIStepInnerAdjointProblem* problem_ptr)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  MRIStepInnerAdjointProblem problem = NULL;
  N_Vector lambda_f, mu_f, lambda, parameters;
  N_Vector* parameter_vectors = NULL;
  N_Vector outer_vectors[2];
  N_Vector view_vectors[2];
  int retval;

  if (!problem_ptr)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The output problem pointer is NULL");
    return ARK_ILL_INPUT;
  }
  *problem_ptr = NULL;

  retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }
  if (ark_mem->step_init != mriStep_Init)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The ARKODE memory was not created by MRIStep");
    return ARK_ILL_INPUT;
  }
  if (!adj_f)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The adjoint right-hand side function is NULL");
    return ARK_ILL_INPUT;
  }
  if (!sf || N_VGetVectorID(sf) != SUNDIALS_NVEC_MANYVECTOR ||
      N_VGetNumSubvectors_ManyVector(sf) != 2)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The terminal state must be a two-subvector ManyVector");
    return ARK_ILL_INPUT;
  }
  if (step_mem->ninner_forcing < 0)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The number of inner forcing vectors is invalid");
    return ARK_ILL_INPUT;
  }

  lambda_f = N_VGetSubvector_ManyVector(sf, 0);
  mu_f     = N_VGetSubvector_ManyVector(sf, 1);
  if (!lambda_f || !mu_f)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The terminal state contains a NULL subvector");
    return ARK_ILL_INPUT;
  }

  problem = (MRIStepInnerAdjointProblem)
    calloc(1, sizeof(struct MRIStepInnerAdjointProblem_));
  if (!problem) { return ARK_MEM_FAIL; }

  problem->adj_f     = adj_f;
  problem->user_data = user_data;
  problem->tshift    = step_mem->inner_tshift;
  problem->tscale    = step_mem->inner_tscale;
  problem->nforcing  = step_mem->ninner_forcing;

  parameter_vectors = (N_Vector*)calloc((size_t)problem->nforcing + 1,
                                        sizeof(N_Vector));
  if (!parameter_vectors)
  {
    MRIStepInnerAdjointProblem_Free(problem);
    return ARK_MEM_FAIL;
  }
  parameter_vectors[0] = N_VClone(mu_f);
  if (!parameter_vectors[0])
  {
    N_VDestroyVectorArray(parameter_vectors, problem->nforcing + 1);
    MRIStepInnerAdjointProblem_Free(problem);
    return ARK_MEM_FAIL;
  }
  N_VScale(SUN_RCONST(1.0), mu_f, parameter_vectors[0]);

  for (int k = 0; k < problem->nforcing; k++)
  {
    parameter_vectors[k + 1] = N_VClone(lambda_f);
    if (!parameter_vectors[k + 1])
    {
      N_VDestroyVectorArray(parameter_vectors, problem->nforcing + 1);
      MRIStepInnerAdjointProblem_Free(problem);
      return ARK_MEM_FAIL;
    }
    N_VConst(SUN_RCONST(0.0), parameter_vectors[k + 1]);
  }

  parameters = N_VNew_ManyVector(problem->nforcing + 1, parameter_vectors,
                                 ark_mem->sunctx);
  if (!parameters)
  {
    N_VDestroyVectorArray(parameter_vectors, problem->nforcing + 1);
    MRIStepInnerAdjointProblem_Free(problem);
    return ARK_MEM_FAIL;
  }
  ((N_VectorContent_ManyVector)parameters->content)->own_data = SUNTRUE;
  free(parameter_vectors);

  lambda = N_VClone(lambda_f);
  if (!lambda)
  {
    N_VDestroy(parameters);
    MRIStepInnerAdjointProblem_Free(problem);
    return ARK_MEM_FAIL;
  }
  N_VScale(SUN_RCONST(1.0), lambda_f, lambda);

  outer_vectors[0] = lambda;
  outer_vectors[1] = parameters;
  problem->sf      = N_VNew_ManyVector(2, outer_vectors, ark_mem->sunctx);
  if (!problem->sf)
  {
    N_VDestroy(parameters);
    N_VDestroy(lambda);
    MRIStepInnerAdjointProblem_Free(problem);
    return ARK_MEM_FAIL;
  }
  ((N_VectorContent_ManyVector)problem->sf->content)->own_data = SUNTRUE;

  view_vectors[0]        = lambda;
  view_vectors[1]        = N_VGetSubvector_ManyVector(parameters, 0);
  problem->sens_view     = N_VNew_ManyVector(2, view_vectors, ark_mem->sunctx);
  problem->sens_dot_view = N_VNew_ManyVector(2, view_vectors, ark_mem->sunctx);
  if (!problem->sens_view || !problem->sens_dot_view)
  {
    MRIStepInnerAdjointProblem_Free(problem);
    return ARK_MEM_FAIL;
  }

  *problem_ptr = problem;
  return ARK_SUCCESS;
}

int MRIStepInnerAdjointProblem_GetAdjRhsFn(MRIStepInnerAdjointProblem problem,
                                           SUNAdjRhsFn* adj_f)
{
  if (!problem || !adj_f) { return ARK_ILL_INPUT; }
  *adj_f = mriStepInnerAdjointProblem_Rhs;
  return ARK_SUCCESS;
}

int MRIStepInnerAdjointProblem_GetTerminalState(MRIStepInnerAdjointProblem problem,
                                                N_Vector* sf)
{
  if (!problem || !sf) { return ARK_ILL_INPUT; }
  *sf = problem->sf;
  return ARK_SUCCESS;
}

int MRIStepInnerAdjointProblem_GetUserData(MRIStepInnerAdjointProblem problem,
                                           void** user_data)
{
  if (!problem || !user_data) { return ARK_ILL_INPUT; }
  *user_data = problem;
  return ARK_SUCCESS;
}

void MRIStepInnerAdjointProblem_Free(MRIStepInnerAdjointProblem problem)
{
  if (!problem) { return; }

  N_VDestroy(problem->sens_view);
  N_VDestroy(problem->sens_dot_view);
  N_VDestroy(problem->sf);

  free(problem);
}
