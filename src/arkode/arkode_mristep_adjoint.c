/* ---------------------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 * ---------------------------------------------------------------------------
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
 * ---------------------------------------------------------------------------
 * MRIStep adjoint implementation.
 * ---------------------------------------------------------------------------*/

#include <nvector/nvector_manyvector.h>
#include <stdlib.h>
#include <sundials/sundials_adjointcheckpointscheme.h>
#include <sundials/sundials_math.h>

#include "arkode_impl.h"
#include "arkode_mristep_impl.h"
#include "sundials_adjointstepper_impl.h"

typedef struct
{
  N_Vector lambda;
  N_Vector mu;
  sunindextype count;
} MRIStepAdjointState;

typedef struct
{
  N_Vector lambda;
  N_Vector mu;
  N_Vector* omega;
} MRIStepInnerAdjointState;

static int mriStepAdjointGetOuterState(N_Vector state, MRIStepAdjointState* view)
{
  if (!state || !view || N_VGetVectorID(state) != SUNDIALS_NVEC_MANYVECTOR)
  {
    return ARK_ILL_INPUT;
  }

  view->count = N_VGetNumSubvectors_ManyVector(state);
  if (view->count != 1 && view->count != 2) { return ARK_ILL_INPUT; }

  view->lambda = N_VGetSubvector_ManyVector(state, 0);
  view->mu = (view->count == 2) ? N_VGetSubvector_ManyVector(state, 1) : NULL;
  if (!view->lambda || (view->count == 2 && !view->mu))
  {
    return ARK_ILL_INPUT;
  }
  return ARK_SUCCESS;
}

static int mriStepAdjointGetInnerState(N_Vector state, sunindextype nsens,
                                       int nforcing,
                                       MRIStepInnerAdjointState* view)
{
  N_VectorContent_ManyVector parameters_content;
  N_Vector parameters;
  int nparameters;

  if (!state || !view || (nsens != 1 && nsens != 2) ||
      N_VGetVectorID(state) != SUNDIALS_NVEC_MANYVECTOR ||
      N_VGetNumSubvectors_ManyVector(state) != 2)
  {
    return ARK_ILL_INPUT;
  }

  view->lambda = N_VGetSubvector_ManyVector(state, 0);
  parameters   = N_VGetSubvector_ManyVector(state, 1);
  nparameters  = nforcing + (int)nsens - 1;
  if (!view->lambda || !parameters ||
      N_VGetVectorID(parameters) != SUNDIALS_NVEC_MANYVECTOR ||
      N_VGetNumSubvectors_ManyVector(parameters) != nparameters)
  {
    return ARK_ILL_INPUT;
  }

  parameters_content = (N_VectorContent_ManyVector)parameters->content;
  if (!parameters_content || !parameters_content->subvec_array)
  {
    return ARK_ILL_INPUT;
  }
  view->omega = parameters_content->subvec_array;
  view->mu = (nsens == 2) ? parameters_content->subvec_array[nforcing] : NULL;
  if (nsens == 2 && !view->mu) { return ARK_ILL_INPUT; }
  for (int k = 0; k < nforcing; k++)
  {
    if (!view->omega[k]) { return ARK_ILL_INPUT; }
  }
  return ARK_SUCCESS;
}

static int mriStepInnerAdjointProblem_SetView(N_Vector view, N_Vector lambda,
                                              N_Vector mu)
{
  N_VectorContent_ManyVector content;

  if (!view || N_VGetVectorID(view) != SUNDIALS_NVEC_MANYVECTOR)
  {
    return ARK_ILL_INPUT;
  }

  content = (N_VectorContent_ManyVector)view->content;
  if (!content || (content->num_subvectors != 1 && content->num_subvectors != 2) ||
      !lambda || (content->num_subvectors == 2 && !mu))
  {
    return ARK_ILL_INPUT;
  }

  content->subvec_array[0] = lambda;
  if (content->num_subvectors == 2) { content->subvec_array[1] = mu; }
  return ARK_SUCCESS;
}

static int mriStepInnerAdjointProblem_Rhs(sunrealtype t, N_Vector y,
                                          N_Vector sens, N_Vector sens_dot,
                                          void* user_data)
{
  MRIStepInnerAdjointProblem problem = (MRIStepInnerAdjointProblem)user_data;
  MRIStepInnerAdjointState state, state_dot;
  sunrealtype tau       = SUN_RCONST(0.0);
  sunrealtype tau_power = SUN_RCONST(1.0);
  sunindextype nsens;
  int retval;

  if (!problem || !sens || !sens_dot) { return ARK_ILL_INPUT; }
  nsens  = N_VGetNumSubvectors_ManyVector(problem->sens_view);
  retval = mriStepAdjointGetInnerState(sens, nsens, problem->nforcing, &state);
  if (retval != ARK_SUCCESS) { return retval; }
  retval = mriStepAdjointGetInnerState(sens_dot, nsens, problem->nforcing,
                                       &state_dot);
  if (retval != ARK_SUCCESS) { return retval; }

  if (problem->nforcing > 1)
  {
    if (problem->tscale == SUN_RCONST(0.0)) { return ARK_ILL_INPUT; }
    tau = (t - problem->tshift) / problem->tscale;
  }

  retval = mriStepInnerAdjointProblem_SetView(problem->sens_view, state.lambda,
                                              state.mu);
  if (retval != ARK_SUCCESS) { return retval; }
  retval = mriStepInnerAdjointProblem_SetView(problem->sens_dot_view,
                                              state_dot.lambda, state_dot.mu);
  if (retval != ARK_SUCCESS) { return retval; }

  retval = problem->adj_f(t, y, problem->sens_view, problem->sens_dot_view,
                          problem->user_data);
  if (retval != ARK_SUCCESS) { return retval; }

  for (int k = 0; k < problem->nforcing; k++)
  {
    N_VScale(tau_power, state.lambda, state_dot.omega[k]);
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
  sunindextype nsens;
  int nparameters;
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
  if (!sf || N_VGetVectorID(sf) != SUNDIALS_NVEC_MANYVECTOR)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The terminal state must be a ManyVector");
    return ARK_ILL_INPUT;
  }
  nsens = N_VGetNumSubvectors_ManyVector(sf);
  if (nsens != 1 && nsens != 2)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "The terminal state must contain lambda and optional mu subvectors");
    return ARK_ILL_INPUT;
  }
  if (step_mem->ninner_forcing < 0)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The number of inner forcing vectors is invalid");
    return ARK_ILL_INPUT;
  }

  lambda_f = N_VGetSubvector_ManyVector(sf, 0);
  mu_f     = (nsens == 2) ? N_VGetSubvector_ManyVector(sf, 1) : NULL;
  if (!lambda_f || (nsens == 2 && !mu_f))
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
  nparameters        = problem->nforcing + (int)nsens - 1;

  parameter_vectors = (N_Vector*)calloc((size_t)nparameters, sizeof(N_Vector));
  if (!parameter_vectors)
  {
    MRIStepInnerAdjointProblem_Free(problem);
    return ARK_MEM_FAIL;
  }
  for (int k = 0; k < problem->nforcing; k++)
  {
    parameter_vectors[k] = N_VClone(lambda_f);
    if (!parameter_vectors[k])
    {
      N_VDestroyVectorArray(parameter_vectors, nparameters);
      MRIStepInnerAdjointProblem_Free(problem);
      return ARK_MEM_FAIL;
    }
    N_VConst(SUN_RCONST(0.0), parameter_vectors[k]);
  }
  if (nsens == 2)
  {
    parameter_vectors[problem->nforcing] = N_VClone(mu_f);
    if (!parameter_vectors[problem->nforcing])
    {
      N_VDestroyVectorArray(parameter_vectors, nparameters);
      MRIStepInnerAdjointProblem_Free(problem);
      return ARK_MEM_FAIL;
    }
    N_VScale(SUN_RCONST(1.0), mu_f, parameter_vectors[problem->nforcing]);
  }

  parameters = N_VNew_ManyVector(nparameters, parameter_vectors, ark_mem->sunctx);
  if (!parameters)
  {
    N_VDestroyVectorArray(parameter_vectors, nparameters);
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

  view_vectors[0] = lambda;
  if (nsens == 2)
  {
    view_vectors[1] = N_VGetSubvector_ManyVector(parameters, problem->nforcing);
  }
  problem->sens_view = N_VNew_ManyVector(nsens, view_vectors, ark_mem->sunctx);
  problem->sens_dot_view = N_VNew_ManyVector(nsens, view_vectors,
                                             ark_mem->sunctx);
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

static int mriStep_fse_Adj(sunrealtype t, N_Vector sens_partial_stage,
                           N_Vector sens_complete_stage, void* content)
{
  SUNAdjointStepper adj_stepper = (SUNAdjointStepper)content;
  if (!adj_stepper || !adj_stepper->adj_sunstepper) { return ARK_MEM_NULL; }

  ARKodeMem ark_mem = (ARKodeMem)adj_stepper->adj_sunstepper->content;
  if (!ark_mem || !ark_mem->step_mem) { return ARK_MEM_NULL; }

  ARKodeMRIStepMem step_mem   = (ARKodeMRIStepMem)ark_mem->step_mem;
  MRIStepAdjointData adj_data = step_mem->adj_data;
  if (!adj_data || !adj_data->fse || !adj_data->stage_states ||
      ark_mem->adj_stage_idx < 0 ||
      ark_mem->adj_stage_idx >= step_mem->stages + 1)
  {
    return ARK_ILL_INPUT;
  }

  return adj_data->fse(t, adj_data->stage_states[ark_mem->adj_stage_idx],
                       sens_partial_stage, sens_complete_stage,
                       adj_stepper->user_data);
}

static SUNErrCode mriStepAdjointSUNStepperReInit(SUNStepper stepper,
                                                 sunrealtype t0, N_Vector y0)
{
  void* arkode_mem = NULL;
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  if (SUNStepper_GetContent(stepper, &arkode_mem) != SUN_SUCCESS)
  {
    return SUN_ERR_OP_FAIL;
  }
  retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return SUN_ERR_OP_FAIL; }

  stepper->last_flag = MRIStepReInit(arkode_mem, step_mem->fse, step_mem->fsi,
                                     t0, y0);
  return (stepper->last_flag == ARK_SUCCESS) ? SUN_SUCCESS : SUN_ERR_OP_FAIL;
}

void mriStep_FreeAdjointData(ARKodeMRIStepMem step_mem)
{
  if (!step_mem || !step_mem->adj_data) { return; }

  MRIStepAdjointData adj_data = step_mem->adj_data;
  int nstages                 = step_mem->stages + 1;

  N_VDestroyVectorArray(adj_data->stage_states, nstages);
  N_VDestroyVectorArray(adj_data->slow_accum, nstages);
  free(adj_data->stage_times);
  N_VDestroy(adj_data->inner_state);
  free(adj_data);
  step_mem->adj_data = NULL;
}

static int mriStepAdjointCreateData(ARKodeMRIStepMem step_mem,
                                    SUNAdjRhsFn adj_fse,
                                    SUNAdjointStepper inner_stepper,
                                    N_Vector state_template)
{
  int nstages = step_mem->stages + 1;

  MRIStepAdjointData adj_data = (MRIStepAdjointData)calloc(1, sizeof(*adj_data));
  if (!adj_data) { return ARK_MEM_FAIL; }
  step_mem->adj_data = adj_data;

  adj_data->fse                  = adj_fse;
  adj_data->inner_stepper        = inner_stepper;
  adj_data->inner_checkpoint_idx = 0;
  adj_data->stage_states = (N_Vector*)calloc((size_t)nstages, sizeof(N_Vector));
  adj_data->slow_accum   = (N_Vector*)calloc((size_t)nstages, sizeof(N_Vector));
  adj_data->stage_times  = (sunrealtype*)calloc((size_t)nstages,
                                                sizeof(sunrealtype));
  if (!adj_data->stage_states || !adj_data->slow_accum || !adj_data->stage_times)
  {
    mriStep_FreeAdjointData(step_mem);
    return ARK_MEM_FAIL;
  }

  for (int i = 0; i < nstages; i++)
  {
    adj_data->stage_states[i] = N_VClone(state_template);
    adj_data->slow_accum[i]   = N_VClone(state_template);
    if (!adj_data->stage_states[i] || !adj_data->slow_accum[i])
    {
      mriStep_FreeAdjointData(step_mem);
      return ARK_MEM_FAIL;
    }
  }
  adj_data->inner_state = N_VClone(state_template);
  if (!adj_data->inner_state)
  {
    mriStep_FreeAdjointData(step_mem);
    return ARK_MEM_FAIL;
  }
  return ARK_SUCCESS;
}

static int mriStepAdjointLoadStages(ARKodeMem ark_mem,
                                    MRIStepAdjointData adj_data,
                                    ARKodeMRIStepMem fwd_step_mem,
                                    suncountertype step_idx,
                                    sunbooleantype* missing)
{
  ARKodeMRIStepMem step_mem = (ARKodeMRIStepMem)ark_mem->step_mem;
  suncountertype slot       = 0;

  *missing = SUNFALSE;
  for (int stage = 0; stage < step_mem->stages; stage++)
  {
    SUNErrCode err =
      SUNAdjointCheckpointScheme_LoadVector(ark_mem->checkpoint_scheme,
                                            step_idx, slot++, SUNTRUE,
                                            &adj_data->stage_states[stage],
                                            &adj_data->stage_times[stage]);
    if (err == SUN_ERR_CHECKPOINT_NOT_FOUND)
    {
      *missing = SUNTRUE;
      return ARK_SUCCESS;
    }
    if (err != SUN_SUCCESS) { return ARK_ADJ_CHECKPOINT_FAIL; }

    if (mriStep_AdjointNeedsRhsCheckpoint(step_mem, stage))
    {
      sunrealtype rhs_time;
      int map = fwd_step_mem->stage_map[stage];
      err = SUNAdjointCheckpointScheme_LoadVector(ark_mem->checkpoint_scheme,
                                                  step_idx, slot++, SUNTRUE,
                                                  &fwd_step_mem->Fse[map],
                                                  &rhs_time);
      if (err == SUN_ERR_CHECKPOINT_NOT_FOUND)
      {
        *missing = SUNTRUE;
        return ARK_SUCCESS;
      }
      if (err != SUN_SUCCESS) { return ARK_ADJ_CHECKPOINT_FAIL; }
    }
  }

  SUNErrCode err =
    SUNAdjointCheckpointScheme_LoadVector(ark_mem->checkpoint_scheme, step_idx,
                                          slot, SUNTRUE,
                                          &adj_data->stage_states[step_mem->stages],
                                          &adj_data->stage_times[step_mem->stages]);
  if (err == SUN_ERR_CHECKPOINT_NOT_FOUND)
  {
    *missing = SUNTRUE;
    return ARK_SUCCESS;
  }
  if (err != SUN_SUCCESS) { return ARK_ADJ_CHECKPOINT_FAIL; }

  return ARK_SUCCESS;
}

static int mriStepAdjointRecomputeOuter(ARKodeMem ark_mem,
                                        MRIStepAdjointData adj_data,
                                        SUNAdjointStepper adj_stepper,
                                        suncountertype curr_step)
{
  ARKodeMRIStepMem step_mem    = (ARKodeMRIStepMem)ark_mem->step_mem;
  N_Vector checkpoint          = adj_data->stage_states[0];
  suncountertype terminal_slot = 0;
  SUNAdjointCheckpointScheme inner_scheme =
    adj_data->inner_stepper->checkpoint_scheme;

  for (int stage = 0; stage < step_mem->stages; stage++)
  {
    terminal_slot++;
    if (mriStep_AdjointNeedsRhsCheckpoint(step_mem, stage)) { terminal_slot++; }
  }

  for (suncountertype start_step = curr_step; start_step >= 0; start_step--)
  {
    sunrealtype checkpoint_t = SUN_RCONST(0.0);
    SUNErrCode err =
      SUNAdjointCheckpointScheme_LoadVector(ark_mem->checkpoint_scheme,
                                            start_step, terminal_slot, SUNTRUE,
                                            &checkpoint, &checkpoint_t);
    if (err == SUN_ERR_CHECKPOINT_NOT_FOUND) { continue; }
    if (err != SUN_SUCCESS) { return ARK_ADJ_CHECKPOINT_FAIL; }

    err = SUNAdjointCheckpointScheme_Enable(inner_scheme, SUNFALSE);
    if (err != SUN_SUCCESS) { return ARK_ADJ_CHECKPOINT_FAIL; }
    err = SUNAdjointStepper_RecomputeFwd(adj_stepper, start_step + 1,
                                         checkpoint_t, checkpoint, ark_mem->tn);
    SUNErrCode enable_err = SUNAdjointCheckpointScheme_Enable(inner_scheme,
                                                              SUNTRUE);
    if (err != SUN_SUCCESS || enable_err != SUN_SUCCESS)
    {
      return ARK_ADJ_RECOMPUTE_FAIL;
    }
    return ARK_SUCCESS;
  }
  return ARK_ADJ_RECOMPUTE_FAIL;
}

static int mriStepAdjointConsumeStages(ARKodeMem ark_mem,
                                       MRIStepAdjointData adj_data,
                                       suncountertype step_idx)
{
  ARKodeMRIStepMem step_mem = (ARKodeMRIStepMem)ark_mem->step_mem;
  sunrealtype t;
  suncountertype slot = 1; /* include the terminal state */
  for (int stage = 0; stage < step_mem->stages; stage++)
  {
    slot++;
    if (mriStep_AdjointNeedsRhsCheckpoint(step_mem, stage)) { slot++; }
  }
  while (slot > 0)
  {
    slot--;
    SUNErrCode err =
      SUNAdjointCheckpointScheme_LoadVector(ark_mem->checkpoint_scheme,
                                            step_idx, slot, SUNFALSE,
                                            &adj_data->inner_state, &t);
    if (err != SUN_SUCCESS) { return ARK_ADJ_CHECKPOINT_FAIL; }
  }
  return ARK_SUCCESS;
}

static int mriStepAdjointSolveInner(
  ARKodeMem fwd_ark_mem, MRIStepAdjointData adj_data,
  ARKodeMRIStepMem fwd_step_mem, MRIStepInnerAdjointProblem inner_problem,
  MRIStepAdjointState* sens, MRIStepInnerAdjointState* inner_sens, int stage)
{
  sunrealtype t0 = adj_data->stage_times[stage - 1];
  sunrealtype tf = adj_data->stage_times[stage];
  int retval = mriStep_ComputeInnerForcing(fwd_ark_mem, fwd_step_mem, stage, t0,
                                           tf);
  if (retval != ARK_SUCCESS) { return retval; }

  inner_problem->tshift = t0;
  inner_problem->tscale = tf - t0;
  N_VScale(SUN_RCONST(1.0), sens->lambda, inner_sens->lambda);
  if (sens->mu) { N_VScale(SUN_RCONST(1.0), sens->mu, inner_sens->mu); }
  for (int k = 0; k < inner_problem->nforcing; k++)
  {
    N_VConst(SUN_RCONST(0.0), inner_sens->omega[k]);
  }

  SUNErrCode err = SUNStepper_SetForcing(adj_data->inner_stepper->fwd_sunstepper,
                                         t0, tf - t0, fwd_step_mem->inner_forcing,
                                         fwd_step_mem->ninner_forcing);
  if (err != SUN_SUCCESS) { return ARK_SUNSTEPPER_ERR; }

  N_VScale(SUN_RCONST(1.0), adj_data->stage_states[stage - 1],
           adj_data->inner_state);
  suncountertype nst_before, nst_after;
  err = SUNStepper_GetNumSteps(adj_data->inner_stepper->fwd_sunstepper,
                               &nst_before);
  if (err == SUN_SUCCESS)
  {
    err = SUNAdjointStepper_RecomputeFwd(adj_data->inner_stepper,
                                         adj_data->inner_checkpoint_idx, t0,
                                         adj_data->inner_state, tf);
  }
  if (err == SUN_SUCCESS)
  {
    err = SUNStepper_GetNumSteps(adj_data->inner_stepper->fwd_sunstepper,
                                 &nst_after);
  }
  SUNErrCode clear_err =
    SUNStepper_SetForcing(adj_data->inner_stepper->fwd_sunstepper,
                          SUN_RCONST(0.0), SUN_RCONST(1.0), NULL, 0);
  if (err != SUN_SUCCESS || clear_err != SUN_SUCCESS)
  {
    return ARK_ADJ_RECOMPUTE_FAIL;
  }

  suncountertype interval_steps = nst_after - nst_before;
  if (interval_steps <= 0) { return ARK_ADJ_RECOMPUTE_FAIL; }
  suncountertype final_inner_idx = adj_data->inner_checkpoint_idx +
                                   interval_steps - 1;
  adj_data->inner_checkpoint_idx += interval_steps;

  err = SUNAdjointStepper_ReInit(adj_data->inner_stepper, tf, inner_problem->sf,
                                 final_inner_idx);
  if (err != SUN_SUCCESS) { return ARK_SUNADJSTEPPER_ERR; }
  err = SUNStepper_SetStopTime(adj_data->inner_stepper->adj_sunstepper, t0);
  if (err != SUN_SUCCESS) { return ARK_SUNADJSTEPPER_ERR; }

  sunrealtype tret = tf;
  for (suncountertype i = 0; i < interval_steps; i++)
  {
    err = SUNAdjointStepper_OneStep(adj_data->inner_stepper, t0,
                                    inner_problem->sf, &tret);
    if (err != SUN_SUCCESS) { return ARK_SUNADJSTEPPER_ERR; }
  }
  return ARK_SUCCESS;
}

static int mriStepAdjointAccumulateStage(
  ARKodeMem ark_mem, ARKodeMRIStepMem step_mem, ARKodeMem fwd_ark_mem,
  ARKodeMRIStepMem fwd_step_mem, MRIStepAdjointData adj_data,
  MRIStepAdjointState* sens, MRIStepAdjointState* slow_input,
  MRIStepAdjointState* slow_output, MRIStepInnerAdjointState* inner_sens,
  int stage)
{
  sunrealtype t0     = adj_data->stage_times[stage - 1];
  sunrealtype tf     = adj_data->stage_times[stage];
  sunrealtype rcdiff = fwd_ark_mem->h / (tf - t0);

  for (int j = 0; j < stage; j++)
  {
    if (fwd_step_mem->stage_map[j] < 0) { continue; }
    for (int k = 0; k < fwd_step_mem->MRIC->nmat; k++)
    {
      sunrealtype coefficient = rcdiff * fwd_step_mem->MRIC->W[k][stage][j];
      N_VLinearSum(SUN_RCONST(1.0), adj_data->slow_accum[j], coefficient,
                   inner_sens->omega[k], adj_data->slow_accum[j]);
    }
  }

  int current_stage = stage - 1;
  N_VScale(SUN_RCONST(1.0), adj_data->slow_accum[current_stage],
           slow_input->lambda);
  if (sens->mu) { N_VConst(SUN_RCONST(0.0), slow_input->mu); }

  ark_mem->adj_stage_idx = current_stage;
  int retval             = step_mem->fse(adj_data->stage_times[current_stage],
                                         ark_mem->tempv2, ark_mem->tempv3,
                                         ark_mem->user_data);
  if (retval != 0) { return ARK_RHSFUNC_FAIL; }

  N_VLinearSum(SUN_RCONST(1.0), inner_sens->lambda, SUN_RCONST(1.0),
               slow_output->lambda, sens->lambda);
  if (sens->mu)
  {
    N_VLinearSum(SUN_RCONST(1.0), inner_sens->mu, SUN_RCONST(1.0),
                 slow_output->mu, sens->mu);
  }
  return ARK_SUCCESS;
}

int mriStep_TakeStep_Adjoint(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  ARKodeMRIStepMem step_mem;
  MRIStepAdjointData adj_data;
  ARKodeMem fwd_ark_mem;
  ARKodeMRIStepMem fwd_step_mem;
  SUNAdjointStepper adj_stepper;
  MRIStepInnerAdjointProblem inner_problem;
  MRIStepAdjointState sens, slow_input, slow_output;
  MRIStepInnerAdjointState inner_sens;
  suncountertype step_idx;
  int retval;

  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }
  adj_data    = step_mem->adj_data;
  adj_stepper = (SUNAdjointStepper)ark_mem->user_data;
  if (!adj_data || !adj_data->inner_stepper) { return ARK_MEM_NULL; }
  inner_problem = (MRIStepInnerAdjointProblem)adj_data->inner_stepper->user_data;
  if (!adj_stepper) { return ARK_MEM_NULL; }
  if (!inner_problem) { return ARK_MEM_NULL; }
  void* fwd_arkode_mem = NULL;
  SUNErrCode err       = SUNStepper_GetContent(adj_stepper->fwd_sunstepper,
                                               &fwd_arkode_mem);
  if (err != SUN_SUCCESS) { return ARK_SUNSTEPPER_ERR; }
  retval = mriStep_AccessARKODEStepMem(fwd_arkode_mem, __func__, &fwd_ark_mem,
                                       &fwd_step_mem);
  if (retval != ARK_SUCCESS) { return retval; }
  retval = mriStepAdjointGetOuterState(ark_mem->ycur, &sens);
  if (retval != ARK_SUCCESS) { return retval; }
  retval = mriStepAdjointGetOuterState(ark_mem->tempv2, &slow_input);
  if (retval != ARK_SUCCESS) { return retval; }
  retval = mriStepAdjointGetOuterState(ark_mem->tempv3, &slow_output);
  if (retval != ARK_SUCCESS || slow_input.count != sens.count ||
      slow_output.count != sens.count)
  {
    return ARK_ILL_INPUT;
  }
  retval = mriStepAdjointGetInnerState(inner_problem->sf, sens.count,
                                       inner_problem->nforcing, &inner_sens);
  if (retval != ARK_SUCCESS) { return retval; }

  step_idx  = adj_stepper->final_step_idx - ark_mem->nst;
  *dsmPtr   = SUN_RCONST(0.0);
  *nflagPtr = ARK_SUCCESS;

  sunbooleantype missing;
  retval = mriStepAdjointLoadStages(ark_mem, adj_data, fwd_step_mem, step_idx,
                                    &missing);
  if (retval != ARK_SUCCESS) { return retval; }
  if (missing)
  {
    retval = mriStepAdjointRecomputeOuter(ark_mem, adj_data, adj_stepper,
                                          step_idx);
    if (retval != ARK_SUCCESS) { return retval; }
    retval = mriStepAdjointLoadStages(ark_mem, adj_data, fwd_step_mem, step_idx,
                                      &missing);
    if (retval != ARK_SUCCESS || missing) { return ARK_ADJ_CHECKPOINT_FAIL; }
  }

  N_VScale(SUN_RCONST(1.0), ark_mem->yn, ark_mem->ycur);
  for (int i = 0; i < step_mem->stages; i++)
  {
    N_VConst(SUN_RCONST(0.0), adj_data->slow_accum[i]);
  }

  for (int stage = step_mem->stages - 1; stage >= 1; stage--)
  {
    retval = mriStepAdjointSolveInner(fwd_ark_mem, adj_data, fwd_step_mem,
                                      inner_problem, &sens, &inner_sens, stage);
    if (retval != ARK_SUCCESS) { return retval; }
    retval = mriStepAdjointAccumulateStage(ark_mem, step_mem, fwd_ark_mem,
                                           fwd_step_mem, adj_data, &sens,
                                           &slow_input, &slow_output,
                                           &inner_sens, stage);
    if (retval != ARK_SUCCESS) { return retval; }
  }

  retval = mriStepAdjointConsumeStages(ark_mem, adj_data, step_idx);
  return retval;
}

static int mriStepCompatibleWithAdjointSolver(ARKodeMem ark_mem,
                                              ARKodeMRIStepMem step_mem,
                                              int lineno, const char* fname,
                                              const char* filename)
{
  if (!ark_mem->fixedstep)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, lineno, fname, filename,
                    "MRIStep must be using a fixed step to work with "
                    "SUNAdjointStepper");
    return ARK_ILL_INPUT;
  }
  if (ark_mem->checkpoint_scheme == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, lineno, fname,
                    filename, "MRIStep requires a checkpoint scheme to work with SUNAdjointStepper");
    return ARK_ILL_INPUT;
  }
  if (!step_mem->MRIC || step_mem->MRIC->type != MRISTEP_EXPLICIT)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, lineno, fname,
                    filename, "SUNAdjointStepper only supports explicit MRI-GARK methods");
    return ARK_ILL_INPUT;
  }
  if (ark_mem->relax_enabled)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, lineno, fname, filename,
                    "SUNAdjointStepper is not compatible with relaxation");
    return ARK_ILL_INPUT;
  }
  if (ark_mem->constraints)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, lineno, fname, filename,
                    "SUNAdjointStepper is not compatible with constraints");
    return ARK_ILL_INPUT;
  }
  return ARK_SUCCESS;
}

static int mriStepAdjointCreateARKodeMem(ARKodeMem fwd_ark_mem,
                                         ARKodeMRIStepMem fwd_step_mem,
                                         SUNAdjointStepper inner_stepper,
                                         SUNAdjRhsFn adj_fse, sunrealtype tf,
                                         N_Vector sf, long int nst,
                                         void** adj_mem_ptr)
{
  *adj_mem_ptr = MRIStepCreate(mriStep_fse_Adj, NULL, tf, sf,
                               inner_stepper->fwd_sunstepper,
                               fwd_ark_mem->sunctx);
  if (!*adj_mem_ptr) { return ARK_MEM_FAIL; }

  ARKodeMem adj_ark_mem         = (ARKodeMem)*adj_mem_ptr;
  ARKodeMRIStepMem adj_step_mem = (ARKodeMRIStepMem)adj_ark_mem->step_mem;
  int retval = MRIStepSetCoupling(*adj_mem_ptr, fwd_step_mem->MRIC);
  if (retval != ARK_SUCCESS) { return retval; }

  retval = ARKodeSetFixedStep(*adj_mem_ptr, -fwd_ark_mem->hin);
  if (retval != ARK_SUCCESS) { return retval; }

  retval = ARKodeSetMaxNumSteps(*adj_mem_ptr, nst);
  if (retval != ARK_SUCCESS) { return retval; }

  /* Although the adjoint stepper should not be using interpolation, if we use
   * the default Hermite interpolant, the adjoint stepper will forward a
   * full_RHS call to the inner stepper which uses a different NVector format,
   * even if the output time aligns with a step. Therefore, we use Lagrange
   * interpolation instead.
   */
  retval = ARKodeSetInterpolantType(*adj_mem_ptr, ARK_INTERP_LAGRANGE);
  if (retval != ARK_SUCCESS) { return retval; }

  retval = ARKodeSetAdjointCheckpointScheme(*adj_mem_ptr,
                                            fwd_ark_mem->checkpoint_scheme);
  if (retval != ARK_SUCCESS) { return retval; }

  retval = mriStepAdjointCreateData(adj_step_mem, adj_fse, inner_stepper,
                                    fwd_ark_mem->yn);
  if (retval != ARK_SUCCESS) { return retval; }

  adj_ark_mem->do_adjoint = SUNTRUE;
  return ARK_SUCCESS;
}

static int mriStepAdjointCreateSUNSteppers(
  ARKodeMem fwd_ark_mem, void* adj_mem, SUNAdjointStepper inner_stepper,
  MRIStepInnerAdjointProblem inner_problem, sunrealtype tf, N_Vector sf,
  suncountertype final_step_idx, SUNContext sunctx,
  SUNAdjointStepper* adj_stepper_ptr)
{
  SUNStepper fwd_sunstepper = NULL;
  SUNStepper adj_sunstepper = NULL;
  int retval = ARKodeCreateSUNStepper(fwd_ark_mem, &fwd_sunstepper);
  if (retval != ARK_SUCCESS)
  {
    ARKodeFree(&adj_mem);
    return retval;
  }

  SUNErrCode err = SUNStepper_SetReInitFn(fwd_sunstepper,
                                          mriStepAdjointSUNStepperReInit);
  if (err != SUN_SUCCESS)
  {
    SUNStepper_Destroy(&fwd_sunstepper);
    ARKodeFree(&adj_mem);
    return ARK_SUNSTEPPER_ERR;
  }

  retval = ARKodeCreateSUNStepper(adj_mem, &adj_sunstepper);
  if (retval != ARK_SUCCESS)
  {
    SUNStepper_Destroy(&fwd_sunstepper);
    SUNStepper_Destroy(&adj_sunstepper);
    ARKodeFree(&adj_mem);
    return retval;
  }

  err = SUNStepper_SetReInitFn(adj_sunstepper, mriStepAdjointSUNStepperReInit);
  if (err != SUN_SUCCESS)
  {
    SUNStepper_Destroy(&fwd_sunstepper);
    SUNStepper_Destroy(&adj_sunstepper);
    ARKodeFree(&adj_mem);
    return ARK_SUNSTEPPER_ERR;
  }

  err = SUNStepper_SetDestroyFn(adj_sunstepper, arkSUNStepperSelfDestruct);
  if (err != SUN_SUCCESS)
  {
    SUNStepper_Destroy(&fwd_sunstepper);
    SUNStepper_Destroy(&adj_sunstepper);
    ARKodeFree(&adj_mem);
    return ARK_SUNSTEPPER_ERR;
  }

  err = SUNAdjointStepper_Create(fwd_sunstepper, SUNTRUE, adj_sunstepper,
                                 SUNTRUE, final_step_idx, tf, sf,
                                 fwd_ark_mem->checkpoint_scheme, sunctx,
                                 adj_stepper_ptr);
  if (err != SUN_SUCCESS)
  {
    SUNStepper_Destroy(&fwd_sunstepper);
    SUNStepper_Destroy(&adj_sunstepper);
    return ARK_SUNADJSTEPPER_ERR;
  }

  err = SUNAdjointStepper_SetUserData(*adj_stepper_ptr, fwd_ark_mem->user_data);
  if (err == SUN_SUCCESS)
  {
    err = SUNAdjointStepper_SetUserData(inner_stepper, inner_problem);
  }
  if (err == SUN_SUCCESS)
  {
    retval = ARKodeSetUserData(adj_mem, *adj_stepper_ptr);
  }
  if (err != SUN_SUCCESS || retval != ARK_SUCCESS)
  {
    SUNAdjointStepper_Destroy(adj_stepper_ptr);
    return (err != SUN_SUCCESS) ? ARK_SUNADJSTEPPER_ERR : retval;
  }
  return ARK_SUCCESS;
}

int MRIStepCreateAdjointStepper(void* arkode_mem, SUNAdjointStepper inner_stepper,
                                SUNAdjRhsFn adj_fse, SUNAdjRhsFn adj_fsi,
                                MRIStepInnerAdjointProblem inner_problem,
                                sunrealtype tf, N_Vector sf, SUNContext sunctx,
                                SUNAdjointStepper* adj_stepper_ptr)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;

  if (!adj_stepper_ptr)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The output adjoint stepper pointer is NULL");
    return ARK_ILL_INPUT;
  }
  *adj_stepper_ptr = NULL;

  int retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                           &step_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The MRIStep memory pointer is NULL");
    return ARK_ILL_INPUT;
  }

  if (mriStepCompatibleWithAdjointSolver(ark_mem, step_mem, __LINE__, __func__,
                                         __FILE__))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "ark_mem provided is not compatible with adjoint calculation");
    return ARK_ILL_INPUT;
  }

  if (adj_fse == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "adj_fse cannot be NULL");
    return ARK_ILL_INPUT;
  }

  if (adj_fsi != NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "Implicit methods are not yet supported by the adjoint stepper.");
    return ARK_ILL_INPUT;
  }

  if (inner_stepper == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "inner_stepper cannot be NULL");
    return ARK_ILL_INPUT;
  }

  if (inner_problem == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "inner_problem cannot be NULL");
    return ARK_ILL_INPUT;
  }

  if (inner_problem->nforcing != step_mem->MRIC->nmat)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "inner_problem has an incompatible number of forcing vectors");
    return ARK_ILL_INPUT;
  }
  // TODO(SBR): fix this
  MRIStepAdjointState outer_state;
  if (mriStepAdjointGetOuterState(sf, &outer_state) != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "Incompatible vector type provided for adjoint calculation");
    return ARK_ILL_INPUT;
  }
  if (!inner_problem->sens_view ||
      N_VGetVectorID(inner_problem->sens_view) != SUNDIALS_NVEC_MANYVECTOR ||
      N_VGetNumSubvectors_ManyVector(inner_problem->sens_view) !=
        outer_state.count)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "The outer and inner terminal states must both include or omit mu");
    return ARK_ILL_INPUT;
  }

  MRIStepInnerAdjointState inner_state;
  if (mriStepAdjointGetInnerState(inner_problem->sf, outer_state.count,
                                  inner_problem->nforcing,
                                  &inner_state) != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The augmented inner terminal state is incompatible");
    return ARK_ILL_INPUT;
  }

  long int nst = 0;
  retval       = ARKodeGetNumSteps(arkode_mem, &nst);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "ARKodeGetNumSteps failed");
    return retval;
  }

  void* adj_mem = NULL;
  retval = mriStepAdjointCreateARKodeMem(ark_mem, step_mem, inner_stepper,
                                         adj_fse, tf, sf, nst, &adj_mem);
  if (retval != ARK_SUCCESS)
  {
    ARKodeFree(&adj_mem);
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to create the adjoint MRIStep memory");
    return retval;
  }

  retval = mriStepAdjointCreateSUNSteppers(ark_mem, adj_mem, inner_stepper,
                                           inner_problem, tf, sf,
                                           (suncountertype)nst - 1, sunctx,
                                           adj_stepper_ptr);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to create the MRIStep SUNAdjointStepper");
    return retval;
  }

  return ARK_SUCCESS;
}
