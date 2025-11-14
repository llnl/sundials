/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
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
 * SUNAdjointCheckpointScheme_Fixed class definition.
 * ----------------------------------------------------------------*/

#include <sunadjointcheckpointscheme/sunadjointcheckpointscheme_fixed.h>
#include <sundials/sundials_adjointcheckpointscheme.h>
#include <sundials/sundials_core.h>

#include "sundatanode/sundatanode_inmem.h"
#include "sundials_adjointcheckpointscheme_impl.h"
#include "sundials_datanode.h"
#include "sundials_logger_impl.h"
#include "sundials_macros.h"
#include "sundials_utils.h"

struct SUNAdjointCheckpointScheme_Fixed_Content_
{
  suncountertype backup_interval;
  suncountertype interval;
  suncountertype step_num_of_current_insert;
  suncountertype step_num_of_current_load;
  SUNMemoryHelper mem_helper;
  SUNDataNode root_node;
  SUNDataNode current_insert_step_node;
  SUNDataNode current_load_step_node;
  SUNDataIOMode io_mode;
  sunbooleantype keep;
};

typedef struct SUNAdjointCheckpointScheme_Fixed_Content_*
  SUNAdjointCheckpointScheme_Fixed_Content;

#define GET_CONTENT(S)       ((SUNAdjointCheckpointScheme_Fixed_Content)S->content)
#define IMPL_MEMBER(S, prop) (GET_CONTENT(S)->prop)

SUNErrCode SUNAdjointCheckpointScheme_Create_Fixed(
  SUNDataIOMode io_mode, SUNMemoryHelper mem_helper, suncountertype interval,
  suncountertype estimate, sunbooleantype keep, SUNContext sunctx,
  SUNAdjointCheckpointScheme* check_scheme_ptr)
{
  SUNFunctionBegin(sunctx);

  SUNAdjointCheckpointScheme check_scheme = NULL;
  SUNCheckCall(SUNAdjointCheckpointScheme_NewEmpty(sunctx, &check_scheme));

  check_scheme->ops->needssaving = SUNAdjointCheckpointScheme_NeedsSaving_Fixed;
  check_scheme->ops->insertvector = SUNAdjointCheckpointScheme_InsertVector_Fixed;
  check_scheme->ops->loadvector  = SUNAdjointCheckpointScheme_LoadVector_Fixed;
  check_scheme->ops->enableDense = SUNAdjointCheckpointScheme_EnableDense_Fixed;
  check_scheme->ops->destroy     = SUNAdjointCheckpointScheme_Destroy_Fixed;

  SUNAdjointCheckpointScheme_Fixed_Content content = NULL;

  content = malloc(sizeof(*content));
  SUNAssert(content, SUN_ERR_MALLOC_FAIL);

  content->mem_helper                 = mem_helper;
  content->interval                   = interval;
  content->keep                       = keep;
  content->root_node                  = NULL;
  content->current_insert_step_node   = NULL;
  content->step_num_of_current_insert = -2;
  content->current_load_step_node     = NULL;
  content->step_num_of_current_load   = -2;
  content->io_mode                    = io_mode;

  SUNCheckCall(
    SUNDataNode_CreateObject(io_mode, estimate, sunctx, &content->root_node));

  check_scheme->content = content;
  *check_scheme_ptr     = check_scheme;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_NeedsSaving_Fixed(
  SUNAdjointCheckpointScheme self, suncountertype step_num,
  SUNDIALS_MAYBE_UNUSED suncountertype stage_num,
  SUNDIALS_MAYBE_UNUSED sunrealtype t, sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);

  if (!(step_num % IMPL_MEMBER(self, interval))) { *yes_or_no = SUNTRUE; }
  else { *yes_or_no = SUNFALSE; }

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_InsertVector_Fixed(
  SUNAdjointCheckpointScheme self, suncountertype step_num,
  SUNDIALS_MAYBE_UNUSED suncountertype stage_num, sunrealtype t, N_Vector y)
{
  SUNFunctionBegin(self->sunctx);

  /* If this is the first state for a step, then we need to create a
     list node first to store the step and all stage solutions in.
     We keep a pointer to the list node until this step is over for
     fast access when inserting stages. */
  SUNDataNode step_data_node = NULL;
  if (step_num != IMPL_MEMBER(self, step_num_of_current_insert))
  {
    SUNCheckCall(SUNDataNode_CreateList(IMPL_MEMBER(self, io_mode), 0, SUNCTX_,
                                        &step_data_node));
    IMPL_MEMBER(self, current_insert_step_node)   = step_data_node;
    IMPL_MEMBER(self, step_num_of_current_insert) = step_num;

    /* Store the step node in the root node object. */
    char* key = sunSignedToString(step_num);
    SUNLogExtraDebug(SUNCTX_->logger, "insert-new-step",
                     "step_num = %d, key = %s", step_num, key);
    SUNCheckCall(SUNDataNode_AddNamedChild(IMPL_MEMBER(self, root_node), key,
                                           step_data_node));
    free(key);
  }
  else { step_data_node = IMPL_MEMBER(self, current_insert_step_node); }

  /* Add the state data as a leaf node in the step node's list of children. */
  SUNDataNode solution_node = NULL;
  SUNCheckCall(SUNDataNode_CreateLeaf(IMPL_MEMBER(self, io_mode),
                                      IMPL_MEMBER(self, mem_helper), SUNCTX_,
                                      &solution_node));
  SUNCheckCall(SUNDataNode_SetDataNvector(solution_node, y, t));

  SUNLogExtraDebug(SUNCTX_->logger, "insert-stage",
                   "step_num = %d, stage_num = %d, t = " SUN_FORMAT_G, step_num,
                   stage_num, t);
  SUNCheckCall(SUNDataNode_AddChild(step_data_node, solution_node));

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_LoadVector_Fixed(
  SUNAdjointCheckpointScheme self, suncountertype step_num,
  suncountertype stage_num, sunbooleantype peek, N_Vector* yout,
  sunrealtype* tout)
{
  SUNFunctionBegin(self->sunctx);

  SUNErrCode errcode = SUN_SUCCESS;

  /* If we are trying to load the step solution, we need to load the list which holds
     the step and stage solutions. We keep a pointer to the list node until
     this step is over for fast access when loading stages. */
  SUNDataNode step_data_node = NULL;
  if (step_num != IMPL_MEMBER(self, step_num_of_current_load))
  {
    char* key = sunSignedToString(step_num);
    SUNLogExtraDebug(SUNCTX_->logger, "try-load-new-step",
                     "step_num = %d, stage_num = %d", step_num, stage_num);
    errcode = SUNDataNode_GetNamedChild(IMPL_MEMBER(self, root_node), key,
                                        &step_data_node);
    free(key);
    if (errcode == SUN_SUCCESS)
    {
      IMPL_MEMBER(self, current_load_step_node)   = step_data_node;
      IMPL_MEMBER(self, step_num_of_current_load) = step_num;
    }
    else if (errcode == SUN_ERR_DATANODE_NODENOTFOUND)
    {
      step_data_node = NULL;
    }
    else { SUNCheckCall(errcode); }
  }
  else { step_data_node = IMPL_MEMBER(self, current_load_step_node); }

  if (!step_data_node)
  {
    SUNLogExtraDebug(SUNCTX_->logger, "step-not-found",
                     "step_num = %d, stage_num = %d", step_num, stage_num);
    return SUN_ERR_CHECKPOINT_NOT_FOUND;
  }

  SUNLogExtraDebug(SUNCTX_->logger, "step-loaded",
                   "step_num = %d, stage_num = %d", step_num, stage_num);

  SUNDataNode solution_node = NULL;
  if (IMPL_MEMBER(self, keep) || peek)
  {
    SUNLogExtraDebug(SUNCTX_->logger, "try-load-stage",
                     "keep = 1, step_num = %d, stage_num = %d", step_num,
                     stage_num);
    errcode = SUNDataNode_GetChild(step_data_node, stage_num, &solution_node);
    if (errcode == SUN_ERR_DATANODE_NODENOTFOUND) { solution_node = NULL; }
    else { SUNCheckCall(errcode); }
  }
  else
  {
    sunbooleantype has_children = SUNFALSE;
    SUNCheckCall(SUNDataNode_HasChildren(step_data_node, &has_children));

    if (has_children)
    {
      SUNLogExtraDebug(SUNCTX_->logger, "try-load-stage",
                       "keep = 0, step_num = %d, stage_num = %d", step_num,
                       stage_num);
      errcode = SUNDataNode_RemoveChild(step_data_node, stage_num,
                                        &solution_node);
      if (errcode == SUN_ERR_DATANODE_NODENOTFOUND) { solution_node = NULL; }
      else { SUNCheckCall(errcode); }
    }

    /* If we just removed the last stage (so has_children==false),
       then we should remove the step too. */
    SUNCheckCall(SUNDataNode_HasChildren(step_data_node, &has_children));
    if (!has_children)
    {
      char* key = sunSignedToString(step_num);
      SUNLogExtraDebug(SUNCTX_->logger, "remove-step", "step_num = %d", step_num);
      SUNCheckCall(SUNDataNode_RemoveNamedChild(IMPL_MEMBER(self, root_node),
                                                key, &step_data_node));
      free(key);
      SUNCheckCall(SUNDataNode_Destroy(&step_data_node));
    }
  }

  if (!solution_node)
  {
    SUNLogExtraDebug(SUNCTX_->logger, "stage-not-found",
                     "step_num = %d, stage_num = %d", step_num, stage_num);
    return SUN_ERR_CHECKPOINT_NOT_FOUND;
  }

  SUNCheckCall(SUNDataNode_GetDataNvector(solution_node, *yout, tout));
  SUNLogExtraDebug(SUNCTX_->logger, "stage-loaded",
                   "step_num = %d, stage_num = %d, t = " SUN_FORMAT_G, step_num,
                   stage_num, *tout);

  /* Cleanup the checkpoint memory if need be */
  if (!(IMPL_MEMBER(self, keep) || peek))
  {
    SUNCheckCall(SUNDataNode_Destroy(&solution_node));
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_Destroy_Fixed(
  SUNAdjointCheckpointScheme* self_ptr)
{
  SUNFunctionBegin((*self_ptr)->sunctx);

  SUNAdjointCheckpointScheme self = *self_ptr;

  SUNCheckCall(SUNDataNode_Destroy(&IMPL_MEMBER(self, root_node)));

  free(self->content);
  free(self->ops);
  free(self);

  *self_ptr = NULL;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_EnableDense_Fixed(
  SUNAdjointCheckpointScheme check_scheme, sunbooleantype on_or_off)
{
  SUNFunctionBegin(check_scheme->sunctx);

  if (on_or_off)
  {
    IMPL_MEMBER(check_scheme, backup_interval) = IMPL_MEMBER(check_scheme,
                                                             interval);
    IMPL_MEMBER(check_scheme, interval)        = 1;
  }
  else
  {
    IMPL_MEMBER(check_scheme, interval) = IMPL_MEMBER(check_scheme,
                                                      backup_interval);
  }

  return SUN_SUCCESS;
}
