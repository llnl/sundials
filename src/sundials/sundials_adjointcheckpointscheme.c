/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * SUNAdjointCheckpointScheme class definition.
 * ----------------------------------------------------------------*/

#include <stdint.h>
#include <sundials/sundials_adjointcheckpointscheme.h>
#include <sundials/sundials_core.h>

#include "sundials/priv/sundials_errors_impl.h"
#include "sundials/sundials_errors.h"
#include "sundials/sundials_types.h"
#include "sundials_adjointcheckpointscheme_impl.h"

SUNErrCode SUNAdjointCheckpointScheme_NewEmpty(
  SUNContext sunctx, SUNAdjointCheckpointScheme* check_scheme_ptr)
{
  SUNFunctionBegin(sunctx);

  SUNAdjointCheckpointScheme self = NULL;
  self                            = malloc(sizeof(*self));
  SUNAssert(self, SUN_ERR_MALLOC_FAIL);

  self->sunctx  = sunctx;
  self->content = NULL;
  self->ops     = NULL;

  SUNAdjointCheckpointScheme_Ops ops = NULL;
  ops                                = malloc(sizeof(*ops));
  SUNAssert(ops, SUN_ERR_MALLOC_FAIL);

  ops->needssaving  = NULL;
  ops->insertvector = NULL;
  ops->loadvector   = NULL;
  ops->enableDense  = NULL;
  ops->destroy      = NULL;

  self->ops         = ops;
  *check_scheme_ptr = self;

  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_NeedsSaving(SUNAdjointCheckpointScheme self,
                                                  suncountertype step_num,
                                                  suncountertype stage_num,
                                                  sunrealtype t,
                                                  sunbooleantype* yes_or_no)
{
  SUNFunctionBegin(self->sunctx);
  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);

  if (self->ops->needssaving)
  {
    SUNErrCode err = self->ops->needssaving(self, step_num, stage_num, t,
                                            yes_or_no);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }
  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNAdjointCheckpointScheme_InsertVector(SUNAdjointCheckpointScheme self,
                                                   suncountertype step_num,
                                                   suncountertype stage_num,
                                                   sunrealtype t, N_Vector state)
{
  SUNFunctionBegin(self->sunctx);
  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);
  if (self->ops->insertvector)
  {
    SUNErrCode err = self->ops->insertvector(self, step_num, stage_num, t, state);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }
  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNAdjointCheckpointScheme_LoadVector(SUNAdjointCheckpointScheme self,
                                                 suncountertype step_num,
                                                 suncountertype stage_num,
                                                 sunbooleantype peek,
                                                 N_Vector* out, sunrealtype* tout)
{
  SUNFunctionBegin(self->sunctx);
  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);
  if (self->ops->loadvector)
  {
    SUNErrCode err = self->ops->loadvector(self, step_num, stage_num, peek, out,
                                           tout);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }
  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNAdjointCheckpointScheme_Destroy(
  SUNAdjointCheckpointScheme* check_scheme_ptr)
{
  SUNFunctionBegin((*check_scheme_ptr)->sunctx);
  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);
  if ((*check_scheme_ptr)->ops->destroy)
  {
    SUNErrCode err = (*check_scheme_ptr)->ops->destroy(check_scheme_ptr);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }
  else if (*check_scheme_ptr)
  {
    free((*check_scheme_ptr)->ops);
    free(*check_scheme_ptr);
  }
  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_EnableDense(SUNAdjointCheckpointScheme self,
                                                  sunbooleantype on_or_off)
{
  SUNFunctionBegin(self->sunctx);
  SUNDIALS_MARK_FUNCTION_BEGIN(SUNCTX_->profiler);
  if (self->ops->enableDense)
  {
    SUNErrCode err = self->ops->enableDense(self, on_or_off);
    SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
    return err;
  }
  SUNDIALS_MARK_FUNCTION_END(SUNCTX_->profiler);
  return SUN_ERR_NOT_IMPLEMENTED;
}

SUNErrCode SUNAdjointCheckpointScheme_SetContent(SUNAdjointCheckpointScheme self,
                                                 void* content)
{
  SUNFunctionBegin(self->sunctx);
  self->content = content;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_GetContent(SUNAdjointCheckpointScheme self,
                                                 void** content)
{
  SUNFunctionBegin(self->sunctx);
  *content = self->content;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_SetNeedsSavingFn(
  SUNAdjointCheckpointScheme self, SUNAdjointCheckpointSchemeNeedsSavingFn fn)
{
  SUNFunctionBegin(self->sunctx);
  self->ops->needssaving = fn;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_SetInsertVectorFn(
  SUNAdjointCheckpointScheme self, SUNAdjointCheckpointSchemeInsertVectorFn fn)
{
  SUNFunctionBegin(self->sunctx);
  self->ops->insertvector = fn;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_SetLoadVectorFn(
  SUNAdjointCheckpointScheme self, SUNAdjointCheckpointSchemeLoadVectorFn fn)
{
  SUNFunctionBegin(self->sunctx);
  self->ops->loadvector = fn;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_SetDestroyFn(
  SUNAdjointCheckpointScheme self, SUNAdjointCheckpointSchemeDestroyFn fn)
{
  SUNFunctionBegin(self->sunctx);
  self->ops->destroy = fn;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdjointCheckpointScheme_SetEnableDenseFn(
  SUNAdjointCheckpointScheme self, SUNAdjointCheckpointSchemeEnableDenseFn fn)
{
  SUNFunctionBegin(self->sunctx);
  self->ops->enableDense = fn;
  return SUN_SUCCESS;
}
