/* -----------------------------------------------------------------
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
 * SUNAdjointCheckpointScheme class declaration.
 * ----------------------------------------------------------------*/

#ifndef _SUNADJOINT_CHECKPOINTSCHEME_H
#define _SUNADJOINT_CHECKPOINTSCHEME_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef _SUNDIALS_STRUCT_ SUNAdjointCheckpointScheme_* SUNAdjointCheckpointScheme;

typedef SUNErrCode (*SUNAdjointCheckpointSchemeNeedsSavingFn)(
  SUNAdjointCheckpointScheme check_scheme, suncountertype step_num,
  suncountertype stage_num, sunrealtype t, sunbooleantype* yes_or_no);

typedef SUNErrCode (*SUNAdjointCheckpointSchemeInsertVectorFn)(
  SUNAdjointCheckpointScheme check_scheme, suncountertype step_num,
  suncountertype stage_num, sunrealtype t, N_Vector y);

typedef SUNErrCode (*SUNAdjointCheckpointSchemeLoadVectorFn)(
  SUNAdjointCheckpointScheme check_scheme, suncountertype step_num,
  suncountertype stage_num, sunbooleantype peek, N_Vector* yout,
  sunrealtype* tout);

typedef SUNErrCode (*SUNAdjointCheckpointSchemeDestroyFn)(
  SUNAdjointCheckpointScheme* check_scheme);

typedef SUNErrCode (*SUNAdjointCheckpointSchemeEnableDenseFn)(
  SUNAdjointCheckpointScheme check_scheme, sunbooleantype on_or_off);

/*
 * "static" base class methods
 */

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_NewEmpty(SUNContext sunctx,
                                               SUNAdjointCheckpointScheme*);

/*
 * Base class methods
 */

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_SetNeedsSavingFn(
  SUNAdjointCheckpointScheme check_scheme,
  SUNAdjointCheckpointSchemeNeedsSavingFn);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_SetInsertVectorFn(
  SUNAdjointCheckpointScheme check_scheme,
  SUNAdjointCheckpointSchemeInsertVectorFn);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_SetLoadVectorFn(
  SUNAdjointCheckpointScheme check_scheme,
  SUNAdjointCheckpointSchemeLoadVectorFn);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_SetDestroyFn(
  SUNAdjointCheckpointScheme check_scheme, SUNAdjointCheckpointSchemeDestroyFn);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_SetEnableDenseFn(
  SUNAdjointCheckpointScheme check_scheme,
  SUNAdjointCheckpointSchemeEnableDenseFn);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_SetContent(
  SUNAdjointCheckpointScheme check_scheme, void* content);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_GetContent(
  SUNAdjointCheckpointScheme check_scheme, void** content);

/*
 * Virtual (overridable) base class methods
 */

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_NeedsSaving(
  SUNAdjointCheckpointScheme check_scheme, suncountertype step_num,
  suncountertype stage_num, sunrealtype t, sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_InsertVector(
  SUNAdjointCheckpointScheme check_scheme, suncountertype step_num,
  suncountertype stage_num, sunrealtype t, N_Vector state);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_LoadVector(
  SUNAdjointCheckpointScheme check_scheme, suncountertype step_num,
  suncountertype stage_num, sunbooleantype peek, N_Vector* out,
  sunrealtype* tout); // nb::rv_policy::reference

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Destroy(SUNAdjointCheckpointScheme*);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_EnableDense(
  SUNAdjointCheckpointScheme check_scheme, sunbooleantype on_or_off);

#ifdef __cplusplus
}
#endif

#endif /*_SUNADJOINT_CHECKPOINTSCHEME_H*/
