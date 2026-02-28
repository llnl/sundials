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
 * -----------------------------------------------------------------
 * SUNAdjointCheckpointScheme_Fixed class declaration.
 * ----------------------------------------------------------------*/

#ifndef _SUNADJOINTCHECKPOINTSCHEME_FIXED_H
#define _SUNADJOINTCHECKPOINTSCHEME_FIXED_H

#include <sundials/sundials_adjointcheckpointscheme.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_export.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Create_Fixed(
  SUNDataIOMode io_mode, SUNMemoryHelper mem_helper, suncountertype interval,
  suncountertype estimate, sunbooleantype keep, SUNContext sunctx,
  SUNAdjointCheckpointScheme* check_scheme_ptr);

SUNDIALS_EXPORT SUNErrCode SUNAdjointCheckpointScheme_NeedsSaving_Fixed(
  SUNAdjointCheckpointScheme check_scheme, suncountertype step_num,
  suncountertype stage_num, sunrealtype t, sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_InsertVector_Fixed(
  SUNAdjointCheckpointScheme check_scheme, suncountertype step_num,
  suncountertype stage_num, sunrealtype t, N_Vector state);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_LoadVector_Fixed(
  SUNAdjointCheckpointScheme check_scheme, suncountertype step_num,
  suncountertype stage_num, sunbooleantype peek, N_Vector* out,
  sunrealtype* tout);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Destroy_Fixed(
  SUNAdjointCheckpointScheme* check_scheme_ptr);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_EnableDense_Fixed(
  SUNAdjointCheckpointScheme check_scheme, sunbooleantype on_or_off);

#ifdef __cplusplus
}
#endif

#endif /* _SUNADJOINTCHECKPOINTSCHEME_FIXED_H */
