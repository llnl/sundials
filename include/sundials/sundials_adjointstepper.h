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
 * SUNAdjointStepper class definition.
 * ----------------------------------------------------------------*/

#ifndef _SUNADJOINT_STEPPER_H
#define _SUNADJOINT_STEPPER_H

#include <sundials/sundials_adjointcheckpointscheme.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_stepper.h>
#include "sundials/sundials_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef int (*SUNAdjRhsFn)(sunrealtype t, N_Vector y, N_Vector sens,
                           N_Vector sens_dot, void* user_data);

typedef struct SUNAdjointStepper_* SUNAdjointStepper;

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_Create(
  SUNStepper fwd_sunstepper, sunbooleantype own_fwd, SUNStepper adj_sunstepper,
  sunbooleantype own_adj, suncountertype final_step_idx, sunrealtype tf,
  N_Vector sf, SUNAdjointCheckpointScheme checkpoint_scheme, SUNContext sunctx,
  SUNAdjointStepper* adj_stepper);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_ReInit(SUNAdjointStepper adj, sunrealtype t0,
                                    N_Vector y0, sunrealtype tf, N_Vector sf);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_Evolve(SUNAdjointStepper adj_stepper,
                                    sunrealtype tout, N_Vector sens,
                                    sunrealtype* tret);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_OneStep(SUNAdjointStepper adj_stepper,
                                     sunrealtype tout, N_Vector sens,
                                     sunrealtype* tret);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_RecomputeFwd(SUNAdjointStepper adj_stepper,
                                          suncountertype start_idx,
                                          sunrealtype t0, N_Vector y0,
                                          sunrealtype tf);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_SetUserData(SUNAdjointStepper, void* user_data);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_GetNumSteps(SUNAdjointStepper adj_stepper,
                                         suncountertype* num_steps);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_GetNumRecompute(SUNAdjointStepper adj_stepper,
                                             suncountertype* num_recompute);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_PrintAllStats(SUNAdjointStepper adj_stepper,
                                           FILE* outfile, SUNOutputFormat fmt);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointStepper_Destroy(SUNAdjointStepper*);

#ifdef __cplusplus
}
#endif
#endif /* _SUNADJOINT_STEPPER_H */
