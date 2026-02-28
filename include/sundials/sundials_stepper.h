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

#ifndef _SUNDIALS_STEPPER_H
#define _SUNDIALS_STEPPER_H

#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus
extern "C" {
#endif

enum SUNFullRhsMode
{
  SUN_FULLRHS_START,
  SUN_FULLRHS_END,
  SUN_FULLRHS_OTHER
};

#ifndef SWIG
typedef enum SUNFullRhsMode SUNFullRhsMode;
#endif

typedef int (*SUNRhsJacFn)(sunrealtype t, N_Vector y, N_Vector fy,
                           SUNMatrix Jac, void* user_data, N_Vector tmp1,
                           N_Vector tmp2, N_Vector tmp3);

typedef int (*SUNRhsJacTimesFn)(N_Vector v, N_Vector Jv, sunrealtype t,
                                N_Vector y, N_Vector fy, void* user_data,
                                N_Vector tmp);

typedef _SUNDIALS_STRUCT_ SUNStepper_* SUNStepper;

typedef SUNErrCode (*SUNStepperEvolveFn)(SUNStepper stepper, sunrealtype tout,
                                         N_Vector vret, sunrealtype* tret);

typedef SUNErrCode (*SUNStepperOneStepFn)(SUNStepper stepper, sunrealtype tout,
                                          N_Vector vret, sunrealtype* tret);

typedef SUNErrCode (*SUNStepperFullRhsFn)(SUNStepper stepper, sunrealtype t,
                                          N_Vector v, N_Vector f,
                                          SUNFullRhsMode mode);

typedef SUNErrCode (*SUNStepperReInitFn)(SUNStepper stepper, sunrealtype t0,
                                         N_Vector v0);

typedef SUNErrCode (*SUNStepperResetFn)(SUNStepper stepper, sunrealtype tR,
                                        N_Vector vR);

typedef SUNErrCode (*SUNStepperResetCheckpointIndexFn)(SUNStepper stepper,
                                                       suncountertype ckptIdxR);

typedef SUNErrCode (*SUNStepperSetStopTimeFn)(SUNStepper stepper,
                                              sunrealtype tstop);

typedef SUNErrCode (*SUNStepperSetStepDirectionFn)(SUNStepper stepper,
                                                   sunrealtype stepdir);

typedef SUNErrCode (*SUNStepperSetForcingFn)(SUNStepper stepper,
                                             sunrealtype tshift,
                                             sunrealtype tscale,
                                             N_Vector* forcing_1d, int nforcing);

typedef SUNErrCode (*SUNStepperGetNumStepsFn)(SUNStepper stepper,
                                              suncountertype* nst);

typedef SUNErrCode (*SUNStepperDestroyFn)(SUNStepper stepper);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_Create(SUNContext sunctx, SUNStepper* stepper);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_Destroy(SUNStepper* stepper);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_Evolve(SUNStepper stepper, sunrealtype tout,
                             N_Vector vret, sunrealtype* tret);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_OneStep(SUNStepper stepper, sunrealtype tout,
                              N_Vector vret, sunrealtype* tret);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_FullRhs(SUNStepper stepper, sunrealtype t, N_Vector v,
                              N_Vector f, SUNFullRhsMode mode);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_ReInit(SUNStepper stepper, sunrealtype t0, N_Vector v0);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_Reset(SUNStepper stepper, sunrealtype tR, N_Vector vR);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_ResetCheckpointIndex(SUNStepper stepper,
                                           suncountertype ckptIdxR);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetStopTime(SUNStepper stepper, sunrealtype tstop);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetStepDirection(SUNStepper stepper, sunrealtype stepdir);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetForcing(SUNStepper stepper, sunrealtype tshift,
                                 sunrealtype tscale, N_Vector* forcing_1d,
                                 int nforcing);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetContent(SUNStepper stepper, void* content);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_GetContent(SUNStepper stepper, void** content);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetLastFlag(SUNStepper stepper, int last_flag);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_GetLastFlag(SUNStepper stepper, int* last_flag);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_GetNumSteps(SUNStepper stepper, suncountertype* nst);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetEvolveFn(SUNStepper stepper, SUNStepperEvolveFn fn);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetOneStepFn(SUNStepper stepper, SUNStepperOneStepFn fn);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetFullRhsFn(SUNStepper stepper, SUNStepperFullRhsFn fn);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetReInitFn(SUNStepper stepper, SUNStepperReInitFn fn);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetResetFn(SUNStepper stepper, SUNStepperResetFn fn);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetResetCheckpointIndexFn(
  SUNStepper stepper, SUNStepperResetCheckpointIndexFn fn);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetStopTimeFn(SUNStepper stepper,
                                    SUNStepperSetStopTimeFn fn);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetStepDirectionFn(SUNStepper stepper,
                                         SUNStepperSetStepDirectionFn fn);

SUNDIALS_EXPORT SUNErrCode SUNStepper_SetForcingFn(SUNStepper stepper,
                                                   SUNStepperSetForcingFn fn);

SUNDIALS_EXPORT
SUNErrCode SUNStepper_SetGetNumStepsFn(SUNStepper stepper,
                                       SUNStepperGetNumStepsFn fn);

SUNDIALS_EXPORT SUNErrCode SUNStepper_SetDestroyFn(SUNStepper stepper,
                                                   SUNStepperDestroyFn fn);

#ifdef __cplusplus
}
#endif

#endif /* _SUNDIALS_STEPPER_H */
