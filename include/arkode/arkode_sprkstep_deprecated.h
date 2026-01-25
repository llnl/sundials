/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
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
 * ---------------------------------------------------------------------------*/

#ifndef _SPRKSTEP_DEPRECATED_H
#define _SPRKSTEP_DEPRECATED_H

#include <arkode/arkode.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* --------------------------------------------------------------------------
 * Deprecated Functions -- all are superseded by shared ARKODE-level routines
 * -------------------------------------------------------------------------- */

SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeReset instead")
int SPRKStepReset(void* arkode_mem, sunrealtype tR, N_Vector yR);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeRootInit instead")
int SPRKStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRootDirection instead")
int SPRKStepSetRootDirection(void* arkode_mem, int* rootdir);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNoInactiveRootWarn instead")
int SPRKStepSetNoInactiveRootWarn(void* arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetDefaults instead")
int SPRKStepSetDefaults(void* arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetOrder instead")
int SPRKStepSetOrder(void* arkode_mem, int maxord);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInterpolantType instead")
int SPRKStepSetInterpolantType(void* arkode_mem, int itype);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInterpolantDegree instead")
int SPRKStepSetInterpolantDegree(void* arkode_mem, int degree);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxNumSteps instead")
int SPRKStepSetMaxNumSteps(void* arkode_mem, long int mxsteps);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetStopTime instead")
int SPRKStepSetStopTime(void* arkode_mem, sunrealtype tstop);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetFixedStep instead")
int SPRKStepSetFixedStep(void* arkode_mem, sunrealtype hfixed);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetUserData instead")
int SPRKStepSetUserData(void* arkode_mem, void* user_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetPostprocessStepFn instead")
int SPRKStepSetPostprocessStepFn(void* arkode_mem, ARKPostProcessFn ProcessStep);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetPostprocessStageFn instead")
int SPRKStepSetPostprocessStageFn(void* arkode_mem,
                                  ARKPostProcessFn ProcessStage);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeEvolve instead")
int SPRKStepEvolve(void* arkode_mem, sunrealtype tout, N_Vector yout,
                   sunrealtype* tret, int itask);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetDky instead")
int SPRKStepGetDky(void* arkode_mem, sunrealtype t, int k, N_Vector dky);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetReturnFlagName instead")
char* SPRKStepGetReturnFlagName(long int flag);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetCurrentState instead")
int SPRKStepGetCurrentState(void* arkode_mem, N_Vector* state);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetCurrentStep instead")
int SPRKStepGetCurrentStep(void* arkode_mem, sunrealtype* hcur);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetCurrentTime instead")
int SPRKStepGetCurrentTime(void* arkode_mem, sunrealtype* tcur);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetLastStep instead")
int SPRKStepGetLastStep(void* arkode_mem, sunrealtype* hlast);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumStepAttempts instead")
int SPRKStepGetNumStepAttempts(void* arkode_mem, long int* step_attempts);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumSteps instead")
int SPRKStepGetNumSteps(void* arkode_mem, long int* nsteps);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetRootInfo instead")
int SPRKStepGetRootInfo(void* arkode_mem, int* rootsfound);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetUserData instead")
int SPRKStepGetUserData(void* arkode_mem, void** user_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodePrintAllStats instead")
int SPRKStepPrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeWriteParameters instead")
int SPRKStepWriteParameters(void* arkode_mem, FILE* fp);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetStepStats instead")
int SPRKStepGetStepStats(void* arkode_mem, long int* nsteps,
                         sunrealtype* hinused, sunrealtype* hlast,
                         sunrealtype* hcur, sunrealtype* tcur);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeFree instead")
void SPRKStepFree(void** arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumRhsEvals instead")
int SPRKStepGetNumRhsEvals(void* arkode_mem, long int* nf1, long int* nf2);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetUseCompensatedSums instead")
int SPRKStepSetUseCompensatedSums(void* arkode_mem, sunbooleantype onoff);

#ifdef __cplusplus
}
#endif

#endif
