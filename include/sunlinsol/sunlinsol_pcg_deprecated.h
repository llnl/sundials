/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds, Ashley Crawford @ UMBC
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
 * Temporary header file for deprecated functions.
 * TODO:(SBR) remove in version 8.1.0
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_SUNLINSOL_SUNLINSOL_PCG_DEPRECATED_H
#define SUNDIALS_SUNLINSOL_SUNLINSOL_PCG_DEPRECATED_H

#include <sunlinsol/sunlinsol_pcg.h>

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolGetType instead; will be removed in version 8.1.0")
SUNLinearSolver_Type SUNLinSolGetType_PCG(SUNLinearSolver S);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolGetID instead; will be removed in version 8.1.0")
SUNLinearSolver_ID SUNLinSolGetID_PCG(SUNLinearSolver S);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolInitialize instead; will be removed in version 8.1.0")
SUNErrCode SUNLinSolInitialize_PCG(SUNLinearSolver S);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolSetATimes instead; will be removed in version 8.1.0")
SUNErrCode SUNLinSolSetATimes_PCG(SUNLinearSolver S, void* A_data,
                                  SUNATimesFn ATimes);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolSetPreconditioner instead; will be removed in version 8.1.0")
SUNErrCode SUNLinSolSetPreconditioner_PCG(SUNLinearSolver S, void* P_data,
                                          SUNPSetupFn Pset, SUNPSolveFn Psol);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolSetScalingVectors instead; will be removed in version 8.1.0")
SUNErrCode SUNLinSolSetScalingVectors_PCG(SUNLinearSolver S, N_Vector s,
                                          N_Vector nul);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolSetZeroGuess instead; will be removed in version 8.1.0")
SUNErrCode SUNLinSolSetZeroGuess_PCG(SUNLinearSolver S, sunbooleantype onoff);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolSetup instead; will be removed in version 8.1.0")
int SUNLinSolSetup_PCG(SUNLinearSolver S, SUNMatrix nul);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolSolve instead; will be removed in version 8.1.0")
int SUNLinSolSolve_PCG(SUNLinearSolver S, SUNMatrix nul, N_Vector x, N_Vector b,
                       sunrealtype tol);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolNumIters instead; will be removed in version 8.1.0")
int SUNLinSolNumIters_PCG(SUNLinearSolver S);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolResNorm instead; will be removed in version 8.1.0")
sunrealtype SUNLinSolResNorm_PCG(SUNLinearSolver S);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolResid instead; will be removed in version 8.1.0")
N_Vector SUNLinSolResid_PCG(SUNLinearSolver S);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolLastFlag instead; will be removed in version 8.1.0")
sunindextype SUNLinSolLastFlag_PCG(SUNLinearSolver S);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolFree instead; will be removed in version 8.1.0")
SUNErrCode SUNLinSolFree_PCG(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

#endif
