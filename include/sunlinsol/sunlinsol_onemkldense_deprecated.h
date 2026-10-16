/* ---------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * -----------------------------------------------------------------
 * Temporary header file for deprecated functions.
 * TODO:(SBR) remove in version 8.1.0
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_SUNLINSOL_SUNLINSOL_ONEMKLDENSE_DEPRECATED_H
#define SUNDIALS_SUNLINSOL_SUNLINSOL_ONEMKLDENSE_DEPRECATED_H

#include <sunlinsol/sunlinsol_onemkldense.h>

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_DEPRECATED_MSG(
  "use SUNLinSolGetType instead; will be removed in version 8.1.0")

static inline SUNLinearSolver_Type SUNLinSolGetType_OneMklDense(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_DIRECT;
}

SUNDIALS_DEPRECATED_MSG(
  "use SUNLinSolGetID instead; will be removed in version 8.1.0")

static inline SUNLinearSolver_ID SUNLinSolGetID_OneMklDense(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_ONEMKLDENSE;
}

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolInitialize instead; will be removed in version 8.1.0")
SUNErrCode SUNLinSolInitialize_OneMklDense(SUNLinearSolver S);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolSetup instead; will be removed in version 8.1.0")
int SUNLinSolSetup_OneMklDense(SUNLinearSolver S, SUNMatrix A);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolSolve instead; will be removed in version 8.1.0")
int SUNLinSolSolve_OneMklDense(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                               N_Vector b, sunrealtype tol);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolLastFlag instead; will be removed in version 8.1.0")
sunindextype SUNLinSolLastFlag_OneMklDense(SUNLinearSolver S);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNLinSolFree instead; will be removed in version 8.1.0")
SUNErrCode SUNLinSolFree_OneMklDense(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

#endif
