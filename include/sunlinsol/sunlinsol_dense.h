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
 * This is the header file for the dense implementation of the
 * SUNLINSOL module, SUNLINSOL_DENSE.
 *
 * Notes:
 *   - The definition of the generic SUNLinearSolver structure can
 *     be found in the header file sundials_linearsolver.h.
 *   - The definition of the type 'sunrealtype' can be found in the
 *     header file sundials_types.h, and it may be changed (at the
 *     configuration stage) according to the user's needs.
 *     The sundials_types.h file also contains the definition
 *     for the type 'sunbooleantype' and 'indextype'.
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_DENSE_H
#define _SUNLINSOL_DENSE_H

#include <sundials/sundials_dense.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sunmatrix/sunmatrix_dense.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ----------------------------------------
 * Dense Implementation of SUNLinearSolver
 * ---------------------------------------- */

struct _SUNLinearSolverContent_Dense
{
  sunindextype N;
  sunindextype* pivots;
  sunindextype last_flag;
};

typedef struct _SUNLinearSolverContent_Dense* SUNLinearSolverContent_Dense;

/* ----------------------------------------
 * Exported Functions for SUNLINSOL_DENSE
 * ---------------------------------------- */

SUNDIALS_EXPORT
SUNLinearSolver SUNLinSol_Dense(N_Vector y, SUNMatrix A, SUNContext sunctx);

SUNDIALS_EXPORT
SUNLinearSolver_Type SUNLinSolGetType_Dense(SUNLinearSolver S);

SUNDIALS_EXPORT
SUNLinearSolver_ID SUNLinSolGetID_Dense(SUNLinearSolver S);

SUNDIALS_EXPORT
SUNErrCode SUNLinSolInitialize_Dense(SUNLinearSolver S);

SUNDIALS_EXPORT
int SUNLinSolSetup_Dense(SUNLinearSolver S, SUNMatrix A);

SUNDIALS_EXPORT
int SUNLinSolSolve_Dense(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b,
                         sunrealtype tol);

SUNDIALS_EXPORT
sunindextype SUNLinSolLastFlag_Dense(SUNLinearSolver S);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
SUNErrCode SUNLinSolSpace_Dense(SUNLinearSolver S, long int* lenrwLS,
                                long int* leniwLS);

SUNDIALS_EXPORT
SUNErrCode SUNLinSolFree_Dense(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

#endif
