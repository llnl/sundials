/*
 * -----------------------------------------------------------------
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
 * This is the header file for the MAGMA dense implementation of the
 * SUNLINSOL module, SUNLINSOL_MAGMADENSE.
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_MAGMADENSE_H
#define _SUNLINSOL_MAGMADENSE_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_memory.h>
#include <sundials/sundials_nvector.h>

#if defined(SUNDIALS_MAGMA_BACKENDS_CUDA)
#define HAVE_CUBLAS
#elif defined(SUNDIALS_MAGMA_BACKENDS_HIP)
#define HAVE_HIP
#endif
#include <magma_v2.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------
 * MAGMA dense implementation of SUNLinearSolver
 * ----------------------------------------------- */

struct _SUNLinearSolverContent_MagmaDense
{
  int last_flag;
  sunbooleantype async;
  sunindextype N;
  SUNMemory pivots;
  SUNMemory pivotsarr;
  SUNMemory dpivotsarr;
  SUNMemory infoarr;
  SUNMemory rhsarr;
  SUNMemoryHelper memhelp;
  magma_queue_t q;
};

typedef struct _SUNLinearSolverContent_MagmaDense* SUNLinearSolverContent_MagmaDense;

SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_MagmaDense(N_Vector y, SUNMatrix A,
                                                     SUNContext sunctx);

SUNDIALS_EXPORT int SUNLinSol_MagmaDense_SetAsync(SUNLinearSolver S,
                                                  sunbooleantype onoff);

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_MagmaDense(SUNLinearSolver S);
SUNDIALS_EXPORT SUNLinearSolver_ID SUNLinSolGetID_MagmaDense(SUNLinearSolver S);
SUNDIALS_EXPORT SUNErrCode SUNLinSolInitialize_MagmaDense(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetup_MagmaDense(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_MagmaDense(SUNLinearSolver S, SUNMatrix A,
                                              N_Vector x, N_Vector b,
                                              sunrealtype tol);
SUNDIALS_EXPORT sunindextype SUNLinSolLastFlag_MagmaDense(SUNLinearSolver S);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
SUNErrCode SUNLinSolSpace_MagmaDense(SUNLinearSolver S, long int* lenrwLS,
                                     long int* leniwLS);
SUNDIALS_EXPORT SUNErrCode SUNLinSolFree_MagmaDense(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

#endif
