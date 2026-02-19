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
 * This is the header file for the PCG implementation of the
 * SUNLINSOL module, SUNLINSOL_PCG.  The PCG algorithm is based
 * on the Preconditioned Conjugate Gradient.
 *
 * Note:
 *   - The definition of the generic SUNLinearSolver structure can
 *     be found in the header file sundials_linearsolver.h.
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_PCG_H
#define _SUNLINSOL_PCG_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default PCG solver parameters */
#define SUNPCG_MAXL_DEFAULT 5

/* --------------------------------------
 * PCG Implementation of SUNLinearSolver
 * -------------------------------------- */

struct _SUNLinearSolverContent_PCG
{
  int maxl;
  int pretype;
  sunbooleantype zeroguess;
  int numiters;
  sunrealtype resnorm;
  int last_flag;

  SUNATimesFn ATimes;
  void* ATData;
  SUNPSetupFn Psetup;
  SUNPSolveFn Psolve;
  void* PData;

  N_Vector s;
  N_Vector r;
  N_Vector p;
  N_Vector z;
  N_Vector Ap;
};

typedef struct _SUNLinearSolverContent_PCG* SUNLinearSolverContent_PCG;

/* -------------------------------------
 * Exported Functions for SUNLINSOL_PCG
 * ------------------------------------- */

SUNDIALS_EXPORT
SUNLinearSolver SUNLinSol_PCG(N_Vector y, int pretype, int maxl,
                              SUNContext sunctx);

SUNDIALS_EXPORT
SUNErrCode SUNLinSol_PCGSetPrecType(SUNLinearSolver S, int pretype);

SUNDIALS_EXPORT
SUNErrCode SUNLinSol_PCGSetMaxl(SUNLinearSolver S, int maxl);

SUNDIALS_EXPORT
SUNLinearSolver_Type SUNLinSolGetType_PCG(SUNLinearSolver S);

SUNDIALS_EXPORT
SUNLinearSolver_ID SUNLinSolGetID_PCG(SUNLinearSolver S);

SUNDIALS_EXPORT
SUNErrCode SUNLinSolInitialize_PCG(SUNLinearSolver S);

SUNDIALS_EXPORT
SUNErrCode SUNLinSolSetATimes_PCG(SUNLinearSolver S, void* A_data,
                                  SUNATimesFn ATimes);

SUNDIALS_EXPORT
SUNErrCode SUNLinSolSetPreconditioner_PCG(SUNLinearSolver S, void* P_data,
                                          SUNPSetupFn Pset, SUNPSolveFn Psol);

SUNDIALS_EXPORT
SUNErrCode SUNLinSolSetScalingVectors_PCG(SUNLinearSolver S, N_Vector s,
                                          N_Vector nul);

SUNDIALS_EXPORT
SUNErrCode SUNLinSolSetZeroGuess_PCG(SUNLinearSolver S, sunbooleantype onoff);

SUNDIALS_EXPORT
int SUNLinSolSetup_PCG(SUNLinearSolver S, SUNMatrix nul);

SUNDIALS_EXPORT
int SUNLinSolSolve_PCG(SUNLinearSolver S, SUNMatrix nul, N_Vector x, N_Vector b,
                       sunrealtype tol);

SUNDIALS_EXPORT
int SUNLinSolNumIters_PCG(SUNLinearSolver S);

SUNDIALS_EXPORT
sunrealtype SUNLinSolResNorm_PCG(SUNLinearSolver S);

SUNDIALS_EXPORT
N_Vector SUNLinSolResid_PCG(SUNLinearSolver S);

SUNDIALS_EXPORT
sunindextype SUNLinSolLastFlag_PCG(SUNLinearSolver S);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
SUNErrCode SUNLinSolSpace_PCG(SUNLinearSolver S, long int* lenrwLS,
                              long int* leniwLS);

SUNDIALS_EXPORT
SUNErrCode SUNLinSolFree_PCG(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

#endif
