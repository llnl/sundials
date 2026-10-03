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

#ifndef SUNDIALS_SUNLINSOL_PCG_H
#define SUNDIALS_SUNLINSOL_PCG_H

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

#ifdef __cplusplus
}
#endif

#endif
