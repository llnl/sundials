/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ UMBC
 * Based on code sundials_spbcgs.h by: Peter Brown and
 *     Aaron Collier @ LLNL
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
 * This is the header file for the SPBCGS implementation of the
 * SUNLINSOL module, SUNLINSOL_SPBCGS.  The SPBCGS algorithm is based
 * on the Scaled Preconditioned Bi-CG-Stabilized method.
 *
 * Note:
 *   - The definition of the generic SUNLinearSolver structure can
 *     be found in the header file sundials_linearsolver.h.
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_SPBCGS_H
#define _SUNLINSOL_SPBCGS_H

#include <stdio.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default SPBCGS solver parameters */
#define SUNSPBCGS_MAXL_DEFAULT 5

/* -----------------------------------------
 * SPBCGS Implementation of SUNLinearSolver
 * ---------------------------------------- */

struct _SUNLinearSolverContent_SPBCGS
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

  N_Vector s1;
  N_Vector s2;
  N_Vector r;
  N_Vector r_star;
  N_Vector p;
  N_Vector q;
  N_Vector u;
  N_Vector Ap;
  N_Vector vtemp;
};

typedef struct _SUNLinearSolverContent_SPBCGS* SUNLinearSolverContent_SPBCGS;

/* ---------------------------------------
 *Exported Functions for SUNLINSOL_SPBCGS
 * --------------------------------------- */

SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_SPBCGS(N_Vector y, int pretype,
                                                 int maxl, SUNContext sunctx);
SUNDIALS_EXPORT SUNErrCode SUNLinSol_SPBCGSSetPrecType(SUNLinearSolver S,
                                                       int pretype);
SUNDIALS_EXPORT SUNErrCode SUNLinSol_SPBCGSSetMaxl(SUNLinearSolver S, int maxl);
SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_SPBCGS(SUNLinearSolver S);
SUNDIALS_EXPORT SUNLinearSolver_ID SUNLinSolGetID_SPBCGS(SUNLinearSolver S);
SUNDIALS_EXPORT SUNErrCode SUNLinSolInitialize_SPBCGS(SUNLinearSolver S);
SUNDIALS_EXPORT SUNErrCode SUNLinSolSetATimes_SPBCGS(SUNLinearSolver S,
                                                     void* A_data,
                                                     SUNATimesFn ATimes);
SUNDIALS_EXPORT SUNErrCode SUNLinSolSetPreconditioner_SPBCGS(SUNLinearSolver S,
                                                             void* P_data,
                                                             SUNPSetupFn Pset,
                                                             SUNPSolveFn Psol);
SUNDIALS_EXPORT SUNErrCode SUNLinSolSetScalingVectors_SPBCGS(SUNLinearSolver S,
                                                             N_Vector s1,
                                                             N_Vector s2);
SUNDIALS_EXPORT SUNErrCode SUNLinSolSetZeroGuess_SPBCGS(SUNLinearSolver S,
                                                        sunbooleantype onoff);
SUNDIALS_EXPORT int SUNLinSolSetup_SPBCGS(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_SPBCGS(SUNLinearSolver S, SUNMatrix A,
                                          N_Vector x, N_Vector b,
                                          sunrealtype tol);
SUNDIALS_EXPORT int SUNLinSolNumIters_SPBCGS(SUNLinearSolver S);
SUNDIALS_EXPORT sunrealtype SUNLinSolResNorm_SPBCGS(SUNLinearSolver S);
SUNDIALS_EXPORT N_Vector SUNLinSolResid_SPBCGS(SUNLinearSolver S);
SUNDIALS_EXPORT sunindextype SUNLinSolLastFlag_SPBCGS(SUNLinearSolver S);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
SUNErrCode SUNLinSolSpace_SPBCGS(SUNLinearSolver S, long int* lenrwLS,
                                 long int* leniwLS);
SUNDIALS_EXPORT SUNErrCode SUNLinSolFree_SPBCGS(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

#endif
