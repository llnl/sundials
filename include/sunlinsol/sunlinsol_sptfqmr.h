/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ UMBC
 * Based on code sundials_sptfqmr.h by: Aaron Collier @ LLNL
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
 * This is the header file for the SPTFQMR implementation of the
 * SUNLINSOL module, SUNLINSOL_SPTFQMR.  The SPTFQMR algorithm is
 * based on the Scaled Preconditioned Transpose-free Quasi-Minimum
 * Residual method.
 *
 * Note:
 *   - The definition of the generic SUNLinearSolver structure can
 *     be found in the header file sundials_linearsolver.h.
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_SPTFQMR_H
#define _SUNLINSOL_SPTFQMR_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default SPTFQMR solver parameters */
#define SUNSPTFQMR_MAXL_DEFAULT 5

/* ------------------------------------------
 * SPTFQMR Implementation of SUNLinearSolver
 * ------------------------------------------ */

struct _SUNLinearSolverContent_SPTFQMR
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
  N_Vector r_star;
  N_Vector q;
  N_Vector d;
  N_Vector v;
  N_Vector p;
  N_Vector* r;
  N_Vector u;
  N_Vector vtemp1;
  N_Vector vtemp2;
  N_Vector vtemp3;
};

typedef struct _SUNLinearSolverContent_SPTFQMR* SUNLinearSolverContent_SPTFQMR;

/* -------------------------------------
 * Exported Functions SUNLINSOL_SPTFQMR
 * -------------------------------------- */

SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_SPTFQMR(N_Vector y, int pretype,
                                                  int maxl, SUNContext sunctx);
SUNDIALS_EXPORT SUNErrCode SUNLinSol_SPTFQMRSetPrecType(SUNLinearSolver S,
                                                        int pretype);
SUNDIALS_EXPORT SUNErrCode SUNLinSol_SPTFQMRSetMaxl(SUNLinearSolver S, int maxl);
SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_SPTFQMR(SUNLinearSolver S);
SUNDIALS_EXPORT SUNLinearSolver_ID SUNLinSolGetID_SPTFQMR(SUNLinearSolver S);
SUNDIALS_EXPORT SUNErrCode SUNLinSolInitialize_SPTFQMR(SUNLinearSolver S);
SUNDIALS_EXPORT SUNErrCode SUNLinSolSetATimes_SPTFQMR(SUNLinearSolver S,
                                                      void* A_data,
                                                      SUNATimesFn ATimes);
SUNDIALS_EXPORT SUNErrCode SUNLinSolSetPreconditioner_SPTFQMR(SUNLinearSolver S,
                                                              void* P_data,
                                                              SUNPSetupFn Pset,
                                                              SUNPSolveFn Psol);
SUNDIALS_EXPORT SUNErrCode SUNLinSolSetScalingVectors_SPTFQMR(SUNLinearSolver S,
                                                              N_Vector s1,
                                                              N_Vector s2);
SUNDIALS_EXPORT SUNErrCode SUNLinSolSetZeroGuess_SPTFQMR(SUNLinearSolver S,
                                                         sunbooleantype onoff);
SUNDIALS_EXPORT int SUNLinSolSetup_SPTFQMR(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_SPTFQMR(SUNLinearSolver S, SUNMatrix A,
                                           N_Vector x, N_Vector b,
                                           sunrealtype tol);
SUNDIALS_EXPORT int SUNLinSolNumIters_SPTFQMR(SUNLinearSolver S);
SUNDIALS_EXPORT sunrealtype SUNLinSolResNorm_SPTFQMR(SUNLinearSolver S);
SUNDIALS_EXPORT N_Vector SUNLinSolResid_SPTFQMR(SUNLinearSolver S);
SUNDIALS_EXPORT sunindextype SUNLinSolLastFlag_SPTFQMR(SUNLinearSolver S);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
SUNErrCode SUNLinSolSpace_SPTFQMR(SUNLinearSolver S, long int* lenrwLS,
                                  long int* leniwLS);
SUNDIALS_EXPORT SUNErrCode SUNLinSolFree_SPTFQMR(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

#endif
