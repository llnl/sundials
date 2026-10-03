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

#ifndef SUNDIALS_SUNLINSOL_SPBCGS_H
#define SUNDIALS_SUNLINSOL_SPBCGS_H

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

struct SUNLinearSolverContent_SPBCGS_
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

typedef struct SUNLinearSolverContent_SPBCGS_* SUNLinearSolverContent_SPBCGS;

/* ---------------------------------------
 *Exported Functions for SUNLINSOL_SPBCGS
 * --------------------------------------- */

SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_SPBCGS(N_Vector y, int pretype,
                                                 int maxl, SUNContext sunctx);
SUNDIALS_EXPORT SUNErrCode SUNLinSol_SPBCGSSetPrecType(SUNLinearSolver S,
                                                       int pretype);
SUNDIALS_EXPORT SUNErrCode SUNLinSol_SPBCGSSetMaxl(SUNLinearSolver S, int maxl);

#ifdef __cplusplus
}
#endif

#endif
