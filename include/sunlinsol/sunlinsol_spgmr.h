/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ UMBC
 * Based on code sundials_spgmr.h by: Scott D. Cohen,
 *      Alan C. Hindmarsh and Radu Serban @ LLNL
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
 * This is the header file for the SPGMR implementation of the
 * SUNLINSOL module, SUNLINSOL_SPGMR.  The SPGMR algorithm is based
 * on the Scaled Preconditioned GMRES (Generalized Minimal Residual)
 * method.
 *
 * Note:
 *   - The definition of the generic SUNLinearSolver structure can
 *     be found in the header file sundials_linearsolver.h.
 * -----------------------------------------------------------------
 */

#ifndef SUNDIALS_SUNLINSOL_SPGMR_H
#define SUNDIALS_SUNLINSOL_SPGMR_H

#include <stdio.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default SPGMR solver parameters */
#define SUNSPGMR_MAXL_DEFAULT   5
#define SUNSPGMR_MAXRS_DEFAULT  0
#define SUNSPGMR_GSTYPE_DEFAULT SUN_MODIFIED_GS

/* ----------------------------------------
 * SPGMR Implementation of SUNLinearSolver
 * ---------------------------------------- */

struct _SUNLinearSolverContent_SPGMR
{
  int maxl;
  int pretype;
  int gstype;
  int max_restarts;
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
  N_Vector* V;
  sunrealtype** Hes;
  sunrealtype* givens;
  N_Vector xcor;
  sunrealtype* yg;
  N_Vector vtemp;

  sunrealtype* cv;
  N_Vector* Xv;
};

typedef struct _SUNLinearSolverContent_SPGMR* SUNLinearSolverContent_SPGMR;

/* ---------------------------------------
 * Exported Functions for SUNLINSOL_SPGMR
 * --------------------------------------- */

SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_SPGMR(N_Vector y, int pretype,
                                                int maxl, SUNContext sunctx);
SUNDIALS_EXPORT SUNErrCode SUNLinSol_SPGMRSetPrecType(SUNLinearSolver S,
                                                      int pretype);
SUNDIALS_EXPORT SUNErrCode SUNLinSol_SPGMRSetGSType(SUNLinearSolver S,
                                                    int gstype);
SUNDIALS_EXPORT SUNErrCode SUNLinSol_SPGMRSetMaxRestarts(SUNLinearSolver S,
                                                         int maxrs);

#ifdef __cplusplus
}
#endif

#endif
