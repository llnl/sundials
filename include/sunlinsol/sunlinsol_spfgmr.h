/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ UMBC
 * Based on code sundials_spfgmr.h by: Daniel R. Reynolds and
 *    Hilari C. Tiedeman @ UMBC
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
 * This is the header file for the SPFGMR implementation of the
 * SUNLINSOL module, SUNLINSOL_SPFGMR.  The SPFGMR algorithm is based
 * on the Scaled Preconditioned FGMRES (Flexible Generalized Minimal
 * Residual) method [Y. Saad, SIAM J. Sci. Comput., 1993].
 *
 * Note:
 *   - The definition of the generic SUNLinearSolver structure can
 *     be found in the header file sundials_linearsolver.h.
 * -----------------------------------------------------------------
 */

#ifndef SUNDIALS_SUNLINSOL_SPFGMR_H
#define SUNDIALS_SUNLINSOL_SPFGMR_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default SPFGMR solver parameters */
#define SUNSPFGMR_MAXL_DEFAULT   5
#define SUNSPFGMR_MAXRS_DEFAULT  0
#define SUNSPFGMR_GSTYPE_DEFAULT SUN_MODIFIED_GS

/* -----------------------------------------
 * SPFGMR Implementation of SUNLinearSolver
 * ----------------------------------------- */

struct SUNLinearSolverContent_SPFGMR_
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
  N_Vector* Z;
  sunrealtype** Hes;
  sunrealtype* givens;
  N_Vector xcor;
  sunrealtype* yg;
  N_Vector vtemp;

  sunrealtype* cv;
  N_Vector* Xv;
};

typedef struct SUNLinearSolverContent_SPFGMR_* SUNLinearSolverContent_SPFGMR;

/* ----------------------------------------
 * Exported Functions for SUNLINSOL_SPFGMR
 * ---------------------------------------- */

SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_SPFGMR(N_Vector y, int pretype,
                                                 int maxl, SUNContext sunctx);
SUNDIALS_EXPORT SUNErrCode SUNLinSol_SPFGMRSetPrecType(SUNLinearSolver S,
                                                       int pretype);
SUNDIALS_EXPORT SUNErrCode SUNLinSol_SPFGMRSetGSType(SUNLinearSolver S,
                                                     int gstype);
SUNDIALS_EXPORT SUNErrCode SUNLinSol_SPFGMRSetMaxRestarts(SUNLinearSolver S,
                                                          int maxrs);

#ifdef __cplusplus
}
#endif

#endif
