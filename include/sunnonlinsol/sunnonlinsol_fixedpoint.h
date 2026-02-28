/* ---------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
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
 * ---------------------------------------------------------------------------
 * This is the header file for the SUNNonlinearSolver module implementation of
 * the Anderson-accelerated fixed-point method.
 *
 * Part I defines the solver-specific content structure.
 *
 * Part II contains prototypes for the solver constructor and operations.
 * ---------------------------------------------------------------------------*/

#ifndef _SUNNONLINSOL_FIXEDPOINT_H
#define _SUNNONLINSOL_FIXEDPOINT_H

#include <sundials/sundials_core.h>

#include "sundials/sundials_types.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*-----------------------------------------------------------------------------
  I. Content structure
  ---------------------------------------------------------------------------*/

struct _SUNNonlinearSolverContent_FixedPoint
{
  /* functions provided by the integrator */
  SUNNonlinSolSysFn Sys;        /* fixed-point iteration function */
  SUNNonlinSolConvTestFn CTest; /* convergence test function      */

  /* nonlinear solver variables */
  int m;                  /* number of acceleration vectors to use          */
  int* imap;              /* array of length m                              */
  sunbooleantype damping; /* flag to apply dampling in acceleration         */
  sunrealtype beta;       /* damping parameter                               */
  sunrealtype* R;         /* array of length m*m                            */
  sunrealtype* gamma;     /* array of length m                              */
  sunrealtype* cvals;     /* array of length m+1 for fused vector op        */
  N_Vector* df;           /* vector array of length m                       */
  N_Vector* dg;           /* vector array of length m                       */
  N_Vector* q;            /* vector array of length m                       */
  N_Vector* Xvecs;        /* array of length m+1 for fused vector op        */
  N_Vector yprev;         /* temporary vectors for performing solve         */
  N_Vector gy;
  N_Vector fold;
  N_Vector gold;
  N_Vector delta;      /* correction vector (change between 2 iterates)  */
  int curiter;         /* current iteration number in a solve attempt    */
  int maxiters;        /* maximum number of iterations per solve attempt */
  long int niters;     /* total number of iterations across all solves   */
  long int nconvfails; /* total number of convergence failures           */
  void* ctest_data;    /* data to pass to convergence test function      */
};

typedef struct _SUNNonlinearSolverContent_FixedPoint* SUNNonlinearSolverContent_FixedPoint;

/* -----------------------------------------------------------------------------
   II: Exported functions
   ---------------------------------------------------------------------------*/

/* Constructor to create solver and allocates memory */
SUNDIALS_EXPORT
SUNNonlinearSolver SUNNonlinSol_FixedPoint(N_Vector y, int m, SUNContext sunctx);

SUNDIALS_EXPORT
SUNNonlinearSolver SUNNonlinSol_FixedPointSens(int count, N_Vector y, int m,
                                               SUNContext sunctx);

/* core functions */
SUNDIALS_EXPORT
SUNNonlinearSolver_Type SUNNonlinSolGetType_FixedPoint(SUNNonlinearSolver NLS);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolInitialize_FixedPoint(SUNNonlinearSolver NLS);

SUNDIALS_EXPORT
int SUNNonlinSolSolve_FixedPoint(SUNNonlinearSolver NLS, N_Vector y0,
                                 N_Vector y, N_Vector w, sunrealtype tol,
                                 sunbooleantype callSetup, void* mem);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolFree_FixedPoint(SUNNonlinearSolver NLS);

/* set functions */
SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetSysFn_FixedPoint(SUNNonlinearSolver NLS,
                                           SUNNonlinSolSysFn SysFn);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetConvTestFn_FixedPoint(SUNNonlinearSolver NLS,
                                                SUNNonlinSolConvTestFn CTestFn,
                                                void* ctest_data);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetMaxIters_FixedPoint(SUNNonlinearSolver NLS,
                                              int maxiters);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetDamping_FixedPoint(SUNNonlinearSolver NLS,
                                             sunrealtype beta);

/* get functions */
SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetNumIters_FixedPoint(SUNNonlinearSolver NLS,
                                              long int* niters);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetCurIter_FixedPoint(SUNNonlinearSolver NLS, int* iter);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetNumConvFails_FixedPoint(SUNNonlinearSolver NLS,
                                                  long int* nconvfails);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetSysFn_FixedPoint(SUNNonlinearSolver NLS,
                                           SUNNonlinSolSysFn* SysFn);

#ifdef __cplusplus
}
#endif

#endif
