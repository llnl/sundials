/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
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
 * This is the header file contains the prototypes for functions to
 * test SUNLinearSolver module implementations.
 * -----------------------------------------------------------------
 */

#include <math.h>
#include <sundials/sundials_core.h>

/* define constants */
#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif
/* Forward declarations for implementation specific utility functions */
int check_vector(N_Vector expected, N_Vector computed, sunrealtype tol);
void sync_device(void);

/* Test function declarations */
int Test_SUNLinSolGetType(SUNLinearSolver S, SUNLinearSolver_Type suntype,
                          int myid);
int Test_SUNLinSolGetID(SUNLinearSolver S, SUNLinearSolver_ID sunid, int myid);
int Test_SUNLinSolLastFlag(SUNLinearSolver S, int myid);
int Test_SUNLinSolSpace(SUNLinearSolver S, int myid);
int Test_SUNLinSolNumIters(SUNLinearSolver S, int myid);
int Test_SUNLinSolResNorm(SUNLinearSolver S, int myid);
int Test_SUNLinSolResid(SUNLinearSolver S, int myid);
int Test_SUNLinSolSetATimes(SUNLinearSolver S, void* ATdata, SUNATimesFn ATimes,
                            int myid);
int Test_SUNLinSolSetPreconditioner(SUNLinearSolver S, void* Pdata,
                                    SUNPSetupFn PSetup, SUNPSolveFn PSolve,
                                    int myid);
int Test_SUNLinSolSetScalingVectors(SUNLinearSolver S, N_Vector s1, N_Vector s2,
                                    int myid);
int Test_SUNLinSolSetZeroGuess(SUNLinearSolver S, int myid);
int Test_SUNLinSolInitialize(SUNLinearSolver S, int myid);
int Test_SUNLinSolSetup(SUNLinearSolver S, SUNMatrix A, int myid);
int Test_SUNLinSolSolve(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b,
                        sunrealtype tol, sunbooleantype zeroguess, int myid);

/* Timing function */
void SetTiming(int onoff);

#ifdef __cplusplus
}
#endif
