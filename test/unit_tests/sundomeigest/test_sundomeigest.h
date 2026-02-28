/*
 * -----------------------------------------------------------------
 * Programmer(s): Mustafa Aggul @ SMU
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
 * test SUNDomEigEstimator module implementations.
 * -----------------------------------------------------------------
 */

#include <sundials/sundials_core.h>

/* define constants */
#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Test function declarations */
int Test_SUNDomEigEstimator_SetATimes(SUNDomEigEstimator DEE, void* ATdata,
                                      SUNATimesFn ATimes, int myid);
int Test_SUNDomEigEstimator_SetMaxIters(SUNDomEigEstimator DEE,
                                        long int max_iters, int myid);
int Test_SUNDomEigEstimator_SetNumPreprocessIters(SUNDomEigEstimator DEE,
                                                  int num_warmups, int myid);
int Test_SUNDomEigEstimator_SetRelTol(SUNDomEigEstimator DEE, sunrealtype tol,
                                      int myid);
int Test_SUNDomEigEstimator_SetInitialGuess(SUNDomEigEstimator DEE, N_Vector q,
                                            int myid);
int Test_SUNDomEigEstimator_Initialize(SUNDomEigEstimator DEE, int myid);
int Test_SUNDomEigEstimator_Estimate(SUNDomEigEstimator DEE, sunrealtype* lambdaR,
                                     sunrealtype* lambdaI, int myid);
int Test_SUNDomEigEstimator_GetRes(SUNDomEigEstimator DEE, sunrealtype* cur_res,
                                   int myid);
int Test_SUNDomEigEstimator_GetNumIters(SUNDomEigEstimator DEE,
                                        long int* curniter, int myid);
int Test_SUNDomEigEstimator_GetNumATimesCalls(SUNDomEigEstimator DEE,
                                              long int* num_ATimes, int myid);
int Test_SUNDomEigEstimator_Write(SUNDomEigEstimator DEE, int myid);

/* Timing function */
void SetTiming(int onoff);

#ifdef __cplusplus
}
#endif
