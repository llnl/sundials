/* -----------------------------------------------------------------
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
 * This is the header file for a generic SUNDomEigEst package.
 * -----------------------------------------------------------------*/

#ifndef _SUNDOMEIGEST_H
#define _SUNDOMEIGEST_H

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------------------
 * Generic definition of SUNDomEigEstimator (DEE)
 * ----------------------------------------------------------------- */

/* Forward reference for pointer to SUNDomEigEstimator_Ops object */
typedef _SUNDIALS_STRUCT_ SUNDomEigEstimator_Ops_* SUNDomEigEstimator_Ops;

/* Forward reference for pointer to SUNDomEigEstimator object */
typedef _SUNDIALS_STRUCT_ SUNDomEigEstimator_* SUNDomEigEstimator;

/* Structure containing function pointers to estimator operations */
struct SUNDomEigEstimator_Ops_
{
  SUNErrCode (*setatimes)(SUNDomEigEstimator, void*, SUNATimesFn);
  SUNErrCode (*setoptions)(SUNDomEigEstimator DEE, const char* Did,
                           const char* file_name, int argc, char* argv[]);
  SUNErrCode (*setmaxiters)(SUNDomEigEstimator, long int);
  SUNErrCode (*setnumpreprocessiters)(SUNDomEigEstimator, int);
  SUNErrCode (*setreltol)(SUNDomEigEstimator, sunrealtype);
  SUNErrCode (*setinitialguess)(SUNDomEigEstimator, N_Vector);
  SUNErrCode (*initialize)(SUNDomEigEstimator);
  SUNErrCode (*estimate)(SUNDomEigEstimator, sunrealtype*, sunrealtype*);
  SUNErrCode (*getres)(SUNDomEigEstimator, sunrealtype*);
  SUNErrCode (*getnumiters)(SUNDomEigEstimator, long int*);
  SUNErrCode (*getnumatimescalls)(SUNDomEigEstimator, long int*);
  SUNErrCode (*write)(SUNDomEigEstimator, FILE*);
  SUNErrCode (*destroy)(SUNDomEigEstimator*);
};

/* An estimator is a structure with an implementation-dependent
   'content' field, and a pointer to a structure of estimator
   operations corresponding to that implementation. */
struct SUNDomEigEstimator_
{
  void* content;
  void* python;
  SUNDomEigEstimator_Ops ops;
  SUNContext sunctx;
};

/* -----------------------------------------------------------------
 * Functions exported by SUNDomEigEstimator module
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT
SUNDomEigEstimator SUNDomEigEstimator_NewEmpty(SUNContext sunctx);

SUNDIALS_EXPORT
void SUNDomEigEstimator_FreeEmpty(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_SetATimes(SUNDomEigEstimator DEE, void* A_data,
                                        SUNATimesFn ATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_SetOptions(SUNDomEigEstimator DEE,
                                         const char* Did, const char* file_name,
                                         int argc, char* argv[]);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_SetMaxIters(SUNDomEigEstimator DEE,
                                          long int max_iters);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_SetNumPreprocessIters(SUNDomEigEstimator DEE,
                                                    int num_iters);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_SetRelTol(SUNDomEigEstimator DEE, sunrealtype tol);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_SetInitialGuess(SUNDomEigEstimator DEE, N_Vector q);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_Initialize(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_Estimate(SUNDomEigEstimator DEE,
                                       sunrealtype* lambdaR,
                                       sunrealtype* lambdaI);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_GetRes(SUNDomEigEstimator DEE, sunrealtype* res);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_GetNumIters(SUNDomEigEstimator DEE,
                                          long int* num_iters);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_GetNumATimesCalls(SUNDomEigEstimator DEE,
                                                long int* num_ATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_Write(SUNDomEigEstimator DEE, FILE* outfile);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_Destroy(SUNDomEigEstimator* DEEptr);

#ifdef __cplusplus
}
#endif

#endif
