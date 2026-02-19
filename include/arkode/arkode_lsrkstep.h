/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds and Mustafa Aggul @ UMBC
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
 * This is the header file for the ARKODE LSRKStep module.
 * -----------------------------------------------------------------*/

#ifndef _LSRKSTEP_H
#define _LSRKSTEP_H

#include <arkode/arkode.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef int (*ARKDomEigFn)(sunrealtype t, N_Vector y, N_Vector fn,
                           sunrealtype* lambdaR, sunrealtype* lambdaI,
                           void* user_data, N_Vector temp1, N_Vector temp2,
                           N_Vector temp3);

/* ------------------
 * LSRKStep Constants
 * ------------------ */

enum ARKODE_LSRKMethodType
{
  ARKODE_LSRK_RKC_2,
  ARKODE_LSRK_RKL_2,
  ARKODE_LSRK_SSP_S_2,
  ARKODE_LSRK_SSP_S_3,
  ARKODE_LSRK_SSP_10_4
};

#ifndef SWIG
typedef enum ARKODE_LSRKMethodType ARKODE_LSRKMethodType;
#endif

/* -------------------
 * Exported Functions
 * ------------------- */

/* Creation and Reinitialization functions */

SUNDIALS_EXPORT void* LSRKStepCreateSTS(ARKRhsFn rhs, sunrealtype t0,
                                        N_Vector y0, SUNContext sunctx);

SUNDIALS_EXPORT void* LSRKStepCreateSSP(ARKRhsFn rhs, sunrealtype t0,
                                        N_Vector y0, SUNContext sunctx);

SUNDIALS_EXPORT int LSRKStepReInitSTS(void* arkode_mem, ARKRhsFn rhs,
                                      sunrealtype t0, N_Vector y0);

SUNDIALS_EXPORT int LSRKStepReInitSSP(void* arkode_mem, ARKRhsFn rhs,
                                      sunrealtype t0, N_Vector y0);

/* Optional input functions -- must be called AFTER a creation routine above */

SUNDIALS_EXPORT int LSRKStepSetSTSMethod(void* arkode_mem,
                                         ARKODE_LSRKMethodType method);

SUNDIALS_EXPORT int LSRKStepSetSSPMethod(void* arkode_mem,
                                         ARKODE_LSRKMethodType method);

SUNDIALS_EXPORT int LSRKStepSetSTSMethodByName(void* arkode_mem,
                                               const char* emethod);

SUNDIALS_EXPORT int LSRKStepSetSSPMethodByName(void* arkode_mem,
                                               const char* emethod);

SUNDIALS_EXPORT int LSRKStepSetDomEigFn(void* arkode_mem, ARKDomEigFn dom_eig);

SUNDIALS_EXPORT int LSRKStepSetDomEigEstimator(void* arkode_mem,
                                               SUNDomEigEstimator DEE);

SUNDIALS_EXPORT int LSRKStepSetDomEigFrequency(void* arkode_mem, long int nsteps);

SUNDIALS_EXPORT int LSRKStepSetMaxNumStages(void* arkode_mem,
                                            int stage_max_limit);

SUNDIALS_EXPORT int LSRKStepSetDomEigSafetyFactor(void* arkode_mem,
                                                  sunrealtype dom_eig_safety);

SUNDIALS_EXPORT int LSRKStepSetNumDomEigEstInitPreprocessIters(void* arkode_mem,
                                                               int num_iters);

SUNDIALS_EXPORT int LSRKStepSetNumDomEigEstPreprocessIters(void* arkode_mem,
                                                           int num_iters);

SUNDIALS_EXPORT int LSRKStepSetNumSSPStages(void* arkode_mem, int num_of_stages);

/* Optional output functions */

SUNDIALS_EXPORT int LSRKStepGetNumDomEigUpdates(void* arkode_mem,
                                                long int* dom_eig_num_evals);

SUNDIALS_EXPORT int LSRKStepGetMaxNumStages(void* arkode_mem, int* stage_max);

SUNDIALS_EXPORT int LSRKStepGetNumDomEigEstRhsEvals(void* arkode_mem,
                                                    long int* nfeDQ);

SUNDIALS_EXPORT int LSRKStepGetNumDomEigEstIters(void* arkode_mem,
                                                 long int* num_iters);

#ifdef __cplusplus
}
#endif

#endif
