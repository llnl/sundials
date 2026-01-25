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
 * This is the header file for the Arnoldi Iteration
 * implementation of the SUNDomEigEst package.
 *
 * Note:
 *   - The definition of the generic SUNDomEigEstimator structure can
 *     be found in the header file sundials_domeigestimator.h.
 * -----------------------------------------------------------------
 */

#ifndef _SUNDOMEIGEST_ARNOLDI_H
#define _SUNDOMEIGEST_ARNOLDI_H

#include <sundials/sundials_domeigestimator.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------
 * Arnoldi Iteration Implementation of SUNDomEigEstimator
 * ----------------------------------------------------- */

struct SUNDomEigEstimatorContent_Arnoldi_
{
  SUNATimesFn ATimes; /* User provided ATimes function */
  void* ATdata;       /* ATimes function data*/

  /* Krylov subspace vectors */
  N_Vector* V;
  N_Vector q;

  int kry_dim;        /* Krylov subspace dimension */
  int num_warmups;    /* Number of preprocessing iterations */
  long int num_iters; /* Number of iterations in last Estimate call */

  long int num_ATimes; /* Number of ATimes calls */

  sunrealtype* LAPACK_A; /* The vector which holds rows of the Hessenberg matrix in the given order */
  sunrealtype* LAPACK_wr;    /* Real parts of eigenvalues */
  sunrealtype* LAPACK_wi;    /* Imaginary parts of eigenvalues */
  sunrealtype* LAPACK_work;  /* Workspace array */
  sunindextype LAPACK_lwork; /* Dimension of the workspace array */
  sunrealtype** LAPACK_arr;  /* an array to sort eigenvalues*/

  sunrealtype** Hes; /* Hessenberg matrix Hes */
};

typedef struct SUNDomEigEstimatorContent_Arnoldi_* SUNDomEigEstimatorContent_Arnoldi;

/* ---------------------------------------
 * Exported Functions for SUNDOMEIGEST_Arnoldi
 * --------------------------------------- */

SUNDIALS_EXPORT
SUNDomEigEstimator SUNDomEigEstimator_Arnoldi(N_Vector q, int kry_dim,
                                              SUNContext sunctx);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_SetATimes_Arnoldi(SUNDomEigEstimator DEE,
                                                void* A_data, SUNATimesFn ATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_SetNumPreprocessIters_Arnoldi(SUNDomEigEstimator DEE,
                                                            int num_iters);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_SetInitialGuess_Arnoldi(SUNDomEigEstimator DEE,
                                                      N_Vector q);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_Initialize_Arnoldi(SUNDomEigEstimator DEE);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_Estimate_Arnoldi(SUNDomEigEstimator DEE,
                                               sunrealtype* lambdaR,
                                               sunrealtype* lambdaI);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_GetNumIters_Arnoldi(SUNDomEigEstimator DEE,
                                                  long int* num_iters);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_GetNumATimesCalls_Arnoldi(SUNDomEigEstimator DEE,
                                                        long int* num_ATimes);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_Write_Arnoldi(SUNDomEigEstimator DEE,
                                            FILE* outfile);

SUNDIALS_EXPORT
SUNErrCode SUNDomEigEstimator_Destroy_Arnoldi(SUNDomEigEstimator* DEEptr);

#ifdef __cplusplus
}
#endif

#endif
