/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * evaluate the performance of an NVECTOR module implementation.
 * -----------------------------------------------------------------*/

#include <math.h>

/* define constants */
#define NEG_ONE SUN_RCONST(-1.0)
#define ZERO    SUN_RCONST(0.0)
#define ONE     SUN_RCONST(1.0)
#define TWO     SUN_RCONST(2.0)
#define TEN     SUN_RCONST(10.0)

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Forward declarations for implementation specific utility functions */
void collect_times(N_Vector X, double* times, int ntimes);
void sync_device(N_Vector x);
void ClearCache(void);
void N_VRand(N_Vector Xvec, sunindextype Xlen, sunrealtype lower,
             sunrealtype upper);
void N_VRandZeroOne(N_Vector Xvec, sunindextype Xlen);
void N_VRandConstraints(N_Vector Xvec, sunindextype Xlen);

/* Standard vector operation tests */
int Test_N_VLinearSum(N_Vector X, sunindextype local_length, int ntests);
int Test_N_VConst(N_Vector X, sunindextype local_length, int ntests);
int Test_N_VProd(N_Vector X, sunindextype local_length, int ntests);
int Test_N_VDiv(N_Vector X, sunindextype local_length, int ntests);
int Test_N_VScale(N_Vector X, sunindextype local_length, int ntests);
int Test_N_VAbs(N_Vector X, sunindextype local_length, int ntests);
int Test_N_VInv(N_Vector X, sunindextype local_length, int ntests);
int Test_N_VAddConst(N_Vector X, sunindextype local_length, int ntests);
int Test_N_VDotProd(N_Vector X, sunindextype local_length, int ntests);
int Test_N_VMaxNorm(N_Vector X, sunindextype local_length, int ntests);
int Test_N_VWrmsNorm(N_Vector X, sunindextype local_length, int ntests);
int Test_N_VWrmsNormMask(N_Vector X, sunindextype local_length, int ntests);
int Test_N_VMin(N_Vector X, sunindextype local_length, int ntests);
int Test_N_VWL2Norm(N_Vector X, sunindextype local_length, int ntests);
int Test_N_VL1Norm(N_Vector X, sunindextype local_length, int ntests);
int Test_N_VCompare(N_Vector X, sunindextype local_length, int ntests);
int Test_N_VInvTest(N_Vector X, sunindextype local_length, int ntests);
int Test_N_VConstrMask(N_Vector X, sunindextype local_length, int ntests);
int Test_N_VMinQuotient(N_Vector X, sunindextype local_length, int ntests);

/* Fused vector operation tests */
int Test_N_VLinearCombination(N_Vector X, sunindextype local_length, int nvecs,
                              int ntests);
int Test_N_VScaleAddMulti(N_Vector X, sunindextype local_length, int nvecs,
                          int ntests);
int Test_N_VDotProdMulti(N_Vector X, sunindextype local_length, int nvecs,
                         int ntests);

/* Vector array operation tests */
int Test_N_VLinearSumVectorArray(N_Vector X, sunindextype local_length,
                                 int nvecs, int ntests);
int Test_N_VScaleVectorArray(N_Vector X, sunindextype local_length, int nvecs,
                             int ntests);
int Test_N_VConstVectorArray(N_Vector X, sunindextype local_length, int nvecs,
                             int ntests);
int Test_N_VWrmsNormVectorArray(N_Vector X, sunindextype local_length,
                                int nvecs, int ntests);
int Test_N_VWrmsNormMaskVectorArray(N_Vector X, sunindextype local_length,
                                    int nvecs, int ntests);
int Test_N_VScaleAddMultiVectorArray(N_Vector X, sunindextype local_length,
                                     int nvecs, int nsums, int ntests);
int Test_N_VLinearCombinationVectorArray(N_Vector X, sunindextype local_length,
                                         int nvecs, int nsums, int tests);
/* Turn timing on/off */
void SetTiming(int onoff, int myid);

/* Set the number of warmup runs (1 by default) */
void SetNumWarmups(int num_warmups);

/* Print output table header */
void PrintTableHeader(int type);

/* Fill sunrealtype arrays with random data */
void rand_realtype(sunrealtype* data, sunindextype len, sunrealtype lower,
                   sunrealtype upper);
void rand_realtype_zero_one(sunrealtype* data, sunindextype len);
void rand_realtype_constraints(sunrealtype* data, sunindextype len);

#ifdef __cplusplus
}
#endif
