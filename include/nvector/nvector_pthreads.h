/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * Acknowledgements: This NVECTOR module is based on the NVECTOR
 *                   Serial module by Scott D. Cohen, Alan C.
 *                   Hindmarsh, Radu Serban, and Aaron Collier
 *                   @ LLNL
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
 * This is the header file for the POSIX Threads (Pthreads)
 * implementation of the NVECTOR module using LOCAL data structs
 * to share data between threads.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be found
 *     in the header file sundials_nvector.h.
 *
 *   - The definition of the type 'sunrealtype' can be found in the
 *     header file sundials_types.h, and it may be changed (at the
 *     configuration stage) according to the user's needs.
 *     The sundials_types.h file also contains the definition
 *     for the type 'sunbooleantype'.
 *
 *   - N_Vector arguments to arithmetic vector operations need not
 *     be distinct. For example, the following call:
 *
 *       N_VLinearSum_Pthreads(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_PTHREADS_H
#define _NVECTOR_PTHREADS_H

#include <pthread.h>
#include <stdio.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Pthreads implementation of N_Vector
 * -----------------------------------------------------------------
 */

struct _N_VectorContent_Pthreads
{
  sunindextype length;     /* vector length           */
  sunbooleantype own_data; /* data ownership flag     */
  sunrealtype* data;       /* data array              */
  int num_threads;         /* number of POSIX threads */
};

typedef struct _N_VectorContent_Pthreads* N_VectorContent_Pthreads;

/* Structure to hold parallelization information for each thread when
   calling "companion" functions to compute vector operations. The
   start and end vector (loop) indices are unique to each thread, the
   sunrealtype variables are the same for each thread, and the mutex
   variable is used to lock variables in reductions. */

struct _Pthreads_Data
{
  sunindextype start;            /* starting index for loop  */
  sunindextype end;              /* ending index for loop    */
  sunrealtype c1, c2;            /* scalar values            */
  sunrealtype *v1, *v2, *v3;     /* vector data              */
  sunrealtype* global_val;       /* shared global variable   */
  pthread_mutex_t* global_mutex; /* lock for shared variable */

  int nvec; /* number of vectors in fused op */
  int nsum; /* number of sums in fused op    */

  sunrealtype* cvals; /* scalar values in fused op */

  N_Vector x1; /* vector array in fused op */
  N_Vector x2; /* vector array in fused op */
  N_Vector x3; /* vector array in fused op */

  N_Vector* Y1; /* vector array in fused op */
  N_Vector* Y2; /* vector array in fused op */
  N_Vector* Y3; /* vector array in fused op */

  N_Vector** ZZ1; /* array of vector arrays in fused op */
  N_Vector** ZZ2; /* array of vector arrays in fused op */
};

typedef struct _Pthreads_Data Pthreads_Data;

/*
 * -----------------------------------------------------------------
 * Macros NV_CONTENT_PT, NV_DATA_PT, NV_OWN_DATA_PT,
 *        NV_LENGTH_PT, and NV_Ith_PT
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_PT(v) ((N_VectorContent_Pthreads)(v->content))

#define NV_LENGTH_PT(v) (NV_CONTENT_PT(v)->length)

#define NV_NUM_THREADS_PT(v) (NV_CONTENT_PT(v)->num_threads)

#define NV_OWN_DATA_PT(v) (NV_CONTENT_PT(v)->own_data)

#define NV_DATA_PT(v) (NV_CONTENT_PT(v)->data)

#define NV_Ith_PT(v, i) (NV_DATA_PT(v)[i])

/*
 * -----------------------------------------------------------------
 * Functions exported by nvector_Pthreads
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
N_Vector N_VNew_Pthreads(sunindextype vec_length, int n_threads,
                         SUNContext sunctx);

SUNDIALS_EXPORT
N_Vector N_VNewEmpty_Pthreads(sunindextype vec_length, int n_threads,
                              SUNContext sunctx);

SUNDIALS_EXPORT
N_Vector N_VMake_Pthreads(sunindextype vec_length, int n_threads,
                          sunrealtype* v_data, SUNContext sunctx);

SUNDIALS_EXPORT
sunindextype N_VGetLength_Pthreads(N_Vector v);

SUNDIALS_EXPORT
void N_VPrint_Pthreads(N_Vector v);

SUNDIALS_EXPORT
void N_VPrintFile_Pthreads(N_Vector v, FILE* outfile);

SUNDIALS_EXPORT
N_Vector_ID N_VGetVectorID_Pthreads(N_Vector v);

SUNDIALS_EXPORT
N_Vector N_VCloneEmpty_Pthreads(N_Vector w);

SUNDIALS_EXPORT
N_Vector N_VClone_Pthreads(N_Vector w);

SUNDIALS_EXPORT
void N_VDestroy_Pthreads(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
void N_VSpace_Pthreads(N_Vector v, sunindextype* lrw, sunindextype* liw);

SUNDIALS_EXPORT
sunrealtype* N_VGetArrayPointer_Pthreads(N_Vector v);

SUNDIALS_EXPORT
void N_VSetArrayPointer_Pthreads(sunrealtype* v_data, N_Vector v);

/* standard vector operations */
SUNDIALS_EXPORT
void N_VLinearSum_Pthreads(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                           N_Vector z);

SUNDIALS_EXPORT
void N_VConst_Pthreads(sunrealtype c, N_Vector z);

SUNDIALS_EXPORT
void N_VProd_Pthreads(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_EXPORT
void N_VDiv_Pthreads(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_EXPORT
void N_VScale_Pthreads(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_EXPORT
void N_VAbs_Pthreads(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
void N_VInv_Pthreads(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
void N_VAddConst_Pthreads(N_Vector x, sunrealtype b, N_Vector z);

SUNDIALS_EXPORT
sunrealtype N_VDotProd_Pthreads(N_Vector x, N_Vector y);

SUNDIALS_EXPORT
sunrealtype N_VMaxNorm_Pthreads(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VWrmsNorm_Pthreads(N_Vector x, N_Vector w);

SUNDIALS_EXPORT
sunrealtype N_VWrmsNormMask_Pthreads(N_Vector x, N_Vector w, N_Vector id);

SUNDIALS_EXPORT
sunrealtype N_VMin_Pthreads(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VWL2Norm_Pthreads(N_Vector x, N_Vector w);

SUNDIALS_EXPORT
sunrealtype N_VL1Norm_Pthreads(N_Vector x);

SUNDIALS_EXPORT
void N_VCompare_Pthreads(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_EXPORT
sunbooleantype N_VInvTest_Pthreads(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
sunbooleantype N_VConstrMask_Pthreads(N_Vector c, N_Vector x, N_Vector m);

SUNDIALS_EXPORT
sunrealtype N_VMinQuotient_Pthreads(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT
SUNErrCode N_VLinearCombination_Pthreads(int nvec, sunrealtype* c, N_Vector* X,
                                         N_Vector z);

SUNDIALS_EXPORT
SUNErrCode N_VScaleAddMulti_Pthreads(int nvec, sunrealtype* a, N_Vector x,
                                     N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT
SUNErrCode N_VDotProdMulti_Pthreads(int nvec, N_Vector x, N_Vector* Y,
                                    sunrealtype* dotprods);

/* vector array operations */
SUNDIALS_EXPORT
SUNErrCode N_VLinearSumVectorArray_Pthreads(int nvec, sunrealtype a,
                                            N_Vector* X, sunrealtype b,
                                            N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT
SUNErrCode N_VScaleVectorArray_Pthreads(int nvec, sunrealtype* c, N_Vector* X,
                                        N_Vector* Z);
SUNDIALS_EXPORT
SUNErrCode N_VConstVectorArray_Pthreads(int nvec, sunrealtype c, N_Vector* Z);

SUNDIALS_EXPORT
SUNErrCode N_VWrmsNormVectorArray_Pthreads(int nvec, N_Vector* X, N_Vector* W,
                                           sunrealtype* nrm);
SUNDIALS_EXPORT
SUNErrCode N_VWrmsNormMaskVectorArray_Pthreads(int nvec, N_Vector* X, N_Vector* W,
                                               N_Vector id, sunrealtype* nrm);
SUNDIALS_EXPORT
SUNErrCode N_VScaleAddMultiVectorArray_Pthreads(int nvec, int nsum,
                                                sunrealtype* a, N_Vector* X,
                                                N_Vector** Y, N_Vector** Z);
SUNDIALS_EXPORT
SUNErrCode N_VLinearCombinationVectorArray_Pthreads(int nvec, int nsum,
                                                    sunrealtype* c,
                                                    N_Vector** X, N_Vector* Z);

/* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT
sunrealtype N_VWSqrSumLocal_Pthreads(N_Vector x, N_Vector w);

SUNDIALS_EXPORT
sunrealtype N_VWSqrSumMaskLocal_Pthreads(N_Vector x, N_Vector w, N_Vector id);

/* OPTIONAL XBraid interface operations */
SUNDIALS_EXPORT
SUNErrCode N_VBufSize_Pthreads(N_Vector x, sunindextype* size);

SUNDIALS_EXPORT
SUNErrCode N_VBufPack_Pthreads(N_Vector x, void* buf);

SUNDIALS_EXPORT
SUNErrCode N_VBufUnpack_Pthreads(N_Vector x, void* buf);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
SUNErrCode N_VEnableFusedOps_Pthreads(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombination_Pthreads(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMulti_Pthreads(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableDotProdMulti_Pthreads(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearSumVectorArray_Pthreads(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleVectorArray_Pthreads(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableConstVectorArray_Pthreads(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormVectorArray_Pthreads(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormMaskVectorArray_Pthreads(N_Vector v,
                                                     sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMultiVectorArray_Pthreads(N_Vector v,
                                                      sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombinationVectorArray_Pthreads(N_Vector v,
                                                          sunbooleantype tf);

#ifdef __cplusplus
}
#endif

#endif
