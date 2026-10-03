/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner and Carol S. Woodward @ LLNL
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
 * This is the header file for the OpenMP implementation of the
 * NVECTOR module.
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
 *       N_VLinearSum(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_NVECTOR_OPENMP_H
#define SUNDIALS_NVECTOR_OPENMP_H

#include <stdio.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * OpenMP implementation of N_Vector
 * -----------------------------------------------------------------
 */

struct N_VectorContent_OpenMP_
{
  sunindextype length;     /* vector length            */
  sunbooleantype own_data; /* data ownership flag      */
  sunrealtype* data;       /* data array               */
  int num_threads;         /* number of OpenMP threads */
};

typedef struct N_VectorContent_OpenMP_* N_VectorContent_OpenMP;

/*
 * -----------------------------------------------------------------
 * Macros NV_CONTENT_OMP, NV_DATA_OMP, NV_OWN_DATA_OMP,
 *        NV_LENGTH_OMP, and NV_Ith_OMP
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_OMP(v) ((N_VectorContent_OpenMP)(v->content))

#define NV_LENGTH_OMP(v) (NV_CONTENT_OMP(v)->length)

#define NV_NUM_THREADS_OMP(v) (NV_CONTENT_OMP(v)->num_threads)

#define NV_OWN_DATA_OMP(v) (NV_CONTENT_OMP(v)->own_data)

#define NV_DATA_OMP(v) (NV_CONTENT_OMP(v)->data)

#define NV_Ith_OMP(v, i) (NV_DATA_OMP(v)[i])

/*
 * -----------------------------------------------------------------
 * Functions exported by nvector_openmp
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
N_Vector N_VNew_OpenMP(sunindextype vec_length, int num_threads,
                       SUNContext sunctx);

SUNDIALS_EXPORT
N_Vector N_VNewEmpty_OpenMP(sunindextype vec_length, int num_threads,
                            SUNContext sunctx);

SUNDIALS_EXPORT
N_Vector N_VMake_OpenMP(sunindextype vec_length, sunrealtype* v_data_1d,
                        int num_threads, SUNContext sunctx);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
SUNErrCode N_VEnableFusedOps_OpenMP(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombination_OpenMP(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMulti_OpenMP(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableDotProdMulti_OpenMP(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearSumVectorArray_OpenMP(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleVectorArray_OpenMP(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableConstVectorArray_OpenMP(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormVectorArray_OpenMP(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormMaskVectorArray_OpenMP(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMultiVectorArray_OpenMP(N_Vector v,
                                                    sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombinationVectorArray_OpenMP(N_Vector v,
                                                        sunbooleantype tf);

#ifdef __cplusplus
}
#endif

#endif
