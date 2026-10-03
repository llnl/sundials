/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                and Aaron Collier @ LLNL
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
 * This is the main header file for the MPI-enabled implementation
 * of the NVECTOR module.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be
 *     found in the header file sundials_nvector.h.
 *
 *   - The definition of the type sunrealtype can be found in the
 *     header file sundials_types.h, and it may be changed (at the
 *     configuration stage) according to the user's needs.
 *     The sundials_types.h file also contains the definition
 *     for the type sunbooleantype.
 *
 *   - N_Vector arguments to arithmetic vector operations need not
 *     be distinct. For example, the following call:
 *
 *        N_VLinearSum(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_NVECTOR_PARALLEL_H
#define SUNDIALS_NVECTOR_PARALLEL_H

#include <mpi.h>
#include <stdio.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_mpi_types.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Parallel implementation of N_Vector
 * -----------------------------------------------------------------
 */

struct N_VectorContent_Parallel_
{
  sunindextype local_length;  /* local vector length         */
  sunindextype global_length; /* global vector length        */
  sunbooleantype own_data;    /* ownership of data           */
  sunrealtype* data;          /* local data array            */
  MPI_Comm comm;              /* pointer to MPI communicator */
};

typedef struct N_VectorContent_Parallel_* N_VectorContent_Parallel;

/*
 * -----------------------------------------------------------------
 * Macros NV_CONTENT_P, NV_DATA_P, NV_OWN_DATA_P,
 *        NV_LOCLENGTH_P, NV_GLOBLENGTH_P,NV_COMM_P, and NV_Ith_P
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_P(v) ((N_VectorContent_Parallel)(v->content))

#define NV_LOCLENGTH_P(v) (NV_CONTENT_P(v)->local_length)

#define NV_GLOBLENGTH_P(v) (NV_CONTENT_P(v)->global_length)

#define NV_OWN_DATA_P(v) (NV_CONTENT_P(v)->own_data)

#define NV_DATA_P(v) (NV_CONTENT_P(v)->data)

#define NV_COMM_P(v) (NV_CONTENT_P(v)->comm)

#define NV_Ith_P(v, i) (NV_DATA_P(v)[i])

/*
 * -----------------------------------------------------------------
 * Functions exported by nvector_parallel
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
N_Vector N_VNew_Parallel(MPI_Comm comm, sunindextype local_length,
                         sunindextype global_length, SUNContext sunctx);

SUNDIALS_EXPORT
N_Vector N_VNewEmpty_Parallel(MPI_Comm comm, sunindextype local_length,
                              sunindextype global_length, SUNContext sunctx);

SUNDIALS_EXPORT
N_Vector N_VMake_Parallel(MPI_Comm comm, sunindextype local_length,
                          sunindextype global_length, sunrealtype* v_data_1d,
                          SUNContext sunctx);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
SUNErrCode N_VEnableFusedOps_Parallel(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombination_Parallel(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMulti_Parallel(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableDotProdMulti_Parallel(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearSumVectorArray_Parallel(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleVectorArray_Parallel(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableConstVectorArray_Parallel(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormVectorArray_Parallel(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormMaskVectorArray_Parallel(N_Vector v,
                                                     sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMultiVectorArray_Parallel(N_Vector v,
                                                      sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombinationVectorArray_Parallel(N_Vector v,
                                                          sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableDotProdMultiLocal_Parallel(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableDotProdMultiLocal_Parallel(N_Vector v, sunbooleantype tf);

#ifdef __cplusplus
}
#endif

#endif
