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
 *        N_VLinearSum_Parallel(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_PARALLEL_H
#define _NVECTOR_PARALLEL_H

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

struct _N_VectorContent_Parallel
{
  sunindextype local_length;  /* local vector length         */
  sunindextype global_length; /* global vector length        */
  sunbooleantype own_data;    /* ownership of data           */
  sunrealtype* data;          /* local data array            */
  MPI_Comm comm;              /* pointer to MPI communicator */
};

typedef struct _N_VectorContent_Parallel* N_VectorContent_Parallel;

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

SUNDIALS_EXPORT
sunindextype N_VGetLength_Parallel(N_Vector v);

SUNDIALS_EXPORT
sunindextype N_VGetLocalLength_Parallel(N_Vector v);

SUNDIALS_EXPORT
void N_VPrint_Parallel(N_Vector v);

SUNDIALS_EXPORT
void N_VPrintFile_Parallel(N_Vector v, FILE* outfile);

SUNDIALS_EXPORT
N_Vector_ID N_VGetVectorID_Parallel(N_Vector v);

SUNDIALS_EXPORT
N_Vector N_VCloneEmpty_Parallel(N_Vector w);

SUNDIALS_EXPORT
N_Vector N_VClone_Parallel(N_Vector w);

SUNDIALS_EXPORT
void N_VDestroy_Parallel(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
void N_VSpace_Parallel(N_Vector v, sunindextype* lrw, sunindextype* liw);

SUNDIALS_EXPORT
sunrealtype* N_VGetArrayPointer_Parallel(N_Vector v);

SUNDIALS_EXPORT
void N_VSetArrayPointer_Parallel(sunrealtype* v_data, N_Vector v);

SUNDIALS_EXPORT
MPI_Comm N_VGetCommunicator_Parallel(N_Vector v);

/* standard vector operations */

SUNDIALS_EXPORT
void N_VLinearSum_Parallel(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                           N_Vector z);

SUNDIALS_EXPORT
void N_VConst_Parallel(sunrealtype c, N_Vector z);

SUNDIALS_EXPORT
void N_VProd_Parallel(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_EXPORT
void N_VDiv_Parallel(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_EXPORT
void N_VScale_Parallel(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_EXPORT
void N_VAbs_Parallel(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
void N_VInv_Parallel(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
void N_VAddConst_Parallel(N_Vector x, sunrealtype b, N_Vector z);

SUNDIALS_EXPORT
sunrealtype N_VDotProd_Parallel(N_Vector x, N_Vector y);

SUNDIALS_EXPORT
sunrealtype N_VMaxNorm_Parallel(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VWrmsNorm_Parallel(N_Vector x, N_Vector w);

SUNDIALS_EXPORT
sunrealtype N_VWrmsNormMask_Parallel(N_Vector x, N_Vector w, N_Vector id);

SUNDIALS_EXPORT
sunrealtype N_VMin_Parallel(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VWL2Norm_Parallel(N_Vector x, N_Vector w);

SUNDIALS_EXPORT
sunrealtype N_VL1Norm_Parallel(N_Vector x);

SUNDIALS_EXPORT
void N_VCompare_Parallel(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_EXPORT
sunbooleantype N_VInvTest_Parallel(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
sunbooleantype N_VConstrMask_Parallel(N_Vector c, N_Vector x, N_Vector m);

SUNDIALS_EXPORT
sunrealtype N_VMinQuotient_Parallel(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT
SUNErrCode N_VLinearCombination_Parallel(int nvec, sunrealtype* c, N_Vector* V,
                                         N_Vector z);

SUNDIALS_EXPORT
SUNErrCode N_VScaleAddMulti_Parallel(int nvec, sunrealtype* a, N_Vector x,
                                     N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT
SUNErrCode N_VDotProdMulti_Parallel(int nvec, N_Vector x, N_Vector* Y,
                                    sunrealtype* dotprods);

/* vector array operations */
SUNDIALS_EXPORT
SUNErrCode N_VLinearSumVectorArray_Parallel(int nvec, sunrealtype a,
                                            N_Vector* X, sunrealtype b,
                                            N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT
SUNErrCode N_VScaleVectorArray_Parallel(int nvec, sunrealtype* c, N_Vector* X,
                                        N_Vector* Z);
SUNDIALS_EXPORT
SUNErrCode N_VConstVectorArray_Parallel(int nvecs, sunrealtype c, N_Vector* Z);

SUNDIALS_EXPORT
SUNErrCode N_VWrmsNormVectorArray_Parallel(int nvecs, N_Vector* X, N_Vector* W,
                                           sunrealtype* nrm);
SUNDIALS_EXPORT
SUNErrCode N_VWrmsNormMaskVectorArray_Parallel(int nvec, N_Vector* X, N_Vector* W,
                                               N_Vector id, sunrealtype* nrm);
SUNDIALS_EXPORT
SUNErrCode N_VScaleAddMultiVectorArray_Parallel(int nvec, int nsum,
                                                sunrealtype* a, N_Vector* X,
                                                N_Vector** Y, N_Vector** Z);
SUNDIALS_EXPORT
SUNErrCode N_VLinearCombinationVectorArray_Parallel(int nvec, int nsum,
                                                    sunrealtype* c,
                                                    N_Vector** X, N_Vector* Z);

/* OPTIONAL local reduction kernels (no parallel communication) */

SUNDIALS_EXPORT
sunrealtype N_VDotProdLocal_Parallel(N_Vector x, N_Vector y);

SUNDIALS_EXPORT
sunrealtype N_VMaxNormLocal_Parallel(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VMinLocal_Parallel(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VL1NormLocal_Parallel(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VWSqrSumLocal_Parallel(N_Vector x, N_Vector w);

SUNDIALS_EXPORT
sunrealtype N_VWSqrSumMaskLocal_Parallel(N_Vector x, N_Vector w, N_Vector id);

SUNDIALS_EXPORT
sunbooleantype N_VInvTestLocal_Parallel(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
sunbooleantype N_VConstrMaskLocal_Parallel(N_Vector c, N_Vector x, N_Vector m);

SUNDIALS_EXPORT
sunrealtype N_VMinQuotientLocal_Parallel(N_Vector num, N_Vector denom);

/* OPTIONAL single buffer reduction operations */

SUNDIALS_EXPORT
SUNErrCode N_VDotProdMultiLocal_Parallel(int nvec, N_Vector x, N_Vector* Y,
                                         sunrealtype* dotprods);

SUNDIALS_EXPORT
SUNErrCode N_VDotProdMultiAllReduce_Parallel(int nvec_total, N_Vector x,
                                             sunrealtype* dotprods);

/* OPTIONAL XBraid interface operations */

SUNDIALS_EXPORT
SUNErrCode N_VBufSize_Parallel(N_Vector x, sunindextype* size);

SUNDIALS_EXPORT
SUNErrCode N_VBufPack_Parallel(N_Vector x, void* buf);

SUNDIALS_EXPORT
SUNErrCode N_VBufUnpack_Parallel(N_Vector x, void* buf);

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
