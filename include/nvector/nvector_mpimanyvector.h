/* -----------------------------------------------------------------
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
 * This is the main header file for the "MPIManyVector" implementation
 * of the NVECTOR module.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be found
 *     in the header file sundials_nvector.h.
 *
 *   - The definitions of the types 'sunrealtype' and 'sunindextype' can
 *     be found in the header file sundials_types.h, and it may be
 *     changed (at the configuration stage) according to the user's needs.
 *     The sundials_types.h file also contains the definition
 *     for the type 'sunbooleantype'.
 *
 *   - N_Vector arguments to arithmetic vector operations need not
 *     be distinct. For example, the following call:
 *
 *       N_VLinearSum_MPIManyVector(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_MPI_MANY_VECTOR_H
#define _NVECTOR_MPI_MANY_VECTOR_H

#include <mpi.h>
#include <stdio.h>
#include <sundials/sundials_mpi_types.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------------------
   ManyVector implementation of N_Vector
   ----------------------------------------------------------------- */

struct _N_VectorContent_MPIManyVector
{
  MPI_Comm comm;               /* overall MPI communicator        */
  sunindextype num_subvectors; /* number of vectors attached       */
  sunindextype global_length;  /* overall global manyvector length */
  N_Vector* subvec_array;      /* pointer to N_Vector array        */
  sunbooleantype own_data;     /* flag indicating data ownership   */
};

typedef struct _N_VectorContent_MPIManyVector* N_VectorContent_MPIManyVector;

/* -----------------------------------------------------------------
   functions exported by ManyVector
   ----------------------------------------------------------------- */

SUNDIALS_EXPORT
N_Vector N_VMake_MPIManyVector(MPI_Comm comm, sunindextype num_subvectors,
                               N_Vector* vec_array_1d, SUNContext sunctx);

SUNDIALS_EXPORT
N_Vector N_VNew_MPIManyVector(sunindextype num_subvectors,
                              N_Vector* vec_array_1d, SUNContext sunctx);

SUNDIALS_EXPORT
N_Vector N_VGetSubvector_MPIManyVector(N_Vector v, sunindextype vec_num);

SUNDIALS_EXPORT
sunrealtype* N_VGetSubvectorArrayPointer_MPIManyVector(N_Vector v,
                                                       sunindextype vec_num);

SUNDIALS_EXPORT
SUNErrCode N_VSetSubvectorArrayPointer_MPIManyVector(sunrealtype* v_data_1d,
                                                     N_Vector v,
                                                     sunindextype vec_num);

SUNDIALS_EXPORT
sunindextype N_VGetNumSubvectors_MPIManyVector(N_Vector v);

/* standard vector operations */
SUNDIALS_EXPORT
N_Vector_ID N_VGetVectorID_MPIManyVector(N_Vector v);

SUNDIALS_EXPORT
void N_VPrint_MPIManyVector(N_Vector v);

SUNDIALS_EXPORT
void N_VPrintFile_MPIManyVector(N_Vector v, FILE* outfile);

SUNDIALS_EXPORT
N_Vector N_VCloneEmpty_MPIManyVector(N_Vector w);

SUNDIALS_EXPORT
N_Vector N_VClone_MPIManyVector(N_Vector w);

SUNDIALS_EXPORT
void N_VDestroy_MPIManyVector(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
void N_VSpace_MPIManyVector(N_Vector v, sunindextype* lrw, sunindextype* liw);

SUNDIALS_EXPORT
MPI_Comm N_VGetCommunicator_MPIManyVector(N_Vector v);

SUNDIALS_EXPORT
sunindextype N_VGetLength_MPIManyVector(N_Vector v);

SUNDIALS_EXPORT
sunindextype N_VGetSubvectorLocalLength_MPIManyVector(N_Vector v,
                                                      sunindextype vec_num);

SUNDIALS_EXPORT
void N_VLinearSum_MPIManyVector(sunrealtype a, N_Vector x, sunrealtype b,
                                N_Vector y, N_Vector z);
SUNDIALS_EXPORT
void N_VConst_MPIManyVector(sunrealtype c, N_Vector z);

SUNDIALS_EXPORT
void N_VProd_MPIManyVector(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_EXPORT
void N_VDiv_MPIManyVector(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_EXPORT
void N_VScale_MPIManyVector(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_EXPORT
void N_VAbs_MPIManyVector(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
void N_VInv_MPIManyVector(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
void N_VAddConst_MPIManyVector(N_Vector x, sunrealtype b, N_Vector z);

SUNDIALS_EXPORT
sunrealtype N_VDotProd_MPIManyVector(N_Vector x, N_Vector y);

SUNDIALS_EXPORT
sunrealtype N_VMaxNorm_MPIManyVector(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VWrmsNorm_MPIManyVector(N_Vector x, N_Vector w);

SUNDIALS_EXPORT
sunrealtype N_VWrmsNormMask_MPIManyVector(N_Vector x, N_Vector w, N_Vector id);

SUNDIALS_EXPORT
sunrealtype N_VMin_MPIManyVector(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VWL2Norm_MPIManyVector(N_Vector x, N_Vector w);

SUNDIALS_EXPORT
sunrealtype N_VL1Norm_MPIManyVector(N_Vector x);

SUNDIALS_EXPORT
void N_VCompare_MPIManyVector(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_EXPORT
sunbooleantype N_VInvTest_MPIManyVector(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
sunbooleantype N_VConstrMask_MPIManyVector(N_Vector c, N_Vector x, N_Vector m);

SUNDIALS_EXPORT
sunrealtype N_VMinQuotient_MPIManyVector(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT
SUNErrCode N_VLinearCombination_MPIManyVector(int nvec, sunrealtype* c,
                                              N_Vector* V, N_Vector z);
SUNDIALS_EXPORT
SUNErrCode N_VScaleAddMulti_MPIManyVector(int nvec, sunrealtype* a, N_Vector x,
                                          N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT
SUNErrCode N_VDotProdMulti_MPIManyVector(int nvec, N_Vector x, N_Vector* Y,
                                         sunrealtype* dotprods);

/* single buffer reduction operations */
SUNDIALS_EXPORT
SUNErrCode N_VDotProdMultiLocal_MPIManyVector(int nvec, N_Vector x, N_Vector* Y,
                                              sunrealtype* dotprods);

SUNDIALS_EXPORT
SUNErrCode N_VDotProdMultiAllReduce_MPIManyVector(int nvec_total, N_Vector x,
                                                  sunrealtype* sum);

/* vector array operations */
SUNDIALS_EXPORT
SUNErrCode N_VLinearSumVectorArray_MPIManyVector(int nvec, sunrealtype a,
                                                 N_Vector* X, sunrealtype b,
                                                 N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT
SUNErrCode N_VScaleVectorArray_MPIManyVector(int nvec, sunrealtype* c,
                                             N_Vector* X, N_Vector* Z);

SUNDIALS_EXPORT
SUNErrCode N_VConstVectorArray_MPIManyVector(int nvecs, sunrealtype c,
                                             N_Vector* Z);

SUNDIALS_EXPORT
SUNErrCode N_VWrmsNormVectorArray_MPIManyVector(int nvecs, N_Vector* X,
                                                N_Vector* W, sunrealtype* nrm);

SUNDIALS_EXPORT
SUNErrCode N_VWrmsNormMaskVectorArray_MPIManyVector(int nvec, N_Vector* X,
                                                    N_Vector* W, N_Vector id,
                                                    sunrealtype* nrm);

/* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT
sunrealtype N_VDotProdLocal_MPIManyVector(N_Vector x, N_Vector y);

SUNDIALS_EXPORT
sunrealtype N_VMaxNormLocal_MPIManyVector(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VMinLocal_MPIManyVector(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VL1NormLocal_MPIManyVector(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VWSqrSumLocal_MPIManyVector(N_Vector x, N_Vector w);

SUNDIALS_EXPORT
sunrealtype N_VWSqrSumMaskLocal_MPIManyVector(N_Vector x, N_Vector w,
                                              N_Vector id);

SUNDIALS_EXPORT
sunbooleantype N_VInvTestLocal_MPIManyVector(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
sunbooleantype N_VConstrMaskLocal_MPIManyVector(N_Vector c, N_Vector x,
                                                N_Vector m);

SUNDIALS_EXPORT
sunrealtype N_VMinQuotientLocal_MPIManyVector(N_Vector num, N_Vector denom);

/* OPTIONAL XBraid interface operations */
SUNDIALS_EXPORT
SUNErrCode N_VBufSize_MPIManyVector(N_Vector x, sunindextype* size);

SUNDIALS_EXPORT
SUNErrCode N_VBufPack_MPIManyVector(N_Vector x, void* buf);

SUNDIALS_EXPORT
SUNErrCode N_VBufUnpack_MPIManyVector(N_Vector x, void* buf);

/* -----------------------------------------------------------------
   Enable / disable fused vector operations
   ----------------------------------------------------------------- */

SUNDIALS_EXPORT
SUNErrCode N_VEnableFusedOps_MPIManyVector(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombination_MPIManyVector(N_Vector v,
                                                    sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMulti_MPIManyVector(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableDotProdMulti_MPIManyVector(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearSumVectorArray_MPIManyVector(N_Vector v,
                                                       sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleVectorArray_MPIManyVector(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableConstVectorArray_MPIManyVector(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormVectorArray_MPIManyVector(N_Vector v,
                                                      sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormMaskVectorArray_MPIManyVector(N_Vector v,
                                                          sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableDotProdMultiLocal_MPIManyVector(N_Vector v,
                                                    sunbooleantype tf);

#ifdef __cplusplus
}
#endif
#endif
