/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
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
 * This is the main header file for the "ManyVector" implementation
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
 *       N_VLinearSum_ManyVector(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_MANY_VECTOR_H
#define _NVECTOR_MANY_VECTOR_H

#include <stdio.h>
#include <sundials/sundials_core.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------------------
   ManyVector implementation of N_Vector
   ----------------------------------------------------------------- */

struct _N_VectorContent_ManyVector
{
  sunindextype num_subvectors; /* number of vectors attached       */
  sunindextype global_length;  /* overall global manyvector length */
  N_Vector* subvec_array;      /* pointer to N_Vector array        */
  sunbooleantype own_data;     /* flag indicating data ownership   */
};

typedef struct _N_VectorContent_ManyVector* N_VectorContent_ManyVector;

/* -----------------------------------------------------------------
   functions exported by ManyVector
   ----------------------------------------------------------------- */

SUNDIALS_EXPORT
N_Vector N_VNew_ManyVector(sunindextype num_subvectors, N_Vector* vec_array_1d,
                           SUNContext sunctx);

SUNDIALS_EXPORT
N_Vector N_VGetSubvector_ManyVector(N_Vector v,
                                    sunindextype vec_num); // nb::rv_policy::reference

SUNDIALS_EXPORT
sunrealtype* N_VGetSubvectorArrayPointer_ManyVector(N_Vector v,
                                                    sunindextype vec_num);

SUNDIALS_EXPORT
SUNErrCode N_VSetSubvectorArrayPointer_ManyVector(sunrealtype* v_data_1d,
                                                  N_Vector v,
                                                  sunindextype vec_num);

SUNDIALS_EXPORT
sunindextype N_VGetNumSubvectors_ManyVector(N_Vector v);

/* standard vector operations */

SUNDIALS_EXPORT
N_Vector_ID N_VGetVectorID_ManyVector(N_Vector v);

SUNDIALS_EXPORT
void N_VPrint_ManyVector(N_Vector v);

SUNDIALS_EXPORT
void N_VPrintFile_ManyVector(N_Vector v, FILE* outfile);

SUNDIALS_EXPORT
N_Vector N_VCloneEmpty_ManyVector(N_Vector w);

SUNDIALS_EXPORT
N_Vector N_VClone_ManyVector(N_Vector w);

SUNDIALS_EXPORT
void N_VDestroy_ManyVector(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
void N_VSpace_ManyVector(N_Vector v, sunindextype* lrw, sunindextype* liw);

SUNDIALS_EXPORT
sunindextype N_VGetLength_ManyVector(N_Vector v);

SUNDIALS_EXPORT
sunindextype N_VGetSubvectorLocalLength_ManyVector(N_Vector v,
                                                   sunindextype vec_num);

SUNDIALS_EXPORT
void N_VLinearSum_ManyVector(sunrealtype a, N_Vector x, sunrealtype b,
                             N_Vector y, N_Vector z);

SUNDIALS_EXPORT
void N_VConst_ManyVector(sunrealtype c, N_Vector z);

SUNDIALS_EXPORT
void N_VProd_ManyVector(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_EXPORT
void N_VDiv_ManyVector(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_EXPORT
void N_VScale_ManyVector(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_EXPORT
void N_VAbs_ManyVector(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
void N_VInv_ManyVector(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
void N_VAddConst_ManyVector(N_Vector x, sunrealtype b, N_Vector z);

SUNDIALS_EXPORT
sunrealtype N_VWrmsNorm_ManyVector(N_Vector x, N_Vector w);

SUNDIALS_EXPORT
sunrealtype N_VWrmsNormMask_ManyVector(N_Vector x, N_Vector w, N_Vector id);

SUNDIALS_EXPORT
sunrealtype N_VWL2Norm_ManyVector(N_Vector x, N_Vector w);

SUNDIALS_EXPORT
void N_VCompare_ManyVector(sunrealtype c, N_Vector x, N_Vector z);

/* fused vector operations */

SUNDIALS_EXPORT
SUNErrCode N_VLinearCombination_ManyVector(int nvec, sunrealtype* c,
                                           N_Vector* V, N_Vector z);

SUNDIALS_EXPORT
SUNErrCode N_VScaleAddMulti_ManyVector(int nvec, sunrealtype* a, N_Vector x,
                                       N_Vector* Y, N_Vector* Z);

SUNDIALS_EXPORT
SUNErrCode N_VDotProdMulti_ManyVector(int nvec, N_Vector x, N_Vector* Y,
                                      sunrealtype* dotprods);

/* vector array operations */

SUNDIALS_EXPORT
SUNErrCode N_VLinearSumVectorArray_ManyVector(int nvec, sunrealtype a,
                                              N_Vector* X, sunrealtype b,
                                              N_Vector* Y, N_Vector* Z);

SUNDIALS_EXPORT
SUNErrCode N_VScaleVectorArray_ManyVector(int nvec, sunrealtype* c, N_Vector* X,
                                          N_Vector* Z);

SUNDIALS_EXPORT
SUNErrCode N_VConstVectorArray_ManyVector(int nvecs, sunrealtype c, N_Vector* Z);

SUNDIALS_EXPORT
SUNErrCode N_VWrmsNormVectorArray_ManyVector(int nvecs, N_Vector* X,
                                             N_Vector* W, sunrealtype* nrm);

SUNDIALS_EXPORT
SUNErrCode N_VWrmsNormMaskVectorArray_ManyVector(int nvec, N_Vector* X,
                                                 N_Vector* W, N_Vector id,
                                                 sunrealtype* nrm);

/* OPTIONAL local reduction kernels (no parallel communication) */

SUNDIALS_EXPORT
sunrealtype N_VDotProdLocal_ManyVector(N_Vector x, N_Vector y);

SUNDIALS_EXPORT
sunrealtype N_VMaxNormLocal_ManyVector(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VMinLocal_ManyVector(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VL1NormLocal_ManyVector(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VWSqrSumLocal_ManyVector(N_Vector x, N_Vector w);

SUNDIALS_EXPORT
sunrealtype N_VWSqrSumMaskLocal_ManyVector(N_Vector x, N_Vector w, N_Vector id);

SUNDIALS_EXPORT
sunbooleantype N_VInvTestLocal_ManyVector(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
sunbooleantype N_VConstrMaskLocal_ManyVector(N_Vector c, N_Vector x, N_Vector m);

SUNDIALS_EXPORT
sunrealtype N_VMinQuotientLocal_ManyVector(N_Vector num, N_Vector denom);

/* OPTIONAL single buffer reduction operations */

SUNDIALS_EXPORT
SUNErrCode N_VDotProdMultiLocal_ManyVector(int nvec, N_Vector x, N_Vector* Y,
                                           sunrealtype* dotprods);

/* OPTIONAL XBraid interface operations */

SUNDIALS_EXPORT
SUNErrCode N_VBufSize_ManyVector(N_Vector x, sunindextype* size);

SUNDIALS_EXPORT
SUNErrCode N_VBufPack_ManyVector(N_Vector x, void* buf);

SUNDIALS_EXPORT
SUNErrCode N_VBufUnpack_ManyVector(N_Vector x, void* buf);

/* -----------------------------------------------------------------
   Enable / disable fused vector operations
   ----------------------------------------------------------------- */

SUNDIALS_EXPORT
SUNErrCode N_VEnableFusedOps_ManyVector(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombination_ManyVector(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMulti_ManyVector(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableDotProdMulti_ManyVector(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearSumVectorArray_ManyVector(N_Vector v,
                                                    sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleVectorArray_ManyVector(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableConstVectorArray_ManyVector(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormVectorArray_ManyVector(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormMaskVectorArray_ManyVector(N_Vector v,
                                                       sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableDotProdMultiLocal_ManyVector(N_Vector v, sunbooleantype tf);

#ifdef __cplusplus
}
#endif
#endif
