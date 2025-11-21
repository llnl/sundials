/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                and Aaron Collier @ LLNL
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
 * This is the header file for the serial implementation of the
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
 *       N_VLinearSum_Serial(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_SERIAL_H
#define _NVECTOR_SERIAL_H

#include <stdio.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * SERIAL implementation of N_Vector
 * -----------------------------------------------------------------
 */

struct _N_VectorContent_Serial
{
  sunindextype length;     /* vector length       */
  sunbooleantype own_data; /* data ownership flag */
  sunrealtype* data;       /* data array          */
};

typedef struct _N_VectorContent_Serial* N_VectorContent_Serial;

/*
 * -----------------------------------------------------------------
 * Macros NV_CONTENT_S, NV_DATA_S, NV_OWN_DATA_S,
 *        NV_LENGTH_S, and NV_Ith_S
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_S(v) ((N_VectorContent_Serial)(v->content))

#define NV_LENGTH_S(v) (NV_CONTENT_S(v)->length)

#define NV_OWN_DATA_S(v) (NV_CONTENT_S(v)->own_data)

#define NV_DATA_S(v) (NV_CONTENT_S(v)->data)

#define NV_Ith_S(v, i) (NV_DATA_S(v)[i])

/*
 * -----------------------------------------------------------------
 * Functions exported by nvector_serial
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
N_Vector N_VNewEmpty_Serial(sunindextype vec_length, SUNContext sunctx);

SUNDIALS_EXPORT
N_Vector N_VNew_Serial(sunindextype vec_length, SUNContext sunctx);

SUNDIALS_EXPORT
N_Vector N_VMake_Serial(sunindextype vec_length, sunrealtype* v_data_1d,
                        SUNContext sunctx);

SUNDIALS_EXPORT
sunindextype N_VGetLength_Serial(N_Vector v);

SUNDIALS_EXPORT
void N_VPrint_Serial(N_Vector v);

SUNDIALS_EXPORT
void N_VPrintFile_Serial(N_Vector v, FILE* outfile);

SUNDIALS_EXPORT
N_Vector_ID N_VGetVectorID_Serial(N_Vector v);

SUNDIALS_EXPORT
N_Vector N_VCloneEmpty_Serial(N_Vector w);

SUNDIALS_EXPORT
N_Vector N_VClone_Serial(N_Vector w);

SUNDIALS_EXPORT
void N_VDestroy_Serial(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
void N_VSpace_Serial(N_Vector v, sunindextype* lrw, sunindextype* liw);

SUNDIALS_EXPORT
sunrealtype* N_VGetArrayPointer_Serial(N_Vector v);

SUNDIALS_EXPORT
void N_VSetArrayPointer_Serial(sunrealtype* v_data, N_Vector v);

/* standard vector operations */
SUNDIALS_EXPORT
void N_VLinearSum_Serial(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                         N_Vector z);
SUNDIALS_EXPORT
void N_VConst_Serial(sunrealtype c, N_Vector z);

SUNDIALS_EXPORT
void N_VProd_Serial(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_EXPORT
void N_VDiv_Serial(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_EXPORT
void N_VScale_Serial(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_EXPORT
void N_VAbs_Serial(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
void N_VInv_Serial(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
void N_VAddConst_Serial(N_Vector x, sunrealtype b, N_Vector z);

SUNDIALS_EXPORT
sunrealtype N_VDotProd_Serial(N_Vector x, N_Vector y);

SUNDIALS_EXPORT
sunrealtype N_VMaxNorm_Serial(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VWrmsNorm_Serial(N_Vector x, N_Vector w);

SUNDIALS_EXPORT
sunrealtype N_VWrmsNormMask_Serial(N_Vector x, N_Vector w, N_Vector id);

SUNDIALS_EXPORT
sunrealtype N_VMin_Serial(N_Vector x);

SUNDIALS_EXPORT
sunrealtype N_VWL2Norm_Serial(N_Vector x, N_Vector w);

SUNDIALS_EXPORT
sunrealtype N_VL1Norm_Serial(N_Vector x);

SUNDIALS_EXPORT
void N_VCompare_Serial(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_EXPORT
sunbooleantype N_VInvTest_Serial(N_Vector x, N_Vector z);

SUNDIALS_EXPORT
sunbooleantype N_VConstrMask_Serial(N_Vector c, N_Vector x, N_Vector m);

SUNDIALS_EXPORT
sunrealtype N_VMinQuotient_Serial(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT
SUNErrCode N_VLinearCombination_Serial(int nvec, sunrealtype* c, N_Vector* V,
                                       N_Vector z);
SUNDIALS_EXPORT
SUNErrCode N_VScaleAddMulti_Serial(int nvec, sunrealtype* a, N_Vector x,
                                   N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT
SUNErrCode N_VDotProdMulti_Serial(int nvec, N_Vector x, N_Vector* Y,
                                  sunrealtype* dotprods);

/* vector array operations */
SUNDIALS_EXPORT
SUNErrCode N_VLinearSumVectorArray_Serial(int nvec, sunrealtype a, N_Vector* X,
                                          sunrealtype b, N_Vector* Y,
                                          N_Vector* Z);

SUNDIALS_EXPORT
SUNErrCode N_VScaleVectorArray_Serial(int nvec, sunrealtype* c, N_Vector* X,
                                      N_Vector* Z);

SUNDIALS_EXPORT
SUNErrCode N_VConstVectorArray_Serial(int nvecs, sunrealtype c, N_Vector* Z);

SUNDIALS_EXPORT
SUNErrCode N_VWrmsNormVectorArray_Serial(int nvecs, N_Vector* X, N_Vector* W,
                                         sunrealtype* nrm);

SUNDIALS_EXPORT
SUNErrCode N_VWrmsNormMaskVectorArray_Serial(int nvecs, N_Vector* X, N_Vector* W,
                                             N_Vector id, sunrealtype* nrm);

SUNDIALS_EXPORT
SUNErrCode N_VScaleAddMultiVectorArray_Serial(int nvec, int nsum,
                                              sunrealtype* a, N_Vector* X,
                                              N_Vector** Y, N_Vector** Z);

SUNDIALS_EXPORT
SUNErrCode N_VLinearCombinationVectorArray_Serial(int nvec, int nsum,
                                                  sunrealtype* c, N_Vector** X,
                                                  N_Vector* Z);

/* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT
sunrealtype N_VWSqrSumLocal_Serial(N_Vector x, N_Vector w);

SUNDIALS_EXPORT
sunrealtype N_VWSqrSumMaskLocal_Serial(N_Vector x, N_Vector w, N_Vector id);

/* OPTIONAL XBraid interface operations */
SUNDIALS_EXPORT
SUNErrCode N_VBufSize_Serial(N_Vector x, sunindextype* size);

SUNDIALS_EXPORT
SUNErrCode N_VBufPack_Serial(N_Vector x, void* buf);

SUNDIALS_EXPORT
SUNErrCode N_VBufUnpack_Serial(N_Vector x, void* buf);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
SUNErrCode N_VEnableFusedOps_Serial(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombination_Serial(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMulti_Serial(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableDotProdMulti_Serial(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearSumVectorArray_Serial(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleVectorArray_Serial(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableConstVectorArray_Serial(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormVectorArray_Serial(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormMaskVectorArray_Serial(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMultiVectorArray_Serial(N_Vector v,
                                                    sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombinationVectorArray_Serial(N_Vector v,
                                                        sunbooleantype tf);

#ifdef __cplusplus
}
#endif

#endif
