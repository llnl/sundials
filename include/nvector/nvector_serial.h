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
 *       N_VLinearSum(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_NVECTOR_SERIAL_H
#define SUNDIALS_NVECTOR_SERIAL_H

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

struct N_VectorContent_Serial_
{
  sunindextype length;     /* vector length       */
  sunbooleantype own_data; /* data ownership flag */
  sunrealtype* data;       /* data array          */
};

typedef struct N_VectorContent_Serial_* N_VectorContent_Serial;

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
