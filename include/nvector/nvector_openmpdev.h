/* -------------------------------------------------------------------
 * Programmer(s): David J. Gardner and Shelby Lockhart @ LLNL
 * -------------------------------------------------------------------
 * Acknowledgements: This NVECTOR module is based on the NVECTOR
 *                   Serial module by Scott D. Cohen, Alan C.
 *                   Hindmarsh, Radu Serban, and Aaron Collier
 *                   @ LLNL
 * -------------------------------------------------------------------
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
 * This is the header file for the OpenMP 4.5+ implementation of the
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

#ifndef SUNDIALS_NVECTOR_OPENMPDEV_H
#define SUNDIALS_NVECTOR_OPENMPDEV_H

#include <stdio.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * OpenMPDEV implementation of N_Vector
 * -----------------------------------------------------------------
 */

struct _N_VectorContent_OpenMPDEV
{
  sunindextype length;     /* vector length       */
  sunbooleantype own_data; /* data ownership flag */
  sunrealtype* host_data;  /* host data array     */
  sunrealtype* dev_data;   /* device data array   */
};

typedef struct _N_VectorContent_OpenMPDEV* N_VectorContent_OpenMPDEV;

/*
 * -----------------------------------------------------------------
 * Macros NV_CONTENT_OMPDEV, NV_DATA_HOST_OMPDEV, NV_OWN_DATA_OMPDEV,
 *        NV_LENGTH_OMPDEV, and NV_Ith_OMPDEV
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_OMPDEV(v) ((N_VectorContent_OpenMPDEV)(v->content))

#define NV_LENGTH_OMPDEV(v) (NV_CONTENT_OMPDEV(v)->length)

#define NV_OWN_DATA_OMPDEV(v) (NV_CONTENT_OMPDEV(v)->own_data)

#define NV_DATA_HOST_OMPDEV(v) (NV_CONTENT_OMPDEV(v)->host_data)

#define NV_DATA_DEV_OMPDEV(v) (NV_CONTENT_OMPDEV(v)->dev_data)

/*
 * -----------------------------------------------------------------
 * Functions exported by nvector_openmpdev
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNew_OpenMPDEV(sunindextype vec_length,
                                          SUNContext sunctx);

SUNDIALS_EXPORT N_Vector N_VNewEmpty_OpenMPDEV(sunindextype vec_length,
                                               SUNContext sunctx);

SUNDIALS_EXPORT N_Vector N_VMake_OpenMPDEV(sunindextype vec_length,
                                           sunrealtype* h_data_1d,
                                           sunrealtype* v_data_1d,
                                           SUNContext sunctx);

SUNDIALS_EXPORT void N_VCopyToDevice_OpenMPDEV(N_Vector v);

SUNDIALS_EXPORT void N_VCopyFromDevice_OpenMPDEV(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNErrCode N_VEnableFusedOps_OpenMPDEV(N_Vector v,
                                                       sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableLinearCombination_OpenMPDEV(N_Vector v,
                                                                sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleAddMulti_OpenMPDEV(N_Vector v,
                                                            sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableDotProdMulti_OpenMPDEV(N_Vector v,
                                                           sunbooleantype tf);

SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearSumVectorArray_OpenMPDEV(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleVectorArray_OpenMPDEV(N_Vector v,
                                                               sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableConstVectorArray_OpenMPDEV(N_Vector v,
                                                               sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormVectorArray_OpenMPDEV(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormMaskVectorArray_OpenMPDEV(N_Vector v,
                                                      sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMultiVectorArray_OpenMPDEV(N_Vector v,
                                                       sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombinationVectorArray_OpenMPDEV(N_Vector v,
                                                           sunbooleantype tf);

#ifdef __cplusplus
}
#endif

#endif
