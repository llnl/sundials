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
 *       N_VLinearSum(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_NVECTOR_MANYVECTOR_H
#define SUNDIALS_NVECTOR_MANYVECTOR_H

#include <stdio.h>
#include <sundials/sundials_core.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------------------
   ManyVector implementation of N_Vector
   ----------------------------------------------------------------- */

struct N_VectorContent_ManyVector_
{
  sunindextype num_subvectors; /* number of vectors attached       */
  sunindextype global_length;  /* overall global manyvector length */
  N_Vector* subvec_array;      /* pointer to N_Vector array        */
  sunbooleantype own_data;     /* flag indicating data ownership   */
};

typedef struct N_VectorContent_ManyVector_* N_VectorContent_ManyVector;

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
sunindextype N_VGetSubvectorLocalLength_ManyVector(N_Vector v,
                                                   sunindextype vec_num);

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
