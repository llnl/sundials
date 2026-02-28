/* -----------------------------------------------------------------
 * Programmer(s): Jean M. Sexton @ UMBC
 *                Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * Based on work by: Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                   and Aaron Collier @ LLNL
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
 * This is the main header file for the ParHyp implementation
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
 *        N_VLinearSum_ParHyp(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_PARHYP_H
#define _NVECTOR_PARHYP_H

#include <mpi.h>
#include <stdio.h>
#include <sundials/sundials_mpi_types.h>
#include <sundials/sundials_nvector.h>

/* hypre header files */
#include <_hypre_parcsr_mv.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * ParHyp implementation of N_Vector
 * -----------------------------------------------------------------
 */

struct _N_VectorContent_ParHyp
{
  sunindextype local_length;    /* local vector length         */
  sunindextype global_length;   /* global vector length        */
  sunbooleantype own_parvector; /* ownership of HYPRE vector   */
  MPI_Comm comm;                /* pointer to MPI communicator */

  HYPRE_ParVector x; /* the actual HYPRE_ParVector object */
};

typedef struct _N_VectorContent_ParHyp* N_VectorContent_ParHyp;

/*
 * -----------------------------------------------------------------
 * Functions exported by nvector_parhyp
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_ParHyp(MPI_Comm comm,
                                            sunindextype local_length,
                                            sunindextype global_length,
                                            SUNContext sunctx);

SUNDIALS_EXPORT N_Vector N_VMake_ParHyp(HYPRE_ParVector x, SUNContext sunctx);

SUNDIALS_EXPORT HYPRE_ParVector N_VGetVector_ParHyp(N_Vector v);

SUNDIALS_EXPORT void N_VPrint_ParHyp(N_Vector v);

SUNDIALS_EXPORT void N_VPrintFile_ParHyp(N_Vector v, FILE* outfile);

SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID_ParHyp(N_Vector v);
SUNDIALS_EXPORT N_Vector N_VCloneEmpty_ParHyp(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_ParHyp(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_ParHyp(N_Vector v);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
void N_VSpace_ParHyp(N_Vector v, sunindextype* lrw, sunindextype* liw);
SUNDIALS_EXPORT MPI_Comm N_VGetCommunicator_ParHyp(N_Vector v);
SUNDIALS_EXPORT sunindextype N_VGetLength_ParHyp(N_Vector v);

/* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum_ParHyp(sunrealtype a, N_Vector x,
                                         sunrealtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_ParHyp(sunrealtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_ParHyp(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_ParHyp(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_ParHyp(sunrealtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_ParHyp(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_ParHyp(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_ParHyp(N_Vector x, sunrealtype b, N_Vector z);
SUNDIALS_EXPORT sunrealtype N_VDotProd_ParHyp(N_Vector x, N_Vector y);
SUNDIALS_EXPORT sunrealtype N_VMaxNorm_ParHyp(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWrmsNorm_ParHyp(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VWrmsNormMask_ParHyp(N_Vector x, N_Vector w,
                                                   N_Vector id);
SUNDIALS_EXPORT sunrealtype N_VMin_ParHyp(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWL2Norm_ParHyp(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VL1Norm_ParHyp(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_ParHyp(sunrealtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VInvTest_ParHyp(N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VConstrMask_ParHyp(N_Vector c, N_Vector x,
                                                    N_Vector m);
SUNDIALS_EXPORT sunrealtype N_VMinQuotient_ParHyp(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT SUNErrCode N_VLinearCombination_ParHyp(int nvec, sunrealtype* c,
                                                       N_Vector* X, N_Vector z);
SUNDIALS_EXPORT SUNErrCode N_VScaleAddMulti_ParHyp(int nvec, sunrealtype* a,
                                                   N_Vector x, N_Vector* Y,
                                                   N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VDotProdMulti_ParHyp(int nvec, N_Vector x,
                                                  N_Vector* Y,
                                                  sunrealtype* dotprods);

/* vector array operations */
SUNDIALS_EXPORT SUNErrCode N_VLinearSumVectorArray_ParHyp(
  int nvec, sunrealtype a, N_Vector* X, sunrealtype b, N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VScaleVectorArray_ParHyp(int nvec, sunrealtype* c,
                                                      N_Vector* X, N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VConstVectorArray_ParHyp(int nvecs, sunrealtype c,
                                                      N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VWrmsNormVectorArray_ParHyp(int nvecs, N_Vector* X,
                                                         N_Vector* W,
                                                         sunrealtype* nrm);
SUNDIALS_EXPORT SUNErrCode N_VWrmsNormMaskVectorArray_ParHyp(
  int nvec, N_Vector* X, N_Vector* W, N_Vector id, sunrealtype* nrm);
SUNDIALS_EXPORT SUNErrCode N_VScaleAddMultiVectorArray_ParHyp(
  int nvec, int nsum, sunrealtype* a, N_Vector* X, N_Vector** Y, N_Vector** Z);
SUNDIALS_EXPORT SUNErrCode N_VLinearCombinationVectorArray_ParHyp(
  int nvec, int nsum, sunrealtype* c, N_Vector** X, N_Vector* Z);

/* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT sunrealtype N_VDotProdLocal_ParHyp(N_Vector x, N_Vector y);
SUNDIALS_EXPORT sunrealtype N_VMaxNormLocal_ParHyp(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VMinLocal_ParHyp(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VL1NormLocal_ParHyp(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWSqrSumLocal_ParHyp(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VWSqrSumMaskLocal_ParHyp(N_Vector x, N_Vector w,
                                                       N_Vector id);
SUNDIALS_EXPORT sunbooleantype N_VInvTestLocal_ParHyp(N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VConstrMaskLocal_ParHyp(N_Vector c, N_Vector x,
                                                         N_Vector m);
SUNDIALS_EXPORT sunrealtype N_VMinQuotientLocal_ParHyp(N_Vector num,
                                                       N_Vector denom);

/* OPTIONAL single buffer reduction operations */
SUNDIALS_EXPORT SUNErrCode N_VDotProdMultiLocal_ParHyp(int nvec, N_Vector x,
                                                       N_Vector* Y,
                                                       sunrealtype* dotprods);
SUNDIALS_EXPORT SUNErrCode N_VDotProdMultiAllReduce_ParHyp(int nvec, N_Vector x,
                                                           sunrealtype* sum);

/* OPTIONAL XBraid interface operations */
SUNDIALS_EXPORT SUNErrCode N_VBufSize_ParHyp(N_Vector x, sunindextype* size);
SUNDIALS_EXPORT SUNErrCode N_VBufPack_ParHyp(N_Vector x, void* buf);
SUNDIALS_EXPORT SUNErrCode N_VBufUnpack_ParHyp(N_Vector x, void* buf);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNErrCode N_VEnableFusedOps_ParHyp(N_Vector v,
                                                    sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableLinearCombination_ParHyp(N_Vector v,
                                                             sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleAddMulti_ParHyp(N_Vector v,
                                                         sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableDotProdMulti_ParHyp(N_Vector v,
                                                        sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableLinearSumVectorArray_ParHyp(N_Vector v,
                                                                sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleVectorArray_ParHyp(N_Vector v,
                                                            sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableConstVectorArray_ParHyp(N_Vector v,
                                                            sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableWrmsNormVectorArray_ParHyp(N_Vector v,
                                                               sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormMaskVectorArray_ParHyp(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMultiVectorArray_ParHyp(N_Vector v,
                                                    sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombinationVectorArray_ParHyp(N_Vector v,
                                                        sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableDotProdMultiLocal_ParHyp(N_Vector v,
                                                             sunbooleantype tf);

SUNDIALS_DEPRECATED_EXPORT_MSG("Not supported, use N_VGetVector_ParHyp")
sunrealtype* N_VGetArrayPointer_ParHyp(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG("Not supported, use N_VGetVector_ParHyp")
void N_VSetArrayPointer_ParHyp(sunrealtype* v_data, N_Vector v);

#ifdef __cplusplus
}
#endif

#endif
