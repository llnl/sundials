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
 * Temporary header file for deprecated functions.
 * TODO:(SBR) remove in version 8.1.0
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_NVECTOR_NVECTOR_PARHYP_DEPRECATED_H
#define SUNDIALS_NVECTOR_NVECTOR_PARHYP_DEPRECATED_H

#include <nvector/nvector_parhyp.h>

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VGetVectorID instead; will be removed in version 8.1.0")
N_Vector_ID N_VGetVectorID_ParHyp(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VCloneEmpty instead; will be removed in version 8.1.0")
N_Vector N_VCloneEmpty_ParHyp(N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VClone instead; will be removed in version 8.1.0")
N_Vector N_VClone_ParHyp(N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDestroy instead; will be removed in version 8.1.0")
void N_VDestroy_ParHyp(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VGetCommunicator instead; will be removed in version 8.1.0")
MPI_Comm N_VGetCommunicator_ParHyp(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VGetLength instead; will be removed in version 8.1.0")
sunindextype N_VGetLength_ParHyp(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VLinearSum instead; will be removed in version 8.1.0")
void N_VLinearSum_ParHyp(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                         N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VConst instead; will be removed in version 8.1.0")
void N_VConst_ParHyp(sunrealtype c, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VProd instead; will be removed in version 8.1.0")
void N_VProd_ParHyp(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDiv instead; will be removed in version 8.1.0")
void N_VDiv_ParHyp(N_Vector x, N_Vector y, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VScale instead; will be removed in version 8.1.0")
void N_VScale_ParHyp(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VAbs instead; will be removed in version 8.1.0")
void N_VAbs_ParHyp(N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VInv instead; will be removed in version 8.1.0")
void N_VInv_ParHyp(N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VAddConst instead; will be removed in version 8.1.0")
void N_VAddConst_ParHyp(N_Vector x, sunrealtype b, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDotProd instead; will be removed in version 8.1.0")
sunrealtype N_VDotProd_ParHyp(N_Vector x, N_Vector y);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMaxNorm instead; will be removed in version 8.1.0")
sunrealtype N_VMaxNorm_ParHyp(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWrmsNorm instead; will be removed in version 8.1.0")
sunrealtype N_VWrmsNorm_ParHyp(N_Vector x, N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWrmsNormMask instead; will be removed in version 8.1.0")
sunrealtype N_VWrmsNormMask_ParHyp(N_Vector x, N_Vector w, N_Vector id);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMin instead; will be removed in version 8.1.0")
sunrealtype N_VMin_ParHyp(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWL2Norm instead; will be removed in version 8.1.0")
sunrealtype N_VWL2Norm_ParHyp(N_Vector x, N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VL1Norm instead; will be removed in version 8.1.0")
sunrealtype N_VL1Norm_ParHyp(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VCompare instead; will be removed in version 8.1.0")
void N_VCompare_ParHyp(sunrealtype c, N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VInvTest instead; will be removed in version 8.1.0")
sunbooleantype N_VInvTest_ParHyp(N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VConstrMask instead; will be removed in version 8.1.0")
sunbooleantype N_VConstrMask_ParHyp(N_Vector c, N_Vector x, N_Vector m);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMinQuotient instead; will be removed in version 8.1.0")
sunrealtype N_VMinQuotient_ParHyp(N_Vector num, N_Vector denom);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VLinearCombination instead; will be removed in version 8.1.0")
SUNErrCode N_VLinearCombination_ParHyp(int nvec, sunrealtype* c, N_Vector* X,
                                       N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VScaleAddMulti instead; will be removed in version 8.1.0")
SUNErrCode N_VScaleAddMulti_ParHyp(int nvec, sunrealtype* a, N_Vector x,
                                   N_Vector* Y, N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDotProdMulti instead; will be removed in version 8.1.0")
SUNErrCode N_VDotProdMulti_ParHyp(int nvec, N_Vector x, N_Vector* Y,
                                  sunrealtype* dotprods);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VLinearSumVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VLinearSumVectorArray_ParHyp(int nvec, sunrealtype a, N_Vector* X,
                                          sunrealtype b, N_Vector* Y,
                                          N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VScaleVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VScaleVectorArray_ParHyp(int nvec, sunrealtype* c, N_Vector* X,
                                      N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VConstVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VConstVectorArray_ParHyp(int nvecs, sunrealtype c, N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWrmsNormVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VWrmsNormVectorArray_ParHyp(int nvecs, N_Vector* X, N_Vector* W,
                                         sunrealtype* nrm);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWrmsNormMaskVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VWrmsNormMaskVectorArray_ParHyp(int nvec, N_Vector* X, N_Vector* W,
                                             N_Vector id, sunrealtype* nrm);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VScaleAddMultiVectorArray instead; will be removed in version 8.1.0")
SUNErrCode N_VScaleAddMultiVectorArray_ParHyp(int nvec, int nsum,
                                              sunrealtype* a, N_Vector* X,
                                              N_Vector** Y, N_Vector** Z);

SUNDIALS_DEPRECATED_EXPORT_MSG("use N_VLinearCombinationVectorArray instead; "
                               "will be removed in version 8.1.0")
SUNErrCode N_VLinearCombinationVectorArray_ParHyp(int nvec, int nsum,
                                                  sunrealtype* c, N_Vector** X,
                                                  N_Vector* Z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDotProdLocal instead; will be removed in version 8.1.0")
sunrealtype N_VDotProdLocal_ParHyp(N_Vector x, N_Vector y);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMaxNormLocal instead; will be removed in version 8.1.0")
sunrealtype N_VMaxNormLocal_ParHyp(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMinLocal instead; will be removed in version 8.1.0")
sunrealtype N_VMinLocal_ParHyp(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VL1NormLocal instead; will be removed in version 8.1.0")
sunrealtype N_VL1NormLocal_ParHyp(N_Vector x);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWSqrSumLocal instead; will be removed in version 8.1.0")
sunrealtype N_VWSqrSumLocal_ParHyp(N_Vector x, N_Vector w);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VWSqrSumMaskLocal instead; will be removed in version 8.1.0")
sunrealtype N_VWSqrSumMaskLocal_ParHyp(N_Vector x, N_Vector w, N_Vector id);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VInvTestLocal instead; will be removed in version 8.1.0")
sunbooleantype N_VInvTestLocal_ParHyp(N_Vector x, N_Vector z);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VConstrMaskLocal instead; will be removed in version 8.1.0")
sunbooleantype N_VConstrMaskLocal_ParHyp(N_Vector c, N_Vector x, N_Vector m);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VMinQuotientLocal instead; will be removed in version 8.1.0")
sunrealtype N_VMinQuotientLocal_ParHyp(N_Vector num, N_Vector denom);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDotProdMultiLocal instead; will be removed in version 8.1.0")
SUNErrCode N_VDotProdMultiLocal_ParHyp(int nvec, N_Vector x, N_Vector* Y,
                                       sunrealtype* dotprods);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VDotProdMultiAllReduce instead; will be removed in version 8.1.0")
SUNErrCode N_VDotProdMultiAllReduce_ParHyp(int nvec, N_Vector x,
                                           sunrealtype* sum);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VBufSize instead; will be removed in version 8.1.0")
SUNErrCode N_VBufSize_ParHyp(N_Vector x, sunindextype* size);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VBufPack instead; will be removed in version 8.1.0")
SUNErrCode N_VBufPack_ParHyp(N_Vector x, void* buf);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VBufUnpack instead; will be removed in version 8.1.0")
SUNErrCode N_VBufUnpack_ParHyp(N_Vector x, void* buf);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VPrint instead; will be removed in version 8.1.0")
void N_VPrint_ParHyp(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use N_VPrintFile instead; will be removed in version 8.1.0")
void N_VPrintFile_ParHyp(N_Vector v, FILE* outfile);

#ifdef __cplusplus
}
#endif

#endif
