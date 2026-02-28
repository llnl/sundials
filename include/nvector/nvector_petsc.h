/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
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
 * This is the main header file for the PETSc vector wrapper
 * for NVECTOR module.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be
 *     found in the header file sundials_nvector.h.
 *
 *   - The definition of the type sunrealtype can be found in the
 *     header file sundials_types.h, and it may be changed (at the
 *     build configuration stage) according to the user's needs.
 *     The sundials_types.h file also contains the definition
 *     for the type sunbooleantype.
 *
 *   - N_Vector arguments to arithmetic vector operations need not
 *     be distinct. For example, the following call:
 *
 *        N_VLinearSum_Petsc(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_PETSC_H
#define _NVECTOR_PETSC_H

#include <mpi.h>
#include <petscvec.h>
#include <sundials/sundials_mpi_types.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * PETSc implementation of N_Vector
 * -----------------------------------------------------------------
 */

struct _N_VectorContent_Petsc
{
  sunindextype local_length;  /* copy of local vector length  */
  sunindextype global_length; /* copy of global vector length */
  sunbooleantype own_data;    /* ownership of data            */
  Vec pvec;                   /* the PETSc Vec object         */
  MPI_Comm comm;              /* copy of MPI communicator     */
};

typedef struct _N_VectorContent_Petsc* N_VectorContent_Petsc;

/*
 * -----------------------------------------------------------------
 * Functions exported by nvector_petsc
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_Petsc(MPI_Comm comm,
                                           sunindextype local_length,
                                           sunindextype global_length,
                                           SUNContext sunctx);

SUNDIALS_EXPORT N_Vector N_VMake_Petsc(Vec v, SUNContext sunctx);

SUNDIALS_EXPORT Vec N_VGetVector_Petsc(N_Vector v);

SUNDIALS_EXPORT void N_VSetVector_Petsc(N_Vector v, Vec p);

SUNDIALS_EXPORT void N_VPrint_Petsc(N_Vector v);

SUNDIALS_EXPORT void N_VPrintFile_Petsc(N_Vector v, const char fname[]);

/* nvector API functions */
SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID_Petsc(N_Vector v);
SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Petsc(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_Petsc(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_Petsc(N_Vector v);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
void N_VSpace_Petsc(N_Vector v, sunindextype* lrw, sunindextype* liw);
SUNDIALS_EXPORT MPI_Comm N_VGetCommunicator_Petsc(N_Vector v);
SUNDIALS_EXPORT sunindextype N_VGetLength_Petsc(N_Vector v);

/* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum_Petsc(sunrealtype a, N_Vector x,
                                        sunrealtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_Petsc(sunrealtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_Petsc(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_Petsc(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_Petsc(sunrealtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_Petsc(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_Petsc(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_Petsc(N_Vector x, sunrealtype b, N_Vector z);
SUNDIALS_EXPORT sunrealtype N_VDotProd_Petsc(N_Vector x, N_Vector y);
SUNDIALS_EXPORT sunrealtype N_VMaxNorm_Petsc(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWrmsNorm_Petsc(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VWrmsNormMask_Petsc(N_Vector x, N_Vector w,
                                                  N_Vector id);
SUNDIALS_EXPORT sunrealtype N_VMin_Petsc(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWL2Norm_Petsc(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VL1Norm_Petsc(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_Petsc(sunrealtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VInvTest_Petsc(N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VConstrMask_Petsc(N_Vector c, N_Vector x,
                                                   N_Vector m);
SUNDIALS_EXPORT sunrealtype N_VMinQuotient_Petsc(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT SUNErrCode N_VLinearCombination_Petsc(int nvec, sunrealtype* c,
                                                      N_Vector* X, N_Vector z);
SUNDIALS_EXPORT SUNErrCode N_VScaleAddMulti_Petsc(int nvec, sunrealtype* a,
                                                  N_Vector x, N_Vector* Y,
                                                  N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VDotProdMulti_Petsc(int nvec, N_Vector x,
                                                 N_Vector* Y,
                                                 sunrealtype* dotprods);

/* vector array operations */
SUNDIALS_EXPORT SUNErrCode N_VLinearSumVectorArray_Petsc(
  int nvec, sunrealtype a, N_Vector* X, sunrealtype b, N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VScaleVectorArray_Petsc(int nvec, sunrealtype* c,
                                                     N_Vector* X, N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VConstVectorArray_Petsc(int nvecs, sunrealtype c,
                                                     N_Vector* Z);
SUNDIALS_EXPORT SUNErrCode N_VWrmsNormVectorArray_Petsc(int nvecs, N_Vector* X,
                                                        N_Vector* W,
                                                        sunrealtype* nrm);
SUNDIALS_EXPORT SUNErrCode N_VWrmsNormMaskVectorArray_Petsc(
  int nvec, N_Vector* X, N_Vector* W, N_Vector id, sunrealtype* nrm);
SUNDIALS_EXPORT SUNErrCode N_VScaleAddMultiVectorArray_Petsc(
  int nvec, int nsum, sunrealtype* a, N_Vector* X, N_Vector** Y, N_Vector** Z);
SUNDIALS_EXPORT SUNErrCode N_VLinearCombinationVectorArray_Petsc(
  int nvec, int nsum, sunrealtype* c, N_Vector** X, N_Vector* Z);

/* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT sunrealtype N_VDotProdLocal_Petsc(N_Vector x, N_Vector y);
SUNDIALS_EXPORT sunrealtype N_VMaxNormLocal_Petsc(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VMinLocal_Petsc(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VL1NormLocal_Petsc(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWSqrSumLocal_Petsc(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VWSqrSumMaskLocal_Petsc(N_Vector x, N_Vector w,
                                                      N_Vector id);
SUNDIALS_EXPORT sunbooleantype N_VInvTestLocal_Petsc(N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VConstrMaskLocal_Petsc(N_Vector c, N_Vector x,
                                                        N_Vector m);
SUNDIALS_EXPORT sunrealtype N_VMinQuotientLocal_Petsc(N_Vector num,
                                                      N_Vector denom);

/* OPTIONAL single buffer reduction operations */
SUNDIALS_EXPORT SUNErrCode N_VDotProdMultiLocal_Petsc(int nvec, N_Vector x,
                                                      N_Vector* Y,
                                                      sunrealtype* dotprods);
SUNDIALS_EXPORT SUNErrCode N_VDotProdMultiAllReduce_Petsc(int nvec, N_Vector x,
                                                          sunrealtype* sum);

/* OPTIONAL XBraid interface operations */
SUNDIALS_EXPORT SUNErrCode N_VBufSize_Petsc(N_Vector x, sunindextype* size);
SUNDIALS_EXPORT SUNErrCode N_VBufPack_Petsc(N_Vector x, void* buf);
SUNDIALS_EXPORT SUNErrCode N_VBufUnpack_Petsc(N_Vector x, void* buf);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNErrCode N_VEnableFusedOps_Petsc(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableLinearCombination_Petsc(N_Vector v,
                                                            sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleAddMulti_Petsc(N_Vector v,
                                                        sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableDotProdMulti_Petsc(N_Vector v,
                                                       sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableLinearSumVectorArray_Petsc(N_Vector v,
                                                               sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleVectorArray_Petsc(N_Vector v,
                                                           sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableConstVectorArray_Petsc(N_Vector v,
                                                           sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableWrmsNormVectorArray_Petsc(N_Vector v,
                                                              sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormMaskVectorArray_Petsc(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMultiVectorArray_Petsc(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombinationVectorArray_Petsc(N_Vector v,
                                                       sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableDotProdMultiLocal_Petsc(N_Vector v,
                                                            sunbooleantype tf);

SUNDIALS_DEPRECATED_EXPORT_MSG("Not supported, use N_VGetVector_Petsc")
sunrealtype* N_VGetArrayPointer_Petsc(N_Vector v);

SUNDIALS_DEPRECATED_EXPORT_MSG("Not supported, use N_VGetVector_Petsc")
void N_VSetArrayPointer_Petsc(sunrealtype* v_data, N_Vector v);

#ifdef __cplusplus
}
#endif

#endif /* _NVECTOR_PETSC_H */
