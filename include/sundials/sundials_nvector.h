/* -----------------------------------------------------------------
 * Programmer(s): Radu Serban and Aaron Collier @ LLNL
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
 * This is the header file for a generic NVECTOR package.
 * It defines the N_Vector structure (_generic_N_Vector) which
 * contains the following fields:
 *   - an implementation-dependent 'content' field which contains
 *     the description and actual data of the vector
 *   - an 'ops' filed which contains a structure listing operations
 *     acting on such vectors
 * -----------------------------------------------------------------
 * This header file contains:
 *   - enumeration constants for all SUNDIALS-defined vector types,
 *     as well as a generic type for user-supplied vector types,
 *   - type declarations for the _generic_N_Vector and
 *     _generic_N_Vector_Ops structures, as well as references to
 *     pointers to such structures (N_Vector), and
 *   - prototypes for the vector functions which operate on
 *     N_Vector objects.
 * -----------------------------------------------------------------
 * At a minimum, a particular implementation of an NVECTOR must
 * do the following:
 *   - specify the 'content' field of N_Vector,
 *   - implement the operations on those N_Vector objects,
 *   - provide a constructor routine for new N_Vector objects
 *
 * Additionally, an NVECTOR implementation may provide the following:
 *   - macros to access the underlying N_Vector data
 *   - a constructor for an array of N_Vectors
 *   - a constructor for an empty N_Vector (i.e., a new N_Vector with
 *     a NULL data pointer).
 *   - a routine to print the content of an N_Vector
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_H
#define _NVECTOR_H

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------------------
 * Implemented N_Vector types
 * ----------------------------------------------------------------- */

enum N_Vector_ID
{
  SUNDIALS_NVEC_SERIAL,
  SUNDIALS_NVEC_PARALLEL,
  SUNDIALS_NVEC_OPENMP,
  SUNDIALS_NVEC_PTHREADS,
  SUNDIALS_NVEC_PARHYP,
  SUNDIALS_NVEC_PETSC,
  SUNDIALS_NVEC_CUDA,
  SUNDIALS_NVEC_HIP,
  SUNDIALS_NVEC_SYCL,
  SUNDIALS_NVEC_RAJA,
  SUNDIALS_NVEC_KOKKOS,
  SUNDIALS_NVEC_OPENMPDEV,
  SUNDIALS_NVEC_TRILINOS,
  SUNDIALS_NVEC_MANYVECTOR,
  SUNDIALS_NVEC_MPIMANYVECTOR,
  SUNDIALS_NVEC_MPIPLUSX,
  SUNDIALS_NVEC_CUSTOM
};

#ifndef SWIG
typedef enum N_Vector_ID N_Vector_ID;
#endif

/* -----------------------------------------------------------------
 * Generic definition of N_Vector
 * ----------------------------------------------------------------- */

/* Forward reference for pointer to N_Vector_Ops object */
typedef _SUNDIALS_STRUCT_ _generic_N_Vector_Ops* N_Vector_Ops;

/* Forward reference for pointer to N_Vector object */
typedef _SUNDIALS_STRUCT_ _generic_N_Vector* N_Vector;

/* Define array of N_Vectors */
/* DEPRECATION NOTICE: N_Vector_S will be removed in the next major release */
typedef N_Vector* N_Vector_S;

/* Structure containing function pointers to vector operations  */
struct _generic_N_Vector_Ops
{
  /*
   * REQUIRED operations.
   *
   * These must be implemented by derivations of the generic N_Vector.
   */

  /* constructors, destructors, and utility operations */
  N_Vector_ID (*nvgetvectorid)(N_Vector);
  N_Vector (*nvclone)(N_Vector);
  N_Vector (*nvcloneempty)(N_Vector);
  void (*nvdestroy)(N_Vector);
  void (*nvspace)(N_Vector, sunindextype*, sunindextype*);
  sunrealtype* (*nvgetarraypointer)(N_Vector);
  sunrealtype* (*nvgetdevicearraypointer)(N_Vector);
  void (*nvsetarraypointer)(sunrealtype*, N_Vector);
  SUNComm (*nvgetcommunicator)(N_Vector);
  sunindextype (*nvgetlength)(N_Vector);
  sunindextype (*nvgetlocallength)(N_Vector);

  /* standard vector operations */
  void (*nvlinearsum)(sunrealtype, N_Vector, sunrealtype, N_Vector, N_Vector);
  void (*nvconst)(sunrealtype, N_Vector);
  void (*nvprod)(N_Vector, N_Vector, N_Vector);
  void (*nvdiv)(N_Vector, N_Vector, N_Vector);
  void (*nvscale)(sunrealtype, N_Vector, N_Vector);
  void (*nvabs)(N_Vector, N_Vector);
  void (*nvinv)(N_Vector, N_Vector);
  void (*nvaddconst)(N_Vector, sunrealtype, N_Vector);
  sunrealtype (*nvdotprod)(N_Vector, N_Vector);
  sunrealtype (*nvmaxnorm)(N_Vector);
  sunrealtype (*nvwrmsnorm)(N_Vector, N_Vector);
  sunrealtype (*nvwrmsnormmask)(N_Vector, N_Vector, N_Vector);
  sunrealtype (*nvmin)(N_Vector);
  sunrealtype (*nvwl2norm)(N_Vector, N_Vector);
  sunrealtype (*nvl1norm)(N_Vector);
  void (*nvcompare)(sunrealtype, N_Vector, N_Vector);
  sunbooleantype (*nvinvtest)(N_Vector, N_Vector);
  sunbooleantype (*nvconstrmask)(N_Vector, N_Vector, N_Vector);
  sunrealtype (*nvminquotient)(N_Vector, N_Vector);

  /*
   * OPTIONAL operations.
   *
   * These operations provide default implementations that may be overridden.
   */

  /* OPTIONAL fused vector operations */
  SUNErrCode (*nvlinearcombination)(int, sunrealtype*, N_Vector*, N_Vector);
  SUNErrCode (*nvscaleaddmulti)(int, sunrealtype*, N_Vector, N_Vector*,
                                N_Vector*);
  SUNErrCode (*nvdotprodmulti)(int, N_Vector, N_Vector*, sunrealtype*);

  /* OPTIONAL vector array operations */
  SUNErrCode (*nvlinearsumvectorarray)(int, sunrealtype, N_Vector*, sunrealtype,
                                       N_Vector*, N_Vector*);
  SUNErrCode (*nvscalevectorarray)(int, sunrealtype*, N_Vector*, N_Vector*);
  SUNErrCode (*nvconstvectorarray)(int, sunrealtype, N_Vector*);
  SUNErrCode (*nvwrmsnormvectorarray)(int, N_Vector*, N_Vector*, sunrealtype*);
  SUNErrCode (*nvwrmsnormmaskvectorarray)(int, N_Vector*, N_Vector*, N_Vector,
                                          sunrealtype*);
  SUNErrCode (*nvscaleaddmultivectorarray)(int, int, sunrealtype*, N_Vector*,
                                           N_Vector**, N_Vector**);
  SUNErrCode (*nvlinearcombinationvectorarray)(int, int, sunrealtype*,
                                               N_Vector**, N_Vector*);

  /*
   * OPTIONAL operations with no default implementation.
   */

  /* Local reduction kernels (no parallel communication) */
  sunrealtype (*nvdotprodlocal)(N_Vector, N_Vector);
  sunrealtype (*nvmaxnormlocal)(N_Vector);
  sunrealtype (*nvminlocal)(N_Vector);
  sunrealtype (*nvl1normlocal)(N_Vector);
  sunbooleantype (*nvinvtestlocal)(N_Vector, N_Vector);
  sunbooleantype (*nvconstrmasklocal)(N_Vector, N_Vector, N_Vector);
  sunrealtype (*nvminquotientlocal)(N_Vector, N_Vector);
  sunrealtype (*nvwsqrsumlocal)(N_Vector, N_Vector);
  sunrealtype (*nvwsqrsummasklocal)(N_Vector, N_Vector, N_Vector);

  /* Single buffer reduction operations */
  SUNErrCode (*nvdotprodmultilocal)(int, N_Vector, N_Vector*, sunrealtype*);
  SUNErrCode (*nvdotprodmultiallreduce)(int, N_Vector, sunrealtype*);

  /* XBraid interface operations */
  SUNErrCode (*nvbufsize)(N_Vector, sunindextype*);
  SUNErrCode (*nvbufpack)(N_Vector, void*);
  SUNErrCode (*nvbufunpack)(N_Vector, void*);

  /* Debugging functions (called when SUNDIALS_DEBUG_PRINTVEC is defined). */
  void (*nvprint)(N_Vector);
  void (*nvprintfile)(N_Vector, FILE*);
};

/* A vector is a structure with an implementation-dependent
   'content' field, and a pointer to a structure of vector
   operations corresponding to that implementation. */
struct _generic_N_Vector
{
  void* content;
  N_Vector_Ops ops;
  SUNContext sunctx;
};

/* -----------------------------------------------------------------
 * Functions exported by NVECTOR module
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT N_Vector N_VNewEmpty(SUNContext sunctx);
SUNDIALS_EXPORT void N_VFreeEmpty(N_Vector v);
SUNDIALS_EXPORT SUNErrCode N_VCopyOps(N_Vector w, N_Vector v);

/*
 * Required operations.
 */

SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VCloneEmpty(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy(N_Vector v);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
void N_VSpace(N_Vector v, sunindextype* lrw, sunindextype* liw);
SUNDIALS_EXPORT sunrealtype* N_VGetArrayPointer(N_Vector v);
SUNDIALS_EXPORT sunrealtype* N_VGetDeviceArrayPointer(N_Vector v);
SUNDIALS_EXPORT void N_VSetArrayPointer(sunrealtype* v_data_1d, N_Vector v);
SUNDIALS_EXPORT SUNComm N_VGetCommunicator(N_Vector v);
SUNDIALS_EXPORT sunindextype N_VGetLength(N_Vector v);
SUNDIALS_EXPORT sunindextype N_VGetLocalLength(N_Vector v);

/* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum(sunrealtype a, N_Vector x, sunrealtype b,
                                  N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst(sunrealtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale(sunrealtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst(N_Vector x, sunrealtype b, N_Vector z);
SUNDIALS_EXPORT sunrealtype N_VDotProd(N_Vector x, N_Vector y);
SUNDIALS_EXPORT sunrealtype N_VMaxNorm(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWrmsNorm(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VWrmsNormMask(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT sunrealtype N_VMin(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWL2Norm(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VL1Norm(N_Vector x);
SUNDIALS_EXPORT void N_VCompare(sunrealtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VInvTest(N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VConstrMask(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT sunrealtype N_VMinQuotient(N_Vector num, N_Vector denom);

/*
 * OPTIONAL operations with default implementations.
 */

/* fused vector operations */
SUNDIALS_EXPORT
SUNErrCode N_VLinearCombination(int nvec, sunrealtype* c_1d, N_Vector* X_1d,
                                N_Vector z);

SUNDIALS_EXPORT
SUNErrCode N_VScaleAddMulti(int nvec, sunrealtype* a_1d, N_Vector x,
                            N_Vector* Y_1d, N_Vector* Z_1d);

SUNDIALS_EXPORT
SUNErrCode N_VDotProdMulti(int nvec, N_Vector x, N_Vector* Y_1d,
                           sunrealtype* dotprods_1d);

/* vector array operations */
SUNDIALS_EXPORT
SUNErrCode N_VLinearSumVectorArray(int nvec, sunrealtype a, N_Vector* X_1d,
                                   sunrealtype b, N_Vector* Y_1d, N_Vector* Z_1d);

SUNDIALS_EXPORT
SUNErrCode N_VScaleVectorArray(int nvec, sunrealtype* c_1d, N_Vector* X_1d,
                               N_Vector* Z_1d);

SUNDIALS_EXPORT
SUNErrCode N_VConstVectorArray(int nvec, sunrealtype c, N_Vector* Z_1d);

SUNDIALS_EXPORT
SUNErrCode N_VWrmsNormVectorArray(int nvec, N_Vector* X_1d, N_Vector* W_1d,
                                  sunrealtype* nrm_1d);

SUNDIALS_EXPORT
SUNErrCode N_VWrmsNormMaskVectorArray(int nvec, N_Vector* X_1d, N_Vector* W_1d,
                                      N_Vector id, sunrealtype* nrm_1d);

SUNDIALS_EXPORT
SUNErrCode N_VScaleAddMultiVectorArray(int nvec, int nsum, sunrealtype* a_1d,
                                       N_Vector* X_1d, N_Vector** Y_2d,
                                       N_Vector** Z_2d);

SUNDIALS_EXPORT
SUNErrCode N_VLinearCombinationVectorArray(int nvec, int nsum, sunrealtype* c_1d,
                                           N_Vector** X_2d, N_Vector* Z_1d);

/*
 * OPTIONAL operations with no default implementation.
 */

/* local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT sunrealtype N_VDotProdLocal(N_Vector x, N_Vector y);
SUNDIALS_EXPORT sunrealtype N_VMaxNormLocal(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VMinLocal(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VL1NormLocal(N_Vector x);
SUNDIALS_EXPORT sunrealtype N_VWSqrSumLocal(N_Vector x, N_Vector w);
SUNDIALS_EXPORT sunrealtype N_VWSqrSumMaskLocal(N_Vector x, N_Vector w,
                                                N_Vector id);
SUNDIALS_EXPORT sunbooleantype N_VInvTestLocal(N_Vector x, N_Vector z);
SUNDIALS_EXPORT sunbooleantype N_VConstrMaskLocal(N_Vector c, N_Vector x,
                                                  N_Vector m);
SUNDIALS_EXPORT sunrealtype N_VMinQuotientLocal(N_Vector num, N_Vector denom);

/* single buffer reduction operations */
SUNDIALS_EXPORT SUNErrCode N_VDotProdMultiLocal(int nvec, N_Vector x,
                                                N_Vector* Y_1d,
                                                sunrealtype* dotprods_1d);
SUNDIALS_EXPORT SUNErrCode N_VDotProdMultiAllReduce(int nvec_total, N_Vector x,
                                                    sunrealtype* sum_1d);

/* XBraid interface operations */
SUNDIALS_EXPORT SUNErrCode N_VBufSize(N_Vector x, sunindextype* size);
SUNDIALS_EXPORT SUNErrCode N_VBufPack(N_Vector x, void* buf);
SUNDIALS_EXPORT SUNErrCode N_VBufUnpack(N_Vector x, void* buf);

/* -----------------------------------------------------------------
 * Additional functions exported by NVECTOR module
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT N_Vector* N_VNewVectorArray(int count, SUNContext sunctx);
SUNDIALS_EXPORT N_Vector* N_VCloneEmptyVectorArray(int count, N_Vector w);
SUNDIALS_EXPORT N_Vector* N_VCloneVectorArray(int count, N_Vector w);
SUNDIALS_EXPORT void N_VDestroyVectorArray(N_Vector* vs_1d, int count);

/* These function are really only for users of the Fortran interface */
SUNDIALS_EXPORT N_Vector N_VGetVecAtIndexVectorArray(N_Vector* vs_1d, int index);
SUNDIALS_EXPORT void N_VSetVecAtIndexVectorArray(N_Vector* vs_1d, int index,
                                                 N_Vector w);

/* -----------------------------------------------------------------
 * Debugging functions
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT void N_VPrint(N_Vector v);
SUNDIALS_EXPORT void N_VPrintFile(N_Vector v, FILE* outfile);

#ifdef __cplusplus
}
#endif

#endif
