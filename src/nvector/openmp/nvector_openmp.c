/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner and Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * Acknowledgements: This NVECTOR module is based on the NVECTOR
 *                   Serial module by Scott D. Cohen, Alan C.
 *                   Hindmarsh, Radu Serban, and Aaron Collier
 *                   @ LLNL
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
 * This is the implementation file for an OpenMP implementation
 * of the NVECTOR module.
 * -----------------------------------------------------------------*/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_openmp_deprecated.h>
#include <sundials/priv/sundials_context_impl.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_errors.h>

#include "sundials_macros.h"

#define ZERO   SUN_RCONST(0.0)
#define HALF   SUN_RCONST(0.5)
#define ONE    SUN_RCONST(1.0)
#define ONEPT5 SUN_RCONST(1.5)

/* Functions attached to the N_Vector */
static void nvAbs_OpenMP(N_Vector x, N_Vector z);
static void nvAddConst_OpenMP(N_Vector x, sunrealtype b, N_Vector z);
static SUNErrCode nvBufPack_OpenMP(N_Vector x, void* buf);
static SUNErrCode nvBufSize_OpenMP(N_Vector x, sunindextype* size);
static SUNErrCode nvBufUnpack_OpenMP(N_Vector x, void* buf);
static N_Vector nvCloneEmpty_OpenMP(N_Vector w);
static N_Vector nvClone_OpenMP(N_Vector w);
static void nvCompare_OpenMP(sunrealtype c, N_Vector x, N_Vector z);
static SUNErrCode nvConstVectorArray_OpenMP(int nvecs, sunrealtype c,
                                            N_Vector* Z);
static void nvConst_OpenMP(sunrealtype c, N_Vector z);
static sunbooleantype nvConstrMask_OpenMP(N_Vector c, N_Vector x, N_Vector m);
static void nvDestroy_OpenMP(N_Vector v);
static void nvDiv_OpenMP(N_Vector x, N_Vector y, N_Vector z);
static SUNErrCode nvDotProdMulti_OpenMP(int nvec, N_Vector x, N_Vector* Y,
                                        sunrealtype* dotprods);
static sunrealtype nvDotProd_OpenMP(N_Vector x, N_Vector y);
static sunrealtype* nvGetArrayPointer_OpenMP(N_Vector v);
static sunindextype nvGetLength_OpenMP(N_Vector v);
static N_Vector_ID nvGetVectorID_OpenMP(N_Vector v);
static sunbooleantype nvInvTest_OpenMP(N_Vector x, N_Vector z);
static void nvInv_OpenMP(N_Vector x, N_Vector z);
static sunrealtype nvL1Norm_OpenMP(N_Vector x);
static SUNErrCode nvLinearCombinationVectorArray_OpenMP(int nvec, int nsum,
                                                        sunrealtype* c,
                                                        N_Vector** X,
                                                        N_Vector* Z);
static SUNErrCode nvLinearCombination_OpenMP(int nvec, sunrealtype* c,
                                             N_Vector* V, N_Vector z);
static SUNErrCode nvLinearSumVectorArray_OpenMP(int nvec, sunrealtype a,
                                                N_Vector* X, sunrealtype b,
                                                N_Vector* Y, N_Vector* Z);
static void nvLinearSum_OpenMP(sunrealtype a, N_Vector x, sunrealtype b,
                               N_Vector y, N_Vector z);
static sunrealtype nvMaxNorm_OpenMP(N_Vector x);
static sunrealtype nvMinQuotient_OpenMP(N_Vector num, N_Vector denom);
static sunrealtype nvMin_OpenMP(N_Vector x);
static void nvPrintFile_OpenMP(N_Vector v, FILE* outfile);
static void nvPrint_OpenMP(N_Vector v);
static void nvProd_OpenMP(N_Vector x, N_Vector y, N_Vector z);
static SUNErrCode nvScaleAddMultiVectorArray_OpenMP(int nvec, int nsum,
                                                    sunrealtype* a, N_Vector* X,
                                                    N_Vector** Y, N_Vector** Z);
static SUNErrCode nvScaleAddMulti_OpenMP(int nvec, sunrealtype* a, N_Vector x,
                                         N_Vector* Y, N_Vector* Z);
static SUNErrCode nvScaleVectorArray_OpenMP(int nvec, sunrealtype* c,
                                            N_Vector* X, N_Vector* Z);
static void nvScale_OpenMP(sunrealtype c, N_Vector x, N_Vector z);
static void nvSetArrayPointer_OpenMP(sunrealtype* v_data, N_Vector v);
static sunrealtype nvWL2Norm_OpenMP(N_Vector x, N_Vector w);
static sunrealtype nvWSqrSumLocal_OpenMP(N_Vector x, N_Vector w);
static sunrealtype nvWSqrSumMaskLocal_OpenMP(N_Vector x, N_Vector w, N_Vector id);
static SUNErrCode nvWrmsNormMaskVectorArray_OpenMP(int nvecs, N_Vector* X,
                                                   N_Vector* W, N_Vector id,
                                                   sunrealtype* nrm);
static sunrealtype nvWrmsNormMask_OpenMP(N_Vector x, N_Vector w, N_Vector id);
static SUNErrCode nvWrmsNormVectorArray_OpenMP(int nvecs, N_Vector* X,
                                               N_Vector* W, sunrealtype* nrm);
static sunrealtype nvWrmsNorm_OpenMP(N_Vector x, N_Vector w);

/* Private functions for special cases of vector operations */
static void VCopy_OpenMP(N_Vector x, N_Vector z);             /* z=x */
static void VSum_OpenMP(N_Vector x, N_Vector y, N_Vector z);  /* z=x+y     */
static void VDiff_OpenMP(N_Vector x, N_Vector y, N_Vector z); /* z=x-y     */
static void VNeg_OpenMP(N_Vector x, N_Vector z);              /* z=-x */
static void VScaleSum_OpenMP(sunrealtype c, N_Vector x, N_Vector y,
                             N_Vector z); /* z=c(x+y)  */
static void VScaleDiff_OpenMP(sunrealtype c, N_Vector x, N_Vector y,
                              N_Vector z); /* z=c(x-y)  */
static void VLin1_OpenMP(sunrealtype a, N_Vector x, N_Vector y,
                         N_Vector z); /* z=ax+y    */
static void VLin2_OpenMP(sunrealtype a, N_Vector x, N_Vector y,
                         N_Vector z);                            /* z=ax-y    */
static void Vaxpy_OpenMP(sunrealtype a, N_Vector x, N_Vector y); /* y <- ax+y */
static void VScaleBy_OpenMP(sunrealtype a, N_Vector x);          /* x <-
                                                                          ax */

/* Private functions for special cases of vector array operations */
static void VSumVectorArray_OpenMP(int nvec, N_Vector* X, N_Vector* Y,
                                   N_Vector* Z); /* Z=X+Y */
static void VDiffVectorArray_OpenMP(int nvec, N_Vector* X, N_Vector* Y,
                                    N_Vector* Z); /* Z=X-Y */
static void VScaleSumVectorArray_OpenMP(int nvec, sunrealtype c, N_Vector* X,
                                        N_Vector* Y, N_Vector* Z); /* Z=c(X+Y)
                                                                         */
static void VScaleDiffVectorArray_OpenMP(int nvec, sunrealtype c, N_Vector* X,
                                         N_Vector* Y, N_Vector* Z); /* Z=c(X-Y)
                                                                          */
static void VLin1VectorArray_OpenMP(int nvec, sunrealtype a, N_Vector* X,
                                    N_Vector* Y, N_Vector* Z); /* Z=aX+Y */
static void VLin2VectorArray_OpenMP(int nvec, sunrealtype a, N_Vector* X,
                                    N_Vector* Y, N_Vector* Z); /* Z=aX-Y */
static void VaxpyVectorArray_OpenMP(int nvec, sunrealtype a, N_Vector* X,
                                    N_Vector* Y); /* Y <- aX+Y
                                                                     */

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID nvGetVectorID_OpenMP(SUNDIALS_MAYBE_UNUSED N_Vector v)
{
  return SUNDIALS_NVEC_OPENMP;
}

/* ----------------------------------------------------------------------------
 * Function to create a new empty vector
 */

N_Vector N_VNewEmpty_OpenMP(sunindextype length, int num_threads,
                            SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  N_Vector v;
  N_VectorContent_OpenMP content;

  SUNAssertNull(length >= 0, SUN_ERR_ARG_OUTOFRANGE);

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = nvGetVectorID_OpenMP;
  v->ops->nvclone           = nvClone_OpenMP;
  v->ops->nvcloneempty      = nvCloneEmpty_OpenMP;
  v->ops->nvdestroy         = nvDestroy_OpenMP;
  v->ops->nvgetarraypointer = nvGetArrayPointer_OpenMP;
  v->ops->nvsetarraypointer = nvSetArrayPointer_OpenMP;
  v->ops->nvgetlength       = nvGetLength_OpenMP;
  v->ops->nvgetlocallength  = nvGetLength_OpenMP;

  /* standard vector operations */
  v->ops->nvlinearsum    = nvLinearSum_OpenMP;
  v->ops->nvconst        = nvConst_OpenMP;
  v->ops->nvprod         = nvProd_OpenMP;
  v->ops->nvdiv          = nvDiv_OpenMP;
  v->ops->nvscale        = nvScale_OpenMP;
  v->ops->nvabs          = nvAbs_OpenMP;
  v->ops->nvinv          = nvInv_OpenMP;
  v->ops->nvaddconst     = nvAddConst_OpenMP;
  v->ops->nvdotprod      = nvDotProd_OpenMP;
  v->ops->nvmaxnorm      = nvMaxNorm_OpenMP;
  v->ops->nvwrmsnormmask = nvWrmsNormMask_OpenMP;
  v->ops->nvwrmsnorm     = nvWrmsNorm_OpenMP;
  v->ops->nvmin          = nvMin_OpenMP;
  v->ops->nvwl2norm      = nvWL2Norm_OpenMP;
  v->ops->nvl1norm       = nvL1Norm_OpenMP;
  v->ops->nvcompare      = nvCompare_OpenMP;
  v->ops->nvinvtest      = nvInvTest_OpenMP;
  v->ops->nvconstrmask   = nvConstrMask_OpenMP;
  v->ops->nvminquotient  = nvMinQuotient_OpenMP;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction kernels */
  v->ops->nvdotprodlocal     = nvDotProd_OpenMP;
  v->ops->nvmaxnormlocal     = nvMaxNorm_OpenMP;
  v->ops->nvminlocal         = nvMin_OpenMP;
  v->ops->nvl1normlocal      = nvL1Norm_OpenMP;
  v->ops->nvinvtestlocal     = nvInvTest_OpenMP;
  v->ops->nvconstrmasklocal  = nvConstrMask_OpenMP;
  v->ops->nvminquotientlocal = nvMinQuotient_OpenMP;
  v->ops->nvwsqrsumlocal     = nvWSqrSumLocal_OpenMP;
  v->ops->nvwsqrsummasklocal = nvWSqrSumMaskLocal_OpenMP;

  /* single buffer reduction operations */
  v->ops->nvdotprodmultilocal = nvDotProdMulti_OpenMP;

  /* XBraid interface operations */
  v->ops->nvbufsize   = nvBufSize_OpenMP;
  v->ops->nvbufpack   = nvBufPack_OpenMP;
  v->ops->nvbufunpack = nvBufUnpack_OpenMP;

  /* debugging functions */
  v->ops->nvprint     = nvPrint_OpenMP;
  v->ops->nvprintfile = nvPrintFile_OpenMP;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_OpenMP)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  v->content = content;

  /* Initialize content */
  content->length      = length;
  content->num_threads = num_threads;
  content->own_data    = SUNFALSE;
  content->data        = NULL;

  return (v);
}

/* ----------------------------------------------------------------------------
 * Function to create a new vector
 */

N_Vector N_VNew_OpenMP(sunindextype length, int num_threads, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  N_Vector v;
  sunrealtype* data;

  SUNAssertNull(length >= 0, SUN_ERR_ARG_OUTOFRANGE);

  v = NULL;
  v = N_VNewEmpty_OpenMP(length, num_threads, sunctx);
  SUNCheckLastErrNull();

  /* Create data */
  data = NULL;
  if (length > 0)
  {
    data = (sunrealtype*)malloc(length * sizeof(sunrealtype));
    SUNAssertNull(data, SUN_ERR_MALLOC_FAIL);

    /* Attach data */
    NV_OWN_DATA_OMP(v) = SUNTRUE;
    NV_DATA_OMP(v)     = data;
  }

  return (v);
}

/* ----------------------------------------------------------------------------
 * Function to create a vector with user data component
 */

N_Vector N_VMake_OpenMP(sunindextype length, sunrealtype* v_data,
                        int num_threads, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  N_Vector v;

  SUNAssertNull(length >= 0, SUN_ERR_ARG_OUTOFRANGE);

  v = NULL;
  v = N_VNewEmpty_OpenMP(length, num_threads, sunctx);
  SUNCheckLastErrNull();

  if (length > 0)
  {
    /* Attach data */
    NV_OWN_DATA_OMP(v) = SUNFALSE;
    NV_DATA_OMP(v)     = v_data;
  }

  return (v);
}

/* ----------------------------------------------------------------------------
 * Function to return number of vector elements
 */
sunindextype nvGetLength_OpenMP(N_Vector v) { return NV_LENGTH_OMP(v); }

/* ----------------------------------------------------------------------------
 * Function to print a vector to stdout
 */

void nvPrint_OpenMP(N_Vector x) { nvPrintFile_OpenMP(x, stdout); }

/* ----------------------------------------------------------------------------
 * Function to print a vector to outfile
 */

void nvPrintFile_OpenMP(N_Vector x, FILE* outfile)
{
  sunindextype i, N;
  sunrealtype* xd;

  xd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);

  for (i = 0; i < N; i++) { fprintf(outfile, SUN_FORMAT_E "\n", xd[i]); }

  return;
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Create new vector from existing vector without attaching data
 */

N_Vector nvCloneEmpty_OpenMP(N_Vector w)
{
  N_Vector v;
  N_VectorContent_OpenMP content;

  SUNFunctionBegin(w->sunctx);

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty(w->sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  SUNCheckCallNull(N_VCopyOps(w, v));

  /* Create content */
  content = NULL;
  content = (N_VectorContent_OpenMP)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  v->content = content;

  /* Initialize content */
  content->length      = NV_LENGTH_OMP(w);
  content->num_threads = NV_NUM_THREADS_OMP(w);
  content->own_data    = SUNFALSE;
  content->data        = NULL;

  return (v);
}

/* ----------------------------------------------------------------------------
 * Create new vector from existing vector and attach data
 */

N_Vector nvClone_OpenMP(N_Vector w)
{
  SUNFunctionBegin(w->sunctx);
  N_Vector v;
  sunrealtype* data;
  sunindextype length;

  v = NULL;
  v = nvCloneEmpty_OpenMP(w);
  SUNCheckLastErrNull();

  length = NV_LENGTH_OMP(w);

  /* Create data */
  data = NULL;
  if (length > 0)
  {
    data = (sunrealtype*)malloc(length * sizeof(sunrealtype));
    SUNAssertNull(data, SUN_ERR_MALLOC_FAIL);
  }

  /* Attach data */
  NV_OWN_DATA_OMP(v) = SUNTRUE;
  NV_DATA_OMP(v)     = data;

  return (v);
}

/* ----------------------------------------------------------------------------
 * Destroy vector and free vector memory
 */

void nvDestroy_OpenMP(N_Vector v)
{
  if (v == NULL) { return; }

  /* free content */
  if (v->content != NULL)
  {
    /* free data array if it's owned by the vector */
    if (NV_OWN_DATA_OMP(v) && NV_DATA_OMP(v) != NULL)
    {
      free(NV_DATA_OMP(v));
      NV_DATA_OMP(v) = NULL;
    }
    free(v->content);
    v->content = NULL;
  }

  /* free ops and vector */
  if (v->ops != NULL)
  {
    free(v->ops);
    v->ops = NULL;
  }
  free(v);
  v = NULL;

  return;
}

/* ----------------------------------------------------------------------------
 * Get storage requirement for N_Vector
 */

/* ----------------------------------------------------------------------------
 * Get vector data pointer
 */

sunrealtype* nvGetArrayPointer_OpenMP(N_Vector v)
{
  return ((sunrealtype*)NV_DATA_OMP(v));
}

/* ----------------------------------------------------------------------------
 * Set vector data pointer
 */

void nvSetArrayPointer_OpenMP(sunrealtype* v_data, N_Vector v)
{
  if (NV_LENGTH_OMP(v) > 0) { NV_DATA_OMP(v) = v_data; }

  return;
}

/* ----------------------------------------------------------------------------
 * Compute linear combination z[i] = a*x[i]+b*y[i]
 */

void nvLinearSum_OpenMP(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                        N_Vector z)
{
  sunindextype i, N;
  sunrealtype c, *xd, *yd, *zd;
  N_Vector v1, v2;
  sunbooleantype test;

  i  = 0; /* initialize to suppress clang warning */
  xd = yd = zd = NULL;

  if ((b == ONE) && (z == y))
  { /* BLAS usage: axpy y <- ax+y */
    Vaxpy_OpenMP(a, x, y);
    return;
  }

  if ((a == ONE) && (z == x))
  { /* BLAS usage: axpy x <- by+x */
    Vaxpy_OpenMP(b, y, x);
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE))
  {
    VSum_OpenMP(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE)))
  {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff_OpenMP(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE))
  {
    c  = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1_OpenMP(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE))
  {
    c  = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2_OpenMP(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b)
  {
    VScaleSum_OpenMP(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b)
  {
    VScaleDiff_OpenMP(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N, a, b, xd, yd, zd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { zd[i] = (a * xd[i]) + (b * yd[i]); }

  return;
}

/* ----------------------------------------------------------------------------
 * Assigns constant value to all vector elements, z[i] = c
 */

void nvConst_OpenMP(sunrealtype c, N_Vector z)
{
  sunindextype i, N;
  sunrealtype* zd;

  i  = 0; /* initialize to suppress clang warning */
  zd = NULL;

  N  = NV_LENGTH_OMP(z);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N, c, zd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(z))
  for (i = 0; i < N; i++) { zd[i] = c; }

  return;
}

/* ----------------------------------------------------------------------------
 * Compute componentwise product z[i] = x[i]*y[i]
 */

void nvProd_OpenMP(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  i  = 0; /* initialize to suppress clang warning */
  xd = yd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N, xd, yd, zd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { zd[i] = xd[i] * yd[i]; }

  return;
}

/* ----------------------------------------------------------------------------
 * Compute componentwise division z[i] = x[i]/y[i]
 */

void nvDiv_OpenMP(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  i  = 0; /* initialize to suppress clang warning */
  xd = yd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N, xd, yd, zd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { zd[i] = xd[i] / yd[i]; }

  return;
}

/* ----------------------------------------------------------------------------
 * Compute scaler multiplication z[i] = c*x[i]
 */

void nvScale_OpenMP(sunrealtype c, N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  i  = 0; /* initialize to suppress clang warning */
  xd = zd = NULL;

  if (z == x)
  { /* BLAS usage: scale x <- cx */
    VScaleBy_OpenMP(c, x);
    return;
  }

  if (c == ONE) { VCopy_OpenMP(x, z); }
  else if (c == -ONE) { VNeg_OpenMP(x, z); }
  else
  {
    N  = NV_LENGTH_OMP(x);
    xd = NV_DATA_OMP(x);
    zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N, c, xd, zd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
    for (i = 0; i < N; i++) { zd[i] = c * xd[i]; }
  }

  return;
}

/* ----------------------------------------------------------------------------
 * Compute absolute value of vector components z[i] = SUNRabs(x[i])
 */

void nvAbs_OpenMP(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  i  = 0; /* initialize to suppress clang warning */
  xd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { zd[i] = SUNRabs(xd[i]); }

  return;
}

/* ----------------------------------------------------------------------------
 * Compute componentwise inverse z[i] = 1 / x[i]
 */

void nvInv_OpenMP(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  i  = 0; /* initialize to suppress clang warning */
  xd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N, xd, zd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { zd[i] = ONE / xd[i]; }

  return;
}

/* ----------------------------------------------------------------------------
 * Compute componentwise addition of a scaler to a vector z[i] = x[i] + b
 */

void nvAddConst_OpenMP(N_Vector x, sunrealtype b, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  i  = 0; /* initialize to suppress clang warning */
  xd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N, b, xd, zd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { zd[i] = xd[i] + b; }

  return;
}

/* ----------------------------------------------------------------------------
 * Computes the dot product of two vectors, a = sum(x[i]*y[i])
 */

sunrealtype nvDotProd_OpenMP(N_Vector x, N_Vector y)
{
  sunindextype i, N;
  sunrealtype sum, *xd, *yd;

  i   = 0; /* initialize to suppress clang warning */
  sum = ZERO;
  xd = yd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);

#pragma omp parallel for default(none) private(i) shared(N, xd, yd) \
  reduction(+ : sum) schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { sum += xd[i] * yd[i]; }

  return (sum);
}

/* ----------------------------------------------------------------------------
 * Computes max norm of a vector
 */

sunrealtype nvMaxNorm_OpenMP(N_Vector x)
{
  sunindextype i, N;
  sunrealtype tmax, max, *xd;

  i   = 0; /* initialize to suppress clang warning */
  max = ZERO;
  xd  = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);

#pragma omp parallel default(none) private(i, tmax) shared(N, max, xd) \
  num_threads(NV_NUM_THREADS_OMP(x))
  {
    tmax = ZERO;
#pragma omp for schedule(static)
    for (i = 0; i < N; i++)
    {
      if (SUNRabs(xd[i]) > tmax) { tmax = SUNRabs(xd[i]); }
    }
#pragma omp critical
    {
      if (tmax > max) { max = tmax; }
    }
  }
  return (max);
}

/* ----------------------------------------------------------------------------
 * Computes weighted root mean square norm of a vector
 */

sunrealtype nvWrmsNorm_OpenMP(N_Vector x, N_Vector w)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype sqr_sum = nvWSqrSumLocal_OpenMP(x, w);
  SUNCheckLastErrNoRet();
  return (SUNRsqrt(sqr_sum / (NV_LENGTH_OMP(x))));
}

/* ----------------------------------------------------------------------------
 * Computes weighted root mean square norm of a masked vector
 */

sunrealtype nvWrmsNormMask_OpenMP(N_Vector x, N_Vector w, N_Vector id)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype sqr_sum = nvWSqrSumMaskLocal_OpenMP(x, w, id);
  SUNCheckLastErrNoRet();
  return (SUNRsqrt(sqr_sum / (NV_LENGTH_OMP(x))));
}

/* ----------------------------------------------------------------------------
 * Finds the minimum component of a vector
 */

sunrealtype nvMin_OpenMP(N_Vector x)
{
  sunindextype i, N;
  sunrealtype min, *xd;
  sunrealtype tmin;

  i  = 0; /* initialize to suppress clang warning */
  xd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);

  min = xd[0];

#pragma omp parallel default(none) private(i, tmin) shared(N, min, xd) \
  num_threads(NV_NUM_THREADS_OMP(x))
  {
    tmin = xd[0];
#pragma omp for schedule(static)
    for (i = 1; i < N; i++)
    {
      if (xd[i] < tmin) { tmin = xd[i]; }
    }
    if (tmin < min)
    {
#pragma omp critical
      {
        if (tmin < min) { min = tmin; }
      }
    }
  }

  return (min);
}

/* ----------------------------------------------------------------------------
 * Computes weighted L2 norm of a vector
 */

sunrealtype nvWL2Norm_OpenMP(N_Vector x, N_Vector w)
{
  sunindextype i, N;
  sunrealtype sum, *xd, *wd;

  i   = 0; /* initialize to suppress clang warning */
  sum = ZERO;
  xd = wd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  wd = NV_DATA_OMP(w);

#pragma omp parallel for default(none) private(i) shared(N, xd, wd) \
  reduction(+ : sum) schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { sum += SUNSQR(xd[i] * wd[i]); }

  return (SUNRsqrt(sum));
}

/* ----------------------------------------------------------------------------
 * Computes L1 norm of a vector
 */

sunrealtype nvL1Norm_OpenMP(N_Vector x)
{
  sunindextype i, N;
  sunrealtype sum, *xd;

  i   = 0; /* initialize to suppress clang warning */
  sum = ZERO;
  xd  = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);

#pragma omp parallel for default(none) private(i) shared(N, xd) \
  reduction(+ : sum) schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { sum += SUNRabs(xd[i]); }

  return (sum);
}

/* ----------------------------------------------------------------------------
 * Compare vector component values to a scaler
 */

void nvCompare_OpenMP(sunrealtype c, N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  i  = 0; /* initialize to suppress clang warning */
  xd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N, c, xd, zd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { zd[i] = (SUNRabs(xd[i]) >= c) ? ONE : ZERO; }

  return;
}

/* ----------------------------------------------------------------------------
 * Compute componentwise inverse z[i] = ONE/x[i] and checks if x[i] == ZERO
 */

sunbooleantype nvInvTest_OpenMP(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd, val;

  i  = 0; /* initialize to suppress clang warning */
  xd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  zd = NV_DATA_OMP(z);

  val = ZERO;

#pragma omp parallel for default(none) private(i) shared(N, val, xd, zd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++)
  {
    if (xd[i] == ZERO) { val = ONE; }
    else { zd[i] = ONE / xd[i]; }
  }

  if (val > ZERO) { return (SUNFALSE); }
  else { return (SUNTRUE); }
}

/* ----------------------------------------------------------------------------
 * Compute constraint mask of a vector
 */

sunbooleantype nvConstrMask_OpenMP(N_Vector c, N_Vector x, N_Vector m)
{
  sunindextype i, N;
  sunrealtype temp;
  sunrealtype *cd, *xd, *md;
  sunbooleantype test;

  i  = 0; /* initialize to suppress clang warning */
  cd = xd = md = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  cd = NV_DATA_OMP(c);
  md = NV_DATA_OMP(m);

  temp = ZERO;

#pragma omp parallel for default(none) private(i, test) \
  shared(N, xd, cd, md, temp) schedule(static)          \
  num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++)
  {
    md[i] = ZERO;

    /* Continue if no constraints were set for the variable */
    if (cd[i] == ZERO) { continue; }

    /* Check if a set constraint has been violated */
    test = (SUNRabs(cd[i]) > ONEPT5 && xd[i] * cd[i] <= ZERO) ||
           (SUNRabs(cd[i]) > HALF && xd[i] * cd[i] < ZERO);
    if (test) { temp = md[i] = ONE; /* Here is a race to write to temp */ }
  }
  /* Return false if any constraint was violated */
  return (temp == ONE) ? SUNFALSE : SUNTRUE;
}

/* ----------------------------------------------------------------------------
 * Compute minimum componentwise quotient
 */

sunrealtype nvMinQuotient_OpenMP(N_Vector num, N_Vector denom)
{
  sunindextype i, N;
  sunrealtype *nd, *dd, min, tmin, val;

  i  = 0; /* initialize to suppress clang warning */
  nd = dd = NULL;

  N  = NV_LENGTH_OMP(num);
  nd = NV_DATA_OMP(num);
  dd = NV_DATA_OMP(denom);

  min = SUN_BIG_REAL;

#pragma omp parallel default(none) private(i, tmin, val) \
  shared(N, min, nd, dd) num_threads(NV_NUM_THREADS_OMP(num))
  {
    tmin = SUN_BIG_REAL;
#pragma omp for schedule(static)
    for (i = 0; i < N; i++)
    {
      if (dd[i] != ZERO)
      {
        val = nd[i] / dd[i];
        if (val < tmin) { tmin = val; }
      }
    }
    if (tmin < min)
    {
#pragma omp critical
      {
        if (tmin < min) { min = tmin; }
      }
    }
  }

  return (min);
}

/* ----------------------------------------------------------------------------
 * Computes weighted square sum of a vector
 */

sunrealtype nvWSqrSumLocal_OpenMP(N_Vector x, N_Vector w)
{
  sunindextype i, N;
  sunrealtype sum, *xd, *wd;

  i   = 0; /* initialize to suppress clang warning */
  sum = ZERO;
  xd = wd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  wd = NV_DATA_OMP(w);

#pragma omp parallel for default(none) private(i) shared(N, xd, wd) \
  reduction(+ : sum) schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { sum += SUNSQR(xd[i] * wd[i]); }

  return (sum);
}

/* ----------------------------------------------------------------------------
 * Computes weighted square sum of a masked vector
 */

sunrealtype nvWSqrSumMaskLocal_OpenMP(N_Vector x, N_Vector w, N_Vector id)
{
  sunindextype i, N;
  sunrealtype sum, *xd, *wd, *idd;

  i   = 0; /* initialize to suppress clang warning */
  sum = ZERO;
  xd = wd = idd = NULL;

  N   = NV_LENGTH_OMP(x);
  xd  = NV_DATA_OMP(x);
  wd  = NV_DATA_OMP(w);
  idd = NV_DATA_OMP(id);

#pragma omp parallel for default(none) private(i) shared(N, xd, wd, idd) \
  reduction(+ : sum) schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++)
  {
    if (idd[i] > ZERO) { sum += SUNSQR(xd[i] * wd[i]); }
  }

  return (sum);
}

/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */

SUNErrCode nvLinearCombination_OpenMP(int nvec, sunrealtype* c, N_Vector* X,
                                      N_Vector z)
{
  SUNFunctionBegin(X[0]->sunctx);

  int i;
  sunindextype j, N;
  sunrealtype* zd = NULL;
  sunrealtype* xd = NULL;

  i = 0; /* initialize to suppress clang warning */
  j = 0;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VScale */
  if (nvec == 1)
  {
    nvScale_OpenMP(c[0], X[0], z);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* should have called N_VLinearSum */
  if (nvec == 2)
  {
    nvLinearSum_OpenMP(c[0], X[0], c[1], X[1], z);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N  = NV_LENGTH_OMP(z);
  zd = NV_DATA_OMP(z);

  /*
   * X[0] += c[i]*X[i], i = 1,...,nvec-1
   */
  if ((X[0] == z) && (c[0] == ONE))
  {
#pragma omp parallel default(none) private(i, j, xd) shared(nvec, X, N, c, zd) \
  num_threads(NV_NUM_THREADS_OMP(z))
    {
      for (i = 1; i < nvec; i++)
      {
        xd = NV_DATA_OMP(X[i]);
#pragma omp for schedule(static)
        for (j = 0; j < N; j++) { zd[j] += c[i] * xd[j]; }
      }
    }
    return SUN_SUCCESS;
  }

  /*
   * X[0] = c[0] * X[0] + sum{ c[i] * X[i] }, i = 1,...,nvec-1
   */
  if (X[0] == z)
  {
#pragma omp parallel default(none) private(i, j, xd) shared(nvec, X, N, c, zd) \
  num_threads(NV_NUM_THREADS_OMP(z))
    {
#pragma omp for schedule(static)
      for (j = 0; j < N; j++) { zd[j] *= c[0]; }

      for (i = 1; i < nvec; i++)
      {
        xd = NV_DATA_OMP(X[i]);
#pragma omp for schedule(static)
        for (j = 0; j < N; j++) { zd[j] += c[i] * xd[j]; }
      }
    }
    return SUN_SUCCESS;
  }

  /*
   * z = sum{ c[i] * X[i] }, i = 0,...,nvec-1
   */
#pragma omp parallel default(none) private(i, j, xd) shared(nvec, X, N, c, zd) \
  num_threads(NV_NUM_THREADS_OMP(z))
  {
    xd = NV_DATA_OMP(X[0]);
#pragma omp for schedule(static)
    for (j = 0; j < N; j++) { zd[j] = c[0] * xd[j]; }

    for (i = 1; i < nvec; i++)
    {
      xd = NV_DATA_OMP(X[i]);
#pragma omp for schedule(static)
      for (j = 0; j < N; j++) { zd[j] += c[i] * xd[j]; }
    }
  }
  return SUN_SUCCESS;
}

SUNErrCode nvScaleAddMulti_OpenMP(int nvec, sunrealtype* a, N_Vector x,
                                  N_Vector* Y, N_Vector* Z)
{
  SUNFunctionBegin(x->sunctx);

  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  i = 0; /* initialize to suppress clang warning */
  j = 0;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VLinearSum */
  if (nvec == 1)
  {
    nvLinearSum_OpenMP(a[0], x, ONE, Y[0], Z[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z)
  {
#pragma omp parallel default(none) private(i, j, yd) shared(nvec, Y, N, a, xd) \
  num_threads(NV_NUM_THREADS_OMP(x))
    {
      for (i = 0; i < nvec; i++)
      {
        yd = NV_DATA_OMP(Y[i]);
#pragma omp for schedule(static)
        for (j = 0; j < N; j++) { yd[j] += a[i] * xd[j]; }
      }
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
#pragma omp parallel default(none) private(i, j, yd, zd) \
  shared(nvec, Y, Z, N, a, xd) num_threads(NV_NUM_THREADS_OMP(x))
  {
    for (i = 0; i < nvec; i++)
    {
      yd = NV_DATA_OMP(Y[i]);
      zd = NV_DATA_OMP(Z[i]);
#pragma omp for schedule(static)
      for (j = 0; j < N; j++) { zd[j] = a[i] * xd[j] + yd[j]; }
    }
  }
  return SUN_SUCCESS;
}

SUNErrCode nvDotProdMulti_OpenMP(int nvec, N_Vector x, N_Vector* Y,
                                 sunrealtype* dotprods)
{
  SUNFunctionBegin(x->sunctx);

  int i;
  sunindextype j, N;
  sunrealtype sum;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;

  i = 0; /* initialize to suppress clang warning */
  j = 0;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VDotProd */
  if (nvec == 1)
  {
    dotprods[0] = nvDotProd_OpenMP(x, Y[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);

  /* initialize dot products */
  for (i = 0; i < nvec; i++) { dotprods[i] = ZERO; }

  /* compute multiple dot products */
#pragma omp parallel default(none) private(i, j, yd, sum) \
  shared(nvec, Y, N, xd, dotprods) num_threads(NV_NUM_THREADS_OMP(x))
  {
    for (i = 0; i < nvec; i++)
    {
      yd  = NV_DATA_OMP(Y[i]);
      sum = ZERO;
#pragma omp for schedule(static)
      for (j = 0; j < N; j++) { sum += xd[j] * yd[j]; }
#pragma omp critical
      {
        dotprods[i] += sum;
      }
    }
  }

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------
 */

SUNErrCode nvLinearSumVectorArray_OpenMP(int nvec, sunrealtype a, N_Vector* X,
                                         sunrealtype b, N_Vector* Y, N_Vector* Z)
{
  SUNFunctionBegin(X[0]->sunctx);

  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;
  sunrealtype c;
  N_Vector* V1;
  N_Vector* V2;
  sunbooleantype test;

  i = 0; /* initialize to suppress clang warning */
  j = 0;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VLinearSum */
  if (nvec == 1)
  {
    nvLinearSum_OpenMP(a, X[0], b, Y[0], Z[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* BLAS usage: axpy y <- ax+y */
  if ((b == ONE) && (Z == Y))
  {
    VaxpyVectorArray_OpenMP(nvec, a, X, Y);
    return SUN_SUCCESS;
  }

  /* BLAS usage: axpy x <- by+x */
  if ((a == ONE) && (Z == X))
  {
    VaxpyVectorArray_OpenMP(nvec, b, Y, X);
    return SUN_SUCCESS;
  }

  /* Case: a == b == 1.0 */
  if ((a == ONE) && (b == ONE))
  {
    VSumVectorArray_OpenMP(nvec, X, Y, Z);
    return SUN_SUCCESS;
  }

  /* Cases:                    */
  /*   (1) a == 1.0, b = -1.0, */
  /*   (2) a == -1.0, b == 1.0 */
  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE)))
  {
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    VDiffVectorArray_OpenMP(nvec, V2, V1, Z);
    return SUN_SUCCESS;
  }

  /* Cases:                                                  */
  /*   (1) a == 1.0, b == other or 0.0,                      */
  /*   (2) a == other or 0.0, b == 1.0                       */
  /* if a or b is 0.0, then user should have called N_VScale */
  if ((test = (a == ONE)) || (b == ONE))
  {
    c  = test ? b : a;
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    VLin1VectorArray_OpenMP(nvec, c, V1, V2, Z);
    return SUN_SUCCESS;
  }

  /* Cases:                     */
  /*   (1) a == -1.0, b != 1.0, */
  /*   (2) a != 1.0, b == -1.0  */
  if ((test = (a == -ONE)) || (b == -ONE))
  {
    c  = test ? b : a;
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    VLin2VectorArray_OpenMP(nvec, c, V1, V2, Z);
    return SUN_SUCCESS;
  }

  /* Case: a == b                                                         */
  /* catches case both a and b are 0.0 - user should have called N_VConst */
  if (a == b)
  {
    VScaleSumVectorArray_OpenMP(nvec, a, X, Y, Z);
    return SUN_SUCCESS;
  }

  /* Case: a == -b */
  if (a == -b)
  {
    VScaleDiffVectorArray_OpenMP(nvec, a, X, Y, Z);
    return SUN_SUCCESS;
  }

  /* Do all cases not handled above:                               */
  /*   (1) a == other, b == 0.0 - user should have called N_VScale */
  /*   (2) a == 0.0, b == other - user should have called N_VScale */
  /*   (3) a,b == other, a !=b, a != -b                            */

  /* get vector length */
  N = NV_LENGTH_OMP(Z[0]);

  /* compute linear sum for each vector pair in vector arrays */
#pragma omp parallel default(none) private(i, j, xd, yd, zd) \
  shared(nvec, X, Y, Z, N, a, b) num_threads(NV_NUM_THREADS_OMP(Z[0]))
  {
    for (i = 0; i < nvec; i++)
    {
      xd = NV_DATA_OMP(X[i]);
      yd = NV_DATA_OMP(Y[i]);
      zd = NV_DATA_OMP(Z[i]);
#pragma omp for schedule(static)
      for (j = 0; j < N; j++) { zd[j] = a * xd[j] + b * yd[j]; }
    }
  }

  return SUN_SUCCESS;
}

SUNErrCode nvScaleVectorArray_OpenMP(int nvec, sunrealtype* c, N_Vector* X,
                                     N_Vector* Z)
{
  SUNFunctionBegin(X[0]->sunctx);

  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* zd = NULL;

  i = 0; /* initialize to suppress clang warning */
  j = 0;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VScale */
  if (nvec == 1)
  {
    nvScale_OpenMP(c[0], X[0], Z[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length */
  N = NV_LENGTH_OMP(Z[0]);

  /*
   * X[i] *= c[i]
   */
  if (X == Z)
  {
#pragma omp parallel default(none) private(i, j, xd) shared(nvec, X, N, c) \
  num_threads(NV_NUM_THREADS_OMP(Z[0]))
    {
      for (i = 0; i < nvec; i++)
      {
        xd = NV_DATA_OMP(X[i]);
#pragma omp for schedule(static)
        for (j = 0; j < N; j++) { xd[j] *= c[i]; }
      }
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[i] = c[i] * X[i]
   */
#pragma omp parallel default(none) private(i, j, xd, zd) \
  shared(nvec, X, Z, N, c) num_threads(NV_NUM_THREADS_OMP(Z[0]))
  {
    for (i = 0; i < nvec; i++)
    {
      xd = NV_DATA_OMP(X[i]);
      zd = NV_DATA_OMP(Z[i]);
#pragma omp for schedule(static)
      for (j = 0; j < N; j++) { zd[j] = c[i] * xd[j]; }
    }
  }
  return SUN_SUCCESS;
}

SUNErrCode nvConstVectorArray_OpenMP(int nvec, sunrealtype c, N_Vector* Z)
{
  SUNFunctionBegin(Z[0]->sunctx);

  int i;
  sunindextype j, N;
  sunrealtype* zd = NULL;

  i = 0; /* initialize to suppress clang warning */
  j = 0;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VConst */
  if (nvec == 1)
  {
    nvConst_OpenMP(c, Z[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length */
  N = NV_LENGTH_OMP(Z[0]);

  /* set each vector in the vector array to a constant */
#pragma omp parallel default(none) private(i, j, zd) shared(nvec, Z, N, c) \
  num_threads(NV_NUM_THREADS_OMP(Z[0]))
  {
    for (i = 0; i < nvec; i++)
    {
      zd = NV_DATA_OMP(Z[i]);
#pragma omp for schedule(static)
      for (j = 0; j < N; j++) { zd[j] = c; }
    }
  }

  return SUN_SUCCESS;
}

SUNErrCode nvWrmsNormVectorArray_OpenMP(int nvec, N_Vector* X, N_Vector* W,
                                        sunrealtype* nrm)
{
  SUNFunctionBegin(X[0]->sunctx);

  int i;
  sunindextype j, N;
  sunrealtype sum;
  sunrealtype* wd = NULL;
  sunrealtype* xd = NULL;

  i = 0; /* initialize to suppress clang warning */
  j = 0;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VWrmsNorm */
  if (nvec == 1)
  {
    nrm[0] = nvWrmsNorm_OpenMP(X[0], W[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length */
  N = NV_LENGTH_OMP(X[0]);

  /* initialize norms */
  for (i = 0; i < nvec; i++) { nrm[i] = ZERO; }

  /* compute the WRMS norm for each vector in the vector array */
#pragma omp parallel default(none) private(i, j, xd, wd, sum) \
  shared(nvec, X, W, N, nrm) num_threads(NV_NUM_THREADS_OMP(X[0]))
  {
    for (i = 0; i < nvec; i++)
    {
      xd  = NV_DATA_OMP(X[i]);
      wd  = NV_DATA_OMP(W[i]);
      sum = ZERO;
#pragma omp for schedule(static)
      for (j = 0; j < N; j++) { sum += SUNSQR(xd[j] * wd[j]); }
#pragma omp critical
      {
        nrm[i] += sum;
      }
    }
  }

  for (i = 0; i < nvec; i++) { nrm[i] = SUNRsqrt(nrm[i] / N); }

  return SUN_SUCCESS;
}

SUNErrCode nvWrmsNormMaskVectorArray_OpenMP(int nvec, N_Vector* X, N_Vector* W,
                                            N_Vector id, sunrealtype* nrm)
{
  SUNFunctionBegin(X[0]->sunctx);

  int i;
  sunindextype j, N;
  sunrealtype sum;
  sunrealtype* wd  = NULL;
  sunrealtype* xd  = NULL;
  sunrealtype* idd = NULL;

  i = 0; /* initialize to suppress clang warning */
  j = 0;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VWrmsNorm */
  if (nvec == 1)
  {
    nrm[0] = nvWrmsNormMask_OpenMP(X[0], W[0], id);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length and mask data array */
  N   = NV_LENGTH_OMP(X[0]);
  idd = NV_DATA_OMP(id);

  /* initialize norms */
  for (i = 0; i < nvec; i++) { nrm[i] = ZERO; }

  /* compute the WRMS norm for each vector in the vector array */
#pragma omp parallel default(none) private(i, j, xd, wd, sum) \
  shared(nvec, X, W, N, idd, nrm) num_threads(NV_NUM_THREADS_OMP(X[0]))
  {
    for (i = 0; i < nvec; i++)
    {
      xd  = NV_DATA_OMP(X[i]);
      wd  = NV_DATA_OMP(W[i]);
      sum = ZERO;
#pragma omp for schedule(static)
      for (j = 0; j < N; j++)
      {
        if (idd[j] > ZERO) { sum += SUNSQR(xd[j] * wd[j]); }
      }
#pragma omp critical
      {
        nrm[i] += sum;
      }
    }
  }

  for (i = 0; i < nvec; i++) { nrm[i] = SUNRsqrt(nrm[i] / N); }

  return SUN_SUCCESS;
}

SUNErrCode nvScaleAddMultiVectorArray_OpenMP(int nvec, int nsum, sunrealtype* a,
                                             N_Vector* X, N_Vector** Y,
                                             N_Vector** Z)
{
  SUNFunctionBegin(X[0]->sunctx);

  int i, j;
  sunindextype k, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  N_Vector* YY;
  N_Vector* ZZ;

  i = 0; /* initialize to suppress clang warning */
  k = 0;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1 && nsum >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1)
  {
    /* should have called N_VLinearSum */
    if (nsum == 1)
    {
      nvLinearSum_OpenMP(a[0], X[0], ONE, Y[0][0], Z[0][0]);
      SUNCheckLastErr();
      return SUN_SUCCESS;
    }

    /* should have called N_VScaleAddMulti */
    YY = (N_Vector*)malloc(nsum * sizeof(N_Vector));
    SUNAssert(YY, SUN_ERR_MALLOC_FAIL);
    ZZ = (N_Vector*)malloc(nsum * sizeof(N_Vector));
    SUNAssert(ZZ, SUN_ERR_MALLOC_FAIL);

    for (j = 0; j < nsum; j++)
    {
      YY[j] = Y[j][0];
      ZZ[j] = Z[j][0];
    }

    SUNCheckCall(nvScaleAddMulti_OpenMP(nsum, a, X[0], YY, ZZ));

    free(YY);
    free(ZZ);
    return SUN_SUCCESS;
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 1)
  {
    SUNCheckCall(nvLinearSumVectorArray_OpenMP(nvec, a[0], X, ONE, Y[0], Z[0]));
    return SUN_SUCCESS;
  }

  /* ----------------------------
   * Compute multiple linear sums
   * ---------------------------- */

  /* get vector length */
  N = NV_LENGTH_OMP(X[0]);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z)
  {
#pragma omp parallel default(none) private(i, j, k, xd, yd) \
  shared(nvec, nsum, X, Y, N, a) num_threads(NV_NUM_THREADS_OMP(X[0]))
    {
      for (i = 0; i < nvec; i++)
      {
        xd = NV_DATA_OMP(X[i]);
        for (j = 0; j < nsum; j++)
        {
          yd = NV_DATA_OMP(Y[j][i]);
#pragma omp for schedule(static)
          for (k = 0; k < N; k++) { yd[k] += a[j] * xd[k]; }
        }
      }
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
#pragma omp parallel default(none) private(i, j, k, xd, yd, zd) \
  shared(nvec, nsum, X, Y, Z, N, a) num_threads(NV_NUM_THREADS_OMP(X[0]))
  {
    for (i = 0; i < nvec; i++)
    {
      xd = NV_DATA_OMP(X[i]);
      for (j = 0; j < nsum; j++)
      {
        yd = NV_DATA_OMP(Y[j][i]);
        zd = NV_DATA_OMP(Z[j][i]);
#pragma omp for schedule(static)
        for (k = 0; k < N; k++) { zd[k] = a[j] * xd[k] + yd[k]; }
      }
    }
  }
  return SUN_SUCCESS;
}

SUNErrCode nvLinearCombinationVectorArray_OpenMP(int nvec, int nsum,
                                                 sunrealtype* c, N_Vector** X,
                                                 N_Vector* Z)
{
  SUNFunctionBegin(X[0][0]->sunctx);

  int i;          /* vector arrays index in summation [0,nsum) */
  int j;          /* vector index in vector array     [0,nvec) */
  sunindextype k; /* element index in vector          [0,N)    */
  sunindextype N;
  sunrealtype* zd = NULL;
  sunrealtype* xd = NULL;

  sunrealtype* ctmp;
  N_Vector* Y;

  i = 0; /* initialize to suppress clang warning */
  k = 0;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1 && nsum >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1)
  {
    /* should have called N_VScale */
    if (nsum == 1)
    {
      nvScale_OpenMP(c[0], X[0][0], Z[0]);
      SUNCheckLastErr();
      return SUN_SUCCESS;
    }

    /* should have called N_VLinearSum */
    if (nsum == 2)
    {
      nvLinearSum_OpenMP(c[0], X[0][0], c[1], X[1][0], Z[0]);
      SUNCheckLastErr();
      return SUN_SUCCESS;
    }

    /* should have called N_VLinearCombination */
    Y = (N_Vector*)malloc(nsum * sizeof(N_Vector));
    SUNAssert(Y, SUN_ERR_MALLOC_FAIL);

    for (i = 0; i < nsum; i++) { Y[i] = X[i][0]; }

    SUNCheckCall(nvLinearCombination_OpenMP(nsum, c, Y, Z[0]));

    free(Y);
    return SUN_SUCCESS;
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VScaleVectorArray */
  if (nsum == 1)
  {
    ctmp = (sunrealtype*)malloc(nvec * sizeof(sunrealtype));

    for (j = 0; j < nvec; j++) { ctmp[j] = c[0]; }

    SUNCheckCall(nvScaleVectorArray_OpenMP(nvec, ctmp, X[0], Z));

    free(ctmp);
    return SUN_SUCCESS;
  }

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 2)
  {
    SUNCheckCall(nvLinearSumVectorArray_OpenMP(nvec, c[0], X[0], c[1], X[1], Z));
    return SUN_SUCCESS;
  }

  /* --------------------------
   * Compute linear combination
   * -------------------------- */

  /* get vector length */
  N = NV_LENGTH_OMP(Z[0]);

  /*
   * X[0][j] += c[i]*X[i][j], i = 1,...,nvec-1
   */
  if ((X[0] == Z) && (c[0] == ONE))
  {
#pragma omp parallel default(none) private(i, j, k, xd, zd) \
  shared(nvec, nsum, X, Z, N, c) num_threads(NV_NUM_THREADS_OMP(Z[0]))
    {
      for (j = 0; j < nvec; j++)
      {
        zd = NV_DATA_OMP(Z[j]);
        for (i = 1; i < nsum; i++)
        {
          xd = NV_DATA_OMP(X[i][j]);
#pragma omp for schedule(static)
          for (k = 0; k < N; k++) { zd[k] += c[i] * xd[k]; }
        }
      }
    }
    return SUN_SUCCESS;
  }

  /*
   * X[0][j] = c[0] * X[0][j] + sum{ c[i] * X[i][j] }, i = 1,...,nvec-1
   */
  if (X[0] == Z)
  {
#pragma omp parallel default(none) private(i, j, k, xd, zd) \
  shared(nvec, nsum, X, Z, N, c) num_threads(NV_NUM_THREADS_OMP(Z[0]))
    {
      for (j = 0; j < nvec; j++)
      {
        zd = NV_DATA_OMP(Z[j]);
#pragma omp for schedule(static)
        for (k = 0; k < N; k++) { zd[k] *= c[0]; }
        for (i = 1; i < nsum; i++)
        {
          xd = NV_DATA_OMP(X[i][j]);
#pragma omp for schedule(static)
          for (k = 0; k < N; k++) { zd[k] += c[i] * xd[k]; }
        }
      }
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[j] = sum{ c[i] * X[i][j] }, i = 0,...,nvec-1
   */
#pragma omp parallel default(none) private(i, j, k, xd, zd) \
  shared(nvec, nsum, X, Z, N, c) num_threads(NV_NUM_THREADS_OMP(Z[0]))
  {
    for (j = 0; j < nvec; j++)
    {
      /* scale first vector in the sum into the output vector */
      xd = NV_DATA_OMP(X[0][j]);
      zd = NV_DATA_OMP(Z[j]);
#pragma omp for schedule(static)
      for (k = 0; k < N; k++) { zd[k] = c[0] * xd[k]; }
      /* scale and sum remaining vectors into the output vector */
      for (i = 1; i < nsum; i++)
      {
        xd = NV_DATA_OMP(X[i][j]);
#pragma omp for schedule(static)
        for (k = 0; k < N; k++) { zd[k] += c[i] * xd[k]; }
      }
    }
  }
  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * OPTIONAL XBraid interface operations
 * -----------------------------------------------------------------
 */

SUNErrCode nvBufSize_OpenMP(N_Vector x, sunindextype* size)
{
  *size = NV_LENGTH_OMP(x) * ((sunindextype)sizeof(sunrealtype));
  return SUN_SUCCESS;
}

SUNErrCode nvBufPack_OpenMP(N_Vector x, void* buf)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype i, N;
  sunrealtype* xd = NULL;
  sunrealtype* bd = NULL;

  SUNAssert(buf, SUN_ERR_ARG_CORRUPT);

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  bd = (sunrealtype*)buf;

#pragma omp for schedule(static)
  for (i = 0; i < N; i++) { bd[i] = xd[i]; }

  return SUN_SUCCESS;
}

SUNErrCode nvBufUnpack_OpenMP(N_Vector x, void* buf)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype i, N;
  sunrealtype* xd = NULL;
  sunrealtype* bd = NULL;

  SUNAssert(buf, SUN_ERR_ARG_CORRUPT);

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  bd = (sunrealtype*)buf;

#pragma omp for schedule(static)
  for (i = 0; i < N; i++) { xd[i] = bd[i]; }

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * private functions for special cases of vector operations
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Copy vector components into a second vector
 */

static void VCopy_OpenMP(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  i  = 0; /* initialize to suppress clang warning */
  xd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N, xd, zd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { zd[i] = xd[i]; }

  return;
}

/* ----------------------------------------------------------------------------
 * Compute vector sum
 */

static void VSum_OpenMP(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  i  = 0; /* initialize to suppress clang warning */
  xd = yd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N, xd, yd, zd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { zd[i] = xd[i] + yd[i]; }

  return;
}

/* ----------------------------------------------------------------------------
 * Compute vector difference
 */

static void VDiff_OpenMP(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  i  = 0; /* initialize to suppress clang warning */
  xd = yd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N, xd, yd, zd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { zd[i] = xd[i] - yd[i]; }

  return;
}

/* ----------------------------------------------------------------------------
 * Compute the negative of a vector
 */

static void VNeg_OpenMP(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  i  = 0; /* initialize to suppress clang warning */
  xd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N, xd, zd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { zd[i] = -xd[i]; }

  return;
}

/* ----------------------------------------------------------------------------
 * Compute scaled vector sum
 */

static void VScaleSum_OpenMP(sunrealtype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  i  = 0; /* initialize to suppress clang warning */
  xd = yd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N, c, xd, yd, zd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { zd[i] = c * (xd[i] + yd[i]); }

  return;
}

/* ----------------------------------------------------------------------------
 * Compute scaled vector difference
 */

static void VScaleDiff_OpenMP(sunrealtype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  i  = 0; /* initialize to suppress clang warning */
  xd = yd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N, c, xd, yd, zd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { zd[i] = c * (xd[i] - yd[i]); }

  return;
}

/* ----------------------------------------------------------------------------
 * Compute vector sum z[i] = a*x[i]+y[i]
 */

static void VLin1_OpenMP(sunrealtype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  i  = 0; /* initialize to suppress clang warning */
  xd = yd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N, a, xd, yd, zd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { zd[i] = (a * xd[i]) + yd[i]; }

  return;
}

/* ----------------------------------------------------------------------------
 * Compute vector difference z[i] = a*x[i]-y[i]
 */

static void VLin2_OpenMP(sunrealtype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  i  = 0; /* initialize to suppress clang warning */
  xd = yd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N, a, xd, yd, zd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { zd[i] = (a * xd[i]) - yd[i]; }

  return;
}

/* ----------------------------------------------------------------------------
 * Compute special cases of linear sum
 */

static void Vaxpy_OpenMP(sunrealtype a, N_Vector x, N_Vector y)
{
  sunindextype i, N;
  sunrealtype *xd, *yd;

  i  = 0; /* initialize to suppress clang warning */
  xd = yd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);

  if (a == ONE)
  {
#pragma omp parallel for default(none) private(i) shared(N, xd, yd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
    for (i = 0; i < N; i++) { yd[i] += xd[i]; }
    return;
  }

  if (a == -ONE)
  {
#pragma omp parallel for default(none) private(i) shared(N, xd, yd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
    for (i = 0; i < N; i++) { yd[i] -= xd[i]; }
    return;
  }

#pragma omp parallel for default(none) private(i) shared(N, a, xd, yd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { yd[i] += a * xd[i]; }

  return;
}

/* ----------------------------------------------------------------------------
 * Compute scaled vector x[i] = a*x[i]
 */

static void VScaleBy_OpenMP(sunrealtype a, N_Vector x)
{
  sunindextype i, N;
  sunrealtype* xd;

  i  = 0; /* initialize to suppress clang warning */
  xd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);

#pragma omp parallel for default(none) private(i) shared(N, a, xd) \
  schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) { xd[i] *= a; }

  return;
}

/*
 * -----------------------------------------------------------------
 * private functions for special cases of vector array operations
 * -----------------------------------------------------------------
 */

static void VSumVectorArray_OpenMP(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  i = 0; /* initialize to suppress clang warning */
  j = 0;

  N = NV_LENGTH_OMP(X[0]);

#pragma omp parallel default(none) private(i, j, xd, yd, zd) \
  shared(nvec, X, Y, Z, N) num_threads(NV_NUM_THREADS_OMP(X[0]))
  {
    for (i = 0; i < nvec; i++)
    {
      xd = NV_DATA_OMP(X[i]);
      yd = NV_DATA_OMP(Y[i]);
      zd = NV_DATA_OMP(Z[i]);
#pragma omp for schedule(static)
      for (j = 0; j < N; j++) { zd[j] = xd[j] + yd[j]; }
    }
  }
}

static void VDiffVectorArray_OpenMP(int nvec, N_Vector* X, N_Vector* Y,
                                    N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  i = 0; /* initialize to suppress clang warning */
  j = 0;

  N = NV_LENGTH_OMP(X[0]);

#pragma omp parallel default(none) private(i, j, xd, yd, zd) \
  shared(nvec, X, Y, Z, N) num_threads(NV_NUM_THREADS_OMP(X[0]))
  {
    for (i = 0; i < nvec; i++)
    {
      xd = NV_DATA_OMP(X[i]);
      yd = NV_DATA_OMP(Y[i]);
      zd = NV_DATA_OMP(Z[i]);
#pragma omp for schedule(static)
      for (j = 0; j < N; j++) { zd[j] = xd[j] - yd[j]; }
    }
  }
}

static void VScaleSumVectorArray_OpenMP(int nvec, sunrealtype c, N_Vector* X,
                                        N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  i = 0; /* initialize to suppress clang warning */
  j = 0;

  N = NV_LENGTH_OMP(X[0]);

#pragma omp parallel default(none) private(i, j, xd, yd, zd) \
  shared(nvec, X, Y, Z, N, c) num_threads(NV_NUM_THREADS_OMP(X[0]))
  {
    for (i = 0; i < nvec; i++)
    {
      xd = NV_DATA_OMP(X[i]);
      yd = NV_DATA_OMP(Y[i]);
      zd = NV_DATA_OMP(Z[i]);
#pragma omp for schedule(static)
      for (j = 0; j < N; j++) { zd[j] = c * (xd[j] + yd[j]); }
    }
  }
}

static void VScaleDiffVectorArray_OpenMP(int nvec, sunrealtype c, N_Vector* X,
                                         N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  i = 0; /* initialize to suppress clang warning */
  j = 0;

  N = NV_LENGTH_OMP(X[0]);

#pragma omp parallel default(none) private(i, j, xd, yd, zd) \
  shared(nvec, X, Y, Z, N, c) num_threads(NV_NUM_THREADS_OMP(X[0]))
  {
    for (i = 0; i < nvec; i++)
    {
      xd = NV_DATA_OMP(X[i]);
      yd = NV_DATA_OMP(Y[i]);
      zd = NV_DATA_OMP(Z[i]);
#pragma omp for schedule(static)
      for (j = 0; j < N; j++) { zd[j] = c * (xd[j] - yd[j]); }
    }
  }
}

static void VLin1VectorArray_OpenMP(int nvec, sunrealtype a, N_Vector* X,
                                    N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  i = 0; /* initialize to suppress clang warning */
  j = 0;

  N = NV_LENGTH_OMP(X[0]);

#pragma omp parallel default(none) private(i, j, xd, yd, zd) \
  shared(nvec, X, Y, Z, N, a) num_threads(NV_NUM_THREADS_OMP(X[0]))
  {
    for (i = 0; i < nvec; i++)
    {
      xd = NV_DATA_OMP(X[i]);
      yd = NV_DATA_OMP(Y[i]);
      zd = NV_DATA_OMP(Z[i]);
#pragma omp for schedule(static)
      for (j = 0; j < N; j++) { zd[j] = (a * xd[j]) + yd[j]; }
    }
  }
}

static void VLin2VectorArray_OpenMP(int nvec, sunrealtype a, N_Vector* X,
                                    N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  i = 0; /* initialize to suppress clang warning */
  j = 0;

  N = NV_LENGTH_OMP(X[0]);

#pragma omp parallel default(none) private(i, j, xd, yd, zd) \
  shared(nvec, X, Y, Z, N, a) num_threads(NV_NUM_THREADS_OMP(X[0]))
  {
    for (i = 0; i < nvec; i++)
    {
      xd = NV_DATA_OMP(X[i]);
      yd = NV_DATA_OMP(Y[i]);
      zd = NV_DATA_OMP(Z[i]);
#pragma omp for schedule(static)
      for (j = 0; j < N; j++) { zd[j] = (a * xd[j]) - yd[j]; }
    }
  }
}

static void VaxpyVectorArray_OpenMP(int nvec, sunrealtype a, N_Vector* X,
                                    N_Vector* Y)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;

  i = 0; /* initialize to suppress clang warning */
  j = 0;

  N = NV_LENGTH_OMP(X[0]);

  if (a == ONE)
  {
#pragma omp parallel default(none) private(i, j, xd, yd) \
  shared(nvec, X, Y, N, a) num_threads(NV_NUM_THREADS_OMP(X[0]))
    {
      for (i = 0; i < nvec; i++)
      {
        xd = NV_DATA_OMP(X[i]);
        yd = NV_DATA_OMP(Y[i]);
#pragma omp for schedule(static)
        for (j = 0; j < N; j++) { yd[j] += xd[j]; }
      }
    }
    return;
  }

  if (a == -ONE)
  {
#pragma omp parallel default(none) private(i, j, xd, yd) \
  shared(nvec, X, Y, N, a) num_threads(NV_NUM_THREADS_OMP(X[0]))
    {
      for (i = 0; i < nvec; i++)
      {
        xd = NV_DATA_OMP(X[i]);
        yd = NV_DATA_OMP(Y[i]);
#pragma omp for schedule(static)
        for (j = 0; j < N; j++) { yd[j] -= xd[j]; }
      }
    }
    return;
  }

#pragma omp parallel default(none) private(i, j, xd, yd) \
  shared(nvec, X, Y, N, a) num_threads(NV_NUM_THREADS_OMP(X[0]))
  {
    for (i = 0; i < nvec; i++)
    {
      xd = NV_DATA_OMP(X[i]);
      yd = NV_DATA_OMP(Y[i]);
#pragma omp for schedule(static)
      for (j = 0; j < N; j++) { yd[j] += a * xd[j]; }
    }
  }
}

/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

SUNErrCode N_VEnableFusedOps_OpenMP(N_Vector v, sunbooleantype tf)
{
  if (tf)
  {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = nvLinearCombination_OpenMP;
    v->ops->nvscaleaddmulti     = nvScaleAddMulti_OpenMP;
    v->ops->nvdotprodmulti      = nvDotProdMulti_OpenMP;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray     = nvLinearSumVectorArray_OpenMP;
    v->ops->nvscalevectorarray         = nvScaleVectorArray_OpenMP;
    v->ops->nvconstvectorarray         = nvConstVectorArray_OpenMP;
    v->ops->nvwrmsnormvectorarray      = nvWrmsNormVectorArray_OpenMP;
    v->ops->nvwrmsnormmaskvectorarray  = nvWrmsNormMaskVectorArray_OpenMP;
    v->ops->nvscaleaddmultivectorarray = nvScaleAddMultiVectorArray_OpenMP;
    v->ops->nvlinearcombinationvectorarray = nvLinearCombinationVectorArray_OpenMP;
    /* enable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = nvDotProdMulti_OpenMP;
  }
  else
  {
    /* disable all fused vector operations */
    v->ops->nvlinearcombination = NULL;
    v->ops->nvscaleaddmulti     = NULL;
    v->ops->nvdotprodmulti      = NULL;
    /* disable all vector array operations */
    v->ops->nvlinearsumvectorarray         = NULL;
    v->ops->nvscalevectorarray             = NULL;
    v->ops->nvconstvectorarray             = NULL;
    v->ops->nvwrmsnormvectorarray          = NULL;
    v->ops->nvwrmsnormmaskvectorarray      = NULL;
    v->ops->nvscaleaddmultivectorarray     = NULL;
    v->ops->nvlinearcombinationvectorarray = NULL;
    /* disable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = NULL;
  }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearCombination_OpenMP(N_Vector v, sunbooleantype tf)
{
  v->ops->nvlinearcombination = tf ? nvLinearCombination_OpenMP : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleAddMulti_OpenMP(N_Vector v, sunbooleantype tf)
{
  v->ops->nvscaleaddmulti = tf ? nvScaleAddMulti_OpenMP : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableDotProdMulti_OpenMP(N_Vector v, sunbooleantype tf)
{
  v->ops->nvdotprodmulti      = tf ? nvDotProdMulti_OpenMP : NULL;
  v->ops->nvdotprodmultilocal = tf ? nvDotProdMulti_OpenMP : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearSumVectorArray_OpenMP(N_Vector v, sunbooleantype tf)
{
  v->ops->nvlinearsumvectorarray = tf ? nvLinearSumVectorArray_OpenMP : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleVectorArray_OpenMP(N_Vector v, sunbooleantype tf)
{
  v->ops->nvscalevectorarray = tf ? nvScaleVectorArray_OpenMP : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableConstVectorArray_OpenMP(N_Vector v, sunbooleantype tf)
{
  v->ops->nvconstvectorarray = tf ? nvConstVectorArray_OpenMP : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableWrmsNormVectorArray_OpenMP(N_Vector v, sunbooleantype tf)
{
  v->ops->nvwrmsnormvectorarray = tf ? nvWrmsNormVectorArray_OpenMP : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableWrmsNormMaskVectorArray_OpenMP(N_Vector v, sunbooleantype tf)
{
  v->ops->nvwrmsnormmaskvectorarray = tf ? nvWrmsNormMaskVectorArray_OpenMP
                                         : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleAddMultiVectorArray_OpenMP(N_Vector v, sunbooleantype tf)
{
  v->ops->nvscaleaddmultivectorarray = tf ? nvScaleAddMultiVectorArray_OpenMP
                                          : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearCombinationVectorArray_OpenMP(N_Vector v,
                                                        sunbooleantype tf)
{
  v->ops->nvlinearcombinationvectorarray =
    tf ? nvLinearCombinationVectorArray_OpenMP : NULL;
  return SUN_SUCCESS;
}

/* Deprecated concrete operation wrappers */

void N_VAbs_OpenMP(N_Vector x, N_Vector z) { nvAbs_OpenMP(x, z); }

void N_VAddConst_OpenMP(N_Vector x, sunrealtype b, N_Vector z)
{
  nvAddConst_OpenMP(x, b, z);
}

SUNErrCode N_VBufPack_OpenMP(N_Vector x, void* buf)
{
  return nvBufPack_OpenMP(x, buf);
}

SUNErrCode N_VBufSize_OpenMP(N_Vector x, sunindextype* size)
{
  return nvBufSize_OpenMP(x, size);
}

SUNErrCode N_VBufUnpack_OpenMP(N_Vector x, void* buf)
{
  return nvBufUnpack_OpenMP(x, buf);
}

N_Vector N_VCloneEmpty_OpenMP(N_Vector w) { return nvCloneEmpty_OpenMP(w); }

N_Vector N_VClone_OpenMP(N_Vector w) { return nvClone_OpenMP(w); }

void N_VCompare_OpenMP(sunrealtype c, N_Vector x, N_Vector z)
{
  nvCompare_OpenMP(c, x, z);
}

SUNErrCode N_VConstVectorArray_OpenMP(int nvecs, sunrealtype c, N_Vector* Z)
{
  return nvConstVectorArray_OpenMP(nvecs, c, Z);
}

void N_VConst_OpenMP(sunrealtype c, N_Vector z) { nvConst_OpenMP(c, z); }

sunbooleantype N_VConstrMask_OpenMP(N_Vector c, N_Vector x, N_Vector m)
{
  return nvConstrMask_OpenMP(c, x, m);
}

void N_VDestroy_OpenMP(N_Vector v) { nvDestroy_OpenMP(v); }

void N_VDiv_OpenMP(N_Vector x, N_Vector y, N_Vector z)
{
  nvDiv_OpenMP(x, y, z);
}

SUNErrCode N_VDotProdMulti_OpenMP(int nvec, N_Vector x, N_Vector* Y,
                                  sunrealtype* dotprods)
{
  return nvDotProdMulti_OpenMP(nvec, x, Y, dotprods);
}

sunrealtype N_VDotProd_OpenMP(N_Vector x, N_Vector y)
{
  return nvDotProd_OpenMP(x, y);
}

sunrealtype* N_VGetArrayPointer_OpenMP(N_Vector v)
{
  return nvGetArrayPointer_OpenMP(v);
}

sunindextype N_VGetLength_OpenMP(N_Vector v) { return nvGetLength_OpenMP(v); }

N_Vector_ID N_VGetVectorID_OpenMP(N_Vector v)
{
  return nvGetVectorID_OpenMP(v);
}

sunbooleantype N_VInvTest_OpenMP(N_Vector x, N_Vector z)
{
  return nvInvTest_OpenMP(x, z);
}

void N_VInv_OpenMP(N_Vector x, N_Vector z) { nvInv_OpenMP(x, z); }

sunrealtype N_VL1Norm_OpenMP(N_Vector x) { return nvL1Norm_OpenMP(x); }

SUNErrCode N_VLinearCombinationVectorArray_OpenMP(int nvec, int nsum,
                                                  sunrealtype* c, N_Vector** X,
                                                  N_Vector* Z)
{
  return nvLinearCombinationVectorArray_OpenMP(nvec, nsum, c, X, Z);
}

SUNErrCode N_VLinearCombination_OpenMP(int nvec, sunrealtype* c, N_Vector* V,
                                       N_Vector z)
{
  return nvLinearCombination_OpenMP(nvec, c, V, z);
}

SUNErrCode N_VLinearSumVectorArray_OpenMP(int nvec, sunrealtype a, N_Vector* X,
                                          sunrealtype b, N_Vector* Y, N_Vector* Z)
{
  return nvLinearSumVectorArray_OpenMP(nvec, a, X, b, Y, Z);
}

void N_VLinearSum_OpenMP(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                         N_Vector z)
{
  nvLinearSum_OpenMP(a, x, b, y, z);
}

sunrealtype N_VMaxNorm_OpenMP(N_Vector x) { return nvMaxNorm_OpenMP(x); }

sunrealtype N_VMinQuotient_OpenMP(N_Vector num, N_Vector denom)
{
  return nvMinQuotient_OpenMP(num, denom);
}

sunrealtype N_VMin_OpenMP(N_Vector x) { return nvMin_OpenMP(x); }

void N_VPrintFile_OpenMP(N_Vector v, FILE* outfile)
{
  nvPrintFile_OpenMP(v, outfile);
}

void N_VPrint_OpenMP(N_Vector v) { nvPrint_OpenMP(v); }

void N_VProd_OpenMP(N_Vector x, N_Vector y, N_Vector z)
{
  nvProd_OpenMP(x, y, z);
}

SUNErrCode N_VScaleAddMultiVectorArray_OpenMP(int nvec, int nsum,
                                              sunrealtype* a, N_Vector* X,
                                              N_Vector** Y, N_Vector** Z)
{
  return nvScaleAddMultiVectorArray_OpenMP(nvec, nsum, a, X, Y, Z);
}

SUNErrCode N_VScaleAddMulti_OpenMP(int nvec, sunrealtype* a, N_Vector x,
                                   N_Vector* Y, N_Vector* Z)
{
  return nvScaleAddMulti_OpenMP(nvec, a, x, Y, Z);
}

SUNErrCode N_VScaleVectorArray_OpenMP(int nvec, sunrealtype* c, N_Vector* X,
                                      N_Vector* Z)
{
  return nvScaleVectorArray_OpenMP(nvec, c, X, Z);
}

void N_VScale_OpenMP(sunrealtype c, N_Vector x, N_Vector z)
{
  nvScale_OpenMP(c, x, z);
}

void N_VSetArrayPointer_OpenMP(sunrealtype* v_data, N_Vector v)
{
  nvSetArrayPointer_OpenMP(v_data, v);
}

sunrealtype N_VWL2Norm_OpenMP(N_Vector x, N_Vector w)
{
  return nvWL2Norm_OpenMP(x, w);
}

sunrealtype N_VWSqrSumLocal_OpenMP(N_Vector x, N_Vector w)
{
  return nvWSqrSumLocal_OpenMP(x, w);
}

sunrealtype N_VWSqrSumMaskLocal_OpenMP(N_Vector x, N_Vector w, N_Vector id)
{
  return nvWSqrSumMaskLocal_OpenMP(x, w, id);
}

SUNErrCode N_VWrmsNormMaskVectorArray_OpenMP(int nvecs, N_Vector* X, N_Vector* W,
                                             N_Vector id, sunrealtype* nrm)
{
  return nvWrmsNormMaskVectorArray_OpenMP(nvecs, X, W, id, nrm);
}

sunrealtype N_VWrmsNormMask_OpenMP(N_Vector x, N_Vector w, N_Vector id)
{
  return nvWrmsNormMask_OpenMP(x, w, id);
}

SUNErrCode N_VWrmsNormVectorArray_OpenMP(int nvecs, N_Vector* X, N_Vector* W,
                                         sunrealtype* nrm)
{
  return nvWrmsNormVectorArray_OpenMP(nvecs, X, W, nrm);
}

sunrealtype N_VWrmsNorm_OpenMP(N_Vector x, N_Vector w)
{
  return nvWrmsNorm_OpenMP(x, w);
}
