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
 * This is the implementation file for a serial implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_serial_deprecated.h>
#include <sundials/priv/sundials_context_impl.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_errors.h>

#include "sundials_macros.h"

#define ZERO   SUN_RCONST(0.0)
#define HALF   SUN_RCONST(0.5)
#define ONE    SUN_RCONST(1.0)
#define ONEPT5 SUN_RCONST(1.5)

/* Functions attached to the N_Vector */
static void nvAbs_Serial(N_Vector x, N_Vector z);
static void nvAddConst_Serial(N_Vector x, sunrealtype b, N_Vector z);
static SUNErrCode nvBufPack_Serial(N_Vector x, void* buf);
static SUNErrCode nvBufSize_Serial(N_Vector x, sunindextype* size);
static SUNErrCode nvBufUnpack_Serial(N_Vector x, void* buf);
static N_Vector nvCloneEmpty_Serial(N_Vector w);
static N_Vector nvClone_Serial(N_Vector w);
static void nvCompare_Serial(sunrealtype c, N_Vector x, N_Vector z);
static SUNErrCode nvConstVectorArray_Serial(int nvecs, sunrealtype c,
                                            N_Vector* Z);
static void nvConst_Serial(sunrealtype c, N_Vector z);
static sunbooleantype nvConstrMask_Serial(N_Vector c, N_Vector x, N_Vector m);
static void nvDestroy_Serial(N_Vector v);
static void nvDiv_Serial(N_Vector x, N_Vector y, N_Vector z);
static SUNErrCode nvDotProdMulti_Serial(int nvec, N_Vector x, N_Vector* Y,
                                        sunrealtype* dotprods);
static sunrealtype nvDotProd_Serial(N_Vector x, N_Vector y);
static sunrealtype* nvGetArrayPointer_Serial(N_Vector v);
static sunindextype nvGetLength_Serial(N_Vector v);
static N_Vector_ID nvGetVectorID_Serial(N_Vector v);
static sunbooleantype nvInvTest_Serial(N_Vector x, N_Vector z);
static void nvInv_Serial(N_Vector x, N_Vector z);
static sunrealtype nvL1Norm_Serial(N_Vector x);
static SUNErrCode nvLinearCombinationVectorArray_Serial(int nvec, int nsum,
                                                        sunrealtype* c,
                                                        N_Vector** X,
                                                        N_Vector* Z);
static SUNErrCode nvLinearCombination_Serial(int nvec, sunrealtype* c,
                                             N_Vector* V, N_Vector z);
static SUNErrCode nvLinearSumVectorArray_Serial(int nvec, sunrealtype a,
                                                N_Vector* X, sunrealtype b,
                                                N_Vector* Y, N_Vector* Z);
static void nvLinearSum_Serial(sunrealtype a, N_Vector x, sunrealtype b,
                               N_Vector y, N_Vector z);
static sunrealtype nvMaxNorm_Serial(N_Vector x);
static sunrealtype nvMinQuotient_Serial(N_Vector num, N_Vector denom);
static sunrealtype nvMin_Serial(N_Vector x);
static void nvPrintFile_Serial(N_Vector v, FILE* outfile);
static void nvPrint_Serial(N_Vector v);
static void nvProd_Serial(N_Vector x, N_Vector y, N_Vector z);
static SUNErrCode nvScaleAddMultiVectorArray_Serial(int nvec, int nsum,
                                                    sunrealtype* a, N_Vector* X,
                                                    N_Vector** Y, N_Vector** Z);
static SUNErrCode nvScaleAddMulti_Serial(int nvec, sunrealtype* a, N_Vector x,
                                         N_Vector* Y, N_Vector* Z);
static SUNErrCode nvScaleVectorArray_Serial(int nvec, sunrealtype* c,
                                            N_Vector* X, N_Vector* Z);
static void nvScale_Serial(sunrealtype c, N_Vector x, N_Vector z);
static void nvSetArrayPointer_Serial(sunrealtype* v_data, N_Vector v);
static sunrealtype nvWL2Norm_Serial(N_Vector x, N_Vector w);
static sunrealtype nvWSqrSumLocal_Serial(N_Vector x, N_Vector w);
static sunrealtype nvWSqrSumMaskLocal_Serial(N_Vector x, N_Vector w, N_Vector id);
static SUNErrCode nvWrmsNormMaskVectorArray_Serial(int nvecs, N_Vector* X,
                                                   N_Vector* W, N_Vector id,
                                                   sunrealtype* nrm);
static sunrealtype nvWrmsNormMask_Serial(N_Vector x, N_Vector w, N_Vector id);
static SUNErrCode nvWrmsNormVectorArray_Serial(int nvecs, N_Vector* X,
                                               N_Vector* W, sunrealtype* nrm);
static sunrealtype nvWrmsNorm_Serial(N_Vector x, N_Vector w);

/* Private functions for special cases of vector operations */
static void VCopy_Serial(N_Vector x, N_Vector z);             /* z=x       */
static void VSum_Serial(N_Vector x, N_Vector y, N_Vector z);  /* z=x+y     */
static void VDiff_Serial(N_Vector x, N_Vector y, N_Vector z); /* z=x-y     */
static void VNeg_Serial(N_Vector x, N_Vector z);              /* z=-x      */
static void VScaleSum_Serial(sunrealtype c, N_Vector x, N_Vector y,
                             N_Vector z); /* z=c(x+y)  */
static void VScaleDiff_Serial(sunrealtype c, N_Vector x, N_Vector y,
                              N_Vector z); /* z=c(x-y)  */
static void VLin1_Serial(sunrealtype a, N_Vector x, N_Vector y,
                         N_Vector z); /* z=ax+y    */
static void VLin2_Serial(sunrealtype a, N_Vector x, N_Vector y,
                         N_Vector z);                            /* z=ax-y    */
static void Vaxpy_Serial(sunrealtype a, N_Vector x, N_Vector y); /* y <- ax+y */
static void VScaleBy_Serial(sunrealtype a, N_Vector x);          /* x <- ax   */

/* Private functions for special cases of vector array operations */
static void VSumVectorArray_Serial(int nvec, N_Vector* X, N_Vector* Y,
                                   N_Vector* Z); /* Z=X+Y     */
static void VDiffVectorArray_Serial(int nvec, N_Vector* X, N_Vector* Y,
                                    N_Vector* Z); /* Z=X-Y     */
static void VScaleSumVectorArray_Serial(int nvec, sunrealtype c, N_Vector* X,
                                        N_Vector* Y, N_Vector* Z); /* Z=c(X+Y)  */
static void VScaleDiffVectorArray_Serial(int nvec, sunrealtype c, N_Vector* X,
                                         N_Vector* Y, N_Vector* Z); /* Z=c(X-Y)  */
static void VLin1VectorArray_Serial(int nvec, sunrealtype a, N_Vector* X,
                                    N_Vector* Y, N_Vector* Z); /* Z=aX+Y    */
static void VLin2VectorArray_Serial(int nvec, sunrealtype a, N_Vector* X,
                                    N_Vector* Y, N_Vector* Z); /* Z=aX-Y    */
static void VaxpyVectorArray_Serial(int nvec, sunrealtype a, N_Vector* X,
                                    N_Vector* Y); /* Y <- aX+Y */

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new empty serial vector
 */

N_Vector N_VNewEmpty_Serial(sunindextype length, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  N_Vector v;
  N_VectorContent_Serial content;

  SUNAssertNull(length >= 0, SUN_ERR_ARG_OUTOFRANGE);

  /* Create an empty vector object */
  v = NULL;
  v = N_VNewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = nvGetVectorID_Serial;
  v->ops->nvclone           = nvClone_Serial;
  v->ops->nvcloneempty      = nvCloneEmpty_Serial;
  v->ops->nvdestroy         = nvDestroy_Serial;
  v->ops->nvgetarraypointer = nvGetArrayPointer_Serial;
  v->ops->nvsetarraypointer = nvSetArrayPointer_Serial;
  v->ops->nvgetlength       = nvGetLength_Serial;
  v->ops->nvgetlocallength  = nvGetLength_Serial;

  /* standard vector operations */
  v->ops->nvlinearsum    = nvLinearSum_Serial;
  v->ops->nvconst        = nvConst_Serial;
  v->ops->nvprod         = nvProd_Serial;
  v->ops->nvdiv          = nvDiv_Serial;
  v->ops->nvscale        = nvScale_Serial;
  v->ops->nvabs          = nvAbs_Serial;
  v->ops->nvinv          = nvInv_Serial;
  v->ops->nvaddconst     = nvAddConst_Serial;
  v->ops->nvdotprod      = nvDotProd_Serial;
  v->ops->nvmaxnorm      = nvMaxNorm_Serial;
  v->ops->nvwrmsnormmask = nvWrmsNormMask_Serial;
  v->ops->nvwrmsnorm     = nvWrmsNorm_Serial;
  v->ops->nvmin          = nvMin_Serial;
  v->ops->nvwl2norm      = nvWL2Norm_Serial;
  v->ops->nvl1norm       = nvL1Norm_Serial;
  v->ops->nvcompare      = nvCompare_Serial;
  v->ops->nvinvtest      = nvInvTest_Serial;
  v->ops->nvconstrmask   = nvConstrMask_Serial;
  v->ops->nvminquotient  = nvMinQuotient_Serial;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvdotprodlocal     = nvDotProd_Serial;
  v->ops->nvmaxnormlocal     = nvMaxNorm_Serial;
  v->ops->nvminlocal         = nvMin_Serial;
  v->ops->nvl1normlocal      = nvL1Norm_Serial;
  v->ops->nvinvtestlocal     = nvInvTest_Serial;
  v->ops->nvconstrmasklocal  = nvConstrMask_Serial;
  v->ops->nvminquotientlocal = nvMinQuotient_Serial;
  v->ops->nvwsqrsumlocal     = nvWSqrSumLocal_Serial;
  v->ops->nvwsqrsummasklocal = nvWSqrSumMaskLocal_Serial;

  /* single buffer reduction operations */
  v->ops->nvdotprodmultilocal = nvDotProdMulti_Serial;

  /* XBraid interface operations */
  v->ops->nvbufsize   = nvBufSize_Serial;
  v->ops->nvbufpack   = nvBufPack_Serial;
  v->ops->nvbufunpack = nvBufUnpack_Serial;

  /* debugging functions */
  v->ops->nvprint     = nvPrint_Serial;
  v->ops->nvprintfile = nvPrintFile_Serial;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Serial)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  v->content = content;

  /* Initialize content */
  content->length   = length;
  content->own_data = SUNFALSE;
  content->data     = NULL;

  return (v);
}

/* ----------------------------------------------------------------------------
 * Function to create a new serial vector
 */

N_Vector N_VNew_Serial(sunindextype length, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  N_Vector v;
  sunrealtype* data;

  SUNAssertNull(length >= 0, SUN_ERR_ARG_OUTOFRANGE);

  v = NULL;
  v = N_VNewEmpty_Serial(length, sunctx);
  SUNCheckLastErrNull();

  /* Create data */
  data = NULL;
  if (length > 0)
  {
    data = (sunrealtype*)malloc(length * sizeof(sunrealtype));
    SUNAssertNull(data, SUN_ERR_MALLOC_FAIL);
  }

  /* Attach data */
  NV_OWN_DATA_S(v) = SUNTRUE;
  NV_DATA_S(v)     = data;

  return (v);
}

/* ----------------------------------------------------------------------------
 * Function to create a serial N_Vector with user data component
 */

N_Vector N_VMake_Serial(sunindextype length, sunrealtype* v_data,
                        SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  N_Vector v;

  SUNAssertNull(length >= 0, SUN_ERR_ARG_OUTOFRANGE);

  v = NULL;
  v = N_VNewEmpty_Serial(length, sunctx);
  SUNCheckLastErrNull();

  if (length > 0)
  {
    /* Attach data */
    NV_OWN_DATA_S(v) = SUNFALSE;
    NV_DATA_S(v)     = v_data;
  }

  return (v);
}

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID nvGetVectorID_Serial(SUNDIALS_MAYBE_UNUSED N_Vector v)
{
  return SUNDIALS_NVEC_SERIAL;
}

/* ----------------------------------------------------------------------------
 * Function to return number of vector elements
 */
sunindextype nvGetLength_Serial(N_Vector v) { return NV_LENGTH_S(v); }

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector to stdout
 */

void nvPrint_Serial(N_Vector x)
{
  SUNFunctionBegin(x->sunctx);
  nvPrintFile_Serial(x, stdout);
  SUNCheckLastErrVoid();
}

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector to outfile
 */

void nvPrintFile_Serial(N_Vector x, FILE* outfile)
{
  sunindextype i, N;
  sunrealtype* xd;

  xd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);

  for (i = 0; i < N; i++) { fprintf(outfile, SUN_FORMAT_E "\n", xd[i]); }

  return;
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector nvCloneEmpty_Serial(N_Vector w)
{
  SUNFunctionBegin(w->sunctx);
  N_Vector v;
  N_VectorContent_Serial content;

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty(w->sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  SUNCheckCallNull(N_VCopyOps(w, v));

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Serial)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  v->content = content;

  /* Initialize content */
  content->length   = NV_LENGTH_S(w);
  content->own_data = SUNFALSE;
  content->data     = NULL;

  return (v);
}

N_Vector nvClone_Serial(N_Vector w)
{
  SUNFunctionBegin(w->sunctx);
  N_Vector v;
  sunrealtype* data;
  sunindextype length;

  v = NULL;
  v = nvCloneEmpty_Serial(w);
  SUNCheckLastErrNull();

  length = NV_LENGTH_S(w);

  /* Create data */
  data = NULL;
  if (length > 0)
  {
    data = (sunrealtype*)malloc(length * sizeof(sunrealtype));
    SUNAssertNull(data, SUN_ERR_MALLOC_FAIL);

    /* Attach data */
    NV_OWN_DATA_S(v) = SUNTRUE;
    NV_DATA_S(v)     = data;
  }

  return (v);
}

void nvDestroy_Serial(N_Vector v)
{
  if (v == NULL) { return; }

  /* free content */
  if (v->content != NULL)
  {
    /* free data array if it's owned by the vector */
    if (NV_OWN_DATA_S(v) && NV_DATA_S(v) != NULL)
    {
      free(NV_DATA_S(v));
      NV_DATA_S(v) = NULL;
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

sunrealtype* nvGetArrayPointer_Serial(N_Vector v)
{
  return ((sunrealtype*)NV_DATA_S(v));
}

void nvSetArrayPointer_Serial(sunrealtype* v_data, N_Vector v)
{
  if (NV_LENGTH_S(v) > 0) { NV_DATA_S(v) = v_data; }

  return;
}

void nvLinearSum_Serial(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                        N_Vector z)
{
  sunindextype i, N;
  sunrealtype c, *xd, *yd, *zd;
  N_Vector v1, v2;
  sunbooleantype test;

  xd = yd = zd = NULL;

  if ((b == ONE) && (z == y))
  { /* BLAS usage: axpy y <- ax+y */
    Vaxpy_Serial(a, x, y);
    return;
  }

  if ((a == ONE) && (z == x))
  { /* BLAS usage: axpy x <- by+x */
    Vaxpy_Serial(b, y, x);
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE))
  {
    VSum_Serial(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE)))
  {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff_Serial(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE))
  {
    c  = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1_Serial(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE))
  {
    c  = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2_Serial(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b)
  {
    VScaleSum_Serial(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b)
  {
    VScaleDiff_Serial(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) { zd[i] = (a * xd[i]) + (b * yd[i]); }

  return;
}

void nvConst_Serial(sunrealtype c, N_Vector z)
{
  sunindextype i, N;
  sunrealtype* zd;

  zd = NULL;

  N  = NV_LENGTH_S(z);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) { zd[i] = c; }

  return;
}

void nvProd_Serial(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i] * yd[i]; }

  return;
}

void nvDiv_Serial(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i] / yd[i]; }

  return;
}

void nvScale_Serial(sunrealtype c, N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  xd = zd = NULL;

  if (z == x)
  { /* BLAS usage: scale x <- cx */
    VScaleBy_Serial(c, x);
    return;
  }

  if (c == ONE) { VCopy_Serial(x, z); }
  else if (c == -ONE) { VNeg_Serial(x, z); }
  else
  {
    N  = NV_LENGTH_S(x);
    xd = NV_DATA_S(x);
    zd = NV_DATA_S(z);
    for (i = 0; i < N; i++) { zd[i] = c * xd[i]; }
  }

  return;
}

void nvAbs_Serial(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) { zd[i] = SUNRabs(xd[i]); }

  return;
}

void nvInv_Serial(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) { zd[i] = ONE / xd[i]; }

  return;
}

void nvAddConst_Serial(N_Vector x, sunrealtype b, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i] + b; }

  return;
}

sunrealtype nvDotProd_Serial(N_Vector x, N_Vector y)
{
  sunindextype i, N;
  sunrealtype sum, *xd, *yd;

  sum = ZERO;
  xd = yd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);

  for (i = 0; i < N; i++) { sum += xd[i] * yd[i]; }

  return (sum);
}

sunrealtype nvMaxNorm_Serial(N_Vector x)
{
  sunindextype i, N;
  sunrealtype max, *xd;

  max = ZERO;
  xd  = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);

  for (i = 0; i < N; i++)
  {
    if (SUNRabs(xd[i]) > max) { max = SUNRabs(xd[i]); }
  }

  return (max);
}

sunrealtype nvWrmsNorm_Serial(N_Vector x, N_Vector w)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype norm = nvWSqrSumLocal_Serial(x, w);
  SUNCheckLastErrNoRet();
  norm = SUNRsqrt(norm / NV_LENGTH_S(x));
  return norm;
}

sunrealtype nvWSqrSumLocal_Serial(N_Vector x, N_Vector w)
{
  sunindextype i, N;
  sunrealtype sum, prodi, *xd, *wd;

  sum = ZERO;
  xd = wd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  wd = NV_DATA_S(w);

  for (i = 0; i < N; i++)
  {
    prodi = xd[i] * wd[i];
    sum += SUNSQR(prodi);
  }

  return (sum);
}

sunrealtype nvWrmsNormMask_Serial(N_Vector x, N_Vector w, N_Vector id)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype norm = nvWSqrSumMaskLocal_Serial(x, w, id);
  SUNCheckLastErrNoRet();
  norm = SUNRsqrt(norm / NV_LENGTH_S(x));
  return norm;
}

sunrealtype nvWSqrSumMaskLocal_Serial(N_Vector x, N_Vector w, N_Vector id)
{
  sunindextype i, N;
  sunrealtype sum, prodi, *xd, *wd, *idd;

  sum = ZERO;
  xd = wd = idd = NULL;

  N   = NV_LENGTH_S(x);
  xd  = NV_DATA_S(x);
  wd  = NV_DATA_S(w);
  idd = NV_DATA_S(id);

  for (i = 0; i < N; i++)
  {
    if (idd[i] > ZERO)
    {
      prodi = xd[i] * wd[i];
      sum += SUNSQR(prodi);
    }
  }

  return (sum);
}

sunrealtype nvMin_Serial(N_Vector x)
{
  sunindextype i, N;
  sunrealtype min, *xd;

  xd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);

  min = xd[0];

  for (i = 1; i < N; i++)
  {
    if (xd[i] < min) { min = xd[i]; }
  }

  return (min);
}

sunrealtype nvWL2Norm_Serial(N_Vector x, N_Vector w)
{
  sunindextype i, N;
  sunrealtype sum, prodi, *xd, *wd;

  sum = ZERO;
  xd = wd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  wd = NV_DATA_S(w);

  for (i = 0; i < N; i++)
  {
    prodi = xd[i] * wd[i];
    sum += SUNSQR(prodi);
  }

  return (SUNRsqrt(sum));
}

sunrealtype nvL1Norm_Serial(N_Vector x)
{
  sunindextype i, N;
  sunrealtype sum, *xd;

  sum = ZERO;
  xd  = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);

  for (i = 0; i < N; i++) { sum += SUNRabs(xd[i]); }

  return (sum);
}

void nvCompare_Serial(sunrealtype c, N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) { zd[i] = (SUNRabs(xd[i]) >= c) ? ONE : ZERO; }

  return;
}

sunbooleantype nvInvTest_Serial(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;
  sunbooleantype no_zero_found;

  xd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  no_zero_found = SUNTRUE;
  for (i = 0; i < N; i++)
  {
    if (xd[i] == ZERO) { no_zero_found = SUNFALSE; }
    else { zd[i] = ONE / xd[i]; }
  }

  return no_zero_found;
}

sunbooleantype nvConstrMask_Serial(N_Vector c, N_Vector x, N_Vector m)
{
  sunindextype i, N;
  sunrealtype temp;
  sunrealtype *cd, *xd, *md;
  sunbooleantype test;

  cd = xd = md = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  cd = NV_DATA_S(c);
  md = NV_DATA_S(m);

  temp = ZERO;

  for (i = 0; i < N; i++)
  {
    md[i] = ZERO;

    /* Continue if no constraints were set for the variable */
    if (cd[i] == ZERO) { continue; }

    /* Check if a set constraint has been violated */
    test = (SUNRabs(cd[i]) > ONEPT5 && xd[i] * cd[i] <= ZERO) ||
           (SUNRabs(cd[i]) > HALF && xd[i] * cd[i] < ZERO);
    if (test) { temp = md[i] = ONE; }
  }

  /* Return false if any constraint was violated */
  return (temp == ONE) ? SUNFALSE : SUNTRUE;
}

sunrealtype nvMinQuotient_Serial(N_Vector num, N_Vector denom)
{
  sunbooleantype notEvenOnce;
  sunindextype i, N;
  sunrealtype *nd, *dd, min;

  nd = dd = NULL;

  N  = NV_LENGTH_S(num);
  nd = NV_DATA_S(num);
  dd = NV_DATA_S(denom);

  notEvenOnce = SUNTRUE;
  min         = SUN_BIG_REAL;

  for (i = 0; i < N; i++)
  {
    if (dd[i] == ZERO) { continue; }
    else
    {
      if (!notEvenOnce) { min = SUNMIN(min, nd[i] / dd[i]); }
      else
      {
        min         = nd[i] / dd[i];
        notEvenOnce = SUNFALSE;
      }
    }
  }

  return (min);
}

/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */

SUNErrCode nvLinearCombination_Serial(int nvec, sunrealtype* c, N_Vector* X,
                                      N_Vector z)
{
  SUNFunctionBegin(X[0]->sunctx);

  int i;
  sunindextype j, N;
  sunrealtype* zd = NULL;
  sunrealtype* xd = NULL;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VScale */
  if (nvec == 1)
  {
    nvScale_Serial(c[0], X[0], z);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* should have called N_VLinearSum */
  if (nvec == 2)
  {
    nvLinearSum_Serial(c[0], X[0], c[1], X[1], z);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N  = NV_LENGTH_S(z);
  zd = NV_DATA_S(z);

  /*
   * X[0] += c[i]*X[i], i = 1,...,nvec-1
   */
  if ((X[0] == z) && (c[0] == ONE))
  {
    for (i = 1; i < nvec; i++)
    {
      xd = NV_DATA_S(X[i]);
      for (j = 0; j < N; j++) { zd[j] += c[i] * xd[j]; }
    }
    return SUN_SUCCESS;
  }

  /*
   * X[0] = c[0] * X[0] + sum{ c[i] * X[i] }, i = 1,...,nvec-1
   */
  if (X[0] == z)
  {
    for (j = 0; j < N; j++) { zd[j] *= c[0]; }
    for (i = 1; i < nvec; i++)
    {
      xd = NV_DATA_S(X[i]);
      for (j = 0; j < N; j++) { zd[j] += c[i] * xd[j]; }
    }
    return SUN_SUCCESS;
  }

  /*
   * z = sum{ c[i] * X[i] }, i = 0,...,nvec-1
   */
  xd = NV_DATA_S(X[0]);
  for (j = 0; j < N; j++) { zd[j] = c[0] * xd[j]; }
  for (i = 1; i < nvec; i++)
  {
    xd = NV_DATA_S(X[i]);
    for (j = 0; j < N; j++) { zd[j] += c[i] * xd[j]; }
  }
  return SUN_SUCCESS;
}

SUNErrCode nvScaleAddMulti_Serial(int nvec, sunrealtype* a, N_Vector x,
                                  N_Vector* Y, N_Vector* Z)
{
  SUNFunctionBegin(x->sunctx);
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VLinearSum */
  if (nvec == 1)
  {
    nvLinearSum_Serial(a[0], x, ONE, Y[0], Z[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z)
  {
    for (i = 0; i < nvec; i++)
    {
      yd = NV_DATA_S(Y[i]);
      for (j = 0; j < N; j++) { yd[j] += a[i] * xd[j]; }
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i = 0; i < nvec; i++)
  {
    yd = NV_DATA_S(Y[i]);
    zd = NV_DATA_S(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = a[i] * xd[j] + yd[j]; }
  }
  return SUN_SUCCESS;
}

SUNErrCode nvDotProdMulti_Serial(int nvec, N_Vector x, N_Vector* Y,
                                 sunrealtype* dotprods)
{
  SUNFunctionBegin(x->sunctx);
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VDotProd */
  if (nvec == 1)
  {
    dotprods[0] = nvDotProd_Serial(x, Y[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);

  /* compute multiple dot products */
  for (i = 0; i < nvec; i++)
  {
    yd          = NV_DATA_S(Y[i]);
    dotprods[i] = ZERO;
    for (j = 0; j < N; j++) { dotprods[i] += xd[j] * yd[j]; }
  }

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------
 */

SUNErrCode nvLinearSumVectorArray_Serial(int nvec, sunrealtype a, N_Vector* X,
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

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VLinearSum */
  if (nvec == 1)
  {
    nvLinearSum_Serial(a, X[0], b, Y[0], Z[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* BLAS usage: axpy y <- ax+y */
  if ((b == ONE) && (Z == Y))
  {
    VaxpyVectorArray_Serial(nvec, a, X, Y);
    return SUN_SUCCESS;
  }

  /* BLAS usage: axpy x <- by+x */
  if ((a == ONE) && (Z == X))
  {
    VaxpyVectorArray_Serial(nvec, b, Y, X);
    return SUN_SUCCESS;
  }

  /* Case: a == b == 1.0 */
  if ((a == ONE) && (b == ONE))
  {
    VSumVectorArray_Serial(nvec, X, Y, Z);
    return SUN_SUCCESS;
  }

  /* Cases:                    */
  /*   (1) a == 1.0, b = -1.0, */
  /*   (2) a == -1.0, b == 1.0 */
  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE)))
  {
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    VDiffVectorArray_Serial(nvec, V2, V1, Z);
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
    VLin1VectorArray_Serial(nvec, c, V1, V2, Z);
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
    VLin2VectorArray_Serial(nvec, c, V1, V2, Z);
    return SUN_SUCCESS;
  }

  /* Case: a == b                                                         */
  /* catches case both a and b are 0.0 - user should have called N_VConst */
  if (a == b)
  {
    VScaleSumVectorArray_Serial(nvec, a, X, Y, Z);
    return SUN_SUCCESS;
  }

  /* Case: a == -b */
  if (a == -b)
  {
    VScaleDiffVectorArray_Serial(nvec, a, X, Y, Z);
    return SUN_SUCCESS;
  }

  /* Do all cases not handled above:                               */
  /*   (1) a == other, b == 0.0 - user should have called N_VScale */
  /*   (2) a == 0.0, b == other - user should have called N_VScale */
  /*   (3) a,b == other, a !=b, a != -b                            */

  /* get vector length */
  N = NV_LENGTH_S(Z[0]);

  /* compute linear sum for each vector pair in vector arrays */
  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_S(X[i]);
    yd = NV_DATA_S(Y[i]);
    zd = NV_DATA_S(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = a * xd[j] + b * yd[j]; }
  }

  return SUN_SUCCESS;
}

SUNErrCode nvScaleVectorArray_Serial(int nvec, sunrealtype* c, N_Vector* X,
                                     N_Vector* Z)
{
  SUNFunctionBegin(X[0]->sunctx);
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* zd = NULL;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VScale */
  if (nvec == 1)
  {
    nvScale_Serial(c[0], X[0], Z[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length */
  N = NV_LENGTH_S(Z[0]);

  /*
   * X[i] *= c[i]
   */
  if (X == Z)
  {
    for (i = 0; i < nvec; i++)
    {
      xd = NV_DATA_S(X[i]);
      for (j = 0; j < N; j++) { xd[j] *= c[i]; }
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[i] = c[i] * X[i]
   */
  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_S(X[i]);
    zd = NV_DATA_S(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = c[i] * xd[j]; }
  }
  return SUN_SUCCESS;
}

SUNErrCode nvConstVectorArray_Serial(int nvec, sunrealtype c, N_Vector* Z)
{
  SUNFunctionBegin(Z[0]->sunctx);
  int i;
  sunindextype j, N;
  sunrealtype* zd = NULL;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VConst */
  if (nvec == 1)
  {
    nvConst_Serial(c, Z[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length */
  N = NV_LENGTH_S(Z[0]);

  /* set each vector in the vector array to a constant */
  for (i = 0; i < nvec; i++)
  {
    zd = NV_DATA_S(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = c; }
  }

  return SUN_SUCCESS;
}

SUNErrCode nvWrmsNormVectorArray_Serial(int nvec, N_Vector* X, N_Vector* W,
                                        sunrealtype* nrm)
{
  SUNFunctionBegin(X[0]->sunctx);
  int i;
  sunindextype j, N;
  sunrealtype* wd = NULL;
  sunrealtype* xd = NULL;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VWrmsNorm */
  if (nvec == 1)
  {
    nrm[0] = nvWrmsNorm_Serial(X[0], W[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length */
  N = NV_LENGTH_S(X[0]);

  /* compute the WRMS norm for each vector in the vector array */
  for (i = 0; i < nvec; i++)
  {
    xd     = NV_DATA_S(X[i]);
    wd     = NV_DATA_S(W[i]);
    nrm[i] = ZERO;
    for (j = 0; j < N; j++) { nrm[i] += SUNSQR(xd[j] * wd[j]); }
    nrm[i] = SUNRsqrt(nrm[i] / N);
  }

  return SUN_SUCCESS;
}

SUNErrCode nvWrmsNormMaskVectorArray_Serial(int nvec, N_Vector* X, N_Vector* W,
                                            N_Vector id, sunrealtype* nrm)
{
  SUNFunctionBegin(X[0]->sunctx);
  int i;
  sunindextype j, N;
  sunrealtype* wd  = NULL;
  sunrealtype* xd  = NULL;
  sunrealtype* idd = NULL;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VWrmsNorm */
  if (nvec == 1)
  {
    nrm[0] = nvWrmsNormMask_Serial(X[0], W[0], id);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length and mask data array */
  N   = NV_LENGTH_S(X[0]);
  idd = NV_DATA_S(id);

  /* compute the WRMS norm for each vector in the vector array */
  for (i = 0; i < nvec; i++)
  {
    xd     = NV_DATA_S(X[i]);
    wd     = NV_DATA_S(W[i]);
    nrm[i] = ZERO;
    for (j = 0; j < N; j++)
    {
      if (idd[j] > ZERO) { nrm[i] += SUNSQR(xd[j] * wd[j]); }
    }
    nrm[i] = SUNRsqrt(nrm[i] / N);
  }

  return SUN_SUCCESS;
}

SUNErrCode nvScaleAddMultiVectorArray_Serial(int nvec, int nsum, sunrealtype* a,
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
      nvLinearSum_Serial(a[0], X[0], ONE, Y[0][0], Z[0][0]);
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

    SUNCheckCall(nvScaleAddMulti_Serial(nsum, a, X[0], YY, ZZ));

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
    SUNCheckCall(nvLinearSumVectorArray_Serial(nvec, a[0], X, ONE, Y[0], Z[0]));
    return SUN_SUCCESS;
  }

  /* ----------------------------
   * Compute multiple linear sums
   * ---------------------------- */

  /* get vector length */
  N = NV_LENGTH_S(X[0]);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z)
  {
    for (i = 0; i < nvec; i++)
    {
      xd = NV_DATA_S(X[i]);
      for (j = 0; j < nsum; j++)
      {
        yd = NV_DATA_S(Y[j][i]);
        for (k = 0; k < N; k++) { yd[k] += a[j] * xd[k]; }
      }
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_S(X[i]);
    for (j = 0; j < nsum; j++)
    {
      yd = NV_DATA_S(Y[j][i]);
      zd = NV_DATA_S(Z[j][i]);
      for (k = 0; k < N; k++) { zd[k] = a[j] * xd[k] + yd[k]; }
    }
  }
  return SUN_SUCCESS;
}

SUNErrCode nvLinearCombinationVectorArray_Serial(int nvec, int nsum,
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
      nvScale_Serial(c[0], X[0][0], Z[0]);
      SUNCheckLastErr();
      return SUN_SUCCESS;
    }

    /* should have called N_VLinearSum */
    if (nsum == 2)
    {
      nvLinearSum_Serial(c[0], X[0][0], c[1], X[1][0], Z[0]);
      SUNCheckLastErr();
      return SUN_SUCCESS;
    }

    /* should have called N_VLinearCombination */
    Y = (N_Vector*)malloc(nsum * sizeof(N_Vector));
    SUNAssert(Y, SUN_ERR_MALLOC_FAIL);

    for (i = 0; i < nsum; i++) { Y[i] = X[i][0]; }

    SUNCheckCall(nvLinearCombination_Serial(nsum, c, Y, Z[0]));

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
    SUNAssert(ctmp, SUN_ERR_MALLOC_FAIL);

    for (j = 0; j < nvec; j++) { ctmp[j] = c[0]; }

    SUNCheckCall(nvScaleVectorArray_Serial(nvec, ctmp, X[0], Z));

    free(ctmp);
    return SUN_SUCCESS;
  }

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 2)
  {
    SUNCheckCall(nvLinearSumVectorArray_Serial(nvec, c[0], X[0], c[1], X[1], Z));
    return SUN_SUCCESS;
  }

  /* --------------------------
   * Compute linear combination
   * -------------------------- */

  /* get vector length */
  N = NV_LENGTH_S(Z[0]);

  /*
   * X[0][j] += c[i]*X[i][j], i = 1,...,nvec-1
   */
  if ((X[0] == Z) && (c[0] == ONE))
  {
    for (j = 0; j < nvec; j++)
    {
      zd = NV_DATA_S(Z[j]);
      for (i = 1; i < nsum; i++)
      {
        xd = NV_DATA_S(X[i][j]);
        for (k = 0; k < N; k++) { zd[k] += c[i] * xd[k]; }
      }
    }
    return SUN_SUCCESS;
  }

  /*
   * X[0][j] = c[0] * X[0][j] + sum{ c[i] * X[i][j] }, i = 1,...,nvec-1
   */
  if (X[0] == Z)
  {
    for (j = 0; j < nvec; j++)
    {
      zd = NV_DATA_S(Z[j]);
      for (k = 0; k < N; k++) { zd[k] *= c[0]; }
      for (i = 1; i < nsum; i++)
      {
        xd = NV_DATA_S(X[i][j]);
        for (k = 0; k < N; k++) { zd[k] += c[i] * xd[k]; }
      }
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[j] = sum{ c[i] * X[i][j] }, i = 0,...,nvec-1
   */
  for (j = 0; j < nvec; j++)
  {
    xd = NV_DATA_S(X[0][j]);
    zd = NV_DATA_S(Z[j]);
    for (k = 0; k < N; k++) { zd[k] = c[0] * xd[k]; }
    for (i = 1; i < nsum; i++)
    {
      xd = NV_DATA_S(X[i][j]);
      for (k = 0; k < N; k++) { zd[k] += c[i] * xd[k]; }
    }
  }
  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * OPTIONAL XBraid interface operations
 * -----------------------------------------------------------------
 */

SUNErrCode nvBufSize_Serial(N_Vector x, sunindextype* size)
{
  *size = NV_LENGTH_S(x) * ((sunindextype)sizeof(sunrealtype));
  return SUN_SUCCESS;
}

SUNErrCode nvBufPack_Serial(N_Vector x, void* buf)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i, N;
  sunrealtype* xd = NULL;
  sunrealtype* bd = NULL;

  SUNAssert(buf, SUN_ERR_ARG_CORRUPT);

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  bd = (sunrealtype*)buf;

  for (i = 0; i < N; i++) { bd[i] = xd[i]; }

  return SUN_SUCCESS;
}

SUNErrCode nvBufUnpack_Serial(N_Vector x, void* buf)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i, N;
  sunrealtype* xd = NULL;
  sunrealtype* bd = NULL;

  SUNAssert(buf, SUN_ERR_ARG_CORRUPT);

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  bd = (sunrealtype*)buf;

  for (i = 0; i < N; i++) { xd[i] = bd[i]; }

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * private functions for special cases of vector operations
 * -----------------------------------------------------------------
 */

static void VCopy_Serial(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i]; }

  return;
}

static void VSum_Serial(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i] + yd[i]; }

  return;
}

static void VDiff_Serial(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i] - yd[i]; }

  return;
}

static void VNeg_Serial(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) { zd[i] = -xd[i]; }

  return;
}

static void VScaleSum_Serial(sunrealtype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) { zd[i] = c * (xd[i] + yd[i]); }

  return;
}

static void VScaleDiff_Serial(sunrealtype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) { zd[i] = c * (xd[i] - yd[i]); }

  return;
}

static void VLin1_Serial(sunrealtype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) { zd[i] = (a * xd[i]) + yd[i]; }

  return;
}

static void VLin2_Serial(sunrealtype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) { zd[i] = (a * xd[i]) - yd[i]; }

  return;
}

static void Vaxpy_Serial(sunrealtype a, N_Vector x, N_Vector y)
{
  sunindextype i, N;
  sunrealtype *xd, *yd;

  xd = yd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);

  if (a == ONE)
  {
    for (i = 0; i < N; i++) { yd[i] += xd[i]; }
    return;
  }

  if (a == -ONE)
  {
    for (i = 0; i < N; i++) { yd[i] -= xd[i]; }
    return;
  }

  for (i = 0; i < N; i++) { yd[i] += a * xd[i]; }

  return;
}

static void VScaleBy_Serial(sunrealtype a, N_Vector x)
{
  sunindextype i, N;
  sunrealtype* xd;

  xd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);

  for (i = 0; i < N; i++) { xd[i] *= a; }

  return;
}

/*
 * -----------------------------------------------------------------
 * private functions for special cases of vector array operations
 * -----------------------------------------------------------------
 */

static void VSumVectorArray_Serial(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  N = NV_LENGTH_S(X[0]);

  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_S(X[i]);
    yd = NV_DATA_S(Y[i]);
    zd = NV_DATA_S(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = xd[j] + yd[j]; }
  }
}

static void VDiffVectorArray_Serial(int nvec, N_Vector* X, N_Vector* Y,
                                    N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  N = NV_LENGTH_S(X[0]);

  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_S(X[i]);
    yd = NV_DATA_S(Y[i]);
    zd = NV_DATA_S(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = xd[j] - yd[j]; }
  }
}

static void VScaleSumVectorArray_Serial(int nvec, sunrealtype c, N_Vector* X,
                                        N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  N = NV_LENGTH_S(X[0]);

  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_S(X[i]);
    yd = NV_DATA_S(Y[i]);
    zd = NV_DATA_S(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = c * (xd[j] + yd[j]); }
  }
}

static void VScaleDiffVectorArray_Serial(int nvec, sunrealtype c, N_Vector* X,
                                         N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  N = NV_LENGTH_S(X[0]);

  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_S(X[i]);
    yd = NV_DATA_S(Y[i]);
    zd = NV_DATA_S(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = c * (xd[j] - yd[j]); }
  }
}

static void VLin1VectorArray_Serial(int nvec, sunrealtype a, N_Vector* X,
                                    N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  N = NV_LENGTH_S(X[0]);

  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_S(X[i]);
    yd = NV_DATA_S(Y[i]);
    zd = NV_DATA_S(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = (a * xd[j]) + yd[j]; }
  }
}

static void VLin2VectorArray_Serial(int nvec, sunrealtype a, N_Vector* X,
                                    N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  N = NV_LENGTH_S(X[0]);

  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_S(X[i]);
    yd = NV_DATA_S(Y[i]);
    zd = NV_DATA_S(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = (a * xd[j]) - yd[j]; }
  }
}

static void VaxpyVectorArray_Serial(int nvec, sunrealtype a, N_Vector* X,
                                    N_Vector* Y)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;

  N = NV_LENGTH_S(X[0]);

  if (a == ONE)
  {
    for (i = 0; i < nvec; i++)
    {
      xd = NV_DATA_S(X[i]);
      yd = NV_DATA_S(Y[i]);
      for (j = 0; j < N; j++) { yd[j] += xd[j]; }
    }
    return;
  }

  if (a == -ONE)
  {
    for (i = 0; i < nvec; i++)
    {
      xd = NV_DATA_S(X[i]);
      yd = NV_DATA_S(Y[i]);
      for (j = 0; j < N; j++) { yd[j] -= xd[j]; }
    }
    return;
  }

  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_S(X[i]);
    yd = NV_DATA_S(Y[i]);
    for (j = 0; j < N; j++) { yd[j] += a * xd[j]; }
  }
}

/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

SUNErrCode N_VEnableFusedOps_Serial(N_Vector v, sunbooleantype tf)
{
  if (tf)
  {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = nvLinearCombination_Serial;
    v->ops->nvscaleaddmulti     = nvScaleAddMulti_Serial;
    v->ops->nvdotprodmulti      = nvDotProdMulti_Serial;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray     = nvLinearSumVectorArray_Serial;
    v->ops->nvscalevectorarray         = nvScaleVectorArray_Serial;
    v->ops->nvconstvectorarray         = nvConstVectorArray_Serial;
    v->ops->nvwrmsnormvectorarray      = nvWrmsNormVectorArray_Serial;
    v->ops->nvwrmsnormmaskvectorarray  = nvWrmsNormMaskVectorArray_Serial;
    v->ops->nvscaleaddmultivectorarray = nvScaleAddMultiVectorArray_Serial;
    v->ops->nvlinearcombinationvectorarray = nvLinearCombinationVectorArray_Serial;
    /* enable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = nvDotProdMulti_Serial;
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

SUNErrCode N_VEnableLinearCombination_Serial(N_Vector v, sunbooleantype tf)
{
  v->ops->nvlinearcombination = tf ? nvLinearCombination_Serial : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleAddMulti_Serial(N_Vector v, sunbooleantype tf)
{
  v->ops->nvscaleaddmulti = tf ? nvScaleAddMulti_Serial : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableDotProdMulti_Serial(N_Vector v, sunbooleantype tf)
{
  v->ops->nvdotprodmulti      = tf ? nvDotProdMulti_Serial : NULL;
  v->ops->nvdotprodmultilocal = tf ? nvDotProdMulti_Serial : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearSumVectorArray_Serial(N_Vector v, sunbooleantype tf)
{
  v->ops->nvlinearsumvectorarray = tf ? nvLinearSumVectorArray_Serial : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleVectorArray_Serial(N_Vector v, sunbooleantype tf)
{
  v->ops->nvscalevectorarray = tf ? nvScaleVectorArray_Serial : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableConstVectorArray_Serial(N_Vector v, sunbooleantype tf)
{
  v->ops->nvconstvectorarray = tf ? nvConstVectorArray_Serial : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableWrmsNormVectorArray_Serial(N_Vector v, sunbooleantype tf)
{
  v->ops->nvwrmsnormvectorarray = tf ? nvWrmsNormVectorArray_Serial : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableWrmsNormMaskVectorArray_Serial(N_Vector v, sunbooleantype tf)
{
  v->ops->nvwrmsnormmaskvectorarray = tf ? nvWrmsNormMaskVectorArray_Serial
                                         : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleAddMultiVectorArray_Serial(N_Vector v, sunbooleantype tf)
{
  v->ops->nvscaleaddmultivectorarray = tf ? nvScaleAddMultiVectorArray_Serial
                                          : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearCombinationVectorArray_Serial(N_Vector v,
                                                        sunbooleantype tf)
{
  v->ops->nvlinearcombinationvectorarray =
    tf ? nvLinearCombinationVectorArray_Serial : NULL;
  return SUN_SUCCESS;
}

/* Deprecated concrete operation wrappers */

void N_VAbs_Serial(N_Vector x, N_Vector z) { nvAbs_Serial(x, z); }

void N_VAddConst_Serial(N_Vector x, sunrealtype b, N_Vector z)
{
  nvAddConst_Serial(x, b, z);
}

SUNErrCode N_VBufPack_Serial(N_Vector x, void* buf)
{
  return nvBufPack_Serial(x, buf);
}

SUNErrCode N_VBufSize_Serial(N_Vector x, sunindextype* size)
{
  return nvBufSize_Serial(x, size);
}

SUNErrCode N_VBufUnpack_Serial(N_Vector x, void* buf)
{
  return nvBufUnpack_Serial(x, buf);
}

N_Vector N_VCloneEmpty_Serial(N_Vector w) { return nvCloneEmpty_Serial(w); }

N_Vector N_VClone_Serial(N_Vector w) { return nvClone_Serial(w); }

void N_VCompare_Serial(sunrealtype c, N_Vector x, N_Vector z)
{
  nvCompare_Serial(c, x, z);
}

SUNErrCode N_VConstVectorArray_Serial(int nvecs, sunrealtype c, N_Vector* Z)
{
  return nvConstVectorArray_Serial(nvecs, c, Z);
}

void N_VConst_Serial(sunrealtype c, N_Vector z) { nvConst_Serial(c, z); }

sunbooleantype N_VConstrMask_Serial(N_Vector c, N_Vector x, N_Vector m)
{
  return nvConstrMask_Serial(c, x, m);
}

void N_VDestroy_Serial(N_Vector v) { nvDestroy_Serial(v); }

void N_VDiv_Serial(N_Vector x, N_Vector y, N_Vector z)
{
  nvDiv_Serial(x, y, z);
}

SUNErrCode N_VDotProdMulti_Serial(int nvec, N_Vector x, N_Vector* Y,
                                  sunrealtype* dotprods)
{
  return nvDotProdMulti_Serial(nvec, x, Y, dotprods);
}

sunrealtype N_VDotProd_Serial(N_Vector x, N_Vector y)
{
  return nvDotProd_Serial(x, y);
}

sunrealtype* N_VGetArrayPointer_Serial(N_Vector v)
{
  return nvGetArrayPointer_Serial(v);
}

sunindextype N_VGetLength_Serial(N_Vector v) { return nvGetLength_Serial(v); }

N_Vector_ID N_VGetVectorID_Serial(N_Vector v)
{
  return nvGetVectorID_Serial(v);
}

sunbooleantype N_VInvTest_Serial(N_Vector x, N_Vector z)
{
  return nvInvTest_Serial(x, z);
}

void N_VInv_Serial(N_Vector x, N_Vector z) { nvInv_Serial(x, z); }

sunrealtype N_VL1Norm_Serial(N_Vector x) { return nvL1Norm_Serial(x); }

SUNErrCode N_VLinearCombinationVectorArray_Serial(int nvec, int nsum,
                                                  sunrealtype* c, N_Vector** X,
                                                  N_Vector* Z)
{
  return nvLinearCombinationVectorArray_Serial(nvec, nsum, c, X, Z);
}

SUNErrCode N_VLinearCombination_Serial(int nvec, sunrealtype* c, N_Vector* V,
                                       N_Vector z)
{
  return nvLinearCombination_Serial(nvec, c, V, z);
}

SUNErrCode N_VLinearSumVectorArray_Serial(int nvec, sunrealtype a, N_Vector* X,
                                          sunrealtype b, N_Vector* Y, N_Vector* Z)
{
  return nvLinearSumVectorArray_Serial(nvec, a, X, b, Y, Z);
}

void N_VLinearSum_Serial(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                         N_Vector z)
{
  nvLinearSum_Serial(a, x, b, y, z);
}

sunrealtype N_VMaxNorm_Serial(N_Vector x) { return nvMaxNorm_Serial(x); }

sunrealtype N_VMinQuotient_Serial(N_Vector num, N_Vector denom)
{
  return nvMinQuotient_Serial(num, denom);
}

sunrealtype N_VMin_Serial(N_Vector x) { return nvMin_Serial(x); }

void N_VPrintFile_Serial(N_Vector v, FILE* outfile)
{
  nvPrintFile_Serial(v, outfile);
}

void N_VPrint_Serial(N_Vector v) { nvPrint_Serial(v); }

void N_VProd_Serial(N_Vector x, N_Vector y, N_Vector z)
{
  nvProd_Serial(x, y, z);
}

SUNErrCode N_VScaleAddMultiVectorArray_Serial(int nvec, int nsum,
                                              sunrealtype* a, N_Vector* X,
                                              N_Vector** Y, N_Vector** Z)
{
  return nvScaleAddMultiVectorArray_Serial(nvec, nsum, a, X, Y, Z);
}

SUNErrCode N_VScaleAddMulti_Serial(int nvec, sunrealtype* a, N_Vector x,
                                   N_Vector* Y, N_Vector* Z)
{
  return nvScaleAddMulti_Serial(nvec, a, x, Y, Z);
}

SUNErrCode N_VScaleVectorArray_Serial(int nvec, sunrealtype* c, N_Vector* X,
                                      N_Vector* Z)
{
  return nvScaleVectorArray_Serial(nvec, c, X, Z);
}

void N_VScale_Serial(sunrealtype c, N_Vector x, N_Vector z)
{
  nvScale_Serial(c, x, z);
}

void N_VSetArrayPointer_Serial(sunrealtype* v_data, N_Vector v)
{
  nvSetArrayPointer_Serial(v_data, v);
}

sunrealtype N_VWL2Norm_Serial(N_Vector x, N_Vector w)
{
  return nvWL2Norm_Serial(x, w);
}

sunrealtype N_VWSqrSumLocal_Serial(N_Vector x, N_Vector w)
{
  return nvWSqrSumLocal_Serial(x, w);
}

sunrealtype N_VWSqrSumMaskLocal_Serial(N_Vector x, N_Vector w, N_Vector id)
{
  return nvWSqrSumMaskLocal_Serial(x, w, id);
}

SUNErrCode N_VWrmsNormMaskVectorArray_Serial(int nvecs, N_Vector* X, N_Vector* W,
                                             N_Vector id, sunrealtype* nrm)
{
  return nvWrmsNormMaskVectorArray_Serial(nvecs, X, W, id, nrm);
}

sunrealtype N_VWrmsNormMask_Serial(N_Vector x, N_Vector w, N_Vector id)
{
  return nvWrmsNormMask_Serial(x, w, id);
}

SUNErrCode N_VWrmsNormVectorArray_Serial(int nvecs, N_Vector* X, N_Vector* W,
                                         sunrealtype* nrm)
{
  return nvWrmsNormVectorArray_Serial(nvecs, X, W, nrm);
}

sunrealtype N_VWrmsNorm_Serial(N_Vector x, N_Vector w)
{
  return nvWrmsNorm_Serial(x, w);
}
