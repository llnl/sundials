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
 * This is the implementation file for a parallel MPI implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_parallel.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/priv/sundials_mpi_errors_impl.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_types.h>

#include "sundials_macros.h"

#define ZERO   SUN_RCONST(0.0)
#define HALF   SUN_RCONST(0.5)
#define ONE    SUN_RCONST(1.0)
#define ONEPT5 SUN_RCONST(1.5)

/* Private functions for special cases of vector operations */
static void VCopy_Parallel(N_Vector x, N_Vector z);             /* z=x       */
static void VSum_Parallel(N_Vector x, N_Vector y, N_Vector z);  /* z=x+y     */
static void VDiff_Parallel(N_Vector x, N_Vector y, N_Vector z); /* z=x-y     */
static void VNeg_Parallel(N_Vector x, N_Vector z);              /* z=-x      */
static void VScaleSum_Parallel(sunrealtype c, N_Vector x, N_Vector y,
                               N_Vector z); /* z=c(x+y)  */
static void VScaleDiff_Parallel(sunrealtype c, N_Vector x, N_Vector y,
                                N_Vector z); /* z=c(x-y)  */
static void VLin1_Parallel(sunrealtype a, N_Vector x, N_Vector y,
                           N_Vector z); /* z=ax+y    */
static void VLin2_Parallel(sunrealtype a, N_Vector x, N_Vector y,
                           N_Vector z); /* z=ax-y    */
static void Vaxpy_Parallel(sunrealtype a, N_Vector x, N_Vector y); /* y <- ax+y */
static void VScaleBy_Parallel(sunrealtype a, N_Vector x); /* x <- ax   */

/* Private functions for special cases of vector array operations */
static void VSumVectorArray_Parallel(int nvec, N_Vector* X, N_Vector* Y,
                                     N_Vector* Z); /* Z=X+Y     */
static void VDiffVectorArray_Parallel(int nvec, N_Vector* X, N_Vector* Y,
                                      N_Vector* Z); /* Z=X-Y     */
static void VScaleSumVectorArray_Parallel(int nvec, sunrealtype c, N_Vector* X,
                                          N_Vector* Y, N_Vector* Z); /* Z=c(X+Y)  */
static void VScaleDiffVectorArray_Parallel(int nvec, sunrealtype c, N_Vector* X,
                                           N_Vector* Y,
                                           N_Vector* Z); /* Z=c(X-Y)  */
static void VLin1VectorArray_Parallel(int nvec, sunrealtype a, N_Vector* X,
                                      N_Vector* Y, N_Vector* Z); /* Z=aX+Y    */
static void VLin2VectorArray_Parallel(int nvec, sunrealtype a, N_Vector* X,
                                      N_Vector* Y, N_Vector* Z); /* Z=aX-Y    */
static void VaxpyVectorArray_Parallel(int nvec, sunrealtype a, N_Vector* X,
                                      N_Vector* Y); /* Y <- aX+Y */

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------
 * Function to create a new parallel vector with empty data array
 */

N_Vector N_VNewEmpty_Parallel(MPI_Comm comm, sunindextype local_length,
                              sunindextype global_length, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  N_Vector v;
  N_VectorContent_Parallel content;
  sunindextype n, Nsum;

  SUNAssertNull(local_length >= 0, SUN_ERR_ARG_OUTOFRANGE);
  SUNAssertNull(global_length >= 0, SUN_ERR_ARG_OUTOFRANGE);

  /* Compute global length as sum of local lengths */
  n = local_length;
  SUNCheckMPICallNull(
    MPI_Allreduce(&n, &Nsum, 1, MPI_SUNINDEXTYPE, MPI_SUM, comm));
  if (Nsum != global_length)
  {
    SUNHandleErrWithMsg(__LINE__, __func__,
                        __FILE__, "global_length does not equal the computed global length",
                        SUN_ERR_ARG_INCOMPATIBLE, SUNCTX_);
    return NULL;
  }

  /* Create an empty vector object */
  v = NULL;
  v = N_VNewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = N_VGetVectorID_Parallel;
  v->ops->nvclone           = N_VClone_Parallel;
  v->ops->nvcloneempty      = N_VCloneEmpty_Parallel;
  v->ops->nvdestroy         = N_VDestroy_Parallel;
  v->ops->nvspace           = N_VSpace_Parallel;
  v->ops->nvgetarraypointer = N_VGetArrayPointer_Parallel;
  v->ops->nvsetarraypointer = N_VSetArrayPointer_Parallel;
  v->ops->nvgetcommunicator = N_VGetCommunicator_Parallel;
  v->ops->nvgetlength       = N_VGetLength_Parallel;
  v->ops->nvgetlocallength  = N_VGetLocalLength_Parallel;

  /* standard vector operations */
  v->ops->nvlinearsum    = N_VLinearSum_Parallel;
  v->ops->nvconst        = N_VConst_Parallel;
  v->ops->nvprod         = N_VProd_Parallel;
  v->ops->nvdiv          = N_VDiv_Parallel;
  v->ops->nvscale        = N_VScale_Parallel;
  v->ops->nvabs          = N_VAbs_Parallel;
  v->ops->nvinv          = N_VInv_Parallel;
  v->ops->nvaddconst     = N_VAddConst_Parallel;
  v->ops->nvdotprod      = N_VDotProd_Parallel;
  v->ops->nvmaxnorm      = N_VMaxNorm_Parallel;
  v->ops->nvwrmsnormmask = N_VWrmsNormMask_Parallel;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_Parallel;
  v->ops->nvmin          = N_VMin_Parallel;
  v->ops->nvwl2norm      = N_VWL2Norm_Parallel;
  v->ops->nvl1norm       = N_VL1Norm_Parallel;
  v->ops->nvcompare      = N_VCompare_Parallel;
  v->ops->nvinvtest      = N_VInvTest_Parallel;
  v->ops->nvconstrmask   = N_VConstrMask_Parallel;
  v->ops->nvminquotient  = N_VMinQuotient_Parallel;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvdotprodlocal     = N_VDotProdLocal_Parallel;
  v->ops->nvmaxnormlocal     = N_VMaxNormLocal_Parallel;
  v->ops->nvminlocal         = N_VMinLocal_Parallel;
  v->ops->nvl1normlocal      = N_VL1NormLocal_Parallel;
  v->ops->nvinvtestlocal     = N_VInvTestLocal_Parallel;
  v->ops->nvconstrmasklocal  = N_VConstrMaskLocal_Parallel;
  v->ops->nvminquotientlocal = N_VMinQuotientLocal_Parallel;
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Parallel;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Parallel;

  /* single buffer reduction operations */
  v->ops->nvdotprodmultilocal     = N_VDotProdMultiLocal_Parallel;
  v->ops->nvdotprodmultiallreduce = N_VDotProdMultiAllReduce_Parallel;

  /* XBraid interface operations */
  v->ops->nvbufsize   = N_VBufSize_Parallel;
  v->ops->nvbufpack   = N_VBufPack_Parallel;
  v->ops->nvbufunpack = N_VBufUnpack_Parallel;

  /* debugging functions */
  v->ops->nvprint     = N_VPrint_Parallel;
  v->ops->nvprintfile = N_VPrintFile_Parallel;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Parallel)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  v->content = content;

  /* Initialize content */
  content->local_length  = local_length;
  content->global_length = global_length;
  content->comm          = comm;
  content->own_data      = SUNFALSE;
  content->data          = NULL;

  return (v);
}

/* ----------------------------------------------------------------
 * Function to create a new parallel vector
 */

N_Vector N_VNew_Parallel(MPI_Comm comm, sunindextype local_length,
                         sunindextype global_length, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  N_Vector v;
  sunrealtype* data;

  SUNAssertNull(local_length >= 0, SUN_ERR_ARG_OUTOFRANGE);
  SUNAssertNull(global_length >= 0, SUN_ERR_ARG_OUTOFRANGE);

  v = NULL;
  v = N_VNewEmpty_Parallel(comm, local_length, global_length, sunctx);
  SUNCheckLastErrNull();

  /* Create data */
  if (local_length > 0)
  {
    /* Allocate memory */
    data = NULL;
    data = (sunrealtype*)malloc(local_length * sizeof(sunrealtype));
    SUNAssertNull(data, SUN_ERR_MALLOC_FAIL);

    /* Attach data */
    NV_OWN_DATA_P(v) = SUNTRUE;
    NV_DATA_P(v)     = data;
  }

  return (v);
}

/* ----------------------------------------------------------------
 * Function to create a parallel N_Vector with user data component
 */

N_Vector N_VMake_Parallel(MPI_Comm comm, sunindextype local_length,
                          sunindextype global_length, sunrealtype* v_data,
                          SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  N_Vector v;

  SUNAssertNull(local_length >= 0, SUN_ERR_ARG_OUTOFRANGE);
  SUNAssertNull(global_length >= 0, SUN_ERR_ARG_OUTOFRANGE);

  v = NULL;
  v = N_VNewEmpty_Parallel(comm, local_length, global_length, sunctx);
  SUNCheckLastErrNull();

  if (local_length > 0)
  {
    /* Attach data */
    NV_OWN_DATA_P(v) = SUNFALSE;
    NV_DATA_P(v)     = v_data;
  }

  return (v);
}

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_Parallel(SUNDIALS_MAYBE_UNUSED N_Vector v)
{
  return SUNDIALS_NVEC_PARALLEL;
}

/* ----------------------------------------------------------------
 * Function to return global vector length
 */

sunindextype N_VGetLength_Parallel(N_Vector v) { return NV_GLOBLENGTH_P(v); }

/* ----------------------------------------------------------------
 * Function to return local vector length
 */

sunindextype N_VGetLocalLength_Parallel(N_Vector v)
{
  return NV_LOCLENGTH_P(v);
}

/* ----------------------------------------------------------------
 * Function to print the local data in a parallel vector to stdout
 */

void N_VPrint_Parallel(N_Vector x) { N_VPrintFile_Parallel(x, stdout); }

/* ----------------------------------------------------------------
 * Function to print the local data in a parallel vector to outfile
 */

void N_VPrintFile_Parallel(N_Vector x, FILE* outfile)
{
  sunindextype i, N;
  sunrealtype* xd;

  xd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);

  for (i = 0; i < N; i++) { fprintf(outfile, SUN_FORMAT_E "\n", xd[i]); }

  return;
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_Parallel(N_Vector w)
{
  SUNFunctionBegin(w->sunctx);
  N_Vector v;
  N_VectorContent_Parallel content;

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty(w->sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  SUNCheckCallNull(N_VCopyOps(w, v));

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Parallel)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  v->content = content;

  /* Initialize content */
  content->local_length  = NV_LOCLENGTH_P(w);
  content->global_length = NV_GLOBLENGTH_P(w);
  content->comm          = NV_COMM_P(w);
  content->own_data      = SUNFALSE;
  content->data          = NULL;

  return (v);
}

N_Vector N_VClone_Parallel(N_Vector w)
{
  SUNFunctionBegin(w->sunctx);
  N_Vector v;
  sunrealtype* data;
  sunindextype local_length;

  v = NULL;
  v = N_VCloneEmpty_Parallel(w);
  SUNCheckLastErrNull();

  local_length = NV_LOCLENGTH_P(w);

  /* Create data */
  if (local_length > 0)
  {
    /* Allocate memory */
    data = NULL;
    data = (sunrealtype*)malloc(local_length * sizeof(sunrealtype));
    SUNAssertNull(data, SUN_ERR_MALLOC_FAIL);

    /* Attach data */
    NV_OWN_DATA_P(v) = SUNTRUE;
    NV_DATA_P(v)     = data;
  }

  return (v);
}

void N_VDestroy_Parallel(N_Vector v)
{
  if (v == NULL) { return; }

  /* free content */
  if (v->content != NULL)
  {
    if (NV_OWN_DATA_P(v) && NV_DATA_P(v) != NULL)
    {
      free(NV_DATA_P(v));
      NV_DATA_P(v) = NULL;
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

void N_VSpace_Parallel(N_Vector v, sunindextype* lrw, sunindextype* liw)
{
  SUNFunctionBegin(v->sunctx);

  MPI_Comm comm;
  int npes;

  SUNAssertVoid(lrw, SUN_ERR_ARG_CORRUPT);
  SUNAssertVoid(liw, SUN_ERR_ARG_CORRUPT);

  comm = NV_COMM_P(v);
  SUNCheckMPICallVoid(MPI_Comm_size(comm, &npes));

  *lrw = NV_GLOBLENGTH_P(v);
  *liw = 2 * npes;

  return;
}

sunrealtype* N_VGetArrayPointer_Parallel(N_Vector v)
{
  return ((sunrealtype*)NV_DATA_P(v));
}

void N_VSetArrayPointer_Parallel(sunrealtype* v_data, N_Vector v)
{
  if (NV_LOCLENGTH_P(v) > 0) { NV_DATA_P(v) = v_data; }

  return;
}

MPI_Comm N_VGetCommunicator_Parallel(N_Vector v) { return NV_COMM_P(v); }

void N_VLinearSum_Parallel(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                           N_Vector z)
{
  sunindextype i, N;
  sunrealtype c, *xd, *yd, *zd;
  N_Vector v1, v2;
  sunbooleantype test;

  xd = yd = zd = NULL;

  if ((b == ONE) && (z == y))
  { /* BLAS usage: axpy y <- ax+y */
    Vaxpy_Parallel(a, x, y);
    return;
  }

  if ((a == ONE) && (z == x))
  { /* BLAS usage: axpy x <- by+x */
    Vaxpy_Parallel(b, y, x);
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE))
  {
    VSum_Parallel(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE)))
  {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff_Parallel(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE))
  {
    c  = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1_Parallel(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE))
  {
    c  = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2_Parallel(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b)
  {
    VScaleSum_Parallel(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b)
  {
    VScaleDiff_Parallel(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) { zd[i] = (a * xd[i]) + (b * yd[i]); }

  return;
}

void N_VConst_Parallel(sunrealtype c, N_Vector z)
{
  sunindextype i, N;
  sunrealtype* zd;

  zd = NULL;

  N  = NV_LOCLENGTH_P(z);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) { zd[i] = c; }

  return;
}

void N_VProd_Parallel(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i] * yd[i]; }

  return;
}

void N_VDiv_Parallel(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i] / yd[i]; }

  return;
}

void N_VScale_Parallel(sunrealtype c, N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  xd = zd = NULL;

  if (z == x)
  { /* BLAS usage: scale x <- cx */
    VScaleBy_Parallel(c, x);
    return;
  }

  if (c == ONE) { VCopy_Parallel(x, z); }
  else if (c == -ONE) { VNeg_Parallel(x, z); }
  else
  {
    N  = NV_LOCLENGTH_P(x);
    xd = NV_DATA_P(x);
    zd = NV_DATA_P(z);
    for (i = 0; i < N; i++) { zd[i] = c * xd[i]; }
  }

  return;
}

void N_VAbs_Parallel(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) { zd[i] = SUNRabs(xd[i]); }

  return;
}

void N_VInv_Parallel(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) { zd[i] = ONE / xd[i]; }

  return;
}

void N_VAddConst_Parallel(N_Vector x, sunrealtype b, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i] + b; }

  return;
}

sunrealtype N_VDotProdLocal_Parallel(N_Vector x, N_Vector y)
{
  sunindextype i, N;
  sunrealtype sum, *xd, *yd;

  sum = ZERO;
  xd = yd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);

  for (i = 0; i < N; i++) { sum += xd[i] * yd[i]; }
  return (sum);
}

sunrealtype N_VDotProd_Parallel(N_Vector x, N_Vector y)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype lsum, gsum;
  lsum = N_VDotProdLocal_Parallel(x, y);
  SUNCheckLastErrNoRet();
  SUNCheckMPICallNoRet(
    MPI_Allreduce(&lsum, &gsum, 1, MPI_SUNREALTYPE, MPI_SUM, NV_COMM_P(x)));
  return (gsum);
}

sunrealtype N_VMaxNormLocal_Parallel(N_Vector x)
{
  sunindextype i, N;
  sunrealtype max, *xd;

  xd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);

  max = ZERO;

  for (i = 0; i < N; i++)
  {
    if (SUNRabs(xd[i]) > max) { max = SUNRabs(xd[i]); }
  }

  return (max);
}

sunrealtype N_VMaxNorm_Parallel(N_Vector x)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype lmax, gmax;
  lmax = N_VMaxNormLocal_Parallel(x);
  SUNCheckLastErrNoRet();
  SUNCheckMPICallNoRet(
    MPI_Allreduce(&lmax, &gmax, 1, MPI_SUNREALTYPE, MPI_MAX, NV_COMM_P(x)));
  return (gmax);
}

sunrealtype N_VWSqrSumLocal_Parallel(N_Vector x, N_Vector w)
{
  sunindextype i, N;
  sunrealtype sum, prodi, *xd, *wd;

  sum = ZERO;
  xd = wd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  wd = NV_DATA_P(w);

  for (i = 0; i < N; i++)
  {
    prodi = xd[i] * wd[i];
    sum += SUNSQR(prodi);
  }

  return (sum);
}

sunrealtype N_VWrmsNorm_Parallel(N_Vector x, N_Vector w)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype lsum, gsum;
  lsum = N_VWSqrSumLocal_Parallel(x, w);
  SUNCheckLastErrNoRet();
  SUNCheckMPICallNoRet(
    MPI_Allreduce(&lsum, &gsum, 1, MPI_SUNREALTYPE, MPI_SUM, NV_COMM_P(x)));
  return (SUNRsqrt(gsum / (NV_GLOBLENGTH_P(x))));
}

sunrealtype N_VWSqrSumMaskLocal_Parallel(N_Vector x, N_Vector w, N_Vector id)
{
  sunindextype i, N;
  sunrealtype sum, prodi, *xd, *wd, *idd;

  sum = ZERO;
  xd = wd = idd = NULL;

  N   = NV_LOCLENGTH_P(x);
  xd  = NV_DATA_P(x);
  wd  = NV_DATA_P(w);
  idd = NV_DATA_P(id);

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

sunrealtype N_VWrmsNormMask_Parallel(N_Vector x, N_Vector w, N_Vector id)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype lsum, gsum;
  lsum = N_VWSqrSumMaskLocal_Parallel(x, w, id);
  SUNCheckLastErrNoRet();
  SUNCheckMPICallNoRet(
    MPI_Allreduce(&lsum, &gsum, 1, MPI_SUNREALTYPE, MPI_SUM, NV_COMM_P(x)));
  return (SUNRsqrt(gsum / (NV_GLOBLENGTH_P(x))));
}

sunrealtype N_VMinLocal_Parallel(N_Vector x)
{
  sunindextype i, N;
  sunrealtype min, *xd;

  xd  = NULL;
  N   = NV_LOCLENGTH_P(x);
  min = SUN_BIG_REAL;

  if (N > 0)
  {
    xd  = NV_DATA_P(x);
    min = xd[0];
    for (i = 1; i < N; i++)
    {
      if (xd[i] < min) { min = xd[i]; }
    }
  }
  return (min);
}

sunrealtype N_VMin_Parallel(N_Vector x)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype lmin, gmin;
  lmin = N_VMinLocal_Parallel(x);
  SUNCheckLastErrNoRet();
  SUNCheckMPICallNoRet(
    MPI_Allreduce(&lmin, &gmin, 1, MPI_SUNREALTYPE, MPI_MIN, NV_COMM_P(x)));
  return (gmin);
}

sunrealtype N_VWL2Norm_Parallel(N_Vector x, N_Vector w)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype lsum, gsum;
  lsum = N_VWSqrSumLocal_Parallel(x, w);
  SUNCheckLastErrNoRet();
  SUNCheckMPICallNoRet(
    MPI_Allreduce(&lsum, &gsum, 1, MPI_SUNREALTYPE, MPI_SUM, NV_COMM_P(x)));
  return (SUNRsqrt(gsum));
}

sunrealtype N_VL1NormLocal_Parallel(N_Vector x)
{
  sunindextype i, N;
  sunrealtype sum, *xd;

  sum = ZERO;
  xd  = NULL;
  N   = NV_LOCLENGTH_P(x);
  xd  = NV_DATA_P(x);

  for (i = 0; i < N; i++) { sum += SUNRabs(xd[i]); }

  return (sum);
}

sunrealtype N_VL1Norm_Parallel(N_Vector x)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype lsum, gsum;
  lsum = N_VL1NormLocal_Parallel(x);
  SUNCheckLastErrNoRet();
  SUNCheckMPICallNoRet(
    MPI_Allreduce(&lsum, &gsum, 1, MPI_SUNREALTYPE, MPI_SUM, NV_COMM_P(x)));
  return (gsum);
}

void N_VCompare_Parallel(sunrealtype c, N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) { zd[i] = (SUNRabs(xd[i]) >= c) ? ONE : ZERO; }

  return;
}

sunbooleantype N_VInvTestLocal_Parallel(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd, val;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  val = ONE;
  for (i = 0; i < N; i++)
  {
    if (xd[i] == ZERO) { val = ZERO; }
    else { zd[i] = ONE / xd[i]; }
  }

  if (val == ZERO) { return (SUNFALSE); }
  else { return (SUNTRUE); }
}

sunbooleantype N_VInvTest_Parallel(N_Vector x, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);
  sunbooleantype itest = N_VInvTestLocal_Parallel(x, z);
  SUNCheckLastErrNoRet();
  sunrealtype val  = (itest) ? ONE : ZERO;
  sunrealtype gval = ZERO;
  SUNCheckMPICallNoRet(
    MPI_Allreduce(&val, &gval, 1, MPI_SUNREALTYPE, MPI_MIN, NV_COMM_P(x)));
  return (gval != ZERO);
}

sunbooleantype N_VConstrMaskLocal_Parallel(N_Vector c, N_Vector x, N_Vector m)
{
  sunindextype i, N;
  sunrealtype temp;
  sunrealtype *cd, *xd, *md;
  sunbooleantype test;

  cd = xd = md = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  cd = NV_DATA_P(c);
  md = NV_DATA_P(m);

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

sunbooleantype N_VConstrMask_Parallel(N_Vector c, N_Vector x, N_Vector m)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype temp, temp2;
  sunbooleantype test = N_VConstrMaskLocal_Parallel(c, x, m);
  SUNCheckLastErrNoRet();
  temp = test ? ZERO : ONE;
  SUNCheckMPICallNoRet(
    MPI_Allreduce(&temp, &temp2, 1, MPI_SUNREALTYPE, MPI_MAX, NV_COMM_P(x)));
  return (temp2 == ONE) ? SUNFALSE : SUNTRUE;
}

sunrealtype N_VMinQuotientLocal_Parallel(N_Vector num, N_Vector denom)
{
  sunbooleantype notEvenOnce;
  sunindextype i, N;
  sunrealtype *nd, *dd, min;

  nd = dd = NULL;

  N  = NV_LOCLENGTH_P(num);
  nd = NV_DATA_P(num);
  dd = NV_DATA_P(denom);

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

sunrealtype N_VMinQuotient_Parallel(N_Vector num, N_Vector denom)
{
  SUNFunctionBegin(num->sunctx);
  sunrealtype lmin, gmin;
  lmin = N_VMinQuotientLocal_Parallel(num, denom);
  SUNCheckLastErrNoRet();
  SUNCheckMPICallNoRet(
    MPI_Allreduce(&lmin, &gmin, 1, MPI_SUNREALTYPE, MPI_MIN, NV_COMM_P(num)));
  return (gmin);
}

/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */

SUNErrCode N_VLinearCombination_Parallel(int nvec, sunrealtype* c, N_Vector* X,
                                         N_Vector z)
{
  SUNFunctionBegin(X[0]->sunctx);

  int i;
  sunindextype j, N;
  sunrealtype* zd = NULL;
  sunrealtype* xd = NULL;

  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VScale */
  if (nvec == 1)
  {
    N_VScale_Parallel(c[0], X[0], z);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* should have called N_VLinearSum */
  if (nvec == 2)
  {
    N_VLinearSum_Parallel(c[0], X[0], c[1], X[1], z);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N  = NV_LOCLENGTH_P(z);
  zd = NV_DATA_P(z);

  /*
   * X[0] += c[i]*X[i], i = 1,...,nvec-1
   */
  if ((X[0] == z) && (c[0] == ONE))
  {
    for (i = 1; i < nvec; i++)
    {
      xd = NV_DATA_P(X[i]);
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
      xd = NV_DATA_P(X[i]);
      for (j = 0; j < N; j++) { zd[j] += c[i] * xd[j]; }
    }
    return SUN_SUCCESS;
  }

  /*
   * z = sum{ c[i] * X[i] }, i = 0,...,nvec-1
   */
  xd = NV_DATA_P(X[0]);
  for (j = 0; j < N; j++) { zd[j] = c[0] * xd[j]; }
  for (i = 1; i < nvec; i++)
  {
    xd = NV_DATA_P(X[i]);
    for (j = 0; j < N; j++) { zd[j] += c[i] * xd[j]; }
  }
  return SUN_SUCCESS;
}

SUNErrCode N_VScaleAddMulti_Parallel(int nvec, sunrealtype* a, N_Vector x,
                                     N_Vector* Y, N_Vector* Z)
{
  SUNFunctionBegin(x->sunctx);

  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VLinearSum */
  if (nvec == 1)
  {
    N_VLinearSum_Parallel(a[0], x, ONE, Y[0], Z[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z)
  {
    for (i = 0; i < nvec; i++)
    {
      yd = NV_DATA_P(Y[i]);
      for (j = 0; j < N; j++) { yd[j] += a[i] * xd[j]; }
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i = 0; i < nvec; i++)
  {
    yd = NV_DATA_P(Y[i]);
    zd = NV_DATA_P(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = a[i] * xd[j] + yd[j]; }
  }
  return SUN_SUCCESS;
}

SUNErrCode N_VDotProdMulti_Parallel(int nvec, N_Vector x, N_Vector* Y,
                                    sunrealtype* dotprods)
{
  SUNFunctionBegin(x->sunctx);

  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  MPI_Comm comm;

  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VDotProd */
  if (nvec == 1)
  {
    dotprods[0] = N_VDotProd_Parallel(x, Y[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length, data array, and communicator */
  N    = NV_LOCLENGTH_P(x);
  xd   = NV_DATA_P(x);
  comm = NV_COMM_P(x);

  /* compute multiple dot products */
  for (i = 0; i < nvec; i++)
  {
    yd          = NV_DATA_P(Y[i]);
    dotprods[i] = ZERO;
    for (j = 0; j < N; j++) { dotprods[i] += xd[j] * yd[j]; }
  }

  SUNCheckMPICall(MPI_Allreduce(MPI_IN_PLACE, dotprods, nvec, MPI_SUNREALTYPE,
                                MPI_SUM, comm));

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * single buffer reduction operations
 * -----------------------------------------------------------------
 */

SUNErrCode N_VDotProdMultiLocal_Parallel(int nvec, N_Vector x, N_Vector* Y,
                                         sunrealtype* dotprods)
{
  SUNFunctionBegin(x->sunctx);

  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;

  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* get vector length and data array */
  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);

  /* compute multiple dot products */
  for (i = 0; i < nvec; i++)
  {
    yd          = NV_DATA_P(Y[i]);
    dotprods[i] = ZERO;
    for (j = 0; j < N; j++) { dotprods[i] += xd[j] * yd[j]; }
  }

  return SUN_SUCCESS;
}

SUNErrCode N_VDotProdMultiAllReduce_Parallel(int nvec, N_Vector x,
                                             sunrealtype* sum)
{
  SUNFunctionBegin(x->sunctx);

  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* get communicator */
  MPI_Comm comm = NV_COMM_P(x);

  /* perform reduction */
  SUNCheckMPICall(
    MPI_Allreduce(MPI_IN_PLACE, sum, nvec, MPI_SUNREALTYPE, MPI_SUM, comm));

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------
 */

SUNErrCode N_VLinearSumVectorArray_Parallel(int nvec, sunrealtype a,
                                            N_Vector* X, sunrealtype b,
                                            N_Vector* Y, N_Vector* Z)
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

  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VLinearSum */
  if (nvec == 1)
  {
    N_VLinearSum_Parallel(a, X[0], b, Y[0], Z[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* BLAS usage: axpy y <- ax+y */
  if ((b == ONE) && (Z == Y))
  {
    VaxpyVectorArray_Parallel(nvec, a, X, Y);
    return SUN_SUCCESS;
  }

  /* BLAS usage: axpy x <- by+x */
  if ((a == ONE) && (Z == X))
  {
    VaxpyVectorArray_Parallel(nvec, b, Y, X);
    return SUN_SUCCESS;
  }

  /* Case: a == b == 1.0 */
  if ((a == ONE) && (b == ONE))
  {
    VSumVectorArray_Parallel(nvec, X, Y, Z);
    return SUN_SUCCESS;
  }

  /* Cases:                    */
  /*   (1) a == 1.0, b = -1.0, */
  /*   (2) a == -1.0, b == 1.0 */
  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE)))
  {
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    VDiffVectorArray_Parallel(nvec, V2, V1, Z);
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
    VLin1VectorArray_Parallel(nvec, c, V1, V2, Z);
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
    VLin2VectorArray_Parallel(nvec, c, V1, V2, Z);
    return SUN_SUCCESS;
  }

  /* Case: a == b                                                         */
  /* catches case both a and b are 0.0 - user should have called N_VConst */
  if (a == b)
  {
    VScaleSumVectorArray_Parallel(nvec, a, X, Y, Z);
    return SUN_SUCCESS;
  }

  /* Case: a == -b */
  if (a == -b)
  {
    VScaleDiffVectorArray_Parallel(nvec, a, X, Y, Z);
    return SUN_SUCCESS;
  }

  /* Do all cases not handled above:                               */
  /*   (1) a == other, b == 0.0 - user should have called N_VScale */
  /*   (2) a == 0.0, b == other - user should have called N_VScale */
  /*   (3) a,b == other, a !=b, a != -b                            */

  /* get vector length */
  N = NV_LOCLENGTH_P(Z[0]);

  /* compute linear sum for each vector pair in vector arrays */
  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_P(X[i]);
    yd = NV_DATA_P(Y[i]);
    zd = NV_DATA_P(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = a * xd[j] + b * yd[j]; }
  }

  return SUN_SUCCESS;
}

SUNErrCode N_VScaleVectorArray_Parallel(int nvec, sunrealtype* c, N_Vector* X,
                                        N_Vector* Z)
{
  SUNFunctionBegin(X[0]->sunctx);

  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* zd = NULL;

  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VScale */
  if (nvec == 1)
  {
    N_VScale_Parallel(c[0], X[0], Z[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length */
  N = NV_LOCLENGTH_P(Z[0]);

  /*
   * X[i] *= c[i]
   */
  if (X == Z)
  {
    for (i = 0; i < nvec; i++)
    {
      xd = NV_DATA_P(X[i]);
      for (j = 0; j < N; j++) { xd[j] *= c[i]; }
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[i] = c[i] * X[i]
   */
  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_P(X[i]);
    zd = NV_DATA_P(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = c[i] * xd[j]; }
  }
  return SUN_SUCCESS;
}

SUNErrCode N_VConstVectorArray_Parallel(int nvec, sunrealtype c, N_Vector* Z)
{
  SUNFunctionBegin(Z[0]->sunctx);

  int i;
  sunindextype j, N;
  sunrealtype* zd = NULL;

  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VConst */
  if (nvec == 1)
  {
    N_VConst_Parallel(c, Z[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length */
  N = NV_LOCLENGTH_P(Z[0]);

  /* set each vector in the vector array to a constant */
  for (i = 0; i < nvec; i++)
  {
    zd = NV_DATA_P(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = c; }
  }

  return SUN_SUCCESS;
}

SUNErrCode N_VWrmsNormVectorArray_Parallel(int nvec, N_Vector* X, N_Vector* W,
                                           sunrealtype* nrm)
{
  SUNFunctionBegin(X[0]->sunctx);

  sunindextype j, Nl, Ng;
  sunrealtype* wd = NULL;
  sunrealtype* xd = NULL;
  MPI_Comm comm;

  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VWrmsNorm */
  if (nvec == 1)
  {
    nrm[0] = N_VWrmsNorm_Parallel(X[0], W[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector lengths and communicator */
  Nl   = NV_LOCLENGTH_P(X[0]);
  Ng   = NV_GLOBLENGTH_P(X[0]);
  comm = NV_COMM_P(X[0]);

  /* compute the WRMS norm for each vector in the vector array */
  for (int i = 0; i < nvec; i++)
  {
    xd     = NV_DATA_P(X[i]);
    wd     = NV_DATA_P(W[i]);
    nrm[i] = ZERO;
    for (j = 0; j < Nl; j++) { nrm[i] += SUNSQR(xd[j] * wd[j]); }
  }
  SUNCheckMPICall(
    MPI_Allreduce(MPI_IN_PLACE, nrm, nvec, MPI_SUNREALTYPE, MPI_SUM, comm));

  for (int i = 0; i < nvec; i++) { nrm[i] = SUNRsqrt(nrm[i] / Ng); }

  return SUN_SUCCESS;
}

SUNErrCode N_VWrmsNormMaskVectorArray_Parallel(int nvec, N_Vector* X, N_Vector* W,
                                               N_Vector id, sunrealtype* nrm)
{
  SUNFunctionBegin(X[0]->sunctx);

  sunindextype j, Nl, Ng;
  sunrealtype* wd  = NULL;
  sunrealtype* xd  = NULL;
  sunrealtype* idd = NULL;
  MPI_Comm comm;

  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VWrmsNorm */
  if (nvec == 1)
  {
    nrm[0] = N_VWrmsNormMask_Parallel(X[0], W[0], id);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector lengths, communicator, and mask data */
  Nl   = NV_LOCLENGTH_P(X[0]);
  Ng   = NV_GLOBLENGTH_P(X[0]);
  comm = NV_COMM_P(X[0]);
  idd  = NV_DATA_P(id);

  /* compute the WRMS norm for each vector in the vector array */
  for (int i = 0; i < nvec; i++)
  {
    xd     = NV_DATA_P(X[i]);
    wd     = NV_DATA_P(W[i]);
    nrm[i] = ZERO;
    for (j = 0; j < Nl; j++)
    {
      if (idd[j] > ZERO) { nrm[i] += SUNSQR(xd[j] * wd[j]); }
    }
  }
  SUNCheckMPICall(
    MPI_Allreduce(MPI_IN_PLACE, nrm, nvec, MPI_SUNREALTYPE, MPI_SUM, comm));

  for (int i = 0; i < nvec; i++) { nrm[i] = SUNRsqrt(nrm[i] / Ng); }

  return SUN_SUCCESS;
}

SUNErrCode N_VScaleAddMultiVectorArray_Parallel(int nvec, int nsum,
                                                sunrealtype* a, N_Vector* X,
                                                N_Vector** Y, N_Vector** Z)
{
  int i, j;
  sunindextype k, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  N_Vector* YY;
  N_Vector* ZZ;

  SUNFunctionBegin(X[0]->sunctx);

  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);
  SUNAssert(nsum >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1)
  {
    /* should have called N_VLinearSum */
    if (nsum == 1)
    {
      N_VLinearSum_Parallel(a[0], X[0], ONE, Y[0][0], Z[0][0]);
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

    SUNCheckCall(N_VScaleAddMulti_Parallel(nsum, a, X[0], YY, ZZ));

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
    SUNCheckCall(N_VLinearSumVectorArray_Parallel(nvec, a[0], X, ONE, Y[0], Z[0]));
    return SUN_SUCCESS;
  }

  /* ----------------------------
   * Compute multiple linear sums
   * ---------------------------- */

  /* get vector length */
  N = NV_LOCLENGTH_P(X[0]);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z)
  {
    for (i = 0; i < nvec; i++)
    {
      xd = NV_DATA_P(X[i]);
      for (j = 0; j < nsum; j++)
      {
        yd = NV_DATA_P(Y[j][i]);
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
    xd = NV_DATA_P(X[i]);
    for (j = 0; j < nsum; j++)
    {
      yd = NV_DATA_P(Y[j][i]);
      zd = NV_DATA_P(Z[j][i]);
      for (k = 0; k < N; k++) { zd[k] = a[j] * xd[k] + yd[k]; }
    }
  }
  return SUN_SUCCESS;
}

SUNErrCode N_VLinearCombinationVectorArray_Parallel(int nvec, int nsum,
                                                    sunrealtype* c,
                                                    N_Vector** X, N_Vector* Z)
{
  int i;          /* vector arrays index in summation [0,nsum) */
  int j;          /* vector index in vector array     [0,nvec) */
  sunindextype k; /* element index in vector          [0,N)    */
  sunindextype N;
  sunrealtype* zd = NULL;
  sunrealtype* xd = NULL;

  sunrealtype* ctmp;
  N_Vector* Y;

  SUNFunctionBegin(Z[0]->sunctx);

  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);
  SUNAssert(nsum >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1)
  {
    /* should have called N_VScale */
    if (nsum == 1)
    {
      N_VScale_Parallel(c[0], X[0][0], Z[0]);
      SUNCheckLastErr();
      return SUN_SUCCESS;
    }

    /* should have called N_VLinearSum */
    if (nsum == 2)
    {
      N_VLinearSum_Parallel(c[0], X[0][0], c[1], X[1][0], Z[0]);
      SUNCheckLastErr();
      return SUN_SUCCESS;
    }

    /* should have called N_VLinearCombination */
    Y = (N_Vector*)malloc(nsum * sizeof(N_Vector));

    for (i = 0; i < nsum; i++) { Y[i] = X[i][0]; }

    SUNCheckCall(N_VLinearCombination_Parallel(nsum, c, Y, Z[0]));

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

    SUNCheckCall(N_VScaleVectorArray_Parallel(nvec, ctmp, X[0], Z));

    free(ctmp);
    return SUN_SUCCESS;
  }

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 2)
  {
    SUNCheckCall(
      N_VLinearSumVectorArray_Parallel(nvec, c[0], X[0], c[1], X[1], Z));
    return SUN_SUCCESS;
  }

  /* --------------------------
   * Compute linear combination
   * -------------------------- */

  /* get vector length */
  N = NV_LOCLENGTH_P(Z[0]);

  /*
   * X[0][j] += c[i]*X[i][j], i = 1,...,nvec-1
   */
  if ((X[0] == Z) && (c[0] == ONE))
  {
    for (j = 0; j < nvec; j++)
    {
      zd = NV_DATA_P(Z[j]);
      for (i = 1; i < nsum; i++)
      {
        xd = NV_DATA_P(X[i][j]);
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
      zd = NV_DATA_P(Z[j]);
      for (k = 0; k < N; k++) { zd[k] *= c[0]; }
      for (i = 1; i < nsum; i++)
      {
        xd = NV_DATA_P(X[i][j]);
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
    xd = NV_DATA_P(X[0][j]);
    zd = NV_DATA_P(Z[j]);
    for (k = 0; k < N; k++) { zd[k] = c[0] * xd[k]; }
    for (i = 1; i < nsum; i++)
    {
      xd = NV_DATA_P(X[i][j]);
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

SUNErrCode N_VBufSize_Parallel(N_Vector x, sunindextype* size)
{
  *size = NV_LOCLENGTH_P(x) * ((sunindextype)sizeof(sunrealtype));
  return SUN_SUCCESS;
}

SUNErrCode N_VBufPack_Parallel(N_Vector x, void* buf)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i, N;
  sunrealtype* xd = NULL;
  sunrealtype* bd = NULL;

  SUNAssert(buf, SUN_ERR_ARG_CORRUPT);

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  bd = (sunrealtype*)buf;

  for (i = 0; i < N; i++) { bd[i] = xd[i]; }

  return SUN_SUCCESS;
}

SUNErrCode N_VBufUnpack_Parallel(N_Vector x, void* buf)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i, N;
  sunrealtype* xd = NULL;
  sunrealtype* bd = NULL;

  SUNAssert(buf, SUN_ERR_ARG_CORRUPT);

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  bd = (sunrealtype*)buf;

  for (i = 0; i < N; i++) { xd[i] = bd[i]; }

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

static void VCopy_Parallel(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i]; }

  return;
}

static void VSum_Parallel(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i] + yd[i]; }

  return;
}

static void VDiff_Parallel(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i] - yd[i]; }

  return;
}

static void VNeg_Parallel(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) { zd[i] = -xd[i]; }

  return;
}

static void VScaleSum_Parallel(sunrealtype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) { zd[i] = c * (xd[i] + yd[i]); }

  return;
}

static void VScaleDiff_Parallel(sunrealtype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) { zd[i] = c * (xd[i] - yd[i]); }

  return;
}

static void VLin1_Parallel(sunrealtype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) { zd[i] = (a * xd[i]) + yd[i]; }

  return;
}

static void VLin2_Parallel(sunrealtype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) { zd[i] = (a * xd[i]) - yd[i]; }

  return;
}

static void Vaxpy_Parallel(sunrealtype a, N_Vector x, N_Vector y)
{
  sunindextype i, N;
  sunrealtype *xd, *yd;

  xd = yd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);

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

static void VScaleBy_Parallel(sunrealtype a, N_Vector x)
{
  sunindextype i, N;
  sunrealtype* xd;

  xd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);

  for (i = 0; i < N; i++) { xd[i] *= a; }

  return;
}

/*
 * -----------------------------------------------------------------
 * private functions for special cases of vector array operations
 * -----------------------------------------------------------------
 */

static void VSumVectorArray_Parallel(int nvec, N_Vector* X, N_Vector* Y,
                                     N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  N = NV_LOCLENGTH_P(X[0]);

  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_P(X[i]);
    yd = NV_DATA_P(Y[i]);
    zd = NV_DATA_P(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = xd[j] + yd[j]; }
  }
}

static void VDiffVectorArray_Parallel(int nvec, N_Vector* X, N_Vector* Y,
                                      N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  N = NV_LOCLENGTH_P(X[0]);

  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_P(X[i]);
    yd = NV_DATA_P(Y[i]);
    zd = NV_DATA_P(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = xd[j] - yd[j]; }
  }
}

static void VScaleSumVectorArray_Parallel(int nvec, sunrealtype c, N_Vector* X,
                                          N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  N = NV_LOCLENGTH_P(X[0]);

  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_P(X[i]);
    yd = NV_DATA_P(Y[i]);
    zd = NV_DATA_P(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = c * (xd[j] + yd[j]); }
  }
}

static void VScaleDiffVectorArray_Parallel(int nvec, sunrealtype c, N_Vector* X,
                                           N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  N = NV_LOCLENGTH_P(X[0]);

  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_P(X[i]);
    yd = NV_DATA_P(Y[i]);
    zd = NV_DATA_P(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = c * (xd[j] - yd[j]); }
  }
}

static void VLin1VectorArray_Parallel(int nvec, sunrealtype a, N_Vector* X,
                                      N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  N = NV_LOCLENGTH_P(X[0]);

  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_P(X[i]);
    yd = NV_DATA_P(Y[i]);
    zd = NV_DATA_P(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = (a * xd[j]) + yd[j]; }
  }
}

static void VLin2VectorArray_Parallel(int nvec, sunrealtype a, N_Vector* X,
                                      N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  N = NV_LOCLENGTH_P(X[0]);

  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_P(X[i]);
    yd = NV_DATA_P(Y[i]);
    zd = NV_DATA_P(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = (a * xd[j]) - yd[j]; }
  }
}

static void VaxpyVectorArray_Parallel(int nvec, sunrealtype a, N_Vector* X,
                                      N_Vector* Y)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;

  N = NV_LOCLENGTH_P(X[0]);

  if (a == ONE)
  {
    for (i = 0; i < nvec; i++)
    {
      xd = NV_DATA_P(X[i]);
      yd = NV_DATA_P(Y[i]);
      for (j = 0; j < N; j++) { yd[j] += xd[j]; }
    }
    return;
  }

  if (a == -ONE)
  {
    for (i = 0; i < nvec; i++)
    {
      xd = NV_DATA_P(X[i]);
      yd = NV_DATA_P(Y[i]);
      for (j = 0; j < N; j++) { yd[j] -= xd[j]; }
    }

    return;
  }

  for (i = 0; i < nvec; i++)
  {
    xd = NV_DATA_P(X[i]);
    yd = NV_DATA_P(Y[i]);
    for (j = 0; j < N; j++) { yd[j] += a * xd[j]; }
  }
}

/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

SUNErrCode N_VEnableFusedOps_Parallel(N_Vector v, sunbooleantype tf)
{
  SUNFunctionBegin(v->sunctx);

  if (tf)
  {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_Parallel;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_Parallel;
    v->ops->nvdotprodmulti      = N_VDotProdMulti_Parallel;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray     = N_VLinearSumVectorArray_Parallel;
    v->ops->nvscalevectorarray         = N_VScaleVectorArray_Parallel;
    v->ops->nvconstvectorarray         = N_VConstVectorArray_Parallel;
    v->ops->nvwrmsnormvectorarray      = N_VWrmsNormVectorArray_Parallel;
    v->ops->nvwrmsnormmaskvectorarray  = N_VWrmsNormMaskVectorArray_Parallel;
    v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_Parallel;
    v->ops->nvlinearcombinationvectorarray =
      N_VLinearCombinationVectorArray_Parallel;
    /* enable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = N_VDotProdMultiLocal_Parallel;
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

  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearCombination_Parallel(N_Vector v, sunbooleantype tf)
{
  SUNFunctionBegin(v->sunctx);

  /* enable/disable operation */
  if (tf) { v->ops->nvlinearcombination = N_VLinearCombination_Parallel; }
  else { v->ops->nvlinearcombination = NULL; }

  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleAddMulti_Parallel(N_Vector v, sunbooleantype tf)
{
  SUNFunctionBegin(v->sunctx);

  /* enable/disable operation */
  if (tf) { v->ops->nvscaleaddmulti = N_VScaleAddMulti_Parallel; }
  else { v->ops->nvscaleaddmulti = NULL; }

  return SUN_SUCCESS;
}

SUNErrCode N_VEnableDotProdMulti_Parallel(N_Vector v, sunbooleantype tf)
{
  SUNFunctionBegin(v->sunctx);

  /* enable/disable operation */
  if (tf) { v->ops->nvdotprodmulti = N_VDotProdMulti_Parallel; }
  else { v->ops->nvdotprodmulti = NULL; }

  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearSumVectorArray_Parallel(N_Vector v, sunbooleantype tf)
{
  SUNFunctionBegin(v->sunctx);

  /* enable/disable operation */
  if (tf) { v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_Parallel; }
  else { v->ops->nvlinearsumvectorarray = NULL; }

  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleVectorArray_Parallel(N_Vector v, sunbooleantype tf)
{
  SUNFunctionBegin(v->sunctx);

  /* enable/disable operation */
  if (tf) { v->ops->nvscalevectorarray = N_VScaleVectorArray_Parallel; }
  else { v->ops->nvscalevectorarray = NULL; }

  return SUN_SUCCESS;
}

SUNErrCode N_VEnableConstVectorArray_Parallel(N_Vector v, sunbooleantype tf)
{
  SUNFunctionBegin(v->sunctx);

  /* enable/disable operation */
  if (tf) { v->ops->nvconstvectorarray = N_VConstVectorArray_Parallel; }
  else { v->ops->nvconstvectorarray = NULL; }

  return SUN_SUCCESS;
}

SUNErrCode N_VEnableWrmsNormVectorArray_Parallel(N_Vector v, sunbooleantype tf)
{
  SUNFunctionBegin(v->sunctx);

  /* enable/disable operation */
  if (tf) { v->ops->nvwrmsnormvectorarray = N_VWrmsNormVectorArray_Parallel; }
  else { v->ops->nvwrmsnormvectorarray = NULL; }

  return SUN_SUCCESS;
}

SUNErrCode N_VEnableWrmsNormMaskVectorArray_Parallel(N_Vector v, sunbooleantype tf)
{
  SUNFunctionBegin(v->sunctx);

  /* enable/disable operation */
  if (tf)
  {
    v->ops->nvwrmsnormmaskvectorarray = N_VWrmsNormMaskVectorArray_Parallel;
  }
  else { v->ops->nvwrmsnormmaskvectorarray = NULL; }

  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleAddMultiVectorArray_Parallel(N_Vector v,
                                                      sunbooleantype tf)
{
  SUNFunctionBegin(v->sunctx);

  /* enable/disable operation */
  if (tf)
  {
    v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_Parallel;
  }
  else { v->ops->nvscaleaddmultivectorarray = NULL; }

  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearCombinationVectorArray_Parallel(N_Vector v,
                                                          sunbooleantype tf)
{
  SUNFunctionBegin(v->sunctx);

  /* enable/disable operation */
  if (tf)
  {
    v->ops->nvlinearcombinationvectorarray =
      N_VLinearCombinationVectorArray_Parallel;
  }
  else { v->ops->nvlinearcombinationvectorarray = NULL; }

  return SUN_SUCCESS;
}

SUNErrCode N_VEnableDotProdMultiLocal_Parallel(N_Vector v, sunbooleantype tf)
{
  SUNFunctionBegin(v->sunctx);

  /* enable/disable operation */
  if (tf) { v->ops->nvdotprodmultilocal = N_VDotProdMultiLocal_Parallel; }
  else { v->ops->nvdotprodmultilocal = NULL; }

  return SUN_SUCCESS;
}
