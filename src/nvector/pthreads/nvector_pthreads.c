/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * This is the implementation file for a POSIX Threads (Pthreads)
 * implementation of the NVECTOR package using a LOCAL array of
 * structures to pass data to threads.
 * -----------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_pthreads_deprecated.h>
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
static void nvAbs_Pthreads(N_Vector x, N_Vector z);
static void nvAddConst_Pthreads(N_Vector x, sunrealtype b, N_Vector z);
static SUNErrCode nvBufPack_Pthreads(N_Vector x, void* buf);
static SUNErrCode nvBufSize_Pthreads(N_Vector x, sunindextype* size);
static SUNErrCode nvBufUnpack_Pthreads(N_Vector x, void* buf);
static N_Vector nvCloneEmpty_Pthreads(N_Vector w);
static N_Vector nvClone_Pthreads(N_Vector w);
static void nvCompare_Pthreads(sunrealtype c, N_Vector x, N_Vector z);
static SUNErrCode nvConstVectorArray_Pthreads(int nvec, sunrealtype c,
                                              N_Vector* Z);
static void nvConst_Pthreads(sunrealtype c, N_Vector z);
static sunbooleantype nvConstrMask_Pthreads(N_Vector c, N_Vector x, N_Vector m);
static void nvDestroy_Pthreads(N_Vector v);
static void nvDiv_Pthreads(N_Vector x, N_Vector y, N_Vector z);
static SUNErrCode nvDotProdMulti_Pthreads(int nvec, N_Vector x, N_Vector* Y,
                                          sunrealtype* dotprods);
static sunrealtype nvDotProd_Pthreads(N_Vector x, N_Vector y);
static sunrealtype* nvGetArrayPointer_Pthreads(N_Vector v);
static sunindextype nvGetLength_Pthreads(N_Vector v);
static N_Vector_ID nvGetVectorID_Pthreads(N_Vector v);
static sunbooleantype nvInvTest_Pthreads(N_Vector x, N_Vector z);
static void nvInv_Pthreads(N_Vector x, N_Vector z);
static sunrealtype nvL1Norm_Pthreads(N_Vector x);
static SUNErrCode nvLinearCombinationVectorArray_Pthreads(int nvec, int nsum,
                                                          sunrealtype* c,
                                                          N_Vector** X,
                                                          N_Vector* Z);
static SUNErrCode nvLinearCombination_Pthreads(int nvec, sunrealtype* c,
                                               N_Vector* X, N_Vector z);
static SUNErrCode nvLinearSumVectorArray_Pthreads(int nvec, sunrealtype a,
                                                  N_Vector* X, sunrealtype b,
                                                  N_Vector* Y, N_Vector* Z);
static void nvLinearSum_Pthreads(sunrealtype a, N_Vector x, sunrealtype b,
                                 N_Vector y, N_Vector z);
static sunrealtype nvMaxNorm_Pthreads(N_Vector x);
static sunrealtype nvMinQuotient_Pthreads(N_Vector num, N_Vector denom);
static sunrealtype nvMin_Pthreads(N_Vector x);
static void nvPrintFile_Pthreads(N_Vector v, FILE* outfile);
static void nvPrint_Pthreads(N_Vector v);
static void nvProd_Pthreads(N_Vector x, N_Vector y, N_Vector z);
static SUNErrCode nvScaleAddMultiVectorArray_Pthreads(int nvec, int nsum,
                                                      sunrealtype* a,
                                                      N_Vector* X, N_Vector** Y,
                                                      N_Vector** Z);
static SUNErrCode nvScaleAddMulti_Pthreads(int nvec, sunrealtype* a, N_Vector x,
                                           N_Vector* Y, N_Vector* Z);
static SUNErrCode nvScaleVectorArray_Pthreads(int nvec, sunrealtype* c,
                                              N_Vector* X, N_Vector* Z);
static void nvScale_Pthreads(sunrealtype c, N_Vector x, N_Vector z);
static void nvSetArrayPointer_Pthreads(sunrealtype* v_data, N_Vector v);
static sunrealtype nvWL2Norm_Pthreads(N_Vector x, N_Vector w);
static sunrealtype nvWSqrSumLocal_Pthreads(N_Vector x, N_Vector w);
static sunrealtype nvWSqrSumMaskLocal_Pthreads(N_Vector x, N_Vector w,
                                               N_Vector id);
static SUNErrCode nvWrmsNormMaskVectorArray_Pthreads(int nvec, N_Vector* X,
                                                     N_Vector* W, N_Vector id,
                                                     sunrealtype* nrm);
static sunrealtype nvWrmsNormMask_Pthreads(N_Vector x, N_Vector w, N_Vector id);
static SUNErrCode nvWrmsNormVectorArray_Pthreads(int nvec, N_Vector* X,
                                                 N_Vector* W, sunrealtype* nrm);
static sunrealtype nvWrmsNorm_Pthreads(N_Vector x, N_Vector w);

/* Private functions for special cases of vector operations */
static void VCopy_Pthreads(N_Vector x, N_Vector z);             /* z=x       */
static void VSum_Pthreads(N_Vector x, N_Vector y, N_Vector z);  /* z=x+y     */
static void VDiff_Pthreads(N_Vector x, N_Vector y, N_Vector z); /* z=x-y     */
static void VNeg_Pthreads(N_Vector x, N_Vector z);              /* z=-x      */
static void VScaleSum_Pthreads(sunrealtype c, N_Vector x, N_Vector y,
                               N_Vector z); /* z=c(x+y)  */
static void VScaleDiff_Pthreads(sunrealtype c, N_Vector x, N_Vector y,
                                N_Vector z); /* z=c(x-y)  */
static void VLin1_Pthreads(sunrealtype a, N_Vector x, N_Vector y,
                           N_Vector z); /* z=ax+y    */
static void VLin2_Pthreads(sunrealtype a, N_Vector x, N_Vector y,
                           N_Vector z); /* z=ax-y    */
static void Vaxpy_Pthreads(sunrealtype a, N_Vector x, N_Vector y); /* y <- ax+y */
static void VScaleBy_Pthreads(sunrealtype a, N_Vector x); /* x <- ax   */

/* Private functions for special cases of vector array operations */
/* Z=X+Y */
static void VSumVectorArray_Pthreads(int nvec, N_Vector* X, N_Vector* Y,
                                     N_Vector* Z);
/* Z=X-Y */
static void VDiffVectorArray_Pthreads(int nvec, N_Vector* X, N_Vector* Y,
                                      N_Vector* Z);
/* Z=c(X+Y) */
static void VScaleSumVectorArray_Pthreads(int nvec, sunrealtype c, N_Vector* X,
                                          N_Vector* Y, N_Vector* Z);
/* Z=c(X-Y) */
static void VScaleDiffVectorArray_Pthreads(int nvec, sunrealtype c, N_Vector* X,
                                           N_Vector* Y, N_Vector* Z);
/* Z=aX+Y */
static void VLin1VectorArray_Pthreads(int nvec, sunrealtype a, N_Vector* X,
                                      N_Vector* Y, N_Vector* Z);
/* Z=aX-Y */
static void VLin2VectorArray_Pthreads(int nvec, sunrealtype a, N_Vector* X,
                                      N_Vector* Y, N_Vector* Z);
/* Y <- aX+Y */
static void VaxpyVectorArray_Pthreads(int nvec, sunrealtype a, N_Vector* X,
                                      N_Vector* Y);

/* Pthread companion functions for vector operations */
static void* nvConstPt(void* thread_data);
static void* nvLinearSumPt(void* thread_data);
static void* nvProdPt(void* thread_data);
static void* nvDivPt(void* thread_data);
static void* nvScalePt(void* thread_data);
static void* nvAbsPt(void* thread_data);
static void* nvInvPt(void* thread_data);
static void* nvAddConstPt(void* thread_data);
static void* nvComparePt(void* thread_data);
static void* nvDotProdPt(void* thread_data);
static void* nvMaxNormPt(void* thread_data);
static void* nvWSqrSumPt(void* thread_data);
static void* nvMinPt(void* thread_data);
static void* nvWL2NormPt(void* thread_data);
static void* nvL1NormPt(void* thread_data);
static void* nvInvTestPt(void* thread_data);
static void* nvWSqrSumMaskPt(void* thread_data);
static void* nvConstrMaskPt(void* thread_data);
static void* nvMinQuotientPt(void* thread_data);

/* Pthread companion functions special cases of vector operations */
static void* VCopy_PT(void* thread_data);
static void* VSum_PT(void* thread_data);
static void* VDiff_PT(void* thread_data);
static void* VNeg_PT(void* thread_data);
static void* VScaleSum_PT(void* thread_data);
static void* VScaleDiff_PT(void* thread_data);
static void* VLin1_PT(void* thread_data);
static void* VLin2_PT(void* thread_data);
static void* VScaleBy_PT(void* thread_data);
static void* Vaxpy_PT(void* thread_data);

/* Pthread companion functions for fused vector operations */
static void* nvLinearCombinationPt(void* thread_data);
static void* nvScaleAddMultiPt(void* thread_data);
static void* nvDotProdMultiPt(void* thread_data);

/* Pthread companion functions for vector array operations */
static void* nvLinearSumVectorArrayPt(void* thread_data);
static void* nvScaleVectorArrayPt(void* thread_data);
static void* nvConstVectorArrayPt(void* thread_data);
static void* nvWrmsNormVectorArrayPt(void* thread_data);
static void* nvWrmsNormMaskVectorArrayPt(void* thread_data);
static void* nvScaleAddMultiVectorArrayPt(void* thread_data);
static void* nvLinearCombinationVectorArrayPt(void* thread_data);

/* Pthread companion functions special cases of vector array operations */
static void* VSumVectorArray_PT(void* thread_data);
static void* VDiffVectorArray_PT(void* thread_data);
static void* VScaleSumVectorArray_PT(void* thread_data);
static void* VScaleDiffVectorArray_PT(void* thread_data);
static void* VLin1VectorArray_PT(void* thread_data);
static void* VLin2VectorArray_PT(void* thread_data);
static void* VaxpyVectorArray_PT(void* thread_data);

/* Pthread companion functions for XBraid interface operations */
static void* VBufPack_PT(void* thread_data);
static void* VBufUnpack_PT(void* thread_data);

/* Function to determine loop values for threads */
static void nvSplitLoop(int myid, int* nthreads, sunindextype* N,
                        sunindextype* start, sunindextype* end);

/* Function to initialize thread data */
static void nvInitThreadData(Pthreads_Data* thread_data);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID nvGetVectorID_Pthreads(SUNDIALS_MAYBE_UNUSED N_Vector v)
{
  return SUNDIALS_NVEC_PTHREADS;
}

/* ----------------------------------------------------------------------------
 * Function to create a new empty vector
 */

N_Vector N_VNewEmpty_Pthreads(sunindextype length, int num_threads,
                              SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  N_Vector v;
  N_VectorContent_Pthreads content;

  SUNAssertNull(length >= 0, SUN_ERR_ARG_OUTOFRANGE);

  /* Create an empty vector object */
  v = NULL;
  v = N_VNewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = nvGetVectorID_Pthreads;
  v->ops->nvclone           = nvClone_Pthreads;
  v->ops->nvcloneempty      = nvCloneEmpty_Pthreads;
  v->ops->nvdestroy         = nvDestroy_Pthreads;
  v->ops->nvgetarraypointer = nvGetArrayPointer_Pthreads;
  v->ops->nvsetarraypointer = nvSetArrayPointer_Pthreads;
  v->ops->nvgetlength       = nvGetLength_Pthreads;
  v->ops->nvgetlocallength  = nvGetLength_Pthreads;

  /* standard vector operations */
  v->ops->nvlinearsum    = nvLinearSum_Pthreads;
  v->ops->nvconst        = nvConst_Pthreads;
  v->ops->nvprod         = nvProd_Pthreads;
  v->ops->nvdiv          = nvDiv_Pthreads;
  v->ops->nvscale        = nvScale_Pthreads;
  v->ops->nvabs          = nvAbs_Pthreads;
  v->ops->nvinv          = nvInv_Pthreads;
  v->ops->nvaddconst     = nvAddConst_Pthreads;
  v->ops->nvdotprod      = nvDotProd_Pthreads;
  v->ops->nvmaxnorm      = nvMaxNorm_Pthreads;
  v->ops->nvwrmsnormmask = nvWrmsNormMask_Pthreads;
  v->ops->nvwrmsnorm     = nvWrmsNorm_Pthreads;
  v->ops->nvmin          = nvMin_Pthreads;
  v->ops->nvwl2norm      = nvWL2Norm_Pthreads;
  v->ops->nvl1norm       = nvL1Norm_Pthreads;
  v->ops->nvcompare      = nvCompare_Pthreads;
  v->ops->nvinvtest      = nvInvTest_Pthreads;
  v->ops->nvconstrmask   = nvConstrMask_Pthreads;
  v->ops->nvminquotient  = nvMinQuotient_Pthreads;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvdotprodlocal     = nvDotProd_Pthreads;
  v->ops->nvmaxnormlocal     = nvMaxNorm_Pthreads;
  v->ops->nvminlocal         = nvMin_Pthreads;
  v->ops->nvl1normlocal      = nvL1Norm_Pthreads;
  v->ops->nvinvtestlocal     = nvInvTest_Pthreads;
  v->ops->nvconstrmasklocal  = nvConstrMask_Pthreads;
  v->ops->nvminquotientlocal = nvMinQuotient_Pthreads;
  v->ops->nvwsqrsumlocal     = nvWSqrSumLocal_Pthreads;
  v->ops->nvwsqrsummasklocal = nvWSqrSumMaskLocal_Pthreads;

  /* single buffer reduction operations */
  v->ops->nvdotprodmultilocal = nvDotProdMulti_Pthreads;

  /* XBraid interface operations */
  v->ops->nvbufsize   = nvBufSize_Pthreads;
  v->ops->nvbufpack   = nvBufPack_Pthreads;
  v->ops->nvbufunpack = nvBufUnpack_Pthreads;

  /* debugging functions */
  v->ops->nvprint     = nvPrint_Pthreads;
  v->ops->nvprintfile = nvPrintFile_Pthreads;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Pthreads)malloc(sizeof *content);
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

N_Vector N_VNew_Pthreads(sunindextype length, int num_threads, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  N_Vector v;
  sunrealtype* data;

  SUNAssertNull(length >= 0, SUN_ERR_ARG_OUTOFRANGE);

  v = NULL;
  v = N_VNewEmpty_Pthreads(length, num_threads, sunctx);
  SUNCheckLastErrNull();

  /* Create data */
  if (length > 0)
  {
    /* Allocate memory */
    data = NULL;
    data = (sunrealtype*)malloc(length * sizeof(sunrealtype));
    SUNAssertNull(data, SUN_ERR_MALLOC_FAIL);

    /* Attach data */
    NV_OWN_DATA_PT(v) = SUNTRUE;
    NV_DATA_PT(v)     = data;
  }

  return (v);
}

/* ----------------------------------------------------------------------------
 * Function to create a vector with user data component
 */

N_Vector N_VMake_Pthreads(sunindextype length, int num_threads,
                          sunrealtype* v_data, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  N_Vector v;

  SUNAssertNull(length >= 0, SUN_ERR_ARG_OUTOFRANGE);

  v = NULL;
  v = N_VNewEmpty_Pthreads(length, num_threads, sunctx);
  SUNCheckLastErrNull();

  if (length > 0)
  {
    /* Attach data */
    NV_OWN_DATA_PT(v) = SUNFALSE;
    NV_DATA_PT(v)     = v_data;
  }

  return (v);
}

/* ----------------------------------------------------------------------------
 * Function to return number of vector elements
 */
sunindextype nvGetLength_Pthreads(N_Vector v) { return NV_LENGTH_PT(v); }

/* ----------------------------------------------------------------------------
 * Function to print a vector to stdout
 */

void nvPrint_Pthreads(N_Vector x) { nvPrintFile_Pthreads(x, stdout); }

/* ----------------------------------------------------------------------------
 * Function to print a vector to outfile
 */

void nvPrintFile_Pthreads(N_Vector x, FILE* outfile)
{
  sunindextype i, N;
  sunrealtype* xd;

  xd = NULL;

  N  = NV_LENGTH_PT(x);
  xd = NV_DATA_PT(x);

  for (i = 0; i < N; i++) { fprintf(outfile, SUN_FORMAT_E "\n", xd[i]); }
  fprintf(outfile, "\n");

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

N_Vector nvCloneEmpty_Pthreads(N_Vector w)
{
  SUNFunctionBegin(w->sunctx);

  N_Vector v;
  N_VectorContent_Pthreads content;

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty(w->sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  SUNCheckCallNull(N_VCopyOps(w, v));

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Pthreads)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  v->content = content;

  /* Initialize content */
  content->length      = NV_LENGTH_PT(w);
  content->num_threads = NV_NUM_THREADS_PT(w);
  content->own_data    = SUNFALSE;
  content->data        = NULL;

  return (v);
}

/* ----------------------------------------------------------------------------
 * Create new vector from existing vector and attach data
 */

N_Vector nvClone_Pthreads(N_Vector w)
{
  SUNFunctionBegin(w->sunctx);

  N_Vector v;
  sunrealtype* data;
  sunindextype length;

  v = NULL;
  v = nvCloneEmpty_Pthreads(w);
  SUNCheckLastErrNull();

  length = NV_LENGTH_PT(w);

  /* Create data */
  if (length > 0)
  {
    /* Allocate memory */
    data = NULL;
    data = (sunrealtype*)malloc(length * sizeof(sunrealtype));
    SUNAssertNull(data, SUN_ERR_MALLOC_FAIL);

    /* Attach data */
    NV_OWN_DATA_PT(v) = SUNTRUE;
    NV_DATA_PT(v)     = data;
  }

  return (v);
}

/* ----------------------------------------------------------------------------
 * Destroy vector and free vector memory
 */

void nvDestroy_Pthreads(N_Vector v)
{
  if (v == NULL) { return; }

  /* free content */
  if (v->content != NULL)
  {
    if (NV_OWN_DATA_PT(v) && NV_DATA_PT(v) != NULL)
    {
      free(NV_DATA_PT(v));
      NV_DATA_PT(v) = NULL;
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
 * Get storage requirement for vector
 */

/* ----------------------------------------------------------------------------
 * Get vector data pointer
 */

sunrealtype* nvGetArrayPointer_Pthreads(N_Vector v)
{
  return ((sunrealtype*)NV_DATA_PT(v));
}

/* ----------------------------------------------------------------------------
 * Set vector data pointer
 */

void nvSetArrayPointer_Pthreads(sunrealtype* v_data, N_Vector v)
{
  if (NV_LENGTH_PT(v) > 0) { NV_DATA_PT(v) = v_data; }

  return;
}

/* ----------------------------------------------------------------------------
 * Compute linear sum z[i] = a*x[i]+b*y[i]
 */

void nvLinearSum_Pthreads(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                          N_Vector z)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  sunrealtype c;
  N_Vector v1, v2;
  sunbooleantype test;

  if ((b == ONE) && (z == y))
  { /* BLAS usage: axpy y <- ax+y */
    Vaxpy_Pthreads(a, x, y);
    return;
  }

  if ((a == ONE) && (z == x))
  { /* BLAS usage: axpy x <- by+x */
    Vaxpy_Pthreads(b, y, x);
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE))
  {
    VSum_Pthreads(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE)))
  {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff_Pthreads(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE))
  {
    c  = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1_Pthreads(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE))
  {
    c  = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2_Pthreads(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b)
  {
    VScaleSum_Pthreads(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b)
  {
    VScaleDiff_Pthreads(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = (pthread_t*)malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = a;
    thread_data[i].c2 = b;
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].v3 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvLinearSumPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VLinearSum
 */

static void* nvLinearSumPt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype a, b;
  sunrealtype *xd, *yd, *zd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  a  = my_data->c1;
  b  = my_data->c2;
  xd = my_data->v1;
  yd = my_data->v2;
  zd = my_data->v3;

  start = my_data->start;
  end   = my_data->end;

  /* compute linear sum */
  for (i = start; i < end; i++) { zd[i] = (a * xd[i]) + (b * yd[i]); }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Assigns constant value to all vector elements, z[i] = c
 */

void nvConst_Pthreads(sunrealtype c, N_Vector z)
{
  SUNFunctionBegin(z->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(z);
  nthreads = NV_NUM_THREADS_PT(z);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = c;
    thread_data[i].v1 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvConstPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VConst
 */

static void* nvConstPt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype c;
  sunrealtype* zd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  c  = my_data->c1;
  zd = my_data->v1;

  start = my_data->start;
  end   = my_data->end;

  /* assign constant values */
  for (i = start; i < end; i++) { zd[i] = c; }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute componentwise product z[i] = x[i]*y[i]
 */

void nvProd_Pthreads(N_Vector x, N_Vector y, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].v3 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvProdPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and exit */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VProd
 */

static void* nvProdPt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype *xd, *yd, *zd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  xd = my_data->v1;
  yd = my_data->v2;
  zd = my_data->v3;

  start = my_data->start;
  end   = my_data->end;

  /* compute componentwise product */
  for (i = start; i < end; i++) { zd[i] = xd[i] * yd[i]; }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute componentwise division z[i] = x[i]/y[i]
 */

void nvDiv_Pthreads(N_Vector x, N_Vector y, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].v3 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvDivPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VDiv
 */

static void* nvDivPt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype *xd, *yd, *zd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  xd = my_data->v1;
  yd = my_data->v2;
  zd = my_data->v3;

  start = my_data->start;
  end   = my_data->end;

  /* compute componentwise division */
  for (i = start; i < end; i++) { zd[i] = xd[i] / yd[i]; }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute scaler multiplication z[i] = c*x[i]
 */

void nvScale_Pthreads(sunrealtype c, N_Vector x, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  if (z == x)
  { /* BLAS usage: scale x <- cx */
    VScaleBy_Pthreads(c, x);
    return;
  }

  if (c == ONE)
  {
    VCopy_Pthreads(x, z);
    return;
  }
  else if (c == -ONE)
  {
    VNeg_Pthreads(x, z);
    return;
  }
  else
  {
    /* allocate threads and thread data structs */
    N        = NV_LENGTH_PT(x);
    nthreads = NV_NUM_THREADS_PT(x);
    threads  = malloc(nthreads * sizeof(pthread_t));
    SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
    thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
    SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

    /* set thread attributes */
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for (i = 0; i < nthreads; i++)
    {
      /* initialize thread data */
      nvInitThreadData(&thread_data[i]);

      /* compute start and end loop index for thread */
      nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

      /* pack thread data */
      thread_data[i].c1 = c;
      thread_data[i].v1 = NV_DATA_PT(x);
      thread_data[i].v2 = NV_DATA_PT(z);

      /* create threads and call pthread companion function */
      pthread_create(&threads[i], &attr, nvScalePt, (void*)&thread_data[i]);
    }

    /* wait for all threads to finish */
    for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

    /* clean up */
    pthread_attr_destroy(&attr);
    free(threads);
    free(thread_data);
  }

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VScale
 */

static void* nvScalePt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype c;
  sunrealtype *xd, *zd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  c  = my_data->c1;
  xd = my_data->v1;
  zd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* compute scaler multiplication */
  for (i = start; i < end; i++) { zd[i] = c * xd[i]; }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute absolute value of vector components z[i] = SUNRabs(x[i])
 */

void nvAbs_Pthreads(N_Vector x, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvAbsPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VAbs
 */

static void* nvAbsPt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype *xd, *zd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  xd = my_data->v1;
  zd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* compute absolute value of components */
  for (i = start; i < end; i++) { zd[i] = SUNRabs(xd[i]); }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute componentwise inverse z[i] = 1 / x[i]
 */

void nvInv_Pthreads(N_Vector x, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvInvPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VInv
 */

static void* nvInvPt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype *xd, *zd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  xd = my_data->v1;
  zd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* compute componentwise inverse */
  for (i = start; i < end; i++) { zd[i] = ONE / xd[i]; }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute componentwise addition of a scaler to a vector z[i] = x[i] + b
 */

void nvAddConst_Pthreads(N_Vector x, sunrealtype b, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = b;
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvAddConstPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VAddConst
 */

static void* nvAddConstPt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype b;
  sunrealtype *xd, *zd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  b  = my_data->c1;
  xd = my_data->v1;
  zd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* compute componentwise constant addition */
  for (i = start; i < end; i++) { zd[i] = xd[i] + b; }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Computes the dot product of two vectors, a = sum(x[i]*y[i])
 */

sunrealtype nvDotProd_Pthreads(N_Vector x, N_Vector y)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;
  sunrealtype sum = ZERO;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1           = NV_DATA_PT(x);
    thread_data[i].v2           = NV_DATA_PT(y);
    thread_data[i].global_val   = &sum;
    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvDotProdPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return (sum);
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VDotProd
 */

static void* nvDotProdPt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype *xd, *yd;
  sunrealtype local_sum, *global_sum;
  Pthreads_Data* my_data;
  pthread_mutex_t* global_mutex;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  xd = my_data->v1;
  yd = my_data->v2;

  global_sum   = my_data->global_val;
  global_mutex = my_data->global_mutex;

  start = my_data->start;
  end   = my_data->end;

  /* compute dot product */
  local_sum = ZERO;
  for (i = start; i < end; i++) { local_sum += xd[i] * yd[i]; }

  /* update global sum */
  pthread_mutex_lock(global_mutex);
  *global_sum += local_sum;
  pthread_mutex_unlock(global_mutex);

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Computes max norm of the vector
 */

sunrealtype nvMaxNorm_Pthreads(N_Vector x)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;
  sunrealtype max = ZERO;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1           = NV_DATA_PT(x);
    thread_data[i].global_val   = &max;
    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvMaxNormPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return (max);
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VMaxNorm
 */

static void* nvMaxNormPt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype* xd;
  sunrealtype local_max, *global_max;
  Pthreads_Data* my_data;
  pthread_mutex_t* global_mutex;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  xd = my_data->v1;

  global_max   = my_data->global_val;
  global_mutex = my_data->global_mutex;

  start = my_data->start;
  end   = my_data->end;

  /* find local max */
  local_max = ZERO;
  for (i = start; i < end; i++)
  {
    if (SUNRabs(xd[i]) > local_max) { local_max = SUNRabs(xd[i]); }
  }

  /* update global max */
  pthread_mutex_lock(global_mutex);
  if (local_max > *global_max) { *global_max = local_max; }
  pthread_mutex_unlock(global_mutex);

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Computes weighted root mean square norm of a vector
 */

sunrealtype nvWrmsNorm_Pthreads(N_Vector x, N_Vector w)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype sqrsum = nvWSqrSumLocal_Pthreads(x, w);
  SUNCheckLastErrNoRet();
  return (SUNRsqrt(sqrsum / (NV_LENGTH_PT(x))));
}

/* ----------------------------------------------------------------------------
 * Computes weighted square sum of a vector
 */

sunrealtype nvWSqrSumLocal_Pthreads(N_Vector x, N_Vector w)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;
  sunrealtype sum = ZERO;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1           = NV_DATA_PT(x);
    thread_data[i].v2           = NV_DATA_PT(w);
    thread_data[i].global_val   = &sum;
    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvWSqrSumPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return (sum);
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VWrmsNorm
 */

static void* nvWSqrSumPt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype *xd, *wd;
  sunrealtype local_sum, *global_sum;
  Pthreads_Data* my_data;
  pthread_mutex_t* global_mutex;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  xd = my_data->v1;
  wd = my_data->v2;

  global_sum   = my_data->global_val;
  global_mutex = my_data->global_mutex;

  start = my_data->start;
  end   = my_data->end;

  /* compute wrms norm */
  local_sum = ZERO;
  for (i = start; i < end; i++) { local_sum += SUNSQR(xd[i] * wd[i]); }

  /* update global sum */
  pthread_mutex_lock(global_mutex);
  *global_sum += local_sum;
  pthread_mutex_unlock(global_mutex);

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Computes weighted root mean square norm of a masked vector
 */

sunrealtype nvWrmsNormMask_Pthreads(N_Vector x, N_Vector w, N_Vector id)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype sqrsummask = nvWSqrSumMaskLocal_Pthreads(x, w, id);
  SUNCheckLastErrNoRet();
  return (SUNRsqrt(sqrsummask / (NV_LENGTH_PT(x))));
}

/* ----------------------------------------------------------------------------
 * Computes weighted square sum of a masked vector
 */

sunrealtype nvWSqrSumMaskLocal_Pthreads(N_Vector x, N_Vector w, N_Vector id)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;
  sunrealtype sum = ZERO;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1           = NV_DATA_PT(x);
    thread_data[i].v2           = NV_DATA_PT(w);
    thread_data[i].v3           = NV_DATA_PT(id);
    thread_data[i].global_val   = &sum;
    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvWSqrSumMaskPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return (sum);
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VWSqrSumMask
 */

static void* nvWSqrSumMaskPt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype *xd, *wd, *idd;
  sunrealtype local_sum, *global_sum;
  Pthreads_Data* my_data;
  pthread_mutex_t* global_mutex;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  xd  = my_data->v1;
  wd  = my_data->v2;
  idd = my_data->v3;

  global_sum   = my_data->global_val;
  global_mutex = my_data->global_mutex;

  start = my_data->start;
  end   = my_data->end;

  /* compute wrms norm with mask */
  local_sum = ZERO;
  for (i = start; i < end; i++)
  {
    if (idd[i] > ZERO) { local_sum += SUNSQR(xd[i] * wd[i]); }
  }

  /* update global sum */
  pthread_mutex_lock(global_mutex);
  *global_sum += local_sum;
  pthread_mutex_unlock(global_mutex);

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Finds the minimum component of a vector
 */

sunrealtype nvMin_Pthreads(N_Vector x)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;
  sunrealtype min;

  /* initialize global min */
  min = NV_Ith_PT(x, 0);

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1           = NV_DATA_PT(x);
    thread_data[i].global_val   = &min;
    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvMinPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return (min);
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VMin
 */

static void* nvMinPt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype* xd;
  sunrealtype local_min, *global_min;
  Pthreads_Data* my_data;
  pthread_mutex_t* global_mutex;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  xd = my_data->v1;

  global_min   = my_data->global_val;
  global_mutex = my_data->global_mutex;

  start = my_data->start;
  end   = my_data->end;

  /* find local min */
  local_min = *global_min;
  for (i = start; i < end; i++)
  {
    if (xd[i] < local_min) { local_min = xd[i]; }
  }

  /* update global min */
  pthread_mutex_lock(global_mutex);
  if (local_min < *global_min) { *global_min = local_min; }
  pthread_mutex_unlock(global_mutex);

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Computes weighted L2 norm of a vector
 */

sunrealtype nvWL2Norm_Pthreads(N_Vector x, N_Vector w)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;
  sunrealtype sum = ZERO;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1           = NV_DATA_PT(x);
    thread_data[i].v2           = NV_DATA_PT(w);
    thread_data[i].global_val   = &sum;
    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvWL2NormPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return (SUNRsqrt(sum));
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VWL2Norm
 */

static void* nvWL2NormPt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype *xd, *wd;
  sunrealtype local_sum, *global_sum;
  Pthreads_Data* my_data;
  pthread_mutex_t* global_mutex;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  xd = my_data->v1;
  wd = my_data->v2;

  global_sum   = my_data->global_val;
  global_mutex = my_data->global_mutex;

  start = my_data->start;
  end   = my_data->end;

  /* compute WL2 norm */
  local_sum = ZERO;
  for (i = start; i < end; i++) { local_sum += SUNSQR(xd[i] * wd[i]); }

  /* update global sum */
  pthread_mutex_lock(global_mutex);
  *global_sum += local_sum;
  pthread_mutex_unlock(global_mutex);

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Computes L1 norm of a vector
 */

sunrealtype nvL1Norm_Pthreads(N_Vector x)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;
  sunrealtype sum = ZERO;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1           = NV_DATA_PT(x);
    thread_data[i].global_val   = &sum;
    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvL1NormPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return (sum);
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VL1Norm
 */

static void* nvL1NormPt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype* xd;
  sunrealtype local_sum, *global_sum;
  Pthreads_Data* my_data;
  pthread_mutex_t* global_mutex;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  xd = my_data->v1;

  global_sum   = my_data->global_val;
  global_mutex = my_data->global_mutex;

  start = my_data->start;
  end   = my_data->end;

  /* compute L1 norm */
  local_sum = ZERO;
  for (i = start; i < end; i++) { local_sum += SUNRabs(xd[i]); }

  /* update global sum */
  pthread_mutex_lock(global_mutex);
  *global_sum += local_sum;
  pthread_mutex_unlock(global_mutex);

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compare vector component values to a scaler
 */

void nvCompare_Pthreads(sunrealtype c, N_Vector x, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = c;
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvComparePt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VCompare
 */

static void* nvComparePt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype c;
  sunrealtype *xd, *zd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  c  = my_data->c1;
  xd = my_data->v1;
  zd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* compare component to scaler */
  for (i = start; i < end; i++) { zd[i] = (SUNRabs(xd[i]) >= c) ? ONE : ZERO; }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute componentwise inverse z[i] = ONE/x[i] and check if x[i] == ZERO
 */

sunbooleantype nvInvTest_Pthreads(N_Vector x, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  sunrealtype val = ZERO;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1         = NV_DATA_PT(x);
    thread_data[i].v2         = NV_DATA_PT(z);
    thread_data[i].global_val = &val;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvInvTestPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  if (val > ZERO) { return (SUNFALSE); }
  else { return (SUNTRUE); }
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VInvTest
 */

static void* nvInvTestPt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype *xd, *zd;
  sunrealtype local_val, *global_val;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  xd = my_data->v1;
  zd = my_data->v2;

  global_val = my_data->global_val;

  start = my_data->start;
  end   = my_data->end;

  /* compute inverse with check for divide by ZERO */
  local_val = ZERO;
  for (i = start; i < end; i++)
  {
    if (xd[i] == ZERO) { local_val = ONE; }
    else { zd[i] = ONE / xd[i]; }
  }

  /* update global val */
  if (local_val > ZERO) { *global_val = local_val; }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute constraint mask of a vector
 */

sunbooleantype nvConstrMask_Pthreads(N_Vector c, N_Vector x, N_Vector m)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  sunrealtype val = ZERO;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1         = NV_DATA_PT(c);
    thread_data[i].v2         = NV_DATA_PT(x);
    thread_data[i].v3         = NV_DATA_PT(m);
    thread_data[i].global_val = &val;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvConstrMaskPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  if (val > ZERO) { return (SUNFALSE); }
  else { return (SUNTRUE); }
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VConstrMask
 */

static void* nvConstrMaskPt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype *cd, *xd, *md;
  sunrealtype local_val, *global_val;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  cd = my_data->v1;
  xd = my_data->v2;
  md = my_data->v3;

  global_val = my_data->global_val;

  start = my_data->start;
  end   = my_data->end;

  /* compute constraint mask */
  local_val = ZERO;
  for (i = start; i < end; i++)
  {
    md[i] = ZERO;

    /* Continue if no constraints were set for the variable */
    if (cd[i] == ZERO) { continue; }

    /* Check if a set constraint has been violated */
    if ((SUNRabs(cd[i]) > ONEPT5 && xd[i] * cd[i] <= ZERO) ||
        (SUNRabs(cd[i]) > HALF && xd[i] * cd[i] < ZERO))
    {
      local_val = md[i] = ONE;
    }
  }

  /* update global val */
  if (local_val > ZERO) { *global_val = local_val; }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute minimum componentwise quotient
 */

sunrealtype nvMinQuotient_Pthreads(N_Vector num, N_Vector denom)
{
  SUNFunctionBegin(num->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;
  sunrealtype min = SUN_BIG_REAL;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(num);
  nthreads = NV_NUM_THREADS_PT(num);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1           = NV_DATA_PT(num);
    thread_data[i].v2           = NV_DATA_PT(denom);
    thread_data[i].global_val   = &min;
    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvMinQuotientPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return (min);
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VConstrMask
 */

static void* nvMinQuotientPt(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype *nd, *dd;
  sunrealtype local_min, *global_min;
  Pthreads_Data* my_data;
  pthread_mutex_t* global_mutex;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  nd = my_data->v1;
  dd = my_data->v2;

  global_min   = my_data->global_val;
  global_mutex = my_data->global_mutex;

  start = my_data->start;
  end   = my_data->end;

  /* compute minimum quotient */
  local_min = SUN_BIG_REAL;
  for (i = start; i < end; i++)
  {
    if (dd[i] == ZERO) { continue; }
    local_min = SUNMIN(local_min, nd[i] / dd[i]);
  }

  /* update global min */
  pthread_mutex_lock(global_mutex);
  if (local_min < *global_min) { *global_min = local_min; }
  pthread_mutex_unlock(global_mutex);

  /* exit */
  pthread_exit(NULL);
}

/*
 * -----------------------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------
 * Compute the linear combination z = c[i]*X[i]
 */

SUNErrCode nvLinearCombination_Pthreads(int nvec, sunrealtype* c, N_Vector* X,
                                        N_Vector z)
{
  SUNFunctionBegin(X[0]->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VScale */
  if (nvec == 1)
  {
    nvScale_Pthreads(c[0], X[0], z);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* should have called N_VLinearSum */
  if (nvec == 2)
  {
    nvLinearSum_Pthreads(c[0], X[0], c[1], X[1], z);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N        = NV_LENGTH_PT(z);
  nthreads = NV_NUM_THREADS_PT(z);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec  = nvec;
    thread_data[i].cvals = c;
    thread_data[i].Y1    = X;
    thread_data[i].x1    = z;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvLinearCombinationPt,
                   (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VLinearCombination
 */

static void* nvLinearCombinationPt(void* thread_data)
{
  Pthreads_Data* my_data;
  sunindextype j, start, end;

  int i;
  sunrealtype* c  = NULL;
  sunrealtype* xd = NULL;
  sunrealtype* zd = NULL;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  start = my_data->start;
  end   = my_data->end;

  c  = my_data->cvals;
  zd = NV_DATA_PT(my_data->x1);

  /*
   * X[0] += c[i]*X[i], i = 1,...,nvec-1
   */
  if ((my_data->Y1[0] == my_data->x1) && (c[0] == ONE))
  {
    for (i = 1; i < my_data->nvec; i++)
    {
      xd = NV_DATA_PT(my_data->Y1[i]);
      for (j = start; j < end; j++) { zd[j] += c[i] * xd[j]; }
    }
    pthread_exit(NULL);
  }

  /*
   * X[0] = c[0] * X[0] + sum{ c[i] * X[i] }, i = 1,...,nvec-1
   */
  if (my_data->Y1[0] == my_data->x1)
  {
    for (j = start; j < end; j++) { zd[j] *= c[0]; }
    for (i = 1; i < my_data->nvec; i++)
    {
      xd = NV_DATA_PT(my_data->Y1[i]);
      for (j = start; j < end; j++) { zd[j] += c[i] * xd[j]; }
    }
    pthread_exit(NULL);
  }

  /*
   * z = sum{ c[i] * X[i] }, i = 0,...,nvec-1
   */
  xd = NV_DATA_PT(my_data->Y1[0]);
  for (j = start; j < end; j++) { zd[j] = c[0] * xd[j]; }
  for (i = 1; i < my_data->nvec; i++)
  {
    xd = NV_DATA_PT(my_data->Y1[i]);
    for (j = start; j < end; j++) { zd[j] += c[i] * xd[j]; }
  }
  pthread_exit(NULL);
}

/* -----------------------------------------------------------------------------
 * Compute multiple linear sums Z[i] = Y[i] + a*x
 */

SUNErrCode nvScaleAddMulti_Pthreads(int nvec, sunrealtype* a, N_Vector x,
                                    N_Vector* Y, N_Vector* Z)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VLinearSum */
  if (nvec == 1)
  {
    nvLinearSum_Pthreads(a[0], x, ONE, Y[0], Z[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec  = nvec;
    thread_data[i].cvals = a;
    thread_data[i].x1    = x;
    thread_data[i].Y1    = Y;
    thread_data[i].Y2    = Z;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvScaleAddMultiPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VScaleAddMulti
 */

static void* nvScaleAddMultiPt(void* thread_data)
{
  Pthreads_Data* my_data;
  sunindextype j, start, end;

  int i;
  sunrealtype* a  = NULL;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  start = my_data->start;
  end   = my_data->end;

  a  = my_data->cvals;
  xd = NV_DATA_PT(my_data->x1);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (my_data->Y1 == my_data->Y2)
  {
    for (i = 0; i < my_data->nvec; i++)
    {
      yd = NV_DATA_PT(my_data->Y1[i]);
      for (j = start; j < end; j++) { yd[j] += a[i] * xd[j]; }
    }
    pthread_exit(NULL);
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i = 0; i < my_data->nvec; i++)
  {
    yd = NV_DATA_PT(my_data->Y1[i]);
    zd = NV_DATA_PT(my_data->Y2[i]);
    for (j = start; j < end; j++) { zd[j] = a[i] * xd[j] + yd[j]; }
  }
  pthread_exit(NULL);
}

/* -----------------------------------------------------------------------------
 * Compute the dot product of a vector with multiple vectors, a[i] = sum(x*Y[i])
 */

SUNErrCode nvDotProdMulti_Pthreads(int nvec, N_Vector x, N_Vector* Y,
                                   sunrealtype* dotprods)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VDotProd */
  if (nvec == 1)
  {
    dotprods[0] = nvDotProd_Pthreads(x, Y[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* initialize output array */
  for (i = 0; i < nvec; i++) { dotprods[i] = ZERO; }

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec  = nvec;
    thread_data[i].x1    = x;
    thread_data[i].Y1    = Y;
    thread_data[i].cvals = dotprods;

    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvDotProdMultiPt, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VDotProdMulti
 */

static void* nvDotProdMultiPt(void* thread_data)
{
  Pthreads_Data* my_data;
  sunindextype j, start, end;
  pthread_mutex_t* lock;

  int i;
  sunrealtype sum;
  sunrealtype* dotprods = NULL;
  sunrealtype* xd       = NULL;
  sunrealtype* yd       = NULL;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  start = my_data->start;
  end   = my_data->end;
  lock  = my_data->global_mutex;

  xd       = NV_DATA_PT(my_data->x1);
  dotprods = my_data->cvals;

  /* compute multiple dot products */
  for (i = 0; i < my_data->nvec; i++)
  {
    yd  = NV_DATA_PT(my_data->Y1[i]);
    sum = ZERO;
    for (j = start; j < end; j++) { sum += xd[j] * yd[j]; }
    /* update global sum */
    pthread_mutex_lock(lock);
    dotprods[i] += sum;
    pthread_mutex_unlock(lock);
  }

  /* exit */
  pthread_exit(NULL);
}

/*
 * -----------------------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------
 * Compute multiple linear sums Z[i] = a*X[i] + b*Y[i]
 */

SUNErrCode nvLinearSumVectorArray_Pthreads(int nvec, sunrealtype a, N_Vector* X,
                                           sunrealtype b, N_Vector* Y,
                                           N_Vector* Z)
{
  SUNFunctionBegin(X[0]->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  sunrealtype c;
  N_Vector* V1;
  N_Vector* V2;
  sunbooleantype test;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VLinearSum */
  if (nvec == 1)
  {
    nvLinearSum_Pthreads(a, X[0], b, Y[0], Z[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* BLAS usage: axpy y <- ax+y */
  if ((b == ONE) && (Z == Y))
  {
    VaxpyVectorArray_Pthreads(nvec, a, X, Y);
    return SUN_SUCCESS;
  }

  /* BLAS usage: axpy x <- by+x */
  if ((a == ONE) && (Z == X))
  {
    VaxpyVectorArray_Pthreads(nvec, b, Y, X);
    return SUN_SUCCESS;
  }

  /* Case: a == b == 1.0 */
  if ((a == ONE) && (b == ONE))
  {
    VSumVectorArray_Pthreads(nvec, X, Y, Z);
    return SUN_SUCCESS;
  }

  /* Cases:                    */
  /*   (1) a == 1.0, b = -1.0, */
  /*   (2) a == -1.0, b == 1.0 */
  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE)))
  {
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    VDiffVectorArray_Pthreads(nvec, V2, V1, Z);
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
    VLin1VectorArray_Pthreads(nvec, c, V1, V2, Z);
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
    VLin2VectorArray_Pthreads(nvec, c, V1, V2, Z);
    return SUN_SUCCESS;
  }

  /* Case: a == b                                                         */
  /* catches case both a and b are 0.0 - user should have called N_VConst */
  if (a == b)
  {
    VScaleSumVectorArray_Pthreads(nvec, a, X, Y, Z);
    return SUN_SUCCESS;
  }

  /* Case: a == -b */
  if (a == -b)
  {
    VScaleDiffVectorArray_Pthreads(nvec, a, X, Y, Z);
    return SUN_SUCCESS;
  }

  /* Do all cases not handled above:                               */
  /*   (1) a == other, b == 0.0 - user should have called N_VScale */
  /*   (2) a == 0.0, b == other - user should have called N_VScale */
  /*   (3) a,b == other, a !=b, a != -b                            */

  /* get vector length and data array */
  N        = NV_LENGTH_PT(Z[0]);
  nthreads = NV_NUM_THREADS_PT(Z[0]);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec = nvec;
    thread_data[i].c1   = a;
    thread_data[i].c2   = b;
    thread_data[i].Y1   = X;
    thread_data[i].Y2   = Y;
    thread_data[i].Y3   = Z;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvLinearSumVectorArrayPt,
                   (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VLinearSumVectorArray
 */

static void* nvLinearSumVectorArrayPt(void* thread_data)
{
  Pthreads_Data* my_data;
  sunindextype j, start, end;

  int i;
  sunrealtype a, b;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  start = my_data->start;
  end   = my_data->end;

  a = my_data->c1;
  b = my_data->c2;

  /* compute linear sum for each vector pair in vector arrays */
  for (i = 0; i < my_data->nvec; i++)
  {
    xd = NV_DATA_PT(my_data->Y1[i]);
    yd = NV_DATA_PT(my_data->Y2[i]);
    zd = NV_DATA_PT(my_data->Y3[i]);
    for (j = start; j < end; j++) { zd[j] = a * xd[j] + b * yd[j]; }
  }

  /* exit */
  pthread_exit(NULL);
}

/* -----------------------------------------------------------------------------
 * Scale multiple vectors Z[i] = c[i]*X[i]
 */

SUNErrCode nvScaleVectorArray_Pthreads(int nvec, sunrealtype* c, N_Vector* X,
                                       N_Vector* Z)
{
  SUNFunctionBegin(X[0]->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VScale */
  if (nvec == 1)
  {
    nvScale_Pthreads(c[0], X[0], Z[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N        = NV_LENGTH_PT(Z[0]);
  nthreads = NV_NUM_THREADS_PT(Z[0]);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec  = nvec;
    thread_data[i].cvals = c;
    thread_data[i].Y1    = X;
    thread_data[i].Y2    = Z;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvScaleVectorArrayPt,
                   (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VScaleVectorArray
 */

static void* nvScaleVectorArrayPt(void* thread_data)
{
  Pthreads_Data* my_data;
  sunindextype j, start, end;

  int i;
  sunrealtype* c;
  sunrealtype* xd = NULL;
  sunrealtype* zd = NULL;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  start = my_data->start;
  end   = my_data->end;

  c = my_data->cvals;

  /*
   * X[i] *= c[i]
   */
  if (my_data->Y1 == my_data->Y2)
  {
    for (i = 0; i < my_data->nvec; i++)
    {
      xd = NV_DATA_PT(my_data->Y1[i]);
      for (j = start; j < end; j++) { xd[j] *= c[i]; }
    }
    pthread_exit(NULL);
  }

  /*
   * Z[i] = c[i] * X[i]
   */
  for (i = 0; i < my_data->nvec; i++)
  {
    xd = NV_DATA_PT(my_data->Y1[i]);
    zd = NV_DATA_PT(my_data->Y2[i]);
    for (j = start; j < end; j++) { zd[j] = c[i] * xd[j]; }
  }
  pthread_exit(NULL);
}

/* -----------------------------------------------------------------------------
 * Set multiple vectors to a constant value Z[i] = c
 */

SUNErrCode nvConstVectorArray_Pthreads(int nvec, sunrealtype c, N_Vector* Z)
{
  SUNFunctionBegin(Z[0]->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VConst */
  if (nvec == 1)
  {
    nvConst_Pthreads(c, Z[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N        = NV_LENGTH_PT(Z[0]);
  nthreads = NV_NUM_THREADS_PT(Z[0]);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec = nvec;
    thread_data[i].c1   = c;
    thread_data[i].Y1   = Z;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvConstVectorArrayPt,
                   (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VConstVectorArray
 */

static void* nvConstVectorArrayPt(void* thread_data)
{
  Pthreads_Data* my_data;
  sunindextype j, start, end;

  int i;
  sunrealtype* zd = NULL;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  start = my_data->start;
  end   = my_data->end;

  /* set each vector in the vector array to a constant */
  for (i = 0; i < my_data->nvec; i++)
  {
    zd = NV_DATA_PT(my_data->Y1[i]);
    for (j = start; j < end; j++) { zd[j] = my_data->c1; }
  }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute the weighted root mean square norm of multiple vectors
 */

SUNErrCode nvWrmsNormVectorArray_Pthreads(int nvec, N_Vector* X, N_Vector* W,
                                          sunrealtype* nrm)
{
  SUNFunctionBegin(X[0]->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VWrmsNorm */
  if (nvec == 1)
  {
    nrm[0] = nvWrmsNorm_Pthreads(X[0], W[0]);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* initialize output array */
  for (i = 0; i < nvec; i++) { nrm[i] = ZERO; }

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(X[0]);
  nthreads = NV_NUM_THREADS_PT(X[0]);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec  = nvec;
    thread_data[i].Y1    = X;
    thread_data[i].Y2    = W;
    thread_data[i].cvals = nrm;

    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvWrmsNormVectorArrayPt,
                   (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* finalize wrms calculation */
  for (i = 0; i < nvec; i++) { nrm[i] = SUNRsqrt(nrm[i] / N); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VWrmsNormVectorArray
 */

static void* nvWrmsNormVectorArrayPt(void* thread_data)
{
  Pthreads_Data* my_data;
  sunindextype j, start, end;
  pthread_mutex_t* lock;

  int i;
  sunrealtype sum;
  sunrealtype* nrm = NULL;
  sunrealtype* xd  = NULL;
  sunrealtype* wd  = NULL;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  start = my_data->start;
  end   = my_data->end;
  lock  = my_data->global_mutex;

  nrm = my_data->cvals;

  /* compute the WRMS norm for each vector in the vector array */
  for (i = 0; i < my_data->nvec; i++)
  {
    xd  = NV_DATA_PT(my_data->Y1[i]);
    wd  = NV_DATA_PT(my_data->Y2[i]);
    sum = ZERO;
    for (j = start; j < end; j++) { sum += SUNSQR(xd[j] * wd[j]); }
    /* update global sum */
    pthread_mutex_lock(lock);
    nrm[i] += sum;
    pthread_mutex_unlock(lock);
  }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute the weighted root mean square norm of multiple vectors
 */

SUNErrCode nvWrmsNormMaskVectorArray_Pthreads(int nvec, N_Vector* X, N_Vector* W,
                                              N_Vector id, sunrealtype* nrm)
{
  SUNFunctionBegin(X[0]->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;

  /* invalid number of vectors */
  SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

  /* should have called N_VWrmsNorm */
  if (nvec == 1)
  {
    nrm[0] = nvWrmsNormMask_Pthreads(X[0], W[0], id);
    SUNCheckLastErr();
    return SUN_SUCCESS;
  }

  /* initialize output array */
  for (i = 0; i < nvec; i++) { nrm[i] = ZERO; }

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(X[0]);
  nthreads = NV_NUM_THREADS_PT(X[0]);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec  = nvec;
    thread_data[i].Y1    = X;
    thread_data[i].Y2    = W;
    thread_data[i].x1    = id;
    thread_data[i].cvals = nrm;

    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvWrmsNormMaskVectorArrayPt,
                   (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* finalize wrms calculation */
  for (i = 0; i < nvec; i++) { nrm[i] = SUNRsqrt(nrm[i] / N); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VWrmsNormVectorArray
 */

static void* nvWrmsNormMaskVectorArrayPt(void* thread_data)
{
  Pthreads_Data* my_data;
  sunindextype j, start, end;
  pthread_mutex_t* lock;

  int i;
  sunrealtype sum;
  sunrealtype* nrm = NULL;
  sunrealtype* xd  = NULL;
  sunrealtype* wd  = NULL;
  sunrealtype* idd = NULL;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  start = my_data->start;
  end   = my_data->end;
  lock  = my_data->global_mutex;

  nrm = my_data->cvals;
  idd = NV_DATA_PT(my_data->x1);

  /* compute the WRMS norm for each vector in the vector array */
  for (i = 0; i < my_data->nvec; i++)
  {
    xd  = NV_DATA_PT(my_data->Y1[i]);
    wd  = NV_DATA_PT(my_data->Y2[i]);
    sum = ZERO;
    for (j = start; j < end; j++)
    {
      if (idd[j] > ZERO) { sum += SUNSQR(xd[j] * wd[j]); }
    }
    /* update global sum */
    pthread_mutex_lock(lock);
    nrm[i] += sum;
    pthread_mutex_unlock(lock);
  }

  /* exit */
  pthread_exit(NULL);
}

/* -----------------------------------------------------------------------------
 * Scale and add a vector to multiple vectors Z = Y + a*X
 */

SUNErrCode nvScaleAddMultiVectorArray_Pthreads(int nvec, int nsum,
                                               sunrealtype* a, N_Vector* X,
                                               N_Vector** Y, N_Vector** Z)
{
  SUNFunctionBegin(X[0]->sunctx);

  sunindextype N;
  int i, j, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

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
      nvLinearSum_Pthreads(a[0], X[0], ONE, Y[0][0], Z[0][0]);
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

    SUNCheckCall(nvScaleAddMulti_Pthreads(nsum, a, X[0], YY, ZZ));

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
    SUNCheckCall(nvLinearSumVectorArray_Pthreads(nvec, a[0], X, ONE, Y[0], Z[0]));
    return SUN_SUCCESS;
  }

  /* ----------------------------
   * Compute multiple linear sums
   * ---------------------------- */

  /* get vector length and data array */
  N        = NV_LENGTH_PT(X[0]);
  nthreads = NV_NUM_THREADS_PT(X[0]);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec  = nvec;
    thread_data[i].nsum  = nsum;
    thread_data[i].cvals = a;
    thread_data[i].Y1    = X;
    thread_data[i].ZZ1   = Y;
    thread_data[i].ZZ2   = Z;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvScaleAddMultiVectorArrayPt,
                   (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VScaleAddMultiVectorArray
 */

static void* nvScaleAddMultiVectorArrayPt(void* thread_data)
{
  Pthreads_Data* my_data;
  sunindextype k, start, end;

  int i, j;
  sunrealtype* a  = NULL;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  start = my_data->start;
  end   = my_data->end;

  a = my_data->cvals;

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (my_data->ZZ1 == my_data->ZZ2)
  {
    for (i = 0; i < my_data->nvec; i++)
    {
      xd = NV_DATA_PT(my_data->Y1[i]);
      for (j = 0; j < my_data->nsum; j++)
      {
        yd = NV_DATA_PT(my_data->ZZ1[j][i]);
        for (k = start; k < end; k++) { yd[k] += a[j] * xd[k]; }
      }
    }
    pthread_exit(NULL);
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i = 0; i < my_data->nvec; i++)
  {
    xd = NV_DATA_PT(my_data->Y1[i]);
    for (j = 0; j < my_data->nsum; j++)
    {
      yd = NV_DATA_PT(my_data->ZZ1[j][i]);
      zd = NV_DATA_PT(my_data->ZZ2[j][i]);
      for (k = start; k < end; k++) { zd[k] = a[j] * xd[k] + yd[k]; }
    }
  }
  pthread_exit(NULL);
}

/* -----------------------------------------------------------------------------
 * Compute a linear combination for multiple vectors
 */

SUNErrCode nvLinearCombinationVectorArray_Pthreads(int nvec, int nsum,
                                                   sunrealtype* c, N_Vector** X,
                                                   N_Vector* Z)
{
  SUNFunctionBegin(X[0][0]->sunctx);

  sunindextype N;
  int i, j, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

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
      nvScale_Pthreads(c[0], X[0][0], Z[0]);
      SUNCheckLastErr();
      return SUN_SUCCESS;
    }

    /* should have called N_VLinearSum */
    if (nsum == 2)
    {
      nvLinearSum_Pthreads(c[0], X[0][0], c[1], X[1][0], Z[0]);
      SUNCheckLastErr();
      return SUN_SUCCESS;
    }

    /* should have called N_VLinearCombination */
    Y = (N_Vector*)malloc(nsum * sizeof(N_Vector));
    SUNAssert(Y, SUN_ERR_MALLOC_FAIL);

    for (i = 0; i < nsum; i++) { Y[i] = X[i][0]; }

    SUNCheckCall(nvLinearCombination_Pthreads(nsum, c, Y, Z[0]));

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

    SUNCheckCall(nvScaleVectorArray_Pthreads(nvec, ctmp, X[0], Z));

    free(ctmp);
    return SUN_SUCCESS;
  }

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 2)
  {
    SUNCheckCall(nvLinearSumVectorArray_Pthreads(nvec, c[0], X[0], c[1], X[1], Z));
    return SUN_SUCCESS;
  }

  /* --------------------------
   * Compute linear combination
   * -------------------------- */

  /* get vector length and data array */
  N        = NV_LENGTH_PT(Z[0]);
  nthreads = NV_NUM_THREADS_PT(Z[0]);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec  = nvec;
    thread_data[i].nsum  = nsum;
    thread_data[i].cvals = c;
    thread_data[i].ZZ1   = X;
    thread_data[i].Y1    = Z;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, nvLinearCombinationVectorArrayPt,
                   (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VLinearCombinationVectorArray
 */

static void* nvLinearCombinationVectorArrayPt(void* thread_data)
{
  Pthreads_Data* my_data;
  sunindextype k, start, end;

  int i; /* vector arrays index in summation [0,nsum) */
  int j; /* vector index in vector array     [0,nvec) */
  sunrealtype* c  = NULL;
  sunrealtype* zd = NULL;
  sunrealtype* xd = NULL;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  start = my_data->start;
  end   = my_data->end;

  c = my_data->cvals;

  /*
   * X[0][j] += c[i]*X[i][j], i = 1,...,nvec-1
   */
  if ((my_data->ZZ1[0] == my_data->Y1) && (c[0] == ONE))
  {
    for (j = 0; j < my_data->nvec; j++)
    {
      zd = NV_DATA_PT(my_data->Y1[j]);
      for (i = 1; i < my_data->nsum; i++)
      {
        xd = NV_DATA_PT(my_data->ZZ1[i][j]);
        for (k = start; k < end; k++) { zd[k] += c[i] * xd[k]; }
      }
    }
    pthread_exit(NULL);
  }

  /*
   * X[0][j] = c[0] * X[0][j] + sum{ c[i] * X[i][j] }, i = 1,...,nvec-1
   */
  if (my_data->ZZ1[0] == my_data->Y1)
  {
    for (j = 0; j < my_data->nvec; j++)
    {
      zd = NV_DATA_PT(my_data->Y1[j]);
      for (k = start; k < end; k++) { zd[k] *= c[0]; }
      for (i = 1; i < my_data->nsum; i++)
      {
        xd = NV_DATA_PT(my_data->ZZ1[i][j]);
        for (k = start; k < end; k++) { zd[k] += c[i] * xd[k]; }
      }
    }
    pthread_exit(NULL);
  }

  /*
   * Z[j] = sum{ c[i] * X[i][j] }, i = 0,...,nvec-1
   */
  for (j = 0; j < my_data->nvec; j++)
  {
    xd = NV_DATA_PT(my_data->ZZ1[0][j]);
    zd = NV_DATA_PT(my_data->Y1[j]);
    for (k = start; k < end; k++) { zd[k] = c[0] * xd[k]; }
    for (i = 1; i < my_data->nsum; i++)
    {
      xd = NV_DATA_PT(my_data->ZZ1[i][j]);
      for (k = start; k < end; k++) { zd[k] += c[i] * xd[k]; }
    }
  }
  pthread_exit(NULL);
}

/*
 * -----------------------------------------------------------------
 * OPTIONAL XBraid interface operations
 * -----------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------
 * Set buffer size
 */

SUNErrCode nvBufSize_Pthreads(N_Vector x, sunindextype* size)
{
  *size = NV_LENGTH_PT(x) * ((sunindextype)sizeof(sunrealtype));
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Pack butter
 */

SUNErrCode nvBufPack_Pthreads(N_Vector x, void* buf)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  SUNAssert(buf, SUN_ERR_ARG_CORRUPT);

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = (sunrealtype*)buf;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VBufPack_PT, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VBufPack
 */

static void* VBufPack_PT(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype *xd, *bd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  xd = my_data->v1;
  bd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* pack the buffer */
  for (i = start; i < end; i++) { bd[i] = xd[i]; }

  /* exit */
  pthread_exit(NULL);
}

/* -----------------------------------------------------------------------------
 * Unpack butter
 */

SUNErrCode nvBufUnpack_Pthreads(N_Vector x, void* buf)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  SUNAssert(buf, SUN_ERR_ARG_CORRUPT);

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssert(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssert(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = (sunrealtype*)buf;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VBufUnpack_PT, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VBufUnpack
 */

static void* VBufUnpack_PT(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype *xd, *bd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  xd = my_data->v1;
  bd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* unpack the buffer */
  for (i = start; i < end; i++) { xd[i] = bd[i]; }

  /* exit */
  pthread_exit(NULL);
}

/*
 * -----------------------------------------------------------------
 * private functions for special cases of vector operations
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Copy vector components into second vector
 */

static void VCopy_Pthreads(N_Vector x, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VCopy_PT, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to VCopy
 */

static void* VCopy_PT(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype *xd, *zd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  xd = my_data->v1;
  zd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* copy vector components */
  for (i = start; i < end; i++) { zd[i] = xd[i]; }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute vector sum
 */

static void VSum_Pthreads(N_Vector x, N_Vector y, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].v3 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VSum_PT, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to VSum
 */

static void* VSum_PT(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype *xd, *yd, *zd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  xd = my_data->v1;
  yd = my_data->v2;
  zd = my_data->v3;

  start = my_data->start;
  end   = my_data->end;

  /* compute vector sum */
  for (i = start; i < end; i++) { zd[i] = xd[i] + yd[i]; }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute vector difference
 */

static void VDiff_Pthreads(N_Vector x, N_Vector y, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].v3 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VDiff_PT, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to VDiff
 */

static void* VDiff_PT(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype *xd, *yd, *zd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  xd = my_data->v1;
  yd = my_data->v2;
  zd = my_data->v3;

  start = my_data->start;
  end   = my_data->end;

  /* compute vector difference */
  for (i = start; i < end; i++) { zd[i] = xd[i] - yd[i]; }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute the negative of a vector
 */

static void VNeg_Pthreads(N_Vector x, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VNeg_PT, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to VNeg
 */

static void* VNeg_PT(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype *xd, *zd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  xd = my_data->v1;
  zd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* compute negative of vector */
  for (i = start; i < end; i++) { zd[i] = -xd[i]; }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute scaled vector sum
 */

static void VScaleSum_Pthreads(sunrealtype c, N_Vector x, N_Vector y, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = c;
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].v3 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VScaleSum_PT, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to VScaleSum
 */

static void* VScaleSum_PT(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype c;
  sunrealtype *xd, *yd, *zd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  c  = my_data->c1;
  xd = my_data->v1;
  yd = my_data->v2;
  zd = my_data->v3;

  start = my_data->start;
  end   = my_data->end;

  /* compute scaled vector sum */
  for (i = start; i < end; i++) { zd[i] = c * (xd[i] + yd[i]); }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute scaled vector difference
 */

static void VScaleDiff_Pthreads(sunrealtype c, N_Vector x, N_Vector y, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = c;
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].v3 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VScaleDiff_PT, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to VScaleDiff
 */

static void* VScaleDiff_PT(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype c;
  sunrealtype *xd, *yd, *zd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  c  = my_data->c1;
  xd = my_data->v1;
  yd = my_data->v2;
  zd = my_data->v3;

  start = my_data->start;
  end   = my_data->end;

  /* compute scaled vector difference */
  for (i = start; i < end; i++) { zd[i] = c * (xd[i] - yd[i]); }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute vector sum z[i] = a*x[i]+y[i]
 */

static void VLin1_Pthreads(sunrealtype a, N_Vector x, N_Vector y, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = a;
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].v3 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VLin1_PT, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to VLin1
 */

static void* VLin1_PT(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype a;
  sunrealtype *xd, *yd, *zd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  a  = my_data->c1;
  xd = my_data->v1;
  yd = my_data->v2;
  zd = my_data->v3;

  start = my_data->start;
  end   = my_data->end;

  /* compute vector sum */
  for (i = start; i < end; i++) { zd[i] = (a * xd[i]) + yd[i]; }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute vector difference z[i] = a*x[i]-y[i]
 */

static void VLin2_Pthreads(sunrealtype a, N_Vector x, N_Vector y, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = a;
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].v3 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VLin2_PT, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to VLin2
 */

static void* VLin2_PT(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype a;
  sunrealtype *xd, *yd, *zd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  a  = my_data->c1;
  xd = my_data->v1;
  yd = my_data->v2;
  zd = my_data->v3;

  start = my_data->start;
  end   = my_data->end;

  /* compute vector difference */
  for (i = start; i < end; i++) { zd[i] = (a * xd[i]) - yd[i]; }

  /* exit */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute special cases of linear sum
 */

static void Vaxpy_Pthreads(sunrealtype a, N_Vector x, N_Vector y)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = a;
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, Vaxpy_PT, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to Vaxpy
 */

static void* Vaxpy_PT(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype a;
  sunrealtype *xd, *yd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  a  = my_data->c1;
  xd = my_data->v1;
  yd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* compute axpy */
  if (a == ONE)
  {
    for (i = start; i < end; i++) { yd[i] += xd[i]; }

    /* exit */
    pthread_exit(NULL);
  }

  if (a == -ONE)
  {
    for (i = start; i < end; i++) { yd[i] -= xd[i]; }

    /* exit */
    pthread_exit(NULL);
  }

  for (i = start; i < end; i++) { yd[i] += a * xd[i]; }

  /* return */
  pthread_exit(NULL);
}

/* ----------------------------------------------------------------------------
 * Compute scaled vector
 */

static void VScaleBy_Pthreads(sunrealtype a, N_Vector x)
{
  SUNFunctionBegin(x->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(x);
  nthreads = NV_NUM_THREADS_PT(x);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i = 0; i < nthreads; i++)
  {
    /* initialize thread data */
    nvInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = a;
    thread_data[i].v1 = NV_DATA_PT(x);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VScaleBy_PT, (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to VScaleBy
 */

static void* VScaleBy_PT(void* thread_data)
{
  sunindextype i, start, end;
  sunrealtype a;
  sunrealtype* xd;
  Pthreads_Data* my_data;

  /* extract thread data */
  my_data = (Pthreads_Data*)thread_data;

  a  = my_data->c1;
  xd = my_data->v1;

  start = my_data->start;
  end   = my_data->end;

  /* compute scaled vector */
  for (i = start; i < end; i++) { xd[i] *= a; }

  /* exit */
  pthread_exit(NULL);
}

/*
 * -----------------------------------------------------------------------------
 * private functions for special cases of vector array operations
 * -----------------------------------------------------------------------------
 */

static void VSumVectorArray_Pthreads(int nvec, N_Vector* X, N_Vector* Y,
                                     N_Vector* Z)
{
  SUNFunctionBegin(X[0]->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(X[0]);
  nthreads = NV_NUM_THREADS_PT(X[0]);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* pack thread data, distribute loop indices, and create threads/call kernel */
  for (i = 0; i < nthreads; i++)
  {
    nvInitThreadData(&thread_data[i]);

    thread_data[i].nvec = nvec;
    thread_data[i].Y1   = X;
    thread_data[i].Y2   = Y;
    thread_data[i].Y3   = Z;

    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    pthread_create(&threads[i], &attr, VSumVectorArray_PT,
                   (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);
}

static void* VSumVectorArray_PT(void* thread_data)
{
  Pthreads_Data* my_data;
  sunindextype j, start, end;

  int i;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  my_data = (Pthreads_Data*)thread_data;
  start   = my_data->start;
  end     = my_data->end;

  for (i = 0; i < my_data->nvec; i++)
  {
    xd = NV_DATA_PT(my_data->Y1[i]);
    yd = NV_DATA_PT(my_data->Y2[i]);
    zd = NV_DATA_PT(my_data->Y3[i]);
    for (j = start; j < end; j++) { zd[j] = xd[j] + yd[j]; }
  }

  pthread_exit(NULL);
}

static void VDiffVectorArray_Pthreads(int nvec, N_Vector* X, N_Vector* Y,
                                      N_Vector* Z)
{
  SUNFunctionBegin(X[0]->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(X[0]);
  nthreads = NV_NUM_THREADS_PT(X[0]);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* pack thread data, distribute loop indices, and create threads/call kernel */
  for (i = 0; i < nthreads; i++)
  {
    nvInitThreadData(&thread_data[i]);

    thread_data[i].nvec = nvec;
    thread_data[i].Y1   = X;
    thread_data[i].Y2   = Y;
    thread_data[i].Y3   = Z;

    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    pthread_create(&threads[i], &attr, VDiffVectorArray_PT,
                   (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);
}

static void* VDiffVectorArray_PT(void* thread_data)
{
  Pthreads_Data* my_data;
  sunindextype j, start, end;

  int i;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  my_data = (Pthreads_Data*)thread_data;
  start   = my_data->start;
  end     = my_data->end;

  for (i = 0; i < my_data->nvec; i++)
  {
    xd = NV_DATA_PT(my_data->Y1[i]);
    yd = NV_DATA_PT(my_data->Y2[i]);
    zd = NV_DATA_PT(my_data->Y3[i]);
    for (j = start; j < end; j++) { zd[j] = xd[j] - yd[j]; }
  }

  pthread_exit(NULL);
}

static void VScaleSumVectorArray_Pthreads(int nvec, sunrealtype c, N_Vector* X,
                                          N_Vector* Y, N_Vector* Z)
{
  SUNFunctionBegin(X[0]->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(X[0]);
  nthreads = NV_NUM_THREADS_PT(X[0]);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* pack thread data, distribute loop indices, and create threads/call kernel */
  for (i = 0; i < nthreads; i++)
  {
    nvInitThreadData(&thread_data[i]);

    thread_data[i].nvec = nvec;
    thread_data[i].c1   = c;
    thread_data[i].Y1   = X;
    thread_data[i].Y2   = Y;
    thread_data[i].Y3   = Z;

    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    pthread_create(&threads[i], &attr, VScaleSumVectorArray_PT,
                   (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);
}

static void* VScaleSumVectorArray_PT(void* thread_data)
{
  Pthreads_Data* my_data;
  sunindextype j, start, end;

  int i;
  sunrealtype c;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  my_data = (Pthreads_Data*)thread_data;
  start   = my_data->start;
  end     = my_data->end;
  c       = my_data->c1;

  for (i = 0; i < my_data->nvec; i++)
  {
    xd = NV_DATA_PT(my_data->Y1[i]);
    yd = NV_DATA_PT(my_data->Y2[i]);
    zd = NV_DATA_PT(my_data->Y3[i]);
    for (j = start; j < end; j++) { zd[j] = c * (xd[j] + yd[j]); }
  }

  pthread_exit(NULL);
}

static void VScaleDiffVectorArray_Pthreads(int nvec, sunrealtype c, N_Vector* X,
                                           N_Vector* Y, N_Vector* Z)
{
  SUNFunctionBegin(X[0]->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(X[0]);
  nthreads = NV_NUM_THREADS_PT(X[0]);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* pack thread data, distribute loop indices, and create threads/call kernel */
  for (i = 0; i < nthreads; i++)
  {
    nvInitThreadData(&thread_data[i]);

    thread_data[i].nvec = nvec;
    thread_data[i].c1   = c;
    thread_data[i].Y1   = X;
    thread_data[i].Y2   = Y;
    thread_data[i].Y3   = Z;

    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    pthread_create(&threads[i], &attr, VScaleDiffVectorArray_PT,
                   (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);
}

static void* VScaleDiffVectorArray_PT(void* thread_data)
{
  Pthreads_Data* my_data;
  sunindextype j, start, end;

  int i;
  sunrealtype c;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  my_data = (Pthreads_Data*)thread_data;
  start   = my_data->start;
  end     = my_data->end;
  c       = my_data->c1;

  for (i = 0; i < my_data->nvec; i++)
  {
    xd = NV_DATA_PT(my_data->Y1[i]);
    yd = NV_DATA_PT(my_data->Y2[i]);
    zd = NV_DATA_PT(my_data->Y3[i]);
    for (j = start; j < end; j++) { zd[j] = c * (xd[j] - yd[j]); }
  }

  pthread_exit(NULL);
}

static void VLin1VectorArray_Pthreads(int nvec, sunrealtype a, N_Vector* X,
                                      N_Vector* Y, N_Vector* Z)
{
  SUNFunctionBegin(X[0]->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(X[0]);
  nthreads = NV_NUM_THREADS_PT(X[0]);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* pack thread data, distribute loop indices, and create threads/call kernel */
  for (i = 0; i < nthreads; i++)
  {
    nvInitThreadData(&thread_data[i]);

    thread_data[i].nvec = nvec;
    thread_data[i].c1   = a;
    thread_data[i].Y1   = X;
    thread_data[i].Y2   = Y;
    thread_data[i].Y3   = Z;

    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    pthread_create(&threads[i], &attr, VLin1VectorArray_PT,
                   (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);
}

static void* VLin1VectorArray_PT(void* thread_data)
{
  Pthreads_Data* my_data;
  sunindextype j, start, end;

  int i;
  sunrealtype a;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  my_data = (Pthreads_Data*)thread_data;
  start   = my_data->start;
  end     = my_data->end;
  a       = my_data->c1;

  for (i = 0; i < my_data->nvec; i++)
  {
    xd = NV_DATA_PT(my_data->Y1[i]);
    yd = NV_DATA_PT(my_data->Y2[i]);
    zd = NV_DATA_PT(my_data->Y3[i]);
    for (j = start; j < end; j++) { zd[j] = (a * xd[j]) + yd[j]; }
  }

  pthread_exit(NULL);
}

static void VLin2VectorArray_Pthreads(int nvec, sunrealtype a, N_Vector* X,
                                      N_Vector* Y, N_Vector* Z)
{
  SUNFunctionBegin(X[0]->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(X[0]);
  nthreads = NV_NUM_THREADS_PT(X[0]);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* pack thread data, distribute loop indices, and create threads/call kernel */
  for (i = 0; i < nthreads; i++)
  {
    nvInitThreadData(&thread_data[i]);

    thread_data[i].nvec = nvec;
    thread_data[i].c1   = a;
    thread_data[i].Y1   = X;
    thread_data[i].Y2   = Y;
    thread_data[i].Y3   = Z;

    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    pthread_create(&threads[i], &attr, VLin2VectorArray_PT,
                   (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);
}

static void* VLin2VectorArray_PT(void* thread_data)
{
  Pthreads_Data* my_data;
  sunindextype j, start, end;

  int i;
  sunrealtype a;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  my_data = (Pthreads_Data*)thread_data;
  start   = my_data->start;
  end     = my_data->end;
  a       = my_data->c1;

  for (i = 0; i < my_data->nvec; i++)
  {
    xd = NV_DATA_PT(my_data->Y1[i]);
    yd = NV_DATA_PT(my_data->Y2[i]);
    zd = NV_DATA_PT(my_data->Y3[i]);
    for (j = start; j < end; j++) { zd[j] = (a * xd[j]) - yd[j]; }
  }

  pthread_exit(NULL);
}

static void VaxpyVectorArray_Pthreads(int nvec, sunrealtype a, N_Vector* X,
                                      N_Vector* Y)
{
  SUNFunctionBegin(X[0]->sunctx);

  sunindextype N;
  int i, nthreads;
  pthread_t* threads;
  Pthreads_Data* thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N        = NV_LENGTH_PT(X[0]);
  nthreads = NV_NUM_THREADS_PT(X[0]);
  threads  = malloc(nthreads * sizeof(pthread_t));
  SUNAssertVoid(threads, SUN_ERR_MALLOC_FAIL);
  thread_data = (Pthreads_Data*)malloc(nthreads * sizeof(struct Pthreads_Data_));
  SUNAssertVoid(thread_data, SUN_ERR_MALLOC_FAIL);

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* pack thread data, distribute loop indices, and create threads/call kernel */
  for (i = 0; i < nthreads; i++)
  {
    nvInitThreadData(&thread_data[i]);

    thread_data[i].nvec = nvec;
    thread_data[i].c1   = a;
    thread_data[i].Y1   = X;
    thread_data[i].Y2   = Y;

    nvSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    pthread_create(&threads[i], &attr, VaxpyVectorArray_PT,
                   (void*)&thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i = 0; i < nthreads; i++) { pthread_join(threads[i], NULL); }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);
}

static void* VaxpyVectorArray_PT(void* thread_data)
{
  Pthreads_Data* my_data;
  sunindextype j, start, end;

  int i;
  sunrealtype a;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;

  my_data = (Pthreads_Data*)thread_data;
  start   = my_data->start;
  end     = my_data->end;
  a       = my_data->c1;

  if (a == ONE)
  {
    for (i = 0; i < my_data->nvec; i++)
    {
      xd = NV_DATA_PT(my_data->Y1[i]);
      yd = NV_DATA_PT(my_data->Y2[i]);
      for (j = start; j < end; j++) { yd[j] += xd[j]; }
    }
    pthread_exit(NULL);
  }

  if (a == -ONE)
  {
    for (i = 0; i < my_data->nvec; i++)
    {
      xd = NV_DATA_PT(my_data->Y1[i]);
      yd = NV_DATA_PT(my_data->Y2[i]);
      for (j = start; j < end; j++) { yd[j] -= xd[j]; }
    }
    pthread_exit(NULL);
  }

  for (i = 0; i < my_data->nvec; i++)
  {
    xd = NV_DATA_PT(my_data->Y1[i]);
    yd = NV_DATA_PT(my_data->Y2[i]);
    for (j = start; j < end; j++) { yd[j] += a * xd[j]; }
  }
  pthread_exit(NULL);
}

/*
 * -----------------------------------------------------------------
 * private utility functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Determine loop indices for a thread
 */

static void nvSplitLoop(int myid, int* nthreads, sunindextype* N,
                        sunindextype* start, sunindextype* end)
{
  sunindextype q, r; /* quotient and remainder */

  /* work per thread and leftover work */
  q = *N / *nthreads;
  r = *N % *nthreads;

  /* assign work */
  if (myid < r)
  {
    *start = myid * q + myid;
    *end   = *start + q + 1;
  }
  else
  {
    *start = myid * q + r;
    *end   = *start + q;
  }
}

/* ----------------------------------------------------------------------------
 * Initialize values of local thread data struct
 */

static void nvInitThreadData(Pthreads_Data* thread_data)
{
  thread_data->start = -1;
  thread_data->end   = -1;

#ifdef NAN
  thread_data->c1 = NAN;
  thread_data->c2 = NAN;
#else
  thread_data->c1 = ZERO;
  thread_data->c2 = ZERO;
#endif

  thread_data->v1           = NULL;
  thread_data->v2           = NULL;
  thread_data->v3           = NULL;
  thread_data->global_val   = NULL;
  thread_data->global_mutex = NULL;

  thread_data->nvec  = ZERO;
  thread_data->nsum  = ZERO;
  thread_data->cvals = NULL;
  thread_data->Y1    = NULL;
  thread_data->Y2    = NULL;
  thread_data->Y3    = NULL;
}

/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

SUNErrCode N_VEnableFusedOps_Pthreads(N_Vector v, sunbooleantype tf)
{
  if (tf)
  {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = nvLinearCombination_Pthreads;
    v->ops->nvscaleaddmulti     = nvScaleAddMulti_Pthreads;
    v->ops->nvdotprodmulti      = nvDotProdMulti_Pthreads;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray     = nvLinearSumVectorArray_Pthreads;
    v->ops->nvscalevectorarray         = nvScaleVectorArray_Pthreads;
    v->ops->nvconstvectorarray         = nvConstVectorArray_Pthreads;
    v->ops->nvwrmsnormvectorarray      = nvWrmsNormVectorArray_Pthreads;
    v->ops->nvwrmsnormmaskvectorarray  = nvWrmsNormMaskVectorArray_Pthreads;
    v->ops->nvscaleaddmultivectorarray = nvScaleAddMultiVectorArray_Pthreads;
    v->ops->nvlinearcombinationvectorarray =
      nvLinearCombinationVectorArray_Pthreads;
    /* enable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = nvDotProdMulti_Pthreads;
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
    /* enable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = NULL;
  }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearCombination_Pthreads(N_Vector v, sunbooleantype tf)
{
  v->ops->nvlinearcombination = tf ? nvLinearCombination_Pthreads : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleAddMulti_Pthreads(N_Vector v, sunbooleantype tf)
{
  v->ops->nvscaleaddmulti = tf ? nvScaleAddMulti_Pthreads : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableDotProdMulti_Pthreads(N_Vector v, sunbooleantype tf)
{
  v->ops->nvdotprodmulti      = tf ? nvDotProdMulti_Pthreads : NULL;
  v->ops->nvdotprodmultilocal = tf ? nvDotProdMulti_Pthreads : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearSumVectorArray_Pthreads(N_Vector v, sunbooleantype tf)
{
  v->ops->nvlinearsumvectorarray = tf ? nvLinearSumVectorArray_Pthreads : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleVectorArray_Pthreads(N_Vector v, sunbooleantype tf)
{
  v->ops->nvscalevectorarray = tf ? nvScaleVectorArray_Pthreads : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableConstVectorArray_Pthreads(N_Vector v, sunbooleantype tf)
{
  v->ops->nvconstvectorarray = tf ? nvConstVectorArray_Pthreads : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableWrmsNormVectorArray_Pthreads(N_Vector v, sunbooleantype tf)
{
  v->ops->nvwrmsnormvectorarray = tf ? nvWrmsNormVectorArray_Pthreads : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableWrmsNormMaskVectorArray_Pthreads(N_Vector v, sunbooleantype tf)
{
  v->ops->nvwrmsnormmaskvectorarray = tf ? nvWrmsNormMaskVectorArray_Pthreads
                                         : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleAddMultiVectorArray_Pthreads(N_Vector v,
                                                      sunbooleantype tf)
{
  v->ops->nvscaleaddmultivectorarray = tf ? nvScaleAddMultiVectorArray_Pthreads
                                          : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearCombinationVectorArray_Pthreads(N_Vector v,
                                                          sunbooleantype tf)
{
  v->ops->nvlinearcombinationvectorarray =
    tf ? nvLinearCombinationVectorArray_Pthreads : NULL;
  return SUN_SUCCESS;
}

/* Deprecated concrete operation wrappers */

void N_VAbs_Pthreads(N_Vector x, N_Vector z) { nvAbs_Pthreads(x, z); }

void N_VAddConst_Pthreads(N_Vector x, sunrealtype b, N_Vector z)
{
  nvAddConst_Pthreads(x, b, z);
}

SUNErrCode N_VBufPack_Pthreads(N_Vector x, void* buf)
{
  return nvBufPack_Pthreads(x, buf);
}

SUNErrCode N_VBufSize_Pthreads(N_Vector x, sunindextype* size)
{
  return nvBufSize_Pthreads(x, size);
}

SUNErrCode N_VBufUnpack_Pthreads(N_Vector x, void* buf)
{
  return nvBufUnpack_Pthreads(x, buf);
}

N_Vector N_VCloneEmpty_Pthreads(N_Vector w) { return nvCloneEmpty_Pthreads(w); }

N_Vector N_VClone_Pthreads(N_Vector w) { return nvClone_Pthreads(w); }

void N_VCompare_Pthreads(sunrealtype c, N_Vector x, N_Vector z)
{
  nvCompare_Pthreads(c, x, z);
}

SUNErrCode N_VConstVectorArray_Pthreads(int nvec, sunrealtype c, N_Vector* Z)
{
  return nvConstVectorArray_Pthreads(nvec, c, Z);
}

void N_VConst_Pthreads(sunrealtype c, N_Vector z) { nvConst_Pthreads(c, z); }

sunbooleantype N_VConstrMask_Pthreads(N_Vector c, N_Vector x, N_Vector m)
{
  return nvConstrMask_Pthreads(c, x, m);
}

void N_VDestroy_Pthreads(N_Vector v) { nvDestroy_Pthreads(v); }

void N_VDiv_Pthreads(N_Vector x, N_Vector y, N_Vector z)
{
  nvDiv_Pthreads(x, y, z);
}

SUNErrCode N_VDotProdMulti_Pthreads(int nvec, N_Vector x, N_Vector* Y,
                                    sunrealtype* dotprods)
{
  return nvDotProdMulti_Pthreads(nvec, x, Y, dotprods);
}

sunrealtype N_VDotProd_Pthreads(N_Vector x, N_Vector y)
{
  return nvDotProd_Pthreads(x, y);
}

sunrealtype* N_VGetArrayPointer_Pthreads(N_Vector v)
{
  return nvGetArrayPointer_Pthreads(v);
}

sunindextype N_VGetLength_Pthreads(N_Vector v)
{
  return nvGetLength_Pthreads(v);
}

N_Vector_ID N_VGetVectorID_Pthreads(N_Vector v)
{
  return nvGetVectorID_Pthreads(v);
}

sunbooleantype N_VInvTest_Pthreads(N_Vector x, N_Vector z)
{
  return nvInvTest_Pthreads(x, z);
}

void N_VInv_Pthreads(N_Vector x, N_Vector z) { nvInv_Pthreads(x, z); }

sunrealtype N_VL1Norm_Pthreads(N_Vector x) { return nvL1Norm_Pthreads(x); }

SUNErrCode N_VLinearCombinationVectorArray_Pthreads(int nvec, int nsum,
                                                    sunrealtype* c,
                                                    N_Vector** X, N_Vector* Z)
{
  return nvLinearCombinationVectorArray_Pthreads(nvec, nsum, c, X, Z);
}

SUNErrCode N_VLinearCombination_Pthreads(int nvec, sunrealtype* c, N_Vector* X,
                                         N_Vector z)
{
  return nvLinearCombination_Pthreads(nvec, c, X, z);
}

SUNErrCode N_VLinearSumVectorArray_Pthreads(int nvec, sunrealtype a,
                                            N_Vector* X, sunrealtype b,
                                            N_Vector* Y, N_Vector* Z)
{
  return nvLinearSumVectorArray_Pthreads(nvec, a, X, b, Y, Z);
}

void N_VLinearSum_Pthreads(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                           N_Vector z)
{
  nvLinearSum_Pthreads(a, x, b, y, z);
}

sunrealtype N_VMaxNorm_Pthreads(N_Vector x) { return nvMaxNorm_Pthreads(x); }

sunrealtype N_VMinQuotient_Pthreads(N_Vector num, N_Vector denom)
{
  return nvMinQuotient_Pthreads(num, denom);
}

sunrealtype N_VMin_Pthreads(N_Vector x) { return nvMin_Pthreads(x); }

void N_VPrintFile_Pthreads(N_Vector v, FILE* outfile)
{
  nvPrintFile_Pthreads(v, outfile);
}

void N_VPrint_Pthreads(N_Vector v) { nvPrint_Pthreads(v); }

void N_VProd_Pthreads(N_Vector x, N_Vector y, N_Vector z)
{
  nvProd_Pthreads(x, y, z);
}

SUNErrCode N_VScaleAddMultiVectorArray_Pthreads(int nvec, int nsum,
                                                sunrealtype* a, N_Vector* X,
                                                N_Vector** Y, N_Vector** Z)
{
  return nvScaleAddMultiVectorArray_Pthreads(nvec, nsum, a, X, Y, Z);
}

SUNErrCode N_VScaleAddMulti_Pthreads(int nvec, sunrealtype* a, N_Vector x,
                                     N_Vector* Y, N_Vector* Z)
{
  return nvScaleAddMulti_Pthreads(nvec, a, x, Y, Z);
}

SUNErrCode N_VScaleVectorArray_Pthreads(int nvec, sunrealtype* c, N_Vector* X,
                                        N_Vector* Z)
{
  return nvScaleVectorArray_Pthreads(nvec, c, X, Z);
}

void N_VScale_Pthreads(sunrealtype c, N_Vector x, N_Vector z)
{
  nvScale_Pthreads(c, x, z);
}

void N_VSetArrayPointer_Pthreads(sunrealtype* v_data, N_Vector v)
{
  nvSetArrayPointer_Pthreads(v_data, v);
}

sunrealtype N_VWL2Norm_Pthreads(N_Vector x, N_Vector w)
{
  return nvWL2Norm_Pthreads(x, w);
}

sunrealtype N_VWSqrSumLocal_Pthreads(N_Vector x, N_Vector w)
{
  return nvWSqrSumLocal_Pthreads(x, w);
}

sunrealtype N_VWSqrSumMaskLocal_Pthreads(N_Vector x, N_Vector w, N_Vector id)
{
  return nvWSqrSumMaskLocal_Pthreads(x, w, id);
}

SUNErrCode N_VWrmsNormMaskVectorArray_Pthreads(int nvec, N_Vector* X, N_Vector* W,
                                               N_Vector id, sunrealtype* nrm)
{
  return nvWrmsNormMaskVectorArray_Pthreads(nvec, X, W, id, nrm);
}

sunrealtype N_VWrmsNormMask_Pthreads(N_Vector x, N_Vector w, N_Vector id)
{
  return nvWrmsNormMask_Pthreads(x, w, id);
}

SUNErrCode N_VWrmsNormVectorArray_Pthreads(int nvec, N_Vector* X, N_Vector* W,
                                           sunrealtype* nrm)
{
  return nvWrmsNormVectorArray_Pthreads(nvec, X, W, nrm);
}

sunrealtype N_VWrmsNorm_Pthreads(N_Vector x, N_Vector w)
{
  return nvWrmsNorm_Pthreads(x, w);
}
