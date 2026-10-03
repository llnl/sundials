/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * This is the implementation file for a SYCL implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sycl/sycl.hpp>

/* SUNDIALS public headers */
#include <nvector/nvector_sycl_deprecated.h>
#include <sunmemory/sunmemory_sycl.h>
/* SUNDIALS private headers */
#include "sundials/sundials_errors.h"
#include "sundials_debug.h"
#include "sundials_sycl.h"

#define ZERO   SUN_RCONST(0.0)
#define HALF   SUN_RCONST(0.5)
#define ONE    SUN_RCONST(1.0)
#define ONEPT5 SUN_RCONST(1.5)

extern "C" {

using namespace sundials;
using namespace sundials::sycl;

/* --------------------------------------------------------------------------
 * Helpful macros
 * -------------------------------------------------------------------------- */

/* Macros to access vector content */
#define NVEC_SYCL_CONTENT(x) ((N_VectorContent_Sycl)(x->content))
#define NVEC_SYCL_LENGTH(x)  (NVEC_SYCL_CONTENT(x)->length)
#define NVEC_SYCL_MEMHELP(x) (NVEC_SYCL_CONTENT(x)->mem_helper)
#define NVEC_SYCL_MEMSIZE(x) \
  (NVEC_SYCL_CONTENT(x)->length * sizeof(sunrealtype))
#define NVEC_SYCL_HDATAp(x) ((sunrealtype*)NVEC_SYCL_CONTENT(x)->host_data->ptr)
#define NVEC_SYCL_DDATAp(x) \
  ((sunrealtype*)NVEC_SYCL_CONTENT(x)->device_data->ptr)
#define NVEC_SYCL_QUEUE(x) (NVEC_SYCL_CONTENT(x)->queue)

/* Macros to access vector private content */
#define NVEC_SYCL_PRIVATE(x) \
  ((N_PrivateVectorContent_Sycl)(NVEC_SYCL_CONTENT(x)->priv))
#define NVEC_SYCL_HBUFFERp(x) \
  ((sunrealtype*)NVEC_SYCL_PRIVATE(x)->reduce_buffer_host->ptr)
#define NVEC_SYCL_DBUFFERp(x) \
  ((sunrealtype*)NVEC_SYCL_PRIVATE(x)->reduce_buffer_dev->ptr)

/* --------------------------------------------------------------------------
 * Private structure definition
 * -------------------------------------------------------------------------- */

struct _N_PrivateVectorContent_Sycl
{
  sunbooleantype use_managed_mem; /* do data pointers use managed memory */

  /* reduction workspace */
  SUNMemory reduce_buffer_dev;  /* device memory for reductions      */
  SUNMemory reduce_buffer_host; /* host memory for reductions        */
  size_t reduce_buffer_bytes;   /* current size of reduction buffers */

  /* fused op workspace */
  SUNMemory fused_buffer_dev;  /* device memory for fused ops    */
  SUNMemory fused_buffer_host; /* host memory for fused ops      */
  size_t fused_buffer_bytes;   /* current size of the buffers    */
  size_t fused_buffer_offset;  /* current offset into the buffer */
};

typedef struct _N_PrivateVectorContent_Sycl* N_PrivateVectorContent_Sycl;

/* --------------------------------------------------------------------------
 * Utility functions
 * -------------------------------------------------------------------------- */

/* Functions attached to the N_Vector */
static void nvAbs_Sycl(N_Vector x, N_Vector z);
static void nvAddConst_Sycl(N_Vector x, sunrealtype b, N_Vector z);
static SUNErrCode nvBufPack_Sycl(N_Vector x, void* buf);
static SUNErrCode nvBufSize_Sycl(N_Vector x, sunindextype* size);
static SUNErrCode nvBufUnpack_Sycl(N_Vector x, void* buf);
static N_Vector nvCloneEmpty_Sycl(N_Vector w);
static N_Vector nvClone_Sycl(N_Vector w);
static void nvCompare_Sycl(sunrealtype c, N_Vector x, N_Vector z);
static SUNErrCode nvConstVectorArray_Sycl(int nvec, sunrealtype c, N_Vector* Z);
static void nvConst_Sycl(sunrealtype c, N_Vector z);
static sunbooleantype nvConstrMask_Sycl(N_Vector c, N_Vector x, N_Vector m);
static void nvDestroy_Sycl(N_Vector v);
static void nvDiv_Sycl(N_Vector x, N_Vector y, N_Vector z);
static sunrealtype nvDotProd_Sycl(N_Vector x, N_Vector y);

static inline sunrealtype* nvGetDeviceArrayPointer_Sycl(N_Vector x)
{
  N_VectorContent_Sycl content = (N_VectorContent_Sycl)x->content;
  return (content->device_data == NULL ? NULL
                                       : (sunrealtype*)content->device_data->ptr);
}

static inline sunrealtype* nvGetHostArrayPointer_Sycl(N_Vector x)
{
  N_VectorContent_Sycl content = (N_VectorContent_Sycl)x->content;
  return (content->host_data == NULL ? NULL
                                     : (sunrealtype*)content->host_data->ptr);
}

static inline sunindextype nvGetLength_Sycl(N_Vector x)
{
  N_VectorContent_Sycl content = (N_VectorContent_Sycl)x->content;
  return content->length;
}

static inline N_Vector_ID nvGetVectorID_Sycl(N_Vector)
{
  return SUNDIALS_NVEC_SYCL;
}

static sunbooleantype nvInvTest_Sycl(N_Vector x, N_Vector z);
static void nvInv_Sycl(N_Vector x, N_Vector z);
static sunrealtype nvL1Norm_Sycl(N_Vector x);
static SUNErrCode nvLinearCombinationVectorArray_Sycl(int nvec, int nsum,
                                                      sunrealtype* c,
                                                      N_Vector** X, N_Vector* Z);
static SUNErrCode nvLinearCombination_Sycl(int nvec, sunrealtype* c,
                                           N_Vector* X, N_Vector Z);
static SUNErrCode nvLinearSumVectorArray_Sycl(int nvec, sunrealtype a,
                                              N_Vector* X, sunrealtype b,
                                              N_Vector* Y, N_Vector* Z);
static void nvLinearSum_Sycl(sunrealtype a, N_Vector x, sunrealtype b,
                             N_Vector y, N_Vector z);
static sunrealtype nvMaxNorm_Sycl(N_Vector x);
static sunrealtype nvMinQuotient_Sycl(N_Vector num, N_Vector denom);
static sunrealtype nvMin_Sycl(N_Vector x);
static void nvPrintFile_Sycl(N_Vector v, FILE* outfile);
static void nvPrint_Sycl(N_Vector v);
static void nvProd_Sycl(N_Vector x, N_Vector y, N_Vector z);
static SUNErrCode nvScaleAddMultiVectorArray_Sycl(int nvec, int nsum,
                                                  sunrealtype* a, N_Vector* X,
                                                  N_Vector** Y, N_Vector** Z);
static SUNErrCode nvScaleAddMulti_Sycl(int nvec, sunrealtype* c, N_Vector X,
                                       N_Vector* Y, N_Vector* Z);
static SUNErrCode nvScaleVectorArray_Sycl(int nvec, sunrealtype* c, N_Vector* X,
                                          N_Vector* Z);
static void nvScale_Sycl(sunrealtype c, N_Vector x, N_Vector z);
static void nvSetDeviceArrayPointer_Sycl(sunrealtype* d_vdata_1d, N_Vector v);
static void nvSetHostArrayPointer_Sycl(sunrealtype* h_vdata_1d, N_Vector v);
static sunrealtype nvWL2Norm_Sycl(N_Vector x, N_Vector w);
static sunrealtype nvWSqrSumLocal_Sycl(N_Vector x, N_Vector w);
static sunrealtype nvWSqrSumMaskLocal_Sycl(N_Vector x, N_Vector w, N_Vector id);
static sunrealtype nvWrmsNormMask_Sycl(N_Vector x, N_Vector w, N_Vector id);
static sunrealtype nvWrmsNorm_Sycl(N_Vector x, N_Vector w);

/* Allocate vector data */
static int AllocateData(N_Vector v);

/* Reduction buffer functions */
static int InitializeReductionBuffer(N_Vector v, const sunrealtype value,
                                     size_t n = 1);
static void FreeReductionBuffer(N_Vector v);
static int CopyReductionBufferFromDevice(N_Vector v, size_t n = 1);

/* Fused operation buffer functions */
static int FusedBuffer_Init(N_Vector v, int nreal, int nptr);
static int FusedBuffer_CopyRealArray(N_Vector v, sunrealtype* r_data, int nval,
                                     sunrealtype** shortcut);
static int FusedBuffer_CopyPtrArray1D(N_Vector v, N_Vector* X, int nvec,
                                      sunrealtype*** shortcut);
static int FusedBuffer_CopyPtrArray2D(N_Vector v, N_Vector** X, int nvec,
                                      int nsum, sunrealtype*** shortcut);
static int FusedBuffer_CopyToDevice(N_Vector v);
static int FusedBuffer_Free(N_Vector v);

/* Kernel launch parameters */
static int GetKernelParameters(N_Vector v, sunbooleantype reduction,
                               size_t& nthreads_total,
                               size_t& nthreads_per_block);

/* --------------------------------------------------------------------------
 * Constructors
 * -------------------------------------------------------------------------- */

N_Vector N_VNewEmpty_Sycl(SUNContext sunctx)
{
  /* Create an empty vector object */
  N_Vector v = N_VNewEmpty(sunctx);
  if (v == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VNewEmpty_Sycl: N_VNewEmpty returned NULL\n");
    return NULL;
  }

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid           = nvGetVectorID_Sycl;
  v->ops->nvclone                 = nvClone_Sycl;
  v->ops->nvcloneempty            = nvCloneEmpty_Sycl;
  v->ops->nvdestroy               = nvDestroy_Sycl;
  v->ops->nvgetlength             = nvGetLength_Sycl;
  v->ops->nvgetarraypointer       = nvGetHostArrayPointer_Sycl;
  v->ops->nvgetdevicearraypointer = nvGetDeviceArrayPointer_Sycl;
  v->ops->nvsetarraypointer       = nvSetHostArrayPointer_Sycl;
  v->ops->nvsetdevicearraypointer = nvSetDeviceArrayPointer_Sycl;

  /* standard vector operations */
  v->ops->nvlinearsum    = nvLinearSum_Sycl;
  v->ops->nvconst        = nvConst_Sycl;
  v->ops->nvprod         = nvProd_Sycl;
  v->ops->nvdiv          = nvDiv_Sycl;
  v->ops->nvscale        = nvScale_Sycl;
  v->ops->nvabs          = nvAbs_Sycl;
  v->ops->nvinv          = nvInv_Sycl;
  v->ops->nvaddconst     = nvAddConst_Sycl;
  v->ops->nvdotprod      = nvDotProd_Sycl;
  v->ops->nvmaxnorm      = nvMaxNorm_Sycl;
  v->ops->nvmin          = nvMin_Sycl;
  v->ops->nvl1norm       = nvL1Norm_Sycl;
  v->ops->nvinvtest      = nvInvTest_Sycl;
  v->ops->nvconstrmask   = nvConstrMask_Sycl;
  v->ops->nvminquotient  = nvMinQuotient_Sycl;
  v->ops->nvwrmsnormmask = nvWrmsNormMask_Sycl;
  v->ops->nvwrmsnorm     = nvWrmsNorm_Sycl;
  v->ops->nvwl2norm      = nvWL2Norm_Sycl;
  v->ops->nvcompare      = nvCompare_Sycl;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvwsqrsumlocal     = nvWSqrSumLocal_Sycl;
  v->ops->nvwsqrsummasklocal = nvWSqrSumMaskLocal_Sycl;
  v->ops->nvdotprodlocal     = nvDotProd_Sycl;
  v->ops->nvmaxnormlocal     = nvMaxNorm_Sycl;
  v->ops->nvminlocal         = nvMin_Sycl;
  v->ops->nvl1normlocal      = nvL1Norm_Sycl;
  v->ops->nvinvtestlocal     = nvInvTest_Sycl;
  v->ops->nvconstrmasklocal  = nvConstrMask_Sycl;
  v->ops->nvminquotientlocal = nvMinQuotient_Sycl;

  /* XBraid interface operations */
  v->ops->nvbufsize   = nvBufSize_Sycl;
  v->ops->nvbufpack   = nvBufPack_Sycl;
  v->ops->nvbufunpack = nvBufUnpack_Sycl;

  /* print operation for debugging */
  v->ops->nvprint     = nvPrint_Sycl;
  v->ops->nvprintfile = nvPrintFile_Sycl;

  /* Allocate content structure */
  v->content = (N_VectorContent_Sycl)malloc(sizeof(_N_VectorContent_Sycl));
  if (v->content == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VNewEmpty_Sycl: content allocation failed\n");
    N_VDestroy(v);
    return NULL;
  }

  /* Allocate private content structure */
  NVEC_SYCL_CONTENT(v)->priv = NULL;
  NVEC_SYCL_CONTENT(v)->priv = malloc(sizeof(_N_PrivateVectorContent_Sycl));
  if (NVEC_SYCL_CONTENT(v)->priv == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VNewEmpty_Sycl: private content allocation failed\n");
    N_VDestroy(v);
    return NULL;
  }

  /* Initialize content */
  NVEC_SYCL_CONTENT(v)->length             = 0;
  NVEC_SYCL_CONTENT(v)->own_helper         = SUNFALSE;
  NVEC_SYCL_CONTENT(v)->host_data          = NULL;
  NVEC_SYCL_CONTENT(v)->device_data        = NULL;
  NVEC_SYCL_CONTENT(v)->stream_exec_policy = NULL;
  NVEC_SYCL_CONTENT(v)->reduce_exec_policy = NULL;
  NVEC_SYCL_CONTENT(v)->mem_helper         = NULL;
  NVEC_SYCL_CONTENT(v)->queue              = NULL;

  /* Initialize private content */
  NVEC_SYCL_PRIVATE(v)->use_managed_mem     = SUNFALSE;
  NVEC_SYCL_PRIVATE(v)->reduce_buffer_dev   = NULL;
  NVEC_SYCL_PRIVATE(v)->reduce_buffer_host  = NULL;
  NVEC_SYCL_PRIVATE(v)->reduce_buffer_bytes = 0;
  NVEC_SYCL_PRIVATE(v)->fused_buffer_dev    = NULL;
  NVEC_SYCL_PRIVATE(v)->fused_buffer_host   = NULL;
  NVEC_SYCL_PRIVATE(v)->fused_buffer_bytes  = 0;
  NVEC_SYCL_PRIVATE(v)->fused_buffer_offset = 0;

  return v;
}

N_Vector N_VNew_Sycl(sunindextype length, ::sycl::queue* Q, SUNContext sunctx)
{
  /* Check inputs */
  if (Q == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNew_Sycl: queue is NULL\n");
    return NULL;
  }

  if (!(Q->is_in_order()))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNew_Sycl: queue is not in-order\n");
    return NULL;
  }

  /* Create vector with empty content */
  N_Vector v = N_VNewEmpty_Sycl(sunctx);
  if (v == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VNew_Sycl: N_VNewEmpty_Sycl returned NULL\n");
    return NULL;
  }

  /* Fill content */
  NVEC_SYCL_CONTENT(v)->length     = length;
  NVEC_SYCL_CONTENT(v)->own_helper = SUNTRUE;
  NVEC_SYCL_CONTENT(v)->stream_exec_policy =
    new ThreadDirectExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->reduce_exec_policy =
    new BlockReduceExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->mem_helper      = SUNMemoryHelper_Sycl(sunctx);
  NVEC_SYCL_CONTENT(v)->queue           = Q;
  NVEC_SYCL_PRIVATE(v)->use_managed_mem = SUNFALSE;

  if (NVEC_SYCL_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNew_Sycl: memory helper is NULL\n");
    N_VDestroy(v);
    return NULL;
  }

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VNew_Sycl: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return NULL;
  }

  return v;
}

N_Vector N_VNewWithMemHelp_Sycl(sunindextype length,
                                sunbooleantype use_managed_mem,
                                SUNMemoryHelper helper, ::sycl::queue* Q,
                                SUNContext sunctx)
{
  /* Check inputs */
  if (Q == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Sycl: queue is NULL\n");
    return NULL;
  }

  if (!(Q->is_in_order()))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VNewWithMemHelp_Sycl: queue is not in-order\n");
    return NULL;
  }

  if (helper == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Sycl: helper is NULL\n");
    return NULL;
  }

  if (!SUNMemoryHelper_ImplementsRequiredOps(helper))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Sycl: helper doesn't "
                         "implement all required ops\n");
    return NULL;
  }

  /* Create vector with empty content */
  N_Vector v = N_VNewEmpty_Sycl(sunctx);
  if (v == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VNewWithMemHelp_Sycl: N_VNewEmpty_Sycl returned NULL\n");
    return NULL;
  }

  /* Fill content */
  NVEC_SYCL_CONTENT(v)->length     = length;
  NVEC_SYCL_CONTENT(v)->own_helper = SUNFALSE;
  NVEC_SYCL_CONTENT(v)->stream_exec_policy =
    new ThreadDirectExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->reduce_exec_policy =
    new BlockReduceExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->mem_helper      = helper;
  NVEC_SYCL_CONTENT(v)->queue           = Q;
  NVEC_SYCL_PRIVATE(v)->use_managed_mem = use_managed_mem;

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VNewWithMemHelp_Sycl: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return NULL;
  }

  return v;
}

N_Vector N_VNewManaged_Sycl(sunindextype length, ::sycl::queue* Q,
                            SUNContext sunctx)
{
  /* Check inputs */
  if (Q == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewManaged_Sycl: queue is NULL\n");
    return NULL;
  }

  if (!(Q->is_in_order()))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VNewManaged_Sycl: queue is not in-order\n");
    return NULL;
  }

  /* Create vector with empty content */
  N_Vector v = N_VNewEmpty_Sycl(sunctx);
  if (v == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VNewManaged_Sycl: N_VNewEmpty_Sycl returned NULL\n");
    return NULL;
  }

  /* Fill content */
  NVEC_SYCL_CONTENT(v)->length     = length;
  NVEC_SYCL_CONTENT(v)->own_helper = SUNTRUE;
  NVEC_SYCL_CONTENT(v)->stream_exec_policy =
    new ThreadDirectExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->reduce_exec_policy =
    new BlockReduceExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->mem_helper      = SUNMemoryHelper_Sycl(sunctx);
  NVEC_SYCL_CONTENT(v)->queue           = Q;
  NVEC_SYCL_PRIVATE(v)->use_managed_mem = SUNTRUE;

  if (NVEC_SYCL_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VNewManaged_Sycl: memory helper is NULL\n");
    N_VDestroy(v);
    return NULL;
  }

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VNewManaged_Sycl: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return NULL;
  }

  return v;
}

N_Vector N_VMake_Sycl(sunindextype length, sunrealtype* h_vdata,
                      sunrealtype* d_vdata, ::sycl::queue* Q, SUNContext sunctx)
{
  /* Check inputs */
  if (Q == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMake_Sycl: queue is NULL\n");
    return NULL;
  }

  if (!(Q->is_in_order()))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMake_Sycl: queue is not in-order\n");
    return NULL;
  }

  if (h_vdata == NULL || d_vdata == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VMake_Sycl: host or device data is NULL\n");
    return NULL;
  }

  /* Create vector with empty content */
  N_Vector v = N_VNewEmpty_Sycl(sunctx);
  if (v == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMake_Sycl: N_VMake_Sycl returned NULL\n");
    return NULL;
  }

  /* Fill content */
  NVEC_SYCL_CONTENT(v)->length     = length;
  NVEC_SYCL_CONTENT(v)->own_helper = SUNTRUE;
  NVEC_SYCL_CONTENT(v)->mem_helper = SUNMemoryHelper_Sycl(sunctx);
  NVEC_SYCL_CONTENT(v)->host_data =
    SUNMemoryHelper_Wrap(NVEC_SYCL_MEMHELP(v), h_vdata, SUNMEMTYPE_HOST);
  NVEC_SYCL_CONTENT(v)->device_data =
    SUNMemoryHelper_Wrap(NVEC_SYCL_MEMHELP(v), d_vdata, SUNMEMTYPE_DEVICE);
  NVEC_SYCL_CONTENT(v)->stream_exec_policy =
    new ThreadDirectExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->reduce_exec_policy =
    new BlockReduceExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->queue           = Q;
  NVEC_SYCL_PRIVATE(v)->use_managed_mem = SUNFALSE;

  if (NVEC_SYCL_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMake_Sycl: memory helper is NULL\n");
    N_VDestroy(v);
    return NULL;
  }

  if (NVEC_SYCL_CONTENT(v)->device_data == NULL ||
      NVEC_SYCL_CONTENT(v)->host_data == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VMake_Sycl: SUNMemoryHelper_Wrap returned NULL\n");
    N_VDestroy(v);
    return NULL;
  }

  return v;
}

N_Vector N_VMakeManaged_Sycl(sunindextype length, sunrealtype* vdata,
                             ::sycl::queue* Q, SUNContext sunctx)
{
  /* Check inputs */
  if (Q == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMakeManaged_Sycl: queue is NULL\n");
    return NULL;
  }

  if (!(Q->is_in_order()))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VMakeManaged_Sycl: queue is not in-order\n");
    return NULL;
  }

  if (vdata == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VMakeManaged_Sycl: host or device data is NULL\n");
    return NULL;
  }

  /* Create vector with empty content */
  N_Vector v = N_VNewEmpty_Sycl(sunctx);
  if (v == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VMakeManaged_Sycl: N_VNewEmpty_Sycl returned NULL\n");
    return NULL;
  }

  /* Fill content */
  NVEC_SYCL_CONTENT(v)->length     = length;
  NVEC_SYCL_CONTENT(v)->mem_helper = SUNMemoryHelper_Sycl(sunctx);
  NVEC_SYCL_CONTENT(v)->own_helper = SUNTRUE;
  NVEC_SYCL_CONTENT(v)->host_data  = SUNMemoryHelper_Wrap(NVEC_SYCL_MEMHELP(v),
                                                          vdata, SUNMEMTYPE_UVM);
  NVEC_SYCL_CONTENT(v)->device_data =
    SUNMemoryHelper_Alias(NVEC_SYCL_MEMHELP(v), NVEC_SYCL_CONTENT(v)->host_data);
  NVEC_SYCL_CONTENT(v)->stream_exec_policy =
    new ThreadDirectExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->reduce_exec_policy =
    new BlockReduceExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->queue           = Q;
  NVEC_SYCL_PRIVATE(v)->use_managed_mem = SUNTRUE;

  if (NVEC_SYCL_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VMakeManaged_Sycl: memory helper is NULL\n");
    N_VDestroy(v);
    return NULL;
  }

  if (NVEC_SYCL_CONTENT(v)->device_data == NULL ||
      NVEC_SYCL_CONTENT(v)->host_data == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VMakeManaged_Sycl: SUNMemoryHelper_Wrap returned NULL\n");
    N_VDestroy(v);
    return NULL;
  }

  return v;
}

/* --------------------------------------------------------------------------
 * Vector get, set, and utility functions
 * -------------------------------------------------------------------------- */

/* Function to return the global length of the vector. This is defined as an
 * inline function in nvector_sycl.h, so we just mark it as extern here. */
extern sunindextype nvGetLength_Sycl(N_Vector v);

/* Return pointer to the raw host data. This is defined as an inline function in
 * nvector_sycl.h, so we just mark it as extern here. */
extern sunrealtype* nvGetHostArrayPointer_Sycl(N_Vector x);

/* Return pointer to the raw device data. This is defined as an inline function
 * in nvector_sycl.h, so we just mark it as extern here. */
extern sunrealtype* nvGetDeviceArrayPointer_Sycl(N_Vector x);

/* Set pointer to the raw host data. Does not free the existing pointer. */
void nvSetHostArrayPointer_Sycl(sunrealtype* h_vdata, N_Vector v)
{
  if (N_VIsManagedMemory_Sycl(v))
  {
    if (NVEC_SYCL_CONTENT(v)->host_data)
    {
      NVEC_SYCL_CONTENT(v)->host_data->ptr   = (void*)h_vdata;
      NVEC_SYCL_CONTENT(v)->device_data->ptr = (void*)h_vdata;
    }
    else
    {
      NVEC_SYCL_CONTENT(v)->host_data =
        SUNMemoryHelper_Wrap(NVEC_SYCL_MEMHELP(v), (void*)h_vdata,
                             SUNMEMTYPE_UVM);
      NVEC_SYCL_CONTENT(v)->device_data =
        SUNMemoryHelper_Alias(NVEC_SYCL_MEMHELP(v),
                              NVEC_SYCL_CONTENT(v)->host_data);
    }
  }
  else
  {
    if (NVEC_SYCL_CONTENT(v)->host_data)
    {
      NVEC_SYCL_CONTENT(v)->host_data->ptr = (void*)h_vdata;
    }
    else
    {
      NVEC_SYCL_CONTENT(v)->host_data =
        SUNMemoryHelper_Wrap(NVEC_SYCL_MEMHELP(v), (void*)h_vdata,
                             SUNMEMTYPE_HOST);
    }
  }
}

/* Set pointer to the raw device data */
void nvSetDeviceArrayPointer_Sycl(sunrealtype* d_vdata, N_Vector v)
{
  if (N_VIsManagedMemory_Sycl(v))
  {
    if (NVEC_SYCL_CONTENT(v)->device_data)
    {
      NVEC_SYCL_CONTENT(v)->device_data->ptr = (void*)d_vdata;
      NVEC_SYCL_CONTENT(v)->host_data->ptr   = (void*)d_vdata;
    }
    else
    {
      NVEC_SYCL_CONTENT(v)->device_data =
        SUNMemoryHelper_Wrap(NVEC_SYCL_MEMHELP(v), (void*)d_vdata,
                             SUNMEMTYPE_UVM);
      NVEC_SYCL_CONTENT(v)->host_data =
        SUNMemoryHelper_Alias(NVEC_SYCL_MEMHELP(v),
                              NVEC_SYCL_CONTENT(v)->device_data);
    }
  }
  else
  {
    if (NVEC_SYCL_CONTENT(v)->device_data)
    {
      NVEC_SYCL_CONTENT(v)->device_data->ptr = (void*)d_vdata;
    }
    else
    {
      NVEC_SYCL_CONTENT(v)->device_data =
        SUNMemoryHelper_Wrap(NVEC_SYCL_MEMHELP(v), (void*)d_vdata,
                             SUNMEMTYPE_DEVICE);
    }
  }
}

/* Return a flag indicating if the memory for the vector data is managed */
sunbooleantype N_VIsManagedMemory_Sycl(N_Vector x)
{
  return NVEC_SYCL_PRIVATE(x)->use_managed_mem;
}

SUNErrCode N_VSetKernelExecPolicy_Sycl(N_Vector x,
                                       SUNSyclExecPolicy* stream_exec_policy,
                                       SUNSyclExecPolicy* reduce_exec_policy)
{
  if (!x)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VSetKernelExecPolicy_Sycl: The input vector is NULL\n");
    return SUN_ERR_GENERIC;
  }

  /* Delete the old policies */
  delete NVEC_SYCL_CONTENT(x)->stream_exec_policy;
  delete NVEC_SYCL_CONTENT(x)->reduce_exec_policy;

  /* Clone input policy or use default */
  if (stream_exec_policy)
  {
    NVEC_SYCL_CONTENT(x)->stream_exec_policy = stream_exec_policy->clone();
  }
  else
    NVEC_SYCL_CONTENT(x)->stream_exec_policy =
      new ThreadDirectExecPolicy(SYCL_BLOCKDIM(NVEC_SYCL_QUEUE(x)));

  if (reduce_exec_policy)
  {
    NVEC_SYCL_CONTENT(x)->reduce_exec_policy = reduce_exec_policy->clone();
  }
  else
    NVEC_SYCL_CONTENT(x)->reduce_exec_policy =
      new BlockReduceExecPolicy(SYCL_BLOCKDIM(NVEC_SYCL_QUEUE(x)));

  return SUN_SUCCESS;
}

/* Copy vector data to the device */
void N_VCopyToDevice_Sycl(N_Vector x)
{
  int copy_fail;

  copy_fail = SUNMemoryHelper_Copy(NVEC_SYCL_MEMHELP(x),
                                   NVEC_SYCL_CONTENT(x)->device_data,
                                   NVEC_SYCL_CONTENT(x)->host_data,
                                   NVEC_SYCL_MEMSIZE(x), NVEC_SYCL_QUEUE(x));

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VCopyToDevice_Sycl: SUNMemoryHelper_Copy returned nonzero\n");
  }

  /* synchronize with the host */
  NVEC_SYCL_QUEUE(x)->wait_and_throw();
}

/* Copy vector data from the device to the host */
void N_VCopyFromDevice_Sycl(N_Vector x)
{
  int copy_fail;

  copy_fail = SUNMemoryHelper_Copy(NVEC_SYCL_MEMHELP(x),
                                   NVEC_SYCL_CONTENT(x)->host_data,
                                   NVEC_SYCL_CONTENT(x)->device_data,
                                   NVEC_SYCL_MEMSIZE(x), NVEC_SYCL_QUEUE(x));

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VCopyFromDevice_Sycl: "
                         "SUNMemoryHelper_Copy returned nonzero\n");
  }

  /* synchronize with the host */
  NVEC_SYCL_QUEUE(x)->wait_and_throw();
}

/* Function to print the a serial vector to stdout */
void nvPrint_Sycl(N_Vector X) { nvPrintFile_Sycl(X, stdout); }

/* Function to print the a serial vector to outfile */
void nvPrintFile_Sycl(N_Vector X, FILE* outfile)
{
  sunindextype i;

  for (i = 0; i < NVEC_SYCL_CONTENT(X)->length; i++)
  {
    fprintf(outfile, SUN_FORMAT_E "\n", NVEC_SYCL_HDATAp(X)[i]);
  }

  return;
}

/* --------------------------------------------------------------------------
 * Vector operations
 * -------------------------------------------------------------------------- */

N_Vector nvCloneEmpty_Sycl(N_Vector w)
{
  /* Check input */
  if (w == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvCloneEmpty_Sycl: input vector is NULL\n");
    return NULL;
  }

  /* Create vector */
  N_Vector v = N_VNewEmpty_Sycl(w->sunctx);
  if (v == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvCloneEmpty_Sycl: N_VNewEmpty returned NULL\n");
    return NULL;
  }

  /* Attach operations */
  if (N_VCopyOps(w, v))
  {
    N_VDestroy(v);
    SUNDIALS_DEBUG_PRINT("ERROR in nvCloneEmpty_Sycl: Error in N_VCopyOps\n");
    return NULL;
  }

  /* Copy some content */
  NVEC_SYCL_CONTENT(v)->length          = NVEC_SYCL_CONTENT(w)->length;
  NVEC_SYCL_CONTENT(v)->queue           = NVEC_SYCL_CONTENT(w)->queue;
  NVEC_SYCL_PRIVATE(v)->use_managed_mem = NVEC_SYCL_PRIVATE(w)->use_managed_mem;

  return v;
}

N_Vector nvClone_Sycl(N_Vector w)
{
  /* Check inputs */
  if (w == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvClone_Sycl: vector is NULL\n");
    return NULL;
  }

  /* Create an empty clone vector */
  N_Vector v = nvCloneEmpty_Sycl(w);
  if (v == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvClone_Sycl: nvCloneEmpty_Sycl returned NULL\n");
    return NULL;
  }

  /* Allocate content */
  NVEC_SYCL_CONTENT(v)->own_helper = SUNTRUE;
  NVEC_SYCL_CONTENT(v)->stream_exec_policy =
    NVEC_SYCL_CONTENT(w)->stream_exec_policy->clone();
  NVEC_SYCL_CONTENT(v)->reduce_exec_policy =
    NVEC_SYCL_CONTENT(w)->reduce_exec_policy->clone();
  NVEC_SYCL_CONTENT(v)->mem_helper = SUNMemoryHelper_Clone(NVEC_SYCL_MEMHELP(w));

  if (NVEC_SYCL_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvClone_Sycl: memory helper is NULL\n");
    N_VDestroy(v);
    return NULL;
  }

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvClone_Sycl: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return NULL;
  }

  return v;
}

void nvDestroy_Sycl(N_Vector v)
{
  N_VectorContent_Sycl vc;
  N_PrivateVectorContent_Sycl vcp;

  if (v == NULL) { return; }

  /* free ops structure */
  if (v->ops != NULL)
  {
    free(v->ops);
    v->ops = NULL;
  }

  /* extract content */
  vc = NVEC_SYCL_CONTENT(v);
  if (vc == NULL)
  {
    free(v);
    v = NULL;
    return;
  }

  /* free private content */
  vcp = (N_PrivateVectorContent_Sycl)vc->priv;
  if (vcp != NULL)
  {
    /* free items in private content */
    FreeReductionBuffer(v);
    FusedBuffer_Free(v);
    free(vcp);
    vc->priv = NULL;
  }

  /* free items in content */
  if (NVEC_SYCL_MEMHELP(v))
  {
    SUNMemoryHelper_Dealloc(NVEC_SYCL_MEMHELP(v), vc->host_data,
                            NVEC_SYCL_QUEUE(v));
    vc->host_data = NULL;
    SUNMemoryHelper_Dealloc(NVEC_SYCL_MEMHELP(v), vc->device_data,
                            NVEC_SYCL_QUEUE(v));
    vc->device_data = NULL;
    if (vc->own_helper) { SUNMemoryHelper_Destroy(vc->mem_helper); }
    vc->mem_helper = NULL;
  }
  else
  {
    SUNDIALS_DEBUG_PRINT(
      "WARNING in nvDestroy_Sycl: mem_helper was NULL when trying to dealloc "
      "data, this could result in a memory leak\n");
  }

  delete vc->stream_exec_policy;
  delete vc->reduce_exec_policy;

  vc->stream_exec_policy = nullptr;
  vc->reduce_exec_policy = nullptr;

  /* free content struct and vector */
  free(vc);
  free(v);

  return;
}

void nvConst_Sycl(sunrealtype c, N_Vector z)
{
  const sunindextype N = NVEC_SYCL_LENGTH(z);
  sunrealtype* zdata   = NVEC_SYCL_DDATAp(z);
  ::sycl::queue* Q     = NVEC_SYCL_QUEUE(z);
  size_t nthreads_total, nthreads_per_block;

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvConst_Sycl: GetKernelParameters returned nonzero\n");
  }

  SYCL_FOR(
    Q, nthreads_total, nthreads_per_block, item,
    GRID_STRIDE_XLOOP(item, i, N) { zdata[i] = c; });
}

void nvLinearSum_Sycl(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                      N_Vector z)
{
  const sunindextype N     = NVEC_SYCL_LENGTH(z);
  const sunrealtype* xdata = NVEC_SYCL_DDATAp(x);
  const sunrealtype* ydata = NVEC_SYCL_DDATAp(y);
  sunrealtype* zdata       = NVEC_SYCL_DDATAp(z);
  ::sycl::queue* Q         = NVEC_SYCL_QUEUE(z);
  size_t nthreads_total, nthreads_per_block;

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvLinearSum_Sycl: GetKernelParameters returned nonzero\n");
  }

  SYCL_FOR(
    Q, nthreads_total, nthreads_per_block, item, GRID_STRIDE_XLOOP(item, i, N) {
      zdata[i] = (a * xdata[i]) + (b * ydata[i]);
    });
}

void nvProd_Sycl(N_Vector x, N_Vector y, N_Vector z)
{
  const sunindextype N     = NVEC_SYCL_LENGTH(z);
  const sunrealtype* xdata = NVEC_SYCL_DDATAp(x);
  const sunrealtype* ydata = NVEC_SYCL_DDATAp(y);
  sunrealtype* zdata       = NVEC_SYCL_DDATAp(z);
  ::sycl::queue* Q         = NVEC_SYCL_QUEUE(z);
  size_t nthreads_total, nthreads_per_block;

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvProd_Sycl: GetKernelParameters returned nonzero\n");
  }

  SYCL_FOR(
    Q, nthreads_total, nthreads_per_block, item,
    GRID_STRIDE_XLOOP(item, i, N) { zdata[i] = xdata[i] * ydata[i]; });
}

void nvDiv_Sycl(N_Vector x, N_Vector y, N_Vector z)
{
  const sunindextype N     = NVEC_SYCL_LENGTH(z);
  const sunrealtype* xdata = NVEC_SYCL_DDATAp(x);
  const sunrealtype* ydata = NVEC_SYCL_DDATAp(y);
  sunrealtype* zdata       = NVEC_SYCL_DDATAp(z);
  ::sycl::queue* Q         = NVEC_SYCL_QUEUE(z);
  size_t nthreads_total, nthreads_per_block;

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvDiv_Sycl: GetKernelParameters returned nonzero\n");
  }

  SYCL_FOR(
    Q, nthreads_total, nthreads_per_block, item,
    GRID_STRIDE_XLOOP(item, i, N) { zdata[i] = xdata[i] / ydata[i]; });
}

void nvScale_Sycl(sunrealtype c, N_Vector x, N_Vector z)
{
  const sunindextype N     = NVEC_SYCL_LENGTH(z);
  const sunrealtype* xdata = NVEC_SYCL_DDATAp(x);
  sunrealtype* zdata       = NVEC_SYCL_DDATAp(z);
  ::sycl::queue* Q         = NVEC_SYCL_QUEUE(z);
  size_t nthreads_total, nthreads_per_block;

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvScale_Sycl: GetKernelParameters returned nonzero\n");
  }

  SYCL_FOR(
    Q, nthreads_total, nthreads_per_block, item,
    GRID_STRIDE_XLOOP(item, i, N) { zdata[i] = c * xdata[i]; });
}

void nvAbs_Sycl(N_Vector x, N_Vector z)
{
  const sunindextype N     = NVEC_SYCL_LENGTH(z);
  const sunrealtype* xdata = NVEC_SYCL_DDATAp(x);
  sunrealtype* zdata       = NVEC_SYCL_DDATAp(z);
  ::sycl::queue* Q         = NVEC_SYCL_QUEUE(z);
  size_t nthreads_total, nthreads_per_block;

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvAbs_Sycl: GetKernelParameters returned nonzero\n");
  }

  SYCL_FOR(
    Q, nthreads_total, nthreads_per_block, item,
    GRID_STRIDE_XLOOP(item, i, N) { zdata[i] = abs(xdata[i]); });
}

void nvInv_Sycl(N_Vector x, N_Vector z)
{
  const sunindextype N     = NVEC_SYCL_LENGTH(z);
  const sunrealtype* xdata = NVEC_SYCL_DDATAp(x);
  sunrealtype* zdata       = NVEC_SYCL_DDATAp(z);
  ::sycl::queue* Q         = NVEC_SYCL_QUEUE(z);
  size_t nthreads_total, nthreads_per_block;

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvInv_Sycl: GetKernelParameters returned nonzero\n");
  }

  SYCL_FOR(
    Q, nthreads_total, nthreads_per_block, item,
    GRID_STRIDE_XLOOP(item, i, N) { zdata[i] = ONE / xdata[i]; });
}

void nvAddConst_Sycl(N_Vector x, sunrealtype b, N_Vector z)
{
  const sunindextype N     = NVEC_SYCL_LENGTH(z);
  const sunrealtype* xdata = NVEC_SYCL_DDATAp(x);
  sunrealtype* zdata       = NVEC_SYCL_DDATAp(z);
  ::sycl::queue* Q         = NVEC_SYCL_QUEUE(z);
  size_t nthreads_total, nthreads_per_block;

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvAddConst_Sycl: GetKernelParameters returned nonzero\n");
  }

  SYCL_FOR(
    Q, nthreads_total, nthreads_per_block, item,
    GRID_STRIDE_XLOOP(item, i, N) { zdata[i] = xdata[i] + b; });
}

sunrealtype nvDotProd_Sycl(N_Vector x, N_Vector y)
{
  const sunindextype N     = NVEC_SYCL_LENGTH(x);
  const sunrealtype* xdata = NVEC_SYCL_DDATAp(x);
  const sunrealtype* ydata = NVEC_SYCL_DDATAp(y);
  ::sycl::queue* Q         = NVEC_SYCL_QUEUE(x);
  size_t nthreads_total, nthreads_per_block;

  if (InitializeReductionBuffer(x, ZERO))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvDotProd_Sycl: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(x, SUNTRUE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvDotProd_Sycl: GetKernelParameters returned nonzero\n");
  }

  /* Shortcut to the reduction buffer */
  sunrealtype* sum = NVEC_SYCL_DBUFFERp(x);

  SYCL_FOR_REDUCE(
    Q, nthreads_total, nthreads_per_block, item, sum, ::sycl::plus<sunrealtype>(),
    GRID_STRIDE_XLOOP(item, i, N) { sum += xdata[i] * ydata[i]; });

  if (CopyReductionBufferFromDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvDotProd_Sycl: "
                         "CopyReductionBufferFromDevice returned nonzero\n");
  }

  return NVEC_SYCL_HBUFFERp(x)[0];
}

sunrealtype nvMaxNorm_Sycl(N_Vector x)
{
  const sunindextype N     = NVEC_SYCL_LENGTH(x);
  const sunrealtype* xdata = NVEC_SYCL_DDATAp(x);
  ::sycl::queue* Q         = NVEC_SYCL_QUEUE(x);
  size_t nthreads_total, nthreads_per_block;

  if (InitializeReductionBuffer(x, ZERO))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvMaxNorm_Sycl: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(x, SUNTRUE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvMaxNorm_Sycl: GetKernelParameters returned nonzero\n");
  }

  /* Shortcut to the reduction buffer */
  sunrealtype* max = NVEC_SYCL_DBUFFERp(x);

  SYCL_FOR_REDUCE(
    Q, nthreads_total, nthreads_per_block, item, max,
    ::sycl::maximum<sunrealtype>(),
    GRID_STRIDE_XLOOP(item, i, N) { max.combine(abs(xdata[i])); });

  if (CopyReductionBufferFromDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvMaxNorm_Sycl: "
                         "CopyReductionBufferFromDevice returned nonzero\n");
  }

  return NVEC_SYCL_HBUFFERp(x)[0];
}

sunrealtype nvWSqrSumLocal_Sycl(N_Vector x, N_Vector w)
{
  const sunindextype N     = NVEC_SYCL_LENGTH(x);
  const sunrealtype* xdata = NVEC_SYCL_DDATAp(x);
  const sunrealtype* wdata = NVEC_SYCL_DDATAp(w);
  ::sycl::queue* Q         = NVEC_SYCL_QUEUE(x);
  size_t nthreads_total, nthreads_per_block;

  if (InitializeReductionBuffer(x, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWSqrSumLocal_Sycl: "
                         "InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(x, SUNTRUE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvWSqrSumLocal_Sycl: GetKernelParameters returned nonzero\n");
  }

  /* Shortcut to the reduction buffer */
  sunrealtype* sum = NVEC_SYCL_DBUFFERp(x);

  SYCL_FOR_REDUCE(
    Q, nthreads_total, nthreads_per_block, item, sum,
    ::sycl::plus<sunrealtype>(), GRID_STRIDE_XLOOP(item, i, N) {
      sum += xdata[i] * wdata[i] * xdata[i] * wdata[i];
    });

  if (CopyReductionBufferFromDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWSqrSumLocal_Sycl: "
                         "CopyReductionBufferFromDevice returned nonzero\n");
  }

  return NVEC_SYCL_HBUFFERp(x)[0];
}

sunrealtype nvWrmsNorm_Sycl(N_Vector x, N_Vector w)
{
  const sunindextype N  = NVEC_SYCL_LENGTH(x);
  const sunrealtype sum = nvWSqrSumLocal_Sycl(x, w);
  return ::sycl::sqrt(sum / N);
}

sunrealtype nvWSqrSumMaskLocal_Sycl(N_Vector x, N_Vector w, N_Vector id)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(x);
  const sunrealtype* xdata  = NVEC_SYCL_DDATAp(x);
  const sunrealtype* wdata  = NVEC_SYCL_DDATAp(w);
  const sunrealtype* iddata = NVEC_SYCL_DDATAp(id);
  ::sycl::queue* Q          = NVEC_SYCL_QUEUE(x);
  size_t nthreads_total, nthreads_per_block;

  if (InitializeReductionBuffer(x, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWSqrSumMaskLocal_Sycl: "
                         "InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(x, SUNTRUE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWSqrSumMaskLocal_Sycl: "
                         "GetKernelParameters returned nonzero\n");
  }

  /* Shortcut to the reduction buffer */
  sunrealtype* sum = NVEC_SYCL_DBUFFERp(x);

  SYCL_FOR_REDUCE(
    Q, nthreads_total, nthreads_per_block, item, sum,
    ::sycl::plus<sunrealtype>(), GRID_STRIDE_XLOOP(item, i, N) {
      if (iddata[i] > ZERO) sum += xdata[i] * wdata[i] * xdata[i] * wdata[i];
    });

  if (CopyReductionBufferFromDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWSqrSumMaskLocal_Sycl: "
                         "CopyReductionBufferFromDevice returned nonzero\n");
  }

  return NVEC_SYCL_HBUFFERp(x)[0];
}

sunrealtype nvWrmsNormMask_Sycl(N_Vector x, N_Vector w, N_Vector id)
{
  const sunindextype N  = NVEC_SYCL_LENGTH(x);
  const sunrealtype sum = nvWSqrSumMaskLocal_Sycl(x, w, id);
  return ::sycl::sqrt(sum / N);
}

sunrealtype nvMin_Sycl(N_Vector x)
{
  const sunindextype N     = NVEC_SYCL_LENGTH(x);
  const sunrealtype* xdata = NVEC_SYCL_DDATAp(x);
  ::sycl::queue* Q         = NVEC_SYCL_QUEUE(x);
  size_t nthreads_total, nthreads_per_block;

  if (InitializeReductionBuffer(x, std::numeric_limits<sunrealtype>::max()))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvMin_Sycl: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(x, SUNTRUE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvMin_Sycl: GetKernelParameters returned nonzero\n");
  }

  /* Shortcut to the reduction buffer */
  sunrealtype* min = NVEC_SYCL_DBUFFERp(x);

  SYCL_FOR_REDUCE(
    Q, nthreads_total, nthreads_per_block, item, min,
    ::sycl::minimum<sunrealtype>(),
    GRID_STRIDE_XLOOP(item, i, N) { min.combine(xdata[i]); });

  if (CopyReductionBufferFromDevice(x))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvMin_Sycl: CopyReductionBufferFromDevice returned nonzero\n");
  }

  return NVEC_SYCL_HBUFFERp(x)[0];
}

sunrealtype nvWL2Norm_Sycl(N_Vector x, N_Vector w)
{
  return ::sycl::sqrt(nvWSqrSumLocal_Sycl(x, w));
}

sunrealtype nvL1Norm_Sycl(N_Vector x)
{
  const sunindextype N     = NVEC_SYCL_LENGTH(x);
  const sunrealtype* xdata = NVEC_SYCL_DDATAp(x);
  ::sycl::queue* Q         = NVEC_SYCL_QUEUE(x);
  size_t nthreads_total, nthreads_per_block;

  if (InitializeReductionBuffer(x, ZERO))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvL1Norm_Sycl: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(x, SUNTRUE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvL1Norm_Sycl: GetKernelParameters returned nonzero\n");
  }

  /* Shortcut to the reduction buffer */
  sunrealtype* sum = NVEC_SYCL_DBUFFERp(x);

  SYCL_FOR_REDUCE(
    Q, nthreads_total, nthreads_per_block, item, sum, ::sycl::plus<sunrealtype>(),
    GRID_STRIDE_XLOOP(item, i, N) { sum += abs(xdata[i]); });

  if (CopyReductionBufferFromDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvL1Norm_Sycl: "
                         "CopyReductionBufferFromDevice returned nonzero\n");
  }

  return NVEC_SYCL_HBUFFERp(x)[0];
}

void nvCompare_Sycl(sunrealtype c, N_Vector x, N_Vector z)
{
  const sunindextype N     = NVEC_SYCL_LENGTH(z);
  const sunrealtype* xdata = NVEC_SYCL_DDATAp(x);
  sunrealtype* zdata       = NVEC_SYCL_DDATAp(z);
  ::sycl::queue* Q         = NVEC_SYCL_QUEUE(z);
  size_t nthreads_total, nthreads_per_block;

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvCompare_Sycl: GetKernelParameters returned nonzero\n");
  }

  SYCL_FOR(
    Q, nthreads_total, nthreads_per_block, item, GRID_STRIDE_XLOOP(item, i, N) {
      zdata[i] = abs(xdata[i]) >= c ? ONE : ZERO;
    });
}

sunbooleantype nvInvTest_Sycl(N_Vector x, N_Vector z)
{
  const sunindextype N     = NVEC_SYCL_LENGTH(z);
  const sunrealtype* xdata = NVEC_SYCL_DDATAp(x);
  sunrealtype* zdata       = NVEC_SYCL_DDATAp(z);
  ::sycl::queue* Q         = NVEC_SYCL_QUEUE(z);
  size_t nthreads_total, nthreads_per_block;

  if (InitializeReductionBuffer(x, ZERO))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvInvTest_Sycl: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(x, SUNTRUE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvInvTest_Sycl: GetKernelParameters returned nonzero\n");
  }

  /* Shortcut to the reduction buffer */
  sunrealtype* sum = NVEC_SYCL_DBUFFERp(x);

  SYCL_FOR_REDUCE(
    Q, nthreads_total, nthreads_per_block, item, sum,
    ::sycl::plus<sunrealtype>(), GRID_STRIDE_XLOOP(item, i, N) {
      if (xdata[i] == ZERO) { sum += ONE; }
      else { zdata[i] = ONE / xdata[i]; }
    });

  if (CopyReductionBufferFromDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvInvTest_Sycl: "
                         "CopyReductionBufferFromDevice returned nonzero\n");
  }

  return (NVEC_SYCL_HBUFFERp(x)[0] < HALF);
}

sunbooleantype nvConstrMask_Sycl(N_Vector c, N_Vector x, N_Vector m)
{
  const sunindextype N     = NVEC_SYCL_LENGTH(x);
  const sunrealtype* cdata = NVEC_SYCL_DDATAp(c);
  const sunrealtype* xdata = NVEC_SYCL_DDATAp(x);
  sunrealtype* mdata       = NVEC_SYCL_DDATAp(m);
  ::sycl::queue* Q         = NVEC_SYCL_QUEUE(x);
  size_t nthreads_total, nthreads_per_block;

  if (InitializeReductionBuffer(x, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvConstrMask_Sycl: "
                         "InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(x, SUNTRUE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvConstrMask_Sycl: GetKernelParameters returned nonzero\n");
  }

  /* Shortcut to the reduction buffer */
  sunrealtype* sum = NVEC_SYCL_DBUFFERp(x);

  SYCL_FOR_REDUCE(
    Q, nthreads_total, nthreads_per_block, item, sum,
    ::sycl::plus<sunrealtype>(), GRID_STRIDE_XLOOP(item, i, N) {
      bool test = (abs(cdata[i]) > ONEPT5 && cdata[i] * xdata[i] <= ZERO) ||
                  (abs(cdata[i]) > HALF && cdata[i] * xdata[i] < ZERO);
      mdata[i] = test ? ONE : ZERO;
      sum += mdata[i];
    });

  if (CopyReductionBufferFromDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvConstrMask_Sycl: "
                         "CopyReductionBufferFromDevice returned nonzero\n");
  }

  return (NVEC_SYCL_HBUFFERp(x)[0] < HALF);
}

sunrealtype nvMinQuotient_Sycl(N_Vector num, N_Vector denom)
{
  const sunindextype N     = NVEC_SYCL_LENGTH(num);
  const sunrealtype* ndata = NVEC_SYCL_DDATAp(num);
  const sunrealtype* ddata = NVEC_SYCL_DDATAp(denom);
  ::sycl::queue* Q         = NVEC_SYCL_QUEUE(num);
  size_t nthreads_total, nthreads_per_block;

  if (InitializeReductionBuffer(num, std::numeric_limits<sunrealtype>::max()))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvMinQuotient_Sycl: "
                         "InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(num, SUNTRUE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvMinQuotient_Sycl: GetKernelParameters returned nonzero\n");
  }

  /* Shortcut to the reduction buffer */
  sunrealtype* min = NVEC_SYCL_DBUFFERp(num);

  SYCL_FOR_REDUCE(
    Q, nthreads_total, nthreads_per_block, item, min,
    ::sycl::minimum<sunrealtype>(), GRID_STRIDE_XLOOP(item, i, N) {
      if (ddata[i] != ZERO) min.combine(ndata[i] / ddata[i]);
    });

  if (CopyReductionBufferFromDevice(num))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvMinQuotient_Sycl: "
                         "CopyReductionBufferFromDevice returned nonzero\n");
  }

  return NVEC_SYCL_HBUFFERp(num)[0];
}

/* --------------------------------------------------------------------------
 * fused vector operations
 * -------------------------------------------------------------------------- */

SUNErrCode nvLinearCombination_Sycl(int nvec, sunrealtype* c, N_Vector* X,
                                    N_Vector z)
{
  const sunindextype N = NVEC_SYCL_LENGTH(z);
  sunrealtype* zdata   = NVEC_SYCL_DDATAp(z);
  ::sycl::queue* Q     = NVEC_SYCL_QUEUE(z);
  size_t nthreads_total, nthreads_per_block;

  /* Fused op workspace shortcuts */
  sunrealtype* cdata  = NULL;
  sunrealtype** xdata = NULL;

  /* Setup the fused op workspace */
  if (FusedBuffer_Init(z, nvec, nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombination_Sycl: FusedBuffer_Init "
                         "returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyRealArray(z, c, nvec, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombination_Sycl: "
                         "FusedBuffer_CopyRealArray returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(z, X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombination_Sycl: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyToDevice(z))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombination_Sycl: "
                         "FusedBuffer_CopyToDevice returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombination_Sycl: "
                         "GetKernelParameters returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  SYCL_FOR(
    Q, nthreads_total, nthreads_per_block, item, GRID_STRIDE_XLOOP(item, i, N) {
      zdata[i] = cdata[0] * xdata[0][i];
      for (int j = 1; j < nvec; j++) { zdata[i] += cdata[j] * xdata[j][i]; }
    });

  return SUN_SUCCESS;
}

SUNErrCode nvScaleAddMulti_Sycl(int nvec, sunrealtype* c, N_Vector x,
                                N_Vector* Y, N_Vector* Z)
{
  const sunindextype N     = NVEC_SYCL_LENGTH(x);
  const sunrealtype* xdata = NVEC_SYCL_DDATAp(x);
  ::sycl::queue* Q         = NVEC_SYCL_QUEUE(x);
  size_t nthreads_total, nthreads_per_block;

  /* Shortcuts to the fused op workspace */
  sunrealtype* cdata  = NULL;
  sunrealtype** ydata = NULL;
  sunrealtype** zdata = NULL;

  /* Setup the fused op workspace */
  if (FusedBuffer_Init(x, nvec, 2 * nvec))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvScaleAddMulti_Sycl: FusedBuffer_Init returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyRealArray(x, c, nvec, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleAddMulti_Sycl: "
                         "FusedBuffer_CopyRealArray returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(x, Y, nvec, &ydata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleAddMulti_Sycl: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(x, Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleAddMulti_Sycl: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyToDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleAddMulti_Sycl: "
                         "FusedBuffer_CopyToDevice returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (GetKernelParameters(x, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvScaleAddMulti_Sycl: GetKernelParameters returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  SYCL_FOR(
    Q, nthreads_total, nthreads_per_block, item, GRID_STRIDE_XLOOP(item, i, N) {
      for (int j = 0; j < nvec; j++)
      {
        zdata[j][i] = cdata[j] * xdata[i] + ydata[j][i];
      }
    });

  return SUN_SUCCESS;
}

/* --------------------------------------------------------------------------
 * vector array operations
 * -------------------------------------------------------------------------- */

SUNErrCode nvLinearSumVectorArray_Sycl(int nvec, sunrealtype a, N_Vector* X,
                                       sunrealtype b, N_Vector* Y, N_Vector* Z)
{
  const sunindextype N = NVEC_SYCL_LENGTH(Z[0]);
  ::sycl::queue* Q     = NVEC_SYCL_QUEUE(Z[0]);
  size_t nthreads_total, nthreads_per_block;

  /* Shortcuts to the fused op workspace */
  sunrealtype** xdata = NULL;
  sunrealtype** ydata = NULL;
  sunrealtype** zdata = NULL;

  /* Setup the fused op workspace */
  if (FusedBuffer_Init(Z[0], 0, 3 * nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearSumVectorArray_Sycl: "
                         "FusedBuffer_Init returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearSumVectorArray_Sycl: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Y, nvec, &ydata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearSumVectorArray_Sycl: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearSumVectorArray_Sycl: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyToDevice(Z[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinaerSumVectorArray_Sycl: "
                         "FusedBuffer_CopyToDevice returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (GetKernelParameters(Z[0], SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinaerSumVectorArray_Sycl: "
                         "GetKernelParameters returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  SYCL_FOR(
    Q, nthreads_total, nthreads_per_block, item, GRID_STRIDE_XLOOP(item, i, N) {
      for (int j = 0; j < nvec; j++)
      {
        zdata[j][i] = a * xdata[j][i] + b * ydata[j][i];
      }
    });

  return SUN_SUCCESS;
}

SUNErrCode nvScaleVectorArray_Sycl(int nvec, sunrealtype* c, N_Vector* X,
                                   N_Vector* Z)
{
  const sunindextype N = NVEC_SYCL_LENGTH(Z[0]);
  ::sycl::queue* Q     = NVEC_SYCL_QUEUE(Z[0]);
  size_t nthreads_total, nthreads_per_block;

  /* Shortcuts to the fused op workspace arrays */
  sunrealtype* cdata  = NULL;
  sunrealtype** xdata = NULL;
  sunrealtype** zdata = NULL;

  /* Setup the fused op workspace */
  if (FusedBuffer_Init(Z[0], nvec, 2 * nvec))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvScaleVectorArray_Sycl: FusedBuffer_Init returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyRealArray(Z[0], c, nvec, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleVectorArray_Sycl: "
                         "FusedBuffer_CopyReadArray returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleVectorArray_Sycl: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleVectorArray_Sycl: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyToDevice(Z[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleVectorArray_Sycl: "
                         "FusedBuffer_CopyToDevice returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (GetKernelParameters(Z[0], SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleVectorArray_Sycl: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  SYCL_FOR(
    Q, nthreads_total, nthreads_per_block, item, GRID_STRIDE_XLOOP(item, i, N) {
      for (int j = 0; j < nvec; j++) { zdata[j][i] = cdata[j] * xdata[j][i]; }
    });

  return SUN_SUCCESS;
}

SUNErrCode nvConstVectorArray_Sycl(int nvec, sunrealtype c, N_Vector* Z)
{
  const sunindextype N = NVEC_SYCL_LENGTH(Z[0]);
  ::sycl::queue* Q     = NVEC_SYCL_QUEUE(Z[0]);
  size_t nthreads_total, nthreads_per_block;

  /* Shortcuts to the fused op workspace arrays */
  sunrealtype** zdata = NULL;

  /* Setup the fused op workspace */
  if (FusedBuffer_Init(Z[0], 0, nvec))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvConstVectorArray_Sycl: FusedBuffer_Init returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvConstVectorArray_Sycl: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyToDevice(Z[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvConstVectorArray_Sycl: "
                         "FusedBuffer_CopyToDevice returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (GetKernelParameters(Z[0], SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvConstVectorArray_Sycl: "
                         "GetKernelParameters returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  SYCL_FOR(
    Q, nthreads_total, nthreads_per_block, item, GRID_STRIDE_XLOOP(item, i, N) {
      for (int j = 0; j < nvec; j++) { zdata[j][i] = c; }
    });

  return SUN_SUCCESS;
}

SUNErrCode nvScaleAddMultiVectorArray_Sycl(int nvec, int nsum, sunrealtype* c,
                                           N_Vector* X, N_Vector** Y,
                                           N_Vector** Z)
{
  const sunindextype N = NVEC_SYCL_LENGTH(X[0]);
  ::sycl::queue* Q     = NVEC_SYCL_QUEUE(X[0]);
  size_t nthreads_total, nthreads_per_block;

  /* Shortcuts to the fused op workspace */
  sunrealtype* cdata  = NULL;
  sunrealtype** xdata = NULL;
  sunrealtype** ydata = NULL;
  sunrealtype** zdata = NULL;

  /* Setup the fused op workspace */
  if (FusedBuffer_Init(X[0], nsum, nvec + 2 * nvec * nsum))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMultiArray_Sycl: "
                         "FusedBuffer_Init returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyRealArray(X[0], c, nsum, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMultiArray_Sycl: "
                         "FusedBuffer_CopyRealArray returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(X[0], X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleAddMultiVectorArray_Sycl: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray2D(X[0], Y, nvec, nsum, &ydata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleAddMultiVectorArray_Sycl: "
                         "FusedBuffer_CopyPtrArray2D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray2D(X[0], Z, nvec, nsum, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleAddMultiVectorArray_Sycl: "
                         "FusedBuffer_CopyPtrArray2D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyToDevice(X[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleVectorArray_Sycl: "
                         "FusedBuffer_CopyToDevice returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (GetKernelParameters(X[0], SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleAddMultiVectorArray_Sycl: "
                         "GetKernelParameters returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  SYCL_FOR(
    Q, nthreads_total, nthreads_per_block, item, GRID_STRIDE_XLOOP(item, i, N) {
      for (int j = 0; j < nvec; j++)
      {
        for (int k = 0; k < nsum; k++)
        {
          zdata[j * nsum + k][i] = cdata[k] * xdata[j][i] +
                                   ydata[j * nsum + k][i];
        }
      }
    });

  return SUN_SUCCESS;
}

SUNErrCode nvLinearCombinationVectorArray_Sycl(int nvec, int nsum, sunrealtype* c,
                                               N_Vector** X, N_Vector* Z)
{
  const sunindextype N = NVEC_SYCL_LENGTH(Z[0]);
  ::sycl::queue* Q     = NVEC_SYCL_QUEUE(Z[0]);
  size_t nthreads_total, nthreads_per_block;

  /* Shortcuts to the fused op workspace arrays */
  sunrealtype* cdata  = NULL;
  sunrealtype** xdata = NULL;
  sunrealtype** zdata = NULL;

  /* Setup the fused op workspace */
  if (FusedBuffer_Init(Z[0], nsum, nvec + nvec * nsum))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombinationVectorArray_Sycl: "
                         "FusedBuffer_Init returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyRealArray(Z[0], c, nsum, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombinationVectorArray_Sycl: "
                         "FusedBuffer_CopyRealArray returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray2D(Z[0], X, nvec, nsum, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombinationVectorArray_Sycl: "
                         "FusedBuffer_CopyPtrArray2D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombinationVectorArray_Sycl: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyToDevice(Z[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombinationVectorArray_Sycl: "
                         "FusedBuffer_CopyToDevice returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (GetKernelParameters(Z[0], SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombinationVectorArray_Sycl: "
                         "GetKernelParameters returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  SYCL_FOR(
    Q, nthreads_total, nthreads_per_block, item, GRID_STRIDE_XLOOP(item, i, N) {
      for (int j = 0; j < nvec; j++)
      {
        zdata[j][i] = cdata[0] * xdata[j * nsum][i];
        for (int k = 1; k < nsum; k++)
        {
          zdata[j][i] += cdata[k] * xdata[j * nsum + k][i];
        }
      }
    });

  return SUN_SUCCESS;
}

/* --------------------------------------------------------------------------
 * OPTIONAL XBraid interface operations
 * -------------------------------------------------------------------------- */

SUNErrCode nvBufSize_Sycl(N_Vector x, sunindextype* size)
{
  if (x == NULL) { return SUN_ERR_GENERIC; }
  *size = (sunindextype)NVEC_SYCL_MEMSIZE(x);
  return SUN_SUCCESS;
}

SUNErrCode nvBufPack_Sycl(N_Vector x, void* buf)
{
  int copy_fail = 0;

  if (x == NULL || buf == NULL) { return SUN_ERR_GENERIC; }

  SUNMemory buf_mem = SUNMemoryHelper_Wrap(NVEC_SYCL_MEMHELP(x), buf,
                                           SUNMEMTYPE_HOST);
  if (buf_mem == NULL) { return SUN_ERR_GENERIC; }

  copy_fail = SUNMemoryHelper_Copy(NVEC_SYCL_MEMHELP(x), buf_mem,
                                   NVEC_SYCL_CONTENT(x)->device_data,
                                   NVEC_SYCL_MEMSIZE(x), NVEC_SYCL_QUEUE(x));

  /* synchronize with the host */
  NVEC_SYCL_QUEUE(x)->wait_and_throw();

  SUNMemoryHelper_Dealloc(NVEC_SYCL_MEMHELP(x), buf_mem, NVEC_SYCL_QUEUE(x));

  return (copy_fail ? SUN_ERR_GENERIC : SUN_SUCCESS);
}

SUNErrCode nvBufUnpack_Sycl(N_Vector x, void* buf)
{
  int copy_fail = 0;

  if (x == NULL || buf == NULL) { return SUN_ERR_GENERIC; }

  SUNMemory buf_mem = SUNMemoryHelper_Wrap(NVEC_SYCL_MEMHELP(x), buf,
                                           SUNMEMTYPE_HOST);
  if (buf_mem == NULL) { return SUN_ERR_GENERIC; }

  copy_fail = SUNMemoryHelper_Copy(NVEC_SYCL_MEMHELP(x),
                                   NVEC_SYCL_CONTENT(x)->device_data, buf_mem,
                                   NVEC_SYCL_MEMSIZE(x), NVEC_SYCL_QUEUE(x));

  /* synchronize with the host */
  NVEC_SYCL_QUEUE(x)->wait_and_throw();

  SUNMemoryHelper_Dealloc(NVEC_SYCL_MEMHELP(x), buf_mem, NVEC_SYCL_QUEUE(x));

  return (copy_fail ? SUN_ERR_GENERIC : SUN_SUCCESS);
}

/* --------------------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -------------------------------------------------------------------------- */

SUNErrCode N_VEnableFusedOps_Sycl(N_Vector v, sunbooleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) { return SUN_ERR_GENERIC; }

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }

  if (tf)
  {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = nvLinearCombination_Sycl;
    v->ops->nvscaleaddmulti     = nvScaleAddMulti_Sycl;
    v->ops->nvdotprodmulti      = NULL;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray     = nvLinearSumVectorArray_Sycl;
    v->ops->nvscalevectorarray         = nvScaleVectorArray_Sycl;
    v->ops->nvconstvectorarray         = nvConstVectorArray_Sycl;
    v->ops->nvwrmsnormvectorarray      = NULL;
    v->ops->nvwrmsnormmaskvectorarray  = NULL;
    v->ops->nvscaleaddmultivectorarray = nvScaleAddMultiVectorArray_Sycl;
    v->ops->nvlinearcombinationvectorarray = nvLinearCombinationVectorArray_Sycl;
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
  }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearCombination_Sycl(N_Vector v, sunbooleantype tf)
{
  if (v == NULL) { return SUN_ERR_GENERIC; }
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }
  v->ops->nvlinearcombination = tf ? nvLinearCombination_Sycl : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleAddMulti_Sycl(N_Vector v, sunbooleantype tf)
{
  if (v == NULL) { return SUN_ERR_GENERIC; }
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }
  v->ops->nvscaleaddmulti = tf ? nvScaleAddMulti_Sycl : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearSumVectorArray_Sycl(N_Vector v, sunbooleantype tf)
{
  if (v == NULL) { return SUN_ERR_GENERIC; }
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }
  v->ops->nvlinearsumvectorarray = tf ? nvLinearSumVectorArray_Sycl : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleVectorArray_Sycl(N_Vector v, sunbooleantype tf)
{
  if (v == NULL) { return SUN_ERR_GENERIC; }
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }
  v->ops->nvscalevectorarray = tf ? nvScaleVectorArray_Sycl : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableConstVectorArray_Sycl(N_Vector v, sunbooleantype tf)
{
  if (v == NULL) { return SUN_ERR_GENERIC; }
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }
  v->ops->nvconstvectorarray = tf ? nvConstVectorArray_Sycl : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleAddMultiVectorArray_Sycl(N_Vector v, sunbooleantype tf)
{
  if (v == NULL) { return SUN_ERR_GENERIC; }
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }
  v->ops->nvscaleaddmultivectorarray = tf ? nvScaleAddMultiVectorArray_Sycl
                                          : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearCombinationVectorArray_Sycl(N_Vector v,
                                                      sunbooleantype tf)
{
  if (v == NULL) { return SUN_ERR_GENERIC; }
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }
  v->ops->nvlinearcombinationvectorarray =
    tf ? nvLinearCombinationVectorArray_Sycl : NULL;
  return SUN_SUCCESS;
}

/* --------------------------------------------------------------------------
 * Private utility functions
 * -------------------------------------------------------------------------- */

static int AllocateData(N_Vector v)
{
  int alloc_fail                  = 0;
  N_VectorContent_Sycl vc         = NVEC_SYCL_CONTENT(v);
  N_PrivateVectorContent_Sycl vcp = NVEC_SYCL_PRIVATE(v);

  if (nvGetLength_Sycl(v) == 0) { return SUN_SUCCESS; }

  if (vcp->use_managed_mem)
  {
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_SYCL_MEMHELP(v), &(vc->device_data),
                                       NVEC_SYCL_MEMSIZE(v), SUNMEMTYPE_UVM,
                                       NVEC_SYCL_QUEUE(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc "
                           "failed for SUNMEMTYPE_UVM\n");
    }
    vc->host_data = SUNMemoryHelper_Alias(NVEC_SYCL_MEMHELP(v), vc->device_data);
  }
  else
  {
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_SYCL_MEMHELP(v), &(vc->host_data),
                                       NVEC_SYCL_MEMSIZE(v), SUNMEMTYPE_HOST,
                                       NVEC_SYCL_QUEUE(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc "
                           "failed to alloc SUNMEMTYPE_HOST\n");
    }

    alloc_fail = SUNMemoryHelper_Alloc(NVEC_SYCL_MEMHELP(v), &(vc->device_data),
                                       NVEC_SYCL_MEMSIZE(v), SUNMEMTYPE_DEVICE,
                                       NVEC_SYCL_QUEUE(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc "
                           "failed to alloc SUNMEMTYPE_DEVICE\n");
    }
  }

  return (alloc_fail ? SUN_ERR_GENERIC : SUN_SUCCESS);
}

/* Allocate and initializes the internal memory used for reductions */
static int InitializeReductionBuffer(N_Vector v, const sunrealtype value, size_t n)
{
  int alloc_fail           = 0;
  int copy_fail            = 0;
  sunbooleantype alloc_mem = SUNFALSE;
  size_t bytes             = n * sizeof(sunrealtype);

  /* Get the vector private memory structure */
  N_PrivateVectorContent_Sycl vcp = NVEC_SYCL_PRIVATE(v);

  /* Wrap the initial value as SUNMemory object */
  SUNMemory value_mem = SUNMemoryHelper_Wrap(NVEC_SYCL_MEMHELP(v),
                                             (void*)&value, SUNMEMTYPE_HOST);

  /* check if the existing reduction memory is not large enough */
  if (vcp->reduce_buffer_bytes < bytes)
  {
    FreeReductionBuffer(v);
    alloc_mem = SUNTRUE;
  }

  if (alloc_mem)
  {
    /* allocate pinned memory on the host */
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_SYCL_MEMHELP(v),
                                       &(vcp->reduce_buffer_host), bytes,
                                       SUNMEMTYPE_PINNED, NVEC_SYCL_QUEUE(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT(
        "WARNING in InitializeReductionBuffer: SUNMemoryHelper_Alloc failed to "
        "alloc SUNMEMTYPE_PINNED, using SUNMEMTYPE_HOST instead\n");

      /* if pinned alloc failed, allocate plain host memory */
      alloc_fail = SUNMemoryHelper_Alloc(NVEC_SYCL_MEMHELP(v),
                                         &(vcp->reduce_buffer_host), bytes,
                                         SUNMEMTYPE_HOST, NVEC_SYCL_QUEUE(v));
      if (alloc_fail)
      {
        SUNDIALS_DEBUG_PRINT(
          "ERROR in InitializeReductionBuffer: SUNMemoryHelper_Alloc failed to "
          "alloc SUNMEMTYPE_HOST\n");
      }
    }

    /* allocate device memory */
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_SYCL_MEMHELP(v),
                                       &(vcp->reduce_buffer_dev), bytes,
                                       SUNMEMTYPE_DEVICE, NVEC_SYCL_QUEUE(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in InitializeReductionBuffer: SUNMemoryHelper_Alloc failed to "
        "alloc SUNMEMTYPE_DEVICE\n");
    }
  }

  if (!alloc_fail)
  {
    /* store the size of the reduction memory */
    vcp->reduce_buffer_bytes = bytes;

    /* initialize the memory with the value */
    copy_fail = SUNMemoryHelper_CopyAsync(NVEC_SYCL_MEMHELP(v),
                                          vcp->reduce_buffer_dev, value_mem,
                                          bytes, (void*)NVEC_SYCL_QUEUE(v));

    /* wait for copy to finish (possible bug in minimum reduction object) */
    NVEC_SYCL_QUEUE(v)->wait_and_throw();

    if (copy_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in InitializeReductionBuffer: "
                           "SUNMemoryHelper_CopyAsync failed\n");
    }
  }

  /* deallocate the wrapper */
  SUNMemoryHelper_Dealloc(NVEC_SYCL_MEMHELP(v), value_mem, NVEC_SYCL_QUEUE(v));

  return ((alloc_fail || copy_fail) ? SUN_ERR_GENERIC : SUN_SUCCESS);
}

/* Free the reduction memory */
static void FreeReductionBuffer(N_Vector v)
{
  N_PrivateVectorContent_Sycl vcp = NVEC_SYCL_PRIVATE(v);

  if (vcp == NULL) { return; }

  /* free device mem */
  if (vcp->reduce_buffer_dev != NULL)
  {
    SUNMemoryHelper_Dealloc(NVEC_SYCL_MEMHELP(v), vcp->reduce_buffer_dev,
                            NVEC_SYCL_QUEUE(v));
  }
  vcp->reduce_buffer_dev = NULL;

  /* free host mem */
  if (vcp->reduce_buffer_host != NULL)
  {
    SUNMemoryHelper_Dealloc(NVEC_SYCL_MEMHELP(v), vcp->reduce_buffer_host,
                            NVEC_SYCL_QUEUE(v));
  }
  vcp->reduce_buffer_host = NULL;

  /* reset allocated memory size */
  vcp->reduce_buffer_bytes = 0;
}

/* Copy the reduction memory from the device to the host. */
static int CopyReductionBufferFromDevice(N_Vector v, size_t n)
{
  int copy_fail;

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_SYCL_MEMHELP(v),
                                        NVEC_SYCL_PRIVATE(v)->reduce_buffer_host,
                                        NVEC_SYCL_PRIVATE(v)->reduce_buffer_dev,
                                        n * sizeof(sunrealtype),
                                        (void*)NVEC_SYCL_QUEUE(v));

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in CopyReductionBufferFromDevice: "
                         "SUNMemoryHelper_CopyAsync returned nonzero\n");
  }

  /* synchronize with respect to the host */
  NVEC_SYCL_QUEUE(v)->wait_and_throw();

  return (copy_fail ? SUN_ERR_GENERIC : SUN_SUCCESS);
}

static int FusedBuffer_Init(N_Vector v, int nreal, int nptr)
{
  int alloc_fail           = 0;
  sunbooleantype alloc_mem = SUNFALSE;
  size_t bytes = nreal * sizeof(sunrealtype) + nptr * sizeof(sunrealtype*);

  /* Get the vector private memory structure */
  N_PrivateVectorContent_Sycl vcp = NVEC_SYCL_PRIVATE(v);

  /* Check if the existing memory is not large enough */
  if (vcp->fused_buffer_bytes < bytes)
  {
    FusedBuffer_Free(v);
    alloc_mem = SUNTRUE;
  }

  if (alloc_mem)
  {
    /* allocate pinned memory on the host */
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_SYCL_MEMHELP(v),
                                       &(vcp->fused_buffer_host), bytes,
                                       SUNMEMTYPE_PINNED, NVEC_SYCL_QUEUE(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT(
        "WARNING in FusedBuffer_Init: SUNMemoryHelper_Alloc failed to alloc "
        "SUNMEMTYPE_PINNED, using SUNMEMTYPE_HOST instead\n");

      /* if pinned alloc failed, allocate plain host memory */
      alloc_fail = SUNMemoryHelper_Alloc(NVEC_SYCL_MEMHELP(v),
                                         &(vcp->fused_buffer_host), bytes,
                                         SUNMEMTYPE_HOST, NVEC_SYCL_QUEUE(v));
      if (alloc_fail)
      {
        SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_Init: SUNMemoryHelper_Alloc "
                             "failed to alloc SUNMEMTYPE_HOST\n");
        return SUN_ERR_GENERIC;
      }
    }

    /* allocate device memory */
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_SYCL_MEMHELP(v),
                                       &(vcp->fused_buffer_dev), bytes,
                                       SUNMEMTYPE_DEVICE, NVEC_SYCL_QUEUE(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_Init: SUNMemoryHelper_Alloc "
                           "failed to alloc SUNMEMTYPE_DEVICE\n");
      return SUN_ERR_GENERIC;
    }

    /* Store the size of the fused op buffer */
    vcp->fused_buffer_bytes = bytes;
  }

  /* Reset the buffer offset */
  vcp->fused_buffer_offset = 0;

  return SUN_SUCCESS;
}

static int FusedBuffer_CopyRealArray(N_Vector v, sunrealtype* rdata, int nval,
                                     sunrealtype** shortcut)
{
  /* Get the vector private memory structure */
  N_PrivateVectorContent_Sycl vcp = NVEC_SYCL_PRIVATE(v);

  /* Check buffer space and fill the host buffer */
  if (vcp->fused_buffer_offset >= vcp->fused_buffer_bytes)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_CopyRealArray: Buffer offset is "
                         "exceedes the buffer size\n");
    return SUN_ERR_GENERIC;
  }

  sunrealtype* h_buffer = reinterpret_cast<sunrealtype*>(
    (char*)(vcp->fused_buffer_host->ptr) + vcp->fused_buffer_offset);

  for (int j = 0; j < nval; j++) { h_buffer[j] = rdata[j]; }

  /* Set shortcut to the device buffer and update offset*/
  *shortcut = reinterpret_cast<sunrealtype*>(
    (char*)(vcp->fused_buffer_dev->ptr) + vcp->fused_buffer_offset);

  vcp->fused_buffer_offset += nval * sizeof(sunrealtype);

  return SUN_SUCCESS;
}

static int FusedBuffer_CopyPtrArray1D(N_Vector v, N_Vector* X, int nvec,
                                      sunrealtype*** shortcut)
{
  /* Get the vector private memory structure */
  N_PrivateVectorContent_Sycl vcp = NVEC_SYCL_PRIVATE(v);

  /* Check buffer space and fill the host buffer */
  if (vcp->fused_buffer_offset >= vcp->fused_buffer_bytes)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_CopyPtrArray1D: Buffer offset "
                         "is exceedes the buffer size\n");
    return SUN_ERR_GENERIC;
  }

  sunrealtype** h_buffer = reinterpret_cast<sunrealtype**>(
    (char*)(vcp->fused_buffer_host->ptr) + vcp->fused_buffer_offset);

  for (int j = 0; j < nvec; j++) { h_buffer[j] = NVEC_SYCL_DDATAp(X[j]); }

  /* Set shortcut to the device buffer and update offset*/
  *shortcut = reinterpret_cast<sunrealtype**>(
    (char*)(vcp->fused_buffer_dev->ptr) + vcp->fused_buffer_offset);

  vcp->fused_buffer_offset += nvec * sizeof(sunrealtype*);

  return SUN_SUCCESS;
}

static int FusedBuffer_CopyPtrArray2D(N_Vector v, N_Vector** X, int nvec,
                                      int nsum, sunrealtype*** shortcut)
{
  /* Get the vector private memory structure */
  N_PrivateVectorContent_Sycl vcp = NVEC_SYCL_PRIVATE(v);

  /* Check buffer space and fill the host buffer */
  if (vcp->fused_buffer_offset >= vcp->fused_buffer_bytes)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_CopyPtrArray2D: Buffer offset "
                         "is exceedes the buffer size\n");
    return SUN_ERR_GENERIC;
  }

  sunrealtype** h_buffer = reinterpret_cast<sunrealtype**>(
    (char*)(vcp->fused_buffer_host->ptr) + vcp->fused_buffer_offset);

  for (int j = 0; j < nvec; j++)
  {
    for (int k = 0; k < nsum; k++)
    {
      h_buffer[j * nsum + k] = NVEC_SYCL_DDATAp(X[k][j]);
    }
  }

  /* Set shortcut to the device buffer and update offset*/
  *shortcut = reinterpret_cast<sunrealtype**>(
    (char*)(vcp->fused_buffer_dev->ptr) + vcp->fused_buffer_offset);

  /* Update the offset */
  vcp->fused_buffer_offset += nvec * nsum * sizeof(sunrealtype*);

  return SUN_SUCCESS;
}

static int FusedBuffer_CopyToDevice(N_Vector v)
{
  /* Get the vector private memory structure */
  N_PrivateVectorContent_Sycl vcp = NVEC_SYCL_PRIVATE(v);

  /* Copy the fused buffer to the device */
  int copy_fail =
    SUNMemoryHelper_CopyAsync(NVEC_SYCL_MEMHELP(v), vcp->fused_buffer_dev,
                              vcp->fused_buffer_host, vcp->fused_buffer_offset,
                              (void*)NVEC_SYCL_QUEUE(v));
  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in FusedBuffer_CopyToDevice: SUNMemoryHelper_CopyAsync failed\n");
    return SUN_ERR_GENERIC;
  }

  return SUN_SUCCESS;
}

static int FusedBuffer_Free(N_Vector v)
{
  N_PrivateVectorContent_Sycl vcp = NVEC_SYCL_PRIVATE(v);

  if (vcp == NULL) { return SUN_SUCCESS; }

  if (vcp->fused_buffer_host)
  {
    SUNMemoryHelper_Dealloc(NVEC_SYCL_MEMHELP(v), vcp->fused_buffer_host,
                            NVEC_SYCL_QUEUE(v));
    vcp->fused_buffer_host = NULL;
  }

  if (vcp->fused_buffer_dev)
  {
    SUNMemoryHelper_Dealloc(NVEC_SYCL_MEMHELP(v), vcp->fused_buffer_dev,
                            NVEC_SYCL_QUEUE(v));
    vcp->fused_buffer_dev = NULL;
  }

  vcp->fused_buffer_bytes  = 0;
  vcp->fused_buffer_offset = 0;

  return SUN_SUCCESS;
}

/* Get the kernel launch parameters based on the kernel type (reduction or not),
 * using the appropriate kernel execution policy. */
static int GetKernelParameters(N_Vector v, sunbooleantype reduction,
                               size_t& nthreads_total, size_t& nthreads_per_block)
{
  /* Get the execution policy */
  SUNSyclExecPolicy* exec_policy = NULL;
  exec_policy = reduction ? NVEC_SYCL_CONTENT(v)->reduce_exec_policy
                          : NVEC_SYCL_CONTENT(v)->stream_exec_policy;

  if (exec_policy == NULL)
  {
    SUNDIALS_DEBUG_ERROR("The execution policy is NULL\n");
    return SUN_ERR_GENERIC;
  }

  /* Get the number of threads per block and total number threads */
  nthreads_per_block = exec_policy->blockSize();
  nthreads_total     = nthreads_per_block *
                   exec_policy->gridSize(NVEC_SYCL_LENGTH(v));

  if (nthreads_per_block == 0)
  {
    SUNDIALS_DEBUG_ERROR("The number of threads per block must be > 0\n");
    return SUN_ERR_GENERIC;
  }

  if (nthreads_total == 0)
  {
    SUNDIALS_DEBUG_ERROR("the total number of threads must be > 0\n");
    return SUN_ERR_GENERIC;
  }

  return SUN_SUCCESS;
}

} /* extern "C" */

/* Deprecated concrete operation wrappers */

extern "C" {

void N_VAbs_Sycl(N_Vector x, N_Vector z) { nvAbs_Sycl(x, z); }

void N_VAddConst_Sycl(N_Vector x, sunrealtype b, N_Vector z)
{
  nvAddConst_Sycl(x, b, z);
}

SUNErrCode N_VBufPack_Sycl(N_Vector x, void* buf)
{
  return nvBufPack_Sycl(x, buf);
}

SUNErrCode N_VBufSize_Sycl(N_Vector x, sunindextype* size)
{
  return nvBufSize_Sycl(x, size);
}

SUNErrCode N_VBufUnpack_Sycl(N_Vector x, void* buf)
{
  return nvBufUnpack_Sycl(x, buf);
}

N_Vector N_VCloneEmpty_Sycl(N_Vector w) { return nvCloneEmpty_Sycl(w); }

N_Vector N_VClone_Sycl(N_Vector w) { return nvClone_Sycl(w); }

void N_VCompare_Sycl(sunrealtype c, N_Vector x, N_Vector z)
{
  nvCompare_Sycl(c, x, z);
}

SUNErrCode N_VConstVectorArray_Sycl(int nvec, sunrealtype c, N_Vector* Z)
{
  return nvConstVectorArray_Sycl(nvec, c, Z);
}

void N_VConst_Sycl(sunrealtype c, N_Vector z) { nvConst_Sycl(c, z); }

sunbooleantype N_VConstrMask_Sycl(N_Vector c, N_Vector x, N_Vector m)
{
  return nvConstrMask_Sycl(c, x, m);
}

void N_VDestroy_Sycl(N_Vector v) { nvDestroy_Sycl(v); }

void N_VDiv_Sycl(N_Vector x, N_Vector y, N_Vector z) { nvDiv_Sycl(x, y, z); }

sunrealtype N_VDotProd_Sycl(N_Vector x, N_Vector y)
{
  return nvDotProd_Sycl(x, y);
}

sunbooleantype N_VInvTest_Sycl(N_Vector x, N_Vector z)
{
  return nvInvTest_Sycl(x, z);
}

void N_VInv_Sycl(N_Vector x, N_Vector z) { nvInv_Sycl(x, z); }

sunrealtype N_VL1Norm_Sycl(N_Vector x) { return nvL1Norm_Sycl(x); }

SUNErrCode N_VLinearCombinationVectorArray_Sycl(int nvec, int nsum,
                                                sunrealtype* c, N_Vector** X,
                                                N_Vector* Z)
{
  return nvLinearCombinationVectorArray_Sycl(nvec, nsum, c, X, Z);
}

SUNErrCode N_VLinearCombination_Sycl(int nvec, sunrealtype* c, N_Vector* X,
                                     N_Vector Z)
{
  return nvLinearCombination_Sycl(nvec, c, X, Z);
}

SUNErrCode N_VLinearSumVectorArray_Sycl(int nvec, sunrealtype a, N_Vector* X,
                                        sunrealtype b, N_Vector* Y, N_Vector* Z)
{
  return nvLinearSumVectorArray_Sycl(nvec, a, X, b, Y, Z);
}

void N_VLinearSum_Sycl(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                       N_Vector z)
{
  nvLinearSum_Sycl(a, x, b, y, z);
}

sunrealtype N_VMaxNorm_Sycl(N_Vector x) { return nvMaxNorm_Sycl(x); }

sunrealtype N_VMinQuotient_Sycl(N_Vector num, N_Vector denom)
{
  return nvMinQuotient_Sycl(num, denom);
}

sunrealtype N_VMin_Sycl(N_Vector x) { return nvMin_Sycl(x); }

void N_VPrintFile_Sycl(N_Vector v, FILE* outfile)
{
  nvPrintFile_Sycl(v, outfile);
}

void N_VPrint_Sycl(N_Vector v) { nvPrint_Sycl(v); }

void N_VProd_Sycl(N_Vector x, N_Vector y, N_Vector z) { nvProd_Sycl(x, y, z); }

SUNErrCode N_VScaleAddMultiVectorArray_Sycl(int nvec, int nsum, sunrealtype* a,
                                            N_Vector* X, N_Vector** Y,
                                            N_Vector** Z)
{
  return nvScaleAddMultiVectorArray_Sycl(nvec, nsum, a, X, Y, Z);
}

SUNErrCode N_VScaleAddMulti_Sycl(int nvec, sunrealtype* c, N_Vector X,
                                 N_Vector* Y, N_Vector* Z)
{
  return nvScaleAddMulti_Sycl(nvec, c, X, Y, Z);
}

SUNErrCode N_VScaleVectorArray_Sycl(int nvec, sunrealtype* c, N_Vector* X,
                                    N_Vector* Z)
{
  return nvScaleVectorArray_Sycl(nvec, c, X, Z);
}

void N_VScale_Sycl(sunrealtype c, N_Vector x, N_Vector z)
{
  nvScale_Sycl(c, x, z);
}

void N_VSetDeviceArrayPointer_Sycl(sunrealtype* d_vdata_1d, N_Vector v)
{
  nvSetDeviceArrayPointer_Sycl(d_vdata_1d, v);
}

void N_VSetHostArrayPointer_Sycl(sunrealtype* h_vdata_1d, N_Vector v)
{
  nvSetHostArrayPointer_Sycl(h_vdata_1d, v);
}

sunrealtype N_VWL2Norm_Sycl(N_Vector x, N_Vector w)
{
  return nvWL2Norm_Sycl(x, w);
}

sunrealtype N_VWSqrSumLocal_Sycl(N_Vector x, N_Vector w)
{
  return nvWSqrSumLocal_Sycl(x, w);
}

sunrealtype N_VWSqrSumMaskLocal_Sycl(N_Vector x, N_Vector w, N_Vector id)
{
  return nvWSqrSumMaskLocal_Sycl(x, w, id);
}

sunrealtype N_VWrmsNormMask_Sycl(N_Vector x, N_Vector w, N_Vector id)
{
  return nvWrmsNormMask_Sycl(x, w, id);
}

sunrealtype N_VWrmsNorm_Sycl(N_Vector x, N_Vector w)
{
  return nvWrmsNorm_Sycl(x, w);
}

} // extern "C"
