/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles, and Cody J. Balos @ LLNL
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
 * This is the implementation file for a CUDA implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <nvector/nvector_cuda_deprecated.h>
#include <sundials/priv/sundials_context_impl.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>
#include "VectorArrayKernels.cuh"
#include "VectorKernels.cuh"
#include "sundials/sundials_errors.h"
#include "sundials_cuda.h"
#include "sundials_debug.h"

#define ZERO SUN_RCONST(0.0)
#define HALF SUN_RCONST(0.5)

using namespace sundials;
using namespace sundials::cuda;
using namespace sundials::cuda::impl;
/*
 * Private function definitions
 */

// Functions attached to the N_Vector
static void nvAbs_Cuda(N_Vector x, N_Vector z);
static void nvAddConst_Cuda(N_Vector x, sunrealtype b, N_Vector z);
static SUNErrCode nvBufPack_Cuda(N_Vector x, void* buf);
static SUNErrCode nvBufSize_Cuda(N_Vector x, sunindextype* size);
static SUNErrCode nvBufUnpack_Cuda(N_Vector x, void* buf);
static N_Vector nvCloneEmpty_Cuda(N_Vector w);
static N_Vector nvClone_Cuda(N_Vector w);
static void nvCompare_Cuda(sunrealtype c, N_Vector x, N_Vector z);
static SUNErrCode nvConstVectorArray_Cuda(int nvec, sunrealtype c, N_Vector* Z);
static void nvConst_Cuda(sunrealtype c, N_Vector z);
static sunbooleantype nvConstrMask_Cuda(N_Vector c, N_Vector x, N_Vector m);
static void nvDestroy_Cuda(N_Vector v);
static void nvDiv_Cuda(N_Vector x, N_Vector y, N_Vector z);
static SUNErrCode nvDotProdMulti_Cuda(int nvec, N_Vector x, N_Vector* Y,
                                      sunrealtype* dotprods);
static sunrealtype nvDotProd_Cuda(N_Vector x, N_Vector y);

static inline sunrealtype* nvGetDeviceArrayPointer_Cuda(N_Vector x)
{
  N_VectorContent_Cuda content = (N_VectorContent_Cuda)x->content;
  return (content->device_data == NULL ? NULL
                                       : (sunrealtype*)content->device_data->ptr);
}

static inline sunrealtype* nvGetHostArrayPointer_Cuda(N_Vector x)
{
  N_VectorContent_Cuda content = (N_VectorContent_Cuda)x->content;
  return (content->host_data == NULL ? NULL
                                     : (sunrealtype*)content->host_data->ptr);
}

static inline sunindextype nvGetLength_Cuda(N_Vector x)
{
  N_VectorContent_Cuda content = (N_VectorContent_Cuda)x->content;
  return content->length;
}

static inline N_Vector_ID nvGetVectorID_Cuda(N_Vector /*v*/)
{
  return SUNDIALS_NVEC_CUDA;
}

static sunbooleantype nvInvTest_Cuda(N_Vector x, N_Vector z);
static void nvInv_Cuda(N_Vector x, N_Vector z);
static sunrealtype nvL1Norm_Cuda(N_Vector x);
static SUNErrCode nvLinearCombinationVectorArray_Cuda(int nvec, int nsum,
                                                      sunrealtype* c,
                                                      N_Vector** X, N_Vector* Z);
static SUNErrCode nvLinearCombination_Cuda(int nvec, sunrealtype* c,
                                           N_Vector* X, N_Vector Z);
static SUNErrCode nvLinearSumVectorArray_Cuda(int nvec, sunrealtype a,
                                              N_Vector* X, sunrealtype b,
                                              N_Vector* Y, N_Vector* Z);
static void nvLinearSum_Cuda(sunrealtype a, N_Vector x, sunrealtype b,
                             N_Vector y, N_Vector z);
static sunrealtype nvMaxNorm_Cuda(N_Vector x);
static sunrealtype nvMinQuotient_Cuda(N_Vector num, N_Vector denom);
static sunrealtype nvMin_Cuda(N_Vector x);
static void nvPrintFile_Cuda(N_Vector v, FILE* outfile);
static void nvPrint_Cuda(N_Vector v);
static void nvProd_Cuda(N_Vector x, N_Vector y, N_Vector z);
static SUNErrCode nvScaleAddMultiVectorArray_Cuda(int nvec, int nsum,
                                                  sunrealtype* a, N_Vector* X,
                                                  N_Vector** Y, N_Vector** Z);
static SUNErrCode nvScaleAddMulti_Cuda(int nvec, sunrealtype* c, N_Vector X,
                                       N_Vector* Y, N_Vector* Z);
static SUNErrCode nvScaleVectorArray_Cuda(int nvec, sunrealtype* c, N_Vector* X,
                                          N_Vector* Z);
static void nvScale_Cuda(sunrealtype c, N_Vector x, N_Vector z);
static void nvSetDeviceArrayPointer_Cuda(sunrealtype* d_vdata_1d, N_Vector v);
static void nvSetHostArrayPointer_Cuda(sunrealtype* h_vdata_1d, N_Vector v);
static sunrealtype nvWL2Norm_Cuda(N_Vector x, N_Vector w);
static sunrealtype nvWSqrSumLocal_Cuda(N_Vector x, N_Vector w);
static sunrealtype nvWSqrSumMaskLocal_Cuda(N_Vector x, N_Vector w, N_Vector id);
static SUNErrCode nvWrmsNormMaskVectorArray_Cuda(int nvec, N_Vector* X,
                                                 N_Vector* W, N_Vector id,
                                                 sunrealtype* nrm);
static sunrealtype nvWrmsNormMask_Cuda(N_Vector x, N_Vector w, N_Vector id);
static SUNErrCode nvWrmsNormVectorArray_Cuda(int nvec, N_Vector* X, N_Vector* W,
                                             sunrealtype* nrm);
static sunrealtype nvWrmsNorm_Cuda(N_Vector x, N_Vector w);

// Allocate vector data
static int AllocateData(N_Vector v);

// Reduction buffer functions
static int InitializeDeviceCounter(N_Vector v);
static int FreeDeviceCounter(N_Vector v);
static int InitializeReductionBuffer(N_Vector v, sunrealtype value, size_t n = 1);
static void FreeReductionBuffer(N_Vector v);
static int CopyReductionBufferFromDevice(N_Vector v, size_t n = 1);

// Fused operation buffer functions
static int FusedBuffer_Init(N_Vector v, int nreal, int nptr);
static int FusedBuffer_CopyRealArray(N_Vector v, sunrealtype* r_data, int nval,
                                     sunrealtype** shortcut);
static int FusedBuffer_CopyPtrArray1D(N_Vector v, N_Vector* X, int nvec,
                                      sunrealtype*** shortcut);
static int FusedBuffer_CopyPtrArray2D(N_Vector v, N_Vector** X, int nvec,
                                      int nsum, sunrealtype*** shortcut);
static int FusedBuffer_CopyToDevice(N_Vector v);
static int FusedBuffer_Free(N_Vector v);

// Kernel launch parameters
static int GetKernelParameters(N_Vector v, sunbooleantype reduction,
                               size_t& grid, size_t& block, size_t& shMemSize,
                               cudaStream_t& stream, size_t n = 0);
static int GetKernelParameters(N_Vector v, sunbooleantype reduction,
                               size_t& grid, size_t& block, size_t& shMemSize,
                               cudaStream_t& stream, bool& atomic, size_t n = 0);
static void PostKernelLaunch();

/*
 * Macro definitions
 */

// Macros to access vector content
#define NVEC_CUDA_CONTENT(x) ((N_VectorContent_Cuda)(x->content))
#define NVEC_CUDA_MEMSIZE(x) \
  (NVEC_CUDA_CONTENT(x)->length * sizeof(sunrealtype))
#define NVEC_CUDA_MEMHELP(x) (NVEC_CUDA_CONTENT(x)->mem_helper)
#define NVEC_CUDA_HDATAp(x)  ((sunrealtype*)NVEC_CUDA_CONTENT(x)->host_data->ptr)
#define NVEC_CUDA_DDATAp(x) \
  ((sunrealtype*)NVEC_CUDA_CONTENT(x)->device_data->ptr)
#define NVEC_CUDA_STREAM(x) (NVEC_CUDA_CONTENT(x)->stream_exec_policy->stream())

// Macros to access vector private content
#define NVEC_CUDA_PRIVATE(x) \
  ((N_PrivateVectorContent_Cuda)(NVEC_CUDA_CONTENT(x)->priv))
#define NVEC_CUDA_HBUFFERp(x) \
  ((sunrealtype*)NVEC_CUDA_PRIVATE(x)->reduce_buffer_host->ptr)
#define NVEC_CUDA_DBUFFERp(x) \
  ((sunrealtype*)NVEC_CUDA_PRIVATE(x)->reduce_buffer_dev->ptr)
#define NVEC_CUDA_DCOUNTERp(x) \
  ((unsigned int*)NVEC_CUDA_PRIVATE(x)->device_counter->ptr)

/*
 * Private structure definition
 */

struct _N_PrivateVectorContent_Cuda
{
  sunbooleantype use_managed_mem; /* do data pointers use managed memory */

  // reduction workspace
  SUNMemory device_counter; // device memory for a counter (used in LDS reductions)
  SUNMemory reduce_buffer_dev;  // device memory for reductions
  SUNMemory reduce_buffer_host; // host memory for reductions
  size_t reduce_buffer_bytes;   // current size of reduction buffers

  // fused op workspace
  SUNMemory fused_buffer_dev;  // device memory for fused ops
  SUNMemory fused_buffer_host; // host memory for fused ops
  size_t fused_buffer_bytes;   // current size of the buffers
  size_t fused_buffer_offset;  // current offset into the buffer
};

typedef struct _N_PrivateVectorContent_Cuda* N_PrivateVectorContent_Cuda;

/* Default policies to clone */
ThreadDirectExecPolicy DEFAULT_STREAMING_EXECPOLICY(256);
BlockReduceAtomicExecPolicy DEFAULT_REDUCTION_EXECPOLICY(256);

extern "C" {

N_Vector N_VNewEmpty_Cuda(SUNContext sunctx)
{
  N_Vector v;

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty(sunctx);
  if (v == NULL) { return (NULL); }

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid           = nvGetVectorID_Cuda;
  v->ops->nvclone                 = nvClone_Cuda;
  v->ops->nvcloneempty            = nvCloneEmpty_Cuda;
  v->ops->nvdestroy               = nvDestroy_Cuda;
  v->ops->nvgetlength             = nvGetLength_Cuda;
  v->ops->nvgetarraypointer       = nvGetHostArrayPointer_Cuda;
  v->ops->nvgetdevicearraypointer = nvGetDeviceArrayPointer_Cuda;
  v->ops->nvsetarraypointer       = nvSetHostArrayPointer_Cuda;
  v->ops->nvsetdevicearraypointer = nvSetDeviceArrayPointer_Cuda;

  /* standard vector operations */
  v->ops->nvlinearsum    = nvLinearSum_Cuda;
  v->ops->nvconst        = nvConst_Cuda;
  v->ops->nvprod         = nvProd_Cuda;
  v->ops->nvdiv          = nvDiv_Cuda;
  v->ops->nvscale        = nvScale_Cuda;
  v->ops->nvabs          = nvAbs_Cuda;
  v->ops->nvinv          = nvInv_Cuda;
  v->ops->nvaddconst     = nvAddConst_Cuda;
  v->ops->nvdotprod      = nvDotProd_Cuda;
  v->ops->nvmaxnorm      = nvMaxNorm_Cuda;
  v->ops->nvmin          = nvMin_Cuda;
  v->ops->nvl1norm       = nvL1Norm_Cuda;
  v->ops->nvinvtest      = nvInvTest_Cuda;
  v->ops->nvconstrmask   = nvConstrMask_Cuda;
  v->ops->nvminquotient  = nvMinQuotient_Cuda;
  v->ops->nvwrmsnormmask = nvWrmsNormMask_Cuda;
  v->ops->nvwrmsnorm     = nvWrmsNorm_Cuda;
  v->ops->nvwl2norm      = nvWL2Norm_Cuda;
  v->ops->nvcompare      = nvCompare_Cuda;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvdotprodlocal     = nvDotProd_Cuda;
  v->ops->nvmaxnormlocal     = nvMaxNorm_Cuda;
  v->ops->nvminlocal         = nvMin_Cuda;
  v->ops->nvl1normlocal      = nvL1Norm_Cuda;
  v->ops->nvinvtestlocal     = nvInvTest_Cuda;
  v->ops->nvconstrmasklocal  = nvConstrMask_Cuda;
  v->ops->nvminquotientlocal = nvMinQuotient_Cuda;
  v->ops->nvwsqrsumlocal     = nvWSqrSumLocal_Cuda;
  v->ops->nvwsqrsummasklocal = nvWSqrSumMaskLocal_Cuda;

  /* single buffer reduction operations */
  v->ops->nvdotprodmultilocal = nvDotProdMulti_Cuda;

  /* XBraid interface operations */
  v->ops->nvbufsize   = nvBufSize_Cuda;
  v->ops->nvbufpack   = nvBufPack_Cuda;
  v->ops->nvbufunpack = nvBufUnpack_Cuda;

  /* print operation for debugging */
  v->ops->nvprint     = nvPrint_Cuda;
  v->ops->nvprintfile = nvPrintFile_Cuda;

  /* Create content */

  v->content = (N_VectorContent_Cuda)malloc(sizeof(_N_VectorContent_Cuda));
  if (v->content == NULL)
  {
    N_VDestroy(v);
    return (NULL);
  }

  NVEC_CUDA_CONTENT(v)->priv = malloc(sizeof(_N_PrivateVectorContent_Cuda));
  if (NVEC_CUDA_CONTENT(v)->priv == NULL)
  {
    N_VDestroy(v);
    return (NULL);
  }

  // Initialize content
  NVEC_CUDA_CONTENT(v)->length             = 0;
  NVEC_CUDA_CONTENT(v)->host_data          = NULL;
  NVEC_CUDA_CONTENT(v)->device_data        = NULL;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy = NULL;
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy = NULL;
  NVEC_CUDA_CONTENT(v)->mem_helper         = NULL;
  NVEC_CUDA_CONTENT(v)->own_helper         = SUNFALSE;

  // Initialize private content
  NVEC_CUDA_PRIVATE(v)->use_managed_mem     = SUNFALSE;
  NVEC_CUDA_PRIVATE(v)->device_counter      = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev   = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_host  = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_bytes = 0;
  NVEC_CUDA_PRIVATE(v)->fused_buffer_dev    = NULL;
  NVEC_CUDA_PRIVATE(v)->fused_buffer_host   = NULL;
  NVEC_CUDA_PRIVATE(v)->fused_buffer_bytes  = 0;
  NVEC_CUDA_PRIVATE(v)->fused_buffer_offset = 0;

  return (v);
}

N_Vector N_VNew_Cuda(sunindextype length, SUNContext sunctx)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Cuda(sunctx);
  if (v == NULL) { return (NULL); }

  NVEC_CUDA_CONTENT(v)->length     = length;
  NVEC_CUDA_CONTENT(v)->mem_helper = SUNMemoryHelper_Cuda(sunctx);
  NVEC_CUDA_CONTENT(v)->stream_exec_policy = DEFAULT_STREAMING_EXECPOLICY.clone();
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy = DEFAULT_REDUCTION_EXECPOLICY.clone();
  NVEC_CUDA_CONTENT(v)->own_helper      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem = SUNFALSE;

  if (NVEC_CUDA_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNew_Cuda: memory helper is NULL\n");
    N_VDestroy(v);
    return (NULL);
  }

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VNew_Cuda: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return (NULL);
  }

  return (v);
}

N_Vector N_VNewWithMemHelp_Cuda(sunindextype length,
                                sunbooleantype use_managed_mem,
                                SUNMemoryHelper helper, SUNContext sunctx)
{
  N_Vector v;

  if (helper == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Cuda: helper is NULL\n");
    return (NULL);
  }

  if (!SUNMemoryHelper_ImplementsRequiredOps(helper))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Cuda: helper doesn't "
                         "implement all required ops\n");
    return (NULL);
  }

  v = NULL;
  v = N_VNewEmpty_Cuda(sunctx);
  if (v == NULL) { return (NULL); }

  NVEC_CUDA_CONTENT(v)->length     = length;
  NVEC_CUDA_CONTENT(v)->mem_helper = helper;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy = DEFAULT_STREAMING_EXECPOLICY.clone();
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy = DEFAULT_REDUCTION_EXECPOLICY.clone();
  NVEC_CUDA_CONTENT(v)->own_helper      = SUNFALSE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem = use_managed_mem;

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VNewWithMemHelp_Cuda: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return (NULL);
  }

  return (v);
}

N_Vector N_VNewManaged_Cuda(sunindextype length, SUNContext sunctx)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Cuda(sunctx);
  if (v == NULL) { return (NULL); }

  NVEC_CUDA_CONTENT(v)->length = length;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy = DEFAULT_STREAMING_EXECPOLICY.clone();
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy = DEFAULT_REDUCTION_EXECPOLICY.clone();
  NVEC_CUDA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Cuda(sunctx);
  NVEC_CUDA_CONTENT(v)->own_helper      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem = SUNTRUE;

  if (NVEC_CUDA_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VNewManaged_Cuda: memory helper is NULL\n");
    N_VDestroy(v);
    return (NULL);
  }

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VNewManaged_Cuda: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return (NULL);
  }

  return (v);
}

N_Vector N_VMake_Cuda(sunindextype length, sunrealtype* h_vdata,
                      sunrealtype* d_vdata, SUNContext sunctx)
{
  N_Vector v;

  if (h_vdata == NULL || d_vdata == NULL) { return (NULL); }

  v = NULL;
  v = N_VNewEmpty_Cuda(sunctx);
  if (v == NULL) { return (NULL); }

  NVEC_CUDA_CONTENT(v)->length     = length;
  NVEC_CUDA_CONTENT(v)->mem_helper = SUNMemoryHelper_Cuda(sunctx);
  NVEC_CUDA_CONTENT(v)->host_data =
    SUNMemoryHelper_Wrap(NVEC_CUDA_MEMHELP(v), h_vdata, SUNMEMTYPE_HOST);
  NVEC_CUDA_CONTENT(v)->device_data =
    SUNMemoryHelper_Wrap(NVEC_CUDA_MEMHELP(v), d_vdata, SUNMEMTYPE_DEVICE);
  NVEC_CUDA_CONTENT(v)->stream_exec_policy = DEFAULT_STREAMING_EXECPOLICY.clone();
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy = DEFAULT_REDUCTION_EXECPOLICY.clone();
  NVEC_CUDA_CONTENT(v)->own_helper      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem = SUNFALSE;

  if (NVEC_CUDA_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMake_Cuda: memory helper is NULL\n");
    N_VDestroy(v);
    return (NULL);
  }

  if (NVEC_CUDA_CONTENT(v)->device_data == NULL ||
      NVEC_CUDA_CONTENT(v)->host_data == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VMake_Cuda: SUNMemoryHelper_Wrap returned NULL\n");
    N_VDestroy(v);
    return (NULL);
  }

  return (v);
}

N_Vector N_VMakeManaged_Cuda(sunindextype length, sunrealtype* vdata,
                             SUNContext sunctx)
{
  N_Vector v;

  if (vdata == NULL) { return (NULL); }

  v = NULL;
  v = N_VNewEmpty_Cuda(sunctx);
  if (v == NULL) { return (NULL); }

  NVEC_CUDA_CONTENT(v)->length     = length;
  NVEC_CUDA_CONTENT(v)->mem_helper = SUNMemoryHelper_Cuda(sunctx);
  NVEC_CUDA_CONTENT(v)->host_data  = SUNMemoryHelper_Wrap(NVEC_CUDA_MEMHELP(v),
                                                          vdata, SUNMEMTYPE_UVM);
  NVEC_CUDA_CONTENT(v)->device_data =
    SUNMemoryHelper_Alias(NVEC_CUDA_MEMHELP(v), NVEC_CUDA_CONTENT(v)->host_data);
  NVEC_CUDA_CONTENT(v)->stream_exec_policy = DEFAULT_STREAMING_EXECPOLICY.clone();
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy = DEFAULT_REDUCTION_EXECPOLICY.clone();
  NVEC_CUDA_CONTENT(v)->own_helper      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem = SUNTRUE;

  if (NVEC_CUDA_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VMakeManaged_Cuda: memory helper is NULL\n");
    N_VDestroy(v);
    return (NULL);
  }

  if (NVEC_CUDA_CONTENT(v)->device_data == NULL ||
      NVEC_CUDA_CONTENT(v)->host_data == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in N_VMakeManaged_Cuda: SUNMemoryHelper_Wrap returned NULL\n");
    N_VDestroy(v);
    return (NULL);
  }

  return (v);
}

/* ----------------------------------------------------------------------------
 * Set pointer to the raw host data. Does not free the existing pointer.
 */

void nvSetHostArrayPointer_Cuda(sunrealtype* h_vdata, N_Vector v)
{
  if (N_VIsManagedMemory_Cuda(v))
  {
    if (NVEC_CUDA_CONTENT(v)->host_data)
    {
      NVEC_CUDA_CONTENT(v)->host_data->ptr   = (void*)h_vdata;
      NVEC_CUDA_CONTENT(v)->device_data->ptr = (void*)h_vdata;
    }
    else
    {
      NVEC_CUDA_CONTENT(v)->host_data =
        SUNMemoryHelper_Wrap(NVEC_CUDA_MEMHELP(v), (void*)h_vdata,
                             SUNMEMTYPE_UVM);
      NVEC_CUDA_CONTENT(v)->device_data =
        SUNMemoryHelper_Alias(NVEC_CUDA_MEMHELP(v),
                              NVEC_CUDA_CONTENT(v)->host_data);
    }
  }
  else
  {
    if (NVEC_CUDA_CONTENT(v)->host_data)
    {
      NVEC_CUDA_CONTENT(v)->host_data->ptr = (void*)h_vdata;
    }
    else
    {
      NVEC_CUDA_CONTENT(v)->host_data =
        SUNMemoryHelper_Wrap(NVEC_CUDA_MEMHELP(v), (void*)h_vdata,
                             SUNMEMTYPE_HOST);
    }
  }
}

/* ----------------------------------------------------------------------------
 * Set pointer to the raw device data
 */

void nvSetDeviceArrayPointer_Cuda(sunrealtype* d_vdata, N_Vector v)
{
  if (N_VIsManagedMemory_Cuda(v))
  {
    if (NVEC_CUDA_CONTENT(v)->device_data)
    {
      NVEC_CUDA_CONTENT(v)->device_data->ptr = (void*)d_vdata;
      NVEC_CUDA_CONTENT(v)->host_data->ptr   = (void*)d_vdata;
    }
    else
    {
      NVEC_CUDA_CONTENT(v)->device_data =
        SUNMemoryHelper_Wrap(NVEC_CUDA_MEMHELP(v), (void*)d_vdata,
                             SUNMEMTYPE_UVM);
      NVEC_CUDA_CONTENT(v)->host_data =
        SUNMemoryHelper_Alias(NVEC_CUDA_MEMHELP(v),
                              NVEC_CUDA_CONTENT(v)->device_data);
    }
  }
  else
  {
    if (NVEC_CUDA_CONTENT(v)->device_data)
    {
      NVEC_CUDA_CONTENT(v)->device_data->ptr = (void*)d_vdata;
    }
    else
    {
      NVEC_CUDA_CONTENT(v)->device_data =
        SUNMemoryHelper_Wrap(NVEC_CUDA_MEMHELP(v), (void*)d_vdata,
                             SUNMEMTYPE_DEVICE);
    }
  }
}

/* ----------------------------------------------------------------------------
 * Return a flag indicating if the memory for the vector data is managed
 */

sunbooleantype N_VIsManagedMemory_Cuda(N_Vector x)
{
  return NVEC_CUDA_PRIVATE(x)->use_managed_mem;
}

SUNErrCode N_VSetKernelExecPolicy_Cuda(N_Vector x,
                                       SUNCudaExecPolicy* stream_exec_policy,
                                       SUNCudaExecPolicy* reduce_exec_policy)
{
  if (x == NULL) { return SUN_ERR_GENERIC; }

  /* Delete the old policies */
  delete NVEC_CUDA_CONTENT(x)->stream_exec_policy;
  delete NVEC_CUDA_CONTENT(x)->reduce_exec_policy;

  /* Reset the policy if it is null */

  if (stream_exec_policy == NULL)
  {
    NVEC_CUDA_CONTENT(x)->stream_exec_policy =
      DEFAULT_STREAMING_EXECPOLICY.clone();
  }
  else
  {
    NVEC_CUDA_CONTENT(x)->stream_exec_policy = stream_exec_policy->clone();
  }

  if (reduce_exec_policy == NULL)
  {
    NVEC_CUDA_CONTENT(x)->reduce_exec_policy =
      DEFAULT_REDUCTION_EXECPOLICY.clone();
  }
  else
  {
    NVEC_CUDA_CONTENT(x)->reduce_exec_policy = reduce_exec_policy->clone();
  }

  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Copy vector data to the device
 */

void N_VCopyToDevice_Cuda(N_Vector x)
{
  int copy_fail;

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_CUDA_MEMHELP(x),
                                        NVEC_CUDA_CONTENT(x)->device_data,
                                        NVEC_CUDA_CONTENT(x)->host_data,
                                        NVEC_CUDA_MEMSIZE(x),
                                        (void*)NVEC_CUDA_STREAM(x));

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VCopyToDevice_Cuda: "
                         "SUNMemoryHelper_CopyAsync returned nonzero\n");
  }

  /* we synchronize with respect to the host, but only in this stream */
  SUNDIALS_CUDA_VERIFY(cudaStreamSynchronize(*NVEC_CUDA_STREAM(x)));
}

/* ----------------------------------------------------------------------------
 * Copy vector data from the device to the host
 */

void N_VCopyFromDevice_Cuda(N_Vector x)
{
  int copy_fail;

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_CUDA_MEMHELP(x),
                                        NVEC_CUDA_CONTENT(x)->host_data,
                                        NVEC_CUDA_CONTENT(x)->device_data,
                                        NVEC_CUDA_MEMSIZE(x),
                                        (void*)NVEC_CUDA_STREAM(x));

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VCopyFromDevice_Cuda: "
                         "SUNMemoryHelper_CopyAsync returned nonzero\n");
  }

  /* we synchronize with respect to the host, but only in this stream */
  SUNDIALS_CUDA_VERIFY(cudaStreamSynchronize(*NVEC_CUDA_STREAM(x)));
}

/* ----------------------------------------------------------------------------
 * Function to print the a CUDA-based vector to stdout
 */

void nvPrint_Cuda(N_Vector x) { nvPrintFile_Cuda(x, stdout); }

/* ----------------------------------------------------------------------------
 * Function to print the a CUDA-based vector to outfile
 */

void nvPrintFile_Cuda(N_Vector x, FILE* outfile)
{
  sunindextype i;

#ifdef SUNDIALS_DEBUG_PRINTVEC
  N_VCopyFromDevice_Cuda(x);
#endif

  for (i = 0; i < NVEC_CUDA_CONTENT(x)->length; i++)
  {
    fprintf(outfile, SUN_FORMAT_E "\n", NVEC_CUDA_HDATAp(x)[i]);
  }

  return;
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector nvCloneEmpty_Cuda(N_Vector w)
{
  N_Vector v;

  if (w == NULL) { return (NULL); }

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty_Cuda(w->sunctx);
  if (v == NULL) { return (NULL); }

  /* Attach operations */
  if (N_VCopyOps(w, v))
  {
    N_VDestroy(v);
    return (NULL);
  }

  /* Set content */
  NVEC_CUDA_CONTENT(v)->length          = NVEC_CUDA_CONTENT(w)->length;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem = NVEC_CUDA_PRIVATE(w)->use_managed_mem;

  return (v);
}

N_Vector nvClone_Cuda(N_Vector w)
{
  N_Vector v;

  v = NULL;
  v = nvCloneEmpty_Cuda(w);
  if (v == NULL) { return (NULL); }

  NVEC_CUDA_MEMHELP(v) = SUNMemoryHelper_Clone(NVEC_CUDA_MEMHELP(w));
  NVEC_CUDA_CONTENT(v)->own_helper = SUNTRUE;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy =
    NVEC_CUDA_CONTENT(w)->stream_exec_policy->clone();
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy =
    NVEC_CUDA_CONTENT(w)->reduce_exec_policy->clone();

  if (NVEC_CUDA_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvClone_Cuda: SUNMemoryHelper_Clone returned NULL\n");
    N_VDestroy(v);
    return (NULL);
  }

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvClone_Cuda: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return (NULL);
  }

  return (v);
}

void nvDestroy_Cuda(N_Vector v)
{
  N_VectorContent_Cuda vc;
  N_PrivateVectorContent_Cuda vcp;

  if (v == NULL) { return; }

  /* free ops structure */
  if (v->ops != NULL)
  {
    free(v->ops);
    v->ops = NULL;
  }

  /* extract content */
  vc = NVEC_CUDA_CONTENT(v);
  if (vc == NULL)
  {
    free(v);
    v = NULL;
    return;
  }

  /* free private content */
  vcp = (N_PrivateVectorContent_Cuda)vc->priv;
  if (vcp != NULL)
  {
    /* free items in private content */
    FreeDeviceCounter(v);
    FreeReductionBuffer(v);
    FusedBuffer_Free(v);
    free(vcp);
    vc->priv = NULL;
  }

  /* free items in content */
  if (NVEC_CUDA_MEMHELP(v))
  {
    SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v), vc->host_data,
                            (void*)NVEC_CUDA_STREAM(v));
    vc->host_data = NULL;
    SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v), vc->device_data,
                            (void*)NVEC_CUDA_STREAM(v));
    vc->device_data = NULL;
    if (vc->own_helper) { SUNMemoryHelper_Destroy(vc->mem_helper); }
    vc->mem_helper = NULL;
  }

  /* we can delete the exec policies now that we are done with the streams */
  delete vc->stream_exec_policy;
  delete vc->reduce_exec_policy;

  /* free content struct */
  free(vc);

  /* free vector */
  free(v);

  return;
}

void nvConst_Cuda(sunrealtype a, N_Vector X)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvConst_Cuda: GetKernelParameters returned nonzero\n");
  }

  setConstKernel<<<grid, block, shMemSize, stream>>>(a, NVEC_CUDA_DDATAp(X),
                                                     NVEC_CUDA_CONTENT(X)->length);
  PostKernelLaunch();
}

void nvLinearSum_Cuda(sunrealtype a, N_Vector X, sunrealtype b, N_Vector Y,
                      N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvLinearSum_Cuda: GetKernelParameters returned nonzero\n");
  }

  linearSumKernel<<<grid, block, shMemSize, stream>>>(a, NVEC_CUDA_DDATAp(X), b,
                                                      NVEC_CUDA_DDATAp(Y),
                                                      NVEC_CUDA_DDATAp(Z),
                                                      NVEC_CUDA_CONTENT(X)->length);
  PostKernelLaunch();
}

void nvProd_Cuda(N_Vector X, N_Vector Y, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvProd_Cuda: GetKernelParameters returned nonzero\n");
  }

  prodKernel<<<grid, block, shMemSize, stream>>>(NVEC_CUDA_DDATAp(X),
                                                 NVEC_CUDA_DDATAp(Y),
                                                 NVEC_CUDA_DDATAp(Z),
                                                 NVEC_CUDA_CONTENT(X)->length);
  PostKernelLaunch();
}

void nvDiv_Cuda(N_Vector X, N_Vector Y, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvDiv_Cuda: GetKernelParameters returned nonzero\n");
  }

  divKernel<<<grid, block, shMemSize, stream>>>(NVEC_CUDA_DDATAp(X),
                                                NVEC_CUDA_DDATAp(Y),
                                                NVEC_CUDA_DDATAp(Z),
                                                NVEC_CUDA_CONTENT(X)->length);
  PostKernelLaunch();
}

void nvScale_Cuda(sunrealtype a, N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvScale_Cuda: GetKernelParameters returned nonzero\n");
  }

  scaleKernel<<<grid, block, shMemSize, stream>>>(a, NVEC_CUDA_DDATAp(X),
                                                  NVEC_CUDA_DDATAp(Z),
                                                  NVEC_CUDA_CONTENT(X)->length);
  PostKernelLaunch();
}

void nvAbs_Cuda(N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvAbs_Cuda: GetKernelParameters returned nonzero\n");
  }

  absKernel<<<grid, block, shMemSize, stream>>>(NVEC_CUDA_DDATAp(X),
                                                NVEC_CUDA_DDATAp(Z),
                                                NVEC_CUDA_CONTENT(X)->length);
  PostKernelLaunch();
}

void nvInv_Cuda(N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvInv_Cuda: GetKernelParameters returned nonzero\n");
  }

  invKernel<<<grid, block, shMemSize, stream>>>(NVEC_CUDA_DDATAp(X),
                                                NVEC_CUDA_DDATAp(Z),
                                                NVEC_CUDA_CONTENT(X)->length);
  PostKernelLaunch();
}

void nvAddConst_Cuda(N_Vector X, sunrealtype b, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvAddConst_Cuda: GetKernelParameters returned nonzero\n");
  }

  addConstKernel<<<grid, block, shMemSize, stream>>>(b, NVEC_CUDA_DDATAp(X),
                                                     NVEC_CUDA_DDATAp(Z),
                                                     NVEC_CUDA_CONTENT(X)->length);
  PostKernelLaunch();
}

sunrealtype nvDotProd_Cuda(N_Vector X, N_Vector Y)
{
  bool atomic;
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  sunrealtype gpu_result = ZERO;

  if (GetKernelParameters(X, true, grid, block, shMemSize, stream, atomic))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvDotProd_Cuda: GetKernelParameters returned nonzero\n");
  }

  // When using atomic reductions, we only need one output value
  const size_t buffer_size = atomic ? 1 : grid;
  if (InitializeReductionBuffer(X, gpu_result, buffer_size))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvDotProd_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  if (atomic)
  {
    dotProdKernel<sunrealtype, sunindextype, GridReducerAtomic>
      <<<grid, block, shMemSize, stream>>>(NVEC_CUDA_DDATAp(X),
                                           NVEC_CUDA_DDATAp(Y),
                                           NVEC_CUDA_DBUFFERp(X),
                                           NVEC_CUDA_CONTENT(X)->length, nullptr);
  }
  else
  {
    dotProdKernel<sunrealtype, sunindextype, GridReducerLDS>
      <<<grid, block, shMemSize, stream>>>(NVEC_CUDA_DDATAp(X),
                                           NVEC_CUDA_DDATAp(Y),
                                           NVEC_CUDA_DBUFFERp(X),
                                           NVEC_CUDA_CONTENT(X)->length,
                                           NVEC_CUDA_DCOUNTERp(X));
  }
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return gpu_result;
}

sunrealtype nvMaxNorm_Cuda(N_Vector X)
{
  bool atomic;
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  sunrealtype gpu_result = ZERO;

  if (GetKernelParameters(X, true, grid, block, shMemSize, stream, atomic))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvMaxNorm_Cuda: GetKernelParameters returned nonzero\n");
  }

  // When using atomic reductions, we only need one output value
  const size_t buffer_size = atomic ? 1 : grid;
  if (InitializeReductionBuffer(X, gpu_result, buffer_size))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvMaxNorm_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  if (atomic)
  {
    maxNormKernel<sunrealtype, sunindextype, GridReducerAtomic>
      <<<grid, block, shMemSize, stream>>>(NVEC_CUDA_DDATAp(X),
                                           NVEC_CUDA_DBUFFERp(X),
                                           NVEC_CUDA_CONTENT(X)->length, nullptr);
  }
  else
  {
    maxNormKernel<sunrealtype, sunindextype, GridReducerLDS>
      <<<grid, block, shMemSize, stream>>>(NVEC_CUDA_DDATAp(X),
                                           NVEC_CUDA_DBUFFERp(X),
                                           NVEC_CUDA_CONTENT(X)->length,
                                           NVEC_CUDA_DCOUNTERp(X));
  }

  PostKernelLaunch();

  // Finish reduction on CPU if there are less than two blocks of data left.
  CopyReductionBufferFromDevice(X);
  gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return gpu_result;
}

sunrealtype nvWSqrSumLocal_Cuda(N_Vector X, N_Vector W)
{
  bool atomic;
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  sunrealtype gpu_result = ZERO;

  if (GetKernelParameters(X, true, grid, block, shMemSize, stream, atomic))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvWSqrSumLocal_Cuda: GetKernelParameters returned nonzero\n");
  }

  const size_t buffer_size = atomic ? 1 : grid;
  if (InitializeReductionBuffer(X, gpu_result, buffer_size))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWSqrSumLocal_Cuda: "
                         "InitializeReductionBuffer returned nonzero\n");
  }

  if (atomic)
  {
    wL2NormSquareKernel<sunrealtype, sunindextype, GridReducerAtomic>
      <<<grid, block, shMemSize, stream>>>(NVEC_CUDA_DDATAp(X),
                                           NVEC_CUDA_DDATAp(W),
                                           NVEC_CUDA_DBUFFERp(X),
                                           NVEC_CUDA_CONTENT(X)->length, nullptr);
  }
  else
  {
    wL2NormSquareKernel<sunrealtype, sunindextype, GridReducerLDS>
      <<<grid, block, shMemSize, stream>>>(NVEC_CUDA_DDATAp(X),
                                           NVEC_CUDA_DDATAp(W),
                                           NVEC_CUDA_DBUFFERp(X),
                                           NVEC_CUDA_CONTENT(X)->length,
                                           NVEC_CUDA_DCOUNTERp(X));
  }

  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return gpu_result;
}

sunrealtype nvWrmsNorm_Cuda(N_Vector X, N_Vector W)
{
  const sunrealtype sum = nvWSqrSumLocal_Cuda(X, W);
  return std::sqrt(sum / NVEC_CUDA_CONTENT(X)->length);
}

sunrealtype nvWSqrSumMaskLocal_Cuda(N_Vector X, N_Vector W, N_Vector Id)
{
  bool atomic;
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  sunrealtype gpu_result = ZERO;

  if (GetKernelParameters(X, true, grid, block, shMemSize, stream, atomic))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWSqrSumMaskLocal_Cuda: "
                         "GetKernelParameters returned nonzero\n");
  }

  const size_t buffer_size = atomic ? 1 : grid;
  if (InitializeReductionBuffer(X, gpu_result, buffer_size))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWSqrSumMaskLocal_Cuda: "
                         "InitializeReductionBuffer returned nonzero\n");
  }

  if (atomic)
  {
    wL2NormSquareMaskKernel<sunrealtype, sunindextype, GridReducerAtomic>
      <<<grid, block, shMemSize, stream>>>(NVEC_CUDA_DDATAp(X),
                                           NVEC_CUDA_DDATAp(W),
                                           NVEC_CUDA_DDATAp(Id),
                                           NVEC_CUDA_DBUFFERp(X),
                                           NVEC_CUDA_CONTENT(X)->length, nullptr);
  }
  else
  {
    wL2NormSquareMaskKernel<sunrealtype, sunindextype, GridReducerLDS>
      <<<grid, block, shMemSize, stream>>>(NVEC_CUDA_DDATAp(X),
                                           NVEC_CUDA_DDATAp(W),
                                           NVEC_CUDA_DDATAp(Id),
                                           NVEC_CUDA_DBUFFERp(X),
                                           NVEC_CUDA_CONTENT(X)->length,
                                           NVEC_CUDA_DCOUNTERp(X));
  }

  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return gpu_result;
}

sunrealtype nvWrmsNormMask_Cuda(N_Vector X, N_Vector W, N_Vector Id)
{
  const sunrealtype sum = nvWSqrSumMaskLocal_Cuda(X, W, Id);
  return std::sqrt(sum / NVEC_CUDA_CONTENT(X)->length);
}

sunrealtype nvMin_Cuda(N_Vector X)
{
  bool atomic;
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  sunrealtype gpu_result = std::numeric_limits<sunrealtype>::max();

  if (GetKernelParameters(X, true, grid, block, shMemSize, stream, atomic))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvMin_Cuda: GetKernelParameters returned nonzero\n");
  }

  const size_t buffer_size = atomic ? 1 : grid;
  if (InitializeReductionBuffer(X, gpu_result, buffer_size))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvMin_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  if (atomic)
  {
    findMinKernel<sunrealtype, sunindextype, GridReducerAtomic>
      <<<grid, block, shMemSize, stream>>>(gpu_result, NVEC_CUDA_DDATAp(X),
                                           NVEC_CUDA_DBUFFERp(X),
                                           NVEC_CUDA_CONTENT(X)->length, nullptr);
  }
  else
  {
    findMinKernel<sunrealtype, sunindextype, GridReducerLDS>
      <<<grid, block, shMemSize, stream>>>(gpu_result, NVEC_CUDA_DDATAp(X),
                                           NVEC_CUDA_DBUFFERp(X),
                                           NVEC_CUDA_CONTENT(X)->length,
                                           NVEC_CUDA_DCOUNTERp(X));
  }

  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return gpu_result;
}

sunrealtype nvWL2Norm_Cuda(N_Vector X, N_Vector W)
{
  const sunrealtype sum = nvWSqrSumLocal_Cuda(X, W);
  return std::sqrt(sum);
}

sunrealtype nvL1Norm_Cuda(N_Vector X)
{
  bool atomic;
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  sunrealtype gpu_result = ZERO;

  if (GetKernelParameters(X, true, grid, block, shMemSize, stream, atomic))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvL1Norm_Cuda: GetKernelParameters returned nonzero\n");
  }

  const size_t buffer_size = atomic ? 1 : grid;
  if (InitializeReductionBuffer(X, gpu_result, buffer_size))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvL1Norm_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  if (atomic)
  {
    L1NormKernel<sunrealtype, sunindextype, GridReducerAtomic>
      <<<grid, block, shMemSize, stream>>>(NVEC_CUDA_DDATAp(X),
                                           NVEC_CUDA_DBUFFERp(X),
                                           NVEC_CUDA_CONTENT(X)->length, nullptr);
  }
  else
  {
    L1NormKernel<sunrealtype, sunindextype, GridReducerLDS>
      <<<grid, block, shMemSize, stream>>>(NVEC_CUDA_DDATAp(X),
                                           NVEC_CUDA_DBUFFERp(X),
                                           NVEC_CUDA_CONTENT(X)->length,
                                           NVEC_CUDA_DCOUNTERp(X));
  }

  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return gpu_result;
}

void nvCompare_Cuda(sunrealtype c, N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvCompare_Cuda: GetKernelParameters returned nonzero\n");
  }

  compareKernel<<<grid, block, shMemSize, stream>>>(c, NVEC_CUDA_DDATAp(X),
                                                    NVEC_CUDA_DDATAp(Z),
                                                    NVEC_CUDA_CONTENT(X)->length);
  PostKernelLaunch();
}

sunbooleantype nvInvTest_Cuda(N_Vector X, N_Vector Z)
{
  bool atomic;
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  sunrealtype gpu_result = ZERO;

  if (GetKernelParameters(X, true, grid, block, shMemSize, stream, atomic))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvInvTest_Cuda: GetKernelParameters returned nonzero\n");
  }

  const size_t buffer_size = atomic ? 1 : grid;
  if (InitializeReductionBuffer(X, gpu_result, buffer_size))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvInvTest_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  if (atomic)
  {
    invTestKernel<sunrealtype, sunindextype, GridReducerAtomic>
      <<<grid, block, shMemSize, stream>>>(NVEC_CUDA_DDATAp(X),
                                           NVEC_CUDA_DDATAp(Z),
                                           NVEC_CUDA_DBUFFERp(X),
                                           NVEC_CUDA_CONTENT(X)->length, nullptr);
  }
  else
  {
    invTestKernel<sunrealtype, sunindextype, GridReducerLDS>
      <<<grid, block, shMemSize, stream>>>(NVEC_CUDA_DDATAp(X),
                                           NVEC_CUDA_DDATAp(Z),
                                           NVEC_CUDA_DBUFFERp(X),
                                           NVEC_CUDA_CONTENT(X)->length,
                                           NVEC_CUDA_DCOUNTERp(X));
  }

  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return (gpu_result < HALF);
}

sunbooleantype nvConstrMask_Cuda(N_Vector C, N_Vector X, N_Vector M)
{
  bool atomic;
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  sunrealtype gpu_result = ZERO;

  if (GetKernelParameters(X, true, grid, block, shMemSize, stream, atomic))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvConstrMask_Cuda: GetKernelParameters returned nonzero\n");
  }

  const size_t buffer_size = atomic ? 1 : grid;
  if (InitializeReductionBuffer(X, gpu_result, buffer_size))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvConstrMask_Cuda: "
                         "InitializeReductionBuffer returned nonzero\n");
  }

  if (atomic)
  {
    constrMaskKernel<sunrealtype, sunindextype, GridReducerAtomic>
      <<<grid, block, shMemSize, stream>>>(NVEC_CUDA_DDATAp(C),
                                           NVEC_CUDA_DDATAp(X),
                                           NVEC_CUDA_DDATAp(M),
                                           NVEC_CUDA_DBUFFERp(X),
                                           NVEC_CUDA_CONTENT(X)->length, nullptr);
  }
  else
  {
    constrMaskKernel<sunrealtype, sunindextype, GridReducerLDS>
      <<<grid, block, shMemSize, stream>>>(NVEC_CUDA_DDATAp(C),
                                           NVEC_CUDA_DDATAp(X),
                                           NVEC_CUDA_DDATAp(M),
                                           NVEC_CUDA_DBUFFERp(X),
                                           NVEC_CUDA_CONTENT(X)->length,
                                           NVEC_CUDA_DCOUNTERp(X));
  }

  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return (gpu_result < HALF);
}

sunrealtype nvMinQuotient_Cuda(N_Vector num, N_Vector denom)
{
  bool atomic;
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  sunrealtype gpu_result = std::numeric_limits<sunrealtype>::max();
  ;

  if (GetKernelParameters(num, true, grid, block, shMemSize, stream, atomic))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvMinQuotient_Cuda: GetKernelParameters returned nonzero\n");
  }

  const size_t buffer_size = atomic ? 1 : grid;
  if (InitializeReductionBuffer(num, gpu_result, buffer_size))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvMinQuotient_Cuda: "
                         "InitializeReductionBuffer returned nonzero\n");
  }

  if (atomic)
  {
    minQuotientKernel<sunrealtype, sunindextype, GridReducerAtomic>
      <<<grid, block, shMemSize, stream>>>(gpu_result, NVEC_CUDA_DDATAp(num),
                                           NVEC_CUDA_DDATAp(denom),
                                           NVEC_CUDA_DBUFFERp(num),
                                           NVEC_CUDA_CONTENT(num)->length,
                                           nullptr);
  }
  else
  {
    minQuotientKernel<sunrealtype, sunindextype, GridReducerLDS>
      <<<grid, block, shMemSize, stream>>>(gpu_result, NVEC_CUDA_DDATAp(num),
                                           NVEC_CUDA_DDATAp(denom),
                                           NVEC_CUDA_DBUFFERp(num),
                                           NVEC_CUDA_CONTENT(num)->length,
                                           NVEC_CUDA_DCOUNTERp(num));
  }

  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(num);
  gpu_result = NVEC_CUDA_HBUFFERp(num)[0];

  return gpu_result;
}

/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */

SUNErrCode nvLinearCombination_Cuda(int nvec, sunrealtype* c, N_Vector* X,
                                    N_Vector z)
{
  // Fused op workspace shortcuts
  sunrealtype* cdata  = NULL;
  sunrealtype** xdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(z, nvec, nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombination_Cuda: FusedBuffer_Init "
                         "returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyRealArray(z, c, nvec, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombination_Cuda: "
                         "FusedBuffer_CopyRealArray returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(z, X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombination_Cuda: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyToDevice(z))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombination_Cuda: "
                         "FusedBuffer_CopyToDevice returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  // Set kernel parameters and launch
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X[0], false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombination_Cuda: "
                         "GetKernelParameters returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  linearCombinationKernel<<<grid, block, shMemSize, stream>>>(nvec, cdata, xdata,
                                                              NVEC_CUDA_DDATAp(z),
                                                              NVEC_CUDA_CONTENT(z)
                                                                ->length);
  PostKernelLaunch();

  return SUN_SUCCESS;
}

SUNErrCode nvScaleAddMulti_Cuda(int nvec, sunrealtype* c, N_Vector x,
                                N_Vector* Y, N_Vector* Z)
{
  // Shortcuts to the fused op workspace
  sunrealtype* cdata  = NULL;
  sunrealtype** ydata = NULL;
  sunrealtype** zdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(x, nvec, 2 * nvec))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvScaleAddMulti_Cuda: FusedBuffer_Init returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyRealArray(x, c, nvec, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleAddMulti_Cuda: "
                         "FusedBuffer_CopyRealArray returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(x, Y, nvec, &ydata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleAddMulti_Cuda: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(x, Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleAddMulti_Cuda: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyToDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleAddMulti_Cuda: "
                         "FusedBuffer_CopyToDevice returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(x, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvScaleAddMulti_Cuda: GetKernelParameters returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  scaleAddMultiKernel<<<grid, block, shMemSize, stream>>>(nvec, cdata,
                                                          NVEC_CUDA_DDATAp(x),
                                                          ydata, zdata,
                                                          NVEC_CUDA_CONTENT(x)
                                                            ->length);
  PostKernelLaunch();

  return SUN_SUCCESS;
}

SUNErrCode nvDotProdMulti_Cuda(int nvec, N_Vector x, N_Vector* Y,
                               sunrealtype* dots)
{
  // Fused op workspace shortcuts
  sunrealtype** ydata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(x, 0, nvec))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvDotProdMulti_Cuda: FusedBuffer_Init returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(x, Y, nvec, &ydata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvDotProdMulti_Cuda: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyToDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvDotProdMulti_Cuda: "
                         "FusedBuffer_CopyToDevice returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(x, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvDotProdMulti_Cuda: GetKernelParameters returned nonzero\n");
    return SUN_ERR_GENERIC;
  }
  grid = nvec;

  if (InitializeReductionBuffer(x, ZERO, nvec))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvDotProd_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  dotProdMultiKernel<sunrealtype, sunindextype, GridReducerAtomic>
    <<<grid, block, shMemSize, stream>>>(nvec, NVEC_CUDA_DDATAp(x), ydata,
                                         NVEC_CUDA_DBUFFERp(x),
                                         NVEC_CUDA_CONTENT(x)->length);

  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(x, nvec);
  for (int i = 0; i < nvec; ++i) { dots[i] = NVEC_CUDA_HBUFFERp(x)[i]; }

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------------------
 */

SUNErrCode nvLinearSumVectorArray_Cuda(int nvec, sunrealtype a, N_Vector* X,
                                       sunrealtype b, N_Vector* Y, N_Vector* Z)
{
  // Shortcuts to the fused op workspace
  sunrealtype** xdata = NULL;
  sunrealtype** ydata = NULL;
  sunrealtype** zdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(Z[0], 0, 3 * nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearSumVectorArray_Cuda: "
                         "FusedBuffer_Init returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearSumVectorArray_Cuda: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Y, nvec, &ydata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearSumVectorArray_Cuda: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearSumVectorArray_Cuda: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyToDevice(Z[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinaerSumVectorArray_Cuda: "
                         "FusedBuffer_CopyToDevice returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(Z[0], false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearSumVectorArray_Cuda: "
                         "GetKernelParameters returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  linearSumVectorArrayKernel<<<grid, block, shMemSize, stream>>>(nvec, a, xdata,
                                                                 b, ydata, zdata,
                                                                 NVEC_CUDA_CONTENT(
                                                                   Z[0])
                                                                   ->length);
  PostKernelLaunch();

  return SUN_SUCCESS;
}

SUNErrCode nvScaleVectorArray_Cuda(int nvec, sunrealtype* c, N_Vector* X,
                                   N_Vector* Z)
{
  // Shortcuts to the fused op workspace arrays
  sunrealtype* cdata  = NULL;
  sunrealtype** xdata = NULL;
  sunrealtype** zdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(Z[0], nvec, 2 * nvec))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvScaleVectorArray_Cuda: FusedBuffer_Init returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyRealArray(Z[0], c, nvec, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleVectorArray_Cuda: "
                         "FusedBuffer_CopyRealArray returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleVectorArray_Cuda: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleVectorArray_Cuda: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyToDevice(Z[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleVectorArray_Cuda: "
                         "FusedBuffer_CopyToDevice returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(Z[0], false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleVectorArray_Cuda: "
                         "GetKernelParameters returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  scaleVectorArrayKernel<<<grid, block, shMemSize, stream>>>(nvec, cdata, xdata,
                                                             zdata,
                                                             NVEC_CUDA_CONTENT(
                                                               Z[0])
                                                               ->length);
  PostKernelLaunch();

  return SUN_SUCCESS;
}

SUNErrCode nvConstVectorArray_Cuda(int nvec, sunrealtype c, N_Vector* Z)
{
  // Shortcuts to the fused op workspace arrays
  sunrealtype** zdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(Z[0], 0, nvec))
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in nvConstVectorArray_Cuda: FusedBuffer_Init returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvConstVectorArray_Cuda: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyToDevice(Z[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvConstVectorArray_Cuda: "
                         "FusedBuffer_CopyToDevice returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(Z[0], false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvConstVectorArray_Cuda: "
                         "GetKernelParameters returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  constVectorArrayKernel<<<grid, block, shMemSize, stream>>>(nvec, c, zdata,
                                                             NVEC_CUDA_CONTENT(
                                                               Z[0])
                                                               ->length);
  PostKernelLaunch();

  return SUN_SUCCESS;
}

SUNErrCode nvWrmsNormVectorArray_Cuda(int nvec, N_Vector* X, N_Vector* W,
                                      sunrealtype* norms)
{
  // Fused op workspace shortcuts
  sunrealtype** xdata = NULL;
  sunrealtype** wdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(W[0], 0, 2 * nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWrmsNormVectorArray_Cuda: "
                         "FusedBuffer_Init returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(W[0], X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWrmsNormVectorArray_Cuda: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(W[0], W, nvec, &wdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWrmsNormVectorArray_Cuda: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyToDevice(W[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWrmsNormVectorArray_Cuda: "
                         "FusedBuffer_CopyToDevice returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (InitializeReductionBuffer(W[0], ZERO, nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWrmsNormVectorArray_Cuda: "
                         "InitializeReductionBuffer returned nonzero\n");
  }

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(W[0], true, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWrmsNormVectorArray_Cuda: "
                         "GetKernelParameters returned nonzero\n");
    return SUN_ERR_GENERIC;
  }
  grid = nvec;

  wL2NormSquareVectorArrayKernel<sunrealtype, sunindextype, GridReducerAtomic>
    <<<grid, block, shMemSize, stream>>>(nvec, xdata, wdata,
                                         NVEC_CUDA_DBUFFERp(W[0]),
                                         NVEC_CUDA_CONTENT(W[0])->length);
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(W[0], nvec);
  for (int i = 0; i < nvec; ++i)
  {
    norms[i] =
      std::sqrt(NVEC_CUDA_HBUFFERp(W[0])[i] / NVEC_CUDA_CONTENT(W[0])->length);
  }

  return SUN_SUCCESS;
}

SUNErrCode nvWrmsNormMaskVectorArray_Cuda(int nvec, N_Vector* X, N_Vector* W,
                                          N_Vector id, sunrealtype* norms)
{
  // Fused op workspace shortcuts
  sunrealtype** xdata = NULL;
  sunrealtype** wdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(W[0], 0, 2 * nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWrmsNormVectorArray_Cuda: "
                         "FusedBuffer_Init returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(W[0], X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWrmsNormVectorArray_Cuda: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(W[0], W, nvec, &wdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWrmsNormVectorArray_Cuda: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyToDevice(W[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWrmsNormVectorArray_Cuda: "
                         "FusedBuffer_CopyToDevice returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (InitializeReductionBuffer(W[0], ZERO, nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWrmsNormVectorArray_Cuda: "
                         "InitializeReductionBuffer returned nonzero\n");
  }

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(W[0], true, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvWrmsNormMaskVectorArray_Cuda: "
                         "GetKernelParameters returned nonzero\n");
    return SUN_ERR_GENERIC;
  }
  grid = nvec;

  wL2NormSquareMaskVectorArrayKernel<sunrealtype, sunindextype, GridReducerAtomic>
    <<<grid, block, shMemSize, stream>>>(nvec, xdata, wdata, NVEC_CUDA_DDATAp(id),
                                         NVEC_CUDA_DBUFFERp(W[0]),
                                         NVEC_CUDA_CONTENT(W[0])->length);
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(W[0], nvec);
  for (int i = 0; i < nvec; ++i)
  {
    norms[i] =
      std::sqrt(NVEC_CUDA_HBUFFERp(W[0])[i] / NVEC_CUDA_CONTENT(W[0])->length);
  }

  return SUN_SUCCESS;
}

SUNErrCode nvScaleAddMultiVectorArray_Cuda(int nvec, int nsum, sunrealtype* c,
                                           N_Vector* X, N_Vector** Y,
                                           N_Vector** Z)
{
  // Shortcuts to the fused op workspace
  sunrealtype* cdata  = NULL;
  sunrealtype** xdata = NULL;
  sunrealtype** ydata = NULL;
  sunrealtype** zdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(X[0], nsum, nvec + 2 * nvec * nsum))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMultiArray_Cuda: "
                         "FusedBuffer_Init returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyRealArray(X[0], c, nsum, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMultiArray_Cuda: "
                         "FusedBuffer_CopyRealArray returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(X[0], X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleAddMultiVectorArray_Cuda: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray2D(X[0], Y, nvec, nsum, &ydata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleAddMultiVectorArray_Cuda: "
                         "FusedBuffer_CopyPtrArray2D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray2D(X[0], Z, nvec, nsum, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleAddMultiVectorArray_Cuda: "
                         "FusedBuffer_CopyPtrArray2D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyToDevice(X[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleVectorArray_Cuda: "
                         "FusedBuffer_CopyToDevice returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X[0], false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvScaleAddMultiVectorArray_Cuda: "
                         "GetKernelParameters returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  scaleAddMultiVectorArrayKernel<<<grid, block, shMemSize, stream>>>(nvec, nsum,
                                                                     cdata, xdata,
                                                                     ydata, zdata,
                                                                     NVEC_CUDA_CONTENT(
                                                                       X[0])
                                                                       ->length);
  PostKernelLaunch();

  return SUN_SUCCESS;
}

SUNErrCode nvLinearCombinationVectorArray_Cuda(int nvec, int nsum, sunrealtype* c,
                                               N_Vector** X, N_Vector* Z)
{
  // Shortcuts to the fused op workspace arrays
  sunrealtype* cdata  = NULL;
  sunrealtype** xdata = NULL;
  sunrealtype** zdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(Z[0], nsum, nvec + nvec * nsum))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombinationVectorArray_Cuda: "
                         "FusedBuffer_Init returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyRealArray(Z[0], c, nsum, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombinationVectorArray_Cuda: "
                         "FusedBuffer_CopyRealArray returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray2D(Z[0], X, nvec, nsum, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombinationVectorArray_Cuda: "
                         "FusedBuffer_CopyPtrArray2D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombinationVectorArray_Cuda: "
                         "FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  if (FusedBuffer_CopyToDevice(Z[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombinationVectorArray_Cuda: "
                         "FusedBuffer_CopyToDevice returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(Z[0], false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in nvLinearCombinationVectorArray_Cuda: "
                         "GetKernelParameters returned nonzero\n");
    return SUN_ERR_GENERIC;
  }

  linearCombinationVectorArrayKernel<<<grid, block, shMemSize, stream>>>(nvec,
                                                                         nsum,
                                                                         cdata,
                                                                         xdata,
                                                                         zdata,
                                                                         NVEC_CUDA_CONTENT(
                                                                           Z[0])
                                                                           ->length);
  PostKernelLaunch();

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * OPTIONAL XBraid interface operations
 * -----------------------------------------------------------------
 */

SUNErrCode nvBufSize_Cuda(N_Vector x, sunindextype* size)
{
  if (x == NULL) { return SUN_ERR_GENERIC; }
  *size = (sunindextype)NVEC_CUDA_MEMSIZE(x);
  return SUN_SUCCESS;
}

SUNErrCode nvBufPack_Cuda(N_Vector x, void* buf)
{
  int copy_fail = 0;
  cudaError_t cuerr;

  if (x == NULL || buf == NULL) { return SUN_ERR_GENERIC; }

  SUNMemory buf_mem = SUNMemoryHelper_Wrap(NVEC_CUDA_MEMHELP(x), buf,
                                           SUNMEMTYPE_HOST);
  if (buf_mem == NULL) { return SUN_ERR_GENERIC; }

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_CUDA_MEMHELP(x), buf_mem,
                                        NVEC_CUDA_CONTENT(x)->device_data,
                                        NVEC_CUDA_MEMSIZE(x),
                                        (void*)NVEC_CUDA_STREAM(x));

  /* we synchronize with respect to the host, but only in this stream */
  cuerr = cudaStreamSynchronize(*NVEC_CUDA_STREAM(x));

  SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(x), buf_mem,
                          (void*)NVEC_CUDA_STREAM(x));

  if (!SUNDIALS_CUDA_VERIFY(cuerr) || copy_fail) { return SUN_ERR_GENERIC; }
  else { return SUN_SUCCESS; }
}

SUNErrCode nvBufUnpack_Cuda(N_Vector x, void* buf)
{
  int copy_fail = 0;
  cudaError_t cuerr;

  if (x == NULL || buf == NULL) { return SUN_ERR_GENERIC; }

  SUNMemory buf_mem = SUNMemoryHelper_Wrap(NVEC_CUDA_MEMHELP(x), buf,
                                           SUNMEMTYPE_HOST);
  if (buf_mem == NULL) { return SUN_ERR_GENERIC; }

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_CUDA_MEMHELP(x),
                                        NVEC_CUDA_CONTENT(x)->device_data,
                                        buf_mem, NVEC_CUDA_MEMSIZE(x),
                                        (void*)NVEC_CUDA_STREAM(x));

  /* we synchronize with respect to the host, but only in this stream */
  cuerr = cudaStreamSynchronize(*NVEC_CUDA_STREAM(x));

  SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(x), buf_mem,
                          (void*)NVEC_CUDA_STREAM(x));

  if (!SUNDIALS_CUDA_VERIFY(cuerr) || copy_fail) { return SUN_ERR_GENERIC; }
  else { return SUN_SUCCESS; }
}

/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

SUNErrCode N_VEnableFusedOps_Cuda(N_Vector v, sunbooleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) { return SUN_ERR_GENERIC; }

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }

  if (tf)
  {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = nvLinearCombination_Cuda;
    v->ops->nvscaleaddmulti     = nvScaleAddMulti_Cuda;
    v->ops->nvdotprodmulti      = nvDotProdMulti_Cuda;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray     = nvLinearSumVectorArray_Cuda;
    v->ops->nvscalevectorarray         = nvScaleVectorArray_Cuda;
    v->ops->nvconstvectorarray         = nvConstVectorArray_Cuda;
    v->ops->nvwrmsnormvectorarray      = nvWrmsNormVectorArray_Cuda;
    v->ops->nvwrmsnormmaskvectorarray  = nvWrmsNormMaskVectorArray_Cuda;
    v->ops->nvscaleaddmultivectorarray = nvScaleAddMultiVectorArray_Cuda;
    v->ops->nvlinearcombinationvectorarray = nvLinearCombinationVectorArray_Cuda;
    /* enable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = nvDotProdMulti_Cuda;
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

SUNErrCode N_VEnableLinearCombination_Cuda(N_Vector v, sunbooleantype tf)
{
  if (v == NULL) { return SUN_ERR_GENERIC; }
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }
  v->ops->nvlinearcombination = tf ? nvLinearCombination_Cuda : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleAddMulti_Cuda(N_Vector v, sunbooleantype tf)
{
  if (v == NULL) { return SUN_ERR_GENERIC; }
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }
  v->ops->nvscaleaddmulti = tf ? nvScaleAddMulti_Cuda : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableDotProdMulti_Cuda(N_Vector v, sunbooleantype tf)
{
  if (v == NULL) { return SUN_ERR_GENERIC; }
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }
  v->ops->nvdotprodmulti      = tf ? nvDotProdMulti_Cuda : NULL;
  v->ops->nvdotprodmultilocal = tf ? nvDotProdMulti_Cuda : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearSumVectorArray_Cuda(N_Vector v, sunbooleantype tf)
{
  if (v == NULL) { return SUN_ERR_GENERIC; }
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }
  v->ops->nvlinearsumvectorarray = tf ? nvLinearSumVectorArray_Cuda : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleVectorArray_Cuda(N_Vector v, sunbooleantype tf)
{
  if (v == NULL) { return SUN_ERR_GENERIC; }
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }
  v->ops->nvscalevectorarray = tf ? nvScaleVectorArray_Cuda : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableConstVectorArray_Cuda(N_Vector v, sunbooleantype tf)
{
  if (v == NULL) { return SUN_ERR_GENERIC; }
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }
  v->ops->nvconstvectorarray = tf ? nvConstVectorArray_Cuda : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableWrmsNormVectorArray_Cuda(N_Vector v, sunbooleantype tf)
{
  if (v == NULL) { return SUN_ERR_GENERIC; }
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }
  v->ops->nvwrmsnormvectorarray = tf ? nvWrmsNormVectorArray_Cuda : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableWrmsNormMaskVectorArray_Cuda(N_Vector v, sunbooleantype tf)
{
  if (v == NULL) { return SUN_ERR_GENERIC; }
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }
  v->ops->nvwrmsnormmaskvectorarray = tf ? nvWrmsNormMaskVectorArray_Cuda : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleAddMultiVectorArray_Cuda(N_Vector v, sunbooleantype tf)
{
  if (v == NULL) { return SUN_ERR_GENERIC; }
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }
  v->ops->nvscaleaddmultivectorarray = tf ? nvScaleAddMultiVectorArray_Cuda
                                          : NULL;
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearCombinationVectorArray_Cuda(N_Vector v,
                                                      sunbooleantype tf)
{
  if (v == NULL) { return SUN_ERR_GENERIC; }
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }
  v->ops->nvlinearcombinationvectorarray =
    tf ? nvLinearCombinationVectorArray_Cuda : NULL;
  return SUN_SUCCESS;
}

} // extern "C"

/*
 * Private helper functions.
 */

static int AllocateData(N_Vector v)
{
  int alloc_fail                  = 0;
  N_VectorContent_Cuda vc         = NVEC_CUDA_CONTENT(v);
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  if (nvGetLength_Cuda(v) == 0) { return SUN_SUCCESS; }

  if (vcp->use_managed_mem)
  {
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v), &(vc->device_data),
                                       NVEC_CUDA_MEMSIZE(v), SUNMEMTYPE_UVM,
                                       (void*)NVEC_CUDA_STREAM(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc "
                           "failed for SUNMEMTYPE_UVM\n");
    }
    vc->host_data = SUNMemoryHelper_Alias(NVEC_CUDA_MEMHELP(v), vc->device_data);
  }
  else
  {
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v), &(vc->host_data),
                                       NVEC_CUDA_MEMSIZE(v), SUNMEMTYPE_HOST,
                                       (void*)NVEC_CUDA_STREAM(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc "
                           "failed to alloc SUNMEMTYPE_HOST\n");
    }

    alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v), &(vc->device_data),
                                       NVEC_CUDA_MEMSIZE(v), SUNMEMTYPE_DEVICE,
                                       (void*)NVEC_CUDA_STREAM(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc "
                           "failed to alloc SUNMEMTYPE_DEVICE\n");
    }
  }

  return (alloc_fail ? SUN_ERR_GENERIC : SUN_SUCCESS);
}

/*
 * Initializes the internal buffer used for reductions.
 * If the buffer is already allocated, it will only be reallocated
 * if it is no longer large enough. This may occur if the length
 * of the vector is increased. The buffer is initialized to the
 * value given.
 */
static int InitializeReductionBuffer(N_Vector v, sunrealtype value, size_t n)
{
  int alloc_fail           = 0;
  int copy_fail            = 0;
  sunbooleantype alloc_mem = SUNFALSE;
  size_t bytes             = n * sizeof(sunrealtype);

  // Get the vector private memory structure
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  // Check if the existing reduction memory is not large enough
  if (vcp->reduce_buffer_bytes < bytes)
  {
    FreeReductionBuffer(v);
    alloc_mem = SUNTRUE;
  }

  if (alloc_mem)
  {
    // Allocate pinned memory on the host
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v),
                                       &(vcp->reduce_buffer_host), bytes,
                                       SUNMEMTYPE_PINNED,
                                       (void*)NVEC_CUDA_STREAM(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT(
        "WARNING in InitializeReductionBuffer: SUNMemoryHelper_Alloc failed to "
        "alloc SUNMEMTYPE_PINNED, using SUNMEMTYPE_HOST instead\n");

      // If pinned alloc failed, allocate plain host memory
      alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v),
                                         &(vcp->reduce_buffer_host), bytes,
                                         SUNMEMTYPE_HOST,
                                         (void*)NVEC_CUDA_STREAM(v));
      if (alloc_fail)
      {
        SUNDIALS_DEBUG_PRINT(
          "ERROR in InitializeReductionBuffer: SUNMemoryHelper_Alloc failed to "
          "alloc SUNMEMTYPE_HOST\n");
      }
    }

    // Allocate device memory
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v),
                                       &(vcp->reduce_buffer_dev), bytes,
                                       SUNMEMTYPE_DEVICE,
                                       (void*)NVEC_CUDA_STREAM(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT(
        "ERROR in InitializeReductionBuffer: SUNMemoryHelper_Alloc failed to "
        "alloc SUNMEMTYPE_DEVICE\n");
    }
  }

  if (!alloc_fail)
  {
    // Store the size of the reduction memory buffer
    vcp->reduce_buffer_bytes = bytes;

    // Initialize the host memory with the value
    for (int i = 0; i < n; ++i)
    {
      ((sunrealtype*)vcp->reduce_buffer_host->ptr)[i] = value;
    }

    // Initialize the device memory with the value
    copy_fail = SUNMemoryHelper_CopyAsync(NVEC_CUDA_MEMHELP(v),
                                          vcp->reduce_buffer_dev,
                                          vcp->reduce_buffer_host, bytes,
                                          (void*)NVEC_CUDA_STREAM(v));

    if (copy_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in InitializeReductionBuffer: "
                           "SUNMemoryHelper_CopyAsync failed\n");
    }
  }

  return ((alloc_fail || copy_fail) ? SUN_ERR_GENERIC : SUN_SUCCESS);
}

/* Free the reduction buffer
 */
static void FreeReductionBuffer(N_Vector v)
{
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  if (vcp == NULL) { return; }

  // Free device mem
  if (vcp->reduce_buffer_dev != NULL)
  {
    SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v), vcp->reduce_buffer_dev,
                            (void*)NVEC_CUDA_STREAM(v));
  }
  vcp->reduce_buffer_dev = NULL;

  // Free host mem
  if (vcp->reduce_buffer_host != NULL)
  {
    SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v), vcp->reduce_buffer_host,
                            (void*)NVEC_CUDA_STREAM(v));
  }
  vcp->reduce_buffer_host = NULL;

  // Reset allocated memory size
  vcp->reduce_buffer_bytes = 0;
}

/* Copy the reduction buffer from the device to the host.
 */
static int CopyReductionBufferFromDevice(N_Vector v, size_t n)
{
  int copy_fail;
  cudaError_t cuerr;

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_CUDA_MEMHELP(v),
                                        NVEC_CUDA_PRIVATE(v)->reduce_buffer_host,
                                        NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev,
                                        n * sizeof(sunrealtype),
                                        (void*)NVEC_CUDA_STREAM(v));

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in CopyReductionBufferFromDevice: "
                         "SUNMemoryHelper_CopyAsync returned nonzero\n");
  }

  /* we synchronize with respect to the host, but only in this stream */
  cuerr = cudaStreamSynchronize(*NVEC_CUDA_STREAM(v));
  if (!SUNDIALS_CUDA_VERIFY(cuerr) || copy_fail) { return SUN_ERR_GENERIC; }
  else { return SUN_SUCCESS; }
}

static int FusedBuffer_Init(N_Vector v, int nreal, int nptr)
{
  int alloc_fail           = 0;
  sunbooleantype alloc_mem = SUNFALSE;

  // pad buffer with single precision data
#if defined(SUNDIALS_SINGLE_PRECISION)
  size_t bytes = nreal * 2 * sizeof(sunrealtype) + nptr * sizeof(sunrealtype*);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  size_t bytes = nreal * sizeof(sunrealtype) + nptr * sizeof(sunrealtype*);
#else
#error Incompatible precision for CUDA
#endif

  // Get the vector private memory structure
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  // Check if the existing memory is not large enough
  if (vcp->fused_buffer_bytes < bytes)
  {
    FusedBuffer_Free(v);
    alloc_mem = SUNTRUE;
  }

  if (alloc_mem)
  {
    // Allocate pinned memory on the host
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v),
                                       &(vcp->fused_buffer_host), bytes,
                                       SUNMEMTYPE_PINNED,
                                       (void*)NVEC_CUDA_STREAM(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT(
        "WARNING in FusedBuffer_Init: SUNMemoryHelper_Alloc failed to alloc "
        "SUNMEMTYPE_PINNED, using SUNMEMTYPE_HOST instead\n");

      // If pinned alloc failed, allocate plain host memory
      alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v),
                                         &(vcp->fused_buffer_host), bytes,
                                         SUNMEMTYPE_HOST,
                                         (void*)NVEC_CUDA_STREAM(v));
      if (alloc_fail)
      {
        SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_Init: SUNMemoryHelper_Alloc "
                             "failed to alloc SUNMEMTYPE_HOST\n");
        return SUN_ERR_GENERIC;
      }
    }

    // Allocate device memory
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v),
                                       &(vcp->fused_buffer_dev), bytes,
                                       SUNMEMTYPE_DEVICE,
                                       (void*)NVEC_CUDA_STREAM(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_Init: SUNMemoryHelper_Alloc "
                           "failed to alloc SUNMEMTYPE_DEVICE\n");
      return SUN_ERR_GENERIC;
    }

    // Store the size of the fused op buffer
    vcp->fused_buffer_bytes = bytes;
  }

  // Reset the buffer offset
  vcp->fused_buffer_offset = 0;

  return SUN_SUCCESS;
}

static int FusedBuffer_CopyRealArray(N_Vector v, sunrealtype* rdata, int nval,
                                     sunrealtype** shortcut)
{
  // Get the vector private memory structure
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  // Check buffer space and fill the host buffer
  if (vcp->fused_buffer_offset >= vcp->fused_buffer_bytes)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_CopyRealArray: Buffer offset is "
                         "exceedes the buffer size\n");
    return SUN_ERR_GENERIC;
  }

  sunrealtype* h_buffer = (sunrealtype*)((char*)(vcp->fused_buffer_host->ptr) +
                                         vcp->fused_buffer_offset);

  for (int j = 0; j < nval; j++) { h_buffer[j] = rdata[j]; }

  // Set shortcut to the device buffer and update offset
  *shortcut = (sunrealtype*)((char*)(vcp->fused_buffer_dev->ptr) +
                             vcp->fused_buffer_offset);

  // accounting for buffer padding
#if defined(SUNDIALS_SINGLE_PRECISION)
  vcp->fused_buffer_offset += nval * 2 * sizeof(sunrealtype);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  vcp->fused_buffer_offset += nval * sizeof(sunrealtype);
#else
#error Incompatible precision for CUDA
#endif

  return SUN_SUCCESS;
}

static int FusedBuffer_CopyPtrArray1D(N_Vector v, N_Vector* X, int nvec,
                                      sunrealtype*** shortcut)
{
  // Get the vector private memory structure
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  // Check buffer space and fill the host buffer
  if (vcp->fused_buffer_offset >= vcp->fused_buffer_bytes)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_CopyPtrArray1D: Buffer offset "
                         "is exceedes the buffer size\n");
    return SUN_ERR_GENERIC;
  }

  sunrealtype** h_buffer = (sunrealtype**)((char*)(vcp->fused_buffer_host->ptr) +
                                           vcp->fused_buffer_offset);

  for (int j = 0; j < nvec; j++) { h_buffer[j] = NVEC_CUDA_DDATAp(X[j]); }

  // Set shortcut to the device buffer and update offset
  *shortcut = (sunrealtype**)((char*)(vcp->fused_buffer_dev->ptr) +
                              vcp->fused_buffer_offset);

  vcp->fused_buffer_offset += nvec * sizeof(sunrealtype*);

  return SUN_SUCCESS;
}

static int FusedBuffer_CopyPtrArray2D(N_Vector v, N_Vector** X, int nvec,
                                      int nsum, sunrealtype*** shortcut)
{
  // Get the vector private memory structure
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  // Check buffer space and fill the host buffer
  if (vcp->fused_buffer_offset >= vcp->fused_buffer_bytes)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_CopyPtrArray2D: Buffer offset "
                         "is exceedes the buffer size\n");
    return SUN_ERR_GENERIC;
  }

  sunrealtype** h_buffer = (sunrealtype**)((char*)(vcp->fused_buffer_host->ptr) +
                                           vcp->fused_buffer_offset);

  for (int j = 0; j < nvec; j++)
  {
    for (int k = 0; k < nsum; k++)
    {
      h_buffer[j * nsum + k] = NVEC_CUDA_DDATAp(X[k][j]);
    }
  }

  // Set shortcut to the device buffer and update offset
  *shortcut = (sunrealtype**)((char*)(vcp->fused_buffer_dev->ptr) +
                              vcp->fused_buffer_offset);

  // Update the offset
  vcp->fused_buffer_offset += nvec * nsum * sizeof(sunrealtype*);

  return SUN_SUCCESS;
}

static int FusedBuffer_CopyToDevice(N_Vector v)
{
  // Get the vector private memory structure
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  // Copy the fused buffer to the device
  int copy_fail =
    SUNMemoryHelper_CopyAsync(NVEC_CUDA_MEMHELP(v), vcp->fused_buffer_dev,
                              vcp->fused_buffer_host, vcp->fused_buffer_offset,
                              (void*)NVEC_CUDA_STREAM(v));
  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT(
      "ERROR in FusedBuffer_CopyToDevice: SUNMemoryHelper_CopyAsync failed\n");
    return SUN_ERR_GENERIC;
  }

  // Synchronize with respect to the host, but only in this stream
  SUNDIALS_CUDA_VERIFY(cudaStreamSynchronize(*NVEC_CUDA_STREAM(v)));

  return SUN_SUCCESS;
}

static int FusedBuffer_Free(N_Vector v)
{
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  if (vcp == NULL) { return SUN_SUCCESS; }

  if (vcp->fused_buffer_host)
  {
    SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v), vcp->fused_buffer_host,
                            (void*)NVEC_CUDA_STREAM(v));
    vcp->fused_buffer_host = NULL;
  }

  if (vcp->fused_buffer_dev)
  {
    SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v), vcp->fused_buffer_dev,
                            (void*)NVEC_CUDA_STREAM(v));
    vcp->fused_buffer_dev = NULL;
  }

  vcp->fused_buffer_bytes  = 0;
  vcp->fused_buffer_offset = 0;

  return SUN_SUCCESS;
}

static int InitializeDeviceCounter(N_Vector v)
{
  int retval = 0;
  if (NVEC_CUDA_PRIVATE(v)->device_counter == NULL)
  {
    retval = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v),
                                   &(NVEC_CUDA_PRIVATE(v)->device_counter),
                                   sizeof(unsigned int), SUNMEMTYPE_DEVICE,
                                   (void*)NVEC_CUDA_STREAM(v));
  }
  cudaMemsetAsync(NVEC_CUDA_DCOUNTERp(v), 0, sizeof(unsigned int),
                  *NVEC_CUDA_STREAM(v));
  return retval;
}

static int FreeDeviceCounter(N_Vector v)
{
  int retval = 0;
  if (NVEC_CUDA_PRIVATE(v)->device_counter)
  {
    retval = SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v),
                                     NVEC_CUDA_PRIVATE(v)->device_counter,
                                     (void*)NVEC_CUDA_STREAM(v));
  }
  return retval;
}

/* Get the kernel launch parameters based on the kernel type (reduction or not),
 * using the appropriate kernel execution policy.
 */
static int GetKernelParameters(N_Vector v, sunbooleantype reduction,
                               size_t& grid, size_t& block, size_t& shMemSize,
                               cudaStream_t& stream, bool& atomic, size_t n)
{
  n = (n == 0) ? NVEC_CUDA_CONTENT(v)->length : n;
  if (reduction)
  {
    SUNCudaExecPolicy* reduce_exec_policy =
      NVEC_CUDA_CONTENT(v)->reduce_exec_policy;
    grid      = reduce_exec_policy->gridSize(n);
    block     = reduce_exec_policy->blockSize();
    shMemSize = 0;
    stream    = *(reduce_exec_policy->stream());
    atomic    = reduce_exec_policy->atomic();

    if (!atomic)
    {
      if (InitializeDeviceCounter(v))
      {
#ifdef SUNDIALS_DEBUG
        throw std::runtime_error("SUNMemoryHelper_Alloc returned nonzero\n");
#endif
        return SUN_ERR_GENERIC;
      }
    }

    if (block % sundials::cuda::WARP_SIZE)
    {
#ifdef SUNDIALS_DEBUG
      throw std::runtime_error(
        "the block size must be a multiple must be of the CUDA warp size");
#endif
      return SUN_ERR_GENERIC;
    }
  }
  else
  {
    SUNCudaExecPolicy* stream_exec_policy =
      NVEC_CUDA_CONTENT(v)->stream_exec_policy;
    grid      = stream_exec_policy->gridSize(n);
    block     = stream_exec_policy->blockSize();
    shMemSize = 0;
    stream    = *(stream_exec_policy->stream());
    atomic    = false;
  }

  if (grid == 0)
  {
#ifdef SUNDIALS_DEBUG
    throw std::runtime_error("the grid size must be > 0");
#endif
    return SUN_ERR_GENERIC;
  }
  if (block == 0)
  {
#ifdef SUNDIALS_DEBUG
    throw std::runtime_error("the block size must be > 0");
#endif
    return SUN_ERR_GENERIC;
  }

  return SUN_SUCCESS;
}

static int GetKernelParameters(N_Vector v, sunbooleantype reduction,
                               size_t& grid, size_t& block, size_t& shMemSize,
                               cudaStream_t& stream, size_t n)
{
  bool atomic;
  return GetKernelParameters(v, reduction, grid, block, shMemSize, stream,
                             atomic, n);
}

/* Should be called after a kernel launch.
 * If SUNDIALS_DEBUG_CUDA_LASTERROR is not defined, then the function does nothing.
 * If it is defined, the function will synchronize and check the last CUDA error.
 */
static void PostKernelLaunch()
{
#ifdef SUNDIALS_DEBUG_CUDA_LASTERROR
  cudaDeviceSynchronize();
  SUNDIALS_CUDA_VERIFY(cudaGetLastError());
#endif
}

/* Deprecated concrete operation wrappers */

extern "C" {

void N_VAbs_Cuda(N_Vector x, N_Vector z) { nvAbs_Cuda(x, z); }

void N_VAddConst_Cuda(N_Vector x, sunrealtype b, N_Vector z)
{
  nvAddConst_Cuda(x, b, z);
}

SUNErrCode N_VBufPack_Cuda(N_Vector x, void* buf)
{
  return nvBufPack_Cuda(x, buf);
}

SUNErrCode N_VBufSize_Cuda(N_Vector x, sunindextype* size)
{
  return nvBufSize_Cuda(x, size);
}

SUNErrCode N_VBufUnpack_Cuda(N_Vector x, void* buf)
{
  return nvBufUnpack_Cuda(x, buf);
}

N_Vector N_VCloneEmpty_Cuda(N_Vector w) { return nvCloneEmpty_Cuda(w); }

N_Vector N_VClone_Cuda(N_Vector w) { return nvClone_Cuda(w); }

void N_VCompare_Cuda(sunrealtype c, N_Vector x, N_Vector z)
{
  nvCompare_Cuda(c, x, z);
}

SUNErrCode N_VConstVectorArray_Cuda(int nvec, sunrealtype c, N_Vector* Z)
{
  return nvConstVectorArray_Cuda(nvec, c, Z);
}

void N_VConst_Cuda(sunrealtype c, N_Vector z) { nvConst_Cuda(c, z); }

sunbooleantype N_VConstrMask_Cuda(N_Vector c, N_Vector x, N_Vector m)
{
  return nvConstrMask_Cuda(c, x, m);
}

void N_VDestroy_Cuda(N_Vector v) { nvDestroy_Cuda(v); }

void N_VDiv_Cuda(N_Vector x, N_Vector y, N_Vector z) { nvDiv_Cuda(x, y, z); }

SUNErrCode N_VDotProdMulti_Cuda(int nvec, N_Vector x, N_Vector* Y,
                                sunrealtype* dotprods)
{
  return nvDotProdMulti_Cuda(nvec, x, Y, dotprods);
}

sunrealtype N_VDotProd_Cuda(N_Vector x, N_Vector y)
{
  return nvDotProd_Cuda(x, y);
}

sunbooleantype N_VInvTest_Cuda(N_Vector x, N_Vector z)
{
  return nvInvTest_Cuda(x, z);
}

void N_VInv_Cuda(N_Vector x, N_Vector z) { nvInv_Cuda(x, z); }

sunrealtype N_VL1Norm_Cuda(N_Vector x) { return nvL1Norm_Cuda(x); }

SUNErrCode N_VLinearCombinationVectorArray_Cuda(int nvec, int nsum,
                                                sunrealtype* c, N_Vector** X,
                                                N_Vector* Z)
{
  return nvLinearCombinationVectorArray_Cuda(nvec, nsum, c, X, Z);
}

SUNErrCode N_VLinearCombination_Cuda(int nvec, sunrealtype* c, N_Vector* X,
                                     N_Vector Z)
{
  return nvLinearCombination_Cuda(nvec, c, X, Z);
}

SUNErrCode N_VLinearSumVectorArray_Cuda(int nvec, sunrealtype a, N_Vector* X,
                                        sunrealtype b, N_Vector* Y, N_Vector* Z)
{
  return nvLinearSumVectorArray_Cuda(nvec, a, X, b, Y, Z);
}

void N_VLinearSum_Cuda(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                       N_Vector z)
{
  nvLinearSum_Cuda(a, x, b, y, z);
}

sunrealtype N_VMaxNorm_Cuda(N_Vector x) { return nvMaxNorm_Cuda(x); }

sunrealtype N_VMinQuotient_Cuda(N_Vector num, N_Vector denom)
{
  return nvMinQuotient_Cuda(num, denom);
}

sunrealtype N_VMin_Cuda(N_Vector x) { return nvMin_Cuda(x); }

void N_VPrintFile_Cuda(N_Vector v, FILE* outfile)
{
  nvPrintFile_Cuda(v, outfile);
}

void N_VPrint_Cuda(N_Vector v) { nvPrint_Cuda(v); }

void N_VProd_Cuda(N_Vector x, N_Vector y, N_Vector z) { nvProd_Cuda(x, y, z); }

SUNErrCode N_VScaleAddMultiVectorArray_Cuda(int nvec, int nsum, sunrealtype* a,
                                            N_Vector* X, N_Vector** Y,
                                            N_Vector** Z)
{
  return nvScaleAddMultiVectorArray_Cuda(nvec, nsum, a, X, Y, Z);
}

SUNErrCode N_VScaleAddMulti_Cuda(int nvec, sunrealtype* c, N_Vector X,
                                 N_Vector* Y, N_Vector* Z)
{
  return nvScaleAddMulti_Cuda(nvec, c, X, Y, Z);
}

SUNErrCode N_VScaleVectorArray_Cuda(int nvec, sunrealtype* c, N_Vector* X,
                                    N_Vector* Z)
{
  return nvScaleVectorArray_Cuda(nvec, c, X, Z);
}

void N_VScale_Cuda(sunrealtype c, N_Vector x, N_Vector z)
{
  nvScale_Cuda(c, x, z);
}

void N_VSetDeviceArrayPointer_Cuda(sunrealtype* d_vdata_1d, N_Vector v)
{
  nvSetDeviceArrayPointer_Cuda(d_vdata_1d, v);
}

void N_VSetHostArrayPointer_Cuda(sunrealtype* h_vdata_1d, N_Vector v)
{
  nvSetHostArrayPointer_Cuda(h_vdata_1d, v);
}

sunrealtype N_VWL2Norm_Cuda(N_Vector x, N_Vector w)
{
  return nvWL2Norm_Cuda(x, w);
}

sunrealtype N_VWSqrSumLocal_Cuda(N_Vector x, N_Vector w)
{
  return nvWSqrSumLocal_Cuda(x, w);
}

sunrealtype N_VWSqrSumMaskLocal_Cuda(N_Vector x, N_Vector w, N_Vector id)
{
  return nvWSqrSumMaskLocal_Cuda(x, w, id);
}

SUNErrCode N_VWrmsNormMaskVectorArray_Cuda(int nvec, N_Vector* X, N_Vector* W,
                                           N_Vector id, sunrealtype* nrm)
{
  return nvWrmsNormMaskVectorArray_Cuda(nvec, X, W, id, nrm);
}

sunrealtype N_VWrmsNormMask_Cuda(N_Vector x, N_Vector w, N_Vector id)
{
  return nvWrmsNormMask_Cuda(x, w, id);
}

SUNErrCode N_VWrmsNormVectorArray_Cuda(int nvec, N_Vector* X, N_Vector* W,
                                       sunrealtype* nrm)
{
  return nvWrmsNormVectorArray_Cuda(nvec, X, W, nrm);
}

sunrealtype N_VWrmsNorm_Cuda(N_Vector x, N_Vector w)
{
  return nvWrmsNorm_Cuda(x, w);
}

} // extern "C"
