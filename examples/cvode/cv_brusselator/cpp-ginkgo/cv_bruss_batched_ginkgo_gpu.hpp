/* ----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * ----------------------------------------------------------------------------
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
 * ----------------------------------------------------------------------------*/

#pragma once

#include <cstdio>

#include <sunlinsol/sunlinsol_ginkgobatch.hpp>

#include "cv_bruss_batched_ginkgo.hpp"

#if defined(USE_CUDA)
#include <nvector/nvector_cuda.h>
#include <sunmemory/sunmemory_cuda.h>
#define GPU_DEVICE_SYNCHRONIZE cudaDeviceSynchronize
#define GPU_GET_LAST_ERROR     cudaGetLastError
#define GPU_SUCCESS            cudaSuccess
#define GPU_GET_ERROR_NAME     cudaGetErrorName
#elif defined(USE_HIP)
#include <nvector/nvector_hip.h>
#include <sunmemory/sunmemory_hip.h>
#define GPU_DEVICE_SYNCHRONIZE hipDeviceSynchronize
#define GPU_GET_LAST_ERROR     hipGetLastError
#define GPU_SUCCESS            hipSuccess
#define GPU_GET_ERROR_NAME     hipGetErrorName
#endif

using GkoBatchMatrixTypeGPU = gko::batch::matrix::Csr<sunrealtype, sunindextype>;
using SUNGkoMatrixTypeGPU = sundials::ginkgo::BatchMatrix<GkoBatchMatrixTypeGPU>;

N_Vector create_vector(sunindextype length, SUNContext sunctx)
{
#if defined(USE_CUDA)
  return N_VNewManaged_Cuda(length, sunctx);
#else
  return N_VNewManaged_Hip(length, sunctx);
#endif
}

SUNMemoryHelper create_memory_helper(SUNContext sunctx)
{
#if defined(USE_CUDA)
  return SUNMemoryHelper_Cuda(sunctx);
#else
  return SUNMemoryHelper_Hip(sunctx);
#endif
}

/* Right hand side function evaluation GPU kernel.
   This kernel needs total number of threads >= num_batches. */
__global__ void f_kernel(sunrealtype t, sunrealtype* ydata, sunrealtype* ydotdata,
                         sunrealtype* A, sunrealtype* B, sunrealtype* Ep,
                         int neq, int num_batches, int batch_size)
{
  sunrealtype u, v, w, a, b, ep;

  int batchj = blockIdx.x * blockDim.x + threadIdx.x;

  if (batchj < num_batches)
  {
    a = A[batchj];
    b = B[batchj], ep = Ep[batchj];

    u = ydata[batchj * batch_size];
    v = ydata[batchj * batch_size + 1];
    w = ydata[batchj * batch_size + 2];

    ydotdata[batchj * batch_size]     = a - (w + 1.0) * u + v * u * u;
    ydotdata[batchj * batch_size + 1] = w * u - v * u * u;
    ydotdata[batchj * batch_size + 2] = (b - w) / ep - w * u;
  }
}

/* Jacobian evaluation GPU kernel
   This kernel needs total number of threads >= num_batches. */
__global__ void j_kernel(sunrealtype* ydata, sunrealtype* Jdata, sunrealtype* A,
                         sunrealtype* B, sunrealtype* Ep, int neq,
                         int num_batches, int batch_size, int nnzper)
{
  sunrealtype u, v, w, ep;

  int batchj = blockIdx.x * blockDim.x + threadIdx.x;

  if (batchj < num_batches)
  {
    ep = Ep[batchj];

    /* get y values */
    u = ydata[batch_size * batchj];
    v = ydata[batch_size * batchj + 1];
    w = ydata[batch_size * batchj + 2];

    /* first row of batch */
    Jdata[nnzper * batchj]     = -(w + 1.0) + 2.0 * u * v;
    Jdata[nnzper * batchj + 1] = u * u;
    Jdata[nnzper * batchj + 2] = -u;

    /* second row of batch */
    Jdata[nnzper * batchj + 3] = w - 2.0 * u * v;
    Jdata[nnzper * batchj + 4] = -u * u;
    Jdata[nnzper * batchj + 5] = u;

    /* third row of batch */
    Jdata[nnzper * batchj + 6] = -w;
    Jdata[nnzper * batchj + 7] = 0.0;
    Jdata[nnzper * batchj + 8] = -1.0 / ep - u;
  }
}

/* Jacobian-vector product GPU kernel
   This kernel needs total number of threads >= num_batches. */
__global__ void jv_kernel(sunrealtype* vdata, sunrealtype* Jvdata,
                          sunrealtype* ydata, sunrealtype* A, sunrealtype* B,
                          sunrealtype* Ep, int neq, int num_batches,
                          int batch_size, int nnzper)
{
  sunrealtype u, v, w, v0, v1, v2, ep;

  int batchj = blockIdx.x * blockDim.x + threadIdx.x;

  if (batchj < num_batches)
  {
    ep = Ep[batchj];

    /* get y values */
    u = ydata[batch_size * batchj];
    v = ydata[batch_size * batchj + 1];
    w = ydata[batch_size * batchj + 2];

    /* get v values */
    v0 = vdata[batch_size * batchj];
    v1 = vdata[batch_size * batchj + 1];
    v2 = vdata[batch_size * batchj + 2];

    /* initialize Jv to zero */
    Jvdata[batch_size * batchj]     = 0.0;
    Jvdata[batch_size * batchj + 1] = 0.0;
    Jvdata[batch_size * batchj + 2] = 0.0;

    // add J[:,0]*v[0]
    Jvdata[batch_size * batchj] += (-(w + 1.0) + 2.0 * u * v) * v0;
    Jvdata[batch_size * batchj + 1] += (w - 2.0 * u * v) * v0;
    Jvdata[batch_size * batchj + 2] += -w * v0;

    // add J[:,1]*v[1]
    Jvdata[batch_size * batchj] += (u * u) * v1;
    Jvdata[batch_size * batchj + 1] += (-u * u) * v1;
    Jvdata[batch_size * batchj + 2] += 0.0;

    // add J[:,2]*v[2]
    Jvdata[batch_size * batchj] += -u * v2;
    Jvdata[batch_size * batchj + 1] += u * v2;
    Jvdata[batch_size * batchj + 2] += (-1.0 / ep - u) * v2;
  }
}

int f_gpu(sunrealtype t, N_Vector y, N_Vector ydot, UserData* udata)
{
  sunrealtype* ydata    = N_VGetArrayPointer(y);
  sunrealtype* ydotdata = N_VGetArrayPointer(ydot);

  unsigned threads_per_block = 256;
  unsigned num_blocks        = (udata->num_batches + threads_per_block - 1) /
                        threads_per_block;
  f_kernel<<<num_blocks, threads_per_block>>>(t, ydata, ydotdata, udata->a.get(),
                                              udata->b.get(), udata->ep.get(),
                                              udata->neq, udata->num_batches,
                                              udata->batch_size);

  GPU_DEVICE_SYNCHRONIZE();
  auto gpu_err = GPU_GET_LAST_ERROR();
  if (gpu_err != GPU_SUCCESS)
  {
    fprintf(stderr, ">>> ERROR in f: GetLastError returned %s\n",
            GPU_GET_ERROR_NAME(gpu_err));
    return -1;
  }

  return 0;
}

int Jac_gpu(N_Vector y, SUNMatrix J, UserData* udata)
{
  auto Jgko = static_cast<SUNGkoMatrixTypeGPU*>(J->content)->GkoMtx();

  sunrealtype* Jdata = Jgko->get_values();
  sunrealtype* ydata = N_VGetArrayPointer(y);

  unsigned threads_per_block = 256;
  unsigned num_blocks        = (udata->num_batches + threads_per_block - 1) /
                        threads_per_block;
  j_kernel<<<num_blocks, threads_per_block>>>(ydata, Jdata, udata->a.get(),
                                              udata->b.get(), udata->ep.get(),
                                              udata->neq, udata->num_batches,
                                              udata->batch_size, udata->nnzper);

  GPU_DEVICE_SYNCHRONIZE();
  auto gpu_err = GPU_GET_LAST_ERROR();
  if (gpu_err != GPU_SUCCESS)
  {
    fprintf(stderr, ">>> ERROR in Jac: GetLastError returned %s\n",
            GPU_GET_ERROR_NAME(gpu_err));
    return -1;
  }

  return 0;
}

int JacVec_gpu(N_Vector v, N_Vector Jv, N_Vector y, UserData* udata)
{
  sunrealtype* vdata  = N_VGetArrayPointer(v);
  sunrealtype* Jvdata = N_VGetArrayPointer(Jv);
  sunrealtype* ydata  = N_VGetArrayPointer(y);

  unsigned threads_per_block = 256;
  unsigned num_blocks        = (udata->num_batches + threads_per_block - 1) /
                        threads_per_block;
  jv_kernel<<<num_blocks, threads_per_block>>>(vdata, Jvdata, ydata,
                                               udata->a.get(), udata->b.get(),
                                               udata->ep.get(), udata->neq,
                                               udata->num_batches,
                                               udata->batch_size, udata->nnzper);

  GPU_DEVICE_SYNCHRONIZE();
  auto gpu_err = GPU_GET_LAST_ERROR();
  if (gpu_err != GPU_SUCCESS)
  {
    fprintf(stderr, ">>> ERROR in JacVec: GetLastError returned %s\n",
            GPU_GET_ERROR_NAME(gpu_err));
    return -1;
  }

  return 0;
}

#undef GPU_DEVICE_SYNCHRONIZE
#undef GPU_GET_LAST_ERROR
#undef GPU_SUCCESS
#undef GPU_GET_ERROR_NAME
