/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 * CUDA/HIP functions for exchange buffers
 * ---------------------------------------------------------------------------*/

#include "diffusion_2D.hpp"

#if defined(USE_HIP)
#include <hip/hip_runtime.h>
#define BLOCK_SIZE   256
#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16
#elif defined(USE_CUDA)
#include <cuda_runtime.h>
#define BLOCK_SIZE   256
#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16
#else
#error Define USE_CUDA or USE_HIP
#endif

namespace
{
int copy_to_host(vector<sunrealtype>& host, const sunrealtype* device)
{
  if (host.empty()) { return 0; }
#if defined(USE_HIP)
  hipError_t flag = hipMemcpy(host.data(), device,
                             host.size() * sizeof(sunrealtype),
                             hipMemcpyDeviceToHost);
  if (flag != hipSuccess)
  {
    cerr << "hipMemcpy device to host failed: " << hipGetErrorString(flag)
         << endl;
    return -1;
  }
#else
  cudaError_t flag = cudaMemcpy(host.data(), device,
                               host.size() * sizeof(sunrealtype),
                               cudaMemcpyDeviceToHost);
  if (flag != cudaSuccess)
  {
    cerr << "cudaMemcpy device to host failed: " << cudaGetErrorString(flag)
         << endl;
    return -1;
  }
#endif
  return 0;
}

int copy_to_device(sunrealtype* device, const vector<sunrealtype>& host)
{
  if (host.empty()) { return 0; }
#if defined(USE_HIP)
  hipError_t flag = hipMemcpy(device, host.data(),
                             host.size() * sizeof(sunrealtype),
                             hipMemcpyHostToDevice);
  if (flag != hipSuccess)
  {
    cerr << "hipMemcpy host to device failed: " << hipGetErrorString(flag)
         << endl;
    return -1;
  }
#else
  cudaError_t flag = cudaMemcpy(device, host.data(),
                               host.size() * sizeof(sunrealtype),
                               cudaMemcpyHostToDevice);
  if (flag != cudaSuccess)
  {
    cerr << "cudaMemcpy host to device failed: " << cudaGetErrorString(flag)
         << endl;
    return -1;
  }
#endif
  return 0;
}
} // namespace

// Pack exchange buffers kernel
__global__ void pack_buffers_kernel(const sunindextype nx_loc,
                                    const sunindextype ny_loc,
                                    const sunrealtype* u, sunrealtype* wbuf,
                                    sunrealtype* ebuf, sunrealtype* sbuf,
                                    sunrealtype* nbuf)
{
  // Thread ID
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  // West and East faces
  if (i < ny_loc)
  {
    if (wbuf) { wbuf[i] = u[i * nx_loc]; }
    if (ebuf) { ebuf[i] = u[(i + 1) * nx_loc - 1]; }
  }

  // South and North faces
  if (i < nx_loc)
  {
    if (sbuf) { sbuf[i] = u[i]; }
    if (nbuf) { nbuf[i] = u[i + (ny_loc - 1) * nx_loc]; }
  }
}

// Pack exchange buffers
int UserData::pack_buffers(const N_Vector u)
{
  // Access data array
  const sunrealtype* uarray =
    N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(u));
  if (check_flag((void*)uarray, "N_VGetDeviceArrayPointer", 0)) return -1;

  sunindextype maxdim = max(nx_loc, ny_loc);
  dim3 block(BLOCK_SIZE);
  dim3 grid(ICEIL(maxdim, BLOCK_SIZE));

  pack_buffers_kernel<<<grid, block>>>(nx_loc, ny_loc, uarray, Wsend, Esend,
                                       Ssend, Nsend);

  return 0;
}

int UserData::stage_send_buffers()
{
  if (copy_to_host(Wsend_host, Wsend) || copy_to_host(Esend_host, Esend) ||
      copy_to_host(Ssend_host, Ssend) || copy_to_host(Nsend_host, Nsend))
  {
    return -1;
  }
  return 0;
}

int UserData::stage_recv_buffers()
{
  if (copy_to_device(Wrecv, Wrecv_host) || copy_to_device(Erecv, Erecv_host) ||
      copy_to_device(Srecv, Srecv_host) || copy_to_device(Nrecv, Nrecv_host))
  {
    return -1;
  }
  return 0;
}

// Allocate exchange buffers
int UserData::allocate_buffers()
{
  if (HaveNbrW)
  {
#if defined(USE_HIP)
    hipMalloc(&Wrecv, ny_loc * sizeof(sunrealtype));
    hipMalloc(&Wsend, ny_loc * sizeof(sunrealtype));
#elif defined(USE_CUDA)
    cudaMalloc(&Wrecv, ny_loc * sizeof(sunrealtype));
    cudaMalloc(&Wsend, ny_loc * sizeof(sunrealtype));
#else
#error Define USE_CUDA or USE_HIP
#endif
  }

  if (HaveNbrE)
  {
#if defined(USE_HIP)
    hipMalloc(&Erecv, ny_loc * sizeof(sunrealtype));
    hipMalloc(&Esend, ny_loc * sizeof(sunrealtype));
#elif defined(USE_CUDA)
    cudaMalloc(&Erecv, ny_loc * sizeof(sunrealtype));
    cudaMalloc(&Esend, ny_loc * sizeof(sunrealtype));
#else
#error Define USE_CUDA or USE_HIP
#endif
  }

  if (HaveNbrS)
  {
#if defined(USE_HIP)
    hipMalloc(&Srecv, nx_loc * sizeof(sunrealtype));
    hipMalloc(&Ssend, nx_loc * sizeof(sunrealtype));
#elif defined(USE_CUDA)
    cudaMalloc(&Srecv, nx_loc * sizeof(sunrealtype));
    cudaMalloc(&Ssend, nx_loc * sizeof(sunrealtype));
#else
#error Define USE_CUDA or USE_HIP
#endif
  }

  if (HaveNbrN)
  {
#if defined(USE_HIP)
    hipMalloc(&Nrecv, nx_loc * sizeof(sunrealtype));
    hipMalloc(&Nsend, nx_loc * sizeof(sunrealtype));
#elif defined(USE_CUDA)
    cudaMalloc(&Nrecv, nx_loc * sizeof(sunrealtype));
    cudaMalloc(&Nsend, nx_loc * sizeof(sunrealtype));
#else
#error Define USE_CUDA or USE_HIP
#endif
  }

  if (host_staging)
  {
    if (HaveNbrW)
    {
      Wrecv_host.resize(ny_loc);
      Wsend_host.resize(ny_loc);
    }
    if (HaveNbrE)
    {
      Erecv_host.resize(ny_loc);
      Esend_host.resize(ny_loc);
    }
    if (HaveNbrS)
    {
      Srecv_host.resize(nx_loc);
      Ssend_host.resize(nx_loc);
    }
    if (HaveNbrN)
    {
      Nrecv_host.resize(nx_loc);
      Nsend_host.resize(nx_loc);
    }
  }

  return 0;
}

// Free exchange buffers
int UserData::free_buffers()
{
#if defined(USE_HIP)
  if (Wrecv != NULL) hipFree(Wrecv);
  if (Wsend != NULL) hipFree(Wsend);
  if (Erecv != NULL) hipFree(Erecv);
  if (Esend != NULL) hipFree(Esend);
  if (Srecv != NULL) hipFree(Srecv);
  if (Ssend != NULL) hipFree(Ssend);
  if (Nrecv != NULL) hipFree(Nrecv);
  if (Nsend != NULL) hipFree(Nsend);
#elif defined(USE_CUDA)
  if (Wrecv != NULL) cudaFree(Wrecv);
  if (Wsend != NULL) cudaFree(Wsend);
  if (Erecv != NULL) cudaFree(Erecv);
  if (Esend != NULL) cudaFree(Esend);
  if (Srecv != NULL) cudaFree(Srecv);
  if (Ssend != NULL) cudaFree(Ssend);
  if (Nrecv != NULL) cudaFree(Nrecv);
  if (Nsend != NULL) cudaFree(Nsend);
#else
#error Define USE_CUDA or USE_HIP
#endif

  return 0;
}
