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
 * CUDA/HIP diffusion functions
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

// Forcing device function
__device__ void add_forcing(const sunrealtype t, const sunrealtype x,
                            const sunrealtype y, const sunrealtype kx,
                            const sunrealtype ky, const sunindextype c,
                            sunrealtype* udot)
{
  sunrealtype sin_sqr_x = sin(PI * x) * sin(PI * x);
  sunrealtype sin_sqr_y = sin(PI * y) * sin(PI * y);

  sunrealtype cos_sqr_x = cos(PI * x) * cos(PI * x);
  sunrealtype cos_sqr_y = cos(PI * y) * cos(PI * y);

  sunrealtype sin_t_cos_t = sin(PI * t) * cos(PI * t);
  sunrealtype cos_sqr_t   = cos(PI * t) * cos(PI * t);

  sunrealtype bx = kx * TWO * PI * PI;
  sunrealtype by = ky * TWO * PI * PI;

  udot[c] += -TWO * PI * sin_sqr_x * sin_sqr_y * sin_t_cos_t -
             bx * (cos_sqr_x - sin_sqr_x) * sin_sqr_y * cos_sqr_t -
             by * (cos_sqr_y - sin_sqr_y) * sin_sqr_x * cos_sqr_t;
}

// Interior diffusion kernel
__global__ void diffusion_interior_kernel(
  const sunrealtype t, const sunrealtype* u, sunrealtype* udot,
  const sunindextype is, const sunindextype js, const sunindextype nx_loc,
  const sunindextype ny_loc, const sunrealtype dx, const sunrealtype dy,
  const sunrealtype kx, const sunrealtype ky, const bool forcing)
{
  // Thread location in the local grid
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  // Only update the interior points
  bool interior = (i > 0 && i < nx_loc - 1 && j > 0 && j < ny_loc - 1);
  if (interior)
  {
    // 1D array index for center, west, east, south, and north nodes
    int c = i + j * nx_loc;
    int w = c - 1;
    int e = c + 1;
    int s = c - nx_loc;
    int n = c + nx_loc;

    // Set diffusion term
    sunrealtype cx = kx / (dx * dx);
    sunrealtype cy = ky / (dy * dy);
    sunrealtype cc = -TWO * (cx + cy);

    udot[c] = cc * u[c] + cx * (u[w] + u[e]) + cy * (u[s] + u[n]);

    if (forcing)
    {
      sunrealtype x = (is + i) * dx;
      sunrealtype y = (js + j) * dy;

      add_forcing(t, x, y, kx, ky, c, udot);
    }
  }
}

// Interior boundary kernel
__global__ void diffusion_boundary_kernel(
  const sunrealtype t, const sunrealtype* u, sunrealtype* udot,
  const sunindextype is, const sunindextype js, const sunindextype nx_loc,
  const sunindextype ny_loc, const sunrealtype dx, const sunrealtype dy,
  const sunrealtype kx, const sunrealtype ky, const bool forcing,
  const sunrealtype* wbuf, const sunrealtype* ebuf, const sunrealtype* sbuf,
  const sunrealtype* nbuf)
{
  // Thread ID
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  // Set diffusion term
  sunrealtype cx = kx / (dx * dx);
  sunrealtype cy = ky / (dy * dy);
  sunrealtype cc = -TWO * (cx + cy);

  // West and East faces excluding corners
  if (i > 0 && i < ny_loc - 1)
  {
    // West face: 1D array index for center, west, east, south, and north nodes
    int c = i * nx_loc;
    int w = i;
    int e = c + 1;
    int s = c - nx_loc;
    int n = c + nx_loc;

    if (wbuf)
    {
      // West processor boundary
      udot[c] = cc * u[c] + cx * (wbuf[w] + u[e]) + cy * (u[s] + u[n]);

      if (forcing)
      {
        sunrealtype x = is * dx;
        sunrealtype y = (js + i) * dy;

        add_forcing(t, x, y, kx, ky, c, udot);
      }
    }
    else
    {
      // West physical boundary
      udot[c] = ZERO;
    }

    // East face: 1D array index for center, west, east, south, and north nodes
    c = (i + 1) * nx_loc - 1;
    w = c - 1;
    e = i;
    s = c - nx_loc;
    n = c + nx_loc;

    if (ebuf)
    {
      // East processor boundary
      udot[c] = cc * u[c] + cx * (u[w] + ebuf[e]) + cy * (u[s] + u[n]);

      if (forcing)
      {
        sunrealtype x = (is + nx_loc - 1) * dx;
        sunrealtype y = (js + i) * dy;

        add_forcing(t, x, y, kx, ky, c, udot);
      }
    }
    else
    {
      // East physical boundary
      udot[c] = ZERO;
    }
  }

  __syncthreads();

  // South and North faces excluding corners
  if (i > 0 && i < nx_loc - 1)
  {
    // South face: 1D array index for center, west, east, south, and north nodes
    int c = i;
    int w = c - 1;
    int e = c + 1;
    int s = i;
    int n = c + nx_loc;

    if (sbuf)
    {
      // South processor boundary
      udot[c] = cc * u[c] + cx * (u[w] + u[e]) + cy * (sbuf[s] + u[n]);

      if (forcing)
      {
        sunrealtype x = (is + i) * dx;
        sunrealtype y = js * dy;

        add_forcing(t, x, y, kx, ky, c, udot);
      }
    }
    else
    {
      // South physical boundary
      udot[c] = ZERO;
    }

    // North face: 1D array index for center, west, east, south, and north nodes
    c = i + (ny_loc - 1) * nx_loc;
    w = c - 1;
    e = c + 1;
    s = c - nx_loc;
    n = i;

    if (nbuf)
    {
      // North processor boundary
      udot[c] = cc * u[c] + cx * (u[w] + u[e]) + cy * (u[s] + nbuf[n]);

      if (forcing)
      {
        sunrealtype x = (is + i) * dx;
        sunrealtype y = (js + ny_loc - 1) * dy;

        add_forcing(t, x, y, kx, ky, c, udot);
      }
    }
    else
    {
      // North physical boundary
      udot[c] = ZERO;
    }
  }

  __syncthreads();

  // Corners
  if (i == 0)
  {
    // South-West
    int c = 0;
    int w = 0;
    int e = c + 1;
    int s = 0;
    int n = c + nx_loc;

    if (sbuf && wbuf)
    {
      // South-West processor boundary
      udot[c] = cc * u[c] + cx * (wbuf[w] + u[e]) + cy * (sbuf[s] + u[n]);

      if (forcing)
      {
        sunrealtype x = is * dx;
        sunrealtype y = js * dy;

        add_forcing(t, x, y, kx, ky, c, udot);
      }
    }
    else
    {
      // South-West physical boundary
      udot[c] = ZERO;
    }

    // South-East
    c = nx_loc - 1;
    w = c - 1;
    e = 0;
    s = nx_loc - 1;
    n = c + nx_loc;

    if (sbuf && ebuf)
    {
      // South-East processor boundary
      udot[c] = cc * u[c] + cx * (u[w] + ebuf[e]) + cy * (sbuf[s] + u[n]);

      if (forcing)
      {
        sunrealtype x = (is + nx_loc - 1) * dx;
        sunrealtype y = js * dy;

        add_forcing(t, x, y, kx, ky, c, udot);
      }
    }
    else
    {
      // South-East physical boundary
      udot[c] = ZERO;
    }

    // North-West
    c = (ny_loc - 1) * nx_loc;
    w = ny_loc - 1;
    e = c + 1;
    s = c - nx_loc;
    n = 0;

    if (nbuf && wbuf)
    {
      // North-West processor boundary
      udot[c] = cc * u[c] + cx * (wbuf[w] + u[e]) + cy * (u[s] + nbuf[n]);

      if (forcing)
      {
        sunrealtype x = is * dx;
        sunrealtype y = (js + ny_loc - 1) * dy;

        add_forcing(t, x, y, kx, ky, c, udot);
      }
    }
    else
    {
      // North-West physical boundary
      udot[c] = ZERO;
    }

    // North-East
    c = nx_loc * ny_loc - 1;
    w = c - 1;
    e = ny_loc - 1;
    s = c - nx_loc;
    n = nx_loc - 1;

    if (nbuf && ebuf)
    {
      // North-East processor boundary
      udot[c] = cc * u[c] + cx * (u[w] + ebuf[e]) + cy * (u[s] + nbuf[n]);

      if (forcing)
      {
        sunrealtype x = (is + nx_loc - 1) * dx;
        sunrealtype y = (js + ny_loc - 1) * dy;

        add_forcing(t, x, y, kx, ky, c, udot);
      }
    }
    else
    {
      // North-East physical boundary
      udot[c] = ZERO;
    }
  }
}

// Diffusion function
int laplacian(sunrealtype t, N_Vector u, N_Vector f, UserData* udata)
{
  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  int flag;

  // Start exchange
  flag = udata->start_exchange(u);
  if (check_flag(&flag, "UserData::start_exchange", 1)) { return -1; }

  // Extract needed constants from user data
  const sunindextype is     = udata->is;
  const sunindextype js     = udata->js;
  const sunindextype nx_loc = udata->nx_loc;
  const sunindextype ny_loc = udata->ny_loc;
  const sunrealtype dx      = udata->dx;
  const sunrealtype dy      = udata->dy;
  const sunrealtype kx      = udata->kx;
  const sunrealtype ky      = udata->ky;
  const bool forcing        = udata->forcing;

  // Access data arrays
  const sunrealtype* uarray =
    N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(u));
  if (check_flag((void*)uarray, "N_VGetDeviceArrayPointer", 0)) { return -1; }

  sunrealtype* farray = N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(f));
  if (check_flag((void*)farray, "N_VGetDeviceArrayPointer", 0)) { return -1; }

  // Update subdomain interior
  dim3 iblock(BLOCK_SIZE_X, BLOCK_SIZE_Y);
  dim3 igrid(ICEIL(nx_loc, BLOCK_SIZE_X), ICEIL(ny_loc, BLOCK_SIZE_Y));

  diffusion_interior_kernel<<<igrid, iblock>>>(t, uarray, farray, is, js, nx_loc,
                                               ny_loc, dx, dy, kx, ky, forcing);

  // Wait for exchange receives
  flag = udata->end_exchange();
  if (check_flag(&flag, "UserData::end_exchagne", 1)) { return -1; }

  // Update subdomain boundary
  const sunrealtype* Warray = (udata->HaveNbrW) ? udata->Wrecv : NULL;
  const sunrealtype* Earray = (udata->HaveNbrE) ? udata->Erecv : NULL;
  const sunrealtype* Sarray = (udata->HaveNbrS) ? udata->Srecv : NULL;
  const sunrealtype* Narray = (udata->HaveNbrN) ? udata->Nrecv : NULL;

  sunindextype maxdim = max(nx_loc, ny_loc);
  dim3 bblock(BLOCK_SIZE);
  dim3 bgrid(ICEIL(maxdim, BLOCK_SIZE));

  diffusion_boundary_kernel<<<bgrid, bblock>>>(t, uarray, farray, is, js, nx_loc,
                                               ny_loc, dx, dy, kx, ky, forcing,
                                               Warray, Earray, Sarray, Narray);

  // Return success
  return 0;
}
