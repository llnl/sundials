/*
 * -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos, and Daniel McGreer @ LLNL
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
 * This header files defines internal utility functions and macros
 * for working with HIP.
 * -----------------------------------------------------------------
 */

#include <hip/hip_runtime.h>
#include <stdio.h>
#include <sundials/sundials_types.h>

#ifndef _SUNDIALS_HIP_H
#define _SUNDIALS_HIP_H

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ---------------------------------------------------------------------------
 * Utility macros
 * ---------------------------------------------------------------------------*/

#define SUNDIALS_HIP_VERIFY(hiperr) \
  SUNDIALS_HIP_Assert(hiperr, __FILE__, __LINE__)

#define SUNDIALS_KERNEL_NAME(...) __VA_ARGS__
#ifndef SUNDIALS_DEBUG_HIP_LASTERROR
#define SUNDIALS_LAUNCH_KERNEL(kernel, gridDim, blockDim, shMem, stream, ...) \
  {                                                                           \
    kernel<<<gridDim, blockDim, shMem, stream>>>(__VA_ARGS__);                \
  }
#else
#define SUNDIALS_LAUNCH_KERNEL(kernel, gridDim, blockDim, shMem, stream, ...) \
  {                                                                           \
    kernel<<<gridDim, blockDim, shMem, stream>>>(__VA_ARGS__);                \
    hipDeviceSynchronize();                                                   \
    SUNDIALS_HIP_VERIFY(hipGetLastError());                                   \
  }
#endif

/* ---------------------------------------------------------------------------
 * Utility functions
 * ---------------------------------------------------------------------------*/
inline sunbooleantype SUNDIALS_HIP_Assert(hipError_t hiperr, const char* file,
                                          int line)
{
  if (hiperr != hipSuccess)
  {
#ifdef SUNDIALS_DEBUG
    fprintf(stderr, "ERROR in HIP runtime operation: %s %s:%d\n",
            hipGetErrorString(hiperr), file, line);
#endif
    return SUNFALSE; /* Assert failed */
  }
  return SUNTRUE; /* Assert OK */
}

#ifdef __cplusplus /* wrapper to enable C++ usage */
}
#endif

#endif /* _SUNDIALS_HIP_H */
