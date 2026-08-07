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
 * ----------------------------------------------------------------------------
 * Header file for batched Brusselator problem using Ginkgo linear solvers.
 * See cv_bruss_batched_ginkgo.cpp for more information.
 * ---------------------------------------------------------------------------*/

#include <cstdio>
#include <memory>
#include <random>

#include <sundials/sundials_core.hpp>

#if defined(USE_CUDA)
#include <nvector/nvector_cuda.h>
#include <sunmemory/sunmemory_cuda.h>
#define SERIAL_CUDA_OR_HIP(a, b, c) b
constexpr auto N_VNew = N_VNewManaged_Cuda;
#elif defined(USE_HIP)
#include <nvector/nvector_hip.h>
#include <sunmemory/sunmemory_hip.h>
#define SERIAL_CUDA_OR_HIP(a, b, c) c
constexpr auto N_VNew = N_VNewManaged_Hip;
#else
#include <nvector/nvector_serial.h>
#include <sunmemory/sunmemory_system.h>
#define SERIAL_CUDA_OR_HIP(a, b, c) a
constexpr auto N_VNew = N_VNew_Serial;
#endif

/* Helper class for Arrays. Built on top of the SUNMemoryHelper. */
template<typename T, typename I>
class Array
{
public:
  Array(I size, SUNMemoryHelper helper) : mem_(nullptr), helper_(helper)
  {
    SUNMemoryHelper_Alloc(helper, &mem_, size * sizeof(T), SUNMEMTYPE_UVM, NULL);
  }

  Array(SUNMemory mem, SUNMemoryHelper helper) : mem_(mem), helper_(helper) {}

  T& operator[](int index) { return static_cast<T*>(mem_->ptr)[index]; }

  T* get() { return static_cast<T*>(mem_->ptr); }

private:
  SUNMemory mem_;
  SUNMemoryHelper helper_;
};

using RealArray = Array<sunrealtype, int>;

/* User data structure. This will be available
   in SUNDIALS callback functions. */
struct UserData
{
  UserData(int num_batches_in, int batch_size_in, int nnzper_in,
           SUNMemoryHelper h_in)
    : num_batches(num_batches_in),
      batch_size(batch_size_in),
      nnzper(nnzper_in),
      neq(batch_size_in * num_batches_in),
      u0{num_batches_in, h_in},
      v0{num_batches_in, h_in},
      w0{num_batches_in, h_in},
      a{num_batches_in, h_in},
      b{num_batches_in, h_in},
      ep{num_batches_in, h_in}
  {}

  int num_batches;      /* number of chemical networks  */
  int batch_size;       /* size of each network         */
  int nnzper;           /* number of nonzeros per batch */
  int neq;              /* total number of equations    */
  RealArray u0, v0, w0; /* initial conditions */
  RealArray a, b;       /* chemical concentrations that are constant */
  RealArray ep;
};
