/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * Shared helpers for exposing N_Vector data to Python array backends.
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS4PY_NVECTOR_ARRAY_HELPERS_HPP
#define SUNDIALS4PY_NVECTOR_ARRAY_HELPERS_HPP

#include "sundials4py.hpp"

#include <sundials/sundials_nvector.h>

namespace sundials4py {
namespace nvector_detail {

enum class ArrayDevice
{
  Cpu,
  Cuda
};

inline bool object_is_none(nanobind::object obj)
{
  return obj.ptr() == Py_None;
}

inline std::string optional_string(nanobind::object value)
{
  if (object_is_none(value)) { return ""; }
  try
  {
    return nanobind::cast<std::string>(value);
  }
  catch (const nanobind::cast_error&)
  {
    throw sundials4py::error_returned("Expected a string or None");
  }
}

inline bool is_cuda_nvector(N_Vector v)
{
  return N_VGetVectorID(v) == SUNDIALS_NVEC_CUDA;
}

inline ArrayDevice parse_device(nanobind::object device, N_Vector v)
{
  auto value = optional_string(device);
  if (value.empty())
  {
    return is_cuda_nvector(v) ? ArrayDevice::Cuda : ArrayDevice::Cpu;
  }
  if (value == "cpu" || value == "host") { return ArrayDevice::Cpu; }
  if (value == "cuda") { return ArrayDevice::Cuda; }

  throw sundials4py::error_returned(
    "device must be 'cpu', 'host', 'cuda', or None");
}

inline sundials4py::Array1d host_array(N_Vector v)
{
  auto ptr = N_VGetArrayPointer(v);
  if (!ptr)
  {
    throw sundials4py::error_returned("Failed to get array pointer");
  }
  auto owner = nanobind::find(v);
  size_t shape[1]{static_cast<size_t>(N_VGetLength(v))};
  return sundials4py::Array1d(ptr, 1, shape, owner);
}

} // namespace nvector_detail
} // namespace sundials4py

#endif
