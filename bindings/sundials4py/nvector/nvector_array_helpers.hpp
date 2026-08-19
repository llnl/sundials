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

enum class MemoryType
{
  Cpu,
  Cuda
};

enum class CopyFrom
{
  None,
  Cpu,
  Device
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

inline MemoryType parse_memory_type(nanobind::object memory_type, N_Vector v)
{
  auto value = optional_string(memory_type);
  if (value.empty())
  {
    return is_cuda_nvector(v) ? MemoryType::Cuda : MemoryType::Cpu;
  }
  if (value == "cpu" || value == "host") { return MemoryType::Cpu; }
  if (value == "cuda") { return MemoryType::Cuda; }

  throw sundials4py::error_returned(
    "memory_type must be 'cpu', 'host', 'cuda', or None");
}

inline CopyFrom parse_copy_from(nanobind::object copy_from)
{
  auto value = optional_string(copy_from);
  if (value.empty()) { return CopyFrom::None; }
  if (value == "cpu") { return CopyFrom::Cpu; }
  if (value == "device") { return CopyFrom::Device; }

  throw sundials4py::error_returned(
    "copy_from must be 'cpu', 'device', or None");
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
