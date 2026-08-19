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

#include <unordered_map>

#include <sundials/sundials_nvector.h>

namespace sundials4py {
namespace nvector_detail {

enum class ArrayDevice
{
  Cpu,
  Cuda
};

using NVectorDestroyFn = void (*)(N_Vector);

struct PythonArrayOwners
{
  nanobind::object host_array;
  nanobind::object device_array;
  NVectorDestroyFn destroy = nullptr;
};

inline std::unordered_map<N_Vector, PythonArrayOwners>& python_array_owners()
{
  // Intentionally keep this map alive until process exit. Destroying a static
  // map of Python objects after Python finalization is unsafe.
  static auto* owners = new std::unordered_map<N_Vector, PythonArrayOwners>;
  return *owners;
}

inline void destroy_nvector_with_python_arrays(N_Vector v)
{
  NVectorDestroyFn destroy = nullptr;
  nanobind::object host_array;
  nanobind::object device_array;

  {
    nanobind::gil_scoped_acquire gil;
    auto& owners = python_array_owners();
    auto it      = owners.find(v);
    if (it != owners.end())
    {
      destroy      = it->second.destroy;
      host_array   = std::move(it->second.host_array);
      device_array = std::move(it->second.device_array);
      owners.erase(it);
    }
  }

  if (destroy) { destroy(v); }

  // Release Python references only while holding the GIL and after SUNDIALS
  // has finished destroying its non-owning memory wrappers.
  if (host_array.ptr() || device_array.ptr())
  {
    nanobind::gil_scoped_acquire gil;
    host_array   = nanobind::object();
    device_array = nanobind::object();
  }
}

inline PythonArrayOwners& prepare_python_array_owners(N_Vector v)
{
  auto& owners        = python_array_owners();
  auto [it, inserted] = owners.try_emplace(v);
  if (inserted)
  {
    it->second.destroy = v->ops->nvdestroy;
    v->ops->nvdestroy  = destroy_nvector_with_python_arrays;
  }
  return it->second;
}

inline void retain_python_host_array(N_Vector v, nanobind::object array)
{
  prepare_python_array_owners(v).host_array = std::move(array);
}

inline void retain_python_device_array(N_Vector v, nanobind::object array)
{
  prepare_python_array_owners(v).device_array = std::move(array);
}

inline void release_python_array(N_Vector v, ArrayDevice device)
{
  auto& owners = python_array_owners();
  auto it      = owners.find(v);
  if (it == owners.end()) { return; }

  if (device == ArrayDevice::Cpu)
  {
    it->second.host_array = nanobind::object();
  }
  else { it->second.device_array = nanobind::object(); }

  if (it->second.host_array.ptr() || it->second.device_array.ptr()) { return; }

  if (v && v->ops && v->ops->nvdestroy == destroy_nvector_with_python_arrays)
  {
    v->ops->nvdestroy = it->second.destroy;
  }
  owners.erase(it);
}

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

#ifdef SUNDIALS_NVECTOR_CUDA
void replace_cuda_array_pointer(N_Vector v, sunrealtype* ptr, ArrayDevice device,
                                nanobind::object owner = nanobind::object());
#endif

} // namespace nvector_detail
} // namespace sundials4py

#endif
