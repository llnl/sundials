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

#include <memory>
#include <type_traits>
#include <unordered_map>

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_nvector.hpp>

namespace sundials4py {
namespace nvector_detail {

inline std::shared_ptr<std::remove_pointer_t<N_Vector>> wrap_nvector(N_Vector v)
{
  return sundials::experimental::our_make_shared<
    std::remove_pointer_t<N_Vector>, sundials::experimental::N_VectorDeleter>(v);
}

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
  // Intentionally process-lifetime: destroying Python objects after interpreter
  // finalization is unsafe. Entries are still erased when each N_Vector dies.
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
  // Post-construction zero-copy replacement cannot use return-value keep_alive;
  // retain the JAX Ref here to prevent early free without affecting hot paths.
  prepare_python_array_owners(v).host_array = std::move(array);
}

inline void retain_python_device_array(N_Vector v, nanobind::object array)
{
  prepare_python_array_owners(v).device_array = std::move(array);
}

inline bool object_is_none(nanobind::object obj)
{
  return obj.ptr() == Py_None;
}

inline bool is_jax_object(nanobind::handle obj)
{
  auto module = nanobind::cast<std::string>(obj.type().attr("__module__"));
  return module == "jax" || module.rfind("jax.", 0) == 0 ||
         module == "jaxlib" || module.rfind("jaxlib.", 0) == 0;
}

inline bool is_jax_ref(nanobind::handle obj)
{
  if (!is_jax_object(obj)) { return false; }
  auto jax = nanobind::module_::import_("jax");
  return nanobind::isinstance(obj, jax.attr("ref").attr("Ref"));
}

inline bool is_jax_array(nanobind::handle obj)
{
  // Some JAX versions make Ref satisfy the Array predicate. Refs have
  // different ownership and mutability semantics, so keep them separate.
  if (!is_jax_object(obj) || is_jax_ref(obj)) { return false; }
  auto jax = nanobind::module_::import_("jax");
  return nanobind::isinstance(obj, jax.attr("Array"));
}

inline std::string read_optional_string(nanobind::object value)
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
  auto value = read_optional_string(device);
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

template<typename GetNumpyArray, typename GetJaxArray, typename GetCupyArray,
         typename GetTorchTensor>
inline void bind_nvector_array_accessors(nanobind::module_& m,
                                         GetNumpyArray get_numpy_array,
                                         GetJaxArray get_jax_array,
                                         GetCupyArray get_cupy_array,
                                         GetTorchTensor get_torch_tensor)
{
  m.def("N_VGetNumpyArray", get_numpy_array, nanobind::arg("v"),
        "Return a NumPy view of an N_Vector's host data.\n\n"
        "The returned one-dimensional array shares storage with the N_Vector, "
        "so modifying the array modifies the vector. This function is not "
        "supported for CUDA N_Vectors; use a backend-specific CUDA accessor "
        "instead.");

  m.def("N_VGetJaxArray", get_jax_array, nanobind::arg("v"),
        "Return a JAX array view of an N_Vector's data.\n\n"
        "The returned array uses the N_Vector's native device and shares "
        "storage with the N_Vector. Use an explicit JAX or NumPy conversion "
        "when a host array is needed. JAX arrays are immutable; use "
        "N_VSetJaxArray to copy values into an N_Vector or to attach mutable "
        "JAX Ref storage.");

  m.def("N_VGetCupyArray", get_cupy_array, nanobind::arg("v"),
        "Return a CuPy view of a CUDA N_Vector's device data.\n\n"
        "This function requires a CUDA N_Vector and a CUDA-enabled "
        "sundials4py build. The returned array shares device storage with "
        "the N_Vector.");

  m.def("N_VGetTorchTensor", get_torch_tensor, nanobind::arg("v"),
        "Return a PyTorch tensor view of an N_Vector's data.\n\n"
        "The returned tensor uses the N_Vector's native device and shares "
        "storage with the N_Vector. Call .cpu() or another PyTorch transfer "
        "operation explicitly when a host copy is needed.");
}

template<typename Array>
inline void require_vector_length(sunindextype vec_length, const Array& array,
                                  const char* array_name = "Array")
{
  if (array.shape(0) != static_cast<size_t>(vec_length))
  {
    throw sundials4py::error_returned(std::string(array_name) +
                                      " shape does not match vector length");
  }
}

inline void require_vector_length(sunindextype vec_length, nanobind::object array,
                                  const char* array_name = "Device array")
{
  if (!PyObject_HasAttrString(array.ptr(), "shape"))
  {
    throw sundials4py::error_returned(std::string(array_name) +
                                      " does not have a shape");
  }

  auto shape = nanobind::cast<nanobind::tuple>(array.attr("shape"));
  if (shape.size() != 1 ||
      nanobind::cast<size_t>(shape[0]) != static_cast<size_t>(vec_length))
  {
    throw sundials4py::error_returned(std::string(array_name) +
                                      " shape does not match vector length");
  }
}

#ifdef SUNDIALS_NVECTOR_CUDA
void replace_cuda_array_pointer(N_Vector v, sunrealtype* ptr,
                                ArrayDevice device, nanobind::object owner);
void copy_to_cuda_array_pointer(N_Vector v, const sunrealtype* ptr,
                                ArrayDevice source_device);
#endif

} // namespace nvector_detail
} // namespace sundials4py

#endif
