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
 * This file is the entrypoint for the Python binding code for the
 * SUNDIALS N_Vector class. It contains hand-written code for functions
 * that require special treatment, and includes the generated code
 * produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <optional>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_nvector.hpp>
#include "sundials/sundials_config.h"

#include "../nvector/nvector_array_helpers.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

namespace {

void replace_host_array_pointer(N_Vector v, sunrealtype* ptr, nb::object owner)
{
  auto vector_id = N_VGetVectorID(v);
  if (vector_id == SUNDIALS_NVEC_SERIAL)
  {
    auto content = static_cast<N_VectorContent_Serial>(v->content);
    if (content->length == 0 || content->data == ptr) { return; }

    nvector_detail::prepare_python_array_owners(v);
    if (content->own_data && content->data) { std::free(content->data); }
    content->data     = ptr;
    content->own_data = SUNFALSE;
    nvector_detail::retain_python_host_array(v, std::move(owner));
    return;
  }

#ifdef SUNDIALS_NVECTOR_CUDA
  if (vector_id == SUNDIALS_NVEC_CUDA)
  {
    nvector_detail::replace_cuda_array_pointer(v, ptr,
                                               nvector_detail::ArrayDevice::Cpu,
                                               std::move(owner));
    return;
  }
#endif

  throw sundials4py::error_returned(
    "Ownership-safe host pointer replacement is not supported for this "
    "N_Vector implementation");
}

nvector_detail::ArrayDevice check_jax_array(N_Vector v, nb::object data, bool copy)
{
  if (copy && !nvector_detail::is_jax_array(data))
  {
    throw sundials4py::error_returned(
      "N_VSetJaxArray with copy=True requires a jax.Array");
  }
  if (!copy && !nvector_detail::is_jax_ref(data))
  {
    throw sundials4py::error_returned(
      "N_VSetJaxArray with copy=False requires a jax.ref.Ref");
  }

  std::string platform;
  if (copy)
  {
    auto devices = nb::cast<nb::list>(
      nb::module_::import_("builtins").attr("list")(data.attr("devices")()));
    if (devices.size() != 1)
    {
      throw sundials4py::error_returned(
        "N_VSetJaxArray requires a single-device jax.Array");
    }
    platform = nb::cast<std::string>(devices[0].attr("platform"));
  }
  else
  {
    auto devices = nb::cast<nb::list>(
      nb::module_::import_("builtins")
        .attr("list")(data.attr("sharding").attr("device_set")));
    if (devices.size() != 1)
    {
      throw sundials4py::error_returned(
        "N_Vector storage requires a single-device JAX Ref");
    }
    platform = nb::cast<std::string>(devices[0].attr("platform"));
  }

  nvector_detail::ArrayDevice device;
  if (platform == "cpu") { device = nvector_detail::ArrayDevice::Cpu; }
  else if (platform == "cuda" || platform == "gpu")
  {
    device = nvector_detail::ArrayDevice::Cuda;
  }
  else
  {
    throw sundials4py::error_returned(
      "N_VSetJaxArray only supports JAX CPU and CUDA devices");
  }

  auto shape = nb::cast<nb::tuple>(data.attr("shape"));
  if (shape.size() != 1)
  {
    throw sundials4py::error_returned(
      "N_VSetJaxArray requires a one-dimensional array");
  }
  nvector_detail::require_vector_length(N_VGetLength(v), data, "JAX array");

  auto dtype = data.attr("dtype");
  if (nb::cast<std::string>(dtype.attr("kind")) != "f")
  {
    throw sundials4py::error_returned(
      "N_VSetJaxArray requires a floating-point dtype");
  }
  if (nb::cast<size_t>(dtype.attr("itemsize")) != sizeof(sunrealtype))
  {
    throw sundials4py::error_returned(
      "JAX array dtype does not match SUNDIALS precision");
  }

  return device;
}

void set_jax_array(nb::object data, N_Vector v, std::optional<bool> copy)
{
  if (!copy.has_value())
  {
    if (nvector_detail::is_jax_ref(data)) { copy = false; }
    else if (nvector_detail::is_jax_array(data)) { copy = true; }
    else
    {
      throw sundials4py::error_returned(
        "N_VSetJaxArray requires a jax.Array or jax.ref.Ref when "
        "copy=None");
    }
  }

  const bool copy_array = *copy;
  auto device           = check_jax_array(v, data, copy_array);
  auto jax              = nb::module_::import_("jax");
  if (copy_array) { data.attr("block_until_ready")(); }
  else { jax.attr("effects_barrier")(); }

  auto source_address =
    nb::cast<std::uintptr_t>(data.attr("unsafe_buffer_pointer")());
  auto source    = reinterpret_cast<sunrealtype*>(source_address);
  auto vector_id = N_VGetVectorID(v);

  if (vector_id == SUNDIALS_NVEC_SERIAL)
  {
    if (device != nvector_detail::ArrayDevice::Cpu)
    {
      throw sundials4py::error_returned(
        "A serial N_Vector requires a JAX array on the CPU");
    }

    if (copy_array)
    {
      auto destination = N_VGetArrayPointer(v);
      auto length      = N_VGetLength(v);
      if (length > 0 && !destination)
      {
        throw sundials4py::error_returned(
          "Serial N_Vector does not have host storage");
      }
      if (length > 0 && !source)
      {
        throw sundials4py::error_returned("JAX array pointer must not be null");
      }
      if (length > 0)
      {
        std::memcpy(destination, source, length * sizeof(sunrealtype));
      }
    }
    else { replace_host_array_pointer(v, source, std::move(data)); }
    return;
  }

#ifdef SUNDIALS_NVECTOR_CUDA
  if (vector_id == SUNDIALS_NVEC_CUDA)
  {
    if (copy_array)
    {
      nvector_detail::copy_to_cuda_array_pointer(v, source, device);
    }
    else
    {
      if (device != nvector_detail::ArrayDevice::Cuda)
      {
        throw sundials4py::error_returned(
          "A CUDA N_Vector requires a JAX Ref on a CUDA device");
      }
      nvector_detail::replace_cuda_array_pointer(v, source,
                                                 nvector_detail::ArrayDevice::Cuda,
                                                 std::move(data));
    }
    return;
  }
#endif

  throw sundials4py::error_returned(
    "N_VSetJaxArray only supports serial and CUDA N_Vectors");
}

#ifndef SUNDIALS_NVECTOR_CUDA
using nvector_detail::ArrayDevice;
using nvector_detail::host_array;

nb::object get_numpy_array(N_Vector v) { return nb::cast(host_array(v)); }

nb::object get_torch_tensor(N_Vector v)
{
  return nb::module_::import_("torch").attr("from_numpy")(nb::cast(host_array(v)));
}

nb::object get_cupy_array(N_Vector v)
{
  (void)v;
  throw sundials4py::error_returned(
    "CUDA array access requires a CUDA-enabled sundials4py build");
}

nb::object get_jax_array(N_Vector v)
{
  return nb::module_::import_("jax").attr("dlpack").attr(
    "from_dlpack")(nb::cast(host_array(v)), nb::none(), nb::none());
}

#endif

} // namespace

void bind_nvector(nb::module_& m)
{
#include "sundials_nvector_generated.hpp"

  m.def(
    "N_VGetArrayPointer",
    [](N_Vector v) { return nvector_detail::host_array(v); },
    nb::rv_policy::reference);

#ifndef SUNDIALS_NVECTOR_CUDA
  m.def(
    "N_VGetDeviceArrayPointer",
    [](N_Vector v) -> std::uintptr_t
    {
      auto ptr = N_VGetDeviceArrayPointer(v);
      if (!ptr)
      {
        throw sundials4py::error_returned("Failed to get device array pointer");
      }
      return reinterpret_cast<std::uintptr_t>(ptr);
    },
    nb::arg("v"));

  nvector_detail::bind_nvector_array_accessors(m, &get_numpy_array,
                                               &get_jax_array, &get_cupy_array,
                                               &get_torch_tensor);
#endif

  m.def("N_VSetJaxArray", &set_jax_array, nb::arg("array"), nb::arg("v"),
        nb::kw_only(), nb::arg("copy").none() = nb::none(),
        "Set an N_Vector from a JAX array.\n\n"
        "If copy is None, its value is inferred from array: jax.Array uses "
        "copy=True and jax.ref.Ref uses copy=False. With copy=True, array "
        "must be a one-dimensional jax.Array whose "
        "length and floating-point dtype match the N_Vector. Its values are "
        "copied into the existing N_Vector storage. With copy=False, array "
        "must be a one-dimensional jax.ref.Ref; its storage is attached to "
        "the N_Vector without copying, and the Ref is retained as its owner. "
        "Only CPU storage is supported for serial N_Vectors, and only CUDA "
        "storage is supported for CUDA N_Vectors.");

  m.def(
    "N_VScaleAddMultiVectorArray",
    [](int nvec, int nsum, sundials4py::Array1d c_1d,
       std::vector<N_Vector> X_1d, std::vector<std::vector<N_Vector>> Y_2d,
       std::vector<std::vector<N_Vector>> Z_2d) -> SUNErrCode
    {
      sunrealtype* c_1d_ptr = reinterpret_cast<sunrealtype*>(c_1d.data());
      N_Vector* X_1d_ptr =
        reinterpret_cast<N_Vector*>(X_1d.empty() ? nullptr : X_1d.data());

      // Convert Y_2d and Z_2d to N_Vector**
      std::vector<N_Vector*> Y_2d_ptrs, Z_2d_ptrs;
      for (auto& row : Y_2d) { Y_2d_ptrs.push_back(row.data()); }
      for (auto& row : Z_2d) { Z_2d_ptrs.push_back(row.data()); }

      N_Vector** Y_2d_ptr = Y_2d_ptrs.data();
      N_Vector** Z_2d_ptr = Z_2d_ptrs.data();

      auto lambda_result = N_VScaleAddMultiVectorArray(nvec, nsum, c_1d_ptr,
                                                       X_1d_ptr, Y_2d_ptr,
                                                       Z_2d_ptr);
      return lambda_result;
    },
    nb::arg("nvec"), nb::arg("nsum"), nb::arg("c_1d"), nb::arg("X_1d"),
    nb::arg("Y_2d"), nb::arg("Z_2d"));

  m.def(
    "N_VLinearCombinationVectorArray",
    [](int nvec, int nsum, sundials4py::Array1d c_1d,
       std::vector<std::vector<N_Vector>> X_2d,
       std::vector<N_Vector> Z_1d) -> SUNErrCode
    {
      sunrealtype* c_1d_ptr = reinterpret_cast<sunrealtype*>(c_1d.data());

      // Convert X_2d to N_Vector**
      std::vector<N_Vector*> X_2d_ptrs;
      for (auto& row : X_2d) { X_2d_ptrs.push_back(row.data()); }
      N_Vector** X_2d_ptr = X_2d_ptrs.data();

      N_Vector* Z_1d_ptr =
        reinterpret_cast<N_Vector*>(Z_1d.empty() ? nullptr : Z_1d.data());

      auto lambda_result = N_VLinearCombinationVectorArray(nvec, nsum, c_1d_ptr,
                                                           X_2d_ptr, Z_1d_ptr);
      return lambda_result;
    },
    nb::arg("nvec"), nb::arg("nsum"), nb::arg("c_1d"), nb::arg("X_2d"),
    nb::arg("Z_1d"));
}

} // namespace sundials4py
