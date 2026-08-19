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

void check_host_array_length(N_Vector v, size_t length)
{
  if (length != static_cast<size_t>(N_VGetLength(v)))
  {
    throw sundials4py::error_returned(
      "Array shape does not match vector length");
  }
}

void set_jax_host_array(nb::object data, N_Vector v)
{
  if (!PyObject_HasAttrString(data.ptr(), "unsafe_buffer_pointer") ||
      !PyObject_HasAttrString(data.ptr(), "device"))
  {
    throw sundials4py::error_returned(
      "Host array must be a NumPy array or a JAX CPU array");
  }

  auto platform = nb::cast<std::string>(data.attr("device").attr("platform"));
  if (platform != "cpu")
  {
    throw sundials4py::error_returned(
      "N_VSetArrayPointer requires a host-accessible array");
  }

  auto numpy_view = nb::module_::import_("numpy").attr("asarray")(data);
  auto array      = nb::cast<sundials4py::Array1d>(numpy_view);
  check_host_array_length(v, array.shape(0));

  auto source_ptr =
    nb::cast<std::uintptr_t>(data.attr("unsafe_buffer_pointer")());
  if (source_ptr != reinterpret_cast<std::uintptr_t>(array.data()))
  {
    throw sundials4py::error_returned(
      "Converting the host array would require a copy");
  }

  replace_host_array_pointer(v, array.data(), std::move(numpy_view));
}

#ifndef SUNDIALS_NVECTOR_CUDA
using nvector_detail::ArrayDevice;
using nvector_detail::host_array;
using nvector_detail::parse_device;

nb::object get_numpy_array(N_Vector v) { return nb::cast(host_array(v)); }

nb::object get_torch_tensor(N_Vector v, nb::object device)
{
  if (parse_device(device, v) == ArrayDevice::Cuda)
  {
    throw sundials4py::error_returned(
      "CUDA tensor access requires a CUDA-enabled sundials4py build");
  }
  return nb::module_::import_("torch").attr("from_numpy")(nb::cast(host_array(v)));
}

nb::object get_cupy_array(N_Vector v)
{
  (void)v;
  throw sundials4py::error_returned(
    "CUDA array access requires a CUDA-enabled sundials4py build");
}

nb::object get_jax_array(N_Vector v, nb::object device)
{
  if (parse_device(device, v) == ArrayDevice::Cuda)
  {
    throw sundials4py::error_returned(
      "CUDA array access requires a CUDA-enabled sundials4py build");
  }
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
    [](N_Vector v)
    {
      auto ptr = N_VGetArrayPointer(v);
      if (!ptr)
      {
        throw sundials4py::error_returned("Failed to get array pointer");
      }
      auto owner = nb::find(v);
      size_t shape[1]{static_cast<size_t>(N_VGetLength(v))};
      return sundials4py::Array1d(ptr, 1, shape, owner);
    },
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

  // Provide a generic N_VSetDeviceArrayPointer binding for non-CUDA builds
  // that accepts a host-side Array1d. CUDA-aware builds provide richer
  // overloads in nvector_cuda.cpp.
  m.def(
    "N_VSetDeviceArrayPointer",
    [](sundials4py::Array1d /* d_vdata_1d */, N_Vector /* v */)
    {
      throw sundials4py::error_returned(
        "Device pointer replacement requires a CUDA-enabled sundials4py "
        "build");
    },
    nb::arg("d_vdata_1d"), nb::arg("v"));

  m.def("N_VGetNumpyArray", &get_numpy_array, nb::arg("v"));

  m.def("N_VGetJaxArray", &get_jax_array, nb::arg("v"),
        nb::arg("device").none() = nb::none());

  m.def("N_VGetCupyArray", &get_cupy_array, nb::arg("v"));

  m.def("N_VGetTorchTensor", &get_torch_tensor, nb::arg("v"),
        nb::arg("device").none() = nb::none());
#endif

  m.def(
    "N_VSetArrayPointer",
    [](sundials4py::Array1d arr, N_Vector v)
    {
      check_host_array_length(v, arr.shape(0));
      replace_host_array_pointer(v, arr.data(), arr.cast());
    },
    nb::arg("arr"), nb::arg("v"),
    "Replace owned host storage with a borrowed NumPy array and retain its "
    "owner.");

  m.def("N_VSetArrayPointer", &set_jax_host_array, nb::arg("arr"), nb::arg("v"),
        "Replace owned host storage with a borrowed JAX CPU array and retain "
        "its owner.");

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
