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

#include <sundials/sundials_nvector.hpp>
#include "sundials/sundials_config.h"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

namespace {

#ifndef SUNDIALS_NVECTOR_CUDA
enum class ArrayDevice
{
  Cpu,
  Cuda
};

bool object_is_none(nb::object obj) { return obj.ptr() == Py_None; }

std::string optional_string(nb::object value, const char* name)
{
  if (object_is_none(value)) { return ""; }
  try
  {
    return nb::cast<std::string>(value);
  }
  catch (const nb::cast_error&)
  {
    throw sundials4py::error_returned(std::string(name) +
                                      " must be a string or None");
  }
}

ArrayDevice parse_device(nb::object device, N_Vector /*v*/)
{
  auto value = optional_string(device, "device");
  if (value.empty() || value == "cpu" || value == "host")
  {
    return ArrayDevice::Cpu;
  }
  if (value == "cuda") { return ArrayDevice::Cuda; }

  throw sundials4py::error_returned(
    "device must be 'cpu', 'host', 'cuda', or None");
}

void validate_copy_from(nb::object copy_from)
{
  auto value = optional_string(copy_from, "copy_from");
  if (value.empty() || value == "cpu") { return; }
  if (value == "device")
  {
    throw sundials4py::error_returned(
      "copy_from='device' requires a CUDA N_Vector");
  }

  throw sundials4py::error_returned(
    "copy_from must be 'cpu', 'device', or None");
}

sundials4py::Array1d host_array(N_Vector v)
{
  auto ptr = N_VGetArrayPointer(v);
  if (!ptr)
  {
    throw sundials4py::error_returned("Failed to get array pointer");
  }
  auto owner = nb::find(v);
  size_t shape[1]{static_cast<size_t>(N_VGetLength(v))};
  return sundials4py::Array1d(ptr, 1, shape, owner);
}

nb::object get_numpy_array(N_Vector v, nb::object device, nb::object copy_from)
{
  if (parse_device(device, v) == ArrayDevice::Cuda)
  {
    throw sundials4py::error_returned(
      "N_VGetNumpyArray does not support device='cuda'");
  }
  validate_copy_from(copy_from);
  return nb::cast(host_array(v));
}

nb::object get_torch_tensor(N_Vector v, nb::object device, nb::object copy_from)
{
  if (parse_device(device, v) == ArrayDevice::Cuda)
  {
    throw sundials4py::error_returned(
      "CUDA tensor access requires a CUDA-enabled sundials4py build");
  }
  validate_copy_from(copy_from);
  return nb::module_::import_("torch").attr("from_numpy")(nb::cast(host_array(v)));
}

nb::object get_cupy_array(N_Vector v, nb::object device, nb::object copy_from)
{
  (void)v;
  if (parse_device(device, v) == ArrayDevice::Cpu)
  {
    throw sundials4py::error_returned(
      "N_VGetCupyArray requires a CUDA N_Vector");
  }
  validate_copy_from(copy_from);
  throw sundials4py::error_returned(
    "CUDA array access requires a CUDA-enabled sundials4py build");
}

nb::object get_jax_array(N_Vector v, nb::object device, nb::object copy_from)
{
  if (parse_device(device, v) == ArrayDevice::Cuda)
  {
    throw sundials4py::error_returned(
      "CUDA array access requires a CUDA-enabled sundials4py build");
  }
  validate_copy_from(copy_from);
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

  m.def("N_VGetNumpyArray", &get_numpy_array, nb::arg("v"),
        nb::arg("device").none()    = nb::none(),
        nb::arg("copy_from").none() = nb::none());

  m.def("N_VGetJaxArray", &get_jax_array, nb::arg("v"),
        nb::arg("device").none()    = nb::none(),
        nb::arg("copy_from").none() = nb::none());

  m.def("N_VGetCupyArray", &get_cupy_array, nb::arg("v"),
        nb::arg("device").none()    = nb::none(),
        nb::arg("copy_from").none() = nb::none());

  m.def("N_VGetTorchTensor", &get_torch_tensor, nb::arg("v"),
        nb::arg("device").none()    = nb::none(),
        nb::arg("copy_from").none() = nb::none());
#endif

  m.def("N_VSetArrayPointer",
        [](sundials4py::Array1d arr, N_Vector v)
        {
          if (arr.shape(0) != N_VGetLength(v))
          {
            throw sundials4py::error_returned(
              "Array shape does not match vector length");
          }
          N_VSetArrayPointer(arr.data(), v);
        });

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
