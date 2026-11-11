/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * All rights reserved.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Initial HIP N_Vector bindings. Array ownership follows the CUDA bindings.
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <cstdint>
#include <nvector/nvector_hip.h>
#include <string>
#include <utility>

#include "nvector_array_helpers.hpp"

namespace nb = nanobind;

namespace sundials4py {

using HipArray1d =
  nb::ndarray<sunrealtype, nb::device::rocm, nb::ndim<1>, nb::c_contig>;

namespace {

void require_hip_nvector(N_Vector v)
{
  if (N_VGetVectorID(v) != SUNDIALS_NVEC_HIP)
  {
    throw sundials4py::error_returned("Expected a HIP N_Vector");
  }
}

hipStream_t vector_stream(N_Vector v)
{
  auto content = static_cast<N_VectorContent_Hip>(v->content);
  return content->stream_exec_policy ? *content->stream_exec_policy->stream()
                                     : hipStream_t{};
}

void sync_device(N_Vector v)
{
  require_hip_nvector(v);
  const auto status = hipStreamSynchronize(vector_stream(v));
  if (status != hipSuccess)
  {
    (void)hipGetLastError();
    throw sundials4py::error_returned("Failed to synchronize HIP N_Vector (" +
                                      std::string(hipGetErrorName(status)) +
                                      ": " + hipGetErrorString(status) + ")");
  }
}

int device_id_for_pointer(void* ptr)
{
  hipPointerAttribute_t attrs;
  if (hipPointerGetAttributes(&attrs, ptr) == hipSuccess)
  {
    return attrs.device;
  }
  (void)hipGetLastError();
  throw sundials4py::error_returned("HIP device pointer is not accessible");
}

struct HipDeviceArrayPointer
{
  N_Vector x;
  nb::object owner;
  int device_id;

  explicit HipDeviceArrayPointer(N_Vector v)
    : HipDeviceArrayPointer(v, nb::find(v))
  {}

  // owner is the Python object that must keep the data alive. It differs from
  // v when v is the local vector of an MPIPlusX vector (the owner is then the
  // parent MPIPlusX wrapper).
  HipDeviceArrayPointer(N_Vector v, nb::handle owner)
    : x(v), owner(nb::borrow(owner)), device_id(0)
  {
    require_hip_nvector(v);
    auto ptr = N_VGetDeviceArrayPointer_Hip(v);
    if (!ptr)
    {
      throw sundials4py::error_returned("HIP N_Vector has no device storage");
    }
    device_id = device_id_for_pointer(ptr);
  }

  HipArray1d dlpack() const
  {
    auto ptr = N_VGetDeviceArrayPointer_Hip(x);
    if (!ptr)
    {
      throw sundials4py::error_returned("HIP N_Vector has no device storage");
    }
    size_t shape[1]{static_cast<size_t>(N_VGetLength(x))};
    // Pass the device id explicitly: the default (0) would mislabel buffers
    // that live on other devices.
    return HipArray1d(ptr, 1, shape, owner, /*strides=*/nullptr,
                      nb::dtype<sunrealtype>(), nb::device::rocm::value,
                      device_id);
  }
};

nb::object device_pointer_object(N_Vector vec)
{
  // MPIPlusX vectors resolve to their local vector (MPIManyVector raises).
  auto target = nvector_detail::resolve_array_target(vec);
  require_hip_nvector(target.v);
  sync_device(target.v);
  return nb::cast(HipDeviceArrayPointer(target.v, target.owner));
}

nb::object get_cupy_array(N_Vector v)
{
  auto pointer = device_pointer_object(v);
  return nb::module_::import_("cupy").attr("from_dlpack")(pointer);
}

nb::object get_torch_tensor(N_Vector v)
{
  auto pointer = device_pointer_object(v);
  return nb::module_::import_("torch").attr("utils").attr("dlpack").attr(
    "from_dlpack")(pointer);
}

nb::object get_jax_array(N_Vector v)
{
  auto pointer = device_pointer_object(v);
  return nb::module_::import_("jax").attr("dlpack").attr(
    "from_dlpack")(pointer, nb::none(), false);
}

void require_hip_pointer(const sunrealtype* ptr)
{
  hipPointerAttribute_t attrs;
  if (hipPointerGetAttributes(&attrs, ptr) != hipSuccess ||
      (attrs.type != hipMemoryTypeDevice &&
       attrs.type != hipMemoryTypeManaged && attrs.type != hipMemoryTypeUnified))
  {
    (void)hipGetLastError();
    throw sundials4py::error_returned("Array pointer is not HIP-accessible");
  }
}

void set_jax_array(nb::object data, N_Vector v, bool copy)
{
  // MPIPlusX vectors act on their local vector (MPIManyVector raises).
  v = nvector_detail::resolve_array_target(v).v;
  require_hip_nvector(v);
  if (copy ? !nvector_detail::is_jax_array(data)
           : !nvector_detail::is_jax_ref(data))
  {
    throw sundials4py::error_returned(
      "copy=True requires a JAX array; copy=False requires a JAX Ref");
  }
  nvector_detail::require_vector_length(N_VGetLength(v), data, "JAX array");
  auto dtype = data.attr("dtype");
  if (nb::cast<std::string>(dtype.attr("kind")) != "f" ||
      nb::cast<size_t>(dtype.attr("itemsize")) != sizeof(sunrealtype))
  {
    throw sundials4py::error_returned(
      "JAX array dtype does not match SUNDIALS precision");
  }

  auto jax = nb::module_::import_("jax");
  if (copy) { data.attr("block_until_ready")(); }
  else { jax.attr("effects_barrier")(); }
  auto address = nb::cast<std::uintptr_t>(data.attr("unsafe_buffer_pointer")());
  auto source  = reinterpret_cast<sunrealtype*>(address);
  if (N_VGetLength(v) == 0) { return; }
  require_hip_pointer(source);

  sync_device(v);
  if (copy)
  {
    auto destination = N_VGetDeviceArrayPointer_Hip(v);
    if (!destination ||
        hipMemcpy(destination, source, N_VGetLength(v) * sizeof(sunrealtype),
                  hipMemcpyDeviceToDevice) != hipSuccess)
    {
      (void)hipGetLastError();
      throw sundials4py::error_returned(
        "Failed to copy JAX array to HIP N_Vector");
    }
    return;
  }

  if (N_VIsManagedMemory_Hip(v))
  {
    throw sundials4py::error_returned(
      "JAX Ref attachment requires a non-managed HIP N_Vector");
  }
  auto content = static_cast<N_VectorContent_Hip>(v->content);
  if (content->device_data && content->device_data->ptr == source) { return; }
  auto replacement = SUNMemoryHelper_Wrap(content->mem_helper, source,
                                          SUNMEMTYPE_DEVICE);
  if (!replacement)
  {
    throw sundials4py::error_returned("Failed to wrap JAX Ref storage");
  }
  nvector_detail::prepare_python_array_owners(v);
  auto previous        = content->device_data;
  content->device_data = replacement;
  nvector_detail::retain_python_device_array(v, std::move(data));
  if (SUNMemoryHelper_Dealloc(content->mem_helper, previous, nullptr) !=
      SUN_SUCCESS)
  {
    throw sundials4py::error_returned(
      "Failed to release old HIP device storage");
  }
}

} // namespace

void bind_nvector_hip(nb::module_& m)
{
#include "nvector_hip_generated.hpp"

  nb::class_<HipDeviceArrayPointer>(m, "_HipDeviceArrayPointer")
    .def(
      "__dlpack__",
      [](const HipDeviceArrayPointer& self, nb::object /* stream */)
      { return self.dlpack(); },
      nb::arg("stream") = nb::none())
    .def("__dlpack_device__", [](const HipDeviceArrayPointer& self)
         { return nb::make_tuple(10, self.device_id); });

  m.def("N_VGetCupyArray_Hip", &get_cupy_array, nb::arg("v"));
  m.def("N_VGetTorchTensor_Hip", &get_torch_tensor, nb::arg("v"));
  m.def("N_VGetJaxArray_Hip", &get_jax_array, nb::arg("v"));
  m.def("N_VSetJaxArray_Hip", &set_jax_array, nb::arg("array"), nb::arg("v"),
        nb::kw_only(), nb::arg("copy") = true);

  m.def(
    "N_VNewEmpty_Hip",
    [](SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    { return nvector_detail::wrap_nvector(N_VNewEmpty_Hip(sunctx)); },
    nb::arg("sunctx"), nb::keep_alive<0, 1>());

  m.def(
    "N_VNew_Hip",
    [](sunindextype vec_length,
       SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    { return nvector_detail::wrap_nvector(N_VNew_Hip(vec_length, sunctx)); },
    nb::arg("vec_length"), nb::arg("sunctx"), nb::keep_alive<0, 2>());

  m.def(
    "N_VNewManaged_Hip",
    [](sunindextype vec_length,
       SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>> {
      return nvector_detail::wrap_nvector(N_VNewManaged_Hip(vec_length, sunctx));
    },
    nb::arg("vec_length"), nb::arg("sunctx"), nb::keep_alive<0, 2>());

  m.def(
    "N_VNewWithMemHelp_Hip",
    [](sunindextype vec_length, sunbooleantype use_managed_mem,
       SUNMemoryHelper helper,
       SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      return nvector_detail::wrap_nvector(
        N_VNewWithMemHelp_Hip(vec_length, use_managed_mem, helper, sunctx));
    },
    nb::arg("vec_length"), nb::arg("use_managed_mem"), nb::arg("helper"),
    nb::arg("sunctx"), nb::keep_alive<0, 3>(), nb::keep_alive<0, 4>());

  m.def(
    "N_VMake_Hip",
    [](sunindextype vec_length, sundials4py::Array1d h_vdata_1d,
       HipArray1d d_vdata_1d,
       SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      nvector_detail::require_vector_length(vec_length, h_vdata_1d);
      nvector_detail::require_vector_length(vec_length, d_vdata_1d, "HIP array");
      auto v = N_VMake_Hip(vec_length, h_vdata_1d.data(), d_vdata_1d.data(),
                           sunctx);
      return nvector_detail::wrap_nvector(v);
    },
    nb::arg("vec_length"), nb::arg("h_vdata_1d").noconvert(),
    nb::arg("d_vdata_1d").noconvert(), nb::arg("sunctx"),
    nb::keep_alive<0, 2>(), nb::keep_alive<0, 3>(), nb::keep_alive<0, 4>());

  m.def(
    "N_VMake_Hip",
    [](sunindextype vec_length, sundials4py::Array1d h_vdata_1d,
       nb::object d_vdata_1d,
       SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      nvector_detail::require_vector_length(vec_length, h_vdata_1d);
      nvector_detail::require_vector_length(vec_length, d_vdata_1d);
      if (!PyObject_HasAttrString(d_vdata_1d.ptr(), "unsafe_buffer_pointer"))
      {
        throw sundials4py::error_returned(
          "HIP device array must provide unsafe_buffer_pointer()");
      }
      auto address =
        nb::cast<std::uintptr_t>(d_vdata_1d.attr("unsafe_buffer_pointer")());
      auto v = N_VMake_Hip(vec_length, h_vdata_1d.data(),
                           reinterpret_cast<sunrealtype*>(address), sunctx);
      return nvector_detail::wrap_nvector(v);
    },
    nb::arg("vec_length"), nb::arg("h_vdata_1d").noconvert(),
    nb::arg("d_vdata_1d"), nb::arg("sunctx"), nb::keep_alive<0, 2>(),
    nb::keep_alive<0, 3>(), nb::keep_alive<0, 4>());

  m.def(
    "N_VMakeManaged_Hip",
    [](sunindextype vec_length, std::uintptr_t vdata,
       SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      return nvector_detail::wrap_nvector(
        N_VMakeManaged_Hip(vec_length, reinterpret_cast<sunrealtype*>(vdata),
                           sunctx));
    },
    nb::arg("vec_length"), nb::arg("vdata"), nb::arg("sunctx"),
    nb::keep_alive<0, 3>());

  m.def(
    "N_VGetDeviceArrayPointer_Hip",
    [](N_Vector vec) -> std::uintptr_t
    {
      // MPIPlusX vectors resolve to their local vector (MPIManyVector raises).
      N_Vector v = nvector_detail::resolve_array_target(vec).v;
      if (N_VGetVectorID(v) != SUNDIALS_NVEC_HIP)
      {
        throw sundials4py::error_returned(
          "HIP device array access requires a HIP N_Vector");
      }
      auto ptr = N_VGetDeviceArrayPointer_Hip(v);
      if (!ptr)
      {
        throw sundials4py::error_returned(
          "Failed to get HIP device array pointer");
      }
      return reinterpret_cast<std::uintptr_t>(ptr);
    },
    nb::arg("v"));

  m.def(
    "N_VCopyToDevice_Hip",
    [](N_Vector v)
    {
      if (N_VGetVectorID(v) != SUNDIALS_NVEC_HIP)
      {
        throw sundials4py::error_returned("Expected a HIP N_Vector");
      }
      N_VCopyToDevice_Hip(v);
    },
    nb::arg("v"));

  m.def(
    "N_VCopyFromDevice_Hip",
    [](N_Vector v)
    {
      if (N_VGetVectorID(v) != SUNDIALS_NVEC_HIP)
      {
        throw sundials4py::error_returned("Expected a HIP N_Vector");
      }
      N_VCopyFromDevice_Hip(v);
    },
    nb::arg("v"));
}

} // namespace sundials4py
