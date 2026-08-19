#include "sundials4py.hpp"

#include <cstdint>

#include <nvector/nvector_cuda.h>
#include <sundials/sundials_nvector.hpp>

#include "nvector_array_helpers.hpp"
#include "sundials/sundials_classview.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

using CudaArray1d =
  nb::ndarray<sunrealtype, nb::device::cuda, nb::ndim<1>, nb::c_contig>;

namespace nvector_detail {

namespace {

void require_cuda_nvector(N_Vector v)
{
  if (!is_cuda_nvector(v))
  {
    throw sundials4py::error_returned(
      "CUDA array access requires a CUDA N_Vector");
  }
}

void* cuda_stream(N_Vector v)
{
  auto content = static_cast<N_VectorContent_Cuda>(v->content);
  return content->stream_exec_policy
           ? static_cast<void*>(
               const_cast<cudaStream_t*>(content->stream_exec_policy->stream()))
           : nullptr;
}

void require_compatible_pointer(N_Vector v, sunrealtype* ptr, ArrayDevice device)
{
  if (!ptr && N_VGetLength(v) == 0) { return; }
  if (!ptr)
  {
    throw sundials4py::error_returned("Array pointer must not be null");
  }

  cudaPointerAttributes attrs;
  auto status = cudaPointerGetAttributes(&attrs, ptr);
  if (N_VIsManagedMemory_Cuda(v))
  {
    if (status != cudaSuccess || attrs.type != cudaMemoryTypeManaged)
    {
      cudaGetLastError();
      throw sundials4py::error_returned(
        "Managed CUDA N_Vector replacement requires managed memory");
    }
    return;
  }

  if (device == ArrayDevice::Cuda &&
      (status != cudaSuccess || (attrs.type != cudaMemoryTypeDevice &&
                                 attrs.type != cudaMemoryTypeManaged)))
  {
    cudaGetLastError();
    throw sundials4py::error_returned(
      "Device array pointer is not CUDA-accessible");
  }

  if (status != cudaSuccess) { cudaGetLastError(); }
}

} // namespace

void replace_cuda_array_pointer(N_Vector v, sunrealtype* ptr,
                                ArrayDevice device, nb::object owner)
{
  require_cuda_nvector(v);
  require_compatible_pointer(v, ptr, device);

  auto content = static_cast<N_VectorContent_Cuda>(v->content);
  auto current = device == ArrayDevice::Cpu ? content->host_data
                                            : content->device_data;
  if (current && current->ptr == ptr) { return; }

  const bool managed = N_VIsManagedMemory_Cuda(v);
  auto mem_type      = managed                      ? SUNMEMTYPE_UVM
                       : device == ArrayDevice::Cpu ? SUNMEMTYPE_HOST
                                                    : SUNMEMTYPE_DEVICE;
  auto replacement   = SUNMemoryHelper_Wrap(content->mem_helper, ptr, mem_type);
  if (!replacement)
  {
    throw sundials4py::error_returned(
      "Failed to wrap replacement CUDA array pointer");
  }

  SUNMemory replacement_alias = nullptr;
  if (managed)
  {
    replacement_alias = SUNMemoryHelper_Alias(content->mem_helper, replacement);
    if (!replacement_alias)
    {
      SUNMemoryHelper_Dealloc(content->mem_helper, replacement, cuda_stream(v));
      throw sundials4py::error_returned(
        "Failed to create managed-memory host alias");
    }
  }

  if (owner.ptr()) { prepare_python_array_owners(v); }

  auto stream = cuda_stream(v);
  auto status = cudaStreamSynchronize(
    stream ? *static_cast<const cudaStream_t*>(stream) : cudaStream_t{});
  if (status != cudaSuccess)
  {
    SUNMemoryHelper_Dealloc(content->mem_helper, replacement_alias, stream);
    SUNMemoryHelper_Dealloc(content->mem_helper, replacement, stream);
    throw sundials4py::error_returned(
      "Failed to synchronize before replacing CUDA array pointer");
  }

  SUNMemory old_host   = content->host_data;
  SUNMemory old_device = content->device_data;

  if (managed)
  {
    content->device_data = replacement;
    content->host_data   = replacement_alias;
  }
  else if (device == ArrayDevice::Cpu) { content->host_data = replacement; }
  else { content->device_data = replacement; }

  if (managed)
  {
    if (owner.ptr())
    {
      auto& owners        = prepare_python_array_owners(v);
      owners.host_array   = nb::object();
      owners.device_array = nb::object();
      if (device == ArrayDevice::Cpu) { owners.host_array = std::move(owner); }
      else { owners.device_array = std::move(owner); }
    }
    else
    {
      release_python_array(v, ArrayDevice::Cpu);
      release_python_array(v, ArrayDevice::Cuda);
    }
  }
  else if (owner.ptr())
  {
    if (device == ArrayDevice::Cpu)
    {
      retain_python_host_array(v, std::move(owner));
    }
    else { retain_python_device_array(v, std::move(owner)); }
  }
  else { release_python_array(v, device); }

  SUNErrCode dealloc_status = SUN_SUCCESS;
  if (managed)
  {
    dealloc_status     = SUNMemoryHelper_Dealloc(content->mem_helper, old_host,
                                                 stream);
    auto device_status = SUNMemoryHelper_Dealloc(content->mem_helper,
                                                 old_device, stream);
    if (dealloc_status == SUN_SUCCESS) { dealloc_status = device_status; }
  }
  else
  {
    auto old       = device == ArrayDevice::Cpu ? old_host : old_device;
    dealloc_status = SUNMemoryHelper_Dealloc(content->mem_helper, old, stream);
  }

  if (dealloc_status != SUN_SUCCESS)
  {
    throw sundials4py::error_returned(
      "Failed to release the previous CUDA array pointer");
  }
}

} // namespace nvector_detail

namespace {

using nvector_detail::ArrayDevice;
using nvector_detail::host_array;
using nvector_detail::is_cuda_nvector;
using nvector_detail::parse_device;
using nvector_detail::replace_cuda_array_pointer;
using nvector_detail::require_cuda_nvector;
using nvector_detail::retain_python_device_array;
using nvector_detail::retain_python_host_array;

void check_length(sunindextype vec_length, sundials4py::Array1d data)
{
  if (data.shape(0) != static_cast<size_t>(vec_length))
  {
    throw sundials4py::error_returned(
      "Array shape does not match vector length");
  }
}

void check_length(sunindextype vec_length, CudaArray1d data)
{
  if (data.shape(0) != static_cast<size_t>(vec_length))
  {
    throw sundials4py::error_returned(
      "CUDA array shape does not match vector length");
  }
}

void check_length(sunindextype vec_length, nb::object data)
{
  if (!PyObject_HasAttrString(data.ptr(), "shape"))
  {
    throw sundials4py::error_returned("Device array does not have a shape");
  }

  auto shape = nb::cast<nb::tuple>(data.attr("shape"));
  if (shape.size() != 1 ||
      nb::cast<size_t>(shape[0]) != static_cast<size_t>(vec_length))
  {
    throw sundials4py::error_returned(
      "Device array shape does not match vector length");
  }
}

std::uintptr_t device_array_pointer(nb::object data)
{
  if (!PyObject_HasAttrString(data.ptr(), "unsafe_buffer_pointer"))
  {
    throw sundials4py::error_returned(
      "Device array does not provide unsafe_buffer_pointer()");
  }

  return nb::cast<std::uintptr_t>(data.attr("unsafe_buffer_pointer")());
}

int device_id_for_pointer(void* ptr)
{
  cudaPointerAttributes attrs;
  auto cuda_status = cudaPointerGetAttributes(&attrs, ptr);
  if (cuda_status == cudaSuccess && attrs.device >= 0) { return attrs.device; }

  // Clear the failed cudaPointerGetAttributes status before falling back.
  cudaGetLastError();

  int device_id = 0;
  cuda_status   = cudaGetDevice(&device_id);
  if (cuda_status != cudaSuccess)
  {
    cudaGetLastError();
    return 0;
  }
  return device_id;
}

void sync_to_host(N_Vector v)
{
  if (is_cuda_nvector(v)) { N_VCopyFromDevice_Cuda(v); }
}

struct CudaDeviceArrayPointer
{
  N_Vector x = nullptr;
  nb::object owner;
  int device_id = 0;

  explicit CudaDeviceArrayPointer(N_Vector x) : x(x), owner(nb::find(x))
  {
    auto ptr = N_VGetDeviceArrayPointer(x);
    if (!ptr)
    {
      throw sundials4py::error_returned("Failed to get device array pointer");
    }
    device_id = device_id_for_pointer(ptr);
  }

  std::uintptr_t address() const
  {
    auto ptr = N_VGetDeviceArrayPointer(x);
    if (!ptr)
    {
      throw sundials4py::error_returned("Failed to get device array pointer");
    }
    return reinterpret_cast<std::uintptr_t>(ptr);
  }

  CudaArray1d dlpack() const
  {
    auto ptr = N_VGetDeviceArrayPointer(x);
    if (!ptr)
    {
      throw sundials4py::error_returned("Failed to get device array pointer");
    }
    size_t shape[1]{static_cast<size_t>(N_VGetLength(x))};
    return CudaArray1d(ptr, 1, shape, owner);
  }
};

nb::object device_pointer_object(N_Vector v)
{
  require_cuda_nvector(v);
  return nb::cast(CudaDeviceArrayPointer(v));
}

nb::object get_numpy_array(N_Vector v)
{
  sync_to_host(v);
  return nb::cast(host_array(v));
}

nb::object get_cupy_array(N_Vector v)
{
  require_cuda_nvector(v);
  return nb::module_::import_("cupy").attr("from_dlpack")(
    device_pointer_object(v));
}

nb::object get_torch_tensor(N_Vector v, nb::object device)
{
  auto target = parse_device(device, v);

  auto torch = nb::module_::import_("torch");
  if (target == ArrayDevice::Cpu)
  {
    sync_to_host(v);
    return torch.attr("from_numpy")(nb::cast(host_array(v)));
  }

  require_cuda_nvector(v);
  return torch.attr("utils").attr("dlpack").attr("from_dlpack")(
    device_pointer_object(v));
}

nb::object get_jax_array(N_Vector v, nb::object device)
{
  auto target = parse_device(device, v);

  if (target == ArrayDevice::Cpu)
  {
    sync_to_host(v);
    return nb::module_::import_("jax").attr("dlpack").attr(
      "from_dlpack")(nb::cast(host_array(v)), nb::none(), nb::none());
  }

  require_cuda_nvector(v);
  return nb::module_::import_("jax").attr("dlpack").attr(
    "from_dlpack")(device_pointer_object(v), nb::none(), false);
}

std::shared_ptr<std::remove_pointer_t<N_Vector>> wrap_nvector(N_Vector v)
{
  return our_make_shared<std::remove_pointer_t<N_Vector>, N_VectorDeleter>(v);
}

} // namespace

void bind_nvector_cuda(nb::module_& m)
{
#include "nvector_cuda_generated.hpp"

  nb::class_<CudaDeviceArrayPointer>(m, "_CudaDeviceArrayPointer")
    .def("__int__",
         [](const CudaDeviceArrayPointer& self) { return self.address(); })
    .def("__index__",
         [](const CudaDeviceArrayPointer& self) { return self.address(); })
    .def(
      "__dlpack__",
      [](const CudaDeviceArrayPointer& self, nb::object /* stream */)
      { return self.dlpack(); },
      nb::arg("stream") = nb::none())
    .def("__dlpack_device__", [](const CudaDeviceArrayPointer& self)
         { return nb::make_tuple(2, self.device_id); });

  m.def(
    "N_VNewEmpty_Cuda",
    [](SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    { return wrap_nvector(N_VNewEmpty_Cuda(sunctx)); },
    nb::arg("sunctx"), nb::keep_alive<0, 1>());

  m.def(
    "N_VNew_Cuda",
    [](sunindextype vec_length,
       SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    { return wrap_nvector(N_VNew_Cuda(vec_length, sunctx)); },
    nb::arg("vec_length"), nb::arg("sunctx"), nb::keep_alive<0, 2>());

  m.def(
    "N_VNewManaged_Cuda",
    [](sunindextype vec_length,
       SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    { return wrap_nvector(N_VNewManaged_Cuda(vec_length, sunctx)); },
    nb::arg("vec_length"), nb::arg("sunctx"), nb::keep_alive<0, 2>());

  m.def(
    "N_VNewWithMemHelp_Cuda",
    [](sunindextype vec_length, sunbooleantype use_managed_mem,
       SUNMemoryHelper helper,
       SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      return wrap_nvector(
        N_VNewWithMemHelp_Cuda(vec_length, use_managed_mem, helper, sunctx));
    },
    nb::arg("vec_length"), nb::arg("use_managed_mem"), nb::arg("helper"),
    nb::arg("sunctx"), nb::keep_alive<0, 3>(), nb::keep_alive<0, 4>());

  m.def(
    "N_VMake_Cuda",
    [](sunindextype vec_length, sundials4py::Array1d h_vdata_1d,
       std::uintptr_t d_vdata,
       SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      check_length(vec_length, h_vdata_1d);
      auto v = N_VMake_Cuda(vec_length, h_vdata_1d.data(),
                            reinterpret_cast<sunrealtype*>(d_vdata), sunctx);
      retain_python_host_array(v, h_vdata_1d.cast());
      return wrap_nvector(v);
    },
    nb::arg("vec_length"), nb::arg("h_vdata_1d"), nb::arg("d_vdata"),
    nb::arg("sunctx"), nb::keep_alive<0, 4>());

  m.def(
    "N_VMake_Cuda",
    [](sunindextype vec_length, sundials4py::Array1d h_vdata_1d,
       CudaArray1d d_vdata_1d,
       SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      check_length(vec_length, h_vdata_1d);
      check_length(vec_length, d_vdata_1d);
      auto v = N_VMake_Cuda(vec_length, h_vdata_1d.data(), d_vdata_1d.data(),
                            sunctx);
      retain_python_host_array(v, h_vdata_1d.cast());
      retain_python_device_array(v, d_vdata_1d.cast());
      return wrap_nvector(v);
    },
    nb::arg("vec_length"), nb::arg("h_vdata_1d").noconvert(),
    nb::arg("d_vdata_1d").noconvert(), nb::arg("sunctx"), nb::keep_alive<0, 4>());

  m.def(
    "N_VMake_Cuda",
    [](sunindextype vec_length, sundials4py::Array1d h_vdata_1d,
       nb::object d_vdata_1d,
       SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      check_length(vec_length, h_vdata_1d);
      check_length(vec_length, d_vdata_1d);
      auto v = N_VMake_Cuda(vec_length, h_vdata_1d.data(),
                            reinterpret_cast<sunrealtype*>(
                              device_array_pointer(d_vdata_1d)),
                            sunctx);
      retain_python_host_array(v, h_vdata_1d.cast());
      retain_python_device_array(v, std::move(d_vdata_1d));
      return wrap_nvector(v);
    },
    nb::arg("vec_length"), nb::arg("h_vdata_1d").noconvert(),
    nb::arg("d_vdata_1d"), nb::arg("sunctx"), nb::keep_alive<0, 4>());

  m.def(
    "N_VMakeManaged_Cuda",
    [](sunindextype vec_length, std::uintptr_t vdata,
       SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      return wrap_nvector(
        N_VMakeManaged_Cuda(vec_length, reinterpret_cast<sunrealtype*>(vdata),
                            sunctx));
    },
    nb::arg("vec_length"), nb::arg("vdata"), nb::arg("sunctx"),
    nb::keep_alive<0, 3>());

  m.def(
    "N_VSetHostArrayPointer_Cuda",
    [](sundials4py::Array1d h_vdata_1d, N_Vector v)
    {
      check_length(N_VGetLength(v), h_vdata_1d);
      replace_cuda_array_pointer(v, h_vdata_1d.data(), ArrayDevice::Cpu,
                                 h_vdata_1d.cast());
    },
    nb::arg("h_vdata_1d"), nb::arg("v"),
    "Replace the CUDA N_Vector host pointer with borrowed Python-owned "
    "storage.");

  // Expose the generic N_VSetDeviceArrayPointer to Python with CUDA-aware
  // overloads. This replaces the previous CUDA-specific Python binding
  // N_VSetDeviceArrayPointer_Cuda which is removed intentionally.
  m.def(
    "N_VSetDeviceArrayPointer",
    [](std::uintptr_t d_vdata, N_Vector v)
    {
      replace_cuda_array_pointer(v, reinterpret_cast<sunrealtype*>(d_vdata),
                                 ArrayDevice::Cuda);
    },
    nb::arg("d_vdata"), nb::arg("v"),
    "Replace owned device storage with a borrowed raw pointer. The caller "
    "must keep the allocation alive.");

  m.def(
    "N_VSetDeviceArrayPointer",
    [](CudaArray1d d_vdata_1d, N_Vector v)
    {
      check_length(N_VGetLength(v), d_vdata_1d);
      replace_cuda_array_pointer(v, d_vdata_1d.data(), ArrayDevice::Cuda,
                                 d_vdata_1d.cast());
    },
    nb::arg("d_vdata_1d").noconvert(), nb::arg("v"),
    "Replace owned device storage with a borrowed Python array and retain "
    "its owner.");

  m.def(
    "N_VSetDeviceArrayPointer",
    [](nb::object d_vdata_1d, N_Vector v)
    {
      check_length(N_VGetLength(v), d_vdata_1d);
      replace_cuda_array_pointer(v,
                                 reinterpret_cast<sunrealtype*>(
                                   device_array_pointer(d_vdata_1d)),
                                 ArrayDevice::Cuda, std::move(d_vdata_1d));
    },
    nb::arg("d_vdata_1d"), nb::arg("v"),
    "Replace owned device storage with a borrowed Python array and retain "
    "its owner.");

  m.def(
    "N_VGetHostArrayPointer_Cuda",
    [](N_Vector x)
    {
      auto ptr = N_VGetHostArrayPointer_Cuda(x);
      if (!ptr)
      {
        throw sundials4py::error_returned("Failed to get host array pointer");
      }
      auto owner = nb::find(x);
      size_t shape[1]{static_cast<size_t>(N_VGetLength(x))};
      return sundials4py::Array1d(ptr, 1, shape, owner);
    },
    nb::rv_policy::reference);

  m.def(
    "N_VGetDeviceArrayPointer",
    [](N_Vector x) { return CudaDeviceArrayPointer(x); }, nb::arg("x"));

  m.def("N_VGetNumpyArray", &get_numpy_array, nb::arg("v"));

  m.def("N_VGetJaxArray", &get_jax_array, nb::arg("v"),
        nb::arg("device").none() = nb::none());

  m.def("N_VGetCupyArray", &get_cupy_array, nb::arg("v"));

  m.def("N_VGetTorchTensor", &get_torch_tensor, nb::arg("v"),
        nb::arg("device").none() = nb::none());
}

} // namespace sundials4py
