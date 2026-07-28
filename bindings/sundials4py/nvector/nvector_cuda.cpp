#include "sundials4py.hpp"

#include <cstdint>
#include <unordered_map>

#include <nvector/nvector_cuda.h>
#include <sundials/sundials_nvector.hpp>

#include "sundials/sundials_classview.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

namespace {

using CudaArray1d =
  nb::ndarray<sunrealtype, nb::device::cuda, nb::ndim<1>, nb::c_contig>;

using NVectorDestroyFn = void (*)(N_Vector);

struct PythonDeviceArrayOwner
{
  nb::object device_array;
  NVectorDestroyFn destroy = nullptr;
};

std::unordered_map<N_Vector, PythonDeviceArrayOwner>& python_device_array_owners()
{
  static auto* owners = new std::unordered_map<N_Vector, PythonDeviceArrayOwner>;
  return *owners;
}

void destroy_cuda_nvector_with_python_device_array(N_Vector v)
{
  NVectorDestroyFn destroy = N_VDestroy_Cuda;
  nb::object device_array;

  {
    nb::gil_scoped_acquire gil;
    auto& owners = python_device_array_owners();
    auto it      = owners.find(v);
    if (it != owners.end())
    {
      destroy      = it->second.destroy;
      device_array = std::move(it->second.device_array);
      owners.erase(it);
    }
  }

  destroy(v);

  if (device_array.ptr())
  {
    nb::gil_scoped_acquire gil;
    device_array = nb::object();
  }
}

void retain_python_device_array(N_Vector v, nb::object device_array)
{
  if (!v || !device_array.ptr()) { return; }

  auto& owners        = python_device_array_owners();
  auto [it, inserted] = owners.try_emplace(v);
  if (inserted)
  {
    it->second.destroy = v->ops->nvdestroy;
    v->ops->nvdestroy  = destroy_cuda_nvector_with_python_device_array;
  }
  it->second.device_array = std::move(device_array);
}

void release_python_device_array(N_Vector v)
{
  auto& owners = python_device_array_owners();
  auto it      = owners.find(v);
  if (it == owners.end()) { return; }

  if (v && v->ops &&
      v->ops->nvdestroy == destroy_cuda_nvector_with_python_device_array)
  {
    v->ops->nvdestroy = it->second.destroy;
  }
  owners.erase(it);
}

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

enum class ArrayDevice
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

bool object_is_none(nb::object obj) { return obj.ptr() == Py_None; }

bool is_cuda_nvector(N_Vector v)
{
  return N_VGetVectorID(v) == SUNDIALS_NVEC_CUDA;
}

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

ArrayDevice parse_device(nb::object device, N_Vector v)
{
  auto value = optional_string(device, "device");
  if (value.empty())
  {
    return is_cuda_nvector(v) ? ArrayDevice::Cuda : ArrayDevice::Cpu;
  }
  if (value == "cpu" || value == "host") { return ArrayDevice::Cpu; }
  if (value == "cuda") { return ArrayDevice::Cuda; }

  throw sundials4py::error_returned(
    "device must be 'cpu', 'host', 'cuda', or None");
}

CopyFrom parse_copy_from(nb::object copy_from)
{
  auto value = optional_string(copy_from, "copy_from");
  if (value.empty()) { return CopyFrom::None; }
  if (value == "cpu") { return CopyFrom::Cpu; }
  if (value == "device") { return CopyFrom::Device; }

  throw sundials4py::error_returned(
    "copy_from must be 'cpu', 'device', or None");
}

void require_cuda_nvector(N_Vector v)
{
  if (!is_cuda_nvector(v))
  {
    throw sundials4py::error_returned(
      "CUDA array access requires a CUDA N_Vector");
  }
}

void apply_copy_from(N_Vector v, ArrayDevice device, CopyFrom copy_from)
{
  if (copy_from == CopyFrom::None) { return; }

  if (copy_from == CopyFrom::Cpu)
  {
    if (device == ArrayDevice::Cpu) { return; }
    require_cuda_nvector(v);
    N_VCopyToDevice_Cuda(v);
    return;
  }

  if (device == ArrayDevice::Cuda) { return; }
  require_cuda_nvector(v);
  N_VCopyFromDevice_Cuda(v);
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

nb::object get_numpy_array(N_Vector v, nb::object device, nb::object copy_from)
{
  auto target = parse_device(device, v);
  if (target == ArrayDevice::Cuda)
  {
    throw sundials4py::error_returned(
      "N_VGetNumpyArray does not support device='cuda'");
  }
  apply_copy_from(v, target, parse_copy_from(copy_from));
  return nb::cast(host_array(v));
}

nb::object get_cupy_array(N_Vector v, nb::object device, nb::object copy_from)
{
  auto target = parse_device(device, v);
  if (target == ArrayDevice::Cpu)
  {
    throw sundials4py::error_returned(
      "N_VGetCupyArray does not support device='cpu'");
  }
  apply_copy_from(v, target, parse_copy_from(copy_from));
  return nb::module_::import_("cupy").attr("from_dlpack")(
    device_pointer_object(v));
}

nb::object get_torch_tensor(N_Vector v, nb::object device, nb::object copy_from)
{
  auto target = parse_device(device, v);
  apply_copy_from(v, target, parse_copy_from(copy_from));

  auto torch = nb::module_::import_("torch");
  if (target == ArrayDevice::Cpu)
  {
    return torch.attr("from_numpy")(nb::cast(host_array(v)));
  }

  return torch.attr("utils").attr("dlpack").attr("from_dlpack")(
    device_pointer_object(v));
}

nb::object get_jax_array(N_Vector v, nb::object device, nb::object copy_from)
{
  auto target = parse_device(device, v);
  apply_copy_from(v, target, parse_copy_from(copy_from));

  if (target == ArrayDevice::Cpu)
  {
    return nb::module_::import_("jax").attr("dlpack").attr(
      "from_dlpack")(nb::cast(host_array(v)), nb::none(), nb::none());
  }

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
      return wrap_nvector(N_VMake_Cuda(vec_length, h_vdata_1d.data(),
                                       reinterpret_cast<sunrealtype*>(d_vdata),
                                       sunctx));
    },
    nb::arg("vec_length"), nb::arg("h_vdata_1d"), nb::arg("d_vdata"),
    nb::arg("sunctx"), nb::keep_alive<0, 2>(), nb::keep_alive<0, 4>());

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
      retain_python_device_array(v, d_vdata_1d.cast());
      return wrap_nvector(v);
    },
    nb::arg("vec_length"), nb::arg("h_vdata_1d").noconvert(),
    nb::arg("d_vdata_1d").noconvert(), nb::arg("sunctx"),
    nb::keep_alive<0, 2>(), nb::keep_alive<0, 3>(), nb::keep_alive<0, 4>());

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
      retain_python_device_array(v, std::move(d_vdata_1d));
      return wrap_nvector(v);
    },
    nb::arg("vec_length"), nb::arg("h_vdata_1d").noconvert(),
    nb::arg("d_vdata_1d"), nb::arg("sunctx"), nb::keep_alive<0, 2>(),
    nb::keep_alive<0, 3>(), nb::keep_alive<0, 4>());

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
      N_VSetHostArrayPointer_Cuda(h_vdata_1d.data(), v);
    },
    nb::arg("h_vdata_1d"), nb::arg("v"));

  m.def(
    "N_VSetDeviceArrayPointer_Cuda",
    [](std::uintptr_t d_vdata, N_Vector v)
    {
      release_python_device_array(v);
      N_VSetDeviceArrayPointer_Cuda(reinterpret_cast<sunrealtype*>(d_vdata), v);
    },
    nb::arg("d_vdata"), nb::arg("v"));

  m.def(
    "N_VSetDeviceArrayPointer_Cuda",
    [](CudaArray1d d_vdata_1d, N_Vector v)
    {
      check_length(N_VGetLength(v), d_vdata_1d);
      N_VSetDeviceArrayPointer_Cuda(d_vdata_1d.data(), v);
      retain_python_device_array(v, d_vdata_1d.cast());
    },
    nb::arg("d_vdata_1d").noconvert(), nb::arg("v"), nb::keep_alive<2, 1>());

  m.def(
    "N_VSetDeviceArrayPointer_Cuda",
    [](nb::object d_vdata_1d, N_Vector v)
    {
      check_length(N_VGetLength(v), d_vdata_1d);
      N_VSetDeviceArrayPointer_Cuda(reinterpret_cast<sunrealtype*>(
                                      device_array_pointer(d_vdata_1d)),
                                    v);
      retain_python_device_array(v, std::move(d_vdata_1d));
    },
    nb::arg("d_vdata_1d"), nb::arg("v"), nb::keep_alive<2, 1>());

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
}

} // namespace sundials4py
