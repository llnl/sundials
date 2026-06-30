#include "sundials4py.hpp"

#include <cstdint>

#include <nvector/nvector_cuda.h>
#include <sundials/sundials_nvector.hpp>

#include "sundials/sundials_classview.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

namespace {

using CudaArray1d =
  nb::ndarray<sunrealtype, nb::device::cuda, nb::ndim<1>, nb::c_contig>;

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

std::shared_ptr<std::remove_pointer_t<N_Vector>> wrap_nvector(N_Vector v)
{
  return our_make_shared<std::remove_pointer_t<N_Vector>, N_VectorDeleter>(v);
}

} // namespace

void bind_nvector_cuda(nb::module_& m)
{
#include "nvector_cuda_generated.hpp"

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
      return wrap_nvector(
        N_VMake_Cuda(vec_length, h_vdata_1d.data(), d_vdata_1d.data(), sunctx));
    },
    nb::arg("vec_length"), nb::arg("h_vdata_1d").noconvert(),
    nb::arg("d_vdata_1d").noconvert(), nb::arg("sunctx"),
    nb::keep_alive<0, 2>(), nb::keep_alive<0, 3>(), nb::keep_alive<0, 4>());

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
      check_length(N_VGetLength_Cuda(v), h_vdata_1d);
      N_VSetHostArrayPointer_Cuda(h_vdata_1d.data(), v);
    },
    nb::arg("h_vdata_1d"), nb::arg("v"));

  m.def(
    "N_VSetDeviceArrayPointer_Cuda",
    [](std::uintptr_t d_vdata, N_Vector v) {
      N_VSetDeviceArrayPointer_Cuda(reinterpret_cast<sunrealtype*>(d_vdata), v);
    },
    nb::arg("d_vdata"), nb::arg("v"));

  m.def(
    "N_VSetDeviceArrayPointer_Cuda",
    [](CudaArray1d d_vdata_1d, N_Vector v)
    {
      check_length(N_VGetLength_Cuda(v), d_vdata_1d);
      N_VSetDeviceArrayPointer_Cuda(d_vdata_1d.data(), v);
    },
    nb::arg("d_vdata_1d").noconvert(), nb::arg("v"), nb::keep_alive<2, 1>());

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
      size_t shape[1]{static_cast<size_t>(N_VGetLength_Cuda(x))};
      return sundials4py::Array1d(ptr, 1, shape, owner);
    },
    nb::rv_policy::reference);

  m.def(
    "N_VGetDeviceArrayPointer_Cuda",
    [](N_Vector x) -> std::uintptr_t
    {
      auto ptr = N_VGetDeviceArrayPointer_Cuda(x);
      if (!ptr)
      {
        throw sundials4py::error_returned("Failed to get device array pointer");
      }
      return reinterpret_cast<std::uintptr_t>(ptr);
    },
    nb::arg("x"));

  m.def(
    "N_VGetDeviceArray_Cuda",
    [](N_Vector x)
    {
      auto ptr = N_VGetDeviceArrayPointer_Cuda(x);
      if (!ptr)
      {
        throw sundials4py::error_returned("Failed to get device array pointer");
      }
      auto owner = nb::find(x);
      size_t shape[1]{static_cast<size_t>(N_VGetLength_Cuda(x))};
      return CudaArray1d(ptr, 1, shape, owner);
    },
    nb::rv_policy::reference);
}

} // namespace sundials4py
