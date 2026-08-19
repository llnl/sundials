#include "sundials4py.hpp"

#include <nvector/nvector_serial.h>
#include <sundials/sundials_nvector.hpp>

#include "nvector_array_helpers.hpp"
#include "sundials/sundials_classview.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

namespace {

std::shared_ptr<std::remove_pointer_t<N_Vector>> wrap_nvector(N_Vector v)
{
  return our_make_shared<std::remove_pointer_t<N_Vector>, N_VectorDeleter>(v);
}

} // namespace

void bind_nvector_serial(nb::module_& m)
{
#include "nvector_serial_generated.hpp"

  m.def(
    "N_VMake_Serial",
    [](sunindextype vec_length, sundials4py::Array1d data,
       SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      if (data.shape(0) != static_cast<size_t>(vec_length))
      {
        throw sundials4py::error_returned(
          "Array shape does not match vector length");
      }
      auto v = N_VMake_Serial(vec_length, data.data(), sunctx);
      nvector_detail::retain_python_host_array(v, data.cast());
      return wrap_nvector(v);
    },
    nb::arg("vec_length"), nb::arg("data"), nb::arg("sunctx"),
    nb::keep_alive<0, 3>());
}

} // namespace sundials4py
