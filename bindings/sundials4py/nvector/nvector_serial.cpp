#include "sundials4py.hpp"

#include <nvector/nvector_serial.h>
#include <sundials/sundials_nvector.hpp>

#include "nvector_array_helpers.hpp"
#include "sundials/sundials_classview.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_nvector_serial(nb::module_& m)
{
#include "nvector_serial_generated.hpp"

  m.def(
    "N_VMake_Serial",
    [](sunindextype vec_length, sundials4py::Array1d data,
       SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      nvector_detail::require_vector_length(vec_length, data);
      auto v = N_VMake_Serial(vec_length, data.data(), sunctx);
      return nvector_detail::wrap_nvector(v);
    },
    nb::arg("vec_length"), nb::arg("data"), nb::arg("sunctx"),
    // Constructor arguments have the same lifetime as the returned wrapper, so
    // keep_alive avoids the explicit owner registry used by later replacements.
    nb::keep_alive<0, 2>(), nb::keep_alive<0, 3>());
}

} // namespace sundials4py
