#include "sundials4py.hpp"

#include <nvector/nvector_serial.h>
#include <sundials/sundials_nvector.hpp>

#include "sundials/sundials_classview.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_nvector_serial(nb::module_& m)
{
#include "nvector_serial_generated.hpp"
}

} // namespace sundials4py
