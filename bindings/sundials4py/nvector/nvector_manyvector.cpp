#include "sundials4py.hpp"

#include <nvector/nvector_manyvector.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_nvector.hpp>

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_nvector_manyvector(nb::module_& m)
{
#include "nvector_manyvector_generated.hpp"

  m.def(
    "N_VNew_ManyVector",
    [](sunindextype num_subvectors, std::vector<N_Vector> vec_array_1d,
       SUNContext sunctx) -> std::shared_ptr<std::remove_pointer_t<N_Vector>>
    {
      N_Vector* vec_array_1d_ptr = vec_array_1d.empty() ? nullptr
                                                        : vec_array_1d.data();
      return our_make_shared<std::remove_pointer_t<N_Vector>, N_VectorDeleter>(
        N_VNew_ManyVector(num_subvectors, vec_array_1d_ptr, sunctx));
    },
    nb::arg("num_subvectors"), nb::arg("vec_array_1d"), nb::arg("sunctx"),
    nb::keep_alive<0, 3>() /* keep the SUNContext alive as long as the N_Vector is */,
    nb::keep_alive<0, 2>() /* keep the list, and thus the elements, alive as long as the N_Vector is */);
}

} // namespace sundials4py
