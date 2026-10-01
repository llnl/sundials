/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *                Daniel R. Reynolds @ UMBC
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
 * SUNDIALS SUNMatrix class. It contains hand-written code for functions
 * that require special treatment, and includes the generated code
 * produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <sundials/sundials_matrix.hpp>

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_sunmatrix(nb::module_& m)
{
  // CustomSUNMatrix is a Python-side base class. Its default methods raise
  // NotImplementedError, while subclasses provide the implementation that the
  // native SUNMatrix vtable calls through sundials4py_types.hpp.
  nb::class_<CustomSUNMatrix>(m, "CustomSUNMatrix", nb::dynamic_attr())
    .def(nb::init<std::shared_ptr<std::remove_pointer_t<SUNContext>>>(),
         nb::arg("sunctx"))
    // This is intentionally private-facing: tests use it to verify that
    // transparent conversion materializes the native handle exactly once.
    .def("_materialization_count", &CustomSUNMatrix::_materialization_count)
    .def_prop_ro("sunctx", &CustomSUNMatrix::sunctx,
                 nb::sig("def sunctx(self) -> object"),
                 "The SUNDIALS context owned by this object.")
    .def("clone", [](CustomSUNMatrix&)
         { return CustomSUNMatrix::base_method_status("clone"); })
    .def("zero", [](CustomSUNMatrix&)
         { return CustomSUNMatrix::base_method_status("zero"); })
    .def("copy", [](CustomSUNMatrix&, nb::object)
         { return CustomSUNMatrix::base_method_status("copy"); })
    .def("scaleadd", [](CustomSUNMatrix&, sunrealtype, nb::object)
         { return CustomSUNMatrix::base_method_status("scaleadd"); })
    .def("scaleaddi", [](CustomSUNMatrix&, sunrealtype)
         { return CustomSUNMatrix::base_method_status("scaleaddi"); })
    .def("matvec", [](CustomSUNMatrix&, nb::object, nb::object)
         { return CustomSUNMatrix::base_method_status("matvec"); })
    .def("matvecsetup", [](CustomSUNMatrix&)
         { return CustomSUNMatrix::base_method_status("matvecsetup"); })
    .def("hermitian_transpose_matvec",
         [](CustomSUNMatrix&, nb::object, nb::object) {
           return CustomSUNMatrix::base_method_status(
             "hermitian_transpose_matvec");
         });

  m.def(
    "SUNMatClone",
    [](SUNMatrix A) -> nb::object
    {
      SUNMatrix clone = SUNMatClone(A);
      if (!clone) { return nb::none(); }
      return nb::cast(
        our_make_shared<std::remove_pointer_t<SUNMatrix>, SUNMatrixDeleter>(
          clone));
    },
    nb::arg("A"));

#include "sundials_matrix_generated.hpp"
}

} // namespace sundials4py
