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
 * SUNDIALS SUNAdaptController class. It contains hand-written code for
 * functions that require special treatment, and includes the generated
 * code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include "sundials/sundials_adaptcontroller.h"
#include "sundials4py.hpp"

#include <sundials/sundials_adaptcontroller.hpp>

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_sunadaptcontroller(nb::module_& m)
{
  // The shared base class owns the lazy native-handle machinery. The concrete
  // H and MRI subclasses below select the required estimate callback shape.
  nb::class_<CustomSUNAdaptController>(m, "CustomSUNAdaptController",
                                       nb::dynamic_attr())
    .def("_materialization_count",
         &CustomSUNAdaptController::_materialization_count)
    .def_prop_ro("sunctx", &CustomSUNAdaptController::sunctx,
                 nb::sig("def sunctx(self) -> object"),
                 "The SUNDIALS context owned by this object.")
    .def("estimate_step",
         [](CustomSUNAdaptController&, sunrealtype, int, sunrealtype) {
           return CustomSUNAdaptController::base_method_status("estimate_step");
         })
    .def("estimate_step_tol",
         [](CustomSUNAdaptController&, sunrealtype, sunrealtype, int,
            sunrealtype, sunrealtype) {
           return CustomSUNAdaptController::base_method_status(
             "estimate_step_tol");
         })
    .def("reset", [](CustomSUNAdaptController&)
         { return CustomSUNAdaptController::base_method_status("reset"); })
    .def("set_defaults", [](CustomSUNAdaptController&)
         { return CustomSUNAdaptController::base_method_status("set_defaults"); })
    .def("set_error_bias",
         [](CustomSUNAdaptController&, sunrealtype) {
           return CustomSUNAdaptController::base_method_status(
             "set_error_bias");
         })
    .def("update_h", [](CustomSUNAdaptController&, sunrealtype, sunrealtype)
         { return CustomSUNAdaptController::base_method_status("update_h"); })
    .def("update_mri_h_tol",
         [](CustomSUNAdaptController&, sunrealtype, sunrealtype, sunrealtype,
            sunrealtype) {
           return CustomSUNAdaptController::base_method_status(
             "update_mri_h_tol");
         });

  nb::class_<CustomSUNHController, CustomSUNAdaptController>(m, "CustomSUNHController",
                                                             nb::dynamic_attr())
    .def(nb::init<std::shared_ptr<std::remove_pointer_t<SUNContext>>>(),
         nb::arg("sunctx"));

  nb::class_<CustomSUNMRIController, CustomSUNAdaptController>(m, "CustomSUNMRIController",
                                                               nb::dynamic_attr())
    .def(nb::init<std::shared_ptr<std::remove_pointer_t<SUNContext>>>(),
         nb::arg("sunctx"));

#include "sundials_adaptcontroller_generated.hpp"

  m.def(
    "SUNAdaptController_SetOptions",
    [](SUNAdaptController self, const std::string& id,
       const std::string& file_name, int argc,
       const std::vector<std::string>& args)
    {
      std::vector<char*> argv;
      argv.reserve(args.size());

      for (const auto& arg : args)
      {
        // We need a non-const char*, so we use data() and an explicit cast.
        // This is safe as long as the underlying std::string is not modified.
        argv.push_back(const_cast<char*>(arg.data()));
      }

      return SUNAdaptController_SetOptions(self,
                                           id.empty() ? nullptr : id.c_str(),
                                           file_name.empty() ? nullptr
                                                             : file_name.c_str(),
                                           argc, argv.data());
    },
    nb::arg("self"), nb::arg("id"), nb::arg("file_name"), nb::arg("argc"),
    nb::arg("args"));
}

} // namespace sundials4py
