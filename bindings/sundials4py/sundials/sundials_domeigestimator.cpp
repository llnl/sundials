/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * SUNDIALS N_Vector class. It contains hand-written code for functions
 * that require special treatment, and includes the generated code
 * produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include <cstdlib>
#include <cstring>

#include "sundials4py.hpp"

#include <sundials/sundials_classview.hpp>
#include <sundials/sundials_domeigestimator.hpp>

#include "sundials_domeigestimator_usersupplied.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_sundomeigestimator(nb::module_& m)
{
  nb::class_<CustomSUNDomEigEstimator>(m, "CustomSUNDomEigEstimator",
                                       nb::dynamic_attr(),
                                       nb::is_weak_referenceable())
    .def(nb::init<std::shared_ptr<std::remove_pointer_t<SUNContext>>>(),
         nb::arg("sunctx"))
    .def("_is_materialized", &CustomSUNDomEigEstimator::_is_materialized)
    .def("validate",
         [](CustomSUNDomEigEstimator& self) { self.validate(nb::find(&self)); })
    .def_prop_ro("sunctx", &CustomSUNDomEigEstimator::sunctx,
                 nb::sig("def sunctx(self) -> object"))
    .def("set_atimes", [](CustomSUNDomEigEstimator&, nb::object)
         { return CustomSUNDomEigEstimator::base_method_status("set_atimes"); })
    .def("set_rhs", [](CustomSUNDomEigEstimator&, nb::object)
         { return CustomSUNDomEigEstimator::base_method_status("set_rhs"); })
    .def("set_rhs_linearization_point",
         [](CustomSUNDomEigEstimator&, sunrealtype, N_Vector)
         {
           return CustomSUNDomEigEstimator::base_method_status(
             "set_rhs_linearization_point");
         })
    .def("set_options", [](CustomSUNDomEigEstimator&, const std::string&,
                           const std::string&, const std::vector<std::string>&)
         { return CustomSUNDomEigEstimator::base_method_status("set_options"); })
    .def("set_max_iters",
         [](CustomSUNDomEigEstimator&, long int) {
           return CustomSUNDomEigEstimator::base_method_status("set_max_iters");
         })
    .def("set_num_preprocess_iters",
         [](CustomSUNDomEigEstimator&, int)
         {
           return CustomSUNDomEigEstimator::base_method_status(
             "set_num_preprocess_iters");
         })
    .def("set_rel_tol", [](CustomSUNDomEigEstimator&, sunrealtype)
         { return CustomSUNDomEigEstimator::base_method_status("set_rel_tol"); })
    .def("set_initial_guess",
         [](CustomSUNDomEigEstimator&, N_Vector) {
           return CustomSUNDomEigEstimator::base_method_status(
             "set_initial_guess");
         })
    .def("initialize", [](CustomSUNDomEigEstimator&)
         { return CustomSUNDomEigEstimator::base_method_status("initialize"); })
    .def("estimate", [](CustomSUNDomEigEstimator&)
         { return CustomSUNDomEigEstimator::base_method_status("estimate"); })
    .def("get_res", [](CustomSUNDomEigEstimator&)
         { return CustomSUNDomEigEstimator::base_method_status("get_res"); })
    .def("get_num_iters",
         [](CustomSUNDomEigEstimator&) {
           return CustomSUNDomEigEstimator::base_method_status("get_num_iters");
         })
    .def("get_num_rhs_evals",
         [](CustomSUNDomEigEstimator&) {
           return CustomSUNDomEigEstimator::base_method_status(
             "get_num_rhs_evals");
         })
    .def("get_num_atimes_calls",
         [](CustomSUNDomEigEstimator&)
         {
           return CustomSUNDomEigEstimator::base_method_status(
             "get_num_atimes_calls");
         })
    .def("write", [](CustomSUNDomEigEstimator&, std::uintptr_t)
         { return CustomSUNDomEigEstimator::base_method_status("write"); });

#include "sundials_domeigestimator_generated.hpp"

  sundials4py::scoped_def(
    m, "SUNDomEigEstimator_SetOptions",
    [](SUNDomEigEstimator self, const std::string& id,
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

      return SUNDomEigEstimator_SetOptions(self,
                                           id.empty() ? nullptr : id.c_str(),
                                           file_name.empty() ? nullptr
                                                             : file_name.c_str(),
                                           argc, argv.data());
    },
    nb::arg("self"), nb::arg("id"), nb::arg("file_name"), nb::arg("argc"),
    nb::arg("args"));

  sundials4py::scoped_def(
    m, "SUNDomEigEstimator_SetATimes",
    [](SUNDomEigEstimator dee,
       std::function<std::remove_pointer_t<SUNATimesFn>> ATimes) -> SUNErrCode
    {
      auto fntable      = domeigestimator_function_table(dee);
      nb::object old_fn = fntable->atimes;
      fntable->atimes   = nb::cast(ATimes);
      SUNErrCode status;
      if (ATimes)
      {
        status = SUNDomEigEstimator_SetATimes(dee, fntable,
                                              sundomeigestimator_atimes_wrapper);
      }
      else { status = SUNDomEigEstimator_SetATimes(dee, fntable, nullptr); }
      if (status != SUN_SUCCESS) { fntable->atimes = std::move(old_fn); }
      return status;
    },
    nb::arg("DEE"), nb::arg("ATimes").none());

  sundials4py::scoped_def(
    m, "SUNDomEigEstimator_SetRhs",
    [](SUNDomEigEstimator DEE,
       std::function<std::remove_pointer_t<SUNRhsFn>> RHSfn) -> SUNErrCode
    {
      auto fntable      = domeigestimator_function_table(DEE);
      nb::object old_fn = fntable->deerhs;
      fntable->deerhs   = nb::cast(RHSfn);
      SUNErrCode status;
      if (RHSfn)
      {
        status = SUNDomEigEstimator_SetRhs(DEE, fntable,
                                           sundomeigestimator_setrhs_wrapper);
      }
      else { status = SUNDomEigEstimator_SetRhs(DEE, fntable, nullptr); }
      if (status != SUN_SUCCESS) { fntable->deerhs = std::move(old_fn); }
      return status;
    },
    nb::arg("DEE"), nb::arg("RHSfn").none());
}

} // namespace sundials4py

extern "C" void SUNDomEigEstimatorFunctionTable_Destroy(void* ptr)
{
  sundials4py::shutdown_safe_delete(
    static_cast<SUNDomEigEstimatorFunctionTable*>(ptr));
}
