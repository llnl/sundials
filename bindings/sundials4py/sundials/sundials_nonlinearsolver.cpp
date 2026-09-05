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
 * SUNDIALS SUNNonlinearSolver class. It contains hand-written code 
 * for functions that require special treatment, and includes the
 * generated code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include "sundials/sundials_nonlinearsolver.h"
#include "sundials4py.hpp"

#include <sundials/sundials_nonlinearsolver.hpp>

#include "sundials_nonlinearsolver_usersupplied.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_sunnonlinearsolver(nb::module_& m)
{
  // The abstract Python methods below define the names used by the native
  // operation trampolines; subclasses override only the SUNDIALS operations they
  // support, with solve() validated as mandatory at materialization time.
  nb::class_<CustomSUNNonlinearSolver>(m, "CustomSUNNonlinearSolver",
                                       nb::dynamic_attr())
    .def(nb::init<std::shared_ptr<std::remove_pointer_t<SUNContext>>,
                  SUNNonlinearSolver_Type>(),
         nb::arg("sunctx"), nb::arg("solver_type"))
    .def("_materialization_count",
         &CustomSUNNonlinearSolver::_materialization_count)
    .def_prop_ro("sunctx", &CustomSUNNonlinearSolver::sunctx,
                 nb::sig("def sunctx(self) -> object"),
                 "The SUNDIALS context owned by this object.")
    .def("initialize", [](CustomSUNNonlinearSolver&)
         { return CustomSUNNonlinearSolver::base_method_int("initialize"); })
    .def("setup", [](CustomSUNNonlinearSolver&, nb::object)
         { return CustomSUNNonlinearSolver::base_method_int("setup"); })
    .def("solve", [](CustomSUNNonlinearSolver&, nb::object, nb::object,
                     nb::object, sunrealtype, sunbooleantype)
         { return CustomSUNNonlinearSolver::base_method_int("solve"); })
    .def("set_max_iters", [](CustomSUNNonlinearSolver&, int)
         { return CustomSUNNonlinearSolver::base_method_int("set_max_iters"); })
    .def("get_num_iters", [](CustomSUNNonlinearSolver&)
         { return CustomSUNNonlinearSolver::base_method_int("get_num_iters"); })
    .def("get_cur_iter", [](CustomSUNNonlinearSolver&)
         { return CustomSUNNonlinearSolver::base_method_int("get_cur_iter"); })
    .def("get_num_conv_fails",
         [](CustomSUNNonlinearSolver&) {
           return CustomSUNNonlinearSolver::base_method_int(
             "get_num_conv_fails");
         })
    // The callback setters below must be bound even though no subclass is
    // required to override them: the override check compares a subclass's
    // descriptor against the one registered here, so an unbound name would make
    // that check raise AttributeError and materialization fail outright.
    .def("set_sys_fn", [](CustomSUNNonlinearSolver&, nb::object)
         { return CustomSUNNonlinearSolver::base_method_int("set_sys_fn"); })
    .def("set_sys_fns", [](CustomSUNNonlinearSolver&, nb::object, nb::object)
         { return CustomSUNNonlinearSolver::base_method_int("set_sys_fns"); })
    .def("set_lsetup_fn", [](CustomSUNNonlinearSolver&, nb::object)
         { return CustomSUNNonlinearSolver::base_method_int("set_lsetup_fn"); })
    .def("set_lsolve_fn", [](CustomSUNNonlinearSolver&, nb::object)
         { return CustomSUNNonlinearSolver::base_method_int("set_lsolve_fn"); })
    .def("set_conv_test_fn",
         [](CustomSUNNonlinearSolver&, nb::object) {
           return CustomSUNNonlinearSolver::base_method_int("set_conv_test_fn");
         })
    .def("set_norm_fn", [](CustomSUNNonlinearSolver&, nb::object)
         { return CustomSUNNonlinearSolver::base_method_int("set_norm_fn"); })
    .def("set_get_update_norm_fn",
         [](CustomSUNNonlinearSolver&, nb::object) {
           return CustomSUNNonlinearSolver::base_method_int(
             "set_get_update_norm_fn");
         })
    .def("set_get_conv_rate_fn",
         [](CustomSUNNonlinearSolver&, nb::object) {
           return CustomSUNNonlinearSolver::base_method_int(
             "set_get_conv_rate_fn");
         })
    .def("set_options", [](CustomSUNNonlinearSolver&, const std::string&,
                           const std::string&, const std::vector<std::string>&)
         { return CustomSUNNonlinearSolver::base_method_int("set_options"); });

#include "sundials_nonlinearsolver_generated.hpp"

  m.def(
    "SUNNonlinSolSetOptions",
    [](SUNNonlinearSolver self, const std::string& id,
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

      return SUNNonlinSolSetOptions(self, id.empty() ? nullptr : id.c_str(),
                                    file_name.empty() ? nullptr
                                                      : file_name.c_str(),
                                    argc, argv.data());
    },
    nb::arg("self"), nb::arg("id"), nb::arg("file_name"), nb::arg("argc"),
    nb::arg("args"));

  m.def(
    "SUNNonlinSolSetup",
    [](SUNNonlinearSolver NLS, N_Vector y)
    {
      // Announce that a hand-written binding -- not a SUNDIALS package -- is
      // driving this call, so a custom solver's trampoline records the function
      // table as direct-binding memory rather than mistaking it for integrator
      // memory.
      DirectBindingScope direct;
      return SUNNonlinSolSetup(NLS, y, sunnonlinearsolver_function_table(NLS));
    },
    nb::arg("NLS"), nb::arg("y"));

  m.def(
    "SUNNonlinSolSolve",
    [](SUNNonlinearSolver NLS, N_Vector y0, N_Vector y, N_Vector w,
       sunrealtype tol, sunbooleantype callLSetup)
    {
      DirectBindingScope direct;
      return SUNNonlinSolSolve(NLS, y0, y, w, tol, callLSetup,
                               sunnonlinearsolver_function_table(NLS));
    },
    nb::arg("NLS"), nb::arg("y0"), nb::arg("y"), nb::arg("w"), nb::arg("tol"),
    nb::arg("callLSetup"));

  m.def(
    "SUNNonlinSolSetSysFn",
    [](SUNNonlinearSolver NLS,
       std::function<std::remove_pointer_t<SUNNonlinSolSysFn>> SysFn) -> SUNErrCode
    {
      DirectBindingScope direct;
      auto fntable   = sunnonlinearsolver_function_table(NLS);
      fntable->sysfn = nb::cast(SysFn);
      if (SysFn)
      {
        return SUNNonlinSolSetSysFn(NLS, sunnonlinearsolver_sysfn_wrapper);
      }
      else { return SUNNonlinSolSetSysFn(NLS, nullptr); }
    },
    nb::arg("NLS"), nb::arg("SysFn").none());

  m.def(
    "SUNNonlinSolSetSysFns",
    [](SUNNonlinearSolver NLS,
       std::function<std::remove_pointer_t<SUNNonlinSolSysFn>> RootFn,
       std::function<std::remove_pointer_t<SUNNonlinSolSysFn>> FixedPointFn) -> SUNErrCode
    {
      DirectBindingScope direct;
      auto fntable             = sunnonlinearsolver_function_table(NLS);
      fntable->rootsysfn       = nb::cast(RootFn);
      fntable->fixedpointsysfn = nb::cast(FixedPointFn);
      return SUNNonlinSolSetSysFns(NLS,
                                   RootFn
                                     ? static_cast<SUNNonlinSolSysFn>(
                                         sunnonlinearsolver_rootsysfn_wrapper)
                                     : nullptr,
                                   FixedPointFn
                                     ? static_cast<SUNNonlinSolSysFn>(
                                         sunnonlinearsolver_fixedpointsysfn_wrapper)
                                     : nullptr);
    },
    nb::arg("NLS"), nb::arg("RootFn").none(), nb::arg("FixedPointFn").none());

  m.def(
    "SUNNonlinSolSetLSetupFn",
    [](SUNNonlinearSolver NLS,
       std::function<SUNNonlinSolLSetupStdFn> SetupFn) -> SUNErrCode
    {
      DirectBindingScope direct;
      auto fntable      = sunnonlinearsolver_function_table(NLS);
      fntable->lsetupfn = nb::cast(SetupFn);
      if (SetupFn)
      {
        return SUNNonlinSolSetLSetupFn(NLS, sunnonlinearsolver_lsetupfn_wrapper);
      }
      else { return SUNNonlinSolSetLSetupFn(NLS, nullptr); }
    },
    nb::arg("NLS"), nb::arg("SetupFn").none());

  m.def(
    "SUNNonlinSolSetLSolveFn",
    [](SUNNonlinearSolver NLS,
       std::function<std::remove_pointer_t<SUNNonlinSolLSolveFn>> SolveFn) -> SUNErrCode
    {
      DirectBindingScope direct;
      auto fntable      = sunnonlinearsolver_function_table(NLS);
      fntable->lsolvefn = nb::cast(SolveFn);
      if (SolveFn)
      {
        return SUNNonlinSolSetLSolveFn(NLS, sunnonlinearsolver_lsolvefn_wrapper);
      }
      else { return SUNNonlinSolSetLSolveFn(NLS, nullptr); }
    },
    nb::arg("NLS"), nb::arg("SolveFn").none());

  m.def(
    "SUNNonlinSolSetNormFn",
    [](SUNNonlinearSolver NLS,
       std::function<SUNNonlinSolNormStdFn> NormFn) -> SUNErrCode
    {
      DirectBindingScope direct;
      auto fntable    = sunnonlinearsolver_function_table(NLS);
      fntable->normfn = nb::cast(NormFn);
      if (NormFn)
      {
        return SUNNonlinSolSetNormFn(NLS, sunnonlinearsolver_normfn_wrapper,
                                     fntable);
      }
      else { return SUNNonlinSolSetNormFn(NLS, nullptr, nullptr); }
    },
    nb::arg("NLS"), nb::arg("NormFn").none());

  m.def(
    "SUNNonlinSolSetGetUpdateNormFn",
    [](SUNNonlinearSolver NLS,
       std::function<SUNNonlinSolGetUpdateNormStdFn> GetUpdateNormFn) -> SUNErrCode
    {
      DirectBindingScope direct;
      auto fntable             = sunnonlinearsolver_function_table(NLS);
      fntable->getupdatenormfn = nb::cast(GetUpdateNormFn);
      if (GetUpdateNormFn)
      {
        return SUNNonlinSolSetGetUpdateNormFn(NLS,
                                              sunnonlinearsolver_getupdatenormfn_wrapper,
                                              fntable);
      }
      else { return SUNNonlinSolSetGetUpdateNormFn(NLS, nullptr, nullptr); }
    },
    nb::arg("NLS"), nb::arg("GetUpdateNormFn").none());

  m.def(
    "SUNNonlinSolSetConvTestFn",
    [](SUNNonlinearSolver NLS,
       std::function<std::remove_pointer_t<SUNNonlinSolConvTestFn>> CTestFn) -> SUNErrCode
    {
      DirectBindingScope direct;
      auto fntable        = sunnonlinearsolver_function_table(NLS);
      fntable->convtestfn = nb::cast(CTestFn);
      if (CTestFn)
      {
        return SUNNonlinSolSetConvTestFn(NLS,
                                         sunnonlinearsolver_convtestfn_wrapper,
                                         fntable);
      }
      else { return SUNNonlinSolSetConvTestFn(NLS, nullptr, nullptr); }
    },
    nb::arg("NLS"), nb::arg("CTestFn").none());

  m.def(
    "SUNNonlinSolSetGetConvRateFn",
    [](SUNNonlinearSolver NLS,
       std::function<SUNNonlinSolGetConvRateStdFn> GetConvRateFn) -> SUNErrCode
    {
      DirectBindingScope direct;
      auto fntable           = sunnonlinearsolver_function_table(NLS);
      fntable->getconvratefn = nb::cast(GetConvRateFn);
      if (GetConvRateFn)
      {
        return SUNNonlinSolSetGetConvRateFn(NLS,
                                            sunnonlinearsolver_getconvratefn_wrapper,
                                            fntable);
      }
      else { return SUNNonlinSolSetGetConvRateFn(NLS, nullptr, nullptr); }
    },
    nb::arg("NLS"), nb::arg("GetConvRateFn").none());
}

} // namespace sundials4py

extern "C" void SUNNonlinearSolverFunctionTable_Destroy(void* ptr)
{
  sundials4py::shutdown_safe_delete(
    static_cast<SUNNonlinearSolverFunctionTable*>(ptr));
}
