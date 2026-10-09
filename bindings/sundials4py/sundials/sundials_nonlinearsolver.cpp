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
                                       nb::dynamic_attr(),
                                       nb::is_weak_referenceable())
    .def(nb::init<std::shared_ptr<std::remove_pointer_t<SUNContext>>, int>(),
         nb::arg("sunctx"), nb::arg("solver_type"))
    .def("_is_materialized", &CustomSUNNonlinearSolver::_is_materialized)
    .def("validate",
         [](CustomSUNNonlinearSolver& self) { self.validate(nb::find(&self)); })
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

  sundials4py::scoped_def(
    m, "SUNNonlinSolSetOptions",
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

  sundials4py::scoped_def(
    m, "SUNNonlinSolSetup",
    [](SUNNonlinearSolver NLS, N_Vector y)
    {
      if (nb::object impl = CustomSUNNonlinearSolver::_python_object_for(NLS);
          impl.is_valid())
      {
        return nb::cast<int>(
          impl.attr("setup")(nb::cast(y, nb::rv_policy::reference)));
      }
      return SUNNonlinSolSetup(NLS, y, sunnonlinearsolver_function_table(NLS));
    },
    nb::arg("NLS"), nb::arg("y"));

  sundials4py::scoped_def(
    m, "SUNNonlinSolSolve",
    [](SUNNonlinearSolver NLS, N_Vector y0, N_Vector y, N_Vector w,
       sunrealtype tol, sunbooleantype callLSetup)
    {
      if (nb::object impl = CustomSUNNonlinearSolver::_python_object_for(NLS);
          impl.is_valid())
      {
        return nb::cast<int>(impl.attr(
          "solve")(nb::cast(y0, nb::rv_policy::reference),
                   nb::cast(y, nb::rv_policy::reference),
                   nb::cast(w, nb::rv_policy::reference), tol, callLSetup));
      }
      return SUNNonlinSolSolve(NLS, y0, y, w, tol, callLSetup,
                               sunnonlinearsolver_function_table(NLS));
    },
    nb::arg("NLS"), nb::arg("y0"), nb::arg("y"), nb::arg("w"), nb::arg("tol"),
    nb::arg("callLSetup"));

  sundials4py::scoped_def(
    m, "SUNNonlinSolSetSysFn",
    [](SUNNonlinearSolver NLS,
       std::function<std::remove_pointer_t<SUNNonlinSolSysFn>> SysFn) -> SUNErrCode
    {
      if (nb::object impl = CustomSUNNonlinearSolver::_python_object_for(NLS);
          impl.is_valid())
      {
        nb::object adapter = nb::none();
        if (SysFn)
        {
          adapter = nb::cpp_function([SysFn](N_Vector y, N_Vector F)
                                     { return SysFn(y, F, nullptr); },
                                     nb::arg("y"), nb::arg("F"));
        }
        return nb::cast<int>(impl.attr("set_sys_fn")(adapter));
      }
      auto fntable      = sunnonlinearsolver_function_table(NLS);
      nb::object old_fn = fntable->sysfn;
      fntable->sysfn    = nb::cast(SysFn);
      SUNErrCode status;
      if (SysFn)
      {
        status = SUNNonlinSolSetSysFn(NLS, sunnonlinearsolver_sysfn_wrapper);
      }
      else { status = SUNNonlinSolSetSysFn(NLS, nullptr); }
      if (status != SUN_SUCCESS) { fntable->sysfn = std::move(old_fn); }
      return status;
    },
    nb::arg("NLS"), nb::arg("SysFn").none());

  sundials4py::scoped_def(
    m, "SUNNonlinSolSetSysFns",
    [](SUNNonlinearSolver NLS,
       std::function<std::remove_pointer_t<SUNNonlinSolSysFn>> RootFn,
       std::function<std::remove_pointer_t<SUNNonlinSolSysFn>> FixedPointFn) -> SUNErrCode
    {
      if (nb::object impl = CustomSUNNonlinearSolver::_python_object_for(NLS);
          impl.is_valid())
      {
        nb::object root        = nb::none();
        nb::object fixed_point = nb::none();
        if (RootFn)
        {
          root = nb::cpp_function([RootFn](N_Vector y, N_Vector F)
                                  { return RootFn(y, F, nullptr); },
                                  nb::arg("y"), nb::arg("F"));
        }
        if (FixedPointFn)
        {
          fixed_point = nb::cpp_function([FixedPointFn](N_Vector y, N_Vector F)
                                         { return FixedPointFn(y, F, nullptr); },
                                         nb::arg("y"), nb::arg("F"));
        }
        return nb::cast<int>(impl.attr("set_sys_fns")(root, fixed_point));
      }
      auto fntable             = sunnonlinearsolver_function_table(NLS);
      nb::object old_root      = fntable->rootsysfn;
      nb::object old_fixed     = fntable->fixedpointsysfn;
      fntable->rootsysfn       = nb::cast(RootFn);
      fntable->fixedpointsysfn = nb::cast(FixedPointFn);
      SUNErrCode status =
        SUNNonlinSolSetSysFns(NLS,
                              RootFn ? static_cast<SUNNonlinSolSysFn>(
                                         sunnonlinearsolver_rootsysfn_wrapper)
                                     : nullptr,
                              FixedPointFn
                                ? static_cast<SUNNonlinSolSysFn>(
                                    sunnonlinearsolver_fixedpointsysfn_wrapper)
                                : nullptr);
      if (status != SUN_SUCCESS)
      {
        fntable->rootsysfn       = std::move(old_root);
        fntable->fixedpointsysfn = std::move(old_fixed);
      }
      return status;
    },
    nb::arg("NLS"), nb::arg("RootFn").none(), nb::arg("FixedPointFn").none());

  sundials4py::scoped_def(
    m, "SUNNonlinSolSetLSetupFn",
    [](SUNNonlinearSolver NLS,
       std::function<SUNNonlinSolLSetupStdFn> SetupFn) -> SUNErrCode
    {
      if (nb::object impl = CustomSUNNonlinearSolver::_python_object_for(NLS);
          impl.is_valid())
      {
        nb::object adapter = nb::none();
        if (SetupFn)
        {
          adapter = nb::cpp_function([SetupFn](sunbooleantype jbad)
                                     { return SetupFn(jbad, nullptr); },
                                     nb::arg("jbad"));
        }
        return nb::cast<int>(impl.attr("set_lsetup_fn")(adapter));
      }
      auto fntable      = sunnonlinearsolver_function_table(NLS);
      nb::object old_fn = fntable->lsetupfn;
      fntable->lsetupfn = nb::cast(SetupFn);
      SUNErrCode status;
      if (SetupFn)
      {
        status = SUNNonlinSolSetLSetupFn(NLS,
                                         sunnonlinearsolver_lsetupfn_wrapper);
      }
      else { status = SUNNonlinSolSetLSetupFn(NLS, nullptr); }
      if (status != SUN_SUCCESS) { fntable->lsetupfn = std::move(old_fn); }
      return status;
    },
    nb::arg("NLS"), nb::arg("SetupFn").none());

  sundials4py::scoped_def(
    m, "SUNNonlinSolSetLSolveFn",
    [](SUNNonlinearSolver NLS,
       std::function<std::remove_pointer_t<SUNNonlinSolLSolveFn>> SolveFn) -> SUNErrCode
    {
      if (nb::object impl = CustomSUNNonlinearSolver::_python_object_for(NLS);
          impl.is_valid())
      {
        nb::object adapter = nb::none();
        if (SolveFn)
        {
          adapter = nb::cpp_function([SolveFn](N_Vector b)
                                     { return SolveFn(b, nullptr); },
                                     nb::arg("b"));
        }
        return nb::cast<int>(impl.attr("set_lsolve_fn")(adapter));
      }
      auto fntable      = sunnonlinearsolver_function_table(NLS);
      nb::object old_fn = fntable->lsolvefn;
      fntable->lsolvefn = nb::cast(SolveFn);
      SUNErrCode status;
      if (SolveFn)
      {
        status = SUNNonlinSolSetLSolveFn(NLS,
                                         sunnonlinearsolver_lsolvefn_wrapper);
      }
      else { status = SUNNonlinSolSetLSolveFn(NLS, nullptr); }
      if (status != SUN_SUCCESS) { fntable->lsolvefn = std::move(old_fn); }
      return status;
    },
    nb::arg("NLS"), nb::arg("SolveFn").none());

  sundials4py::scoped_def(
    m, "SUNNonlinSolSetNormFn",
    [](SUNNonlinearSolver NLS,
       std::function<SUNNonlinSolNormStdFn> NormFn) -> SUNErrCode
    {
      if (nb::object impl = CustomSUNNonlinearSolver::_python_object_for(NLS);
          impl.is_valid())
      {
        nb::object adapter = nb::none();
        if (NormFn)
        {
          adapter = nb::cpp_function([NormFn](N_Vector delta, N_Vector w)
                                     { return NormFn(delta, w, nullptr); },
                                     nb::arg("delta"), nb::arg("w"));
        }
        return nb::cast<int>(impl.attr("set_norm_fn")(adapter));
      }
      auto fntable      = sunnonlinearsolver_function_table(NLS);
      nb::object old_fn = fntable->normfn;
      fntable->normfn   = nb::cast(NormFn);
      SUNErrCode status;
      if (NormFn)
      {
        status = SUNNonlinSolSetNormFn(NLS, sunnonlinearsolver_normfn_wrapper,
                                       fntable);
      }
      else { status = SUNNonlinSolSetNormFn(NLS, nullptr, nullptr); }
      if (status != SUN_SUCCESS) { fntable->normfn = std::move(old_fn); }
      return status;
    },
    nb::arg("NLS"), nb::arg("NormFn").none());

  sundials4py::scoped_def(
    m, "SUNNonlinSolSetGetUpdateNormFn",
    [](SUNNonlinearSolver NLS,
       std::function<SUNNonlinSolGetUpdateNormStdFn> GetUpdateNormFn) -> SUNErrCode
    {
      if (nb::object impl = CustomSUNNonlinearSolver::_python_object_for(NLS);
          impl.is_valid())
      {
        nb::object adapter = nb::none();
        if (GetUpdateNormFn)
        {
          adapter = nb::cpp_function([GetUpdateNormFn]()
                                     { return GetUpdateNormFn(nullptr); });
        }
        return nb::cast<int>(impl.attr("set_get_update_norm_fn")(adapter));
      }
      auto fntable             = sunnonlinearsolver_function_table(NLS);
      nb::object old_fn        = fntable->getupdatenormfn;
      fntable->getupdatenormfn = nb::cast(GetUpdateNormFn);
      SUNErrCode status;
      if (GetUpdateNormFn)
      {
        status =
          SUNNonlinSolSetGetUpdateNormFn(NLS,
                                         sunnonlinearsolver_getupdatenormfn_wrapper,
                                         fntable);
      }
      else { status = SUNNonlinSolSetGetUpdateNormFn(NLS, nullptr, nullptr); }
      if (status != SUN_SUCCESS)
      {
        fntable->getupdatenormfn = std::move(old_fn);
      }
      return status;
    },
    nb::arg("NLS"), nb::arg("GetUpdateNormFn").none());

  sundials4py::scoped_def(
    m, "SUNNonlinSolSetConvTestFn",
    [](SUNNonlinearSolver NLS,
       std::function<std::remove_pointer_t<SUNNonlinSolConvTestFn>> CTestFn) -> SUNErrCode
    {
      if (nb::object impl = CustomSUNNonlinearSolver::_python_object_for(NLS);
          impl.is_valid())
      {
        nb::object adapter = nb::none();
        if (CTestFn)
        {
          adapter = nb::
            cpp_function([NLS, CTestFn](N_Vector y, N_Vector delta,
                                        sunrealtype tol, N_Vector ewt)
                         { return CTestFn(NLS, y, delta, tol, ewt, nullptr); },
                         nb::arg("y"), nb::arg("delta"), nb::arg("tol"),
                         nb::arg("ewt"));
        }
        return nb::cast<int>(impl.attr("set_conv_test_fn")(adapter));
      }
      auto fntable        = sunnonlinearsolver_function_table(NLS);
      nb::object old_fn   = fntable->convtestfn;
      fntable->convtestfn = nb::cast(CTestFn);
      SUNErrCode status;
      if (CTestFn)
      {
        status = SUNNonlinSolSetConvTestFn(NLS,
                                           sunnonlinearsolver_convtestfn_wrapper,
                                           fntable);
      }
      else { status = SUNNonlinSolSetConvTestFn(NLS, nullptr, nullptr); }
      if (status != SUN_SUCCESS) { fntable->convtestfn = std::move(old_fn); }
      return status;
    },
    nb::arg("NLS"), nb::arg("CTestFn").none());

  sundials4py::scoped_def(
    m, "SUNNonlinSolSetGetConvRateFn",
    [](SUNNonlinearSolver NLS,
       std::function<SUNNonlinSolGetConvRateStdFn> GetConvRateFn) -> SUNErrCode
    {
      if (nb::object impl = CustomSUNNonlinearSolver::_python_object_for(NLS);
          impl.is_valid())
      {
        nb::object adapter = nb::none();
        if (GetConvRateFn)
        {
          adapter = nb::cpp_function([GetConvRateFn]()
                                     { return GetConvRateFn(nullptr); });
        }
        return nb::cast<int>(impl.attr("set_get_conv_rate_fn")(adapter));
      }
      auto fntable           = sunnonlinearsolver_function_table(NLS);
      nb::object old_fn      = fntable->getconvratefn;
      fntable->getconvratefn = nb::cast(GetConvRateFn);
      SUNErrCode status;
      if (GetConvRateFn)
      {
        status =
          SUNNonlinSolSetGetConvRateFn(NLS,
                                       sunnonlinearsolver_getconvratefn_wrapper,
                                       fntable);
      }
      else { status = SUNNonlinSolSetGetConvRateFn(NLS, nullptr, nullptr); }
      if (status != SUN_SUCCESS) { fntable->getconvratefn = std::move(old_fn); }
      return status;
    },
    nb::arg("NLS"), nb::arg("GetConvRateFn").none());
}

} // namespace sundials4py

extern "C" void SUNNonlinearSolverFunctionTable_Destroy(void* ptr)
{
  sundials4py::shutdown_safe_delete(
    static_cast<SUNNonlinearSolverFunctionTable*>(ptr));
}
