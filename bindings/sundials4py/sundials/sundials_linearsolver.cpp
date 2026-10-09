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
 * SUNDIALS SUNLinearSolver class. It contains hand-written code for
 * functions that require special treatment, and includes the generated
 * code produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include "sundials/sundials_linearsolver.h"
#include "sundials/sundials_iterative.h"
#include "sundials4py.hpp"

#include <sundials/sundials_linearsolver.hpp>
#include <sundials/sundials_nvector.hpp>

namespace nb = nanobind;
using namespace sundials::experimental;

#include "sundials_linearsolver_usersupplied.hpp"

namespace sundials4py {

void bind_sunlinearsolver(nb::module_& m)
{
  // CustomSUNLinearSolver exposes the Python subclassing surface. The C handle
  // is created lazily by the custom nanobind caster when a SUNLinearSolver is
  // required by a generated wrapper or hand-written binding.
  nb::class_<CustomSUNLinearSolver>(m, "CustomSUNLinearSolver",
                                    nb::dynamic_attr(),
                                    nb::is_weak_referenceable())
    .def(nb::init<std::shared_ptr<std::remove_pointer_t<SUNContext>>, int>(),
         nb::arg("sunctx"), nb::arg("solver_type"))
    .def("_is_materialized", &CustomSUNLinearSolver::_is_materialized)
    .def("validate",
         [](CustomSUNLinearSolver& self) { self.validate(nb::find(&self)); })
    .def_prop_ro("sunctx", &CustomSUNLinearSolver::sunctx,
                 nb::sig("def sunctx(self) -> object"),
                 "The SUNDIALS context owned by this object.")
    .def("set_atimes", [](CustomSUNLinearSolver&, nb::object)
         { return CustomSUNLinearSolver::base_method_int("set_atimes"); })
    .def("set_preconditioner", [](CustomSUNLinearSolver&, nb::object, nb::object)
         { return CustomSUNLinearSolver::base_method_int("set_preconditioner"); })
    .def("set_scaling_vectors",
         [](CustomSUNLinearSolver&, nb::object, nb::object) {
           return CustomSUNLinearSolver::base_method_int("set_scaling_vectors");
         })
    .def("set_zero_guess", [](CustomSUNLinearSolver&, sunbooleantype)
         { return CustomSUNLinearSolver::base_method_int("set_zero_guess"); })
    .def("set_options", [](CustomSUNLinearSolver&, const std::string&,
                           const std::string&, const std::vector<std::string>&)
         { return CustomSUNLinearSolver::base_method_int("set_options"); })
    .def("initialize", [](CustomSUNLinearSolver&)
         { return CustomSUNLinearSolver::base_method_int("initialize"); })
    .def("setup", [](CustomSUNLinearSolver&, nb::object)
         { return CustomSUNLinearSolver::base_method_int("setup"); })
    .def("solve", [](CustomSUNLinearSolver&, nb::object, nb::object, nb::object,
                     sunrealtype)
         { return CustomSUNLinearSolver::base_method_int("solve"); })
    .def("num_iters", [](CustomSUNLinearSolver&)
         { return CustomSUNLinearSolver::base_method_int("num_iters"); })
    .def("res_norm", [](CustomSUNLinearSolver&)
         { return CustomSUNLinearSolver::base_method_real("res_norm"); })
    .def("resid", [](CustomSUNLinearSolver&)
         { return CustomSUNLinearSolver::base_method_int("resid"); });

#include "sundials_linearsolver_generated.hpp"

  sundials4py::scoped_def(
    m, "SUNLinSolSetOptions",
    [](SUNLinearSolver self, const std::string& id, const std::string& file_name,
       int argc, const std::vector<std::string>& args)
    {
      std::vector<char*> argv;
      argv.reserve(args.size());

      for (const auto& arg : args)
      {
        // We need a non-const char*, so we use data() and an explicit cast.
        // This is safe as long as the underlying std::string is not modified.
        argv.push_back(const_cast<char*>(arg.data()));
      }

      return SUNLinSolSetOptions(self, id.empty() ? nullptr : id.c_str(),
                                 file_name.empty() ? nullptr : file_name.c_str(),
                                 argc, argv.data());
    },
    nb::arg("self"), nb::arg("id"), nb::arg("file_name"), nb::arg("argc"),
    nb::arg("args"));

  sundials4py::scoped_def(m, "SUNLinSolSolve", SUNLinSolSolve, nb::arg("S"),
                          nb::arg("A").none(), nb::arg("x"), nb::arg("b"),
                          nb::arg("tol"));

  sundials4py::scoped_def(
    m, "SUNLinSolSetATimes",
    [](SUNLinearSolver LS,
       std::function<std::remove_pointer_t<SUNATimesFn>> ATimesFn) -> SUNErrCode
    {
      auto fn_table       = sunlinearsolver_function_table(LS);
      nb::object previous = fn_table->ATimesFn;
      fn_table->ATimesFn  = nb::cast(ATimesFn);
      SUNErrCode status;
      if (ATimesFn)
      {
        status = SUNLinSolSetATimes(LS, fn_table,
                                    sunlinearsolver_atimesfn_wrapper);
      }
      else { status = SUNLinSolSetATimes(LS, nullptr, nullptr); }
      if (status != SUN_SUCCESS) { fn_table->ATimesFn = std::move(previous); }
      return status;
    },
    nb::arg("LS"), nb::arg("ATimes").none());

  sundials4py::scoped_def(
    m, "SUNLinSolSetPreconditioner",
    [](SUNLinearSolver LS,
       std::function<std::remove_pointer_t<SUNPSetupFn>> PSetupFn,
       std::function<std::remove_pointer_t<SUNPSolveFn>> PSolveFn) -> SUNErrCode
    {
      auto fn_table        = sunlinearsolver_function_table(LS);
      nb::object old_setup = fn_table->PSetupFn;
      nb::object old_solve = fn_table->PSolveFn;
      fn_table->PSetupFn   = nb::cast(PSetupFn);
      fn_table->PSolveFn   = nb::cast(PSolveFn);
      SUNErrCode status;
      if (PSetupFn && !PSolveFn)
      {
        status = SUNLinSolSetPreconditioner(LS, fn_table,
                                            sunlinearsolver_psetupfn_wrapper,
                                            nullptr);
      }
      else if (!PSetupFn && PSolveFn)
      {
        status = SUNLinSolSetPreconditioner(LS, fn_table, nullptr,
                                            sunlinearsolver_psolvefn_wrapper);
      }
      else if (PSetupFn && PSolveFn)
      {
        status = SUNLinSolSetPreconditioner(LS, fn_table,
                                            sunlinearsolver_psetupfn_wrapper,
                                            sunlinearsolver_psolvefn_wrapper);
      }
      else
      {
        status = SUNLinSolSetPreconditioner(LS, nullptr, nullptr, nullptr);
      }
      if (status != SUN_SUCCESS)
      {
        fn_table->PSetupFn = std::move(old_setup);
        fn_table->PSolveFn = std::move(old_solve);
      }
      return status;
    },
    nb::arg("LS"), nb::arg("PSetupFn").none(), nb::arg("PSolveFn").none());
}

} // namespace sundials4py

extern "C" void SUNLinearSolverFunctionTable_Destroy(void* ptr)
{
  sundials4py::shutdown_safe_delete(
    static_cast<SUNLinearSolverFunctionTable*>(ptr));
}
