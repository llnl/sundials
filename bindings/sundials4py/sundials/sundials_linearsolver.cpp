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
#include "sundials_linearsolver_generated.hpp"

  m.def(
    "SUNLinSolSetOptions",
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

  m.def("SUNLinSolSolve", SUNLinSolSolve, nb::arg("S"), nb::arg("A").none(),
        nb::arg("x"), nb::arg("b"), nb::arg("tol"));

  m.def(
    "SUNLinSolSetATimes",
    [](SUNLinearSolver LS,
       std::function<std::remove_pointer_t<SUNATimesFn>> ATimesFn) -> SUNErrCode
    {
      if (!LS->python) { LS->python = new SUNLinearSolverFunctionTable; }
      auto fn_table = static_cast<SUNLinearSolverFunctionTable*>(LS->python);
      fn_table->ATimesFn = nb::cast(ATimesFn);
      if (ATimesFn)
      {
        return SUNLinSolSetATimes(LS, LS->python,
                                  sunlinearsolver_atimesfn_wrapper);
      }
      else { return SUNLinSolSetATimes(LS, nullptr, nullptr); }
    },
    nb::arg("LS"), nb::arg("ATimes").none());

  m.def(
    "SUNLinSolSetPreconditioner",
    [](SUNLinearSolver LS,
       std::function<std::remove_pointer_t<SUNPSetupFn>> PSetupFn,
       std::function<std::remove_pointer_t<SUNPSetupFn>> PSolveFn) -> SUNErrCode
    {
      if (!LS->python) { LS->python = new SUNLinearSolverFunctionTable; }
      auto fn_table = static_cast<SUNLinearSolverFunctionTable*>(LS->python);
      fn_table->PSetupFn = nb::cast(PSetupFn);
      fn_table->PSolveFn = nb::cast(PSolveFn);
      if (!PSetupFn && PSolveFn)
      {
        return SUNLinSolSetPreconditioner(LS, LS->python, nullptr,
                                          sunlinearsolver_psolvefn_wrapper);
      }
      else if (PSetupFn && PSolveFn)
      {
        return SUNLinSolSetPreconditioner(LS, LS->python,
                                          sunlinearsolver_psetupfn_wrapper,
                                          sunlinearsolver_psolvefn_wrapper);
      }
      else { return SUNLinSolSetPreconditioner(LS, nullptr, nullptr, nullptr); }
    },
    nb::arg("LS"), nb::arg("PSetupFn").none(), nb::arg("PSolveFn"));
}

} // namespace sundials4py

extern "C" void SUNLinearSolverFunctionTable_Destroy(void* ptr)
{
  delete static_cast<SUNLinearSolverFunctionTable*>(ptr);
}
