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
 * This file defines the sundials4py test module, which contains
 * utility functions for testing purposes.
 * -----------------------------------------------------------------*/

#include "sundials4py.hpp"

#include <cstdint>
#include <sundials/sundials_context.hpp>
#include <sundials/sundials_errors.h>
#include <tuple>
#include <utility>

#include <sundials/priv/sundials_errors_impl.h>

#include <thread>

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

namespace {

int test_conv_callback_one(SUNNonlinearSolver, N_Vector, N_Vector, sunrealtype,
                           N_Vector, void*)
{ return SUN_SUCCESS; }

int test_conv_callback_two(SUNNonlinearSolver, N_Vector, N_Vector, sunrealtype,
                           N_Vector, void*)
{ return SUN_ERR_EXT_FAIL; }

} // namespace

void bind_test(nb::module_& m)
{
  sundials4py::scoped_def(
    m, "SUNContext_TestErrHandler",
    [](SUNContext sunctx)
    {
      SUNHandleErrWithMsg(__LINE__, __func__, __FILE__,
                          "create an error to test the error handlers",
                          SUN_ERR_ARG_CORRUPT, sunctx);
    },
    "This function is for testing purposes and should not be called.");

  sundials4py::scoped_def(
    m, "SUNMat_TestDestroyCloneOnNativeThread",
    [](SUNMatrix matrix)
    {
      SUNMatrix clone = SUNMatClone(matrix);
      if (!clone) { throw sundials4py::error_returned("SUNMatClone failed"); }

      nb::object implementation = CustomSUNMatrix::_python_object_for(clone);
      if (!implementation.is_valid())
      {
        SUNMatDestroy(clone);
        throw nb::type_error("matrix is not a custom SUNMatrix");
      }

      nb::weakref implementation_ref(implementation);
      implementation.reset();
      {
        nb::gil_scoped_release release;
        std::thread([clone] { SUNMatDestroy(clone); }).join();
      }
      return implementation_ref;
    },
    nb::arg("matrix"),
    "Clone a custom matrix and destroy its native handle on a native thread.");

  sundials4py::scoped_def(
    m, "SUNMat_TestHandleAddress",
    [](SUNMatrix matrix) { return reinterpret_cast<std::uintptr_t>(matrix); },
    nb::arg("matrix"), "Return a SUNMatrix address for clone identity tests.");

  sundials4py::scoped_def(
    m, "SUNMat_TestStrongClone",
    [](SUNMatrix matrix) -> std::tuple<nb::object, nb::capsule>
    {
      SUNMatrix clone = SUNMatClone(matrix);
      if (!clone) { throw sundials4py::error_returned("SUNMatClone failed"); }
      nb::object implementation = CustomSUNMatrix::_python_object_for(clone);
      if (!implementation.is_valid())
      {
        SUNMatDestroy(clone);
        throw nb::type_error("matrix is not a custom SUNMatrix");
      }
      nb::capsule keepalive(clone, [](void* ptr) noexcept
                            { SUNMatDestroy(reinterpret_cast<SUNMatrix>(ptr)); });
      return std::make_tuple(implementation, std::move(keepalive));
    },
    nb::arg("matrix"), "Return a SUNDIALS-created custom matrix clone and keep it alive for tests.");

  sundials4py::scoped_def(
    m, "N_VTestHandleAddress",
    [](N_Vector vector) { return reinterpret_cast<std::uintptr_t>(vector); },
    nb::arg("vector"), "Return an N_Vector address for clone identity tests.");

  sundials4py::scoped_def(
    m, "N_VTestStrongClone",
    [](N_Vector vector) -> std::tuple<nb::object, nb::capsule>
    {
      N_Vector clone = N_VClone(vector);
      if (!clone) { throw sundials4py::error_returned("N_VClone failed"); }
      nb::object implementation = CustomNVector::_python_object_for(clone);
      if (!implementation.is_valid())
      {
        N_VDestroy(clone);
        throw nb::type_error("vector is not a custom N_Vector");
      }
      nb::capsule keepalive(clone, [](void* ptr) noexcept
                            { N_VDestroy(reinterpret_cast<N_Vector>(ptr)); });
      return std::make_tuple(implementation, std::move(keepalive));
    },
    nb::arg("vector"), "Return a SUNDIALS-created custom N_Vector clone and keep it alive for tests.");

  sundials4py::scoped_def(
    m, "SUNNonlinSol_TestSetupFromPackage",
    [](SUNNonlinearSolver solver, N_Vector y)
    { return solver->ops->setup(solver, y, solver); }, nb::arg("solver"),
    nb::arg("y"));

  sundials4py::scoped_def(
    m, "SUNNonlinSol_TestSolveFromPackage",
    [](SUNNonlinearSolver solver, N_Vector y0, N_Vector y, N_Vector w,
       sunrealtype tol, sunbooleantype call_lsetup)
    { return solver->ops->solve(solver, y0, y, w, tol, call_lsetup, solver); },
    nb::arg("solver"), nb::arg("y0"), nb::arg("y"), nb::arg("w"),
    nb::arg("tol"), nb::arg("call_lsetup"));

  sundials4py::scoped_def(
    m, "SUNNonlinSol_TestSetConvTestFnFromPackage",
    [](SUNNonlinearSolver solver, int which)
    {
      SUNNonlinSolConvTestFn fn = which == 1 ? test_conv_callback_one
                                             : test_conv_callback_two;
      return solver->ops->setctestfn(solver, fn, nullptr);
    },
    nb::arg("solver"), nb::arg("which"));
}

} // namespace sundials4py
