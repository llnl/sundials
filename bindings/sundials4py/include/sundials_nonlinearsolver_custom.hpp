/*------------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *------------------------------------------------------------------------------
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
 *------------------------------------------------------------------------------
 * CustomSUNNonlinearSolver: the base class Python code subclasses to implement a
 * SUNNonlinearSolver, including scoped handling of package memory pointers.
 *----------------------------------------------------------------------------*/

#ifndef SUNDIALS4PY_NONLINEARSOLVER_CUSTOM_HPP
#define SUNDIALS4PY_NONLINEARSOLVER_CUSTOM_HPP

#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <sundials/sundials_classview.hpp>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_nonlinearsolver.h>
#include <sundials/sundials_nonlinearsolver.hpp>

#include "sundials4py_core_types.hpp"
#include "sundials4py_custom_object.hpp"

namespace sundials4py {

/*
 * Python-owned SUNNonlinearSolver implementation.
 *
 * This follows the same lazy-handle design as the matrix and linear solver
 * wrappers, but nonlinear solver callbacks pass N_Vector arguments directly
 * through nanobind's existing native vector casters.
 *
 * The distinguishing complication is the opaque `mem` pointer. SUNDIALS hands it
 * to setup() and solve(), and the system/linear-setup/linear-solve callbacks
 * need it, but it is not available when those callbacks are installed. The
 * pointer is therefore recorded for the duration of each package-driven
 * setup()/solve() call. Handwritten Python entry points call a Python custom
 * solver directly and do not route their callback table through this scope.
 */
class CustomSUNNonlinearSolver
  : public CustomObjectBase<CustomSUNNonlinearSolver,
                            std::remove_pointer_t<SUNNonlinearSolver>>
{
public:
  CustomSUNNonlinearSolver(std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx,
                           int solver_type)
    : CustomObjectBase(std::move(sunctx), "CustomSUNNonlinearSolver"),
      solver_type_(checked_type(solver_type))
  {}

  void validate(nb::handle self) const { validate_required_methods(self); }

  static int base_method_int(const char* name)
  {
    PyErr_SetString(PyExc_NotImplementedError, name);
    nb::raise_python_error();
    return SUN_ERR_EXT_FAIL;
  }

private:
  static SUNNonlinearSolver_Type checked_type(int solver_type)
  {
    switch (solver_type)
    {
    case SUNNONLINEARSOLVER_ROOTFIND:
    case SUNNONLINEARSOLVER_FIXEDPOINT:
    case SUNNONLINEARSOLVER_HYBRID:
      return static_cast<SUNNonlinearSolver_Type>(solver_type);
    default: throw nb::value_error("invalid SUNNonlinearSolver_Type");
    }
  }

  struct Content : CustomContentBase
  {
    // Keep the solver kind in content so gettype works from the opaque C handle
    // without ever calling into Python.
    SUNNonlinearSolver_Type solver_type{SUNNONLINEARSOLVER_ROOTFIND};

    // The pointer SUNDIALS supplied for the setup()/solve() call currently in
    // progress.
    void* active_mem{nullptr};
    std::thread::id active_mem_owner{};
  };

  using MemScope = ActiveMemScope<Content>;

  static constexpr const char* label = "CustomSUNNonlinearSolver";

  friend class CustomObjectBase<CustomSUNNonlinearSolver,
                                std::remove_pointer_t<SUNNonlinearSolver>>;

  std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>> make_handle(
    nb::handle impl);

  static std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>> make_handle(
    nb::handle impl,
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
    SUNNonlinearSolver_Type solver_type);

  static Content* get_content(SUNNonlinearSolver NLS)
  {
    if (!NLS || !NLS->ops || NLS->ops->free != custom_nls_free || !NLS->content)
    {
      return nullptr;
    }
    return static_cast<Content*>(NLS->content);
  }

public:
  static nb::object _python_object_for(SUNNonlinearSolver NLS) noexcept
  {
    try
    {
      Content* content = get_content(NLS);
      if (!content) { return nb::object(); }
      nb::object impl = content->weak_impl();
      return impl.is_none() ? nb::object() : impl;
    }
    catch (...)
    {
      return nb::object();
    }
  }

private:
  static nb::object get_impl(SUNNonlinearSolver NLS)
  { return custom_content_impl(get_content(NLS), label); }

  static bool method_overridden(nb::handle impl, const char* name)
  { return custom_method_overridden<CustomSUNNonlinearSolver>(impl, name); }

  static void validate_required_methods(nb::handle impl)
  {
    if (!method_overridden(impl, "solve"))
    {
      throw nb::type_error(
        "CustomSUNNonlinearSolver subclass must override solve()");
    }
  }

  static SUNNonlinearSolver_Type custom_nls_gettype(SUNNonlinearSolver NLS);

  static SUNErrCode custom_nls_initialize(SUNNonlinearSolver NLS);

  static int custom_nls_setup(SUNNonlinearSolver NLS, N_Vector y, void* mem);

  static int custom_nls_solve(SUNNonlinearSolver NLS, N_Vector y0, N_Vector y,
                              N_Vector w, sunrealtype tol,
                              sunbooleantype call_lsetup, void* mem);

  /*
   * Wrap a native system function as a plain Python callable f(y, F) -> status.
   *
   * The `mem` argument the native function needs is not known here, so it is
   * fetched from the active package scope at call time.
   */
  static nb::object make_sys_fn(
    Content* content, SUNNonlinSolSysFn SysFn, const char* what,
    std::shared_ptr<NativeCallbackState<SUNNonlinSolSysFn>>& state)
  {
    state = NativeCallbackRegistry::prepare(SysFn, nullptr);
    if (!state) { return nb::none(); }

    return nb::cpp_function(
      [content, state, what](N_Vector y, N_Vector F) -> int
      {
        require_valid_callback(state.get(), what);
        void* mem = require_active_mem(content, what);
        return state->fn(y, F, mem);
      },
      nb::arg("y"), nb::arg("F"));
  }

  static SUNErrCode custom_nls_setsysfn(SUNNonlinearSolver NLS,
                                        SUNNonlinSolSysFn SysFn);

  static SUNErrCode custom_nls_setsysfns(SUNNonlinearSolver NLS,
                                         SUNNonlinSolSysFn root_fn,
                                         SUNNonlinSolSysFn fixed_point_fn);

  static SUNErrCode custom_nls_setlsetupfn(SUNNonlinearSolver NLS,
                                           SUNNonlinSolLSetupFn SetupFn);

  static SUNErrCode custom_nls_setlsolvefn(SUNNonlinearSolver NLS,
                                           SUNNonlinSolLSolveFn SolveFn);

  static SUNErrCode custom_nls_setctestfn(SUNNonlinearSolver NLS,
                                          SUNNonlinSolConvTestFn CTestFn,
                                          void* ctest_data);

  static SUNErrCode custom_nls_setnormfn(SUNNonlinearSolver NLS,
                                         SUNNonlinSolNormFn NormFn,
                                         void* norm_fn_data);

  static SUNErrCode custom_nls_setgetupdatenormfn(
    SUNNonlinearSolver NLS, SUNNonlinSolGetUpdateNormFn GetUpdateNormFn,
    void* getupdatenorm_data);

  static SUNErrCode custom_nls_setgetconvratefn(
    SUNNonlinearSolver NLS, SUNNonlinSolGetConvRateFn GetConvRateFn,
    void* getconvrate_data);

  static SUNErrCode custom_nls_setoptions(SUNNonlinearSolver NLS,
                                          const char* NLSid,
                                          const char* file_name, int argc,
                                          char* argv[]);

  static SUNErrCode custom_nls_setmaxiters(SUNNonlinearSolver NLS, int maxiters);

  template<typename Value>
  static SUNErrCode call_tuple_getter(SUNNonlinearSolver NLS, const char* name,
                                      Value* out, const char* operation)
  {
    // Getter wrappers return (status, value) in Python because the C API uses
    // an output pointer plus a status code.
    *out = Value{};
    try
    {
      nb::gil_scoped_acquire gil;
      if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
      auto result = nb::cast<std::tuple<int, Value>>(get_impl(NLS).attr(name)());
      *out = std::get<1>(result);
      return static_cast<SUNErrCode>(std::get<0>(result));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr, operation,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_nls_getnumiters(SUNNonlinearSolver NLS,
                                           long int* niters);

  static SUNErrCode custom_nls_getcuriter(SUNNonlinearSolver NLS, int* iter);

  static SUNErrCode custom_nls_getnumconvfails(SUNNonlinearSolver NLS,
                                               long int* nconvfails);

  static SUNErrCode custom_nls_free(SUNNonlinearSolver NLS);

  SUNNonlinearSolver_Type solver_type_;
};

} // namespace sundials4py

#endif // SUNDIALS4PY_NONLINEARSOLVER_CUSTOM_HPP
