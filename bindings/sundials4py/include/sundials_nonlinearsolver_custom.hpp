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
 * SUNNonlinearSolver, including explicit handling of active-memory modes.
 *----------------------------------------------------------------------------*/

#ifndef _SUNDIALS4PY_NONLINEARSOLVER_CUSTOM_HPP
#define _SUNDIALS4PY_NONLINEARSOLVER_CUSTOM_HPP

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
 * pointer is therefore recorded in an explicit scope for the duration of each
 * setup()/solve() call, tagged with which provider supplied it, and each adapter
 * refuses to run unless the provider matches what it was installed for.
 */
class CustomSUNNonlinearSolver
{
public:
  enum class HandleState
  {
    unmaterialized,
    materializing,
    materialized
  };

  CustomSUNNonlinearSolver(std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx,
                           SUNNonlinearSolver_Type solver_type)
    : sunctx_owner_(std::move(sunctx)), solver_type_(solver_type)
  {}

  std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>> _get_sundials_handle(
    nb::handle self)
  {
    if (!sunctx_owner_)
    {
      throw nb::type_error(
        "CustomSUNNonlinearSolver base constructor was not initialized");
    }
    if (state_ == HandleState::materializing)
    {
      throw std::runtime_error(
        "reentrant CustomSUNNonlinearSolver native handle materialization");
    }
    if (state_ == HandleState::materialized) { return solver_; }

    // Materialize once on demand and cache the native shell for future
    // generated-wrapper calls.
    state_ = HandleState::materializing;
    try
    {
      solver_ = make_handle(self, sunctx_owner_, solver_type_);
      materialization_count_++;
      state_ = HandleState::materialized;
    }
    catch (...)
    {
      solver_.reset();
      state_ = HandleState::unmaterialized;
      throw;
    }
    return solver_;
  }

  int _materialization_count() const { return materialization_count_; }

  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx() const
  {
    return sunctx_owner_;
  }

  static int base_method_int(const char* name)
  {
    PyErr_SetString(PyExc_NotImplementedError, name);
    nb::raise_python_error();
    return SUN_ERR_EXT_FAIL;
  }

private:
  struct Content
  {
    SUNDIALS4PY_CUSTOM_CONTENT_MEMBERS(0x53554e4e4c535031ULL);

    // Keep the solver kind in content so gettype works from the opaque C handle
    // without ever calling into Python.
    SUNNonlinearSolver_Type solver_type{SUNNONLINEARSOLVER_ROOTFIND};

    // The pointer SUNDIALS supplied for the setup()/solve() call
    // currently in progress, and which provider supplied it.
    ActiveMemMode active_mem_mode{ActiveMemMode::none};
    void* active_mem{nullptr};
  };

  using MemScope = ActiveMemScope<Content>;

  static constexpr const char* label = "CustomSUNNonlinearSolver";

  static std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>> make_handle(
    nb::handle impl,
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
    SUNNonlinearSolver_Type solver_type)
  {
    validate_required_methods(impl);

    // Resolve every optional operation before allocating the native shell, so
    // an exception from Python's attribute machinery cannot leak that shell.
    const bool has_set_sys_fn       = method_overridden(impl, "set_sys_fn");
    const bool has_set_sys_fns      = method_overridden(impl, "set_sys_fns");
    const bool has_set_lsetup_fn    = method_overridden(impl, "set_lsetup_fn");
    const bool has_set_lsolve_fn    = method_overridden(impl, "set_lsolve_fn");
    const bool has_set_conv_test_fn = method_overridden(impl,
                                                        "set_conv_test_fn");
    const bool has_set_norm_fn      = method_overridden(impl, "set_norm_fn");
    const bool has_set_get_update_norm_fn =
      method_overridden(impl, "set_get_update_norm_fn");
    const bool has_set_get_conv_rate_fn =
      method_overridden(impl, "set_get_conv_rate_fn");
    const bool has_set_options   = method_overridden(impl, "set_options");
    const bool has_initialize    = method_overridden(impl, "initialize");
    const bool has_setup         = method_overridden(impl, "setup");
    const bool has_set_max_iters = method_overridden(impl, "set_max_iters");
    const bool has_get_num_iters = method_overridden(impl, "get_num_iters");
    const bool has_get_cur_iter  = method_overridden(impl, "get_cur_iter");
    const bool has_get_num_conv_fails = method_overridden(impl,
                                                          "get_num_conv_fails");

    // Install solve unconditionally after validation and attach optional
    // operations only when the Python subclass provides concrete overrides.
    std::unique_ptr<Content> content(new Content);
    content->weak_impl    = nb::weakref(impl);
    content->sunctx_owner = std::move(sunctx_owner);
    content->solver_type  = solver_type;

    SUNNonlinearSolver NLS = SUNNonlinSolNewEmpty(content->sunctx_owner.get());
    if (!NLS) { throw error_returned("SUNNonlinSolNewEmpty failed"); }

    NLS->content      = content.release();
    NLS->ops->gettype = custom_nls_gettype;
    NLS->ops->solve   = custom_nls_solve;
    NLS->ops->free    = custom_nls_free;

    if (has_set_sys_fn) { NLS->ops->setsysfn = custom_nls_setsysfn; }
    if (has_set_sys_fns) { NLS->ops->setsysfns = custom_nls_setsysfns; }
    if (has_set_lsetup_fn) { NLS->ops->setlsetupfn = custom_nls_setlsetupfn; }
    if (has_set_lsolve_fn) { NLS->ops->setlsolvefn = custom_nls_setlsolvefn; }
    if (has_set_conv_test_fn) { NLS->ops->setctestfn = custom_nls_setctestfn; }
    if (has_set_norm_fn) { NLS->ops->setnormfn = custom_nls_setnormfn; }
    if (has_set_get_update_norm_fn)
    {
      NLS->ops->setgetupdatenormfn = custom_nls_setgetupdatenormfn;
    }
    if (has_set_get_conv_rate_fn)
    {
      NLS->ops->setgetconvratefn = custom_nls_setgetconvratefn;
    }
    if (has_set_options) { NLS->ops->setoptions = custom_nls_setoptions; }
    if (has_initialize) { NLS->ops->initialize = custom_nls_initialize; }
    if (has_setup) { NLS->ops->setup = custom_nls_setup; }
    if (has_set_max_iters) { NLS->ops->setmaxiters = custom_nls_setmaxiters; }
    if (has_get_num_iters) { NLS->ops->getnumiters = custom_nls_getnumiters; }
    if (has_get_cur_iter) { NLS->ops->getcuriter = custom_nls_getcuriter; }
    if (has_get_num_conv_fails)
    {
      NLS->ops->getnumconvfails = custom_nls_getnumconvfails;
    }

    return sundials::experimental::our_make_shared<
      std::remove_pointer_t<SUNNonlinearSolver>,
      sundials::experimental::SUNNonlinearSolverDeleter>(NLS);
  }

  static Content* get_content(SUNNonlinearSolver NLS)
  {
    return custom_content_cast<Content>(NLS);
  }

  static nb::object get_impl(SUNNonlinearSolver NLS)
  {
    return custom_content_impl(get_content(NLS), label);
  }

  static bool method_overridden(nb::handle impl, const char* name)
  {
    return custom_method_overridden<CustomSUNNonlinearSolver>(impl, name);
  }

  static void validate_required_methods(nb::handle impl)
  {
    if (!method_overridden(impl, "solve"))
    {
      throw nb::type_error(
        "CustomSUNNonlinearSolver subclass must override solve()");
    }
  }

  /*
   * Determine which provider is driving the vtable right now. A hand-written
   * sundials4py wrapper announces itself with DirectBindingScope; anything else
   * reaching setup()/solve() is a SUNDIALS package.
   */
  static ActiveMemMode entry_mode(void* mem)
  {
    if (!mem) { return ActiveMemMode::none; }
    return DirectBindingScope::active() ? ActiveMemMode::direct_binding
                                        : ActiveMemMode::integrator;
  }

  static SUNNonlinearSolver_Type custom_nls_gettype(SUNNonlinearSolver NLS)
  {
    Content* content = get_content(NLS);
    return content ? content->solver_type : SUNNONLINEARSOLVER_ROOTFIND;
  }

  static SUNErrCode custom_nls_initialize(SUNNonlinearSolver NLS)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      return static_cast<SUNErrCode>(
        nb::cast<int>(get_impl(NLS).attr("initialize")()));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static int custom_nls_setup(SUNNonlinearSolver NLS, N_Vector y, void* mem)
  {
    Content* content = get_content(NLS);
    try
    {
      nb::gil_scoped_acquire gil;
      // The scope restores the previous active memory on every exit path, so a
      // package that calls setup() from inside solve() nests correctly and a
      // Python exception cannot leave a stale pointer behind.
      MemScope scope(content, entry_mode(mem), mem);
      return nb::cast<int>(
        get_impl(NLS).attr("setup")(nb::cast(y, nb::rv_policy::reference)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static int custom_nls_solve(SUNNonlinearSolver NLS, N_Vector y0, N_Vector y,
                              N_Vector w, sunrealtype tol,
                              sunbooleantype call_lsetup, void* mem)
  {
    Content* content = get_content(NLS);
    try
    {
      nb::gil_scoped_acquire gil;
      MemScope scope(content, entry_mode(mem), mem);
      return nb::cast<int>(get_impl(NLS).attr(
        "solve")(nb::cast(y0, nb::rv_policy::reference),
                 nb::cast(y, nb::rv_policy::reference),
                 nb::cast(w, nb::rv_policy::reference), tol, call_lsetup));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  /*
   * Wrap a native system function as a plain Python callable f(y, F) -> status.
   *
   * The `mem` argument the native function needs is not known here, so it is
   * fetched from the active scope at call time and checked against the provider
   * that installed this callback.
   */
  static nb::object make_sys_fn(SUNNonlinearSolver NLS, Content* content,
                                const char* slot, SUNNonlinSolSysFn SysFn,
                                const char* what)
  {
    auto state = content ? content->callbacks.install(slot, SysFn, nullptr)
                         : nullptr;
    if (!state) { return nb::none(); }

    return nb::cpp_function(
      [content, state, what](N_Vector y, N_Vector F) -> int
      {
        require_valid_callback(state.get(), what);
        void* mem = require_active_mem(content, state->required_mode, what);
        return state->fn(y, F, mem);
      },
      nb::arg("y"), nb::arg("F"));
  }

  static SUNErrCode custom_nls_setsysfn(SUNNonlinearSolver NLS,
                                        SUNNonlinSolSysFn SysFn)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      Content* content = get_content(NLS);
      nb::object impl  = get_impl(NLS);
      return static_cast<SUNErrCode>(nb::cast<int>(impl.attr("set_sys_fn")(
        make_sys_fn(NLS, content, "sysfn", SysFn, "nonlinear system"))));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_nls_setsysfns(SUNNonlinearSolver NLS,
                                         SUNNonlinSolSysFn root_fn,
                                         SUNNonlinSolSysFn fixed_point_fn)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      Content* content = get_content(NLS);
      nb::object impl  = get_impl(NLS);
      // Hybrid solvers receive both forms and choose per iteration.
      nb::object root        = make_sys_fn(NLS, content, "rootsysfn", root_fn,
                                           "root-find nonlinear system");
      nb::object fixed_point = make_sys_fn(NLS, content, "fixedpointsysfn",
                                           fixed_point_fn,
                                           "fixed-point nonlinear system");
      return static_cast<SUNErrCode>(
        nb::cast<int>(impl.attr("set_sys_fns")(root, fixed_point)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_nls_setlsetupfn(SUNNonlinearSolver NLS,
                                           SUNNonlinSolLSetupFn SetupFn)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      Content* content = get_content(NLS);
      nb::object impl  = get_impl(NLS);

      auto state = content
                     ? content->callbacks.install("lsetupfn", SetupFn, nullptr)
                     : nullptr;

      nb::object setup = nb::none();
      if (state)
      {
        // The C signature returns the updated Jacobian status through a pointer,
        // which becomes the second element of a Python tuple.
        setup = nb::cpp_function(
          [content, state](sunbooleantype jbad) -> std::tuple<int, sunbooleantype>
          {
            require_valid_callback(state.get(), "linear solver setup");
            void* mem = require_active_mem(content, state->required_mode,
                                           "linear solver setup");
            sunbooleantype jcur = SUNFALSE;
            int status          = state->fn(jbad, &jcur, mem);
            return std::make_tuple(status, jcur);
          },
          nb::arg("jbad"));
      }
      return static_cast<SUNErrCode>(
        nb::cast<int>(impl.attr("set_lsetup_fn")(setup)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_nls_setlsolvefn(SUNNonlinearSolver NLS,
                                           SUNNonlinSolLSolveFn SolveFn)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      Content* content = get_content(NLS);
      nb::object impl  = get_impl(NLS);

      auto state = content
                     ? content->callbacks.install("lsolvefn", SolveFn, nullptr)
                     : nullptr;

      nb::object solve = nb::none();
      if (state)
      {
        solve = nb::cpp_function(
          [content, state](N_Vector b) -> int
          {
            require_valid_callback(state.get(), "linear solver solve");
            void* mem = require_active_mem(content, state->required_mode,
                                           "linear solver solve");
            return state->fn(b, mem);
          },
          nb::arg("b"));
      }
      return static_cast<SUNErrCode>(
        nb::cast<int>(impl.attr("set_lsolve_fn")(solve)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_nls_setctestfn(SUNNonlinearSolver NLS,
                                          SUNNonlinSolConvTestFn CTestFn,
                                          void* ctest_data)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      Content* content = get_content(NLS);
      nb::object impl  = get_impl(NLS);

      // This callback carries its own data pointer, so it does not consult the
      // active-memory scope; it only needs to be revocable.
      auto state = content
                     ? content->callbacks.install("ctestfn", CTestFn, ctest_data)
                     : nullptr;

      nb::object ctest = nb::none();
      if (state)
      {
        // `delta` rather than `del`, which is a Python keyword and so could not
        // be passed by name.
        ctest = nb::cpp_function(
          [NLS, state](N_Vector y, N_Vector delta, sunrealtype tol,
                       N_Vector ewt) -> int
          {
            require_valid_callback(state.get(), "convergence test");
            return state->fn(NLS, y, delta, tol, ewt, state->data);
          },
          nb::arg("y"), nb::arg("delta"), nb::arg("tol"), nb::arg("ewt"));
      }
      return static_cast<SUNErrCode>(
        nb::cast<int>(impl.attr("set_conv_test_fn")(ctest)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_nls_setnormfn(SUNNonlinearSolver NLS,
                                         SUNNonlinSolNormFn NormFn,
                                         void* norm_fn_data)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      Content* content = get_content(NLS);
      nb::object impl  = get_impl(NLS);

      auto state = content
                     ? content->callbacks.install("normfn", NormFn, norm_fn_data)
                     : nullptr;

      nb::object norm = nb::none();
      if (state)
      {
        norm = nb::cpp_function(
          [state](N_Vector delta, N_Vector w) -> std::tuple<SUNErrCode, sunrealtype>
          {
            require_valid_callback(state.get(), "convergence-test norm");
            sunrealtype delnrm = SUN_RCONST(0.0);
            SUNErrCode status  = state->fn(delta, w, &delnrm, state->data);
            return std::make_tuple(status, delnrm);
          },
          nb::arg("delta"), nb::arg("w"));
      }
      return static_cast<SUNErrCode>(
        nb::cast<int>(impl.attr("set_norm_fn")(norm)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_nls_setgetupdatenormfn(
    SUNNonlinearSolver NLS, SUNNonlinSolGetUpdateNormFn GetUpdateNormFn,
    void* getupdatenorm_data)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      Content* content = get_content(NLS);
      nb::object impl  = get_impl(NLS);

      auto state = content ? content->callbacks.install("getupdatenormfn",
                                                        GetUpdateNormFn,
                                                        getupdatenorm_data)
                           : nullptr;

      nb::object get_update_norm = nb::none();
      if (state)
      {
        get_update_norm = nb::cpp_function(
          [state]() -> std::tuple<SUNErrCode, sunrealtype>
          {
            require_valid_callback(state.get(), "update-norm getter");
            sunrealtype delnrm = SUN_RCONST(0.0);
            SUNErrCode status  = state->fn(&delnrm, state->data);
            return std::make_tuple(status, delnrm);
          });
      }
      return static_cast<SUNErrCode>(
        nb::cast<int>(impl.attr("set_get_update_norm_fn")(get_update_norm)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_nls_setgetconvratefn(
    SUNNonlinearSolver NLS, SUNNonlinSolGetConvRateFn GetConvRateFn,
    void* getconvrate_data)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      Content* content = get_content(NLS);
      nb::object impl  = get_impl(NLS);

      auto state = content
                     ? content->callbacks.install("getconvratefn", GetConvRateFn,
                                                  getconvrate_data)
                     : nullptr;

      nb::object get_conv_rate = nb::none();
      if (state)
      {
        get_conv_rate = nb::cpp_function(
          [state]() -> std::tuple<SUNErrCode, sunrealtype>
          {
            require_valid_callback(state.get(), "convergence-rate getter");
            sunrealtype crate = SUN_RCONST(0.0);
            SUNErrCode status = state->fn(&crate, state->data);
            return std::make_tuple(status, crate);
          });
      }
      return static_cast<SUNErrCode>(
        nb::cast<int>(impl.attr("set_get_conv_rate_fn")(get_conv_rate)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_nls_setoptions(SUNNonlinearSolver NLS,
                                          const char* NLSid,
                                          const char* file_name, int argc,
                                          char* argv[])
  {
    try
    {
      nb::gil_scoped_acquire gil;
      std::vector<std::string> args;
      args.reserve(static_cast<size_t>(argc));
      for (int i = 0; i < argc; i++) { args.emplace_back(argv[i]); }
      return static_cast<SUNErrCode>(nb::cast<int>(get_impl(NLS).attr(
        "set_options")(NLSid ? NLSid : "", file_name ? file_name : "", args)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_nls_setmaxiters(SUNNonlinearSolver NLS, int maxiters)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      return static_cast<SUNErrCode>(
        nb::cast<int>(get_impl(NLS).attr("set_max_iters")(maxiters)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  template<typename Value>
  static SUNErrCode call_tuple_getter(SUNNonlinearSolver NLS, const char* name,
                                      Value* out)
  {
    // Getter wrappers return (status, value) in Python because the C API uses
    // an output pointer plus a status code.
    *out = Value{};
    try
    {
      nb::gil_scoped_acquire gil;
      auto result = nb::cast<std::tuple<int, Value>>(get_impl(NLS).attr(name)());
      *out = std::get<1>(result);
      return static_cast<SUNErrCode>(std::get<0>(result));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_nls_getnumiters(SUNNonlinearSolver NLS,
                                           long int* niters)
  {
    return call_tuple_getter(NLS, "get_num_iters", niters);
  }

  static SUNErrCode custom_nls_getcuriter(SUNNonlinearSolver NLS, int* iter)
  {
    return call_tuple_getter(NLS, "get_cur_iter", iter);
  }

  static SUNErrCode custom_nls_getnumconvfails(SUNNonlinearSolver NLS,
                                               long int* nconvfails)
  {
    return call_tuple_getter(NLS, "get_num_conv_fails", nconvfails);
  }

  static SUNErrCode custom_nls_free(SUNNonlinearSolver NLS)
  {
    if (!NLS) { return SUN_SUCCESS; }
    Content* content = get_content(NLS);
    custom_content_destroy(content);
    NLS->content = nullptr;
    SUNNonlinSolFreeEmpty(NLS);
    return SUN_SUCCESS;
  }

  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner_;
  SUNNonlinearSolver_Type solver_type_;
  std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>> solver_;
  HandleState state_{HandleState::unmaterialized};
  int materialization_count_{0};
};

} // namespace sundials4py

#endif // _SUNDIALS4PY_NONLINEARSOLVER_CUSTOM_HPP
