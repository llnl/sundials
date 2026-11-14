/*-----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *----------------------------------------------------------------------------*/

#include "sundials4py.hpp"

namespace sundials4py {

std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>> CustomSUNNonlinearSolver::make_handle(
  nb::handle impl)
{
  return make_handle(impl, sunctx_owner_, solver_type_);
}

std::shared_ptr<std::remove_pointer_t<SUNNonlinearSolver>> CustomSUNNonlinearSolver::make_handle(
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
  const bool has_set_conv_test_fn = method_overridden(impl, "set_conv_test_fn");
  const bool has_set_norm_fn      = method_overridden(impl, "set_norm_fn");
  const bool has_set_get_update_norm_fn =
    method_overridden(impl, "set_get_update_norm_fn");
  const bool has_set_get_conv_rate_fn =
    method_overridden(impl, "set_get_conv_rate_fn");
  const bool has_set_options        = method_overridden(impl, "set_options");
  const bool has_initialize         = method_overridden(impl, "initialize");
  const bool has_setup              = method_overridden(impl, "setup");
  const bool has_set_max_iters      = method_overridden(impl, "set_max_iters");
  const bool has_get_num_iters      = method_overridden(impl, "get_num_iters");
  const bool has_get_cur_iter       = method_overridden(impl, "get_cur_iter");
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

SUNNonlinearSolver_Type CustomSUNNonlinearSolver::custom_nls_gettype(
  SUNNonlinearSolver NLS)
{
  Content* content = get_content(NLS);
  return content ? content->solver_type : SUNNONLINEARSOLVER_ROOTFIND;
}

SUNErrCode CustomSUNNonlinearSolver::custom_nls_initialize(SUNNonlinearSolver NLS)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    return static_cast<SUNErrCode>(
      nb::cast<int>(get_impl(NLS).attr("initialize")()));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr,
                               "CustomSUNNonlinearSolver.initialize",
                               SUN_ERR_EXT_FAIL)
}

int CustomSUNNonlinearSolver::custom_nls_setup(SUNNonlinearSolver NLS,
                                               N_Vector y, void* mem)
{
  Content* content = get_content(NLS);
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    // The scope restores the previous active memory on every exit path, so a
    // package that calls setup() from inside solve() nests correctly and a
    // Python exception cannot leave a stale pointer behind.
    MemScope scope(content, mem);
    return nb::cast<int>(
      get_impl(NLS).attr("setup")(nb::cast(y, nb::rv_policy::reference)));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr,
                               "CustomSUNNonlinearSolver.setup", SUN_ERR_EXT_FAIL)
}

int CustomSUNNonlinearSolver::custom_nls_solve(SUNNonlinearSolver NLS,
                                               N_Vector y0, N_Vector y,
                                               N_Vector w, sunrealtype tol,
                                               sunbooleantype call_lsetup,
                                               void* mem)
{
  Content* content = get_content(NLS);
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    MemScope scope(content, mem);
    return nb::cast<int>(get_impl(NLS).attr(
      "solve")(nb::cast(y0, nb::rv_policy::reference),
               nb::cast(y, nb::rv_policy::reference),
               nb::cast(w, nb::rv_policy::reference), tol, call_lsetup));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr,
                               "CustomSUNNonlinearSolver.solve", SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNNonlinearSolver::custom_nls_setsysfn(SUNNonlinearSolver NLS,
                                                         SUNNonlinSolSysFn SysFn)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    Content* content = get_content(NLS);
    nb::object impl  = get_impl(NLS);
    std::shared_ptr<NativeCallbackState<SUNNonlinSolSysFn>> state;
    nb::object callback = make_sys_fn(content, SysFn, "nonlinear system", state);
    PreparedCallback prepared(state);
    SUNErrCode status =
      static_cast<SUNErrCode>(nb::cast<int>(impl.attr("set_sys_fn")(callback)));
    if (status == SUN_SUCCESS)
    {
      prepared.commit(content->callbacks, NativeCallbackSlot::sysfn);
    }
    return status;
  }
  SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr,
                               "CustomSUNNonlinearSolver.set_sys_fn",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNNonlinearSolver::custom_nls_setsysfns(
  SUNNonlinearSolver NLS, SUNNonlinSolSysFn root_fn,
  SUNNonlinSolSysFn fixed_point_fn)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    Content* content = get_content(NLS);
    nb::object impl  = get_impl(NLS);
    // Hybrid solvers receive both forms and choose per iteration.
    std::shared_ptr<NativeCallbackState<SUNNonlinSolSysFn>> root_state;
    std::shared_ptr<NativeCallbackState<SUNNonlinSolSysFn>> fixed_state;
    nb::object root        = make_sys_fn(content, root_fn,
                                         "root-find nonlinear system", root_state);
    nb::object fixed_point = make_sys_fn(content, fixed_point_fn,
                                         "fixed-point nonlinear system",
                                         fixed_state);
    PreparedCallback prepared_root(root_state);
    PreparedCallback prepared_fixed(fixed_state);
    SUNErrCode status = static_cast<SUNErrCode>(
      nb::cast<int>(impl.attr("set_sys_fns")(root, fixed_point)));
    if (status == SUN_SUCCESS)
    {
      prepared_root.commit(content->callbacks, NativeCallbackSlot::rootsysfn);
      prepared_fixed.commit(content->callbacks,
                            NativeCallbackSlot::fixedpointsysfn);
    }
    return status;
  }
  SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr,
                               "CustomSUNNonlinearSolver.set_sys_fns",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNNonlinearSolver::custom_nls_setlsetupfn(
  SUNNonlinearSolver NLS, SUNNonlinSolLSetupFn SetupFn)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    Content* content = get_content(NLS);
    nb::object impl  = get_impl(NLS);

    auto state = NativeCallbackRegistry::prepare(SetupFn, nullptr);
    PreparedCallback prepared(state);

    nb::object setup = nb::none();
    if (state)
    {
      // The C signature returns the updated Jacobian status through a pointer,
      // which becomes the second element of a Python tuple.
      setup =
        nb::cpp_function(sundials4py::scoped(
                           [content, state](sunbooleantype jbad)
                             -> std::tuple<int, sunbooleantype>
                           {
                             require_valid_callback(state.get(),
                                                    "linear solver setup");
                             void* mem =
                               require_active_mem(content,
                                                  "linear solver setup");
                             sunbooleantype jcur = SUNFALSE;
                             int status          = state->fn(jbad, &jcur, mem);
                             return std::make_tuple(status, jcur);
                           }),
                         nb::arg("jbad"));
    }
    SUNErrCode status =
      static_cast<SUNErrCode>(nb::cast<int>(impl.attr("set_lsetup_fn")(setup)));
    if (status == SUN_SUCCESS)
    {
      prepared.commit(content->callbacks, NativeCallbackSlot::lsetupfn);
    }
    return status;
  }
  SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr,
                               "CustomSUNNonlinearSolver.set_lsetup_fn",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNNonlinearSolver::custom_nls_setlsolvefn(
  SUNNonlinearSolver NLS, SUNNonlinSolLSolveFn SolveFn)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    Content* content = get_content(NLS);
    nb::object impl  = get_impl(NLS);

    auto state = NativeCallbackRegistry::prepare(SolveFn, nullptr);
    PreparedCallback prepared(state);

    nb::object solve = nb::none();
    if (state)
    {
      solve =
        nb::cpp_function(sundials4py::scoped(
                           [content, state](N_Vector b) -> int
                           {
                             require_valid_callback(state.get(),
                                                    "linear solver solve");
                             void* mem =
                               require_active_mem(content,
                                                  "linear solver solve");
                             return state->fn(b, mem);
                           }),
                         nb::arg("b"));
    }
    SUNErrCode status =
      static_cast<SUNErrCode>(nb::cast<int>(impl.attr("set_lsolve_fn")(solve)));
    if (status == SUN_SUCCESS)
    {
      prepared.commit(content->callbacks, NativeCallbackSlot::lsolvefn);
    }
    return status;
  }
  SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr,
                               "CustomSUNNonlinearSolver.set_lsolve_fn",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNNonlinearSolver::custom_nls_setctestfn(
  SUNNonlinearSolver NLS, SUNNonlinSolConvTestFn CTestFn, void* ctest_data)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    Content* content = get_content(NLS);
    nb::object impl  = get_impl(NLS);

    // This callback carries its own data pointer, so it does not consult the
    // active-memory scope; it only needs to be revocable.
    auto state = NativeCallbackRegistry::prepare(CTestFn, ctest_data);
    PreparedCallback prepared(state);

    nb::object ctest = nb::none();
    if (state)
    {
      // `delta` rather than `del`, which is a Python keyword and so could not
      // be passed by name.
      ctest = nb::cpp_function(sundials4py::scoped(
                                 [NLS, state](N_Vector y, N_Vector delta,
                                              sunrealtype tol, N_Vector ewt) -> int
                                 {
                                   require_valid_callback(state.get(),
                                                          "convergence test");
                                   return state->fn(NLS, y, delta, tol, ewt,
                                                    state->data);
                                 }),
                               nb::arg("y"), nb::arg("delta"), nb::arg("tol"),
                               nb::arg("ewt"));
    }
    SUNErrCode status = static_cast<SUNErrCode>(
      nb::cast<int>(impl.attr("set_conv_test_fn")(ctest)));
    if (status == SUN_SUCCESS)
    {
      prepared.commit(content->callbacks, NativeCallbackSlot::ctestfn);
    }
    return status;
  }
  SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr,
                               "CustomSUNNonlinearSolver.set_conv_test_fn",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNNonlinearSolver::custom_nls_setnormfn(
  SUNNonlinearSolver NLS, SUNNonlinSolNormFn NormFn, void* norm_fn_data)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    Content* content = get_content(NLS);
    nb::object impl  = get_impl(NLS);

    auto state = NativeCallbackRegistry::prepare(NormFn, norm_fn_data);
    PreparedCallback prepared(state);

    nb::object norm = nb::none();
    if (state)
    {
      norm =
        nb::cpp_function(sundials4py::scoped(
                           [state](N_Vector delta, N_Vector w)
                             -> std::tuple<SUNErrCode, sunrealtype>
                           {
                             require_valid_callback(state.get(),
                                                    "convergence-test norm");
                             sunrealtype delnrm = SUN_RCONST(0.0);
                             SUNErrCode status  = state->fn(delta, w, &delnrm,
                                                            state->data);
                             return std::make_tuple(status, delnrm);
                           }),
                         nb::arg("delta"), nb::arg("w"));
    }
    SUNErrCode status =
      static_cast<SUNErrCode>(nb::cast<int>(impl.attr("set_norm_fn")(norm)));
    if (status == SUN_SUCCESS)
    {
      prepared.commit(content->callbacks, NativeCallbackSlot::normfn);
    }
    return status;
  }
  SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr,
                               "CustomSUNNonlinearSolver.set_norm_fn",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNNonlinearSolver::custom_nls_setgetupdatenormfn(
  SUNNonlinearSolver NLS, SUNNonlinSolGetUpdateNormFn GetUpdateNormFn,
  void* getupdatenorm_data)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    Content* content = get_content(NLS);
    nb::object impl  = get_impl(NLS);

    auto state = NativeCallbackRegistry::prepare(GetUpdateNormFn,
                                                 getupdatenorm_data);
    PreparedCallback prepared(state);

    nb::object get_update_norm = nb::none();
    if (state)
    {
      get_update_norm = nb::cpp_function(sundials4py::scoped(
        [state]() -> std::tuple<SUNErrCode, sunrealtype>
        {
          require_valid_callback(state.get(), "update-norm getter");
          sunrealtype delnrm = SUN_RCONST(0.0);
          SUNErrCode status  = state->fn(&delnrm, state->data);
          return std::make_tuple(status, delnrm);
        }));
    }
    SUNErrCode status = static_cast<SUNErrCode>(
      nb::cast<int>(impl.attr("set_get_update_norm_fn")(get_update_norm)));
    if (status == SUN_SUCCESS)
    {
      prepared.commit(content->callbacks, NativeCallbackSlot::getupdatenormfn);
    }
    return status;
  }
  SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr, "CustomSUNNonlinearSolver.set_get_update_norm_fn",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNNonlinearSolver::custom_nls_setgetconvratefn(
  SUNNonlinearSolver NLS, SUNNonlinSolGetConvRateFn GetConvRateFn,
  void* getconvrate_data)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    Content* content = get_content(NLS);
    nb::object impl  = get_impl(NLS);

    auto state = NativeCallbackRegistry::prepare(GetConvRateFn, getconvrate_data);
    PreparedCallback prepared(state);

    nb::object get_conv_rate = nb::none();
    if (state)
    {
      get_conv_rate = nb::cpp_function(sundials4py::scoped(
        [state]() -> std::tuple<SUNErrCode, sunrealtype>
        {
          require_valid_callback(state.get(), "convergence-rate getter");
          sunrealtype crate = SUN_RCONST(0.0);
          SUNErrCode status = state->fn(&crate, state->data);
          return std::make_tuple(status, crate);
        }));
    }
    SUNErrCode status = static_cast<SUNErrCode>(
      nb::cast<int>(impl.attr("set_get_conv_rate_fn")(get_conv_rate)));
    if (status == SUN_SUCCESS)
    {
      prepared.commit(content->callbacks, NativeCallbackSlot::getconvratefn);
    }
    return status;
  }
  SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr,
                               "CustomSUNNonlinearSolver.set_get_conv_rate_fn",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNNonlinearSolver::custom_nls_setoptions(SUNNonlinearSolver NLS,
                                                           const char* NLSid,
                                                           const char* file_name,
                                                           int argc, char* argv[])
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    std::vector<std::string> args;
    args.reserve(static_cast<size_t>(argc));
    for (int i = 0; i < argc; i++) { args.emplace_back(argv[i]); }
    return static_cast<SUNErrCode>(nb::cast<int>(get_impl(NLS).attr(
      "set_options")(NLSid ? NLSid : "", file_name ? file_name : "", args)));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr,
                               "CustomSUNNonlinearSolver.set_options",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNNonlinearSolver::custom_nls_setmaxiters(SUNNonlinearSolver NLS,
                                                            int maxiters)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    return static_cast<SUNErrCode>(
      nb::cast<int>(get_impl(NLS).attr("set_max_iters")(maxiters)));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(NLS ? NLS->sunctx : nullptr,
                               "CustomSUNNonlinearSolver.set_max_iters",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNNonlinearSolver::custom_nls_getnumiters(SUNNonlinearSolver NLS,
                                                            long int* niters)
{
  return call_tuple_getter(NLS, "get_num_iters", niters,
                           "CustomSUNNonlinearSolver.get_num_iters");
}

SUNErrCode CustomSUNNonlinearSolver::custom_nls_getcuriter(SUNNonlinearSolver NLS,
                                                           int* iter)
{
  return call_tuple_getter(NLS, "get_cur_iter", iter,
                           "CustomSUNNonlinearSolver.get_cur_iter");
}

SUNErrCode CustomSUNNonlinearSolver::custom_nls_getnumconvfails(
  SUNNonlinearSolver NLS, long int* nconvfails)
{
  return call_tuple_getter(NLS, "get_num_conv_fails", nconvfails,
                           "CustomSUNNonlinearSolver.get_num_conv_fails");
}

SUNErrCode CustomSUNNonlinearSolver::custom_nls_free(SUNNonlinearSolver NLS)
{
  if (!NLS) { return SUN_SUCCESS; }
  Content* content = get_content(NLS);
  custom_content_destroy(content);
  NLS->content = nullptr;
  SUNNonlinSolFreeEmpty(NLS);
  return SUN_SUCCESS;
}

} // namespace sundials4py
