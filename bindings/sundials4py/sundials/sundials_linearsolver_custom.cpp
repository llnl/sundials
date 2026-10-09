/*-----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *----------------------------------------------------------------------------*/

#include <limits>

#include "sundials4py.hpp"

namespace sundials4py {

std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>> CustomSUNLinearSolver::make_handle(
  nb::handle impl)
{
  return make_handle(impl, sunctx_owner_, solver_type_);
}

std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>> CustomSUNLinearSolver::make_handle(
  nb::handle impl,
  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
  SUNLinearSolver_Type solver_type)
{
  validate_required_methods(impl);

  // Resolve every optional operation before allocating the native shell, so
  // an exception from Python's attribute machinery cannot leak that shell.
  const bool has_set_atimes          = method_overridden(impl, "set_atimes");
  const bool has_set_preconditioner  = method_overridden(impl,
                                                         "set_preconditioner");
  const bool has_set_scaling_vectors = method_overridden(impl,
                                                         "set_scaling_vectors");
  const bool has_set_zero_guess = method_overridden(impl, "set_zero_guess");
  const bool has_set_options    = method_overridden(impl, "set_options");
  const bool has_initialize     = method_overridden(impl, "initialize");
  const bool has_setup          = method_overridden(impl, "setup");
  const bool has_num_iters      = method_overridden(impl, "num_iters");
  const bool has_res_norm       = method_overridden(impl, "res_norm");
  const bool has_resid          = method_overridden(impl, "resid");

  // Build one native shell and populate only the operations Python actually
  // implements. A NULL operation pointer preserves native SUNDIALS semantics.
  std::unique_ptr<Content> content(new Content);
  content->weak_impl    = nb::weakref(impl);
  content->sunctx_owner = std::move(sunctx_owner);
  content->solver_type  = solver_type;

  SUNLinearSolver S = SUNLinSolNewEmpty(content->sunctx_owner.get());
  if (!S) { throw error_returned("SUNLinSolNewEmpty failed"); }

  S->content      = content.release();
  S->ops->gettype = custom_linsol_gettype;
  S->ops->getid   = custom_linsol_getid;
  S->ops->solve   = custom_linsol_solve;
  S->ops->free    = custom_linsol_free;

  if (has_set_atimes) { S->ops->setatimes = custom_linsol_setatimes; }
  if (has_set_preconditioner)
  {
    S->ops->setpreconditioner = custom_linsol_setpreconditioner;
  }
  if (has_set_scaling_vectors)
  {
    S->ops->setscalingvectors = custom_linsol_setscalingvectors;
  }
  if (has_set_zero_guess) { S->ops->setzeroguess = custom_linsol_setzeroguess; }
  if (has_set_options) { S->ops->setoptions = custom_linsol_setoptions; }
  if (has_initialize) { S->ops->initialize = custom_linsol_initialize; }
  if (has_setup) { S->ops->setup = custom_linsol_setup; }
  if (has_num_iters) { S->ops->numiters = custom_linsol_numiters; }
  if (has_res_norm) { S->ops->resnorm = custom_linsol_resnorm; }
  if (has_resid) { S->ops->resid = custom_linsol_resid; }

  return sundials::experimental::our_make_shared<
    std::remove_pointer_t<SUNLinearSolver>,
    sundials::experimental::SUNLinearSolverDeleter>(S);
}

SUNLinearSolver_Type CustomSUNLinearSolver::custom_linsol_gettype(SUNLinearSolver S)
{
  Content* content = get_content(S);
  return content ? content->solver_type : SUNLINEARSOLVER_DIRECT;
}

SUNLinearSolver_ID CustomSUNLinearSolver::custom_linsol_getid(SUNLinearSolver)
{
  return SUNLINEARSOLVER_CUSTOM;
}

SUNErrCode CustomSUNLinearSolver::custom_linsol_setoptions(SUNLinearSolver S,
                                                           const char* LSid,
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
    return static_cast<SUNErrCode>(nb::cast<int>(get_impl(S).attr(
      "set_options")(LSid ? LSid : "", file_name ? file_name : "", args)));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr,
                               "CustomSUNLinearSolver.set_options",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNLinearSolver::custom_linsol_setatimes(SUNLinearSolver S,
                                                          void* A_data,
                                                          SUNATimesFn ATimes)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    Content* content = get_content(S);
    nb::object impl  = get_impl(S);

    // Register revocable state instead of capturing the raw
    // (function, data) pair, and invalidate whatever occupied this slot
    // before. Converting to a plain Python callable keeps SUNDIALS' opaque
    // A_data out of the subclass's interface entirely.
    auto state = NativeCallbackRegistry::prepare(ATimes, A_data);
    PreparedCallback prepared(state);

    nb::object atimes = nb::none();
    if (state)
    {
      atimes = nb::cpp_function(sundials4py::scoped(
                                  [state](N_Vector x, N_Vector y) -> int
                                  {
                                    require_valid_callback(state.get(),
                                                           "ATimes");
                                    return state->fn(state->data, x, y);
                                  }),
                                nb::arg("x"), nb::arg("y"));
    }
    SUNErrCode status =
      static_cast<SUNErrCode>(nb::cast<int>(impl.attr("set_atimes")(atimes)));
    if (status == SUN_SUCCESS)
    {
      prepared.commit(content->callbacks, NativeCallbackSlot::atimes);
    }
    return status;
  }
  SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr,
                               "CustomSUNLinearSolver.set_atimes",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNLinearSolver::custom_linsol_setpreconditioner(
  SUNLinearSolver S, void* P_data, SUNPSetupFn Pset, SUNPSolveFn Psol)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    Content* content = get_content(S);
    nb::object impl  = get_impl(S);

    // SUNDIALS sets both preconditioner halves in one call, so both slots are
    // replaced (and any previous adapters revoked) together.
    auto setup_state = NativeCallbackRegistry::prepare(Pset, P_data);
    auto solve_state = NativeCallbackRegistry::prepare(Psol, P_data);
    PreparedCallback prepared_setup(setup_state);
    PreparedCallback prepared_solve(solve_state);

    nb::object psetup = nb::none();
    nb::object psolve = nb::none();
    if (setup_state)
    {
      psetup = nb::cpp_function(sundials4py::scoped(
        [setup_state]() -> int
        {
          require_valid_callback(setup_state.get(), "preconditioner setup");
          return setup_state->fn(setup_state->data);
        }));
    }
    if (solve_state)
    {
      psolve =
        nb::cpp_function(sundials4py::scoped(
                           [solve_state](N_Vector r, N_Vector z,
                                         sunrealtype tol, int lr) -> int
                           {
                             require_valid_callback(solve_state.get(),
                                                    "preconditioner solve");
                             return solve_state->fn(solve_state->data, r, z,
                                                    tol, lr);
                           }),
                         nb::arg("r"), nb::arg("z"), nb::arg("tol"),
                         nb::arg("lr"));
    }
    SUNErrCode status = static_cast<SUNErrCode>(
      nb::cast<int>(impl.attr("set_preconditioner")(psetup, psolve)));
    if (status == SUN_SUCCESS)
    {
      prepared_setup.commit(content->callbacks, NativeCallbackSlot::psetup);
      prepared_solve.commit(content->callbacks, NativeCallbackSlot::psolve);
    }
    return status;
  }
  SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr,
                               "CustomSUNLinearSolver.set_preconditioner",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNLinearSolver::custom_linsol_setscalingvectors(
  SUNLinearSolver S, N_Vector s1, N_Vector s2)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    Content* content = get_content(S);
    nb::object o1    = s1 ? nb::cast(s1, nb::rv_policy::reference) : nb::none();
    nb::object o2    = s2 ? nb::cast(s2, nb::rv_policy::reference) : nb::none();
    auto status      = static_cast<SUNErrCode>(
      nb::cast<int>(get_impl(S).attr("set_scaling_vectors")(o1, o2)));
    if (status == SUN_SUCCESS && content)
    {
      content->scaling_owners[0] = std::move(o1);
      content->scaling_owners[1] = std::move(o2);
    }
    return status;
  }
  SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr,
                               "CustomSUNLinearSolver.set_scaling_vectors",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNLinearSolver::custom_linsol_setzeroguess(SUNLinearSolver S,
                                                             sunbooleantype onoff)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    return static_cast<SUNErrCode>(
      nb::cast<int>(get_impl(S).attr("set_zero_guess")(onoff)));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr,
                               "CustomSUNLinearSolver.set_zero_guess",
                               SUN_ERR_EXT_FAIL)
}

SUNErrCode CustomSUNLinearSolver::custom_linsol_initialize(SUNLinearSolver S)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    return static_cast<SUNErrCode>(
      nb::cast<int>(get_impl(S).attr("initialize")()));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr,
                               "CustomSUNLinearSolver.initialize",
                               SUN_ERR_EXT_FAIL)
}

int CustomSUNLinearSolver::custom_linsol_setup(SUNLinearSolver S, SUNMatrix A)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    return nb::cast<int>(get_impl(S).attr("setup")(matrix_arg(A)));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr,
                               "CustomSUNLinearSolver.setup", SUN_ERR_EXT_FAIL)
}

int CustomSUNLinearSolver::custom_linsol_solve(SUNLinearSolver S, SUNMatrix A,
                                               N_Vector x, N_Vector b,
                                               sunrealtype tol)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return SUN_ERR_EXT_FAIL; }
    return nb::cast<int>(get_impl(S).attr(
      "solve")(matrix_arg(A), nb::cast(x, nb::rv_policy::reference),
               nb::cast(b, nb::rv_policy::reference), tol));
  }
  SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr,
                               "CustomSUNLinearSolver.solve", SUN_ERR_EXT_FAIL)
}

int CustomSUNLinearSolver::custom_linsol_numiters(SUNLinearSolver S)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return 0; }
    return nb::cast<int>(get_impl(S).attr("num_iters")());
  }
  SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr,
                               "CustomSUNLinearSolver.num_iters", 0)
}

sunrealtype CustomSUNLinearSolver::custom_linsol_resnorm(SUNLinearSolver S)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending())
    {
      return std::numeric_limits<sunrealtype>::infinity();
    }
    return nb::cast<sunrealtype>(get_impl(S).attr("res_norm")());
  }
  SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr,
                               "CustomSUNLinearSolver.res_norm",
                               std::numeric_limits<sunrealtype>::infinity())
}

N_Vector CustomSUNLinearSolver::custom_linsol_resid(SUNLinearSolver S)
{
  try
  {
    nb::gil_scoped_acquire gil;
    if (custom_exception_pending()) { return nullptr; }
    Content* content = get_content(S);
    nb::object owner = get_impl(S).attr("resid")();
    N_Vector resid   = nb::cast<N_Vector>(owner);

    // SUNLinSolResid() returns a borrowed vector that the caller
    // reads after this trampoline unwinds. Were the Python object holding it
    // a temporary, it would be collected on return and the caller would read
    // freed memory, so retain it here. The reference is dropped on the next
    // resid() call or when the solver is destroyed, which matches the native
    // contract that the residual vector belongs to the solver.
    if (content) { content->resid_owner = std::move(owner); }
    return resid;
  }
  SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr,
                               "CustomSUNLinearSolver.resid", nullptr)
}

SUNErrCode CustomSUNLinearSolver::custom_linsol_free(SUNLinearSolver S)
{
  if (!S) { return SUN_SUCCESS; }
  Content* content = get_content(S);
  custom_content_destroy(content);
  S->content = nullptr;
  SUNLinSolFreeEmpty(S);
  return SUN_SUCCESS;
}

} // namespace sundials4py
