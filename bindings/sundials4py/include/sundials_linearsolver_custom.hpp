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
 * CustomSUNLinearSolver: the base class Python code subclasses to implement a
 * SUNLinearSolver.
 *----------------------------------------------------------------------------*/

#ifndef _SUNDIALS4PY_LINEARSOLVER_CUSTOM_HPP
#define _SUNDIALS4PY_LINEARSOLVER_CUSTOM_HPP

#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include <sundials/sundials_classview.hpp>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_linearsolver.hpp>

#include "sundials4py_core_types.hpp"
#include "sundials4py_custom_object.hpp"

namespace sundials4py {

/*
 * Python-owned SUNLinearSolver implementation.
 *
 * The vtable mirrors SUNDIALS' optional-operation model: solve() is required,
 * while setup, callbacks, and diagnostics are attached only if the subclass
 * overrides the named Python method. A NULL operation pointer is meaningful to
 * SUNDIALS, so installing a trampoline for an operation the subclass does not
 * implement would change solver behavior rather than merely fail later.
 */
class CustomSUNLinearSolver
{
public:
  enum class HandleState
  {
    unmaterialized,
    materializing,
    materialized
  };

  CustomSUNLinearSolver(std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx,
                        SUNLinearSolver_Type solver_type)
    : sunctx_owner_(std::move(sunctx)), solver_type_(solver_type)
  {}

  std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>> _get_sundials_handle(
    nb::handle self)
  {
    if (!sunctx_owner_)
    {
      throw nb::type_error(
        "CustomSUNLinearSolver base constructor was not initialized");
    }

    if (state_ == HandleState::materializing)
    {
      throw std::runtime_error(
        "reentrant CustomSUNLinearSolver native handle materialization");
    }

    if (state_ == HandleState::materialized) { return solver_; }

    // See CustomSUNMatrix for the lazy materialization pattern; the same cache
    // guarantees that repeated transparent conversions reuse one C handle.
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

  static sunrealtype base_method_real(const char* name)
  {
    PyErr_SetString(PyExc_NotImplementedError, name);
    nb::raise_python_error();
    return SUN_RCONST(0.0);
  }

private:
  struct Content
  {
    SUNDIALS4PY_CUSTOM_CONTENT_MEMBERS(0x53554e4c53505931ULL);

    SUNLinearSolver_Type solver_type{SUNLINEARSOLVER_DIRECT};

    // Keeps the Python object that owns the vector returned by resid() alive,
    // because SUNLinSolResid() hands back a borrowed pointer its caller may
    // read after the trampoline has returned.
    nb::object resid_owner;
  };

  static constexpr const char* label = "CustomSUNLinearSolver";

  static std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>> make_handle(
    nb::handle impl,
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
    SUNLinearSolver_Type solver_type)
  {
    validate_required_methods(impl);

    // Resolve every optional operation before allocating the native shell, so
    // an exception from Python's attribute machinery cannot leak that shell.
    const bool has_set_atimes         = method_overridden(impl, "set_atimes");
    const bool has_set_preconditioner = method_overridden(impl,
                                                          "set_preconditioner");
    const bool has_set_scaling_vectors =
      method_overridden(impl, "set_scaling_vectors");
    const bool has_set_zero_guess = method_overridden(impl, "set_zero_guess");
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
    if (has_set_zero_guess)
    {
      S->ops->setzeroguess = custom_linsol_setzeroguess;
    }
    if (has_initialize) { S->ops->initialize = custom_linsol_initialize; }
    if (has_setup) { S->ops->setup = custom_linsol_setup; }
    if (has_num_iters) { S->ops->numiters = custom_linsol_numiters; }
    if (has_res_norm) { S->ops->resnorm = custom_linsol_resnorm; }
    if (has_resid) { S->ops->resid = custom_linsol_resid; }

    return sundials::experimental::our_make_shared<
      std::remove_pointer_t<SUNLinearSolver>,
      sundials::experimental::SUNLinearSolverDeleter>(S);
  }

  static Content* get_content(SUNLinearSolver S)
  {
    return custom_content_cast<Content>(S);
  }

  static nb::object get_impl(SUNLinearSolver S)
  {
    return custom_content_impl(get_content(S), label);
  }

  static bool method_overridden(nb::handle impl, const char* name)
  {
    return custom_method_overridden<CustomSUNLinearSolver>(impl, name);
  }

  static void validate_required_methods(nb::handle impl)
  {
    if (!method_overridden(impl, "solve"))
    {
      throw nb::type_error(
        "CustomSUNLinearSolver subclass must override solve()");
    }
  }

  static nb::object matrix_arg(SUNMatrix A)
  {
    if (!A) { return nb::none(); }
    return nb::cast(A, nb::rv_policy::reference);
  }

  static SUNLinearSolver_Type custom_linsol_gettype(SUNLinearSolver S)
  {
    Content* content = get_content(S);
    return content ? content->solver_type : SUNLINEARSOLVER_DIRECT;
  }

  static SUNLinearSolver_ID custom_linsol_getid(SUNLinearSolver)
  {
    return SUNLINEARSOLVER_CUSTOM;
  }

  static SUNErrCode custom_linsol_setatimes(SUNLinearSolver S, void* A_data,
                                            SUNATimesFn ATimes)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      Content* content = get_content(S);
      nb::object impl  = get_impl(S);

      // Register revocable state instead of capturing the raw
      // (function, data) pair, and invalidate whatever occupied this slot
      // before. Converting to a plain Python callable keeps SUNDIALS' opaque
      // A_data out of the subclass's interface entirely.
      auto state = content ? content->callbacks.install("atimes", ATimes, A_data)
                           : nullptr;

      nb::object atimes = nb::none();
      if (state)
      {
        atimes = nb::cpp_function(
          [state](N_Vector x, N_Vector y) -> int
          {
            require_valid_callback(state.get(), "ATimes");
            return state->fn(state->data, x, y);
          },
          nb::arg("x"), nb::arg("y"));
      }
      return static_cast<SUNErrCode>(
        nb::cast<int>(impl.attr("set_atimes")(atimes)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_linsol_setpreconditioner(SUNLinearSolver S,
                                                    void* P_data,
                                                    SUNPSetupFn Pset,
                                                    SUNPSolveFn Psol)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      Content* content = get_content(S);
      nb::object impl  = get_impl(S);

      // SUNDIALS sets both preconditioner halves in one call, so both slots are
      // replaced (and any previous adapters revoked) together.
      auto setup_state =
        content ? content->callbacks.install("psetup", Pset, P_data) : nullptr;
      auto solve_state =
        content ? content->callbacks.install("psolve", Psol, P_data) : nullptr;

      nb::object psetup = nb::none();
      nb::object psolve = nb::none();
      if (setup_state)
      {
        psetup = nb::cpp_function(
          [setup_state]() -> int
          {
            require_valid_callback(setup_state.get(), "preconditioner setup");
            return setup_state->fn(setup_state->data);
          });
      }
      if (solve_state)
      {
        psolve = nb::cpp_function(
          [solve_state](N_Vector r, N_Vector z, sunrealtype tol, int lr) -> int
          {
            require_valid_callback(solve_state.get(), "preconditioner solve");
            return solve_state->fn(solve_state->data, r, z, tol, lr);
          },
          nb::arg("r"), nb::arg("z"), nb::arg("tol"), nb::arg("lr"));
      }
      return static_cast<SUNErrCode>(
        nb::cast<int>(impl.attr("set_preconditioner")(psetup, psolve)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_linsol_setscalingvectors(SUNLinearSolver S,
                                                    N_Vector s1, N_Vector s2)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      return static_cast<SUNErrCode>(nb::cast<int>(get_impl(S).attr(
        "set_scaling_vectors")(nb::cast(s1, nb::rv_policy::reference),
                               nb::cast(s2, nb::rv_policy::reference))));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_linsol_setzeroguess(SUNLinearSolver S,
                                               sunbooleantype onoff)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      return static_cast<SUNErrCode>(
        nb::cast<int>(get_impl(S).attr("set_zero_guess")(onoff)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static SUNErrCode custom_linsol_initialize(SUNLinearSolver S)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      return static_cast<SUNErrCode>(
        nb::cast<int>(get_impl(S).attr("initialize")()));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static int custom_linsol_setup(SUNLinearSolver S, SUNMatrix A)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      return nb::cast<int>(get_impl(S).attr("setup")(matrix_arg(A)));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static int custom_linsol_solve(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                                 N_Vector b, sunrealtype tol)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      return nb::cast<int>(get_impl(S).attr(
        "solve")(matrix_arg(A), nb::cast(x, nb::rv_policy::reference),
                 nb::cast(b, nb::rv_policy::reference), tol));
    }
    SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static int custom_linsol_numiters(SUNLinearSolver S)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      return nb::cast<int>(get_impl(S).attr("num_iters")());
    }
    SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr, __func__,
                                 SUN_ERR_EXT_FAIL)
  }

  static sunrealtype custom_linsol_resnorm(SUNLinearSolver S)
  {
    try
    {
      nb::gil_scoped_acquire gil;
      return nb::cast<sunrealtype>(get_impl(S).attr("res_norm")());
    }
    SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr, __func__,
                                 SUN_RCONST(0.0))
  }

  static N_Vector custom_linsol_resid(SUNLinearSolver S)
  {
    try
    {
      nb::gil_scoped_acquire gil;
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
    SUNDIALS4PY_CATCH_AND_REPORT(S ? S->sunctx : nullptr, __func__, nullptr)
  }

  static SUNErrCode custom_linsol_free(SUNLinearSolver S)
  {
    if (!S) { return SUN_SUCCESS; }
    Content* content = get_content(S);
    custom_content_destroy(content);
    S->content = nullptr;
    SUNLinSolFreeEmpty(S);
    return SUN_SUCCESS;
  }

  std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner_;
  SUNLinearSolver_Type solver_type_;
  std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>> solver_;
  HandleState state_{HandleState::unmaterialized};
  int materialization_count_{0};
};

} // namespace sundials4py

#endif // _SUNDIALS4PY_LINEARSOLVER_CUSTOM_HPP
