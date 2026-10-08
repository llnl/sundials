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

#ifndef SUNDIALS4PY_LINEARSOLVER_CUSTOM_HPP
#define SUNDIALS4PY_LINEARSOLVER_CUSTOM_HPP

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
  : public CustomObjectBase<CustomSUNLinearSolver, std::remove_pointer_t<SUNLinearSolver>>
{
public:
  CustomSUNLinearSolver(std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx,
                        int solver_type)
    : CustomObjectBase(std::move(sunctx), "CustomSUNLinearSolver"),
      solver_type_(checked_type(solver_type))
  {}

  void validate(nb::handle self) const { validate_required_methods(self); }

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
  static SUNLinearSolver_Type checked_type(int solver_type)
  {
    switch (solver_type)
    {
    case SUNLINEARSOLVER_DIRECT:
    case SUNLINEARSOLVER_ITERATIVE:
    case SUNLINEARSOLVER_MATRIX_ITERATIVE:
    case SUNLINEARSOLVER_MATRIX_EMBEDDED:
      return static_cast<SUNLinearSolver_Type>(solver_type);
    default: throw nb::value_error("invalid SUNLinearSolver_Type");
    }
  }

  struct Content : CustomContentBase
  {
    SUNLinearSolver_Type solver_type{SUNLINEARSOLVER_DIRECT};

    // Keeps the Python object that owns the vector returned by resid() alive,
    // because SUNLinSolResid() hands back a borrowed pointer its caller may
    // read after the trampoline has returned.
    nb::object resid_owner;
  };

  static constexpr const char* label = "CustomSUNLinearSolver";

  friend class CustomObjectBase<CustomSUNLinearSolver,
                                std::remove_pointer_t<SUNLinearSolver>>;

  std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>> make_handle(
    nb::handle impl);

  static std::shared_ptr<std::remove_pointer_t<SUNLinearSolver>> make_handle(
    nb::handle impl,
    std::shared_ptr<std::remove_pointer_t<SUNContext>> sunctx_owner,
    SUNLinearSolver_Type solver_type);

  static Content* get_content(SUNLinearSolver S)
  {
    if (!S || !S->ops || S->ops->free != custom_linsol_free || !S->content)
    {
      return nullptr;
    }
    return static_cast<Content*>(S->content);
  }

public:
  static nb::object _python_object_for(SUNLinearSolver S) noexcept
  {
    try
    {
      Content* content = get_content(S);
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
  static nb::object get_impl(SUNLinearSolver S)
  { return custom_content_impl(get_content(S), label); }

  static bool method_overridden(nb::handle impl, const char* name)
  { return custom_method_overridden<CustomSUNLinearSolver>(impl, name); }

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

  static SUNLinearSolver_Type custom_linsol_gettype(SUNLinearSolver S);

  static SUNLinearSolver_ID custom_linsol_getid(SUNLinearSolver);

  static SUNErrCode custom_linsol_setoptions(SUNLinearSolver S, const char* LSid,
                                             const char* file_name, int argc,
                                             char* argv[]);

  static SUNErrCode custom_linsol_setatimes(SUNLinearSolver S, void* A_data,
                                            SUNATimesFn ATimes);

  static SUNErrCode custom_linsol_setpreconditioner(SUNLinearSolver S,
                                                    void* P_data,
                                                    SUNPSetupFn Pset,
                                                    SUNPSolveFn Psol);

  static SUNErrCode custom_linsol_setscalingvectors(SUNLinearSolver S,
                                                    N_Vector s1, N_Vector s2);

  static SUNErrCode custom_linsol_setzeroguess(SUNLinearSolver S,
                                               sunbooleantype onoff);

  static SUNErrCode custom_linsol_initialize(SUNLinearSolver S);

  static int custom_linsol_setup(SUNLinearSolver S, SUNMatrix A);

  static int custom_linsol_solve(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                                 N_Vector b, sunrealtype tol);

  static int custom_linsol_numiters(SUNLinearSolver S);

  static sunrealtype custom_linsol_resnorm(SUNLinearSolver S);

  static N_Vector custom_linsol_resid(SUNLinearSolver S);

  static SUNErrCode custom_linsol_free(SUNLinearSolver S);

  SUNLinearSolver_Type solver_type_;
};

} // namespace sundials4py

#endif // SUNDIALS4PY_LINEARSOLVER_CUSTOM_HPP
