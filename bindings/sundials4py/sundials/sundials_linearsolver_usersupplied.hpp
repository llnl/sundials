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
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS4PY_LINEARSOLVER_USERSUPPLIED_HPP
#define _SUNDIALS4PY_LINEARSOLVER_USERSUPPLIED_HPP

#include <cstdlib>
#include <cstring>

#include "sundials/sundials_iterative.h"
#include "sundials4py.hpp"

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_linearsolver.hpp>

#include "sundials4py_helpers.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

struct SUNLinearSolverFunctionTable
{
  nb::object ATimesFn;
  nb::object PSetupFn;
  nb::object PSolveFn;

  /* The solver's own free operation, displaced by the interposed one below.
     NULL means the solver had no free operation at all. */
  SUNErrCode (*original_free)(SUNLinearSolver){nullptr};
};

/*
 * Free operation installed in place of a solver's own whenever a Python function
 * table is attached to it.
 *
 * SUNLinSolFree() delegates to ops->free and returns immediately; only its
 * fallback path (used when there is no free operation) destroys S->python. No
 * native SUNDIALS linear solver's free operation knows anything about S->python,
 * so without this interposition every function table attached to a native solver
 * would leak, along with every Python callable it holds.
 */
inline SUNErrCode sunlinearsolver_interposed_free(SUNLinearSolver S)
{
  if (S == nullptr) { return SUN_SUCCESS; }

  auto* table = static_cast<SUNLinearSolverFunctionTable*>(S->python);
  SUNErrCode (*original)(SUNLinearSolver) = table ? table->original_free
                                                  : nullptr;

  /* Detach and destroy the table, and put the displaced operation back, before
     doing anything else. Whatever runs below then sees a solver with no Python
     state, so a nested FunctionTable_Destroy(NULL) is a harmless no-op. */
  S->python = nullptr;
  sundials4py::shutdown_safe_delete(table);
  if (S->ops != nullptr) { S->ops->free = original; }

  if (original != nullptr) { return original(S); }

  /* The solver had no free operation of its own, so SUNLinSolFree() would have
     taken its generic fallback path. Since interposing replaced the branch that
     selects that path, reproduce it here. Note that SUNLinSolFreeEmpty frees ops
     and the solver but not content. */
  free(S->content);
  S->content = nullptr;
  SUNLinSolFreeEmpty(S);
  return SUN_SUCCESS;
}

/*
 * Return the solver's Python function table, creating and attaching it (with the
 * free interposition above) on first use.
 *
 * Every hand-written setter that stores a Python callable goes through this, so
 * the attach-and-interpose policy is stated exactly once.
 */
inline SUNLinearSolverFunctionTable* sunlinearsolver_function_table(SUNLinearSolver S)
{
  if (S->python == nullptr)
  {
    auto* table = new SUNLinearSolverFunctionTable;
    if (S->ops != nullptr)
    {
      table->original_free = S->ops->free;
      S->ops->free         = sunlinearsolver_interposed_free;
    }
    S->python = table;
  }
  return static_cast<SUNLinearSolverFunctionTable*>(S->python);
}

template<typename... Args>
inline int sunlinearsolver_atimesfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNATimesFn>, SUNLinearSolverFunctionTable,
    3>(&SUNLinearSolverFunctionTable::ATimesFn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int sunlinearsolver_psetupfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNPSetupFn>, SUNLinearSolverFunctionTable,
    1>(&SUNLinearSolverFunctionTable::PSetupFn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int sunlinearsolver_psolvefn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNPSolveFn>, SUNLinearSolverFunctionTable,
    5>(&SUNLinearSolverFunctionTable::PSolveFn, std::forward<Args>(args)...);
}

#endif // _SUNDIALS4PY_LINEARSOLVER_USERSUPPLIED_HPP
