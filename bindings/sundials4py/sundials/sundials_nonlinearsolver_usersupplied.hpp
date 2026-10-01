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

#ifndef _SUNDIALS4PY_NONLINEARSOLVER_USERSUPPLIED_HPP
#define _SUNDIALS4PY_NONLINEARSOLVER_USERSUPPLIED_HPP

#include <cstdlib>
#include <cstring>

#include "sundials4py.hpp"

#include <sundials/sundials_nonlinearsolver.h>
#include <sundials/sundials_nonlinearsolver.hpp>

#include "sundials4py_helpers.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

struct SUNNonlinearSolverFunctionTable
{
  nb::object sysfn;
  nb::object rootsysfn;
  nb::object fixedpointsysfn;
  nb::object lsetupfn;
  nb::object lsolvefn;
  nb::object convtestfn;
  nb::object normfn;
  nb::object getupdatenormfn;
  nb::object getconvratefn;

  /* The solver's own free operation, displaced by the interposed one below.
     NULL means the solver had no free operation at all. */
  SUNErrCode (*original_free)(SUNNonlinearSolver){nullptr};
};

/*
 * Free operation installed in place of a solver's own whenever a Python function
 * table is attached to it.
 *
 * SUNNonlinSolFree() delegates to ops->free and returns immediately; only its
 * fallback path (used when there is no free operation) destroys NLS->python. No
 * native SUNDIALS solver's free operation knows anything about NLS->python, so
 * without this interposition every function table attached to a native solver
 * would leak, along with every Python callable it holds.
 */
inline SUNErrCode sunnonlinearsolver_interposed_free(SUNNonlinearSolver NLS)
{
  if (NLS == nullptr) { return SUN_SUCCESS; }

  auto* table = static_cast<SUNNonlinearSolverFunctionTable*>(NLS->python);
  SUNErrCode (*original)(SUNNonlinearSolver) = table ? table->original_free
                                                     : nullptr;

  /* Detach and destroy the table, and put the displaced operation back, before
     doing anything else. Whatever runs below then sees a solver with no Python
     state, so a nested FunctionTable_Destroy(NULL) is a harmless no-op. */
  NLS->python = nullptr;
  sundials4py::shutdown_safe_delete(table);
  if (NLS->ops != nullptr) { NLS->ops->free = original; }

  if (original != nullptr) { return original(NLS); }

  /* The solver had no free operation of its own, so SUNNonlinSolFree() would
     have taken its generic fallback path. Since interposing replaced the branch
     that selects that path, reproduce it here. Note that SUNNonlinSolFreeEmpty
     frees ops and the solver but not content. */
  free(NLS->content);
  NLS->content = nullptr;
  SUNNonlinSolFreeEmpty(NLS);
  return SUN_SUCCESS;
}

/*
 * Return the solver's Python function table, creating and attaching it (with the
 * free interposition above) on first use.
 *
 * Every hand-written setter that stores a Python callable goes through this, so
 * the attach-and-interpose policy is stated exactly once.
 */
inline SUNNonlinearSolverFunctionTable* sunnonlinearsolver_function_table(
  SUNNonlinearSolver NLS)
{
  if (NLS->python == nullptr)
  {
    auto* table = new SUNNonlinearSolverFunctionTable;
    if (NLS->ops != nullptr)
    {
      table->original_free = NLS->ops->free;
      NLS->ops->free       = sunnonlinearsolver_interposed_free;
    }
    NLS->python = table;
  }
  return static_cast<SUNNonlinearSolverFunctionTable*>(NLS->python);
}

template<typename... Args>
inline int sunnonlinearsolver_sysfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNNonlinSolSysFn>, SUNNonlinearSolverFunctionTable,
    1>(&SUNNonlinearSolverFunctionTable::sysfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int sunnonlinearsolver_rootsysfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNNonlinSolSysFn>, SUNNonlinearSolverFunctionTable,
    1>(&SUNNonlinearSolverFunctionTable::rootsysfn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int sunnonlinearsolver_fixedpointsysfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNNonlinSolSysFn>, SUNNonlinearSolverFunctionTable,
    1>(&SUNNonlinearSolverFunctionTable::fixedpointsysfn,
       std::forward<Args>(args)...);
}

using SUNNonlinSolLSetupStdFn = std::tuple<int, sunbooleantype>(sunbooleantype jbad,
                                                                void* mem);

inline int sunnonlinearsolver_lsetupfn_wrapper(sunbooleantype jbad,
                                               sunbooleantype* jcur, void* mem)
{
  auto fn_table = static_cast<SUNNonlinearSolverFunctionTable*>(mem);
  auto fn = nb::cast<std::function<SUNNonlinSolLSetupStdFn>>(fn_table->lsetupfn);

  auto result = fn(jbad, nullptr);

  *jcur = std::get<1>(result);

  return std::get<0>(result);
}

using SUNNonlinSolConvRateStdFn = std::tuple<SUNErrCode, sunrealtype>(void* mem);

using SUNNonlinSolNormStdFn = std::tuple<SUNErrCode, sunrealtype>(N_Vector del,
                                                                  N_Vector w,
                                                                  void* mem);

inline int sunnonlinearsolver_normfn_wrapper(N_Vector del, N_Vector w,
                                             sunrealtype* delnrm, void* mem)
{
  auto fn_table = static_cast<SUNNonlinearSolverFunctionTable*>(mem);
  auto fn = nb::cast<std::function<SUNNonlinSolNormStdFn>>(fn_table->normfn);

  auto result = fn(del, w, nullptr);

  *delnrm = std::get<1>(result);

  return std::get<0>(result);
}

using SUNNonlinSolGetUpdateNormStdFn =
  std::tuple<SUNErrCode, sunrealtype>(void* mem);

inline int sunnonlinearsolver_getupdatenormfn_wrapper(sunrealtype* delnrm,
                                                      void* mem)
{
  auto fn_table = static_cast<SUNNonlinearSolverFunctionTable*>(mem);
  auto fn       = nb::cast<std::function<SUNNonlinSolGetUpdateNormStdFn>>(
    fn_table->getupdatenormfn);

  auto result = fn(nullptr);

  *delnrm = std::get<1>(result);

  return std::get<0>(result);
}

using SUNNonlinSolGetConvRateStdFn = std::tuple<SUNErrCode, sunrealtype>(void* mem);

inline int sunnonlinearsolver_getconvratefn_wrapper(sunrealtype* crate, void* mem)
{
  auto fn_table = static_cast<SUNNonlinearSolverFunctionTable*>(mem);
  auto fn       = nb::cast<std::function<SUNNonlinSolGetConvRateStdFn>>(
    fn_table->getconvratefn);

  auto result = fn(nullptr);

  *crate = std::get<1>(result);

  return std::get<0>(result);
}

template<typename... Args>
inline int sunnonlinearsolver_lsolvefn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNNonlinSolLSolveFn>, SUNNonlinearSolverFunctionTable,
    1>(&SUNNonlinearSolverFunctionTable::lsolvefn, std::forward<Args>(args)...);
}

template<typename... Args>
inline int sunnonlinearsolver_convtestfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNNonlinSolConvTestFn>, SUNNonlinearSolverFunctionTable,
    1>(&SUNNonlinearSolverFunctionTable::convtestfn, std::forward<Args>(args)...);
}

#endif // _SUNDIALS4PY_NONLINEARSOLVER_USERSUPPLIED_HPP
