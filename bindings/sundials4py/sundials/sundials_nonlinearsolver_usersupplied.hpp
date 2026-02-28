/* -----------------------------------------------------------------
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
  nb::object lsetupfn;
  nb::object lsolvefn;
  nb::object convtestfn;
};

template<typename... Args>
inline int sunnonlinearsolver_sysfn_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNNonlinSolSysFn>, SUNNonlinearSolverFunctionTable,
    1>(&SUNNonlinearSolverFunctionTable::sysfn, std::forward<Args>(args)...);
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
