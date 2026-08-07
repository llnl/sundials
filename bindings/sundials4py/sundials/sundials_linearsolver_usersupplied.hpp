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
};

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
