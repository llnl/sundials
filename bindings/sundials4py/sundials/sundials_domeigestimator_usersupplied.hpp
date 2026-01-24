/*------------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *------------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2025, Lawrence Livermore National Security,
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
 *----------------------------------------------------------------------------*/

#ifndef _SUNDIALS4PY_SUNDOMEIGESTIMATOR_USERSUPPLIED_HPP
#define _SUNDIALS4PY_SUNDOMEIGESTIMATOR_USERSUPPLIED_HPP

#include <cstdlib>
#include <cstring>

#include "sundials4py.hpp"

#include <sundials/sundials_domeigestimator.hpp>

#include "sundials4py_helpers.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

struct SUNDomEigEstimatorFunctionTable
{
  nb::object atimes;
};

template<typename... Args>
SUNErrCode sundomeigestimator_atimes_wrapper(Args... args)
{
  return sundials4py::user_supplied_fn_caller<
    std::remove_pointer_t<SUNATimesFn>, SUNDomEigEstimatorFunctionTable,
    3>(&SUNDomEigEstimatorFunctionTable::atimes, std::forward<Args>(args)...);
}

#endif