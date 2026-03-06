/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 * C++ view of SUNDIALS SUNNonlinearSolver
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_NONLINEARSOLVER_HPP
#define _SUNDIALS_NONLINEARSOLVER_HPP

#include <utility>

#include <sundials/sundials_base.hpp>
#include <sundials/sundials_classview.hpp>
#include <sundials/sundials_nonlinearsolver.h>

namespace sundials {
namespace impl {
using BaseNonlinearSolver =
  BaseObject<_generic_SUNNonlinearSolver, _generic_SUNNonlinearSolver_Ops>;
} // namespace impl

namespace experimental {
struct SUNNonlinearSolverDeleter
{
  void operator()(SUNNonlinearSolver NLS)
  {
    if (NLS) { SUNNonlinSolFree(NLS); }
  }
};

using SUNNonlinearSolverView =
  ClassView<SUNNonlinearSolver, SUNNonlinearSolverDeleter>;
} // namespace experimental
} // namespace sundials

#endif
