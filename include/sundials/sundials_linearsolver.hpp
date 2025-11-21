/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
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
 * C++ view of SUNDIALS SUNLinaerSolver
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_LINEARSOLVER_HPP
#define _SUNDIALS_LINEARSOLVER_HPP

#include <utility>

#include <sundials/sundials_base.hpp>
#include <sundials/sundials_classview.hpp>
#include <sundials/sundials_linearsolver.h>

namespace sundials {
namespace impl {
using BaseLinearSolver =
  BaseObject<_generic_SUNLinearSolver, _generic_SUNLinearSolver_Ops>;
} // namespace impl

namespace experimental {
struct SUNLinearSolverDeleter
{
  void operator()(SUNLinearSolver LS) { SUNLinSolFree(LS); }
};

using SUNLinearSolverView = ClassView<SUNLinearSolver, SUNLinearSolverDeleter>;

} // namespace experimental
} // namespace sundials

#endif
