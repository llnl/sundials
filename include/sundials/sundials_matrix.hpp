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
 * C++ view of SUNDIALS SUNMatrix
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_MATRIX_HPP
#define _SUNDIALS_MATRIX_HPP

#include <utility>

#include <sundials/sundials_base.hpp>
#include <sundials/sundials_classview.hpp>
#include <sundials/sundials_matrix.h>

namespace sundials {
namespace impl {
using BaseMatrix = BaseObject<_generic_SUNMatrix, _generic_SUNMatrix_Ops>;
} // namespace impl

namespace experimental {
struct SUNMatrixDeleter
{
  void operator()(SUNMatrix A) { SUNMatDestroy(A); }
};

using SUNMatrixView = ClassView<SUNMatrix, SUNMatrixDeleter>;
} // namespace experimental
} // namespace sundials

#endif
