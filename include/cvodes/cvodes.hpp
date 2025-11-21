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
 * C++ specific CVODES definitions.
 * ---------------------------------------------------------------------------*/

#ifndef _CVODES_HPP
#define _CVODES_HPP

#include <sundials/sundials_classview.hpp>

#include <cvodes/cvodes.h>

namespace sundials {
namespace experimental {

struct CVodeDeleter
{
  void operator()(void* v) { CVodeFree(&v); }
};

using CVodeView = ClassView<void*, CVodeDeleter>;

} // namespace experimental
} // namespace sundials

#endif
