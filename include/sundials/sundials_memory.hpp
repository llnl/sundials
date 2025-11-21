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
 * C++ view of SUNDIALS SUNMemoryHelper
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_MEMORY_HPP
#define _SUNDIALS_MEMORY_HPP

#include <utility>

#include <sundials/sundials_base.hpp>
#include <sundials/sundials_classview.hpp>
#include <sundials/sundials_memory.h>

namespace sundials {
namespace impl {
using BaseMemoryHelper = BaseObject<SUNMemoryHelper_, SUNMemoryHelper_Ops_>;
} // namespace impl

namespace experimental {
struct SUNMemoryHelperDeleter
{
  void operator()(SUNMemoryHelper helper) { SUNMemoryHelper_Destroy(helper); }
};

using SUNMemoryHelperView = ClassView<SUNMemoryHelper, SUNMemoryHelperDeleter>;
} // namespace experimental
} // namespace sundials

#endif
