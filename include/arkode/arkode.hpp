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
 * C++ specific ARKODE definitions.
 * ---------------------------------------------------------------------------*/

#ifndef _ARKODE_HPP
#define _ARKODE_HPP

#include <sundials/sundials_classview.hpp>

#include <arkode/arkode.h>

namespace sundials {
namespace experimental {

struct ARKodeDeleter
{
  void operator()(void* v) { ARKodeFree(&v); }
};

using ARKodeView = ClassView<void*, ARKodeDeleter>;

class ARKodeBorrowedView : public sundials::ConvertibleTo<void*>
{
public:
  explicit ARKodeBorrowedView(void* object = nullptr) noexcept : object_(object)
  {}

  void* get() noexcept override { return object_; }

  void* get() const noexcept override { return object_; }

  operator void*() noexcept override { return object_; }

  operator void*() const noexcept override { return object_; }

private:
  void* object_;
};

struct ARKodeButcherTableDeleter
{
  void operator()(ARKodeButcherTable t) { ARKodeButcherTable_Free(t); }
};

} // namespace experimental
} // namespace sundials

#endif
