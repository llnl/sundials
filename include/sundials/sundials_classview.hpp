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
 * Base classes for C++ implementations wrappers (views) of SUNDIALS objects.
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_CLASSVIEW_HPP
#define _SUNDIALS_CLASSVIEW_HPP

#include <memory>
#include <type_traits>
#include <utility>

#include <sundials/sundials_convertibleto.hpp>

namespace sundials {
namespace experimental {

template<class T, class Deleter>
std::shared_ptr<T> our_make_shared(T* ptr)
{
  return std::shared_ptr<T>(ptr, Deleter{});
}

template<class T, class Deleter>
class ClassView : public sundials::ConvertibleTo<T>
{
public:
  static_assert(std::is_pointer<T>::value, "ClassView type must be a pointer");

  ClassView(T object = nullptr) noexcept : object_(object) {}

  ClassView(const ClassView&) = delete;

  ClassView(ClassView&& other) = default;

  ClassView& operator=(const ClassView&) = delete;

  ClassView& operator=(ClassView&& rhs) = default;

  ~ClassView() = default;

  // Override ConvertibleTo functions
  T get() noexcept override { return object_.get(); }

  T get() const noexcept override { return object_.get(); }

  operator T() noexcept override { return object_.get(); }

  operator T() const noexcept override { return object_.get(); }

protected:
  std::unique_ptr<std::remove_pointer_t<T>, Deleter> object_;
};

} // namespace experimental
} // namespace sundials

#endif
