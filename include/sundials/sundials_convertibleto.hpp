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
 * Base class for converting C++ wrappers (views) to SUNDIALS objects
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_CONVERTIBLETO_HPP
#define _SUNDIALS_CONVERTIBLETO_HPP

namespace sundials {

template<class T>
class ConvertibleTo
{
public:
  // Explicit conversion to the underlying type
  virtual T get() noexcept       = 0;
  virtual T get() const noexcept = 0;

  [[deprecated("This function will be removed in the next major release, use "
               "get() instead.")]] T
  Convert() noexcept
  {
    return get();
  }

  [[deprecated("This function will be removed in the next major release, use "
               "get() instead.")]] T
  Convert() const noexcept
  {
    return get();
  }

  // Implicit conversion to the underlying type
  virtual operator T()       = 0;
  virtual operator T() const = 0;

  virtual ~ConvertibleTo() = default;
};

} // namespace sundials

#endif // _SUNDIALS_CONVERTIBLETO_HPP
