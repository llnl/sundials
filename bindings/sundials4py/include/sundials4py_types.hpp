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

#ifndef _SUNDIALS4PY_TYPES_HPP
#define _SUNDIALS4PY_TYPES_HPP

#include <stdexcept>

#include <sundials/sundials_types.h>

#include "sundials4py.hpp"

namespace nb = nanobind;

namespace sundials4py {

using Array1d = nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig>;

class error_returned : public std::runtime_error
{
public:
  explicit error_returned(const char* message)
    : std::runtime_error(base_message + message)
  {}

  // Constructor that takes a std::string message
  explicit error_returned(const std::string& message)
    : std::runtime_error(base_message + message)
  {}

private:
  inline static const std::string base_message =
    "[sundials4py] a SUNDIALS function returned a code indicating an error, "
    "details are given below:\n\t";
};

class illegal_value : public std::runtime_error
{
public:
  explicit illegal_value(const char* message)
    : std::runtime_error(base_message + message)
  {}

  // Constructor that takes a std::string message
  explicit illegal_value(const std::string& message)
    : std::runtime_error(base_message + message)
  {}

private:
  inline static const std::string base_message =
    "[sundials4py] an illegal value was given, "
    "details are given below:\n\t";
};

class null_function_table : public std::runtime_error
{
public:
  explicit null_function_table(const char* message)
    : std::runtime_error(base_message + message)
  {}

  // Constructor that takes a std::string message
  explicit null_function_table(const std::string& message)
    : std::runtime_error(base_message + message)
  {}

private:
  inline static const std::string base_message =
    "[sundials4py] the python function table was null:\n\t";
};

} // namespace sundials4py

#endif
