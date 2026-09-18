/*------------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *                Daniel R. Reynolds @ UMBC
 *------------------------------------------------------------------------------
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
 *------------------------------------------------------------------------------
 * Basic types and exceptions used throughout sundials4py.
 *
 * This is the leaf of the sundials4py header graph: it depends only on nanobind
 * and SUNDIALS, so the per-family custom object headers can include it without
 * creating a cycle with sundials4py_types.hpp (which pulls in all of them).
 *----------------------------------------------------------------------------*/

#ifndef _SUNDIALS4PY_CORE_TYPES_HPP
#define _SUNDIALS4PY_CORE_TYPES_HPP

#include <stdexcept>
#include <string>

#include <sundials/sundials_types.h>

/* nanobind is included directly rather than through sundials4py.hpp: that header
   has no include guard and itself pulls in sundials4py_types.hpp, which pulls in
   this file, so depending on it here would make the header graph order-sensitive.
   This file must remain a leaf. */
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;

namespace sundials4py {

using Array1d = nb::ndarray<sunrealtype, nb::numpy, nb::ndim<1>, nb::c_contig>;
using IntArray1d  = nb::ndarray<int, nb::numpy, nb::ndim<1>, nb::c_contig>;
using LongArray1d = nb::ndarray<long, nb::numpy, nb::ndim<1>, nb::c_contig>;

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

#endif // _SUNDIALS4PY_CORE_TYPES_HPP
