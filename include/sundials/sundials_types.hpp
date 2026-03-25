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
 * C++ view of SUNDIALS SUNStepper
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_TYPES_HPP
#define _SUNDIALS_TYPES_HPP

#include <sundials/sundials_types.h>

#ifdef SUNDIALS_FLOAT128_PRECISION

#include <cstdio>
#include <iomanip>
#include <iostream>
#include <quadmath.h>

/* This defines an output stream operator for the `__float128` type.*/
inline std::ostream& operator<<(std::ostream& os, __float128 value)
{
  // Get current stream formatting state
  const long int width = os.width();         // Width set by std::setw
  const long int precision = os.precision(); // Precision set by std::setprecision
  const std::ios_base::fmtflags flags =
    os.flags(); // Format flags (e.g., scientific notation)

  // Determine format specifier based on stream flags (e/f/g)
  char format_specifier = 'g';
  if (flags & std::ios_base::scientific) { format_specifier = 'e'; }
  else if (flags & std::ios_base::fixed) { format_specifier = 'f'; }

  // Dynamically generate format string (e.g., "%20.15Qe")
  char format_buffer[64];
  std::snprintf(format_buffer, sizeof(format_buffer),
                "%%%d.%dQ%c", // Format template: %[width].[precision]Q[e/f/g]
                static_cast<int>(width),     // Width from setw
                static_cast<int>(precision), // Precision from setprecision
                format_specifier);

  // Format __float128 to string
  char value_buffer[128];
  const auto n = static_cast<long unsigned int>(
    quadmath_snprintf(value_buffer, sizeof(value_buffer),
                      format_buffer, // Dynamically generated format (e.g., "%20.15Qe")
                      value));

  // Write to output stream
  if (n < sizeof(value_buffer)) { os << value_buffer; }
  else { os << "[FORMAT ERROR]"; }

  // Reset stream width (setw has one-time effect)
  os.width(0);

  return os;
}

#endif

#endif
