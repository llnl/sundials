/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * ODE test problem from D. Estep, et al., "An a posteriori-a priori analysis of
 * multiscale operator splitting," SIAM Journal on Numerical Analysis, 2008:
 *
 *   y' = -lambda * y + y^2
 *
 * where lambda is positive and the initial condition is y(0) = 1. The system
 * as the analytic solution:
 *
 *   y(t) = lambda * y(0) / (y(0) - (y(0) - lambda) * exp(lambda * t))
 *
 * Right-hand side (RHS) functions are defined for the full problem as well as
 * an additive partitioning:
 *
 *   Full RHS    = ode_rhs   = -lambda * y + y^2
 *   Partition 1 = ode_rhs_1 = -lambda * y
 *   Partition 2 = ode_rhs_2 = y^2
 * ---------------------------------------------------------------------------*/

#ifndef ESTEP_HPP_
#define ESTEP_HPP_

#include <sundials/sundials_core.hpp>

namespace problems {
namespace estep {

// problem constants
constexpr sunrealtype zero = SUN_RCONST(0.0);
constexpr sunrealtype one  = SUN_RCONST(1.0);
constexpr sunrealtype two  = SUN_RCONST(2.0);

// initial condition
constexpr sunrealtype y0 = one;

// lambda constant
sunrealtype problem_data = two;

inline int initial_condition(N_Vector y_vec)
{
  if (y_vec == nullptr) { return 1; }

  sunrealtype* y_data = N_VGetArrayPointer(y_vec);
  if (y_data == nullptr) { return 1; }

  y_data[0] = y0;

  return 0;
}

inline int true_solution(sunrealtype t, sunrealtype lambda, N_Vector y_vec)
{
  if (y_vec == nullptr) { return 1; }

  sunrealtype* y_data = N_VGetArrayPointer(y_vec);
  if (y_data == nullptr) { return 1; }

  y_data[0] = lambda * y0 / (y0 - (y0 - lambda) * std::exp(lambda * t));

  return 0;
}

inline int ode_rhs(sunrealtype t, N_Vector y_vec, N_Vector f_vec, void* user_data)
{
  if (y_vec == nullptr || f_vec == nullptr || user_data == nullptr)
  {
    return 1;
  }

  sunrealtype* y_data = N_VGetArrayPointer(y_vec);
  sunrealtype* f_data = N_VGetArrayPointer(f_vec);
  if (y_data == nullptr || f_data == nullptr) { return 1; }

  const sunrealtype lambda = *static_cast<sunrealtype*>(user_data);

  f_data[0] = -lambda * y_data[0] + y_data[0] * y_data[0];

  return 0;
}

inline int ode_rhs_1(sunrealtype t, N_Vector y_vec, N_Vector f_vec,
                     void* user_data)
{
  if (y_vec == nullptr || f_vec == nullptr || user_data == nullptr)
  {
    return 1;
  }

  sunrealtype* y_data = N_VGetArrayPointer(y_vec);
  sunrealtype* f_data = N_VGetArrayPointer(f_vec);
  if (y_data == nullptr || f_data == nullptr) { return 1; }

  const sunrealtype lambda = *static_cast<sunrealtype*>(user_data);

  f_data[0] = -lambda * y_data[0];

  return 0;
}

inline int ode_rhs_2(sunrealtype t, N_Vector y_vec, N_Vector f_vec,
                     void* user_data)
{
  if (y_vec == nullptr || f_vec == nullptr) { return 1; }

  sunrealtype* y_data = N_VGetArrayPointer(y_vec);
  sunrealtype* f_data = N_VGetArrayPointer(f_vec);
  if (y_data == nullptr || f_data == nullptr) { return 1; }

  f_data[0] = y_data[0] * y_data[0];

  return 0;
}

} // namespace estep
} // namespace problems

#endif
