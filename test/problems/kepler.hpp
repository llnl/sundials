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
 * Kepler ODE test problem
 *
 *   q1' = p1
 *   q2' = p2
 *   p1' = -q1 / (q1^2 + q2^2)^(3/2)
 *   p2' = -q2 / (q1^2 + q2^2)^(3/2)
 *
 * with the initial condition
 *
 *   q1(0) = 1 - e
 *   q2(0) = 0
 *   p1(0) = 0
 *   p2(0) = sqrt((1 + e) / (1 - e))
 *
 * where e is the eccentricity. The Hamiltonian for the system is
 *
 *   H(p,q) = 1/2 * (p1^2 + p2^2) - 1 / sqrt(q1^2 + q2^2)
 *
 * is conserved as well as the angular momentum,
 *
 *   L(p,q) = q1 * p2 - q2 * p1.
 *
 * Right-hand side (RHS) functions are defined for the full problem as well as
 * a component partitioning:
 *
 *   Full RHS     = ode_rhs          = RHS terms for q' and p'
 *   Velocity RHS = ode_rhs_velocity = RHS terms for q'
 *   Force RHS    = ode_rhs_force    = RHS terms for p'
 * ---------------------------------------------------------------------------*/

#ifndef KEPLER_HPP_
#define KEPLER_HPP_

#include <cmath>
#include <sundials/sundials_core.hpp>

namespace problems {
namespace kepler {

// Problem constants
constexpr sunrealtype zero = SUN_RCONST(0.0);
constexpr sunrealtype half = SUN_RCONST(0.5);
constexpr sunrealtype one  = SUN_RCONST(1.0);

// eccentricity
sunrealtype eccentricity = SUN_RCONST(0.6);

inline int initial_condition(N_Vector y_vec, sunrealtype ecc)
{
  if (y_vec == nullptr) { return 1; }

  sunrealtype* y_data = N_VGetArrayPointer(y_vec);
  if (y_data == nullptr) { return 1; }

  y_data[0] = one - ecc;
  y_data[1] = zero;
  y_data[2] = zero;
  y_data[3] = std::sqrt((one + ecc) / (one - ecc));

  return 0;
}

inline int hamiltonian(N_Vector y_vec, sunrealtype* H)
{
  if (y_vec == nullptr || H == nullptr) { return 1; }

  sunrealtype* y_data = N_VGetArrayPointer(y_vec);
  if (y_data == nullptr) { return 1; }

  const sunrealtype q1 = y_data[0];
  const sunrealtype q2 = y_data[1];
  const sunrealtype p1 = y_data[2];
  const sunrealtype p2 = y_data[3];

  const sunrealtype qTq = q1 * q1 + q2 * q2;
  const sunrealtype pTp = p1 * p1 + p2 * p2;

  *H = half * pTp - one / std::sqrt(qTq);

  return 0;
}

inline int angular_momentum(N_Vector y_vec, sunrealtype* L)
{
  if (y_vec == nullptr || L == nullptr) { return 1; }

  sunrealtype* y_data = N_VGetArrayPointer(y_vec);
  if (y_data == nullptr) { return 1; }

  const sunrealtype q1 = y_data[0];
  const sunrealtype q2 = y_data[1];
  const sunrealtype p1 = y_data[2];
  const sunrealtype p2 = y_data[3];

  *L = q1 * p2 - q2 * p1;

  return 0;
}

inline int ode_rhs_velocity(sunrealtype t, N_Vector y_vec, N_Vector f_vec,
                            void* user_data)
{
  if (y_vec == nullptr || f_vec == nullptr) { return 1; }

  sunrealtype* y_data = N_VGetArrayPointer(y_vec);
  sunrealtype* f_data = N_VGetArrayPointer(f_vec);
  if (y_data == nullptr || f_data == nullptr) { return 1; }

  const sunrealtype p1 = y_data[2];
  const sunrealtype p2 = y_data[3];

  f_data[0] = p1;
  f_data[1] = p2;

  return 0;
}

inline int ode_rhs_force(sunrealtype t, N_Vector y_vec, N_Vector f_vec,
                         void* user_data)
{
  if (y_vec == nullptr || f_vec == nullptr) { return 1; }

  sunrealtype* y_data = N_VGetArrayPointer(y_vec);
  sunrealtype* f_data = N_VGetArrayPointer(f_vec);
  if (y_data == nullptr || f_data == nullptr) { return 1; }

  const sunrealtype q1 = y_data[0];
  const sunrealtype q2 = y_data[1];

  const sunrealtype sqrt_qTq = std::sqrt(q1 * q1 + q2 * q2);

  f_data[2] = -q1 / std::pow(sqrt_qTq, 3);
  f_data[3] = -q2 / std::pow(sqrt_qTq, 3);

  return 0;
}

inline int ode_rhs(sunrealtype t, N_Vector y_vec, N_Vector f_vec, void* user_data)
{
  int retval = 0;

  retval += ode_rhs_velocity(t, y_vec, f_vec, user_data);
  retval += ode_rhs_force(t, y_vec, f_vec, user_data);

  return retval;
}

} // namespace kepler
} // namespace problems

#endif
