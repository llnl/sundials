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
 * Prothero-Robinson (PR) ODE test problem:
 *
 *   y' = lambda * (y - phi(t)) + phi'(t),
 *
 * where phi(t) = atan(t) and phi'(t) = 1 / (1 + t^2). This problem has the
 * analytic solution:
 *
 *   y(t) = atan(t).
 *
 * The stiffness of the problem depends on the value of lambda. For a well-posed
 * problem, the value of lambda should be negative.
 * ---------------------------------------------------------------------------*/

#ifndef PR_HPP_
#define PR_HPP_

#include <cmath>
#include <sundials/sundials_core.hpp>
#include <sunmatrix/sunmatrix_band.h>

namespace problems {
namespace pr {

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// Compute phi(t)
inline sunrealtype phi(sunrealtype t) { return std::atan(t); }

// Compute phi'(t)
inline sunrealtype phi_prime(sunrealtype t)
{
  return SUN_RCONST(1.0) / (SUN_RCONST(1.0) + t * t);
}

// Compute the true solution
inline int true_solution(sunrealtype t, N_Vector y_vec)
{
  if (y_vec == nullptr) { return -1; }

  sunrealtype* y_data = N_VGetArrayPointer(y_vec);
  if (y_data == nullptr) { return -1; }

  const sunrealtype phi_t   = phi(t);
  const sunindextype length = N_VGetLength(y_vec);

  for (sunindextype i = 0; i < length; i++) { y_data[i] = phi_t; }

  return 0;
}

// -----------------------------------------------------------------------------
// Problem functions
// -----------------------------------------------------------------------------

// ODE RHS function
inline int ode_rhs(sunrealtype t, N_Vector y_vec, N_Vector f_vec, void* user_data)
{
  if (y_vec == nullptr || f_vec == nullptr || user_data == nullptr)
  {
    return -1;
  }

  const sunrealtype lambda      = *static_cast<sunrealtype*>(user_data);
  const sunrealtype phi_t       = phi(t);
  const sunrealtype phi_prime_t = phi_prime(t);
  const sunindextype length     = N_VGetLength(y_vec);

  sunrealtype* y_data = N_VGetArrayPointer(y_vec);
  sunrealtype* f_data = N_VGetArrayPointer(f_vec);
  if (y_data == nullptr || f_data == nullptr) { return -1; }

  for (sunindextype i = 0; i < length; i++)
  {
    f_data[i] = lambda * (y_data[i] - phi_t) + phi_prime_t;
  }

  return 0;
}

// ODE RHS Jacobian function
inline int ode_rhs_jac(sunrealtype t, N_Vector y_vec, N_Vector fy_vec,
                       SUNMatrix J_mat, void* user_data, N_Vector tmp1,
                       N_Vector tmp2, N_Vector tmp3)
{
  if (y_vec == nullptr || J_mat == nullptr || user_data == nullptr)
  {
    return -1;
  }

  const sunrealtype lambda  = *static_cast<sunrealtype*>(user_data);
  const sunindextype length = N_VGetLength(y_vec);

  sunrealtype* J_data = SUNBandMatrix_Data(J_mat);
  if (J_data == nullptr) { return -1; }

  for (sunindextype i = 0; i < length; i++) { J_data[i] = lambda; }

  return 0;
}

} // namespace pr
} // namespace problems

#endif
