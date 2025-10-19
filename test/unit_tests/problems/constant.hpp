/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * Test problem with constant derivatives
 *
 *   y' = a
 *
 * with the initial condition
 *
 *   y(t_0) = y_0
 *
 * This problem has analytical solution
 *
 *   y(t) = a * (t - t_0) + y_0
 *
 * ---------------------------------------------------------------------------*/

#ifndef CONSTANT_HPP_
#define CONSTANT_HPP_

#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>

namespace problems {
namespace constant {

struct UserData
{
  sunrealtype t0_;
  N_Vector y0_;
  N_Vector rhs_;

  UserData(sunrealtype t0, N_Vector y0, N_Vector rhs)
    : t0_(t0), y0_(y0), rhs_(rhs){};
};

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// Compute the true solution
inline int true_sol(sunrealtype t, N_Vector y, UserData& user_data)
{
  N_VLinearSum((t - user_data.t0_), user_data.rhs_, SUN_RCONST(1.0),
               user_data.y0_, y);
  return 0;
}

// Compute the true solution derivative
inline int true_sol_p(sunrealtype t, N_Vector yp, UserData& user_data)
{
  N_VScale(SUN_RCONST(1.0), user_data.rhs_, yp);
  return 0;
}

// ODE RHS function, f(t,y) = a
inline int ode_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  auto udata = static_cast<UserData*>(user_data);
  N_VScale(SUN_RCONST(1.0), udata->rhs_, ydot);
  return 0;
}

// ODE RHS Jacobin, J = df/dy = 0
inline int ode_rhs_jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                       void* user_data, N_Vector tmp1, N_Vector tmp2,
                       N_Vector tmp3)
{
  SUNMatZero(J);
  return 0;
}

// DAE residual function, F(t, y, y') = y' - f(t,y) = y' - a
inline int dae_res(sunrealtype t, N_Vector y, N_Vector yp, N_Vector rr,
                   void* user_data)
{
  auto udata = static_cast<UserData*>(user_data);
  N_VLinearSum(SUN_RCONST(1.0), yp, SUN_RCONST(-1.0), udata->rhs_, rr);
  return 0;
}

// DAE residual Jacobian, J = dF/dy + alpha dF/dy' = alpha I
inline int dae_res_jac(sunrealtype t, sunrealtype cj, N_Vector y, N_Vector yp,
                       N_Vector rr, SUNMatrix J, void* user_data,
                       N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
  const sunindextype length = N_VGetLength(y);
  sunrealtype* jac_data     = SUNDenseMatrix_Data(J);

  SUNMatZero(J);
  for (sunindextype i = 0; i < length; ++i) { jac_data[i * length + i] = cj; }

  return 0;
}

} // namespace constant
} // namespace problems

#endif // CONSTANT_HPP_
