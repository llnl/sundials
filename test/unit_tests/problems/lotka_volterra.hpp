/* -----------------------------------------------------------------------------
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
 * This header provides right-hand-side and related functions (e.g., Jacobian)
 * for the four parameter Lotka-Volterra problem,
 *
 *     u = [dx/dt] = [ p_0*x - p_1*x*y  ]
 *         [dy/dt]   [ -p_2*y + p_3*x*y ].
 *
 * with parameters p.
 * ---------------------------------------------------------------------------*/

#ifndef _LOTKA_VOLTERRA_HPP
#define _LOTKA_VOLTERRA_HPP

#include <sundials/sundials_core.hpp>
#include <sunmatrix/sunmatrix_dense.h>

namespace problems {
namespace lotka_volterra {

inline int ode_rhs(sunrealtype t, N_Vector uvec, N_Vector udotvec, void* user_data)
{
  sunrealtype* p    = (sunrealtype*)user_data;
  sunrealtype* u    = N_VGetArrayPointer(uvec);
  sunrealtype* udot = N_VGetArrayPointer(udotvec);

  udot[0] = p[0] * u[0] - p[1] * u[0] * u[1];
  udot[1] = -p[2] * u[1] + p[3] * u[0] * u[1];

  return 0;
}

inline int ode_jac(sunrealtype t, N_Vector uvec, N_Vector udotvec, SUNMatrix Jac,
                   void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* p = (sunrealtype*)user_data;
  sunrealtype* u = N_VGetArrayPointer(uvec);
  sunrealtype* J = SUNDenseMatrix_Data(Jac);

  J[0] = p[0] - p[1] * u[1];
  J[2] = -p[1] * u[0];
  J[1] = p[3] * u[1];
  J[3] = p[3] * u[0] - p[2];

  return 0;
}

inline int ode_vjp(N_Vector vvec, N_Vector Jvvec, sunrealtype t, N_Vector uvec,
                   N_Vector udotvec, void* user_data, N_Vector tmp)
{
  sunrealtype* p  = (sunrealtype*)user_data;
  sunrealtype* u  = N_VGetArrayPointer(uvec);
  sunrealtype* v  = N_VGetArrayPointer(vvec);
  sunrealtype* Jv = N_VGetArrayPointer(Jvvec);

  Jv[0] = (p[0] - p[1] * u[1]) * v[0] + p[3] * u[1] * v[1];
  Jv[1] = -p[1] * u[0] * v[0] + (-p[2] + p[3] * u[0]) * v[1];

  return 0;
}

inline int parameter_jacobian(sunrealtype t, N_Vector uvec, N_Vector udotvec,
                              SUNMatrix Jac, void* user_data, N_Vector tmp1,
                              N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* u = N_VGetArrayPointer(uvec);
  sunrealtype* J = SUNDenseMatrix_Data(Jac);

  J[0] = u[0];
  J[1] = SUN_RCONST(0.0);
  J[2] = -u[0] * u[1];
  J[3] = SUN_RCONST(0.0);
  J[4] = SUN_RCONST(0.0);
  J[5] = -u[1];
  J[6] = SUN_RCONST(0.0);
  J[7] = u[0] * u[1];

  return 0;
}

inline int parameter_vjp(N_Vector vvec, N_Vector Jvvec, sunrealtype t,
                         N_Vector uvec, N_Vector udotvec, void* user_data,
                         N_Vector tmp)
{
  sunrealtype* u  = N_VGetArrayPointer(uvec);
  sunrealtype* v  = N_VGetArrayPointer(vvec);
  sunrealtype* Jv = N_VGetArrayPointer(Jvvec);

  Jv[0] = u[0] * v[0];
  Jv[1] = -u[0] * u[1] * v[0];
  Jv[2] = -u[1] * v[1];
  Jv[3] = u[0] * u[1] * v[1];

  return 0;
}

} // namespace lotka_volterra
} // namespace problems

#endif
