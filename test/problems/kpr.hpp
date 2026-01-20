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
 * Kvaerno-Prothero-Robinson (KPR) ODE test problem:
 *
 *   [u]' = [ a  b ] [ (-1 + u^2 - rho(t))   / (2u) ] + [ rho'(t) / (2u) ]
 *   [v]    [ c  d ] [ (-2 + v^2 - sigma(t)) / (2v) ]   [ sigma'(t) / (2v) ]
 *
 * This problem has analytical solution given by
 *
 *   u(t) = sqrt(1 + rho(t))
 *   v(t) = sqrt(2 + sigma(t))
 *
 * where, in this test, we use the functions
 *
 *   rho(t) = 0.5 * cos(t)
 *   sigma(t) = cos(2t)
 *
 * For ImEx methods, the first term is treated implicitly while the second term
 * is the treated explicitly.
 *
 * For MRI methods the u equation is considered slow (potentially with the same
 * ImEx splitting as above) while the v equation is considered fast.
 * ---------------------------------------------------------------------------*/

#ifndef KPR_HPP_
#define KPR_HPP_

#include <cmath>
#include <sunmatrix/sunmatrix_dense.h>

namespace problems {
namespace kpr {

// Macros for problem constants
constexpr sunrealtype zero   = SUN_RCONST(0.0);
constexpr sunrealtype half   = SUN_RCONST(0.5);
constexpr sunrealtype one    = SUN_RCONST(1.0);
constexpr sunrealtype two    = SUN_RCONST(2.0);
constexpr sunrealtype twenty = SUN_RCONST(20.0);

// constants a, b, c, and d
sunrealtype problem_data[4] = {-two, half, half, -one};

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// Compute rho(t)
inline sunrealtype rho(sunrealtype t) { return half * std::cos(t); }

// Compute rho'(t)
inline sunrealtype rho_p(sunrealtype t) { return -half * std::sin(t); }

// Compute sigma(t)
inline sunrealtype sig(sunrealtype t) { return cos(twenty * t); }

// Compute sigma'(t)
inline sunrealtype sig_p(sunrealtype t)
{
  return -twenty * std::sin(twenty * t);
}

// Compute the true solution
inline int true_sol(sunrealtype t, sunrealtype* u, sunrealtype* v)
{
  *u = std::sqrt(one + rho(t));
  *v = std::sqrt(two + sig(t));

  return 0;
}

// Compute the true solution derivative
inline int true_sol_p(sunrealtype t, sunrealtype* up, sunrealtype* vp)
{
  *up = rho_p(t) / (two * std::sqrt(one + rho(t)));
  *vp = sig_p(t) / (two * std::sqrt(two + sig(t)));

  return 0;
}

/* -----------------------------------------------------------------------------
 * ODE RHS function:
 *   [a  b] * [ (-1 + u^2 - r(t)) / (2u) ] + [ r'(t) / (2u) ]
 *   [c  d]   [ (-2 + v^2 - s(t)) / (2v) ]   [ s'(t) / (2v) ]
 * ---------------------------------------------------------------------------*/
inline int ode_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* udata  = (sunrealtype*)user_data;
  const sunrealtype a = udata[0];
  const sunrealtype b = udata[1];
  const sunrealtype c = udata[2];
  const sunrealtype d = udata[3];

  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* fdata = N_VGetArrayPointer(ydot);

  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];

  const sunrealtype tmp1 = (-one + u * u - rho(t)) / (two * u);
  const sunrealtype tmp2 = (-two + v * v - sig(t)) / (two * v);

  fdata[0] = a * tmp1 + b * tmp2 + rho_p(t) / (two * u);
  fdata[1] = c * tmp1 + d * tmp2 + sig_p(t) / (two * v);

  return 0;
}

/* -----------------------------------------------------------------------------
 * ODE RHS Jacobin:
 *   [a/2 + (a(1+r(t))-rdot(t))/(2u^2)     b/2 + b*(2+s(t))/(2*v^2)         ]
 *   [c/2 + c(1+r(t))/(2u^2)               d/2 + (d(2+s(t))-sdot(t))/(2u^2) ]
 * ---------------------------------------------------------------------------*/
inline int ode_rhs_jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                       void* user_data, N_Vector tmp1, N_Vector tmp2,
                       N_Vector tmp3)
{
  sunrealtype* udata  = (sunrealtype*)user_data;
  const sunrealtype a = udata[0];
  const sunrealtype b = udata[1];
  const sunrealtype c = udata[2];
  const sunrealtype d = udata[3];

  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* Jdata = SUNDenseMatrix_Data(J);

  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];

  Jdata[0] = a / two + (a * (one + rho(t)) - rho_p(t)) / (two * u * u);
  Jdata[1] = c / two + c * (one + rho(t)) / (two * u * u);
  Jdata[2] = b / two + b * (two + sig(t)) / (two * v * v);
  Jdata[3] = d / two + (d * (two + sig(t)) - sig_p(t)) / (two * v * v);

  return 0;
}

/* -----------------------------------------------------------------------------
 * ODE Explicit RHS function:
 *   [ r'(t) / (2u) ]
 *   [ s'(t) / (2v) ]
 * ---------------------------------------------------------------------------*/
inline int ode_rhs_ex(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* fdata = N_VGetArrayPointer(ydot);

  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];

  fdata[0] = rho_p(t) / (two * u);
  fdata[1] = sig_p(t) / (two * v);

  return 0;
}

/* -----------------------------------------------------------------------------
 * ODE Implicit RHS function:
 *   [a  b] * [ (-1 + u^2 - r(t)) / (2u) ]
 *   [c  d]   [ (-2 + v^2 - s(t)) / (2v) ]
 * ---------------------------------------------------------------------------*/
inline int ode_rhs_im(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* udata  = (sunrealtype*)user_data;
  const sunrealtype a = udata[0];
  const sunrealtype b = udata[1];
  const sunrealtype c = udata[2];
  const sunrealtype d = udata[3];

  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* fdata = N_VGetArrayPointer(ydot);

  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];

  const sunrealtype tmp1 = (-one + u * u - rho(t)) / (two * u);
  const sunrealtype tmp2 = (-two + v * v - sig(t)) / (two * v);

  fdata[0] = a * tmp1 + b * tmp2;
  fdata[1] = c * tmp1 + d * tmp2;

  return 0;
}

/* -----------------------------------------------------------------------------
 * ODE Implicit RHS Jacobin:
 *   [a/2 + (a(1+r(t)))/(2u^2)     b/2 + b*(2+s(t))/(2v^2)         ]
 *   [c/2 + c(1+r(t))/(2u^2)       d/2 + (d(2+s(t)))/(2u^2) ]
 * ---------------------------------------------------------------------------*/
inline int ode_rhs_jac_im(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                          void* user_data, N_Vector tmp1, N_Vector tmp2,
                          N_Vector tmp3)
{
  sunrealtype* udata  = (sunrealtype*)user_data;
  const sunrealtype a = udata[0];
  const sunrealtype b = udata[1];
  const sunrealtype c = udata[2];
  const sunrealtype d = udata[3];

  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* Jdata = SUNDenseMatrix_Data(J);

  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];

  Jdata[0] = a / two + (a * (one + rho(t))) / (two * u * u);
  Jdata[1] = c / two + c * (one + rho(t)) / (two * u * u);
  Jdata[2] = b / two + b * (two + sig(t)) / (two * v * v);
  Jdata[3] = d / two + (d * (two + sig(t))) / (two * v * v);

  return 0;
}

/* -----------------------------------------------------------------------------
 * ODE Slow RHS function:
 *   [a  b] * [ (-1 + u^2 - r(t)) / (2u) ] + [ r'(t) / (2u) ]
 *   [c  d]   [ 0                        ]   [ 0            ]
 * ---------------------------------------------------------------------------*/
inline int ode_rhs_s(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* udata  = (sunrealtype*)user_data;
  const sunrealtype a = udata[0];
  const sunrealtype b = udata[1];

  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* fdata = N_VGetArrayPointer(ydot);

  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];

  const sunrealtype tmp1 = (-one + u * u - rho(t)) / (two * u);
  const sunrealtype tmp2 = (-two + v * v - sig(t)) / (two * v);

  fdata[0] = a * tmp1 + b * tmp2 + rho_p(t) / (two * u);
  fdata[1] = zero;

  return 0;
}

/* -----------------------------------------------------------------------------
 * ODE Slow Explicit RHS function:
 *   [ r'(t) / (2u) ]
 *   [ 0            ]
 * ---------------------------------------------------------------------------*/
inline int ode_rhs_se(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* fdata = N_VGetArrayPointer(ydot);

  const sunrealtype u = ydata[0];

  fdata[0] = rho_p(t) / (two * u);
  fdata[1] = zero;

  return 0;
}

/* -----------------------------------------------------------------------------
 * ODE Slow Implicit RHS function:
 *   [a  b] * [ (-1 + u^2 - r(t)) / (2u) ]
 *   [c  d]   [ 0                        ]
 * ---------------------------------------------------------------------------*/
inline int ode_rhs_si(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* udata  = (sunrealtype*)user_data;
  const sunrealtype a = udata[0];
  const sunrealtype b = udata[1];

  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* fdata = N_VGetArrayPointer(ydot);

  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];

  const sunrealtype tmp1 = (-one + u * u - rho(t)) / (two * u);
  const sunrealtype tmp2 = (-two + v * v - sig(t)) / (two * v);

  fdata[0] = a * tmp1 + b * tmp2;
  fdata[1] = zero;

  return 0;
}

/* -----------------------------------------------------------------------------
 * ODE Fast RHS function:
 *   [a  b] * [ 0                        ] + [ 0            ]
 *   [c  d]   [ (-2 + v^2 - s(t)) / (2v) ]   [ s'(t) / (2v) ]
 * ---------------------------------------------------------------------------*/
inline int ode_rhs_ff(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* udata  = (sunrealtype*)user_data;
  const sunrealtype c = udata[2];
  const sunrealtype d = udata[3];

  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* fdata = N_VGetArrayPointer(ydot);

  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];

  const sunrealtype tmp1 = (-one + u * u - rho(t)) / (two * u);
  const sunrealtype tmp2 = (-two + v * v - sig(t)) / (two * v);

  fdata[0] = zero;
  fdata[1] = c * tmp1 + d * tmp2 + sig_p(t) / (two * v);

  return 0;
}

/* -----------------------------------------------------------------------------
 * DAE residual function:
 *   ru = [a  b] * [ (-1 + u^2 - r(t)) / (2*u) ] + [ r'(t) / (2u) ] - [u']
 *   rv = [c  d]   [ (-2 + v^2 - s(t)) / (2*v) ]   [ s'(t) / (2v) ] - [v']
 * ---------------------------------------------------------------------------*/
inline int dae_res(sunrealtype t, N_Vector y, N_Vector yp, N_Vector rr,
                   void* user_data)
{
  sunrealtype* udata  = (sunrealtype*)user_data;
  const sunrealtype a = udata[0];
  const sunrealtype b = udata[1];
  const sunrealtype c = udata[2];
  const sunrealtype d = udata[3];

  sunrealtype* ydata  = N_VGetArrayPointer(y);
  sunrealtype* ypdata = N_VGetArrayPointer(yp);
  sunrealtype* rdata  = N_VGetArrayPointer(rr);

  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];

  const sunrealtype up = ypdata[0];
  const sunrealtype vp = ypdata[1];

  const sunrealtype tmp1 = (-one + u * u - rho(t)) / (two * u);
  const sunrealtype tmp2 = (-two + v * v - sig(t)) / (two * v);

  rdata[0] = (a * tmp1 + b * tmp2 + rho_p(t) / (two * u)) - up;
  rdata[1] = (c * tmp1 + d * tmp2 + sig_p(t) / (two * v)) - vp;

  return 0;
}

/* -----------------------------------------------------------------------------
 * DAE residual Jacobian:
 *   [a/2 + (a(1+r(t))-r'(t))/(2u^2) - cj  b/2 + b*(2+s(t))/(2*v^2)            ]
 *   [c/2 + c(1+r(t))/(2u^2)               d/2 + (d(2+s(t))-s'(t))/(2u^2) - cj ]
 * ---------------------------------------------------------------------------*/
inline int dae_res_jac(sunrealtype t, sunrealtype cj, N_Vector y, N_Vector yp,
                       N_Vector rr, SUNMatrix J, void* user_data,
                       N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
  sunrealtype* udata  = (sunrealtype*)user_data;
  const sunrealtype a = udata[0];
  const sunrealtype b = udata[1];
  const sunrealtype c = udata[2];
  const sunrealtype d = udata[3];

  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* Jdata = SUNDenseMatrix_Data(J);

  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];

  Jdata[0] = (a / two + (a * (one + rho(t)) - rho_p(t)) / (two * u * u)) - cj;
  Jdata[1] = c / two + c * (one + rho(t)) / (two * u * u);
  Jdata[2] = b / two + b * (two + sig(t)) / (two * v * v);
  Jdata[3] = (d / two + (d * (two + sig(t)) - sig_p(t)) / (two * v * v)) - cj;

  return 0;
}

} // namespace kpr
} // namespace problems

#endif // KPR_HPP_
