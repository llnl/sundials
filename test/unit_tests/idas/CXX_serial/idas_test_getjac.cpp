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
 * Kvaerno-Prothero-Robinson test problem:
 *
 *   r[0] = [ a  b ] [ (-1 + u^2 - r(t)) ] + [ r'(t) ] - [2 u u']
 *   r[1]   [ c  d ] [ (-2 + v^2 - s(t)) ]   [ s'(t) ] - [2 v v']
 *
 * This problem has analytical solution given by
 *
 *   u(t) = sqrt(1 + r(t))
 *   v(t) = sqrt(2 + s(t))
 *
 * where, in this test, we use the functions
 *
 *   r(t) = 0.5 * cos(t)    r'(t) = -0.5 * sin(t)
 *   s(t) = cos(20 t)       s'(t) = -20 sin(20 t)
 *
 * Note: The setup function in dense linear solver is disabled and the solve
 * function is overridden with a custom function that does the setup and solve
 * using a copy of the input matrix. This is done do the internally computed
 * Jacobian matrix is not overwritten with the factorization.
 * ---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

// Include desired integrators, vectors, linear solvers, and nonlinear solvers
#include "idas/idas.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_core.hpp"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

// Macros for problem constants
#define ZERO   SUN_RCONST(0.0)
#define HALF   SUN_RCONST(0.5)
#define ONE    SUN_RCONST(1.0)
#define TWO    SUN_RCONST(2.0)
#define TWENTY SUN_RCONST(20.0)

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

// Compute r(t)
static sunrealtype r(sunrealtype t) { return HALF * cos(t); }

// Compute the derivative of r(t)
static sunrealtype rdot(sunrealtype t) { return -HALF * sin(t); }

// Compute s(t)
static sunrealtype s(sunrealtype t) { return cos(TWENTY * t); }

// Compute the derivative of s(t)
static sunrealtype sdot(sunrealtype t) { return -TWENTY * sin(TWENTY * t); }

// Compute the true solution
static int ytrue(sunrealtype t, N_Vector y)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);

  ydata[0] = sqrt(ONE + r(t));
  ydata[1] = sqrt(TWO + s(t));

  return 0;
}

// Compute the derivative of the true solution
static int yptrue(sunrealtype t, N_Vector yp)
{
  sunrealtype* ypdata = N_VGetArrayPointer(yp);

  ypdata[0] = HALF * rdot(t) / sqrt(ONE + r(t));
  ypdata[1] = HALF * sdot(t) / sqrt(TWO + s(t));

  return 0;
}

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS integrator
// -----------------------------------------------------------------------------

/* -----------------------------------------------------------------------------
 * Compute the DAE residual function:
 *   [a  b] * [ (-1 + u^2 - r(t)) ] + [ r'(t) ] - [ 2 u u'] = 0
 *   [c  d]   [ (-2 + v^2 - s(t)) ]   [ s'(t) ] - [ 2 v v'] = 0
 * ---------------------------------------------------------------------------*/
static int res(sunrealtype t, N_Vector y, N_Vector yp, N_Vector res,
               void* user_data)
{
  sunrealtype* udata  = (sunrealtype*)user_data;
  const sunrealtype a = udata[0];
  const sunrealtype b = udata[1];
  const sunrealtype c = udata[2];
  const sunrealtype d = udata[3];

  sunrealtype* ydata  = N_VGetArrayPointer(y);
  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];

  sunrealtype* ypdata  = N_VGetArrayPointer(yp);
  const sunrealtype up = ypdata[0];
  const sunrealtype vp = ypdata[1];

  const sunrealtype tmp1 = -ONE + u * u - r(t);
  const sunrealtype tmp2 = -TWO + v * v - s(t);

  sunrealtype* rdata = N_VGetArrayPointer(res);

  rdata[0] = a * tmp1 + b * tmp2 + rdot(t) - TWO * u * up;
  rdata[1] = c * tmp1 + d * tmp2 + sdot(t) - TWO * v * vp;

  return 0;
}

/* -----------------------------------------------------------------------------
 * Compute the DAE residual Jacobin:
 *   [2 a u - 2 u' - 2 cj u   2 b v                 ]
 *   [2 c u                   2 d v - 2 v' - 2 cj v ]
 * ---------------------------------------------------------------------------*/
static int J(sunrealtype t, sunrealtype cj, N_Vector y, N_Vector yp,
             N_Vector res, SUNMatrix J, void* user_data, N_Vector tmp1,
             N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* udata  = (sunrealtype*)user_data;
  const sunrealtype a = udata[0];
  const sunrealtype b = udata[1];
  const sunrealtype c = udata[2];
  const sunrealtype d = udata[3];

  sunrealtype* ydata  = N_VGetArrayPointer(y);
  sunrealtype* ypdata = N_VGetArrayPointer(yp);
  sunrealtype* Jdata  = SUNDenseMatrix_Data(J);

  const sunrealtype u  = ydata[0];
  const sunrealtype v  = ydata[1];
  const sunrealtype up = ypdata[0];
  const sunrealtype vp = ypdata[1];

  Jdata[0] = TWO * u * (a - cj) - TWO * up;
  Jdata[1] = TWO * c * u;
  Jdata[2] = TWO * b * v;
  Jdata[3] = TWO * v * (d - cj) - TWO * vp;

  return 0;
}

// -----------------------------------------------------------------------------
// Custom linear solver solve function
// -----------------------------------------------------------------------------

static int DenseSetupAndSolve(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                              N_Vector b, sunrealtype tol)
{
  // Create a copy of the matrix for factorization
  SUNMatrix Acpy = SUNMatClone(A);

  // Copy the input matrix
  int flag = SUNMatCopy_Dense(A, Acpy);
  if (flag) { return flag; }

  // Factor the matrix
  flag = SUNLinSolSetup_Dense(S, Acpy);
  if (flag) { return flag; }

  // Solve the system
  flag = SUNLinSolSolve_Dense(S, A, x, b, tol);
  if (flag) { return flag; }

  // Destroy matrix copy
  SUNMatDestroy(Acpy);

  return 0;
}

// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Check function return flag
static int check_flag(int flag, const std::string funcname)
{
  if (!flag) { return 0; }
  if (flag < 0) { std::cerr << "ERROR: "; }
  if (flag > 0) { std::cerr << "WARNING: "; }
  std::cerr << funcname << " returned " << flag << std::endl;
  return 1;
}

// Check if a function returned a NULL pointer
static int check_ptr(void* ptr, const std::string funcname)
{
  if (ptr) { return 0; }
  std::cerr << "ERROR: " << funcname << " returned NULL" << std::endl;
  return 1;
}

// -----------------------------------------------------------------------------
// Test main
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Comparison tolerance
  sunrealtype tol = 100 * std::sqrt(SUN_UNIT_ROUNDOFF);
  if (argc > 1)
  {
#if defined(SUNDIALS_SINGLE_PRECISION)
    tol = std::stof(argv[1]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    tol = std::stod(argv[1]);
#elif defined(SUNDIALS_EXTENDED_PRECISION)
    tol = std::stold(argv[1]);
#else
#error "SUNDIALS precision macro not defined"
#endif
    if (tol <= ZERO)
    {
      std::cerr << "ERROR: Invalid tolerance, tol = " << tol << std::endl;
      return 1;
    }
  }

  // Integration tolerances
  const sunrealtype atol = 100 * SUN_UNIT_ROUNDOFF;
#if defined(SUNDIALS_SINGLE_PRECISION)
  const sunrealtype rtol = SUN_RCONST(1.0e-3);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  const sunrealtype rtol = SUN_RCONST(1.0e-6);
#elif defined(SUNDIALS_EXTENDED_PRECISION)
  const sunrealtype rtol = SUN_RCONST(1.0e-9);
#else
#error "SUNDIALS precision macro not defined"
#endif

  // Create initial condition
  N_Vector y = N_VNew_Serial(2, sunctx);
  if (check_ptr(y, "N_VNew_Serial")) { return 1; }

  N_Vector yp = N_VClone(y);
  if (check_ptr(y, "N_VClone")) { return 1; }

  int flag = ytrue(ZERO, y);
  if (check_flag(flag, "ytrue")) { return 1; }

  flag = yptrue(ZERO, yp);
  if (check_flag(flag, "yptrue")) { return 1; }

  // Create IDA memory structure
  void* ida_mem = IDACreate(sunctx);
  if (check_ptr(ida_mem, "IDACreate")) { return 1; }

  flag = IDAInit(ida_mem, res, ZERO, y, yp);
  if (check_flag(flag, "IDAInit")) { return 1; }

  flag = IDASStolerances(ida_mem, rtol, atol);
  if (check_flag(flag, "IDASStolerances")) { return 1; }

  SUNMatrix A = SUNDenseMatrix(2, 2, sunctx);
  if (check_ptr(A, "SUNDenseMatrix")) { return 1; }

  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_ptr(LS, "SUNLinSol_Dense")) { return 1; }

  // Disable the linear solver setup function and attach custom solve function
  LS->ops->setup = nullptr;
  LS->ops->solve = DenseSetupAndSolve;

  flag = IDASetLinearSolver(ida_mem, LS, A);
  if (check_flag(flag, "IDASetLinearSolver")) { return 1; }

  sunrealtype udata[4] = {-TWO, HALF, HALF, -ONE};

  flag = IDASetUserData(ida_mem, udata);
  if (check_flag(flag, "IDASetUserData")) { return 1; }

  // Initial time and fist output time
  sunrealtype tret = ZERO;
  sunrealtype tout = tret + SUN_RCONST(0.1);

  // Advance one step in time
  flag = IDASolve(ida_mem, tout, &tret, y, yp, IDA_ONE_STEP);
  if (check_flag(flag, "IDASolve")) { return 1; }

  // Get the internal finite difference approximation to J
  SUNMatrix Jdq;
  flag = IDAGetJac(ida_mem, &Jdq);
  if (check_flag(flag, "IDAGetJac")) { return 1; }

  // Get the step, time, and cj for the Jacobian approximation
  long int nst_Jdq;
  flag = IDAGetJacNumSteps(ida_mem, &nst_Jdq);
  if (check_flag(flag, "IDAGetJacNumSteps")) { return 1; }

  sunrealtype t_Jdq;
  flag = IDAGetJacTime(ida_mem, &t_Jdq);
  if (check_flag(flag, "IDAGetJacTime")) { return 1; }

  sunrealtype cj_Jdq;
  flag = IDAGetJacCj(ida_mem, &cj_Jdq);
  if (check_flag(flag, "IDAGetJacCj")) { return 1; }

  // Compute the true Jacobian
  SUNMatrix Jtrue = SUNDenseMatrix(2, 2, sunctx);
  if (check_ptr(Jtrue, "SUNDenseMatrix")) { return 1; }

  flag = ytrue(t_Jdq, y);
  if (check_flag(flag, "ytrue")) { return 1; }

  flag = yptrue(t_Jdq, yp);
  if (check_flag(flag, "yptrue")) { return 1; }

  flag = J(t_Jdq, cj_Jdq, y, yp, nullptr, Jtrue, &udata, nullptr, nullptr,
           nullptr);
  if (check_flag(flag, "J")) { return 1; }

  // Compare finite difference and true Jacobian
  sunrealtype* Jdq_data = SUNDenseMatrix_Data(Jdq);
  if (check_ptr(Jdq_data, "SUNDenseMatrix_Data")) { return 1; }

  sunrealtype* Jtrue_data = SUNDenseMatrix_Data(Jtrue);
  if (check_ptr(Jtrue_data, "SUNDenseMatrix_Data")) { return 1; }

  // Output Jacobian data
  std::cout << std::scientific;
  std::cout << std::setprecision(std::numeric_limits<sunrealtype>::digits10);
  std::cout << "Jac nst = " << nst_Jdq << std::endl;
  std::cout << "Jac t   = " << t_Jdq << std::endl;
  std::cout << "Jac cj  = " << cj_Jdq << std::endl;
  std::cout << std::endl;
  std::cout << std::setw(8) << std::right << "Index" << std::setw(25)
            << std::right << "J DQ" << std::setw(25) << std::right << "J true"
            << std::setw(25) << std::right << "absolute difference"
            << std::setw(25) << std::right << "relative difference" << std::endl;
  for (int i = 0; i < 4 * 25 + 8; i++) { std::cout << "-"; }
  std::cout << std::endl;

  int result         = 0;
  sunindextype ldata = SUNDenseMatrix_LData(Jtrue);
  for (sunindextype i = 0; i < ldata; i++)
  {
    std::cout << std::setw(8) << std::right << i << std::setw(25) << std::right
              << Jdq_data[i] << std::setw(25) << std::right << Jtrue_data[i]
              << std::setw(25) << std::right
              << std::abs(Jdq_data[i] - Jtrue_data[i]) << std::setw(25)
              << std::right
              << std::abs(Jdq_data[i] - Jtrue_data[i]) / Jtrue_data[i]
              << std::endl;
    result += SUNRCompareTol(Jdq_data[i], Jtrue_data[i], tol);
  }

  // Clean up and return with successful completion
  N_VDestroy(y);
  N_VDestroy(yp);
  SUNMatDestroy(A);
  SUNMatDestroy(Jtrue);
  SUNLinSolFree(LS);
  IDAFree(&ida_mem);

  return result;
}

/*---- end of file ----*/
