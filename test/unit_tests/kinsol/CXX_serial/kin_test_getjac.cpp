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
 * Nonlinear system:
 *
 *   3x - cos((y-1)z) - 1/2 = 0
 *   x^2 - 81(y-0.9)^2 + sin(z) + 1.06 = 0
 *   exp(-x(y-1)) + 20z + (10pi-3) / 3 = 0
 *
 * Jacobian:
 *
 *  [ 3                   sin((y-1)z)z    sin((y-1)z)(y-1) ]
 *  [ 2x                  -162(y-0.9)     cos(z)           ]
 *  [ exp(-x(y-1))(1-y)   -exp(-x(y-1))x  20               ]
 *
 * This system has the analytic solution x = 1/2, y = 1, z = -pi/6.
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

// Include KINSOL,vectors, and linear solvers
#include "kinsol/kinsol.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_core.hpp"
#include "sundials/sundials_math.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

// Macros for problem constants
#define ZERO         SUN_RCONST(0.0)
#define PTONE        SUN_RCONST(0.1)
#define HALF         SUN_RCONST(0.5)
#define PTNINE       SUN_RCONST(0.9)
#define ONE          SUN_RCONST(1.0)
#define ONEPTZEROSIX SUN_RCONST(1.06)
#define TWO          SUN_RCONST(2.0)
#define THREE        SUN_RCONST(3.0)
#define TEN          SUN_RCONST(10.0)
#define TWENTY       SUN_RCONST(20.0)
#define EIGHTYONE    SUN_RCONST(81.0)
#define PI           SUN_RCONST(3.141592653589793238462643383279502884197169)

#define SQR(x) ((x) * (x))

// -----------------------------------------------------------------------------
// Functions provided to the SUNDIALS solver
// -----------------------------------------------------------------------------

/* -----------------------------------------------------------------------------
 * Compute the nonlinear residual function:
 *
 * 3x - cos((y-1)z) - 1/2 = 0
 * x^2 - 81(y-0.9)^2 + sin(z) + 1.06 = 0
 * exp(-x(y-1)) + 20z + (10 pi - 3)/3 = 0
 * ---------------------------------------------------------------------------*/
static int res(N_Vector uu, N_Vector fuu, void* user_data)
{
  /* Get vector data arrays */
  sunrealtype* udata = N_VGetArrayPointer(uu);
  sunrealtype* fdata = N_VGetArrayPointer(fuu);

  const sunrealtype x = udata[0];
  const sunrealtype y = udata[1];
  const sunrealtype z = udata[2];

  fdata[0] = THREE * x - std::cos((y - ONE) * z) - HALF;
  fdata[1] = SQR(x) - EIGHTYONE * SQR((y - PTNINE)) + std::sin(z) + ONEPTZEROSIX;
  fdata[2] = std::exp(-x * (y - ONE)) + TWENTY * z + (TEN * PI - THREE) / THREE;

  return 0;
}

/* -----------------------------------------------------------------------------
 * Compute the nonlinear system Jacobian:
 *
 *  [ 3                  sin((y-1)z)z    sin((y-1)z)(y-1) ]
 *  [ 2x                 -162(y-0.9)     cos(z)           ]
 *  [ exp(-x(y-1))(1-y)   -exp(-x(y-1))x  20              ]
 * ---------------------------------------------------------------------------*/

static int J(N_Vector uu, N_Vector fuu, SUNMatrix J, void* user_data,
             N_Vector tmp1, N_Vector tmp2)
{
  sunrealtype* udata = N_VGetArrayPointer(uu);
  sunrealtype* Jdata = SUNDenseMatrix_Data(J);

  const sunrealtype x = udata[0];
  const sunrealtype y = udata[1];
  const sunrealtype z = udata[2];

  // First column
  Jdata[0] = THREE;
  Jdata[1] = TWO * x;
  Jdata[2] = std::exp(-x * (y - ONE)) * (ONE - y);

  // Second column
  Jdata[3] = std::sin((y - ONE) * z) * z;
  Jdata[4] = -TWO * EIGHTYONE * (y - PTNINE);
  Jdata[5] = -std::exp(-x * (y - ONE)) * x;

  // Third column
  Jdata[6] = std::sin((y - ONE) * z) * (y - ONE);
  Jdata[7] = cos(z);
  Jdata[8] = TWENTY;

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
  const sunrealtype ftol = 10 * std::sqrt(SUN_UNIT_ROUNDOFF);

  // Create initial guess and scaling vectors
  N_Vector uu = N_VNew_Serial(3, sunctx);
  if (check_ptr(uu, "N_VNew_Serial")) { return 1; }

  sunrealtype* udata = N_VGetArrayPointer(uu);
  if (check_ptr(uu, "N_VGetArrayPointer")) { return 1; }

  udata[0] = PTONE;
  udata[1] = PTONE;
  udata[2] = -PTONE;

  N_Vector scale = N_VClone(uu);
  if (check_ptr(scale, "N_VNew_Serial")) { return 1; }

  N_VConst(ONE, scale);

  // Create KINSOL memory structure
  void* kin_mem = KINCreate(sunctx);
  if (check_ptr(kin_mem, "KINCreate")) { return 1; }

  int flag = KINInit(kin_mem, res, uu);
  if (check_flag(flag, "KINInit")) { return 1; }

  SUNMatrix A = SUNDenseMatrix(3, 3, sunctx);
  if (check_ptr(A, "SUNDenseMatrix")) { return 1; }

  SUNLinearSolver LS = SUNLinSol_Dense(uu, A, sunctx);
  if (check_ptr(LS, "SUNLinSol_Dense")) { return 1; }

  // Disable the linear solver setup function and attach custom solve function
  LS->ops->setup = nullptr;
  LS->ops->solve = DenseSetupAndSolve;

  flag = KINSetLinearSolver(kin_mem, LS, A);
  if (check_flag(flag, "KINSetLinearSolver")) { return 1; }

  // Disable residual monitoring
  flag = KINSetNoResMon(kin_mem, SUNTRUE);
  if (check_flag(flag, "KINSetNoResMon")) { return 1; }

  // Update Jacobian every iteration
  flag = KINSetMaxSetupCalls(kin_mem, 1);
  if (check_flag(flag, "KINSetMaxSetupCalls")) { return 1; }

  // Set function norm tolerance
  flag = KINSetFuncNormTol(kin_mem, ftol);
  if (check_flag(flag, "KINSetFuncNormTol")) { return 1; }

  // Call main solver
  flag = KINSol(kin_mem, uu, KIN_NONE, scale, scale);
  if (check_flag(flag, "KINSol")) { return 1; }

  long int nni;
  flag = KINGetNumNonlinSolvIters(kin_mem, &nni);
  if (check_flag(flag, "KINGetNumNonlinSolvIters")) { return 1; }

  // Get the internal finite difference approximation to J
  SUNMatrix Jdq;
  flag = KINGetJac(kin_mem, &Jdq);
  if (check_flag(flag, "KINGetJac")) { return 1; }

  // Get the iteration the Jacobian was computed
  long int nni_Jdq;
  flag = KINGetJacNumIters(kin_mem, &nni_Jdq);
  if (check_flag(flag, "KINGetJacNumSteps")) { return 1; }

  // Compute the Jacobian at the returned solution
  SUNMatrix Jtrue = SUNDenseMatrix(3, 3, sunctx);
  if (check_ptr(Jtrue, "SUNDenseMatrix")) { return 1; }

  flag = J(uu, nullptr, Jtrue, nullptr, nullptr, nullptr);
  if (check_flag(flag, "J")) { return 1; }

  // Compare finite difference and true Jacobian
  sunrealtype* Jdq_data = SUNDenseMatrix_Data(Jdq);
  if (check_ptr(Jdq_data, "SUNDenseMatrix_Data")) { return 1; }

  sunrealtype* Jtrue_data = SUNDenseMatrix_Data(Jtrue);
  if (check_ptr(Jtrue_data, "SUNDenseMatrix_Data")) { return 1; }

  // Output Jacobian data
  std::cout << std::scientific;
  std::cout << std::setprecision(std::numeric_limits<sunrealtype>::digits10);
  std::cout << "nni = " << nni << std::endl;
  std::cout << "Jac nni = " << nni_Jdq << std::endl;
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
  N_VDestroy(uu);
  N_VDestroy(scale);
  SUNMatDestroy(A);
  SUNMatDestroy(Jtrue);
  SUNLinSolFree(LS);
  KINFree(&kin_mem);

  return result;
}

/*---- end of file ----*/
