// -----------------------------------------------------------------------------
// Programmer(s): David J. Gardner @ LLNL
// -----------------------------------------------------------------------------
// SUNDIALS Copyright Start
// Copyright (c) 2025-2026, Lawrence Livermore National Security,
// University of Maryland Baltimore County, and the SUNDIALS contributors.
// Copyright (c) 2013-2025, Lawrence Livermore National Security
// and Southern Methodist University.
// Copyright (c) 2002-2013, Lawrence Livermore National Security.
// All rights reserved.
//
// See the top-level LICENSE and NOTICE files for details.
//
// SPDX-License-Identifier: BSD-3-Clause
// SUNDIALS Copyright End
// -----------------------------------------------------------------------------
// Solving a nonlinear index-1 DAE with an IMEX Runge-Kutta method
//
// y1' = -y1 * y3  +  sin(t)
// y2' = -y2 * y3  +  cos(t)
//  0  =  y3 - y1^2 - y2^2
//
// M y' = f_e(t,y) + f_i(t,y)
//
// [ 1 0 0 ] [ y1' ]   [ sin(t) ]   [ -y1 * y3 ]
// [ 0 1 0 ] [ y2' ] = [ cos(t) ] + [ -y2 * y3 ]
// [ 0 0 0 ] [ y3' ]   [   0    ]   [ y3 - y1^2 - y2^2 ]
//
//                 [ -y3     0     -y1  ]
// J_i = df_i/dy = [  0     -y3    -y2  ]
//                 [ -2*y1  -2*y2    1  ]
//
// y1(0) = 1
// y2(0) = 0
// y3(0) = 1 consistent: y3 = y1^2 + y2^2 = 1 + 0 = 1
// -----------------------------------------------------------------------------

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

// SUNDIALS headers
#include <arkode/arkode_arkstep.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.hpp>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "arkode/arkode.h"
#include "arkode/arkode_butcher.h"
#include "arkode/arkode_ls.h"
#include "sundials/sundials_nvector.h"

using namespace std;

class DAEProblem
{
public:
  // Get number of equations
  int getNumEquations() const { return 3; }

  // Set initial conditions
  void setInitialConditions(N_Vector y)
  {
    sunrealtype* ydata = N_VGetArrayPointer(y);
    ydata[0]           = 1.0;
    ydata[1]           = 0.0;
    ydata[2]           = 1.0;
  }

  // Explicit right-hand side function
  int explicitRHS(sunrealtype t, N_Vector y, N_Vector ydot)
  {
    sunrealtype* dydata = N_VGetArrayPointer(ydot);

    dydata[0] = sin(t);
    dydata[1] = cos(t);
    dydata[2] = 0.0;

    return 0;
  }

  // Static wrapper for RHS (to pass to SUNDIALS)
  static int explicitRHSWrapper(sunrealtype t, N_Vector y, N_Vector ydot,
                                void* user_data)
  {
    DAEProblem* problem = static_cast<DAEProblem*>(user_data);
    return problem->explicitRHS(t, y, ydot);
  }

  // Implicit right-hand side function
  int implicitRHS(sunrealtype t, N_Vector y, N_Vector ydot)
  {
    sunrealtype* ydata  = N_VGetArrayPointer(y);
    sunrealtype* dydata = N_VGetArrayPointer(ydot);

    sunrealtype y1 = ydata[0];
    sunrealtype y2 = ydata[1];
    sunrealtype y3 = ydata[2];

    dydata[0] = -y1 * y3;
    dydata[1] = -y2 * y3;
    dydata[2] = y3 - y1 * y1 - y2 * y2;

    return 0; // Success
  }

  // Static wrapper for RHS (to pass to SUNDIALS)
  static int implicitRHSWrapper(sunrealtype t, N_Vector y, N_Vector ydot,
                                void* user_data)
  {
    DAEProblem* problem = static_cast<DAEProblem*>(user_data);
    return problem->implicitRHS(t, y, ydot);
  }

  // Jacobian function (instance method)
  int computeJac(sunrealtype t, N_Vector y, SUNMatrix J)
  {
    sunrealtype* ydata = N_VGetArrayPointer(y);
    sunrealtype* Jdata = SUNDenseMatrix_Data(J);

    sunrealtype y1 = ydata[0];
    sunrealtype y2 = ydata[1];
    sunrealtype y3 = ydata[2];

    // Fill Jacobian (column-major order for SUNDIALS dense matrices)

    // Column 0
    Jdata[0] = -y3;
    Jdata[1] = 0.0;
    Jdata[2] = -2 * y1;

    // Column 1
    Jdata[3] = 0.0;
    Jdata[4] = -y3;
    Jdata[5] = -2 * y2;

    // Column 2
    Jdata[6] = -y1;
    Jdata[7] = -y2;
    Jdata[8] = 1.0;

    return 0; // Success
  }

  // Static wrapper for Jacobian (to pass to SUNDIALS)
  static int jacWrapper(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                        void* user_data, N_Vector tmp1, N_Vector tmp2,
                        N_Vector tmp3)
  {
    DAEProblem* problem = static_cast<DAEProblem*>(user_data);
    return problem->computeJac(t, y, J);
  }

  // Jacobian function (instance method)
  int computeMass(sunrealtype t, SUNMatrix M)
  {
    sunrealtype* Mdata = SUNDenseMatrix_Data(M);

    // Fill Jacobian (column-major order for SUNDIALS dense matrices)

    // Column 0
    Mdata[0] = 1.0;
    Mdata[1] = 0.0;
    Mdata[2] = 0.0;

    // Column 1
    Mdata[3] = 0.0;
    Mdata[4] = 1.0;
    Mdata[5] = 0.0;

    // Column 2
    Mdata[6] = 0.0;
    Mdata[7] = 0.0;
    Mdata[8] = 0.0;

    return 0; // Success
  }

  // Static wrapper for Jacobian (to pass to SUNDIALS)
  static int massWrapper(sunrealtype t, SUNMatrix M, void* user_data,
                         N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
  {
    DAEProblem* problem = static_cast<DAEProblem*>(user_data);
    return problem->computeMass(t, M);
  }
};

int main(void)
{
  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Default problem parameters
  sunrealtype t0 = 0.0;
  sunrealtype tf = 1.0;
  sunrealtype dt = 0.01;
  string output_file = "data.txt";
  int ierr           = 0;

  // Create problem instance
  DAEProblem problem;

  // Create SUNDIALS vector for initial conditions
  N_Vector y = N_VNew_Serial(problem.getNumEquations(), sunctx);
  if (y == nullptr)
  {
    cerr << "Error creating N_Vector" << endl;
    return 1;
  }

  // Set initial conditions
  problem.setInitialConditions(y);

  // Create ARKODE memory structure based on method choice
  void* arkode_mem = ARKStepCreate(DAEProblem::explicitRHSWrapper,
                                   DAEProblem::implicitRHSWrapper, t0, y, sunctx);
  if (arkode_mem == nullptr)
  {
    cerr << "Error creating ARKODE memory" << endl;
    N_VDestroy(y);
    return 1;
  }

  // Set user data (pointer to problem instance)
  ierr = ARKodeSetUserData(arkode_mem, &problem);
  if (ierr)
  {
    cerr << "Error setting user data" << endl;
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    return 1;
  }

  // Attach IMEX Euler method
  ARKodeButcherTable Be = ARKodeButcherTable_Alloc(2, SUNFALSE);
  if (Be == nullptr)
  {
    cerr << "Error creating Butcher table" << endl;
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    return 1;
  }

  Be->A[1][0] = 1.0;
  Be->b[0]    = 1.0;
  Be->c[1]    = 1.0;
  Be->q       = 1;

  ARKodeButcherTable Bi = ARKodeButcherTable_Alloc(2, SUNFALSE);
  if (Bi == nullptr)
  {
    cerr << "Error creating Butcher table" << endl;
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    return 1;
  }

  Bi->A[1][1] = 1.0;
  Bi->b[1]    = 1.0;
  Bi->c[1]    = 1.0;
  Bi->q       = 1;

  ierr = ARKStepSetTables(arkode_mem, 1, 0, Bi, Be);
  if (ierr)
  {
    cerr << "Error setting Butcher tables" << endl;
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    return 1;
  }

  ARKodeButcherTable_Free(Be);
  ARKodeButcherTable_Free(Bi);

  // Set fixed step size
  ierr = ARKodeSetFixedStep(arkode_mem, dt);
  if (ierr)
  {
    cerr << "Error setting fixed step size" << endl;
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    return 1;
  }

  // Set tolerances
  sunrealtype reltol = 1e-6;
  sunrealtype abstol = 1e-8;
  ierr               = ARKodeSStolerances(arkode_mem, reltol, abstol);
  if (ierr)
  {
    cerr << "Error setting tolerances" << endl;
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    return 1;
  }

  // Create and attach mass matrix
  SUNMatrix M = SUNDenseMatrix(problem.getNumEquations(),
                               problem.getNumEquations(), sunctx);
  if (M == nullptr)
  {
    cerr << "Error creating mass matrix" << endl;
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    return 1;
  }

  ierr = ARKodeSetMassMatrix(arkode_mem, M, DAEProblem::massWrapper,
                             false, true);
  if (ierr)
  {
    cerr << "Error setting mass matrix" << endl;
    SUNMatDestroy(M);
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    return 1;
  }

  // Create dense matrix and linear solver
  SUNMatrix J = SUNDenseMatrix(problem.getNumEquations(),
                               problem.getNumEquations(), sunctx);
  if (J == nullptr)
  {
    cerr << "Error creating Jacobian matrix" << endl;
    SUNMatDestroy(M);
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    return 1;
  }

  SUNLinearSolver LS = SUNLinSol_Dense(y, J, sunctx);
  if (LS == nullptr)
  {
    cerr << "Error creating linear solver" << endl;
    SUNMatDestroy(J);
    SUNMatDestroy(M);
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    return 1;
  }

  // Attach linear solver
  ierr = ARKodeSetLinearSolver(arkode_mem, LS, J);
  if (ierr)
  {
    cerr << "Error setting linear solver" << endl;
    SUNMatDestroy(J);
    SUNMatDestroy(M);
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    return 1;
  }

  // Set Jacobian function
  ierr = ARKodeSetJacFn(arkode_mem, DAEProblem::jacWrapper);
  if (ierr)
  {
    cerr << "Error setting Jacobian function" << endl;
    SUNMatDestroy(J);
    SUNMatDestroy(M);
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    return 1;
  }

  ierr = ARKodeSetInterpolantType(arkode_mem, ARK_INTERP_LAGRANGE);
  if (ierr)
  {
    cerr << "Error setting interpolation type" << endl;
    SUNMatDestroy(J);
    SUNMatDestroy(M);
    ARKodeFree(&arkode_mem);
    N_VDestroy(y);
    return 1;
  }

  // Open file for writing data
  ofstream datafile(output_file);
  datafile << setprecision(17) << scientific;
  datafile << "# t, y1, y2, y3\n";

  // Initial output
  sunrealtype* ydata = N_VGetArrayPointer(y);

  datafile << setw(26) << t0 << setw(26) << ydata[0] << setw(26) << ydata[1]
           << setw(26) << ydata[2] << endl;

  // Time integration loop
  sunrealtype t    = t0;
  sunrealtype tout = t0 + dt;

  while (t < tf)
  {
    int flag = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_ONE_STEP);

    if (flag < 0)
    {
      cerr << "ARKODE error, flag = " << flag << endl;
      break;
    }

    // Step output
    datafile << setw(26) << t << setw(26) << ydata[0] << setw(26) << ydata[1]
             << setw(26) << ydata[2] << endl;

    tout += dt;
    if (tout > tf) tout = tf;
  }

  datafile.close();

  // Print solver statistics
  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  // Clean up
  ARKodeFree(&arkode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(J);
  SUNMatDestroy(M);
  N_VDestroy(y);

  return 0;
}
