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
// Solving a nonlinear index-1 DAE with BDF methods
//
// y1' = -y1 * y3  +  sin(t)
// y2' = -y2 * y3  +  cos(t)
//  0  =  y3 - y1^2 - y2^2
//
// F(t, y, y') = 0
//
// [ -y1 * y3 + sin(t) ]   [ y1' ]   [ 0 ]
// [ -y2 * y3 + cos(t) ] - [ y2' ] = [ 0 ]
// [  y3 - y1^2 - y2^2 ]   [  0  ]   [ 0 ]
//
//                            [ -y3    0    -y1 ]         [ 1 0 0 ]
// J = dF/dy + alpha dF/dy' = [  0    -y3   -y2 ] - alpha [ 0 1 0 ]
//                            [ -2*y1 -2*y2   1 ]         [ 0 0 0 ]
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
#include <ida/ida.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.hpp>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "sundials/sundials_nvector.h"

using namespace std;

class DAEProblem
{
public:
  // Get number of equations
  int getNumEquations() const { return 3; }

  // Set initial conditions
  void setInitialConditions(N_Vector y, N_Vector yp)
  {
    sunrealtype* ydata = N_VGetArrayPointer(y);
    ydata[0]           = 1.0;
    ydata[1]           = 0.0;
    ydata[2]           = 1.0;

    sunrealtype* ypdata = N_VGetArrayPointer(yp);
    ypdata[0]           = -1.0;
    ypdata[1]           = 1.0;
    ypdata[2]           = 0.0;
  }

  // Explicit right-hand side function
  int computeRes(sunrealtype t, N_Vector y, N_Vector ydot, N_Vector res)
  {
    sunrealtype* ydata   = N_VGetArrayPointer(y);
    sunrealtype* dydata  = N_VGetArrayPointer(ydot);
    sunrealtype* resdata = N_VGetArrayPointer(res);

    sunrealtype y1 = ydata[0];
    sunrealtype y2 = ydata[1];
    sunrealtype y3 = ydata[2];

    sunrealtype dy1 = dydata[0];
    sunrealtype dy2 = dydata[1];

    resdata[0] = -y1 * y3 + sin(t) - dy1;
    resdata[1] = -y2 * y3 + cos(t) - dy2;
    resdata[2] = y3 - y1 * y1 - y2 * y2;

    return 0;
  }

  // Static wrapper for RHS (to pass to SUNDIALS)
  static int resWrapper(sunrealtype t, N_Vector y, N_Vector ydot, N_Vector res,
                        void* user_data)
  {
    DAEProblem* problem = static_cast<DAEProblem*>(user_data);
    return problem->computeRes(t, y, ydot, res);
  }

  // Jacobian function (instance method)
  int computeJac(sunrealtype t, sunrealtype cj, N_Vector y, SUNMatrix J)
  {
    sunrealtype* ydata = N_VGetArrayPointer(y);
    sunrealtype* Jdata = SUNDenseMatrix_Data(J);

    sunrealtype y1 = ydata[0];
    sunrealtype y2 = ydata[1];
    sunrealtype y3 = ydata[2];

    // Fill Jacobian (column-major order for SUNDIALS dense matrices)

    // Column 0
    Jdata[0] = -y3 - cj;
    Jdata[1] = 0.0;
    Jdata[2] = -2 * y1;

    // Column 1
    Jdata[3] = 0.0;
    Jdata[4] = -y3 - cj;
    Jdata[5] = -2 * y2;

    // Column 2
    Jdata[6] = -y1;
    Jdata[7] = -y2;
    Jdata[8] = 1.0;

    return 0; // Success
  }

  // Static wrapper for Jacobian (to pass to SUNDIALS)
  static int jacWrapper(sunrealtype t, sunrealtype cj, N_Vector y, N_Vector yp,
                        N_Vector r, SUNMatrix J, void* user_data, N_Vector tmp1,
                        N_Vector tmp2, N_Vector tmp3)
  {
    DAEProblem* problem = static_cast<DAEProblem*>(user_data);
    return problem->computeJac(t, cj, y, J);
  }
};

int main(void)
{
  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Default problem parameters
  sunrealtype t0     = 0.0;
  sunrealtype tf     = 1.0;
  sunrealtype dt     = 0.01;
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

  N_Vector yp = N_VClone(y);
  if (y == nullptr)
  {
    cerr << "Error creating N_Vector" << endl;
    N_VDestroy(y);
    return 1;
  }

  // Set initial conditions
  problem.setInitialConditions(y, yp);

  // Create IDA memory structure based on method choice
  void* ida_mem = IDACreate(sunctx);
  if (ida_mem == nullptr)
  {
    cerr << "Error creating IDA memory" << endl;
    N_VDestroy(y);
    N_VDestroy(yp);
    return 1;
  }

  ierr = IDAInit(ida_mem, DAEProblem::resWrapper, t0, y, yp);

  // Set user data (pointer to problem instance)
  ierr = IDASetUserData(ida_mem, &problem);
  if (ierr)
  {
    cerr << "Error setting user data" << endl;
    IDAFree(&ida_mem);
    N_VDestroy(y);
    N_VDestroy(yp);
    return 1;
  }

  // Set tolerances
  sunrealtype reltol = 1e-6;
  sunrealtype abstol = 1e-8;
  ierr               = IDASStolerances(ida_mem, reltol, abstol);
  if (ierr)
  {
    cerr << "Error setting tolerances" << endl;
    IDAFree(&ida_mem);
    N_VDestroy(y);
    N_VDestroy(yp);
    return 1;
  }

  // Create dense matrix and linear solver
  SUNMatrix J = SUNDenseMatrix(problem.getNumEquations(),
                               problem.getNumEquations(), sunctx);
  if (J == nullptr)
  {
    cerr << "Error creating Jacobian matrix" << endl;
    IDAFree(&ida_mem);
    N_VDestroy(y);
    N_VDestroy(yp);
    return 1;
  }

  SUNLinearSolver LS = SUNLinSol_Dense(y, J, sunctx);
  if (LS == nullptr)
  {
    cerr << "Error creating linear solver" << endl;
    SUNMatDestroy(J);
    IDAFree(&ida_mem);
    N_VDestroy(y);
    N_VDestroy(yp);
    return 1;
  }

  // Attach linear solver
  ierr = IDASetLinearSolver(ida_mem, LS, J);
  if (ierr)
  {
    cerr << "Error setting linear solver" << endl;
    SUNMatDestroy(J);
    IDAFree(&ida_mem);
    N_VDestroy(y);
    N_VDestroy(yp);
    return 1;
  }

  // Set Jacobian function
  ierr = IDASetJacFn(ida_mem, DAEProblem::jacWrapper);
  if (ierr)
  {
    cerr << "Error setting Jacobian function" << endl;
    SUNMatDestroy(J);
    IDAFree(&ida_mem);
    N_VDestroy(y);
    N_VDestroy(yp);
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
    int flag = IDASolve(ida_mem, tout, &t, y, yp, IDA_NORMAL);

    if (flag < 0)
    {
      cerr << "IDA error, flag = " << flag << endl;
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
  IDAPrintAllStats(ida_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  // Clean up
  IDAFree(&ida_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(J);
  N_VDestroy(y);

  return 0;
}
