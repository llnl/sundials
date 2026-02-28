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
 * Kvaerno-Prothero-Robinson ODE test problem:
 *
 *   [u]' = [ a  b ] [ (-1 + u^2 - r(t)) / (2u) ] + [ r'(t) / (2u) ]
 *   [v]    [ c  d ] [ (-2 + v^2 - s(t)) / (2v) ]   [ s'(t) / (2v) ]
 *
 * This problem has analytical solution given by
 *
 *   u(t) = sqrt(1 + r(t))
 *   v(t) = sqrt(2 + s(t))
 *
 * where, in this test, we use the functions
 *
 *   r(t) = 0.5 * cos(t)
 *   s(t) = cos(2t)
 * ---------------------------------------------------------------------------*/

#include "ida_test_kpr.hpp"

int main(int argc, char* argv[])
{
  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Read input options
  TestOptions opts;

  vector<string> args(argv + 1, argv + argc);

  int flag = ReadInputs(args, opts, sunctx);
  if (check_flag(flag, "ReadInputs")) { return 1; }

  // Create initial condition
  N_Vector y = N_VNew_Serial(2, sunctx);
  if (check_ptr(y, "N_VNew_Serial")) { return 1; }

  sunrealtype utrue, vtrue;
  flag = true_sol(ZERO, &utrue, &vtrue);
  if (check_flag(flag, "true_sol")) { return 1; }

  sunrealtype* ydata = N_VGetArrayPointer(y);
  ydata[0]           = utrue;
  ydata[1]           = vtrue;

  N_Vector yp = N_VNew_Serial(2, sunctx);
  if (check_ptr(y, "N_VNew_Serial")) { return 1; }

  sunrealtype uptrue, vptrue;
  flag = true_sol_p(ZERO, &uptrue, &vptrue);
  if (check_flag(flag, "true_sol")) { return 1; }

  sunrealtype* ypdata = N_VGetArrayPointer(yp);
  ypdata[0]           = uptrue;
  ypdata[1]           = vptrue;

  // Create IDA memory structure
  void* ida_mem = IDACreate(sunctx);
  if (check_ptr(ida_mem, "IDACreate")) { return 1; }

  flag = IDAInit(ida_mem, res, ZERO, y, yp);
  if (check_flag(flag, "IDAInit")) { return 1; }

  flag = IDASStolerances(ida_mem, opts.rtol, opts.atol);
  if (check_flag(flag, "IDASStolerances")) { return 1; }

  SUNMatrix A = SUNDenseMatrix(2, 2, sunctx);
  if (check_ptr(A, "SUNDenseMatrix")) { return 1; }

  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_ptr(LS, "SUNLinSol_Dense")) { return 1; }

  flag = IDASetLinearSolver(ida_mem, LS, A);
  if (check_flag(flag, "IDASetLinearSolver")) { return 1; }

  flag = IDASetJacFn(ida_mem, J);
  if (check_flag(flag, "IDASetJacFn")) { return 1; }

  sunrealtype udata[4] = {-TWO, HALF, HALF, -ONE};
  flag                 = IDASetUserData(ida_mem, udata);
  if (check_flag(flag, "IDASetUserData")) { return 1; }

  flag = IDASetEtaFixedStepBounds(ida_mem, opts.eta_min_fx, opts.eta_max_fx);
  if (check_flag(flag, "IDASetEtaFixeStepBounds")) { return 1; }

  flag = IDASetEtaMax(ida_mem, opts.eta_max);
  if (check_flag(flag, "IDASetEtaMax")) { return 1; }

  flag = IDASetEtaMin(ida_mem, opts.eta_min);
  if (check_flag(flag, "IDASetEtaMin")) { return 1; }

  flag = IDASetEtaMinErrFail(ida_mem, opts.eta_min_ef);
  if (check_flag(flag, "IDASetEtaMinErrFail")) { return 1; }

  flag = IDASetEtaConvFail(ida_mem, opts.eta_cf);
  if (check_flag(flag, "IDASetEtaConvFail")) { return 1; }

  flag = IDASetDeltaCjLSetup(ida_mem, opts.dcj);
  if (check_flag(flag, "IDASetDeltaCjLSetup")) { return 1; }

  // Initial time and fist output time
  sunrealtype tret = ZERO;
  sunrealtype tout = tret + opts.dtout;

  // Output initial contion
  cout << scientific;
  cout << setprecision(numeric_limits<sunrealtype>::digits10);
  cout << "           t              ";
  cout << "          u              ";
  cout << "          v              ";
  cout << "        u err            ";
  cout << "        v err            " << endl;
  for (int i = 0; i < 9; i++) { cout << "--------------"; }
  cout << endl;

  cout << setw(22) << tret << setw(25) << ydata[0] << setw(25) << ydata[1]
       << setw(25) << abs(ydata[0] - utrue) << setw(25) << abs(ydata[1] - vtrue)
       << endl;

  // Advance in time
  for (int i = 0; i < opts.nout; i++)
  {
    flag = IDASolve(ida_mem, tout, &tret, y, yp, IDA_NORMAL);
    if (check_flag(flag, "IDA")) { return 1; }

    flag = true_sol(tret, &utrue, &vtrue);
    if (check_flag(flag, "true_sol")) { return 1; }

    cout << setw(22) << tret << setw(25) << ydata[0] << setw(25) << ydata[1]
         << setw(25) << abs(ydata[0] - utrue) << setw(25)
         << abs(ydata[1] - vtrue) << endl;

    // update output time
    tout += opts.dtout;
  }
  for (int i = 0; i < 9; i++) { cout << "--------------"; }
  cout << endl;

  // Print some final statistics
  flag = IDAPrintAllStats(ida_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(flag, "IDAPrintAllStats")) { return 1; }

  // Clean up and return with successful completion
  N_VDestroy(y);
  N_VDestroy(yp);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  IDAFree(&ida_mem);

  return 0;
}

/* -----------------------------------------------------------------------------
 * Compute the DAE residual function:
 *   ru = [a  b] * [ (-1 + u^2 - r(t)) / (2*u) ] + [ r'(t) / (2u) ] - [u']
 *   rv = [c  d]   [ (-2 + v^2 - s(t)) / (2*v) ]   [ s'(t) / (2v) ] - [v']
 * ---------------------------------------------------------------------------*/
int res(sunrealtype t, N_Vector y, N_Vector yp, N_Vector rr, void* user_data)
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

  const sunrealtype tmp1 = (-ONE + u * u - r(t)) / (TWO * u);
  const sunrealtype tmp2 = (-TWO + v * v - s(t)) / (TWO * v);

  sunrealtype* rdata = N_VGetArrayPointer(rr);
  rdata[0]           = (a * tmp1 + b * tmp2 + rdot(t) / (TWO * u)) - up;
  rdata[1]           = (c * tmp1 + d * tmp2 + sdot(t) / (TWO * v)) - vp;

  return 0;
}

/* -----------------------------------------------------------------------------
 * Compute the DAE residual Jacobian:
 *   [a/2 + (a(1+r(t))-r'(t))/(2u^2) - cj  b/2 + b*(2+s(t))/(2*v^2)            ]
 *   [c/2 + c(1+r(t))/(2u^2)               d/2 + (d(2+s(t))-s'(t))/(2u^2) - cj ]
 * ---------------------------------------------------------------------------*/
int J(sunrealtype t, sunrealtype cj, N_Vector y, N_Vector yp, N_Vector rr,
      SUNMatrix J, void* user_data, N_Vector tempv1, N_Vector tempv2,
      N_Vector tempv3)
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

  Jdata[0] = (a / TWO + (a * (ONE + r(t)) - rdot(t)) / (TWO * u * u)) - cj;
  Jdata[1] = c / TWO + c * (ONE + r(t)) / (TWO * u * u);
  Jdata[2] = b / TWO + b * (TWO + s(t)) / (TWO * v * v);
  Jdata[3] = (d / TWO + (d * (TWO + s(t)) - sdot(t)) / (TWO * v * v)) - cj;

  return 0;
}

/*---- end of file ----*/
