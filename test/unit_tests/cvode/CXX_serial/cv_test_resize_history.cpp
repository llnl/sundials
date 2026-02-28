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
 * Test resizing the CVODE history array
 * ---------------------------------------------------------------------------*/

#include <cstdio>
#include <iomanip>
#include <iostream>

#include "cvode/cvode.h"
#include "cvode/cvode_impl.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_context.hpp"
#include "sundials/sundials_linearsolver.h"
#include "sundials/sundials_matrix.h"
#include "sundials/sundials_nvector.h"
#include "sunlinsol/sunlinsol_band.h"
#include "sunmatrix/sunmatrix_band.h"
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"

#include "problems/pr.hpp"
#include "utilities/check_return.hpp"

using namespace std;
using namespace problems::pr;

constexpr sunrealtype zero = SUN_RCONST(0.0);
constexpr sunrealtype one  = SUN_RCONST(1.0);

static int save_history(sunrealtype t_n, N_Vector y_n, sunrealtype* t_hist,
                        N_Vector* y_hist, N_Vector* f_hist, int hist_size,
                        int step, void* user_data)
{
  // Move old values back one overwriting the oldest value and insert the most
  // recent time, solution vector, and rhs vector as the first history entry
  int i_start = hist_size - 2;
  if (step < hist_size) { i_start = step - 1; }

  for (int i = i_start; i >= 0; i--) { t_hist[i + 1] = t_hist[i]; }
  t_hist[0] = t_n;

  for (int i = i_start; i >= 0; i--)
  {
    N_VScale(one, y_hist[i], y_hist[i + 1]);
    N_VScale(one, f_hist[i], f_hist[i + 1]);
  }
  N_VScale(one, y_n, y_hist[0]);
  int retval = ode_rhs(t_n, y_n, f_hist[0], user_data);
  if (retval) { return 1; }

  return 0;
}

static int resize_history(sunrealtype* t_hist, N_Vector* y_hist, N_Vector* f_hist,
                          int hist_size, SUNContext sunctx, void* user_data)
{
  // Add one new element to the solution and rhs history vectors
  int new_size = N_VGetLength(y_hist[0]) + 1;

  for (int i = 0; i < hist_size; i++)
  {
    N_Vector new_vec      = N_VNew_Serial(new_size, sunctx);
    sunrealtype* old_data = N_VGetArrayPointer(y_hist[i]);
    sunrealtype* new_data = N_VGetArrayPointer(new_vec);
    for (int j = 0; j < new_size; j++) { new_data[j] = old_data[0]; }
    N_VDestroy(y_hist[i]);
    y_hist[i] = new_vec;
  }

  for (int i = 0; i < hist_size; i++)
  {
    N_VDestroy(f_hist[i]);
    f_hist[i]  = N_VClone(y_hist[i]);
    int retval = ode_rhs(t_hist[i], y_hist[i], f_hist[i], user_data);
    if (retval) { return 1; }
  }

  return 0;
}

static int PrintNordsieckArray(void* cvode_mem, bool print_all)
{
  // Access private CVODE memory to print Nordsieck array
  CVodeMem cv_mem = (CVodeMem)cvode_mem;

  const sunindextype N = N_VGetLength(cv_mem->cv_zn[0]);
  sunrealtype* vdata   = nullptr;

  cout << setw(4) << "idx" << setw(25) << "zn" << endl;

  int i_max = cv_mem->cv_q;
  if (print_all) { i_max = cv_mem->cv_qmax; }

  for (int i = 0; i <= i_max; i++)
  {
    vdata = N_VGetArrayPointer(cv_mem->cv_zn[i]);
    for (sunindextype j = 0; j < N; j++)
    {
      cout << setw(4) << i << setw(25) << vdata[j] << endl;
    }
  }
  cout << endl;

  return 0;
}

int main(int argc, char* argv[])
{
  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // --------------------
  // Handle input options
  // --------------------

  // resize = 0 -- do not resize
  // resize = 1 -- call resize but with the same problem size
  // resize = 2 -- grow the problem by one element each time step
  int resize = 0;
  if (argc > 1) { resize = atoi(argv[1]); }
  if (resize > 2)
  {
    cerr << "Invalid resize value" << endl;
    return 1;
  }

  // Use Adams (1) or BDF methods (2)
  int method = CV_BDF;
  if (argc > 2)
  {
    if (atoi(argv[2]) == 1) { method = CV_ADAMS; }
    else if (atoi(argv[2]) == 2) { method = CV_BDF; }
    else
    {
      cerr << "Invalid method option" << endl;
      return 1;
    }
  }

  // Use fixed-point (0) or Newton (1) nonlinear solver
  int nonlinear_solver = 0;
  if (argc > 3)
  {
    nonlinear_solver = atoi(argv[3]);
    if (nonlinear_solver != 0 && nonlinear_solver != 1)
    {
      cerr << "Invalid nonlinear solver option" << endl;
      return 1;
    }
  }

  // Number of steps to advance
  int max_steps = 60;
  if (argc > 4) { max_steps = atoi(argv[4]); }

  cout << "CVODE Resize History Test\n";
  cout << "Method: ";
  if (method == 1) { cout << "Adams\n"; }
  else { cout << "BDF\n"; }
  cout << "Algebraic solvers: ";
  if (nonlinear_solver == 0) { cout << "Fixed-point\n"; }
  else { cout << "Newton + Band\n"; }
  cout << "Case: ";
  if (resize == 0) { cout << "No resize\n"; }
  else if (resize == 1)
  {
    cout << "Resize but do not change the problem size\n";
  }
  else { cout << "Resize each step with the problem size increased by one\n"; }

  // -------------
  // Setup problem
  // -------------

  // Set lambda
  sunrealtype lambda = -one;

  // Create initial condition
  N_Vector y = N_VNew_Serial(1, sunctx);
  if (check_ptr(y, "N_VNew_Serial")) { return 1; }

  int flag = true_solution(zero, y);
  if (check_flag(flag, "true_solution")) { return 1; }

  // Create CVODE
  void* cvode_mem = CVodeCreate(method, sunctx);
  if (check_ptr(cvode_mem, "CVodeCreate")) { return 1; }

  flag = CVodeInit(cvode_mem, ode_rhs, zero, y);
  if (check_flag(flag, "CVodeInit")) { return 1; }

  // Set relative and absolute tolerances
  const sunrealtype rtol = SUN_RCONST(1.0e-5);
  const sunrealtype atol = SUN_RCONST(1.0e-8);

  flag = CVodeSStolerances(cvode_mem, rtol, atol);
  if (check_flag(flag, "CVodeSStolerances")) { return 1; }

  // Attach user data
  flag = CVodeSetUserData(cvode_mem, &lambda);
  if (check_flag(flag, "CVodeSetUserData")) { return 1; }

  // Attach algebraic solvers
  SUNNonlinearSolver NLS = nullptr;
  SUNLinearSolver LS     = nullptr;
  SUNMatrix A            = nullptr;

  if (nonlinear_solver == 0)
  {
    // Use fixed-point nonlinear solver
    NLS = SUNNonlinSol_FixedPoint(y, 2, sunctx);
    if (check_ptr(NLS, "SUNNonlinSol_FixedPoint")) { return 1; }

    flag = CVodeSetNonlinearSolver(cvode_mem, NLS);
    if (check_flag(flag, "CVodeSetNonlinearSolver")) { return 1; }

    flag = CVodeSetMaxNonlinIters(cvode_mem, 10);
    if (check_flag(flag, "CVodeSetMaxNonlinIters")) { return 1; }
  }
  else
  {
    // Use default Newton solver, attach banded matrix and linear solver
    A = SUNBandMatrix(1, 0, 0, sunctx);
    if (check_ptr(A, "SUNBandMatrix")) { return 1; }

    LS = SUNLinSol_Band(y, A, sunctx);
    if (check_ptr(LS, "SUNLinSol_Band")) { return 1; }

    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (check_flag(flag, "CVodeSetLinearSolver")) { return 1; }

    flag = CVodeSetJacFn(cvode_mem, ode_rhs_jac);
    if (check_flag(flag, "CVodeSetJacFn")) { return 1; }
  }

  // ------------------------
  // Initialize saved history
  // ------------------------

  int hist_size = 0;
  if (method == CV_BDF) { hist_size = 5; }
  else { hist_size = 12; }

  sunrealtype* t_hist = new sunrealtype[hist_size];
  for (int i = 0; i < hist_size; i++) { t_hist[i] = zero; }

  N_Vector* y_hist = N_VCloneVectorArray(hist_size, y);
  if (check_ptr(y_hist, "N_VCloneVectorArray")) { return 1; }
  N_VScale(one, y, y_hist[0]);

  N_Vector* f_hist = N_VCloneVectorArray(hist_size, y);
  if (check_ptr(f_hist, "N_VCloneVectorArray")) { return 1; }
  flag = ode_rhs(t_hist[0], y_hist[0], f_hist[0], &lambda);
  if (check_flag(flag, "ode_rhs")) { return 1; }

  // Set output formatting
  cout << scientific;
  cout << setprecision(16);
  // cout << setprecision(numeric_limits<sunrealtype>::digits10);
  cout << endl;

  // Output initial data
  cout << "t:      " << t_hist[0] << endl;
  cout << "y:      " << N_VGetArrayPointer(y_hist[0])[0] << endl;
  cout << "y_true: " << N_VGetArrayPointer(y)[0] << endl;
  cout << "Error:  " << zero << endl;
  cout << "f:      " << N_VGetArrayPointer(f_hist[0])[0] << endl;

  // ---------------
  // Advance in time
  // ---------------

  // Final integration time
  sunrealtype tf = SUN_RCONST(10.0);

  // CVODE return time
  sunrealtype t_ret;

  for (int i = 1; i <= max_steps; i++)
  {
    N_Vector tmp = N_VClone(y);

    cout << flush;
    cerr << flush;
    cout << "\n========== Start Step " << i << " ==========\n\n";

    flag = CVode(cvode_mem, tf, y, &(t_ret), CV_ONE_STEP);
    if (check_flag(flag, "CVode"))
    {
      N_VDestroy(tmp);
      break;
    }

    int q_last;
    flag = CVodeGetLastOrder(cvode_mem, &q_last);

    sunrealtype h_last;
    flag = CVodeGetLastStep(cvode_mem, &h_last);

    cout << " Step Number: " << setw(3) << i << " | Time: " << setw(21) << t_ret
         << " | Step Size: " << setw(21) << h_last << " | Order: " << q_last
         << endl;

    PrintNordsieckArray(cvode_mem, false);
    if (check_flag(flag, "PrintNordsieckArray")) { return 1; }

    cout << "t:      " << t_ret << endl;
    cout << "y:      " << N_VGetArrayPointer(y)[0] << endl;

    flag = true_solution(t_ret, tmp);
    if (check_flag(flag, "true_solution")) { return 1; }
    cout << "y_true: " << N_VGetArrayPointer(tmp)[0] << endl;

    N_VLinearSum(one, y, -one, tmp, tmp);
    cout << "Error:  " << N_VMaxNorm(tmp) << endl;

    flag = ode_rhs(t_ret, y, tmp, &lambda);
    if (check_flag(flag, "ode_rhs")) { return 1; }
    cout << "f:      " << N_VGetArrayPointer(tmp)[0] << endl;

    cout << "========== End Step " << i << " ==========\n";

    // -------------------------------
    // Update Saved History and Resize
    // -------------------------------

    if (resize == 1)
    {
      // Save history and "resize" but do not change the problem size
      flag = save_history(t_ret, y, t_hist, y_hist, f_hist, hist_size, i,
                          &lambda);

      int n_hist = (i < hist_size) ? i + 1 : hist_size;

      flag = CVodeResizeHistory(cvode_mem, t_hist, y_hist, f_hist, n_hist,
                                n_hist);
      if (check_flag(flag, "CVodeResizeHistory")) { return 1; }
    }
    else if (resize == 2)
    {
      // Save history and increase the problem size
      flag = save_history(t_ret, y, t_hist, y_hist, f_hist, hist_size, i,
                          &lambda);

      flag = resize_history(t_hist, y_hist, f_hist, hist_size, sunctx, &lambda);

      int n_hist = (i < hist_size) ? i + 1 : hist_size;

      flag = CVodeResizeHistory(cvode_mem, t_hist, y_hist, f_hist, n_hist,
                                n_hist);
      if (check_flag(flag, "CVodeResizeHistory")) { return 1; }

      // "Resize" output vector and nonlinear solver
      N_VDestroy(y);
      y = N_VClone(y_hist[0]);

      if (nonlinear_solver == 0)
      {
        // Use fixed-point nonlinear solver
        SUNNonlinSolFree(NLS);

        NLS = SUNNonlinSol_FixedPoint(y, 2, sunctx);
        if (check_ptr(NLS, "SUNNonlinSol_FixedPoint")) { return 1; }

        flag = CVodeSetNonlinearSolver(cvode_mem, NLS);
        if (check_flag(flag, "CVodeSetNonlinearSolver")) { return 1; }

        flag = CVodeSetMaxNonlinIters(cvode_mem, 10);
        if (check_flag(flag, "CVodeSetMaxNonlinIters")) { return 1; }
      }
      else
      {
        SUNMatDestroy(A);
        SUNLinSolFree(LS);

        // Use default Newton solver, attach banded matrix and linear solver
        A = SUNBandMatrix(N_VGetLength(y), 0, 0, sunctx);
        if (check_ptr(A, "SUNBandMatrix")) { return 1; }

        LS = SUNLinSol_Band(y, A, sunctx);
        if (check_ptr(LS, "SUNLinSol_Band")) { return 1; }

        flag = CVodeSetLinearSolver(cvode_mem, LS, A);
        if (check_flag(flag, "CVodeSetLinearSolver")) { return 1; }

        flag = CVodeSetJacFn(cvode_mem, ode_rhs_jac);
        if (check_flag(flag, "CVodeSetJacFn")) { return 1; }
      }
    }

    N_VDestroy(tmp);
  }
  cout << endl;

  // Print final statistics
  flag = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(flag, "CVodePrintAllStats")) { return 1; }

  // Clean up and return with successful completion
  N_VDestroy(y);
  N_VDestroyVectorArray(y_hist, hist_size);
  N_VDestroyVectorArray(f_hist, hist_size);
  SUNNonlinSolFree(NLS);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  CVodeFree(&cvode_mem);
  delete[] t_hist;

  return 0;
}

/*---- end of file ----*/
