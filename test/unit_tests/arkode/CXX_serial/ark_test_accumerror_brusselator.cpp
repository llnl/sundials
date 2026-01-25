/*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *---------------------------------------------------------------
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
 *---------------------------------------------------------------
 * Routine to test the accumulated temporal error estimation
 * approaches from ARKStep.  Uses the
 * "stiff Brusselator" test problem with 3 components,
 *    du/dt = a - (w+1)*u + v*u^2
 *    dv/dt = w*u - v*u^2
 *    dw/dt = (b-w)/ep - w*u
 * for t in the interval [0.0, 10.0], with initial conditions
 * Y0 = [u0,v0,w0].
 *
 * The stiffness of the problem is essentially determined
 * by ep, wherein if the dynamical time step is given by H and
 * the explicit time step size is given by h, then H/h = 1/(100 ep),
 * i.e., the stability-limited step takes over at values of
 * ep < 1e-2.  This file defaults to a moderately stiff setup with
 * ep = 1/2500.
 *
 * We may run the problem in one of 3 different testing scenarios:
 *
 * Test 1:  u0=3.9,  v0=1.1,  w0=2.8,  a=1.2,  b=2.5
 *    Here, all three components exhibit a rapid transient change
 *    during the first 0.2 time units, followed by a slow and
 *    smooth evolution.
 *
 * Test 2 [default]:  u0=1.2,  v0=3.1,  w0=3,  a=1,  b=3.5
 *    Here, w experiences a fast initial transient, jumping 0.5
 *    within a few steps.  All values proceed smoothly until
 *    around t=6.5, when both u and v undergo a sharp transition,
 *    with u increaseing from around 0.5 to 5 and v decreasing
 *    from around 6 to 1 in less than 0.5 time units.  After this
 *    transition, both u and v continue to evolve somewhat
 *    rapidly for another 1.4 time units, and finish off smoothly.
 *
 * Test 3:  u0=3,  v0=3,  w0=3.5,  a=0.5,  b=3
 *    Here, all components undergo very rapid initial transients
 *    during the first 0.3 time units, and all then proceed very
 *    smoothly for the remainder of the simulation.
 *
 * We partition the full time integration interval, 0 < t < 5, into
 * Npart pieces, and run the accumulation test over each.
 *
 * We use either the ARKStep/DIRK/Newton/Dense solver (0) or
 * the ARKStep/ERK solver (1).  Either defaults to using a
 * 4th-order method.
 *
 * By default, all runs use temporal adaptivity; however, if the
 * requested 'ord' command-line input is negative, we run with
 * order |ord|, using fixed step sizes.
 *
 * The program should be run with arguments in the following order:
 *   $ a.out Npart ord method ep test
 * Not all arguments are required, but these must be omitted from
 * end-to-beginning, i.e. any one of
 *   $ a.out Npart ord method ep
 *   $ a.out Npart ord method
 *   $ a.out Npart ord
 *   $ a.out Npart
 *   $ a.out
 * are acceptable.  We require:
 *   * method = {0, 1}
 *   * test = {1, 2, 3}
 *   * ep > 0
 *   * Npart > 0
 *
 * For either temporally adaptive (ord >= 0) or fixed-step (ord < 0)
 * runs, we test a variety of tolerances/step sizes, and compare
 * the error at the end of each partition (computed via a reference
 * solution) against the integrator-reported accumulated error
 * estimate.
 *-----------------------------------------------------------------*/

// Header files
#include <arkode/arkode_arkstep.h>
#include <cmath>
#include <iostream>
#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <string.h>
#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <vector>

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)

using namespace std;

// User data structure
struct UserData
{
  sunrealtype a;
  sunrealtype b;
  sunrealtype ep;
  int Npart;
};

// User-supplied Functions Called by the Solver
static int fn(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Private utility functions
static int adaptive_run(void* arkode_mem, N_Vector y, sunrealtype T0,
                        sunrealtype Tf, int rk_type, int order, N_Vector* yref,
                        UserData& udata);
static int fixed_run(void* arkode_mem, N_Vector y, sunrealtype T0, sunrealtype Tf,
                     int rk_type, int order, N_Vector* yref, UserData& udata);
static int computeErrorWeights(N_Vector ycur, N_Vector weight, sunrealtype rtol,
                               sunrealtype atol, N_Vector vtemp);
static int check_retval(void* returnvalue, const char* funcname, int opt);

// Main Program
int main(int argc, char* argv[])
{
  // general problem parameters
  sunrealtype T0          = SUN_RCONST(0.0);  // initial time
  sunrealtype Tf          = SUN_RCONST(10.0); // final time
  sunindextype NEQ        = 3;                // number of dependent vars.
  int rk_type             = 1;       // type of RK method [DIRK=0, ERK=1]
  int order               = 4;       // order of accuracy for RK method
  int test                = 2;       // test problem to run
  sunbooleantype adaptive = SUNTRUE; // adaptive run vs convergence order
  sunrealtype u0, v0, w0;            // parameters

  // general problem variables
  int retval;              // reusable error-checking flag
  N_Vector y       = NULL; // empty vector for storing solution
  N_Vector* yref   = NULL; // empty vectors for storing reference solution
  void* arkode_mem = NULL; // empty ARKStep memory structure
  void* arkode_ref = NULL; // empty ARKStep memory structure for reference solution
  SUNMatrix A        = NULL; // empty matrix for solver
  SUNLinearSolver LS = NULL; // empty linear solver object
  UserData udata;            // user-data structure
  sunrealtype* ydata = NULL;
  udata.ep           = SUN_RCONST(0.0004); // stiffness parameter
  udata.Npart        = 20;                 // partition size

  //
  // Initialization
  //

  // Retrieve the command-line options:  Npart ord method ep test
  if (argc > 1) udata.Npart = atoi(argv[1]);
  if (argc > 2) order = atoi(argv[2]);
  if (argc > 3) rk_type = atoi(argv[3]);
  if (argc > 4) udata.ep = SUNStrToReal(argv[4]);
  if (argc > 5) test = atoi(argv[5]);

  // Check arguments for validity
  //    method = {0, 1}
  //    test = {1, 2, 3}
  //    ep > 0
  //    Npart > 0
  if ((rk_type < 0) || (rk_type > 1))
  {
    cerr << "ERROR: RK type be an integer in {0,1} \n";
    return (-1);
  }
  if ((test < 1) || (test > 3))
  {
    cerr << "ERROR: test type be an integer in {1,2,3} \n";
    return (-1);
  }
  if (udata.ep <= ZERO)
  {
    cerr << "ERROR: ep must be a positive real number\n";
    return (-1);
  }
  if (udata.Npart < 1)
  {
    cerr << "ERROR: Npart must be a positive integer\n";
    return (-1);
  }

  // Handle adaptive run vs order-of-convergence run
  if (order < 0)
  {
    adaptive = SUNFALSE;
    order    = abs(order);
  }
  if (order == 0) order = 4;

  // set up the test problem according to the desired test
  if (test == 1)
  {
    u0      = SUN_RCONST(3.9);
    v0      = SUN_RCONST(1.1);
    w0      = SUN_RCONST(2.8);
    udata.a = SUN_RCONST(1.2);
    udata.b = SUN_RCONST(2.5);
  }
  else if (test == 3)
  {
    u0      = SUN_RCONST(3.0);
    v0      = SUN_RCONST(3.0);
    w0      = SUN_RCONST(3.5);
    udata.a = SUN_RCONST(0.5);
    udata.b = SUN_RCONST(3.0);
  }
  else
  {
    u0      = SUN_RCONST(1.2);
    v0      = SUN_RCONST(3.1);
    w0      = SUN_RCONST(3.0);
    udata.a = SUN_RCONST(1.0);
    udata.b = SUN_RCONST(3.5);
  }

  // Initial problem output (and set implicit solver tolerances as needed)
  cout
    << "\nAccumulated error estimation test (stiff Brusselator ODE problem):\n";
  cout << "    time domain:  (" << T0 << "," << Tf << "]\n";
  cout << "    partition size = " << udata.Npart << endl;
  cout << "    initial conditions:  u0 = " << u0 << ",  v0 = " << v0
       << ",  w0 = " << w0 << endl;
  cout << "    problem parameters:  a = " << udata.a << ",  b = " << udata.b
       << ",  ep = " << udata.ep << endl;
  if (rk_type == 0) { cout << "    DIRK solver, order = " << order << endl; }
  else if (rk_type == 1)
  {
    cout << "    ERK solver, order = " << order << endl;
  }

  //
  // Problem Setup
  //

  // Create SUNDIALS context
  sundials::Context ctx;

  // Create serial vectors for the solution and reference
  y = N_VNew_Serial(NEQ, ctx);
  if (check_retval((void*)y, "N_VNew_Serial", 0)) return 1;
  yref = N_VCloneVectorArray(udata.Npart + 1, y);
  if (check_retval((void*)yref, "N_VNew_Serial", 0)) return 1;
  ydata = N_VGetArrayPointer(y);
  if (check_retval((void*)ydata, "N_VGetArrayPointer", 0)) return 1;

  // Generate reference solution
  ydata[0]   = u0;
  ydata[1]   = v0;
  ydata[2]   = w0;
  arkode_ref = ARKStepCreate(fn, NULL, T0, y, ctx);
  if (check_retval((void*)arkode_ref, "ARKStepCreate", 0)) return 1;
  retval = ARKodeSetUserData(arkode_ref, (void*)&udata);
  if (check_retval(&retval, "ARKodeSetUserData", 1)) return 1;
  retval = ARKodeSetOrder(arkode_ref, 5);
  if (check_retval(&retval, "ARKodeSetOrder", 1)) return 1;
  retval = ARKodeSStolerances(arkode_ref, SUN_RCONST(1.e-10), SUN_RCONST(1.e-12));
  if (check_retval(&retval, "ARKodeSStolerances", 1)) return 1;
  retval = ARKodeSetMaxNumSteps(arkode_ref, 1000000);
  if (check_retval(&retval, "ARKodeSetMaxNumSteps", 1)) return (1);
  N_VScale(ONE, y, yref[0]);
  sunrealtype hpart = (Tf - T0) / udata.Npart;
  for (int ipart = 0; ipart < udata.Npart; ipart++)
  {
    sunrealtype t = T0 + ipart * hpart;
    retval        = ARKodeSetStopTime(arkode_ref, t + hpart);
    if (check_retval(&retval, "ARKodeSetStopTime", 1)) return 1;
    retval = ARKodeEvolve(arkode_ref, t + hpart, y, &t, ARK_NORMAL);
    if (check_retval(&retval, "ARKodeEvolve", 1)) return 1;
    N_VScale(ONE, y, yref[ipart + 1]);
  }
  ARKodeFree(&arkode_ref);

  // Set up ARKStep integrator
  ydata[0] = u0;
  ydata[1] = v0;
  ydata[2] = w0;
  if (rk_type == 0)
  { // DIRK method
    arkode_mem = ARKStepCreate(NULL, fn, T0, y, ctx);
  }
  else
  { // ERK method
    arkode_mem = ARKStepCreate(fn, NULL, T0, y, ctx);
  }
  if (check_retval((void*)arkode_mem, "ARKStepCreate", 0)) return 1;
  retval = ARKodeSetUserData(arkode_mem, (void*)&udata);
  if (check_retval(&retval, "ARKodeSetUserData", 1)) return 1;
  retval = ARKodeSetOrder(arkode_mem, order);
  if (check_retval(&retval, "ARKodeSetOrder", 1)) return 1;
  retval = ARKodeSStolerances(arkode_mem, SUN_RCONST(1.e-4), SUN_RCONST(1.e-9));
  if (check_retval(&retval, "ARKodeSStolerances", 1)) return 1;
  if (rk_type == 0)
  { // DIRK method
    A = SUNDenseMatrix(NEQ, NEQ, ctx);
    if (check_retval((void*)A, "SUNDenseMatrix", 0)) return 1;
    LS = SUNLinSol_Dense(y, A, ctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return 1;
    retval = ARKodeSetLinearSolver(arkode_mem, LS, A);
    if (check_retval(&retval, "ARKodeSetLinearSolver", 1)) return 1;
    retval = ARKodeSetJacFn(arkode_mem, Jac);
    if (check_retval(&retval, "ARKodeSetJacFn", 1)) return 1;
  }
  else
  { // ERK method
    retval = ARKodeSetMaxNumSteps(arkode_mem, 1000000);
    if (check_retval(&retval, "ARKodeSetMaxNumSteps", 1)) return (1);
  }

  // Integrate ODE, based on run type
  if (adaptive)
  {
    retval = adaptive_run(arkode_mem, y, T0, Tf, rk_type, order, yref, udata);
    if (check_retval(&retval, "adaptive_run", 1)) return 1;
  }
  else
  {
    retval = fixed_run(arkode_mem, y, T0, Tf, rk_type, order, yref, udata);
    if (check_retval(&retval, "fixed_run", 1)) return 1;
  }

  // Clean up and return
  ARKodeFree(&arkode_mem);
  if (LS) { SUNLinSolFree(LS); } // free system linear solver
  if (A) { SUNMatDestroy(A); }   // free system matrix
  N_VDestroy(y);                 // Free y and yref vectors
  N_VDestroyVectorArray(yref, udata.Npart + 1);
  return 0;
}

//------------------------------
// Functions called by the solver
//------------------------------

static int fn(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  UserData* udata     = (UserData*)user_data;
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  sunrealtype* dydata = N_VGetArrayPointer(ydot);
  const sunrealtype u = ydata[0]; // access solution values
  const sunrealtype v = ydata[1];
  const sunrealtype w = ydata[2];

  // fill in the RHS function
  dydata[0] = udata->a - (w + ONE) * u + v * u * u;
  dydata[1] = w * u - v * u * u;
  dydata[2] = (udata->b - w) / udata->ep - w * u;

  // Return with success
  return 0;
}

static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData* udata     = (UserData*)user_data;
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  const sunrealtype u = ydata[0]; // access solution values
  const sunrealtype v = ydata[1];
  const sunrealtype w = ydata[2];

  // fill in the Jacobian
  SM_ELEMENT_D(J, 0, 0) = -(w + ONE) + TWO * u * v;
  SM_ELEMENT_D(J, 0, 1) = u * u;
  SM_ELEMENT_D(J, 0, 2) = -u;

  SM_ELEMENT_D(J, 1, 0) = w - TWO * u * v;
  SM_ELEMENT_D(J, 1, 1) = -u * u;
  SM_ELEMENT_D(J, 1, 2) = u;

  SM_ELEMENT_D(J, 2, 0) = -w;
  SM_ELEMENT_D(J, 2, 1) = ZERO;
  SM_ELEMENT_D(J, 2, 2) = -ONE / udata->ep - u;

  // Return with success
  return 0;
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

//------------------------------
// Private helper functions
//------------------------------

static int adaptive_run(void* arkode_mem, N_Vector y, sunrealtype T0,
                        sunrealtype Tf, int rk_type, int order, N_Vector* yref,
                        UserData& udata)
{
  // Reused variables
  int retval;
  sunrealtype t;
  sunrealtype hpart                 = (Tf - T0) / udata.Npart;
  sunrealtype abstol                = SUN_RCONST(1.e-12);
  vector<sunrealtype> rtols         = {SUN_RCONST(1.e-2), SUN_RCONST(1.e-4),
                                       SUN_RCONST(1.e-6)};
  vector<ARKAccumError> accum_types = {ARK_ACCUMERROR_MAX, ARK_ACCUMERROR_SUM,
                                       ARK_ACCUMERROR_AVG};
  vector<sunrealtype> dsm(udata.Npart);
  vector<sunrealtype> dsm_est(udata.Npart);
  vector<long int> Nsteps(udata.Npart);

  // Loop over tolerances
  cout << "\nAdaptive-step runs:\n";
  for (size_t irtol = 0; irtol < rtols.size(); irtol++)
  {
    // Loop over accumulation types
    for (size_t iaccum = 0; iaccum < accum_types.size(); iaccum++)
    {
      // Loop over partition
      for (int ipart = 0; ipart < udata.Npart; ipart++)
      {
        // Reset integrator for this run, and evolve over partition interval
        t = T0 + ipart * hpart;
        N_VScale(ONE, yref[ipart], y);
        if (rk_type == 0)
        { // DIRK
          retval = ARKStepReInit(arkode_mem, NULL, fn, t, y);
        }
        else
        { // ERK
          retval = ARKStepReInit(arkode_mem, fn, NULL, t, y);
        }
        if (check_retval(&retval, "ARKStepReInit", 1)) return 1;
        retval = ARKodeSetAccumulatedErrorType(arkode_mem, accum_types[iaccum]);
        if (check_retval(&retval, "ARKodeSetAccumulatedErrorType", 1)) return 1;
        retval = ARKodeResetAccumulatedError(arkode_mem);
        if (check_retval(&retval, "ARKodeResetAccumulatedError", 1)) return 1;
        retval = ARKodeSStolerances(arkode_mem, rtols[irtol], abstol);
        if (check_retval(&retval, "ARKodeSStolerances", 1)) return 1;
        retval = ARKodeSetStopTime(arkode_mem, t + hpart);
        if (check_retval(&retval, "ARKodeSetStopTime", 1)) return 1;
        retval = ARKodeSetMaxNumSteps(arkode_mem, 1000000);
        if (check_retval(&retval, "ARKodeSetMaxNumSteps", 1)) return (1);
        retval = ARKodeEvolve(arkode_mem, t + hpart, y, &t, ARK_NORMAL);
        if (check_retval(&retval, "ARKodeEvolve", 1)) break;
        retval = ARKodeGetAccumulatedError(arkode_mem, &(dsm_est[ipart]));
        if (check_retval(&retval, "ARKodeGetAccumulatedError", 1)) break;
        retval = ARKodeGetNumSteps(arkode_mem, &(Nsteps[ipart]));
        if (check_retval(&retval, "ARKodeGetNumSteps", 1)) break;

        // Compute/print solution error
        sunrealtype* ydata     = N_VGetArrayPointer(y);
        sunrealtype* yrefdata  = N_VGetArrayPointer(yref[ipart + 1]);
        const sunrealtype udsm = abs(ydata[0] - yrefdata[0]) /
                                 (abstol + rtols[irtol] * abs(yrefdata[0]));
        const sunrealtype vdsm = abs(ydata[1] - yrefdata[1]) /
                                 (abstol + rtols[irtol] * abs(yrefdata[1]));
        const sunrealtype wdsm = abs(ydata[2] - yrefdata[2]) /
                                 (abstol + rtols[irtol] * abs(yrefdata[2]));
        dsm[ipart] =
          rtols[irtol] *
          sqrt((udsm * udsm + vdsm * vdsm + wdsm * wdsm) / SUN_RCONST(3.0));
        cout << "  rtol " << rtols[irtol] << "  rk_type " << rk_type
             << "  order " << order << "  acc " << accum_types[iaccum] << "  t "
             << t << "  dsm " << dsm[ipart] << "  dsm_est " << dsm_est[ipart]
             << "  nsteps " << Nsteps[ipart] << endl;
      }
    }
  }

  return (0);
}

static int fixed_run(void* arkode_mem, N_Vector y, sunrealtype T0, sunrealtype Tf,
                     int rk_type, int order, N_Vector* yref, UserData& udata)
{
  // Reused variables
  int retval;
  sunrealtype hpart = (Tf - T0) / udata.Npart;
  long int nsteps2;
  sunrealtype t, t2;
  sunrealtype reltol = SUN_RCONST(1.e-9);
  sunrealtype abstol = SUN_RCONST(1.e-12);
  N_Vector y2        = N_VClone(y);
  N_Vector ewt       = N_VClone(y);
  N_Vector vtemp     = N_VClone(y);

  // Set array of fixed step sizes to use, storage for corresponding errors/orders
  sunrealtype hmax = (Tf - T0) / 400;
  if (rk_type == 1) hmax = min(hmax, udata.ep);
  vector<sunrealtype> hvals         = {hmax, hmax / 4, hmax / 16, hmax / 64};
  vector<ARKAccumError> accum_types = {ARK_ACCUMERROR_MAX, ARK_ACCUMERROR_SUM,
                                       ARK_ACCUMERROR_AVG};
  vector<sunrealtype> dsm(udata.Npart);
  vector<sunrealtype> dsm_est(udata.Npart);
  vector<long int> Nsteps(udata.Npart);

  // Loop over step sizes
  cout << "\nFixed-step runs:\n";
  for (size_t ih = 0; ih < hvals.size(); ih++)
  {
    // Loop over built-in accumulation types
    for (size_t iaccum = 0; iaccum < accum_types.size(); iaccum++)
    {
      // Loop over partition
      for (int ipart = 0; ipart < udata.Npart; ipart++)
      {
        // Reset integrator for this run, and evolve to Tf
        t = T0 + ipart * hpart;
        N_VScale(ONE, yref[ipart], y);
        if (rk_type == 0)
        { // DIRK
          retval = ARKStepReInit(arkode_mem, NULL, fn, t, y);
        }
        else
        { // ERK
          retval = ARKStepReInit(arkode_mem, fn, NULL, t, y);
        }
        if (check_retval(&retval, "ARKStepReInit", 1)) return 1;
        retval = ARKodeSetAccumulatedErrorType(arkode_mem, accum_types[iaccum]);
        if (check_retval(&retval, "ARKodeSetAccumulatedErrorType", 1)) return 1;
        retval = ARKodeResetAccumulatedError(arkode_mem);
        if (check_retval(&retval, "ARKodeResetAccumulatedError", 1)) return 1;
        retval = ARKodeSetFixedStep(arkode_mem, hvals[ih]);
        if (check_retval(&retval, "ARKodeSetFixedStep", 1)) return 1;
        retval = ARKodeSetMaxNumSteps(arkode_mem, 1000000);
        if (check_retval(&retval, "ARKodeSetMaxNumSteps", 1)) return (1);
        retval = ARKodeSetStopTime(arkode_mem, t + hpart);
        if (check_retval(&retval, "ARKodeSetStopTime", 1)) return 1;
        retval = ARKodeSStolerances(arkode_mem, reltol, abstol);
        if (check_retval(&retval, "ARKodeSStolerances", 1)) return 1;
        if (rk_type == 0)
        {
          retval = ARKodeSetJacEvalFrequency(arkode_mem, 1);
          if (check_retval(&retval, "ARKodeSetJacEvalFrequency", 1)) return 1;
          retval = ARKodeSetLSetupFrequency(arkode_mem, 1);
          if (check_retval(&retval, "ARKodeSetLSetupFrequency", 1)) return 1;
          retval = ARKodeSetMaxNonlinIters(arkode_mem, 20);
          if (check_retval(&retval, "ARKodeSetMaxNonlinIters", 1)) return 1;
        }
        retval = ARKodeEvolve(arkode_mem, t + hpart, y, &t, ARK_NORMAL);
        if (check_retval(&retval, "ARKodeEvolve", 1)) break;
        retval = ARKodeGetAccumulatedError(arkode_mem, &(dsm_est[ipart]));
        if (check_retval(&retval, "ARKodeGetAccumulatedError", 1)) break;
        retval = ARKodeGetNumSteps(arkode_mem, &(Nsteps[ipart]));
        if (check_retval(&retval, "ARKodeGetNumSteps", 1)) break;

        // Compute/print solution error
        sunrealtype* ydata     = N_VGetArrayPointer(y);
        sunrealtype* yrefdata  = N_VGetArrayPointer(yref[ipart + 1]);
        const sunrealtype udsm = abs(ydata[0] - yrefdata[0]) /
                                 (abstol + reltol * abs(yrefdata[0]));
        const sunrealtype vdsm = abs(ydata[1] - yrefdata[1]) /
                                 (abstol + reltol * abs(yrefdata[1]));
        const sunrealtype wdsm = abs(ydata[2] - yrefdata[2]) /
                                 (abstol + reltol * abs(yrefdata[2]));
        dsm[ipart] = reltol * sqrt((udsm * udsm + vdsm * vdsm + wdsm * wdsm) /
                                   SUN_RCONST(3.0));
        cout << "  h " << hvals[ih] << "  rk_type " << rk_type << "  order "
             << order << "  acc " << accum_types[iaccum] << "  t " << t
             << "  dsm " << dsm[ipart] << "  dsm_est " << dsm_est[ipart]
             << "  nsteps " << Nsteps[ipart] << endl;
      }
    }

    // Test double-step error estimator

    // Loop over partition
    for (int ipart = 0; ipart < udata.Npart; ipart++)
    {
      // Reset integrator for this run, and evolve over partition interval
      t = t2 = T0 + ipart * hpart;
      N_VScale(ONE, yref[ipart], y);
      N_VScale(ONE, yref[ipart], y2);
      if (rk_type == 0)
      { // DIRK
        retval = ARKStepReInit(arkode_mem, NULL, fn, t, y);
      }
      else
      { // ERK
        retval = ARKStepReInit(arkode_mem, fn, NULL, t, y);
      }
      if (check_retval(&retval, "ARKStepReInit", 1)) return 1;
      retval = ARKodeSetAccumulatedErrorType(arkode_mem, ARK_ACCUMERROR_NONE);
      if (check_retval(&retval, "ARKodeSetAccumulatedErrorType", 1)) return 1;
      retval = ARKodeSetFixedStep(arkode_mem, hvals[ih]);
      if (check_retval(&retval, "ARKodeSetFixedStep", 1)) return 1;
      retval = ARKodeSetMaxNumSteps(arkode_mem, 1000000);
      if (check_retval(&retval, "ARKodeSetMaxNumSteps", 1)) return (1);
      retval = ARKodeSetStopTime(arkode_mem, t + hpart);
      if (check_retval(&retval, "ARKodeSetStopTime", 1)) return 1;
      retval = ARKodeSStolerances(arkode_mem, reltol, abstol);
      if (check_retval(&retval, "ARKodeSStolerances", 1)) return 1;
      if (rk_type == 0)
      {
        retval = ARKodeSetJacEvalFrequency(arkode_mem, 1);
        if (check_retval(&retval, "ARKodeSetJacEvalFrequency", 1)) return 1;
        retval = ARKodeSetLSetupFrequency(arkode_mem, 1);
        if (check_retval(&retval, "ARKodeSetLSetupFrequency", 1)) return 1;
        retval = ARKodeSetMaxNonlinIters(arkode_mem, 20);
        if (check_retval(&retval, "ARKodeSetMaxNonlinIters", 1)) return 1;
      }
      retval = ARKodeEvolve(arkode_mem, t + hpart, y, &t, ARK_NORMAL);
      if (check_retval(&retval, "ARKodeEvolve", 1)) break;
      retval = ARKodeGetNumSteps(arkode_mem, &(Nsteps[ipart]));
      if (check_retval(&retval, "ARKodeGetNumSteps", 1)) break;

      if (rk_type == 0)
      { // DIRK
        retval = ARKStepReInit(arkode_mem, NULL, fn, t2, y2);
      }
      else
      { // ERK
        retval = ARKStepReInit(arkode_mem, fn, NULL, t2, y2);
      }
      if (check_retval(&retval, "ARKStepReInit", 1)) return 1;
      retval = ARKodeSetFixedStep(arkode_mem, 2.0 * hvals[ih]);
      if (check_retval(&retval, "ARKodeSetFixedStep", 1)) return 1;
      retval = ARKodeSetStopTime(arkode_mem, t2 + hpart);
      if (check_retval(&retval, "ARKodeSetStopTime", 1)) return 1;
      retval = ARKodeEvolve(arkode_mem, t2 + hpart, y2, &t2, ARK_NORMAL);
      if (check_retval(&retval, "ARKodeEvolve", 1)) break;
      retval = ARKodeGetNumSteps(arkode_mem, &nsteps2);
      if (check_retval(&retval, "ARKodeGetNumSteps", 1)) break;

      retval = computeErrorWeights(y2, ewt, reltol, abstol, vtemp);
      if (check_retval(&retval, "computeErrorWeights", 1)) break;
      N_VLinearSum(ONE, y2, -ONE, y, y2);
      dsm_est[ipart] = reltol * N_VWrmsNorm(y2, ewt);
      Nsteps[ipart] += nsteps2;

      sunrealtype* ydata     = N_VGetArrayPointer(y);
      sunrealtype* yrefdata  = N_VGetArrayPointer(yref[ipart + 1]);
      const sunrealtype udsm = abs(ydata[0] - yrefdata[0]) /
                               (abstol + reltol * abs(yrefdata[0]));
      const sunrealtype vdsm = abs(ydata[1] - yrefdata[1]) /
                               (abstol + reltol * abs(yrefdata[1]));
      const sunrealtype wdsm = abs(ydata[2] - yrefdata[2]) /
                               (abstol + reltol * abs(yrefdata[2]));
      dsm[ipart] = reltol * sqrt((udsm * udsm + vdsm * vdsm + wdsm * wdsm) /
                                 SUN_RCONST(3.0));
      cout << "  h " << hvals[ih] << "  rk_type " << rk_type << "  order "
           << order << "  acc " << 2 << "  t " << t << "  dsm " << dsm[ipart]
           << "  dsm_est " << dsm_est[ipart] << "  nsteps " << Nsteps[ipart]
           << endl;
    }
  }

  N_VDestroy(y2);
  N_VDestroy(ewt);
  N_VDestroy(vtemp);
  return (0);
}

/* Error weight calculation routine (mimics what's in ARKODE already) */
static int computeErrorWeights(N_Vector ycur, N_Vector weight, sunrealtype rtol,
                               sunrealtype atol, N_Vector vtemp)
{
  N_VAbs(ycur, vtemp);
  N_VScale(rtol, vtemp, vtemp);
  N_VAddConst(vtemp, atol, vtemp);
  N_VInv(vtemp, weight);
  return (0);
}

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a retval so check if
             retval >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* retval;

  // Check if SUNDIALS function returned NULL pointer - no memory allocated
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  // Check if retval < 0
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return 1;
    }
  }

  // Check if function returned NULL pointer - no memory allocated
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}

/*---- end of file ----*/
