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
 * Routine to test an MRI method's embedding-based error
 * estimate.  Uses the "stiff Brusselator" test problem with 3
 * components,
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
 * Npart pieces.  We then run a single time step starting at the
 * beginning of each partition, using a variety of slow step sizes,
 * H = {hmax, hmax/4, hmax/16, hmax/64, hmax/256} with
 * hmax=(t_f-t_0)/20/Npart.
 *
 * We place the entire ODE in the "slow" RHS partition.  For IMEX
 * methods, the third row is treated implicitly, and the first two
 * are treated explicitly.  For the fast time scale, all tests use
 * ARKODE's default fifth-order ERK method, with relative and
 * absolute tolerances set to 1e-10 and 1e-12, respectively.
 *
 * We select the slow integrator based on a command-line argument,
 * with the default being ARKODE_MRI_GARK_ERK33a.
 *
 * The program should be run with arguments in the following order:
 *   $ a.out method Npart ep test
 * Not all arguments are required, but these must be omitted from
 * end-to-beginning, i.e. any one of
 *   $ a.out method Npart ep
 *   $ a.out method Npart
 *   $ a.out method
 *   $ a.out
 * are acceptable.  We require:
 *   * method = string corresponding to a valid embedded ARKODE_MRITableID
 *   * Npart > 0
 *   * ep > 0
 *   * test = {1, 2, 3}
 *-----------------------------------------------------------------*/

// Header files
#include <algorithm>
#include <arkode/arkode_erkstep.h>
#include <arkode/arkode_mristep.h>
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
static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fn(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fe(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Jn(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Ji(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Private utility functions
static int computeErrorWeights(N_Vector ycur, N_Vector weight, sunrealtype rtol,
                               sunrealtype atol, N_Vector vtemp);
static int check_retval(void* returnvalue, const char* funcname, int opt);
static int run_test(void* mristep_mem, void* arkode_ref, N_Vector y,
                    sunrealtype T0, sunrealtype Tf, N_Vector* yref,
                    vector<sunrealtype> Hvals, string method,
                    sunrealtype reltol, sunrealtype abstol, UserData& udata);

// Main Program
int main(int argc, char* argv[])
{
  // general problem parameters
  sunrealtype T0   = SUN_RCONST(0.0);      // initial time
  sunrealtype Tf   = SUN_RCONST(10.0);     // final time
  sunindextype NEQ = 3;                    // number of dependent vars.
  string method;                           // MRI method name
  int test = 2;                            // test problem to run
  sunrealtype u0, v0, w0;                  // parameters
  sunrealtype reltol = SUN_RCONST(1.e-10); // fast solver tolerances
  sunrealtype abstol = SUN_RCONST(1.e-12);

  // general problem variables
  int retval;                       // reusable error-checking flag
  UserData udata;                   // user-data structure
  udata.ep    = SUN_RCONST(0.0004); // stiffness parameter
  udata.Npart = 20;                 // partition size

  //
  // Initialization
  //

  // Retrieve the command-line options:  method Npart ep test
  if (argc > 1) { method = argv[1]; }
  else { method = "ARKODE_MRI_GARK_ERK33a"; }
  if (argc > 2) udata.Npart = atoi(argv[2]);
  if (argc > 3) udata.ep = SUNStrToReal(argv[3]);
  if (argc > 4) test = atoi(argv[4]);

  // Check arguments for validity
  //    Npart > 0
  //    ep > 0
  //    test = {1, 2, 3}
  if (udata.Npart < 1)
  {
    cerr << "ERROR: Npart must be a positive integer\n";
    return (-1);
  }
  if (udata.ep <= ZERO)
  {
    cerr << "ERROR: ep must be a positive real number\n";
    return (-1);
  }
  if ((test < 1) || (test > 3))
  {
    cerr << "ERROR: test type be an integer in {1,2,3} \n";
    return (-1);
  }

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

  sunbooleantype implicit = SUNFALSE;
  sunbooleantype imex     = SUNFALSE;
  if ((method == "ARKODE_MRI_GARK_IRK21a") ||
      (method == "ARKODE_MRI_GARK_ESDIRK34a") ||
      (method == "ARKODE_MRI_GARK_ESDIRK46a"))
  {
    implicit = SUNTRUE;
  }
  if ((method == "ARKODE_IMEX_MRI_SR21") ||
      (method == "ARKODE_IMEX_MRI_SR32") || (method == "ARKODE_IMEX_MRI_SR43"))
  {
    imex     = SUNTRUE;
    implicit = SUNTRUE;
  }

  // Initial problem output (and set implicit solver tolerances as needed)
  cout << "\nSlow error estimation test (stiff Brusselator ODE problem):\n";
  cout << "    time domain:  (" << T0 << "," << Tf << "]\n";
  cout << "    partition size = " << udata.Npart << endl;
  cout << "    initial conditions:  u0 = " << u0 << ",  v0 = " << v0
       << ",  w0 = " << w0 << endl;
  cout << "    problem parameters:  a = " << udata.a << ",  b = " << udata.b
       << ",  ep = " << udata.ep << endl;
  cout << "    MRI method: " << method;
  if (imex) { cout << " (ImEx)" << endl; }
  else if (implicit) { cout << " (implicit)" << endl; }
  else { cout << " (explicit)" << endl; }

  //
  // Problem Setup
  //

  // Create SUNDIALS context
  sundials::Context ctx;

  // Create serial vectors for the solution and reference
  N_Vector y = N_VNew_Serial(NEQ, ctx);
  if (check_retval((void*)y, "N_VNew_Serial", 0)) return 1;
  N_Vector* yref = N_VCloneVectorArray(udata.Npart + 1, y);
  if (check_retval((void*)yref, "N_VNew_Serial", 0)) return 1;
  sunrealtype* ydata = N_VGetArrayPointer(y);
  if (check_retval((void*)ydata, "N_VGetArrayPointer", 0)) return 1;
  ydata[0] = u0;
  ydata[1] = v0;
  ydata[2] = w0;

  // Generate reference solution
  void* arkode_ref = ERKStepCreate(fn, T0, y, ctx);
  if (check_retval((void*)arkode_ref, "ERKStepCreate", 0)) return 1;
  retval = ARKodeSetUserData(arkode_ref, (void*)&udata);
  if (check_retval(&retval, "ARKodeSetUserData", 1)) return 1;
  retval = ARKodeSetOrder(arkode_ref, 5);
  if (check_retval(&retval, "ARKodeSetOrder", 1)) return 1;
  retval = ARKodeSStolerances(arkode_ref, reltol, abstol);
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

  // Set up fast ERKStep integrator as fifth-order adaptive method
  ydata[0]               = u0;
  ydata[1]               = v0;
  ydata[2]               = w0;
  void* inner_arkode_mem = ERKStepCreate(f0, T0, y, ctx);
  if (check_retval((void*)inner_arkode_mem, "ERKStepCreate", 0)) return 1;
  retval = ARKodeSetOrder(inner_arkode_mem, 5);
  if (check_retval(&retval, "ARKodeSetOrder", 1)) return 1;
  retval = ARKodeSStolerances(inner_arkode_mem, reltol, abstol);
  if (check_retval(&retval, "ARKodeSStolerances", 1)) return 1;
  retval = ARKodeSetMaxNumSteps(inner_arkode_mem, 1000000);
  if (check_retval(&retval, "ARKodeSetMaxNumSteps", 1)) return (1);

  // Create inner stepper wrapper
  MRIStepInnerStepper inner_stepper = NULL; // inner stepper
  retval = ARKodeCreateMRIStepInnerStepper(inner_arkode_mem, &inner_stepper);
  if (check_retval(&retval, "ARKodeCreateMRIStepInnerStepper", 1)) return 1;

  // Set up slow MRIStep integrator
  void* mristep_mem = NULL;
  if (imex) { mristep_mem = MRIStepCreate(fe, fi, T0, y, inner_stepper, ctx); }
  else if (implicit)
  {
    mristep_mem = MRIStepCreate(NULL, fn, T0, y, inner_stepper, ctx);
  }
  else { mristep_mem = MRIStepCreate(fn, NULL, T0, y, inner_stepper, ctx); }
  if (check_retval((void*)mristep_mem, "MRIStepCreate", 0)) return 1;
  MRIStepCoupling C = MRIStepCoupling_LoadTableByName(method.c_str());
  if (check_retval((void*)C, "MRIStepCoupling_LoadTableByName", 0)) return 1;
  retval = MRIStepSetCoupling(mristep_mem, C);
  if (check_retval(&retval, "MRIStepSetCoupling", 1)) return 1;
  SUNMatrix A        = NULL; // matrix for slow solver
  SUNLinearSolver LS = NULL; // slow linear solver object
  if (implicit)
  {
    A = SUNDenseMatrix(NEQ, NEQ, ctx);
    if (check_retval((void*)A, "SUNDenseMatrix", 0)) return 1;
    LS = SUNLinSol_Dense(y, A, ctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return 1;
    retval = ARKodeSetLinearSolver(mristep_mem, LS, A);
    if (check_retval(&retval, "ARKodeSetLinearSolver", 1)) return 1;
    if (imex) { retval = ARKodeSetJacFn(mristep_mem, Ji); }
    else { retval = ARKodeSetJacFn(mristep_mem, Jn); }
    if (check_retval(&retval, "ARKodeSetJacFn", 1)) return 1;
    retval = ARKodeSetJacEvalFrequency(mristep_mem, 1);
    if (check_retval(&retval, "ARKodeSetJacEvalFrequency", 1)) return 1;
    retval = ARKodeSetLSetupFrequency(mristep_mem, 1);
    if (check_retval(&retval, "ARKodeSetLSetupFrequency", 1)) return 1;
    retval = ARKodeSetMaxNonlinIters(mristep_mem, 50);
    if (check_retval(&retval, "ARKodeSetMaxNonlinIters", 1)) return 1;
  }
  retval = ARKodeSStolerances(mristep_mem, reltol, abstol);
  if (check_retval(&retval, "ARKodeSStolerances", 1)) return 1;
  retval = ARKodeSetUserData(mristep_mem, (void*)&udata);
  if (check_retval(&retval, "ARKodeSetUserData", 1)) return 1;
  retval = ARKodeSetAccumulatedErrorType(mristep_mem, ARK_ACCUMERROR_MAX);
  if (check_retval(&retval, "ARKodeSetAccumulatedErrorType", 1)) return 1;

  // Run test for various H values
  sunrealtype hmax = (Tf - T0) / 20.0 / udata.Npart;
  vector<sunrealtype> Hvals(5);
  for (size_t i = 0; i < Hvals.size(); i++)
  {
    Hvals[i] = hmax / SUNRpowerI(SUN_RCONST(4.0), (int)i);
  }
  retval = run_test(mristep_mem, arkode_ref, y, T0, Tf, yref, Hvals, method,
                    reltol, abstol, udata);
  if (check_retval(&retval, "run_test", 1)) return 1;

  // Clean up and return
  MRIStepCoupling_Free(C);
  ARKodeFree(&arkode_ref);
  ARKodeFree(&inner_arkode_mem);
  MRIStepInnerStepper_Free(&inner_stepper);
  ARKodeFree(&mristep_mem);
  if (LS) { SUNLinSolFree(LS); } // free system linear solver
  if (A) { SUNMatDestroy(A); }   // free system matrix
  N_VDestroy(y);                 // Free y and yref vectors
  N_VDestroyVectorArray(yref, udata.Npart + 1);
  return 0;
}

//------------------------------
// Functions called by the solver
//------------------------------

static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  // fill in the RHS function with zeros and return with success
  N_VConst(ZERO, ydot);
  return 0;
}

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

static int fi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  UserData* udata     = (UserData*)user_data;
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  sunrealtype* dydata = N_VGetArrayPointer(ydot);
  const sunrealtype u = ydata[0]; // access solution values
  const sunrealtype w = ydata[2];

  // fill in the RHS function
  dydata[0] = ZERO;
  dydata[1] = ZERO;
  dydata[2] = (udata->b - w) / udata->ep - w * u;

  // Return with success
  return 0;
}

static int fe(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
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
  dydata[2] = ZERO;

  // Return with success
  return 0;
}

static int Jn(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
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

static int Ji(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData* udata     = (UserData*)user_data;
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  const sunrealtype u = ydata[0]; // access solution values
  const sunrealtype w = ydata[2];

  // fill in the Jacobian
  SM_ELEMENT_D(J, 0, 0) = ZERO;
  SM_ELEMENT_D(J, 0, 1) = ZERO;
  SM_ELEMENT_D(J, 0, 2) = ZERO;

  SM_ELEMENT_D(J, 1, 0) = ZERO;
  SM_ELEMENT_D(J, 1, 1) = ZERO;
  SM_ELEMENT_D(J, 1, 2) = ZERO;

  SM_ELEMENT_D(J, 2, 0) = -w;
  SM_ELEMENT_D(J, 2, 1) = ZERO;
  SM_ELEMENT_D(J, 2, 2) = -ONE / udata->ep - u;

  // Return with success
  return 0;
}

//------------------------------
// Private helper functions
//------------------------------

static int run_test(void* mristep_mem, void* arkode_ref, N_Vector y,
                    sunrealtype T0, sunrealtype Tf, N_Vector* yref,
                    vector<sunrealtype> Hvals, string method,
                    sunrealtype reltol, sunrealtype abstol, UserData& udata)
{
  // Reused variables
  int retval;
  sunrealtype hpart = (Tf - T0) / udata.Npart;
  sunrealtype t, t2;
  N_Vector y2         = N_VClone(y);
  N_Vector ele        = N_VClone(y);
  N_Vector ewt        = N_VClone(y);
  N_Vector vtemp      = N_VClone(y);
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  sunrealtype* y2data = N_VGetArrayPointer(y2);

  // Set storage for errors
  vector<vector<sunrealtype>> dsm(Hvals.size(),
                                  vector<sunrealtype>(udata.Npart, ZERO));
  vector<vector<sunrealtype>> dsm_est(Hvals.size(),
                                      vector<sunrealtype>(udata.Npart, ZERO));

  // Loop over step sizes
  for (size_t iH = 0; iH < Hvals.size(); iH++)
  {
    // Loop over partition
    for (int ipart = 0; ipart < udata.Npart; ipart++)
    {
      // Reset integrators for this run
      t = t2 = T0 + ipart * hpart;
      N_VScale(ONE, yref[ipart], y);
      retval = ARKodeReset(mristep_mem, t, y);
      if (check_retval(&retval, "ARKodeReset", 1)) return 1;
      retval = ARKodeSetFixedStep(mristep_mem, Hvals[iH]);
      if (check_retval(&retval, "ARKodeSetFixedStep", 1)) return 1;
      retval = ARKodeResetAccumulatedError(mristep_mem);
      if (check_retval(&retval, "ARKodeResetAccumulatedError", 1)) return 1;
      N_VScale(ONE, yref[ipart], y2);
      retval = ARKodeReset(arkode_ref, t2, y2);
      if (check_retval(&retval, "ARKodeReset", 1)) return 1;
      retval = ARKodeSetStopTime(arkode_ref, t2 + Hvals[iH]);
      if (check_retval(&retval, "ARKodeSetStopTime", 1)) return 1;

      // Run ERKStep to compute reference solution, and MRIStep to compute one step
      retval = ARKodeEvolve(arkode_ref, t2 + Hvals[iH], y2, &t2, ARK_NORMAL);
      if (check_retval(&retval, "ARKodeEvolve", 1)) return 1;
      retval = ARKodeEvolve(mristep_mem, t + Hvals[iH], y, &t, ARK_ONE_STEP);
      if (check_retval(&retval, "ARKodeEvolve", 1)) return 1;
      retval = ARKodeGetEstLocalErrors(mristep_mem, ele);
      if (check_retval(&retval, "ARKodeGetEstLocalErrors", 1)) return 1;
      retval = computeErrorWeights(y, ewt, reltol, abstol, vtemp);
      if (check_retval(&retval, "computeErrorWeights", 1)) return 1;
      dsm_est[iH][ipart] = N_VWrmsNorm(ewt, ele);

      // Compute/print solution error
      sunrealtype udsm = abs(ydata[0] - y2data[0]) /
                         (abstol + reltol * abs(y2data[0]));
      sunrealtype vdsm = abs(ydata[1] - y2data[1]) /
                         (abstol + reltol * abs(y2data[1]));
      sunrealtype wdsm = abs(ydata[2] - y2data[2]) /
                         (abstol + reltol * abs(y2data[2]));
      dsm[iH][ipart] =
        sqrt((udsm * udsm + vdsm * vdsm + wdsm * wdsm) / SUN_RCONST(3.0));
      cout << "  H " << Hvals[iH] << "  method " << method << "  t " << t
           << "  dsm " << dsm[iH][ipart] << "  dsm_est " << dsm_est[iH][ipart]
           << endl;
    }
  }

  cout << endl << method << " summary:" << endl;
  for (size_t iH = 0; iH < Hvals.size(); iH++)
  {
    cout << "  Stepsize " << Hvals[iH] << "  \tmaxdsm "
         << *max_element(dsm[iH].begin(), dsm[iH].end()) << "  \tmaxdsmest "
         << *max_element(dsm_est[iH].begin(), dsm_est[iH].end()) << endl;
  }

  N_VDestroy(ele);
  N_VDestroy(ewt);
  N_VDestroy(vtemp);
  N_VDestroy(y2);
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
