/* ----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 * ----------------------------------------------------------------
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
 * ----------------------------------------------------------------
 * Routine to test an MRI method's embedding-based error
 * estimate.  Uses a nonlinear Kvaerno-Prothero-Robinson ODE test
 * problem with analytical solution,
 *
 *    [u]' =  [ G  e ] [(u^2-p-1)/(2u)] +  [ p'(t)/(2u) ]
 *    [v]     [ e -1 ] [(v^2-q-2)/(2v)]    [ q'(t)/(2v) ]
 *
 * where p(t) = cos(t), and q(t) = cos(omega*t*(1+exp(-(t-2)^2))).
 * This problem has analytical solution given by
 *    u(t) = sqrt(2+p(t)),  v(t) = sqrt(2+q(t)).
 * We use the parameters: e = 0.1 and G = -10 [default]
 *
 * The stiffness of the problem is essentially determined
 * by G, for |G| > 50 it is "stiff" and ideally suited to an
 * implicit method.
 *
 * Coupling between the two components is determined by e, with
 * coupling strength proportional to |e|.
 *
 * The "fast" variable, v, oscillates at a frequency "omega" times
 * faster than u.
 *
 * We partition the full time integration interval, 0 < t < 5, into
 * Npart pieces.  We then run a single time step starting at the
 * beginning of each partition, using a variety of slow step sizes,
 * H = {hmax, hmax/4, hmax/16, hmax/64, hmax/256} with
 * hmax=(t_f-t_0)/20/Npart.
 *
 * We place the entire ODE in the "slow" RHS partition.  For IMEX
 * methods, the first row is treated implicitly, and the second is
 * treated explicitly.  For the fast time scale, all tests use
 * ARKODE's default fifth-order ERK method, with relative and
 * absolute tolerances set to 1e-10 and 1e-12, respectively.
 *
 * We select the slow integrator based on a command-line argument,
 * with the default being ARKODE_MRI_GARK_ERK33a.
 *
 * The program should be run with arguments in the following order:
 *   $ a.out method Npart G e omega
 * Not all arguments are required, but these must be omitted from
 * end-to-beginning, i.e. any one of
 *   $ a.out method Npart G e
 *   $ a.out method Npart G
 *   $ a.out method Npart
 *   $ a.out method
 *   $ a.out
 * are acceptable.  We require:
 *   * method = string corresponding to a valid embedded ARKODE_MRITableID
 *   * Npart > 0
 *   * G < 0.0
 *   * omega > 0.0
 * ----------------------------------------------------------------*/

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
  sunrealtype G;
  sunrealtype e;
  sunrealtype omega;
  int Npart;
};

// User-supplied functions called by the solver
static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fn(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fe(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Jn(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Ji(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Private utility functions
static sunrealtype p(sunrealtype t);
static sunrealtype q(sunrealtype t, UserData& udata);
static sunrealtype pdot(sunrealtype t);
static sunrealtype qdot(sunrealtype t, UserData& udata);
static sunrealtype utrue(sunrealtype t);
static sunrealtype vtrue(sunrealtype t, UserData& udata);
static int Ytrue(sunrealtype t, N_Vector y, UserData& udata);
static int computeErrorWeights(N_Vector ycur, N_Vector weight, sunrealtype rtol,
                               sunrealtype atol, N_Vector vtemp);
static int check_retval(void* returnvalue, const char* funcname, int opt);
static int run_test(void* mristep_mem, N_Vector y, sunrealtype T0,
                    sunrealtype Tf, vector<sunrealtype> Hvals, string method,
                    sunrealtype reltol, sunrealtype abstol, UserData& udata);

// Main Program
int main(int argc, char* argv[])
{
  // general problem parameters
  sunrealtype T0   = SUN_RCONST(0.0);      // initial time
  sunrealtype Tf   = SUN_RCONST(5.0);      // final time
  sunindextype NEQ = 2;                    // number of dependent vars.
  string method;                           // MRI method name
  sunrealtype reltol = SUN_RCONST(1.e-10); // fast solver tolerances
  sunrealtype abstol = SUN_RCONST(1.e-12);

  // general problem variables
  int retval;                      // reusable error-checking flag
  UserData udata;                  // user-data structure
  udata.G     = SUN_RCONST(-10.0); // stiffness parameter
  udata.e     = SUN_RCONST(0.1);   // coupling strength
  udata.omega = SUN_RCONST(5.0);   // time scale ratio
  udata.Npart = 20;                // partition size

  //
  // Initialization
  //

  // Retrieve the command-line options:  method Npart G e omega
  if (argc > 1) { method = argv[1]; }
  else { method = "ARKODE_MRI_GARK_ERK33a"; }
  if (argc > 2) udata.Npart = atoi(argv[2]);
  if (argc > 3) udata.G = SUNStrToReal(argv[3]);
  if (argc > 4) udata.e = SUNStrToReal(argv[4]);
  if (argc > 5) udata.omega = SUNStrToReal(argv[5]);

  // Check arguments for validity
  //   G < 0.0
  //   omega > 0.0
  //   Npart > 0
  if (udata.G >= ZERO)
  {
    cerr << "ERROR: G must be a negative real number\n";
    return (-1);
  }
  if (udata.omega <= ZERO)
  {
    cerr << "ERROR: omega must be a positive real number\n";
    return (-1);
  }
  if (udata.Npart < 1)
  {
    cerr << "ERROR: Npart must be a positive integer\n";
    return (-1);
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
  cout << "\nSlow error estimation test (Nonlinear Kvaerno-Prothero-Robinson "
          "problem):\n";
  cout << "    time domain:  (" << T0 << "," << Tf << "]\n";
  cout << "    partition size = " << udata.Npart << endl;
  cout << "    problem parameters:  G = " << udata.G << ",  e = " << udata.e
       << ",  omega = " << udata.omega << endl;
  cout << "    MRI method: " << method;
  if (imex) { cout << " (ImEx)" << endl; }
  else if (implicit) { cout << " (implicit)" << endl; }
  else { cout << " (explicit)" << endl; }

  //
  // Problem Setup
  //

  // Create SUNDIALS context
  sundials::Context ctx;

  // Create and initialize serial vector for the solution
  N_Vector y = N_VNew_Serial(NEQ, ctx);
  if (check_retval((void*)y, "N_VNew_Serial", 0)) return 1;
  retval = Ytrue(T0, y, udata);
  if (check_retval(&retval, "Ytrue", 1)) return 1;

  // Set up fast ERKStep integrator as fifth-order adaptive method
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
  sunrealtype hmax = (Tf - T0) / 20 / udata.Npart;
  vector<sunrealtype> Hvals(5);
  for (size_t i = 0; i < Hvals.size(); i++)
  {
    Hvals[i] = hmax / SUNRpowerI(SUN_RCONST(4.0), (int)i);
  }
  retval = run_test(mristep_mem, y, T0, Tf, Hvals, method, reltol, abstol, udata);
  if (check_retval(&retval, "run_test", 1)) return 1;

  // Clean up and return
  MRIStepCoupling_Free(C);
  ARKodeFree(&inner_arkode_mem);
  MRIStepInnerStepper_Free(&inner_stepper);
  ARKodeFree(&mristep_mem);
  if (LS) { SUNLinSolFree(LS); } // free system linear solver
  if (A) { SUNMatDestroy(A); }   // free system matrix
  N_VDestroy(y);                 // Free y vector
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
  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];
  sunrealtype tmp1, tmp2;

  // fill in the RHS function:
  //   [G  e]*[(-2+u^2-p(t))/(2*u)] + [pdot(t)/(2u)]
  //   [e -1] [(-2+v^2-s(t))/(2*v)]   [qdot(t)/(2v)]
  tmp1      = (-TWO + u * u - p(t)) / (TWO * u);
  tmp2      = (-TWO + v * v - q(t, *udata)) / (TWO * v);
  dydata[0] = udata->G * tmp1 + udata->e * tmp2 + pdot(t) / (TWO * u);
  dydata[1] = udata->e * tmp1 - tmp2 + qdot(t, *udata) / (TWO * v);

  // Return with success
  return 0;
}

static int fi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  UserData* udata     = (UserData*)user_data;
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  sunrealtype* dydata = N_VGetArrayPointer(ydot);
  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];
  sunrealtype tmp1, tmp2;

  // fill in the RHS function:
  //   [G  e]*[(-2+u^2-p(t))/(2*u)] + [pdot(t)/(2u)]
  //   [0  0] [(-2+v^2-s(t))/(2*v)]   [     0      ]
  tmp1      = (-TWO + u * u - p(t)) / (TWO * u);
  tmp2      = (-TWO + v * v - q(t, *udata)) / (TWO * v);
  dydata[0] = udata->G * tmp1 + udata->e * tmp2 + pdot(t) / (TWO * u);
  dydata[1] = ZERO;

  // Return with success
  return 0;
}

static int fe(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  UserData* udata     = (UserData*)user_data;
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  sunrealtype* dydata = N_VGetArrayPointer(ydot);
  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];
  sunrealtype tmp1, tmp2;

  // fill in the RHS function:
  //   [0  0]*[(-2+u^2-p(t))/(2*u)] + [     0      ]
  //   [e -1] [(-2+v^2-s(t))/(2*v)]   [qdot(t)/(2v)]
  tmp1      = (-TWO + u * u - p(t)) / (TWO * u);
  tmp2      = (-TWO + v * v - q(t, *udata)) / (TWO * v);
  dydata[0] = ZERO;
  dydata[1] = udata->e * tmp1 - tmp2 + qdot(t, *udata) / (TWO * v);

  // Return with success
  return 0;
}

static int Jn(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData* udata     = (UserData*)user_data;
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];
  sunrealtype t11, t22;

  // fill in the Jacobian:
  //   [G  e]*[1-(u^2-p(t)-2)/(2*u^2),  0] + [-p'(t)/(2*u^2),  0]
  //   [e -1] [0,  1-(v^2-q(t)-2)/(2*v^2)]   [0,  -q'(t)/(2*v^2)]
  t11                   = ONE - (u * u - p(t) - TWO) / (TWO * u * u);
  t22                   = ONE - (v * v - q(t, *udata) - TWO) / (TWO * v * v);
  SM_ELEMENT_D(J, 0, 0) = udata->G * t11 - pdot(t) / (TWO * u * u);
  SM_ELEMENT_D(J, 0, 1) = udata->e * t22;
  SM_ELEMENT_D(J, 1, 0) = udata->e * t11;
  SM_ELEMENT_D(J, 1, 1) = -t22 - qdot(t, *udata) / (TWO * v * v);

  // Return with success
  return 0;
}

static int Ji(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData* udata     = (UserData*)user_data;
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];

  // fill in the Jacobian:
  //   [G/2 + (G*(1+p(t))-pdot(t))/(2*u^2)   e/2+e*(2+s(t))/(2*v^2)]
  //   [                 0                             0           ]
  SM_ELEMENT_D(J, 0, 0) = udata->G / TWO +
                          (udata->G * (ONE + p(t)) - pdot(t)) / (TWO * u * u);
  SM_ELEMENT_D(J, 0, 1) = udata->e / TWO +
                          udata->e * (TWO + q(t, *udata)) / (TWO * v * v);
  SM_ELEMENT_D(J, 1, 0) = ZERO;
  SM_ELEMENT_D(J, 1, 1) = ZERO;

  // Return with success
  return 0;
}

//------------------------------
// Private helper functions
//------------------------------

static int run_test(void* mristep_mem, N_Vector y, sunrealtype T0,
                    sunrealtype Tf, vector<sunrealtype> Hvals, string method,
                    sunrealtype reltol, sunrealtype abstol, UserData& udata)
{
  // Reused variables
  int retval;
  sunrealtype hpart = (Tf - T0) / udata.Npart;
  sunrealtype t;
  N_Vector ele       = N_VClone(y);
  N_Vector ewt       = N_VClone(y);
  N_Vector vtemp     = N_VClone(y);
  sunrealtype* ydata = N_VGetArrayPointer(y);

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
      // Reset integrator for this run
      t      = T0 + ipart * hpart;
      retval = Ytrue(t, y, udata);
      if (check_retval(&retval, "Ytrue", 1)) return 1;
      retval = ARKodeReset(mristep_mem, t, y);
      if (check_retval(&retval, "ARKodeReset", 1)) return 1;
      retval = ARKodeSetFixedStep(mristep_mem, Hvals[iH]);
      if (check_retval(&retval, "ARKodeSetFixedStep", 1)) return 1;
      retval = ARKodeResetAccumulatedError(mristep_mem);
      if (check_retval(&retval, "ARKodeResetAccumulatedError", 1)) return 1;

      // Run ARKodeEvolve to compute one step
      retval = ARKodeEvolve(mristep_mem, t + Hvals[iH], y, &t, ARK_ONE_STEP);
      if (check_retval(&retval, "ARKodeEvolve", 1)) return 1;
      retval = ARKodeGetEstLocalErrors(mristep_mem, ele);
      if (check_retval(&retval, "ARKodeGetEstLocalErrors", 1)) return 1;
      retval = computeErrorWeights(y, ewt, reltol, abstol, vtemp);
      if (check_retval(&retval, "computeErrorWeights", 1)) return 1;
      dsm_est[iH][ipart] = N_VWrmsNorm(ewt, ele);

      // Compute/print solution error
      sunrealtype udsm = abs(ydata[0] - utrue(t)) /
                         (abstol + reltol * abs(utrue(t)));
      sunrealtype vdsm = abs(ydata[1] - vtrue(t, udata)) /
                         (abstol + reltol * abs(vtrue(t, udata)));
      dsm[iH][ipart] = sqrt(0.5 * (udsm * udsm + vdsm * vdsm));
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
  return (0);
}

static sunrealtype p(sunrealtype t) { return (cos(t)); }

static sunrealtype q(sunrealtype t, UserData& udata)
{
  return (cos(udata.omega * t * (ONE + exp(-(t - 2) * (t - 2)))));
}

static sunrealtype pdot(sunrealtype t) { return (-sin(t)); }

static sunrealtype qdot(sunrealtype t, UserData& udata)
{
  return (-sin(udata.omega * t * (ONE + exp(-(t - 2) * (t - 2)))) * udata.omega *
          (ONE + exp(-(t - 2) * (t - 2)) -
           t * 2 * (t - 2) * (exp(-(t - 2) * (t - 2)))));
}

static sunrealtype utrue(sunrealtype t) { return (SUNRsqrt(TWO + p(t))); }

static sunrealtype vtrue(sunrealtype t, UserData& udata)
{
  return (SUNRsqrt(TWO + q(t, udata)));
}

static int Ytrue(sunrealtype t, N_Vector y, UserData& udata)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);
  ydata[0]           = utrue(t);
  ydata[1]           = vtrue(t, udata);
  return (0);
}

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
             retval < 0
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

//---- end of file ----
