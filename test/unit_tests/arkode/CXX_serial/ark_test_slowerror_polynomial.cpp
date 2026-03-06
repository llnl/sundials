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
 * estimate.  Uses a simple polynomial test problem with 1
 * component,
 *    dy/dt = a*t^2 + b*t + c,  y(0) = 1.
 *
 * We run a single time step, using a variety of slow step sizes,
 * H = {1, 1/2, 1/4, 1/8, 1/16}.
 *
 * We place the entire ODE in the "slow" RHS partition.  All tests
 * use ARKODE's default fifth-order ERK method, with relative and
 * absolute tolerances set to 1e-10 and 1e-12, respectively.
 *
 * When using an IMEX method (e.g., IMEX-MRI-SR), the entire slow
 * partition is treated implicitly.
 *
 * We select the slow integrator based on a command-line argument,
 * with the default being ARKODE_MRI_GARK_ERK22a.
 *
 * The program should be run with arguments in the following order:
 *   $ a.out method a b c
 * Not all arguments are required, but these must be omitted from
 * end-to-beginning, i.e. any one of
 *   $ a.out method a b
 *   $ a.out method a
 *   $ a.out method
 *   $ a.out
 * are acceptable.  We require that the "method" argument be a
 * string corresponding to a valid embedded ARKODE_MRITableID.
 *-----------------------------------------------------------------*/

// Header files
#include <algorithm>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_mristep.h>
#include <cmath>
#include <iomanip>
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

using namespace std;

// User data structure
struct UserData
{
  sunrealtype a;
  sunrealtype b;
  sunrealtype c;
};

// User-supplied Functions Called by the Solver
static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fn(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Jn(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Private utility functions
static int Ytrue(sunrealtype t, N_Vector y, UserData& udata);
static int computeErrorWeights(N_Vector ycur, N_Vector weight, sunrealtype rtol,
                               sunrealtype atol, N_Vector vtemp);
static int check_retval(void* returnvalue, const char* funcname, int opt);
static int run_test(void* mristep_mem, N_Vector y, sunrealtype T0,
                    vector<sunrealtype>& Hvals, string method,
                    sunrealtype reltol, sunrealtype abstol, UserData& udata);

// Main Program
int main(int argc, char* argv[])
{
  // general problem parameters
  sunrealtype T0   = ZERO;                 // initial time
  sunindextype NEQ = 1;                    // number of dependent vars.
  string method;                           // MRI method name
  sunrealtype reltol = SUN_RCONST(1.e-10); // fast solver tolerances
  sunrealtype abstol = SUN_RCONST(1.e-12);

  // general problem variables
  int retval;     // reusable error-checking flag
  UserData udata; // user-data structure

  //
  // Initialization
  //

  // Retrieve the command-line options:  method Npart ep test
  if (argc > 1) { method = argv[1]; }
  else { method = "ARKODE_MRI_GARK_ERK22a"; }
  if (argc > 2) { udata.a = SUNStrToReal(argv[2]); }
  else { udata.a = ONE; }
  if (argc > 3) { udata.b = SUNStrToReal(argv[3]); }
  else { udata.b = ONE; }
  if (argc > 4) { udata.c = SUNStrToReal(argv[4]); }
  else { udata.c = ONE; }

  sunbooleantype implicit = SUNFALSE;
  if ((method == "ARKODE_MRI_GARK_IRK21a") ||
      (method == "ARKODE_MRI_GARK_ESDIRK34a") ||
      (method == "ARKODE_MRI_GARK_ESDIRK46a") ||
      (method == "ARKODE_IMEX_MRI_SR21") ||
      (method == "ARKODE_IMEX_MRI_SR32") || (method == "ARKODE_IMEX_MRI_SR43"))
  {
    implicit = SUNTRUE;
  }

  // Initial problem output (and set implicit solver tolerances as needed)
  cout << "\nSlow error estimation test (polynomial ODE problem):\n";
  cout << "    problem parameters:  a = " << udata.a << ",  b = " << udata.b
       << ",  c = " << udata.c << endl;
  cout << "    MRI method: " << method;
  if (implicit) { cout << " (implicit)" << endl; }
  else { cout << " (explicit)" << endl; }

  //
  // Problem Setup
  //

  // Create SUNDIALS context
  sundials::Context ctx;

  // Create serial vectors for the solution and reference
  N_Vector y = N_VNew_Serial(NEQ, ctx);
  if (check_retval((void*)y, "N_VNew_Serial", 0)) return 1;
  retval = Ytrue(T0, y, udata);
  if (check_retval(&retval, "Ytrue", 1)) return 1;

  // Set up fast ARKStep integrator as fifth-order adaptive ERK
  void* inner_arkode_mem = ARKStepCreate(f0, NULL, T0, y, ctx);
  if (check_retval((void*)inner_arkode_mem, "ARKStepCreate", 0)) return 1;
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
  if (implicit)
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
    retval = ARKodeSetJacFn(mristep_mem, Jn);
    if (check_retval(&retval, "ARKodeSetJacFn", 1)) return 1;
    retval = ARKodeSetJacEvalFrequency(mristep_mem, 1);
    if (check_retval(&retval, "ARKodeSetJacEvalFrequency", 1)) return 1;
    retval = ARKodeSetLSetupFrequency(mristep_mem, 1);
    if (check_retval(&retval, "ARKodeSetLSetupFrequency", 1)) return 1;
    retval = ARKodeSetMaxNonlinIters(mristep_mem, 20);
    if (check_retval(&retval, "ARKodeSetMaxNonlinIters", 1)) return 1;
  }
  retval = ARKodeSStolerances(mristep_mem, reltol, abstol);
  if (check_retval(&retval, "ARKodeSStolerances", 1)) return 1;
  retval = ARKodeSetUserData(mristep_mem, (void*)&udata);
  if (check_retval(&retval, "ARKodeSetUserData", 1)) return 1;
  retval = ARKodeSetAccumulatedErrorType(mristep_mem, ARK_ACCUMERROR_MAX);
  if (check_retval(&retval, "ARKodeSetAccumulatedErrorType", 1)) return 1;

  // Run test for various H values
  vector<sunrealtype> Hvals(5);
  for (size_t i = 0; i < Hvals.size(); i++)
  {
    Hvals[i] = SUN_RCONST(0.01) / SUNRpowerI(SUN_RCONST(2.0), (int)i);
  }
  retval = run_test(mristep_mem, y, T0, Hvals, method, reltol, abstol, udata);
  if (check_retval(&retval, "run_test", 1)) return 1;

  // Clean up and return
  MRIStepCoupling_Free(C);
  ARKodeFree(&inner_arkode_mem);
  MRIStepInnerStepper_Free(&inner_stepper);
  ARKodeFree(&mristep_mem);
  if (LS) { SUNLinSolFree(LS); } // free system linear solver
  if (A) { SUNMatDestroy(A); }   // free system matrix
  N_VDestroy(y);                 // Free y and yref vectors
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
  // fill in the RHS function and return with success
  UserData* udata     = (UserData*)user_data;
  sunrealtype* dydata = N_VGetArrayPointer(ydot);
  dydata[0]           = udata->a * t * t + udata->b * t + udata->c;
  return 0;
}

static int Jn(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  // fill in the Jacobian and return with success
  SM_ELEMENT_D(J, 0, 0) = ZERO;
  return 0;
}

//------------------------------
// Private helper functions
//------------------------------

static int run_test(void* mristep_mem, N_Vector y, sunrealtype T0,
                    vector<sunrealtype>& Hvals, string method,
                    sunrealtype reltol, sunrealtype abstol, UserData& udata)
{
  // Reused variables
  int retval;
  sunrealtype t;
  N_Vector vtemp     = N_VClone(y);
  N_Vector ele       = N_VClone(y);
  N_Vector ewt       = N_VClone(y);
  sunrealtype* ydata = N_VGetArrayPointer(y);
  sunrealtype* vdata = N_VGetArrayPointer(vtemp);

  // Set storage for errors
  vector<sunrealtype> dsm(Hvals.size(), ZERO);
  vector<sunrealtype> dsm_est(Hvals.size(), ZERO);

  // Loop over step sizes
  for (size_t iH = 0; iH < Hvals.size(); iH++)
  {
    // Reset integrator for this run
    t      = T0;
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
    dsm_est[iH] = N_VWrmsNorm(ewt, ele);

    // Compute/print solution error
    retval = Ytrue(t, vtemp, udata);
    if (check_retval(&retval, "Ytrue", 1)) return 1;
    dsm[iH] = abs(ydata[0] - vdata[0]) / (abstol + reltol * abs(vdata[0]));
    if (method == "ARKODE_MRI_GARK_ERK22a")
    {
      cout << fixed << setprecision(5);
      cout << setw(10) << "H  " << Hvals[iH];
      cout << scientific << setprecision(2);
      cout << setw(9) << "dsm  " << dsm[iH];
      cout << setw(13) << "dsm_est  " << dsm_est[iH];
      cout << setw(14) << "dsm_anal  "
           << Hvals[iH] * Hvals[iH] * Hvals[iH] * abs(udata.a / 12.0) /
                (abstol + reltol * abs(vdata[0]));
      cout << setw(18) << "dsm_est_anal  "
           << Hvals[iH] * Hvals[iH] *
                abs(udata.a * Hvals[iH] / 4.0 + udata.b / 2.0) /
                (abstol + reltol * abs(vdata[0]));
      cout << endl;
    }
    else if (method == "ARKODE_MRI_GARK_IRK21a")
    {
      cout << fixed << setprecision(5);
      cout << setw(10) << "H  " << Hvals[iH];
      cout << scientific << setprecision(2);
      cout << setw(9) << "dsm  " << dsm[iH];
      cout << setw(13) << "dsm_est  " << dsm_est[iH];
      cout << setw(14) << "dsm_anal  "
           << Hvals[iH] * Hvals[iH] * Hvals[iH] * abs(udata.a / 6.0) /
                (abstol + reltol * abs(vdata[0]));
      cout << setw(18) << "dsm_est_anal  "
           << Hvals[iH] * Hvals[iH] *
                abs(udata.a * Hvals[iH] / 2.0 + udata.b / 2.0) /
                (abstol + reltol * abs(vdata[0]));
      cout << endl;
    }
    else if (method == "ARKODE_IMEX_MRI_SR21")
    {
      cout << fixed << setprecision(5);
      cout << setw(10) << "H  " << Hvals[iH];
      cout << scientific << setprecision(2);
      cout << setw(9) << "dsm  " << dsm[iH];
      cout << setw(13) << "dsm_est  " << dsm_est[iH];
      cout << setw(14) << "dsm_anal  "
           << Hvals[iH] * Hvals[iH] * Hvals[iH] *
                abs(udata.a * 3137.0 / 50370.0) /
                (abstol + reltol * abs(vdata[0]));
      cout << setw(18) << "dsm_est_anal  "
           << Hvals[iH] * Hvals[iH] * Hvals[iH] *
                abs(udata.a * 3137.0 / 50370.0) /
                (abstol + reltol * abs(vdata[0]));
      cout << endl;
    }
    else
    {
      cout << fixed << setprecision(5);
      cout << setw(10) << "H  " << Hvals[iH];
      cout << scientific << setprecision(2);
      cout << setw(9) << "dsm  " << dsm[iH];
      cout << setw(13) << "dsm_est  " << dsm_est[iH];
      cout << endl;
    }
  }

  N_VDestroy(ele);
  N_VDestroy(ewt);
  N_VDestroy(vtemp);
  return (0);
}

static int Ytrue(sunrealtype t, N_Vector y, UserData& udata)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);
  ydata[0]           = udata.a / SUN_RCONST(3.0) * t * t * t +
             udata.b / SUN_RCONST(2.0) * t * t + udata.c * t + ONE;
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
