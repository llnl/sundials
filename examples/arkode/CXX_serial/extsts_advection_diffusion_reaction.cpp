/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
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
 * This example simulates the same problem as in
 * ark_advection_diffusion_reaction.cpp, by David J. Gardner @ LLNL,
 *
 *   u_t = -c u_x + d u_xx + A - (w + 1) * u + v * u^2
 *   v_t = -c v_x + d v_xx + w * u - v * u^2
 *   w_t = -c w_x + d w_xx + (B - w) / eps - w * u
 *
 * where u, v, and w represent the concentrations of chemical species, c = 0.01
 * is the advection speed, d = 0.1 is the diffusion rate, and the species with
 * constant concentration over time are A = 0.6 and B = 2.0.
 *
 * The problem is evolved for t in [0, 3] and x in [0, 1], with initial
 * conditions given by
 *
 *   u(0,x) =  A  + 0.1 * sin(pi * x)
 *   v(0,x) = B/A + 0.1 * sin(pi * x)
 *   w(0,x) =  B  + 0.1 * sin(pi * x)
 *
 * and stationary boundary conditions i.e.,
 *
 *   u_t(t,0) = u_t(t,1) = 0,
 *   v_t(t,0) = v_t(t,1) = 0,
 *   w_t(t,0) = w_t(t,1) = 0.
 *
 * The system is advanced in time an Extended Super Time Stepping (ExtSTS)
 * method, wherein diffusion is treated explicitly using a STS method, advection
 * is treated explicitly using the ExtSTS method, and reaction is treated
 * implicitly using the ExtSTS method.
 *
 * Several command line options are available to change the problem parameters
 * and integrator settings. Use the flag --help for more information.
 * ---------------------------------------------------------------------------*/

#include "extsts_advection_diffusion_reaction.hpp"

int main(int argc, char* argv[])
{
  // SUNDIALS context object for this simulation
  sundials::Context ctx;

  // -----------------
  // Setup the problem
  // -----------------

  UserData udata;
  UserOptions uopts;

  vector<string> args(argv + 1, argv + argc);

  int flag = ReadInputs(args, udata, uopts, ctx);
  if (flag < 0)
  {
    cerr << "ERROR: ReadInputs returned " << flag << endl;
    return 1;
  }
  if (flag > 0) { return 0; }

  flag = PrintSetup(udata, uopts);
  if (check_flag(flag, "PrintSetup")) { return 1; }

  // Create state vector and set initial condition
  N_Vector y = N_VNew_Serial(udata.neq, ctx);
  if (check_ptr(y, "N_VNew_Serial")) { return 1; }

  flag = SetIC(y, udata);
  if (check_flag(flag, "SetIC")) { return 1; }

  // --------------------
  // Setup the integrator
  // --------------------

  // ARKODE memory structures
  void* arkode_mem = nullptr;
  void* arkref_mem = nullptr;

  // Matrix and linear solver for IMEX or ExtSTS integrators
  SUNMatrix A           = nullptr;
  SUNLinearSolver LS    = nullptr;
  SUNMatrix Aref        = nullptr;
  SUNLinearSolver LSref = nullptr;

  // Create integrator
  ARKRhsFn fe_RHS;   // explicit RHS function
  ARKRhsFn fi_RHS;   // implicit RHS function
  ARKLsJacFn Ji_RHS; // implicit RHS Jacobian function

  fe_RHS = (udata.advection) ? f_advection : nullptr;
  fi_RHS = (udata.reaction) ? f_reaction : nullptr;
  Ji_RHS = (udata.reaction) ? J_reaction : nullptr;

  // -------------------------------
  // Setup the ExtSTS integrator
  // -------------------------------

  // Create MRIStep ExtSTS solver memory
  arkode_mem = MRIStepCreateExtSTS(f_diffusion, fe_RHS, fi_RHS, ZERO, y, ctx);
  if (check_ptr(arkode_mem, "MRIStepCreateExtSTS")) { return 1; }

  // Access inner LSRKStep solver
  void* sts_mem = MRIStep_GetSTSStepper(arkode_mem);
  if (check_ptr(sts_mem, "MRIStep_GetSTSStepper")) { return 1; }

  // Attach user data (this attaches to both MRIStep and LSRKStep)
  flag = ARKodeSetUserData(arkode_mem, &udata);
  if (check_flag(flag, "ARKodeSetUserData")) { return 1; }

  // Select STS method
  ARKODE_LSRKMethodType ststype = (uopts.sts_method == 0) ? ARKODE_LSRK_RKC_2
                                                          : ARKODE_LSRK_RKL_2;
  flag = LSRKStepSetSTSMethod(sts_mem, ststype);
  if (check_flag(flag, "LSRKStepSetSTSMethod")) { return 1; }

  // Set dominant eigenvalue function and frequency
  flag = LSRKStepSetDomEigFn(sts_mem, diffusion_domeig);
  if (check_flag(flag, "LSRKStepSetDomEigFn")) { return 1; }
  flag = LSRKStepSetDomEigFrequency(sts_mem, uopts.ls_setup_freq);
  if (check_flag(flag, "LSRKStepSetDomEigFrequency")) { return 1; }

  // Set fixed step size
  if (uopts.fixed_h > ZERO)
  {
    flag = ARKodeSetFixedStep(arkode_mem, uopts.fixed_h);
    if (check_flag(flag, "ARKodeSetFixedStep")) { return 1; }
  }

  // Specify tolerances
  flag = ARKodeSStolerances(arkode_mem, uopts.rtol, uopts.atol);
  if (check_flag(flag, "ARKodeSStolerances")) { return 1; }

  // Attach user data
  flag = ARKodeSetUserData(arkode_mem, &udata);
  if (check_flag(flag, "ARKodeSetUserData")) { return 1; }

  // If implicit, setup solvers
  if (udata.reaction)
  {
    // Create banded matrix
    A = SUNBandMatrix(udata.neq, 2, 2, ctx);
    if (check_ptr(A, "SUNBandMatrix")) { return 1; }

    // Create linear solver
    LS = SUNLinSol_Band(y, A, ctx);
    if (check_ptr(LS, "SUNLinSol_Band")) { return 1; }

    // Attach linear solver
    flag = ARKodeSetLinearSolver(arkode_mem, LS, A);
    if (check_flag(flag, "ARKodeSetLinearSolver")) { return 1; }

    // Attach Jacobian function
    flag = ARKodeSetJacFn(arkode_mem, Ji_RHS);
    if (check_flag(flag, "ARKodeSetJacFn")) { return 1; }

    // Set linear solver setup frequency
    flag = ARKodeSetLSetupFrequency(arkode_mem, uopts.ls_setup_freq);
    if (check_flag(flag, "ARKodeSetLSetupFrequency")) { return 1; }

    // Tighten implicit solver tolerances
    flag = ARKodeSetNonlinConvCoef(arkode_mem, 1.e-1);
    if (check_flag(flag, "ARKodeSetNonlinConvCoef")) { return 1; }
    flag = ARKodeSetEpsLin(arkode_mem, 1.e-1);
    if (check_flag(flag, "ARKodeSetEpsLin")) { return 1; }

    // Use "deduce implicit RHS" option
    flag = ARKodeSetDeduceImplicitRhs(arkode_mem, SUNTRUE);
    if (check_flag(flag, "ARKodeSetDeduceImplicitRhs")) { return 1; }
  }

  // Select ExtSTS method via MRIStepCoupling structure
  MRIStepCoupling C;
  if (uopts.build_mri_table)
  {
    ARKodeButcherTable B = nullptr;
    if (udata.reaction)
    { B = ARKodeButcherTable_LoadDIRKByName(uopts.mri_method.c_str()); }
    else
    { B = ARKodeButcherTable_LoadERKByName(uopts.mri_method.c_str()); }
    if (check_ptr(B, "ARKodeButcherTable_Load")) { return 1; }
    C = MRIStepCoupling_MIStoMRI(B, B->q, B->p);
    if (check_ptr(C, "MRIStepCoupling_MIStoMRI")) { return 1; }
  }
  else
  {
    C = MRIStepCoupling_LoadTableByName(uopts.mri_method.c_str());
    if (check_ptr(C, "MRIStepCoupling_LoadTableByName")) { return 1; }
  }
  flag = MRIStepSetCoupling(arkode_mem, C);
  if (check_flag(flag, "MRIStepSetCoupling")) { return 1; }
  MRIStepCoupling_Free(C);

  // Set max steps between outputs
  flag = ARKodeSetMaxNumSteps(arkode_mem, uopts.maxsteps);
  if (check_flag(flag, "ARKodeSetMaxNumSteps")) { return 1; }

  // Tighten safety factor for time step selection
  flag = ARKodeSetSafetyFactor(arkode_mem, 0.8);
  if (check_flag(flag, "ARKodeSetSafetyFactor")) { return 1; }

  // Set stopping time
  flag = ARKodeSetStopTime(arkode_mem, udata.tf);
  if (check_flag(flag, "ARKodeSetStopTime")) { return 1; }


  // ----------------------
  // Evolve problem in time
  // ----------------------

  // Initial time, time between outputs, output time
  sunrealtype t     = ZERO;
  sunrealtype t2    = ZERO;
  sunrealtype dTout = udata.tf / uopts.nout;
  sunrealtype tout  = dTout;

  // Initial output
  flag = OpenOutput(udata, uopts);
  if (check_flag(flag, "OpenOutput")) { return 1; }

  flag = WriteOutput(t, y, udata, uopts);
  if (check_flag(flag, "WriteOutput")) { return 1; }

  // Loop over output times
  for (int iout = 0; iout < uopts.nout; iout++)
  {
    // Evolve
    flag = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_flag(flag, "ARKodeEvolve")) { return 1; }

    // Output solution
    flag = WriteOutput(t, y, udata, uopts);
    if (check_flag(flag, "WriteOutput")) { return 1; }

    // Update output time
    tout += dTout;
    tout = (tout > udata.tf) ? udata.tf : tout;
  }

  // Close output
  flag = CloseOutput(uopts);
  if (check_flag(flag, "CloseOutput")) { return 1; }

  // ------------
  // Output stats
  // ------------

  if (uopts.output)
  {
    cout << "Final integrator statistics:" << endl;
    cout << fixed << setprecision(6);
    cout << endl << "ExtSTS Integrator:" << endl;
    flag = ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
    if (check_flag(flag, "ARKodePrintAllStats")) { return -1; }
    cout << endl;
    cout << endl << "Inner STS Method:" << endl;
    flag = ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
    if (check_flag(flag, "ARKodePrintAllStats")) { return -1; }
  }

  // --------
  // Clean up
  // --------

  ARKodeFree(&arkode_mem);
  N_VDestroy(y);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);

  return 0;
}


// -----------------------------------------------------------------------------
// Functions called by the integrator
// -----------------------------------------------------------------------------

// Advection RHS function
int f_advection(sunrealtype t, N_Vector y, N_Vector f, void* user_data)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Access data arrays
  sunrealtype* ydata = N_VGetArrayPointer(y);
  if (check_ptr(ydata, "N_VGetArrayPointer")) { return -1; }

  sunrealtype* fdata = N_VGetArrayPointer(f);
  if (check_ptr(fdata, "N_VGetArrayPointer")) { return -1; }

  // Compute advection RHS
  sunrealtype ul, ur;
  sunrealtype vl, vr;
  sunrealtype wl, wr;

  sunrealtype c = -ONE * udata->c / (TWO * udata->dx);

  N_VConst(ZERO, f);
  for (sunindextype i = 1; i < udata->nx - 1; i++)
  {
    ul = ydata[UIDX(i - 1)];
    ur = ydata[UIDX(i + 1)];

    vl = ydata[VIDX(i - 1)];
    vr = ydata[VIDX(i + 1)];

    wl = ydata[WIDX(i - 1)];
    wr = ydata[WIDX(i + 1)];

    fdata[UIDX(i)] = c * (ur - ul);
    fdata[VIDX(i)] = c * (vr - vl);
    fdata[WIDX(i)] = c * (wr - wl);
  }

  return 0;
}

// Diffusion RHS function
int f_diffusion(sunrealtype t, N_Vector y, N_Vector f, void* user_data)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Access data arrays
  sunrealtype* ydata = N_VGetArrayPointer(y);
  if (check_ptr(ydata, "N_VGetArrayPointer")) { return -1; }

  sunrealtype* fdata = N_VGetArrayPointer(f);
  if (check_ptr(fdata, "N_VGetArrayPointer")) { return -1; }

  // Compute diffusion RHS
  sunrealtype ul, uc, ur;
  sunrealtype vl, vc, vr;
  sunrealtype wl, wc, wr;

  sunrealtype d = udata->d / (udata->dx * udata->dx);

  N_VConst(ZERO, f);
  for (sunindextype i = 1; i < udata->nx - 1; i++)
  {
    ul = ydata[UIDX(i - 1)];
    uc = ydata[UIDX(i)];
    ur = ydata[UIDX(i + 1)];

    vl = ydata[VIDX(i - 1)];
    vc = ydata[VIDX(i)];
    vr = ydata[VIDX(i + 1)];

    wl = ydata[WIDX(i - 1)];
    wc = ydata[WIDX(i)];
    wr = ydata[WIDX(i + 1)];

    fdata[UIDX(i)] = d * (ul - TWO * uc + ur);
    fdata[VIDX(i)] = d * (vl - TWO * vc + vr);
    fdata[WIDX(i)] = d * (wl - TWO * wc + wr);
  }

  return 0;
}

// Reaction RHS function
int f_reaction(sunrealtype t, N_Vector y, N_Vector f, void* user_data)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Access data arrays
  sunrealtype* ydata = N_VGetArrayPointer(y);
  if (check_ptr(ydata, "N_VGetArrayPointer")) { return -1; }

  sunrealtype* fdata = N_VGetArrayPointer(f);
  if (check_ptr(fdata, "N_VGetArrayPointer")) { return -1; }

  // Compute reaction RHS
  sunrealtype u, v, w;

  N_VConst(ZERO, f);
  for (sunindextype i = 1; i < udata->nx - 1; i++)
  {
    u = ydata[UIDX(i)];
    v = ydata[VIDX(i)];
    w = ydata[WIDX(i)];

    fdata[UIDX(i)] = udata->A - (w + ONE) * u + v * u * u;
    fdata[VIDX(i)] = w * u - v * u * u;
    fdata[WIDX(i)] = ((udata->B - w) / udata->eps) - w * u;
  }

  return 0;
}

// Reaction Jacobian function
int J_reaction(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Access data array
  sunrealtype* ydata = N_VGetArrayPointer(y);
  if (check_ptr(ydata, "N_VGetArrayPointer")) { return 1; }

  sunrealtype u, v, w;

  SUNMatZero(J);
  for (sunindextype i = 1; i < udata->nx - 1; i++)
  {
    u = ydata[UIDX(i)];
    v = ydata[VIDX(i)];
    w = ydata[WIDX(i)];

    // all vars wrt u
    SM_ELEMENT_B(J, UIDX(i), UIDX(i)) = -(w + ONE) + TWO * u * v;
    SM_ELEMENT_B(J, VIDX(i), UIDX(i)) = w - TWO * u * v;
    SM_ELEMENT_B(J, WIDX(i), UIDX(i)) = -w;

    // all vars wrt v
    SM_ELEMENT_B(J, UIDX(i), VIDX(i)) = u * u;
    SM_ELEMENT_B(J, VIDX(i), VIDX(i)) = -u * u;

    // all vars wrt w
    SM_ELEMENT_B(J, UIDX(i), WIDX(i)) = -u;
    SM_ELEMENT_B(J, VIDX(i), WIDX(i)) = u;
    SM_ELEMENT_B(J, WIDX(i), WIDX(i)) = (-ONE / udata->eps) - u;
  }

  return 0;
}

// Dominant eigenvalue function (for diffusion operator in LSRKStep)
int diffusion_domeig(sunrealtype t, N_Vector y, N_Vector fn,
                     sunrealtype* lambdaR, sunrealtype* lambdaI, void* user_data,
                     N_Vector temp1, N_Vector temp2, N_Vector temp3)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Fill in spectral radius value
  *lambdaR = -SUN_RCONST(4.0) * udata->d / udata->dx / udata->dx;
  *lambdaI = SUN_RCONST(0.0);

  return 0;
}

// Compute the initial condition
int SetIC(N_Vector y, UserData& udata)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);
  if (check_ptr(ydata, "N_VGetArrayPointer")) { return -1; }

  sunrealtype x, p;

  for (sunindextype i = 0; i < udata.nx; i++)
  {
    x              = udata.xl + i * udata.dx;
    p              = SUN_RCONST(0.1) * sin(M_PI * x);
    ydata[UIDX(i)] = udata.A + p;
    ydata[VIDX(i)] = udata.B / udata.A + p;
    ydata[WIDX(i)] = udata.B + p;
  }

  return 0;
}

//---- end of file ----
