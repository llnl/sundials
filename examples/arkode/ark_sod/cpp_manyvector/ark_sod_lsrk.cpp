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
 * This example performs a standard 1D sod shock tube (see section 5.2 of J.A.
 * Greenough and W.J. Rider, "A quantitative comparison of numerical methods for
 * the compressible Euler equations: fifth-order WENO and piecewise-linear
 * Godunov," J. Comput. Phys., 196:259-281, 2004)
 *    [rho, vx, p] = { [1, 0, 1]        if x < 0.5
 *                   { [0.125, 0, 0.1]  if x > 0.5
 *
 * The code solves the 1D compressible Euler equations in conserved variables,
 * over the domain (t,x) in [0, 0.2] x [0, 1].
 *
 * Since the Sod shock tube is specified in terms of primitive variables, we
 * convert between primitive and conserved variables for the initial conditions
 * and accuracy results.
 *
 * This problem should be run with homogeneous Neumann boundary conditions.
 *
 * The system is advanced in time using one of the strong-stability-preserving
 * Runge--Kutta methods from LSRKStep, based on the --integrator METHOD input
 * value. The following options are available:
 *
 *   SSPRK(s,2) -- specified via METHOD = ARKODE_LSRK_SSP_S_2.  The number of
 *                 stages to use defaults to 10.
 *
 *   SSPRK(s,3) -- specified via METHOD = ARKODE_LSRK_SSP_S_3.  The number of
 *                 stages to use defaults to 9.
 *
 *   SSPRK(10,4) -- specified via METHOD = ARKODE_LSRK_SSP_10_4.
 *
 * Both the SSPRK(s,2) and SSPRK(s,3) methods allow specification of a
 * non-default number of stages.  This may be specified using the --stages S
 * input value, where S is an integer.  Note: SSPRK(s,2) requires S at least
 * 2, and SSPRK(s,3) requires S be a perfect square, with S at least 4.
 *
 * Alternately, if METHOD corresponds with a valid ARKODE_ERKTableID then
 * the system will be advanced using that method in ERKStep.
 *
 * Several additional command line options are available to change the
 * and integrator settings. Use the flag --help for more information.
 * ---------------------------------------------------------------------------*/

#include "ark_sod_lsrk.hpp"
using namespace std;

int main(int argc, char* argv[])
{
  // SUNDIALS context object for this simulation
  sundials::Context ctx;

  // -----------------
  // Setup the problem
  // -----------------

  EulerData udata;
  ARKODEParameters uopts;

  vector<string> args(argv + 1, argv + argc);

  int flag = ReadInputs(args, udata, uopts, ctx);
  if (check_flag(flag, "ReadInputs")) { return 1; }
  if (flag > 0) { return 0; }

  flag = PrintSetup(udata, uopts);
  if (check_flag(flag, "PrintSetup")) { return 1; }

  // Create state vector and set initial condition
  N_Vector vecs[NSPECIES];
  for (int i = 0; i < NSPECIES; i++)
  {
    vecs[i] = N_VNew_Serial((sunindextype)udata.nx, ctx); // rho (density)
    if (check_ptr(vecs[i], "N_VNew_Serial")) { return 1; }
  }
  N_Vector y = N_VNew_ManyVector(NSPECIES, vecs, ctx);
  if (check_ptr(y, "N_VNew_ManyVector")) { return 1; }

  flag = SetIC(y, udata);
  if (check_flag(flag, "SetIC")) { return 1; }

  // --------------------
  // Setup the integrator
  // --------------------
  void* arkode_mem = nullptr;

  // Determine type (LSRKStep vs ERKStep)
  bool lsrk = false;
  if (uopts.integrator == "ARKODE_LSRK_SSP_S_2") { lsrk = true; }
  if (uopts.integrator == "ARKODE_LSRK_SSP_S_3") { lsrk = true; }
  if (uopts.integrator == "ARKODE_LSRK_SSP_10_4") { lsrk = true; }

  if (lsrk) // Setup LSRKStep
  {
    // ARKODE memory structure
    arkode_mem = LSRKStepCreateSSP(frhs, udata.t0, y, ctx);
    if (check_ptr(arkode_mem, "LSRKStepCreateSSP")) { return 1; }

    // Select SSPRK method type
    flag = LSRKStepSetSSPMethodByName(arkode_mem, uopts.integrator.c_str());
    if (check_flag(flag, "LSRKStepSetSSPMethodByName")) { return 1; }

    // Select number of SSPRK stages
    if (uopts.stages > 0)
    {
      flag = LSRKStepSetNumSSPStages(arkode_mem, uopts.stages);
      if (check_flag(flag, "LSRKStepSetNumSSPStages")) { return 1; }
    }
  }
  else
  { // Setup ERKStep

    // ARKODE memory structure
    arkode_mem = ERKStepCreate(frhs, udata.t0, y, ctx);
    if (check_ptr(arkode_mem, "ERKStepCreate")) { return 1; }

    // Select ERK method
    flag = ERKStepSetTableName(arkode_mem, uopts.integrator.c_str());
    if (check_flag(flag, "ERKStepSetTableName")) { return 1; }
  }

  // Shared setup

  //   Specify tolerances
  flag = ARKodeSStolerances(arkode_mem, uopts.rtol, uopts.atol);
  if (check_flag(flag, "ARKodeSStolerances")) { return 1; }

  //   Attach user data
  flag = ARKodeSetUserData(arkode_mem, &udata);
  if (check_flag(flag, "ARKodeSetUserData")) { return 1; }

  //   Set fixed step size or adaptivity method
  if (uopts.fixed_h > ZERO)
  {
    flag = ARKodeSetFixedStep(arkode_mem, uopts.fixed_h);
    if (check_flag(flag, "ARKodeSetFixedStep")) { return 1; }
  }

  //   Set max steps between outputs
  flag = ARKodeSetMaxNumSteps(arkode_mem, uopts.maxsteps);
  if (check_flag(flag, "ARKodeSetMaxNumSteps")) { return 1; }

  //   Set stopping time
  flag = ARKodeSetStopTime(arkode_mem, udata.tf);
  if (check_flag(flag, "ARKodeSetStopTime")) { return 1; }

  // ----------------------
  // Evolve problem in time
  // ----------------------

  // Initial time, time between outputs, output time
  sunrealtype t     = ZERO;
  sunrealtype dTout = udata.tf / uopts.nout;
  sunrealtype tout  = dTout;

  // initial output
  flag = OpenOutput(udata, uopts);
  if (check_flag(flag, "OpenOutput")) { return 1; }

  flag = WriteOutput(t, y, udata, uopts);
  if (check_flag(flag, "WriteOutput")) { return 1; }

  // Loop over output times
  for (int iout = 0; iout < uopts.nout; iout++)
  {
    // Evolve
    if (uopts.output == 3)
    {
      // Stop at output time (do not interpolate output)
      flag = ARKodeSetStopTime(arkode_mem, tout);
      if (check_flag(flag, "ARKodeSetStopTime")) { return 1; }
    }

    //   Advance in time
    flag = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_flag(flag, "ARKodeEvolve")) { break; }

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
    flag = ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  }

  // --------
  // Clean up
  // --------

  ARKodeFree(&arkode_mem);
  for (int i = 0; i < NSPECIES; i++) { N_VDestroy(vecs[i]); }
  N_VDestroy(y);

  return 0;
}

// -----------------------------------------------------------------------------
// Functions called by the integrator
// -----------------------------------------------------------------------------

// ODE RHS function
int frhs(sunrealtype t, N_Vector y, N_Vector f, void* user_data)
{
  // Access problem data
  EulerData* udata = (EulerData*)user_data;

  // initialize output to zeros
  N_VConst(ZERO, f);

  // Access data arrays
  sunrealtype* rho = N_VGetSubvectorArrayPointer_ManyVector(y, 0);
  if (check_ptr(rho, "N_VGetSubvectorArrayPointer_ManyVector")) { return -1; }
  sunrealtype* mx = N_VGetSubvectorArrayPointer_ManyVector(y, 1);
  if (check_ptr(mx, "N_VGetSubvectorArrayPointer_ManyVector")) { return -1; }
  sunrealtype* my = N_VGetSubvectorArrayPointer_ManyVector(y, 2);
  if (check_ptr(my, "N_VGetSubvectorArrayPointer_ManyVector")) { return -1; }
  sunrealtype* mz = N_VGetSubvectorArrayPointer_ManyVector(y, 3);
  if (check_ptr(mz, "N_VGetSubvectorArrayPointer_ManyVector")) { return -1; }
  sunrealtype* et = N_VGetSubvectorArrayPointer_ManyVector(y, 4);
  if (check_ptr(et, "N_VGetSubvectorArrayPointer_ManyVector")) { return -1; }

  sunrealtype* rhodot = N_VGetSubvectorArrayPointer_ManyVector(f, 0);
  if (check_ptr(rhodot, "N_VGetSubvectorArrayPointer_ManyVector"))
  {
    return -1;
  }
  sunrealtype* mxdot = N_VGetSubvectorArrayPointer_ManyVector(f, 1);
  if (check_ptr(mxdot, "N_VGetSubvectorArrayPointer_ManyVector")) { return -1; }
  sunrealtype* mydot = N_VGetSubvectorArrayPointer_ManyVector(f, 2);
  if (check_ptr(mydot, "N_VGetSubvectorArrayPointer_ManyVector")) { return -1; }
  sunrealtype* mzdot = N_VGetSubvectorArrayPointer_ManyVector(f, 3);
  if (check_ptr(mzdot, "N_VGetSubvectorArrayPointer_ManyVector")) { return -1; }
  sunrealtype* etdot = N_VGetSubvectorArrayPointer_ManyVector(f, 4);
  if (check_ptr(etdot, "N_VGetSubvectorArrayPointer_ManyVector")) { return -1; }

  // Set shortcut variables
  const long int nx    = udata->nx;
  const sunrealtype dx = udata->dx;
  sunrealtype* flux    = udata->flux;

  // compute face-centered fluxes over domain interior: pack 1D x-directional array
  // of variable shortcuts, and compute flux at lower x-directional face
  for (long int i = 3; i < nx - 2; i++)
  {
    udata->pack1D(rho, mx, my, mz, et, i);
    face_flux(udata->w1d, &(flux[i * NSPECIES]), *udata);
  }

  // compute face-centered fluxes at left boundary
  for (long int i = 0; i < 3; i++)
  {
    udata->pack1D_bdry(rho, mx, my, mz, et, i);
    face_flux(udata->w1d, &(flux[i * NSPECIES]), *udata);
  }

  // compute face-centered fluxes at right boundary
  for (long int i = nx - 2; i <= nx; i++)
  {
    udata->pack1D_bdry(rho, mx, my, mz, et, i);
    face_flux(udata->w1d, &(flux[i * NSPECIES]), *udata);
  }

  // iterate over subdomain, updating RHS
  for (long int i = 0; i < nx; i++)
  {
    rhodot[i] -= (flux[(i + 1) * NSPECIES + 0] - flux[i * NSPECIES + 0]) / dx;
    mxdot[i] -= (flux[(i + 1) * NSPECIES + 1] - flux[i * NSPECIES + 1]) / dx;
    mydot[i] -= (flux[(i + 1) * NSPECIES + 2] - flux[i * NSPECIES + 2]) / dx;
    mzdot[i] -= (flux[(i + 1) * NSPECIES + 3] - flux[i * NSPECIES + 3]) / dx;
    etdot[i] -= (flux[(i + 1) * NSPECIES + 4] - flux[i * NSPECIES + 4]) / dx;
  }

  return 0;
}

// given a 6-point stencil of solution values,
//   w(x_{j-3}) w(x_{j-2}) w(x_{j-1}), w(x_j), w(x_{j+1}), w(x_{j+2}),
// compute the face-centered flux (flux) at the center of the stencil, x_{j-1/2}.
//
// This precisely follows the recipe laid out in:
// Chi-Wang Shu (2003) "High-order Finite Difference and Finite Volume WENO
// Schemes and Discontinuous Galerkin Methods for CFD," International Journal of
// Computational Fluid Dynamics, 17:2, 107-118, DOI: 10.1080/1061856031000104851
// with the only change that since this is 1D, we manually set the y- and
// z-velocities, v and w, to zero.
void face_flux(sunrealtype (&w1d)[STSIZE][NSPECIES], sunrealtype* f_face,
               const EulerData& udata)
{
  // local data
  int i, j;
  sunrealtype rhosqrL, rhosqrR, rhosqrbar, u, v, w, H, qsq, csnd, cinv, gamm,
    alpha, beta1, beta2, beta3, w1, w2, w3, f1, f2, f3;
  sunrealtype RV[5][5], LV[5][5], p[STSIZE], flux[STSIZE][NSPECIES],
    fproj[5][NSPECIES], fs[5][NSPECIES], ff[NSPECIES];
  const sunrealtype bc =
    SUN_RCONST(1.083333333333333333333333333333333333333); // 13/12
  const sunrealtype epsilon = SUN_RCONST(1.0e-6);

  // compute pressures over stencil
  for (i = 0; i < STSIZE; i++)
  {
    p[i] = udata.eos(w1d[i][0], w1d[i][1], w1d[i][2], w1d[i][3], w1d[i][4]);
  }

  // compute Roe-average state at face:
  //   wbar = [sqrt(rho), sqrt(rho)*vx, sqrt(rho)*vy, sqrt(rho)*vz, (e+p)/sqrt(rho)]
  //          [sqrt(rho), mx/sqrt(rho), my/sqrt(rho), mz/sqrt(rho), (e+p)/sqrt(rho)]
  //   u = wbar_2 / wbar_1
  //   v = wbar_3 / wbar_1
  //   w = wbar_4 / wbar_1
  //   H = wbar_5 / wbar_1
  rhosqrL   = sqrt(w1d[2][0]);
  rhosqrR   = sqrt(w1d[3][0]);
  rhosqrbar = HALF * (rhosqrL + rhosqrR);
  u         = HALF * (w1d[2][1] / rhosqrL + w1d[3][1] / rhosqrR) / rhosqrbar;
  v         = HALF * (w1d[2][2] / rhosqrL + w1d[3][2] / rhosqrR) / rhosqrbar;
  w         = HALF * (w1d[2][3] / rhosqrL + w1d[3][3] / rhosqrR) / rhosqrbar;
  H = HALF * ((p[2] + w1d[2][4]) / rhosqrL + (p[3] + w1d[3][4]) / rhosqrR) /
      rhosqrbar;

  // compute eigenvectors at face (note: eigenvectors for tracers are just identity)
  qsq  = u * u + v * v + w * w;
  gamm = udata.gamma - ONE;
  csnd = gamm * (H - HALF * qsq);
  cinv = ONE / csnd;
  for (i = 0; i < 5; i++)
  {
    for (j = 0; j < 5; j++)
    {
      RV[i][j] = ZERO;
      LV[i][j] = ZERO;
    }
  }

  RV[0][0] = ONE;
  RV[0][3] = ONE;
  RV[0][4] = ONE;

  RV[1][0] = u - csnd;
  RV[1][3] = u;
  RV[1][4] = u + csnd;

  RV[2][0] = v;
  RV[2][1] = ONE;
  RV[2][3] = v;
  RV[2][4] = v;

  RV[3][0] = w;
  RV[3][2] = ONE;
  RV[3][3] = w;
  RV[3][4] = w;

  RV[4][0] = H - u * csnd;
  RV[4][1] = v;
  RV[4][2] = w;
  RV[4][3] = HALF * qsq;
  RV[4][4] = H + u * csnd;

  LV[0][0] = HALF * cinv * (u + HALF * gamm * qsq);
  LV[0][1] = -HALF * cinv * (gamm * u + ONE);
  LV[0][2] = -HALF * v * gamm * cinv;
  LV[0][3] = -HALF * w * gamm * cinv;
  LV[0][4] = HALF * gamm * cinv;

  LV[1][0] = -v;
  LV[1][2] = ONE;

  LV[2][0] = -w;
  LV[2][3] = ONE;

  LV[3][0] = -gamm * cinv * (qsq - H);
  LV[3][1] = u * gamm * cinv;
  LV[3][2] = v * gamm * cinv;
  LV[3][3] = w * gamm * cinv;
  LV[3][4] = -gamm * cinv;

  LV[4][0] = -HALF * cinv * (u - HALF * gamm * qsq);
  LV[4][1] = -HALF * cinv * (gamm * u - ONE);
  LV[4][2] = -HALF * v * gamm * cinv;
  LV[4][3] = -HALF * w * gamm * cinv;
  LV[4][4] = HALF * gamm * cinv;

  // compute fluxes and max wave speed over stencil
  alpha = ZERO;
  for (j = 0; j < STSIZE; j++)
  {
    u          = w1d[j][1] / w1d[j][0];  // u = vx = mx/rho
    flux[j][0] = w1d[j][1];              // f_rho = rho*u = mx
    flux[j][1] = u * w1d[j][1] + p[j];   // f_mx = rho*u*u + p = mx*u + p
    flux[j][2] = u * w1d[j][2];          // f_my = rho*v*u = my*u
    flux[j][3] = u * w1d[j][3];          // f_mz = rho*w*u = mz*u
    flux[j][4] = u * (w1d[j][4] + p[j]); // f_et = u*(et + p)
    csnd  = sqrt(udata.gamma * p[j] / w1d[j][0]); // csnd = sqrt(gamma*p/rho)
    alpha = max(alpha, abs(u) + csnd);
  }

  // compute flux from right side of face at x_{i+1/2}:

  //   compute right-shifted Lax-Friedrichs flux over left portion of patch
  for (j = 0; j < 5; j++)
  {
    for (i = 0; i < NSPECIES; i++)
    {
      fs[j][i] = HALF * (flux[j][i] + alpha * w1d[j][i]);
    }
  }

  // compute projected flux for fluid fields
  for (j = 0; j < 5; j++)
  {
    for (i = 0; i < 5; i++)
    {
      fproj[j][i] = LV[i][0] * fs[j][0] + LV[i][1] * fs[j][1] +
                    LV[i][2] * fs[j][2] + LV[i][3] * fs[j][3] +
                    LV[i][4] * fs[j][4];
    }
  }

  //   compute WENO signed flux
  for (i = 0; i < NSPECIES; i++)
  {
    // smoothness indicators
    beta1 = bc * pow(fproj[2][i] - SUN_RCONST(2.0) * fproj[3][i] + fproj[4][i],
                     2) +
            FOURTH * pow(SUN_RCONST(3.0) * fproj[2][i] -
                           SUN_RCONST(4.0) * fproj[3][i] + fproj[4][i],
                         2);
    beta2 = bc * pow(fproj[1][i] - SUN_RCONST(2.0) * fproj[2][i] + fproj[3][i],
                     2) +
            FOURTH * pow(fproj[1][i] - fproj[3][i], 2);
    beta3 = bc * pow(fproj[0][i] - SUN_RCONST(2.0) * fproj[1][i] + fproj[2][i],
                     2) +
            FOURTH * pow(fproj[0][i] - SUN_RCONST(4.0) * fproj[1][i] +
                           SUN_RCONST(3.0) * fproj[2][i],
                         2);
    // nonlinear weights
    w1 = SUN_RCONST(0.3) / ((epsilon + beta1) * (epsilon + beta1));
    w2 = SUN_RCONST(0.6) / ((epsilon + beta2) * (epsilon + beta2));
    w3 = SUN_RCONST(0.1) / ((epsilon + beta3) * (epsilon + beta3));
    // flux stencils
    f1 = SUN_RCONST(0.3333333333333333333333333333333333333333) * fproj[2][i] +
         SUN_RCONST(0.8333333333333333333333333333333333333333) * fproj[3][i] -
         SUN_RCONST(0.1666666666666666666666666666666666666667) * fproj[4][i];
    f2 = -SUN_RCONST(0.1666666666666666666666666666666666666667) * fproj[1][i] +
         SUN_RCONST(0.8333333333333333333333333333333333333333) * fproj[2][i] +
         SUN_RCONST(0.3333333333333333333333333333333333333333) * fproj[3][i];
    f3 = SUN_RCONST(0.3333333333333333333333333333333333333333) * fproj[0][i] -
         SUN_RCONST(1.166666666666666666666666666666666666667) * fproj[1][i] +
         SUN_RCONST(1.833333333333333333333333333333333333333) * fproj[2][i];
    // resulting signed flux at face
    ff[i] = (f1 * w1 + f2 * w2 + f3 * w3) / (w1 + w2 + w3);
  }

  // compute flux from left side of face at x_{i+1/2}:

  //   compute left-shifted Lax-Friedrichs flux over right portion of patch
  for (j = 0; j < 5; j++)
  {
    for (i = 0; i < NSPECIES; i++)
    {
      fs[j][i] = HALF * (flux[j + 1][i] - alpha * w1d[j + 1][i]);
    }
  }

  // compute projected flux for fluid fields
  for (j = 0; j < 5; j++)
  {
    for (i = 0; i < 5; i++)
    {
      fproj[j][i] = LV[i][0] * fs[j][0] + LV[i][1] * fs[j][1] +
                    LV[i][2] * fs[j][2] + LV[i][3] * fs[j][3] +
                    LV[i][4] * fs[j][4];
    }
  }

  //   compute WENO signed fluxes
  for (i = 0; i < NSPECIES; i++)
  {
    // smoothness indicators
    beta1 = bc * pow(fproj[2][i] - SUN_RCONST(2.0) * fproj[3][i] + fproj[4][i],
                     2) +
            FOURTH * pow(SUN_RCONST(3.0) * fproj[2][i] -
                           SUN_RCONST(4.0) * fproj[3][i] + fproj[4][i],
                         2);
    beta2 = bc * pow(fproj[1][i] - SUN_RCONST(2.0) * fproj[2][i] + fproj[3][i],
                     2) +
            FOURTH * pow(fproj[1][i] - fproj[3][i], 2);
    beta3 = bc * pow(fproj[0][i] - SUN_RCONST(2.0) * fproj[1][i] + fproj[2][i],
                     2) +
            FOURTH * pow(fproj[0][i] - SUN_RCONST(4.0) * fproj[1][i] +
                           SUN_RCONST(3.0) * fproj[2][i],
                         2);
    // nonlinear weights
    w1 = SUN_RCONST(0.1) / ((epsilon + beta1) * (epsilon + beta1));
    w2 = SUN_RCONST(0.6) / ((epsilon + beta2) * (epsilon + beta2));
    w3 = SUN_RCONST(0.3) / ((epsilon + beta3) * (epsilon + beta3));
    // flux stencils
    f1 = SUN_RCONST(1.833333333333333333333333333333333333333) * fproj[2][i] -
         SUN_RCONST(1.166666666666666666666666666666666666667) * fproj[3][i] +
         SUN_RCONST(0.3333333333333333333333333333333333333333) * fproj[4][i];
    f2 = SUN_RCONST(0.3333333333333333333333333333333333333333) * fproj[1][i] +
         SUN_RCONST(0.8333333333333333333333333333333333333333) * fproj[2][i] -
         SUN_RCONST(0.1666666666666666666666666666666666666667) * fproj[3][i];
    f3 = -SUN_RCONST(0.1666666666666666666666666666666666666667) * fproj[0][i] +
         SUN_RCONST(0.8333333333333333333333333333333333333333) * fproj[1][i] +
         SUN_RCONST(0.3333333333333333333333333333333333333333) * fproj[2][i];
    // resulting signed flux (add to ff)
    ff[i] += (f1 * w1 + f2 * w2 + f3 * w3) / (w1 + w2 + w3);
  }

  // combine signed fluxes into output, converting back to conserved variables
  for (i = 0; i < NSPECIES; i++)
  {
    f_face[i] = RV[i][0] * ff[0] + RV[i][1] * ff[1] + RV[i][2] * ff[2] +
                RV[i][3] * ff[3] + RV[i][4] * ff[4];
  }
  return;
}

// Compute the initial condition
int SetIC(N_Vector y, EulerData& udata)
{
  sunrealtype* rho = N_VGetSubvectorArrayPointer_ManyVector(y, 0);
  if (check_ptr(rho, "N_VGetSubvectorArrayPointer_ManyVector")) { return -1; }
  sunrealtype* mx = N_VGetSubvectorArrayPointer_ManyVector(y, 1);
  if (check_ptr(mx, "N_VGetSubvectorArrayPointer_ManyVector")) { return -1; }
  sunrealtype* my = N_VGetSubvectorArrayPointer_ManyVector(y, 2);
  if (check_ptr(my, "N_VGetSubvectorArrayPointer_ManyVector")) { return -1; }
  sunrealtype* mz = N_VGetSubvectorArrayPointer_ManyVector(y, 3);
  if (check_ptr(mz, "N_VGetSubvectorArrayPointer_ManyVector")) { return -1; }
  sunrealtype* et = N_VGetSubvectorArrayPointer_ManyVector(y, 4);
  if (check_ptr(et, "N_VGetSubvectorArrayPointer_ManyVector")) { return -1; }

  for (long int i = 0; i < udata.nx; i++)
  {
    sunrealtype xloc = ((sunrealtype)i + HALF) * udata.dx + udata.xl;
    if (xloc < HALF)
    {
      rho[i] = rhoL;
      et[i]  = udata.eos_inv(rhoL, uL, ZERO, ZERO, pL);
      mx[i]  = rhoL * uL;
    }
    else
    {
      rho[i] = rhoR;
      et[i]  = udata.eos_inv(rhoR, uR, ZERO, ZERO, pR);
      mx[i]  = rhoR * uR;
    }
    my[i] = ZERO;
    mz[i] = ZERO;
  }

  return 0;
}

//---- end of file ----
