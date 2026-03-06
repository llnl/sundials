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
 * Nested multirate nonlinear Kvaerno-Prothero-Robinson ODE test
 * problem:
 *
 *    [u]' = [ G   e   e ] [(u^2-p-2)/(2u)] +  [ p'(t)/(2u) ]
 *    [v]    [ e  al  be ] [(v^2-q-2)/(2v)]    [ q'(t)/(2v) ]
 *    [w]    [ e -be  al ] [(w^2-r-2)/(2w)]    [ r'(t)/(2w) ]
 *
 * where p(t) = 0.5*cos(t), q(t) = cos(om*t*(1+exp(-(t-2)^2))),
 * and r(t) = cos(om*om*t*(1+exp(-(t-3)^2))).
 *
 * The first row corresponds to the slowest time scale; when an ImEx method is applied
 * to this time scale, we set
 *
 *    fsi = [ G   e   e ] [(u^2-p-2)/(2u)]
 *          [ 0   0   0 ] [(v^2-q-2)/(2v)]
 *          [ 0   0   0 ] [(w^2-r-2)/(2w)]
 *
 *    fse = [ p'(t)/(2u) ]
 *          [      0     ]
 *          [      0     ]
 *
 * The second row corresponds to the intermediate time scale; when an ImEx method
 * is applied to this time scale, we set
 *
 *    fmi = [ 0   0   0 ] [(u^2-p-2)/(2u)]
 *          [ e  al  be ] [(v^2-q-2)/(2v)]
 *          [ 0   0   0 ] [(w^2-r-2)/(2w)]
 *
 *    fme = [      0     ]
 *          [ q'(t)/(2v) ]
 *          [      0     ]
 *
 * The fast RHS is always just the full third row of the IVP.
 *
 * This problem has analytical solution given by
 *    u(t) = sqrt(2+p(t)), v(t) = sqrt(2+q(t)), w(t) = sqrt(2+r(t)).
 * However, we use a reference solver here to assess performance
 * of the local multirate adaptivity controller.
 *
 * This program allows a number of parameters:
 *   G: stiffness at slow time scale [default = -10]
 *   e: fast/slow coupling strength [default = 0.5]
 *   al,be: oscillatory coupling between v and w [default al = -1, be = 1]
 *   om: time-scale separation factor [default = 50]
 *
 * The stiffness of the slow time scale is essentially determined
 * by G, for |G| > 50 it is 'stiff' and ideally suited to a
 * multirate method that is implicit at the slow time scale.
 *
 * Coupling between the slow and faster components is determined by e, with
 * coupling strength proportional to |e|.  Coupling between the intermediate
 * and fast components is determined by (al,be).
 *
 * The "intermediate" variable, v, oscillates at a frequency "om" times
 * faster than u, and the "fast" variable, w, oscillates at a frequency
 * "om" times faster than v.
 *
 * Additional input options may be used to select between various
 * solver options:
 * - slow fixed/initial step size:  hs [default = 0.01]
 * - intermediate fixed/initial step size:  hm [default = 0.001]
 * - fast fixed/initial step size:  hf [default = 0.0001]
 * - set initial adaptive step size as hs/hm/hf above:  set_h0 [default 0]
 * - relative solution tolerance:  rtol [default = 1e-4]
 * - absolute solution tolerance:  atol [default = 1e-11]
 * - relative solution tolerance for each fast integrator (as compared
 *     to next-slower integrator):  fast_rtol [default = 1e-4]
 * - use p (0) vs q (1) for slow and intermediate adaptivity:  slow_pq [default = 0]
 * - use p (0) vs q (1) for fast adaptivity:  fast_pq [default = 0]
 * - slow stepsize safety factor:  safety [default = 0.96]
 * - "slow" MRI method:  mri_method [default = ARKODE_MRI_GARK_ERK45a]
 * - "intermediate" MRI method:  mid_method [default = ARKODE_MRI_GARK_ERK45a]
 * - "fast" ERKStep method order: fast_order [default 4]
 * - "slow" and "intermediate" MRI temporal adaptivity controllers: scontrol [default = 1]
 *      0:  no controller [fixed time steps]
 *      1:  I controller (as part of MRI-HTOL)
 *      2:  PI controller (as part of MRI-HTOL)
 *      3:  PID controller (as part of MRI-HTOL)
 *      4:  ExpGus controller (as part of MRI-HTOL)
 *      5:  ImpGus controller (as part of MRI-HTOL)
 *      6:  ImExGus controller (as part of MRI-HTOL)
 *      7:  I controller (alone)
 *      8:  PI controller (alone)
 *      9:  PID controller (alone)
 *      10: ExpGus controller (alone)
 *      11: ImpGus controller (alone)
 *      12: ImExGus controller (alone)
 * - "fast" ERKStep temporal adaptivity controller: fcontrol [default = 1]
 *      0:  no controller [fixed time steps]
 *      1:  I controller
 *      2:  PI controller
 *      3:  PID controller
 *      4:  ExpGus controller
 *      5:  ImpGus controller
 *      6:  ImExGus controller
 * - "intermediate" and "fast" accumulated error type: faccum [default = 0]
 *     -1:  no accumulation
 *      0:  maximum accumulation
 *      1:  additive accumulation
 *      2:  average accumulation
 * - controller parameters: (k1s, k2s, k3s, k1f, k2f, k3f,
 *                           bias, htol_relch, htol_minfac, htol_maxfac)
 *     slow and intermediate single-rate controllers: use k1s through k3s,
 *     as appropriate.  fast single-rate controllers: use k1f through k3f,
 *     as appropriate.  MRIHTol controllers: use htol_relch, htol_minfac,
 *     htol_maxfac.  all controllers use bias.
 *     ** if any one of a relevant set are "-1" then the defaults are used
 *
 * Outputs and solution error values are printed at equal intervals
 * of 0.5 and run statistics are printed at the end.
 * ----------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <arkode/arkode_erkstep.h>  // prototypes for ERKStep fcts., consts
#include <arkode/arkode_mristep.h>  // prototypes for MRIStep fcts., consts
#include <nvector/nvector_serial.h> // serial N_Vector type, fcts., macros
#include <sunadaptcontroller/sunadaptcontroller_imexgus.h>
#include <sunadaptcontroller/sunadaptcontroller_mrihtol.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>
#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_dense.h> // dense linear solver
#include <sunmatrix/sunmatrix_dense.h> // dense matrix type, fcts., macros

#include <example_utilities.hpp> // common utility functions

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define ESYM "Le"
#define FSYM "Lf"
#else
#define ESYM "e"
#define FSYM "f"
#endif

#define ZERO  SUN_RCONST(0.0)
#define HALF  SUN_RCONST(0.5)
#define ONE   SUN_RCONST(1.0)
#define TWO   SUN_RCONST(2.0)
#define THREE SUN_RCONST(3.0)

// Problem options
struct Options
{
  // Problem parameters
  sunrealtype e  = SUN_RCONST(0.5);
  sunrealtype G  = SUN_RCONST(-10.0);
  sunrealtype om = SUN_RCONST(50.0);
  sunrealtype al = SUN_RCONST(-1.0);
  sunrealtype be = SUN_RCONST(1.0);

  // Step sizes and tolerances
  int set_h0            = 0;
  sunrealtype hs        = SUN_RCONST(1.0e-2);
  sunrealtype hm        = SUN_RCONST(1.0e-3);
  sunrealtype hf        = SUN_RCONST(1.0e-4);
  sunrealtype rtol      = SUN_RCONST(1.0e-4);
  sunrealtype atol      = SUN_RCONST(1.0e-11);
  sunrealtype fast_rtol = SUN_RCONST(1.0e-4);

  // Method selection
  std::string mri_method = "ARKODE_MRI_GARK_ERK45a";
  std::string mid_method = "ARKODE_MRI_GARK_ERK45a";
  int fast_order         = 4;
  int scontrol           = 1;
  int fcontrol           = 1;
  int faccum             = 0;
  int slow_pq            = 0;
  int fast_pq            = 0;

  // controller parameters
  sunrealtype k1s         = NAN;
  sunrealtype k2s         = NAN;
  sunrealtype k3s         = NAN;
  sunrealtype k1f         = NAN;
  sunrealtype k2f         = NAN;
  sunrealtype k3f         = NAN;
  sunrealtype bias        = NAN;
  sunrealtype htol_relch  = NAN;
  sunrealtype htol_minfac = NAN;
  sunrealtype htol_maxfac = NAN;
  sunrealtype slow_safety = NAN;
};

// User-supplied functions called by the solver
static int fse(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fsi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fmi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fme(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fm(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int ff(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fn(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Js(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Jsi(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Jm(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Jmi(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Utility functions
static void InputHelp();
static int ReadInputs(std::vector<std::string>& args, Options& opts,
                      SUNContext ctx);
static void PrintSlowAdaptivity(Options opts);
static void PrintFastAdaptivity(Options opts);
static sunrealtype p(sunrealtype t, const Options& opts);
static sunrealtype q(sunrealtype t, const Options& opts);
static sunrealtype r(sunrealtype t, const Options& opts);
static sunrealtype pdot(sunrealtype t, const Options& opts);
static sunrealtype qdot(sunrealtype t, const Options& opts);
static sunrealtype rdot(sunrealtype t, const Options& opts);
static sunrealtype utrue(sunrealtype t, const Options& opts);
static sunrealtype vtrue(sunrealtype t, const Options& opts);
static sunrealtype wtrue(sunrealtype t, const Options& opts);
static int Ytrue(sunrealtype t, N_Vector y, const Options& opts);

// Main Program
int main(int argc, char* argv[])
{
  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Read input options
  Options opts;
  std::vector<std::string> args(argv + 1, argv + argc);
  int flag = ReadInputs(args, opts, sunctx);
  if (check_flag(flag, "ReadInputs")) return 1;

  // General problem parameters
  sunrealtype T0   = SUN_RCONST(0.0); // initial time
  sunrealtype Tf   = SUN_RCONST(5.0); // final time
  sunindextype NEQ = 3;               // number of dependent vars.
  int Nt           = 20;              // number of output times

  // Initial problem output
  //    While traversing these, set various function pointers, table constants, and method orders.
  ARKRhsFn f_f, f_me, f_mi, f_se, f_si;
  ARKLsJacFn J_s, J_m;
  int retval;
  sunbooleantype slowimplicit, slowimex, midimplicit, midimex;
  slowimplicit = slowimex = midimplicit = midimex = SUNFALSE;

  f_mi = nullptr;
  f_me = fm;
  f_si = nullptr;
  f_se = fs;
  J_m  = nullptr;
  J_s  = nullptr;
  f_f  = ff;

  if ((opts.mri_method == "ARKODE_MRI_GARK_IRK21a") ||
      (opts.mri_method == "ARKODE_MRI_GARK_ESDIRK34a") ||
      (opts.mri_method == "ARKODE_MRI_GARK_ESDIRK46a"))
  {
    slowimplicit = SUNTRUE;
    f_se         = nullptr;
    f_si         = fs;
    J_s          = Js;
  }
  if ((opts.mri_method == "ARKODE_IMEX_MRI_SR21") ||
      (opts.mri_method == "ARKODE_IMEX_MRI_SR32") ||
      (opts.mri_method == "ARKODE_IMEX_MRI_SR43"))
  {
    slowimex     = SUNTRUE;
    slowimplicit = SUNTRUE;
    f_se         = fse;
    f_si         = fsi;
    J_s          = Jsi;
  }
  if ((opts.mid_method == "ARKODE_MRI_GARK_IRK21a") ||
      (opts.mid_method == "ARKODE_MRI_GARK_ESDIRK34a") ||
      (opts.mid_method == "ARKODE_MRI_GARK_ESDIRK46a"))
  {
    midimplicit = SUNTRUE;
    f_me        = nullptr;
    f_mi        = fm;
    J_m         = Jm;
  }
  if ((opts.mid_method == "ARKODE_IMEX_MRI_SR21") ||
      (opts.mid_method == "ARKODE_IMEX_MRI_SR32") ||
      (opts.mid_method == "ARKODE_IMEX_MRI_SR43"))
  {
    midimex     = SUNTRUE;
    midimplicit = SUNTRUE;
    f_me        = fme;
    f_mi        = fmi;
    J_m         = Jmi;
  }
  std::cout
    << "\nAdaptive nested multirate nonlinear Kvaerno-Prothero-Robinson test "
       "problem:\n";
  std::cout << "    time domain:  (" << T0 << "," << Tf << "]\n";
  std::cout << "    G = " << opts.G << std::endl;
  std::cout << "    e = " << opts.e << std::endl;
  std::cout << "    al = " << opts.al << std::endl;
  std::cout << "    be = " << opts.be << std::endl;
  std::cout << "    om = " << opts.om << std::endl;
  std::cout << "\n  Slow integrator: " << opts.mri_method;
  if (slowimex) { std::cout << " (ImEx)" << std::endl; }
  else if (slowimplicit) { std::cout << " (implicit)" << std::endl; }
  else { std::cout << " (explicit)" << std::endl; }
  std::cout << "\n  Intermediate integrator: " << opts.mid_method;
  if (midimex) { std::cout << " (ImEx)" << std::endl; }
  else if (midimplicit) { std::cout << " (implicit)" << std::endl; }
  else { std::cout << " (explicit)" << std::endl; }
  PrintSlowAdaptivity(opts);
  std::cout << "\n  Fast order " << opts.fast_order << std::endl;
  PrintFastAdaptivity(opts);

  // Create and initialize serial vectors for the solution and reference
  N_Vector y = N_VNew_Serial(NEQ, sunctx);
  if (check_ptr((void*)y, "N_VNew_Serial")) return 1;
  N_Vector yref = N_VClone(y);
  if (check_ptr((void*)yref, "N_VClone")) return 1;

  // Set initial conditions
  retval = Ytrue(T0, y, opts);
  if (check_flag(retval, "Ytrue")) return 1;
  N_VScale(ONE, y, yref);

  // Create and configure reference solver object
  void* arkode_ref = ERKStepCreate(fn, T0, yref, sunctx);
  if (check_ptr((void*)arkode_ref, "ERKStepCreate")) return 1;
  retval = ARKodeSetUserData(arkode_ref, (void*)&opts);
  if (check_flag(retval, "ARKodeSetUserData")) return 1;
  retval = ARKodeSetOrder(arkode_ref, 5);
  if (check_flag(retval, "ARKodeSetOrder")) return 1;
  retval = ARKodeSStolerances(arkode_ref, SUN_RCONST(1.e-10), SUN_RCONST(1.e-12));
  if (check_flag(retval, "ARKodeSStolerances")) return 1;
  retval = ARKodeSetMaxNumSteps(arkode_ref, 10000000);
  if (check_flag(retval, "ARKodeSetMaxNumSteps")) return (1);

  // Create and configure fast controller object
  SUNAdaptController fcontrol = nullptr;
  switch (opts.fcontrol)
  {
  case (1):
    fcontrol = SUNAdaptController_I(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_I")) return 1;
    if (!std::isnan(opts.k1f))
    {
      retval = SUNAdaptController_SetParams_I(fcontrol, opts.k1f);
      if (check_flag(retval, "SUNAdaptController_SetParams_I")) return 1;
    }
    break;
  case (2):
    fcontrol = SUNAdaptController_PI(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_PI")) return 1;
    if (!(std::isnan(opts.k1f) || std::isnan(opts.k2f)))
    {
      retval = SUNAdaptController_SetParams_PI(fcontrol, opts.k1f, opts.k2f);
      if (check_flag(retval, "SUNAdaptController_SetParams_PI")) return 1;
    }
    break;
  case (3):
    fcontrol = SUNAdaptController_PID(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_PID")) return 1;
    if (!(std::isnan(opts.k1f) || std::isnan(opts.k2f) || std::isnan(opts.k3f)))
    {
      retval = SUNAdaptController_SetParams_PID(fcontrol, opts.k1f, opts.k2f,
                                                opts.k3f);
      if (check_flag(retval, "SUNAdaptController_SetParams_PID")) return 1;
    }
    break;
  case (4):
    fcontrol = SUNAdaptController_ExpGus(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_ExpGus")) return 1;
    if (!(std::isnan(opts.k1f) || std::isnan(opts.k2f)))
    {
      retval = SUNAdaptController_SetParams_ExpGus(fcontrol, opts.k1f, opts.k2f);
      if (check_flag(retval, "SUNAdaptController_SetParams_ExpGus")) return 1;
    }
    break;
  case (5):
    fcontrol = SUNAdaptController_ImpGus(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_ImpGus")) return 1;
    if (!(std::isnan(opts.k1f) || std::isnan(opts.k2f)))
    {
      retval = SUNAdaptController_SetParams_ImpGus(fcontrol, opts.k1f, opts.k2f);
      if (check_flag(retval, "SUNAdaptController_SetParams_ImpGus")) return 1;
    }
    break;
  case (6):
    fcontrol = SUNAdaptController_ImExGus(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_ImExGus")) return 1;
    break;
  }
  if (!std::isnan(opts.bias) && (opts.fcontrol > 0))
  {
    retval = SUNAdaptController_SetErrorBias(fcontrol, opts.bias);
    if (check_flag(retval, "SUNAdaptController_SetErrorBias")) return 1;
  }

  // Create ERKStep (fast) integrator
  void* inner_arkode_mem = ERKStepCreate(f_f, T0, y, sunctx);
  if (check_ptr((void*)inner_arkode_mem, "ERKStepCreate")) return 1;
  retval = ARKodeSetOrder(inner_arkode_mem, opts.fast_order);
  if (check_flag(retval, "ARKodeSetOrder")) return 1;
  retval = ARKodeSStolerances(inner_arkode_mem, opts.fast_rtol, opts.atol);
  if (check_flag(retval, "ARKodeSStolerances")) return 1;
  if (opts.fcontrol != 0)
  {
    retval = ARKodeSetAdaptController(inner_arkode_mem, fcontrol);
    if (check_flag(retval, "ARKodeSetAdaptController")) return 1;
    if (opts.set_h0 != 0)
    {
      retval = ARKodeSetInitStep(inner_arkode_mem, opts.hf);
      if (check_flag(retval, "ARKodeSetInitStep")) return 1;
    }
    if (opts.fast_pq == 1)
    {
      retval = ARKodeSetAdaptivityAdjustment(inner_arkode_mem, 0);
      if (check_flag(retval, "ARKodeSetAdaptivityAdjustment")) return 1;
    }
  }
  else
  {
    retval = ARKodeSetFixedStep(inner_arkode_mem, opts.hf);
    if (check_flag(retval, "ARKodeSetFixedStep")) return 1;
  }
  ARKAccumError acc_type = ARK_ACCUMERROR_NONE;
  if (opts.faccum == 0) { acc_type = ARK_ACCUMERROR_MAX; }
  if (opts.faccum == 1) { acc_type = ARK_ACCUMERROR_SUM; }
  if (opts.faccum == 2) { acc_type = ARK_ACCUMERROR_AVG; }
  retval = ARKodeSetAccumulatedErrorType(inner_arkode_mem, acc_type);
  if (check_flag(retval, "ARKodeSetAccumulatedErrorType")) return 1;
  retval = ARKodeSetMaxNumSteps(inner_arkode_mem, 1000000);
  if (check_flag(retval, "ARKodeSetMaxNumSteps")) return 1;
  retval = ARKodeSetUserData(inner_arkode_mem, (void*)&opts);
  if (check_flag(retval, "ARKodeSetUserData")) return 1;

  // Create inner stepper
  MRIStepInnerStepper inner_stepper = nullptr;
  retval = ARKodeCreateMRIStepInnerStepper(inner_arkode_mem, &inner_stepper);
  if (check_flag(retval, "ARKodeCreateMRIStepInnerStepper")) return 1;

  // Create intermediate and slow controller objects, and select orders of accuracy as relevant
  SUNAdaptController scontrol     = nullptr;
  SUNAdaptController scontrol_H   = nullptr;
  SUNAdaptController scontrol_Tol = nullptr;
  SUNAdaptController mcontrol     = nullptr;
  SUNAdaptController mcontrol_H   = nullptr;
  SUNAdaptController mcontrol_Tol = nullptr;
  switch (opts.scontrol)
  {
  case (1):
  {
    scontrol_H = SUNAdaptController_I(sunctx);
    if (check_ptr((void*)scontrol_H, "SUNAdaptController_I (slow H)")) return 1;
    scontrol_Tol = SUNAdaptController_I(sunctx);
    if (check_ptr((void*)scontrol_Tol, "SUNAdaptController_I (slow Tol)"))
      return 1;
    mcontrol_H = SUNAdaptController_I(sunctx);
    if (check_ptr((void*)mcontrol_H, "SUNAdaptController_I (mid H)")) return 1;
    mcontrol_Tol = SUNAdaptController_I(sunctx);
    if (check_ptr((void*)mcontrol_Tol, "SUNAdaptController_I (mid Tol)"))
      return 1;
    if (!std::isnan(opts.k1s))
    {
      retval = SUNAdaptController_SetParams_I(scontrol_H, opts.k1s);
      if (check_flag(retval, "SUNAdaptController_SetParams_I")) return 1;
      retval = SUNAdaptController_SetParams_I(scontrol_Tol, opts.k1s);
      if (check_flag(retval, "SUNAdaptController_SetParams_I")) return 1;
      retval = SUNAdaptController_SetParams_I(mcontrol_H, opts.k1s);
      if (check_flag(retval, "SUNAdaptController_SetParams_I")) return 1;
      retval = SUNAdaptController_SetParams_I(mcontrol_Tol, opts.k1s);
      if (check_flag(retval, "SUNAdaptController_SetParams_I")) return 1;
    }
    scontrol = SUNAdaptController_MRIHTol(scontrol_H, scontrol_Tol, sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    mcontrol = SUNAdaptController_MRIHTol(mcontrol_H, mcontrol_Tol, sunctx);
    if (check_ptr((void*)mcontrol, "SUNAdaptController_MRIHTol")) return 1;
    if (!(std::isnan(opts.htol_relch) || std::isnan(opts.htol_minfac) ||
          std::isnan(opts.htol_maxfac)))
    {
      retval = SUNAdaptController_SetParams_MRIHTol(scontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
      retval = SUNAdaptController_SetParams_MRIHTol(mcontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
    }
    break;
  }
  case (2):
  {
    scontrol_H = SUNAdaptController_PI(sunctx);
    if (check_ptr((void*)scontrol_H, "SUNAdaptController_PI (slow H)"))
      return 1;
    scontrol_Tol = SUNAdaptController_PI(sunctx);
    if (check_ptr((void*)scontrol_Tol, "SUNAdaptController_PI (slow Tol)"))
      return 1;
    mcontrol_H = SUNAdaptController_PI(sunctx);
    if (check_ptr((void*)mcontrol_H, "SUNAdaptController_PI (mid H)")) return 1;
    mcontrol_Tol = SUNAdaptController_PI(sunctx);
    if (check_ptr((void*)mcontrol_Tol, "SUNAdaptController_PI (mid Tol)"))
      return 1;
    if (!(std::isnan(opts.k1s) || std::isnan(opts.k2s)))
    {
      retval = SUNAdaptController_SetParams_PI(scontrol_H, opts.k1s, opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_PI")) return 1;
      retval = SUNAdaptController_SetParams_PI(scontrol_Tol, opts.k1s, opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_PI")) return 1;
      retval = SUNAdaptController_SetParams_PI(mcontrol_H, opts.k1s, opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_PI")) return 1;
      retval = SUNAdaptController_SetParams_PI(mcontrol_Tol, opts.k1s, opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_PI")) return 1;
    }
    scontrol = SUNAdaptController_MRIHTol(scontrol_H, scontrol_Tol, sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    mcontrol = SUNAdaptController_MRIHTol(mcontrol_H, mcontrol_Tol, sunctx);
    if (check_ptr((void*)mcontrol, "SUNAdaptController_MRIHTol")) return 1;
    if (!(std::isnan(opts.htol_relch) || std::isnan(opts.htol_minfac) ||
          std::isnan(opts.htol_maxfac)))
    {
      retval = SUNAdaptController_SetParams_MRIHTol(scontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
      retval = SUNAdaptController_SetParams_MRIHTol(mcontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
    }
    break;
  }
  case (3):
  {
    scontrol_H = SUNAdaptController_PID(sunctx);
    if (check_ptr((void*)scontrol_H, "SUNAdaptController_PID (slow H)"))
      return 1;
    scontrol_Tol = SUNAdaptController_PID(sunctx);
    if (check_ptr((void*)scontrol_Tol, "SUNAdaptController_PID (slow Tol)"))
      return 1;
    mcontrol_H = SUNAdaptController_PID(sunctx);
    if (check_ptr((void*)mcontrol_H, "SUNAdaptController_PID (mid H)"))
      return 1;
    mcontrol_Tol = SUNAdaptController_PID(sunctx);
    if (check_ptr((void*)mcontrol_Tol, "SUNAdaptController_PID (mid Tol)"))
      return 1;
    if (!(std::isnan(opts.k1s) || std::isnan(opts.k2s) || std::isnan(opts.k3s)))
    {
      retval = SUNAdaptController_SetParams_PID(scontrol_H, opts.k1s, opts.k2s,
                                                opts.k3s);
      if (check_flag(retval, "SUNAdaptController_SetParams_PID")) return 1;
      retval = SUNAdaptController_SetParams_PID(scontrol_Tol, opts.k1s,
                                                opts.k2s, opts.k3s);
      if (check_flag(retval, "SUNAdaptController_SetParams_PID")) return 1;
      retval = SUNAdaptController_SetParams_PID(mcontrol_H, opts.k1s, opts.k2s,
                                                opts.k3s);
      if (check_flag(retval, "SUNAdaptController_SetParams_PID")) return 1;
      retval = SUNAdaptController_SetParams_PID(mcontrol_Tol, opts.k1s,
                                                opts.k2s, opts.k3s);
      if (check_flag(retval, "SUNAdaptController_SetParams_PID")) return 1;
    }
    scontrol = SUNAdaptController_MRIHTol(scontrol_H, scontrol_Tol, sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    mcontrol = SUNAdaptController_MRIHTol(mcontrol_H, mcontrol_Tol, sunctx);
    if (check_ptr((void*)mcontrol, "SUNAdaptController_MRIHTol")) return 1;
    if (!(std::isnan(opts.htol_relch) || std::isnan(opts.htol_minfac) ||
          std::isnan(opts.htol_maxfac)))
    {
      retval = SUNAdaptController_SetParams_MRIHTol(scontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
      retval = SUNAdaptController_SetParams_MRIHTol(mcontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
    }
    break;
  }
  case (4):
  {
    scontrol_H = SUNAdaptController_ExpGus(sunctx);
    if (check_ptr((void*)scontrol_H, "SUNAdaptController_ExpGus (slow H)"))
      return 1;
    scontrol_Tol = SUNAdaptController_ExpGus(sunctx);
    if (check_ptr((void*)scontrol_Tol, "SUNAdaptController_ExpGus (slow Tol)"))
      return 1;
    mcontrol_H = SUNAdaptController_ExpGus(sunctx);
    if (check_ptr((void*)mcontrol_H, "SUNAdaptController_ExpGus (mid H)"))
      return 1;
    mcontrol_Tol = SUNAdaptController_ExpGus(sunctx);
    if (check_ptr((void*)mcontrol_Tol, "SUNAdaptController_ExpGus (mid Tol)"))
      return 1;
    if (!(std::isnan(opts.k1s) || std::isnan(opts.k2s)))
    {
      retval = SUNAdaptController_SetParams_ExpGus(scontrol_H, opts.k1s,
                                                   opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_ExpGus")) return 1;
      retval = SUNAdaptController_SetParams_ExpGus(scontrol_Tol, opts.k1s,
                                                   opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_ExpGus")) return 1;
      retval = SUNAdaptController_SetParams_ExpGus(mcontrol_H, opts.k1s,
                                                   opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_ExpGus")) return 1;
      retval = SUNAdaptController_SetParams_ExpGus(mcontrol_Tol, opts.k1s,
                                                   opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_ExpGus")) return 1;
    }
    scontrol = SUNAdaptController_MRIHTol(scontrol_H, scontrol_Tol, sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    mcontrol = SUNAdaptController_MRIHTol(mcontrol_H, mcontrol_Tol, sunctx);
    if (check_ptr((void*)mcontrol, "SUNAdaptController_MRIHTol")) return 1;
    if (!(std::isnan(opts.htol_relch) || std::isnan(opts.htol_minfac) ||
          std::isnan(opts.htol_maxfac)))
    {
      retval = SUNAdaptController_SetParams_MRIHTol(scontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
      retval = SUNAdaptController_SetParams_MRIHTol(mcontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
    }
    break;
  }
  case (5):
  {
    scontrol_H = SUNAdaptController_ImpGus(sunctx);
    if (check_ptr((void*)scontrol_H, "SUNAdaptController_ImpGus (slow H)"))
      return 1;
    scontrol_Tol = SUNAdaptController_ImpGus(sunctx);
    if (check_ptr((void*)scontrol_Tol, "SUNAdaptController_ImpGus (slow Tol)"))
      return 1;
    mcontrol_H = SUNAdaptController_ImpGus(sunctx);
    if (check_ptr((void*)mcontrol_H, "SUNAdaptController_ImpGus (mid H)"))
      return 1;
    mcontrol_Tol = SUNAdaptController_ImpGus(sunctx);
    if (check_ptr((void*)mcontrol_Tol, "SUNAdaptController_ImpGus (mid Tol)"))
      return 1;
    if (!(std::isnan(opts.k1s) || std::isnan(opts.k2s)))
    {
      retval = SUNAdaptController_SetParams_ImpGus(scontrol_H, opts.k1s,
                                                   opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_ImpGus")) return 1;
      retval = SUNAdaptController_SetParams_ImpGus(scontrol_Tol, opts.k1s,
                                                   opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_ImpGus")) return 1;
      retval = SUNAdaptController_SetParams_ImpGus(mcontrol_H, opts.k1s,
                                                   opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_ImpGus")) return 1;
      retval = SUNAdaptController_SetParams_ImpGus(mcontrol_Tol, opts.k1s,
                                                   opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_ImpGus")) return 1;
    }
    scontrol = SUNAdaptController_MRIHTol(scontrol_H, scontrol_Tol, sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    mcontrol = SUNAdaptController_MRIHTol(mcontrol_H, mcontrol_Tol, sunctx);
    if (check_ptr((void*)mcontrol, "SUNAdaptController_MRIHTol")) return 1;
    if (!(std::isnan(opts.htol_relch) || std::isnan(opts.htol_minfac) ||
          std::isnan(opts.htol_maxfac)))
    {
      retval = SUNAdaptController_SetParams_MRIHTol(scontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
      retval = SUNAdaptController_SetParams_MRIHTol(mcontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
    }
    break;
  }
  case (6):
  {
    scontrol_H = SUNAdaptController_ImExGus(sunctx);
    if (check_ptr((void*)scontrol_H, "SUNAdaptController_ImExGus (slow H)"))
      return 1;
    scontrol_Tol = SUNAdaptController_ImExGus(sunctx);
    if (check_ptr((void*)scontrol_Tol, "SUNAdaptController_ImExGus (slow Tol)"))
      return 1;
    mcontrol_H = SUNAdaptController_ImExGus(sunctx);
    if (check_ptr((void*)mcontrol_H, "SUNAdaptController_ImExGus (mid H)"))
      return 1;
    mcontrol_Tol = SUNAdaptController_ImExGus(sunctx);
    if (check_ptr((void*)mcontrol_Tol, "SUNAdaptController_ImExGus (mid Tol)"))
      return 1;
    scontrol = SUNAdaptController_MRIHTol(scontrol_H, scontrol_Tol, sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    mcontrol = SUNAdaptController_MRIHTol(scontrol_H, scontrol_Tol, sunctx);
    if (check_ptr((void*)mcontrol, "SUNAdaptController_MRIHTol")) return 1;
    if (!(std::isnan(opts.htol_relch) || std::isnan(opts.htol_minfac) ||
          std::isnan(opts.htol_maxfac)))
    {
      retval = SUNAdaptController_SetParams_MRIHTol(scontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
      retval = SUNAdaptController_SetParams_MRIHTol(mcontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
    }
    break;
  }
  case (7):
  {
    scontrol = SUNAdaptController_I(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptControllerI (slow)")) return 1;
    mcontrol = SUNAdaptController_I(sunctx);
    if (check_ptr((void*)mcontrol, "SUNAdaptControllerI (mid)")) return 1;
    if (!std::isnan(opts.k1s))
    {
      retval = SUNAdaptController_SetParams_I(scontrol, opts.k1s);
      if (check_flag(retval, "SUNAdaptController_SetParams_I")) return 1;
      retval = SUNAdaptController_SetParams_I(mcontrol, opts.k1s);
      if (check_flag(retval, "SUNAdaptController_SetParams_I")) return 1;
    }
    break;
  }
  case (8):
  {
    scontrol = SUNAdaptController_PI(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_PI (slow)")) return 1;
    mcontrol = SUNAdaptController_PI(sunctx);
    if (check_ptr((void*)mcontrol, "SUNAdaptController_PI (mid)")) return 1;
    if (!(std::isnan(opts.k1s) || std::isnan(opts.k2s)))
    {
      retval = SUNAdaptController_SetParams_PI(scontrol, opts.k1s, opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_PI")) return 1;
      retval = SUNAdaptController_SetParams_PI(mcontrol, opts.k1s, opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_PI")) return 1;
    }
    break;
  }
  case (9):
  {
    scontrol = SUNAdaptController_PID(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_PID (slow)")) return 1;
    mcontrol = SUNAdaptController_PID(sunctx);
    if (check_ptr((void*)mcontrol, "SUNAdaptController_PID (mid)")) return 1;
    if (!(std::isnan(opts.k1s) || std::isnan(opts.k2s) || std::isnan(opts.k3s)))
    {
      retval = SUNAdaptController_SetParams_PID(scontrol, opts.k1s, opts.k2s,
                                                opts.k3s);
      if (check_flag(retval, "SUNAdaptController_SetParams_PID")) return 1;
      retval = SUNAdaptController_SetParams_PID(mcontrol, opts.k1s, opts.k2s,
                                                opts.k3s);
      if (check_flag(retval, "SUNAdaptController_SetParams_PID")) return 1;
    }
    break;
  }
  case (10):
  {
    scontrol = SUNAdaptController_ExpGus(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_ExpGus (slow)"))
      return 1;
    mcontrol = SUNAdaptController_ExpGus(sunctx);
    if (check_ptr((void*)mcontrol, "SUNAdaptController_ExpGus (mid)")) return 1;
    if (!(std::isnan(opts.k1s) || std::isnan(opts.k2s)))
    {
      retval = SUNAdaptController_SetParams_ExpGus(scontrol, opts.k1s, opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_ExpGus")) return 1;
      retval = SUNAdaptController_SetParams_ExpGus(mcontrol, opts.k1s, opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_ExpGus")) return 1;
    }
    break;
  }
  case (11):
  {
    scontrol = SUNAdaptController_ImpGus(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_ImpGus (slow)"))
      return 1;
    mcontrol = SUNAdaptController_ImpGus(sunctx);
    if (check_ptr((void*)mcontrol, "SUNAdaptController_ImpGus (mid)")) return 1;
    if (!(std::isnan(opts.k1s) || std::isnan(opts.k2s)))
    {
      retval = SUNAdaptController_SetParams_ImpGus(scontrol, opts.k1s, opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_ImpGus")) return 1;
      retval = SUNAdaptController_SetParams_ImpGus(mcontrol, opts.k1s, opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_ImpGus")) return 1;
    }
    break;
  }
  case (12):
  {
    scontrol = SUNAdaptController_ImExGus(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_ImExGus (slow)"))
      return 1;
    mcontrol = SUNAdaptController_ImExGus(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_ImExGus (mid)"))
      return 1;
    break;
  }
  }
  if (!std::isnan(opts.bias) && (opts.scontrol > 0))
  {
    retval = SUNAdaptController_SetErrorBias(scontrol, opts.bias);
    if (check_flag(retval, "SUNAdaptController_SetErrorBias")) return 1;
    retval = SUNAdaptController_SetErrorBias(mcontrol, opts.bias);
    if (check_flag(retval, "SUNAdaptController_SetErrorBias")) return 1;
  }

  // Create MRI (intermediate) integrator
  void* mid_arkode_mem = MRIStepCreate(f_me, f_mi, T0, y, inner_stepper, sunctx);
  if (check_ptr((void*)mid_arkode_mem, "MRIStepCreate")) return 1;
  MRIStepCoupling Cm = MRIStepCoupling_LoadTableByName((opts.mid_method).c_str());
  if (check_ptr((void*)Cm, "MRIStepCoupling_LoadTableByName")) return 1;
  retval = MRIStepSetCoupling(mid_arkode_mem, Cm);
  if (check_flag(retval, "MRIStepSetCoupling")) return 1;
  SUNMatrix Am        = nullptr; // matrix for intermediate solver
  SUNLinearSolver LSm = nullptr; // intermediate linear solver object
  if (midimplicit)
  {
    Am = SUNDenseMatrix(NEQ, NEQ, sunctx);
    if (check_ptr((void*)Am, "SUNDenseMatrix")) return 1;
    LSm = SUNLinSol_Dense(y, Am, sunctx);
    if (check_ptr((void*)LSm, "SUNLinSol_Dense")) return 1;
    retval = ARKodeSetLinearSolver(mid_arkode_mem, LSm, Am);
    if (check_flag(retval, "ARKodeSetLinearSolver")) return 1;
    retval = ARKodeSetJacFn(mid_arkode_mem, J_m);
    if (check_flag(retval, "ARKodeSetJacFn")) return 1;
  }
  retval = ARKodeSStolerances(mid_arkode_mem, opts.rtol, opts.atol);
  if (check_flag(retval, "ARKodeSStolerances")) return 1;
  retval = ARKodeSetMaxNumSteps(mid_arkode_mem, 100000);
  if (check_flag(retval, "ARKodeSetMaxNumSteps")) return 1;
  retval = ARKodeSetAccumulatedErrorType(mid_arkode_mem, acc_type);
  if (check_flag(retval, "ARKodeSetAccumulatedErrorType")) return 1;
  retval = ARKodeSetUserData(mid_arkode_mem, (void*)&opts);
  if (check_flag(retval, "ARKodeSetUserData")) return 1;
  if (opts.scontrol != 0)
  {
    retval = ARKodeSetAdaptController(mid_arkode_mem, mcontrol);
    if (check_flag(retval, "ARKodeSetAdaptController")) return 1;
    if (opts.set_h0 != 0)
    {
      retval = ARKodeSetInitStep(mid_arkode_mem, opts.hm);
      if (check_flag(retval, "ARKodeSetInitStep")) return 1;
    }
    if (opts.slow_pq == 1)
    {
      retval = ARKodeSetAdaptivityAdjustment(mid_arkode_mem, 0);
      if (check_flag(retval, "ARKodeSetAdaptivityAdjustment")) return 1;
    }
    if (!std::isnan(opts.slow_safety))
    {
      retval = ARKodeSetSafetyFactor(mid_arkode_mem, opts.slow_safety);
      if (check_flag(retval, "ARKodeSetSafetyFactor")) return 1;
    }
  }
  else
  {
    retval = ARKodeSetFixedStep(mid_arkode_mem, opts.hm);
    if (check_flag(retval, "ARKodeSetFixedStep")) return 1;
  }

  // Create intermediate stepper
  MRIStepInnerStepper intermediate_stepper = nullptr;
  retval = ARKodeCreateMRIStepInnerStepper(mid_arkode_mem, &intermediate_stepper);
  if (check_flag(retval, "ARKodeCreateMRIStepInnerStepper")) return 1;

  // Create MRI (slow) integrator
  void* arkode_mem = MRIStepCreate(f_se, f_si, T0, y, intermediate_stepper,
                                   sunctx);
  if (check_ptr((void*)arkode_mem, "MRIStepCreate")) return 1;
  MRIStepCoupling Cs = MRIStepCoupling_LoadTableByName((opts.mri_method).c_str());
  if (check_ptr((void*)Cs, "MRIStepCoupling_LoadTableByName")) return 1;
  retval = MRIStepSetCoupling(arkode_mem, Cs);
  if (check_flag(retval, "MRIStepSetCoupling")) return 1;
  SUNMatrix As        = nullptr; // matrix for slow solver
  SUNLinearSolver LSs = nullptr; // slow linear solver object
  if (slowimplicit)
  {
    As = SUNDenseMatrix(NEQ, NEQ, sunctx);
    if (check_ptr((void*)As, "SUNDenseMatrix")) return 1;
    LSs = SUNLinSol_Dense(y, As, sunctx);
    if (check_ptr((void*)LSs, "SUNLinSol_Dense")) return 1;
    retval = ARKodeSetLinearSolver(arkode_mem, LSs, As);
    if (check_flag(retval, "ARKodeSetLinearSolver")) return 1;
    retval = ARKodeSetJacFn(arkode_mem, J_s);
    if (check_flag(retval, "ARKodeSetJacFn")) return 1;
  }
  retval = ARKodeSStolerances(arkode_mem, opts.rtol, opts.atol);
  if (check_flag(retval, "ARKodeSStolerances")) return 1;
  retval = ARKodeSetMaxNumSteps(arkode_mem, 100000);
  if (check_flag(retval, "ARKodeSetMaxNumSteps")) return 1;
  retval = ARKodeSetUserData(arkode_mem, (void*)&opts);
  if (check_flag(retval, "ARKodeSetUserData")) return 1;
  if (opts.scontrol != 0)
  {
    retval = ARKodeSetAdaptController(arkode_mem, scontrol);
    if (check_flag(retval, "ARKodeSetAdaptController")) return 1;
    if (opts.set_h0 != 0)
    {
      retval = ARKodeSetInitStep(arkode_mem, opts.hs);
      if (check_flag(retval, "ARKodeSetInitStep")) return 1;
    }
    if (opts.slow_pq == 1)
    {
      retval = ARKodeSetAdaptivityAdjustment(arkode_mem, 0);
      if (check_flag(retval, "ARKodeSetAdaptivityAdjustment")) return 1;
    }
    if (!std::isnan(opts.slow_safety))
    {
      retval = ARKodeSetSafetyFactor(arkode_mem, opts.slow_safety);
      if (check_flag(retval, "ARKodeSetSafetyFactor")) return 1;
    }
  }
  else
  {
    retval = ARKodeSetFixedStep(arkode_mem, opts.hs);
    if (check_flag(retval, "ARKodeSetFixedStep")) return 1;
  }

  //
  // Integrate ODE
  //

  // Main time-stepping loop: calls ARKodeEvolve to perform the
  // integration, then prints results. Stops when the final time
  // has been reached
  sunrealtype t         = T0;
  sunrealtype t2        = T0;
  sunrealtype dTout     = (Tf - T0) / Nt;
  sunrealtype tout      = T0 + dTout;
  sunrealtype* ydata    = N_VGetArrayPointer(y);
  sunrealtype* yrefdata = N_VGetArrayPointer(yref);
  sunrealtype u, v, w, uerr, verr, werr, uerrtot, verrtot, werrtot, errtot,
    accuracy;
  uerr = verr = werr = uerrtot = verrtot = werrtot = errtot = accuracy = ZERO;
  printf("        t          u          v          w        uerr       verr    "
         "   werr\n");
  printf("   "
         "---------------------------------------------------------------------"
         "-------\n");
  printf("  %10.6" FSYM " %10.6" FSYM " %10.6" FSYM " %10.6" FSYM "   %.2" ESYM
         "   %.2" ESYM "   %.2" ESYM "\n",
         t, ydata[0], ydata[1], ydata[2], uerr, verr, werr);
  while (Tf - t > SUN_RCONST(1.0e-8))
  {
    // reset reference solver so that it begins with identical state
    retval = ARKodeReset(arkode_ref, t, y);

    // evolve solution in one-step mode
    retval = ARKodeSetStopTime(arkode_mem, tout);
    if (check_flag(retval, "ARKodeSetStopTime")) return 1;
    retval = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_ONE_STEP);
    if (retval < 0)
    {
      printf("ARKodeEvolve error (%i)\n", retval);
      return 1;
    }

    // evolve reference solver to same time in "normal" mode
    retval = ARKodeSetStopTime(arkode_ref, t);
    if (check_flag(retval, "ARKodeSetStopTime")) return 1;
    retval = ARKodeEvolve(arkode_ref, t, yref, &t2, ARK_NORMAL);
    if (retval < 0)
    {
      printf("ARKodeEvolve reference solution error (%i)\n", retval);
      return 1;
    }

    // access/print solution and error
    u    = ydata[0];
    v    = ydata[1];
    w    = ydata[2];
    uerr = std::abs(yrefdata[0] - u);
    verr = std::abs(yrefdata[1] - v);
    werr = std::abs(yrefdata[2] - w);
    uerrtot += uerr * uerr;
    verrtot += verr * verr;
    werrtot += werr * werr;
    errtot += uerr * uerr + verr * verr + werr * werr;
    accuracy = std::max(accuracy,
                        uerr / std::abs(opts.atol + opts.rtol * yrefdata[0]));
    accuracy = std::max(accuracy,
                        verr / std::abs(opts.atol + opts.rtol * yrefdata[1]));
    accuracy = std::max(accuracy,
                        werr / std::abs(opts.atol + opts.rtol * yrefdata[2]));

    // Periodically output current results to screen
    if (t >= tout)
    {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
      printf("  %10.6" FSYM " %10.6" FSYM " %10.6" FSYM " %10.6" FSYM
             "   %.2" ESYM "   %.2" ESYM "   %.2" ESYM "\n",
             t, u, v, w, uerr, verr, werr);
    }
  }
  printf("   "
         "---------------------------------------------------------------------"
         "-------\n");

  //
  // Finalize
  //

  // Get some slow integrator statistics
  long int nsts, natts, netfs, nfse, nfsi, nifs;
  retval = ARKodeGetNumSteps(arkode_mem, &nsts);
  check_flag(retval, "ARKodeGetNumSteps");
  retval = ARKodeGetNumStepAttempts(arkode_mem, &natts);
  check_flag(retval, "ARKodeGetNumStepAttempts");
  retval = ARKodeGetNumErrTestFails(arkode_mem, &netfs);
  check_flag(retval, "ARKodeGetNumErrTestFails");
  retval = ARKodeGetNumRhsEvals(arkode_mem, 0, &nfse);
  check_flag(retval, "ARKodeGetNumRhsEvals");
  retval = ARKodeGetNumRhsEvals(arkode_mem, 1, &nfsi);
  check_flag(retval, "ARKodeGetNumRhsEvals");
  retval = MRIStepGetNumInnerStepperFails(arkode_mem, &nifs);
  check_flag(retval, "MRIStepGetNumInnerStepperFails");

  // Get some intermediate integrator statistics
  long int nstm, nattm, netfm, nfme, nfmi, nifm;
  retval = ARKodeGetNumSteps(mid_arkode_mem, &nstm);
  check_flag(retval, "ARKodeGetNumSteps");
  retval = ARKodeGetNumStepAttempts(mid_arkode_mem, &nattm);
  check_flag(retval, "ARKodeGetNumStepAttempts");
  retval = ARKodeGetNumErrTestFails(mid_arkode_mem, &netfm);
  check_flag(retval, "ARKodeGetNumErrTestFails");
  retval = ARKodeGetNumRhsEvals(mid_arkode_mem, 0, &nfme);
  check_flag(retval, "ARKodeGetNumRhsEvals");
  retval = ARKodeGetNumRhsEvals(mid_arkode_mem, 1, &nfmi);
  check_flag(retval, "ARKodeGetNumRhsEvals");
  retval = MRIStepGetNumInnerStepperFails(mid_arkode_mem, &nifm);
  check_flag(retval, "MRIStepGetNumInnerStepperFails");

  // Get some fast integrator statistics
  long int nstf, nattf, netff, nff;
  retval = ARKodeGetNumSteps(inner_arkode_mem, &nstf);
  check_flag(retval, "ARKodeGetNumSteps");
  retval = ARKodeGetNumStepAttempts(inner_arkode_mem, &nattf);
  check_flag(retval, "ARKodeGetNumStepAttempts");
  retval = ARKodeGetNumErrTestFails(inner_arkode_mem, &netff);
  check_flag(retval, "ARKodeGetNumErrTestFails");
  retval = ARKodeGetNumRhsEvals(inner_arkode_mem, 0, &nff);
  check_flag(retval, "ARKodeGetNumRhsEvals");

  // Print some final statistics
  uerrtot = std::sqrt(uerrtot / (sunrealtype)nsts);
  verrtot = std::sqrt(verrtot / (sunrealtype)nsts);
  werrtot = std::sqrt(werrtot / (sunrealtype)nsts);
  errtot  = std::sqrt(errtot / SUN_RCONST(3.0) / (sunrealtype)nsts);
  std::cout << "\nFinal Solver Statistics:\n";
  std::cout << "   Slow steps = " << nsts << "  (attempts = " << natts
            << ",  fails = " << netfs << ",  innerfails = " << nifs << ")\n";
  std::cout << "   Intermediate steps = " << nstm << "  (attempts = " << nattm
            << ",  fails = " << netfm << ",  innerfails = " << nifm << ")\n";
  std::cout << "   Fast steps = " << nstf << "  (attempts = " << nattf
            << ",  fails = " << netff << ")\n";
  std::cout << "   u error = " << uerrtot << ", v error = " << verrtot
            << ", total error = " << errtot << std::endl;
  std::cout << "   Relative accuracy = " << accuracy << std::endl;
  std::cout << "   Total RHS evals:  Fse = " << nfse << ", Fsi = " << nfsi
            << ", Fme = " << nfme << ", Fmi = " << nfmi << ", Ff = " << nff
            << std::endl;

  // Get/print slow integrator implicit solver statistics
  if (slowimplicit)
  {
    long int nnis, nncs, njes;
    retval = ARKodeGetNonlinSolvStats(arkode_mem, &nnis, &nncs);
    check_flag(retval, "ARKodeGetNonlinSolvStats");
    retval = ARKodeGetNumJacEvals(arkode_mem, &njes);
    check_flag(retval, "ARKodeGetNumJacEvals");
    std::cout << "   Slow Newton iters = " << nnis << std::endl;
    std::cout << "   Slow Newton iters/attempt = "
              << (sunrealtype)nnis / (sunrealtype)natts << std::endl;
    std::cout << "   Slow Newton conv fails = " << nncs << std::endl;
    std::cout << "   Slow Jacobian evals = " << njes << std::endl;
    std::cout << "   Slow Jacobian evals/Newton = "
              << (sunrealtype)njes / (sunrealtype)nnis << std::endl;
  }

  // Get/print intermediate integrator implicit solver statistics
  if (midimplicit)
  {
    long int nnim, nncm, njem;
    retval = ARKodeGetNonlinSolvStats(mid_arkode_mem, &nnim, &nncm);
    check_flag(retval, "ARKodeGetNonlinSolvStats");
    retval = ARKodeGetNumJacEvals(mid_arkode_mem, &njem);
    check_flag(retval, "ARKodeGetNumJacEvals");
    std::cout << "   Intermediate Newton iters = " << nnim << std::endl;
    std::cout << "   Intermediate Newton iters/attempt = "
              << (sunrealtype)nnim / (sunrealtype)nattm << std::endl;
    std::cout << "   Intermediate Newton conv fails = " << nncm << std::endl;
    std::cout << "   Intermediate Jacobian evals = " << njem << std::endl;
    std::cout << "   Intermediate Jacobian evals/Newton = "
              << (sunrealtype)njem / (sunrealtype)nnim << std::endl;
  }

  // Clean up and return
  N_VDestroy(y);
  N_VDestroy(yref);
  MRIStepCoupling_Free(Cs);
  MRIStepCoupling_Free(Cm);
  SUNMatDestroy(As);
  SUNLinSolFree(LSs);
  SUNMatDestroy(Am);
  SUNLinSolFree(LSm);
  SUNAdaptController_Destroy(scontrol);
  SUNAdaptController_Destroy(scontrol_H);
  SUNAdaptController_Destroy(scontrol_Tol);
  SUNAdaptController_Destroy(mcontrol);
  SUNAdaptController_Destroy(mcontrol_H);
  SUNAdaptController_Destroy(mcontrol_Tol);
  SUNAdaptController_Destroy(fcontrol);
  ARKodeFree(&inner_arkode_mem); // Free fast integrator memory
  ARKodeFree(&mid_arkode_mem);   // Free intermediate integrator memory
  MRIStepInnerStepper_Free(&inner_stepper); // Free inner stepper structures
  MRIStepInnerStepper_Free(&intermediate_stepper);
  ARKodeFree(&arkode_mem); // Free slow integrator memory
  ARKodeFree(&arkode_ref); // Free reference solver memory

  return 0;
}

// ------------------------------
// Functions called by the solver
// -----------------------------

// fn routine to compute the full ODE RHS.
static int fn(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Options* opts         = static_cast<Options*>(user_data);
  sunrealtype* ydata    = N_VGetArrayPointer(y);
  sunrealtype* ydotdata = N_VGetArrayPointer(ydot);
  const sunrealtype u   = ydata[0];
  const sunrealtype v   = ydata[1];
  const sunrealtype w   = ydata[2];
  const sunrealtype G   = opts->G;
  const sunrealtype e   = opts->e;
  const sunrealtype al  = opts->al;
  const sunrealtype be  = opts->be;
  sunrealtype tmp1, tmp2, tmp3;

  // fill in the RHS function:
  //   [ G   e   e ] [(u^2-p-2)/(2u)] +  [ pdot(t)/(2u) ]
  //   [ e  al  be ] [(v^2-q-2)/(2v)]    [ qdot(t)/(2v) ]
  //   [ e -be  al ] [(w^2-r-2)/(2w)]    [ rdot(t)/(2w) ]
  tmp1        = (-TWO + u * u - p(t, *opts)) / (TWO * u);
  tmp2        = (-TWO + v * v - q(t, *opts)) / (TWO * v);
  tmp3        = (-TWO + w * w - r(t, *opts)) / (TWO * w);
  ydotdata[0] = G * tmp1 + e * tmp2 + e * tmp3 + pdot(t, *opts) / (TWO * u);
  ydotdata[1] = e * tmp1 + al * tmp2 + be * tmp3 + qdot(t, *opts) / (TWO * v);
  ydotdata[2] = e * tmp1 - be * tmp2 + al * tmp3 + rdot(t, *opts) / (TWO * w);

  // Return with success
  return 0;
}

// ff routine to compute the fast portion of the ODE RHS.
static int ff(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Options* opts         = static_cast<Options*>(user_data);
  sunrealtype* ydata    = N_VGetArrayPointer(y);
  sunrealtype* ydotdata = N_VGetArrayPointer(ydot);
  const sunrealtype u   = ydata[0];
  const sunrealtype v   = ydata[1];
  const sunrealtype w   = ydata[2];
  const sunrealtype e   = opts->e;
  const sunrealtype al  = opts->al;
  const sunrealtype be  = opts->be;
  sunrealtype tmp1, tmp2, tmp3;

  // fill in the RHS function:
  //   [ 0   0   0 ] [(u^2-p-2)/(2u)] +  [      0       ]
  //   [ 0   0   0 ] [(v^2-q-2)/(2v)]    [      0       ]
  //   [ e -be  al ] [(w^2-r-2)/(2w)]    [ rdot(t)/(2w) ]
  tmp1        = (-TWO + u * u - p(t, *opts)) / (TWO * u);
  tmp2        = (-TWO + v * v - q(t, *opts)) / (TWO * v);
  tmp3        = (-TWO + w * w - r(t, *opts)) / (TWO * w);
  ydotdata[0] = ZERO;
  ydotdata[1] = ZERO;
  ydotdata[2] = e * tmp1 - be * tmp2 + al * tmp3 + rdot(t, *opts) / (TWO * w);

  // Return with success
  return 0;
}

// fm routine to compute the intermediate portion of the ODE RHS.
static int fm(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Options* opts         = static_cast<Options*>(user_data);
  sunrealtype* ydata    = N_VGetArrayPointer(y);
  sunrealtype* ydotdata = N_VGetArrayPointer(ydot);
  const sunrealtype u   = ydata[0];
  const sunrealtype v   = ydata[1];
  const sunrealtype w   = ydata[2];
  const sunrealtype e   = opts->e;
  const sunrealtype al  = opts->al;
  const sunrealtype be  = opts->be;
  sunrealtype tmp1, tmp2, tmp3;

  // fill in the RHS function:
  //   [ 0   0   0 ] [(u^2-p-2)/(2u)] +  [      0       ]
  //   [ e  al  be ] [(v^2-q-2)/(2v)]    [ qdot(t)/(2v) ]
  //   [ 0   0   0 ] [(w^2-r-2)/(2w)]    [      0       ]
  tmp1        = (-TWO + u * u - p(t, *opts)) / (TWO * u);
  tmp2        = (-TWO + v * v - q(t, *opts)) / (TWO * v);
  tmp3        = (-TWO + w * w - r(t, *opts)) / (TWO * w);
  ydotdata[0] = ZERO;
  ydotdata[1] = e * tmp1 + al * tmp2 + be * tmp3 + qdot(t, *opts) / (TWO * v);
  ydotdata[2] = ZERO;

  return 0;
}

// fme routine to compute the explicit intermediate portion of the ODE RHS.
static int fme(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Options* opts         = static_cast<Options*>(user_data);
  sunrealtype* ydata    = N_VGetArrayPointer(y);
  sunrealtype* ydotdata = N_VGetArrayPointer(ydot);
  const sunrealtype v   = ydata[1];

  // fill in the RHS function:
  //   [      0       ]
  //   [ qdot(t)/(2v) ]
  //   [      0       ]
  ydotdata[0] = ZERO;
  ydotdata[1] = qdot(t, *opts) / (TWO * v);
  ydotdata[2] = ZERO;

  return 0;
}

// fmi routine to compute the implicit intermediate portion of the ODE RHS.
static int fmi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Options* opts         = static_cast<Options*>(user_data);
  sunrealtype* ydata    = N_VGetArrayPointer(y);
  sunrealtype* ydotdata = N_VGetArrayPointer(ydot);
  const sunrealtype u   = ydata[0];
  const sunrealtype v   = ydata[1];
  const sunrealtype w   = ydata[2];
  const sunrealtype e   = opts->e;
  const sunrealtype al  = opts->al;
  const sunrealtype be  = opts->be;
  sunrealtype tmp1, tmp2, tmp3;

  // fill in the RHS function:
  //   [ 0   0   0 ] [(u^2-p-2)/(2u)]
  //   [ e  al  be ] [(v^2-q-2)/(2v)]
  //   [ 0   0   0 ] [(w^2-r-2)/(2w)]
  tmp1        = (-TWO + u * u - p(t, *opts)) / (TWO * u);
  tmp2        = (-TWO + v * v - q(t, *opts)) / (TWO * v);
  tmp3        = (-TWO + w * w - r(t, *opts)) / (TWO * w);
  ydotdata[0] = ZERO;
  ydotdata[1] = e * tmp1 + al * tmp2 + be * tmp3;
  ydotdata[2] = ZERO;

  // Return with success
  return 0;
}

// fs routine to compute the slow portion of the ODE RHS.
static int fs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Options* opts         = static_cast<Options*>(user_data);
  sunrealtype* ydata    = N_VGetArrayPointer(y);
  sunrealtype* ydotdata = N_VGetArrayPointer(ydot);
  const sunrealtype u   = ydata[0];
  const sunrealtype v   = ydata[1];
  const sunrealtype w   = ydata[2];
  const sunrealtype G   = opts->G;
  const sunrealtype e   = opts->e;
  sunrealtype tmp1, tmp2, tmp3;

  // fill in the RHS function:
  //   [ G   e   e ] [(u^2-p-2)/(2u)] +  [ pdot(t)/(2u) ]
  //   [ 0   0   0 ] [(v^2-q-2)/(2v)]    [      0       ]
  //   [ 0   0   0 ] [(w^2-r-2)/(2w)]    [      0       ]
  tmp1        = (-TWO + u * u - p(t, *opts)) / (TWO * u);
  tmp2        = (-TWO + v * v - q(t, *opts)) / (TWO * v);
  tmp3        = (-TWO + w * w - r(t, *opts)) / (TWO * w);
  ydotdata[0] = G * tmp1 + e * tmp2 + e * tmp3 + pdot(t, *opts) / (TWO * u);
  ydotdata[1] = ZERO;
  ydotdata[2] = ZERO;

  return 0;
}

// fse routine to compute the explicit slow portion of the ODE RHS.
static int fse(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Options* opts         = static_cast<Options*>(user_data);
  sunrealtype* ydata    = N_VGetArrayPointer(y);
  sunrealtype* ydotdata = N_VGetArrayPointer(ydot);
  const sunrealtype u   = ydata[0];

  // fill in the RHS function:
  //   [ pdot(t)/(2u) ]
  //   [      0       ]
  //   [      0       ]
  ydotdata[0] = pdot(t, *opts) / (TWO * u);
  ydotdata[1] = ZERO;
  ydotdata[2] = ZERO;

  return 0;
}

// fsi routine to compute the implicit slow portion of the ODE RHS.
static int fsi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Options* opts         = static_cast<Options*>(user_data);
  sunrealtype* ydata    = N_VGetArrayPointer(y);
  sunrealtype* ydotdata = N_VGetArrayPointer(ydot);
  const sunrealtype u   = ydata[0];
  const sunrealtype v   = ydata[1];
  const sunrealtype w   = ydata[2];
  const sunrealtype G   = opts->G;
  const sunrealtype e   = opts->e;
  sunrealtype tmp1, tmp2, tmp3;

  // fill in the RHS function:
  //   [ G   e   e ] [(u^2-p-2)/(2u)]
  //   [ 0   0   0 ] [(v^2-q-2)/(2v)]
  //   [ 0   0   0 ] [(w^2-r-2)/(2w)]
  tmp1        = (-TWO + u * u - p(t, *opts)) / (TWO * u);
  tmp2        = (-TWO + v * v - q(t, *opts)) / (TWO * v);
  tmp3        = (-TWO + w * w - r(t, *opts)) / (TWO * w);
  ydotdata[0] = G * tmp1 + e * tmp2 + e * tmp3;
  ydotdata[1] = ZERO;
  ydotdata[2] = ZERO;

  // Return with success
  return 0;
}

// Jacobian of fm
static int Jm(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  Options* opts       = static_cast<Options*>(user_data);
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];
  const sunrealtype w = ydata[2];
  sunrealtype t11, t22, t33;

  // fill in the Jacobian:
  //   [0   0   0]*[1-(u^2-p(t)-2)/(2*u^2),  0,  0] + [0,        0,         0]
  //   [e  al  be] [0,  1-(v^2-q(t)-2)/(2*v^2),  0]   [0,  -q'(t)/(2*v^2),  0]
  //   [0   0   0] [0,  0,  1-(w^2-r(t)-2)/(2*w^2)]   [0,        0,         0]
  t11                   = ONE - (u * u - p(t, *opts) - TWO) / (TWO * u * u);
  t22                   = ONE - (v * v - q(t, *opts) - TWO) / (TWO * v * v);
  t33                   = ONE - (w * w - r(t, *opts) - TWO) / (TWO * w * w);
  SM_ELEMENT_D(J, 0, 0) = ZERO;
  SM_ELEMENT_D(J, 0, 1) = ZERO;
  SM_ELEMENT_D(J, 0, 2) = ZERO;
  SM_ELEMENT_D(J, 1, 0) = opts->e * t11;
  SM_ELEMENT_D(J, 1, 1) = opts->al * t22 - qdot(t, *opts) / (TWO * v * v);
  SM_ELEMENT_D(J, 1, 2) = opts->be * t33;
  SM_ELEMENT_D(J, 2, 0) = ZERO;
  SM_ELEMENT_D(J, 2, 1) = ZERO;
  SM_ELEMENT_D(J, 2, 2) = ZERO;

  // Return with success
  return 0;
}

// Jacobian of fmi
static int Jmi(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  Options* opts       = static_cast<Options*>(user_data);
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];
  const sunrealtype w = ydata[2];
  sunrealtype t11, t22, t33;

  // fill in the Jacobian:
  //   [0   0   0]*[1-(u^2-p(t)-2)/(2*u^2),  0,  0]
  //   [e  al  be] [0,  1-(v^2-q(t)-2)/(2*v^2),  0]
  //   [0   0   0] [0,  0,  1-(w^2-r(t)-2)/(2*w^2)]
  t11                   = ONE - (u * u - p(t, *opts) - TWO) / (TWO * u * u);
  t22                   = ONE - (v * v - q(t, *opts) - TWO) / (TWO * v * v);
  t33                   = ONE - (w * w - r(t, *opts) - TWO) / (TWO * w * w);
  SM_ELEMENT_D(J, 0, 0) = ZERO;
  SM_ELEMENT_D(J, 0, 1) = ZERO;
  SM_ELEMENT_D(J, 0, 2) = ZERO;
  SM_ELEMENT_D(J, 1, 0) = opts->e * t11;
  SM_ELEMENT_D(J, 1, 1) = opts->al * t22;
  SM_ELEMENT_D(J, 1, 2) = opts->be * t33;
  SM_ELEMENT_D(J, 2, 0) = ZERO;
  SM_ELEMENT_D(J, 2, 1) = ZERO;
  SM_ELEMENT_D(J, 2, 2) = ZERO;

  // Return with success
  return 0;
}

// Jacobian of fs
static int Js(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  Options* opts       = static_cast<Options*>(user_data);
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];
  const sunrealtype w = ydata[2];
  sunrealtype t11, t22, t33;

  // fill in the Jacobian:
  //   [G   e   e]*[1-(u^2-p(t)-2)/(2*u^2),  0,  0] + [-p'(t)/(2*u^2),  0,  0]
  //   [0   0   0] [0,  1-(v^2-q(t)-2)/(2*v^2),  0]   [0,               0,  0]
  //   [0   0   0] [0,  0,  1-(w^2-r(t)-2)/(2*w^2)]   [0,               0,  0]
  t11                   = ONE - (u * u - p(t, *opts) - TWO) / (TWO * u * u);
  t22                   = ONE - (v * v - q(t, *opts) - TWO) / (TWO * v * v);
  t33                   = ONE - (w * w - r(t, *opts) - TWO) / (TWO * w * w);
  SM_ELEMENT_D(J, 0, 0) = opts->G * t11 - pdot(t, *opts) / (TWO * u * u);
  SM_ELEMENT_D(J, 0, 1) = opts->e * t22;
  SM_ELEMENT_D(J, 0, 2) = opts->e * t33;
  SM_ELEMENT_D(J, 1, 0) = ZERO;
  SM_ELEMENT_D(J, 1, 1) = ZERO;
  SM_ELEMENT_D(J, 1, 2) = ZERO;
  SM_ELEMENT_D(J, 2, 0) = ZERO;
  SM_ELEMENT_D(J, 2, 1) = ZERO;
  SM_ELEMENT_D(J, 2, 2) = ZERO;

  // Return with success
  return 0;
}

// Jacobian of fsi
static int Jsi(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  Options* opts       = static_cast<Options*>(user_data);
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];
  const sunrealtype w = ydata[2];
  sunrealtype t11, t22, t33;

  // fill in the Jacobian:
  //   [G   e   e]*[1-(u^2-p(t)-2)/(2*u^2),  0,  0]
  //   [0   0   0] [0,  1-(v^2-q(t)-2)/(2*v^2),  0]
  //   [0   0   0] [0,  0,  1-(w^2-r(t)-2)/(2*w^2)]
  t11                   = ONE - (u * u - p(t, *opts) - TWO) / (TWO * u * u);
  t22                   = ONE - (v * v - q(t, *opts) - TWO) / (TWO * v * v);
  t33                   = ONE - (w * w - r(t, *opts) - TWO) / (TWO * w * w);
  SM_ELEMENT_D(J, 0, 0) = opts->G * t11;
  SM_ELEMENT_D(J, 0, 1) = opts->e * t22;
  SM_ELEMENT_D(J, 0, 2) = opts->e * t33;
  SM_ELEMENT_D(J, 1, 0) = ZERO;
  SM_ELEMENT_D(J, 1, 1) = ZERO;
  SM_ELEMENT_D(J, 1, 2) = ZERO;
  SM_ELEMENT_D(J, 2, 0) = ZERO;
  SM_ELEMENT_D(J, 2, 1) = ZERO;
  SM_ELEMENT_D(J, 2, 2) = ZERO;

  // Return with success
  return 0;
}

// ------------------------------
// Private helper functions
// -----------------------------

static sunrealtype p(sunrealtype t, const Options& opts)
{
  return HALF * cos(t);
}

static sunrealtype q(sunrealtype t, const Options& opts)
{
  return (cos(opts.om * t * (ONE + exp(-(t - TWO) * (t - TWO)))));
}

static sunrealtype r(sunrealtype t, const Options& opts)
{
  return (cos(opts.om * opts.om * t * (ONE + exp(-(t - THREE) * (t - THREE)))));
}

static sunrealtype pdot(sunrealtype t, const Options& opts)
{
  return -HALF * sin(t);
}

static sunrealtype qdot(sunrealtype t, const Options& opts)
{
  const sunrealtype tTwo  = t - TWO;
  const sunrealtype eterm = exp(-tTwo * tTwo);
  return (-sin(opts.om * t * (ONE + eterm)) * opts.om *
          (ONE + eterm * (ONE - TWO * t * tTwo)));
}

static sunrealtype rdot(sunrealtype t, const Options& opts)
{
  const sunrealtype tThree = t - THREE;
  const sunrealtype eterm  = exp(-tThree * tThree);
  return (-sin(opts.om * opts.om * t * (ONE + eterm)) * opts.om * opts.om *
          (ONE + eterm * (ONE - TWO * t * tThree)));
}

static sunrealtype utrue(sunrealtype t, const Options& opts)
{
  return (std::sqrt(TWO + p(t, opts)));
}

static sunrealtype vtrue(sunrealtype t, const Options& opts)
{
  return (std::sqrt(TWO + q(t, opts)));
}

static sunrealtype wtrue(sunrealtype t, const Options& opts)
{
  return (std::sqrt(TWO + r(t, opts)));
}

static int Ytrue(sunrealtype t, N_Vector y, const Options& opts)
{
  NV_Ith_S(y, 0) = utrue(t, opts);
  NV_Ith_S(y, 1) = vtrue(t, opts);
  NV_Ith_S(y, 2) = wtrue(t, opts);
  return (0);
}

// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Print command line options
void InputHelp()
{
  std::cout << std::endl;
  std::cout << "Command line options:" << std::endl;
  std::cout << "  --help         : print options and exit\n";
  std::cout << "  --G            : stiffness at slow time scale\n";
  std::cout << "  --e            : fast/slow coupling strength\n";
  std::cout << "  --al           : med/fast oscillation term\n";
  std::cout << "  --be           : med/fast oscillation term\n";
  std::cout << "  --om           : time-scale separation factor\n";
  std::cout << "  --hs           : slow (fixed/initial) step size\n";
  std::cout << "  --hm           : intermediate (fixed/initial) step size\n";
  std::cout << "  --hf           : fast (fixed/initial) step size\n";
  std::cout
    << "  --set_h0       : use hs/hf above to set the initial step size\n";
  std::cout << "  --rtol         : relative solution tolerance\n";
  std::cout << "  --atol         : absolute solution tolerance\n";
  std::cout
    << "  --fast_rtol    : relative solution tolerance for fast method\n";
  std::cout << "  --safety : slow time step safety factor\n";
  std::cout
    << "  --mri_method   : slow MRI method name (valid ARKODE_MRITableID)\n";
  std::cout << "  --mid_method   : intermediate MRI method name (valid "
               "ARKODE_MRITableID)\n";
  std::cout << "  --fast_order   : fast RK method order\n";
  std::cout << "  --scontrol     : slow/intermediate time step controllers, "
               "int in [0,12] "
               "(see source)\n";
  std::cout << "  --fcontrol     : fast time step controller, int in [0,6] "
               "(see source)\n";
  std::cout << "  --faccum       : fast error accumulation type {-1,0,1,2}\n";
  std::cout << "  --slow_pq      : use p (0) vs q (1) for slow/intermediate "
               "adaptivity\n";
  std::cout << "  --fast_pq      : use p (0) vs q (1) for fast adaptivity\n";
  std::cout
    << "  --k1s, --k2s, -k3s : slow/intermediate controller parameters\n";
  std::cout << "  --k1f, --k2f, -k3f : fast controller parameters\n";
  std::cout << "  --bias : slow and fast controller bias factors\n";
  std::cout
    << "  --htol_relch : HTol controller maximum relative tolerance change\n";
  std::cout
    << "  --htol_minfac : HTol controller minimum relative tolerance factor\n";
  std::cout
    << "  --htol_maxfac : HTol controller maximum relative tolerance factor\n";
}

// Read input options
int ReadInputs(std::vector<std::string>& args, Options& opts, SUNContext ctx)
{
  if (find(args.begin(), args.end(), "--help") != args.end())
  {
    InputHelp();
    return 1;
  }

  // Problem options
  find_arg(args, "--G", opts.G);
  find_arg(args, "--e", opts.e);
  find_arg(args, "--al", opts.al);
  find_arg(args, "--be", opts.be);
  find_arg(args, "--om", opts.om);
  find_arg(args, "--hs", opts.hs);
  find_arg(args, "--hm", opts.hm);
  find_arg(args, "--hf", opts.hf);
  find_arg(args, "--set_h0", opts.set_h0);
  find_arg(args, "--rtol", opts.rtol);
  find_arg(args, "--atol", opts.atol);
  find_arg(args, "--fast_rtol", opts.fast_rtol);
  find_arg(args, "--safety", opts.slow_safety);
  find_arg(args, "--mri_method", opts.mri_method);
  find_arg(args, "--mid_method", opts.mid_method);
  find_arg(args, "--fast_order", opts.fast_order);
  find_arg(args, "--scontrol", opts.scontrol);
  find_arg(args, "--fcontrol", opts.fcontrol);
  find_arg(args, "--faccum", opts.faccum);
  find_arg(args, "--slow_pq", opts.slow_pq);
  find_arg(args, "--fast_pq", opts.fast_pq);
  find_arg(args, "--k1s", opts.k1s);
  find_arg(args, "--k2s", opts.k2s);
  find_arg(args, "--k3s", opts.k3s);
  find_arg(args, "--k1f", opts.k1f);
  find_arg(args, "--k2f", opts.k2f);
  find_arg(args, "--k3f", opts.k3f);
  find_arg(args, "--bias", opts.bias);
  find_arg(args, "--htol_relch", opts.htol_relch);
  find_arg(args, "--htol_minfac", opts.htol_minfac);
  find_arg(args, "--htol_maxfac", opts.htol_maxfac);

  // Check inputs for validity
  //   0 < rtol < 1
  if ((opts.rtol <= ZERO) || (opts.rtol >= ONE))
  {
    std::cerr << "ERROR: rtol must be in (0,1), (" << opts.rtol << " input)\n";
    return -1;
  }
  //   0 < atol < 1
  if ((opts.atol <= ZERO) || (opts.atol >= ONE))
  {
    std::cerr << "ERROR: atol must be in (0,1), (" << opts.atol << " input)\n";
    return -1;
  }
  //   0 < fast_rtol < 1
  if ((opts.fast_rtol <= ZERO) || (opts.fast_rtol >= ONE))
  {
    std::cerr << "ERROR: fast_rtol must be in (0,1), (" << opts.fast_rtol
              << " input)\n";
    return -1;
  }
  //   0 < fast_order
  if (opts.fast_order <= 0)
  {
    std::cerr << "ERROR: fast_order must be at least 1, (" << opts.fast_order
              << " input)\n";
    return -1;
  }
  //   slow_pq in {0,1}
  if ((opts.slow_pq < 0) || (opts.slow_pq > 1))
  {
    std::cerr << "ERROR: slow_pq must be in {0,1}, (" << opts.slow_pq
              << " input)\n";
    return -1;
  }
  //   fast_pq in {0,1}
  if ((opts.fast_pq < 0) || (opts.fast_pq > 1))
  {
    std::cerr << "ERROR: fast_pq must be in {0,1}, (" << opts.fast_pq
              << " input)\n";
    return -1;
  }
  //   scontrol in [0,12]
  if ((opts.scontrol < 0) || (opts.scontrol > 12))
  {
    std::cerr << "ERROR: scontrol must be in [0,12], (" << opts.scontrol
              << " input)\n";
    return -1;
  }
  //   fcontrol in [0,6]
  if ((opts.fcontrol < 0) || (opts.fcontrol > 6))
  {
    std::cerr << "ERROR: fcontrol must be in [0,6], (" << opts.fcontrol
              << " input)\n";
    return -1;
  }
  //   hs > 0 and hm > 0 if scontrol == 0
  if (((opts.hs <= 0) || (opts.hm <= 0)) && (opts.scontrol == 0))
  {
    std::cerr << "ERROR: positive hs required with scontrol = 0, (" << opts.hs
              << " input)\n";
    return -1;
  }
  //   hf > 0 if fcontrol == 0
  if ((opts.hf <= 0) && (opts.fcontrol == 0))
  {
    std::cerr << "ERROR: positive hf required with fcontrol = 0, (" << opts.hf
              << " input)\n";
    return -1;
  }
  //   G < 0.0
  if (opts.G >= ZERO)
  {
    std::cerr << "ERROR: G must be a negative real number, (" << opts.G
              << " input)\n";
    return -1;
  }
  //   om >= 1.0
  if (opts.om < ONE)
  {
    std::cerr << "ERROR: w must be >= 1.0, (" << opts.om << " input)\n";
    return -1;
  }

  return 0;
}

static void PrintSlowAdaptivity(Options opts)
{
  switch (opts.scontrol)
  {
  case (0):
  {
    std::cout << "    fixed steps, hs = " << opts.hs << ", hm = " << opts.hm
              << std::endl;
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  }
  case (1):
  {
    std::cout
      << "    MRI-HTOL controller (using I for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    std::cout << "    fast error accumulation strategy = " << opts.faccum << "\n";
    if (!std::isnan(opts.k1s))
    {
      std::cout << "    slow/intermediate controller parameter: " << opts.k1s
                << "\n";
    }
    if (!(std::isnan(opts.htol_relch) || std::isnan(opts.htol_minfac) ||
          std::isnan(opts.htol_maxfac)))
    {
      std::cout << "    HTol controller parameters: " << opts.htol_relch << " "
                << opts.htol_minfac << " " << opts.htol_maxfac << "\n";
    }
    break;
  }
  case (2):
  {
    std::cout
      << "    MRI-HTOL controller (using PI for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    std::cout << "    fast error accumulation strategy = " << opts.faccum << "\n";
    if (!(std::isnan(opts.k1s) || std::isnan(opts.k2s)))
    {
      std::cout << "    slow/intermediate controller parameters: " << opts.k1s
                << " " << opts.k2s << "\n";
    }
    if (!(std::isnan(opts.htol_relch) || std::isnan(opts.htol_minfac) ||
          std::isnan(opts.htol_maxfac)))
    {
      std::cout << "    HTol controller parameters: " << opts.htol_relch << " "
                << opts.htol_minfac << " " << opts.htol_maxfac << "\n";
    }
    break;
  }
  case (3):
  {
    std::cout
      << "    MRI-HTOL controller (using PID for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    std::cout << "    fast error accumulation strategy = " << opts.faccum << "\n";
    if (!(std::isnan(opts.k1s) || std::isnan(opts.k2s) || std::isnan(opts.k3s)))
    {
      std::cout << "    slow/intermediate controller parameters: " << opts.k1s
                << " " << opts.k2s << " " << opts.k3s << "\n";
    }
    if (!(std::isnan(opts.htol_relch) || std::isnan(opts.htol_minfac) ||
          std::isnan(opts.htol_maxfac)))
    {
      std::cout << "    HTol controller parameters: " << opts.htol_relch << " "
                << opts.htol_minfac << " " << opts.htol_maxfac << "\n";
    }
    break;
  }
  case (4):
  {
    std::cout
      << "    MRI-HTOL controller (using ExpGus for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    std::cout << "    fast error accumulation strategy = " << opts.faccum << "\n";
    if (!(std::isnan(opts.k1s) || std::isnan(opts.k2s)))
    {
      std::cout << "    slow/intermediate controller parameters: " << opts.k1s
                << " " << opts.k2s << "\n";
    }
    if (!(std::isnan(opts.htol_relch) || std::isnan(opts.htol_minfac) ||
          std::isnan(opts.htol_maxfac)))
    {
      std::cout << "    HTol controller parameters: " << opts.htol_relch << " "
                << opts.htol_minfac << " " << opts.htol_maxfac << "\n";
    }
    break;
  }
  case (5):
  {
    std::cout
      << "    MRI-HTOL controller (using ImpGus for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    std::cout << "    fast error accumulation strategy = " << opts.faccum << "\n";
    if (!(std::isnan(opts.k1s) || std::isnan(opts.k2s)))
    {
      std::cout << "    slow/intermediate controller parameters: " << opts.k1s
                << " " << opts.k2s << "\n";
    }
    if (!(std::isnan(opts.htol_relch) || std::isnan(opts.htol_minfac) ||
          std::isnan(opts.htol_maxfac)))
    {
      std::cout << "    HTol controller parameters: " << opts.htol_relch << " "
                << opts.htol_minfac << " " << opts.htol_maxfac << "\n";
    }
    break;
  }
  case (6):
  {
    std::cout
      << "    MRI-HTOL controller (using ImExGus for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    std::cout << "    fast error accumulation strategy = " << opts.faccum << "\n";
    break;
  }
  case (7):
  {
    std::cout << "    Decoupled I controller for slow time scale, based on "
                 "order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    if (!std::isnan(opts.k1s))
    {
      std::cout << "    slow/intermediate controller parameter: " << opts.k1s
                << "\n";
    }
    break;
  }
  case (8):
  {
    std::cout << "    Decoupled PI controller for slow time scale, based on "
                 "order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    if (!(std::isnan(opts.k1s) || std::isnan(opts.k2s)))
    {
      std::cout << "    slow/intermediate controller parameters: " << opts.k1s
                << " " << opts.k2s << "\n";
    }
    break;
  }
  case (9):
  {
    std::cout << "    Decoupled PID controller for slow time scale, based on "
                 "order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    if (!(std::isnan(opts.k1s) || std::isnan(opts.k2s) || std::isnan(opts.k3s)))
    {
      std::cout << "    slow/intermediate controller parameters: " << opts.k1s
                << " " << opts.k2s << " " << opts.k3s << "\n";
    }
    break;
  }
  case (10):
  {
    std::cout << "    Decoupled ExpGus controller for slow time scale, based "
                 "on order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    if (!(std::isnan(opts.k1s) || std::isnan(opts.k2s)))
    {
      std::cout << "    slow/intermediate controller parameters: " << opts.k1s
                << " " << opts.k2s << "\n";
    }
    break;
  }
  case (11):
  {
    std::cout << "    Decoupled ImpGus controller for slow time scale, based "
                 "on order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    if (!(std::isnan(opts.k1s) || std::isnan(opts.k2s)))
    {
      std::cout << "    slow/intermediate controller parameters: " << opts.k1s
                << " " << opts.k2s << "\n";
    }
    break;
  }
  case (12):
  {
    std::cout << "    Decoupled ImExGus controller for slow time scale, based "
                 "on order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  }
  }
  if (!std::isnan(opts.bias))
  {
    std::cout << "    controller bias factor: " << opts.bias << "\n";
  }
  if (!std::isnan(opts.slow_safety))
  {
    std::cout << "    slow step safety factor: " << opts.slow_safety << "\n";
  }
}

static void PrintFastAdaptivity(Options opts)
{
  switch (opts.fcontrol)
  {
  case (0):
  {
    std::cout << "    fixed steps, hf = " << opts.hf << std::endl;
    std::cout << "    fast_rtol = " << opts.fast_rtol
              << ", atol = " << opts.atol << "\n";
    break;
  }
  case (1):
  {
    std::cout << "    I controller for fast time scale, based on order of RK "
              << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    fast_rtol = " << opts.fast_rtol
              << ", atol = " << opts.atol << "\n";
    if (!std::isnan(opts.k1f))
    {
      std::cout << "    fast controller parameter: " << opts.k1f << "\n";
    }
    break;
  }
  case (2):
  {
    std::cout << "    PI controller for fast time scale, based on order of RK "
              << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    fast_rtol = " << opts.fast_rtol
              << ", atol = " << opts.atol << "\n";
    if (!(std::isnan(opts.k1f) || std::isnan(opts.k2f)))
    {
      std::cout << "    fast controller parameters: " << opts.k1f << " "
                << opts.k2f << "\n";
    }
    break;
  }
  case (3):
  {
    std::cout << "    PID controller for fast time scale, based on order of RK "
              << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    fast_rtol = " << opts.fast_rtol
              << ", atol = " << opts.atol << "\n";
    if (!(std::isnan(opts.k1f) || std::isnan(opts.k2f) || std::isnan(opts.k3f)))
    {
      std::cout << "    fast controller parameters: " << opts.k1f << " "
                << opts.k2f << " " << opts.k3f << "\n";
    }
    break;
  }
  case (4):
  {
    std::cout
      << "    ExpGus controller for fast time scale, based on order of RK "
      << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    fast_rtol = " << opts.fast_rtol
              << ", atol = " << opts.atol << "\n";
    if (!(std::isnan(opts.k1f) || std::isnan(opts.k2f)))
    {
      std::cout << "    fast controller parameters: " << opts.k1f << " "
                << opts.k2f << "\n";
    }
    break;
  }
  case (5):
  {
    std::cout
      << "    ImpGus controller for fast time scale, based on order of RK "
      << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    fast_rtol = " << opts.fast_rtol
              << ", atol = " << opts.atol << "\n";
    if (!(std::isnan(opts.k1f) || std::isnan(opts.k2f)))
    {
      std::cout << "    fast controller parameters: " << opts.k1f << " "
                << opts.k2f << "\n";
    }
    break;
  }
  case (6):
  {
    std::cout
      << "    ImExGus controller for fast time scale, based on order of RK "
      << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    fast_rtol = " << opts.fast_rtol
              << ", atol = " << opts.atol << "\n";
    break;
  }
  }
}

//---- end of file ----//
