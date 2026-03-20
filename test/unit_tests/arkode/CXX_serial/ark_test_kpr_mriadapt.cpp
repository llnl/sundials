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
 * Multirate nonlinear Kvaerno-Prothero-Robinson ODE test problem:
 *
 *    [u]' = [ G  e ] [(u^2-r-1)/(2u)] +  [ r'(t)/(2u) ]
 *    [v]    [ e -1 ] [(v^2-s-2)/(2v)]    [ s'(t)/(2v) ]
 *         = [ fs(t,u,v) ]
 *           [ ff(t,u,v) ]
 *
 * where r(t) = cos(t), and s(t) = cos(w*t*(1+exp(-(t-2)^2))).
 *
 * This problem has analytical solution given by
 *    u(t) = sqrt(2+r(t)),  v(t) = sqrt(2+s(t)).
 * However, we use a reference solver here to assess performance
 * of the local multirate adaptivity controller.
 *
 * This program allows a number of parameters:
 *   e: fast/slow coupling strength [default = 0.5]
 *   G: stiffness at slow time scale [default = -1e2]
 *   w: time-scale separation factor [default = 100]
 *
 * The stiffness of the slow time scale is essentially determined
 * by G, for |G| > 50 it is 'stiff' and ideally suited to a
 * multirate method that is implicit at the slow time scale.
 *
 * Coupling between the two components is determined by e, with
 * coupling strength proportional to |e|.
 *
 * The "fast" variable, v, oscillates at a frequency "w" times
 * faster than u.
 *
 * Additional input options may be used to select between various
 * solver options:
 * - slow fixed/initial step size:  hs [default = 0.01]
 * - fast fixed/initial step size:  hf [default = 0.0001]
 * - set initial adaptive step size as hs/hf above:  set_h0 [default 0]
 * - relative solution tolerance:  rtol [default = 1e-4]
 * - absolute solution tolerance:  atol [default = 1e-11]
 * - relative solution tolerance for fast integrator:  fast_rtol [default = 1e-4]
 * - use p (0) vs q (1) for slow adaptivity:  slow_pq [default = 0]
 * - use p (0) vs q (1) for fast adaptivity:  fast_pq [default = 0]
 * - slow stepsize safety factor:  safety [default = 0.96]
 * - "slow" MRI method:  mri_method [default = ARKODE_MRI_GARK_ERK45a]
 * - "fast" ERKStep method order: fast_order [default 4]
 *      To put all physics at the slow scale, use "0", otherwise
 *      specify a valid explicit method order.
 * - "slow" MRI temporal adaptivity controller: scontrol [default = 6]
 *      0:  no controller [fixed time steps]
 *      5:  I controller (as part of MRI-HTOL)
 *      6:  I controller (alone)
 *      7:  PI controller (as part of MRI-HTOL)
 *      8:  PI controller (alone)
 *      9:  PID controller (as part of MRI-HTOL)
 *      10: PID controller (alone)
 *      11: ExpGus controller (as part of MRI-HTOL)
 *      12: ExpGus controller (alone)
 *      13: ImpGus controller (as part of MRI-HTOL)
 *      14: ImpGus controller (alone)
 *      15: ImExGus controller (as part of MRI-HTOL)
 *      16: ImExGus controller (alone)
 *      17: H0211 controller (as part of MRI-HTOL)
 *      18: H0211 controller (alone)
 *      19: H0321 controller (as part of MRI-HTOL)
 *      20: H0321 controller (alone)
 *      21: H211 controller (as part of MRI-HTOL)
 *      22: H211 controller (alone)
 *      23: H312 controller (as part of MRI-HTOL)
 *      24: H312 controller (alone)
 * - "fast" ERKStep temporal adaptivity controller: fcontrol [default = 1]
 *   Note that this will only be used for 5 <= scontrol <= 24.
 *      0:   no controller [fixed time steps]
 *      1:   I controller
 *      2:   PI controller
 *      3:   PID controller
 *      4:   ExpGus controller
 *      5:   ImpGus controller
 *      6:   ImExGus controller
 *      7:   H0211 controller
 *      8:   H0321 controller
 *      9:   H211 controller
 *      10:  H312 controller
 * - "fast" ERKStep accumulated error type: faccum [default = 0]
 *   Note that this will only be used for multirate scontrol options
 *     -1:  no accumulation
 *      0:  maximum accumulation
 *      1:  additive accumulation
 *      2:  average accumulation
 * - controller parameters: (k1s, k2s, k3s, k1f, k2f, k3f,
 *                           bias, htol_relch, htol_minfac, htol_maxfac)
 *     slow single-rate controllers: use k1s through k3s, as appropriate.
 *     fast single-rate controllers: use k1f through k3f, as appropriate.
 *     MRIHTol controllers: use htol_relch, htol_minfac, htol_maxfac.
 *     all controllers (fast and slow) use bias.
 *     ** if any one of a relevant set are "-1" then the defaults are used
 *
 * Outputs and solution error values are printed at equal intervals
 * of 0.5 and run statistics are printed at the end.
 * ----------------------------------------------------------------*/

// Header files
#include <arkode/arkode_erkstep.h> // prototypes for ERKStep fcts., consts
#include <arkode/arkode_mristep.h> // prototypes for MRIStep fcts., consts
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <nvector/nvector_serial.h> // serial N_Vector type, fcts., macros
#include <sunadaptcontroller/sunadaptcontroller_imexgus.h>
#include <sunadaptcontroller/sunadaptcontroller_mrihtol.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>
#include <sundials/sundials_core.hpp>
#include <sundials/sundials_logger.h>
#include <sunlinsol/sunlinsol_dense.h> // dense linear solver
#include <sunmatrix/sunmatrix_dense.h> // dense matrix type, fcts., macros
#include "../../utilities/test_utilities.hpp" // common utility functions

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define ESYM "Le"
#define FSYM "Lf"
#else
#define ESYM "e"
#define FSYM "f"
#endif

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)

// Problem options
struct Options
{
  // Problem parameters
  sunrealtype e = SUN_RCONST(0.5);
  sunrealtype G = SUN_RCONST(-100.0);
  sunrealtype w = SUN_RCONST(100.0);

  // Step sizes and tolerances
  int set_h0            = 0;
  sunrealtype hs        = SUN_RCONST(1.0e-2);
  sunrealtype hf        = SUN_RCONST(1.0e-4);
  sunrealtype rtol      = SUN_RCONST(1.0e-4);
  sunrealtype atol      = SUN_RCONST(1.0e-11);
  sunrealtype fast_rtol = SUN_RCONST(1.0e-4);

  // Method selection
  std::string mri_method = "ARKODE_MRI_GARK_ERK45a";
  int fast_order         = 4;
  int scontrol           = 6;
  int fcontrol           = 1;
  int faccum             = 0;
  int slow_pq            = 0;
  int fast_pq            = 0;

  // controller parameters
  sunrealtype k1s         = SUN_RCONST(-1.0);
  sunrealtype k2s         = SUN_RCONST(-1.0);
  sunrealtype k3s         = SUN_RCONST(-1.0);
  sunrealtype k1f         = SUN_RCONST(-1.0);
  sunrealtype k2f         = SUN_RCONST(-1.0);
  sunrealtype k3f         = SUN_RCONST(-1.0);
  sunrealtype bias        = SUN_RCONST(-1.0);
  sunrealtype htol_relch  = SUN_RCONST(-1.0);
  sunrealtype htol_minfac = SUN_RCONST(-1.0);
  sunrealtype htol_maxfac = SUN_RCONST(-1.0);
  sunrealtype slow_safety = SUN_RCONST(-1.0);
};

// User-supplied functions called by the solver
static int fse(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fsi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int ff(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fn(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Js(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Jsi(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Jn(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Utility functions
static void InputHelp();
static int ReadInputs(std::vector<std::string>& args, Options& opts,
                      SUNContext ctx);
static void PrintSlowAdaptivity(Options opts);
static void PrintFastAdaptivity(Options opts);
static sunrealtype r(sunrealtype t, Options* opts);
static sunrealtype s(sunrealtype t, Options* opts);
static sunrealtype rdot(sunrealtype t, Options* opts);
static sunrealtype sdot(sunrealtype t, Options* opts);
static sunrealtype utrue(sunrealtype t, Options* opts);
static sunrealtype vtrue(sunrealtype t, Options* opts);
static int Ytrue(sunrealtype t, N_Vector y, Options* opts);

// Main Program
int main(int argc, char* argv[])
{
  // SUNDIALS context objects
  sundials::Context sunctx; // main solver
  sundials::Context refctx; // reference solver

  // Read input options
  Options opts;
  std::vector<std::string> args(argv + 1, argv + argc);
  int flag = ReadInputs(args, opts, sunctx);
  if (check_flag(flag, "ReadInputs")) return 1;

  // General problem parameters
  sunrealtype T0   = SUN_RCONST(0.0); // initial time
  sunrealtype Tf   = SUN_RCONST(5.0); // final time
  sunindextype NEQ = 2;               // number of dependent vars.
  int Nt           = 20;              // number of output times

  // Initial problem output
  //    While traversing these, set various function pointers, table constants, and method orders.
  ARKRhsFn f_f, f_se, f_si;
  ARKLsJacFn J_s;
  int retval;
  sunbooleantype slowimplicit, slowimex;
  slowimplicit = slowimex = SUNFALSE;
  f_si                    = NULL;
  J_s                     = NULL;
  f_f                     = (opts.fast_order == 0) ? f0 : ff;
  f_se                    = (opts.fast_order == 0) ? fn : fs;
  if ((opts.mri_method == "ARKODE_MRI_GARK_IRK21a") ||
      (opts.mri_method == "ARKODE_MRI_GARK_ESDIRK34a") ||
      (opts.mri_method == "ARKODE_MRI_GARK_ESDIRK46a"))
  {
    slowimplicit = SUNTRUE;
    f_se         = NULL;
    f_si         = (opts.fast_order == 0) ? fn : fs;
    J_s          = (opts.fast_order == 0) ? Jn : Js;
  }
  if ((opts.mri_method == "ARKODE_IMEX_MRI_SR21") ||
      (opts.mri_method == "ARKODE_IMEX_MRI_SR32") ||
      (opts.mri_method == "ARKODE_IMEX_MRI_SR43"))
  {
    slowimex     = SUNTRUE;
    slowimplicit = SUNTRUE;
    f_se         = (opts.fast_order == 0) ? f0 : fse;
    f_si         = (opts.fast_order == 0) ? fn : fsi;
    J_s          = (opts.fast_order == 0) ? Jn : Jsi;
  }
  std::cout << "\nAdaptive multirate nonlinear Kvaerno-Prothero-Robinson test "
               "problem:\n";
  std::cout << "    time domain:  (" << T0 << "," << Tf << "]\n";
  std::cout << "    G = " << opts.G << std::endl;
  std::cout << "    w = " << opts.w << std::endl;
  std::cout << "    e = " << opts.e << std::endl;
  std::cout << "\n  Slow integrator: " << opts.mri_method;
  if (slowimex) { std::cout << " (ImEx)" << std::endl; }
  else if (slowimplicit) { std::cout << " (implicit)" << std::endl; }
  else { std::cout << " (explicit)" << std::endl; }
  PrintSlowAdaptivity(opts);
  if (opts.fast_order == 0) { std::cout << "\n  Fast integrator disabled"; }
  else { std::cout << "\n  Fast order " << opts.fast_order << std::endl; }
  PrintFastAdaptivity(opts);

  // If SUNLogger is enabled, manually disable it for the reference solver
  SUNLogger logger = NULL;
  retval           = SUNLogger_Create(SUN_COMM_NULL, 0, &logger);
  retval           = SUNContext_SetLogger(refctx, logger);
  retval           = SUNLogger_SetErrorFilename(logger, "/dev/null");
  retval           = SUNLogger_SetWarningFilename(logger, "/dev/null");
  retval           = SUNLogger_SetInfoFilename(logger, "/dev/null");
  retval           = SUNLogger_SetDebugFilename(logger, "/dev/null");

  // Create and initialize serial vectors for the solution and reference
  N_Vector y = N_VNew_Serial(NEQ, sunctx);
  if (check_ptr((void*)y, "N_VNew_Serial")) return 1;
  N_Vector yref = N_VNew_Serial(NEQ, refctx);
  if (check_ptr((void*)yref, "N_VNew_Serial")) return 1;
  sunrealtype* ydata = N_VGetArrayPointer(y);
  if (check_ptr((void*)ydata, "N_VGetArrayPointer")) return 1;
  sunrealtype* yrefdata = N_VGetArrayPointer(yref);
  if (check_ptr((void*)yrefdata, "N_VGetArrayPointer")) return 1;

  // Set initial conditions
  retval = Ytrue(T0, y, &opts);
  if (check_flag(retval, "Ytrue")) return 1;
  N_VScale(ONE, y, yref);

  // Create and configure reference solver object
  void* arkode_ref = ERKStepCreate(fn, T0, yref, refctx);
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
  SUNAdaptController fcontrol = NULL;
  switch (opts.fcontrol)
  {
  case (1):
    fcontrol = SUNAdaptController_I(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_I")) return 1;
    if (opts.k1f > -1)
    {
      retval = SUNAdaptController_SetParams_I(fcontrol, opts.k1f);
      if (check_flag(retval, "SUNAdaptController_SetParams_I")) return 1;
    }
    break;
  case (2):
    fcontrol = SUNAdaptController_PI(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_PI")) return 1;
    if (std::min(opts.k1f, opts.k2f) > -1)
    {
      retval = SUNAdaptController_SetParams_PI(fcontrol, opts.k1f, opts.k2f);
      if (check_flag(retval, "SUNAdaptController_SetParams_PI")) return 1;
    }
    break;
  case (3):
    fcontrol = SUNAdaptController_PID(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_PID")) return 1;
    if (std::min(opts.k1f, std::min(opts.k2f, opts.k3f)) > -1)
    {
      retval = SUNAdaptController_SetParams_PID(fcontrol, opts.k1f, opts.k2f,
                                                opts.k3f);
      if (check_flag(retval, "SUNAdaptController_SetParams_PID")) return 1;
    }
    break;
  case (4):
    fcontrol = SUNAdaptController_ExpGus(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_ExpGus")) return 1;
    if (std::min(opts.k1f, opts.k2f) > -1)
    {
      retval = SUNAdaptController_SetParams_ExpGus(fcontrol, opts.k1f, opts.k2f);
      if (check_flag(retval, "SUNAdaptController_SetParams_ExpGus")) return 1;
    }
    break;
  case (5):
    fcontrol = SUNAdaptController_ImpGus(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_ImpGus")) return 1;
    if (std::min(opts.k1f, opts.k2f) > -1)
    {
      retval = SUNAdaptController_SetParams_ImpGus(fcontrol, opts.k1f, opts.k2f);
      if (check_flag(retval, "SUNAdaptController_SetParams_ImpGus")) return 1;
    }
    break;
  case (6):
    fcontrol = SUNAdaptController_ImExGus(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_ImExGus")) return 1;
    break;
  case (7):
    fcontrol = SUNAdaptController_H0211(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_H0211")) return 1;
    break;
  case (8):
    fcontrol = SUNAdaptController_H0321(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_H0321")) return 1;
    break;
  case (9):
    fcontrol = SUNAdaptController_H211(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_H211")) return 1;
    break;
  case (10):
    fcontrol = SUNAdaptController_H312(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_H312")) return 1;
    break;
  }
  if ((opts.bias > -1) && (opts.fcontrol > 0))
  {
    retval = SUNAdaptController_SetErrorBias(fcontrol, opts.bias);
    if (check_flag(retval, "SUNAdaptController_SetErrorBias")) return 1;
  }

  // Create ERKStep (fast) integrator
  void* inner_arkode_mem = NULL; // ARKode memory structure
  inner_arkode_mem       = ERKStepCreate(f_f, T0, y, sunctx);
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
  MRIStepInnerStepper inner_stepper = NULL; // inner stepper
  retval = ARKodeCreateMRIStepInnerStepper(inner_arkode_mem, &inner_stepper);
  if (check_flag(retval, "ARKodeCreateMRIStepInnerStepper")) return 1;

  // Create slow controller object, and select orders of accuracy as relevant
  SUNAdaptController scontrol     = NULL;
  SUNAdaptController scontrol_H   = NULL;
  SUNAdaptController scontrol_Tol = NULL;
  switch (opts.scontrol)
  {
  case (5):
    scontrol_H = SUNAdaptController_I(sunctx);
    if (check_ptr((void*)scontrol_H, "SUNAdaptController_I (slow H)")) return 1;
    scontrol_Tol = SUNAdaptController_I(sunctx);
    if (check_ptr((void*)scontrol_Tol, "SUNAdaptController_I (slow Tol)"))
      return 1;
    if (opts.k1s > -1)
    {
      retval = SUNAdaptController_SetParams_I(scontrol_H, opts.k1s);
      if (check_flag(retval, "SUNAdaptController_SetParams_I")) return 1;
      retval = SUNAdaptController_SetParams_I(scontrol_Tol, opts.k1s);
      if (check_flag(retval, "SUNAdaptController_SetParams_I")) return 1;
    }
    scontrol = SUNAdaptController_MRIHTol(scontrol_H, scontrol_Tol, sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    if (std::min(opts.htol_relch, std::min(opts.htol_minfac, opts.htol_maxfac)) >
        -1)
    {
      retval = SUNAdaptController_SetParams_MRIHTol(scontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
    }
    break;
  case (6):
    scontrol = SUNAdaptController_I(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptControllerI (slow)")) return 1;
    if (opts.k1s > -1)
    {
      retval = SUNAdaptController_SetParams_I(scontrol, opts.k1s);
      if (check_flag(retval, "SUNAdaptController_SetParams_I")) return 1;
    }
    break;
  case (7):
    scontrol_H = SUNAdaptController_PI(sunctx);
    if (check_ptr((void*)scontrol_H, "SUNAdaptController_PI (slow H)"))
      return 1;
    scontrol_Tol = SUNAdaptController_PI(sunctx);
    if (check_ptr((void*)scontrol_Tol, "SUNAdaptController_PI (slow Tol)"))
      return 1;
    if (std::min(opts.k1s, opts.k2s) > -1)
    {
      retval = SUNAdaptController_SetParams_PI(scontrol_H, opts.k1s, opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_PI")) return 1;
      retval = SUNAdaptController_SetParams_PI(scontrol_Tol, opts.k1s, opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_PI")) return 1;
    }
    scontrol = SUNAdaptController_MRIHTol(scontrol_H, scontrol_Tol, sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    if (std::min(opts.htol_relch, std::min(opts.htol_minfac, opts.htol_maxfac)) >
        -1)
    {
      retval = SUNAdaptController_SetParams_MRIHTol(scontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
    }
    break;
  case (8):
    scontrol = SUNAdaptController_PI(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_PI (slow)")) return 1;
    if (std::min(opts.k1s, opts.k2s) > -1)
    {
      retval = SUNAdaptController_SetParams_PI(scontrol, opts.k1s, opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_PI")) return 1;
    }
    break;
  case (9):
    scontrol_H = SUNAdaptController_PID(sunctx);
    if (check_ptr((void*)scontrol_H, "SUNAdaptController_PID (slow H)"))
      return 1;
    scontrol_Tol = SUNAdaptController_PID(sunctx);
    if (check_ptr((void*)scontrol_Tol, "SUNAdaptController_PID (slow Tol)"))
      return 1;
    if (std::min(opts.k1s, std::min(opts.k2s, opts.k3s)) > -1)
    {
      retval = SUNAdaptController_SetParams_PID(scontrol_H, opts.k1s, opts.k2s,
                                                opts.k3s);
      if (check_flag(retval, "SUNAdaptController_SetParams_PID")) return 1;
      retval = SUNAdaptController_SetParams_PID(scontrol_Tol, opts.k1s,
                                                opts.k2s, opts.k3s);
      if (check_flag(retval, "SUNAdaptController_SetParams_PID")) return 1;
    }
    scontrol = SUNAdaptController_MRIHTol(scontrol_H, scontrol_Tol, sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    if (std::min(opts.htol_relch, std::min(opts.htol_minfac, opts.htol_maxfac)) >
        -1)
    {
      retval = SUNAdaptController_SetParams_MRIHTol(scontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
    }
    break;
  case (10):
    scontrol = SUNAdaptController_PID(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_PID (slow)")) return 1;
    if (std::min(opts.k1s, std::min(opts.k2s, opts.k3s)) > -1)
    {
      retval = SUNAdaptController_SetParams_PID(scontrol, opts.k1s, opts.k2s,
                                                opts.k3s);
      if (check_flag(retval, "SUNAdaptController_SetParams_PID")) return 1;
    }
    break;
  case (11):
    scontrol_H = SUNAdaptController_ExpGus(sunctx);
    if (check_ptr((void*)scontrol_H, "SUNAdaptController_ExpGus (slow H)"))
      return 1;
    scontrol_Tol = SUNAdaptController_ExpGus(sunctx);
    if (check_ptr((void*)scontrol_Tol, "SUNAdaptController_ExpGus (slow Tol)"))
      return 1;
    if (std::min(opts.k1s, opts.k2s) > -1)
    {
      retval = SUNAdaptController_SetParams_ExpGus(scontrol_H, opts.k1s,
                                                   opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_ExpGus")) return 1;
      retval = SUNAdaptController_SetParams_ExpGus(scontrol_Tol, opts.k1s,
                                                   opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_ExpGus")) return 1;
    }
    scontrol = SUNAdaptController_MRIHTol(scontrol_H, scontrol_Tol, sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    if (std::min(opts.htol_relch, std::min(opts.htol_minfac, opts.htol_maxfac)) >
        -1)
    {
      retval = SUNAdaptController_SetParams_MRIHTol(scontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
    }
    break;
  case (12):
    scontrol = SUNAdaptController_ExpGus(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_ExpGus (slow)"))
      return 1;
    if (std::min(opts.k1s, opts.k2s) > -1)
    {
      retval = SUNAdaptController_SetParams_ExpGus(scontrol, opts.k1s, opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_ExpGus")) return 1;
    }
    break;
  case (13):
    scontrol_H = SUNAdaptController_ImpGus(sunctx);
    if (check_ptr((void*)scontrol_H, "SUNAdaptController_ImpGus (slow H)"))
      return 1;
    scontrol_Tol = SUNAdaptController_ImpGus(sunctx);
    if (check_ptr((void*)scontrol_Tol, "SUNAdaptController_ImpGus (slow Tol)"))
      return 1;
    if (std::min(opts.k1s, opts.k2s) > -1)
    {
      retval = SUNAdaptController_SetParams_ImpGus(scontrol_H, opts.k1s,
                                                   opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_ImpGus")) return 1;
      retval = SUNAdaptController_SetParams_ImpGus(scontrol_Tol, opts.k1s,
                                                   opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_ImpGus")) return 1;
    }
    scontrol = SUNAdaptController_MRIHTol(scontrol_H, scontrol_Tol, sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    if (std::min(opts.htol_relch, std::min(opts.htol_minfac, opts.htol_maxfac)) >
        -1)
    {
      retval = SUNAdaptController_SetParams_MRIHTol(scontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
    }
    break;
  case (14):
    scontrol = SUNAdaptController_ImpGus(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_ImpGus (slow)"))
      return 1;
    if (std::min(opts.k1s, opts.k2s) > -1)
    {
      retval = SUNAdaptController_SetParams_ImpGus(scontrol, opts.k1s, opts.k2s);
      if (check_flag(retval, "SUNAdaptController_SetParams_ImpGus")) return 1;
    }
    break;
  case (15):
    scontrol_H = SUNAdaptController_ImExGus(sunctx);
    if (check_ptr((void*)scontrol_H, "SUNAdaptController_ImExGus (slow H)"))
      return 1;
    scontrol_Tol = SUNAdaptController_ImExGus(sunctx);
    if (check_ptr((void*)scontrol_Tol, "SUNAdaptController_ImExGus (slow Tol)"))
      return 1;
    scontrol = SUNAdaptController_MRIHTol(scontrol_H, scontrol_Tol, sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    if (std::min(opts.htol_relch, std::min(opts.htol_minfac, opts.htol_maxfac)) >
        -1)
    {
      retval = SUNAdaptController_SetParams_MRIHTol(scontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
    }
    break;
  case (16):
    scontrol = SUNAdaptController_ImExGus(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_ImExGus (slow)"))
      return 1;
    break;
  case (17):
    scontrol_H = SUNAdaptController_H0211(sunctx);
    if (check_ptr((void*)scontrol_H, "SUNAdaptController_H0211 (slow H)"))
      return 1;
    scontrol_Tol = SUNAdaptController_H0211(sunctx);
    if (check_ptr((void*)scontrol_Tol, "SUNAdaptController_H0211 (slow Tol)"))
      return 1;
    scontrol = SUNAdaptController_MRIHTol(scontrol_H, scontrol_Tol, sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    if (std::min(opts.htol_relch, std::min(opts.htol_minfac, opts.htol_maxfac)) >
        -1)
    {
      retval = SUNAdaptController_SetParams_MRIHTol(scontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
    }
    break;
  case (18):
    scontrol = SUNAdaptController_H0211(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_H0211 (slow)")) return 1;
    break;
  case (19):
    scontrol_H = SUNAdaptController_H0321(sunctx);
    if (check_ptr((void*)scontrol_H, "SUNAdaptController_H0321 (slow H)"))
      return 1;
    scontrol_Tol = SUNAdaptController_H0321(sunctx);
    if (check_ptr((void*)scontrol_Tol, "SUNAdaptController_H0321 (slow Tol)"))
      return 1;
    scontrol = SUNAdaptController_MRIHTol(scontrol_H, scontrol_Tol, sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    if (std::min(opts.htol_relch, std::min(opts.htol_minfac, opts.htol_maxfac)) >
        -1)
    {
      retval = SUNAdaptController_SetParams_MRIHTol(scontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
    }
    break;
  case (20):
    scontrol = SUNAdaptController_H0321(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_H0321 (slow)")) return 1;
    break;
  case (21):
    scontrol_H = SUNAdaptController_H211(sunctx);
    if (check_ptr((void*)scontrol_H, "SUNAdaptController_H211 (slow H)"))
      return 1;
    scontrol_Tol = SUNAdaptController_H211(sunctx);
    if (check_ptr((void*)scontrol_Tol, "SUNAdaptController_H211 (slow Tol)"))
      return 1;
    scontrol = SUNAdaptController_MRIHTol(scontrol_H, scontrol_Tol, sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    if (std::min(opts.htol_relch, std::min(opts.htol_minfac, opts.htol_maxfac)) >
        -1)
    {
      retval = SUNAdaptController_SetParams_MRIHTol(scontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
    }
    break;
  case (22):
    scontrol = SUNAdaptController_H211(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_H211 (slow)")) return 1;
    break;
  case (23):
    scontrol_H = SUNAdaptController_H312(sunctx);
    if (check_ptr((void*)scontrol_H, "SUNAdaptController_H312 (slow H)"))
      return 1;
    scontrol_Tol = SUNAdaptController_H312(sunctx);
    if (check_ptr((void*)scontrol_Tol, "SUNAdaptController_H312 (slow Tol)"))
      return 1;
    scontrol = SUNAdaptController_MRIHTol(scontrol_H, scontrol_Tol, sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    if (std::min(opts.htol_relch, std::min(opts.htol_minfac, opts.htol_maxfac)) >
        -1)
    {
      retval = SUNAdaptController_SetParams_MRIHTol(scontrol, opts.htol_relch,
                                                    opts.htol_minfac,
                                                    opts.htol_maxfac);
      if (check_flag(retval, "SUNAdaptController_SetParams_MRIHTol")) return 1;
    }
    break;
  case (24):
    scontrol = SUNAdaptController_H312(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_H312 (slow)")) return 1;
    break;
  }
  if ((opts.bias > -1) && (opts.scontrol > 0))
  {
    retval = SUNAdaptController_SetErrorBias(scontrol, opts.bias);
    if (check_flag(retval, "SUNAdaptController_SetErrorBias")) return 1;
  }

  // Create MRI (slow) integrator
  void* arkode_mem = NULL; // ARKode memory structure
  arkode_mem       = MRIStepCreate(f_se, f_si, T0, y, inner_stepper, sunctx);
  if (check_ptr((void*)arkode_mem, "MRIStepCreate")) return 1;
  MRIStepCoupling C = MRIStepCoupling_LoadTableByName((opts.mri_method).c_str());
  if (check_ptr((void*)C, "MRIStepCoupling_LoadTableByName")) return 1;
  retval = MRIStepSetCoupling(arkode_mem, C);
  if (check_flag(retval, "MRIStepSetCoupling")) return 1;
  SUNMatrix As        = NULL; // matrix for slow solver
  SUNLinearSolver LSs = NULL; // slow linear solver object
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
    if (opts.slow_safety > -1)
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
  sunrealtype t     = T0;
  sunrealtype t2    = T0;
  sunrealtype dTout = (Tf - T0) / Nt;
  sunrealtype tout  = T0 + dTout;
  sunrealtype u, v, uerr, verr, uerrtot, verrtot, errtot, accuracy;
  uerr = verr = uerrtot = verrtot = errtot = accuracy = ZERO;
  printf("        t           u           v       uerr      verr\n");
  printf("   ------------------------------------------------------\n");
  printf("  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM "  %.2" ESYM "  %.2" ESYM
         "\n",
         t, ydata[0], ydata[1], uerr, verr);
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
    uerr = std::abs(yrefdata[0] - u);
    verr = std::abs(yrefdata[1] - v);
    uerrtot += uerr * uerr;
    verrtot += verr * verr;
    errtot += uerr * uerr + verr * verr;
    accuracy = std::max(accuracy,
                        uerr / std::abs(opts.atol + opts.rtol * yrefdata[0]));
    accuracy = std::max(accuracy,
                        verr / std::abs(opts.atol + opts.rtol * yrefdata[1]));

    // Periodically output current results to screen
    if (t >= tout)
    {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
      printf("  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM "  %.2" ESYM
             "  %.2" ESYM "\n",
             t, u, v, uerr, verr);
    }
  }
  printf("   ------------------------------------------------------\n");

  //
  // Finalize
  //

  // Get some slow integrator statistics
  long int nsts, natts, netfs, nfse, nfsi;
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
  errtot  = std::sqrt(errtot / SUN_RCONST(2.0) / (sunrealtype)nsts);
  std::cout << "\nFinal Solver Statistics:\n";
  std::cout << "   Slow steps = " << nsts << "  (attempts = " << natts
            << ",  fails = " << netfs << ")\n";
  std::cout << "   Fast steps = " << nstf << "  (attempts = " << nattf
            << ",  fails = " << netff << ")\n";
  std::cout << "   u error = " << uerrtot << ", v error = " << verrtot
            << ", total error = " << errtot << std::endl;
  std::cout << "   Relative accuracy = " << accuracy << std::endl;
  std::cout << "   Total RHS evals:  Fse = " << nfse << ", Fsi = " << nfsi
            << ", Ff = " << nff << std::endl;

  // Get/print slow integrator decoupled implicit solver statistics
  if (slowimplicit)
  {
    long int nnis, nncs, njes;
    retval = ARKodeGetNonlinSolvStats(arkode_mem, &nnis, &nncs);
    check_flag(retval, "ARKodeGetNonlinSolvStats");
    retval = ARKodeGetNumJacEvals(arkode_mem, &njes);
    check_flag(retval, "ARKodeGetNumJacEvals");
    std::cout << "   Slow Newton iters = " << nnis << std::endl;
    std::cout << "   Slow Newton conv fails = " << nncs << std::endl;
    std::cout << "   Slow Jacobian evals = " << njes << std::endl;
  }

  // Clean up and return
  N_VDestroy(y);
  N_VDestroy(yref);
  MRIStepCoupling_Free(C);
  if (As) { SUNMatDestroy(As); }
  if (LSs) { SUNLinSolFree(LSs); }
  if (scontrol) { SUNAdaptController_Destroy(scontrol); }
  if (scontrol_H) { SUNAdaptController_Destroy(scontrol_H); }
  if (scontrol_Tol) { SUNAdaptController_Destroy(scontrol_Tol); }
  if (fcontrol) { SUNAdaptController_Destroy(fcontrol); }
  ARKodeFree(&inner_arkode_mem);            // Free fast integrator memory
  MRIStepInnerStepper_Free(&inner_stepper); // Free inner stepper structure
  ARKodeFree(&arkode_mem);                  // Free slow integrator memory
  ARKodeFree(&arkode_ref);                  // Free reference solver memory
  SUNLogger_Destroy(&logger);               // Free logger

  return 0;
}

// ------------------------------
// Functions called by the solver
// -----------------------------

// ff routine to compute the fast portion of the ODE RHS.
static int ff(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Options* opts       = (Options*)user_data;
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  sunrealtype* dydata = N_VGetArrayPointer(ydot);
  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];
  sunrealtype tmp1, tmp2;

  // fill in the RHS function:
  //   [0  0]*[(-2+u^2-r(t))/(2*u)] + [     0      ]
  //   [e -1] [(-2+v^2-s(t))/(2*v)]   [sdot(t)/(2v)]
  tmp1      = (-TWO + u * u - r(t, opts)) / (TWO * u);
  tmp2      = (-TWO + v * v - s(t, opts)) / (TWO * v);
  dydata[0] = ZERO;
  dydata[1] = opts->e * tmp1 - tmp2 + sdot(t, opts) / (TWO * v);

  // Return with success
  return 0;
}

// fs routine to compute the slow portion of the ODE RHS.
static int fs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Options* opts       = (Options*)user_data;
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  sunrealtype* dydata = N_VGetArrayPointer(ydot);
  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];
  sunrealtype tmp1, tmp2;

  // fill in the RHS function:
  //   [G  e]*[(-2+u^2-r(t))/(2*u)] + [rdot(t)/(2u)]
  //   [0  0] [(-2+v^2-s(t))/(2*v)]   [      0     ]
  tmp1      = (-TWO + u * u - r(t, opts)) / (TWO * u);
  tmp2      = (-TWO + v * v - s(t, opts)) / (TWO * v);
  dydata[0] = opts->G * tmp1 + opts->e * tmp2 + rdot(t, opts) / (TWO * u);
  dydata[1] = ZERO;

  // Return with success
  return 0;
}

// fse routine to compute the slow portion of the ODE RHS.
static int fse(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Options* opts       = (Options*)user_data;
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  sunrealtype* dydata = N_VGetArrayPointer(ydot);
  const sunrealtype u = ydata[0];

  // fill in the slow explicit RHS function:
  //   [rdot(t)/(2u)]
  //   [      0     ]
  dydata[0] = rdot(t, opts) / (TWO * u);
  dydata[1] = ZERO;

  // Return with success
  return 0;
}

// fsi routine to compute the slow portion of the ODE RHS.(currently same as fse)
static int fsi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Options* opts       = (Options*)user_data;
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  sunrealtype* dydata = N_VGetArrayPointer(ydot);
  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];
  sunrealtype tmp1, tmp2;

  // fill in the slow implicit RHS function:
  //   [G  e]*[(-2+u^2-r(t))/(2*u)]
  //   [0  0] [(-2+v^2-s(t))/(2*v)]
  tmp1      = (-TWO + u * u - r(t, opts)) / (TWO * u);
  tmp2      = (-TWO + v * v - s(t, opts)) / (TWO * v);
  dydata[0] = opts->G * tmp1 + opts->e * tmp2;
  dydata[1] = ZERO;

  // Return with success
  return 0;
}

static int fn(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Options* opts       = (Options*)user_data;
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  sunrealtype* dydata = N_VGetArrayPointer(ydot);
  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];
  sunrealtype tmp1, tmp2;

  // fill in the RHS function:
  //   [G  e]*[(-2+u^2-r(t))/(2*u)] + [rdot(t)/(2u)]
  //   [e -1] [(-2+v^2-s(t))/(2*v)]   [sdot(t)/(2v)]
  tmp1      = (-TWO + u * u - r(t, opts)) / (TWO * u);
  tmp2      = (-TWO + v * v - s(t, opts)) / (TWO * v);
  dydata[0] = opts->G * tmp1 + opts->e * tmp2 + rdot(t, opts) / (TWO * u);
  dydata[1] = opts->e * tmp1 - tmp2 + sdot(t, opts) / (TWO * v);

  // Return with success
  return 0;
}

static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VConst(ZERO, ydot);
  return (0);
}

static int Js(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  Options* opts       = (Options*)user_data;
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];
  sunrealtype t11, t22;

  // fill in the Jacobian:
  //   [G  e]*[1-(u^2-r(t)-2)/(2*u^2),  0] + [-r'(t)/(2*u^2),  0]
  //   [0  0] [0,  1-(v^2-s(t)-2)/(2*v^2)]   [0,               0]
  t11                   = ONE - (u * u - r(t, opts) - TWO) / (TWO * u * u);
  t22                   = ONE - (v * v - s(t, opts) - TWO) / (TWO * v * v);
  SM_ELEMENT_D(J, 0, 0) = opts->G * t11 - rdot(t, opts) / (TWO * u * u);
  SM_ELEMENT_D(J, 0, 1) = opts->e * t22;
  SM_ELEMENT_D(J, 1, 0) = ZERO;
  SM_ELEMENT_D(J, 1, 1) = ZERO;

  // Return with success
  return 0;
}

static int Jsi(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  Options* opts       = (Options*)user_data;
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];
  sunrealtype t11, t22;

  // fill in the Jacobian:
  //   [G  e]*[1-(u^2-r(t)-2)/(2*u^2),  0]
  //   [0  0] [0,  1-(v^2-s(t)-2)/(2*v^2)]
  t11                   = ONE - (u * u - r(t, opts) - TWO) / (TWO * u * u);
  t22                   = ONE - (v * v - s(t, opts) - TWO) / (TWO * v * v);
  SM_ELEMENT_D(J, 0, 0) = opts->G * t11;
  SM_ELEMENT_D(J, 0, 1) = opts->e * t22;
  SM_ELEMENT_D(J, 1, 0) = ZERO;
  SM_ELEMENT_D(J, 1, 1) = ZERO;

  // Return with success
  return 0;
}

static int Jn(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  Options* opts       = (Options*)user_data;
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  const sunrealtype u = ydata[0];
  const sunrealtype v = ydata[1];
  sunrealtype t11, t22;

  // fill in the Jacobian:
  //   [G  e]*[1-(u^2-r(t)-2)/(2*u^2),  0] + [-r'(t)/(2*u^2),  0]
  //   [e -1] [0,  1-(v^2-s(t)-2)/(2*v^2)]   [0,  -s'(t)/(2*v^2)]
  t11                   = ONE - (u * u - r(t, opts) - TWO) / (TWO * u * u);
  t22                   = ONE - (v * v - s(t, opts) - TWO) / (TWO * v * v);
  SM_ELEMENT_D(J, 0, 0) = opts->G * t11 - rdot(t, opts) / (TWO * u * u);
  SM_ELEMENT_D(J, 0, 1) = opts->e * t22;
  SM_ELEMENT_D(J, 1, 0) = opts->e * t11;
  SM_ELEMENT_D(J, 1, 1) = -t22 - sdot(t, opts) / (TWO * v * v);

  // Return with success
  return 0;
}

// ------------------------------
// Private helper functions
// -----------------------------

static sunrealtype r(sunrealtype t, Options* opts) { return (cos(t)); }

static sunrealtype s(sunrealtype t, Options* opts)
{
  return (cos(opts->w * t * (ONE + exp(-(t - TWO) * (t - TWO)))));
}

static sunrealtype rdot(sunrealtype t, Options* opts) { return (-sin(t)); }

static sunrealtype sdot(sunrealtype t, Options* opts)
{
  const sunrealtype tTwo  = t - TWO;
  const sunrealtype eterm = exp(-tTwo * tTwo);
  return (-sin(opts->w * t * (ONE + eterm)) * opts->w *
          (ONE + eterm * (ONE - TWO * t * tTwo)));
}

static sunrealtype utrue(sunrealtype t, Options* opts)
{
  return (std::sqrt(TWO + r(t, opts)));
}

static sunrealtype vtrue(sunrealtype t, Options* opts)
{
  return (std::sqrt(TWO + s(t, opts)));
}

static int Ytrue(sunrealtype t, N_Vector y, Options* opts)
{
  sunrealtype* ydata = N_VGetArrayPointer(y);
  ydata[0]           = utrue(t, opts);
  ydata[1]           = vtrue(t, opts);
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
  std::cout << "  --e            : fast/slow coupling strength\n";
  std::cout << "  --G            : stiffness at slow time scale\n";
  std::cout << "  --w            : time-scale separation factor\n";
  std::cout << "  --hs           : slow (fixed/initial) step size\n";
  std::cout << "  --hf           : fast (fixed/initial) step size\n";
  std::cout
    << "  --set_h0       : use hs/hf above to set the initial step size\n";
  std::cout << "  --rtol         : relative solution tolerance\n";
  std::cout << "  --atol         : absolute solution tolerance\n";
  std::cout
    << "  --fast_rtol    : relative solution tolerance for fast method\n";
  std::cout << "  --mri_method   : MRI method name (valid ARKODE_MRITableID)\n";
  std::cout << "  --fast_order   : fast RK method order\n";
  std::cout << "  --scontrol     : slow time step controller, int in [0,16] "
               "(see source)\n";
  std::cout << "  --fcontrol     : fast time step controller, int in [0,6] "
               "(see source)\n";
  std::cout << "  --faccum       : fast error accumulation type {-1,0,1,2}\n";
  std::cout << "  --slow_pq      : use p (0) vs q (1) for slow adaptivity\n";
  std::cout << "  --fast_pq      : use p (0) vs q (1) for fast adaptivity\n";
  std::cout << "  --k1s, --k2s, ..., -k6s : slow controller parameters\n";
  std::cout << "  --k1f, --k2f, -k3f : fast controller parameters\n";
  std::cout << "  --bias : slow and fast controller bias factors\n";
  std::cout << "  --safety : slow time step safety factor\n";
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
  find_arg(args, "--e", opts.e);
  find_arg(args, "--G", opts.G);
  find_arg(args, "--w", opts.w);
  find_arg(args, "--hs", opts.hs);
  find_arg(args, "--hf", opts.hf);
  find_arg(args, "--set_h0", opts.set_h0);
  find_arg(args, "--rtol", opts.rtol);
  find_arg(args, "--atol", opts.atol);
  find_arg(args, "--fast_rtol", opts.fast_rtol);
  find_arg(args, "--mri_method", opts.mri_method);
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
  find_arg(args, "--safety", opts.slow_safety);
  find_arg(args, "--htol_relch", opts.htol_relch);
  find_arg(args, "--htol_minfac", opts.htol_minfac);
  find_arg(args, "--htol_maxfac", opts.htol_maxfac);

  // Check inputs for validity
  //   0 < rtol < 1
  if ((opts.rtol < ZERO) || (opts.rtol > ONE))
  {
    std::cerr << "ERROR: rtol must be in (0,1), (" << opts.rtol << " input)\n";
    return -1;
  }
  //   0 < atol < 1
  if ((opts.atol < ZERO) || (opts.atol > ONE))
  {
    std::cerr << "ERROR: atol must be in (0,1), (" << opts.atol << " input)\n";
    return -1;
  }
  //   0 < fast_rtol < 1
  if ((opts.fast_rtol < ZERO) || (opts.fast_rtol > ONE))
  {
    std::cerr << "ERROR: fast_rtol must be in (0,1), (" << opts.fast_rtol
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
  //   scontrol in [0,24]
  if ((opts.scontrol < 0) || (opts.scontrol > 24))
  {
    std::cerr << "ERROR: scontrol must be in [0,24], (" << opts.scontrol
              << " input)\n";
    return -1;
  }
  //   fcontrol in [0,10]
  if ((opts.fcontrol < 0) || (opts.fcontrol > 10))
  {
    std::cerr << "ERROR: fcontrol must be in [0,10], (" << opts.fcontrol
              << " input)\n";
    return -1;
  }
  //   hs > 0 if scontrol == 0
  if ((opts.hs <= 0) && (opts.scontrol == 0))
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
  //   w >= 1.0
  if (opts.w < ONE)
  {
    std::cerr << "ERROR: w must be >= 1.0, (" << opts.w << " input)\n";
    return -1;
  }

  return 0;
}

static void PrintSlowAdaptivity(Options opts)
{
  switch (opts.scontrol)
  {
  case (0):
    std::cout << "    fixed steps, hs = " << opts.hs << std::endl;
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (5):
    std::cout
      << "    MRI-HTOL controller (using I for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    std::cout << "    fast error accumulation strategy = " << opts.faccum << "\n";
    if (opts.k1s > -1)
    {
      std::cout << "    slow controller parameter: " << opts.k1s << "\n";
    }
    if (std::min(opts.htol_relch, std::min(opts.htol_minfac, opts.htol_maxfac)) >
        -1)
    {
      std::cout << "    HTol controller parameters: " << opts.htol_relch << " "
                << opts.htol_minfac << " " << opts.htol_maxfac << "\n";
    }
    break;
  case (6):
    std::cout << "    Decoupled I controller for slow time scale, based on "
                 "order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    if (opts.k1s > -1)
    {
      std::cout << "    slow controller parameter: " << opts.k1s << "\n";
    }
    break;
  case (7):
    std::cout
      << "    MRI-HTOL controller (using PI for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    std::cout << "    fast error accumulation strategy = " << opts.faccum << "\n";
    if (std::min(opts.k1s, opts.k2s) > -1)
    {
      std::cout << "    slow controller parameters: " << opts.k1s << " "
                << opts.k2s << "\n";
    }
    if (std::min(opts.htol_relch, std::min(opts.htol_minfac, opts.htol_maxfac)) >
        -1)
    {
      std::cout << "    HTol controller parameters: " << opts.htol_relch << " "
                << opts.htol_minfac << " " << opts.htol_maxfac << "\n";
    }
    break;
  case (8):
    std::cout << "    Decoupled PI controller for slow time scale, based on "
                 "order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    if (std::min(opts.k1s, opts.k2s) > -1)
    {
      std::cout << "    slow controller parameters: " << opts.k1s << " "
                << opts.k2s << "\n";
    }
    break;
  case (9):
    std::cout
      << "    MRI-HTOL controller (using PID for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    std::cout << "    fast error accumulation strategy = " << opts.faccum << "\n";
    if (std::min(opts.k1s, std::min(opts.k2s, opts.k3s)) > -1)
    {
      std::cout << "    slow controller parameters: " << opts.k1s << " "
                << opts.k2s << " " << opts.k3s << "\n";
    }
    if (std::min(opts.htol_relch, std::min(opts.htol_minfac, opts.htol_maxfac)) >
        -1)
    {
      std::cout << "    HTol controller parameters: " << opts.htol_relch << " "
                << opts.htol_minfac << " " << opts.htol_maxfac << "\n";
    }
    break;
  case (10):
    std::cout << "    Decoupled PID controller for slow time scale, based on "
                 "order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    if (std::min(opts.k1s, std::min(opts.k2s, opts.k3s)) > -1)
    {
      std::cout << "    slow controller parameters: " << opts.k1s << " "
                << opts.k2s << " " << opts.k3s << "\n";
    }
    break;
  case (11):
    std::cout
      << "    MRI-HTOL controller (using ExpGus for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    std::cout << "    fast error accumulation strategy = " << opts.faccum << "\n";
    if (std::min(opts.k1s, opts.k2s) > -1)
    {
      std::cout << "    slow controller parameters: " << opts.k1s << " "
                << opts.k2s << "\n";
    }
    if (std::min(opts.htol_relch, std::min(opts.htol_minfac, opts.htol_maxfac)) >
        -1)
    {
      std::cout << "    HTol controller parameters: " << opts.htol_relch << " "
                << opts.htol_minfac << " " << opts.htol_maxfac << "\n";
    }
    break;
  case (12):
    std::cout << "    Decoupled ExpGus controller for slow time scale, based "
                 "on order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    if (std::min(opts.k1s, opts.k2s) > -1)
    {
      std::cout << "    slow controller parameters: " << opts.k1s << " "
                << opts.k2s << "\n";
    }
    break;
  case (13):
    std::cout
      << "    MRI-HTOL controller (using ImpGus for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    std::cout << "    fast error accumulation strategy = " << opts.faccum << "\n";
    if (std::min(opts.k1s, opts.k2s) > -1)
    {
      std::cout << "    slow controller parameters: " << opts.k1s << " "
                << opts.k2s << "\n";
    }
    if (std::min(opts.htol_relch, std::min(opts.htol_minfac, opts.htol_maxfac)) >
        -1)
    {
      std::cout << "    HTol controller parameters: " << opts.htol_relch << " "
                << opts.htol_minfac << " " << opts.htol_maxfac << "\n";
    }
    break;
  case (14):
    std::cout << "    Decoupled ImpGus controller for slow time scale, based "
                 "on order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    if (std::min(opts.k1s, opts.k2s) > -1)
    {
      std::cout << "    slow controller parameters: " << opts.k1s << " "
                << opts.k2s << "\n";
    }
    break;
  case (15):
    std::cout
      << "    MRI-HTOL controller (using ImExGus for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    std::cout << "    fast error accumulation strategy = " << opts.faccum << "\n";
    break;
  case (16):
    std::cout << "    Decoupled ImExGus controller for slow time scale, based "
                 "on order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (17):
    std::cout
      << "    MRI-HTOL controller (using H0211 for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    std::cout << "    fast error accumulation strategy = " << opts.faccum << "\n";
    break;
  case (18):
    std::cout << "    Decoupled H0211 controller for slow time scale, based "
                 "on order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (19):
    std::cout
      << "    MRI-HTOL controller (using H0321 for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    std::cout << "    fast error accumulation strategy = " << opts.faccum << "\n";
    break;
  case (20):
    std::cout << "    Decoupled H0321 controller for slow time scale, based "
                 "on order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (21):
    std::cout
      << "    MRI-HTOL controller (using H211 for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    std::cout << "    fast error accumulation strategy = " << opts.faccum << "\n";
    break;
  case (22):
    std::cout << "    Decoupled H211 controller for slow time scale, based "
                 "on order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (23):
    std::cout
      << "    MRI-HTOL controller (using H312 for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    std::cout << "    fast error accumulation strategy = " << opts.faccum << "\n";
    break;
  case (24):
    std::cout << "    Decoupled H312 controller for slow time scale, based "
                 "on order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  }
  if (opts.bias > -1)
  {
    std::cout << "    controller bias factor: " << opts.bias << "\n";
  }
  if (opts.slow_safety > -1)
  {
    std::cout << "    slow step safety factor: " << opts.slow_safety << "\n";
  }
}

static void PrintFastAdaptivity(Options opts)
{
  switch (opts.fcontrol)
  {
  case (0):
    std::cout << "    fixed steps, hf = " << opts.hf << std::endl;
    std::cout << "    fast_rtol = " << opts.fast_rtol
              << ", atol = " << opts.atol << "\n";
    break;
  case (1):
    std::cout << "    I controller for fast time scale, based on order of RK "
              << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    fast_rtol = " << opts.fast_rtol
              << ", atol = " << opts.atol << "\n";
    if (opts.k1f > -1)
    {
      std::cout << "    fast controller parameter: " << opts.k1f << "\n";
    }
    break;
  case (2):
    std::cout << "    PI controller for fast time scale, based on order of RK "
              << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    fast_rtol = " << opts.fast_rtol
              << ", atol = " << opts.atol << "\n";
    if (std::min(opts.k1f, opts.k2f) > -1)
    {
      std::cout << "    fast controller parameters: " << opts.k1f << " "
                << opts.k2f << "\n";
    }
    break;
  case (3):
    std::cout << "    PID controller for fast time scale, based on order of RK "
              << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    fast_rtol = " << opts.fast_rtol
              << ", atol = " << opts.atol << "\n";
    if (std::min(opts.k1f, std::min(opts.k2f, opts.k3f)) > -1)
    {
      std::cout << "    fast controller parameters: " << opts.k1f << " "
                << opts.k2f << " " << opts.k3f << "\n";
    }
    break;
  case (4):
    std::cout
      << "    ExpGus controller for fast time scale, based on order of RK "
      << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    fast_rtol = " << opts.fast_rtol
              << ", atol = " << opts.atol << "\n";
    if (std::min(opts.k1f, opts.k2f) > -1)
    {
      std::cout << "    fast controller parameters: " << opts.k1f << " "
                << opts.k2f << "\n";
    }
    break;
  case (5):
    std::cout
      << "    ImpGus controller for fast time scale, based on order of RK "
      << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    fast_rtol = " << opts.fast_rtol
              << ", atol = " << opts.atol << "\n";
    if (std::min(opts.k1f, opts.k2f) > -1)
    {
      std::cout << "    fast controller parameters: " << opts.k1f << " "
                << opts.k2f << "\n";
    }
    break;
  case (6):
    std::cout
      << "    ImExGus controller for fast time scale, based on order of RK "
      << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    fast_rtol = " << opts.fast_rtol
              << ", atol = " << opts.atol << "\n";
    break;
  case (7):
    std::cout
      << "    H0211 controller for fast time scale, based on order of RK "
      << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    fast_rtol = " << opts.fast_rtol
              << ", atol = " << opts.atol << "\n";
    break;
  case (8):
    std::cout
      << "    H0321 controller for fast time scale, based on order of RK "
      << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    fast_rtol = " << opts.fast_rtol
              << ", atol = " << opts.atol << "\n";
    break;
  case (9):
    std::cout
      << "    H211 controller for fast time scale, based on order of RK "
      << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    fast_rtol = " << opts.fast_rtol
              << ", atol = " << opts.atol << "\n";
    break;
  case (10):
    std::cout
      << "    H312 controller for fast time scale, based on order of RK "
      << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    fast_rtol = " << opts.fast_rtol
              << ", atol = " << opts.atol << "\n";
    break;
  }
}

//---- end of file ----//
