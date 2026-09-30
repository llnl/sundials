/* ------------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 * ------------------------------------------------------------------
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
 * This example solves the Lotka--Volterra ODE with four parameters,
 *
 *     u = [dx/dt] = [ p_0*x - p_1*x*y  ]
 *         [dy/dt]   [ -p_2*y + p_3*x*y ],
 *
 * using MRIStep. The first equation is the slow partition and the second is
 * the fast partition. The initial condition is u(t_0) = 1 and the parameters
 * are p = [1.5, 1.0, 3.0, 1.0]. After the forward solve, discrete adjoint
 * sensitivity analysis computes the gradient of
 *
 *     g(u(t_f), p) = ||1 - u(t_f, p)||^2 / 2
 *
 * with respect to the initial condition and the parameters.
 *
 * ./ark_lotka_volterra_mri_ASA options:
 * --tf <real>         final simulation time
 * --hs <real>         slow timestep size
 * --hf <real>         fast timestep size
 * --check-freq <int>  how often to checkpoint (in steps)
 * --dont-keep         don't keep checkpoints after loading
 * --help              print these options
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <arkode/arkode.h>
#include <arkode/arkode_erkstep.h>
#include <arkode/arkode_mristep.h>
#include <nvector/nvector_manyvector.h>
#include <nvector/nvector_serial.h>
#include <sunadjointcheckpointscheme/sunadjointcheckpointscheme_fixed.h>
#include <sundials/sundials_adjointstepper.h>
#include <sundials/sundials_core.h>
#include <sunmemory/sunmemory_system.h>

typedef struct
{
  sunrealtype tf;
  sunrealtype hs;
  sunrealtype hf;
  int check_freq;
  sunbooleantype keep_checks;
} ProgramArgs;

static sunrealtype params[4] = {SUN_RCONST(1.5), SUN_RCONST(1.0),
                                SUN_RCONST(3.0), SUN_RCONST(1.0)};

static void parse_args(int argc, char* argv[], ProgramArgs* args);
static void print_help(int argc, char* argv[], int exit_code);
static int check_retval(void* retval_ptr, const char* funcname, int opt);
static int slow_rhs(sunrealtype t, N_Vector u, N_Vector udot, void* user_data);
static int fast_rhs(sunrealtype t, N_Vector u, N_Vector udot, void* user_data);
static int slow_adj_rhs(sunrealtype t, N_Vector y, N_Vector sens,
                        N_Vector sens_dot, void* user_data);
static int fast_adj_rhs(sunrealtype t, N_Vector y, N_Vector sens,
                        N_Vector sens_dot, void* user_data);
static void dgdu(N_Vector u, N_Vector dg);
static void dgdp(N_Vector dg);

int main(int argc, char* argv[])
{
  int retval        = 0;
  SUNContext sunctx = NULL;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  ProgramArgs args;
  args.tf          = SUN_RCONST(10.0);
  args.hs          = SUN_RCONST(1.0) / SUN_RCONST(128.0);
  args.hf          = SUN_RCONST(1.0) / SUN_RCONST(3840.0);
  args.check_freq  = 2;
  args.keep_checks = SUNTRUE;
  parse_args(argc, argv, &args);

  const sunrealtype t0 = SUN_RCONST(0.0);
  const sunrealtype tf = args.tf;
  const int nsteps     = (int)ceil((tf - t0) / args.hs);
  const int nfast      = (int)ceil(args.hs / args.hf);

  N_Vector u = N_VNew_Serial(2, sunctx);
  if (check_retval(u, "N_VNew_Serial", 0)) { return 1; }
  N_VConst(SUN_RCONST(1.0), u);

  /* Create the ERK fast integrator and expose it as a SUNStepper. */
  void* inner_mem = ERKStepCreate(fast_rhs, t0, u, sunctx);
  if (check_retval(inner_mem, "ERKStepCreate", 0)) { return 1; }
  retval = ARKodeSetUserData(inner_mem, params);
  if (check_retval(&retval, "ARKodeSetUserData", 1)) { return 1; }
  retval = ARKodeSetOrder(inner_mem, 4);
  if (check_retval(&retval, "ARKodeSetOrder", 1)) { return 1; }
  retval = ARKodeSetFixedStep(inner_mem, args.hf);
  if (check_retval(&retval, "ARKodeSetFixedStep", 1)) { return 1; }
  retval = ARKodeSetMaxNumSteps(inner_mem, nfast + 1);
  if (check_retval(&retval, "ARKodeSetMaxNumSteps", 1)) { return 1; }

  SUNStepper inner_stepper = NULL;
  retval                   = ARKodeCreateSUNStepper(inner_mem, &inner_stepper);
  if (check_retval(&retval, "ARKodeCreateSUNStepper", 1)) { return 1; }

  /* Create the outer MRIStep integrator with a third-order MRI-GARK method. */
  void* arkode_mem = MRIStepCreate(slow_rhs, NULL, t0, u, inner_stepper, sunctx);
  if (check_retval(arkode_mem, "MRIStepCreate", 0)) { return 1; }
  retval = ARKodeSetUserData(arkode_mem, params);
  if (check_retval(&retval, "ARKodeSetUserData", 1)) { return 1; }
  retval = ARKodeSetFixedStep(arkode_mem, args.hs);
  if (check_retval(&retval, "ARKodeSetFixedStep", 1)) { return 1; }
  retval = ARKodeSetMaxNumSteps(arkode_mem, nsteps + 1);
  if (check_retval(&retval, "ARKodeSetMaxNumSteps", 1)) { return 1; }
  retval = ARKodeSetStopTime(arkode_mem, tf);
  if (check_retval(&retval, "ARKodeSetStopTime", 1)) { return 1; }

  MRIStepCoupling coupling = MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_ERK33a);
  if (check_retval(coupling, "MRIStepCoupling_LoadTable", 0)) { return 1; }
  retval = MRIStepSetCoupling(arkode_mem, coupling);
  MRIStepCoupling_Free(coupling);
  if (check_retval(&retval, "MRIStepSetCoupling", 1)) { return 1; }

  SUNMemoryHelper mem_helper = SUNMemoryHelper_Sys(sunctx);
  if (check_retval(mem_helper, "SUNMemoryHelper_Sys", 0)) { return 1; }

  SUNAdjointCheckpointScheme outer_scheme = NULL;
  retval = SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM,
                                                   mem_helper, args.check_freq,
                                                   3 * nsteps, args.keep_checks,
                                                   sunctx, &outer_scheme);
  if (check_retval(&retval, "SUNAdjointCheckpointScheme_Create_Fixed", 1))
  {
    return 1;
  }
  retval = ARKodeSetAdjointCheckpointScheme(arkode_mem, outer_scheme);
  if (check_retval(&retval, "ARKodeSetAdjointCheckpointScheme", 1))
  {
    return 1;
  }

  printf("Initial condition:\n");
  N_VPrint(u);

  sunrealtype tret = t0;
  retval           = ARKodeEvolve(arkode_mem, tf, u, &tret, ARK_NORMAL);
  if (check_retval(&retval, "ARKodeEvolve", 1)) { return 1; }

  printf("Forward Solution:\n");
  N_VPrint(u);
  printf("MRIStep Stats for Forward Solution:\n");
  retval = ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_retval(&retval, "ARKodePrintAllStats", 1)) { return 1; }
  printf("\n");

  /* Construct the public terminal adjoint state [lambda, mu]. */
  N_Vector lambda = N_VClone(u);
  N_Vector mu     = N_VNew_Serial(4, sunctx);
  if (check_retval(lambda, "N_VClone", 0) || check_retval(mu, "N_VNew_Serial", 0))
  {
    return 1;
  }
  dgdu(u, lambda);
  dgdp(mu);
  N_Vector sens[2] = {lambda, mu};
  N_Vector sf      = N_VNew_ManyVector(2, sens, sunctx);
  if (check_retval(sf, "N_VNew_ManyVector", 0)) { return 1; }

  printf("Adjoint terminal condition:\n");
  N_VPrint(sf);

  /* Augment the fast adjoint problem with the MRI forcing moments. */
  MRIStepInnerAdjointProblem inner_problem = NULL;
  retval = MRIStepInnerAdjointProblem_Create(arkode_mem, fast_adj_rhs, sf,
                                             params, &inner_problem);
  if (check_retval(&retval, "MRIStepInnerAdjointProblem_Create", 1))
  {
    return 1;
  }

  SUNAdjRhsFn inner_adj_rhs = NULL;
  N_Vector inner_sf         = NULL;
  retval = MRIStepInnerAdjointProblem_GetAdjRhsFn(inner_problem, &inner_adj_rhs);
  if (check_retval(&retval, "MRIStepInnerAdjointProblem_GetAdjRhsFn", 1))
  {
    return 1;
  }
  retval = MRIStepInnerAdjointProblem_GetTerminalState(inner_problem, &inner_sf);
  if (check_retval(&retval, "MRIStepInnerAdjointProblem_GetTerminalState", 1))
  {
    return 1;
  }

  /* The inner checkpoints are populated as MRIStep recomputes each fast IVP. */
  SUNAdjointCheckpointScheme inner_scheme = NULL;
  retval = SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM,
                                                   mem_helper, args.check_freq,
                                                   2 * nsteps * (nfast + 1),
                                                   args.keep_checks, sunctx,
                                                   &inner_scheme);
  if (check_retval(&retval, "SUNAdjointCheckpointScheme_Create_Fixed", 1))
  {
    return 1;
  }
  retval = ARKodeSetAdjointCheckpointScheme(inner_mem, inner_scheme);
  if (check_retval(&retval, "ARKodeSetAdjointCheckpointScheme", 1))
  {
    return 1;
  }

  SUNAdjointStepper inner_adjoint = NULL;
  retval = ERKStepCreateAdjointStepper(inner_mem, inner_adj_rhs, tf, inner_sf,
                                       sunctx, &inner_adjoint);
  if (check_retval(&retval, "ERKStepCreateAdjointStepper", 1)) { return 1; }

  SUNAdjointStepper outer_adjoint = NULL;
  retval = MRIStepCreateAdjointStepper(arkode_mem, inner_adjoint, slow_adj_rhs,
                                       NULL, inner_problem, tf, sf, sunctx,
                                       &outer_adjoint);
  if (check_retval(&retval, "MRIStepCreateAdjointStepper", 1)) { return 1; }

  retval = SUNAdjointStepper_Evolve(outer_adjoint, t0, sf, &tret);
  if (check_retval(&retval, "SUNAdjointStepper_Evolve", 1)) { return 1; }

  printf("Adjoint Solution:\n");
  N_VPrint(sf);
  printf("\nSUNAdjointStepper Stats:\n");
  retval = SUNAdjointStepper_PrintAllStats(outer_adjoint, stdout,
                                           SUN_OUTPUTFORMAT_TABLE);
  if (check_retval(&retval, "SUNAdjointStepper_PrintAllStats", 1)) { return 1; }
  printf("\n");

  SUNAdjointStepper_Destroy(&outer_adjoint);
  SUNAdjointStepper_Destroy(&inner_adjoint);
  MRIStepInnerAdjointProblem_Free(inner_problem);
  SUNAdjointCheckpointScheme_Destroy(&inner_scheme);
  SUNAdjointCheckpointScheme_Destroy(&outer_scheme);
  ARKodeFree(&arkode_mem);
  SUNStepper_Destroy(&inner_stepper);
  ARKodeFree(&inner_mem);
  N_VDestroy(sf);
  N_VDestroy(mu);
  N_VDestroy(lambda);
  N_VDestroy(u);
  SUNMemoryHelper_Destroy(mem_helper);
  SUNContext_Free(&sunctx);

  return 0;
}

static int slow_rhs(sunrealtype t, N_Vector uvec, N_Vector udotvec,
                    void* user_data)
{
  sunrealtype* p    = (sunrealtype*)user_data;
  sunrealtype* u    = N_VGetArrayPointer(uvec);
  sunrealtype* udot = N_VGetArrayPointer(udotvec);

  udot[0] = p[0] * u[0] - p[1] * u[0] * u[1];
  udot[1] = SUN_RCONST(0.0);
  return 0;
}

static int fast_rhs(sunrealtype t, N_Vector uvec, N_Vector udotvec,
                    void* user_data)
{
  sunrealtype* p    = (sunrealtype*)user_data;
  sunrealtype* u    = N_VGetArrayPointer(uvec);
  sunrealtype* udot = N_VGetArrayPointer(udotvec);

  udot[0] = SUN_RCONST(0.0);
  udot[1] = -p[2] * u[1] + p[3] * u[0] * u[1];
  return 0;
}

static int slow_adj_rhs(sunrealtype t, N_Vector y, N_Vector sens,
                        N_Vector sens_dot, void* user_data)
{
  sunrealtype* p      = (sunrealtype*)user_data;
  sunrealtype* u      = N_VGetArrayPointer(y);
  sunrealtype* lambda = N_VGetArrayPointer(N_VGetSubvector_ManyVector(sens, 0));
  sunrealtype* lambda_dot =
    N_VGetArrayPointer(N_VGetSubvector_ManyVector(sens_dot, 0));
  sunrealtype* mu_dot =
    N_VGetArrayPointer(N_VGetSubvector_ManyVector(sens_dot, 1));

  lambda_dot[0] = (p[0] - p[1] * u[1]) * lambda[0];
  lambda_dot[1] = -p[1] * u[0] * lambda[0];
  mu_dot[0]     = u[0] * lambda[0];
  mu_dot[1]     = -u[0] * u[1] * lambda[0];
  mu_dot[2]     = SUN_RCONST(0.0);
  mu_dot[3]     = SUN_RCONST(0.0);
  return 0;
}

static int fast_adj_rhs(sunrealtype t, N_Vector y, N_Vector sens,
                        N_Vector sens_dot, void* user_data)
{
  sunrealtype* p      = (sunrealtype*)user_data;
  sunrealtype* u      = N_VGetArrayPointer(y);
  sunrealtype* lambda = N_VGetArrayPointer(N_VGetSubvector_ManyVector(sens, 0));
  sunrealtype* lambda_dot =
    N_VGetArrayPointer(N_VGetSubvector_ManyVector(sens_dot, 0));
  sunrealtype* mu_dot =
    N_VGetArrayPointer(N_VGetSubvector_ManyVector(sens_dot, 1));

  lambda_dot[0] = p[3] * u[1] * lambda[1];
  lambda_dot[1] = (-p[2] + p[3] * u[0]) * lambda[1];
  mu_dot[0]     = SUN_RCONST(0.0);
  mu_dot[1]     = SUN_RCONST(0.0);
  mu_dot[2]     = -u[1] * lambda[1];
  mu_dot[3]     = u[0] * u[1] * lambda[1];
  return 0;
}

static void dgdu(N_Vector uvec, N_Vector dgvec)
{
  sunrealtype* u  = N_VGetArrayPointer(uvec);
  sunrealtype* dg = N_VGetArrayPointer(dgvec);
  dg[0]           = u[0] - SUN_RCONST(1.0);
  dg[1]           = u[1] - SUN_RCONST(1.0);
}

static void dgdp(N_Vector dgvec) { N_VConst(SUN_RCONST(0.0), dgvec); }

static void print_help(int argc, char* argv[], int exit_code)
{
  if (exit_code) { fprintf(stderr, "%s: option not recognized\n", argv[0]); }
  else
  {
    fprintf(stderr, "%s ", argv[0]);
  }
  fprintf(stderr, "options:\n");
  fprintf(stderr, "--tf <real>         final simulation time\n");
  fprintf(stderr, "--hs <real>         slow timestep size\n");
  fprintf(stderr, "--hf <real>         fast timestep size\n");
  fprintf(stderr, "--check-freq <int>  how often to checkpoint (in steps)\n");
  fprintf(stderr, "--dont-keep         don't keep checkpoints after loading\n");
  fprintf(stderr, "--help              print these options\n");
  exit(exit_code);
}

static void parse_args(int argc, char* argv[], ProgramArgs* args)
{
  for (int argi = 1; argi < argc; ++argi)
  {
    const char* arg = argv[argi];
    if (!strcmp(arg, "--tf")) { args->tf = atof(argv[++argi]); }
    else if (!strcmp(arg, "--hs")) { args->hs = atof(argv[++argi]); }
    else if (!strcmp(arg, "--hf")) { args->hf = atof(argv[++argi]); }
    else if (!strcmp(arg, "--check-freq"))
    {
      args->check_freq = atoi(argv[++argi]);
    }
    else if (!strcmp(arg, "--dont-keep")) { args->keep_checks = SUNFALSE; }
    else if (!strcmp(arg, "--help")) { print_help(argc, argv, 0); }
    else
    {
      print_help(argc, argv, 1);
    }
  }
}

static int check_retval(void* retval_ptr, const char* funcname, int opt)
{
  if (opt == 0 && retval_ptr == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }
  if (opt == 1)
  {
    int* retval = (int*)retval_ptr;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return 1;
    }
  }
  return 0;
}
