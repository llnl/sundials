/* ------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * This example solves the Lotka-Volterra ODE with four parameters,
 *
 *     u = [dx/dt] = [ p_0*x - p_1*x*y  ]
 *         [dy/dt]   [ -p_2*y + p_3*x*y ].
 *
 * The initial condition is u(t_0) = 1.0 and we use the parameters
 * p  = [1.5, 1.0, 3.0, 1.0]. The integration interval can be controlled via
 * the --tf command line argument, but by default it is t \in [0, 10.].
 * An explicit Runge--Kutta method is employed via the ARKStep time stepper
 * provided by ARKODE. After solving the forward problem, adjoint sensitivity
 * analysis (ASA) is performed using the discrete adjoint method available with
 * with ARKStep in order to obtain the gradient of the scalar cost function,
 *
 *    g(u(t_f), p) = || 1 - u(t_f, p) ||^2 / 2
 *
 * with respect to the initial condition and the parameters.
 *
 * ./ark_lotka_volterra_adj options:
 * --tf <real>         the final simulation time
 * --dt <real>         the timestep size
 * --order <int>       the order of the RK method
 * --check-freq <int>  how often to checkpoint (in steps)
 * --dont-keep         don't keep checkpoints around after loading
 * --help              print these options
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/sundials_core.h>

#include <nvector/nvector_manyvector.h>
#include <nvector/nvector_serial.h>
#include <sunadjointcheckpointscheme/sunadjointcheckpointscheme_fixed.h>
#include <sundials/sundials_adjointstepper.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmemory/sunmemory_system.h>

#include <arkode/arkode.h>
#include <arkode/arkode_arkstep.h>
#include "sundials/sundials_memory.h"
#include "sundials/sundials_nvector.h"

typedef struct
{
  sunrealtype tf;
  sunrealtype dt;
  int order;
  int check_freq;
  sunbooleantype keep_checks;
} ProgramArgs;

static sunrealtype params[4] = {SUN_RCONST(1.5), SUN_RCONST(1.0),
                                SUN_RCONST(3.0), SUN_RCONST(1.0)};
static void parse_args(int argc, char* argv[], ProgramArgs* args);
static void print_help(int argc, char* argv[], int exit_code);
static int check_retval(void* retval_ptr, const char* funcname, int opt);
static int lotka_volterra(sunrealtype t, N_Vector uvec, N_Vector udotvec,
                          void* user_data);
static int vjp(N_Vector vvec, N_Vector Jvvec, sunrealtype t, N_Vector uvec,
               N_Vector udotvec, void* user_data, N_Vector tmp);
static int parameter_vjp(N_Vector vvec, N_Vector Jvvec, sunrealtype t,
                         N_Vector uvec, N_Vector udotvec, void* user_data,
                         N_Vector tmp);
static int adj_rhs(sunrealtype t, N_Vector y, N_Vector sens, N_Vector sens_dot,
                   void* user_data);
static void dgdu(N_Vector uvec, N_Vector dgvec, const sunrealtype* p);
static void dgdp(N_Vector uvec, N_Vector dgvec, const sunrealtype* p);

int main(int argc, char* argv[])
{
  int retval        = 0;
  SUNContext sunctx = NULL;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  ProgramArgs args;
  args.tf          = SUN_RCONST(10.0);
  args.dt          = SUN_RCONST(1e-3);
  args.order       = 4;
  args.keep_checks = SUNTRUE;
  args.check_freq  = 2;
  parse_args(argc, argv, &args);

  //
  // Create the initial conditions vector
  //

  sunindextype neq = 2;
  N_Vector u       = N_VNew_Serial(neq, sunctx);
  N_Vector u0      = N_VClone(u);
  N_VConst(SUN_RCONST(1.0), u0);
  N_VConst(SUN_RCONST(1.0), u);

  //
  // Create the ARKODE stepper that will be used for the forward evolution.
  //

  const sunrealtype dt = args.dt;
  sunrealtype t0       = SUN_RCONST(0.0);
  sunrealtype tf       = args.tf;
  const int nsteps     = (int)ceil(((tf - t0) / dt));
  const int order      = args.order;
  void* arkode_mem     = ARKStepCreate(lotka_volterra, NULL, t0, u, sunctx);

  retval = ARKodeSetOrder(arkode_mem, order);
  if (check_retval(&retval, "ARKodeSetOrder", 1)) { return 1; }

  // Due to roundoff in the `t` accumulation within the integrator,
  // the integrator may actually use nsteps + 1 time steps to reach tf.
  retval = ARKodeSetMaxNumSteps(arkode_mem, nsteps + 1);
  if (check_retval(&retval, "ARKodeSetMaxNumSteps", 1)) { return 1; }

  // Enable checkpointing during the forward solution.
  const int check_interval                     = args.check_freq;
  const int ncheck                             = nsteps * order;
  const sunbooleantype keep_check              = args.keep_checks;
  SUNAdjointCheckpointScheme checkpoint_scheme = NULL;
  SUNMemoryHelper mem_helper                   = SUNMemoryHelper_Sys(sunctx);

  retval = SUNAdjointCheckpointScheme_Create_Fixed(SUNDATAIOMODE_INMEM,
                                                   mem_helper, check_interval,
                                                   ncheck, keep_check, sunctx,
                                                   &checkpoint_scheme);
  if (check_retval(&retval, "SUNAdjointCheckpointScheme_Create_Fixed", 1))
  {
    return 1;
  }

  retval = ARKodeSetAdjointCheckpointScheme(arkode_mem, checkpoint_scheme);
  if (check_retval(&retval, "ARKodeSetAdjointCheckpointScheme", 1))
  {
    return 1;
  }

  //
  // Compute the forward solution
  //

  printf("Initial condition:\n");
  N_VPrint(u);

  retval = ARKodeSetUserData(arkode_mem, (void*)params);
  if (check_retval(&retval, "ARKodeSetUserData", 1)) { return 1; }

  retval = ARKodeSetFixedStep(arkode_mem, dt);
  if (check_retval(&retval, "ARKodeSetFixedStep", 1)) { return 1; }

  sunrealtype tret = t0;

  retval = ARKodeEvolve(arkode_mem, tf, u, &tret, ARK_NORMAL);
  if (check_retval(&retval, "ARKodeEvolve", 1)) { return 1; }

  printf("Forward Solution:\n");
  N_VPrint(u);

  printf("ARKODE Stats for Forward Solution:\n");
  retval = ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_retval(&retval, "ARKodePrintAllStats", 1)) { return 1; }
  printf("\n");

  //
  // Create the adjoint stepper
  //

  sunindextype num_params = 4;
  N_Vector sensu0         = N_VClone(u);
  N_Vector sensp          = N_VNew_Serial(num_params, sunctx);
  N_Vector sens[2]        = {sensu0, sensp};
  N_Vector sf             = N_VNew_ManyVector(2, sens, sunctx);

  // Set the terminal condition for the adjoint system, which
  // should be the the gradient of our cost function at tf.
  dgdu(u, sensu0, params);
  dgdp(u, sensp, params);

  printf("Adjoint terminal condition:\n");
  N_VPrint(sf);

  SUNAdjointStepper adj_stepper;
  retval = ARKStepCreateAdjointStepper(arkode_mem, adj_rhs, NULL, tf, sf,
                                       sunctx, &adj_stepper);
  if (check_retval(&retval, "ARKStepCreateAdjointStepper", 1)) { return 1; }

  //
  // Now compute the adjoint solution
  //

  retval = SUNAdjointStepper_Evolve(adj_stepper, t0, sf, &tret);
  if (check_retval(&retval, "SUNAdjointStepper_Evolve", 1)) { return 1; }

  printf("Adjoint Solution:\n");
  N_VPrint(sf);

  printf("\nSUNAdjointStepper Stats:\n");
  retval = SUNAdjointStepper_PrintAllStats(adj_stepper, stdout,
                                           SUN_OUTPUTFORMAT_TABLE);
  if (check_retval(&retval, "SUNAdjointStepper_PrintAllStats", 1)) { return 1; }
  printf("\n");

  //
  // Cleanup
  //

  N_VDestroy(sensu0);
  N_VDestroy(sensp);
  N_VDestroy(sf);
  N_VDestroy(u);
  N_VDestroy(u0);
  SUNAdjointCheckpointScheme_Destroy(&checkpoint_scheme);
  SUNAdjointStepper_Destroy(&adj_stepper);
  ARKodeFree(&arkode_mem);
  SUNMemoryHelper_Destroy(mem_helper);
  SUNContext_Free(&sunctx);

  return 0;
}

int lotka_volterra(sunrealtype t, N_Vector uvec, N_Vector udotvec, void* user_data)
{
  sunrealtype* p    = (sunrealtype*)user_data;
  sunrealtype* u    = N_VGetArrayPointer(uvec);
  sunrealtype* udot = N_VGetArrayPointer(udotvec);

  udot[0] = p[0] * u[0] - p[1] * u[0] * u[1];
  udot[1] = -p[2] * u[1] + p[3] * u[0] * u[1];

  return 0;
}

int vjp(N_Vector vvec, N_Vector Jvvec, sunrealtype t, N_Vector uvec,
        N_Vector udotvec, void* user_data, N_Vector tmp)
{
  sunrealtype* p  = (sunrealtype*)user_data;
  sunrealtype* u  = N_VGetArrayPointer(uvec);
  sunrealtype* v  = N_VGetArrayPointer(vvec);
  sunrealtype* Jv = N_VGetArrayPointer(Jvvec);

  Jv[0] = (p[0] - p[1] * u[1]) * v[0] + p[3] * u[1] * v[1];
  Jv[1] = -p[1] * u[0] * v[0] + (-p[2] + p[3] * u[0]) * v[1];

  return 0;
}

int parameter_vjp(N_Vector vvec, N_Vector Jvvec, sunrealtype t, N_Vector uvec,
                  N_Vector udotvec, void* user_data, N_Vector tmp)
{
  if (user_data != params) { return -1; }

  sunrealtype* u  = N_VGetArrayPointer(uvec);
  sunrealtype* v  = N_VGetArrayPointer(vvec);
  sunrealtype* Jv = N_VGetArrayPointer(Jvvec);

  Jv[0] = u[0] * v[0];
  Jv[1] = -u[0] * u[1] * v[0];
  Jv[2] = -u[1] * v[1];
  Jv[3] = u[0] * u[1] * v[1];

  return 0;
}

static int adj_rhs(sunrealtype t, N_Vector y, N_Vector sens, N_Vector sens_dot,
                   void* user_data)
{
  N_Vector Lambda_part = N_VGetSubvector_ManyVector(sens, 0);
  N_Vector Lambda      = N_VGetSubvector_ManyVector(sens_dot, 0);
  N_Vector nu          = N_VGetSubvector_ManyVector(sens_dot, 1);
  vjp(Lambda_part, Lambda, t, y, NULL, user_data, NULL);
  parameter_vjp(Lambda_part, nu, t, y, NULL, user_data, NULL);
  return 0;
}

void dgdu(N_Vector uvec, N_Vector dgvec, const sunrealtype* p)
{
  sunrealtype* u  = N_VGetArrayPointer(uvec);
  sunrealtype* dg = N_VGetArrayPointer(dgvec);

  dg[0] = -SUN_RCONST(1.0) + u[0];
  dg[1] = -SUN_RCONST(1.0) + u[1];
}

void dgdp(N_Vector uvec, N_Vector dgvec, const sunrealtype* p)
{
  sunrealtype* dg = N_VGetArrayPointer(dgvec);

  dg[0] = SUN_RCONST(0.0);
  dg[1] = SUN_RCONST(0.0);
  dg[2] = SUN_RCONST(0.0);
  dg[3] = SUN_RCONST(0.0);
}

void print_help(int argc, char* argv[], int exit_code)
{
  if (exit_code) { fprintf(stderr, "%s: option not recognized\n", argv[0]); }
  else { fprintf(stderr, "%s ", argv[0]); }
  fprintf(stderr, "options:\n");
  fprintf(stderr, "--tf <real>         the final simulation time\n");
  fprintf(stderr, "--dt <real>         the timestep size\n");
  fprintf(stderr, "--order <int>       the order of the RK method\n");
  fprintf(stderr, "--check-freq <int>  how often to checkpoint (in steps)\n");
  fprintf(stderr,
          "--dont-keep         don't keep checkpoints around after loading\n");
  fprintf(stderr, "--help              print these options\n");
  exit(exit_code);
}

void parse_args(int argc, char* argv[], ProgramArgs* args)
{
  for (int argi = 1; argi < argc; ++argi)
  {
    const char* arg = argv[argi];
    if (!strcmp(arg, "--tf")) { args->tf = atof(argv[++argi]); }
    else if (!strcmp(arg, "--dt")) { args->dt = atof(argv[++argi]); }
    else if (!strcmp(arg, "--order")) { args->order = atoi(argv[++argi]); }
    else if (!strcmp(arg, "--check-freq"))
    {
      args->check_freq = atoi(argv[++argi]);
    }
    else if (!strcmp(arg, "--dont-keep")) { args->keep_checks = SUNFALSE; }
    else if (!strcmp(arg, "--help")) { print_help(argc, argv, 0); }
    else { print_help(argc, argv, 1); }
  }
}

int check_retval(void* retval_ptr, const char* funcname, int opt)
{
  int* retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && retval_ptr == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)retval_ptr;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return 1;
    }
  }

  return (0);
}
