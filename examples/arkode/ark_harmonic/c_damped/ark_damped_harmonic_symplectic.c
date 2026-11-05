/* clang-format off */
/* ----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * ----------------------------------------------------------------------------
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
 * ----------------------------------------------------------------------------
 * In this example we consider the time-dependent damped harmonic oscillator
 *    q'(t) = p(t) exp(-F(t))
 *    p'(t) = -(F(t) * p + omega^2(t) * q)
 * With the initial conditions q(0) = 1, p(0) = 0.
 * The Hamiltonian for the system is
 *    H(p,q,t) = (p^2 * exp(-F(t)))/2 + (omega^2(t) * q^2 * exp(F(t)))/2
 * where omega(t) = cos(t/2), F(t) = 0.018*sin(t/pi).
 * We simulate the problem on t = [0, 30] using the symplectic methods in
 * SPRKStep.
 *
 * This is example 7.2 from:
 * Struckmeier, J., & Riedel, C. (2002). Canonical transformations and exact
 * invariants for time‐dependent Hamiltonian systems. Annalen der Physik, 11(1),
 * 15-38.
 *
 * The example has the following command line arguments:
 *   --order <int>               the order of the method to use (default 4)
 *   --dt <Real>                 the fixed-time step size to use (default 0.01)
 *   --nout <int>                the number of output times (default 100)
 *   --disable-tstop             turns off tstop mode
 *   --use-compensated-sums      turns on compensated summation in ARKODE where
 *                               applicable
 * --------------------------------------------------------------------------*/
/* clang-format on */

#include "ark_damped_harmonic_symplectic.h"

#include <arkode/arkode.h>
#include <arkode/arkode_sprk.h>
#include <arkode/arkode_sprkstep.h> /* prototypes for SPRKStep fcts., consts */
#include <math.h>
#include <nvector/nvector_serial.h> /* serial N_Vector type, fcts., macros  */
#include <stdio.h>
#include <string.h>
#include <sundials/sundials_math.h> /* def. math fcns, 'sunrealtype'           */
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

/* RHS functions */
static int pdot(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int qdot(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/* Other functions */
static sunrealtype Hamiltonian(N_Vector yvec, sunrealtype t);

int main(int argc, char* argv[])
{
  ProgramArgs args;
  SUNContext sunctx    = NULL;
  N_Vector y           = NULL;
  sunrealtype* ydata   = NULL;
  void* arkode_mem     = NULL;
  int iout             = 0;
  int retval           = 0;
  int order            = 0;
  int use_compsums     = 0;
  int num_output_times = 0;
  sunrealtype Tf       = SUN_RCONST(0.0);
  sunrealtype dt       = SUN_RCONST(0.0);
  sunrealtype dTout    = SUN_RCONST(0.0);
  const sunrealtype T0 = SUN_RCONST(0.0);

  /* Parse the command line arguments */
  if (ParseArgs(argc, argv, &args)) { return 1; };

  /* Default integrator options */
  order            = args.order;
  use_compsums     = args.use_compsums;
  num_output_times = args.num_output_times;

  /* Default problem parameters */
  Tf    = args.Tf;
  dt    = args.dt;
  dTout = (Tf - T0) / ((sunrealtype)num_output_times);

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return 1; }

  printf("\n   Begin time-dependent damped harmonic oscillator problem\n\n");

  /* Allocate our state vector */
  y = N_VNew_Serial(2, sunctx);

  /* Fill the initial conditions */
  ydata    = N_VGetArrayPointer(y);
  ydata[0] = 0; /* \dot{q} = p */
  ydata[1] = 1; /* \ddot{q} = \dot{p} */

  /* Create SPRKStep integrator */
  arkode_mem = SPRKStepCreate(qdot, pdot, T0, y, sunctx);

  retval = ARKodeSetOrder(arkode_mem, order);
  if (check_retval(&retval, "ARKodeSetOrder", 1)) { return 1; }

  retval = ARKodeSetUseCompensatedSums(arkode_mem, use_compsums);
  if (check_retval(&retval, "ARKodeSetUseCompensatedSums", 1)) { return 1; }

  retval = ARKodeSetFixedStep(arkode_mem, dt);
  if (check_retval(&retval, "ARKodeSetFixedStep", 1)) { return 1; }

  retval = ARKodeSetMaxNumSteps(arkode_mem, ((long int)ceil(Tf / dt)) + 2);
  if (check_retval(&retval, "ARKodeSetMaxNumSteps", 1)) { return 1; }

  /* Print out starting Hamiltonian before integrating */
  sunrealtype tret = T0;
  sunrealtype tout = T0 + dTout;
  /* Output current integration status */
  fprintf(stdout, "t = %.6Lf, q(t) = %.6Lf, H = %.6Lf\n", (long double)tret,
          (long double)ydata[1], (long double)Hamiltonian(y, tret));

  /* Do integration */
  for (iout = 0; iout < num_output_times; iout++)
  {
    if (args.use_tstop) { ARKodeSetStopTime(arkode_mem, tout); }
    retval = ARKodeEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);

    /* Output current integration status */
    fprintf(stdout, "t = %.6Lf, q(t) = %.6Lf, H = %.6Lf\n", (long double)tret,
            (long double)ydata[1], (long double)Hamiltonian(y, tret));

    /* Check if the solve was successful, if so, update the time and continue */
    if (retval >= 0)
    {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
    else
    {
      fprintf(stderr, "Solver failure, stopping integration\n");
      break;
    }
  }

  fprintf(stdout, "\n");
  ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  N_VDestroy(y);
  ARKodeFree(&arkode_mem);
  SUNContext_Free(&sunctx);

  return 0;
}

static sunrealtype omega(sunrealtype t) { return cos(t / SUN_RCONST(2.0)); }

static sunrealtype F(sunrealtype t) { return SUN_RCONST(0.018) * sin(t / PI); }

sunrealtype Hamiltonian(N_Vector yvec, sunrealtype t)
{
  sunrealtype H       = SUN_RCONST(0.0);
  sunrealtype* y      = N_VGetArrayPointer(yvec);
  const sunrealtype p = y[0];
  const sunrealtype q = y[1];

  H = (p * p * exp(-F(t))) / SUN_RCONST(2.0) +
      (omega(t) * omega(t) * q * q * exp(F(t))) / SUN_RCONST(2.0);

  return H;
}

int qdot(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  sunrealtype* y      = N_VGetArrayPointer(yvec);
  sunrealtype* ydot   = N_VGetArrayPointer(ydotvec);
  const sunrealtype p = y[0];

  ydot[1] = p * exp(-F(t));

  return 0;
}

int pdot(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  sunrealtype* y      = N_VGetArrayPointer(yvec);
  sunrealtype* ydot   = N_VGetArrayPointer(ydotvec);
  const sunrealtype p = y[0];
  const sunrealtype q = y[1];

  ydot[0] = -(F(t) * p + omega(t) * omega(t) * q);

  return 0;
}
