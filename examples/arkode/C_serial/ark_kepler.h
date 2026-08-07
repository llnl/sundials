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
 * Utilities for the arkode_kepler example.
 * ---------------------------------------------------------------------------*/

#ifndef ARK_KEPLER_H
#define ARK_KEPLER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

typedef struct
{
  int step_mode;
  int stepper;
  int num_output_times;
  int use_compsums;
  int use_tstop;
  int count_orbits;
  int check_order;
  sunrealtype dt;
  sunrealtype tf;
  const char* method_name;
} ProgramArgs;

static int ComputeConvergence(int num_dt, sunrealtype* orders,
                              sunrealtype expected_order, sunrealtype a11,
                              sunrealtype a12, sunrealtype a21, sunrealtype a22,
                              sunrealtype b1, sunrealtype b2,
                              sunrealtype* ord_avg, sunrealtype* ord_max,
                              sunrealtype* ord_est)
{
  /* Compute/print overall estimated convergence rate */
  int i           = 0;
  sunrealtype det = 0;
  *ord_avg = 0, *ord_max = 0, *ord_est = 0;
  for (i = 1; i < num_dt; i++)
  {
    *ord_avg += orders[i - 1];
    *ord_max = SUNMAX(*ord_max, orders[i - 1]);
  }
  *ord_avg = *ord_avg / ((sunrealtype)num_dt - 1);
  det      = a11 * a22 - a12 * a21;
  *ord_est = (a11 * b2 - a21 * b1) / det;
  return 0;
}

static void PrintHelp(void)
{
  fprintf(stderr, "ark_kepler: an ARKODE example demonstrating the SPRKStep "
                  "time-stepping module solving the Kepler problem\n");
  /* clang-format off */
  fprintf(stderr, "  --step-mode <fixed, adapt>  should we use a fixed time-step or adaptive time-step (default fixed)\n");
  fprintf(stderr, "  --stepper <SPRK, ERK>       should we use SPRKStep or ARKStep with an ERK method (default SPRK)\n");
  fprintf(stderr, "  --method <string>           which method to use (default ARKODE_SPRK_MCLACHLAN_4_4)\n");
  fprintf(stderr, "  --use-compensated-sums      turns on compensated summation in ARKODE where applicable\n");
  fprintf(stderr, "  --disable-tstop             turns off tstop mode\n");
  fprintf(stderr, "  --dt <Real>                 the fixed-time step size to use if fixed time stepping is turned on (default 0.01)\n");
  fprintf(stderr, "  --tf <Real>                 the final time for the simulation (default 100)\n");
  fprintf(stderr, "  --nout <int>                the number of output times (default 100)\n");
  fprintf(stderr, "  --count-orbits              use rootfinding to count the number of completed orbits\n");
  fprintf(stderr, "  --check-order               compute the order of the method used and check if it is within range of the expected\n");
  /* clang-format on */
}

static int ParseArgs(int argc, char* argv[], ProgramArgs* args)
{
  int argi = 0;

  args->step_mode        = 0;
  args->stepper          = 0;
  args->method_name      = NULL;
  args->count_orbits     = 0;
  args->use_compsums     = 0;
  args->use_tstop        = 1;
  args->dt               = SUN_RCONST(1e-2);
  args->tf               = SUN_RCONST(100.);
  args->check_order      = 0;
  args->num_output_times = 50;

  for (argi = 1; argi < argc; argi++)
  {
    if (!strcmp(argv[argi], "--step-mode"))
    {
      argi++;
      if (!strcmp(argv[argi], "fixed")) { args->step_mode = 0; }
      else if (!strcmp(argv[argi], "adapt")) { args->step_mode = 1; }
      else
      {
        fprintf(stderr, "ERROR: --step-mode must be 'fixed' or 'adapt'\n");
        return 1;
      }
    }
    else if (!strcmp(argv[argi], "--stepper"))
    {
      argi++;
      if (!strcmp(argv[argi], "SPRK")) { args->stepper = 0; }
      else if (!strcmp(argv[argi], "ERK")) { args->stepper = 1; }
      else
      {
        fprintf(stderr, "ERROR: --stepper must be 'SPRK' or 'ERK'\n");
        return 1;
      }
    }
    else if (!strcmp(argv[argi], "--method"))
    {
      argi++;
      args->method_name = argv[argi];
    }
    else if (!strcmp(argv[argi], "--dt"))
    {
      argi++;
      args->dt = atof(argv[argi]);
    }
    else if (!strcmp(argv[argi], "--tf"))
    {
      argi++;
      args->tf = atof(argv[argi]);
    }
    else if (!strcmp(argv[argi], "--nout"))
    {
      argi++;
      args->num_output_times = atoi(argv[argi]);
    }
    else if (!strcmp(argv[argi], "--count-orbits")) { args->count_orbits = 1; }
    else if (!strcmp(argv[argi], "--disable-tstop")) { args->use_tstop = 0; }
    else if (!strcmp(argv[argi], "--use-compensated-sums"))
    {
      args->use_compsums = 1;
    }
    else if (!strcmp(argv[argi], "--check-order")) { args->check_order = 1; }
    else if (!strcmp(argv[argi], "--help"))
    {
      PrintHelp();
      return 1;
    }
    else
    {
      fprintf(stderr, "ERROR: unrecognized argument %s\n", argv[argi]);
      PrintHelp();
      return 1;
    }
  }

  if (!args->method_name)
  {
    if (args->stepper == 0) { args->method_name = "ARKODE_SPRK_MCLACHLAN_4_4"; }
    else if (args->stepper == 1)
    {
      args->method_name = "ARKODE_ZONNEVELD_5_3_4";
    }
  }

  return 0;
}

static void PrintArgs(ProgramArgs* args)
{
  fprintf(stdout, "Problem Arguments:\n");
  fprintf(stdout, "  stepper:              %d\n", args->stepper);
  fprintf(stdout, "  step mode:            %d\n", args->step_mode);
  fprintf(stdout, "  use tstop:            %d\n", args->use_tstop);
  fprintf(stdout, "  use compensated sums: %d\n", args->use_compsums);
  fprintf(stdout, "  dt:                   %Lg\n", (long double)args->dt);
  fprintf(stdout, "  Tf:                   %Lg\n", (long double)args->tf);
  fprintf(stdout, "  nout:                 %d\n\n", args->num_output_times);
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

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return 1;
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nERROR: %s() failed with retval = %d\n\n", funcname,
              *retval);
      return 1;
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}

#endif /* ARK_KEPLER_H */
