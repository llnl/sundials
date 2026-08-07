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
 * Utilities for the arkode_damped_harmonic_symplectic example.
 * ---------------------------------------------------------------------------*/

#ifndef ARK_DAMPED_HARMONIC_SYMPLECTIC_H
#define ARK_DAMPED_HARMONIC_SYMPLECTIC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_types.h>

#define PI SUN_RCONST(3.14159265358979323846264338327950)

typedef struct
{
  int order;
  int num_output_times;
  int use_compsums;
  int use_tstop;
  sunrealtype Tf;
  sunrealtype dt;
} ProgramArgs;

static void PrintHelp(void)
{
  fprintf(stderr,
          "ark_damped_harmonic_symplectic: an ARKODE example demonstrating "
          "the SPRKStep time-stepping module solving a time-dependent "
          "damped harmonic oscillator\n");
  /* clang-format off */
  fprintf(stderr, "  --order <int>               the order of the method to use (default 4)\n");
  fprintf(stderr, "  --dt <Real>                 the fixed-time step size to use (default 0.01)\n");
  fprintf(stderr, "  --nout <int>                the number of output times (default 100)\n");
  fprintf(stderr, "  --use-compensated-sums      turns on compensated summation in ARKODE where applicable\n");
  fprintf(stderr, "  --disable-tstop             turns off tstop mode\n");
  /* clang-format on */
}

static int ParseArgs(int argc, char* argv[], ProgramArgs* args)
{
  int argi = 0;

  args->order            = 4;
  args->num_output_times = 8;
  args->use_compsums     = 0;
  args->use_tstop        = 1;
  args->Tf               = SUN_RCONST(10.0) * PI;
  args->dt               = SUN_RCONST(1e-3);

  for (argi = 1; argi < argc; argi++)
  {
    if (!strcmp(argv[argi], "--order"))
    {
      argi++;
      args->order = atoi(argv[argi]);
    }
    else if (!strcmp(argv[argi], "--tf"))
    {
      argi++;
      args->Tf = atof(argv[argi]);
    }
    else if (!strcmp(argv[argi], "--dt"))
    {
      argi++;
      args->dt = atof(argv[argi]);
    }
    else if (!strcmp(argv[argi], "--nout"))
    {
      argi++;
      args->num_output_times = atoi(argv[argi]);
    }
    else if (!strcmp(argv[argi], "--use-compensated-sums"))
    {
      args->use_compsums = 1;
    }
    else if (!strcmp(argv[argi], "--disable-tstop")) { args->use_tstop = 0; }
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

  return 0;
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

#endif /* ARK_DAMPED_HARMONIC_SYMPLECTIC_H */
