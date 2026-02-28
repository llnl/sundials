/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * Unit test for GetUserData functions
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "cvode/cvode.h"
#include "nvector/nvector_serial.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/* Dummy user-supplied function */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  return 0;
}

/* Main program */
int main(int argc, char* argv[])
{
  int retval        = 0;
  SUNContext sunctx = NULL;
  N_Vector y        = NULL;
  void* cvode_mem   = NULL;
  int udata_in      = 1;
  void* udata_out   = NULL;

  /* Create the SUNDIALS context object for this simulation. */
  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (retval)
  {
    fprintf(stderr, "SUNContext_Create returned %i\n", retval);
    return 1;
  }

  /* Create solution vector and initialize to zero */
  y = N_VNew_Serial(1, sunctx);
  if (!y)
  {
    fprintf(stderr, "N_VNew_Serial returned NULL\n");
    return 1;
  }
  N_VConst(ONE, y);

  /* Create CVODE mem structure */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (!cvode_mem)
  {
    fprintf(stderr, "CVodeCreate returned NULL\n");
    return 1;
  }

  retval = CVodeInit(cvode_mem, f, ZERO, y);
  if (retval)
  {
    fprintf(stderr, "CVodeInit returned %i\n", retval);
    return 1;
  }

  /* Set user data */
  retval = CVodeSetUserData(cvode_mem, &udata_in);
  if (retval)
  {
    fprintf(stderr, "CVodeSetUserData returned %i\n", retval);
    return 1;
  }

  /* Get user data */
  retval = CVodeGetUserData(cvode_mem, &udata_out);
  if (retval)
  {
    fprintf(stderr, "CVodeGetUserData returned %i\n", retval);
    return 1;
  }

  if (!udata_out)
  {
    fprintf(stderr, "udata_out is NULL\n");
    return 1;
  }

  if (&udata_in != (int*)udata_out)
  {
    fprintf(stderr, "udata_in != udata_out\n");
    return 1;
  }

  /* Clean up */
  CVodeFree(&cvode_mem);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  printf("SUCCESS\n");

  return 0;
}

/*---- end of file ----*/
