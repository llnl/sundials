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

#include "ida/ida.h"
#include "nvector/nvector_serial.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/* Dummy user-supplied function */
static int r(sunrealtype t, N_Vector y, N_Vector ydot, N_Vector res,
             void* user_data)
{
  return 0;
}

/* Main program */
int main(int argc, char* argv[])
{
  int retval        = 0;
  SUNContext sunctx = NULL;
  N_Vector y        = NULL;
  N_Vector yp       = NULL;
  void* ida_mem     = NULL;
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

  yp = N_VClone(y);
  if (!yp)
  {
    fprintf(stderr, "N_VClone returned NULL\n");
    return 1;
  }
  N_VConst(ONE, yp);

  /* Create IDA mem structure */
  ida_mem = IDACreate(sunctx);
  if (!ida_mem)
  {
    fprintf(stderr, "IDACreate returned NULL\n");
    return 1;
  }

  retval = IDAInit(ida_mem, r, ZERO, y, yp);
  if (retval)
  {
    fprintf(stderr, "IDAInit returned %i\n", retval);
    return 1;
  }

  /* Set user data */
  retval = IDASetUserData(ida_mem, &udata_in);
  if (retval)
  {
    fprintf(stderr, "IDASetUserData returned %i\n", retval);
    return 1;
  }

  /* Get user data */
  retval = IDAGetUserData(ida_mem, &udata_out);
  if (retval)
  {
    fprintf(stderr, "IDAGetUserData returned %i\n", retval);
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
  IDAFree(&ida_mem);
  N_VDestroy(y);
  N_VDestroy(yp);
  SUNContext_Free(&sunctx);

  printf("SUCCESS\n");

  return 0;
}

/*---- end of file ----*/
