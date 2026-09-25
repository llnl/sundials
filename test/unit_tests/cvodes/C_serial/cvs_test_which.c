/* -----------------------------------------------------------------------------
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
 * Unit test for out of range backward problem identifiers
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "cvodes/cvodes.h"
#include "nvector/nvector_serial.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#define ONE SUN_RCONST(1.0)

/* Dummy forward and backward right-hand side functions */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  NV_Ith_S(ydot, 0) = -NV_Ith_S(y, 0);
  return 0;
}

static int fB(sunrealtype t, N_Vector y, N_Vector yB, N_Vector yBdot,
              void* user_dataB)
{
  NV_Ith_S(yBdot, 0) = -NV_Ith_S(yB, 0);
  return 0;
}

/* Main program */
int main(int argc, char* argv[])
{
  /* Identifiers no backward problem can have. The negative values used to walk
     the backward problem list off the end and dereference NULL. */
  const int bad_which[] = {-1, -5, 99};
  const int nbad        = 3;

  int retval        = 0;
  SUNContext sunctx = NULL;
  N_Vector y        = NULL;
  N_Vector yB       = NULL;
  SUNMatrix A       = NULL;
  SUNMatrix AB      = NULL;
  SUNLinearSolver LS  = NULL;
  SUNLinearSolver LSB = NULL;
  void* cvode_mem   = NULL;
  sunrealtype t     = ONE;
  int ncheck        = 0;
  int which         = 0;
  int i;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx))
  {
    printf("ERROR: SUNContext_Create returned nonzero\n");
    return 1;
  }

  y = N_VNew_Serial(1, sunctx);
  if (!y) { return 1; }
  NV_Ith_S(y, 0) = ONE;

  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (!cvode_mem) { return 1; }

  A  = SUNDenseMatrix(1, 1, sunctx);
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (!A || !LS) { return 1; }

  if (CVodeInit(cvode_mem, f, SUN_RCONST(0.0), y) ||
      CVodeSStolerances(cvode_mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-8)) ||
      CVodeSetLinearSolver(cvode_mem, LS, A) ||
      CVodeAdjInit(cvode_mem, 10, CV_HERMITE))
  {
    printf("ERROR: setting up the forward problem failed\n");
    return 1;
  }

  if (CVodeF(cvode_mem, ONE, y, &t, CV_NORMAL, &ncheck) < 0)
  {
    printf("ERROR: CVodeF failed\n");
    return 1;
  }

  yB = N_VNew_Serial(1, sunctx);
  if (!yB) { return 1; }
  NV_Ith_S(yB, 0) = ONE;

  AB  = SUNDenseMatrix(1, 1, sunctx);
  LSB = SUNLinSol_Dense(yB, AB, sunctx);
  if (!AB || !LSB) { return 1; }

  if (CVodeCreateB(cvode_mem, CV_BDF, &which) ||
      CVodeInitB(cvode_mem, which, fB, ONE, yB) ||
      CVodeSStolerancesB(cvode_mem, which, SUN_RCONST(1.0e-6),
                         SUN_RCONST(1.0e-8)) ||
      CVodeSetLinearSolverB(cvode_mem, which, LSB, AB))
  {
    printf("ERROR: setting up the backward problem failed\n");
    return 1;
  }

  for (i = 0; i < nbad; i++)
  {
    sunrealtype tBret = ONE;

    if (CVodeGetB(cvode_mem, bad_which[i], &tBret, yB) != CV_ILL_INPUT)
    {
      printf("ERROR: CVodeGetB did not reject which = %d\n", bad_which[i]);
      retval++;
    }

    if (CVodeSetMaxNumStepsB(cvode_mem, bad_which[i], 10) != CV_ILL_INPUT)
    {
      printf("ERROR: CVodeSetMaxNumStepsB did not reject which = %d\n",
             bad_which[i]);
      retval++;
    }

    if (CVodeSStolerancesB(cvode_mem, bad_which[i], SUN_RCONST(1.0e-6),
                           SUN_RCONST(1.0e-8)) != CV_ILL_INPUT)
    {
      printf("ERROR: CVodeSStolerancesB did not reject which = %d\n",
             bad_which[i]);
      retval++;
    }

    if (CVodeGetAdjCVodeBmem(cvode_mem, bad_which[i]) != NULL)
    {
      printf("ERROR: CVodeGetAdjCVodeBmem did not reject which = %d\n",
             bad_which[i]);
      retval++;
    }
  }

  /* The valid identifier must still work */
  if (CVodeGetB(cvode_mem, which, &t, yB) != CV_SUCCESS)
  {
    printf("ERROR: CVodeGetB rejected the valid which = %d\n", which);
    retval++;
  }

  N_VDestroy(y);
  N_VDestroy(yB);
  SUNLinSolFree(LS);
  SUNLinSolFree(LSB);
  SUNMatDestroy(A);
  SUNMatDestroy(AB);
  CVodeFree(&cvode_mem);
  SUNContext_Free(&sunctx);

  if (retval == 0) { printf("SUCCESS\n"); }

  return retval;
}
