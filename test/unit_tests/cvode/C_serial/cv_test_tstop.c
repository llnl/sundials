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
 * Unit test for setting stop time
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "cvode/cvode.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_matrix.h"
#include "sundials/sundials_nvector.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#else
#define GSYM "g"
#endif

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

static int ode_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* ydot_data = N_VGetArrayPointer(ydot);
  ydot_data[0]           = ONE;
  return 0;
}

static int ode_jac(sunrealtype t, N_Vector y, N_Vector f, SUNMatrix J,
                   void* user_data, N_Vector tempv1, N_Vector tempv2,
                   N_Vector tempv3)
{
  sunrealtype* J_data = SUNDenseMatrix_Data(J);
  J_data[0]           = ZERO;
  return 0;
}

int main(int argc, char* argv[])
{
  SUNContext sunctx  = NULL;
  N_Vector y         = NULL;
  SUNMatrix A        = NULL;
  SUNLinearSolver LS = NULL;
  void* cvode_mem    = NULL;

  int flag             = 0;
  int cvode_flag       = 0;
  int i                = 0;
  sunrealtype tout     = SUN_RCONST(0.10);
  sunrealtype dt_tout  = SUN_RCONST(0.25);
  sunrealtype tstop    = SUN_RCONST(0.30);
  sunrealtype dt_tstop = SUN_RCONST(0.30);
  sunrealtype tret     = ZERO;
  sunrealtype tcur     = ZERO;

  /* --------------
   * Create context
   * -------------- */

  flag = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (flag)
  {
    fprintf(stderr, "SUNContext_Create returned %i\n", flag);
    return 1;
  }

  /* -----------------------
   * Setup initial condition
   * ----------------------- */

  y = N_VNew_Serial(1, sunctx);
  if (!y) { return 1; }
  N_VConst(ZERO, y);

  /* -----------
   * Setup CVODE
   * ----------- */

  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (!cvode_mem) { return 1; }

  flag = CVodeInit(cvode_mem, ode_rhs, ZERO, y);
  if (flag) { return 1; }

  flag = CVodeSStolerances(cvode_mem, SUN_RCONST(1.0e-4), SUN_RCONST(1.0e-8));
  if (flag) { return 1; }

  A = SUNDenseMatrix(1, 1, sunctx);
  if (!A) { return 1; }

  LS = SUNLinSol_Dense(y, A, sunctx);
  if (!LS) { return 1; }

  flag = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (flag) { return 1; }

  flag = CVodeSetJacFn(cvode_mem, ode_jac);
  if (flag) { return 1; }

  flag = CVodeSetMaxOrd(cvode_mem, 1);
  if (flag) { return 1; }

  flag = CVodeSetStopTime(cvode_mem, tstop);
  if (flag) { return 1; }

  /* ---------------
   * Advance in time
   * --------------- */

  printf("0: tout = %" GSYM ", tstop = %" GSYM ", tret = %" GSYM
         ", tcur = %" GSYM "\n",
         tout, tstop, tret, tcur);

  for (i = 1; i <= 6; i++)
  {
    cvode_flag = CVode(cvode_mem, tout, y, &tret, CV_NORMAL);
    if (cvode_flag < 0)
    {
      flag = 1;
      break;
    }

    flag = CVodeGetCurrentTime(cvode_mem, &tcur);
    if (flag) { break; }

    printf("%i: tout = %" GSYM ", tstop = %" GSYM ", tret = %" GSYM
           ", tcur = %" GSYM ", return = %i\n",
           i, tout, tstop, tret, tcur, cvode_flag);

    /* First return: output time < stop time */
    if (i == 1 && cvode_flag != CV_SUCCESS)
    {
      printf("ERROR: Expected output return!\n");
      flag = 1;
      break;
    }

    /* Second return: output time > stop time */
    if (i == 2)
    {
      if (cvode_flag != CV_TSTOP_RETURN)
      {
        printf("ERROR: Expected stop return!\n");
        flag = 1;
        break;
      }

      /* Update stop time */
      tstop += dt_tstop;
      flag = CVodeSetStopTime(cvode_mem, tstop);
      if (flag) { break; }
    }

    /* Third return: output time = stop time */
    if (i == 3)
    {
      if (cvode_flag != CV_TSTOP_RETURN)
      {
        printf("ERROR: Expected stop return!\n");
        flag = 1;
        break;
      }

      /* Update stop time */
      tstop += dt_tstop;
      flag = CVodeSetStopTime(cvode_mem, tstop);
      if (flag) { break; }
    }

    /* Fourth return: output time < stop time but both output time and the stop
       time are overtaken in the same step */
    if (i == 4)
    {
      if (cvode_flag != CV_SUCCESS)
      {
        printf("ERROR: Expected output return!\n");
        flag = 1;
        break;
      }
    }

    /* Fifth return: output time > stop time after step where both output time
       and the stop time were overtaken in the same step */
    if (i == 5)
    {
      if (cvode_flag != CV_TSTOP_RETURN)
      {
        printf("ERROR: Expected stop return!\n");
        flag = 1;
        break;
      }
    }

    /* Sixth return: output time < stop time (not updated) */
    if (i == 6 && cvode_flag != CV_SUCCESS)
    {
      printf("ERROR: Expected output return!\n");
      flag = 1;
      break;
    }

    /* update output time */
    tout += dt_tout;
  }

  /* --------
   * Clean up
   * -------- */

  CVodeFree(&cvode_mem);
  N_VDestroy(y);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  SUNContext_Free(&sunctx);

  if (!flag) { printf("SUCCESS\n"); }

  return flag;
}

/*---- end of file ----*/
