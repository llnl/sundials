/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * We solve the classic Van der Pol problem:
 *   y'' - mu*(1 - y^2)*y' + y = 0,  y(0) = 2,  y'(0) = 0.
 * This second-order ODE is converted to a first-order system by defining
 *   y0 = y,  y1 = y'
 * giving
 *   y0' = y1
 *   y1' = mu*(1 - y0^2)*y1 - y0.
 * We use the SUNNonlinearSolver_Auto module to solve the implicit
 * system. This solver automatically switches between modified Newton
 * iteration and fixed-point iteration using a stiffness metric.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cvode/cvode.h>
#include <sundials/sundials_core.h>

#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunnonlinsol/sunnonlinsol_auto.h>

/* Problem constants */
#define NEQ   2
#define T0    SUN_RCONST(0.0)
#define TF    SUN_RCONST(250.0)
#define DTOUT SUN_RCONST(10.0)

typedef struct
{
  sunrealtype mu;
} UserData;

static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  UserData* ud         = (UserData*)user_data;
  const sunrealtype mu = ud->mu;

  sunrealtype* ydata    = N_VGetArrayPointer(y);
  sunrealtype* ydotdata = N_VGetArrayPointer(ydot);

  const sunrealtype y0 = ydata[0];
  const sunrealtype y1 = ydata[1];

  ydotdata[0] = y1;
  ydotdata[1] = mu * (SUN_RCONST(1.0) - y0 * y0) * y1 - y0;

  return 0;
}

static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData* ud         = (UserData*)user_data;
  const sunrealtype mu = ud->mu;
  sunrealtype* ydata   = N_VGetArrayPointer(y);
  const sunrealtype y0 = ydata[0];
  const sunrealtype y1 = ydata[1];

  SM_ELEMENT_D(J, 0, 0) = SUN_RCONST(0.0);
  SM_ELEMENT_D(J, 0, 1) = SUN_RCONST(1.0);
  SM_ELEMENT_D(J, 1, 0) = -SUN_RCONST(2.0) * mu * y0 * y1 - SUN_RCONST(1.0);
  SM_ELEMENT_D(J, 1, 1) = mu * (SUN_RCONST(1.0) - y0 * y0);

  return 0;
}

static int check_retval(const void* retval, const char* funcname, int opt)
{
  if (opt == 0 && retval == NULL)
  {
    fprintf(stderr, "ERROR: %s() returned NULL\n", funcname);
    return 1;
  }
  if (opt == 1)
  {
    const int err = *((const int*)retval);
    if (err < 0)
    {
      fprintf(stderr, "ERROR: %s() returned %d\n", funcname, err);
      return 1;
    }
  }
  return 0;
}

int main(int argc, char* argv[])
{
  int retval;
  SUNContext sunctx;

  /* Problem setup */
  UserData user_data;
  user_data.mu = SUN_RCONST(100.0);

  const sunrealtype y10    = SUN_RCONST(2.0);
  const sunrealtype y20    = SUN_RCONST(0.0);
  const sunrealtype reltol = SUN_RCONST(1.0e-4);
  const sunrealtype abstol = SUN_RCONST(1.0e-4);

  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return 1; }

  N_Vector y = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)y, "N_VNew_Serial", 0)) { return 1; }

  sunrealtype* ydata = N_VGetArrayPointer(y);
  ydata[0]           = y10;
  ydata[1]           = y20;

  void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval(cvode_mem, "CVodeCreate", 0)) { return 1; }

  retval = CVodeInit(cvode_mem, f, T0, y);
  if (check_retval(&retval, "CVodeInit", 1)) { return 1; }

  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSStolerances", 1)) { return 1; }

  retval = CVodeSetUserData(cvode_mem, &user_data);
  if (check_retval(&retval, "CVodeSetUserData", 1)) { return 1; }

  retval = CVodeSetMaxNumSteps(cvode_mem, 10000);
  if (check_retval(&retval, "CVodeSetMaxNumSteps", 1)) { return 1; }

  /* Create nonlinear solver (auto) */
  SUNNonlinearSolver NLS = SUNNonlinSol_Auto(y, 0, SUNNONLINSOL_AUTO_NEWTON,
                                             sunctx);
  if (check_retval((void*)NLS, "SUNNonlinSol_Auto", 0)) { return 1; }

  retval = CVodeSetNonlinearSolver(cvode_mem, NLS);
  if (check_retval(&retval, "CVodeSetNonlinearSolver", 1)) { return 1; }

  /* Provide Newton linear solver and Jacobian for when Newton is active */
  SUNMatrix A        = SUNDenseMatrix(NEQ, NEQ, sunctx);
  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_retval((void*)A, "SUNDenseMatrix", 0)) { return 1; }
  if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) { return 1; }

  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) { return 1; }

  retval = CVodeSetJacFn(cvode_mem, Jac);
  if (check_retval(&retval, "CVodeSetJacFn", 1)) { return 1; }

  /* Parse any remaining command line arguments */
  retval = CVodeSetOptions(cvode_mem, "", "", argc, argv);
  if (check_retval(&retval, "CVodeSetOptions", 1)) { return 1; }

  printf("\nVan der Pol oscillator (CVODE):\n");
  printf("    initial conditions: y1 = %.6f, y2 = %.6f\n", (double)y10,
         (double)y20);
  printf("    mu = %.6f\n", (double)user_data.mu);
  printf("    reltol = %.2e, abstol = %.2e\n\n", (double)reltol, (double)abstol);
  printf("        t           y1           y2\n");
  printf("   -----------------------------------\n");
  printf("  %10.6f  %10.6f  %10.6f\n", (double)T0, (double)ydata[0],
         (double)ydata[1]);

  const int Nt     = (int)SUNRceil(TF / DTOUT);
  sunrealtype tout = T0 + DTOUT;
  for (int iout = 0; iout < Nt; iout++)
  {
    sunrealtype tret;
    retval = CVode(cvode_mem, tout, y, &tret, CV_NORMAL);
    printf("  %10.6f  %10.6f  %10.6f\n", (double)tret, (double)ydata[0],
           (double)ydata[1]);

    if (retval == CV_SUCCESS)
    {
      tout += DTOUT;
      tout = SUNMIN(tout, TF);
    }
    else
    {
      printf("Solver failure, stopping integration\n");
      break;
    }
  }
  printf("   -----------------------------------\n");

  retval = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_retval(&retval, "CVodePrintAllStats", 1)) { return 1; }

  {
    long int nfp, nnewt;
    retval = SUNNonlinSolGetTotalNumItersByType_Auto(NLS, &nfp, &nnewt);
    if (check_retval(&retval, "SUNNonlinSolGetTotalNumItersByType_Auto", 1))
    {
      return 1;
    }
    printf("   Auto nonlinear solver iteration totals: newton = %ld, "
           "fixed-point = %ld\n",
           nnewt, nfp);
  }

  CVodeFree(&cvode_mem);
  SUNNonlinSolFree(NLS);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  return 0;
}
