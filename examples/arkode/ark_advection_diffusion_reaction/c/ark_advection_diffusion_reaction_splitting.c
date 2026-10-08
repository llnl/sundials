/*---------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 *---------------------------------------------------------------
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
 *---------------------------------------------------------------
 * Example problem:
 *
 * The following test simulates a simple 1D advection-diffusion-
 * reaction equation,
 *    u_t = (a/2)*(u^2)_x + b*u_xx + c*(u - u^3)
 * for t in [0, 1], x in [0, 1], with initial conditions
 *    u(0,x) = u_0
 * and Dirichlet boundary conditions at x=0 and x=1
 *    u(0,t) = u(1,t) = u_0
 *
 * The spatial derivatives are computed using second-order
 * centered differences, with the data distributed over N
 * points (excluding boundary points) on a uniform spatial grid.
 *
 * This program solves the problem with an operator splitting
 * method where advection is treated with a strong stability
 * preserving ERK method, diffusion is treated with a DIRK
 * method, and reaction is treated with a different ERK method.
 *
 * Outputs are printed at equal intervals, and run statistics are
 * printed at the end.
 *---------------------------------------------------------------*/

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <arkode/arkode_splittingstep.h>
#include <math.h>
#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunmatrix/sunmatrix_band.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define FSYM "Lf"
#else
#define GSYM "g"
#define FSYM "f"
#endif

/* user data structure */
typedef struct
{
  sunindextype N; /* number of grid points (excluding boundaries) */
  sunrealtype dx; /* mesh spacing */
  sunrealtype a;  /* advection coefficient */
  sunrealtype b;  /* diffusion coefficient */
  sunrealtype c;  /* reaction coefficient */
  sunrealtype u0; /* initial and boundary values */
} UserData;

/* User-supplied Functions Called by the Solver */
static int f_advection(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int f_diffusion(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int jac_diffusion(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac,
                         void* user_data, N_Vector tmp1, N_Vector tmp2,
                         N_Vector tmp3);
static int f_reaction(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/* Private function to check function return values */
static int check_flag(void* flagvalue, const char* funcname, int opt);

int main(void)
{
  /* Problem parameters */
  sunrealtype T0 = SUN_RCONST(0.0);
  sunrealtype Tf = SUN_RCONST(1.0);
  sunrealtype DT = SUN_RCONST(0.06);
  sunindextype N = 128;
  UserData udata = {.N  = N,
                    .dx = SUN_RCONST(1.0) / (N + 1),
                    .a  = SUN_RCONST(1.0),
                    .b  = SUN_RCONST(0.125),
                    .c  = SUN_RCONST(4.0),
                    .u0 = SUN_RCONST(0.1)};

  printf("\n1D Advection-Diffusion-Reaction PDE test problem:\n");
  printf("  N = %li\n", (long int)udata.N);
  printf("  advection coefficient = %" GSYM "\n", udata.a);
  printf("  diffusion coefficient = %" GSYM "\n", udata.b);
  printf("  reaction coefficient  = %" GSYM "\n\n", udata.c);

  /* Create the SUNDIALS context object for this simulation */
  SUNContext ctx;
  int flag = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_flag(&flag, "SUNContext_Create", 1)) { return 1; }

  /* Initialize vector with initial condition */
  N_Vector y = N_VNew_Serial(udata.N, ctx);
  if (check_flag(y, "N_VNew_Serial", 0)) { return 1; }
  N_VConst(udata.u0, y);

  /* Create advection integrator */
  void* advection_mem = ERKStepCreate(f_advection, T0, y, ctx);
  if (check_flag(advection_mem, "ERKStepCreate", 0)) { return 1; }

  flag = ARKodeSetUserData(advection_mem, &udata);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }

  /* Choose a strong stability preserving method for advecton */
  flag = ERKStepSetTableNum(advection_mem, ARKODE_SHU_OSHER_3_2_3);
  if (check_flag(&flag, "ERKStepSetTableNum", 1)) { return 1; }

  SUNStepper advection_stepper;
  flag = ARKodeCreateSUNStepper(advection_mem, &advection_stepper);
  if (check_flag(&flag, "ARKodeCreateSUNStepper", 1)) { return 1; }

  /* Create diffusion integrator */
  void* diffusion_mem = ARKStepCreate(NULL, f_diffusion, T0, y, ctx);
  if (check_flag(diffusion_mem, "ARKStepCreate", 0)) { return 1; }

  flag = ARKodeSetUserData(diffusion_mem, &udata);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }

  flag = ARKodeSetOrder(diffusion_mem, 3);
  if (check_flag(&flag, "ARKStepSetOrder", 1)) { return 1; }

  SUNMatrix jac_mat = SUNBandMatrix(udata.N, 1, 1, ctx);
  if (check_flag(jac_mat, "SUNBandMatrix", 0)) { return 1; }

  SUNLinearSolver ls = SUNLinSol_Band(y, jac_mat, ctx);
  if (check_flag(ls, "SUNLinSol_Band", 0)) { return 1; }

  flag = ARKodeSetLinearSolver(diffusion_mem, ls, jac_mat);
  if (check_flag(&flag, "ARKStepSetOrder", 1)) { return 1; }

  flag = ARKodeSetJacFn(diffusion_mem, jac_diffusion);
  if (check_flag(&flag, "ARKodeSetJacFn", 1)) { return 1; }

  flag = ARKodeSetLinear(diffusion_mem, SUNFALSE);
  if (check_flag(&flag, "ARKodeSetLinear", 1)) { return 1; }

  SUNStepper diffusion_stepper;
  flag = ARKodeCreateSUNStepper(diffusion_mem, &diffusion_stepper);
  if (check_flag(&flag, "ARKodeCreateSUNStepper", 1)) { return 1; }

  /* Create reaction integrator */
  void* reaction_mem = ERKStepCreate(f_reaction, T0, y, ctx);
  if (check_flag(reaction_mem, "ERKStepCreate", 0)) { return 1; }

  flag = ARKodeSetUserData(reaction_mem, &udata);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }

  flag = ARKodeSetOrder(reaction_mem, 3);
  if (check_flag(&flag, "ARKodeSetOrder", 1)) { return 1; }

  SUNStepper reaction_stepper;
  flag = ARKodeCreateSUNStepper(reaction_mem, &reaction_stepper);
  if (check_flag(&flag, "ARKodeCreateSUNStepper", 1)) { return 1; }

  /* Create operator splitting integrator */
  SUNStepper steppers[] = {advection_stepper, diffusion_stepper,
                           reaction_stepper};
  void* arkode_mem      = SplittingStepCreate(steppers, 3, T0, y, ctx);
  if (check_flag(arkode_mem, "SplittingStepCreate", 0)) { return 1; }

  flag = ARKodeSetFixedStep(arkode_mem, DT);
  if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }

  flag = ARKodeSetStopTime(arkode_mem, Tf);
  if (check_flag(&flag, "ARKodeSetStopTime", 1)) { return 1; }

  /* Evolve solution in time */
  sunrealtype tret = T0;
  printf("        t      ||u||_rms\n");
  printf("   ----------------------\n");
  printf("  %10.6" FSYM "  %10.6f\n", tret, sqrt(N_VDotProd(y, y) / udata.N));
  while (tret < Tf)
  {
    flag = ARKodeEvolve(arkode_mem, Tf, y, &tret, ARK_ONE_STEP);
    if (check_flag(&flag, "ARKodeEvolve", 1)) { return 1; }
    printf("  %10.6" FSYM "  %10.6f\n", tret, sqrt(N_VDotProd(y, y) / udata.N));
  }
  printf("   ----------------------\n");

  /* Print statistics */
  printf("\nSplitting Stepper Statistics:\n");
  flag = ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(&flag, "ARKodePrintAllStats", 1)) { return 1; }

  printf("\nAdvection Stepper Statistics:\n");
  flag = ARKodePrintAllStats(advection_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(&flag, "ARKodePrintAllStats", 1)) { return 1; }

  printf("\nDiffusion Stepper Statistics:\n");
  flag = ARKodePrintAllStats(diffusion_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(&flag, "ARKodePrintAllStats", 1)) { return 1; }

  printf("\nReaction Stepper Statistics:\n");
  flag = ARKodePrintAllStats(reaction_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(&flag, "ARKodePrintAllStats", 1)) { return 1; }

  /* Clean up and return with successful completion */
  N_VDestroy(y);
  ARKodeFree(&advection_mem);
  SUNStepper_Destroy(&advection_stepper);
  ARKodeFree(&diffusion_mem);
  SUNStepper_Destroy(&diffusion_stepper);
  ARKodeFree(&reaction_mem);
  SUNStepper_Destroy(&reaction_stepper);
  ARKodeFree(&arkode_mem);
  SUNLinSolFree(ls);
  SUNMatDestroy(jac_mat);
  SUNContext_Free(&ctx);

  return 0;
}

/* f routine to compute the advection RHS function. */
static int f_advection(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  UserData* udata = (UserData*)user_data;
  sunrealtype* Y  = NULL;
  Y               = N_VGetArrayPointer(y); /* access data arrays */
  if (check_flag((void*)Y, "N_VGetArrayPointer", 0)) { return 1; }

  sunrealtype* Ydot = NULL;
  Ydot              = N_VGetArrayPointer(ydot);
  if (check_flag((void*)Ydot, "N_VGetArrayPointer", 0)) { return 1; }

  sunrealtype coeff  = udata->a / (SUN_RCONST(4.0) * udata->dx);
  sunrealtype u0_sqr = udata->u0 * udata->u0;

  /* Left boundary */
  Ydot[0] = coeff * (Y[1] * Y[1] - u0_sqr);
  /* Interior */
  for (sunindextype i = 1; i < udata->N - 1; i++)
  {
    Ydot[i] = coeff * (Y[i + 1] * Y[i + 1] - Y[i - 1] * Y[i - 1]);
  }
  /* Right boundary */
  Ydot[udata->N - 1] = coeff * (u0_sqr - Y[udata->N - 1] * Y[udata->N - 1]);

  return 0;
}

/* f routine to compute the diffusion RHS function. */
static int f_diffusion(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  UserData* udata = (UserData*)user_data;
  sunrealtype* Y  = NULL;
  Y               = N_VGetArrayPointer(y); /* access data arrays */
  if (check_flag((void*)Y, "N_VGetArrayPointer", 0)) { return 1; }

  sunrealtype* Ydot = NULL;
  Ydot              = N_VGetArrayPointer(ydot);
  if (check_flag((void*)Ydot, "N_VGetArrayPointer", 0)) { return 1; }

  sunrealtype coeff = udata->b / (udata->dx * udata->dx);

  /* Left boundary */
  Ydot[0] = coeff * (udata->u0 - 2 * Y[0] + Y[1]);
  /* Interior */
  for (sunindextype i = 1; i < udata->N - 1; i++)
  {
    Ydot[i] = coeff * (Y[i + 1] - 2 * Y[i] + Y[i - 1]);
  }
  /* Right boundary */
  Ydot[udata->N - 1] = coeff *
                       (Y[udata->N - 2] - 2 * Y[udata->N - 1] + udata->u0);

  return 0;
}

/* Routine to compute the diffusion Jacobian function. */
static int jac_diffusion(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac,
                         void* user_data, N_Vector tmp1, N_Vector tmp2,
                         N_Vector tmp3)
{
  UserData* udata   = (UserData*)user_data;
  sunrealtype coeff = udata->b / (udata->dx * udata->dx);

  SM_ELEMENT_B(Jac, 0, 0) = -2 * coeff;
  for (int i = 1; i < udata->N; i++)
  {
    SM_ELEMENT_B(Jac, i - 1, i) = coeff;
    SM_ELEMENT_B(Jac, i, i)     = -2 * coeff;
    SM_ELEMENT_B(Jac, i, i - 1) = coeff;
  }

  return 0;
}

/* f routine to compute the reaction RHS function. */
static int f_reaction(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  UserData* udata = (UserData*)user_data;
  sunrealtype* Y  = NULL;
  Y               = N_VGetArrayPointer(y); /* access data arrays */
  if (check_flag((void*)Y, "N_VGetArrayPointer", 0)) { return 1; }

  sunrealtype* Ydot = NULL;
  Ydot              = N_VGetArrayPointer(ydot);
  if (check_flag((void*)Ydot, "N_VGetArrayPointer", 0)) { return 1; }

  for (sunindextype i = 0; i < udata->N; i++)
  {
    Ydot[i] = udata->c * Y[i] * (SUN_RCONST(1.0) - Y[i] * Y[i]);
  }

  return 0;
}

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_flag(void* flagvalue, const char* funcname, int opt)
{
  int* errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  /* Check if flag < 0 */
  else if (opt == 1)
  {
    errflag = (int*)flagvalue;
    if (*errflag < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return 1;
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}
