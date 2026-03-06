/*------------------------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 *------------------------------------------------------------------------------
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
 *------------------------------------------------------------------------------
 * We consider the initial value problem
 *    y' + lambda*y = y^2, y(0) = 1
 * proposed in 
 * 
 * Estep, D., et al. "An a posteriori–a priori analysis of multiscale operator
 * splitting." SIAM Journal on Numerical Analysis 46.3 (2008): 1116-1146.
 *
 * The parameter lambda is positive, t is in [0, 1], and the exact solution is
 * 
 *    y(t) = lambda*y / (y(0) - (y(0) - lambda)*exp(lambda*t))
 *
 * This program solves the problem with a splitting or forcing method which can
 * be specified with the command line syntax
 * 
 * ./ark_analytic_partitioned <integrator> <coefficients>
 *    integrator: either 'splitting' or 'forcing'
 *    coefficients (splitting only): the SplittingStepCoefficients to load
 * 
 * The linear term lambda*y and nonlinear term y^2 are treated as the two
 * partitions. The former is integrated using a time step of 5e-3, while the
 * later uses a time step of 1e-3. The overall splitting or forcing integrator
 * uses a time step of 1e-2. Once solved, the program prints the error and
 * statistics.
 *----------------------------------------------------------------------------*/

/* Header files */
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <arkode/arkode_forcingstep.h>
#include <arkode/arkode_splittingstep.h>
#include <nvector/nvector_serial.h>
#include <string.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#else
#define GSYM "g"
#endif

typedef struct
{
  sunrealtype lambda;
} UserData;

/* User-supplied Functions Called by the Solver */
static int f_linear(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int f_nonlinear(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static N_Vector exact_sol(N_Vector y0, sunrealtype tf, UserData* user_data);

/* Private function to check function return values */
static int check_flag(void* flagvalue, const char* funcname, int opt);

int main(int argc, char* argv[])
{
  /* Parse arguments */
  const char* integrator_name = (argc > 1) ? argv[1] : "splitting";
  if (strcmp(integrator_name, "splitting") != 0 &&
      strcmp(integrator_name, "forcing") != 0)
  {
    fprintf(stderr, "Invalid integrator: %s\nMust be 'splitting' or 'forcing'\n",
            integrator_name);
    return 1;
  }
  char* coefficients_name = (argc > 2) ? argv[2] : NULL;

  /* Problem parameters */
  sunrealtype t0           = SUN_RCONST(0.0);  /* initial time */
  sunrealtype tf           = SUN_RCONST(1.0);  /* final time */
  sunrealtype dt           = SUN_RCONST(0.01); /* outer time step */
  sunrealtype dt_linear    = dt / 5;           /* linear integrator time step */
  sunrealtype dt_nonlinear = dt / 10; /* nonlinear integrator time step */

  UserData user_data = {.lambda = SUN_RCONST(2.0)};

  /* Create the SUNDIALS context object for this simulation */
  SUNContext ctx;
  int flag = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_flag(&flag, "SUNContext_Create", 1)) { return 1; }

  /* Initialize vector with initial condition */
  N_Vector y = N_VNew_Serial(1, ctx);
  if (check_flag(y, "N_VNew_Serial", 0)) { return 1; }
  N_VConst(SUN_RCONST(1.0), y);

  N_Vector y_exact = exact_sol(y, tf, &user_data);

  printf("\nAnalytical ODE test problem:\n");
  printf("   integrator = %s method\n", integrator_name);
  if (coefficients_name != NULL)
  {
    printf("   coefficients = %s\n", coefficients_name);
  }
  printf("   lambda     = %" GSYM "\n", user_data.lambda);

  /* Create the integrator for the linear partition */
  void* linear_mem = ERKStepCreate(f_linear, t0, y, ctx);
  if (check_flag(linear_mem, "N_VNew_Serial", 0)) { return 1; }

  flag = ARKodeSetUserData(linear_mem, &user_data);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) { return 1; }

  flag = ARKodeSetFixedStep(linear_mem, dt_linear);
  if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }

  /* Create the integrator for the nonlinear partition */
  void* nonlinear_mem = ARKStepCreate(f_nonlinear, NULL, t0, y, ctx);
  if (check_flag(nonlinear_mem, "N_VNew_Serial", 0)) { return 1; }

  flag = ARKodeSetFixedStep(nonlinear_mem, dt_nonlinear);
  if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }

  /* Create SUNSteppers out of the integrators */
  SUNStepper steppers[2];
  ARKodeCreateSUNStepper(linear_mem, &steppers[0]);
  ARKodeCreateSUNStepper(nonlinear_mem, &steppers[1]);

  /* Create the outer integrator */
  void* arkode_mem;
  if (strcmp(integrator_name, "splitting") == 0)
  {
    arkode_mem = SplittingStepCreate(steppers, 2, t0, y, ctx);
    if (check_flag(arkode_mem, "SplittingStepCreate", 0)) { return 1; }

    if (coefficients_name != NULL)
    {
      SplittingStepCoefficients coefficients =
        SplittingStepCoefficients_LoadCoefficientsByName(coefficients_name);
      if (check_flag(coefficients,
                     "SplittingStepCoefficients_LoadCoefficientsByName", 0))
      {
        return 1;
      }

      flag = SplittingStepSetCoefficients(arkode_mem, coefficients);
      if (check_flag(&flag, "SplittingStepSetCoefficients", 1)) { return 1; }

      SplittingStepCoefficients_Destroy(&coefficients);
    }
  }
  else
  {
    arkode_mem = ForcingStepCreate(steppers[0], steppers[1], t0, y, ctx);
    if (check_flag(arkode_mem, "ForcingStepCreate", 0)) { return 1; }
  }

  flag = ARKodeSetFixedStep(arkode_mem, dt);
  if (check_flag(&flag, "ARKodeSetFixedStep", 1)) { return 1; }

  /* Compute the numerical solution */
  sunrealtype tret;
  flag = ARKodeEvolve(arkode_mem, tf, y, &tret, ARK_NORMAL);
  if (check_flag(&flag, "ARKodeEvolve", 1)) { return 1; }

  /* Print the numerical error and statistics */
  N_Vector y_err = N_VClone(y);
  if (check_flag(y_err, "N_VClone", 0)) { return 1; }
  N_VLinearSum(SUN_RCONST(1.0), y, -SUN_RCONST(1.0), y_exact, y_err);
  printf("\nError: %" GSYM "\n", N_VMaxNorm(y_err));

  printf("\nSplitting Stepper Statistics:\n");
  flag = ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(&flag, "ARKodePrintAllStats", 1)) { return 1; }

  printf("\nLinear Stepper Statistics:\n");
  flag = ARKodePrintAllStats(linear_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(&flag, "ARKodePrintAllStats", 1)) { return 1; }

  printf("\nNonlinear Stepper Statistics:\n");
  flag = ARKodePrintAllStats(nonlinear_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (check_flag(&flag, "ARKodePrintAllStats", 1)) { return 1; }

  /* Free memory */
  N_VDestroy(y);
  N_VDestroy(y_exact);
  N_VDestroy(y_err);
  ARKodeFree(&linear_mem);
  SUNStepper_Destroy(&steppers[0]);
  ARKodeFree(&nonlinear_mem);
  SUNStepper_Destroy(&steppers[1]);
  ARKodeFree(&arkode_mem);
  SUNContext_Free(&ctx);

  return 0;
}

/* RHS for f^1(t, y) = -lambda * y */
static int f_linear(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype lambda = ((UserData*)user_data)->lambda;
  N_VScale(-lambda, y, ydot);
  return 0;
}

/* RHS for f^2(t, y) = y^2 */
static int f_nonlinear(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VProd(y, y, ydot);
  return 0;
}

/* Compute the exact analytic solution */
static N_Vector exact_sol(N_Vector y0, sunrealtype tf, UserData* user_data)
{
  N_Vector sol       = N_VClone(y0);
  sunrealtype y0_val = N_VGetArrayPointer(y0)[0];
  sunrealtype lambda = user_data->lambda;
  N_VGetArrayPointer(sol)[0] =
    lambda * y0_val / (y0_val - (y0_val - lambda) * SUNRexp(lambda * tf));
  return sol;
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
