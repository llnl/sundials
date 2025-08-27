/* -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Example (serial):
 *
 * This example solves a nonlinear system from.
 *
 * Source: "Handbook of Test Problems in Local and Global Optimization",
 *             C.A. Floudas, P.M. Pardalos et al.
 *             Kluwer Academic Publishers, 1999.
 * Test problem 2 from Section 14.1, Chapter 14: Equilibrium Combustion
 *
 * This problem aims to identify the concentrations of the products of a
 * hydrocarbon combustion process at equilibrium (Meintjes and Morgan, 1990).
 *    x1*x2+x1-SUN_RCONST(3.0)*x5=0
 *    SUN_RCONST(2.0)*x1*x2+x1+SUN_RCONST(3.0)*R10*SUNRpowerI(x2,2)
 *       +x2*SUNRpowerI(x3,2)+R7*x2*x3+R9*x2*x4+R8*x2-R*x5=0
 *    SUN_RCONST(2.0)*x2*x3+R7*x2*x3+SUN_RCONST(2.0)*x5*SUNRpowerI(x3,2)
 *       +R6*x3-8*x5=0
 *    R9*x2*x4+SUN_RCONST(2.0)*SUNRpowerI(x4,2)-SUN_RCONST(4.0)*R*x5=0
 *    x1*x2+x1+R10*SUNRpowerI(x2,2)+x2*SUNRpowerI(x3,2)+R7*x2*x3
 *       +R9*x2*x4+R8*x2+R5*SUNRpowerI(x3,2)+R6*x3+SUNRpowerI(x4,2)
 *       -SUN_RCONST(1.0)=0
 * such that
 *    0.0001 <= xi <= 100 i=1,2,...,5
 * coefficient:
 *    R=SUN_RCONST(10.0);
 *    R5 = SUN_RCONST(0.193);
 *    R6 = SUN_RCONST(4.10622e-4);
 *    R7 = SUN_RCONST(5.45177e-4);
 *    R8 = SUN_RCONST(4.4975e-7);
 *    R9 = SUN_RCONST(3.40735e-5);
 *    R10 = SUN_RCONST(9.615e-7);
 *
 * The treatment of the bound constraints on x1, x2, x3, x4 and x5 is done using
 * the additional variables
 *    l1 = x1 - x1_min >= 0
 *    L1 = x1 - x1_max <= 0
 *    l2 = x2 - x2_min >= 0
 *    L2 = x2 - x2_max >= 0
 *    ...
 *    l5 = x5 - x5_min >= 0
 *    L5 = x5 - x5_max >= 0
 *
 * and using the constraint feature in KINSOL to impose
 *    l1 >= 0    l2 >= 0   ...   l5 >= 0
 *    L1 <= 0    L2 <= 0   ...   L5 >= 0
 *
 * The Equilibrium Combustion test problem has one known solutions.
 * The nonlinear system is solved by KINSOL using different
 * combinations of globalization and Jacobian update strategies
 * and with different initial guesses (leading to the known solutions).
 *
 * Constraints are imposed to make all components of the solution
 * positive.
 * -----------------------------------------------------------------
 */

#include <kinsol/kinsol.h> /* access to KINSOL func., consts. */
#include <math.h>
#include <nvector/nvector_serial.h> /* access to serial N_Vector       */
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>   /* defs. of sunrealtype, sunindextype */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix       */

/* Problem Constants */

#define NVAR 5
#define NEQ  3 * NVAR

#define FTOL SUN_RCONST(1.e-5) /* function tolerance */
#define STOL SUN_RCONST(1.e-5) /* step tolerance     */
#define R    SUN_RCONST(10.0)
#define R5   SUN_RCONST(0.193)
#define R6   SUN_RCONST(4.10622e-4)
#define R7   SUN_RCONST(5.45177e-4)
#define R8   SUN_RCONST(4.4975e-7)
#define R9   SUN_RCONST(3.40735e-5)
#define R10  SUN_RCONST(9.615e-7)

#define XMIN  SUN_RCONST(1.0e-4)
#define XMAX  SUN_RCONST(1.0e+2)
#define ONE   SUN_RCONST(1.0)

typedef struct
{
  sunrealtype lb[NVAR];
  sunrealtype ub[NVAR];
}* UserData;

/* Accessor macro */
#define Ith(v, i) NV_Ith_S(v, i - 1)

/* Functions Called by the KINSOL Solver */
static int func(N_Vector u, N_Vector f, void* user_data);

/* Private Helper Functions */
static void SetInitialGuess(N_Vector u, UserData data, sunindextype n);
static int SolveIt(void* kmem, N_Vector u, N_Vector s, int glstr, int mset);
static void PrintHeader(sunrealtype fnormtol, sunrealtype scsteptol);
static void PrintOutput(N_Vector u);
static void PrintFinalStats(void* kmem);
static int check_retval(void* retvalvalue, const char* funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(void)
{
  SUNContext sunctx;
  UserData data;
  sunrealtype fnormtol, scsteptol;
  N_Vector u1, u, s, c;
  int glstr, mset, retval;
  void* kmem;
  SUNMatrix J;
  SUNLinearSolver LS;

  u1 = u = NULL;
  s = c = NULL;
  kmem  = NULL;
  J     = NULL;
  LS    = NULL;
  data  = NULL;

  /* Create the SUNDIALS context that all SUNDIALS objects require */
  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return (1); }

  /* User data */

  data        = (UserData)malloc(sizeof *data);
  data->lb[0] = XMIN;
  data->ub[0] = XMAX;
  data->lb[1] = XMIN;
  data->ub[1] = XMAX;
  data->lb[2] = XMIN;
  data->ub[2] = XMAX;
  data->lb[3] = XMIN;
  data->ub[3] = XMAX;
  data->lb[4] = XMIN;
  data->ub[4] = XMAX;

  /* Create serial vectors of length NEQ */
  u1 = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)u1, "N_VNew_Serial", 0)) { return (1); }

  u = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)u, "N_VNew_Serial", 0)) { return (1); }

  s = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)s, "N_VNew_Serial", 0)) { return (1); }

  c = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)c, "N_VNew_Serial", 0)) { return (1); }

  N_VConst(ONE, s); /* no scaling */

  Ith(c, 1) = SUN_RCONST(0.0); /* no constraint on x1 */
  Ith(c, 2) = SUN_RCONST(0.0); /* no constraint on x2 */
  Ith(c, 3) = SUN_RCONST(0.0); /* no constraint on x3 */
  Ith(c, 4) = SUN_RCONST(0.0); /* no constraint on x4 */
  Ith(c, 5) = SUN_RCONST(0.0); /* no constraint on x5 */
  Ith(c, 6) = SUN_RCONST(1.0);  /* l1 = x1 - XMIN >= 0 */
  Ith(c, 7) = SUN_RCONST(-1.0); /* L1 = x1 - XMAX <= 0 */
  Ith(c, 8) = SUN_RCONST(1.0);  /* l2 = x2 - XMIN >= 0 */
  Ith(c, 9) = SUN_RCONST(-1.0); /* L2 = x2 - XMAX <= 0 */
  Ith(c, 10) = SUN_RCONST(1.0);  /* l3 = x3 - XMIN >= 0 */
  Ith(c, 11) = SUN_RCONST(-1.0); /* L3 = x3 - XMAX <= 0 */
  Ith(c, 12) = SUN_RCONST(1.0);  /* l4 = x4 - XMIN >= 0 */
  Ith(c, 13) = SUN_RCONST(-1.0); /* L4 = x4 - XMAX <= 0 */
  Ith(c, 14) = SUN_RCONST(1.0);  /* l5 = x5 - XMIN >= 0 */
  Ith(c, 15) = SUN_RCONST(-1.0); /* L5 = x5 - XMAX <= 0 */

  fnormtol  = FTOL;
  scsteptol = STOL;

  kmem = KINCreate(sunctx);
  if (check_retval((void*)kmem, "KINCreate", 0)) { return (1); }

  retval = KINSetUserData(kmem, data);
  if (check_retval(&retval, "KINSetUserData", 1)) { return (1); }
  retval = KINSetConstraints(kmem, c);
  if (check_retval(&retval, "KINSetConstraints", 1)) { return (1); }
  retval = KINSetFuncNormTol(kmem, fnormtol);
  if (check_retval(&retval, "KINSetFuncNormTol", 1)) { return (1); }
  retval = KINSetScaledStepTol(kmem, scsteptol);
  if (check_retval(&retval, "KINSetScaledStepTol", 1)) { return (1); }

  retval = KINSetMaxNewtonStep(kmem, SUN_RCONST(100000));
  if (check_retval(&retval, "KINSetMaxNewtonStep", 1)) { return (1); }

  retval = KINInit(kmem, func, u);
  if (check_retval(&retval, "KINInit", 1)) { return (1); }

  /* Create dense SUNMatrix */
  J = SUNDenseMatrix(NEQ, NEQ, sunctx);
  if (check_retval((void*)J, "SUNDenseMatrix", 0)) { return (1); }

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(u, J, sunctx);
  if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) { return (1); }

  /* Attach the matrix and linear solver to KINSOL */
  retval = KINSetLinearSolver(kmem, LS, J);
  if (check_retval(&retval, "KINSetLinearSolver", 1)) { return (1); }

  /* Print out the problem size, solution parameters, initial guess. */
  PrintHeader(fnormtol, scsteptol);

  /* --------------------------- */
  for (sunindextype solutionIndex=1; solutionIndex<=3; solutionIndex++)
  {
    SetInitialGuess(u1, data, solutionIndex);
    printf("\n------------------------------------------\n");
    printf("\nInitial guess on lower bounds\n");
    printf("  [x1,x2] = ");
    PrintOutput(u1);

    N_VScale(ONE, u1, u);
    glstr = KIN_NONE;
    mset  = 10;
    SolveIt(kmem, u, s, glstr, mset);

    /* --------------------------- */
    N_VScale(ONE, u1, u);
    glstr = KIN_LINESEARCH;
    mset  = 10;
    SolveIt(kmem, u, s, glstr, mset);

    /* --------------------------- */
    N_VScale(ONE, u1, u);
    glstr = KIN_NONE;
    mset  = 0;
    SolveIt(kmem, u, s, glstr, mset);

    /* --------------------------- */
    N_VScale(ONE, u1, u);
    glstr = KIN_LINESEARCH;
    mset  = 0;
    SolveIt(kmem, u, s, glstr, mset);
  }

  /* Free memory */
  N_VDestroy(u1);
  N_VDestroy(u);
  N_VDestroy(s);
  N_VDestroy(c);
  KINFree(&kmem);
  SUNLinSolFree(LS);
  SUNMatDestroy(J);
  free(data);
  SUNContext_Free(&sunctx);

  return (0);
}

static int SolveIt(void* kmem, N_Vector u, N_Vector s, int glstr, int mset)
{
  int retval;

  printf("\n");

  if (mset == 1) { printf("Exact Newton"); }
  else { printf("Modified Newton"); }

  if (glstr == KIN_NONE) { printf("\n"); }
  else { printf(" with line search\n"); }

  retval = KINSetMaxSetupCalls(kmem, mset);
  if (check_retval(&retval, "KINSetMaxSetupCalls", 1)) { return (1); }

  retval = KINSol(kmem, u, glstr, s, s);
  if (check_retval(&retval, "KINSol", 1)) { return (1); }

  printf("Solution:\n  [x1,x2] = ");
  PrintOutput(u);

  PrintFinalStats(kmem);

  return (0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY KINSOL
 *--------------------------------------------------------------------
 */

/*
 * System function for predator-prey system
 */

static int func(N_Vector u, N_Vector f, void* user_data)
{
  sunrealtype *udata, *fdata;
  sunrealtype x1, l1, L1, x2, l2, L2, x3, l3, L3, x4, l4, L4, x5, l5, L5;
  sunrealtype *lb, *ub;
  UserData data;

  data = (UserData)user_data;
  lb   = data->lb;
  ub   = data->ub;

  udata = N_VGetArrayPointer(u);
  fdata = N_VGetArrayPointer(f);

  x1 = udata[0];
  x2 = udata[1];
  x3 = udata[2];
  x4 = udata[3];
  x5 = udata[4];

  l1 = udata[5];
  L1 = udata[6];
  l2 = udata[7];
  L2 = udata[8];
  l3 = udata[9];
  L3 = udata[10];
  l4 = udata[11];
  L4 = udata[12];
  l5 = udata[13];
  L5 = udata[14];

 /*    x1*x2+x1-SUN_RCONST(3.0)*x5=0
  *    SUN_RCONST(2.0)*x1*x2+x1+SUN_RCONST(3.0)*R10*SUNRpowerI(x2,2)
  *       +x2*SUNRpowerI(x3,2)+R7*x2*x3+R9*x2*x4+R8*x2-R*x5=0
  *    SUN_RCONST(2.0)*x2*SUNRpowerI(x3,2)+R7*x2*x3+SUN_RCONST(2.0)*R5*SUNRpowerI(x3,2)
  *       +R6*x3-SUN_RCONST(8.0)*x5=0
  *    R9*x2*x4+SUN_RCONST(2.0)*SUNRpowerI(x4,2)-SUN_RCONST(4.0)*R*x5=0
  *    x1*x2+x1+R10*SUNRpowerI(x2,2)+x2*SUNRpowerI(x3,2)+R7*x2*x3
  *       +R9*x2*x4+R8*x2+R5*SUNRpowerI(x3,2)+R6*x3+SUNRpowerI(x4,2)
  *       -SUN_RCONST(1.0)=0
*/
  fdata[0] = x1*x2+x1-SUN_RCONST(3.0)*x5;
  fdata[1] = SUN_RCONST(2.0)*x1*x2+x1+SUN_RCONST(3.0)*R10*SUNRpowerI(x2,2)
            +x2*SUNRpowerI(x3,2)+R7*x2*x3+R9*x2*x4+R8*x2-R*x5;
  fdata[2] = SUN_RCONST(2.0)*x2*SUNRpowerI(x3,2)+R7*x2*x3+SUN_RCONST(2.0)*R5*SUNRpowerI(x3,2)
            +R6*x3-SUN_RCONST(8.0)*x5;
  fdata[3] = R9*x2*x4+SUN_RCONST(2.0)*SUNRpowerI(x4,2)-SUN_RCONST(4.0)*R*x5;
  fdata[4] = x1*x2+x1+R10*SUNRpowerI(x2,2)+x2*SUNRpowerI(x3,2)+R7*x2*x3
            +R9*x2*x4+R8*x2+R5*SUNRpowerI(x3,2)+R6*x3+SUNRpowerI(x4,2)
            -SUN_RCONST(1.0);
  fdata[5] = l1 - x1 + lb[0];
  fdata[6] = L1 - x1 + ub[0];
  fdata[7] = l2 - x2 + lb[1];
  fdata[8] = L2 - x2 + ub[1];
  fdata[9] = l3 - x3 + lb[2];
  fdata[10] = L3 - x3 + ub[2];
  fdata[11] = l4 - x4 + lb[3];
  fdata[12] = L4 - x4 + ub[3];
  fdata[13] = l5 - x5 + lb[4];
  fdata[14] = L5 - x5 + ub[4];

  return (0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Initial guesses
 */

static void SetInitialGuess(N_Vector u, UserData data, sunindextype n)
{
  sunrealtype x1, x2, x3, x4, x5;
  sunrealtype* udata;
  sunrealtype *lb, *ub;

  udata = N_VGetArrayPointer(u);

  lb = data->lb;
  ub = data->ub;

  /* There are nine known solutions for this problem */
  switch (n) {
    case 1:
      /* this init. guess should take us to (-5; -5) */
      x1 = ONE;
      x2 = ONE;
      x3 = ONE;
      x4 = ONE;
      x5 = ONE;
      break;
    case 2:
      /* this init. guess should take us to (-5; -3) */
      x1 = XMAX;
      x2 = XMAX;
      x3 = XMAX;
      x4 = XMAX;
      x5 = XMAX;
      break;
    default:
      /* this init. guess should take us to (-5; 5) */
      x1 = (XMIN+XMAX)/SUN_RCONST(2.0);
      x2 = (XMIN+XMAX)/SUN_RCONST(2.0);
      x3 = (XMIN+XMAX)/SUN_RCONST(2.0);
      x4 = (XMIN+XMAX)/SUN_RCONST(2.0);
      x5 = (XMIN+XMAX)/SUN_RCONST(2.0);
      break;
  }


  udata[0] = x1;
  udata[1] = x2;
  udata[2] = x3;
  udata[3] = x4;
  udata[4] = x5;
  udata[5] = x1 - lb[0];
  udata[6] = x1 - ub[0];
  udata[7] = x2 - lb[1];
  udata[8] = x2 - ub[1];
  udata[9] = x3 - lb[2];
  udata[10] = x3 - ub[2];
  udata[11] = x4 - lb[3];
  udata[12] = x4 - ub[3];
  udata[13] = x5 - lb[4];
  udata[14] = x5 - ub[4];
}

/*
 * Print first lines of output (problem description)
 */

static void PrintHeader(sunrealtype fnormtol, sunrealtype scsteptol)
{
  printf("\nHimmelblau function problem\n");
  printf("Tolerance parameters:\n");
#if defined(SUNDIALS_FLOAT128_PRECISION)
  printf("  fnormtol  = %10.6Qg\n  scsteptol = %10.6Qg\n", fnormtol, scsteptol);
#elif defined(SUNDIALS_EXTENDED_PRECISION)
  printf("  fnormtol  = %10.6Lg\n  scsteptol = %10.6Lg\n", fnormtol, scsteptol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("  fnormtol  = %10.6g\n  scsteptol = %10.6g\n", fnormtol, scsteptol);
#else
  printf("  fnormtol  = %10.6g\n  scsteptol = %10.6g\n", fnormtol, scsteptol);
#endif
}

/*
 * Print solution
 */

static void PrintOutput(N_Vector u)
{
#if defined(SUNDIALS_FLOAT128_PRECISION)
  printf(" %8.6Qg  %8.6Qg  %8.6Qg  %8.6Qg  %8.6Qg\n", Ith(u, 1), Ith(u, 2), Ith(u, 3), Ith(u, 4), Ith(u, 5));
#elif defined(SUNDIALS_EXTENDED_PRECISION)
  printf(" %8.6Lg  %8.6Lg  %8.6Lg  %8.6Lg  %8.6Lg\n", Ith(u, 1), Ith(u, 2), Ith(u, 3), Ith(u, 4), Ith(u, 5));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf(" %8.6g  %8.6g  %8.6g  %8.6g  %8.6g\n", Ith(u, 1), Ith(u, 2), Ith(u, 3), Ith(u, 4), Ith(u, 5));
#else
  printf(" %8.6g  %8.6g  %8.6g  %8.6g  %8.6g\n", Ith(u, 1), Ith(u, 2), Ith(u, 3), Ith(u, 4), Ith(u, 5));
#endif
}

/*
 * Print final statistics contained in iopt
 */

static void PrintFinalStats(void* kmem)
{
  long int nni, nfe, nje, nfeD;
  int retval;

  retval = KINGetNumNonlinSolvIters(kmem, &nni);
  check_retval(&retval, "KINGetNumNonlinSolvIters", 1);
  retval = KINGetNumFuncEvals(kmem, &nfe);
  check_retval(&retval, "KINGetNumFuncEvals", 1);

  retval = KINGetNumJacEvals(kmem, &nje);
  check_retval(&retval, "KINGetNumJacEvals", 1);
  retval = KINGetNumLinFuncEvals(kmem, &nfeD);
  check_retval(&retval, "KINGetNumLinFuncEvals", 1);

  printf("Final Statistics:\n");
  printf("  nni = %5ld    nfe  = %5ld \n", nni, nfe);
  printf("  nje = %5ld    nfeD = %5ld \n", nje, nfeD);
}

/*
 * Check function return value...
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns a retval so check if
 *             retval >= 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer
 */

static int check_retval(void* retvalvalue, const char* funcname, int opt)
{
  int* errretval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && retvalvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    errretval = (int*)retvalvalue;
    if (*errretval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *errretval);
      return (1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && retvalvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  return (0);
}
