/* -----------------------------------------------------------------------------
 * Programmer(s): Sylvia Amihere @ SMU
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This example solves a sequence of NX=100 identical nonlinear systems
 *
 *    4x    - sin(y) - zi     - 1  = 0
 *    -x^2  + 5y     - cos(z) - 2i = 0
 *    -exp(-x) - y   + 6z     - 3  = 0
 *
 * using the Newton-based solvers in KINSOL.  Here, we use either the dense or
 * banded direct linear solvers, leveraging KINSOL's internal difference-
 * quotient Jacobian approximations.  Alternately, we use the SPGMR iterative
 * linear solver, either unpreconditioned, or preconditioned with KINBBDPRE.
 * Control over the linear solver choice is made using the --linear_solver
 * command-line argument:
 *
 *     0 => dense
 *     1 => banded
 *     2 => SPGMR (no preconditioning)
 *     3 => SPGMR (preconditioned with KINBBDPRE)
 *
 * The system has the analytic solution:
 *
 *    x = 0.28443101049565 + 0.27031686078054i
 *    y = 0.16117132843381 + 0.42622240595676i
 *    z = 0.64771494226506 + 0.03754877135588i
 *
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <kinsol/kinsol_bbdpre.h>      /* prototypes for KINBBDPRE module  */
#include <sunlinsol/sunlinsol_band.h>  /* access to banded SUNLinearSolver */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver  */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver  */
#include <sunmatrix/sunmatrix_band.h>  /* access to banded SUNMatrix       */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix        */
#include "kinsol/kinsol.h"             /* access to KINSOL func., consts.  */
#include "nvector/nvector_serial.h"    /* access to serial N_Vector        */

/* precision specific formatting macros */
#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#else
#define GSYM "g"
#endif

/* precision specific math function macros */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define ABS(x)  (fabs((x)))
#define SQRT(x) (sqrt((x)))
#define EXP(x)  (cexp((x)))
#define SIN(x)  (csin((x)))
#define COS(x)  (ccos((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define ABS(x)  (fabsf((x)))
#define SQRT(x) (sqrtf((x)))
#define EXP(x)  (cexpf((x)))
#define SIN(x)  (csinf((x)))
#define COS(x)  (ccosf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define ABS(x)  (fabsl((x)))
#define SQRT(x) (sqrtl((x)))
#define EXP(x)  (cexpl((x)))
#define SIN(x)  (csinl((x)))
#define COS(x)  (ccosl((x)))
#endif

/* problem constants */
#define NX  100    /* number of locations */
#define NEQ 3 * NX /* number of equations */

#define ZERO  SUN_RCONST(0.0)  /* real 0.0  */
#define PTONE SUN_RCONST(0.1)  /* real 0.1  */
#define HALF  SUN_RCONST(0.5)  /* real 0.5  */
#define ONE   SUN_RCONST(1.0)  /* real 1.0  */
#define TWO   SUN_RCONST(2.0)  /* real 2.0  */
#define THREE SUN_RCONST(3.0)  /* real 3.0  */
#define FOUR  SUN_RCONST(4.0)  /* real 4.0  */
#define FIVE  SUN_RCONST(5.0)  /* real 5.0  */
#define SIX   SUN_RCONST(6.0)  /* real 6.0  */
#define TEN   SUN_RCONST(10.0) /* real 10.0 */

/* analytic solution */
#define XTRUE SUN_CCONST(0.28443101049565, 0.27031686078054)
#define YTRUE SUN_CCONST(0.16117132843381, 0.42622240595676)
#define ZTRUE SUN_CCONST(0.64771494226506, 0.03754877135588)

/* problem options */
typedef struct
{
  sunrealtype tol;     /* solve tolerance          */
  long int maxiter;    /* max number of iterations */
  int linear_solver;   /* linear solver type       */
  long int maxliniter; /* max num. lin. iterations */
}* UserOpt;

/* Nonlinear residual function */
static int ResFunction(N_Vector u, N_Vector f, void* user_data);

/* "Local" nonlinear residual function for KINBBDPRE */
static int LocResFunction(sunindextype Nlocal, N_Vector uu, N_Vector gval,
                          void* user_data);

/* Check the system solution */
static int check_ans(N_Vector u, sunrealtype tol);

/* Set default options */
static int SetDefaults(UserOpt* uopt);

/* Read command line inputs */
static int ReadInputs(int* argc, char*** argv, UserOpt uopt);

/* Print command line options */
static void InputHelp(void);

/* Check function return values */
static int check_retval(void* returnvalue, const char* funcname, int opt);

/* -----------------------------------------------------------------------------
  * Main program
  * ---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  SUNContext sunctx;
  int retval     = 0;         /* return value flag   */
  UserOpt uopt   = NULL;      /* user options struct */
  N_Vector u     = NULL;      /* solution vector     */
  N_Vector scale = NULL;      /* scaling vector      */
  FILE* infofp   = NULL;      /* KINSOL log file     */
  int j;                      /* loop counter        */
  sunscalartype* data = NULL; /* vector data array   */
  void* kmem          = NULL; /* KINSOL memory       */
  SUNMatrix A         = NULL; /* template matrix     */
  SUNLinearSolver LS  = NULL; /* lin. solv. object.  */

  /* Set default options */
  retval = SetDefaults(&uopt);
  if (check_retval(&retval, "SetDefaults", 1)) { return (1); }

  retval = ReadInputs(&argc, &argv, uopt);
  if (check_retval(&retval, "ReadInputs", 1))
  {
    free(uopt);
    return (1);
  }

  /* -------------------------
    * Print problem description
    * ------------------------- */

  printf("Solve the nonlinear system in %d bins:\n", NX);
  printf("    4x    - sin(y) - zi     - 1  = 0\n");
  printf("   -x^2   + 5y     - cos(z) - 2i = 0\n");
  printf("   -e^-x  - y      + 6z     - 3  = 0\n");
  printf("Analytic solution in each bin:\n");
  printf("    x = %f  + %fI\n", SUN_CREAL(XTRUE), SUN_CIMAG(XTRUE));
  printf("    y = %f  + %fI\n", SUN_CREAL(YTRUE), SUN_CIMAG(YTRUE));
  printf("    z = %f  + %fI\n", SUN_CREAL(ZTRUE), SUN_CIMAG(ZTRUE));
  printf("Solution method: Newton\n");
  printf("    tolerance     = %" GSYM "\n", uopt->tol);
  printf("    max iters     = %ld\n", uopt->maxiter);
  if (uopt->linear_solver == 0) { printf("    linear solver = DENSE\n"); }
  else if (uopt->linear_solver == 1) { printf("    linear solver = BANDED\n"); }
  else if (uopt->linear_solver == 2)
  {
    printf("    linear solver = SPGMR\n");
    printf("    max lin iters = %ld\n", uopt->maxliniter);
  }
  else if (uopt->linear_solver == 3)
  {
    printf("    linear solver = SPGMR+KINBBDPRE\n");
    printf("    max lin iters = %ld\n", uopt->maxliniter);
  }

  /* Create the SUNDIALS context that all SUNDIALS objects require */
  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return (1); }

  /* --------------------------------------
   * Create vectors for solution and scales
   * -------------------------------------- */

  u = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)u, "N_VNew_Serial", 0)) { return (1); }

  scale = N_VClone(u);
  if (check_retval((void*)scale, "N_VClone", 0)) { return (1); }

  /* -----------------------------------------
   * Initialize and allocate memory for KINSOL
   * ----------------------------------------- */

  kmem = KINCreate(sunctx);
  if (check_retval((void*)kmem, "KINCreate", 0)) { return (1); }

  retval = KINInit(kmem, ResFunction, u);
  if (check_retval(&retval, "KINInit", 1)) { return (1); }

  /* -----------------------------------
   * Initialize and attach linear solver
   * ----------------------------------- */

  /* Initialize matrix data structure and solver */
  if (uopt->linear_solver == 0)
  {
    A = SUNDenseMatrix(NEQ, NEQ, sunctx);
    if (check_retval((void*)A, "SUNDenseMatrix", 0)) { return 1; }
    LS = SUNLinSol_Dense(u, A, sunctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) { return 1; }
  }
  else if (uopt->linear_solver == 1)
  {
    A = SUNBandMatrix(NEQ, 2, 2, sunctx);
    if (check_retval((void*)A, "SUNBandMatrix", 0)) { return 1; }
    LS = SUNLinSol_Band(u, A, sunctx);
    if (check_retval((void*)LS, "SUNLinSol_Band", 0)) { return 1; }
  }
  else if (uopt->linear_solver == 2)
  {
    LS = SUNLinSol_SPGMR(u, SUN_PREC_NONE, uopt->maxliniter, sunctx);
    if (check_retval((void*)LS, "SUNLinSol_SPGMR", 0)) { return 1; }
  }
  else
  {
    LS = SUNLinSol_SPGMR(u, SUN_PREC_RIGHT, uopt->maxliniter, sunctx);
    if (check_retval((void*)LS, "SUNLinSol_SPGMR", 0)) { return 1; }
  }

  retval = KINSetLinearSolver(kmem, LS, A);
  if (check_retval(&retval, "KINInit", 1)) { return (1); }

  if (uopt->linear_solver == 3)
  {
    retval = KINBBDPrecInit(kmem, NEQ, 2, 2, 2, 2, ZERO, LocResFunction, NULL);
    if (check_retval(&retval, "KINBBDPrecInit", 1)) { return (1); }
  }

  /* -------------------
   * Set optional inputs
   * ------------------- */

  /* Specify stopping tolerance based on residual */
  retval = KINSetFuncNormTol(kmem, uopt->tol);
  if (check_retval(&retval, "KINSetFuncNormTol", 1)) { return (1); }

  /* Set maximum number of iterations */
  retval = KINSetNumMaxIters(kmem, uopt->maxiter);
  if (check_retval(&retval, "KINSetNumMaxIters", 1)) { return (1); }

  /* Set info log file and print level */
  infofp = fopen("kinsol.log", "w");
  if (check_retval((void*)infofp, "fopen", 0)) { return (1); }

  /* -------------
    * Initial guess
    * ------------- */

  /* Get vector data array */
  data = N_VGetArrayPointer(u);
  if (check_retval((void*)data, "N_VGetArrayPointer", 0)) { return (1); }

  for (j = 0; j < NX; j++)
  {
    data[3 * j + 0] = SUN_CCONST(0.0, 0.5); /* x */
    data[3 * j + 1] = SUN_CCONST(0.5, 0.0); /* y */
    data[3 * j + 2] = SUN_CCONST(1.0, 0.0); /* z */
  }

  /* ----------------------------
    * Call KINSol to solve problem
    * ---------------------------- */

  /* No scaling used */
  N_VConst(ONE, scale);

  /* Call main solver */
  retval = KINSol(kmem,     /* KINSol memory block */
                  u,        /* initial guess on input; solution vector */
                  KIN_NONE, /* global strategy choice */
                  scale,    /* scaling vector, for the variable cc */
                  scale);   /* scaling vector for function values fval */
  if (check_retval(&retval, "KINSol", 1)) { return (1); }

  /* ------------------------------------
    * Get solver statistics
    * ------------------------------------ */

  /* get solver stats */
  retval = KINPrintAllStats(kmem, stdout, SUN_OUTPUTFORMAT_TABLE);
  check_retval(&retval, "KINPrintAllStats", 1);

  /* ------------------------------------
    * Print solution and check error
    * ------------------------------------ */

  /* check solution */
  retval = check_ans(u, uopt->tol);

  /* -----------
    * Free memory
    * ----------- */

  fclose(infofp);
  N_VDestroy(u);
  N_VDestroy(scale);
  KINFree(&kmem);
  free(uopt);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  SUNContext_Free(&sunctx);

  return (retval);
}

/* -----------------------------------------------------------------------------
  * Nonlinear system (repeated NX times):
  *    4x    - sin(y) - zi     - 1  = 0
  *    -x^2  + 5y     - cos(z) - 2i = 0
  *    -exp(-x) - y   + 6z     - 3  = 0
  * ---------------------------------------------------------------------------*/
int ResFunction(N_Vector u, N_Vector f, void* user_data)
{
  sunscalartype* udata = NULL;
  sunscalartype* fdata = NULL;
  sunscalartype x, y, z;
  sunindextype j;

  /* Get vector data arrays */
  udata = N_VGetArrayPointer(u);
  if (check_retval((void*)udata, "N_VGetArrayPointer", 0)) { return (-1); }

  fdata = N_VGetArrayPointer(f);
  if (check_retval((void*)fdata, "N_VGetArrayPointer", 0)) { return (-1); }

  for (j = 0; j < NX; j++)
  {
    x = udata[3 * j + 0];
    y = udata[3 * j + 1];
    z = udata[3 * j + 2];

    fdata[3 * j + 0] = FOUR * x - SIN(y) - SUN_I * z - ONE;
    fdata[3 * j + 1] = -x * x + FIVE * y - COS(z) - TWO * SUN_I;
    fdata[3 * j + 2] = -EXP(-x) - y + SIX * z - THREE;
  }

  return (0);
}

/* -----------------------------------------------------------------------------
  * "Local" nonlinear system for KINBBDPRE
  * ---------------------------------------------------------------------------*/
int LocResFunction(sunindextype Nlocal, N_Vector uu, N_Vector gval,
                   void* user_data)
{
  return ResFunction(uu, gval, user_data);
}

/* -----------------------------------------------------------------------------
  * Check the solution of the nonlinear system and return PASS or FAIL
  * ---------------------------------------------------------------------------*/
static int check_ans(N_Vector u, sunrealtype tol)
{
  sunscalartype* data = NULL;
  sunrealtype exR, eyR, ezR, exI, eyI, ezI;
  int j;

  /* Get vector data array */
  data = N_VGetArrayPointer(u);
  if (check_retval((void*)data, "N_VGetArrayPointer", 0)) { return (1); }

  for (j = 0; j < NX; j++)
  {
    /* solution error */
    exR = ABS(SUN_CREAL(data[3 * j + 0]) - SUN_CREAL(XTRUE));
    eyR = ABS(SUN_CREAL(data[3 * j + 1]) - SUN_CREAL(YTRUE));
    ezR = ABS(SUN_CREAL(data[3 * j + 2]) - SUN_CREAL(ZTRUE));

    exI = ABS(SUN_CIMAG(data[3 * j + 0]) - SUN_CIMAG(XTRUE));
    eyI = ABS(SUN_CIMAG(data[3 * j + 1]) - SUN_CIMAG(YTRUE));
    ezI = ABS(SUN_CIMAG(data[3 * j + 2]) - SUN_CIMAG(ZTRUE));

    /* print the solution and error in the first cell */
    if (j == 0)
    {
      printf("Computed solution (first cell):\n");
      printf("    x = %f + %fI\n", SUN_CREAL(data[0]), SUN_CIMAG(data[0]));
      printf("    y = %f + %fI\n", SUN_CREAL(data[1]), SUN_CIMAG(data[1]));
      printf("    z = %f + %fI\n", SUN_CREAL(data[2]), SUN_CIMAG(data[2]));

      printf("Solution error (first cell):\n");
      printf("    ex = %e + %eI\n", exR, exI);
      printf("    ey = %e + %eI\n", eyR, eyI);
      printf("    ez = %e + %eI\n", ezR, ezI);
    }

    tol *= TEN;
    if (exR > tol || eyR > tol || ezR > tol || exI > tol || eyI > tol || ezI > tol)
    {
      printf("FAIL (cell %d)\n", j);
      return (1);
    }
  }

  printf("PASS\n");
  return (0);
}

/* -----------------------------------------------------------------------------
  * Set default options
  * ---------------------------------------------------------------------------*/
static int SetDefaults(UserOpt* uopt)
{
  /* Allocate options structure */
  *uopt = NULL;
  *uopt = (UserOpt)malloc(sizeof **uopt);
  if (*uopt == NULL) { return (-1); }

  /* Set default options values */
  (*uopt)->tol           = 100 * SQRT(SUN_UNIT_ROUNDOFF);
  (*uopt)->maxiter       = 30;
  (*uopt)->maxliniter    = 10;
  (*uopt)->linear_solver = 1; /* banded LS */

  return (0);
}

/* -----------------------------------------------------------------------------
  * Read command line inputs
  * ---------------------------------------------------------------------------*/
static int ReadInputs(int* argc, char*** argv, UserOpt uopt)
{
  int arg_index = 1;

  while (arg_index < (*argc))
  {
    if (strcmp((*argv)[arg_index], "--tol") == 0)
    {
      arg_index++;
      uopt->tol = atof((*argv)[arg_index++]);
    }
    else if (strcmp((*argv)[arg_index], "--maxiter") == 0)
    {
      arg_index++;
      uopt->maxiter = atoi((*argv)[arg_index++]);
    }
    else if (strcmp((*argv)[arg_index], "--linear_solver") == 0)
    {
      arg_index++;
      uopt->linear_solver = atoi((*argv)[arg_index++]);
    }
    else if (strcmp((*argv)[arg_index], "--maxliniter") == 0)
    {
      arg_index++;
      uopt->maxliniter = atoi((*argv)[arg_index++]);
    }
    else if (strcmp((*argv)[arg_index], "--help") == 0)
    {
      InputHelp();
      return (-1);
    }
    else
    {
      printf("Error: Invalid command line parameter %s\n", (*argv)[arg_index]);
      InputHelp();
      return (-1);
    }
  }

  return (0);
}

/* -----------------------------------------------------------------------------
  * Print command line options
  * ---------------------------------------------------------------------------*/
static void InputHelp(void)
{
  printf("\n");
  printf(" Command line options:\n");
  printf("   --tol           : nonlinear solver tolerance\n");
  printf("   --maxiter       : max number of nonlinear iterations\n");
  printf("   --linear_solver : linear solver type (0 dense, 1 band, 2 SPGMR, 3 "
         "SPGMR+BBDPRE)\n");
  printf("   --maxliniter    : max number of SPGMR iterations\n");
  return;
}

/* -----------------------------------------------------------------------------
  * Check function return value
  *   opt == 0 check if returned NULL pointer
  *   opt == 1 check if returned a non-zero value
  * ---------------------------------------------------------------------------*/
static int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* errflag;

  /* Check if the function returned a NULL pointer -- no memory allocated */
  if (opt == 0)
  {
    if (returnvalue == NULL)
    {
      fprintf(stderr, "\nERROR: %s() failed -- returned NULL\n\n", funcname);
      return (1);
    }
    else { return (0); }
  }

  /* Check if the function returned a non-zero value -- internal failure */
  if (opt == 1)
  {
    errflag = (int*)returnvalue;
    if (*errflag != 0)
    {
      fprintf(stderr, "\nERROR: %s() failed -- returned %d\n\n", funcname,
              *errflag);
      return (1);
    }
    else { return (0); }
  }

  /* If we make it here then opt was not 0 or 1 */
  fprintf(stderr, "\nERROR: check_retval failed -- Invalid opt value\n\n");

  return (1);
}
