/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * Based on kinLaplace_picard_bnd.c by Carol S. Woodward @ LLNL
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
 * This example solves a 2D elliptic PDE
 *
 *    d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u - 2.0
 *
 * subject to homogeneous Dirichlet boundary conditions.
 * The PDE is discretized on a uniform NX+2 by NY+2 grid with
 * central differencing, and with boundary values eliminated,
 * leaving a system of size NEQ = NX*NY.
 * The nonlinear system is solved by KINSOL using the Picard
 * iteration and the SPGMR linear solver.
 * -----------------------------------------------------------------
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <kinsol/kinsol.h>             /* access to KINSOL func., consts. */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sundials/sundials_math.h>    /* access to SUNRexp               */
#include <sundials/sundials_types.h>   /* defs. of sunrealtype, sunindextype */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver */

/* Problem Constants */

#define NX  31     /* no. of points in x direction */
#define NY  31     /* no. of points in y direction */
#define NEQ NX* NY /* problem dimension */

#define SKIP 3 /* no. of points skipped for printing */

/* function tolerance */
#if defined(SUNDIALS_SINGLE_PRECISION)
#define FTOL SUN_RCONST(1.e-4)
#else
#define FTOL SUN_RCONST(1.e-12)
#endif

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)

/* IJth is defined in order to isolate the translation from the
   mathematical 2-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage.
   IJth(vdata,i,j) references the element in the vdata array for
   u at mesh point (i,j), where 1 <= i <= NX, 1 <= j <= NY.
   The vdata array is obtained via the call vdata = N_VGetArrayPointer(v),
   where v is an N_Vector.
   The variables are ordered by the y index j, then by the x index i. */

#define IJth(vdata, i, j) (vdata[(j - 1) + (i - 1) * NY])

/* Private functions */

static int func(N_Vector u, N_Vector f, void* user_data);
static int jactimes(N_Vector v, N_Vector Jv, N_Vector u, sunbooleantype* new_u,
                    void* user_data);
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
  sunrealtype fnormtol, fnorm;
  N_Vector y, scale;
  int retval;
  void* kmem;
  SUNLinearSolver LS;

  y = scale = NULL;
  kmem      = NULL;
  LS        = NULL;

  /* -------------------------
   * Print problem description
   * ------------------------- */

  printf("\n2D elliptic PDE on unit square\n");
  printf("   d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u + 2.0\n");
  printf(" + homogeneous Dirichlet boundary conditions\n\n");
  printf("Solution method: Anderson accelerated Picard iteration with SPGMR "
         "linear solver.\n");
  printf("Problem size: %2ld x %2ld = %4ld\n", (long int)NX, (long int)NY,
         (long int)NEQ);

  /* Create the SUNDIALS context that all SUNDIALS objects require */
  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return (1); }

  /* --------------------------------------
   * Create vectors for solution and scales
   * -------------------------------------- */

  y = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)y, "N_VNew_Serial", 0)) { return (1); }

  scale = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)scale, "N_VNew_Serial", 0)) { return (1); }

  /* ----------------------------------------------------------------------------------
   * Initialize and allocate memory for KINSOL, set parameters for Anderson acceleration
   * ---------------------------------------------------------------------------------- */

  kmem = KINCreate(sunctx);
  if (check_retval((void*)kmem, "KINCreate", 0)) { return (1); }

  /* y is used as a template */

  retval = KINInit(kmem, func, y);
  if (check_retval(&retval, "KINInit", 1)) { return (1); }

  /* -------------------
   * Set optional inputs
   * ------------------- */

  /* Use acceleration with up to 3 prior residuals */

  retval = KINSetMAA(kmem, 3);
  if (check_retval(&retval, "KINSetMAA", 1)) { return (1); }

  /* Specify stopping tolerance based on residual */

  fnormtol = FTOL;
  retval   = KINSetFuncNormTol(kmem, fnormtol);
  if (check_retval(&retval, "KINSetFuncNormTol", 1)) { return (1); }

  /* ----------------------
   * Create SUNLinearSolver
   * ---------------------- */

  LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 10, sunctx);
  if (check_retval((void*)LS, "SUNLinSol_SPGMR", 0)) { return (1); }

  /* --------------------
   * Attach linear solver
   * -------------------- */

  retval = KINSetLinearSolver(kmem, LS, NULL);
  if (check_retval(&retval, "KINSetLinearSolver", 1)) { return (1); }

  /* ------------------------------------
   * Set Jacobian vector product function
   * ------------------------------------ */

  retval = KINSetJacTimesVecFn(kmem, jactimes);
  if (check_retval(&retval, "KINSetJacTimesVecFn", 1)) { return (1); }

  /* -------------
   * Initial guess
   * ------------- */

  N_VConst(ZERO, y);
  IJth(N_VGetArrayPointer(y), 2, 2) = ONE;

  /* ----------------------------
   * Call KINSol to solve problem
   * ---------------------------- */

  /* No scaling used */
  N_VConst(ONE, scale);

  /* Call main solver */
  retval = KINSol(kmem,       /* KINSol memory block */
                  y,          /* initial guess on input; solution vector */
                  KIN_PICARD, /* global strategy choice */
                  scale,      /* scaling vector, for the variable cc */
                  scale);     /* scaling vector for function values fval */
  if (check_retval(&retval, "KINSol", 1)) { return (1); }

  /* ------------------------------------
   * Print solution and solver statistics
   * ------------------------------------ */

  /* Get scaled norm of the system function */

  retval = KINGetFuncNorm(kmem, &fnorm);
  if (check_retval(&retval, "KINGetfuncNorm", 1)) { return (1); }

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("\nComputed solution (||F|| = %Lg):\n\n", fnorm);
#else
  printf("\nComputed solution (||F|| = %g):\n\n", fnorm);
#endif
  PrintOutput(y);

  PrintFinalStats(kmem);

  /* -----------
   * Free memory
   * ----------- */

  N_VDestroy(y);
  N_VDestroy(scale);
  KINFree(&kmem);
  SUNLinSolFree(LS);
  SUNContext_Free(&sunctx);

  return (0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * System function
 */

static int func(N_Vector u, N_Vector f, void* user_data)
{
  sunrealtype dx, dy, hdiff, vdiff;
  sunrealtype hdc, vdc;
  sunrealtype uij, udn, uup, ult, urt;
  sunrealtype *udata, *fdata;

  int i, j;

  dx  = ONE / (NX + 1);
  dy  = ONE / (NY + 1);
  hdc = ONE / (dx * dx);
  vdc = ONE / (dy * dy);

  udata = N_VGetArrayPointer(u);
  fdata = N_VGetArrayPointer(f);

  for (j = 1; j <= NY; j++)
  {
    for (i = 1; i <= NX; i++)
    {
      /* Extract u at x_i, y_j and four neighboring points */

      uij = IJth(udata, i, j);
      udn = (j == 1) ? ZERO : IJth(udata, i, j - 1);
      uup = (j == NY) ? ZERO : IJth(udata, i, j + 1);
      ult = (i == 1) ? ZERO : IJth(udata, i - 1, j);
      urt = (i == NX) ? ZERO : IJth(udata, i + 1, j);

      /* Evaluate diffusion components */

      hdiff = hdc * (ult - TWO * uij + urt);
      vdiff = vdc * (uup - TWO * uij + udn);

      /* Set residual at x_i, y_j */

      IJth(fdata, i, j) = hdiff + vdiff + uij - uij * uij * uij + 2.0;
    }
  }

  return (0);
}

/*
 * Jacobian vector product function
 */

static int jactimes(N_Vector v, N_Vector Jv, N_Vector u, sunbooleantype* new_u,
                    void* user_data)
{
  sunrealtype dx, dy, hdiff, vdiff;
  sunrealtype hdc, vdc;
  sunrealtype vij, vdn, vup, vlt, vrt;
  sunrealtype *vdata, *Jvdata;

  int i, j;

  dx  = ONE / (NX + 1);
  dy  = ONE / (NY + 1);
  hdc = ONE / (dx * dx);
  vdc = ONE / (dy * dy);

  vdata  = N_VGetArrayPointer(v);
  Jvdata = N_VGetArrayPointer(Jv);

  for (j = 1; j <= NY; j++)
  {
    for (i = 1; i <= NX; i++)
    {
      /* Extract v at x_i, y_j and four neighboring points */

      vij = IJth(vdata, i, j);
      vdn = (j == 1) ? ZERO : IJth(vdata, i, j - 1);
      vup = (j == NY) ? ZERO : IJth(vdata, i, j + 1);
      vlt = (i == 1) ? ZERO : IJth(vdata, i - 1, j);
      vrt = (i == NX) ? ZERO : IJth(vdata, i + 1, j);

      /* Evaluate diffusion components */

      hdiff = hdc * (vlt - TWO * vij + vrt);
      vdiff = vdc * (vup - TWO * vij + vdn);

      /* Set Jv at x_i, y_j */

      IJth(Jvdata, i, j) = hdiff + vdiff;
    }
  }

  return (0);
}

/*
 * Print solution at selected points
 */

static void PrintOutput(N_Vector u)
{
  int i, j;
  sunrealtype dx, dy, x, y;
  sunrealtype* udata;

  dx = ONE / (NX + 1);
  dy = ONE / (NY + 1);

  udata = N_VGetArrayPointer(u);

  printf("            ");
  for (i = 1; i <= NX; i += SKIP)
  {
    x = i * dx;
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%-8.5Lf ", x);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%-8.5f ", x);
#else
    printf("%-8.5f ", x);
#endif
  }
  printf("\n\n");

  for (j = 1; j <= NY; j += SKIP)
  {
    y = j * dy;
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%-8.5Lf    ", y);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%-8.5f    ", y);
#else
    printf("%-8.5f    ", y);
#endif
    for (i = 1; i <= NX; i += SKIP)
    {
#if defined(SUNDIALS_EXTENDED_PRECISION)
      printf("%-8.5Lf ", IJth(udata, i, j));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      printf("%-8.5f ", IJth(udata, i, j));
#else
      printf("%-8.5f ", IJth(udata, i, j));
#endif
    }
    printf("\n");
  }
}

/*
 * Print final statistics
 */

static void PrintFinalStats(void* kmem)
{
  long int nni, nfe, nli, npe, nps, ncfl, nfeLS, njvevals;
  int retval;

  /* Main solver statistics */

  retval = KINGetNumNonlinSolvIters(kmem, &nni);
  check_retval(&retval, "KINGetNumNonlinSolvIters", 1);
  retval = KINGetNumFuncEvals(kmem, &nfe);
  check_retval(&retval, "KINGetNumFuncEvals", 1);

  /* Linear solver statistics */

  retval = KINGetNumLinIters(kmem, &nli);
  check_retval(&retval, "KINGetNumLinIters", 1);
  retval = KINGetNumLinFuncEvals(kmem, &nfeLS);
  check_retval(&retval, "KINGetNumLinFuncEvals", 1);
  retval = KINGetNumLinConvFails(kmem, &ncfl);
  check_retval(&retval, "KINGetNumLinConvFails", 1);
  retval = KINGetNumJtimesEvals(kmem, &njvevals);
  check_retval(&retval, "KINGetNumJtimesEvals", 1);
  retval = KINGetNumPrecEvals(kmem, &npe);
  check_retval(&retval, "KINGetNumPrecEvals", 1);
  retval = KINGetNumPrecSolves(kmem, &nps);
  check_retval(&retval, "KINGetNumPrecSolves", 1);

  printf("\nFinal Statistics.. \n\n");
  printf("nni = %6ld  nli   = %6ld  ncfl = %6ld\n", nni, nli, ncfl);
  printf("nfe = %6ld  nfeLS = %6ld  njt  = %6ld\n", nfe, nfeLS, njvevals);
  printf("npe = %6ld  nps   = %6ld", npe, nps);
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
