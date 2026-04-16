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
 * This example solves the nonlinear system
 *
 * 4x    - sin(y) - zi     - 1  = 0
 * -x^2  + 5y     - cos(z) - 2i = 0
 * -e^-x - y      + 6z     - 3  = 0
 *
 * using the accelerated fixed pointer solver in KINSOL. The nonlinear fixed point function is
 *
 * g1(x,y,z) = (1/4)*(sin(y) + zi + 1)
 * g2(x,y,z) = (1/5)*(x^2 + cos(z) + 2i)
 * g3(x,y,z) = (1/6)*(exp(-x) + y + 3)
 *
 *  This system has the analytic solution:
 *                                        x = 0.28443101049565 + 0.27031686078054i
 *                                        y = 0.16117132843381 + 0.42622240595676i
 *                                        z = 0.64771494226506 + 0.03754877135588i
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kinsol/kinsol.h"          /* access to KINSOL func., consts. */
#include "nvector/nvector_serial.h" /* access to serial N_Vector       */

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
#define NEQ 3 /* number of equations */

#define ZERO  SUN_RCONST(0.0)  /* real 0.0  */
#define PTONE SUN_RCONST(0.1)  /* real 0.1  */
#define HALF  SUN_RCONST(0.5)  /* real 0.5  */
#define ONE   SUN_RCONST(1.0)  /* real 1.0 */
#define TWO   SUN_RCONST(2.0)  /* real 2.0 */
#define THREE SUN_RCONST(3.0)  /* real 3.0  */
#define FOUR  SUN_RCONST(4.0)  /* real 4.0 */
#define FIVE  SUN_RCONST(5.0)  /* real 5.0 */
#define SIX   SUN_RCONST(6.0)  /* real 6.0 */
#define TEN   SUN_RCONST(10.0) /* real 10.0 */

/* analytic solution */
#define XTRUE SUN_CCONST(0.28443101049565, 0.27031686078054)
#define YTRUE SUN_CCONST(0.16117132843381, 0.42622240595676)
#define ZTRUE SUN_CCONST(0.64771494226506, 0.03754877135588)

/* problem options */
typedef struct
{
  sunrealtype tol;               /* solve tolerance                  */
  long int maxiter;              /* max number of iterations         */
  long int m_aa;                 /* number of acceleration vectors   */
  long int delay_aa;             /* number of iterations to delay AA */
  int orth_aa;                   /* orthogonalization method         */
  sunrealtype damping_fp;        /* damping parameter for FP         */
  sunrealtype damping_aa;        /* damping parameter for AA         */
  sunbooleantype use_damping_fn; /* damping function                 */
  sunbooleantype use_depth_fn;   /* depth function                   */
}* UserOpt;

/* Nonlinear fixed point function */
static int FPFunction(N_Vector u, N_Vector f, void* user_data);

static int DampingFn(long int iter, N_Vector u_val, N_Vector g_val,
                     sunrealtype* qt_fn, long int depth, void* user_data,
                     sunrealtype* damping_factor);

static int DepthFn(long int iter, N_Vector u_val, N_Vector g_val,
                   N_Vector f_val, N_Vector* df, sunrealtype* R_mat,
                   long int depth, void* user_data, long int* new_depth,
                   sunbooleantype* remove_index);

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
  int retval     = 0;    /* return value flag   */
  UserOpt uopt   = NULL; /* user options struct */
  N_Vector u     = NULL; /* solution vector     */
  N_Vector scale = NULL; /* scaling vector      */
  FILE* infofp   = NULL; /* KINSOL log file     */
  long int nni, nfe;     /* solver outputs      */
  sunscalartype* data;   /* vector data array   */
  void* kmem;            /* KINSOL memory       */

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

  printf("Solve the nonlinear system:\n");
  printf("    4x    - sin(y) - zi     - 1  = 0\n");
  printf("   -x^2   + 5y     - cos(z) - 2i = 0\n");
  printf("   -e^-x  - y      + 6z     - 3  = 0\n");
  printf("Analytic solution:\n");
  printf("    x = %f  + %fI\n", creal(XTRUE), cimag(XTRUE));
  printf("    y = %f  + %fI\n", creal(YTRUE), cimag(YTRUE));
  printf("    z = %f  + %fI\n", creal(ZTRUE), cimag(ZTRUE));
  printf("Solution method: Anderson accelerated fixed point iteration.\n");
  printf("    tolerance    = %" GSYM "\n", uopt->tol);
  printf("    max iters    = %ld\n", uopt->maxiter);
  printf("    m_aa         = %ld\n", uopt->m_aa);
  printf("    delay_aa     = %ld\n", uopt->delay_aa);
  printf("    damping_aa   = %" GSYM "\n", uopt->damping_aa);
  printf("    damping_fp   = %" GSYM "\n", uopt->damping_fp);
  if (uopt->use_damping_fn) { printf("    damping_fn   = ON\n"); }
  else { printf("    damping_fn   = OFF\n"); }
  if (uopt->use_depth_fn) { printf("    depth_fn     = ON\n"); }
  else { printf("    depth_fn     = OFF\n"); }
  printf("    orth routine = %d\n", uopt->orth_aa);

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

  /* Set number of prior residuals used in Anderson acceleration */
  retval = KINSetMAA(kmem, uopt->m_aa);
  if (check_retval(&retval, "KINSetMAA", 1)) { return (1); }

  /* Set orthogonalization routine used in Anderson acceleration */
  retval = KINSetOrthAA(kmem, uopt->orth_aa);
  if (check_retval(&retval, "KINSetOrthAA", 1)) { return (1); }

  retval = KINInit(kmem, FPFunction, u);
  if (check_retval(&retval, "KINInit", 1)) { return (1); }

  /* -------------------
    * Set optional inputs
    * ------------------- */

  /* Specify stopping tolerance based on residual */
  retval = KINSetFuncNormTol(kmem, uopt->tol);
  if (check_retval(&retval, "KINSetFuncNormTol", 1)) { return (1); }

  /* Set maximum number of iterations */
  retval = KINSetNumMaxIters(kmem, uopt->maxiter);
  if (check_retval(&retval, "KINSetNumMaxItersFuncNormTol", 1)) { return (1); }

  /* Set Fixed point damping parameter */
  if (uopt->m_aa == 0) 
  { 
    retval = KINSetDamping(kmem, uopt->damping_fp);
    if (check_retval(&retval, "KINSetDamping", 1)) { return (1);}
   }

  /* Set Anderson acceleration options */
  if (uopt->m_aa > 0)
  {
    /* Set damping parameter */
    retval = KINSetDampingAA(kmem, uopt->damping_aa);
    if (check_retval(&retval, "KINSetDampingAA", 1)) { return (1); }

    /* Set acceleration delay */
    retval = KINSetDelayAA(kmem, uopt->delay_aa);
    if (check_retval(&retval, "KINSetDelayAA", 1)) { return (1); }
  }

  if (uopt->use_damping_fn)
  {
    /* Attach user defined damping function */
    retval = KINSetDampingFn(kmem, DampingFn);
    if (check_retval(&retval, "KINSetDampingFn", 1)) { return (1); }
  }

  if (uopt->use_depth_fn)
  {
    /* Attach user defined depth function */
    retval = KINSetDepthFn(kmem, DepthFn);
    if (check_retval(&retval, "KINSetDepthFn", 1)) { return (1); }
  }

  /* Set info log file and print level */
  infofp = fopen("kinsol.log", "w");
  if (check_retval((void*)infofp, "fopen", 0)) { return (1); }

  /* -------------
    * Initial guess
    * ------------- */

  /* Get vector data array */
  data = N_VGetArrayPointer(u);
  if (check_retval((void*)data, "N_VGetArrayPointer", 0)) { return (1); }

  data[0] = SUN_CCONST(0.0, 0.0);
  data[1] = SUN_CCONST(0.0, 0.0);
  data[2] = SUN_CCONST(0.0, 0.0);

  /* ----------------------------
    * Call KINSol to solve problem
    * ---------------------------- */

  /* No scaling used */
  N_VConst(ONE, scale);

  /* Call main solver */
  retval = KINSol(kmem,   /* KINSol memory block */
                  u,      /* initial guess on input; solution vector */
                  KIN_FP, /* global strategy choice */
                  scale,  /* scaling vector, for the variable cc */
                  scale); /* scaling vector for function values fval */
  if (check_retval(&retval, "KINSol", 1)) { return (1); }

  /* ------------------------------------
    * Get solver statistics
    * ------------------------------------ */

  /* get solver stats */
  retval = KINGetNumNonlinSolvIters(kmem, &nni);
  check_retval(&retval, "KINGetNumNonlinSolvIters", 1);

  retval = KINGetNumFuncEvals(kmem, &nfe);
  check_retval(&retval, "KINGetNumFuncEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("Number of nonlinear iterations: %6ld\n", nni);
  printf("Number of function evaluations: %6ld\n", nfe);

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
  SUNContext_Free(&sunctx);

  return (retval);
}

/* -----------------------------------------------------------------------------
  * Nonlinear system
  *
  * 4x    - sin(y) - zi     - 1  = 0
  * -x^2  + 5y     - cos(z) - 2i = 0
  * -e^-x - y      + 6z     - 3  = 0
  *
  * Nonlinear fixed point function
  *
  * g1(x,y,z) = (1/4)*(sin(y) + zi + 1)
  * g2(x,y,z) = (1/5)*(x^2 + cos(z) + 2i)
  * g3(x,y,z) = (1/6)*(exp(-x) + y + 3)
  *
  * ---------------------------------------------------------------------------*/
int FPFunction(N_Vector u, N_Vector g, void* user_data)
{
  sunscalartype* udata = NULL;
  sunscalartype* gdata = NULL;
  sunscalartype x, y, z;

  /* Get vector data arrays */
  udata = N_VGetArrayPointer(u);
  if (check_retval((void*)udata, "N_VGetArrayPointer", 0)) { return (-1); }

  gdata = N_VGetArrayPointer(g);
  if (check_retval((void*)gdata, "N_VGetArrayPointer", 0)) { return (-1); }

  x = udata[0];
  y = udata[1];
  z = udata[2];

  gdata[0] = (ONE / FOUR) * (SIN(y) + SUN_CCONST(0.0, 1.0) * z + ONE);
  gdata[1] = (ONE / FIVE) * (x * x + COS(z) + SUN_CCONST(0.0, 2.0));
  gdata[2] = (ONE / SIX) * (EXP(-x) + y + THREE);

  return (0);
}

static int DampingFn(long int iter, N_Vector u_val, N_Vector g_val,
                     sunrealtype* qt_fn, long int depth, void* user_data,
                     sunrealtype* damping_factor)
{
  if (depth == 0) { *damping_factor = 0.5; }
  else
  {
    /* Compute ||Q^T fn||^2 */
    sunrealtype qt_fn_norm_sqr = ZERO;
    for (long int i = 0; i < depth; i++)
    {
      qt_fn_norm_sqr += qt_fn[i] * qt_fn[i];
    }

    /* Compute ||fn||^2 = ||G(u_n) - u_n||^2 */
    sunscalartype* g_data = N_VGetArrayPointer(g_val);
    sunscalartype* u_data = N_VGetArrayPointer(u_val);
    sunscalartype fn[3];
    for (int i = 0; i < 3; i++) { fn[i] = g_data[i] - u_data[i]; }
    sunrealtype fn_norm_sqr = ZERO;
    for (int i = 0; i < 3; i++) { fn_norm_sqr += SUN_CREAL(fn[i] * SUNCONJ(fn[i])); }

    /* Compute the gain = sqrt(1 - ||Q^T fn||^2 / ||fn||^2) */
    sunrealtype gain = SUNRsqrt(ONE - qt_fn_norm_sqr / fn_norm_sqr);

    *damping_factor = 0.9 - 0.5 * gain;
  }

  return 0;
}

static int DepthFn(long int iter, N_Vector u_val, N_Vector g_val,
                   N_Vector f_val, N_Vector* df, sunrealtype* R_mat,
                   long int depth, void* user_data, long int* new_depth,
                   sunbooleantype* remove_index)
{
  if (iter < 2) { *new_depth = 1; }
  else { *new_depth = depth; };

  return 0;
}

/* -----------------------------------------------------------------------------
  * Check the solution of the nonlinear system and return PASS or FAIL
  * ---------------------------------------------------------------------------*/
static int check_ans(N_Vector u, sunrealtype tol)
{
  sunscalartype* data = NULL;
  sunrealtype exR, eyR, ezR, exI, eyI, ezI;

  /* Get vector data array */
  data = N_VGetArrayPointer(u);
  if (check_retval((void*)data, "N_VGetArrayPointer", 0)) { return (1); }

  /* print the solution */
  printf("Computed solution:\n");
  printf("    x = %f + %fI\n", SUN_CREAL(data[0]), SUN_CIMAG(data[0]));
  printf("    y = %f + %fI\n", SUN_CREAL(data[1]), SUN_CIMAG(data[1]));
  printf("    z = %f + %fI\n", SUN_CREAL(data[2]), SUN_CIMAG(data[2]));

  /* solution error */
  exR = ABS(SUN_CREAL(data[0]) - SUN_CREAL(XTRUE));
  eyR = ABS(SUN_CREAL(data[1]) - SUN_CREAL(YTRUE));
  ezR = ABS(SUN_CREAL(data[2]) - SUN_CREAL(ZTRUE));

  exI = ABS(SUN_CIMAG(data[0]) - SUN_CIMAG(XTRUE));
  eyI = ABS(SUN_CIMAG(data[1]) - SUN_CIMAG(YTRUE));
  ezI = ABS(SUN_CIMAG(data[2]) - SUN_CIMAG(ZTRUE));

  /* print the solution error */
  printf("Solution error:\n");
  printf("    ex = %e + %eI\n", exR, exI);
  printf("    ey = %e + %eI\n", eyR, eyI);
  printf("    ez = %e + %eI\n", ezR, ezI);

  tol *= TEN;
  if (exR > tol || eyR > tol || ezR > tol || exI > tol || eyI > tol || ezI > tol)
  {
    printf("FAIL\n");
    return (1);
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
  (*uopt)->tol            = 100 * SQRT(SUN_UNIT_ROUNDOFF);
  (*uopt)->maxiter        = 30;
  (*uopt)->m_aa           = 0;               /* no acceleration */
  (*uopt)->delay_aa       = 0;               /* no delay        */
  (*uopt)->orth_aa        = 0;               /* MGS             */
  (*uopt)->damping_fp     = SUN_RCONST(1.0); /* no FP dampig    */
  (*uopt)->damping_aa     = SUN_RCONST(1.0); /* no AA damping   */
  (*uopt)->use_damping_fn = SUNFALSE;        /* no damping fn   */
  (*uopt)->use_depth_fn   = SUNFALSE;        /* no depth fn     */

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
    else if (strcmp((*argv)[arg_index], "--m_aa") == 0)
    {
      arg_index++;
      uopt->m_aa = atoi((*argv)[arg_index++]);
    }
    else if (strcmp((*argv)[arg_index], "--delay_aa") == 0)
    {
      arg_index++;
      uopt->delay_aa = atoi((*argv)[arg_index++]);
    }
    else if (strcmp((*argv)[arg_index], "--damping_fp") == 0)
    {
      arg_index++;
      uopt->damping_fp = atof((*argv)[arg_index++]);
    }
    else if (strcmp((*argv)[arg_index], "--damping_aa") == 0)
    {
      arg_index++;
      uopt->damping_aa = atof((*argv)[arg_index++]);
    }
    else if (strcmp((*argv)[arg_index], "--damping_fn") == 0)
    {
      arg_index++;
      uopt->use_damping_fn = SUNTRUE;
    }
    else if (strcmp((*argv)[arg_index], "--depth_fn") == 0)
    {
      arg_index++;
      uopt->use_depth_fn = SUNTRUE;
    }
    else if (strcmp((*argv)[arg_index], "--orth_aa") == 0)
    {
      arg_index++;
      uopt->orth_aa = atoi((*argv)[arg_index++]);
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
  printf("   --tol        : nonlinear solver tolerance\n");
  printf("   --maxiter    : max number of nonlinear iterations\n");
  printf("   --m_aa       : number of Anderson acceleration vectors\n");
  printf("   --delay_aa   : Anderson acceleration delay\n");
  printf("   --damping_fp : fixed point damping parameter\n");
  printf("   --damping_aa : Anderson acceleration damping parameter\n");
  printf("   --orth_aa    : Anderson acceleration orthogonalization method\n");
  printf("   --damping_fn : user defined damping function\n");
  printf("   --depth_fn   : user defined depth function\n");

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
