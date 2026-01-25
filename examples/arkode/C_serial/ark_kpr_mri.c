/* ----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *                Rujeko Chinomona @ UMBC
 * ----------------------------------------------------------------
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
 * ----------------------------------------------------------------
 * Multirate nonlinear Kvaerno-Prothero-Robinson ODE test problem:
 *
 *    [u]' = [ G  e ] [(-1+u^2-r)/(2u)] + [      r'(t)/(2u)        ]
 *    [v]    [ e -1 ] [(-2+v^2-s)/(2v)]   [ s'(t)/(2*sqrt(2+s(t))) ]
 *         = [ fs(t,u,v) ]
 *           [ ff(t,u,v) ]
 *
 * where r(t) = 0.5*cos(t),  s(t) = cos(w*t),  0 < t < 5.
 *
 * This problem has analytical solution given by
 *    u(t) = sqrt(1+r(t)),  v(t) = sqrt(2+s(t)).
 *
 * We use the parameters:
 *   e = 0.5 (fast/slow coupling strength) [default]
 *   G = -1e2 (stiffness at slow time scale) [default]
 *   w = 100  (time-scale separation factor) [default]
 *   hs = 0.01 (slow step size) [default]
 *
 * The stiffness of the slow time scale is essentially determined
 * by G, for |G| > 50 it is 'stiff' and ideally suited to a
 * multirate method that is implicit at the slow time scale.
 *
 * We select the MRI method to use based on additional inputs:
 *
 *   slow_type:
 *      0 - none (full problem at fast scale)
 *      1 - ARKODE_MIS_KW3
 *      2 - ARKODE_MRI_GARK_ERK45a
 *      3 - ARKODE_MERK21
 *      4 - ARKODE_MERK32
 *      5 - ARKODE_MERK43
 *      6 - ARKODE_MERK54
 *      7 - ARKODE_MRI_GARK_IRK21a
 *      8 - ARKODE_MRI_GARK_ESDIRK34a
 *      9 - ARKODE_IMEX_MRI_GARK3b
 *     10 - ARKODE_IMEX_MRI_GARK4
 *     11 - ARKODE_IMEX_MRI_SR21
 *     12 - ARKODE_IMEX_MRI_SR32
 *     13 - ARKODE_IMEX_MRI_SR43
 *
 *   fast_type:
 *      0 - none (full problem at slow scale)
 *      1 - esdirk-3-3 (manually entered non-embedded table)
 *      2 - ARKODE_HEUN_EULER_2_1_2
 *      3 - erk-3-3 (manually entered non-embedded table)
 *      4 - erk-4-4 (manually entered non-embeded table)
 *      5 - ARKODE_DORMAND_PRINCE_7_4_5
  *
 * The program should be run with arguments in the following order:
 *   $ ark_kpr_mri slow_type fast_type h G w e deduce_rhs
 * Not all arguments are required, but these must be omitted from
 * end-to-beginning, i.e. any one of
 *   $ ark_kpr_mri slow_type fast_type h G w e deduce_rhs
 *   $ ark_kpr_mri slow_type fast_type h G w e
 *   $ ark_kpr_mri slow_type fast_type h G w
 *   $ ark_kpr_mri slow_type fast_type h G
 *   $ ark_kpr_mri slow_type fast_type h
 *   $ ark_kpr_mri slow_type fast_type
 * are acceptable.  We require:
 *   * 0 <= solve_type <= 9
 *   * 0 < h < 1/|G|
 *   * G < 0.0
 *   * w >= 1.0
 *
 * This program solves the problem with the MRI stepper. Outputs are
 * printed at equal intervals of 0.1 and run statistics are printed
 * at the end.
 * ----------------------------------------------------------------*/

/* Header files */
#include <arkode/arkode_arkstep.h> /* prototypes for ARKStep fcts., consts */
#include <arkode/arkode_mristep.h> /* prototypes for MRIStep fcts., consts */
#include <math.h>
#include <nvector/nvector_serial.h> /* serial N_Vector type, fcts., macros  */
#include <stdio.h>
#include <sundials/sundials_math.h> /* def. math fcns, 'sunrealtype'           */
#include <sunlinsol/sunlinsol_dense.h> /* dense linear solver                  */
#include <sunmatrix/sunmatrix_dense.h> /* dense matrix type, fcts., macros     */

/* Convenience macros for calling precision-specific math functions */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define CEIL(x) (ceil((x)))
#define COS(x)  (cos((x)))
#define FABS(x) (fabs((x)))
#define SIN(x)  (sin((x)))
#define SQRT(x) (sqrt((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define CEIL(x) (ceilf((x)))
#define COS(x)  (cosf((x)))
#define FABS(x) (fabsf((x)))
#define SIN(x)  (sinf((x)))
#define SQRT(x) (sqrtf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define CEIL(x) (ceill((x)))
#define COS(x)  (cosl((x)))
#define FABS(x) (fabsl((x)))
#define SIN(x)  (sinl((x)))
#define SQRT(x) (sqrtl((x)))
#endif

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM ".20Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)

/* User-supplied functions called by the solver */
static int fse(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fsi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int ff(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fn(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Js(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Jsi(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Jn(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Jf(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */
static sunrealtype r(sunrealtype t, void* user_data);
static sunrealtype s(sunrealtype t, void* user_data);
static sunrealtype rdot(sunrealtype t, void* user_data);
static sunrealtype sdot(sunrealtype t, void* user_data);
static sunrealtype utrue(sunrealtype t, void* user_data);
static sunrealtype vtrue(sunrealtype t, void* user_data);
static int Ytrue(sunrealtype t, N_Vector y, void* user_data);
static int check_retval(void* returnvalue, const char* funcname, int opt);

/* Main Program */
int main(int argc, char* argv[])
{
  SUNContext ctx;

  /* general problem parameters */
  sunrealtype T0     = SUN_RCONST(0.0);       /* initial time */
  sunrealtype Tf     = SUN_RCONST(5.0);       /* final time */
  sunrealtype dTout  = SUN_RCONST(0.1);       /* time between outputs */
  sunindextype NEQ   = 2;                     /* number of dependent vars. */
  int Nt             = (int)CEIL(Tf / dTout); /* number of output times */
  int slow_type      = 0;                     /* problem configuration type */
  int fast_type      = 0;                     /* problem configuration type */
  sunrealtype hs     = SUN_RCONST(0.01);      /* slow step size */
  sunrealtype e      = SUN_RCONST(0.5);       /* fast/slow coupling strength */
  sunrealtype G      = SUN_RCONST(-100.0);    /* stiffness at slow time scale */
  sunrealtype w      = SUN_RCONST(100.0);     /* time-scale separation factor */
  sunrealtype reltol = SUN_RCONST(0.01);
  sunrealtype abstol = SUN_RCONST(1.0e-11);

  /* general problem variables */
  int retval;                               /* reusable error-checking flag */
  N_Vector y                        = NULL; /* vector for the solution      */
  void* arkode_mem                  = NULL; /* ARKode memory structure      */
  void* inner_arkode_mem            = NULL; /* ARKode memory structure      */
  MRIStepInnerStepper inner_stepper = NULL; /* inner stepper                */
  ARKodeButcherTable B              = NULL; /* fast method Butcher table    */
  MRIStepCoupling C                 = NULL; /* slow coupling coefficients   */
  SUNMatrix Af                      = NULL; /* matrix for fast solver       */
  SUNLinearSolver LSf               = NULL; /* fast linear solver object    */
  SUNMatrix As                      = NULL; /* matrix for slow solver       */
  SUNLinearSolver LSs               = NULL; /* slow linear solver object    */
  sunbooleantype implicit_slow      = SUNFALSE;
  sunbooleantype imex_slow          = SUNFALSE;
  sunbooleantype explicit_slow      = SUNFALSE;
  sunbooleantype no_slow            = SUNFALSE;
  sunbooleantype implicit_fast      = SUNFALSE;
  sunbooleantype explicit_fast      = SUNFALSE;
  sunbooleantype no_fast            = SUNFALSE;
  sunbooleantype deduce_rhs         = SUNFALSE;
  FILE* UFID;
  sunrealtype hf, gamma, beta, t, tout, rpar[3];
  sunrealtype uerr, verr, uerrtot, verrtot, errtot;
  int iout;
  long int nsts, nstf, nfse, nfsi, nff, nnif, nncf, njef, nnis, nncs, njes;

  /*
   * Initialization
   */

  /* Retrieve the command-line options: slow_type fast_type h G w e deduce_rhs */
  if (argc < 3)
  {
    printf("ERROR: executable requires at least two arguments [slow_type "
           "fast_type]\n");
    printf("Usage:\n");
    printf("  ark_kpr_mri slow_type fast_type h G w e deduce_rhs");
    return (-1);
  }
  slow_type = atoi(argv[1]);
  fast_type = atoi(argv[2]);
  if (argc > 3) { hs = SUNStrToReal(argv[3]); }
  if (argc > 4) { G = SUNStrToReal(argv[4]); }
  if (argc > 5) { w = SUNStrToReal(argv[5]); }
  if (argc > 6) { e = SUNStrToReal(argv[6]); }
  if (argc > 7) { deduce_rhs = (sunbooleantype)atoi(argv[7]); }

  /* Check arguments for validity */
  /*   0 <= slow_type <= 13      */
  /*   0 <= fast_type <= 5       */
  /*   G < 0.0                   */
  /*   h > 0                     */
  /*   h < 1/|G| (explicit slow) */
  /*   w >= 1.0                  */
  if ((slow_type < 0) || (slow_type > 13))
  {
    printf("ERROR: slow_type be an integer in [0,13] \n");
    return (-1);
  }
  if ((fast_type < 0) || (fast_type > 5))
  {
    printf("ERROR: fast_type be an integer in [0,5] \n");
    return (-1);
  }
  if ((slow_type == 0) && (fast_type == 0))
  {
    printf("ERROR: at least one of slow_type and fast_type must be nonzero\n");
    return (-1);
  }
  if ((slow_type >= 9) && (fast_type == 0))
  {
    printf("ERROR: example not configured for ImEx slow solver with no fast "
           "solver\n");
    return (-1);
  }
  if (G >= ZERO)
  {
    printf("ERROR: G must be a negative real number\n");
    return (-1);
  }
  if (hs <= ZERO)
  {
    printf("ERROR: hs must be in positive\n");
    return (-1);
  }
  if ((hs > ONE / FABS(G)) && (!implicit_slow))
  {
    printf("ERROR: hs must be in (0, 1/|G|)\n");
    return (-1);
  }
  if (w < ONE)
  {
    printf("ERROR: w must be >= 1.0\n");
    return (-1);
  }
  rpar[0] = G;
  rpar[1] = w;
  rpar[2] = e;
  hf      = hs / w;

  /* Initial problem output (and set implicit solver tolerances as needed) */
  printf("\nMultirate nonlinear Kvaerno-Prothero-Robinson test problem:\n");
  printf("    time domain:  (%" GSYM ",%" GSYM "]\n", T0, Tf);
  printf("    hs = %" GSYM "\n", hs);
  printf("    hf = %" GSYM "\n", hf);
  printf("    G = %" GSYM "\n", G);
  printf("    w = %" GSYM "\n", w);
  printf("    e = %" GSYM "\n", e);
  switch (slow_type)
  {
  case (0):
    printf("    slow solver: none\n");
    no_slow = SUNTRUE;
    break;
  case (1):
    printf("    slow solver: ARKODE_MIS_KW3\n");
    explicit_slow = SUNTRUE;
    break;
  case (2):
    printf("    slow solver: ARKODE_MRI_GARK_ERK45a\n");
    explicit_slow = SUNTRUE;
    break;
  case (3):
    printf("    slow solver: ARKODE_MERK21\n");
    explicit_slow = SUNTRUE;
    break;
  case (4):
    printf("    slow solver: ARKODE_MERK32\n");
    explicit_slow = SUNTRUE;
    break;
  case (5):
    printf("    slow solver: ARKODE_MERK43\n");
    explicit_slow = SUNTRUE;
    break;
  case (6):
    printf("    slow solver: ARKODE_MERK54\n");
    explicit_slow = SUNTRUE;
    break;
  case (7):
    printf("    slow solver: ARKODE_MRI_GARK_IRK21a\n");
    implicit_slow = SUNTRUE;
    reltol        = SUNMAX(hs * hs, SUN_RCONST(1.0e-10));
    abstol        = SUN_RCONST(1.0e-11);
    printf("      reltol = %.2" ESYM ",  abstol = %.2" ESYM "\n", reltol, abstol);
    break;
  case (8):
    printf("    slow solver: ARKODE_MRI_GARK_ESDIRK34a\n");
    implicit_slow = SUNTRUE;
    reltol        = SUNMAX(hs * hs * hs, SUN_RCONST(1.0e-10));
    abstol        = SUN_RCONST(1.0e-11);
    printf("      reltol = %.2" ESYM ",  abstol = %.2" ESYM "\n", reltol, abstol);
    break;
  case (9):
    printf("    slow solver: ARKODE_IMEX_MRI_GARK3b\n");
    imex_slow = SUNTRUE;
    reltol    = SUNMAX(hs * hs * hs, SUN_RCONST(1.0e-10));
    abstol    = SUN_RCONST(1.0e-11);
    printf("      reltol = %.2" ESYM ",  abstol = %.2" ESYM "\n", reltol, abstol);
    break;
  case (10):
    printf("    slow solver: ARKODE_IMEX_MRI_GARK4\n");
    imex_slow = SUNTRUE;
    reltol    = SUNMAX(hs * hs * hs * hs, SUN_RCONST(1.0e-14));
    abstol    = SUN_RCONST(1.0e-14);
    printf("      reltol = %.2" ESYM ",  abstol = %.2" ESYM "\n", reltol, abstol);
    break;
  case (11):
    printf("    slow solver: ARKODE_IMEX_MRI_SR21\n");
    imex_slow = SUNTRUE;
    reltol    = SUNMAX(hs * hs, SUN_RCONST(1.0e-10));
    abstol    = SUN_RCONST(1.0e-11);
    printf("      reltol = %.2" ESYM ",  abstol = %.2" ESYM "\n", reltol, abstol);
    break;
  case (12):
    printf("    slow solver: ARKODE_IMEX_MRI_SR32\n");
    imex_slow = SUNTRUE;
    reltol    = SUNMAX(hs * hs * hs, SUN_RCONST(1.0e-10));
    abstol    = SUN_RCONST(1.0e-11);
    printf("      reltol = %.2" ESYM ",  abstol = %.2" ESYM "\n", reltol, abstol);
    break;
  case (13):
    printf("    slow solver: ARKODE_IMEX_MRI_SR43\n");
    imex_slow = SUNTRUE;
    reltol    = SUNMAX(hs * hs * hs * hs, SUN_RCONST(1.0e-14));
    abstol    = SUN_RCONST(1.0e-14);
    printf("      reltol = %.2" ESYM ",  abstol = %.2" ESYM "\n", reltol, abstol);
    break;
  }
  switch (fast_type)
  {
  case (0):
    printf("    fast solver: none\n");
    no_fast = SUNTRUE;
    break;
  case (1):
    printf("    fast solver: esdirk-3-3\n");
    implicit_fast = SUNTRUE;
    reltol        = SUNMAX(hs * hs * hs, SUN_RCONST(1.0e-10));
    abstol        = SUN_RCONST(1.0e-11);
    printf("      reltol = %.2" ESYM ",  abstol = %.2" ESYM "\n", reltol, abstol);
    break;
  case (2):
    printf("    fast solver: ARKODE_HEUN_EULER_2_1_2\n");
    explicit_fast = SUNTRUE;
    break;
  case (3):
    printf("    fast solver: erk-3-3\n");
    explicit_fast = SUNTRUE;
    break;
  case (4):
    printf("    fast solver: erk-4-4\n");
    explicit_fast = SUNTRUE;
    break;
  case (5):
    printf("    fast solver: ARKODE_DORMAND_PRINCE_7_4_5\n");
    explicit_fast = SUNTRUE;
    break;
  }

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return 1; }

  /* Create and initialize serial vector for the solution */
  y = N_VNew_Serial(NEQ, ctx);
  if (check_retval((void*)y, "N_VNew_Serial", 0)) { return 1; }
  retval = Ytrue(T0, y, rpar);
  if (check_retval(&retval, "Ytrue", 1)) { return 1; }

  /*
   * Create the fast integrator and set options
   */

  /* Initialize the fast integrator. Specify the fast right-hand side
     function in y'=fs(t,y)+ff(t,y) = fse(t,y)+fsi(t,y)+ff(t,y), the initial time T0,
     and the initial dependent variable vector y.  If the fast scale is implicit,
     set up matrix, linear solver, and Jacobian function */
  if (implicit_fast)
  {
    Af = SUNDenseMatrix(NEQ, NEQ, ctx);
    if (check_retval((void*)Af, "SUNDenseMatrix", 0)) { return 1; }
    LSf = SUNLinSol_Dense(y, Af, ctx);
    if (check_retval((void*)LSf, "SUNLinSol_Dense", 0)) { return 1; }
  }
  if (no_fast)
  {
    inner_arkode_mem = ARKStepCreate(f0, NULL, T0, y, ctx);
    if (check_retval((void*)inner_arkode_mem, "ARKStepCreate", 0)) { return 1; }
  }
  else if (explicit_fast && !no_slow)
  {
    inner_arkode_mem = ARKStepCreate(ff, NULL, T0, y, ctx);
    if (check_retval((void*)inner_arkode_mem, "ARKStepCreate", 0)) { return 1; }
  }
  else if (explicit_fast && no_slow)
  {
    inner_arkode_mem = ARKStepCreate(fn, NULL, T0, y, ctx);
    if (check_retval((void*)inner_arkode_mem, "ARKStepCreate", 0)) { return 1; }
  }
  else if (implicit_fast && no_slow)
  {
    inner_arkode_mem = ARKStepCreate(NULL, fn, T0, y, ctx);
    if (check_retval((void*)inner_arkode_mem, "ARKStepCreate", 0)) { return 1; }
    retval = ARKodeSetLinearSolver(inner_arkode_mem, LSf, Af);
    if (check_retval(&retval, "ARKodeSetLinearSolver", 1)) { return 1; }
    retval = ARKodeSetJacFn(inner_arkode_mem, Jn);
    if (check_retval(&retval, "ARKodeSetJacFn", 1)) { return 1; }
  }
  else if (implicit_fast && !no_slow)
  {
    inner_arkode_mem = ARKStepCreate(NULL, ff, T0, y, ctx);
    if (check_retval((void*)inner_arkode_mem, "ARKStepCreate", 0)) { return 1; }
    retval = ARKodeSetLinearSolver(inner_arkode_mem, LSf, Af);
    if (check_retval(&retval, "ARKodeSetLinearSolver", 1)) { return 1; }
    retval = ARKodeSetJacFn(inner_arkode_mem, Jf);
    if (check_retval(&retval, "ARKodeSetJacFn", 1)) { return 1; }
  }

  /* Set Butcher table for fast integrator */
  switch (fast_type)
  {
  case (0):
    B = ARKodeButcherTable_Alloc(3, SUNTRUE);
    if (check_retval((void*)B, "ARKodeButcherTable_Alloc", 0)) { return 1; }
    B->A[1][0] = SUN_RCONST(0.5);
    B->A[2][0] = -ONE;
    B->A[2][1] = TWO;
    B->b[0]    = ONE / SUN_RCONST(6.0);
    B->b[1]    = TWO / SUN_RCONST(3.0);
    B->b[2]    = ONE / SUN_RCONST(6.0);
    B->d[1]    = ONE;
    B->c[1]    = SUN_RCONST(0.5);
    B->c[2]    = ONE;
    B->q       = 3;
    B->p       = 2;
    retval     = ARKStepSetTables(inner_arkode_mem, 3, 2, NULL, B);
    if (check_retval(&retval, "ARKStepSetTables", 1)) { return 1; }
    break;
  case (1):
    B = ARKodeButcherTable_Alloc(3, SUNFALSE);
    if (check_retval((void*)B, "ARKodeButcherTable_Alloc", 0)) { return 1; }
    beta       = SQRT(SUN_RCONST(3.0)) / SUN_RCONST(6.0) + SUN_RCONST(0.5);
    gamma      = (-ONE / SUN_RCONST(8.0)) * (SQRT(SUN_RCONST(3.0)) + ONE);
    B->A[1][0] = SUN_RCONST(4.0) * gamma + TWO * beta;
    B->A[1][1] = ONE - SUN_RCONST(4.0) * gamma - TWO * beta;
    B->A[2][0] = SUN_RCONST(0.5) - beta - gamma;
    B->A[2][1] = gamma;
    B->A[2][2] = beta;
    B->b[0]    = ONE / SUN_RCONST(6.0);
    B->b[1]    = ONE / SUN_RCONST(6.0);
    B->b[2]    = TWO / SUN_RCONST(3.0);
    B->c[1]    = ONE;
    B->c[2]    = SUN_RCONST(0.5);
    B->q       = 3;
    retval     = ARKStepSetTables(inner_arkode_mem, 3, 0, B, NULL);
    if (check_retval(&retval, "ARKStepSetTables", 1)) { return 1; }
    break;
  case (2):
    B = ARKodeButcherTable_LoadERK(ARKODE_HEUN_EULER_2_1_2);
    if (check_retval((void*)B, "ARKodeButcherTable_LoadERK", 0)) { return 1; }
    retval = ARKStepSetTables(inner_arkode_mem, 2, 1, NULL, B);
    if (check_retval(&retval, "ARKStepSetTables", 1)) { return 1; }
    break;
  case (3):
    B = ARKodeButcherTable_Alloc(3, SUNTRUE);
    if (check_retval((void*)B, "ARKodeButcherTable_Alloc", 0)) { return 1; }
    B->A[1][0] = SUN_RCONST(0.5);
    B->A[2][0] = -ONE;
    B->A[2][1] = TWO;
    B->b[0]    = ONE / SUN_RCONST(6.0);
    B->b[1]    = TWO / SUN_RCONST(3.0);
    B->b[2]    = ONE / SUN_RCONST(6.0);
    B->d[1]    = ONE;
    B->c[1]    = SUN_RCONST(0.5);
    B->c[2]    = ONE;
    B->q       = 3;
    B->p       = 2;
    retval     = ARKStepSetTables(inner_arkode_mem, 3, 2, NULL, B);
    if (check_retval(&retval, "ARKStepSetTables", 1)) { return 1; }
    break;
  case (4):
    B = ARKodeButcherTable_Alloc(4, SUNFALSE);
    if (check_retval((void*)B, "ARKodeButcherTable_Alloc", 0)) { return 1; }
    B->A[1][0] = SUN_RCONST(0.5);
    B->A[2][1] = SUN_RCONST(0.5);
    B->A[3][2] = ONE;
    B->b[0]    = ONE / SUN_RCONST(6.0);
    B->b[1]    = ONE / SUN_RCONST(3.0);
    B->b[2]    = ONE / SUN_RCONST(3.0);
    B->b[3]    = ONE / SUN_RCONST(6.0);
    B->c[1]    = SUN_RCONST(0.5);
    B->c[2]    = SUN_RCONST(0.5);
    B->c[3]    = ONE;
    B->q       = 4;
    retval     = ARKStepSetTables(inner_arkode_mem, 4, 0, NULL, B);
    if (check_retval(&retval, "ARKStepSetTables", 1)) { return 1; }
    break;
  case (5):
    B = ARKodeButcherTable_LoadERK(ARKODE_DORMAND_PRINCE_7_4_5);
    if (check_retval((void*)B, "ARKodeButcherTable_LoadERK", 0)) { return 1; }
    retval = ARKStepSetTables(inner_arkode_mem, 5, 4, NULL, B);
    if (check_retval(&retval, "ARKStepSetTables", 1)) { return 1; }
    break;
  }
  ARKodeButcherTable_Free(B);

  /* Set the tolerances */
  retval = ARKodeSStolerances(inner_arkode_mem, reltol, abstol);
  if (check_retval(&retval, "ARKodeSStolerances", 1)) { return 1; }

  /* Set the user data pointer */
  retval = ARKodeSetUserData(inner_arkode_mem, (void*)rpar);
  if (check_retval(&retval, "ARKodeSetUserData", 1)) { return 1; }

  /* Set the fast step size */
  retval = ARKodeSetFixedStep(inner_arkode_mem, hf);
  if (check_retval(&retval, "ARKodeSetFixedStep", 1)) { return 1; }

  /* Override any current settings with command-line options -- enforce
     the prefix "inner" */
  retval = ARKodeSetOptions(inner_arkode_mem, "inner", "", argc, argv);
  if (check_retval(&retval, "ARKodeSetOptions", 1)) { return 1; }

  /* Create inner stepper */
  retval = ARKodeCreateMRIStepInnerStepper(inner_arkode_mem, &inner_stepper);
  if (check_retval(&retval, "ARKodeCreateMRIStepInnerStepper", 1)) { return 1; }

  /*
   * Create the slow integrator and set options
   */

  /* Initialize the slow integrator. Specify the slow right-hand side
     function in y'=fs(t,y)+ff(t,y) = fse(t,y)+fsi(t,y)+ff(t,y), the initial time
     T0, the initial dependent variable vector y, and the fast integrator.  If
     the slow scale contains an implicit component, set up matrix, linear solver,
     and Jacobian function. */
  if (implicit_slow || imex_slow)
  {
    As = SUNDenseMatrix(NEQ, NEQ, ctx);
    if (check_retval((void*)As, "SUNDenseMatrix", 0)) { return 1; }
    LSs = SUNLinSol_Dense(y, As, ctx);
    if (check_retval((void*)LSs, "SUNLinSol_Dense", 0)) { return 1; }
  }
  if (no_slow)
  {
    arkode_mem = MRIStepCreate(f0, NULL, T0, y, inner_stepper, ctx);
    if (check_retval((void*)arkode_mem, "MRIStepCreate", 0)) { return 1; }
  }
  else if (explicit_slow && !no_fast)
  {
    arkode_mem = MRIStepCreate(fs, NULL, T0, y, inner_stepper, ctx);
    if (check_retval((void*)arkode_mem, "MRIStepCreate", 0)) { return 1; }
  }
  else if (explicit_slow && no_fast)
  {
    arkode_mem = MRIStepCreate(fn, NULL, T0, y, inner_stepper, ctx);
    if (check_retval((void*)arkode_mem, "MRIStepCreate", 0)) { return 1; }
  }
  else if (implicit_slow && !no_fast)
  {
    arkode_mem = MRIStepCreate(NULL, fs, T0, y, inner_stepper, ctx);
    if (check_retval((void*)arkode_mem, "MRIStepCreate", 0)) { return 1; }
    retval = ARKodeSetLinearSolver(arkode_mem, LSs, As);
    if (check_retval(&retval, "ARKodeSetLinearSolver", 1)) { return 1; }
    retval = ARKodeSetJacFn(arkode_mem, Js);
    if (check_retval(&retval, "ARKodeSetJacFn", 1)) { return 1; }
  }
  else if (implicit_slow && no_fast)
  {
    arkode_mem = MRIStepCreate(NULL, fn, T0, y, inner_stepper, ctx);
    if (check_retval((void*)arkode_mem, "MRIStepCreate", 0)) { return 1; }
    retval = ARKodeSetLinearSolver(arkode_mem, LSs, As);
    if (check_retval(&retval, "ARKodeSetLinearSolver", 1)) { return 1; }
    retval = ARKodeSetJacFn(arkode_mem, Jn);
    if (check_retval(&retval, "ARKodeSetJacFn", 1)) { return 1; }
  }
  else if (imex_slow)
  {
    arkode_mem = MRIStepCreate(fse, fsi, T0, y, inner_stepper, ctx);
    if (check_retval((void*)arkode_mem, "MRIStepCreate", 0)) { return 1; }
    retval = ARKodeSetLinearSolver(arkode_mem, LSs, As);
    if (check_retval(&retval, "ARKodeSetLinearSolver", 1)) { return 1; }
    retval = ARKodeSetJacFn(arkode_mem, Jsi);
    if (check_retval(&retval, "ARKodeSetJacFn", 1)) { return 1; }
  }

  /* Set coupling table for slow integrator */
  switch (slow_type)
  {
  case (0): /* no slow dynamics (use ERK-2-2) */
    B = ARKodeButcherTable_Alloc(2, SUNFALSE);
    if (check_retval((void*)B, "ARKodeButcherTable_Alloc", 0)) { return 1; }
    B->A[1][0] = TWO / SUN_RCONST(3.0);
    B->b[0]    = SUN_RCONST(0.25);
    B->b[1]    = SUN_RCONST(0.75);
    B->c[1]    = TWO / SUN_RCONST(3.0);
    B->q       = 2;
    C          = MRIStepCoupling_MIStoMRI(B, 2, 0);
    if (check_retval((void*)C, "MRIStepCoupling_MIStoMRI", 0)) { return 1; }
    ARKodeButcherTable_Free(B);
    break;
  case (1):
    C = MRIStepCoupling_LoadTable(ARKODE_MIS_KW3);
    if (check_retval((void*)C, "MRIStepCoupling_LoadTable", 0)) { return 1; }
    break;
  case (2):
    C = MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_ERK45a);
    if (check_retval((void*)C, "MRIStepCoupling_LoadTable", 1)) { return 1; }
    break;
  case (3):
    C = MRIStepCoupling_LoadTable(ARKODE_MERK21);
    if (check_retval((void*)C, "MRIStepCoupling_LoadTable", 1)) { return 1; }
    break;
  case (4):
    C = MRIStepCoupling_LoadTable(ARKODE_MERK32);
    if (check_retval((void*)C, "MRIStepCoupling_LoadTable", 1)) { return 1; }
    break;
  case (5):
    C = MRIStepCoupling_LoadTable(ARKODE_MERK43);
    if (check_retval((void*)C, "MRIStepCoupling_LoadTable", 1)) { return 1; }
    break;
  case (6):
    C = MRIStepCoupling_LoadTable(ARKODE_MERK54);
    if (check_retval((void*)C, "MRIStepCoupling_LoadTable", 1)) { return 1; }
    break;
  case (7):
    C = MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_IRK21a);
    if (check_retval((void*)C, "MRIStepCoupling_LoadTable", 1)) { return 1; }
    break;
  case (8):
    C = MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_ESDIRK34a);
    if (check_retval((void*)C, "MRIStepCoupling_LoadTable", 1)) { return 1; }
    break;
  case (9):
    C = MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_GARK3b);
    if (check_retval((void*)C, "MRIStepCoupling_LoadTable", 0)) { return 1; }
    break;
  case (10):
    C = MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_GARK4);
    if (check_retval((void*)C, "MRIStepCoupling_LoadTable", 0)) { return 1; }
    break;
  case (11):
    C = MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_SR21);
    if (check_retval((void*)C, "MRIStepCoupling_LoadTable", 0)) { return 1; }
    break;
  case (12):
    C = MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_SR32);
    if (check_retval((void*)C, "MRIStepCoupling_LoadTable", 0)) { return 1; }
    break;
  case (13):
    C = MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_SR43);
    if (check_retval((void*)C, "MRIStepCoupling_LoadTable", 0)) { return 1; }
    break;
  }
  retval = MRIStepSetCoupling(arkode_mem, C);
  if (check_retval(&retval, "MRIStepSetCoupling", 1)) { return 1; }
  MRIStepCoupling_Free(C); /* free coupling coefficients */

  /* Set the tolerances */
  retval = ARKodeSStolerances(arkode_mem, reltol, abstol);
  if (check_retval(&retval, "ARKodeSStolerances", 1)) { return 1; }

  /* Set the user data pointer */
  retval = ARKodeSetUserData(arkode_mem, (void*)rpar);
  if (check_retval(&retval, "ARKodeSetUserData", 1)) { return 1; }

  retval = ARKodeSetDeduceImplicitRhs(arkode_mem, deduce_rhs);
  if (check_retval(&retval, "ARKodeSetDeduceImplicitRhs", 1)) { return 1; }

  /* Set the slow step size */
  retval = ARKodeSetFixedStep(arkode_mem, hs);
  if (check_retval(&retval, "ARKodeSetFixedStep", 1)) { return 1; }

  /* Override any current settings with command-line options -- enforce
     the prefix "outer" */
  retval = ARKodeSetOptions(arkode_mem, "outer", "", argc, argv);
  if (check_retval(&retval, "ARKodeSetOptions", 1)) { return 1; }

  /*
   * Integrate ODE
   */

  /* Open output stream for results, output comment line */
  UFID = fopen("ark_kpr_mri_solution.txt", "w");
  fprintf(UFID, "# t u v uerr verr\n");

  /* output initial condition to disk */
  fprintf(UFID,
          " %.16" ESYM " %.16" ESYM " %.16" ESYM " %.16" ESYM " %.16" ESYM "\n",
          T0, NV_Ith_S(y, 0), NV_Ith_S(y, 1),
          FABS(NV_Ith_S(y, 0) - utrue(T0, rpar)),
          FABS(NV_Ith_S(y, 1) - vtrue(T0, rpar)));

  /* Main time-stepping loop: calls ARKodeEvolve to perform the
     integration, then prints results. Stops when the final time
     has been reached */
  t       = T0;
  tout    = T0 + dTout;
  uerr    = ZERO;
  verr    = ZERO;
  uerrtot = ZERO;
  verrtot = ZERO;
  errtot  = ZERO;
  printf("        t           u           v       uerr      verr\n");
  printf("   ------------------------------------------------------\n");
  printf("  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM "  %.2" ESYM "  %.2" ESYM
         "\n",
         t, NV_Ith_S(y, 0), NV_Ith_S(y, 1), uerr, verr);

  for (iout = 0; iout < Nt; iout++)
  {
    /* call integrator */
    retval = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_retval(&retval, "ARKodeEvolve", 1)) { break; }

    /* access/print solution and error */
    uerr = FABS(NV_Ith_S(y, 0) - utrue(t, rpar));
    verr = FABS(NV_Ith_S(y, 1) - vtrue(t, rpar));
    printf("  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM "  %.2" ESYM
           "  %.2" ESYM "\n",
           t, NV_Ith_S(y, 0), NV_Ith_S(y, 1), uerr, verr);
    fprintf(UFID,
            " %.16" ESYM " %.16" ESYM " %.16" ESYM " %.16" ESYM " %.16" ESYM "\n",
            t, NV_Ith_S(y, 0), NV_Ith_S(y, 1), uerr, verr);
    uerrtot += uerr * uerr;
    verrtot += verr * verr;
    errtot += uerr * uerr + verr * verr;

    /* successful solve: update time */
    tout += dTout;
    tout = (tout > Tf) ? Tf : tout;
  }
  uerrtot = SQRT(uerrtot / Nt);
  verrtot = SQRT(verrtot / Nt);
  errtot  = SQRT(errtot / Nt / 2);
  printf("   ------------------------------------------------------\n");
  fclose(UFID);

  /*
   * Finalize
   */

  /* Get some slow integrator statistics */
  retval = ARKodeGetNumSteps(arkode_mem, &nsts);
  check_retval(&retval, "ARKodeGetNumSteps", 1);
  retval = ARKodeGetNumRhsEvals(arkode_mem, 0, &nfse);
  check_retval(&retval, "ARKodeGetNumRhsEvals", 1);
  retval = ARKodeGetNumRhsEvals(arkode_mem, 1, &nfsi);
  check_retval(&retval, "ARKodeGetNumRhsEvals", 1);

  /* Get some fast integrator statistics */
  retval = ARKodeGetNumSteps(inner_arkode_mem, &nstf);
  check_retval(&retval, "ARKodeGetNumSteps", 1);
  retval = ARKodeGetNumRhsEvals(inner_arkode_mem, 0, &nff);
  check_retval(&retval, "ARKodeGetNumRhsEvals", 1);

  /* Print some final statistics */
  printf("\nFinal Solver Statistics:\n");
  printf("   Steps: nsts = %li, nstf = %li\n", nsts, nstf);
  printf("   u error = %.3" ESYM ", v error = %.3" ESYM
         ", total error = %.3" ESYM "\n",
         uerrtot, verrtot, errtot);
  if (imex_slow)
  {
    printf("   Total RHS evals:  Fse = %li, Fsi = %li,  Ff = %li\n", nfse, nfsi,
           nff);
  }
  else if (implicit_slow)
  {
    printf("   Total RHS evals:  Fs = %li,  Ff = %li\n", nfsi, nff);
  }
  else { printf("   Total RHS evals:  Fs = %li,  Ff = %li\n", nfse, nff); }

  /* Get/print slow integrator decoupled implicit solver statistics */
  if (implicit_slow || imex_slow)
  {
    retval = ARKodeGetNonlinSolvStats(arkode_mem, &nnis, &nncs);
    check_retval(&retval, "ARKodeGetNonlinSolvStats", 1);
    retval = ARKodeGetNumJacEvals(arkode_mem, &njes);
    check_retval(&retval, "ARKodeGetNumJacEvals", 1);
    printf("   Slow Newton iters = %li\n", nnis);
    printf("   Slow Newton conv fails = %li\n", nncs);
    printf("   Slow Jacobian evals = %li\n", njes);
  }

  /* Get/print fast integrator implicit solver statistics */
  if (implicit_fast)
  {
    retval = ARKodeGetNonlinSolvStats(inner_arkode_mem, &nnif, &nncf);
    check_retval(&retval, "ARKodeGetNonlinSolvStats", 1);
    retval = ARKodeGetNumJacEvals(inner_arkode_mem, &njef);
    check_retval(&retval, "ARKodeGetNumJacEvals", 1);
    printf("   Fast Newton iters = %li\n", nnif);
    printf("   Fast Newton conv fails = %li\n", nncf);
    printf("   Fast Jacobian evals = %li\n", njef);
  }

  /* Clean up and return */
  N_VDestroy(y);                            /* Free y vector */
  SUNMatDestroy(Af);                        /* free fast matrix */
  SUNLinSolFree(LSf);                       /* free fast linear solver */
  SUNMatDestroy(As);                        /* free fast matrix */
  SUNLinSolFree(LSs);                       /* free fast linear solver */
  ARKodeFree(&inner_arkode_mem);            /* Free fast integrator memory */
  MRIStepInnerStepper_Free(&inner_stepper); /* Free inner stepper */
  ARKodeFree(&arkode_mem);                  /* Free slow integrator memory */
  SUNContext_Free(&ctx);                    /* Free context */

  return 0;
}

/* ------------------------------
 * Functions called by the solver
 * ------------------------------*/

/* ff routine to compute the fast portion of the ODE RHS. */
static int ff(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* rpar   = (sunrealtype*)user_data;
  const sunrealtype e = rpar[2];
  const sunrealtype u = NV_Ith_S(y, 0);
  const sunrealtype v = NV_Ith_S(y, 1);
  sunrealtype tmp1, tmp2;

  /* fill in the RHS function:
     [0  0]*[(-1+u^2-r(t))/(2*u)] + [         0          ]
     [e -1] [(-2+v^2-s(t))/(2*v)]   [sdot(t)/(2*vtrue(t))] */
  tmp1              = (-ONE + u * u - r(t, rpar)) / (TWO * u);
  tmp2              = (-TWO + v * v - s(t, rpar)) / (TWO * v);
  NV_Ith_S(ydot, 0) = ZERO;
  NV_Ith_S(ydot, 1) = e * tmp1 - tmp2 + sdot(t, rpar) / (TWO * vtrue(t, rpar));

  /* Return with success */
  return 0;
}

/* fs routine to compute the slow portion of the ODE RHS. */
static int fs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* rpar   = (sunrealtype*)user_data;
  const sunrealtype G = rpar[0];
  const sunrealtype e = rpar[2];
  const sunrealtype u = NV_Ith_S(y, 0);
  const sunrealtype v = NV_Ith_S(y, 1);
  sunrealtype tmp1, tmp2;

  /* fill in the RHS function:
     [G e]*[(-1+u^2-r(t))/(2*u))] + [rdot(t)/(2*u)]
     [0 0] [(-2+v^2-s(t))/(2*v)]    [      0      ] */
  tmp1              = (-ONE + u * u - r(t, rpar)) / (TWO * u);
  tmp2              = (-TWO + v * v - s(t, rpar)) / (TWO * v);
  NV_Ith_S(ydot, 0) = G * tmp1 + e * tmp2 + rdot(t, rpar) / (TWO * u);
  NV_Ith_S(ydot, 1) = ZERO;

  /* Return with success */
  return 0;
}

/* fse routine to compute the slow portion of the ODE RHS. */
static int fse(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* rpar   = (sunrealtype*)user_data;
  const sunrealtype u = NV_Ith_S(y, 0);

  /* fill in the slow explicit RHS function:
     [rdot(t)/(2*u)]
     [      0      ] */
  NV_Ith_S(ydot, 0) = rdot(t, rpar) / (TWO * u);
  NV_Ith_S(ydot, 1) = ZERO;

  /* Return with success */
  return 0;
}

/* fsi routine to compute the slow portion of the ODE RHS.(currently same as fse) */
static int fsi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* rpar   = (sunrealtype*)user_data;
  const sunrealtype G = rpar[0];
  const sunrealtype e = rpar[2];
  const sunrealtype u = NV_Ith_S(y, 0);
  const sunrealtype v = NV_Ith_S(y, 1);
  sunrealtype tmp1, tmp2;

  /* fill in the slow implicit RHS function:
     [G e]*[(-1+u^2-r(t))/(2*u))]
     [0 0] [(-2+v^2-s(t))/(2*v)]  */
  tmp1              = (-ONE + u * u - r(t, rpar)) / (TWO * u);
  tmp2              = (-TWO + v * v - s(t, rpar)) / (TWO * v);
  NV_Ith_S(ydot, 0) = G * tmp1 + e * tmp2;
  NV_Ith_S(ydot, 1) = ZERO;

  /* Return with success */
  return 0;
}

static int fn(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* rpar   = (sunrealtype*)user_data;
  const sunrealtype G = rpar[0];
  const sunrealtype e = rpar[2];
  const sunrealtype u = NV_Ith_S(y, 0);
  const sunrealtype v = NV_Ith_S(y, 1);
  sunrealtype tmp1, tmp2;

  /* fill in the RHS function:
     [G e]*[(-1+u^2-r(t))/(2*u))] + [rdot(t)/(2*u)]
     [e -1] [(-2+v^2-s(t))/(2*v)]   [sdot(t)/(2*vtrue(t))] */
  tmp1              = (-ONE + u * u - r(t, rpar)) / (TWO * u);
  tmp2              = (-TWO + v * v - s(t, rpar)) / (TWO * v);
  NV_Ith_S(ydot, 0) = G * tmp1 + e * tmp2 + rdot(t, rpar) / (TWO * u);
  NV_Ith_S(ydot, 1) = e * tmp1 - tmp2 + sdot(t, rpar) / (TWO * vtrue(t, rpar));

  /* Return with success */
  return 0;
}

static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VConst(ZERO, ydot);
  return (0);
}

static int Js(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* rpar   = (sunrealtype*)user_data;
  const sunrealtype G = rpar[0];
  const sunrealtype e = rpar[2];
  const sunrealtype u = NV_Ith_S(y, 0);
  const sunrealtype v = NV_Ith_S(y, 1);

  /* fill in the Jacobian:
     [G/2 + (G*(1+r(t))-rdot(t))/(2*u^2)   e/2+e*(2+s(t))/(2*v^2)]
     [                 0                             0           ] */
  SM_ELEMENT_D(J, 0, 0) = G / TWO + (G * (ONE + r(t, rpar)) - rdot(t, rpar)) /
                                      (2 * u * u);
  SM_ELEMENT_D(J, 0, 1) = e / TWO + e * (TWO + s(t, rpar)) / (TWO * v * v);
  SM_ELEMENT_D(J, 1, 0) = ZERO;
  SM_ELEMENT_D(J, 1, 1) = ZERO;

  /* Return with success */
  return 0;
}

static int Jsi(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* rpar   = (sunrealtype*)user_data;
  const sunrealtype G = rpar[0];
  const sunrealtype e = rpar[2];
  const sunrealtype u = NV_Ith_S(y, 0);
  const sunrealtype v = NV_Ith_S(y, 1);

  /* fill in the Jacobian:
     [G/2 + (G*(1+r(t)))/(2*u^2)   e/2 + e*(2+s(t))/(2*v^2)]
     [                 0                       0           ] */
  SM_ELEMENT_D(J, 0, 0) = G / TWO + (G * (ONE + r(t, rpar))) / (2 * u * u);
  SM_ELEMENT_D(J, 0, 1) = e / TWO + e * (TWO + s(t, rpar)) / (TWO * v * v);
  SM_ELEMENT_D(J, 1, 0) = ZERO;
  SM_ELEMENT_D(J, 1, 1) = ZERO;

  /* Return with success */
  return 0;
}

static int Jn(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* rpar   = (sunrealtype*)user_data;
  const sunrealtype G = rpar[0];
  const sunrealtype e = rpar[2];
  const sunrealtype u = NV_Ith_S(y, 0);
  const sunrealtype v = NV_Ith_S(y, 1);

  /* fill in the Jacobian:
     [G/2 + (G*(1+r(t))-rdot(t))/(2*u^2)     e/2 + e*(2+s(t))/(2*v^2)]
     [e/2+e*(1+r(t))/(2*u^2)                -1/2 - (2+s(t))/(2*v^2)  ] */
  SM_ELEMENT_D(J, 0, 0) = G / TWO + (G * (ONE + r(t, rpar)) - rdot(t, rpar)) /
                                      (2 * u * u);
  SM_ELEMENT_D(J, 0, 1) = e / TWO + e * (TWO + s(t, rpar)) / (TWO * v * v);
  SM_ELEMENT_D(J, 1, 0) = e / TWO + e * (ONE + r(t, rpar)) / (TWO * u * u);
  SM_ELEMENT_D(J, 1, 1) = -ONE / TWO - (TWO + s(t, rpar)) / (TWO * v * v);

  /* Return with success */
  return 0;
}

static int Jf(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* rpar   = (sunrealtype*)user_data;
  const sunrealtype e = rpar[2];
  const sunrealtype u = NV_Ith_S(y, 0);
  const sunrealtype v = NV_Ith_S(y, 1);

  /* fill in the Jacobian:
     [        0                           0        ]
     [e/2+e*(1+r(t))/(2*u^2)  -1/2-(2+s(t))/(2*v^2)] */
  SM_ELEMENT_D(J, 0, 0) = ZERO;
  SM_ELEMENT_D(J, 0, 1) = ZERO;
  SM_ELEMENT_D(J, 1, 0) = e / TWO + e * (ONE + r(t, rpar)) / (TWO * u * u);
  SM_ELEMENT_D(J, 1, 1) = -ONE / TWO - (TWO + s(t, rpar)) / (TWO * v * v);

  /* Return with success */
  return 0;
}

/* ------------------------------
 * Private helper functions
 * ------------------------------*/

static sunrealtype r(sunrealtype t, void* user_data)
{
  return (SUN_RCONST(0.5) * COS(t));
}

static sunrealtype s(sunrealtype t, void* user_data)
{
  sunrealtype* rpar = (sunrealtype*)user_data;
  return (COS(rpar[1] * t));
}

static sunrealtype rdot(sunrealtype t, void* user_data)
{
  return (-SUN_RCONST(0.5) * SIN(t));
}

static sunrealtype sdot(sunrealtype t, void* user_data)
{
  sunrealtype* rpar = (sunrealtype*)user_data;
  return (-rpar[1] * SIN(rpar[1] * t));
}

static sunrealtype utrue(sunrealtype t, void* user_data)
{
  return (SQRT(ONE + r(t, user_data)));
}

static sunrealtype vtrue(sunrealtype t, void* user_data)
{
  return (SQRT(TWO + s(t, user_data)));
}

static int Ytrue(sunrealtype t, N_Vector y, void* user_data)
{
  NV_Ith_S(y, 0) = utrue(t, user_data);
  NV_Ith_S(y, 1) = vtrue(t, user_data);
  return (0);
}

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a retval so check if
             retval < 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return 1;
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}

/*---- end of file ----*/
