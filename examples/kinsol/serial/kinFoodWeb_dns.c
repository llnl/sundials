/* -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
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
 * This example solves a nonlinear system that arises from a system
 * of partial differential equations. The PDE system is a food web
 * population model, with predator-prey interaction and diffusion
 * on the unit square in two dimensions. The dependent variable
 * vector is the following:
 *
 *       1   2         ns
 * c = (c , c ,  ..., c  )     (denoted by the variable cc)
 *
 * and the PDE's are as follows:
 *
 *                    i       i
 *         0 = d(i)*(c     + c    )  +  f  (x,y,c)   (i=1,...,ns)
 *                    xx      yy         i
 *
 *   where
 *
 *                   i             ns         j
 *   f  (x,y,c)  =  c  * (b(i)  + sum a(i,j)*c )
 *    i                           j=1
 *
 * The number of species is ns = 2 * np, with the first np being
 * prey and the last np being predators. The number np is both the
 * number of prey and predator species. The coefficients a(i,j),
 * b(i), d(i) are:
 *
 *   a(i,i) = -AA   (all i)
 *   a(i,j) = -GG   (i <= np , j >  np)
 *   a(i,j) =  EE   (i >  np,  j <= np)
 *   b(i) = BB * (1 + alpha * x * y)   (i <= np)
 *   b(i) =-BB * (1 + alpha * x * y)   (i >  np)
 *   d(i) = DPREY   (i <= np)
 *   d(i) = DPRED   ( i > np)
 *
 * The various scalar parameters are set using define's or in
 * routine InitUserData.
 *
 * The boundary conditions are: normal derivative = 0, and the
 * initial guess is constant in x and y, but the final solution
 * is not.
 *
 * The PDEs are discretized by central differencing on an MX by
 * MY mesh.
 *
 * The nonlinear system is solved by KINSOL using the method
 * specified in the local variable globalstrat.
 *
 * The preconditioner matrix is a block-diagonal matrix based on
 * the partial derivatives of the interaction terms f only.
 *
 * Constraints are imposed to make all components of the solution
 * positive.
 * -----------------------------------------------------------------
 * References:
 *
 * 1. Peter N. Brown and Youcef Saad,
 *    Hybrid Krylov Methods for Nonlinear Systems of Equations
 *    LLNL report UCRL-97645, November 1987.
 *
 * 2. Peter N. Brown and Alan C. Hindmarsh,
 *    Reduced Storage Matrix Methods in Stiff ODE systems,
 *    Lawrence Livermore National Laboratory Report  UCRL-95088,
 *    Rev. 1, June 1987, and  Journal of Applied Mathematics and
 *    Computation, Vol. 31 (May 1989), pp. 40-91. (Presents a
 *    description of the time-dependent version of this test
 *    problem.)
 * -----------------------------------------------------------------
 */

#include <kinsol/kinsol.h> /* access to KINSOL func., consts.      */
#include <math.h>
#include <nvector/nvector_serial.h> /* access to serial N_Vector            */
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_dense.h> /* use generic dense solver in precond. */
#include <sundials/sundials_math.h>  /* access to SUNMAX, SUNRabs, SUNRsqrt  */
#include <sundials/sundials_types.h> /* defs. of sunrealtype, sunindextype      */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver */


/* Problem Constants */

#define NUM_SPECIES \
  6 /* must equal 2*(number of prey or predators)
                              number of prey = number of predators       */

#define MX       8 /* MX = number of x mesh points */
#define MY       8 /* MY = number of y mesh points */
#define NSMX     (NUM_SPECIES * MX)
#define NEQ      (NSMX * MY)        /* number of equations in the system */
#define AA       SUN_RCONST(1.0)    /* value of coefficient AA in above eqns */
#define EE       SUN_RCONST(10000.) /* value of coefficient EE in above eqns */
#define GG       SUN_RCONST(0.5e-6) /* value of coefficient GG in above eqns */
#define BB       SUN_RCONST(1.0)    /* value of coefficient BB in above eqns */
#define DPREY    SUN_RCONST(1.0)    /* value of coefficient dprey above */
#define DPRED    SUN_RCONST(0.5)    /* value of coefficient dpred above */
#define ALPHA    SUN_RCONST(1.0)    /* value of coefficient alpha above */
#define AX       SUN_RCONST(1.0)    /* total range of x variable */
#define AY       SUN_RCONST(1.0)    /* total range of y variable */
#define FTOL     SUN_RCONST(1.e-15)  /* ftol tolerance */
#define STOL     SUN_RCONST(1.e-15) /* stol tolerance */
#define ZERO     SUN_RCONST(0.0)    /* 0. */
#define ONE      SUN_RCONST(1.0)    /* 1. */
#define TWO      SUN_RCONST(2.0)    /* 2. */
#define PREYIN   SUN_RCONST(1.2)    /* initial guess for prey concentrations. */
#define PREDIN   SUN_RCONST(33000.0) /* initial guess for predator concs.      */

/* User-defined vector access macro: IJ_Vptr */

/* IJ_Vptr is defined in order to translate from the underlying 3D structure
   of the dependent variable vector to the 1D storage scheme for an N-vector.
   IJ_Vptr(vv,i,j) returns a pointer to the location in vv corresponding to
   indices is = 0, jx = i, jy = j.    */

#define IJ_Vptr(vv, i, j) (&NV_Ith_S(vv, i * NUM_SPECIES + j * NSMX))

/* Type : UserData
   contains preconditioner blocks, pivot arrays, and problem constants */

typedef struct
{
  sunrealtype** P[MX][MY];
  sunindextype* pivot[MX][MY];
  sunrealtype **acoef, *bcoef;
  N_Vector rates;
  sunrealtype *cox, *coy;
  sunrealtype ax, ay, dx, dy;
  sunrealtype uround, sqruround;
  sunindextype mx, my, ns, np;
}* UserData;

/* Functions Called by the KINSOL Solver */

static int func(N_Vector cc, N_Vector fval, void* user_data);

/* Private Helper Functions */

static UserData AllocUserData(void);
static void InitUserData(UserData data);
static void FreeUserData(UserData data);
static void SetInitialProfiles(N_Vector cc, N_Vector sc);
static void PrintHeader(int globalstrategy,
                        sunrealtype fnormtol, sunrealtype scsteptol);
static void PrintOutput(N_Vector cc);
static void PrintFinalStats(void* kmem);
static void WebRate(sunrealtype xx, sunrealtype yy, sunrealtype* cxy,
                    sunrealtype* ratesxy, void* user_data);
static sunrealtype DotProd(sunindextype size, sunrealtype* x1, sunrealtype* x2);
static int check_retval(void* retvalvalue, const char* funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(void)
{
  SUNContext sunctx;
  int globalstrategy;
  sunrealtype fnormtol, scsteptol;
  N_Vector cc, sc, fval;
  UserData data;
  int retval, maxl, maxlrst;
  void* kmem;
  SUNMatrix J;
  SUNLinearSolver LS;

  cc = sc               = NULL;
  kmem                  = NULL;
  LS                    = NULL;
  data                  = NULL;

  /* Create the SUNDIALS context that all SUNDIALS objects require */
  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return (1); }

  /* Allocate memory, and set problem data, initial values, tolerances */
  globalstrategy = KIN_NONE;

  data = AllocUserData();
  if (check_retval((void*)data, "AllocUserData", 2)) { return (1); }
  InitUserData(data);

  /* Create serial vectors of length NEQ */
  /* Create serial vectors of length NEQ */
  fval = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)fval, "N_VNew_Serial", 0)) { return (1); }
  cc = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)cc, "N_VNew_Serial", 0)) { return (1); }
  sc = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)sc, "N_VNew_Serial", 0)) { return (1); }
  data->rates = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)data->rates, "N_VNew_Serial", 0)) { return (1); }



  SetInitialProfiles(cc, sc);
  //N_VConst(ONE, sc); /* no scaling */
  fnormtol  = FTOL;
  scsteptol = STOL;

  /* Call KINCreate/KINInit to initialize KINSOL.
     A pointer to KINSOL problem memory is returned and stored in kmem. */
  kmem = KINCreate(sunctx);
  if (check_retval((void*)kmem, "KINCreate", 0)) { return (1); }

  /* Vector cc passed as template vector. */
  retval = KINInit(kmem, func, cc);
  if (check_retval(&retval, "KINInit", 1)) { return (1); }

  retval = KINSetUserData(kmem, data);
  if (check_retval(&retval, "KINSetUserData", 1)) { return (1); }

  retval = KINSetFuncNormTol(kmem, fnormtol);
  if (check_retval(&retval, "KINSetFuncNormTol", 1)) { return (1); }
  retval = KINSetScaledStepTol(kmem, scsteptol);
  if (check_retval(&retval, "KINSetScaledStepTol", 1)) { return (1); }



  /* Create dense SUNMatrix */
  J = SUNDenseMatrix(NEQ, NEQ, sunctx);
  if (check_retval((void*)J, "SUNDenseMatrix", 0)) { return (1); }

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(cc, J, sunctx);
  if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) { return (1); }

  /* Attach the matrix and linear solver to KINSOL */
  retval = KINSetLinearSolver(kmem, LS, J);
  if (check_retval(&retval, "KINSetLinearSolver", 1)) { return (1); }

  /* Print out the problem size, solution parameters, initial guess. */
  PrintHeader(globalstrategy, fnormtol, scsteptol);

  retval = func(cc, fval, data);

  /* Call KINSol and print output concentration profile */
  retval = KINSol(kmem,           /* KINSol memory block */
                  cc,             /* initial guess on input; solution vector */
                  globalstrategy, /* global strategy choice */
                  sc,             /* scaling vector for the variable cc */
                  sc);            /* scaling vector for function values fval */
  if (check_retval(&retval, "KINSol", 1)) { return (1); }

  retval = func(cc, fval, data);
  printf("\n\nComputed equilibrium species concentrations:\n");
  PrintOutput(cc);

  /* Print final statistics and free memory */
  PrintFinalStats(kmem);

  N_VDestroy(cc);
  N_VDestroy(sc);
  KINFree(&kmem);
  SUNLinSolFree(LS);
  FreeUserData(data);
  SUNContext_Free(&sunctx);

  return (0);
}

/* Readability definitions used in other routines below */

#define acoef (data->acoef)
#define bcoef (data->bcoef)
#define cox   (data->cox)
#define coy   (data->coy)

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY KINSOL
 *--------------------------------------------------------------------
 */

/*
 * System function for predator-prey system
 */

static int func(N_Vector cc, N_Vector fval, void* user_data)
{
  sunrealtype xx, yy, delx, dely, *cxy, *rxy, *fxy, dcyli, dcyui, dcxli, dcxri;
  sunindextype jx, jy, is, idyu, idyl, idxr, idxl;
  UserData data;

  data = (UserData)user_data;
  delx = data->dx;
  dely = data->dy;

  /* Loop over all mesh points, evaluating rate array at each point*/
  for (jy = 0; jy < MY; jy++)
  {
    yy = dely * jy;

    /* Set lower/upper index shifts, special at boundaries. */
    idyl = (jy != 0) ? NSMX : -NSMX;
    idyu = (jy != MY - 1) ? NSMX : -NSMX;

    for (jx = 0; jx < MX; jx++)
    {
      xx = delx * jx;

      /* Set left/right index shifts, special at boundaries. */
      idxl = (jx != 0) ? NUM_SPECIES : -NUM_SPECIES;
      idxr = (jx != MX - 1) ? NUM_SPECIES : -NUM_SPECIES;

      cxy = IJ_Vptr(cc, jx, jy);
      rxy = IJ_Vptr(data->rates, jx, jy);
      fxy = IJ_Vptr(fval, jx, jy);

      /* Get species interaction rate array at (xx,yy) */
      WebRate(xx, yy, cxy, rxy, user_data);

      for (is = 0; is < NUM_SPECIES; is++)
      {
        /* Differencing in x direction */
          dcyli = *(cxy + is) - *(cxy - idyl + is);
          dcyui = *(cxy + idyu + is) - *(cxy + is);


        /* Differencing in y direction */
          dcxli = *(cxy + is) - *(cxy - idxl + is);
          dcxri = *(cxy + idxr + is) - *(cxy + is);

        /* Compute the total rate value at (xx,yy) */
        fxy[is] = (coy)[is] * (dcyui - dcyli) + (cox)[is] * (dcxri - dcxli) +
                  rxy[is];

      } /* end of is loop */

    } /* end of jx loop */

  } /* end of jy loop */

  return (0);
}
/*
 * Interaction rate function routine
 */

static void WebRate(sunrealtype xx, sunrealtype yy, sunrealtype* cxy,
                    sunrealtype* ratesxy, void* user_data)
{
  sunindextype i;
  sunrealtype fac;
  UserData data;

  data = (UserData)user_data;

  for (i = 0; i < NUM_SPECIES; i++)
  {
    ratesxy[i] = DotProd(NUM_SPECIES, cxy, acoef[i]);
  }

  fac = ONE + ALPHA * xx * yy;

  for (i = 0; i < NUM_SPECIES; i++)
  {
    ratesxy[i] = cxy[i] * (bcoef[i] * fac + ratesxy[i]);
  }
}

/*
 * Dot product routine for sunrealtype arrays
 */

static sunrealtype DotProd(sunindextype size, sunrealtype* x1, sunrealtype* x2)
{
  sunindextype i;
  sunrealtype *xx1, *xx2, temp = ZERO;

  xx1 = x1;
  xx2 = x2;
  for (i = 0; i < size; i++) { temp += (*xx1++) * (*xx2++); }

  return (temp);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Allocate memory for data structure of type UserData
 */

static UserData AllocUserData(void)
{
  int jx, jy;
  UserData data;

  data = (UserData)malloc(sizeof *data);

  for (jx = 0; jx < MX; jx++)
  {
    for (jy = 0; jy < MY; jy++)
    {
      (data->P)[jx][jy]     = SUNDlsMat_newDenseMat(NUM_SPECIES, NUM_SPECIES);
      (data->pivot)[jx][jy] = SUNDlsMat_newIndexArray(NUM_SPECIES);
    }
  }

  acoef = SUNDlsMat_newDenseMat(NUM_SPECIES, NUM_SPECIES);
  bcoef = (sunrealtype*)malloc(NUM_SPECIES * sizeof(sunrealtype));
  cox   = (sunrealtype*)malloc(NUM_SPECIES * sizeof(sunrealtype));
  coy   = (sunrealtype*)malloc(NUM_SPECIES * sizeof(sunrealtype));

  return (data);
}

/*
 * Load problem constants in data
 */

static void InitUserData(UserData data)
{
  sunindextype i, j, np;
  sunrealtype *a1, *a2, *a3, *a4, dx2, dy2;

  data->mx        = MX;
  data->my        = MY;
  data->ns        = NUM_SPECIES;
  data->np        = NUM_SPECIES / 2;
  data->ax        = AX;
  data->ay        = AY;
  data->dx        = (data->ax) / (MX - 1);
  data->dy        = (data->ay) / (MY - 1);
  data->uround    = SUN_UNIT_ROUNDOFF;
  data->sqruround = SUNRsqrt(data->uround);

  /* Set up the coefficients a and b plus others found in the equations */
  np = data->np;

  dx2 = (data->dx) * (data->dx);
  dy2 = (data->dy) * (data->dy);

  for (i = 0; i < np; i++)
  {
    a1 = &(acoef[i][np]);
    a2 = &(acoef[i + np][0]);
    a3 = &(acoef[i][0]);
    a4 = &(acoef[i + np][np]);

    /*  Fill in the portion of acoef in the four quadrants, row by row */
    for (j = 0; j < np; j++)
    {
      *a1++ = -GG;
      *a2++ = EE;
      *a3++ = ZERO;
      *a4++ = ZERO;
    }

    /* and then change the diagonal elements of acoef to -AA */
    acoef[i][i]           = -AA;
    acoef[i + np][i + np] = -AA;

    bcoef[i]      = BB;
    bcoef[i + np] = -BB;

    cox[i]      = DPREY / dx2;
    cox[i + np] = DPRED / dx2;

    coy[i]      = DPREY / dy2;
    coy[i + np] = DPRED / dy2;
  }
}

/*
 * Free data memory
 */

static void FreeUserData(UserData data)
{
  int jx, jy;

  for (jx = 0; jx < MX; jx++)
  {
    for (jy = 0; jy < MY; jy++)
    {
      SUNDlsMat_destroyMat((data->P)[jx][jy]);
      SUNDlsMat_destroyArray((data->pivot)[jx][jy]);
    }
  }

  SUNDlsMat_destroyMat(acoef);
  free(bcoef);
  free(cox);
  free(coy);
  N_VDestroy(data->rates);
  free(data);
}

/*
 * Set initial conditions in cc
 */

static void SetInitialProfiles(N_Vector cc, N_Vector sc)
{
  int i, jx, jy;
  sunrealtype *cloc, *sloc;
  sunrealtype ctemp[NUM_SPECIES], stemp[NUM_SPECIES];

  /* Initialize arrays ctemp and stemp used in the loading process */
  for (i = 0; i < NUM_SPECIES / 2; i++)
  {
    ctemp[i] = PREYIN;
    stemp[i] = ONE;
  }
  for (i = NUM_SPECIES / 2; i < NUM_SPECIES; i++)
  {
    ctemp[i] = PREDIN;
    stemp[i] = SUN_RCONST(1.0);
  }

  /* Load initial profiles into cc and sc vector from ctemp and stemp. */
  for (jy = 0; jy < MY; jy++)
  {
    for (jx = 0; jx < MX; jx++)
    {
      cloc = IJ_Vptr(cc, jx, jy);
      sloc = IJ_Vptr(sc, jx, jy);
      for (i = 0; i < NUM_SPECIES; i++)
      {
        cloc[i] = ctemp[i];
        sloc[i] = stemp[i];
      }
    }
  }
}

/*
 * Print first lines of output (problem description)
 */

static void PrintHeader(int globalstrategy,
                        sunrealtype fnormtol, sunrealtype scsteptol)
{
  printf("\nPredator-prey test problem --  KINSol (serial version)\n\n");
  printf("Mesh dimensions = %d X %d\n", MX, MY);
  printf("Number of species = %d\n", NUM_SPECIES);
  printf("Total system size = %d\n\n", NEQ);
  printf("Flag globalstrategy = %d (0 = None, 1 = Linesearch)\n", globalstrategy);
  printf("Preconditioning uses interaction-only block-diagonal matrix\n");
  printf("Positivity constraints imposed on all components \n");
#if defined(SUNDIALS_FLOAT128_PRECISION)
  printf("Tolerance parameters:  fnormtol = %Qg   scsteptol = %Qg\n", fnormtol,
         scsteptol);
#elif defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters:  fnormtol = %Lg   scsteptol = %Lg\n", fnormtol,
         scsteptol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters:  fnormtol = %g   scsteptol = %g\n", fnormtol,
         scsteptol);
#else
  printf("Tolerance parameters:  fnormtol = %g   scsteptol = %g\n", fnormtol,
         scsteptol);
#endif

  printf("\nInitial profile of concentration\n");
#if defined(SUNDIALS_FLOAT128_PRECISION)
  printf("At all mesh points:  %Qg %Qg %Qg   %Qg %Qg %Qg\n", PREYIN, PREYIN,
         PREYIN, PREDIN, PREDIN, PREDIN);
#elif defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At all mesh points:  %Lg %Lg %Lg   %Lg %Lg %Lg\n", PREYIN, PREYIN,
         PREYIN, PREDIN, PREDIN, PREDIN);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At all mesh points:  %g %g %g   %g %g %g\n", PREYIN, PREYIN, PREYIN,
         PREDIN, PREDIN, PREDIN);
#else
  printf("At all mesh points:  %g %g %g   %g %g %g\n", PREYIN, PREYIN, PREYIN,
         PREDIN, PREDIN, PREDIN);
#endif
}

/*
 * Print sampled values of current cc
 */

static void PrintOutput(N_Vector cc)
{
  int is, jx, jy;
  sunrealtype* ct;

  jy = 0;
  jx = 0;
  ct = IJ_Vptr(cc, jx, jy);
  printf("\nAt bottom left:");

  /* Print out lines with up to 6 values per line */
  for (is = 0; is < NUM_SPECIES; is++)
  {
    if ((is % 6) * 6 == is) { printf("\n"); }
#if defined(SUNDIALS_FLOAT128_PRECISION)
    printf(" %Qg", ct[is]);
#elif defined(SUNDIALS_EXTENDED_PRECISION)
    printf(" %Lg", ct[is]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf(" %g", ct[is]);
#else
    printf(" %g", ct[is]);
#endif
  }

  jy = MY - 1;
  jx = MX - 1;
  ct = IJ_Vptr(cc, jx, jy);
  printf("\n\nAt top right:");

  /* Print out lines with up to 6 values per line */
  for (is = 0; is < NUM_SPECIES; is++)
  {
    if ((is % 6) * 6 == is) { printf("\n"); }
#if defined(SUNDIALS_FLOAT128_PRECISION)
    printf(" %Qg", ct[is]);
#elif defined(SUNDIALS_EXTENDED_PRECISION)
    printf(" %Lg", ct[is]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf(" %g", ct[is]);
#else
    printf(" %g", ct[is]);
#endif
  }
  printf("\n\n");
}

/*
 * Print final statistics contained in iopt
 */

static void PrintFinalStats(void* kmem)
{
  long int nni, nfe, nli, npe, nps, ncfl, nfeSG;
  int retval;

  retval = KINGetNumNonlinSolvIters(kmem, &nni);
  check_retval(&retval, "KINGetNumNonlinSolvIters", 1);
  retval = KINGetNumFuncEvals(kmem, &nfe);
  check_retval(&retval, "KINGetNumFuncEvals", 1);
  retval = KINGetNumLinIters(kmem, &nli);
  check_retval(&retval, "KINGetNumLinIters", 1);
  retval = KINGetNumPrecEvals(kmem, &npe);
  check_retval(&retval, "KINGetNumPrecEvals", 1);
  retval = KINGetNumPrecSolves(kmem, &nps);
  check_retval(&retval, "KINGetNumPrecSolves", 1);
  retval = KINGetNumLinConvFails(kmem, &ncfl);
  check_retval(&retval, "KINGetNumLinConvFails", 1);
  retval = KINGetNumLinFuncEvals(kmem, &nfeSG);
  check_retval(&retval, "KINGetNumLinFuncEvals", 1);

  printf("Final Statistics.. \n");
  printf("nni    = %5ld    nli   = %5ld\n", nni, nli);
  printf("nfe    = %5ld    nfeSG = %5ld\n", nfe, nfeSG);
  printf("nps    = %5ld    npe   = %5ld     ncfl  = %5ld\n", nps, npe, ncfl);
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
