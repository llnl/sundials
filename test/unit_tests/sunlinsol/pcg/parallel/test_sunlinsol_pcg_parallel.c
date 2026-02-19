/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds, Ashley Crawford @ UMBC
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
 * This is the testing routine to check the SUNLinSol PCG module
 * implementation.
 * -----------------------------------------------------------------
 */

#include <nvector/nvector_parallel.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_pcg.h>

#include "mpi.h"
#include "test_sunlinsol.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#else
#define GSYM "g"
#endif

/* constants */
#define FIVE     SUN_RCONST(5.0)
#define THOUSAND SUN_RCONST(1000.0)

/* user data structure */
typedef struct
{
  sunindextype Nloc; /* local problem size */
  N_Vector d;        /* matrix diagonal */
  N_Vector s;        /* scaling vector supplied to PCG */
  MPI_Comm comm;     /* communicator object */
  int myid;          /* MPI process ID */
  int nprocs;        /* total number of MPI processes */
} UserData;

/* private functions */
/*    matrix-vector product  */
int ATimes(void* ProbData, N_Vector v, N_Vector z);
/*    preconditioner setup */
int PSetup(void* ProbData);
/*    preconditioner solve */
int PSolve(void* ProbData, N_Vector r, N_Vector z, sunrealtype tol, int lr);
/*    checks function return values  */
static int check_flag(void* flagvalue, const char* funcname, int opt);
/*    uniform random number generator in [0,1] */
static sunrealtype urand(void);

/* global copy of Nloc (for check_vector routine) */
sunindextype local_problem_size;

/* ----------------------------------------------------------------------
 * SUNOCG Linear Solver Testing Routine
 *
 * We run multiple tests to exercise this solver:
 * 1. simple tridiagonal system (no preconditioning)
 * 2. simple tridiagonal system (Jacobi preconditioning)
 * 3. tridiagonal system w/ scale vector s (no preconditioning)
 * 4. tridiagonal system w/ scale vector s (Jacobi preconditioning)
 *
 * Note: We construct a tridiagonal matrix Ahat, a random solution
 *       xhat, and a corresponding rhs vector bhat = Ahat*xhat, such
 *       that each of these is unit-less.  To test scaling, we use
 *       the matrix
 *             A = (S-inverse) Ahat (S-inverse),
 *       solution vector
 *             x = S xhat;
 *       and construct b = A*x.  Hence the linear system has both rows
 *       and columns scaled by (S-inverse), where S is the diagonal
 *       matrix with entries from the vector s, the 'scaling' vector
 *       supplied to PCG having strictly positive entries.
 *
 *       When this is combined with preconditioning, we construct
 *       P \approx (A-inverse) by taking a unit-less preconditioner
 *       Phat \approx (Ahat-inverse), and constructing the operator
 *       P via
 *             P = S Phat S \approx S (Ahat-inverse) S = A-inverse
 *       We apply this via the steps:
 *             z = Pr = S Phat S r
 *       Since both S and Phat are diagonal matrices, this is
 *       equivalent to
 *             z(i) = s(i)^2 Phat(i) r(i)
 * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  int fails    = 0;    /* counter for test failures */
  int passfail = 0;    /* overall passfail flag     */
  SUNLinearSolver LS;  /* linear solver object      */
  N_Vector xhat, x, b; /* test vectors              */
  UserData ProbData;   /* problem data structure    */
  int maxl, print_timing;
  sunindextype i;
  sunrealtype* vecdata;
  double tol;
  SUNContext sunctx;

  /* Set up MPI environment */
  fails = MPI_Init(&argc, &argv);
  if (check_flag(&fails, "MPI_Init", 1)) { return 1; }
  ProbData.comm = MPI_COMM_WORLD;
  fails         = MPI_Comm_size(ProbData.comm, &(ProbData.nprocs));
  if (check_flag(&fails, "MPI_Comm_size", 1)) { return 1; }
  fails = MPI_Comm_rank(ProbData.comm, &(ProbData.myid));
  if (check_flag(&fails, "MPI_Comm_rank", 1)) { return 1; }

  if (SUNContext_Create(ProbData.comm, &sunctx))
  {
    printf("ERROR: SUNContext_Create failed\n");
    return (-1);
  }

  /* check inputs: local problem size, timing flag */
  if (argc < 5)
  {
    printf("ERROR: FOUR (4) Inputs required:\n");
    printf("  Local problem size should be >0\n");
    printf("  Maximum Krylov subspace dimension should be >0\n");
    printf("  Solver tolerance should be >0\n");
    printf("  timing output flag should be 0 or 1 \n");
    return 1;
  }
  ProbData.Nloc      = (sunindextype)atol(argv[1]);
  local_problem_size = ProbData.Nloc;
  if (ProbData.Nloc <= 0)
  {
    printf("ERROR: local problem size must be a positive integer\n");
    return 1;
  }
  maxl = atoi(argv[2]);
  if (maxl <= 0)
  {
    printf(
      "ERROR: Maximum Krylov subspace dimension must be a positive integer\n");
    return 1;
  }
  tol = atof(argv[3]);
  if (tol <= ZERO)
  {
    printf("ERROR: Solver tolerance must be a positive real number\n");
    return 1;
  }
  print_timing = atoi(argv[4]);
  SetTiming(print_timing);

  if (ProbData.myid == 0)
  {
    printf("\nPCG linear solver test:\n");
    printf("  nprocs = %i\n", ProbData.nprocs);
    printf("  local/global problem sizes = %ld/%lld\n", (long int)ProbData.Nloc,
           (long long int)ProbData.nprocs * ProbData.Nloc);
    printf("  Maximum Krylov subspace dimension = %i\n", maxl);
    printf("  Solver Tolerance = %g\n", tol);
    printf("  timing output flag = %i\n\n", print_timing);
  }

  /* Create vectors */
  x = N_VNew_Parallel(ProbData.comm, ProbData.Nloc,
                      ProbData.nprocs * ProbData.Nloc, sunctx);
  if (check_flag(x, "N_VNew_Parallel", 0)) { return 1; }
  xhat = N_VNew_Parallel(ProbData.comm, ProbData.Nloc,
                         ProbData.nprocs * ProbData.Nloc, sunctx);
  if (check_flag(xhat, "N_VNew_Parallel", 0)) { return 1; }
  b = N_VNew_Parallel(ProbData.comm, ProbData.Nloc,
                      ProbData.nprocs * ProbData.Nloc, sunctx);
  if (check_flag(b, "N_VNew_Parallel", 0)) { return 1; }
  ProbData.d = N_VNew_Parallel(ProbData.comm, ProbData.Nloc,
                               ProbData.nprocs * ProbData.Nloc, sunctx);
  if (check_flag(ProbData.d, "N_VNew_Parallel", 0)) { return 1; }
  ProbData.s = N_VNew_Parallel(ProbData.comm, ProbData.Nloc,
                               ProbData.nprocs * ProbData.Nloc, sunctx);
  if (check_flag(ProbData.s, "N_VNew_Parallel", 0)) { return 1; }

  /* Fill xhat vector with uniform random data in [1,2] */
  vecdata = N_VGetArrayPointer(xhat);
  for (i = 0; i < ProbData.Nloc; i++) { vecdata[i] = ONE + urand(); }

  /* Fill Jacobi vector with matrix diagonal */
  N_VConst(FIVE, ProbData.d);

  /* Create PCG linear solver */
  LS = SUNLinSol_PCG(x, SUN_PREC_RIGHT, maxl, sunctx);
  fails += Test_SUNLinSolGetType(LS, SUNLINEARSOLVER_ITERATIVE, ProbData.myid);
  fails += Test_SUNLinSolGetID(LS, SUNLINEARSOLVER_PCG, ProbData.myid);
  fails += Test_SUNLinSolSetATimes(LS, &ProbData, ATimes, ProbData.myid);
  fails += Test_SUNLinSolSetPreconditioner(LS, &ProbData, PSetup, PSolve,
                                           ProbData.myid);
  fails += Test_SUNLinSolSetScalingVectors(LS, ProbData.s, NULL, ProbData.myid);
  fails += Test_SUNLinSolSetZeroGuess(LS, ProbData.myid);
  fails += Test_SUNLinSolInitialize(LS, ProbData.myid);
  fails += Test_SUNLinSolSpace(LS, ProbData.myid);
  if (fails)
  {
    printf("FAIL: SUNLinSol_PCG module failed %i initialization tests\n\n",
           fails);
    return 1;
  }
  else if (ProbData.myid == 0)
  {
    printf("SUCCESS: SUNLinSol_PCG module passed all initialization tests\n\n");
  }

  /*** Test 1: simple Poisson-like solve (no preconditioning) ***/

  /* set scaling vector */
  N_VConst(ONE, ProbData.s);

  /* Fill x vector with scaled version */
  N_VProd(xhat, ProbData.s, x);

  /* Fill b vector with result of matrix-vector product */
  fails = ATimes(&ProbData, x, b);
  if (check_flag(&fails, "ATimes", 1)) { return 1; }

  /* Run test with this setup */
  fails += SUNLinSol_PCGSetPrecType(LS, SUN_PREC_NONE);
  fails += Test_SUNLinSolSetup(LS, NULL, ProbData.myid);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNTRUE, ProbData.myid);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNFALSE, ProbData.myid);
  fails += Test_SUNLinSolLastFlag(LS, ProbData.myid);
  fails += Test_SUNLinSolNumIters(LS, ProbData.myid);
  fails += Test_SUNLinSolResNorm(LS, ProbData.myid);
  fails += Test_SUNLinSolResid(LS, ProbData.myid);

  /* Print result */
  if (fails)
  {
    printf("FAIL: SUNLinSol_PCG module, problem 1, failed %i tests\n\n", fails);
    passfail += 1;
  }
  else if (ProbData.myid == 0)
  {
    printf("SUCCESS: SUNLinSol_PCG module, problem 1, passed all tests\n\n");
  }

  /*** Test 2: simple Poisson-like solve (Jacobi preconditioning) ***/

  /* set scaling vector */
  N_VConst(ONE, ProbData.s);

  /* Fill x vector with scaled version */
  N_VProd(xhat, ProbData.s, x);

  /* Fill b vector with result of matrix-vector product */
  fails = ATimes(&ProbData, x, b);
  if (check_flag(&fails, "ATimes", 1)) { return 1; }

  /* Run tests with this setup */
  fails += SUNLinSol_PCGSetPrecType(LS, SUN_PREC_RIGHT);
  fails += Test_SUNLinSolSetup(LS, NULL, ProbData.myid);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNTRUE, ProbData.myid);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNFALSE, ProbData.myid);
  fails += Test_SUNLinSolLastFlag(LS, ProbData.myid);
  fails += Test_SUNLinSolNumIters(LS, ProbData.myid);
  fails += Test_SUNLinSolResNorm(LS, ProbData.myid);
  fails += Test_SUNLinSolResid(LS, ProbData.myid);

  /* Print result */
  if (fails)
  {
    printf("FAIL: SUNLinSol_PCG module, problem 2, failed %i tests\n\n", fails);
    passfail += 1;
  }
  else if (ProbData.myid == 0)
  {
    printf("SUCCESS: SUNLinSol_PCG module, problem 2, passed all tests\n\n");
  }

  /*** Test 3: Poisson-like solve w/ scaling (no preconditioning) ***/

  /* set scaling vector */
  vecdata = N_VGetArrayPointer(ProbData.s);
  for (i = 0; i < ProbData.Nloc; i++) { vecdata[i] = ONE + THOUSAND * urand(); }

  /* Fill x vector with scaled version */
  N_VProd(xhat, ProbData.s, x);

  /* Fill b vector with result of matrix-vector product */
  fails = ATimes(&ProbData, x, b);
  if (check_flag(&fails, "ATimes", 1)) { return 1; }

  /* Run tests with this setup */
  fails += SUNLinSol_PCGSetPrecType(LS, SUN_PREC_NONE);
  fails += Test_SUNLinSolSetup(LS, NULL, ProbData.myid);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNTRUE, ProbData.myid);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNFALSE, ProbData.myid);
  fails += Test_SUNLinSolLastFlag(LS, ProbData.myid);
  fails += Test_SUNLinSolNumIters(LS, ProbData.myid);
  fails += Test_SUNLinSolResNorm(LS, ProbData.myid);
  fails += Test_SUNLinSolResid(LS, ProbData.myid);

  /* Print result */
  if (fails)
  {
    printf("FAIL: SUNLinSol_PCG module, problem 3, failed %i tests\n\n", fails);
    passfail += 1;
  }
  else if (ProbData.myid == 0)
  {
    printf("SUCCESS: SUNLinSol_PCG module, problem 3, passed all tests\n\n");
  }

  /*** Test 4: Poisson-like solve w/ scaling (Jacobi preconditioning) ***/

  /* set scaling vectors */
  vecdata = N_VGetArrayPointer(ProbData.s);
  for (i = 0; i < ProbData.Nloc; i++) { vecdata[i] = ONE + THOUSAND * urand(); }

  /* Fill x vector with scaled version */
  N_VProd(xhat, ProbData.s, x);

  /* Fill b vector with result of matrix-vector product */
  fails = ATimes(&ProbData, x, b);
  if (check_flag(&fails, "ATimes", 1)) { return 1; }

  /* Run tests with this setup */
  fails += SUNLinSol_PCGSetPrecType(LS, SUN_PREC_RIGHT);
  fails += Test_SUNLinSolSetup(LS, NULL, ProbData.myid);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNTRUE, ProbData.myid);
  fails += Test_SUNLinSolSolve(LS, NULL, x, b, tol, SUNFALSE, ProbData.myid);
  fails += Test_SUNLinSolLastFlag(LS, ProbData.myid);
  fails += Test_SUNLinSolNumIters(LS, ProbData.myid);
  fails += Test_SUNLinSolResNorm(LS, ProbData.myid);
  fails += Test_SUNLinSolResid(LS, ProbData.myid);

  /* Print result */
  if (fails)
  {
    printf("FAIL: SUNLinSol_PCG module, problem 4, failed %i tests\n\n", fails);
    passfail += 1;
  }
  else if (ProbData.myid == 0)
  {
    printf("SUCCESS: SUNLinSol_PCG module, problem 4, passed all tests\n\n");
  }

  /* check if any other process failed */
  (void)MPI_Allreduce(&passfail, &fails, 1, MPI_INT, MPI_MAX, ProbData.comm);

  /* Free solver and vectors */
  SUNLinSolFree(LS);
  N_VDestroy(x);
  N_VDestroy(xhat);
  N_VDestroy(b);
  N_VDestroy(ProbData.d);
  N_VDestroy(ProbData.s);
  SUNContext_Free(&sunctx);

  MPI_Finalize();
  return (fails);
}

/* ----------------------------------------------------------------------
 * Private helper functions
 * --------------------------------------------------------------------*/

/* matrix-vector product  */
int ATimes(void* Data, N_Vector v_vec, N_Vector z_vec)
{
  /* local variables */
  sunrealtype *v, *z, *s, vL, vR, vsL, vsR;
  sunindextype i, Nloc;
  int ierr;
  UserData* ProbData;
  MPI_Request SendReqL, SendReqR, RecvReqL, RecvReqR;
  MPI_Status stat;

  /* access user data structure and vector data */
  ProbData = (UserData*)Data;
  v        = N_VGetArrayPointer(v_vec);
  if (check_flag(v, "N_VGetArrayPointer", 0)) { return 1; }
  z = N_VGetArrayPointer(z_vec);
  if (check_flag(z, "N_VGetArrayPointer", 0)) { return 1; }
  s = N_VGetArrayPointer(ProbData->s);
  if (check_flag(s, "N_VGetArrayPointer", 0)) { return 1; }
  Nloc = ProbData->Nloc;

  /* send/recv boundary data with neighbors */
  vL = vR = ZERO;
  vsL     = v[0] / s[0];
  vsR     = v[Nloc - 1] / s[Nloc - 1];
  if (ProbData->myid > 0)
  { /* left neighbor exists */
    ierr = MPI_Irecv(&vL, 1, MPI_SUNREALTYPE, ProbData->myid - 1, MPI_ANY_TAG,
                     ProbData->comm, &RecvReqL);
    if (ierr != MPI_SUCCESS) { return 1; }
    ierr = MPI_Isend(&vsL, 1, MPI_SUNREALTYPE, ProbData->myid - 1, 0,
                     ProbData->comm, &SendReqL);
    if (ierr != MPI_SUCCESS) { return 1; }
  }
  if (ProbData->myid < ProbData->nprocs - 1)
  { /* right neighbor exists */
    ierr = MPI_Irecv(&vR, 1, MPI_SUNREALTYPE, ProbData->myid + 1, MPI_ANY_TAG,
                     ProbData->comm, &RecvReqR);
    if (ierr != MPI_SUCCESS) { return 1; }
    ierr = MPI_Isend(&vsR, 1, MPI_SUNREALTYPE, ProbData->myid + 1, 1,
                     ProbData->comm, &SendReqR);
    if (ierr != MPI_SUCCESS) { return 1; }
  }

  /* iterate through interior of local domain, performing product */
  for (i = 1; i < Nloc - 1; i++)
  {
    z[i] = (-v[i - 1] / s[i - 1] + FIVE * v[i] / s[i] - v[i + 1] / s[i + 1]) /
           s[i];
  }

  /* wait on neighbor data to arrive */
  if (ProbData->myid > 0)
  { /* left neighbor exists */
    ierr = MPI_Wait(&RecvReqL, &stat);
    if (ierr != MPI_SUCCESS) { return 1; }
  }
  if (ProbData->myid < ProbData->nprocs - 1)
  { /* right neighbor exists */
    ierr = MPI_Wait(&RecvReqR, &stat);
    if (ierr != MPI_SUCCESS) { return 1; }
  }

  /* perform product at subdomain boundaries (note: vL/vR are zero at boundary)*/
  z[0] = (-vL + FIVE * v[0] / s[0] - v[1] / s[1]) / s[0];
  z[Nloc - 1] =
    (-v[Nloc - 2] / s[Nloc - 2] + FIVE * v[Nloc - 1] / s[Nloc - 1] - vR) /
    s[Nloc - 1];

  /* return with success */
  return 0;
}

/* preconditioner setup -- nothing to do here since everything is already stored */
int PSetup(void* Data) { return 0; }

/* preconditioner solve */
int PSolve(void* Data, N_Vector r_vec, N_Vector z_vec, sunrealtype tol, int lr)
{
  /* local variables */
  sunrealtype *r, *z, *d, *s;
  sunindextype i;
  UserData* ProbData;

  /* access user data structure and vector data */
  ProbData = (UserData*)Data;
  r        = N_VGetArrayPointer(r_vec);
  if (check_flag(r, "N_VGetArrayPointer", 0)) { return 1; }
  z = N_VGetArrayPointer(z_vec);
  if (check_flag(z, "N_VGetArrayPointer", 0)) { return 1; }
  d = N_VGetArrayPointer(ProbData->d);
  if (check_flag(d, "N_VGetArrayPointer", 0)) { return 1; }
  s = N_VGetArrayPointer(ProbData->s);
  if (check_flag(s, "N_VGetArrayPointer", 0)) { return 1; }

  /* iterate through domain, performing Jacobi solve */
  for (i = 0; i < ProbData->Nloc; i++) { z[i] = s[i] * s[i] * r[i] / d[i]; }

  /* return with success */
  return 0;
}

/* uniform random number generator */
static sunrealtype urand(void)
{
  return ((sunrealtype)rand() / (sunrealtype)RAND_MAX);
}

/* Check function return value based on "opt" input:
     0:  function allocates memory so check for NULL pointer
     1:  function returns a flag so check for flag != 0 */
static int check_flag(void* flagvalue, const char* funcname, int opt)
{
  int* errflag;

  /* Check if function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL)
  {
    fprintf(stderr, "\nERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return 1;
  }

  /* Check if flag != 0 */
  if (opt == 1)
  {
    errflag = (int*)flagvalue;
    if (*errflag != 0)
    {
      fprintf(stderr, "\nERROR: %s() failed with flag = %d\n\n", funcname,
              *errflag);
      return 1;
    }
  }

  return 0;
}

/* ----------------------------------------------------------------------
 * Implementation-specific 'check' routines
 * --------------------------------------------------------------------*/
int check_vector(N_Vector X, N_Vector Y, sunrealtype tol)
{
  int failure = 0;
  sunindextype i;
  sunrealtype *Xdata, *Ydata, maxerr;

  Xdata = N_VGetArrayPointer(X);
  Ydata = N_VGetArrayPointer(Y);

  /* check vector data */
  for (i = 0; i < local_problem_size; i++)
  {
    failure += SUNRCompareTol(Xdata[i], Ydata[i], tol);
  }

  if (failure > ZERO)
  {
    maxerr = ZERO;
    for (i = 0; i < local_problem_size; i++)
    {
      maxerr = SUNMAX(SUNRabs(Xdata[i] - Ydata[i]) / SUNRabs(Xdata[i]), maxerr);
    }
    printf("check err failure: maxerr = %" GSYM " (tol = %" GSYM ")\n", maxerr,
           tol);
    return (1);
  }
  else { return (0); }
}

void sync_device(void) {}
