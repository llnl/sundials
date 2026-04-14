/* -----------------------------------------------------------------
 * Programmer: Radu Serban and Cosmin Petra @ LLNL
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
 * Simulation of a slider-crank mechanism modelled with 3 generalized
 * coordinates: crank angle, connecting bar angle, and slider location.
 * The mechanism moves under the action of a constant horizontal
 * force applied to the connecting rod and a spring-damper connecting
 * the crank and connecting rod.
 *
 * The equations of motion are formulated as a system of stabilized
 * index-2 DAEs (Gear-Gupta-Leimkuhler formulation).
 *
 * IDAS also computes sensitivities with respect to the problem
 * parameters k (spring constant) and c (damper constant) of the
 * kinetic energy:
 *   G = int_t0^tend g(t,y,p) dt,
 * where
 *   g(t,y,p) = 0.5*J1*v1^2 + 0.5*J2*v3^2 + 0.5*m2*v2^2
 *
 * -----------------------------------------------------------------*/

#include <idas/idas.h> /* prototypes for IDA fcts., consts.    */
#include <math.h>
#include <nvector/nvector_serial.h> /* access to serial N_Vector            */
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h> /* defs. of sunrealtype, sunindextype      */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */

#define Ith(v, i) NV_Ith_S(v, i - 1) /* i-th vector component i= 1..NEQ */

/* Problem Constants */

#define NEQ 10
#define NP  2

#define TBEGIN SUN_RCONST(0.0)
#define TEND   SUN_RCONST(10.000)

#define RTOLF SUN_RCONST(1.0e-06)
#define ATOLF SUN_RCONST(1.0e-07)

#define RTOLQ SUN_RCONST(1.0e-06)
#define ATOLQ SUN_RCONST(1.0e-08)

#define RTOLFD SUN_RCONST(1.0e-06)
#define ATOLFD SUN_RCONST(1.0e-08)

#define ZERO SUN_RCONST(0.00)
#define HALF SUN_RCONST(0.50)
#define ONE  SUN_RCONST(1.00)
#define TWO  SUN_RCONST(2.00)
#define FOUR SUN_RCONST(4.00)

typedef struct
{
  sunrealtype a;
  sunrealtype J1, J2, m1, m2;
  sunrealtype l0;
  sunrealtype params[2];
  sunrealtype F;
}* UserData;

static int ressc(sunrealtype tres, N_Vector yy, N_Vector yp, N_Vector resval,
                 void* user_data);
static int rhsQ(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector qdot,
                void* user_data);

static int rhsQS(int Ns, sunrealtype t, N_Vector yy, N_Vector yp, N_Vector* yyS,
                 N_Vector* ypS, N_Vector rrQ, N_Vector* rhsvalQS,
                 void* user_data, N_Vector yytmp, N_Vector yptmp, N_Vector tmpQS);

static void setIC(N_Vector yy, N_Vector yp, UserData data);
static void force(N_Vector yy, sunrealtype* Q, UserData data);

static int PrintFinalStats(void* mem);
static int check_retval(void* returnvalue, const char* funcname, int opt);

/*
 *--------------------------------------------------------------------
 * Main Program
 *--------------------------------------------------------------------
 */

int main(void)
{
  UserData data;

  void* mem;
  N_Vector yy, yp, id, q, *yyS, *ypS, *qS;
  sunrealtype tret;
  sunrealtype pbar[2];
  sunrealtype dp, G, Gm[2], Gp[2];
  int retval, is;
  sunrealtype atolS[NP];
  SUNMatrix A;
  SUNLinearSolver LS;
  SUNContext ctx;

  A  = NULL;
  LS = NULL;

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return 1; }

  id = N_VNew_Serial(NEQ, ctx);
  yy = N_VClone(id);
  yp = N_VClone(id);
  q  = N_VNew_Serial(1, ctx);

  yyS = N_VCloneVectorArray(NP, yy);
  ypS = N_VCloneVectorArray(NP, yp);
  qS  = N_VCloneVectorArray(NP, q);

  data = (UserData)malloc(sizeof *data);

  data->a         = 0.5; /* half-length of crank */
  data->J1        = 1.0; /* crank moment of inertia */
  data->m2        = 1.0; /* mass of connecting rod */
  data->m1        = 1.0;
  data->J2        = 2.0; /* moment of inertia of connecting rod */
  data->params[0] = 1.0; /* spring constant */
  data->params[1] = 1.0; /* damper constant */
  data->l0        = 1.0; /* spring free length */
  data->F         = 1.0; /* external constant force */

  N_VConst(ONE, id);
  NV_Ith_S(id, 9) = ZERO;
  NV_Ith_S(id, 8) = ZERO;
  NV_Ith_S(id, 7) = ZERO;
  NV_Ith_S(id, 6) = ZERO;

  printf("\nSlider-Crank example for IDAS:\n");

  /* Consistent IC*/
  setIC(yy, yp, data);

  for (is = 0; is < NP; is++)
  {
    N_VConst(ZERO, yyS[is]);
    N_VConst(ZERO, ypS[is]);
  }

  /* IDA initialization */
  mem    = IDACreate(ctx);
  retval = IDAInit(mem, ressc, TBEGIN, yy, yp);
  retval = IDASStolerances(mem, RTOLF, ATOLF);
  retval = IDASetUserData(mem, data);
  retval = IDASetId(mem, id);
  retval = IDASetSuppressAlg(mem, SUNTRUE);
  retval = IDASetMaxNumSteps(mem, 20000);

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ, ctx);
  if (check_retval((void*)A, "SUNDenseMatrix", 0)) { return (1); }

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(yy, A, ctx);
  if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) { return (1); }

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolver(mem, LS, A);
  if (check_retval(&retval, "IDASetLinearSolver", 1)) { return (1); }

  retval  = IDASensInit(mem, NP, IDA_SIMULTANEOUS, NULL, yyS, ypS);
  pbar[0] = data->params[0];
  pbar[1] = data->params[1];
  retval  = IDASetSensParams(mem, data->params, pbar, NULL);
  retval  = IDASensEEtolerances(mem);
  IDASetSensErrCon(mem, SUNTRUE);

  N_VConst(ZERO, q);
  retval = IDAQuadInit(mem, rhsQ, q);
  retval = IDAQuadSStolerances(mem, RTOLQ, ATOLQ);
  retval = IDASetQuadErrCon(mem, SUNTRUE);

  for (is = 0; is < NP; is++) { N_VConst(ZERO, qS[is]); }
  retval   = IDAQuadSensInit(mem, rhsQS, qS);
  atolS[0] = atolS[1] = ATOLQ;
  retval              = IDAQuadSensSStolerances(mem, RTOLQ, atolS);
  retval              = IDASetQuadSensErrCon(mem, SUNTRUE);

  /* Perform forward run */
  printf("\nForward integration ... ");

  retval = IDASolve(mem, TEND, &tret, yy, yp, IDA_NORMAL);
  if (check_retval(&retval, "IDASolve", 1)) { return (1); }

  printf("done!\n");

  retval = PrintFinalStats(mem);
  if (check_retval(&retval, "PrintFinalStats", 1)) { return (1); }

  IDAGetQuad(mem, &tret, q);
  printf("--------------------------------------------\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("  G = %24.16Lf\n", Ith(q, 1));
#else
  printf("  G = %24.16f\n", Ith(q, 1));
#endif
  printf("--------------------------------------------\n\n");

  IDAGetQuadSens(mem, &tret, qS);
  printf("-------------F O R W A R D------------------\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("   dG/dp:  %12.4Le %12.4Le\n", Ith(qS[0], 1), Ith(qS[1], 1));
#else
  printf("   dG/dp:  %12.4e %12.4e\n", Ith(qS[0], 1), Ith(qS[1], 1));
#endif
  printf("--------------------------------------------\n\n");

  IDAFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);

  /* Finite differences for dG/dp */
  dp              = 0.00001;
  data->params[0] = ONE;
  data->params[1] = ONE;

  mem = IDACreate(ctx);

  setIC(yy, yp, data);
  retval = IDAInit(mem, ressc, TBEGIN, yy, yp);
  retval = IDASStolerances(mem, RTOLFD, ATOLFD);
  retval = IDASetUserData(mem, data);
  retval = IDASetId(mem, id);
  retval = IDASetSuppressAlg(mem, SUNTRUE);

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ, ctx);
  if (check_retval((void*)A, "SUNDenseMatrix", 0)) { return (1); }

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(yy, A, ctx);
  if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) { return (1); }

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolver(mem, LS, A);
  if (check_retval(&retval, "IDASetLinearSolver", 1)) { return (1); }

  N_VConst(ZERO, q);
  IDAQuadInit(mem, rhsQ, q);
  IDAQuadSStolerances(mem, RTOLQ, ATOLQ);
  IDASetQuadErrCon(mem, SUNTRUE);

  IDASolve(mem, TEND, &tret, yy, yp, IDA_NORMAL);

  IDAGetQuad(mem, &tret, q);
  G = Ith(q, 1);
  /*printf("  G  =%12.6e\n", Ith(q,1));*/

  /******************************
  * BACKWARD for k
  ******************************/
  data->params[0] -= dp;
  setIC(yy, yp, data);

  IDAReInit(mem, TBEGIN, yy, yp);

  N_VConst(ZERO, q);
  IDAQuadReInit(mem, q);

  IDASolve(mem, TEND, &tret, yy, yp, IDA_NORMAL);
  IDAGetQuad(mem, &tret, q);
  Gm[0] = Ith(q, 1);
  /*printf("Gm[0]=%12.6e\n", Ith(q,1));*/

  /****************************
  * FORWARD for k *
  ****************************/
  data->params[0] += (TWO * dp);
  setIC(yy, yp, data);
  IDAReInit(mem, TBEGIN, yy, yp);

  N_VConst(ZERO, q);
  IDAQuadReInit(mem, q);

  IDASolve(mem, TEND, &tret, yy, yp, IDA_NORMAL);
  IDAGetQuad(mem, &tret, q);
  Gp[0] = Ith(q, 1);
  /*printf("Gp[0]=%12.6e\n", Ith(q,1));*/

  /* Backward for c */
  data->params[0] = ONE;
  data->params[1] -= dp;
  setIC(yy, yp, data);
  IDAReInit(mem, TBEGIN, yy, yp);

  N_VConst(ZERO, q);
  IDAQuadReInit(mem, q);

  IDASolve(mem, TEND, &tret, yy, yp, IDA_NORMAL);
  IDAGetQuad(mem, &tret, q);
  Gm[1] = Ith(q, 1);

  /* Forward for c */
  data->params[1] += (TWO * dp);
  setIC(yy, yp, data);
  IDAReInit(mem, TBEGIN, yy, yp);

  N_VConst(ZERO, q);
  IDAQuadReInit(mem, q);

  IDASolve(mem, TEND, &tret, yy, yp, IDA_NORMAL);
  IDAGetQuad(mem, &tret, q);
  Gp[1] = Ith(q, 1);

  IDAFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);

  printf("\n\n   Checking using Finite Differences \n\n");

  printf("---------------BACKWARD------------------\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("   dG/dp:  %12.4Le %12.4Le\n", (G - Gm[0]) / dp, (G - Gm[1]) / dp);
#else
  printf("   dG/dp:  %12.4e %12.4e\n", (G - Gm[0]) / dp, (G - Gm[1]) / dp);
#endif
  printf("-----------------------------------------\n\n");

  printf("---------------FORWARD-------------------\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("   dG/dp:  %12.4Le %12.4Le\n", (Gp[0] - G) / dp, (Gp[1] - G) / dp);
#else
  printf("   dG/dp:  %12.4e %12.4e\n", (Gp[0] - G) / dp, (Gp[1] - G) / dp);
#endif
  printf("-----------------------------------------\n\n");

  printf("--------------CENTERED-------------------\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("   dG/dp:  %12.4Le %12.4Le\n", (Gp[0] - Gm[0]) / (TWO * dp),
         (Gp[1] - Gm[1]) / (TWO * dp));
#else
  printf("   dG/dp:  %12.4e %12.4e\n", (Gp[0] - Gm[0]) / (TWO * dp),
         (Gp[1] - Gm[1]) / (TWO * dp));
#endif
  printf("-----------------------------------------\n\n");

  /* Free memory */
  free(data);

  N_VDestroy(id);
  N_VDestroy(yy);
  N_VDestroy(yp);
  N_VDestroy(q);

  N_VDestroyVectorArray(yyS, NP);
  N_VDestroyVectorArray(ypS, NP);
  N_VDestroyVectorArray(qS, NP);

  SUNContext_Free(&ctx);

  return (0);
}

static void setIC(N_Vector yy, N_Vector yp, UserData data)
{
  sunrealtype pi;
  sunrealtype a, J1, m2, J2;
  sunrealtype q, p, x;
  sunrealtype Q[3];

  N_VConst(ZERO, yy);
  N_VConst(ZERO, yp);

  pi = FOUR * atan(ONE);

  a  = data->a;
  J1 = data->J1;
  m2 = data->m2;
  J2 = data->J2;

  q = pi / TWO;
  p = asin(-a);
  x = cos(p);

  NV_Ith_S(yy, 0) = q;
  NV_Ith_S(yy, 1) = x;
  NV_Ith_S(yy, 2) = p;

  force(yy, Q, data);

  NV_Ith_S(yp, 3) = Q[0] / J1;
  NV_Ith_S(yp, 4) = Q[1] / m2;
  NV_Ith_S(yp, 5) = Q[2] / J2;
}

static void force(N_Vector yy, sunrealtype* Q, UserData data)
{
  sunrealtype a, k, c, l0, F;
  sunrealtype q, x, p;
  sunrealtype qd, xd, pd;
  sunrealtype s1, c1, s2, c2, s21, c21;
  sunrealtype l2, l, ld;
  sunrealtype f, fl;

  a  = data->a;
  k  = data->params[0];
  c  = data->params[1];
  l0 = data->l0;
  F  = data->F;

  q = NV_Ith_S(yy, 0);
  x = NV_Ith_S(yy, 1);
  p = NV_Ith_S(yy, 2);

  qd = NV_Ith_S(yy, 3);
  xd = NV_Ith_S(yy, 4);
  pd = NV_Ith_S(yy, 5);

  s1  = sin(q);
  c1  = cos(q);
  s2  = sin(p);
  c2  = cos(p);
  s21 = s2 * c1 - c2 * s1;
  c21 = c2 * c1 + s2 * s1;

  l2 = x * x - x * (c2 + a * c1) + (ONE + a * a) / FOUR + a * c21 / TWO;
  l  = sqrt(l2);
  ld = TWO * x * xd - xd * (c2 + a * c1) + x * (s2 * pd + a * s1 * qd) -
       a * s21 * (pd - qd) / TWO;
  ld /= TWO * l;

  f  = k * (l - l0) + c * ld;
  fl = f / l;

  Q[0] = -fl * a * (s21 / TWO + x * s1) / TWO;
  Q[1] = fl * (c2 / TWO - x + a * c1 / TWO) + F;
  Q[2] = -fl * (x * s2 - a * s21 / TWO) / TWO - F * s2;
}

static int ressc(sunrealtype tres, N_Vector yy, N_Vector yp, N_Vector rr,
                 void* user_data)
{
  UserData data;
  sunrealtype Q[3];
  sunrealtype a, J1, m2, J2;
  sunrealtype *yval, *ypval, *rval;
  sunrealtype q, x, p;
  sunrealtype qd, xd, pd;
  sunrealtype lam1, lam2, mu1, mu2;
  sunrealtype s1, c1, s2, c2;

  data = (UserData)user_data;

  a  = data->a;
  J1 = data->J1;
  m2 = data->m2;
  J2 = data->J2;

  yval  = N_VGetArrayPointer(yy);
  ypval = N_VGetArrayPointer(yp);
  rval  = N_VGetArrayPointer(rr);

  q = yval[0];
  x = yval[1];
  p = yval[2];

  qd = yval[3];
  xd = yval[4];
  pd = yval[5];

  lam1 = yval[6];
  lam2 = yval[7];

  mu1 = yval[8];
  mu2 = yval[9];

  s1 = sin(q);
  c1 = cos(q);
  s2 = sin(p);
  c2 = cos(p);

  force(yy, Q, data);

  rval[0] = ypval[0] - qd + a * s1 * mu1 - a * c1 * mu2;
  rval[1] = ypval[1] - xd + mu1;
  rval[2] = ypval[2] - pd + s2 * mu1 - c2 * mu2;

  rval[3] = J1 * ypval[3] - Q[0] + a * s1 * lam1 - a * c1 * lam2;
  rval[4] = m2 * ypval[4] - Q[1] + lam1;
  rval[5] = J2 * ypval[5] - Q[2] + s2 * lam1 - c2 * lam2;

  rval[6] = x - c2 - a * c1;
  rval[7] = -s2 - a * s1;

  rval[8] = a * s1 * qd + xd + s2 * pd;
  rval[9] = -a * c1 * qd - c2 * pd;

  return (0);
}

static int rhsQ(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector qdot,
                void* user_data)
{
  sunrealtype v1, v2, v3;
  sunrealtype J1, m2, J2;
  UserData data;

  data = (UserData)user_data;
  J1   = data->J1;
  m2   = data->m2;
  J2   = data->J2;

  v1 = Ith(yy, 4);
  v2 = Ith(yy, 5);
  v3 = Ith(yy, 6);

  Ith(qdot, 1) = HALF * (J1 * v1 * v1 + m2 * v2 * v2 + J2 * v3 * v3);

  return (0);
}

static int rhsQS(int Ns, sunrealtype t, N_Vector yy, N_Vector yp, N_Vector* yyS,
                 N_Vector* ypS, N_Vector rrQ, N_Vector* rhsvalQS,
                 void* user_data, N_Vector yytmp, N_Vector yptmp, N_Vector tmpQS)
{
  sunrealtype v1, v2, v3;
  sunrealtype J1, m2, J2;
  UserData data;
  sunrealtype s1, s2, s3;

  data = (UserData)user_data;

  J1 = data->J1;
  m2 = data->m2;
  J2 = data->J2;

  v1 = Ith(yy, 4);
  v2 = Ith(yy, 5);
  v3 = Ith(yy, 6);

  /* Sensitivities of v. */
  s1 = Ith(yyS[0], 4);
  s2 = Ith(yyS[0], 5);
  s3 = Ith(yyS[0], 6);

  Ith(rhsvalQS[0], 1) = J1 * v1 * s1 + m2 * v2 * s2 + J2 * v3 * s3;

  s1 = Ith(yyS[1], 4);
  s2 = Ith(yyS[1], 5);
  s3 = Ith(yyS[1], 6);

  Ith(rhsvalQS[1], 1) = J1 * v1 * s1 + m2 * v2 * s2 + J2 * v3 * s3;

  return (0);
}

static int PrintFinalStats(void* mem)
{
  int retval;
  long int nst, nni, nnf, nje, nre, nreLS, netf, ncfn;

  retval = IDAGetNumSteps(mem, &nst);
  retval = IDAGetNumResEvals(mem, &nre);
  retval = IDAGetNumJacEvals(mem, &nje);
  retval = IDAGetNumNonlinSolvIters(mem, &nni);
  retval = IDAGetNumErrTestFails(mem, &netf);
  retval = IDAGetNumNonlinSolvConvFails(mem, &nnf);
  retval = IDAGetNumStepSolveFails(mem, &ncfn);
  retval = IDAGetNumLinResEvals(mem, &nreLS);

  printf("\nFinal Run Statistics: \n\n");
  printf("Number of steps                    = %ld\n", nst);
  printf("Number of residual evaluations     = %ld\n", nre + nreLS);
  printf("Number of Jacobian evaluations     = %ld\n", nje);
  printf("Number of nonlinear iterations     = %ld\n", nni);
  printf("Number of error test failures      = %ld\n", netf);
  printf("Number of nonlinear conv. failures = %ld\n", nnf);
  printf("Number of step solver failures     = %ld\n", ncfn);

  return (retval);
}

static int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return (1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  return (0);
}
