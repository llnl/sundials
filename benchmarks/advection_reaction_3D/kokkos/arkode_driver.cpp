/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner, Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
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
 * ---------------------------------------------------------------------------*/

#include "advection_reaction_3D.hpp"
#include "arkode/arkode_arkstep.h"
#include "arkode/arkode_erkstep.h"
#include "rhs3D.hpp"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"
#include "sunnonlinsol/sunnonlinsol_newton.h"

/*
 * Definitions for a custom task local SUNNonlinearSolver
 */

typedef struct
{
  int myid;
  int nprocs;
  long int ncnf;
  MPI_Comm comm;
  SUNNonlinearSolver local_nls;
}* TaskLocalNewton_Content;

/* Content accessor macros */
#define GET_NLS_CONTENT(NLS) ((TaskLocalNewton_Content)(NLS->content))
#define LOCAL_NLS(NLS)       (GET_NLS_CONTENT(NLS)->local_nls)

/* SUNNonlinearSolver constructor */
SUNNonlinearSolver TaskLocalNewton(SUNContext ctx, N_Vector y);

/* --------------------------------------------------------------
 * Evolve functions
 * --------------------------------------------------------------*/

/* Setup ARKODE and evolve problem in time with IMEX method */
int EvolveProblemDIRK(N_Vector y, UserData* udata, UserOptions* uopt)
{
  void* arkode_mem       = NULL; /* empty ARKODE memory structure    */
  SUNNonlinearSolver NLS = NULL; /* empty nonlinear solver structure */
  SUNLinearSolver LS     = NULL; /* empty linear solver structure    */

  sunrealtype t, dtout, tout; /* current/output time data     */
  int retval;                 /* reusable error-checking flag */
  int iout;                   /* output counter               */
  long int nst, nst_a, netf;  /* step stats                   */
  long int nfe, nfi;          /* RHS stats                    */
  long int nni, ncnf;         /* nonlinear solver stats       */
  long int nli, npsol;        /* linear solver stats          */

  /* Additively split methods should not add the advection and reaction terms */
  udata->add_reactions = true;

  /* Create the ARK timestepper module */
  arkode_mem = ARKStepCreate(NULL, AdvectionReaction, uopt->t0, y, udata->ctx);
  if (check_retval((void*)arkode_mem, "ARKStepCreate", 0, udata->myid))
  {
    return 1;
  }

  /* Select the method order */
  retval = ARKodeSetOrder(arkode_mem, uopt->order);
  if (check_retval(&retval, "ARKodeSetOrder", 1, udata->myid)) { return 1; }

  /* Attach user data */
  retval = ARKodeSetUserData(arkode_mem, (void*)udata);
  if (check_retval(&retval, "ARKodeSetUserData*", 1, udata->myid)) { return 1; }

  /* Specify tolerances */
  retval = ARKodeSStolerances(arkode_mem, uopt->rtol, uopt->atol);
  if (check_retval(&retval, "ARKodeSStolerances", 1, udata->myid)) { return 1; }

  /* Increase the max number of steps allowed between outputs */
  retval = ARKodeSetMaxNumSteps(arkode_mem, 100000);
  if (check_retval(&retval, "ARKodeSetMaxNumSteps", 1, udata->myid))
  {
    return 1;
  }

  /* Create the (non)linear solver */
  if (uopt->nls == "newton")
  {
    /* Create nonlinear solver */
    NLS = SUNNonlinSol_Newton(y, udata->ctx);
    if (check_retval((void*)NLS, "SUNNonlinSol_Newton", 0, udata->myid))
    {
      return 1;
    }

    /* Attach nonlinear solver */
    retval = ARKodeSetNonlinearSolver(arkode_mem, NLS);
    if (check_retval(&retval, "ARKodeSetNonlinearSolver", 1, udata->myid))
    {
      return 1;
    }

    /* Create linear solver */
    LS = uopt->precond ? SUNLinSol_SPGMR(y, SUN_PREC_LEFT, 0, udata->ctx)
                       : SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, udata->ctx);
    if (check_retval((void*)LS, "SUNLinSol_SPGMR", 0, udata->myid))
    {
      return 1;
    }

    /* Attach linear solver */
    retval = ARKodeSetLinearSolver(arkode_mem, LS, NULL);
    if (check_retval(&retval, "ARKodeSetLinearSolver", 1, udata->myid))
    {
      return 1;
    }

    /* Attach preconditioner */
    retval = ARKodeSetPreconditioner(arkode_mem, NULL, PSolve);
    if (check_retval(&retval, "ARKodeSetPreconditioner", 1, udata->myid))
    {
      return 1;
    }
  }
  else if (uopt->nls == "fixedpoint")
  {
    /* Create nonlinear solver */
    NLS = SUNNonlinSol_FixedPoint(y, uopt->fpaccel, udata->ctx);
    if (check_retval((void*)NLS, "SUNNonlinSol_FixedPoint", 0, udata->myid))
    {
      return 1;
    }

    /* Attach nonlinear solver */
    retval = ARKodeSetNonlinearSolver(arkode_mem, NLS);
    if (check_retval(&retval, "ARKodeSetNonlinearSolver", 1, udata->myid))
    {
      return 1;
    }
  }
  else
  {
    fprintf(stderr, "\nERROR: ARK-DIRK is not compatible with the nls option "
                    "provided\n");
    return 1;
  }

  /* Output initial condition */
  if (uopt->nout > 0)
  {
    if (udata->myid == 0)
    {
      printf("\n          t         ||u||_rms   ||v||_rms   ||w||_rms\n");
      printf("   ----------------------------------------------------\n");
    }
    WriteOutput(uopt->t0, y, udata, uopt);
  }

  /* Integrate to final time */
  t     = uopt->t0;
  dtout = (uopt->tf - uopt->t0);
  if (uopt->nout != 0) { dtout /= uopt->nout; }
  tout = t + dtout;
  iout = 0;

  do {
    /* Integrate to output time */
    retval = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_retval(&retval, "ARKodeEvolve", 1, udata->myid)) { break; }

    /* Output state */
    if (uopt->nout > 0) { WriteOutput(t, y, udata, uopt); }

    /* Update output time */
    tout += dtout;
    tout = (tout > uopt->tf) ? uopt->tf : tout;

    iout++;
  }
  while (iout < uopt->nout);

  /* Get final statistics */
  retval = ARKodeGetNumSteps(arkode_mem, &nst);
  check_retval(&retval, "ARKodeGetNumSteps", 1, udata->myid);
  retval = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  check_retval(&retval, "ARKodeGetNumStepAttempts", 1, udata->myid);
  retval = ARKodeGetNumRhsEvals(arkode_mem, 0, &nfe);
  check_retval(&retval, "ARKodeGetNumRhsEvals", 1, udata->myid);
  retval = ARKodeGetNumRhsEvals(arkode_mem, 1, &nfi);
  check_retval(&retval, "ARKodeGetNumRhsEvals", 1, udata->myid);
  retval = ARKodeGetNumErrTestFails(arkode_mem, &netf);
  check_retval(&retval, "ARKodeGetNumErrTestFails", 1, udata->myid);
  retval = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
  check_retval(&retval, "ARKodeGetNumNonlinSolvIters", 1, udata->myid);
  retval = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncnf);
  check_retval(&retval, "ARKodeGetNumNonlinSolvConvFails", 1, udata->myid);
  if (uopt->nls == "newton")
  {
    retval = ARKodeGetNumLinIters(arkode_mem, &nli);
    check_retval(&retval, "ARKodeGetNumLinIters", 1, udata->myid);
    retval = ARKodeGetNumPrecSolves(arkode_mem, &npsol);
    check_retval(&retval, "ARKodeGetNumPrecSolves", 1, udata->myid);
  }

  /* Print final statistics */
  if (udata->myid == 0)
  {
    printf("\nFinal Solver Statistics (for processor 0):\n");
    printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
    printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi + udata->nnlfi);
    printf("   Total number of error test failures = %li\n", netf);
    printf("   Total number of nonlinear solver convergence failures = %li\n",
           ncnf);
    printf("   Total number of nonlinear iterations = %li\n", nni);
    if (uopt->nls == "newton")
    {
      printf("   Total number of linear iterations = %li\n", nli);
      printf("   Total number of preconditioner solves = %li\n", npsol);
    }
  }

  /* Clean up */
  ARKodeFree(&arkode_mem);
  SUNNonlinSolFree(NLS);
  if (LS) { SUNLinSolFree(LS); }

  /* Return success */
  return (0);
}

/* Setup ARKODE and evolve problem in time with IMEX method */
int EvolveProblemIMEX(N_Vector y, UserData* udata, UserOptions* uopt)
{
  void* arkode_mem       = NULL; /* empty ARKODE memory structure    */
  SUNNonlinearSolver NLS = NULL; /* empty nonlinear solver structure */
  SUNLinearSolver LS     = NULL; /* empty linear solver structure    */

  sunrealtype t, dtout, tout; /* current/output time data     */
  int retval;                 /* reusable error-checking flag */
  int iout;                   /* output counter               */
  long int nst, nst_a, netf;  /* step stats                   */
  long int nfe, nfi;          /* RHS stats                    */
  long int nni, ncnf;         /* nonlinear solver stats       */
  long int nli, npsol;        /* linear solver stats          */

  /* Additively split methods should not add the advection and reaction terms */
  udata->add_reactions = false;

  /* Create the ARK timestepper module */
  arkode_mem = ARKStepCreate(Advection, Reaction, uopt->t0, y, udata->ctx);
  if (check_retval((void*)arkode_mem, "ARKStepCreate", 0, udata->myid))
  {
    return 1;
  }

  /* Select the method order */
  retval = ARKodeSetOrder(arkode_mem, uopt->order);
  if (check_retval(&retval, "ARKodeSetOrder", 1, udata->myid)) { return 1; }

  /* Attach user data */
  retval = ARKodeSetUserData(arkode_mem, (void*)udata);
  if (check_retval(&retval, "ARKodeSetUserData*", 1, udata->myid)) { return 1; }

  /* Specify tolerances */
  retval = ARKodeSStolerances(arkode_mem, uopt->rtol, uopt->atol);
  if (check_retval(&retval, "ARKodeSStolerances", 1, udata->myid)) { return 1; }

  /* Increase the max number of steps allowed between outputs */
  retval = ARKodeSetMaxNumSteps(arkode_mem, 100000);
  if (check_retval(&retval, "ARKodeSetMaxNumSteps", 1, udata->myid))
  {
    return 1;
  }

  /* Create the (non)linear solver */
  if (uopt->nls == "newton")
  {
    /* Create nonlinear solver */
    NLS = SUNNonlinSol_Newton(y, udata->ctx);
    if (check_retval((void*)NLS, "SUNNonlinSol_Newton", 0, udata->myid))
    {
      return 1;
    }

    /* Attach nonlinear solver */
    retval = ARKodeSetNonlinearSolver(arkode_mem, NLS);
    if (check_retval(&retval, "ARKodeSetNonlinearSolver", 1, udata->myid))
    {
      return 1;
    }

    /* Create linear solver */
    LS = SUNLinSol_SPGMR(y, SUN_PREC_LEFT, 0, udata->ctx);
    if (check_retval((void*)LS, "SUNLinSol_SPGMR", 0, udata->myid))
    {
      return 1;
    }

    /* Attach linear solver */
    retval = ARKodeSetLinearSolver(arkode_mem, LS, NULL);
    if (check_retval(&retval, "ARKodeSetLinearSolver", 1, udata->myid))
    {
      return 1;
    }

    /* Attach preconditioner */
    retval = ARKodeSetPreconditioner(arkode_mem, NULL, PSolve);
    if (check_retval(&retval, "ARKodeSetPreconditioner", 1, udata->myid))
    {
      return 1;
    }
  }
  else if (uopt->nls == "tl-newton")
  {
    /* The custom task-local nonlinear solver handles the linear solve
       as well, so we do not need a SUNLinearSolver. */
    NLS = TaskLocalNewton(udata->ctx, y);
    if (check_retval((void*)NLS, "TaskLocalNewton", 0, udata->myid))
    {
      return 1;
    }

    /* Attach nonlinear solver */
    retval = ARKodeSetNonlinearSolver(arkode_mem, NLS);
    if (check_retval(&retval, "ARKodeSetNonlinearSolver", 1, udata->myid))
    {
      return 1;
    }
  }
  else if (uopt->nls == "fixedpoint")
  {
    /* Create nonlinear solver */
    NLS = SUNNonlinSol_FixedPoint(y, uopt->fpaccel, udata->ctx);
    if (check_retval((void*)NLS, "SUNNonlinSol_FixedPoint", 0, udata->myid))
    {
      return 1;
    }

    /* Attach nonlinear solver */
    retval = ARKodeSetNonlinearSolver(arkode_mem, NLS);
    if (check_retval(&retval, "ARKodeSetNonlinearSolver", 1, udata->myid))
    {
      return 1;
    }
  }
  else
  {
    fprintf(stderr, "\nERROR: ARK-IMEX method is not compatible with the nls "
                    "option provided\n");
    return 1;
  }

  /* Output initial condition */
  if (uopt->nout > 0)
  {
    if (udata->myid == 0)
    {
      printf("\n          t         ||u||_rms   ||v||_rms   ||w||_rms\n");
      printf("   ----------------------------------------------------\n");
    }
    WriteOutput(uopt->t0, y, udata, uopt);
  }

  /* Integrate to final time */
  t     = uopt->t0;
  dtout = (uopt->tf - uopt->t0);
  if (uopt->nout != 0) { dtout /= uopt->nout; }
  tout = t + dtout;
  iout = 0;

  do {
    /* Integrate to output time */
    retval = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_retval(&retval, "ARKodeEvolve", 1, udata->myid)) { break; }

    /* Output state */
    if (uopt->nout > 0) { WriteOutput(t, y, udata, uopt); }

    /* Update output time */
    tout += dtout;
    tout = (tout > uopt->tf) ? uopt->tf : tout;

    iout++;
  }
  while (iout < uopt->nout);

  /* Get final statistics */
  retval = ARKodeGetNumSteps(arkode_mem, &nst);
  check_retval(&retval, "ARKodeGetNumSteps", 1, udata->myid);
  retval = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  check_retval(&retval, "ARKodeGetNumStepAttempts", 1, udata->myid);
  retval = ARKodeGetNumRhsEvals(arkode_mem, 0, &nfe);
  check_retval(&retval, "ARKodeGetNumRhsEvals", 1, udata->myid);
  retval = ARKodeGetNumRhsEvals(arkode_mem, 1, &nfi);
  check_retval(&retval, "ARKodeGetNumRhsEvals", 1, udata->myid);
  retval = ARKodeGetNumErrTestFails(arkode_mem, &netf);
  check_retval(&retval, "ARKodeGetNumErrTestFails", 1, udata->myid);
  retval = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
  check_retval(&retval, "ARKodeGetNumNonlinSolvIters", 1, udata->myid);
  retval = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncnf);
  check_retval(&retval, "ARKodeGetNumNonlinSolvConvFails", 1, udata->myid);
  if (uopt->nls == "newton")
  {
    retval = ARKodeGetNumLinIters(arkode_mem, &nli);
    check_retval(&retval, "ARKodeGetNumLinIters", 1, udata->myid);
    retval = ARKodeGetNumPrecSolves(arkode_mem, &npsol);
    check_retval(&retval, "ARKodeGetNumPrecSolves", 1, udata->myid);
  }

  /* Print final statistics */
  if (udata->myid == 0)
  {
    printf("\nFinal Solver Statistics (for processor 0):\n");
    printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
    printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi + udata->nnlfi);
    printf("   Total number of error test failures = %li\n", netf);
    printf("   Total number of nonlinear solver convergence failures = %li\n",
           ncnf);
    printf("   Total number of nonlinear iterations = %li\n", nni);
    if (uopt->nls == "newton")
    {
      printf("   Total number of linear iterations = %li\n", nli);
      printf("   Total number of preconditioner solves = %li\n", npsol);
    }
  }

  /* Clean up */
  ARKodeFree(&arkode_mem);
  if (NLS) { SUNNonlinSolFree(NLS); }
  if (LS) { SUNLinSolFree(LS); }

  /* Return success */
  return (0);
}

/* Setup ARKODE and evolve problem in time explicitly */
int EvolveProblemExplicit(N_Vector y, UserData* udata, UserOptions* uopt)
{
  void* arkode_mem = NULL;    /* empty ARKODE memory structure */
  sunrealtype t, dtout, tout; /* current/output time data      */
  int retval;                 /* reusable error-checking flag  */
  int iout;                   /* output counter                */
  long int nst, nst_a, netf;  /* step stats                    */
  long int nfe;               /* RHS stats                     */

  /* Additively split methods should not add the advection and reaction terms */
  udata->add_reactions = true;

  /* Create the ERK timestepper module */
  arkode_mem = ERKStepCreate(AdvectionReaction, uopt->t0, y, udata->ctx);
  if (check_retval((void*)arkode_mem, "ERKStepCreate", 0, udata->myid))
  {
    return 1;
  }

  /* Select the method order */
  retval = ARKodeSetOrder(arkode_mem, uopt->order);
  if (check_retval(&retval, "ARKodeSetOrder", 1, udata->myid)) { return 1; }

  /* Attach user data */
  retval = ARKodeSetUserData(arkode_mem, (void*)udata);
  if (check_retval(&retval, "ARKodeSetUserData", 1, udata->myid)) { return 1; }

  /* Specify tolerances */
  retval = ARKodeSStolerances(arkode_mem, uopt->rtol, uopt->atol);
  if (check_retval(&retval, "ARKodeSStolerances", 1, udata->myid)) { return 1; }

  /* Increase the max number of steps allowed between outputs */
  retval = ARKodeSetMaxNumSteps(arkode_mem, 1000000);
  if (check_retval(&retval, "ARKodeSetMaxNumSteps", 1, udata->myid))
  {
    return 1;
  }

  /* Set fixed step size */
  retval = ARKodeSetFixedStep(arkode_mem, 1e-5);
  if (check_retval(&retval, "ARKodeSetFixedStep", 1, udata->myid)) { return 1; }

  /* Output initial condition */
  if (uopt->nout > 0)
  {
    if (udata->myid == 0)
    {
      printf("\n          t         ||u||_rms   ||v||_rms   ||w||_rms\n");
      printf("   ----------------------------------------------------\n");
    }
    WriteOutput(uopt->t0, y, udata, uopt);
  }

  /* Integrate to final time */
  t     = uopt->t0;
  dtout = (uopt->tf - uopt->t0);
  if (uopt->nout != 0) { dtout /= uopt->nout; }
  tout = t + dtout;
  iout = 0;

  do {
    /* Integrate to output time */
    retval = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_retval(&retval, "ARKodeEvolve", 1, udata->myid)) { break; }

    /* Output state */
    if (uopt->nout > 0) { WriteOutput(t, y, udata, uopt); }

    /* Update output time */
    tout += dtout;
    tout = (tout > uopt->tf) ? uopt->tf : tout;

    iout++;
  }
  while (iout < uopt->nout);

  /* Get final statistics */
  retval = ARKodeGetNumSteps(arkode_mem, &nst);
  check_retval(&retval, "ARKodeGetNumSteps", 1, udata->myid);
  retval = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  check_retval(&retval, "ARKodeGetNumStepAttempts", 1, udata->myid);
  retval = ARKodeGetNumRhsEvals(arkode_mem, 0, &nfe);
  check_retval(&retval, "ARKodeGetNumRhsEvals", 1, udata->myid);
  retval = ARKodeGetNumErrTestFails(arkode_mem, &netf);
  check_retval(&retval, "ARKodeGetNumErrTestFails", 1, udata->myid);

  /* Print final statistics */
  if (udata->myid == 0)
  {
    printf("\nFinal Solver Statistics (for processor 0):\n");
    printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
    printf("   Total RHS evals:  Fe = %li\n", nfe);
    printf("   Total number of error test failures = %li\n", netf);
  }

  /* Clean up */
  ARKodeFree(&arkode_mem);

  /* Return success */
  return (0);
}

/* --------------------------------------------------------------
 * (Non)linear system functions
 * --------------------------------------------------------------*/

static int TaskLocalNlsResidual(N_Vector ycor, N_Vector F, void* arkode_mem)
{
  /* temporary variables */
  UserData* udata;
  int retval;
  sunrealtype c[3];
  N_Vector X[3];

  /* nonlinear system data */
  N_Vector z, zpred, Fi, sdata;
  sunrealtype tcur, gamma;
  void* user_data;

  ARKodeGetNonlinearSystemData(arkode_mem, &tcur, &zpred, &z, &Fi, &gamma,
                               &sdata, &user_data);
  udata = (UserData*)user_data;

  /* update 'z' value as stored predictor + current corrector */
  N_VLinearSum(1.0, N_VGetLocalVector_MPIPlusX(zpred), 1.0, (ycor),
               N_VGetLocalVector_MPIPlusX(z));

  /* compute implicit RHS and save for later */
  retval = Reaction(tcur, N_VGetLocalVector_MPIPlusX(z),
                    N_VGetLocalVector_MPIPlusX(Fi), user_data);
  udata->nnlfi++; /* count calls to Fi as part of the nonlinear residual */
  if (retval < 0) { return (-1); }
  if (retval > 0) { return (+1); }

  /* update with y, sdata, and gamma * fy */
  X[0] = ycor;
  c[0] = 1.0;
  c[1] = -1.0;
  X[1] = N_VGetLocalVector_MPIPlusX(sdata);
  c[2] = -gamma;
  X[2] = N_VGetLocalVector_MPIPlusX(Fi);

  retval = N_VLinearCombination(3, c, X, F);
  if (retval != 0) { return (-1); }

  return (0);
}

static int TaskLocalLSolve(N_Vector delta, void* arkode_mem)
{
  /* local variables */
  UserData* udata = NULL;
  int retval;

  /* nonlinear system data */
  N_Vector z, zpred, Fi, sdata;
  sunrealtype tcur, gamma;
  void* user_data = NULL;

  ARKodeGetNonlinearSystemData(arkode_mem, &tcur, &zpred, &z, &Fi, &gamma,
                               &sdata, &user_data);
  udata = (UserData*)user_data;

  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  /* set up I - gamma*J and solve */
  retval = SolveReactionLinSys(z, delta, delta, gamma, udata);

  return (retval);
}

static SUNNonlinearSolver_Type TaskLocalNewton_GetType(SUNNonlinearSolver NLS)
{
  return SUNNONLINEARSOLVER_ROOTFIND;
}

static int TaskLocalNewton_Initialize(SUNNonlinearSolver NLS)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL) { return SUN_ERR_ARG_CORRUPT; }

  /* override default system and lsolve functions with local versions */
  SUNNonlinSolSetSysFn(LOCAL_NLS(NLS), TaskLocalNlsResidual);
  SUNNonlinSolSetLSolveFn(LOCAL_NLS(NLS), TaskLocalLSolve);

  return (SUNNonlinSolInitialize(LOCAL_NLS(NLS)));
}

static int TaskLocalNewton_Solve(SUNNonlinearSolver NLS, N_Vector y0,
                                 N_Vector ycor, N_Vector w, sunrealtype tol,
                                 sunbooleantype callLSetup, void* mem)
{
  /* local variables */
  MPI_Comm comm;
  int solve_status, recover, nonrecover;

  /* check that the inputs are non-null */
  if ((NLS == NULL) || (y0 == NULL) || (ycor == NULL) || (w == NULL) ||
      (mem == NULL))
  {
    return SUN_ERR_ARG_CORRUPT;
  }

  /* shortcuts */
  comm = GET_NLS_CONTENT(NLS)->comm;

  /* each tasks solves the local nonlinear system */
  solve_status =
    SUNNonlinSolSolve(LOCAL_NLS(NLS), N_VGetLocalVector_MPIPlusX(y0),
                      N_VGetLocalVector_MPIPlusX(ycor),
                      N_VGetLocalVector_MPIPlusX(w), tol, callLSetup, mem);

  /* if any process had a nonrecoverable failure, return it */
  MPI_Allreduce(&solve_status, &nonrecover, 1, MPI_INT, MPI_MIN, comm);
  if (nonrecover < 0) { return nonrecover; }

  /* check if any process has a recoverable convergence failure */
  MPI_Allreduce(&solve_status, &recover, 1, MPI_INT, MPI_MAX, comm);
  if (recover == SUN_NLS_CONV_RECVR) { GET_NLS_CONTENT(NLS)->ncnf++; }

  /* return success (recover == 0) or a recoverable error code (recover > 0) */
  return recover;
}

static int TaskLocalNewton_Free(SUNNonlinearSolver NLS)
{
  /* return if NLS is already free */
  if (NLS == NULL) { return SUN_SUCCESS; }

  /* free items from contents, then the generic structure */
  if (NLS->content)
  {
    SUNNonlinSolFree(LOCAL_NLS(NLS));
    free(NLS->content);
    NLS->content = NULL;
  }

  /* free the ops structure */
  if (NLS->ops)
  {
    free(NLS->ops);
    NLS->ops = NULL;
  }

  /* free the nonlinear solver */
  free(NLS);

  return SUN_SUCCESS;
}

static int TaskLocalNewton_SetSysFn(SUNNonlinearSolver NLS,
                                    SUNNonlinSolSysFn SysFn)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL) { return SUN_ERR_ARG_CORRUPT; }

  return (SUNNonlinSolSetSysFn(LOCAL_NLS(NLS), SysFn));
}

static int TaskLocalNewton_SetConvTestFn(SUNNonlinearSolver NLS,
                                         SUNNonlinSolConvTestFn CTestFn,
                                         void* ctest_data)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL) { return SUN_ERR_ARG_CORRUPT; }

  return (SUNNonlinSolSetConvTestFn(LOCAL_NLS(NLS), CTestFn, ctest_data));
}

static int TaskLocalNewton_GetNumConvFails(SUNNonlinearSolver NLS,
                                           long int* nconvfails)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL) { return SUN_ERR_ARG_CORRUPT; }

  *nconvfails = GET_NLS_CONTENT(NLS)->ncnf;
  return (0);
}

SUNNonlinearSolver TaskLocalNewton(SUNContext ctx, N_Vector y)
{
  SUNNonlinearSolver NLS;
  TaskLocalNewton_Content content;

  /* Check that the supplied N_Vector is non-NULL */
  if (y == NULL) { return NULL; }

  /* Check that the supplied N_Vector is an MPIPlusX */
  if (N_VGetVectorID(y) != SUNDIALS_NVEC_MPIPLUSX) { return NULL; }

  /* Create an empty nonlinear linear solver object */
  NLS = SUNNonlinSolNewEmpty(ctx);
  if (NLS == NULL) { return NULL; }

  /* Attach operations */
  NLS->ops->gettype         = TaskLocalNewton_GetType;
  NLS->ops->initialize      = TaskLocalNewton_Initialize;
  NLS->ops->solve           = TaskLocalNewton_Solve;
  NLS->ops->free            = TaskLocalNewton_Free;
  NLS->ops->setsysfn        = TaskLocalNewton_SetSysFn;
  NLS->ops->setctestfn      = TaskLocalNewton_SetConvTestFn;
  NLS->ops->getnumconvfails = TaskLocalNewton_GetNumConvFails;

  /* Create content */
  content = NULL;
  content = (TaskLocalNewton_Content)malloc(sizeof *content);
  if (content == NULL)
  {
    SUNNonlinSolFree(NLS);
    return NULL;
  }

  /* Initialize all components of content to 0/NULL */
  memset(content, 0, sizeof(*content));

  /* Attach content */
  NLS->content = content;

  /* Fill general content */
  void* tmpcomm = N_VGetCommunicator(y);
  if (tmpcomm == NULL)
  {
    SUNNonlinSolFree(NLS);
    return NULL;
  }

  MPI_Comm* comm = (MPI_Comm*)tmpcomm;
  if ((*comm) == MPI_COMM_NULL)
  {
    SUNNonlinSolFree(NLS);
    return NULL;
  }

  content->comm = *comm;

  content->local_nls = SUNNonlinSol_Newton(N_VGetLocalVector_MPIPlusX(y), ctx);
  if (content->local_nls == NULL)
  {
    SUNNonlinSolFree(NLS);
    return NULL;
  }

  MPI_Comm_rank(content->comm, &content->myid);
  MPI_Comm_size(content->comm, &content->nprocs);

  content->ncnf = 0;

  return NLS;
}
