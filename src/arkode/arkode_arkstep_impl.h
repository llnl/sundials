/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *---------------------------------------------------------------
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
 *---------------------------------------------------------------
 * Implementation header file for ARKODE's ARK time stepper
 * module.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_ARKSTEP_IMPL_H
#define _ARKODE_ARKSTEP_IMPL_H

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_mristep.h>

#include "arkode_impl.h"
#include "arkode_ls_impl.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  ARK time step module constants
  ===============================================================*/

/* max number of nonlinear iterations */
#define MAXCOR 3
/* constant to estimate the convergence rate for the nonlinear equation */
#define CRDOWN SUN_RCONST(0.3)
/* if |gamma/gammap-1| > DGMAX then call lsetup */
#define DGMAX SUN_RCONST(0.2)
/* declare divergence if ratio del/delp > RDIV */
#define RDIV SUN_RCONST(2.3)
/* max no. of steps between lsetup calls */
#define MSBP 20

/* Default solver tolerance factor */
/* #define NLSCOEF   SUN_RCONST(0.003) */ /* Hairer & Wanner constant */
/* #define NLSCOEF   SUN_RCONST(0.2)   */ /* CVODE constant */
#define NLSCOEF SUN_RCONST(0.1)

/* Mass matrix types */
#define MASS_IDENTITY 0
#define MASS_FIXED    1
#define MASS_TIMEDEP  2

/*===============================================================
  ARK time step module data structure
  ===============================================================*/

/*---------------------------------------------------------------
  Types : struct ARKodeARKStepMemRec, ARKodeARKStepMem
  ---------------------------------------------------------------
  The type ARKodeARKStepMem is type pointer to struct
  ARKodeARKStepMemRec.  This structure contains fields to
  perform an additive Runge-Kutta time step.
  ---------------------------------------------------------------*/
typedef struct ARKodeARKStepMemRec
{
  /* ARK problem specification */
  ARKRhsFn fe; /* My' = fe(t,y) + fi(t,y) */
  ARKRhsFn fi;
  sunbooleantype autonomous;     /* SUNTRUE if fi depends on t     */
  sunbooleantype linear;         /* SUNTRUE if fi is linear        */
  sunbooleantype linear_timedep; /* SUNTRUE if dfi/dy depends on t */
  sunbooleantype explicit;       /* SUNTRUE if fe is enabled       */
  sunbooleantype implicit;       /* SUNTRUE if fi is enabled       */
  sunbooleantype deduce_rhs;     /* SUNTRUE if fi is deduced after
                                   a nonlinear solve               */

  /* Adjoint problem specification */
  SUNAdjRhsFn adj_fe;

  /* ARK method storage and parameters */
  N_Vector* Fe;          /* explicit RHS at each stage */
  N_Vector* Fi;          /* implicit RHS at each stage */
  N_Vector* z;           /* stages (for relaxation)    */
  N_Vector sdata;        /* old stage data in residual */
  N_Vector zpred;        /* predicted stage solution   */
  N_Vector zcor;         /* stage correction           */
  int q;                 /* method order               */
  int p;                 /* embedding order            */
  int istage;            /* current stage              */
  int stages;            /* number of stages           */
  ARKodeButcherTable Be; /* ERK Butcher table          */
  ARKodeButcherTable Bi; /* IRK Butcher table          */

  /* User-supplied stage predictor routine */
  ARKStagePredictFn stage_predict;

  /* (Non)Linear solver parameters & data */
  SUNNonlinearSolver NLS; /* generic SUNNonlinearSolver object     */
  sunbooleantype ownNLS;  /* flag indicating ownership of NLS      */
  ARKRhsFn nls_fi;        /* fi(t,y) used in the nonlinear solver  */
  sunrealtype gamma;      /* gamma = h * A(i,i)                       */
  sunrealtype gammap;     /* gamma at the last setup call             */
  sunrealtype gamrat;     /* gamma / gammap                           */
  sunrealtype dgmax;      /* call lsetup if |gamma/gammap-1| >= dgmax */

  int predictor;       /* implicit prediction method to use        */
  sunrealtype crdown;  /* nonlinear conv rate estimation constant  */
  sunrealtype rdiv;    /* nonlin divergence if del/delp > rdiv     */
  sunrealtype crate;   /* estimated nonlin convergence rate        */
  sunrealtype delp;    /* norm of previous nonlinear solver update */
  sunrealtype eRNrm;   /* estimated residual norm, used in nonlin
                            and linear solver convergence tests      */
  sunrealtype nlscoef; /* coefficient in nonlin. convergence test  */

  int msbp;       /* positive => max # steps between lsetup
                            negative => call at each Newton iter     */
  long int nstlp; /* step number of last setup call           */

  int maxcor; /* max num iterations for solving the
                            nonlinear equation                       */

  int convfail;         /* NLS fail flag (for interface routines)   */
  sunbooleantype jcur;  /* is Jacobian info for lin solver current? */
  N_Vector fn_implicit; /* alias to saved implicit function evaluation */

  /* Linear Solver Data */
  ARKLinsolInitFn linit;
  ARKLinsolSetupFn lsetup;
  ARKLinsolSolveFn lsolve;
  ARKLinsolFreeFn lfree;
  void* lmem;
  SUNLinearSolver_Type lsolve_type;

  /* Mass matrix solver data */
  ARKMassInitFn minit;
  ARKMassSetupFn msetup;
  ARKMassMultFn mmult;
  ARKMassSolveFn msolve;
  ARKMassFreeFn mfree;
  void* mass_mem;
  int mass_type; /* 0=identity, 1=fixed, 2=time-dep */
  SUNLinearSolver_Type msolve_type;

  /* Counters */
  long int nfe;       /* num fe calls               */
  long int nfi;       /* num fi calls               */
  long int nsetups;   /* num setup calls            */
  long int nls_iters; /* num nonlinear solver iters */
  long int nls_fails; /* num nonlinear solver fails */

  /* Reusable arrays for fused vector operations */
  sunrealtype* cvals; /* scalar array for fused ops       */
  N_Vector* Xvecs;    /* array of vectors for fused ops   */
  int nfusedopvecs;   /* length of cvals and Xvecs arrays */

  /* Data for using ARKStep with external polynomial forcing */
  sunbooleantype expforcing; /* add forcing to explicit RHS */
  sunbooleantype impforcing; /* add forcing to implicit RHS */
  sunrealtype tshift;        /* time normalization shift    */
  sunrealtype tscale;        /* time normalization scaling  */
  N_Vector* forcing;         /* array of forcing vectors    */
  int nforcing;              /* number of forcing vectors   */
  sunrealtype* stage_times;  /* workspace for applying forcing */
  sunrealtype* stage_coefs;  /* workspace for applying forcing */

}* ARKodeARKStepMem;

/*===============================================================
  ARK time step module private function prototypes
  ===============================================================*/

/* Interface routines supplied to ARKODE */
int arkStep_AttachLinsol(ARKodeMem ark_mem, ARKLinsolInitFn linit,
                         ARKLinsolSetupFn lsetup, ARKLinsolSolveFn lsolve,
                         ARKLinsolFreeFn lfree,
                         SUNLinearSolver_Type lsolve_type, void* lmem);
int arkStep_AttachMasssol(ARKodeMem ark_mem, ARKMassInitFn minit,
                          ARKMassSetupFn msetup, ARKMassMultFn mmult,
                          ARKMassSolveFn msolve, ARKMassFreeFn lfree,
                          sunbooleantype time_dep,
                          SUNLinearSolver_Type msolve_type, void* mass_mem);
void arkStep_DisableLSetup(ARKodeMem ark_mem);
void arkStep_DisableMSetup(ARKodeMem ark_mem);
int arkStep_Init(ARKodeMem ark_mem, sunrealtype tout, int init_type);
void* arkStep_GetLmem(ARKodeMem ark_mem);
void* arkStep_GetMassMem(ARKodeMem ark_mem);
ARKRhsFn arkStep_GetImplicitRHS(ARKodeMem ark_mem);
int arkStep_GetGammas(ARKodeMem ark_mem, sunrealtype* gamma, sunrealtype* gamrat,
                      sunbooleantype** jcur, sunbooleantype* dgamma_fail);
int arkStep_FullRHS(ARKodeMem ark_mem, sunrealtype t, N_Vector y, N_Vector f,
                    int mode);
int arkStep_TakeStep_ERK_Adjoint(ARKodeMem ark_mem, sunrealtype* dsmPtr,
                                 int* nflagPtr);
int arkStep_TakeStep_Z(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr);
int arkStep_SetOptions(ARKodeMem ark_mem, int* argidx, char* argv[],
                       size_t offset, sunbooleantype* arg_used);
int arkStep_SetUserData(ARKodeMem ark_mem, void* user_data);
int arkStep_SetDefaults(ARKodeMem ark_mem);
int arkStep_SetOrder(ARKodeMem ark_mem, int ord);
int arkStep_SetNonlinearSolver(ARKodeMem ark_mem, SUNNonlinearSolver NLS);
int arkStep_SetNlsRhsFn(ARKodeMem ark_mem, ARKRhsFn nls_fi);
int arkStep_SetLinear(ARKodeMem ark_mem, int timedepend);
int arkStep_SetNonlinear(ARKodeMem ark_mem);
int arkStep_SetAutonomous(ARKodeMem ark_mem, sunbooleantype autonomous);
int arkStep_SetNonlinCRDown(ARKodeMem ark_mem, sunrealtype crdown);
int arkStep_SetNonlinRDiv(ARKodeMem ark_mem, sunrealtype rdiv);
int arkStep_SetDeltaGammaMax(ARKodeMem ark_mem, sunrealtype dgmax);
int arkStep_SetLSetupFrequency(ARKodeMem ark_mem, int msbp);
int arkStep_SetPredictorMethod(ARKodeMem ark_mem, int pred_method);
int arkStep_SetMaxNonlinIters(ARKodeMem ark_mem, int maxcor);
int arkStep_SetNonlinConvCoef(ARKodeMem ark_mem, sunrealtype nlscoef);
int arkStep_SetStagePredictFn(ARKodeMem ark_mem, ARKStagePredictFn PredictStage);
int arkStep_SetDeduceImplicitRhs(ARKodeMem ark_mem, sunbooleantype deduce);
int arkStep_GetNumRhsEvals(ARKodeMem ark_mem, int partition_index,
                           long int* rhs_evals);
int arkStep_GetEstLocalErrors(ARKodeMem ark_mem, N_Vector ele);
int arkStep_GetCurrentGamma(ARKodeMem ark_mem, sunrealtype* gamma);
int arkStep_GetNonlinearSystemData(ARKodeMem ark_mem, sunrealtype* tcur,
                                   N_Vector* zpred, N_Vector* z, N_Vector* Fi,
                                   sunrealtype* gamma, N_Vector* sdata,
                                   void** user_data);
int arkStep_GetNumLinSolvSetups(ARKodeMem ark_mem, long int* nlinsetups);
int arkStep_GetNumNonlinSolvIters(ARKodeMem ark_mem, long int* nniters);
int arkStep_GetNumNonlinSolvConvFails(ARKodeMem ark_mem, long int* nnfails);
int arkStep_GetNonlinSolvStats(ARKodeMem ark_mem, long int* nniters,
                               long int* nnfails);
int arkStep_PrintAllStats(ARKodeMem ark_mem, FILE* outfile, SUNOutputFormat fmt);
int arkStep_WriteParameters(ARKodeMem ark_mem, FILE* fp);
int arkStep_Reset(ARKodeMem ark_mem, sunrealtype tR, N_Vector yR);
int arkStep_Resize(ARKodeMem ark_mem, N_Vector y0, sunrealtype hscale,
                   sunrealtype t0, ARKVecResizeFn resize, void* resize_data);
int arkStep_ComputeState(ARKodeMem ark_mem, N_Vector zcor, N_Vector z);
void arkStep_Free(ARKodeMem ark_mem);
void arkStep_PrintMem(ARKodeMem ark_mem, FILE* outfile);
int arkStep_SetInnerForcing(ARKodeMem ark_mem, sunrealtype tshift,
                            sunrealtype tscale, N_Vector* f, int nvecs);

/* Internal utility routines */
int arkStep_AccessARKODEStepMem(void* arkode_mem, const char* fname,
                                ARKodeMem* ark_mem, ARKodeARKStepMem* step_mem);
int arkStep_AccessStepMem(ARKodeMem ark_mem, const char* fname,
                          ARKodeARKStepMem* step_mem);
int arkStep_SetButcherTables(ARKodeMem ark_mem);
int arkStep_CheckButcherTables(ARKodeMem ark_mem);
int arkStep_Predict(ARKodeMem ark_mem, int istage, N_Vector yguess);
int arkStep_StageSetup(ARKodeMem ark_mem, sunbooleantype implicit);
int arkStep_NlsInit(ARKodeMem ark_mem);
int arkStep_Nls(ARKodeMem ark_mem, int nflag);
int arkStep_SetNlsSysFn(ARKodeMem ark_mem);
int arkStep_ComputeSolutions(ARKodeMem ark_mem, sunrealtype* dsm);
int arkStep_ComputeSolutions_MassFixed(ARKodeMem ark_mem, sunrealtype* dsm);
void arkStep_ApplyForcing(ARKodeARKStepMem step_mem, sunrealtype* stage_times,
                          sunrealtype* stage_coefs, int jmax, int* nvec);

/* private functions passed to nonlinear solver */
int arkStep_NlsResidual_MassIdent(N_Vector zcor, N_Vector r, void* arkode_mem);
int arkStep_NlsResidual_MassIdent_TrivialPredAutonomous(N_Vector zcor, N_Vector r,
                                                        void* arkode_mem);
int arkStep_NlsResidual_MassFixed(N_Vector zcor, N_Vector r, void* arkode_mem);
int arkStep_NlsResidual_MassFixed_TrivialPredAutonomous(N_Vector zcor, N_Vector r,
                                                        void* arkode_mem);
int arkStep_NlsResidual_MassTDep(N_Vector zcor, N_Vector r, void* arkode_mem);
int arkStep_NlsFPFunction_MassIdent(N_Vector zcor, N_Vector g, void* arkode_mem);
int arkStep_NlsFPFunction_MassIdent_TrivialPredAutonomous(N_Vector zcor,
                                                          N_Vector g,
                                                          void* arkode_mem);
int arkStep_NlsFPFunction_MassFixed(N_Vector zcor, N_Vector g, void* arkode_mem);
int arkStep_NlsFPFunction_MassFixed_TrivialPredAutonomous(N_Vector zcor,
                                                          N_Vector g,
                                                          void* arkode_mem);
int arkStep_NlsFPFunction_MassTDep(N_Vector zcor, N_Vector g, void* arkode_mem);
int arkStep_NlsLSetup(sunbooleantype jbad, sunbooleantype* jcur,
                      void* arkode_mem);
int arkStep_NlsLSolve(N_Vector delta, void* arkode_mem);
int arkStep_NlsConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del,
                        sunrealtype tol, N_Vector ewt, void* arkode_mem);

/* private functions for relaxation */
int arkStep_SetRelaxFn(ARKodeMem ark_mem, ARKRelaxFn rfn, ARKRelaxJacFn rjac);
int arkStep_RelaxDeltaE(ARKodeMem ark_mem, ARKRelaxJacFn relax_jac_fn,
                        long int* relax_jac_fn_evals, sunrealtype* delta_e_out);
int arkStep_GetOrder(ARKodeMem ark_mem);

/* private functions for adjoints */
int arkStep_fe_Adj(sunrealtype t, N_Vector sens_partial_stage,
                   N_Vector sens_complete_stage, void* content);

int arkStepCompatibleWithAdjointSolver(ARKodeMem ark_mem,
                                       ARKodeARKStepMem step_mem, int lineno,
                                       const char* fname, const char* filename);

/*===============================================================
  Reusable ARKStep Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
#define MSG_ARKSTEP_NO_MEM "Time step module memory is NULL."
#define MSG_NLS_INIT_FAIL  "The nonlinear solver's init routine failed."

/* Other error messages */
#define MSG_ARK_MISSING_FE                                               \
  "Cannot specify that method is explicit without providing a function " \
  "pointer to fe(t,y)."
#define MSG_ARK_MISSING_FI                                               \
  "Cannot specify that method is implicit without providing a function " \
  "pointer to fi(t,y)."
#define MSG_ARK_MISSING_F                                                      \
  "Cannot specify that method is ImEx without providing function pointers to " \
  "fi(t,y) and fe(t,y)."

#ifdef __cplusplus
}
#endif

#endif
