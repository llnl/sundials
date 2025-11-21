/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ UMBC
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
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
 * -----------------------------------------------------------------------------
 * Implementation header file for ARKODE's MRI time stepper module.
 * ---------------------------------------------------------------------------*/

#ifndef _ARKODE_MRISTEP_IMPL_H
#define _ARKODE_MRISTEP_IMPL_H

/* Public header file */
#include "arkode/arkode_mristep.h"

/* Private header files */
#include "arkode_impl.h"
#include "arkode_ls_impl.h"
#include "arkode_mri_tables_impl.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Stage type identifiers */
#define MRISTAGE_FIRST       -2
#define MRISTAGE_STIFF_ACC   -1
#define MRISTAGE_ERK_FAST    0
#define MRISTAGE_ERK_NOFAST  1
#define MRISTAGE_DIRK_NOFAST 2
#define MRISTAGE_DIRK_FAST   3

/* Implicit solver constants (duplicate from arkode_arkstep_impl.h) */
/*   max number of nonlinear iterations */
#define MAXCOR 3
/*   constant to estimate the convergence rate for the nonlinear equation */
#define CRDOWN SUN_RCONST(0.3)
/*   if |gamma/gammap-1| > DGMAX then call lsetup */
#define DGMAX SUN_RCONST(0.2)
/*   declare divergence if ratio del/delp > RDIV */
#define RDIV SUN_RCONST(2.3)
/*   max no. of steps between lsetup calls */
#define MSBP 20
/*   default solver tolerance factor */
#define NLSCOEF SUN_RCONST(0.1)

/*===============================================================
  MRI time step module data structure
  ===============================================================*/

/*---------------------------------------------------------------
  The type ARKodeMRIStepMem is type pointer to struct
  ARKodeMRIStepMemRec. This structure contains fields to
  perform a MRI time step.
  ---------------------------------------------------------------*/
typedef struct ARKodeMRIStepMemRec
{
  /* MRI problem specification */
  ARKRhsFn fse; /* y' = fse(t,y) + fsi(t,y) + ff(t,y) */
  ARKRhsFn fsi;
  sunbooleantype linear;         /* SUNTRUE if fi is linear        */
  sunbooleantype linear_timedep; /* SUNTRUE if dfi/dy depends on t */
  sunbooleantype explicit_rhs;   /* SUNTRUE if fse is provided     */
  sunbooleantype implicit_rhs;   /* SUNTRUE if fsi is provided     */
  sunbooleantype deduce_rhs;     /* SUNTRUE if fi is deduced after
                                    a nonlinear solve              */

  /* Outer RK method storage and parameters */
  N_Vector* Fse;           /* explicit RHS at each stage               */
  N_Vector* Fsi;           /* implicit RHS at each stage               */
  sunbooleantype unify_Fs; /* Fse and Fsi point at the same memory     */
  sunbooleantype fse_is_current;
  sunbooleantype fsi_is_current;
  MRIStepCoupling MRIC;  /* slow->fast coupling table                */
  int q;                 /* method order                             */
  int p;                 /* embedding order                          */
  int stages;            /* total number of stages                   */
  int nstages_active;    /* number of active stage RHS vectors       */
  int nstages_allocated; /* number of stage RHS vectors allocated    */
  int* stage_map;        /* index map for storing stage RHS vectors  */
  int* stagetypes;       /* type flags for stages                    */
  sunrealtype* Ae_row;   /* equivalent explicit RK coeffs            */
  sunrealtype* Ai_row;   /* equivalent implicit RK coeffs            */

  /* Algebraic solver data and parameters */
  N_Vector sdata;         /* old stage data in residual               */
  N_Vector zpred;         /* predicted stage solution                 */
  N_Vector zcor;          /* stage correction                         */
  int istage;             /* current stage index                      */
  SUNNonlinearSolver NLS; /* generic SUNNonlinearSolver object        */
  sunbooleantype ownNLS;  /* flag indicating ownership of NLS         */
  ARKRhsFn nls_fsi;       /* fsi(t,y) used in the nonlinear solver    */
  sunrealtype gamma;      /* gamma = h * A(i,i)                       */
  sunrealtype gammap;     /* gamma at the last setup call             */
  sunrealtype gamrat;     /* gamma / gammap                           */
  sunrealtype dgmax;      /* call lsetup if |gamma/gammap-1| >= dgmax */
  int predictor;          /* implicit prediction method to use        */
  sunrealtype crdown;     /* nonlinear conv rate estimation constant  */
  sunrealtype rdiv;       /* nonlin divergence if del/delp > rdiv     */
  sunrealtype crate;      /* estimated nonlin convergence rate        */
  sunrealtype delp;       /* norm of previous nonlinear solver update */
  sunrealtype eRNrm;      /* estimated residual norm, used in nonlin
                             and linear solver convergence tests      */
  sunrealtype nlscoef;    /* coefficient in nonlin. convergence test  */

  int msbp;       /* positive => max # steps between lsetup
                     negative => call at each Newton iter     */
  long int nstlp; /* step number of last setup call           */

  int maxcor;          /* max num iterations for solving the
                          nonlinear equation                       */
  int convfail;        /* NLS fail flag (for interface routines)   */
  sunbooleantype jcur; /* is Jacobian info for lin solver current? */
  ARKStagePredictFn stage_predict; /* User-supplied stage predictor   */

  /* Linear Solver Data */
  ARKLinsolInitFn linit;
  ARKLinsolSetupFn lsetup;
  ARKLinsolSolveFn lsolve;
  ARKLinsolFreeFn lfree;
  void* lmem;

  /* Inner stepper */
  MRIStepInnerStepper stepper;

  /* User-supplied pre and post inner evolve functions */
  MRIStepPreInnerFn pre_inner_evolve;
  MRIStepPostInnerFn post_inner_evolve;

  /* MRI adaptivity parameters */
  sunrealtype inner_rtol_factor;     /* prev control parameter */
  sunrealtype inner_dsm;             /* prev inner stepper accumulated error */
  sunrealtype inner_rtol_factor_new; /* upcoming control parameter */

  /* Counters */
  long int nfse;        /* num fse calls                    */
  long int nfsi;        /* num fsi calls                    */
  long int nsetups;     /* num linear solver setup calls    */
  long int nls_iters;   /* num nonlinear solver iters       */
  long int nls_fails;   /* num nonlinear solver fails       */
  long int inner_fails; /* num recov. inner solver fails  */
  int nfusedopvecs;     /* length of cvals and Xvecs arrays */

  /* Data for using MRIStep with external polynomial forcing */
  sunbooleantype expforcing; /* add forcing to explicit RHS */
  sunbooleantype impforcing; /* add forcing to implicit RHS */
  sunrealtype tshift;        /* time normalization shift    */
  sunrealtype tscale;        /* time normalization scaling  */
  N_Vector* forcing;         /* array of forcing vectors    */
  int nforcing;              /* number of forcing vectors   */

  /* Reusable arrays for fused vector operations */
  sunrealtype* cvals;
  N_Vector* Xvecs;

}* ARKodeMRIStepMem;

/*===============================================================
  MRI innter time stepper data structure
  ===============================================================*/

typedef struct _MRIStepInnerStepper_Ops* MRIStepInnerStepper_Ops;

struct _MRIStepInnerStepper_Ops
{
  MRIStepInnerEvolveFn evolve;
  MRIStepInnerFullRhsFn fullrhs;
  MRIStepInnerResetFn reset;
  MRIStepInnerGetAccumulatedError geterror;
  MRIStepInnerResetAccumulatedError reseterror;
  MRIStepInnerSetRTol setrtol;
};

struct _MRIStepInnerStepper
{
  /* stepper specific content and operations */
  void* content;
  void* python;
  MRIStepInnerStepper_Ops ops;

  /* stepper context */
  SUNContext sunctx;

  /* base class data */
  N_Vector* forcing;      /* array of forcing vectors            */
  int nforcing;           /* number of forcing vectors active    */
  int nforcing_allocated; /* number of forcing vectors allocated */
  int last_flag;          /* last stepper return flag            */
  sunrealtype tshift;     /* time normalization shift            */
  sunrealtype tscale;     /* time normalization scaling          */

  /* fused op workspace */
  sunrealtype* vals;
  N_Vector* vecs;

  /* Space requirements */
  sunindextype lrw1; /* no. of sunrealtype words in 1 N_Vector          */
  sunindextype liw1; /* no. of integer words in 1 N_Vector           */
  long int lrw;      /* no. of sunrealtype words in ARKODE work vectors */
  long int liw;      /* no. of integer words in ARKODE work vectors  */
};

/*===============================================================
  MRI time step module private function prototypes
  ===============================================================*/

/* Interface routines supplied to ARKODE */
int mriStep_AttachLinsol(ARKodeMem ark_mem, ARKLinsolInitFn linit,
                         ARKLinsolSetupFn lsetup, ARKLinsolSolveFn lsolve,
                         ARKLinsolFreeFn lfree,
                         SUNLinearSolver_Type lsolve_type, void* lmem);
void mriStep_DisableLSetup(ARKodeMem ark_mem);
int mriStep_Init(ARKodeMem ark_mem, sunrealtype tout, int init_type);
void* mriStep_GetLmem(ARKodeMem ark_mem);
ARKRhsFn mriStep_GetImplicitRHS(ARKodeMem ark_mem);
int mriStep_GetGammas(ARKodeMem ark_mem, sunrealtype* gamma, sunrealtype* gamrat,
                      sunbooleantype** jcur, sunbooleantype* dgamma_fail);
int mriStep_FullRHS(ARKodeMem ark_mem, sunrealtype t, N_Vector y, N_Vector f,
                    int mode);
int mriStep_UpdateF0(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem,
                     sunrealtype t, N_Vector y, int mode);
int mriStep_TakeStepMRIGARK(ARKodeMem ark_mem, sunrealtype* dsmPtr,
                            int* nflagPtr);
int mriStep_TakeStepMRISR(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr);
int mriStep_TakeStepMERK(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr);
int mriStep_SetOptions(ARKodeMem ark_mem, int* argidx, char* argv[],
                       size_t offset, sunbooleantype* arg_used);
int mriStep_SetAdaptController(ARKodeMem ark_mem, SUNAdaptController C);
int mriStep_SetUserData(ARKodeMem ark_mem, void* user_data);
int mriStep_SetDefaults(ARKodeMem ark_mem);
int mriStep_SetOrder(ARKodeMem ark_mem, int ord);
int mriStep_SetNonlinearSolver(ARKodeMem ark_mem, SUNNonlinearSolver NLS);
int mriStep_SetNlsRhsFn(ARKodeMem ark_mem, ARKRhsFn nls_fi);
int mriStep_SetLinear(ARKodeMem ark_mem, int timedepend);
int mriStep_SetNonlinear(ARKodeMem ark_mem);
int mriStep_SetNonlinCRDown(ARKodeMem ark_mem, sunrealtype crdown);
int mriStep_SetNonlinRDiv(ARKodeMem ark_mem, sunrealtype rdiv);
int mriStep_SetDeltaGammaMax(ARKodeMem ark_mem, sunrealtype dgmax);
int mriStep_SetLSetupFrequency(ARKodeMem ark_mem, int msbp);
int mriStep_SetPredictorMethod(ARKodeMem ark_mem, int pred_method);
int mriStep_SetMaxNonlinIters(ARKodeMem ark_mem, int maxcor);
int mriStep_SetNonlinConvCoef(ARKodeMem ark_mem, sunrealtype nlscoef);
int mriStep_SetStagePredictFn(ARKodeMem ark_mem, ARKStagePredictFn PredictStage);
int mriStep_SetDeduceImplicitRhs(ARKodeMem ark_mem, sunbooleantype deduce);
int mriStep_GetEstLocalErrors(ARKodeMem ark_mem, N_Vector ele);
int mriStep_GetNumRhsEvals(ARKodeMem ark_mem, int partition_index,
                           long int* rhs_evals);
int mriStep_GetCurrentGamma(ARKodeMem ark_mem, sunrealtype* gamma);
int mriStep_GetNonlinearSystemData(ARKodeMem ark_mem, sunrealtype* tcur,
                                   N_Vector* zpred, N_Vector* z, N_Vector* Fi,
                                   sunrealtype* gamma, N_Vector* sdata,
                                   void** user_data);
int mriStep_GetNumLinSolvSetups(ARKodeMem ark_mem, long int* nlinsetups);
int mriStep_GetNumNonlinSolvIters(ARKodeMem ark_mem, long int* nniters);
int mriStep_GetNumNonlinSolvConvFails(ARKodeMem ark_mem, long int* nnfails);
int mriStep_GetNonlinSolvStats(ARKodeMem ark_mem, long int* nniters,
                               long int* nnfails);
int mriStep_PrintAllStats(ARKodeMem ark_mem, FILE* outfile, SUNOutputFormat fmt);
int mriStep_WriteParameters(ARKodeMem ark_mem, FILE* fp);
int mriStep_Reset(ARKodeMem ark_mem, sunrealtype tR, N_Vector yR);
int mriStep_Resize(ARKodeMem ark_mem, N_Vector y0, sunrealtype hscale,
                   sunrealtype t0, ARKVecResizeFn resize, void* resize_data);
int mriStep_ComputeState(ARKodeMem ark_mem, N_Vector zcor, N_Vector z);
void mriStep_Free(ARKodeMem ark_mem);
void mriStep_PrintMem(ARKodeMem ark_mem, FILE* outfile);
int mriStep_SetInnerForcing(ARKodeMem ark_mem, sunrealtype tshift,
                            sunrealtype tscale, N_Vector* f, int nvecs);

/* Internal utility routines */
int mriStep_AccessARKODEStepMem(void* arkode_mem, const char* fname,
                                ARKodeMem* ark_mem, ARKodeMRIStepMem* step_mem);
int mriStep_AccessStepMem(ARKodeMem ark_mem, const char* fname,
                          ARKodeMRIStepMem* step_mem);
int mriStep_SetCoupling(ARKodeMem ark_mem);
int mriStep_CheckCoupling(ARKodeMem ark_mem);
int mriStep_StageERKFast(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem,
                         sunrealtype t0, sunrealtype tf, N_Vector ycur,
                         N_Vector ytemp, sunbooleantype get_inner_dsm);
int mriStep_StageERKNoFast(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem, int is);
int mriStep_StageDIRKFast(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem, int is,
                          int* nflagPtr);
int mriStep_StageDIRKNoFast(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem,
                            int is, int* nflagPtr);
int mriStep_Predict(ARKodeMem ark_mem, int istage, N_Vector yguess);
int mriStep_StageSetup(ARKodeMem ark_mem);
int mriStep_NlsInit(ARKodeMem ark_mem);
int mriStep_Nls(ARKodeMem ark_mem, int nflag);
int mriStep_SlowRHS(ARKodeMem ark_mem, sunrealtype t, N_Vector y, N_Vector f,
                    int mode);
int mriStep_Hin(ARKodeMem ark_mem, sunrealtype tcur, sunrealtype tout,
                N_Vector fcur, sunrealtype* h);
void mriStep_ApplyForcing(ARKodeMRIStepMem step_mem, sunrealtype t,
                          sunrealtype s, int* nvec);

/* private functions passed to nonlinear solver */
int mriStep_NlsResidual(N_Vector yy, N_Vector res, void* arkode_mem);
int mriStep_NlsFPFunction(N_Vector yy, N_Vector res, void* arkode_mem);
int mriStep_NlsLSetup(sunbooleantype jbad, sunbooleantype* jcur,
                      void* arkode_mem);
int mriStep_NlsLSolve(N_Vector delta, void* arkode_mem);
int mriStep_NlsConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del,
                        sunrealtype tol, N_Vector ewt, void* arkode_mem);

/* Inner stepper functions */
int mriStepInnerStepper_HasRequiredOps(MRIStepInnerStepper stepper);
sunbooleantype mriStepInnerStepper_SupportsRTolAdaptivity(
  MRIStepInnerStepper stepper);
int mriStepInnerStepper_Evolve(MRIStepInnerStepper stepper, sunrealtype t0,
                               sunrealtype tout, N_Vector y);
int mriStepInnerStepper_EvolveSUNStepper(MRIStepInnerStepper stepper,
                                         sunrealtype t0, sunrealtype tout,
                                         N_Vector y);
int mriStepInnerStepper_FullRhs(MRIStepInnerStepper stepper, sunrealtype t,
                                N_Vector y, N_Vector f, int mode);
int mriStepInnerStepper_FullRhsSUNStepper(MRIStepInnerStepper stepper,
                                          sunrealtype t, N_Vector y, N_Vector f,
                                          int mode);
int mriStepInnerStepper_Reset(MRIStepInnerStepper stepper, sunrealtype tR,
                              N_Vector yR);
int mriStepInnerStepper_GetAccumulatedError(MRIStepInnerStepper stepper,
                                            sunrealtype* accum_error);
int mriStepInnerStepper_ResetAccumulatedError(MRIStepInnerStepper stepper);
int mriStepInnerStepper_SetRTol(MRIStepInnerStepper stepper, sunrealtype rtol);
int mriStepInnerStepper_ResetSUNStepper(MRIStepInnerStepper stepper,
                                        sunrealtype tR, N_Vector yR);
int mriStepInnerStepper_AllocVecs(MRIStepInnerStepper stepper, int count,
                                  N_Vector tmpl);
int mriStepInnerStepper_Resize(MRIStepInnerStepper stepper, ARKVecResizeFn resize,
                               void* resize_data, sunindextype lrw_diff,
                               sunindextype liw_diff, N_Vector tmpl);
int mriStepInnerStepper_FreeVecs(MRIStepInnerStepper stepper);
void mriStepInnerStepper_PrintMem(MRIStepInnerStepper stepper, FILE* outfile);

/* Compute forcing for inner stepper */
int mriStep_ComputeInnerForcing(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem,
                                int stage, sunrealtype t0, sunrealtype tf);

/* Return effective RK coefficients (nofast stage) */
int mriStep_RKCoeffs(MRIStepCoupling MRIC, int is, int* stage_map,
                     sunrealtype* Ae_row, sunrealtype* Ai_row);

/*===============================================================
  MRIStep SUNAdaptController wrapper module -- this is used to
  insert MRIStep in-between ARKODE at the "slow" time scale, and
  the inner-steppers that comprise the MRI time-step evolution.
  Since ARKODE itself only calls single-scale controller
  functions (e.g., EstimateStep and UpdateH), then this serves to
  translate those single-rate controller functions to the multi-
  rate context, leveraging MRIStep-specific knowledge of the
  slow+fast time scale relationship to CALL multi-rate controller
  functions (e.g., EstimateMRISteps, EstimateStepTol, UpdateMRIH,
  and UpdateMRIHTol) provided by the underlying multi-rate
  controller.
  ===============================================================*/

typedef struct _mriStepControlContent
{
  ARKodeMem ark_mem;         /* ARKODE memory pointer */
  ARKodeMRIStepMem step_mem; /* MRIStep memory pointer */
  SUNAdaptController C;      /* attached controller pointer */
}* mriStepControlContent;

#define MRICONTROL_C(C) (((mriStepControlContent)(C->content))->C)
#define MRICONTROL_A(C) (((mriStepControlContent)(C->content))->ark_mem)
#define MRICONTROL_S(C) (((mriStepControlContent)(C->content))->step_mem)

SUNAdaptController SUNAdaptController_MRIStep(ARKodeMem ark_mem,
                                              SUNAdaptController C);
SUNAdaptController_Type SUNAdaptController_GetType_MRIStep(SUNAdaptController C);
SUNErrCode SUNAdaptController_EstimateStep_MRIStep(SUNAdaptController C,
                                                   sunrealtype H, int P,
                                                   sunrealtype DSM,
                                                   sunrealtype* Hnew);
SUNErrCode SUNAdaptController_Reset_MRIStep(SUNAdaptController C);
SUNErrCode SUNAdaptController_Write_MRIStep(SUNAdaptController C, FILE* fptr);
SUNErrCode SUNAdaptController_UpdateH_MRIStep(SUNAdaptController C,
                                              sunrealtype h, sunrealtype dsm);
SUNErrCode SUNAdaptController_Space_MRIStep(SUNAdaptController C,
                                            long int* lenrw, long int* leniw);

/*===============================================================
  Reusable MRIStep Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
#define MSG_MRISTEP_NO_MEM      "Time step module memory is NULL."
#define MSG_NLS_INIT_FAIL       "The nonlinear solver's init routine failed."
#define MSG_MRISTEP_NO_COUPLING "The MRIStepCoupling is NULL."

#ifdef __cplusplus
}
#endif

#endif
