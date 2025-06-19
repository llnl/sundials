/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * Implementation header file for the main ARKODE integrator.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_IMPL_H
#define _ARKODE_IMPL_H

#include <stdarg.h>

#include <arkode/arkode.h>
#include <arkode/arkode_butcher.h>
#include <arkode/arkode_butcher_dirk.h>
#include <arkode/arkode_butcher_erk.h>
#include <arkode/arkode_mristep.h>

#include <sundials/priv/sundials_context_impl.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_adaptcontroller.h>
#include <sundials/sundials_adjointcheckpointscheme.h>
#include <sundials/sundials_adjointstepper.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_types.h>

#include "arkode_adapt_impl.h"
#include "arkode_relaxation_impl.h"
#include "arkode_root_impl.h"
#include "arkode_types_impl.h"
#include "sundials_logger_impl.h"
#include "sundials_macros.h"
#include "sundials_stepper_impl.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  SHORTCUTS
  ===============================================================*/

#define ARK_PROFILER ark_mem->sunctx->profiler
#define ARK_LOGGER   ark_mem->sunctx->logger

/*===============================================================
  MACROS
  ===============================================================*/

/* TODO(DJG): replace with signbit when C99+ is required */
#define DIFFERENT_SIGN(a, b) (((a) < 0 && (b) > 0) || ((a) > 0 && (b) < 0))
#define SAME_SIGN(a, b)      (((a) > 0 && (b) > 0) || ((a) < 0 && (b) < 0))

/*===============================================================
  ARKODE Private Constants
  ===============================================================*/

/* Basic ARKODE defaults */
/*   method order */
#define Q_DEFAULT 4
/*   max steps between returns */
#define MXSTEP_DEFAULT 500
/*   max number of error failures */
#define MAXNEF 7
/*   max number of convergence failures */
#define MAXNCF 10
/*   max number of constraint failures */
#define MAXCONSTRFAILS 10
/*   max number of t+h==h warnings */
#define MXHNIL 10

/* Numeric constants */
#define ZERO  SUN_RCONST(0.0)
#define TINY  SUN_RCONST(1.0e-10)
#define TENTH SUN_RCONST(0.1)
#define HALF  SUN_RCONST(0.5)
#define ONE   SUN_RCONST(1.0)
#define TWO   SUN_RCONST(2.0)
#define THREE SUN_RCONST(3.0)
#define FOUR  SUN_RCONST(4.0)
#define FIVE  SUN_RCONST(5.0)

/* Control constants for tolerances */
#define ARK_SS 0
#define ARK_SV 1
#define ARK_WF 2

/*===============================================================
  ARKODE Routine-Specific Constants
  ===============================================================*/

/*---------------------------------------------------------------
  Initialization types
  ---------------------------------------------------------------*/
#define FIRST_INIT  0 /* first step (re-)initialization */
#define RESET_INIT  1 /* reset initialization           */
#define RESIZE_INIT 2 /* resize initialization          */

/*---------------------------------------------------------------
  Control constants for lower-level time-stepping functions
  ---------------------------------------------------------------*/
#define PREDICT_AGAIN  +3
#define CONV_FAIL      +4
#define TRY_AGAIN      +5
#define FIRST_CALL     +6
#define PREV_CONV_FAIL +7
#define PREV_ERR_FAIL  +8
#define RHSFUNC_RECVR  +9
#define CONSTR_RECVR   +10
#define ARK_RETRY_STEP +11

/*---------------------------------------------------------------
  Return values for lower-level rootfinding functions
  ---------------------------------------------------------------*/
#define RTFOUND +1
#define CLOSERT +3

/*---------------------------------------------------------------
  Algorithmic constants
  ---------------------------------------------------------------
  ARKodeGetDky and arkStep:  FUZZ_FACTOR

  arkHin:  H0_LBFACTOR, H0_UBFACTOR, H0_BIAS and H0_ITERS

  time comparison factors:
     ONEPSM      safety factor for floating point comparisons
     ONEMSM      safety factor for floating point comparisons
  ---------------------------------------------------------------*/
#define FUZZ_FACTOR SUN_RCONST(100.0)

#define H0_LBFACTOR SUN_RCONST(100.0)
#define H0_UBFACTOR SUN_RCONST(0.1)
#define H0_BIAS     HALF
#define H0_ITERS    4

#define ONEPSM SUN_RCONST(1.000001)
#define ONEMSM SUN_RCONST(0.999999)

/*---------------------------------------------------------------
  Input flag to linear solver setup routine:  CONVFAIL

     ARK_NO_FAILURES : Either this is the first lsetup call for
                       this step, or the local error test failed
                       on the previous attempt at this step (but
                       the Newton iteration converged).

     ARK_FAIL_BAD_J  : This value is passed to lsetup if

                       (a) The previous Newton corrector
                           iteration did not converge and the
                           linear solver's setup routine
                           indicated that its Jacobian-related
                           data is not current
                       or

                       (b) During the previous Newton corrector
                           iteration, the linear solver's solve
                           routine failed in a recoverable manner
                           and the linear solver's setup routine
                           indicated that its Jacobian-related
                           data is not current.

     ARK_FAIL_OTHER  : During the current internal step try, the
                       previous Newton iteration failed to
                       converge even though the linear solver was
                       using current Jacobian-related data.
  --------------------------------------------------------------*/
#define ARK_NO_FAILURES 0
#define ARK_FAIL_BAD_J  1
#define ARK_FAIL_OTHER  2

/*===============================================================
  ARKODE Interface function definitions
  ===============================================================*/

/* NOTE: documentation for the purpose of these internal interface
   functions is located at the end of this file */

/* linear solver interface functions */
typedef int (*ARKLinsolInitFn)(ARKodeMem ark_mem);
typedef int (*ARKLinsolSetupFn)(ARKodeMem ark_mem, int convfail,
                                sunrealtype tpred, N_Vector ypred, N_Vector fpred,
                                sunbooleantype* jcurPtr, N_Vector vtemp1,
                                N_Vector vtemp2, N_Vector vtemp3);
typedef int (*ARKLinsolSolveFn)(ARKodeMem ark_mem, N_Vector b, sunrealtype tcur,
                                N_Vector ycur, N_Vector fcur,
                                sunrealtype client_tol, int mnewt);
typedef int (*ARKLinsolFreeFn)(ARKodeMem ark_mem);

/* mass matrix solver interface functions */
typedef int (*ARKMassInitFn)(ARKodeMem ark_mem);
typedef int (*ARKMassSetupFn)(ARKodeMem ark_mem, sunrealtype t, N_Vector vtemp1,
                              N_Vector vtemp2, N_Vector vtemp3);
typedef int (*ARKMassMultFn)(void* arkode_mem, N_Vector v, N_Vector Mv);
typedef int (*ARKMassSolveFn)(ARKodeMem ark_mem, N_Vector b,
                              sunrealtype client_tol);
typedef int (*ARKMassFreeFn)(ARKodeMem ark_mem);

/* time stepper interface functions -- general */
typedef int (*ARKTimestepInitFn)(ARKodeMem ark_mem, sunrealtype tout,
                                 int init_type);
typedef int (*ARKTimestepFullRHSFn)(ARKodeMem ark_mem, sunrealtype t,
                                    N_Vector y, N_Vector f, int mode);
typedef int (*ARKTimestepStepFn)(ARKodeMem ark_mem, sunrealtype* dsm, int* nflag);
typedef int (*ARKTimetepSetUserDataFn)(ARKodeMem ark_mem, void* user_data);
typedef int (*ARKTimestepPrintAllStats)(ARKodeMem ark_mem, FILE* outfile,
                                        SUNOutputFormat fmt);
typedef int (*ARKTimestepWriteParameters)(ARKodeMem ark_mem, FILE* fp);
typedef int (*ARKTimestepResize)(ARKodeMem ark_mem, N_Vector ynew,
                                 sunrealtype hscale, sunrealtype t0,
                                 ARKVecResizeFn resize, void* resize_data);
typedef int (*ARKTimestepReset)(ARKodeMem ark_mem, sunrealtype tR, N_Vector yR);
typedef void (*ARKTimestepFree)(ARKodeMem ark_mem);
typedef void (*ARKTimestepPrintMem)(ARKodeMem ark_mem, FILE* outfile);
typedef int (*ARKTimestepSetDefaults)(ARKodeMem ark_mem);
typedef int (*ARKTimestepSetOrder)(ARKodeMem ark_mem, int maxord);
typedef int (*ARKTimestepGetNumRhsEvals)(ARKodeMem ark_mem, int partition_index,
                                         long int* num_rhs_evals);
typedef int (*ARKTimestepSetStepDirection)(ARKodeMem ark_mem,
                                           sunrealtype stepdir);
typedef int (*ARKTimestepSetUseCompensatedSums)(ARKodeMem ark_mem,
                                                sunbooleantype onoff);

/* time stepper interface functions -- temporal adaptivity */
typedef int (*ARKTimestepGetEstLocalErrors)(ARKodeMem ark_mem, N_Vector ele);
typedef int (*ARKSetAdaptControllerFn)(ARKodeMem ark_mem, SUNAdaptController C);

/* time stepper interface functions -- relaxation */
typedef int (*ARKTimestepSetRelaxFn)(ARKodeMem ark_mem, ARKRelaxFn rfn,
                                     ARKRelaxJacFn rjac);

/* time stepper interface functions -- implicit solvers */
typedef int (*ARKTimestepAttachLinsolFn)(ARKodeMem ark_mem, ARKLinsolInitFn linit,
                                         ARKLinsolSetupFn lsetup,
                                         ARKLinsolSolveFn lsolve,
                                         ARKLinsolFreeFn lfree,
                                         SUNLinearSolver_Type lsolve_type,
                                         void* lmem);
typedef void (*ARKTimestepDisableLSetup)(ARKodeMem ark_mem);
typedef void* (*ARKTimestepGetLinMemFn)(ARKodeMem ark_mem);
typedef ARKRhsFn (*ARKTimestepGetImplicitRHSFn)(ARKodeMem ark_mem);
typedef int (*ARKTimestepGetGammasFn)(ARKodeMem ark_mem, sunrealtype* gamma,
                                      sunrealtype* gamrat, sunbooleantype** jcur,
                                      sunbooleantype* dgamma_fail);
typedef int (*ARKTimestepComputeState)(ARKodeMem ark_mem, N_Vector zcor,
                                       N_Vector z);
typedef int (*ARKTimestepSetNonlinearSolver)(ARKodeMem ark_mem,
                                             SUNNonlinearSolver NLS);
typedef int (*ARKTimestepSetLinear)(ARKodeMem ark_mem, int timedepend);
typedef int (*ARKTimestepSetNonlinear)(ARKodeMem ark_mem);
typedef int (*ARKTimestepSetAutonomous)(ARKodeMem ark_mem,
                                        sunbooleantype autonomous);
typedef int (*ARKTimestepSetNlsRhsFn)(ARKodeMem ark_mem, ARKRhsFn nls_fi);
typedef int (*ARKTimestepSetDeduceImplicitRhs)(ARKodeMem ark_mem,
                                               sunbooleantype deduce);
typedef int (*ARKTimestepSetNonlinCRDown)(ARKodeMem ark_mem, sunrealtype crdown);
typedef int (*ARKTimestepSetNonlinRDiv)(ARKodeMem ark_mem, sunrealtype rdiv);
typedef int (*ARKTimestepSetDeltaGammaMax)(ARKodeMem ark_mem, sunrealtype dgmax);
typedef int (*ARKTimestepSetLSetupFrequency)(ARKodeMem ark_mem, int msbp);
typedef int (*ARKTimestepSetPredictorMethod)(ARKodeMem ark_mem, int method);
typedef int (*ARKTimestepSetMaxNonlinIters)(ARKodeMem ark_mem, int maxcor);
typedef int (*ARKTimestepSetNonlinConvCoef)(ARKodeMem ark_mem,
                                            sunrealtype nlscoef);
typedef int (*ARKTimestepSetStagePredictFn)(ARKodeMem ark_mem,
                                            ARKStagePredictFn PredictStage);
typedef int (*ARKTimestepGetNumLinSolvSetups)(ARKodeMem ark_mem,
                                              long int* nlinsetups);
typedef int (*ARKTimestepGetCurrentGamma)(ARKodeMem ark_mem, sunrealtype* gamma);
typedef int (*ARKTimestepGetNonlinearSystemData)(
  ARKodeMem ark_mem, sunrealtype* tcur, N_Vector* zpred, N_Vector* z,
  N_Vector* Fi, sunrealtype* gamma, N_Vector* sdata, void** user_data);
typedef int (*ARKTimestepGetNumNonlinSolvIters)(ARKodeMem ark_mem,
                                                long int* nniters);
typedef int (*ARKTimestepGetNumNonlinSolvConvFails)(ARKodeMem ark_mem,
                                                    long int* nnfails);
typedef int (*ARKTimestepGetNonlinSolvStats)(ARKodeMem ark_mem, long int* nniters,
                                             long int* nnfails);

/* time stepper interface functions -- non-identity mass matrices */
typedef int (*ARKTimestepAttachMasssolFn)(
  ARKodeMem ark_mem, ARKMassInitFn minit, ARKMassSetupFn msetup,
  ARKMassMultFn mmult, ARKMassSolveFn msolve, ARKMassFreeFn mfree,
  sunbooleantype time_dep, SUNLinearSolver_Type msolve_type, void* mass_mem);
typedef void (*ARKTimestepDisableMSetup)(ARKodeMem ark_mem);
typedef void* (*ARKTimestepGetMassMemFn)(ARKodeMem ark_mem);

/* time stepper interface functions -- forcing */
typedef int (*ARKTimestepSetForcingFn)(ARKodeMem ark_mem, sunrealtype tshift,
                                       sunrealtype tscale, N_Vector* f,
                                       int nvecs);

/*===============================================================
  ARKODE interpolation module definition
  ===============================================================*/

/* Forward reference for pointer to ARKInterp_Ops object */
typedef struct _generic_ARKInterpOps* ARKInterpOps;

/* Forward reference for pointer to ARKInterp object */
typedef struct _generic_ARKInterp* ARKInterp;

/* Structure containing function pointers to interpolation operations  */
struct _generic_ARKInterpOps
{
  int (*resize)(ARKodeMem ark_mem, ARKInterp interp, ARKVecResizeFn resize,
                void* resize_data, sunindextype lrw_diff, sunindextype liw_diff,
                N_Vector tmpl);
  void (*free)(ARKodeMem ark_mem, ARKInterp interp);
  void (*print)(ARKInterp interp, FILE* outfile);
  int (*setdegree)(ARKodeMem ark_mem, ARKInterp interp, int degree);
  int (*init)(ARKodeMem ark_mem, ARKInterp interp, sunrealtype tnew);
  int (*update)(ARKodeMem ark_mem, ARKInterp interp, sunrealtype tnew);
  int (*evaluate)(ARKodeMem ark_mem, ARKInterp interp, sunrealtype tau, int d,
                  int order, N_Vector yout);
};

/* An interpolation module consists of an implementation-dependent 'content'
   structure, and a pointer to a structure of implementation-dependent operations. */
struct _generic_ARKInterp
{
  void* content;
  ARKInterpOps ops;
};

/* ARKInterp module functions */
int arkInterpResize(ARKodeMem ark_mem, ARKInterp interp, ARKVecResizeFn resize,
                    void* resize_data, sunindextype lrw_diff,
                    sunindextype liw_diff, N_Vector tmpl);
void arkInterpFree(ARKodeMem ark_mem, ARKInterp interp);
void arkInterpPrintMem(ARKInterp interp, FILE* outfile);
int arkInterpSetDegree(ARKodeMem ark_mem, ARKInterp interp, int degree);
int arkInterpInit(ARKodeMem ark_mem, ARKInterp interp, sunrealtype tnew);
int arkInterpUpdate(ARKodeMem ark_mem, ARKInterp interp, sunrealtype tnew);
int arkInterpEvaluate(ARKodeMem ark_mem, ARKInterp interp, sunrealtype tau,
                      int d, int order, N_Vector yout);

/*===============================================================
  ARKODE data structures
  ===============================================================*/

/*---------------------------------------------------------------
  Types : struct ARKodeMassMemRec, ARKodeMassMem
  ---------------------------------------------------------------
  The type ARKodeMassMem is type pointer to struct
  ARKodeMassMemRec.  This structure contains data pertaining to
  the use of a non-identity mass matrix.
  ---------------------------------------------------------------*/
typedef struct ARKodeMassMemRec
{
  /* mass matrix linear solver interface function pointers */
  ARKMassInitFn minit;
  ARKMassSetupFn msetup;
  ARKMassMultFn mmult;
  ARKMassSolveFn msolve;
  ARKMassFreeFn mfree;
  void* sol_mem;   /* mass matrix solver interface data */
  int msolve_type; /* mass matrix interface type:
                                  0=iterative; 1=direct; 2=custom */

}* ARKodeMassMem;

/*---------------------------------------------------------------
  Types : struct ARKodeMemRec, ARKodeMem
  ---------------------------------------------------------------
  The type ARKodeMem is type pointer to struct ARKodeMemRec.
  This structure contains fields to keep track of problem state.
  ---------------------------------------------------------------*/
struct ARKodeMemRec
{
  SUNContext sunctx;

  sunrealtype uround; /* machine unit roundoff */

  /* Problem specification data */
  void* user_data;               /* user ptr passed to supplied functions */
  int itol;                      /* itol = ARK_SS (scalar, default),
                                         ARK_SV (vector),
                                         ARK_WF (user weight function)  */
  int ritol;                     /* itol = ARK_SS (scalar, default),
                                         ARK_SV (vector),
                                         ARK_WF (user weight function)  */
  sunrealtype reltol;            /* relative tolerance                    */
  sunrealtype Sabstol;           /* scalar absolute solution tolerance    */
  N_Vector Vabstol;              /* vector absolute solution tolerance    */
  sunbooleantype atolmin0;       /* flag indicating that min(abstol) = 0  */
  sunrealtype SRabstol;          /* scalar absolute residual tolerance    */
  N_Vector VRabstol;             /* vector absolute residual tolerance    */
  sunbooleantype Ratolmin0;      /* flag indicating that min(Rabstol) = 0 */
  sunbooleantype user_efun;      /* SUNTRUE if user sets efun             */
  ARKEwtFn efun;                 /* function to set ewt                   */
  void* e_data;                  /* user pointer passed to efun           */
  sunbooleantype user_rfun;      /* SUNTRUE if user sets rfun             */
  ARKRwtFn rfun;                 /* function to set rwt                   */
  void* r_data;                  /* user pointer passed to rfun           */
  sunbooleantype constraintsSet; /* check inequality constraints          */

  /* Time stepper module -- general */
  void* step_mem;
  ARKTimestepInitFn step_init;
  ARKTimestepFullRHSFn step_fullrhs;
  ARKTimestepStepFn step;
  ARKTimetepSetUserDataFn step_setuserdata;
  ARKTimestepPrintAllStats step_printallstats;
  ARKTimestepWriteParameters step_writeparameters;
  ARKTimestepResize step_resize;
  ARKTimestepReset step_reset;
  ARKTimestepFree step_free;
  ARKTimestepPrintMem step_printmem;
  ARKTimestepSetDefaults step_setdefaults;
  ARKTimestepSetOrder step_setorder;
  ARKTimestepGetNumRhsEvals step_getnumrhsevals;
  ARKTimestepSetStepDirection step_setstepdirection;
  ARKTimestepSetUseCompensatedSums step_setusecompensatedsums;

  /* Time stepper module -- temporal adaptivity */
  sunbooleantype step_supports_adaptive;
  ARKSetAdaptControllerFn step_setadaptcontroller;
  ARKTimestepGetEstLocalErrors step_getestlocalerrors;

  /* Time stepper module -- relaxation */
  sunbooleantype step_supports_relaxation;
  ARKTimestepSetRelaxFn step_setrelaxfn;

  /* Time stepper module -- implicit solvers */
  sunbooleantype step_supports_implicit;
  ARKTimestepAttachLinsolFn step_attachlinsol;
  ARKTimestepDisableLSetup step_disablelsetup;
  ARKTimestepGetLinMemFn step_getlinmem;
  ARKTimestepGetImplicitRHSFn step_getimplicitrhs;
  ARKTimestepGetGammasFn step_getgammas;
  ARKTimestepComputeState step_computestate;
  ARKTimestepSetNonlinearSolver step_setnonlinearsolver;
  ARKTimestepSetLinear step_setlinear;
  ARKTimestepSetAutonomous step_setautonomous;
  ARKTimestepSetNonlinear step_setnonlinear;
  ARKTimestepSetNlsRhsFn step_setnlsrhsfn;
  ARKTimestepSetDeduceImplicitRhs step_setdeduceimplicitrhs;
  ARKTimestepSetNonlinCRDown step_setnonlincrdown;
  ARKTimestepSetNonlinRDiv step_setnonlinrdiv;
  ARKTimestepSetDeltaGammaMax step_setdeltagammamax;
  ARKTimestepSetLSetupFrequency step_setlsetupfrequency;
  ARKTimestepSetPredictorMethod step_setpredictormethod;
  ARKTimestepSetMaxNonlinIters step_setmaxnonliniters;
  ARKTimestepSetNonlinConvCoef step_setnonlinconvcoef;
  ARKTimestepSetStagePredictFn step_setstagepredictfn;
  ARKTimestepGetNumLinSolvSetups step_getnumlinsolvsetups;
  ARKTimestepGetCurrentGamma step_getcurrentgamma;
  ARKTimestepGetNonlinearSystemData step_getnonlinearsystemdata;
  ARKTimestepGetNumNonlinSolvIters step_getnumnonlinsolviters;
  ARKTimestepGetNumNonlinSolvConvFails step_getnumnonlinsolvconvfails;
  ARKTimestepGetNonlinSolvStats step_getnonlinsolvstats;

  /* Time stepper module -- non-identity mass matrices */
  sunbooleantype step_supports_massmatrix;
  ARKTimestepAttachMasssolFn step_attachmasssol;
  ARKTimestepDisableMSetup step_disablemsetup;
  ARKTimestepGetMassMemFn step_getmassmem;
  ARKMassMultFn step_mmult;

  /* Time stepper module -- forcing */
  ARKTimestepSetForcingFn step_setforcing;

  /* N_Vector storage */
  N_Vector ewt;                 /* error weight vector                        */
  N_Vector rwt;                 /* residual weight vector                     */
  sunbooleantype rwt_is_ewt;    /* SUNTRUE if rwt is a pointer to ewt         */
  N_Vector ycur;                /* pointer to user-provided solution memory;
                                   used as evolving solution by the time stepper
                                   modules */
  N_Vector yn;                  /* solution from the last successful step     */
  N_Vector fn;                  /* full IVP right-hand side from last step    */
  sunbooleantype fn_is_current; /* SUNTRUE if fn has been evaluated at yn     */
  N_Vector tempv1;              /* temporary storage vectors (for local use   */
  N_Vector tempv2;              /* and by time-stepping modules)              */
  N_Vector tempv3;
  N_Vector tempv4;
  N_Vector tempv5;

  N_Vector constraints; /* vector of inequality constraint options         */

  /* Temporal interpolation module */
  ARKInterp interp;
  int interp_type;
  int interp_degree;

  /* Tstop information */
  sunbooleantype tstopset;
  sunbooleantype tstopinterp;
  sunrealtype tstop;

  /* Time step data */
  sunrealtype hin;            /* initial step size                        */
  sunrealtype h;              /* current step size                        */
  sunrealtype hmin;           /* |h| >= hmin                              */
  sunrealtype hmax_inv;       /* |h| <= 1/hmax_inv                        */
  sunrealtype hprime;         /* next actual step size to be used         */
  sunrealtype next_h;         /* next dynamical step size (only used in
                                  getCurrentStep); note that this could
                                  overtake tstop */
  sunrealtype eta;            /* eta = hprime / h                         */
  sunrealtype tcur;           /* current internal value of t
                                  (changes with each stage)               */
  sunrealtype tretlast;       /* value of tret last returned by ARKODE    */
  sunbooleantype fixedstep;   /* flag to disable temporal adaptivity      */
  ARKodeHAdaptMem hadapt_mem; /* time step adaptivity structure           */

  /* Limits and various solver parameters */
  long int mxstep;    /* max number of internal steps for one user call */
  int mxhnil;         /* max number of warning messages issued to the
                              user that t+h == t for the next internal step  */
  int maxconstrfails; /* max number of constraint check failures        */
  int maxnef;         /* max error test fails in one step               */
  int maxncf;         /* max num alg. solver conv. fails in one step    */

  /* Counters */
  long int nst_attempts; /* number of attempted steps                  */
  long int nst;          /* number of internal steps taken             */
  int nhnil;             /* number of messages issued to the user that
                             t+h == t for the next iternal step         */
  long int ncfn;         /* num corrector convergence failures         */
  long int netf;         /* num error test failures                    */
  long int nconstrfails; /* number of constraint failures              */

  /* Space requirements for ARKODE */
  sunindextype lrw1; /* no. of sunrealtype words in 1 N_Vector          */
  sunindextype liw1; /* no. of integer words in 1 N_Vector           */
  long int lrw;      /* no. of sunrealtype words in ARKODE work vectors */
  long int liw;      /* no. of integer words in ARKODE work vectors  */

  /* Saved Values */
  sunrealtype h0u;   /* actual initial stepsize                     */
  sunrealtype tn;    /* time of last successful step                */
  sunrealtype terr;  /* error in tn for compensated sums            */
  sunrealtype hold;  /* last successful h value used                */
  sunrealtype tolsf; /* tolerance scale factor (suggestion to user) */
  ARKAccumError AccumErrorType; /* accumulated error estimation type   */
  sunrealtype AccumErrorStart;  /* time of last accumulated error reset */
  sunrealtype AccumError;       /* accumulated error estimate               */
  sunbooleantype VabstolMallocDone;
  sunbooleantype VRabstolMallocDone;
  sunbooleantype MallocDone;
  sunbooleantype initsetup;    /* denotes a call to InitialSetup is needed   */
  int init_type;               /* initialization type (see constants above)  */
  sunbooleantype firststage;   /* denotes first stage in simulation          */
  sunbooleantype initialized;  /* denotes arkInitialSetup has been done      */
  sunbooleantype call_fullrhs; /* denotes the full RHS fn will be called     */

  /* Rootfinding Data */
  ARKodeRootMem root_mem; /* root-finding structure */

  /* Relaxation Data */
  sunbooleantype relax_enabled; /* is relaxation enabled?    */
  ARKodeRelaxMem relax_mem;     /* relaxation data structure */

  /* User-supplied step solution post-processing function */
  ARKPostProcessFn ProcessStep;
  void* ps_data; /* pointer to user_data */

  /* User-supplied stage solution post-processing function */
  ARKPostProcessFn ProcessStage;

  sunbooleantype use_compensated_sums;

  /* Adjoint solver data */
  sunbooleantype load_checkpoint_fail;
  sunbooleantype do_adjoint;
  suncountertype adj_stage_idx; /* current stage index (only valid in adjoint context)*/
  suncountertype adj_step_idx; /* current step index (only valid in adjoint context)*/

  /* Checkpointing data */
  SUNAdjointCheckpointScheme checkpoint_scheme;
  suncountertype checkpoint_step_idx; /* the step number for checkpointing */

  /* XBraid interface variables */
  sunbooleantype force_pass; /* when true the step attempt loop will ignore the
                              return value (kflag) from arkCheckTemporalError
                              and set kflag = ARK_SUCCESS to force the step
                              attempt to always pass (if a solver failure did
                              not occur before the error test). */
  int last_kflag; /* last value of the return flag (kflag) from a call
                              to arkCheckTemporalError. This is only set when
                              force_pass is true and is used by the XBraid
                              interface to determine if a time step passed or
                              failed the time step error test.  */
};

/*===============================================================
  ARKODE PROTOTYPE FUNCTIONS (MAY BE REPLACED BY USER)
  ===============================================================*/

/* Prototype of internal rwtSet function */
int arkRwtSet(N_Vector ycur, N_Vector weight, void* data);

/* Prototype of internal explicit stability estimation function */
int arkExpStab(N_Vector y, sunrealtype t, sunrealtype* hstab, void* user_data);

/*===============================================================
  HIGH LEVEL ERROR HANDLER, USED THROUGHOUT ARKODE
  ===============================================================*/

void arkProcessError(ARKodeMem ark_mem, int error_code, int line,
                     const char* func, const char* file, const char* msgfmt, ...);

/*===============================================================
  ARKODE PRIVATE FUNCTION PROTOTYPES
  ===============================================================*/

ARKodeMem arkCreate(SUNContext sunctx);
int arkInit(ARKodeMem ark_mem, sunrealtype t0, N_Vector y0, int init_type);
sunbooleantype arkAllocVec(ARKodeMem ark_mem, N_Vector tmpl, N_Vector* v);
sunbooleantype arkAllocVecArray(int count, N_Vector tmpl, N_Vector** v,
                                sunindextype lrw1, long int* lrw,
                                sunindextype liw1, long int* liw);
sunbooleantype arkAllocVectors(ARKodeMem ark_mem, N_Vector tmpl);
sunbooleantype arkResizeVectors(ARKodeMem ark_mem, ARKVecResizeFn resize,
                                void* resize_data, sunindextype lrw_diff,
                                sunindextype liw_diff, N_Vector tmpl);
sunbooleantype arkResizeVec(ARKodeMem ark_mem, ARKVecResizeFn resize,
                            void* resize_data, sunindextype lrw_diff,
                            sunindextype liw_diff, N_Vector tmpl, N_Vector* v);
sunbooleantype arkResizeVecArray(ARKVecResizeFn resize, void* resize_data,
                                 int count, N_Vector tmpl, N_Vector** v,
                                 sunindextype lrw_diff, long int* lrw,
                                 sunindextype liw_diff, long int* liw);
void arkFreeVec(ARKodeMem ark_mem, N_Vector* v);
void arkFreeVecArray(int count, N_Vector** v, sunindextype lrw1, long int* lrw,
                     sunindextype liw1, long int* liw);
void arkFreeVectors(ARKodeMem ark_mem);
sunbooleantype arkCheckTimestepper(ARKodeMem ark_mem);
sunbooleantype arkCheckNvectorRequired(N_Vector tmpl);
sunbooleantype arkCheckNvectorOptional(ARKodeMem ark_mem);

int arkInitialSetup(ARKodeMem ark_mem, sunrealtype tout);
int arkStopTests(ARKodeMem ark_mem, sunrealtype tout, N_Vector yout,
                 sunrealtype* tret, int itask, int* ier);
int arkHin(ARKodeMem ark_mem, sunrealtype tout);
sunrealtype arkUpperBoundH0(ARKodeMem ark_mem, sunrealtype tdist);
int arkYddNorm(ARKodeMem ark_mem, sunrealtype hg, sunrealtype* yddnrm);

int arkCompleteStep(ARKodeMem ark_mem, sunrealtype dsm);
int arkHandleFailure(ARKodeMem ark_mem, int flag);

int arkEwtSetSS(N_Vector ycur, N_Vector weight, void* arkode_mem);
int arkEwtSetSV(N_Vector ycur, N_Vector weight, void* arkode_mem);
int arkEwtSetSmallReal(N_Vector ycur, N_Vector weight, void* arkode_mem);
int arkRwtSetSS(ARKodeMem ark_mem, N_Vector My, N_Vector weight);
int arkRwtSetSV(ARKodeMem ark_mem, N_Vector My, N_Vector weight);

int arkPredict_MaximumOrder(ARKodeMem ark_mem, sunrealtype tau, N_Vector yguess);
int arkPredict_VariableOrder(ARKodeMem ark_mem, sunrealtype tau, N_Vector yguess);
int arkPredict_CutoffOrder(ARKodeMem ark_mem, sunrealtype tau, N_Vector yguess);
int arkPredict_Bootstrap(ARKodeMem ark_mem, sunrealtype hj, sunrealtype tau,
                         int nvec, sunrealtype* cvals, N_Vector* Xvecs,
                         N_Vector yguess);
int arkCheckConvergence(ARKodeMem ark_mem, int* nflagPtr, int* ncfPtr);
int arkCheckConstraints(ARKodeMem ark_mem, int* constrfails, int* nflag);
int arkCheckTemporalError(ARKodeMem ark_mem, int* nflagPtr, int* nefPtr,
                          sunrealtype dsm);
int arkAccessHAdaptMem(void* arkode_mem, const char* fname, ARKodeMem* ark_mem,
                       ARKodeHAdaptMem* hadapt_mem);

int arkReplaceAdaptController(ARKodeMem ark_mem, SUNAdaptController C,
                              sunbooleantype take_ownership);
int arkSetAdaptivityMethod(void* arkode_mem, int imethod, int idefault, int pq,
                           sunrealtype adapt_params[3]);
int arkSetAdaptivityFn(void* arkode_mem, ARKAdaptFn hfun, void* h_data);

ARKODE_DIRKTableID arkButcherTableDIRKNameToID(const char* imethod);
ARKODE_ERKTableID arkButcherTableERKNameToID(const char* emethod);

/* utility functions for wrapping ARKODE as an MRIStep inner stepper */
int ark_MRIStepInnerEvolve(MRIStepInnerStepper stepper, sunrealtype t0,
                           sunrealtype tout, N_Vector y);
int ark_MRIStepInnerFullRhs(MRIStepInnerStepper stepper, sunrealtype t,
                            N_Vector y, N_Vector f, int mode);
int ark_MRIStepInnerReset(MRIStepInnerStepper stepper, sunrealtype tR,
                          N_Vector yR);
int ark_MRIStepInnerGetAccumulatedError(MRIStepInnerStepper stepper,
                                        sunrealtype* accum_error);
int ark_MRIStepInnerResetAccumulatedError(MRIStepInnerStepper stepper);
int ark_MRIStepInnerSetRTol(MRIStepInnerStepper stepper, sunrealtype rtol);

/* utility functions for wrapping ARKODE as a SUNStepper */
SUNErrCode arkSUNStepperSelfDestruct(SUNStepper stepper);

/* XBraid interface functions */
int arkSetForcePass(void* arkode_mem, sunbooleantype force_pass);
int arkGetLastKFlag(void* arkode_mem, int* last_kflag);

/*===============================================================
  Reusable ARKODE Error Messages
  ===============================================================*/

#define MSG_TIME   "t = " SUN_FORMAT_G
#define MSG_TIME_H "t = " SUN_FORMAT_G " and h = " SUN_FORMAT_G
#define MSG_TIME_INT                                                \
  "t = " SUN_FORMAT_G " is not between tcur - hold = " SUN_FORMAT_G \
  " and tcur = " SUN_FORMAT_G
#define MSG_TIME_TOUT  "tout = " SUN_FORMAT_G
#define MSG_TIME_TSTOP "tstop = " SUN_FORMAT_G

/* Initialization and I/O error messages */
#define MSG_ARK_NO_MEM         "arkode_mem = NULL illegal."
#define MSG_ARK_ARKMEM_FAIL    "Allocation of arkode_mem failed."
#define MSG_ARK_MEM_FAIL       "A memory request failed."
#define MSG_ARK_NO_MALLOC      "Attempt to call before ARKodeInit."
#define MSG_ARK_BAD_HMIN_HMAX  "Inconsistent step size limits: hmin > hmax."
#define MSG_ARK_BAD_RELTOL     "reltol < 0 illegal."
#define MSG_ARK_BAD_ABSTOL     "abstol has negative component(s) (illegal)."
#define MSG_ARK_NULL_ABSTOL    "abstol = NULL illegal."
#define MSG_ARK_BAD_RABSTOL    "rabstol has negative component(s) (illegal)."
#define MSG_ARK_NULL_RABSTOL   "rabstol = NULL illegal."
#define MSG_ARK_NULL_Y0        "y0 = NULL illegal."
#define MSG_ARK_Y0_FAIL_CONSTR "y0 fails to satisfy constraints."
#define MSG_ARK_NULL_F         "Must specify at least one of fe, fi (both NULL)."
#define MSG_ARK_NULL_G         "g = NULL illegal."
#define MSG_ARK_BAD_NVECTOR    "A required vector operation is not implemented."
#define MSG_ARK_BAD_CONSTR     "Illegal values in constraints vector."
#define MSG_ARK_NULL_DKY       "dky = NULL illegal."
#define MSG_ARK_BAD_T          "Illegal value for t. " MSG_TIME_INT
#define MSG_ARK_NO_ROOT        "Rootfinding was not initialized."

/* ARKODE Error Messages */
#define MSG_ARK_YOUT_NULL "yout = NULL illegal."
#define MSG_ARK_TRET_NULL "tret = NULL illegal."
#define MSG_ARK_BAD_EWT   "Initial ewt has component(s) equal to zero (illegal)."
#define MSG_ARK_EWT_NOW_BAD \
  "At " MSG_TIME ", a component of ewt has become <= 0."
#define MSG_ARK_BAD_RWT "Initial rwt has component(s) equal to zero (illegal)."
#define MSG_ARK_RWT_NOW_BAD \
  "At " MSG_TIME ", a component of rwt has become <= 0."
#define MSG_ARK_BAD_ITASK "Illegal value for itask."
#define MSG_ARK_BAD_H0    "h0 and tout - t0 inconsistent."
#define MSG_ARK_BAD_TOUT                    \
  "Trouble interpolating at " MSG_TIME_TOUT \
  ". tout too far back in direction of integration"
#define MSG_ARK_EWT_FAIL "The user-provide EwtSet function failed."
#define MSG_ARK_EWT_NOW_FAIL \
  "At " MSG_TIME ", the user-provide EwtSet function failed."
#define MSG_ARK_RWT_FAIL "The user-provide RwtSet function failed."
#define MSG_ARK_RWT_NOW_FAIL \
  "At " MSG_TIME ", the user-provide RwtSet function failed."
#define MSG_ARK_LINIT_FAIL "The linear solver's init routine failed."
#define MSG_ARK_HNIL_DONE                                                  \
  "The above warning has been issued mxhnil times and will not be issued " \
  "again for this problem."
#define MSG_ARK_TOO_CLOSE "tout too close to t0 to start integration."
#define MSG_ARK_MAX_STEPS \
  "At " MSG_TIME ", mxstep steps taken before reaching tout."
#define MSG_ARK_TOO_MUCH_ACC "At " MSG_TIME ", too much accuracy requested."
#define MSG_ARK_HNIL                                                       \
  "Internal " MSG_TIME_H " are such that t + h = t on the next step. The " \
  "solver will continue anyway."
#define MSG_ARK_ERR_FAILS \
  "At " MSG_TIME_H ", the error test failed repeatedly or with |h| = hmin."
#define MSG_ARK_CONV_FAILS \
  "At " MSG_TIME_H         \
  ", the solver convergence test failed repeatedly or with |h| = hmin."
#define MSG_ARK_SETUP_FAILED \
  "At " MSG_TIME ", the setup routine failed in an unrecoverable manner."
#define MSG_ARK_SOLVE_FAILED \
  "At " MSG_TIME ", the solve routine failed in an unrecoverable manner."
#define MSG_ARK_FAILED_CONSTR \
  "At " MSG_TIME ", unable to satisfy inequality constraints."
#define MSG_ARK_RHSFUNC_FAILED \
  "At " MSG_TIME               \
  ", the right-hand side routine failed in an unrecoverable manner."
#define MSG_ARK_RHSFUNC_UNREC                                                 \
  "At " MSG_TIME ", the right-hand side failed in a recoverable manner, but " \
  "no recovery is possible."
#define MSG_ARK_RHSFUNC_REPTD \
  "At " MSG_TIME " repeated recoverable right-hand side function errors."
#define MSG_ARK_RTFUNC_FAILED                                            \
  "At " MSG_TIME ", the rootfinding routine failed in an unrecoverable " \
  "manner."
#define MSG_ARK_CLOSE_ROOTS "Root found at and very near " MSG_TIME "."
#define MSG_ARK_BAD_TSTOP                                    \
  "The value " MSG_TIME_TSTOP " is behind current " MSG_TIME \
  " in the direction of integration."
#define MSG_ARK_INACTIVE_ROOTS                                         \
  "At the end of the first step, there are still some root functions " \
  "identically 0. This warning will not be issued again."
#define MSG_ARK_RESIZE_FAIL    "Error in user-supplied resize() function."
#define MSG_ARK_MASSINIT_FAIL  "The mass matrix solver's init routine failed."
#define MSG_ARK_MASSSETUP_FAIL "The mass matrix solver's setup routine failed."
#define MSG_ARK_MASSSOLVE_FAIL "The mass matrix solver failed."
#define MSG_ARK_NLS_FAIL \
  "At " MSG_TIME " the nonlinear solver failed in an unrecoverable manner."
#define MSG_ARK_USER_PREDICT_FAIL \
  "At " MSG_TIME                  \
  " the user-supplied predictor failed in an unrecoverable manner."
#define MSG_ARKADAPT_NO_MEM  "Adaptivity memory structure not allocated."
#define MSG_ARK_VECTOROP_ERR "At " MSG_TIME ", a vector operation failed."
#define MSG_ARK_INNERSTEP_FAILED \
  "At " MSG_TIME ", the inner stepper failed in an unrecoverable manner."
#define MSG_ARK_POSTPROCESS_STEP_FAIL \
  "At " MSG_TIME                      \
  ", the step postprocessing routine failed in an unrecoverable manner."
#define MSG_ARK_POSTPROCESS_STAGE_FAIL \
  "At " MSG_TIME                       \
  ", the stage postprocessing routine failed in an unrecoverable manner."
#define MSG_ARK_NULL_SUNCTX "sunctx = NULL illegal."
#define MSG_ARK_CONTEXT_MISMATCH \
  "Outer and inner steppers have different contexts."
#define MSG_ARK_MISSING_FULLRHS                                          \
  "Time-stepping module missing fullrhs routine (required by requested " \
  "solver configuration)."
#define MSG_ARK_INTERPOLATION_FAIL \
  "At " MSG_TIME ", interpolating the solution failed."
#define MSG_ARK_ADJOINT_BAD_VECTOR                                           \
  "JacPFn or JPvpFn was provided, but the number of subvectors in y is not " \
  "2. To perform ASA w.r.t. parameters, one subvector should be the state "  \
  "vector, and the other should be the parameter vector."

/*===============================================================

  Documentation for internal ARKODE interfaces

  ===============================================================

  Interfaces To Implicit Solvers

  ---------------------------------------------------------------

  ARKLinsolInitFn

  This function should complete initializations for a specific
  ARKODE linear solver interface, such as counters and statistics.
  This should return 0 if it has successfully initialized the
  ARKODE linear solver interface and a negative value otherwise.
  If an error does occur, an appropriate message should be sent
  to the error handler function.

  ---------------------------------------------------------------

  ARKLinsolSetupFn

  This function should prepare the linear solver interface for
  subsequent calls to the ARKLinsolSolveFn routine. It may
  recompute Jacobian-related data is it deems necessary. Its
  parameters are as follows:

  arkode_mem - void* problem memory pointer of type ARKodeMem. See
               the typedef earlier in this file.

  convfail - a flag to indicate any problem that occurred during
             the solution of the nonlinear equation on the
             current time step for which the linear solver is
             being used. This flag can be used to help decide
             whether the Jacobian data kept by a ARKODE linear
             solver needs to be updated or not.
             Its possible values have been documented above.

  tpred - the time for the current ARKODE internal step.

  ypred - the predicted y vector for the current ARKODE internal
          step.

  fpred - f(tpred, ypred).

  jcurPtr - a pointer to a boolean to be filled in by lsetup.
            The function should set *jcurPtr=SUNTRUE if its Jacobian
            data is current after the call and should set
            *jcurPtr=SUNFALSE if its Jacobian data is not current.
            Note: If lsetup calls for re-evaluation of
            Jacobian data (based on convfail and ARKODE state
            data), it should return *jcurPtr=SUNTRUE always;
            otherwise an infinite loop can result.

  vtemp1 - temporary N_Vector provided for use by lsetup.

  vtemp3 - temporary N_Vector provided for use by lsetup.

  vtemp3 - temporary N_Vector provided for use by lsetup.

  This routine should return 0 if successful, a positive value
  for a recoverable error, and a negative value for an
  unrecoverable error.

  ---------------------------------------------------------------

  ARKLinsolSolveFn

  This routine must solve the linear equation P x = b, where
  P is some approximation to (M - gamma J), M is the system mass
  matrix, J = (df/dy)(tcur,ycur), and the RHS vector b is input. The
  N-vector ycur contains the solver's current approximation to
  y(tcur) and the vector fcur contains the N_Vector f(tcur,ycur).
  The input client_tol contains the desired accuracy (in the wrms
  norm) of the routine calling the solver; the direct solvers
  ignore this value and iterative solvers tighten it by the
  factor eplifac.  The input mnewt is the current nonlinear
  iteration index (ignored by direct solvers, used by iterative
  solvers).

  Additional vectors that are set within the ARKODE memory
  structure, and that may be of use within an iterative linear
  solver, include:

  ewt - the error weight vector (scaling for solution vector)

  rwt - the residual weight vector (scaling for rhs vector)

  The solution is to be returned in the vector b.  This should
  return a positive value for a recoverable error and a
  negative value for an unrecoverable error. Success is
  indicated by a 0 return value.

  ---------------------------------------------------------------

  ARKLinsolFreeFn

  This should free up any memory allocated by the linear solver
  interface. This routine is called once a problem has been
  completed and the linear solver is no longer needed.  It should
  return 0 upon success, or a nonzero on failure.

  ===============================================================

  Interfaces For Non-Identity Mass Matrix Support

  ---------------------------------------------------------------

  ARKMassInitFn

  This function should complete initializations for a specific
  mass matrix linear solver interface, such as counters and
  statistics. A function of this type should return 0 if it
  has successfully initialized the mass matrix linear solver and
  a negative value otherwise.  If an error does occur, an
  appropriate message should be sent to the error handler function.

  ---------------------------------------------------------------

  ARKMassSetupFn

  This should prepare the mass matrix solver interface for
  subsequent calls to the ARKMassMultFn and ARKMassSolveFn
  routines. It may recompute mass matrix related data is it deems
  necessary. Its parameters are as follows:

  arkode_mem - void* problem memory pointer of type ARKodeMem. See
           the typedef earlier in this file.
  t - the 'time' at which to setup the mass matrix
  vtemp1, vtemp2, vtemp3 - temporary N_Vectors

  This routine should return 0 if successful, and a negative
  value for an unrecoverable error.

  ---------------------------------------------------------------

  ARKMassMultFn

  This must compute the matrix-vector product, z = M*v, where M is
  the system mass matrix the vector v is input, and the vector z
  is output. The mmult routine returns a positive value for a
  recoverable error and a negative value for an unrecoverable
  error. Success is indicated by a 0 return value.

  ---------------------------------------------------------------

  ARKMassSolveFn

  This must solve the linear equation M x = b, where M is the
  system mass matrix, and the RHS vector b is input. The
  sunrealtype client_tol contains the desired accuracy (in the wrms
  norm) of the routine calling the solver; the ARKDLS solver
  ignore this value and the ARKSPILS solver tightens it by the
  factor eplifac.  The solution is to be returned in the vector b.

  Additional vectors that are set within the ARKODE memory
  structure, and that may be of use within an iterative linear
  solver, include:

  ewt - the error weight vector (scaling for solution vector)

  rwt - the residual weight vector (scaling for rhs vector)

  This routine should return a positive value for a recoverable
  error and a negative value for an unrecoverable error. Success
  is indicated by a 0 return value.

  ---------------------------------------------------------------

  ARKMassFreeFn

  This should free up any memory allocated by the mass matrix
  solver interface. This routine is called once a problem has been
  completed and the solver is no longer needed.  It should return
  0 upon success, or a nonzero on failure.

  ===============================================================

  Internal Interface to Time Steppers -- General

  ---------------------------------------------------------------

  ARKTimestepInitFn

  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup. It should complete initializations for
  a specific ARKODE time stepping module, such as verifying
  compatibility of user-specified linear and nonlinear solver
  objects. The input init_type flag indicates if the call is
  for (re-)initializing, resizing, or resetting the problem.

  This routine should return 0 if it has successfully initialized
  the ARKODE time stepper module and a negative value otherwise.
  If an error does occur, an appropriate message should be sent
  to the error handler function.

  ---------------------------------------------------------------

  ARKTimestepFullRHSFn

  This routine must compute the full ODE right-hand side function
  at the inputs (t,y), and store the result in the N_Vector f.
  Depending on the type of stepper, this may be just the single
  ODE RHS function supplied (e.g. ERK, DIRK), or it may be
  the sum of many ODE RHS functions (e.g. ARK, MRI).  The 'mode'
  indicates where this routine is called:

     ARK_FULLRHS_START -> called in the following circumstances:
                          (a) at the beginning of a simulation
                              i.e., at (tn, yn) = (t0, y0) or
                              (tR, yR), or
                          (b) when transitioning between time
                              steps t_{n-1} \to t_{n} to fill
                              f_{n-1} within the Hermite
                              interpolation module.

                          In each case, the stepper may check the
                          fn_is_current flag to know whether
                          ARKODE believes that the RHS may have
                          already been computed at this (t,y)
                          value, in which case the stepper can
                          copy the RHS data from its own internal
                          storage instead of recomputing. If
                          these values are not current, then the
                          RHS should be recomputed, and the
                          stepper should consider storing the
                          values internally, for potential reuse
                          later.

     ARK_FULLRHS_END   -> called in the following circumstances:
                          (a) when temporal root-finding is
                              enabled, this will be called
                              in-between steps t_{n-1} \to t_{n}
                              to fill f_{n},
                          (b) when high-order dense output is
                              requested from the Hermite
                              interpolation module in-between
                              steps t_{n-1} \to t_{n} to fill
                              f_{n}, or
                          (c) when an implicit predictor is
                              requested from the Hermite
                              interpolation module within the
                              time step t_{n} \to t_{n+1}, in
                              which case f_{n} needs to be
                              filled.

                          Again, the stepper may check the
                          fn_is_current flag to know whether
                          ARKODE believes that the RHS may have
                          already been computed at this (t,y)
                          value, in which case the stepper can
                          copy the RHS data from its own
                          internal storage instead of
                          recomputing. If these values are not
                          current, then the RHS should be
                          recomputed, and the stepper should
                          consider storing the values internally,
                          for potential reuse later.

     ARK_FULLRHS_OTHER -> called in the following circumstances:
                          (a) when estimating the initial time
                              step size,
                          (b) for high-order dense output with the
                              Hermite interpolation module,
                          (c) by an "outer" stepper when ARKODE is
                              used as an inner solver), or
                          (d) when a high-order implicit predictor
                              is requested from the Hermite
                              interpolation module within the time
                              step t_{n} \to t_{n+1}.

                          While instances (a)-(c) will occur
                          in-between calls to the stepper's
                          ARKTimestepStepFn, instance (d) would
                          occur from within the ARKTimestepStepFn
                          (only if it calls an arkPredict_* function
                          that internally calls arkInterpEvaluate).
                          Since the (t,y) input does not correspond
                          to an "official" time step, the RHS
                          functions should always be evaluated, and
                          the values should *not* be stored anywhere
                          that will interfere with other reused
                          stepper data.

  This routine is only *required* to be supplied to ARKODE if:
  * ARKODE's initial time step selection algorithm (arkHin) is used,
  * the user requests temporal root-finding,
  * the Hermite interpolation module is used, or
  * the time-stepping module requests the "bootstrap" implicit predictor.
  Note that any stepper can itself require that this routine
  exist for its own internal business (as in both ERKStep and ARKStep),
  and/or call this routine for its own internal purposes.

  This routine should return 0 if successful, and a negative value
  otherwise.  If an error does occur, an appropriate message
  should be sent to the error handler function.

  ---------------------------------------------------------------

  ARKTimestepStepFn

  This routine serves the primary purpose of any ARKODE
  time-stepping module: it performs a single time step of the
  method (with embedding, if possible).

  It is assumed that this routine uses/modifies general problem
  data directly out of the main ARKodeMem structure, but that all
  method-specific data be stored in the step-module-specific data
  structure.  Relevant items in the ARKodeMem structure for this
  purpose include:
  - tcur -- the current "t" value
  - ycur -- the current "y" value on input; should hold the
    time-evolved solution on output
  - h -- the suggested/maximum "h" value to use; if the step
    eventually completes with a smaller "h" value, then that
    should be stored here
  - tn -- "t" value at end of the last successful step
  - nst -- the counter for overall successful steps
  - user_data -- the (void *) pointer returned to user for
    RHS calls
  - report / diagfp -- if any diagnostic information is
    to be saved to disk, the report flag indicates whether
    this is enabled, and diagfp provides the file pointer
    where this information should be written

  The output variable dsmPtr should contain estimate of the
  weighted local error if an embedding is present; otherwise it
  should be 0.

  The input/output variable nflagPtr is used to gauge convergence
  of any algebraic solvers within the step.  At the start of a new
  time step, this will initially have the value FIRST_CALL.  On
  return from this function, nflagPtr should have a value:
            0 => algebraic solve completed successfully
           >0 => solve did not converge at this step size
                 (but may with a smaller stepsize)
           <0 => solve encountered an unrecoverable failure

  The return value from this routine is:
            0 => step completed successfully
           >0 => step encountered recoverable failure;
                 reduce step and retry (if possible)
           <0 => step encountered unrecoverable failure

  ---------------------------------------------------------------

  ARKTimetepSetUserDataFn

  This optional routine provides the input from ARKodeSetUserData
  to the stepper.

  ---------------------------------------------------------------

  ARKTimestepPrintAllStats

  This optional routine allows the stepper to optionally print
  out any internally-stored solver statistics when
  ARKodePrintAllStats is called.

  ---------------------------------------------------------------

  ARKTimestepWriteParameters

  This optional routine allows the stepper to optionally print
  out any solver parameters when ARKodeWriteParameters is called.

  ---------------------------------------------------------------

  ARKTimestepResize

  This optional routine allows the stepper to optionally resize
  any internal vector storage when ARKodeResize is called.

  ---------------------------------------------------------------

  ARKTimestepReset

  This optional routine allows the stepper to reset any internal
  data when ARKodeReset is called.

  ---------------------------------------------------------------

  ARKTimestepFree

  This optional routine allows the stepper the free any
  interally-stored memory when ARKodeFree is called.

  ---------------------------------------------------------------

  ARKTimestepPrintMem

  This optional routine allows the stepper to output any internal
  memory (typically for debugging purposes) when ARKodePrintMem is
  called.

  ---------------------------------------------------------------

  ARKTimestepSetDefaults

  This optional routine allows the stepper to reset any internal
  solver parameters to their default values, and is called by
  ARKodeSetDefaults.

  ---------------------------------------------------------------

  ARKTimestepSetOrder

  This optional routine allows the stepper to accept any user-
  requested method order parameter that was passed to
  ARKodeSetOrder.

  ===============================================================

  Internal Interface to Time Steppers -- Temporal Adaptivity

  These should only be provided if the stepper supports temporal
  adaptivity, and should be indicated by setting the flag
  "step_supports_adaptive" to SUNTRUE.

  ---------------------------------------------------------------

  ARKTimestepGetEstLocalErrors

  This routine requests the stepper to copy its internal
  estimate of the local truncation error to the output (called by
  ARKodeGetEstLocalErrors).

  ===============================================================

  Internal Interface to Time Steppers -- Relaxation

  These should only be provided if the stepper supports
  "relaxation Runge--Kutta methods" (or similar), and should be
  indicated by setting the flag "step_supports_relaxation" to SUNTRUE.

  ---------------------------------------------------------------

  ARKTimestepSetRelaxFn

  This routine is called by ARKodeSetRelaxFn, and expects the
  stepper to call ARKODE's "arkRelaxCreate" routine with the
  appropriate stepper-specific function pointers.

  ===============================================================

  Internal Interface to Time Steppers -- Implicit Solvers

  These should only be provided if the stepper uses implicit
  linear/nonlinear solvers, and should be indicated by setting
  the flag "step_supports_implicit" to SUNTRUE.

  ---------------------------------------------------------------

  ARKTimestepAttachLinsolFn

  This routine should attach the various set of system linear
  solver interface routines, linear solver interface data
  structure, and system linear solver type to the ARKODE time
  stepping module pointed to in ark_mem->step_mem.  This will
  be called by the ARKODE linear solver interface.

  This routine should return 0 if it has successfully attached
  these items and a negative value otherwise.  If an error does
  occur, an appropriate message should be sent to the ARKODE
  error handler function.

  ---------------------------------------------------------------

  ARKTimestepDisableLSetup

  This routine should NULLify any ARKLinsolSetupFn function
  pointer stored in the ARKODE time stepping module (initially set
  in a call to ARKTimestepAttachLinsolFn).

  This routine has no return value.

  ---------------------------------------------------------------

  ARKTimestepGetLinMemFn

  This routine should return the linear solver memory structure
  used by the ARKODE time stepping module pointed to in
  ark_mem->step_mem.  This will be called by the ARKODE linear
  solver interface.

  This routine should return NULL if no linear solver memory
  structure is attached.

  ---------------------------------------------------------------

  ARKTimestepGetImplicitRHSFn

  This routine should return the implicit RHS function pointer for
  the current nonlinear solve (if there are multiple); it is used
  inside the linear solver interfaces for approximation of
  Jacobian matrix elements and/or matrix-vector products.

  This routine should return NULL if no implicit RHS function is
  active.

  ---------------------------------------------------------------

  ARKTimestepGetGammasFn

  This routine should fill the current value of gamma, the ratio
  of the current gamma value to the gamma value when the
  Jacobian/preconditioner was last updated, a pointer to the
  time step module internal sunbooleantype variable indicating
  whether the preconditioner is current, and a logic value
  indicating whether the gamma value is sufficiently stale
  to cause recomputation of Jacobian/preconditioner data.  Here,
  gamma is the coefficient preceding the RHS Jacobian
  matrix, J, in the full nonlinear system Jacobian,
  A = M - gamma*J.

  The time step module must contain a sunbooleantype variable to
  provide for the boolentype pointer (jcur).  This is only used
  by iterative linear solvers, so could be NULL for time step
  modules that only work with direct linear solvers.  Optionally,
  the value of this parameter could be set to SUNFALSE prior to
  return from the ARKTimestepGetGammasFn to force recalculation
  of preconditioner information.

  The value of the logic flag is used as follows:  if a previous
  Newton iteration failed due to a bad Jacobian/preconditioner,
  and this flag is SUNFALSE, this will trigger recalculation of
  the Jacobian/preconditioner.

  This routine should return 0 if it has successfully attached
  these items, and a negative value otherwise.  If an error does
  occur, an appropriate message should be sent to the ARKODE
  error handler function.

  ---------------------------------------------------------------

  ARKTimestepComputeState

  This routine should combine any stepper-stored prediction with
  the input correction to fill the current state vector within
  an implicit solve (called by ARKodeComputeState).

  ---------------------------------------------------------------

  ARKTimestepSetNonlinearSolver

  This routine is called by ARKodeSetNonlinearSolver, and allows
  the stepper to store the corresponding user input.

  ---------------------------------------------------------------

  ARKTimestepSetLinear

  This routine is called by ARKodeSetLinear, and allows the
  stepper to store the corresponding user input.

  ---------------------------------------------------------------

  ARKTimestepSetAutonomous

  This routine is called by ARKodeSetAutonomous, and allows the
  stepper to store the corresponding user input.

  ---------------------------------------------------------------

  ARKTimestepSetNonlinear

  This routine is called by ARKodeSetNonlinear, and allows the
  stepper to store the corresponding user input.

  ---------------------------------------------------------------

  ARKTimestepSetNlsRhsFn

  This routine is called by ARKodeSetNlsRhsFn, and allows the
  stepper to store the corresponding user input.

  ---------------------------------------------------------------

  ARKTimestepSetDeduceImplicitRhs

  This routine is called by ARKodeSetDeduceImplicitRhs, and
  allows the stepper to store the corresponding user input.

  ---------------------------------------------------------------

  ARKTimestepSetNonlinCRDown

  This routine is called by ARKodeSetNonlinCRDown, and allows
  the stepper to store the corresponding user input.

  ---------------------------------------------------------------

  ARKTimestepSetNonlinRDiv

  This routine is called by ARKodeSetNonlinRDiv, and allows
  the stepper to store the corresponding user input.

  ---------------------------------------------------------------

  ARKTimestepSetDeltaGammaMax

  This routine is called by ARKodeSetDeltaGammaMax, and allows
  the stepper to store the corresponding user input.

  ---------------------------------------------------------------

  ARKTimestepSetLSetupFrequency

  This routine is called by ARKodeSetLSetupFrequency, and allows
  the stepper to store the corresponding user input.

  ---------------------------------------------------------------

  ARKTimestepSetPredictorMethod

  This routine is called by ARKodeSetPredictorMethod, and allows
  the stepper to store the corresponding user input.

  ---------------------------------------------------------------

  ARKTimestepSetMaxNonlinIters

  This routine is called by ARKodeSetMaxNonlinIters, and allows
  the stepper to store the corresponding user input.

  ---------------------------------------------------------------

  ARKTimestepSetNonlinConvCoef

  This routine is called by ARKodeSetNonlinConvCoef, and allows
  the stepper to store the corresponding user input.

  ---------------------------------------------------------------

  ARKTimestepSetStagePredictFn

  This routine is called by ARKodeSetStagePredictFn, and allows
  the stepper to store the corresponding user input.

  ---------------------------------------------------------------

  ARKTimestepGetNumLinSolvSetups

  This routine is called by ARKodeGetNumLinSolvSetups, and
  requests that the stepper return the corresponding output
  value.

  ---------------------------------------------------------------

  ARKTimestepGetCurrentGamma

  This routine is called by ARKodeGetCurrentGamma, and
  requests that the stepper return the corresponding output
  value.

  ---------------------------------------------------------------

  ARKTimestepGetNonlinearSystemData

  This routine is called by ARKodeGetNonlinearSystemData, and
  requests that the stepper return the corresponding output
  values.

  ---------------------------------------------------------------

  ARKTimestepGetNumNonlinSolvIters

  This routine is called by ARKodeGetNumNonlinSolvIters, and
  requests that the stepper return the corresponding output
  value.

  ---------------------------------------------------------------

  ARKTimestepGetNumNonlinSolvConvFails

  This routine is called by ARKodeGetNumNonlinSolvConvFails, and
  requests that the stepper return the corresponding output
  value.

  ---------------------------------------------------------------

  ARKTimestepGetNonlinSolvStats

  This routine is called by ARKodeGetNonlinSolvStats, and
  requests that the stepper return the corresponding output
  values.

  ===============================================================

  Internal Interface to Time Steppers -- Non-identity Mass
  Matrices

  These should only be provided if the stepper supports problems
  with non-identity mass matrices, and should be indicated by
  setting the flag "step_supports_massmatrix" to SUNTRUE.

  ---------------------------------------------------------------

  ARKTimestepAttachMasssolFn

  This routine should attach the various set of mass matrix
  linear solver interface routines, data structure, mass matrix
  type, and solver type to the ARKODE time stepping module
  pointed to in ark_mem->step_mem.  This will be called by the
  ARKODE linear solver interface.

  This routine should return 0 if it has successfully attached
  these items, and a negative value otherwise.  If an error does
  occur, an appropriate message should be sent to the ARKODE
  error handler function.

  ---------------------------------------------------------------

  ARKTimestepDisableMSetup

  This routine should NULLify any ARKMassSetupFn function pointer
  stored in the ARKODE time stepping module (initially set in a
  call to ARKTimestepAttachMasssolFn).

  This routine has no return value.

  ---------------------------------------------------------------

  ARKTimestepGetMassMemFn

  This routine should return the mass matrix linear solver memory
  structure used by the ARKODE time stepping module pointed to in
  ark_mem->step_mem.  This will be called the ARKODE mass matrix
  solver interface.

  This routine should return NULL if no mass matrix solver memory
  structure is attached.

  ===============================================================*/

#ifdef __cplusplus
}
#endif

#endif
