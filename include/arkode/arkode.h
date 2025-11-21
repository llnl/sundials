/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * This is the header file for the main ARKODE infrastructure.
 * -----------------------------------------------------------------
 * ARKODE is used to numerically solve the ordinary initial value
 * problems using one-step methods.  Users do not call ARKODE
 * infrastructure routines directly; they instead interact with
 * one of the time stepping modules built on top of ARKODE.
 * These time step modules define their supported problem types,
 * solver options, etc.
 *
 * This file serves to define constants and provide function
 * prototypes for use across ARKODE-based time integration
 * modules.
 * -----------------------------------------------------------------*/

#ifndef _ARKODE_H
#define _ARKODE_H

#include <arkode/arkode_butcher.h>
#include <stdio.h>
#include <sundials/sundials_adjointstepper.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_stepper.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------
 * ARKODE Constants
 * ----------------- */

/* usage modes (itask) */
#define ARK_NORMAL   1
#define ARK_ONE_STEP 2

/* adaptivity module flags */
#define ARK_ADAPT_CUSTOM   -1
#define ARK_ADAPT_PID      0
#define ARK_ADAPT_PI       1
#define ARK_ADAPT_I        2
#define ARK_ADAPT_EXP_GUS  3
#define ARK_ADAPT_IMP_GUS  4
#define ARK_ADAPT_IMEX_GUS 5

/* Constants for evaluating the full RHS */
#define ARK_FULLRHS_START 0
#define ARK_FULLRHS_END   1
#define ARK_FULLRHS_OTHER 2

/* interpolation module flags */

/*    max allowed degree */
#define ARK_INTERP_MAX_DEGREE 5

/*    interpolation module types */
#define ARK_INTERP_NONE     -1
#define ARK_INTERP_HERMITE  0
#define ARK_INTERP_LAGRANGE 1

/* return values */

#define ARK_SUCCESS      0
#define ARK_TSTOP_RETURN 1
#define ARK_ROOT_RETURN  2

#define ARK_WARNING 99

#define ARK_TOO_MUCH_WORK -1
#define ARK_TOO_MUCH_ACC  -2
#define ARK_ERR_FAILURE   -3
#define ARK_CONV_FAILURE  -4

#define ARK_LINIT_FAIL        -5
#define ARK_LSETUP_FAIL       -6
#define ARK_LSOLVE_FAIL       -7
#define ARK_RHSFUNC_FAIL      -8
#define ARK_FIRST_RHSFUNC_ERR -9
#define ARK_REPTD_RHSFUNC_ERR -10
#define ARK_UNREC_RHSFUNC_ERR -11
#define ARK_RTFUNC_FAIL       -12
#define ARK_LFREE_FAIL        -13
#define ARK_MASSINIT_FAIL     -14
#define ARK_MASSSETUP_FAIL    -15
#define ARK_MASSSOLVE_FAIL    -16
#define ARK_MASSFREE_FAIL     -17
#define ARK_MASSMULT_FAIL     -18

#define ARK_CONSTR_FAIL -19
#define ARK_MEM_FAIL    -20
#define ARK_MEM_NULL    -21
#define ARK_ILL_INPUT   -22
#define ARK_NO_MALLOC   -23
#define ARK_BAD_K       -24
#define ARK_BAD_T       -25
#define ARK_BAD_DKY     -26
#define ARK_TOO_CLOSE   -27

#define ARK_VECTOROP_ERR -28

#define ARK_NLS_INIT_FAIL   -29
#define ARK_NLS_SETUP_FAIL  -30
#define ARK_NLS_SETUP_RECVR -31
#define ARK_NLS_OP_ERR      -32

#define ARK_INNERSTEP_ATTACH_ERR -33
#define ARK_INNERSTEP_FAIL       -34
#define ARK_OUTERTOINNER_FAIL    -35
#define ARK_INNERTOOUTER_FAIL    -36

/* ARK_POSTPROCESS_FAIL equals ARK_POSTPROCESS_STEP_FAIL
   for backwards compatibility */
#define ARK_POSTPROCESS_FAIL       -37
#define ARK_POSTPROCESS_STEP_FAIL  -37
#define ARK_POSTPROCESS_STAGE_FAIL -38

#define ARK_USER_PREDICT_FAIL -39
#define ARK_INTERP_FAIL       -40

#define ARK_INVALID_TABLE -41

#define ARK_CONTEXT_ERR -42

#define ARK_RELAX_FAIL      -43
#define ARK_RELAX_MEM_NULL  -44
#define ARK_RELAX_FUNC_FAIL -45
#define ARK_RELAX_JAC_FAIL  -46

#define ARK_CONTROLLER_ERR -47

#define ARK_STEPPER_UNSUPPORTED -48

#define ARK_DOMEIG_FAIL          -49
#define ARK_MAX_STAGE_LIMIT_FAIL -50

#define ARK_SUNSTEPPER_ERR     -51
#define ARK_STEP_DIRECTION_ERR -52

#define ARK_ADJ_CHECKPOINT_FAIL -53
#define ARK_ADJ_RECOMPUTE_FAIL  -54
#define ARK_SUNADJSTEPPER_ERR   -55

#define ARK_DEE_FAIL -56

#define ARK_UNRECOGNIZED_ERROR -99

/* ------------------------------
 * User-Supplied Function Types
 * ------------------------------ */

typedef int (*ARKRhsFn)(sunrealtype t, N_Vector y, N_Vector ydot,
                        void* user_data);

typedef int (*ARKRootFn)(sunrealtype t, N_Vector y, sunrealtype* gout_1d,
                         void* user_data);

typedef int (*ARKEwtFn)(N_Vector y, N_Vector ewt, void* user_data);

typedef int (*ARKRwtFn)(N_Vector y, N_Vector rwt, void* user_data);

typedef int (*ARKAdaptFn)(N_Vector y, sunrealtype t, sunrealtype h1,
                          sunrealtype h2, sunrealtype h3, sunrealtype e1,
                          sunrealtype e2, sunrealtype e3, int q, int p,
                          sunrealtype* hnew, void* user_data);

typedef int (*ARKExpStabFn)(N_Vector y, sunrealtype t, sunrealtype* hstab,
                            void* user_data);

typedef int (*ARKVecResizeFn)(N_Vector y, N_Vector ytemplate, void* user_data);

typedef int (*ARKPostProcessFn)(sunrealtype t, N_Vector y, void* user_data);

typedef int (*ARKStagePredictFn)(sunrealtype t, N_Vector zpred, void* user_data);

typedef int (*ARKRelaxFn)(N_Vector y, sunrealtype* r, void* user_data);

typedef int (*ARKRelaxJacFn)(N_Vector y, N_Vector J, void* user_data);

/* ------------------------------------------------
 * MRIStep Inner Stepper Type (forward declaration)
 * ------------------------------------------------ */

typedef _SUNDIALS_STRUCT_ _MRIStepInnerStepper* MRIStepInnerStepper;

/* --------------------------
 * Relaxation Solver Options
 * -------------------------- */
enum ARKRelaxSolver
{
  ARK_RELAX_BRENT,
  ARK_RELAX_NEWTON
};

#ifndef SWIG
typedef enum ARKRelaxSolver ARKRelaxSolver;
#endif

/* --------------------------
 * Error Accumulation Options
 * -------------------------- */

enum ARKAccumError
{
  ARK_ACCUMERROR_NONE,
  ARK_ACCUMERROR_MAX,
  ARK_ACCUMERROR_SUM,
  ARK_ACCUMERROR_AVG
};

#ifndef SWIG
typedef enum ARKAccumError ARKAccumError;
#endif

/* --------------------------
 * Shared API routines
 * -------------------------- */

/* Command-line control over ARKODE options */
SUNDIALS_EXPORT int ARKodeSetOptions(void* arkode_mem, const char* arkid,
                                     const char* file_name, int argc,
                                     char* argv[]);

/* Resize and Reset functions */
SUNDIALS_EXPORT int ARKodeResize(void* arkode_mem, N_Vector ynew,
                                 sunrealtype hscale, sunrealtype t0,
                                 ARKVecResizeFn resize, void* resize_data);
SUNDIALS_EXPORT int ARKodeReset(void* arkode_mem, sunrealtype tR, N_Vector yR);

/* Utility to wrap ARKODE as an MRIStepInnerStepper */
SUNDIALS_EXPORT int ARKodeCreateMRIStepInnerStepper(void* arkode_mem,
                                                    MRIStepInnerStepper* stepper);

/* Tolerance input functions */
SUNDIALS_EXPORT int ARKodeSStolerances(void* arkode_mem, sunrealtype reltol,
                                       sunrealtype abstol);
SUNDIALS_EXPORT int ARKodeSVtolerances(void* arkode_mem, sunrealtype reltol,
                                       N_Vector abstol);
SUNDIALS_EXPORT int ARKodeWFtolerances(void* arkode_mem, ARKEwtFn efun);

/* Residual tolerance input functions */
SUNDIALS_EXPORT int ARKodeResStolerance(void* arkode_mem, sunrealtype rabstol);
SUNDIALS_EXPORT int ARKodeResVtolerance(void* arkode_mem, N_Vector rabstol);
SUNDIALS_EXPORT int ARKodeResFtolerance(void* arkode_mem, ARKRwtFn rfun);

/* Rootfinding */
SUNDIALS_EXPORT int ARKodeRootInit(void* arkode_mem, int nrtfn, ARKRootFn g);
SUNDIALS_EXPORT int ARKodeSetRootDirection(void* arkode_mem, int* rootdir);
SUNDIALS_EXPORT int ARKodeSetNoInactiveRootWarn(void* arkode_mem);

/* Optional input functions (general) */
SUNDIALS_EXPORT int ARKodeSetDefaults(void* arkode_mem);
SUNDIALS_EXPORT int ARKodeSetOrder(void* arkode_mem, int maxord);
SUNDIALS_EXPORT int ARKodeSetInterpolantType(void* arkode_mem, int itype);
SUNDIALS_EXPORT int ARKodeSetInterpolantDegree(void* arkode_mem, int degree);
SUNDIALS_EXPORT int ARKodeSetMaxNumSteps(void* arkode_mem, long int mxsteps);
SUNDIALS_EXPORT int ARKodeSetInterpolateStopTime(void* arkode_mem,
                                                 sunbooleantype interp);
SUNDIALS_EXPORT int ARKodeSetStopTime(void* arkode_mem, sunrealtype tstop);
SUNDIALS_EXPORT int ARKodeClearStopTime(void* arkode_mem);
SUNDIALS_EXPORT int ARKodeSetFixedStep(void* arkode_mem, sunrealtype hfixed);
SUNDIALS_EXPORT int ARKodeSetStepDirection(void* arkode_mem, sunrealtype stepdir);
SUNDIALS_EXPORT int ARKodeSetUserData(void* arkode_mem, void* user_data);
SUNDIALS_EXPORT int ARKodeSetPostprocessStepFn(void* arkode_mem,
                                               ARKPostProcessFn ProcessStep);
SUNDIALS_EXPORT int ARKodeSetPostprocessStageFn(void* arkode_mem,
                                                ARKPostProcessFn ProcessStage);

/* Optional input functions (implicit solver) */
SUNDIALS_EXPORT int ARKodeSetNonlinearSolver(void* arkode_mem,
                                             SUNNonlinearSolver NLS);
SUNDIALS_EXPORT int ARKodeSetLinear(void* arkode_mem, int timedepend);
SUNDIALS_EXPORT int ARKodeSetNonlinear(void* arkode_mem);
SUNDIALS_EXPORT int ARKodeSetAutonomous(void* arkode_mem,
                                        sunbooleantype autonomous);
SUNDIALS_EXPORT int ARKodeSetNlsRhsFn(void* arkode_mem, ARKRhsFn nls_fi);
SUNDIALS_EXPORT int ARKodeSetDeduceImplicitRhs(void* arkode_mem,
                                               sunbooleantype deduce);
SUNDIALS_EXPORT int ARKodeSetNonlinCRDown(void* arkode_mem, sunrealtype crdown);
SUNDIALS_EXPORT int ARKodeSetNonlinRDiv(void* arkode_mem, sunrealtype rdiv);
SUNDIALS_EXPORT int ARKodeSetDeltaGammaMax(void* arkode_mem, sunrealtype dgmax);
SUNDIALS_EXPORT int ARKodeSetLSetupFrequency(void* arkode_mem, int msbp);
SUNDIALS_EXPORT int ARKodeSetPredictorMethod(void* arkode_mem, int method);
SUNDIALS_EXPORT int ARKodeSetMaxNonlinIters(void* arkode_mem, int maxcor);
SUNDIALS_EXPORT int ARKodeSetMaxConvFails(void* arkode_mem, int maxncf);
SUNDIALS_EXPORT int ARKodeSetNonlinConvCoef(void* arkode_mem,
                                            sunrealtype nlscoef);
SUNDIALS_EXPORT int ARKodeSetStagePredictFn(void* arkode_mem,
                                            ARKStagePredictFn PredictStage);

/* Optional input functions (temporal adaptivity) */
SUNDIALS_EXPORT int ARKodeSetAdaptController(void* arkode_mem,
                                             SUNAdaptController C);
SUNDIALS_EXPORT int ARKodeSetAdaptControllerByName(void* arkode_mem,
                                                   const char* cname);
SUNDIALS_EXPORT int ARKodeSetAdaptivityAdjustment(void* arkode_mem, int adjust);
SUNDIALS_EXPORT int ARKodeSetCFLFraction(void* arkode_mem, sunrealtype cfl_frac);
SUNDIALS_EXPORT int ARKodeSetErrorBias(void* arkode_mem, sunrealtype bias);
SUNDIALS_EXPORT int ARKodeSetSafetyFactor(void* arkode_mem, sunrealtype safety);
SUNDIALS_EXPORT int ARKodeSetMaxGrowth(void* arkode_mem, sunrealtype mx_growth);
SUNDIALS_EXPORT int ARKodeSetMinReduction(void* arkode_mem, sunrealtype eta_min);
SUNDIALS_EXPORT int ARKodeSetFixedStepBounds(void* arkode_mem, sunrealtype lb,
                                             sunrealtype ub);
SUNDIALS_EXPORT int ARKodeSetMaxFirstGrowth(void* arkode_mem, sunrealtype etamx1);
SUNDIALS_EXPORT int ARKodeSetMaxEFailGrowth(void* arkode_mem, sunrealtype etamxf);
SUNDIALS_EXPORT int ARKodeSetSmallNumEFails(void* arkode_mem, int small_nef);
SUNDIALS_EXPORT int ARKodeSetMaxCFailGrowth(void* arkode_mem, sunrealtype etacf);
SUNDIALS_EXPORT int ARKodeSetStabilityFn(void* arkode_mem, ARKExpStabFn EStab,
                                         void* estab_data);
SUNDIALS_EXPORT int ARKodeSetMaxErrTestFails(void* arkode_mem, int maxnef);
SUNDIALS_EXPORT int ARKodeSetConstraints(void* arkode_mem, N_Vector constraints);
SUNDIALS_EXPORT int ARKodeSetMaxHnilWarns(void* arkode_mem, int mxhnil);
SUNDIALS_EXPORT int ARKodeSetInitStep(void* arkode_mem, sunrealtype hin);
SUNDIALS_EXPORT int ARKodeSetMinStep(void* arkode_mem, sunrealtype hmin);
SUNDIALS_EXPORT int ARKodeSetMaxStep(void* arkode_mem, sunrealtype hmax);
SUNDIALS_EXPORT int ARKodeSetMaxNumConstrFails(void* arkode_mem, int maxfails);
SUNDIALS_EXPORT
int ARKodeSetAdjointCheckpointScheme(void* arkode_mem,
                                     SUNAdjointCheckpointScheme checkpoint_scheme);
SUNDIALS_EXPORT
int ARKodeSetAdjointCheckpointIndex(void* arkode_mem, suncountertype step_index);
SUNDIALS_EXPORT
int ARKodeSetUseCompensatedSums(void* arkode_mem, sunbooleantype onoff);
SUNDIALS_EXPORT int ARKodeSetAccumulatedErrorType(void* arkode_mem,
                                                  ARKAccumError accum_type);
SUNDIALS_EXPORT int ARKodeResetAccumulatedError(void* arkode_mem);

/* Integrate the ODE over an interval in t */
SUNDIALS_EXPORT int ARKodeEvolve(void* arkode_mem, sunrealtype tout,
                                 N_Vector yout, sunrealtype* tret, int itask);

/* Computes the kth derivative of the y function at time t */
SUNDIALS_EXPORT int ARKodeGetDky(void* arkode_mem, sunrealtype t, int k,
                                 N_Vector dky);

/* Utility function to update/compute y based on zcor */
SUNDIALS_EXPORT int ARKodeComputeState(void* arkode_mem, N_Vector zcor,
                                       N_Vector z);

/* Optional output functions (general) */
SUNDIALS_EXPORT int ARKodeGetNumRhsEvals(void* arkode_mem, int partition_index,
                                         long int* num_rhs_evals);
SUNDIALS_EXPORT int ARKodeGetNumStepAttempts(void* arkode_mem,
                                             long int* step_attempts);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
int ARKodeGetWorkSpace(void* arkode_mem, long int* lenrw, long int* leniw);
SUNDIALS_EXPORT int ARKodeGetNumSteps(void* arkode_mem, long int* nsteps);
SUNDIALS_EXPORT int ARKodeGetLastStep(void* arkode_mem, sunrealtype* hlast);
SUNDIALS_EXPORT int ARKodeGetCurrentStep(void* arkode_mem, sunrealtype* hcur);
SUNDIALS_EXPORT int ARKodeGetStepDirection(void* arkode_mem,
                                           sunrealtype* stepdir);
SUNDIALS_EXPORT int ARKodeGetErrWeights(void* arkode_mem, N_Vector eweight);
SUNDIALS_EXPORT int ARKodeGetNumGEvals(void* arkode_mem, long int* ngevals);
SUNDIALS_EXPORT int ARKodeGetRootInfo(void* arkode_mem, int* rootsfound);
SUNDIALS_EXPORT int ARKodeGetUserData(void* arkode_mem, void** user_data);
SUNDIALS_EXPORT int ARKodePrintAllStats(void* arkode_mem, FILE* outfile,
                                        SUNOutputFormat fmt);
SUNDIALS_EXPORT char* ARKodeGetReturnFlagName(long int flag);
SUNDIALS_EXPORT int ARKodeWriteParameters(void* arkode_mem, FILE* fp);

/* Optional output functions (temporal adaptivity) */
SUNDIALS_EXPORT int ARKodeGetNumExpSteps(void* arkode_mem, long int* expsteps);
SUNDIALS_EXPORT int ARKodeGetNumAccSteps(void* arkode_mem, long int* accsteps);
SUNDIALS_EXPORT int ARKodeGetNumErrTestFails(void* arkode_mem,
                                             long int* netfails);
SUNDIALS_EXPORT int ARKodeGetEstLocalErrors(void* arkode_mem, N_Vector ele);
SUNDIALS_EXPORT int ARKodeGetActualInitStep(void* arkode_mem,
                                            sunrealtype* hinused);
SUNDIALS_EXPORT int ARKodeGetTolScaleFactor(void* arkode_mem,
                                            sunrealtype* tolsfac);
SUNDIALS_EXPORT int ARKodeGetNumConstrFails(void* arkode_mem,
                                            long int* nconstrfails);
SUNDIALS_EXPORT int ARKodeGetStepStats(void* arkode_mem, long int* nsteps,
                                       sunrealtype* hinused, sunrealtype* hlast,
                                       sunrealtype* hcur, sunrealtype* tcur);
SUNDIALS_EXPORT int ARKodeGetAccumulatedError(void* arkode_mem,
                                              sunrealtype* accum_error);

/* Optional output functions (implicit solver) */
SUNDIALS_EXPORT int ARKodeGetNumLinSolvSetups(void* arkode_mem,
                                              long int* nlinsetups);
SUNDIALS_EXPORT int ARKodeGetCurrentTime(void* arkode_mem, sunrealtype* tcur);
SUNDIALS_EXPORT int ARKodeGetCurrentState(void* arkode_mem,
                                          N_Vector* state); // nb::rv_policy::reference
SUNDIALS_EXPORT int ARKodeGetCurrentGamma(void* arkode_mem, sunrealtype* gamma);
SUNDIALS_EXPORT int ARKodeGetNonlinearSystemData(
  void* arkode_mem, sunrealtype* tcur, N_Vector* zpred, N_Vector* z,
  N_Vector* Fi, sunrealtype* gamma, N_Vector* sdata, void** user_data);
SUNDIALS_EXPORT int ARKodeGetNumNonlinSolvIters(void* arkode_mem,
                                                long int* nniters);
SUNDIALS_EXPORT int ARKodeGetNumNonlinSolvConvFails(void* arkode_mem,
                                                    long int* nnfails);
SUNDIALS_EXPORT int ARKodeGetNonlinSolvStats(void* arkode_mem, long int* nniters,
                                             long int* nnfails);
SUNDIALS_EXPORT int ARKodeGetNumStepSolveFails(void* arkode_mem,
                                               long int* nncfails);
SUNDIALS_EXPORT int ARKodeGetJac(void* arkode_mem,
                                 SUNMatrix* J); // nb::rv_policy::reference
SUNDIALS_EXPORT int ARKodeGetJacTime(void* arkode_mem, sunrealtype* t_J);
SUNDIALS_EXPORT int ARKodeGetJacNumSteps(void* arkode_mem, long int* nst_J);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
int ARKodeGetLinWorkSpace(void* arkode_mem, long int* lenrwLS, long int* leniwLS);
SUNDIALS_EXPORT int ARKodeGetNumJacEvals(void* arkode_mem, long int* njevals);
SUNDIALS_EXPORT int ARKodeGetNumPrecEvals(void* arkode_mem, long int* npevals);
SUNDIALS_EXPORT int ARKodeGetNumPrecSolves(void* arkode_mem, long int* npsolves);
SUNDIALS_EXPORT int ARKodeGetNumLinIters(void* arkode_mem, long int* nliters);
SUNDIALS_EXPORT int ARKodeGetNumLinConvFails(void* arkode_mem,
                                             long int* nlcfails);
SUNDIALS_EXPORT int ARKodeGetNumJTSetupEvals(void* arkode_mem,
                                             long int* njtsetups);
SUNDIALS_EXPORT int ARKodeGetNumJtimesEvals(void* arkode_mem, long int* njvevals);
SUNDIALS_EXPORT int ARKodeGetNumLinRhsEvals(void* arkode_mem,
                                            long int* nfevalsLS);
SUNDIALS_EXPORT int ARKodeGetLastLinFlag(void* arkode_mem, long int* flag);
SUNDIALS_EXPORT char* ARKodeGetLinReturnFlagName(long int flag);

/* Optional output functions (non-identity mass matrices) */
SUNDIALS_EXPORT int ARKodeGetCurrentMassMatrix(void* arkode_mem,
                                               SUNMatrix* M); // nb::rv_policy::reference
SUNDIALS_EXPORT int ARKodeGetResWeights(void* arkode_mem, N_Vector rweight);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
int ARKodeGetMassWorkSpace(void* arkode_mem, long int* lenrwMLS,
                           long int* leniwMLS);
SUNDIALS_EXPORT int ARKodeGetNumMassSetups(void* arkode_mem, long int* nmsetups);
SUNDIALS_EXPORT int ARKodeGetNumMassMultSetups(void* arkode_mem,
                                               long int* nmvsetups);
SUNDIALS_EXPORT int ARKodeGetNumMassMult(void* arkode_mem, long int* nmvevals);
SUNDIALS_EXPORT int ARKodeGetNumMassSolves(void* arkode_mem, long int* nmsolves);
SUNDIALS_EXPORT int ARKodeGetNumMassPrecEvals(void* arkode_mem,
                                              long int* nmpevals);
SUNDIALS_EXPORT int ARKodeGetNumMassPrecSolves(void* arkode_mem,
                                               long int* nmpsolves);
SUNDIALS_EXPORT int ARKodeGetNumMassIters(void* arkode_mem, long int* nmiters);
SUNDIALS_EXPORT int ARKodeGetNumMassConvFails(void* arkode_mem,
                                              long int* nmcfails);
SUNDIALS_EXPORT int ARKodeGetNumMTSetups(void* arkode_mem, long int* nmtsetups);
SUNDIALS_EXPORT int ARKodeGetLastMassFlag(void* arkode_mem, long int* flag);

/* Free function */
SUNDIALS_EXPORT void ARKodeFree(void** arkode_mem);

/* Output the ARKODE memory structure (useful when debugging) */
SUNDIALS_EXPORT void ARKodePrintMem(void* arkode_mem, FILE* outfile);

/* Relaxation functions */
SUNDIALS_EXPORT int ARKodeSetRelaxFn(void* arkode_mem, ARKRelaxFn rfn,
                                     ARKRelaxJacFn rjac);
SUNDIALS_EXPORT int ARKodeSetRelaxEtaFail(void* arkode_mem, sunrealtype eta_rf);
SUNDIALS_EXPORT int ARKodeSetRelaxLowerBound(void* arkode_mem, sunrealtype lower);
SUNDIALS_EXPORT int ARKodeSetRelaxMaxFails(void* arkode_mem, int max_fails);
SUNDIALS_EXPORT int ARKodeSetRelaxMaxIters(void* arkode_mem, int max_iters);
SUNDIALS_EXPORT int ARKodeSetRelaxSolver(void* arkode_mem, ARKRelaxSolver solver);
SUNDIALS_EXPORT int ARKodeSetRelaxResTol(void* arkode_mem, sunrealtype res_tol);
SUNDIALS_EXPORT int ARKodeSetRelaxTol(void* arkode_mem, sunrealtype rel_tol,
                                      sunrealtype abs_tol);
SUNDIALS_EXPORT int ARKodeSetRelaxUpperBound(void* arkode_mem, sunrealtype upper);
SUNDIALS_EXPORT int ARKodeGetNumRelaxFnEvals(void* arkode_mem, long int* r_evals);
SUNDIALS_EXPORT int ARKodeGetNumRelaxJacEvals(void* arkode_mem,
                                              long int* J_evals);
SUNDIALS_EXPORT int ARKodeGetNumRelaxFails(void* arkode_mem,
                                           long int* relax_fails);
SUNDIALS_EXPORT int ARKodeGetNumRelaxBoundFails(void* arkode_mem,
                                                long int* fails);
SUNDIALS_EXPORT int ARKodeGetNumRelaxSolveFails(void* arkode_mem,
                                                long int* fails);
SUNDIALS_EXPORT int ARKodeGetNumRelaxSolveIters(void* arkode_mem,
                                                long int* iters);

/* SUNStepper functions */
SUNDIALS_EXPORT int ARKodeCreateSUNStepper(void* arkode_mem, SUNStepper* stepper);

#ifdef __cplusplus
}
#endif

#endif
