/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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

#ifndef _ARKSTEP_DEPRECATED_H
#define _ARKSTEP_DEPRECATED_H

#include <arkode/arkode.h>
#include <arkode/arkode_ls.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* --------------------------------------------------------------------------
 * Deprecated Functions -- all are superseded by shared ARKODE-level routines
 * -------------------------------------------------------------------------- */

SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeCreateMRIStepInnerStepper instead")
int ARKStepCreateMRIStepInnerStepper(void* arkode_mem,
                                     MRIStepInnerStepper* stepper);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeResize instead")
int ARKStepResize(void* arkode_mem, N_Vector ynew, sunrealtype hscale,
                  sunrealtype t0, ARKVecResizeFn resize, void* resize_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeReset instead")
int ARKStepReset(void* arkode_mem, sunrealtype tR, N_Vector yR);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSStolerances instead")
int ARKStepSStolerances(void* arkode_mem, sunrealtype reltol, sunrealtype abstol);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSVtolerances instead")
int ARKStepSVtolerances(void* arkode_mem, sunrealtype reltol, N_Vector abstol);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeWFtolerances instead")
int ARKStepWFtolerances(void* arkode_mem, ARKEwtFn efun);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeResStolerance instead")
int ARKStepResStolerance(void* arkode_mem, sunrealtype rabstol);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeResVtolerance instead")
int ARKStepResVtolerance(void* arkode_mem, N_Vector rabstol);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeResFtolerance instead")
int ARKStepResFtolerance(void* arkode_mem, ARKRwtFn rfun);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetLinearSolver instead")
int ARKStepSetLinearSolver(void* arkode_mem, SUNLinearSolver LS, SUNMatrix A);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMassLinearSolver instead")
int ARKStepSetMassLinearSolver(void* arkode_mem, SUNLinearSolver LS,
                               SUNMatrix M, sunbooleantype time_dep);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeRootInit instead")
int ARKStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetDefaults instead")
int ARKStepSetDefaults(void* arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("adjust parameters individually instead")
int ARKStepSetOptimalParams(void* arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetOrder instead")
int ARKStepSetOrder(void* arkode_mem, int maxord);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInterpolantType instead")
int ARKStepSetInterpolantType(void* arkode_mem, int itype);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInterpolantDegree instead")
int ARKStepSetInterpolantDegree(void* arkode_mem, int degree);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInterpolantDegree instead")
int ARKStepSetDenseOrder(void* arkode_mem, int dord);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNonlinearSolver instead")
int ARKStepSetNonlinearSolver(void* arkode_mem, SUNNonlinearSolver NLS);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNlsRhsFn instead")
int ARKStepSetNlsRhsFn(void* arkode_mem, ARKRhsFn nls_fi);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetLinear instead")
int ARKStepSetLinear(void* arkode_mem, int timedepend);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNonlinear instead")
int ARKStepSetNonlinear(void* arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetDeduceImplicitRhs instead")
int ARKStepSetDeduceImplicitRhs(void* arkode_mem, sunbooleantype deduce);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetAdaptController instead")
int ARKStepSetAdaptController(void* arkode_mem, SUNAdaptController C);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetAdaptivityAdjustment instead")
int ARKStepSetAdaptivityAdjustment(void* arkode_mem, int adjust);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetCFLFraction instead")
int ARKStepSetCFLFraction(void* arkode_mem, sunrealtype cfl_frac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetSafetyFactor instead")
int ARKStepSetSafetyFactor(void* arkode_mem, sunrealtype safety);
SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNAdaptController instead")
int ARKStepSetErrorBias(void* arkode_mem, sunrealtype bias);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxGrowth instead")
int ARKStepSetMaxGrowth(void* arkode_mem, sunrealtype mx_growth);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMinReduction instead")
int ARKStepSetMinReduction(void* arkode_mem, sunrealtype eta_min);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetFixedStepBounds instead")
int ARKStepSetFixedStepBounds(void* arkode_mem, sunrealtype lb, sunrealtype ub);
SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNAdaptController instead")
int ARKStepSetAdaptivityMethod(void* arkode_mem, int imethod, int idefault,
                               int pq, sunrealtype adapt_params[3]);
SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNAdaptController instead")
int ARKStepSetAdaptivityFn(void* arkode_mem, ARKAdaptFn hfun, void* h_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxFirstGrowth instead")
int ARKStepSetMaxFirstGrowth(void* arkode_mem, sunrealtype etamx1);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxEFailGrowth instead")
int ARKStepSetMaxEFailGrowth(void* arkode_mem, sunrealtype etamxf);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetSmallNumEFails instead")
int ARKStepSetSmallNumEFails(void* arkode_mem, int small_nef);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxCFailGrowth instead")
int ARKStepSetMaxCFailGrowth(void* arkode_mem, sunrealtype etacf);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNonlinCRDown instead")
int ARKStepSetNonlinCRDown(void* arkode_mem, sunrealtype crdown);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNonlinRDiv instead")
int ARKStepSetNonlinRDiv(void* arkode_mem, sunrealtype rdiv);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetDeltaGammaMax instead")
int ARKStepSetDeltaGammaMax(void* arkode_mem, sunrealtype dgmax);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetLSetupFrequency instead")
int ARKStepSetLSetupFrequency(void* arkode_mem, int msbp);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetPredictorMethod instead")
int ARKStepSetPredictorMethod(void* arkode_mem, int method);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetStabilityFn instead")
int ARKStepSetStabilityFn(void* arkode_mem, ARKExpStabFn EStab, void* estab_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxErrTestFails instead")
int ARKStepSetMaxErrTestFails(void* arkode_mem, int maxnef);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxNonlinIters instead")
int ARKStepSetMaxNonlinIters(void* arkode_mem, int maxcor);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxConvFails instead")
int ARKStepSetMaxConvFails(void* arkode_mem, int maxncf);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNonlinConvCoef instead")
int ARKStepSetNonlinConvCoef(void* arkode_mem, sunrealtype nlscoef);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetConstraints instead")
int ARKStepSetConstraints(void* arkode_mem, N_Vector constraints);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxNumSteps instead")
int ARKStepSetMaxNumSteps(void* arkode_mem, long int mxsteps);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxHnilWarns instead")
int ARKStepSetMaxHnilWarns(void* arkode_mem, int mxhnil);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInitStep instead")
int ARKStepSetInitStep(void* arkode_mem, sunrealtype hin);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMinStep instead")
int ARKStepSetMinStep(void* arkode_mem, sunrealtype hmin);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxStep instead")
int ARKStepSetMaxStep(void* arkode_mem, sunrealtype hmax);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInterpolateStopTime instead")
int ARKStepSetInterpolateStopTime(void* arkode_mem, sunbooleantype interp);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetStopTime instead")
int ARKStepSetStopTime(void* arkode_mem, sunrealtype tstop);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeClearStopTime instead")
int ARKStepClearStopTime(void* arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetFixedStep instead")
int ARKStepSetFixedStep(void* arkode_mem, sunrealtype hfixed);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxNumConstrFails instead")
int ARKStepSetMaxNumConstrFails(void* arkode_mem, int maxfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRootDirection instead")
int ARKStepSetRootDirection(void* arkode_mem, int* rootdir);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNoInactiveRootWarn instead")
int ARKStepSetNoInactiveRootWarn(void* arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetUserData instead")
int ARKStepSetUserData(void* arkode_mem, void* user_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetPostprocessStepFn instead")
int ARKStepSetPostprocessStepFn(void* arkode_mem, ARKPostProcessFn ProcessStep);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetPostprocessStageFn instead")
int ARKStepSetPostprocessStageFn(void* arkode_mem, ARKPostProcessFn ProcessStage);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetStagePredictFn instead")
int ARKStepSetStagePredictFn(void* arkode_mem, ARKStagePredictFn PredictStage);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetJacFn instead")
int ARKStepSetJacFn(void* arkode_mem, ARKLsJacFn jac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMassFn instead")
int ARKStepSetMassFn(void* arkode_mem, ARKLsMassFn mass);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetJacEvalFrequency instead")
int ARKStepSetJacEvalFrequency(void* arkode_mem, long int msbj);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetLinearSolutionScaling instead")
int ARKStepSetLinearSolutionScaling(void* arkode_mem, sunbooleantype onoff);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetEpsLin instead")
int ARKStepSetEpsLin(void* arkode_mem, sunrealtype eplifac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMassEpsLin instead")
int ARKStepSetMassEpsLin(void* arkode_mem, sunrealtype eplifac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetLSNormFactor instead")
int ARKStepSetLSNormFactor(void* arkode_mem, sunrealtype nrmfac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMassLSNormFactor instead")
int ARKStepSetMassLSNormFactor(void* arkode_mem, sunrealtype nrmfac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetPreconditioner instead")
int ARKStepSetPreconditioner(void* arkode_mem, ARKLsPrecSetupFn psetup,
                             ARKLsPrecSolveFn psolve);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMassPreconditioner instead")
int ARKStepSetMassPreconditioner(void* arkode_mem, ARKLsMassPrecSetupFn psetup,
                                 ARKLsMassPrecSolveFn psolve);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetJacTimes instead")
int ARKStepSetJacTimes(void* arkode_mem, ARKLsJacTimesSetupFn jtsetup,
                       ARKLsJacTimesVecFn jtimes);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetJacTimesRhsFn instead")
int ARKStepSetJacTimesRhsFn(void* arkode_mem, ARKRhsFn jtimesRhsFn);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMassTimes instead")
int ARKStepSetMassTimes(void* arkode_mem, ARKLsMassTimesSetupFn msetup,
                        ARKLsMassTimesVecFn mtimes, void* mtimes_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetLinSysFn instead")
int ARKStepSetLinSysFn(void* arkode_mem, ARKLsLinSysFn linsys);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeEvolve instead")
int ARKStepEvolve(void* arkode_mem, sunrealtype tout, N_Vector yout,
                  sunrealtype* tret, int itask);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetDky instead")
int ARKStepGetDky(void* arkode_mem, sunrealtype t, int k, N_Vector dky);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeComputeState instead")
int ARKStepComputeState(void* arkode_mem, N_Vector zcor, N_Vector z);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumExpSteps instead")
int ARKStepGetNumExpSteps(void* arkode_mem, long int* expsteps);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumAccSteps instead")
int ARKStepGetNumAccSteps(void* arkode_mem, long int* accsteps);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumStepAttempts instead")
int ARKStepGetNumStepAttempts(void* arkode_mem, long int* step_attempts);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumLinSolvSetups instead")
int ARKStepGetNumLinSolvSetups(void* arkode_mem, long int* nlinsetups);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumErrTestFails instead")
int ARKStepGetNumErrTestFails(void* arkode_mem, long int* netfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetEstLocalErrors instead")
int ARKStepGetEstLocalErrors(void* arkode_mem, N_Vector ele);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
int ARKStepGetWorkSpace(void* arkode_mem, long int* lenrw, long int* leniw);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumSteps instead")
int ARKStepGetNumSteps(void* arkode_mem, long int* nsteps);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetActualInitStep instead")
int ARKStepGetActualInitStep(void* arkode_mem, sunrealtype* hinused);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetLastStep instead")
int ARKStepGetLastStep(void* arkode_mem, sunrealtype* hlast);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetCurrentStep instead")
int ARKStepGetCurrentStep(void* arkode_mem, sunrealtype* hcur);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetCurrentTime instead")
int ARKStepGetCurrentTime(void* arkode_mem, sunrealtype* tcur);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetCurrentState instead")
int ARKStepGetCurrentState(void* arkode_mem, N_Vector* state);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetCurrentGamma instead")
int ARKStepGetCurrentGamma(void* arkode_mem, sunrealtype* gamma);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetCurrentMassMatrix instead")
int ARKStepGetCurrentMassMatrix(void* arkode_mem, SUNMatrix* M);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetTolScaleFactor instead")
int ARKStepGetTolScaleFactor(void* arkode_mem, sunrealtype* tolsfac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetErrWeights instead")
int ARKStepGetErrWeights(void* arkode_mem, N_Vector eweight);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetResWeights instead")
int ARKStepGetResWeights(void* arkode_mem, N_Vector rweight);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumGEvals instead")
int ARKStepGetNumGEvals(void* arkode_mem, long int* ngevals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetRootInfo instead")
int ARKStepGetRootInfo(void* arkode_mem, int* rootsfound);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumConstrFails instead")
int ARKStepGetNumConstrFails(void* arkode_mem, long int* nconstrfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetUserData instead")
int ARKStepGetUserData(void* arkode_mem, void** user_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodePrintAllStats instead")
int ARKStepPrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetReturnFlagName instead")
char* ARKStepGetReturnFlagName(long int flag);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeWriteParameters instead")
int ARKStepWriteParameters(void* arkode_mem, FILE* fp);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use ARKStepGetCurrentButcherTables and ARKodeButcherTable_Write instead")
int ARKStepWriteButcher(void* arkode_mem, FILE* fp);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetStepStats instead")
int ARKStepGetStepStats(void* arkode_mem, long int* nsteps, sunrealtype* hinused,
                        sunrealtype* hlast, sunrealtype* hcur, sunrealtype* tcur);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNonlinearSystemData instead")
int ARKStepGetNonlinearSystemData(void* arkode_mem, sunrealtype* tcur,
                                  N_Vector* zpred, N_Vector* z, N_Vector* Fi,
                                  sunrealtype* gamma, N_Vector* sdata,
                                  void** user_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumNonlinSolvIters instead")
int ARKStepGetNumNonlinSolvIters(void* arkode_mem, long int* nniters);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumNonlinSolvConvFails instead")
int ARKStepGetNumNonlinSolvConvFails(void* arkode_mem, long int* nnfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNonlinSolvStats instead")
int ARKStepGetNonlinSolvStats(void* arkode_mem, long int* nniters,
                              long int* nnfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumStepSolveFails instead")
int ARKStepGetNumStepSolveFails(void* arkode_mem, long int* nncfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetJac instead")
int ARKStepGetJac(void* arkode_mem, SUNMatrix* J);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetJacTime instead")
int ARKStepGetJacTime(void* arkode_mem, sunrealtype* t_J);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetJacNumSteps instead")
int ARKStepGetJacNumSteps(void* arkode_mem, long int* nst_J);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
int ARKStepGetLinWorkSpace(void* arkode_mem, long int* lenrwLS,
                           long int* leniwLS);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumJacEvals instead")
int ARKStepGetNumJacEvals(void* arkode_mem, long int* njevals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumPrecEvals instead")
int ARKStepGetNumPrecEvals(void* arkode_mem, long int* npevals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumPrecSolves instead")
int ARKStepGetNumPrecSolves(void* arkode_mem, long int* npsolves);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumLinIters instead")
int ARKStepGetNumLinIters(void* arkode_mem, long int* nliters);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumLinConvFails instead")
int ARKStepGetNumLinConvFails(void* arkode_mem, long int* nlcfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumJTSetupEvals instead")
int ARKStepGetNumJTSetupEvals(void* arkode_mem, long int* njtsetups);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumJtimesEvals instead")
int ARKStepGetNumJtimesEvals(void* arkode_mem, long int* njvevals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumLinRhsEvals instead")
int ARKStepGetNumLinRhsEvals(void* arkode_mem, long int* nfevalsLS);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetLastLinFlag instead")
int ARKStepGetLastLinFlag(void* arkode_mem, long int* flag);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
int ARKStepGetMassWorkSpace(void* arkode_mem, long int* lenrwMLS,
                            long int* leniwMLS);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumMassSetups instead")
int ARKStepGetNumMassSetups(void* arkode_mem, long int* nmsetups);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumMassMultSetups instead")
int ARKStepGetNumMassMultSetups(void* arkode_mem, long int* nmvsetups);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumMassMult instead")
int ARKStepGetNumMassMult(void* arkode_mem, long int* nmvevals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumMassSolves instead")
int ARKStepGetNumMassSolves(void* arkode_mem, long int* nmsolves);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumMassPrecEvals instead")
int ARKStepGetNumMassPrecEvals(void* arkode_mem, long int* nmpevals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumMassPrecSolves instead")
int ARKStepGetNumMassPrecSolves(void* arkode_mem, long int* nmpsolves);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumMassIters instead")
int ARKStepGetNumMassIters(void* arkode_mem, long int* nmiters);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumMassConvFails instead")
int ARKStepGetNumMassConvFails(void* arkode_mem, long int* nmcfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumMTSetups instead")
int ARKStepGetNumMTSetups(void* arkode_mem, long int* nmtsetups);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetLastMassFlag instead")
int ARKStepGetLastMassFlag(void* arkode_mem, long int* flag);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetLinReturnFlagName instead")
char* ARKStepGetLinReturnFlagName(long int flag);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeFree instead")
void ARKStepFree(void** arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodePrintMem instead")
void ARKStepPrintMem(void* arkode_mem, FILE* outfile);

/* Relaxation functions */
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRelaxFn instead")
int ARKStepSetRelaxFn(void* arkode_mem, ARKRelaxFn rfn, ARKRelaxJacFn rjac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRelaxEtaFail instead")
int ARKStepSetRelaxEtaFail(void* arkode_mem, sunrealtype eta_rf);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRelaxLowerBound instead")
int ARKStepSetRelaxLowerBound(void* arkode_mem, sunrealtype lower);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRelaxMaxFails instead")
int ARKStepSetRelaxMaxFails(void* arkode_mem, int max_fails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRelaxMaxIters instead")
int ARKStepSetRelaxMaxIters(void* arkode_mem, int max_iters);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRelaxSolver instead")
int ARKStepSetRelaxSolver(void* arkode_mem, ARKRelaxSolver solver);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRelaxResTol instead")
int ARKStepSetRelaxResTol(void* arkode_mem, sunrealtype res_tol);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRelaxTol instead")
int ARKStepSetRelaxTol(void* arkode_mem, sunrealtype rel_tol,
                       sunrealtype abs_tol);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRelaxUpperBound instead")
int ARKStepSetRelaxUpperBound(void* arkode_mem, sunrealtype upper);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumRelaxFnEvals instead")
int ARKStepGetNumRelaxFnEvals(void* arkode_mem, long int* r_evals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumRelaxJacEvals instead")
int ARKStepGetNumRelaxJacEvals(void* arkode_mem, long int* J_evals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumRelaxFails instead")
int ARKStepGetNumRelaxFails(void* arkode_mem, long int* relax_fails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumRelaxBoundFails instead")
int ARKStepGetNumRelaxBoundFails(void* arkode_mem, long int* fails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumRelaxSolveFails instead")
int ARKStepGetNumRelaxSolveFails(void* arkode_mem, long int* fails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumRelaxSolveIters instead")
int ARKStepGetNumRelaxSolveIters(void* arkode_mem, long int* iters);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumRhsEvals instead")
int ARKStepGetNumRhsEvals(void* arkode_mem, long int* nfe_evals,
                          long int* nfi_evals);

#ifdef __cplusplus
}
#endif

#endif /* _ARKSTEP_DEPRECATED_H */
