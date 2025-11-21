/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * ---------------------------------------------------------------------------*/

#ifndef _MRISTEP_DEPRECATED_H
#define _MRISTEP_DEPRECATED_H

#include <arkode/arkode.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* --------------------------------------------------------------------------
 * Deprecated Functions -- all are superseded by shared ARKODE-level routines
 * -------------------------------------------------------------------------- */

SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeResize instead")
int MRIStepResize(void* arkode_mem, N_Vector ynew, sunrealtype t0,
                  ARKVecResizeFn resize, void* resize_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeReset instead")
int MRIStepReset(void* arkode_mem, sunrealtype tR, N_Vector yR);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSStolerances instead")
int MRIStepSStolerances(void* arkode_mem, sunrealtype reltol, sunrealtype abstol);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSVtolerances instead")
int MRIStepSVtolerances(void* arkode_mem, sunrealtype reltol, N_Vector abstol);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeWFtolerances instead")
int MRIStepWFtolerances(void* arkode_mem, ARKEwtFn efun);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetLinearSolver instead")
int MRIStepSetLinearSolver(void* arkode_mem, SUNLinearSolver LS, SUNMatrix A);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeRootInit instead")
int MRIStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetDefaults instead")
int MRIStepSetDefaults(void* arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetOrder instead")
int MRIStepSetOrder(void* arkode_mem, int ord);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInterpolantType instead")
int MRIStepSetInterpolantType(void* arkode_mem, int itype);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInterpolantDegree instead")
int MRIStepSetInterpolantDegree(void* arkode_mem, int degree);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInterpolantDegree instead")
int MRIStepSetDenseOrder(void* arkode_mem, int dord);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNonlinearSolver instead")
int MRIStepSetNonlinearSolver(void* arkode_mem, SUNNonlinearSolver NLS);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNlsRhsFn instead")
int MRIStepSetNlsRhsFn(void* arkode_mem, ARKRhsFn nls_fs);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetLinear instead")
int MRIStepSetLinear(void* arkode_mem, int timedepend);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNonlinear instead")
int MRIStepSetNonlinear(void* arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxNumSteps instead")
int MRIStepSetMaxNumSteps(void* arkode_mem, long int mxsteps);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNonlinCRDown instead")
int MRIStepSetNonlinCRDown(void* arkode_mem, sunrealtype crdown);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNonlinRDiv instead")
int MRIStepSetNonlinRDiv(void* arkode_mem, sunrealtype rdiv);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetDeltaGammaMax instead")
int MRIStepSetDeltaGammaMax(void* arkode_mem, sunrealtype dgmax);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetLSetupFrequency instead")
int MRIStepSetLSetupFrequency(void* arkode_mem, int msbp);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetPredictorMethod instead")
int MRIStepSetPredictorMethod(void* arkode_mem, int method);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxNonlinIters instead")
int MRIStepSetMaxNonlinIters(void* arkode_mem, int maxcor);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNonlinConvCoef instead")
int MRIStepSetNonlinConvCoef(void* arkode_mem, sunrealtype nlscoef);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxHnilWarns instead")
int MRIStepSetMaxHnilWarns(void* arkode_mem, int mxhnil);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInterpolateStopTime instead")
int MRIStepSetInterpolateStopTime(void* arkode_mem, sunbooleantype interp);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetStopTime instead")
int MRIStepSetStopTime(void* arkode_mem, sunrealtype tstop);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeClearStopTime instead")
int MRIStepClearStopTime(void* arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetFixedStep instead")
int MRIStepSetFixedStep(void* arkode_mem, sunrealtype hsfixed);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRootDirection instead")
int MRIStepSetRootDirection(void* arkode_mem, int* rootdir);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNoInactiveRootWarn instead")
int MRIStepSetNoInactiveRootWarn(void* arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetUserData instead")
int MRIStepSetUserData(void* arkode_mem, void* user_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetPostprocessStepFn instead")
int MRIStepSetPostprocessStepFn(void* arkode_mem, ARKPostProcessFn ProcessStep);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetPostprocessStageFn instead")
int MRIStepSetPostprocessStageFn(void* arkode_mem, ARKPostProcessFn ProcessStage);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetStagePredictFn instead")
int MRIStepSetStagePredictFn(void* arkode_mem, ARKStagePredictFn PredictStage);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetDeduceImplicitRhs instead")
int MRIStepSetDeduceImplicitRhs(void* arkode_mem, sunbooleantype deduce);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetJacFn instead")
int MRIStepSetJacFn(void* arkode_mem, ARKLsJacFn jac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetJacEvalFrequency instead")
int MRIStepSetJacEvalFrequency(void* arkode_mem, long int msbj);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetLinearSolutionScaling instead")
int MRIStepSetLinearSolutionScaling(void* arkode_mem, sunbooleantype onoff);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetEpsLin instead")
int MRIStepSetEpsLin(void* arkode_mem, sunrealtype eplifac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetLSNormFactor instead")
int MRIStepSetLSNormFactor(void* arkode_mem, sunrealtype nrmfac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetPreconditioner instead")
int MRIStepSetPreconditioner(void* arkode_mem, ARKLsPrecSetupFn psetup,
                             ARKLsPrecSolveFn psolve);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetJacTimes instead")
int MRIStepSetJacTimes(void* arkode_mem, ARKLsJacTimesSetupFn jtsetup,
                       ARKLsJacTimesVecFn jtimes);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetJacTimesRhsFn instead")
int MRIStepSetJacTimesRhsFn(void* arkode_mem, ARKRhsFn jtimesRhsFn);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetLinSysFn instead")
int MRIStepSetLinSysFn(void* arkode_mem, ARKLsLinSysFn linsys);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeEvolve instead")
int MRIStepEvolve(void* arkode_mem, sunrealtype tout, N_Vector yout,
                  sunrealtype* tret, int itask);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetDky instead")
int MRIStepGetDky(void* arkode_mem, sunrealtype t, int k, N_Vector dky);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeComputeState instead")
int MRIStepComputeState(void* arkode_mem, N_Vector zcor, N_Vector z);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumLinSolvSetups instead")
int MRIStepGetNumLinSolvSetups(void* arkode_mem, long int* nlinsetups);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
int MRIStepGetWorkSpace(void* arkode_mem, long int* lenrw, long int* leniw);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumSteps instead")
int MRIStepGetNumSteps(void* arkode_mem, long int* nssteps);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetLastStep instead")
int MRIStepGetLastStep(void* arkode_mem, sunrealtype* hlast);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetCurrentTime instead")
int MRIStepGetCurrentTime(void* arkode_mem, sunrealtype* tcur);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetCurrentState instead")
int MRIStepGetCurrentState(void* arkode_mem, N_Vector* state);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetCurrentGamma instead")
int MRIStepGetCurrentGamma(void* arkode_mem, sunrealtype* gamma);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetTolScaleFactor instead")
int MRIStepGetTolScaleFactor(void* arkode_mem, sunrealtype* tolsfac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetErrWeights instead")
int MRIStepGetErrWeights(void* arkode_mem, N_Vector eweight);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumGEvals instead")
int MRIStepGetNumGEvals(void* arkode_mem, long int* ngevals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetRootInfo instead")
int MRIStepGetRootInfo(void* arkode_mem, int* rootsfound);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetUserData instead")
int MRIStepGetUserData(void* arkode_mem, void** user_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodePrintAllStats instead")
int MRIStepPrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetReturnFlagName instead")
char* MRIStepGetReturnFlagName(long int flag);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeWriteParameters instead")
int MRIStepWriteParameters(void* arkode_mem, FILE* fp);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use MRIStepGetCurrentCoupling and MRIStepCoupling_Write instead")
int MRIStepWriteCoupling(void* arkode_mem, FILE* fp);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNonlinearSystemData instead")
int MRIStepGetNonlinearSystemData(void* arkode_mem, sunrealtype* tcur,
                                  N_Vector* zpred, N_Vector* z, N_Vector* F,
                                  sunrealtype* gamma, N_Vector* sdata,
                                  void** user_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumNonlinSolvIters instead")
int MRIStepGetNumNonlinSolvIters(void* arkode_mem, long int* nniters);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumNonlinSolvConvFails instead")
int MRIStepGetNumNonlinSolvConvFails(void* arkode_mem, long int* nnfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNonlinSolvStats instead")
int MRIStepGetNonlinSolvStats(void* arkode_mem, long int* nniters,
                              long int* nnfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumStepSolveFails instead")
int MRIStepGetNumStepSolveFails(void* arkode_mem, long int* nncfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetJac instead")
int MRIStepGetJac(void* arkode_mem, SUNMatrix* J);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetJacTime instead")
int MRIStepGetJacTime(void* arkode_mem, sunrealtype* t_J);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetJacNumSteps instead")
int MRIStepGetJacNumSteps(void* arkode_mem, long* nst_J);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
int MRIStepGetLinWorkSpace(void* arkode_mem, long int* lenrwLS,
                           long int* leniwLS);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumJacEvals instead")
int MRIStepGetNumJacEvals(void* arkode_mem, long int* njevals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumPrecEvals instead")
int MRIStepGetNumPrecEvals(void* arkode_mem, long int* npevals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumPrecSolves instead")
int MRIStepGetNumPrecSolves(void* arkode_mem, long int* npsolves);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumLinIters instead")
int MRIStepGetNumLinIters(void* arkode_mem, long int* nliters);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumLinConvFails instead")
int MRIStepGetNumLinConvFails(void* arkode_mem, long int* nlcfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumJTSetupEvals instead")
int MRIStepGetNumJTSetupEvals(void* arkode_mem, long int* njtsetups);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumJtimesEvals instead")
int MRIStepGetNumJtimesEvals(void* arkode_mem, long int* njvevals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumLinRhsEvals instead")
int MRIStepGetNumLinRhsEvals(void* arkode_mem, long int* nfevalsLS);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetLastLinFlag instead")
int MRIStepGetLastLinFlag(void* arkode_mem, long int* flag);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetLinReturnFlagName instead")
char* MRIStepGetLinReturnFlagName(long int flag);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeFree instead")
void MRIStepFree(void** arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodePrintMem instead")
void MRIStepPrintMem(void* arkode_mem, FILE* outfile);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumRhsEvals instead")
int MRIStepGetNumRhsEvals(void* arkode_mem, long int* nfse_evals,
                          long int* nfsi_evals);

#ifdef __cplusplus
}
#endif

#endif
