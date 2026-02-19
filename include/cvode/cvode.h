/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                and Dan Shumaker @ LLNL
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
 * This is the header file for the main CVODE integrator.
 * -----------------------------------------------------------------*/

#ifndef _CVODE_H
#define _CVODE_H

#include <cvode/cvode_ls.h>
#include <cvode/cvode_proj.h>
#include <stdio.h>
#include <sundials/sundials_core.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------
 * CVODE Constants
 * ----------------- */

/* lmm */
#define CV_ADAMS 1
#define CV_BDF   2

/* itask */
#define CV_NORMAL   1
#define CV_ONE_STEP 2

/* return values */

#define CV_SUCCESS      0
#define CV_TSTOP_RETURN 1
#define CV_ROOT_RETURN  2

#define CV_WARNING 99

#define CV_TOO_MUCH_WORK -1
#define CV_TOO_MUCH_ACC  -2
#define CV_ERR_FAILURE   -3
#define CV_CONV_FAILURE  -4

#define CV_LINIT_FAIL        -5
#define CV_LSETUP_FAIL       -6
#define CV_LSOLVE_FAIL       -7
#define CV_RHSFUNC_FAIL      -8
#define CV_FIRST_RHSFUNC_ERR -9
#define CV_REPTD_RHSFUNC_ERR -10
#define CV_UNREC_RHSFUNC_ERR -11
#define CV_RTFUNC_FAIL       -12
#define CV_NLS_INIT_FAIL     -13
#define CV_NLS_SETUP_FAIL    -14
#define CV_CONSTR_FAIL       -15
#define CV_NLS_FAIL          -16

#define CV_MEM_FAIL     -20
#define CV_MEM_NULL     -21
#define CV_ILL_INPUT    -22
#define CV_NO_MALLOC    -23
#define CV_BAD_K        -24
#define CV_BAD_T        -25
#define CV_BAD_DKY      -26
#define CV_TOO_CLOSE    -27
#define CV_VECTOROP_ERR -28

#define CV_PROJ_MEM_NULL      -29
#define CV_PROJFUNC_FAIL      -30
#define CV_REPTD_PROJFUNC_ERR -31

#define CV_CONTEXT_ERR -32

#define CV_UNRECOGNIZED_ERR -99

/* ------------------------------
 * User-Supplied Function Types
 * ------------------------------ */

typedef int (*CVRhsFn)(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

typedef int (*CVRootFn)(sunrealtype t, N_Vector y, sunrealtype* gout_1d,
                        void* user_data);

typedef int (*CVEwtFn)(N_Vector y, N_Vector ewt, void* user_data);

typedef int (*CVMonitorFn)(void* cvode_mem, void* user_data);

/* -------------------
 * Exported Functions
 * ------------------- */

/* Initialization functions */

SUNDIALS_EXPORT void* CVodeCreate(int lmm, SUNContext sunctx);

SUNDIALS_EXPORT int CVodeInit(void* cvode_mem, CVRhsFn f, sunrealtype t0,
                              N_Vector y0);
SUNDIALS_EXPORT int CVodeReInit(void* cvode_mem, sunrealtype t0, N_Vector y0);
SUNDIALS_EXPORT int CVodeResizeHistory(void* cvode_mem, sunrealtype* t_hist_1d,
                                       N_Vector* y_hist_1d, N_Vector* f_hist_1d,
                                       int num_y_hist, int num_f_hist);

/* Tolerance input functions */

SUNDIALS_EXPORT int CVodeSStolerances(void* cvode_mem, sunrealtype reltol,
                                      sunrealtype abstol);
SUNDIALS_EXPORT int CVodeSVtolerances(void* cvode_mem, sunrealtype reltol,
                                      N_Vector abstol);
SUNDIALS_EXPORT int CVodeWFtolerances(void* cvode_mem, CVEwtFn efun);

/* Optional input functions */

SUNDIALS_EXPORT int CVodeSetOptions(void* cvode_mem, const char* cvid,
                                    const char* file_name, int argc,
                                    char* argv[]);
SUNDIALS_EXPORT int CVodeSetConstraints(void* cvode_mem, N_Vector constraints);
SUNDIALS_EXPORT int CVodeSetMaxNumConstraintFails(void* cvode_mem, int max_fails);
SUNDIALS_EXPORT int CVodeSetDeltaGammaMaxLSetup(void* cvode_mem,
                                                sunrealtype dgmax_lsetup);
SUNDIALS_EXPORT int CVodeSetInitStep(void* cvode_mem, sunrealtype hin);
SUNDIALS_EXPORT int CVodeSetLSetupFrequency(void* cvode_mem, long int msbp);
SUNDIALS_EXPORT int CVodeSetMaxConvFails(void* cvode_mem, int maxncf);
SUNDIALS_EXPORT int CVodeSetMaxErrTestFails(void* cvode_mem, int maxnef);
SUNDIALS_EXPORT int CVodeSetMaxHnilWarns(void* cvode_mem, int mxhnil);
SUNDIALS_EXPORT int CVodeSetMaxNonlinIters(void* cvode_mem, int maxcor);
SUNDIALS_EXPORT int CVodeSetMaxNumSteps(void* cvode_mem, long int mxsteps);
SUNDIALS_EXPORT int CVodeSetMaxOrd(void* cvode_mem, int maxord);
SUNDIALS_EXPORT int CVodeSetMaxStep(void* cvode_mem, sunrealtype hmax);
SUNDIALS_EXPORT int CVodeSetMinStep(void* cvode_mem, sunrealtype hmin);
SUNDIALS_EXPORT int CVodeSetMonitorFn(void* cvode_mem, CVMonitorFn fn);
SUNDIALS_EXPORT int CVodeSetMonitorFrequency(void* cvode_mem, long int nst);
SUNDIALS_EXPORT int CVodeSetNlsRhsFn(void* cvode_mem, CVRhsFn f);
SUNDIALS_EXPORT int CVodeSetNonlinConvCoef(void* cvode_mem, sunrealtype nlscoef);
SUNDIALS_EXPORT int CVodeSetNonlinearSolver(void* cvode_mem,
                                            SUNNonlinearSolver NLS);
SUNDIALS_EXPORT int CVodeSetStabLimDet(void* cvode_mem, sunbooleantype stldet);
SUNDIALS_EXPORT int CVodeSetStopTime(void* cvode_mem, sunrealtype tstop);
SUNDIALS_EXPORT int CVodeSetInterpolateStopTime(void* cvode_mem,
                                                sunbooleantype interp);
SUNDIALS_EXPORT int CVodeClearStopTime(void* cvode_mem);
SUNDIALS_EXPORT int CVodeSetUseIntegratorFusedKernels(void* cvode_mem,
                                                      sunbooleantype onoff);
SUNDIALS_EXPORT int CVodeSetUserData(void* cvode_mem, void* user_data);

/* Optional step adaptivity input functions */

SUNDIALS_EXPORT
int CVodeSetEtaFixedStepBounds(void* cvode_mem, sunrealtype eta_min_fx,
                               sunrealtype eta_max_fx);
SUNDIALS_EXPORT
int CVodeSetEtaMaxFirstStep(void* cvode_mem, sunrealtype eta_max_fs);
SUNDIALS_EXPORT
int CVodeSetEtaMaxEarlyStep(void* cvode_mem, sunrealtype eta_max_es);
SUNDIALS_EXPORT
int CVodeSetNumStepsEtaMaxEarlyStep(void* cvode_mem, long int small_nst);
SUNDIALS_EXPORT
int CVodeSetEtaMax(void* cvode_mem, sunrealtype eta_max_gs);
SUNDIALS_EXPORT
int CVodeSetEtaMin(void* cvode_mem, sunrealtype eta_min);
SUNDIALS_EXPORT
int CVodeSetEtaMinErrFail(void* cvode_mem, sunrealtype eta_min_ef);
SUNDIALS_EXPORT
int CVodeSetEtaMaxErrFail(void* cvode_mem, sunrealtype eta_max_ef);
SUNDIALS_EXPORT
int CVodeSetNumFailsEtaMaxErrFail(void* cvode_mem, int small_nef);
SUNDIALS_EXPORT
int CVodeSetEtaConvFail(void* cvode_mem, sunrealtype eta_cf);

/* Rootfinding initialization function */
SUNDIALS_EXPORT int CVodeRootInit(void* cvode_mem, int nrtfn, CVRootFn g);

/* Rootfinding optional input functions */

SUNDIALS_EXPORT int CVodeSetRootDirection(void* cvode_mem, int* rootdir_1d);
SUNDIALS_EXPORT int CVodeSetNoInactiveRootWarn(void* cvode_mem);

/* Solver function */
SUNDIALS_EXPORT int CVode(void* cvode_mem, sunrealtype tout, N_Vector yout,
                          sunrealtype* tret, int itask);

/* Utility functions to update/compute y based on ycor */

SUNDIALS_EXPORT int CVodeComputeState(void* cvode_mem, N_Vector ycor, N_Vector y);

/* Dense output function */
SUNDIALS_EXPORT int CVodeGetDky(void* cvode_mem, sunrealtype t, int k,
                                N_Vector dky);

/* Optional output functions */

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
int CVodeGetWorkSpace(void* cvode_mem, long int* lenrw, long int* leniw);
SUNDIALS_EXPORT int CVodeGetNumSteps(void* cvode_mem, long int* nsteps);
SUNDIALS_EXPORT int CVodeGetNumRhsEvals(void* cvode_mem, long int* nfevals);
SUNDIALS_EXPORT int CVodeGetNumLinSolvSetups(void* cvode_mem,
                                             long int* nlinsetups);
SUNDIALS_EXPORT int CVodeGetNumErrTestFails(void* cvode_mem, long int* netfails);
SUNDIALS_EXPORT int CVodeGetLastOrder(void* cvode_mem, int* qlast);
SUNDIALS_EXPORT int CVodeGetCurrentOrder(void* cvode_mem, int* qcur);
SUNDIALS_EXPORT int CVodeGetCurrentGamma(void* cvode_mem, sunrealtype* gamma);
SUNDIALS_EXPORT int CVodeGetNumStabLimOrderReds(void* cvode_mem,
                                                long int* nslred);
SUNDIALS_EXPORT int CVodeGetActualInitStep(void* cvode_mem, sunrealtype* hinused);
SUNDIALS_EXPORT int CVodeGetLastStep(void* cvode_mem, sunrealtype* hlast);
SUNDIALS_EXPORT int CVodeGetCurrentStep(void* cvode_mem, sunrealtype* hcur);
SUNDIALS_EXPORT int CVodeGetCurrentState(void* cvode_mem,
                                         N_Vector* y); // nb::rv_policy::reference
SUNDIALS_EXPORT int CVodeGetCurrentTime(void* cvode_mem, sunrealtype* tcur);
SUNDIALS_EXPORT int CVodeGetTolScaleFactor(void* cvode_mem, sunrealtype* tolsfac);
SUNDIALS_EXPORT int CVodeGetErrWeights(void* cvode_mem, N_Vector eweight);
SUNDIALS_EXPORT int CVodeGetEstLocalErrors(void* cvode_mem, N_Vector ele);
SUNDIALS_EXPORT int CVodeGetNumGEvals(void* cvode_mem, long int* ngevals);
SUNDIALS_EXPORT int CVodeGetRootInfo(void* cvode_mem, int* rootsfound_1d);
SUNDIALS_EXPORT int CVodeGetIntegratorStats(
  void* cvode_mem, long int* nsteps, long int* nfevals, long int* nlinsetups,
  long int* netfails, int* qlast, int* qcur, sunrealtype* hinused,
  sunrealtype* hlast, sunrealtype* hcur, sunrealtype* tcur);
SUNDIALS_EXPORT int CVodeGetNonlinearSystemData(
  void* cvode_mem, sunrealtype* tcur, N_Vector* ypred, N_Vector* yn,
  N_Vector* fn, sunrealtype* gamma, sunrealtype* rl1, N_Vector* zn1,
  void** user_data); // nb::rv_policy::reference
SUNDIALS_EXPORT int CVodeGetNumNonlinSolvIters(void* cvode_mem,
                                               long int* nniters);
SUNDIALS_EXPORT int CVodeGetNumNonlinSolvConvFails(void* cvode_mem,
                                                   long int* nnfails);
SUNDIALS_EXPORT int CVodeGetNonlinSolvStats(void* cvode_mem, long int* nniters,
                                            long int* nnfails);
SUNDIALS_EXPORT int CVodeGetNumStepSolveFails(void* cvode_mem,
                                              long int* nncfails);
SUNDIALS_EXPORT int CVodeGetNumConstraintFails(void* cvode_mem,
                                               long int* num_fails_out);
SUNDIALS_EXPORT int CVodeGetNumConstraintCorrections(void* cvode_mem,
                                                     long int* num_corrections_out);
SUNDIALS_EXPORT int CVodeGetUserData(void* cvode_mem, void** user_data);
SUNDIALS_EXPORT int CVodePrintAllStats(void* cvode_mem, FILE* outfile,
                                       SUNOutputFormat fmt);
SUNDIALS_EXPORT char* CVodeGetReturnFlagName(long int flag);

/* Free function */
SUNDIALS_EXPORT void CVodeFree(void** cvode_mem);

/* CVLS interface function that depends on CVRhsFn */
SUNDIALS_EXPORT int CVodeSetJacTimesRhsFn(void* cvode_mem, CVRhsFn jtimesRhsFn);

#ifdef __cplusplus
}
#endif

#endif
