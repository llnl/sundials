/* -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier, Shelby Lockhart @ LLNL
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
 * This is the header file for the main KINSOL solver.
 * -----------------------------------------------------------------*/

#ifndef _KINSOL_H
#define _KINSOL_H

#include <kinsol/kinsol_ls.h>
#include <stdio.h>
#include <sundials/sundials_core.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------
 * KINSOL Constants
 * ----------------- */

/*  return values */

#define KIN_SUCCESS          0
#define KIN_INITIAL_GUESS_OK 1
#define KIN_STEP_LT_STPTOL   2

#define KIN_WARNING 99

#define KIN_MEM_NULL            -1
#define KIN_ILL_INPUT           -2
#define KIN_NO_MALLOC           -3
#define KIN_MEM_FAIL            -4
#define KIN_LINESEARCH_NONCONV  -5
#define KIN_MAXITER_REACHED     -6
#define KIN_MXNEWT_5X_EXCEEDED  -7
#define KIN_LINESEARCH_BCFAIL   -8
#define KIN_LINSOLV_NO_RECOVERY -9
#define KIN_LINIT_FAIL          -10
#define KIN_LSETUP_FAIL         -11
#define KIN_LSOLVE_FAIL         -12
#define KIN_SYSFUNC_FAIL        -13
#define KIN_FIRST_SYSFUNC_ERR   -14
#define KIN_REPTD_SYSFUNC_ERR   -15
#define KIN_VECTOROP_ERR        -16
#define KIN_CONTEXT_ERR         -17
#define KIN_DAMPING_FN_ERR      -18
#define KIN_DEPTH_FN_ERR        -19

/* Anderson Acceleration Orthogonalization Choice */
#define KIN_ORTH_MGS   0
#define KIN_ORTH_ICWY  1
#define KIN_ORTH_CGS2  2
#define KIN_ORTH_DCGS2 3

/* Enumeration for eta choice */
#define KIN_ETACHOICE1  1
#define KIN_ETACHOICE2  2
#define KIN_ETACONSTANT 3

/* Enumeration for global strategy */
#define KIN_NONE       0
#define KIN_LINESEARCH 1
#define KIN_PICARD     2
#define KIN_FP         3

/* ------------------------------
 * User-Supplied Function Types
 * ------------------------------ */

typedef int (*KINSysFn)(N_Vector uu, N_Vector fval, void* user_data);

typedef int (*KINDampingFn)(long int iter, N_Vector u_val, N_Vector g_val,
                            sunrealtype* qt_fn_1d, long int depth,
                            void* user_data, sunrealtype* damping_factor_ptr);

typedef int (*KINDepthFn)(long int iter, N_Vector u_val, N_Vector g_val,
                          N_Vector f_val, N_Vector* df_1d,
                          sunrealtype* R_mat_1d, long int depth,
                          void* user_data, long int* new_depth_ptr,
                          sunbooleantype* remove_indices_1d);

/* -------------------
 * Exported Functions
 * ------------------- */

/* Creation function */
SUNDIALS_EXPORT void* KINCreate(SUNContext sunctx);

/* Initialization function */
SUNDIALS_EXPORT int KINInit(void* kinmem, KINSysFn func, N_Vector tmpl);

/* Solver function */
SUNDIALS_EXPORT int KINSol(void* kinmem, N_Vector uu, int strategy,
                           N_Vector u_scale, N_Vector f_scale);

/* Optional input functions */
SUNDIALS_EXPORT int KINSetOptions(void* kinmem, const char* kinid,
                                  const char* file_name, int argc, char* argv[]);
SUNDIALS_EXPORT int KINSetUserData(void* kinmem, void* user_data);
SUNDIALS_EXPORT int KINSetDamping(void* kinmem, sunrealtype beta);
SUNDIALS_EXPORT int KINSetMAA(void* kinmem, long int maa);
SUNDIALS_EXPORT int KINSetOrthAA(void* kinmem, int orthaa);
SUNDIALS_EXPORT int KINSetDelayAA(void* kinmem, long int delay);
SUNDIALS_EXPORT int KINSetDampingAA(void* kinmem, sunrealtype beta);
SUNDIALS_EXPORT int KINSetDampingFn(void* kinmem, KINDampingFn damping_fn);
SUNDIALS_EXPORT int KINSetDepthFn(void* kinmem, KINDepthFn depth_fn);
SUNDIALS_EXPORT int KINSetReturnNewest(void* kinmem, sunbooleantype ret_newest);
SUNDIALS_EXPORT int KINSetNumMaxIters(void* kinmem, long int mxiter);
SUNDIALS_EXPORT int KINSetNoInitSetup(void* kinmem, sunbooleantype noInitSetup);
SUNDIALS_EXPORT int KINSetNoResMon(void* kinmem, sunbooleantype noNNIResMon);
SUNDIALS_EXPORT int KINSetMaxSetupCalls(void* kinmem, long int msbset);
SUNDIALS_EXPORT int KINSetMaxSubSetupCalls(void* kinmem, long int msbsetsub);
SUNDIALS_EXPORT int KINSetEtaForm(void* kinmem, int etachoice);
SUNDIALS_EXPORT int KINSetEtaConstValue(void* kinmem, sunrealtype eta);
SUNDIALS_EXPORT int KINSetEtaParams(void* kinmem, sunrealtype egamma,
                                    sunrealtype ealpha);
SUNDIALS_EXPORT int KINSetResMonParams(void* kinmem, sunrealtype omegamin,
                                       sunrealtype omegamax);
SUNDIALS_EXPORT int KINSetResMonConstValue(void* kinmem, sunrealtype omegaconst);
SUNDIALS_EXPORT int KINSetNoMinEps(void* kinmem, sunbooleantype noMinEps);
SUNDIALS_EXPORT int KINSetMaxNewtonStep(void* kinmem, sunrealtype mxnewtstep);
SUNDIALS_EXPORT int KINSetMaxBetaFails(void* kinmem, long int mxnbcf);
SUNDIALS_EXPORT int KINSetRelErrFunc(void* kinmem, sunrealtype relfunc);
SUNDIALS_EXPORT int KINSetFuncNormTol(void* kinmem, sunrealtype fnormtol);
SUNDIALS_EXPORT int KINSetScaledStepTol(void* kinmem, sunrealtype scsteptol);
SUNDIALS_EXPORT int KINSetConstraints(void* kinmem, N_Vector constraints);
SUNDIALS_EXPORT int KINSetSysFunc(void* kinmem, KINSysFn func);

/* Optional output functions */
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
int KINGetWorkSpace(void* kinmem, long int* lenrw, long int* leniw);
SUNDIALS_EXPORT int KINGetNumNonlinSolvIters(void* kinmem, long int* nniters);
SUNDIALS_EXPORT int KINGetNumFuncEvals(void* kinmem, long int* nfevals);
SUNDIALS_EXPORT int KINGetNumBetaCondFails(void* kinmem, long int* nbcfails);
SUNDIALS_EXPORT int KINGetNumBacktrackOps(void* kinmem, long int* nbacktr);
SUNDIALS_EXPORT int KINGetFuncNorm(void* kinmem, sunrealtype* fnorm);
SUNDIALS_EXPORT int KINGetStepLength(void* kinmem, sunrealtype* steplength);
SUNDIALS_EXPORT int KINGetUserData(void* kinmem, void** user_data);
SUNDIALS_EXPORT int KINPrintAllStats(void* kinmem, FILE* outfile,
                                     SUNOutputFormat fmt);
SUNDIALS_EXPORT char* KINGetReturnFlagName(long int flag);

/* Free function */
SUNDIALS_EXPORT void KINFree(void** kinmem);

/* KINLS interface function that depends on KINSysFn */
SUNDIALS_EXPORT int KINSetJacTimesVecSysFn(void* kinmem, KINSysFn jtimesSysFn);

#ifdef __cplusplus
}
#endif

#endif
