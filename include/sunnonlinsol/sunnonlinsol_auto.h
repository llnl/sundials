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
 * -----------------------------------------------------------------------------
 * This is the header file for the SUNNonlinearSolver module implementation that
 * automatically switches between Newton and fixed-point iterations.
 * ---------------------------------------------------------------------------*/

#ifndef SUNDIALS_NONLINSOLAUTO_H_
#define SUNDIALS_NONLINSOLAUTO_H_

#include <sundials/sundials_context.h>
#include <sundials/sundials_nonlinearsolver.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus
extern "C" {
#endif

enum SUNNonlinSolAutoType
{
  SUNNONLINSOL_AUTO_FIXEDPOINT = 0,
  SUNNONLINSOL_AUTO_NEWTON     = 1
};

#ifndef SWIG
typedef enum SUNNonlinSolAutoType SUNNonlinSolAutoType;
#endif

struct SUNNonlinearSolverContent_Auto_
{
  SUNNonlinSolAutoType active_solver_type;
  SUNNonlinearSolver fp_solver;
  SUNNonlinearSolver newton_solver;
  SUNNonlinSolGetConvRateFn getconvrate_fn;
  void* getconvrate_data;
  long int fp_to_newt_delay;
  long int newt_to_fp_delay;
  long int num_solves_since_switch;
  sunrealtype newt_to_fp_threshold;
  sunrealtype fp_to_newt_threshold;
  long int num_iters;
  long int num_conv_fails;
  long int switch_count;
  long int fp_num_iters_total;
  long int newton_num_iters_total;
  long int fp_num_conv_fails_total;
  long int newton_num_conv_fails_total;
  void* auto_ctest_data;
};

typedef struct SUNNonlinearSolverContent_Auto_* SUNNonlinearSolverContent_Auto;

SUNDIALS_EXPORT
SUNNonlinearSolver SUNNonlinSol_Auto(N_Vector y, int m,
                                     SUNNonlinSolAutoType initial_solver_type,
                                     SUNContext sunctx);

SUNDIALS_EXPORT
SUNNonlinearSolver_Type SUNNonlinSolGetType_Auto(SUNNonlinearSolver NLS);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolInitialize_Auto(SUNNonlinearSolver NLS);

SUNDIALS_EXPORT
int SUNNonlinSolSolve_Auto(SUNNonlinearSolver NLS, N_Vector y0, N_Vector ycor,
                           N_Vector w, sunrealtype tol,
                           sunbooleantype callSetup, void* mem);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolFree_Auto(SUNNonlinearSolver NLS);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetSysFns_Auto(SUNNonlinearSolver NLS,
                                      SUNNonlinSolSysFn root_fn,
                                      SUNNonlinSolSysFn fixed_point_fn);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetConvTestFn_Auto(SUNNonlinearSolver NLS,
                                          SUNNonlinSolConvTestFn CTestFn,
                                          void* ctest_data);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetLSetupFn_Auto(SUNNonlinearSolver NLS,
                                        SUNNonlinSolLSetupFn LSetupFn);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetLSolveFn_Auto(SUNNonlinearSolver NLS,
                                        SUNNonlinSolLSolveFn LSolveFn);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetOptions_Auto(SUNNonlinearSolver NLS,
                                       const char* NLSid, const char* file_name,
                                       int argc, char* argv[]);
SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetNormFn_Auto(SUNNonlinearSolver NLS,
                                      SUNNonlinSolNormFn NormFn,
                                      void* norm_fn_data);
SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetMaxIters_Auto(SUNNonlinearSolver NLS, int maxiters);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetGetConvRateFn_Auto(SUNNonlinearSolver NLS,
                                             SUNNonlinSolGetConvRateFn GetConvRateFn,
                                             void* getconvrate_data);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetGetUpdateNormFn_Auto(
  SUNNonlinearSolver NLS, SUNNonlinSolGetUpdateNormFn GetUpdateNormFn,
  void* getupdatenorm_data);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetSwitchingParameters_Auto(
  SUNNonlinearSolver NLS, sunrealtype newt_to_fp_threshold,
  long int newt_to_fp_delay, sunrealtype fp_to_newt_threshold,
  long int fp_to_newt_delay);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetFixedPointSolver_Auto(
  SUNNonlinearSolver NLS,
  SUNNonlinearSolver* fp_nls); // nb::rv_policy::reference

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetNewtonSolver_Auto(
  SUNNonlinearSolver NLS,
  SUNNonlinearSolver* newton_nls); // nb::rv_policy::reference

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetActiveSolverType_Auto(
  SUNNonlinearSolver NLS, SUNNonlinSolAutoType* active_solver_type);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetSwitchCount_Auto(SUNNonlinearSolver NLS,
                                           long int* switch_count);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetNumIters_Auto(SUNNonlinearSolver NLS, long int* niters);

/* Get the iteration counts for each sub-solver separately. */
SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetTotalNumItersByType_Auto(SUNNonlinearSolver NLS,
                                                   long int* fp_iters,
                                                   long int* newt_iters);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetCurIter_Auto(SUNNonlinearSolver NLS, int* iter);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetNumConvFails_Auto(SUNNonlinearSolver NLS,
                                            long int* nconvfails);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetTotalNumConvFailsByType_Auto(SUNNonlinearSolver NLS,
                                                       long int* fp_nconvfails,
                                                       long int* newt_nconvfails);

#ifdef __cplusplus
}
#endif

#endif /* SUNDIALS_NONLINSOLAUTO_H_ */
