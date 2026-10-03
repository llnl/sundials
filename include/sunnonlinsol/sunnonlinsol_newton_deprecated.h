/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * -----------------------------------------------------------------
 * Temporary header file for deprecated functions.
 * TODO:(SBR) remove in version 8.1.0
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_SUNNONLINSOL_SUNNONLINSOL_NEWTON_DEPRECATED_H
#define SUNDIALS_SUNNONLINSOL_SUNNONLINSOL_NEWTON_DEPRECATED_H

#include <sunnonlinsol/sunnonlinsol_newton.h>

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolGetType instead; will be removed in version 8.1.0")
SUNNonlinearSolver_Type SUNNonlinSolGetType_Newton(SUNNonlinearSolver NLS);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolInitialize instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolInitialize_Newton(SUNNonlinearSolver NLS);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolSolve instead; will be removed in version 8.1.0")
int SUNNonlinSolSolve_Newton(SUNNonlinearSolver NLS, N_Vector y0, N_Vector y,
                             N_Vector w, sunrealtype tol,
                             sunbooleantype callLSetup, void* mem);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolFree instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolFree_Newton(SUNNonlinearSolver NLS);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolSetSysFn instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolSetSysFn_Newton(SUNNonlinearSolver NLS,
                                       SUNNonlinSolSysFn SysFn);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolSetLSetupFn instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolSetLSetupFn_Newton(SUNNonlinearSolver NLS,
                                          SUNNonlinSolLSetupFn LSetupFn);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolSetLSolveFn instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolSetLSolveFn_Newton(SUNNonlinearSolver NLS,
                                          SUNNonlinSolLSolveFn LSolveFn);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolSetConvTestFn instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolSetConvTestFn_Newton(SUNNonlinearSolver NLS,
                                            SUNNonlinSolConvTestFn CTestFn,
                                            void* ctest_data);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolSetMaxIters instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolSetMaxIters_Newton(SUNNonlinearSolver NLS, int maxiters);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNNonlinSolSetGetUpdateNormFn instead; "
                               "will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolSetGetUpdateNormFn_Newton(
  SUNNonlinearSolver NLS, SUNNonlinSolGetUpdateNormFn GetUpdateNormFn,
  void* getupdatenorm_data);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolSetNormFn instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolSetNormFn_Newton(SUNNonlinearSolver NLS,
                                        SUNNonlinSolNormFn NormFn,
                                        void* norm_fn_data);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolGetNumIters instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolGetNumIters_Newton(SUNNonlinearSolver NLS,
                                          long int* niters);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolGetCurIter instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolGetCurIter_Newton(SUNNonlinearSolver NLS, int* iter);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolGetNumConvFails instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolGetNumConvFails_Newton(SUNNonlinearSolver NLS,
                                              long int* nconvfails);

#ifdef __cplusplus
}
#endif

#endif
