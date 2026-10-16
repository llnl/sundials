/* ---------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 * ---------------------------------------------------------------------------
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

#ifndef SUNDIALS_SUNNONLINSOL_SUNNONLINSOL_FIXEDPOINT_DEPRECATED_H
#define SUNDIALS_SUNNONLINSOL_SUNNONLINSOL_FIXEDPOINT_DEPRECATED_H

#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolGetType instead; will be removed in version 8.1.0")
SUNNonlinearSolver_Type SUNNonlinSolGetType_FixedPoint(SUNNonlinearSolver NLS);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolInitialize instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolInitialize_FixedPoint(SUNNonlinearSolver NLS);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolSolve instead; will be removed in version 8.1.0")
int SUNNonlinSolSolve_FixedPoint(SUNNonlinearSolver NLS, N_Vector y0,
                                 N_Vector y, N_Vector w, sunrealtype tol,
                                 sunbooleantype callSetup, void* mem);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolFree instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolFree_FixedPoint(SUNNonlinearSolver NLS);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolSetOptions instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolSetOptions_FixedPoint(SUNNonlinearSolver NLS,
                                             const char* NLSid,
                                             const char* file_name, int argc,
                                             char* argv[]);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolSetSysFn instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolSetSysFn_FixedPoint(SUNNonlinearSolver NLS,
                                           SUNNonlinSolSysFn SysFn);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolSetConvTestFn instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolSetConvTestFn_FixedPoint(SUNNonlinearSolver NLS,
                                                SUNNonlinSolConvTestFn CTestFn,
                                                void* ctest_data);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolSetMaxIters instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolSetMaxIters_FixedPoint(SUNNonlinearSolver NLS,
                                              int maxiters);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolSetNormFn instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolSetNormFn_FixedPoint(SUNNonlinearSolver NLS,
                                            SUNNonlinSolNormFn NormFn,
                                            void* norm_fn_data);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNNonlinSolSetGetUpdateNormFn instead; "
                               "will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolSetGetUpdateNormFn_FixedPoint(
  SUNNonlinearSolver NLS, SUNNonlinSolGetUpdateNormFn GetUpdateNormFn,
  void* getupdatenorm_data);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolGetNumIters instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolGetNumIters_FixedPoint(SUNNonlinearSolver NLS,
                                              long int* niters);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolGetCurIter instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolGetCurIter_FixedPoint(SUNNonlinearSolver NLS, int* iter);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolGetNumConvFails instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolGetNumConvFails_FixedPoint(SUNNonlinearSolver NLS,
                                                  long int* nconvfails);

#ifdef __cplusplus
}
#endif

#endif
