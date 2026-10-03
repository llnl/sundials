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
 * -----------------------------------------------------------------
 * Temporary header file for deprecated functions.
 * TODO:(SBR) remove in version 8.1.0
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_SUNNONLINSOL_SUNNONLINSOL_PETSCSNES_DEPRECATED_H
#define SUNDIALS_SUNNONLINSOL_SUNNONLINSOL_PETSCSNES_DEPRECATED_H

#include <sunnonlinsol/sunnonlinsol_petscsnes.h>

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolGetType instead; will be removed in version 8.1.0")
SUNNonlinearSolver_Type SUNNonlinSolGetType_PetscSNES(SUNNonlinearSolver NLS);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolInitialize instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolInitialize_PetscSNES(SUNNonlinearSolver NLS);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolSolve instead; will be removed in version 8.1.0")
int SUNNonlinSolSolve_PetscSNES(SUNNonlinearSolver NLS, N_Vector y0, N_Vector y,
                                N_Vector w, sunrealtype tol,
                                sunbooleantype callLSetup, void* mem);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolSetSysFn instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolSetSysFn_PetscSNES(SUNNonlinearSolver NLS,
                                          SUNNonlinSolSysFn SysFn);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolGetNumIters instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolGetNumIters_PetscSNES(SUNNonlinearSolver NLS,
                                             long int* nni);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolGetNumConvFails instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolGetNumConvFails_PetscSNES(SUNNonlinearSolver NLS,
                                                 long int* nconvfails);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNNonlinSolFree instead; will be removed in version 8.1.0")
SUNErrCode SUNNonlinSolFree_PetscSNES(SUNNonlinearSolver NLS);

#ifdef __cplusplus
}
#endif

#endif
