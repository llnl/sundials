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
 * This is the header file for the SUNNonlinearSolver module implementation of
 * a wrapper to the PETSc SNES nonlinear solvers.
 *
 * Part I defines the solver-specific content structure.
 *
 * Part II contains prototypes for the solver constructor and operations.
 * ---------------------------------------------------------------------------*/

#ifndef SUNDIALS_SUNNONLINSOL_PETSCSNES_H
#define SUNDIALS_SUNNONLINSOL_PETSCSNES_H

#include <petscsnes.h>

#include "sundials/sundials_nonlinearsolver.h"
#include "sundials/sundials_nvector.h"
#include "sundials/sundials_types.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------------------------------
 * I. Content structure
 * ---------------------------------------------------------------------------*/

struct SUNNonlinearSolverContent_PetscSNES_
{
  int sysfn_last_err; /* last error returned by the system function Sys */
  PetscErrorCode petsc_last_err; /* last error return by PETSc */
  long int nconvfails; /* number of nonlinear converge failures (recoverable or not) */
  long int nni;        /* number of nonlinear iterations */
  void* imem;          /* SUNDIALS integrator memory */
  SNES snes;           /* PETSc SNES context */
  Vec r;               /* nonlinear residual */
  N_Vector y, f; /* wrappers for PETSc vectors in system function */
  /* functions provided by the integrator */
  SUNNonlinSolSysFn Sys; /* nonlinear system function         */
};

typedef struct SUNNonlinearSolverContent_PetscSNES_* SUNNonlinearSolverContent_PetscSNES;

/* -----------------------------------------------------------------------------
 * II: Exported functions
 * ---------------------------------------------------------------------------*/

/* Constructor to create solver and allocates memory */

SUNDIALS_EXPORT
SUNNonlinearSolver SUNNonlinSol_PetscSNES(N_Vector y, SNES snes,
                                          SUNContext sunctx);

/* SUNNonlinearSolver API functions */

/* Implementation specific functions */

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetSNES_PetscSNES(SUNNonlinearSolver NLS, SNES* snes);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetPetscError_PetscSNES(SUNNonlinearSolver NLS,
                                               PetscErrorCode* err);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetSysFn_PetscSNES(SUNNonlinearSolver NLS,
                                          SUNNonlinSolSysFn* SysFn);

#ifdef __cplusplus
}
#endif

#endif
