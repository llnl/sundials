/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ UMBC
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
 * This is the header file for the LAPACK band implementation of the
 * SUNLINSOL module, SUNLINSOL_LAPACKBAND.
 *
 * Note:
 *   - The definition of the generic SUNLinearSolver structure can
 *     be found in the header file sundials_linearsolver.h.
 * -----------------------------------------------------------------
 */

#ifndef SUNDIALS_SUNLINSOL_LAPACKBAND_H
#define SUNDIALS_SUNLINSOL_LAPACKBAND_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sunmatrix/sunmatrix_band.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ----------------------------------------------
 * LAPACK band implementation of SUNLinearSolver
 * ---------------------------------------------- */

struct _SUNLinearSolverContent_LapackBand
{
  sunindextype N;
  sunindextype* pivots;
  sunindextype last_flag;
};

typedef struct _SUNLinearSolverContent_LapackBand* SUNLinearSolverContent_LapackBand;

/* --------------------------------------------
 * Exported Functions for SUNLINSOL_LAPACKBAND
 * -------------------------------------------- */

SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_LapackBand(N_Vector y, SUNMatrix A,
                                                     SUNContext sunctx);

#ifdef __cplusplus
}
#endif

#endif
