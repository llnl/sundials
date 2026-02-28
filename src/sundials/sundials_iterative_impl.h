/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner and Shelby Lockhart @ LLNL
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
 * This is the implementation header file for SUNDIALS functions used by
 * different iterative solvers.
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_ITERATIVE_IMPL_H
#define _SUNDIALS_ITERATIVE_IMPL_H

#include <sundials/sundials_iterative.h>

/* -----------------------------------------------------------------------------
 * Type: SUNQRData
 * -----------------------------------------------------------------------------
 * A SUNQRData struct holds temporary workspace vectors and sunrealtype arrays for
 * a SUNQRAddFn. The N_Vectors and sunrealtype arrays it contains are created by
 * the routine calling a SUNQRAdd function.
 * ---------------------------------------------------------------------------*/

typedef struct _SUNQRData* SUNQRData;

struct _SUNQRData
{
  N_Vector vtemp;
  N_Vector vtemp2;
  sunrealtype* temp_array;
};

#endif
