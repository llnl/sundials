/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *                Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
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
 * This is the header file for the KINBBDPRE module, for a
 * band-block-diagonal preconditioner, i.e. a block-diagonal
 * matrix with banded blocks.
 * -----------------------------------------------------------------*/

#ifndef _KINBBDPRE_H
#define _KINBBDPRE_H

#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* KINBBDPRE return values */

#define KINBBDPRE_SUCCESS      0
#define KINBBDPRE_PDATA_NULL   -11
#define KINBBDPRE_FUNC_UNRECVR -12

/* User-supplied function Types */

typedef int (*KINBBDCommFn)(sunindextype Nlocal, N_Vector u, void* user_data);

typedef int (*KINBBDLocalFn)(sunindextype Nlocal, N_Vector uu, N_Vector gval,
                             void* user_data);

/* Exported Functions */

SUNDIALS_EXPORT int KINBBDPrecInit(void* kinmem, sunindextype Nlocal,
                                   sunindextype mudq, sunindextype mldq,
                                   sunindextype mukeep, sunindextype mlkeep,
                                   sunrealtype dq_rel_uu, KINBBDLocalFn gloc,
                                   KINBBDCommFn gcomm);

/* Optional output functions */

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
int KINBBDPrecGetWorkSpace(void* kinmem, long int* lenrwBBDP,
                           long int* leniwBBDP);

SUNDIALS_EXPORT int KINBBDPrecGetNumGfnEvals(void* kinmem, long int* ngevalsBBDP);

#ifdef __cplusplus
}
#endif

#endif
