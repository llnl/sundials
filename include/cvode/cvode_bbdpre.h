/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *                Michael Wittman, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
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
 * This is the header file for the CVBBDPRE module, for a
 * band-block-diagonal preconditioner, i.e. a block-diagonal
 * matrix with banded blocks.
 * -----------------------------------------------------------------*/

#ifndef _CVBBDPRE_H
#define _CVBBDPRE_H

#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* User-supplied function Types */

typedef int (*CVLocalFn)(sunindextype Nlocal, sunrealtype t, N_Vector y,
                         N_Vector g, void* user_data);

typedef int (*CVCommFn)(sunindextype Nlocal, sunrealtype t, N_Vector y,
                        void* user_data);

/* Exported Functions */

SUNDIALS_EXPORT int CVBBDPrecInit(void* cvode_mem, sunindextype Nlocal,
                                  sunindextype mudq, sunindextype mldq,
                                  sunindextype mukeep, sunindextype mlkeep,
                                  sunrealtype dqrely, CVLocalFn gloc,
                                  CVCommFn cfn);

SUNDIALS_EXPORT int CVBBDPrecReInit(void* cvode_mem, sunindextype mudq,
                                    sunindextype mldq, sunrealtype dqrely);

/* Optional output functions */

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
int CVBBDPrecGetWorkSpace(void* cvode_mem, long int* lenrwBBDP,
                          long int* leniwBBDP);

SUNDIALS_EXPORT int CVBBDPrecGetNumGfnEvals(void* cvode_mem,
                                            long int* ngevalsBBDP);

#ifdef __cplusplus
}
#endif

#endif
