/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
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
 * This is the header file for the ARKBANDPRE module, which provides
 * a banded difference quotient Jacobian-based preconditioner.
 * -----------------------------------------------------------------*/

#ifndef _ARKBANDPRE_H
#define _ARKBANDPRE_H

#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* BandPrec initialization function */

SUNDIALS_EXPORT int ARKBandPrecInit(void* arkode_mem, sunindextype N,
                                    sunindextype mu, sunindextype ml);

/* Optional output functions */
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
int ARKBandPrecGetWorkSpace(void* arkode_mem, long int* lenrwLS,
                            long int* leniwLS);
SUNDIALS_EXPORT int ARKBandPrecGetNumRhsEvals(void* arkode_mem,
                                              long int* nfevalsBP);

#ifdef __cplusplus
}
#endif

#endif
