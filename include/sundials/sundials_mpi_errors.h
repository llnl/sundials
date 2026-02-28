/* -----------------------------------------------------------------
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
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_MPI_ERRORS_H
#define _SUNDIALS_MPI_ERRORS_H

#include <sundials/sundials_errors.h>

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_EXPORT
void SUNMPIAbortErrHandlerFn(int line, const char* func, const char* file,
                             const char* msg, SUNErrCode err_code,
                             void* err_user_data, SUNContext sunctx);

#ifdef __cplusplus
}
#endif

#endif /* _SUNDIALS_MPI_ERRORS_H */
