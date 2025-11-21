/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
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
 * SUNDIALS FILE utility definitions.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_FUTILS_H
#define _SUNDIALS_FUTILS_H

#include <stdio.h>

#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Create a file pointer with the given file name and mode. */
SUNDIALS_EXPORT
SUNErrCode SUNFileOpen(const char* filename, const char* modes, FILE** fp);

SUNDIALS_DEPRECATED_EXPORT_MSG("Use SUNFileOpen")
SUNErrCode SUNDIALSFileOpen(const char* filename, const char* modes, FILE** fp);

/* Close a file pointer with the given file name. */
SUNDIALS_EXPORT
SUNErrCode SUNFileClose(FILE** fp);

SUNDIALS_DEPRECATED_EXPORT_MSG("Use SUNFileClose")
SUNErrCode SUNDIALSFileClose(FILE** fp);

#ifdef __cplusplus
}
#endif

#endif
