/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * This file implements functions for getting SUNDIALS version
 * information.
 * -----------------------------------------------------------------*/

#include <string.h>
#include <sundials/sundials_version.h>

#include "sundials/sundials_errors.h"

/* note strlen does not include terminating null character hence the
   use of >= when checking len below and strncpy copies up to len
   characters including the terminating null character */

/* fill string with SUNDIALS version information */
SUNErrCode SUNDIALSGetVersion(char* version, int len)
{
  if (version == NULL) { return SUN_ERR_ARG_CORRUPT; }
  if (strlen(SUNDIALS_VERSION) >= (size_t)len)
  {
    return SUN_ERR_ARG_OUTOFRANGE;
  }

  strncpy(version, SUNDIALS_VERSION, (size_t)len);

  return SUN_SUCCESS;
}

/* fill integers with SUNDIALS major, minor, and patch release
   numbers and fill a string with the release label */
SUNErrCode SUNDIALSGetVersionNumber(int* major, int* minor, int* patch,
                                    char* label, int len)
{
  if (major == NULL || minor == NULL || patch == NULL || label == NULL)
  {
    return SUN_ERR_ARG_CORRUPT;
  }
  if (strlen(SUNDIALS_VERSION_LABEL) >= (size_t)len)
  {
    return SUN_ERR_ARG_OUTOFRANGE;
  }

  *major = SUNDIALS_VERSION_MAJOR;
  *minor = SUNDIALS_VERSION_MINOR;
  *patch = SUNDIALS_VERSION_PATCH;
  strncpy(label, SUNDIALS_VERSION_LABEL, (size_t)len);

  return SUN_SUCCESS;
}
