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
 * SUNDIALS FILE interface utility implementations.
 * -----------------------------------------------------------------*/

#include <string.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_futils.h>

/* Create a file pointer with the given file name and mode. */
SUNErrCode SUNFileOpen(const char* filename, const char* mode, FILE** fp_out)
{
  SUNErrCode err = SUN_SUCCESS;
  FILE* fp       = *fp_out;

  if (filename)
  {
    if (!strcmp(filename, "stdout")) { fp = stdout; }
    else if (!strcmp(filename, "stderr")) { fp = stderr; }
    else { fp = fopen(filename, mode); }
  }

  if (!fp) { err = SUN_ERR_FILE_OPEN; }

  *fp_out = fp;
  return err;
}

SUNErrCode SUNDIALSFileOpen(const char* filename, const char* mode, FILE** fp_out)
{
  return SUNFileOpen(filename, mode, fp_out);
}

/* Close a file pointer with the given file name. */
SUNErrCode SUNFileClose(FILE** fp_ptr)
{
  if (!fp_ptr) { return SUN_SUCCESS; }
  FILE* fp = *fp_ptr;
  if (fp && (fp != stdout) && (fp != stderr)) { fclose(fp); }
  return SUN_SUCCESS;
}

SUNErrCode SUNDIALSFileClose(FILE** fp_ptr) { return SUNFileClose(fp_ptr); }
