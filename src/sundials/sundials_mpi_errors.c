/* -----------------------------------------------------------------
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
 * -----------------------------------------------------------------*/

#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/priv/sundials_mpi_errors_impl.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_mpi_types.h>

#include "sundials_macros.h"
#include "sundials_utils.h"

void SUNMPIAbortErrHandlerFn(int line, const char* func, const char* file,
                             SUNDIALS_MAYBE_UNUSED const char* msg,
                             SUNErrCode err_code,
                             SUNDIALS_MAYBE_UNUSED void* err_user_data,
                             SUNContext sunctx)
{
  char* file_and_line = sunCombineFileAndLine(line, file);
  SUNLogger_QueueMsg(sunctx->logger, SUN_LOGLEVEL_ERROR, file_and_line, func,
                     "SUNMPIAbortErrHandler: Calling MPI_Abort now, use a "
                     "different "
                     "error handler to avoid program termination.\n");
  free(file_and_line);
  MPI_Abort(sunctx->comm, err_code);
}
