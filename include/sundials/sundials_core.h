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
 * -----------------------------------------------------------------
 * Header file that includes the SUNDIALS core.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_CORE_H
#define _SUNDIALS_CORE_H

#include <sundials/sundials_adaptcontroller.h>
#include <sundials/sundials_adjointcheckpointscheme.h>
#include <sundials/sundials_adjointstepper.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_domeigestimator.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_logger.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_memory.h>
#include <sundials/sundials_nonlinearsolver.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_profiler.h>
#include <sundials/sundials_stepper.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_version.h>

#if SUNDIALS_MPI_ENABLED
#include <sundials/sundials_mpi_errors.h>
#endif

#endif /* _SUNDIALS_CORE_H */
