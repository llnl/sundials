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
 * This header file contains definitions of MPI data types, which
 * are used by MPI parallel vector implementations.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_MPI_TYPES_H
#define _SUNDIALS_MPI_TYPES_H

#include <mpi.h>
#include <sundials/sundials_types.h>

/* define MPI data types */

#if defined(SUNDIALS_SINGLE_PRECISION)
#define MPI_SUNREALTYPE MPI_FLOAT
#elif defined(SUNDIALS_DOUBLE_PRECISION)
#define MPI_SUNREALTYPE MPI_DOUBLE
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define MPI_SUNREALTYPE MPI_LONG_DOUBLE
#endif

#if defined(SUNDIALS_INT64_T)
#define MPI_SUNINDEXTYPE MPI_INT64_T
#elif defined(SUNDIALS_INT32_T)
#define MPI_SUNINDEXTYPE MPI_INT32_T
#endif

#endif /* _SUNDIALS_MPI_TYPES_H */
