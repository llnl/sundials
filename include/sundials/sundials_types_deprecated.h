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
 * This header file exports realtype and booleantype and a few macros
 * related to realtype for backwards compatibility. It is preferable
 * to only use the types defined in sundials_types.h .
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_TYPES_DEPRECATED_H
#define _SUNDIALS_TYPES_DEPRECATED_H

#include <float.h>
#include <stddef.h>
#include <stdint.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 *------------------------------------------------------------------
 * Type realtype
 * Macro RCONST
 * Constants SMALL_REAL, BIG_REAL, and UNIT_ROUNDOFF
 *------------------------------------------------------------------
 */

#if defined(SUNDIALS_SINGLE_PRECISION)

typedef float realtype;
#define RCONST(x)     x##F
#define BIG_REAL      FLT_MAX
#define SMALL_REAL    FLT_MIN
#define UNIT_ROUNDOFF FLT_EPSILON

#elif defined(SUNDIALS_DOUBLE_PRECISION)

typedef double realtype;
#define RCONST(x)     x
#define BIG_REAL      DBL_MAX
#define SMALL_REAL    DBL_MIN
#define UNIT_ROUNDOFF DBL_EPSILON

#elif defined(SUNDIALS_EXTENDED_PRECISION)

typedef long double realtype;
#define RCONST(x)     x##L
#define BIG_REAL      LDBL_MAX
#define SMALL_REAL    LDBL_MIN
#define UNIT_ROUNDOFF LDBL_EPSILON

#endif

#ifndef booleantype
#define booleantype int
#endif

#ifdef __cplusplus
}
#endif

#endif /* _SUNDIALS_TYPES_DEPRECATED_H */
