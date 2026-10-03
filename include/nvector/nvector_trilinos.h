/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
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
 * This is the main header file for the Trilinos vector wrapper
 * for NVECTOR module.
 *
 * Part I contains declarations specific to the Trilinos vector wrapper
 * implementation.
 *
 * Part II contains the prototype for the constructor
 * N_VMake_Trilinos as well as Trilinos-specific prototypes
 * for various useful vector operations.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be
 *     found in the header file sundials_nvector.h.
 *
 *   - The definition of the type sunrealtype can be found in the
 *     header file sundials_types.h, and it may be changed (at the
 *     build configuration stage) according to the user's needs.
 *     The sundials_types.h file also contains the definition
 *     for the type sunbooleantype.
 *
 *   - N_Vector arguments to arithmetic vector operations need not
 *     be distinct. For example, the following call:
 *
 *        N_VLinearSum(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_NVECTOR_TRILINOS_H
#define SUNDIALS_NVECTOR_TRILINOS_H

#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * PART I: N_Vector interface to Trilinos vector
 * -----------------------------------------------------------------
 */

/*
 * Dummy _N_VectorContent_Trilinos structure is used for
 * interfacing C with C++ code
 */

struct _N_VectorContent_Trilinos
{};

typedef struct _N_VectorContent_Trilinos* N_VectorContent_Trilinos;

/*
 * -----------------------------------------------------------------
 * PART II: functions exported by nvector_Trilinos
 *
 * CONSTRUCTORS:
 *    N_VNewEmpty_Trilinos
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : N_VNewEmpty_Trilinos
 * -----------------------------------------------------------------
 * This function creates a new N_Vector wrapper for a Trilinos
 * vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_Trilinos(SUNContext sunctx);

/*
 * -----------------------------------------------------------------
 * Trilinos implementations of the vector operations
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus
}
#endif

#endif /* SUNDIALS_NVECTOR_TRILINOS_H */
