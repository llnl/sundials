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
 * Stack of temporary N_Vectors for reuse.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_VECSTACK_H
#define _SUNDIALS_VECSTACK_H

#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_EXPORT
SUNErrCode SUNVecStack_Create(N_Vector tmpl, int init_size, SUNContext sunctx,
                              SUNVecStack* stack_out);

SUNDIALS_EXPORT
SUNErrCode SUNVecStack_Destroy(SUNVecStack* stack_in);

SUNDIALS_EXPORT
SUNErrCode SUNVecStack_Pop(SUNVecStack stack, N_Vector* vec_out);

SUNDIALS_EXPORT
SUNErrCode SUNVecStack_Push(SUNVecStack stack, N_Vector* vec_in);

SUNDIALS_EXPORT
SUNErrCode SUNVecStack_GetNumVecs(SUNVecStack stack, int64_t* num_vecs);

SUNDIALS_EXPORT
SUNErrCode SUNVecStack_GetNumActiveVecs(SUNVecStack stack,
                                        int64_t* num_active_vecs);

SUNDIALS_EXPORT
SUNErrCode SUNVecStack_GetNumIdleVecs(SUNVecStack stack,
                                      int64_t* num_idle_vecs);

#ifdef __cplusplus
}
#endif

#endif /* _SUNDIALS_VECSTACK_H */
