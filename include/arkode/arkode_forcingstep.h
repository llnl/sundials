/*---------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 *---------------------------------------------------------------
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
 *---------------------------------------------------------------
 * This is the header file for the ARKODE ForcingStep module.
 *--------------------------------------------------------------*/

#ifndef ARKODE_FORCINGINGSTEP_H_
#define ARKODE_FORCINGINGSTEP_H_

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_stepper.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

SUNDIALS_EXPORT void* ForcingStepCreate(SUNStepper stepper1,
                                        SUNStepper stepper2, sunrealtype t0,
                                        N_Vector y0, SUNContext sunctx);

SUNDIALS_EXPORT int ForcingStepReInit(void* arkode_mem, SUNStepper stepper1,
                                      SUNStepper stepper2, sunrealtype t0,
                                      N_Vector y0);

SUNDIALS_EXPORT int ForcingStepGetNumEvolves(void* arkode_mem, int partition,
                                             long int* evolves);

#ifdef __cplusplus
}
#endif

#endif
