/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * This is the header file for the ARKODE SPRKStep module.
 * -----------------------------------------------------------------*/

#ifndef _ARKODE_SPRKSTEP_H
#define _ARKODE_SPRKSTEP_H

#include <arkode/arkode.h>
#include <arkode/arkode_sprk.h>
#include <arkode/arkode_sprkstep_deprecated.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------
 * SPRKStep Constants
 * ----------------- */

static const int SPRKSTEP_DEFAULT_1  = ARKODE_SPRK_EULER_1_1;
static const int SPRKSTEP_DEFAULT_2  = ARKODE_SPRK_LEAPFROG_2_2;
static const int SPRKSTEP_DEFAULT_3  = ARKODE_SPRK_MCLACHLAN_3_3;
static const int SPRKSTEP_DEFAULT_4  = ARKODE_SPRK_MCLACHLAN_4_4;
static const int SPRKSTEP_DEFAULT_5  = ARKODE_SPRK_MCLACHLAN_5_6;
static const int SPRKSTEP_DEFAULT_6  = ARKODE_SPRK_YOSHIDA_6_8;
static const int SPRKSTEP_DEFAULT_8  = ARKODE_SPRK_SUZUKI_UMENO_8_16;
static const int SPRKSTEP_DEFAULT_10 = ARKODE_SPRK_SOFRONIOU_10_36;

/* -------------------
 * Exported Functions
 * ------------------- */

/* Creation and Reinitialization functions */

SUNDIALS_EXPORT void* SPRKStepCreate(ARKRhsFn f1, ARKRhsFn f2, sunrealtype t0,
                                     N_Vector y0, SUNContext sunctx);
SUNDIALS_EXPORT int SPRKStepReInit(void* arkode_mem, ARKRhsFn f1, ARKRhsFn f2,
                                   sunrealtype t0, N_Vector y0);

/* Optional input functions -- must be called AFTER SPRKStepCreate */

SUNDIALS_EXPORT int SPRKStepSetMethod(void* arkode_mem,
                                      ARKodeSPRKTable sprk_storage);
SUNDIALS_EXPORT int SPRKStepSetMethodName(void* arkode_mem, const char* method);

/* Optional output functions */

SUNDIALS_EXPORT int SPRKStepGetCurrentMethod(
  void* arkode_mem,
  ARKodeSPRKTable* sprk_storage); // nb::rv_policy::reference

#ifdef __cplusplus
}
#endif

#endif
