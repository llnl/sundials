/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
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
 * This is the header file for the ARKODE ERKStep module.
 * -----------------------------------------------------------------*/

#ifndef _ERKSTEP_H
#define _ERKSTEP_H

#include <arkode/arkode.h>
#include <arkode/arkode_butcher_erk.h>
#include <arkode/arkode_erkstep_deprecated.h>
#include <sunadaptcontroller/sunadaptcontroller_imexgus.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>
#include <sundials/sundials_adjointstepper.h>
#include <sundials/sundials_stepper.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------
 * ERKStep Constants
 * ----------------- */

/* Default Butcher tables for each order */

static const int ERKSTEP_DEFAULT_1 = ARKODE_FORWARD_EULER_1_1;
static const int ERKSTEP_DEFAULT_2 = ARKODE_RALSTON_3_1_2;
static const int ERKSTEP_DEFAULT_3 = ARKODE_BOGACKI_SHAMPINE_4_2_3;
static const int ERKSTEP_DEFAULT_4 = ARKODE_SOFRONIOU_SPALETTA_5_3_4;
static const int ERKSTEP_DEFAULT_5 = ARKODE_TSITOURAS_7_4_5;
static const int ERKSTEP_DEFAULT_6 = ARKODE_VERNER_9_5_6;
static const int ERKSTEP_DEFAULT_7 = ARKODE_VERNER_10_6_7;
static const int ERKSTEP_DEFAULT_8 = ARKODE_VERNER_13_7_8;
static const int ERKSTEP_DEFAULT_9 = ARKODE_VERNER_16_8_9;

/* -------------------
 * Exported Functions
 * ------------------- */

/* Creation and Reinitialization functions */
SUNDIALS_EXPORT void* ERKStepCreate(ARKRhsFn f, sunrealtype t0, N_Vector y0,
                                    SUNContext sunctx);
SUNDIALS_EXPORT int ERKStepReInit(void* arkode_mem, ARKRhsFn f, sunrealtype t0,
                                  N_Vector y0);

/* Optional input functions -- must be called AFTER ERKStepCreate */
SUNDIALS_EXPORT int ERKStepSetTable(void* arkode_mem, ARKodeButcherTable B);
SUNDIALS_EXPORT int ERKStepSetTableNum(void* arkode_mem,
                                       ARKODE_ERKTableID etable);
SUNDIALS_EXPORT int ERKStepSetTableName(void* arkode_mem, const char* etable);

/* Optional output functions */
SUNDIALS_EXPORT int ERKStepGetCurrentButcherTable(
  void* arkode_mem,
  ARKodeButcherTable* B); // nb::rv_policy::reference

/* Grouped optional output functions */
SUNDIALS_EXPORT int ERKStepGetTimestepperStats(
  void* arkode_mem, long int* expsteps, long int* accsteps,
  long int* step_attempts, long int* nfevals, long int* netfails);

/* Adjoint solver functions */
SUNDIALS_EXPORT
int ERKStepCreateAdjointStepper(void* arkode_mem, SUNAdjRhsFn adj_f,
                                sunrealtype tf, N_Vector sf, SUNContext sunctx,
                                SUNAdjointStepper* adj_stepper_ptr);

#ifdef __cplusplus
}
#endif

#endif
