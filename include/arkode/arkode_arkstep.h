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
 * This is the header file for the ARKODE ARKStep module.
 * -----------------------------------------------------------------*/

#ifndef _ARKSTEP_H
#define _ARKSTEP_H

#include <arkode/arkode.h>
#include <arkode/arkode_arkstep_deprecated.h>
#include <arkode/arkode_butcher_dirk.h>
#include <arkode/arkode_butcher_erk.h>
#include <arkode/arkode_ls.h>
#include <sunadaptcontroller/sunadaptcontroller_imexgus.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>
#include <sundials/sundials_adjointstepper.h>
#include <sundials/sundials_stepper.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------
 * ARKStep Constants
 * ----------------- */

/* Default Butcher tables for each method/order */

/* Ideally these defaults would be declared with types ARKODE_ERKTableID and
 * ARKODE_DIRKTableID, but this causes swig to unnecessarily append `C_INT` to
 * the variable names */

/*    explicit */
static const int ARKSTEP_DEFAULT_ERK_1 = ARKODE_FORWARD_EULER_1_1;
static const int ARKSTEP_DEFAULT_ERK_2 = ARKODE_RALSTON_3_1_2;
static const int ARKSTEP_DEFAULT_ERK_3 = ARKODE_BOGACKI_SHAMPINE_4_2_3;
static const int ARKSTEP_DEFAULT_ERK_4 = ARKODE_SOFRONIOU_SPALETTA_5_3_4;
static const int ARKSTEP_DEFAULT_ERK_5 = ARKODE_TSITOURAS_7_4_5;
static const int ARKSTEP_DEFAULT_ERK_6 = ARKODE_VERNER_9_5_6;
static const int ARKSTEP_DEFAULT_ERK_7 = ARKODE_VERNER_10_6_7;
static const int ARKSTEP_DEFAULT_ERK_8 = ARKODE_VERNER_13_7_8;
static const int ARKSTEP_DEFAULT_ERK_9 = ARKODE_VERNER_16_8_9;

/*    implicit */
static const int ARKSTEP_DEFAULT_DIRK_1 = ARKODE_BACKWARD_EULER_1_1;
static const int ARKSTEP_DEFAULT_DIRK_2 = ARKODE_ARK2_DIRK_3_1_2;
static const int ARKSTEP_DEFAULT_DIRK_3 = ARKODE_ESDIRK325L2SA_5_2_3;
static const int ARKSTEP_DEFAULT_DIRK_4 = ARKODE_ESDIRK436L2SA_6_3_4;
static const int ARKSTEP_DEFAULT_DIRK_5 = ARKODE_ESDIRK547L2SA2_7_4_5;

/*    ImEx */
static const int ARKSTEP_DEFAULT_ARK_ETABLE_2 = ARKODE_ARK2_ERK_3_1_2;
static const int ARKSTEP_DEFAULT_ARK_ETABLE_3 = ARKODE_ARK324L2SA_ERK_4_2_3;
static const int ARKSTEP_DEFAULT_ARK_ETABLE_4 = ARKODE_ARK437L2SA_ERK_7_3_4;
static const int ARKSTEP_DEFAULT_ARK_ETABLE_5 = ARKODE_ARK548L2SAb_ERK_8_4_5;
static const int ARKSTEP_DEFAULT_ARK_ITABLE_2 = ARKODE_ARK2_DIRK_3_1_2;
static const int ARKSTEP_DEFAULT_ARK_ITABLE_3 = ARKODE_ARK324L2SA_DIRK_4_2_3;
static const int ARKSTEP_DEFAULT_ARK_ITABLE_4 = ARKODE_ARK437L2SA_DIRK_7_3_4;
static const int ARKSTEP_DEFAULT_ARK_ITABLE_5 = ARKODE_ARK548L2SAb_DIRK_8_4_5;

/* -------------------
 * Exported Functions
 * ------------------- */

/* Creation and Reinitialization functions */
SUNDIALS_EXPORT void* ARKStepCreate(ARKRhsFn fe, ARKRhsFn fi, sunrealtype t0,
                                    N_Vector y0, SUNContext sunctx);
SUNDIALS_EXPORT int ARKStepReInit(void* arkode_mem, ARKRhsFn fe, ARKRhsFn fi,
                                  sunrealtype t0, N_Vector y0);

/* Optional input functions -- must be called AFTER ARKStepCreate */
SUNDIALS_EXPORT int ARKStepSetExplicit(void* arkode_mem);
SUNDIALS_EXPORT int ARKStepSetImplicit(void* arkode_mem);
SUNDIALS_EXPORT int ARKStepSetImEx(void* arkode_mem);
SUNDIALS_EXPORT int ARKStepSetTables(void* arkode_mem, int q, int p,
                                     ARKodeButcherTable Bi,
                                     ARKodeButcherTable Be);
SUNDIALS_EXPORT int ARKStepSetTableNum(void* arkode_mem,
                                       ARKODE_DIRKTableID itable,
                                       ARKODE_ERKTableID etable);
SUNDIALS_EXPORT int ARKStepSetTableName(void* arkode_mem, const char* itable,
                                        const char* etable);

/* Optional output functions */
SUNDIALS_EXPORT int ARKStepGetCurrentButcherTables(
  void* arkode_mem, ARKodeButcherTable* Bi,
  ARKodeButcherTable* Be); // nb::rv_policy::reference
SUNDIALS_EXPORT int ARKStepGetTimestepperStats(
  void* arkode_mem, long int* expsteps, long int* accsteps,
  long int* step_attempts, long int* nfe_evals, long int* nfi_evals,
  long int* nlinsetups, long int* netfails);

/* Adjoint solver functions */
SUNDIALS_EXPORT
int ARKStepCreateAdjointStepper(void* arkode_mem, SUNAdjRhsFn adj_fe,
                                SUNAdjRhsFn adj_fi, sunrealtype tf, N_Vector sf,
                                SUNContext sunctx,
                                SUNAdjointStepper* adj_stepper_ptr);

#ifdef __cplusplus
}
#endif

#endif
