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
 * This header file defines the ARKodeSPRKTable structure.
 * -----------------------------------------------------------------*/

#ifndef _ARKODE_SPRKTABLE_H
#define _ARKODE_SPRKTABLE_H

#include <arkode/arkode_butcher.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

enum ARKODE_SPRKMethodID
{
  ARKODE_SPRK_NONE = -1, /* ensure enum is signed int */
  /* WARNING:  ARKODE_MIN_SPRK_NUM must come after the first entry, ARKODE_SPRK_EULER_1_1,
     because Python enums will only expose the member that is defined first. Due to
     this and how pybind/nanobind handle the enums, if we defined ARKODE_MRI_NUM first,
     then ARKODE_SPRK_EULER_1_1 would not be usable from the module scope (the MIN/MAX) entries
     will still be usable when accessing through the IntEnum object, but not from module scope. */
  ARKODE_SPRK_EULER_1_1 = 0,
  ARKODE_MIN_SPRK_NUM   = 0,
  ARKODE_SPRK_LEAPFROG_2_2,
  ARKODE_SPRK_PSEUDO_LEAPFROG_2_2,
  ARKODE_SPRK_RUTH_3_3,
  ARKODE_SPRK_MCLACHLAN_2_2,
  ARKODE_SPRK_MCLACHLAN_3_3,
  ARKODE_SPRK_CANDY_ROZMUS_4_4,
  ARKODE_SPRK_MCLACHLAN_4_4,
  ARKODE_SPRK_MCLACHLAN_5_6,
  ARKODE_SPRK_YOSHIDA_6_8,
  ARKODE_SPRK_SUZUKI_UMENO_8_16,
  ARKODE_SPRK_SOFRONIOU_10_36,
  ARKODE_MAX_SPRK_NUM = ARKODE_SPRK_SOFRONIOU_10_36
};

#ifndef SWIG
typedef enum ARKODE_SPRKMethodID ARKODE_SPRKMethodID;
#endif

struct ARKodeSPRKTableMem
{
  /* method order of accuracy */
  int q;
  /* number of stages */
  int stages;
  /* the a_i coefficients generate the explicit Butcher table */
  sunrealtype* a;
  /* the ahat_i coefficients generate the diagonally-implicit Butcher table */
  sunrealtype* ahat;
};

typedef _SUNDIALS_STRUCT_ ARKodeSPRKTableMem* ARKodeSPRKTable;

/* Utility routines to allocate/free/output SPRK structures */
SUNDIALS_EXPORT
ARKodeSPRKTable ARKodeSPRKTable_Create(int s, int q, const sunrealtype* a_1d,
                                       const sunrealtype* ahat_1d);

SUNDIALS_EXPORT
ARKodeSPRKTable ARKodeSPRKTable_Alloc(int stages);

SUNDIALS_EXPORT
ARKodeSPRKTable ARKodeSPRKTable_Load(ARKODE_SPRKMethodID id);

SUNDIALS_EXPORT
ARKodeSPRKTable ARKodeSPRKTable_LoadByName(const char* method);

SUNDIALS_EXPORT
ARKodeSPRKTable ARKodeSPRKTable_Copy(ARKodeSPRKTable that_sprk_storage);

SUNDIALS_EXPORT
void ARKodeSPRKTable_Write(ARKodeSPRKTable sprk_table, FILE* outfile);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
void ARKodeSPRKTable_Space(ARKodeSPRKTable sprk_storage, sunindextype* liw,
                           sunindextype* lrw);
SUNDIALS_EXPORT
void ARKodeSPRKTable_Free(ARKodeSPRKTable sprk_storage);

SUNDIALS_EXPORT
int ARKodeSPRKTable_ToButcher(ARKodeSPRKTable sprk_storage,
                              ARKodeButcherTable* a_ptr,
                              ARKodeButcherTable* b_ptr);

/* Different methods */

#ifdef __cplusplus
}
#endif

#endif
