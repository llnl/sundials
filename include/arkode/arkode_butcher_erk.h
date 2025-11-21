/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
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
 * This is the header file for ARKode's built-in ERK Butcher tables.
 * -----------------------------------------------------------------*/

#ifndef _ARKODE_ERK_TABLES_H
#define _ARKODE_ERK_TABLES_H

#include <arkode/arkode_butcher.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

enum ARKODE_ERKTableID
{
  ARKODE_ERK_NONE = -1, /* ensure enum is signed int */
  /* WARNING:  ARKODE_MIN_ERK_NUM must come after the first entry, ARKODE_HEUN_EULER_2_1_2,
     because Python enums will only expose the member that is defined first. Due to
     this and how pybind/nanobind handle the enums, if we defined ARKODE_MRI_NUM first,
     then ARKODE_HEUN_EULER_2_1_2 would not be usable from the module scope (the MIN/MAX) entries
     will still be usable when accessing through the IntEnum object, but not from module scope. */
  ARKODE_HEUN_EULER_2_1_2 = 0,
  ARKODE_MIN_ERK_NUM      = 0,
  ARKODE_BOGACKI_SHAMPINE_4_2_3,
  ARKODE_ARK324L2SA_ERK_4_2_3,
  ARKODE_ZONNEVELD_5_3_4,
  ARKODE_ARK436L2SA_ERK_6_3_4,
  ARKODE_SAYFY_ABURUB_6_3_4,
  ARKODE_CASH_KARP_6_4_5,
  ARKODE_FEHLBERG_6_4_5,
  ARKODE_DORMAND_PRINCE_7_4_5,
  ARKODE_ARK548L2SA_ERK_8_4_5,
  ARKODE_VERNER_8_5_6,
  ARKODE_FEHLBERG_13_7_8,
  ARKODE_KNOTH_WOLKE_3_3,
  ARKODE_ARK437L2SA_ERK_7_3_4,
  ARKODE_ARK548L2SAb_ERK_8_4_5,
  ARKODE_ARK2_ERK_3_1_2,
  ARKODE_SOFRONIOU_SPALETTA_5_3_4,
  ARKODE_SHU_OSHER_3_2_3,
  ARKODE_VERNER_9_5_6,
  ARKODE_VERNER_10_6_7,
  ARKODE_VERNER_13_7_8,
  ARKODE_VERNER_16_8_9,
  ARKODE_FORWARD_EULER_1_1,
  ARKODE_RALSTON_EULER_2_1_2,
  ARKODE_EXPLICIT_MIDPOINT_EULER_2_1_2,
  ARKODE_RALSTON_3_1_2,
  ARKODE_TSITOURAS_7_4_5,
  ARKODE_MAX_ERK_NUM = ARKODE_TSITOURAS_7_4_5
};

#ifndef SWIG
typedef enum ARKODE_ERKTableID ARKODE_ERKTableID;
#endif

/* Accessor routine to load built-in ERK table */
SUNDIALS_EXPORT ARKodeButcherTable
ARKodeButcherTable_LoadERK(ARKODE_ERKTableID emethod);

SUNDIALS_EXPORT ARKodeButcherTable
ARKodeButcherTable_LoadERKByName(const char* emethod);

SUNDIALS_EXPORT const char* ARKodeButcherTable_ERKIDToName(ARKODE_ERKTableID emethod);

#ifdef __cplusplus
}
#endif

#endif
