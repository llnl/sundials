/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
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
 * This is the implementation file for built-in ERK Butcher
 * tables.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <arkode/arkode_butcher_erk.h>
#include <sundials/sundials_math.h>

#include "arkode_impl.h"

/*---------------------------------------------------------------
  Returns Butcher table structure for pre-set Runge Kutta methods.

  Input:  emethod -- integer key for the desired method
  ---------------------------------------------------------------*/
ARKodeButcherTable ARKodeButcherTable_LoadERK(ARKODE_ERKTableID emethod)
{
  /* Use X-macro to test each method name */
  switch (emethod)
  {
#define ARK_BUTCHER_TABLE(name, coeff) \
  case name: coeff break;
#include "arkode_butcher_erk.def"
#undef ARK_BUTCHER_TABLE

  default:
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Unknown Butcher table");
    return NULL;
  }
}

/*---------------------------------------------------------------
  Returns Butcher table structure for pre-set Runge Kutta methods.

  Input:  emethod -- string key for the desired method
  ---------------------------------------------------------------*/
ARKodeButcherTable ARKodeButcherTable_LoadERKByName(const char* emethod)
{
  return ARKodeButcherTable_LoadERK(arkButcherTableERKNameToID(emethod));
}

/*---------------------------------------------------------------
  Returns the string name for a pre-set Runge Kutta method by its ID.

  Input:  emethod -- integer key for the desired method
  ---------------------------------------------------------------*/
const char* ARKodeButcherTable_ERKIDToName(ARKODE_ERKTableID emethod)
{
  /* Use X-macro to test each method name */
  switch (emethod)
  {
#define ARK_BUTCHER_TABLE(name, coeff) \
  case name: return #name;
#include "arkode_butcher_erk.def"
#undef ARK_BUTCHER_TABLE

  default:
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Unknown Butcher table");
    return NULL;
  }
}

/*---------------------------------------------------------------
  Returns Butcher table ID for pre-set Runge Kutta methods.

  Input:  emethod -- string key for the desired method
  ---------------------------------------------------------------*/
ARKODE_ERKTableID arkButcherTableERKNameToID(const char* emethod)
{
  /* Use X-macro to test each method name */
#define ARK_BUTCHER_TABLE(name, coeff) \
  if (strcmp(#name, emethod) == 0) { return name; }
#include "arkode_butcher_erk.def"
#undef ARK_BUTCHER_TABLE

  arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                  "Unknown Butcher table");

  return ARKODE_ERK_NONE;
}

/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
