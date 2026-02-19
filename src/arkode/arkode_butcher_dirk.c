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
 * This is the implementation file for built-in DIRK Butcher
 * tables.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <arkode/arkode_butcher_dirk.h>
#include <sundials/sundials_math.h>

#include "arkode_impl.h"

/*---------------------------------------------------------------
  Returns Butcher table structure for pre-set DIRK methods.

  Input:  imethod -- integer key for the desired method
  ---------------------------------------------------------------*/
ARKodeButcherTable ARKodeButcherTable_LoadDIRK(ARKODE_DIRKTableID imethod)
{
  /* Use X-macro to test each method name */
  switch (imethod)
  {
#define ARK_BUTCHER_TABLE(name, coeff) \
  case name: coeff break;
#include "arkode_butcher_dirk.def"
#undef ARK_BUTCHER_TABLE

  default:
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Unknown Butcher table");
    return NULL;
  }
}

/*---------------------------------------------------------------
  Returns Butcher table structure for pre-set DIRK methods.

  Input:  method -- string key for the desired method
  ---------------------------------------------------------------*/
ARKodeButcherTable ARKodeButcherTable_LoadDIRKByName(const char* imethod)
{
  return ARKodeButcherTable_LoadDIRK(arkButcherTableDIRKNameToID(imethod));
}

/*---------------------------------------------------------------
  Returns the string name for a pre-set DIRK method by its ID.

  Input:  imethod -- integer key for the desired method
  ---------------------------------------------------------------*/
const char* ARKodeButcherTable_DIRKIDToName(ARKODE_DIRKTableID imethod)
{
  /* Use X-macro to test each method name */
  switch (imethod)
  {
#define ARK_BUTCHER_TABLE(name, coeff) \
  case name: return #name;
#include "arkode_butcher_dirk.def"
#undef ARK_BUTCHER_TABLE

  default:
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Unknown Butcher table");
    return NULL;
  }
}

/*---------------------------------------------------------------
  Returns Butcher table ID for pre-set DIRK methods.

  Input:  method -- string key for the desired method
  ---------------------------------------------------------------*/
ARKODE_DIRKTableID arkButcherTableDIRKNameToID(const char* imethod)
{
  /* Use X-macro to test each method name */
#define ARK_BUTCHER_TABLE(name, coeff) \
  if (strcmp(#name, imethod) == 0) { return name; }
#include "arkode_butcher_dirk.def"
#undef ARK_BUTCHER_TABLE

  arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                  "Unknown Butcher table");

  return ARKODE_DIRK_NONE;
}

/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
