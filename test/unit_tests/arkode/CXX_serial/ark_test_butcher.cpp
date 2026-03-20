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
 *----------------------------------------------------------------
 * Routine to test all built-in Butcher tables to check their
 * order of accuracy.
 *-----------------------------------------------------------------*/

// Header files
#include <arkode/arkode.h>
#include <arkode/arkode_butcher_dirk.h>
#include <arkode/arkode_butcher_erk.h>
#include <iostream>
#include <ostream>
#include <sundials/sundials_types.h>
#include <vector>

#include "arkode/arkode_impl.h"

struct ARK_Table
{
  const char* const name;
  const char* const erk_table;
  const char* const dirk_table;
};

// Main Program
int main()
{
  std::vector<ARK_Table> ark_tables =
    {{"ARKODE_ARK2_3_1_2", "ARKODE_ARK2_ERK_3_1_2", "ARKODE_ARK2_DIRK_3_1_2"},
     {"ARKODE_ARK324L2SA_4_2_3", "ARKODE_ARK324L2SA_ERK_4_2_3",
      "ARKODE_ARK324L2SA_DIRK_4_2_3"},
     {"ARKODE_ARK436L2SA_6_3_4", "ARKODE_ARK436L2SA_ERK_6_3_4",
      "ARKODE_ARK436L2SA_DIRK_6_3_4"},
     {"ARKODE_ARK437L2SA_7_3_4", "ARKODE_ARK437L2SA_ERK_7_3_4",
      "ARKODE_ARK437L2SA_DIRK_7_3_4"},
     {"ARKODE_ARK548L2SA_8_4_5", "ARKODE_ARK548L2SA_ERK_8_4_5",
      "ARKODE_ARK548L2SA_DIRK_8_4_5"},
     {"ARKODE_ARK548L2SAb_8_4_5", "ARKODE_ARK548L2SAb_ERK_8_4_5",
      "ARKODE_ARK548L2SAb_DIRK_8_4_5"}};

  int numfails = 0;

  // loop over individual ERK tables
  std::cout << "\nTesting individual ERK methods:\n\n";
  for (int i = ARKODE_MIN_ERK_NUM; i <= ARKODE_MAX_ERK_NUM; i++)
  {
    ARKODE_ERKTableID id = static_cast<ARKODE_ERKTableID>(i);
    std::cout << "Testing method " << ARKodeButcherTable_ERKIDToName(id) << ":";

    // load Butcher table
    ARKodeButcherTable B = ARKodeButcherTable_LoadERK(id);
    if (B == NULL)
    {
      std::cout << "  error retrieving table, aborting\n";
      return 1;
    }

    // test method, upon questionable results, re-test with statistics printing
    int q, p, retval;
    retval = ARKodeButcherTable_CheckOrder(B, &q, &p, NULL);
    if (retval == 0)
    {
      std::cout << "  table matches predicted method/embedding orders of " << q
                << "/" << p << "\n";
    }
    else if (retval > 0)
    {
      std::cout << "  WARNING:\n";
      retval = ARKodeButcherTable_CheckOrder(B, &q, &p, stdout);
    }
    else
    {
      std::cout << "  ERROR:\n";
      retval = ARKodeButcherTable_CheckOrder(B, &q, &p, stdout);
      numfails++;
    }

    // clean up after this test
    ARKodeButcherTable_Free(B);
  }

  // loop over individual DIRK tables
  std::cout << "\nTesting individual DIRK methods:\n\n";
  for (int i = ARKODE_MIN_DIRK_NUM; i <= ARKODE_MAX_DIRK_NUM; i++)
  {
    ARKODE_DIRKTableID id = static_cast<ARKODE_DIRKTableID>(i);
    std::cout << "Testing method " << ARKodeButcherTable_DIRKIDToName(id) << ":";

    // load Butcher table
    ARKodeButcherTable B = ARKodeButcherTable_LoadDIRK(id);
    if (B == NULL)
    {
      std::cout << "  error retrieving table, aborting\n";
      return 1;
    }

    // test method, upon questionable results, re-test with statistics printing
    int q, p, retval;
    retval = ARKodeButcherTable_CheckOrder(B, &q, &p, NULL);
    if (retval == 0)
    {
      std::cout << "  table matches predicted method/embedding orders of " << q
                << "/" << p << "\n";
    }
    else if (retval > 0)
    {
      std::cout << "  WARNING:\n";
      retval = ARKodeButcherTable_CheckOrder(B, &q, &p, stdout);
    }
    else
    {
      std::cout << "  ERROR:\n";
      retval = ARKodeButcherTable_CheckOrder(B, &q, &p, stdout);
      numfails++;
    }

    // clean up after this test
    ARKodeButcherTable_Free(B);
  }

  // loop over ARK pairs
  std::cout << "\nTesting ARK pairs:\n\n";
  for (ARK_Table& ark_table : ark_tables)
  {
    std::cout << "Testing method " << ark_table.name << ":";

    // load Butcher tables
    ARKodeButcherTable Be = ARKodeButcherTable_LoadERKByName(ark_table.erk_table);
    if (Be == NULL)
    {
      std::cout << "  error retrieving explicit table, aborting\n";
      return 1;
    }
    ARKodeButcherTable Bi =
      ARKodeButcherTable_LoadDIRKByName(ark_table.dirk_table);
    if (Bi == NULL)
    {
      std::cout << "  error retrieving implicit table, aborting";
      return 1;
    }

    // test method, upon questionable results, re-test with statistics printing
    int q, p, retval;
    retval = ARKodeButcherTable_CheckARKOrder(Be, Bi, &q, &p, NULL);
    if (retval == 0)
    {
      std::cout << "  Method/embedding match predicted orders of " << q << "/"
                << p << "\n";
    }
    else if (retval > 0)
    {
      std::cout << "  WARNING:\n";
      retval = ARKodeButcherTable_CheckARKOrder(Be, Bi, &q, &p, stdout);
    }
    else
    {
      std::cout << "  ERROR:\n";
      retval = ARKodeButcherTable_CheckARKOrder(Be, Bi, &q, &p, stdout);
      numfails++;
    }

    // clean up after this test
    ARKodeButcherTable_Free(Be);
    ARKodeButcherTable_Free(Bi);
  }

  // determine overall success/failure and return
  if (numfails == 0)
  {
    std::cout << "All Butcher tables passed\n";
    return 0;
  }
  else
  {
    std::cout << numfails << " Butcher tables failed\n";
    return 1;
  }
}

/*---- end of file ----*/
