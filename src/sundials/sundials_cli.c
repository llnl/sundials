/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
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
 * This file provides command-line control over optional inputs
 * to ARKODE.
 *--------------------------------------------------------------*/

#include "sundials_cli.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

/*===============================================================
  Command-line input utility routines
  ===============================================================*/

int sunCheckAndSetIntArgs(void* mem, int* argidx, char* argv[],
                          const size_t offset,
                          const struct sunKeyIntPair* testpairs, int numpairs,
                          sunbooleantype* arg_used, int* failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    *arg_used = SUNFALSE;
    if (strcmp(argv[*argidx] + offset, testpairs[j].key) == 0)
    {
      (*argidx) += 1;
      int iarg   = atoi(argv[*argidx]);
      int retval = testpairs[j].set(mem, iarg);
      if (retval != SUN_SUCCESS)
      {
        *failedarg = j;
        return retval;
      }
      *arg_used = SUNTRUE;
      return SUN_SUCCESS;
    }
  }
  return SUN_SUCCESS;
}

int sunCheckAndSetTwoIntArgs(void* mem, int* argidx, char* argv[],
                             const size_t offset,
                             const struct sunKeyTwoIntPair* testpairs,
                             int numpairs, sunbooleantype* arg_used,
                             int* failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    *arg_used = SUNFALSE;
    if (strcmp(argv[*argidx] + offset, testpairs[j].key) == 0)
    {
      (*argidx) += 1;
      int iarg1 = atoi(argv[*argidx]);
      (*argidx) += 1;
      int iarg2  = atoi(argv[*argidx]);
      int retval = testpairs[j].set(mem, iarg1, iarg2);
      if (retval != SUN_SUCCESS)
      {
        *failedarg = j;
        return retval;
      }
      *arg_used = SUNTRUE;
      return SUN_SUCCESS;
    }
  }
  return SUN_SUCCESS;
}

int sunCheckAndSetLongArgs(void* mem, int* argidx, char* argv[],
                           const size_t offset,
                           const struct sunKeyLongPair* testpairs, int numpairs,
                           sunbooleantype* arg_used, int* failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    *arg_used = SUNFALSE;
    if (strcmp(argv[*argidx] + offset, testpairs[j].key) == 0)
    {
      (*argidx) += 1;
      long int iarg = atol(argv[*argidx]);
      int retval    = testpairs[j].set(mem, iarg);
      if (retval != SUN_SUCCESS)
      {
        *failedarg = j;
        return retval;
      }
      *arg_used = SUNTRUE;
      return SUN_SUCCESS;
    }
  }
  return SUN_SUCCESS;
}

int sunCheckAndSetIntRealArgs(void* mem, int* argidx, char* argv[],
                              const size_t offset,
                              const struct sunKeyIntRealPair* testpairs,
                              int numpairs, sunbooleantype* arg_used,
                              int* failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    *arg_used = SUNFALSE;
    if (strcmp(argv[*argidx] + offset, testpairs[j].key) == 0)
    {
      (*argidx) += 1;
      int iarg = atoi(argv[*argidx]);
      (*argidx) += 1;
      sunrealtype rarg = SUNStrToReal(argv[*argidx]);
      int retval       = testpairs[j].set(mem, iarg, rarg);
      if (retval != SUN_SUCCESS)
      {
        *failedarg = j;
        return retval;
      }
      *arg_used = SUNTRUE;
      return SUN_SUCCESS;
    }
  }
  return SUN_SUCCESS;
}

int sunCheckAndSetIntRealRealArgs(void* mem, int* argidx, char* argv[],
                                  const size_t offset,
                                  const struct sunKeyIntRealRealPair* testpairs,
                                  int numpairs, sunbooleantype* arg_used,
                                  int* failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    *arg_used = SUNFALSE;
    if (strcmp(argv[*argidx] + offset, testpairs[j].key) == 0)
    {
      (*argidx) += 1;
      int iarg = atoi(argv[*argidx]);
      (*argidx) += 1;
      sunrealtype rarg1 = SUNStrToReal(argv[*argidx]);
      (*argidx) += 1;
      sunrealtype rarg2 = SUNStrToReal(argv[*argidx]);
      int retval        = testpairs[j].set(mem, iarg, rarg1, rarg2);
      if (retval != SUN_SUCCESS)
      {
        *failedarg = j;
        return retval;
      }
      *arg_used = SUNTRUE;
      return SUN_SUCCESS;
    }
  }
  return SUN_SUCCESS;
}

int sunCheckAndSetIntLongArgs(void* mem, int* argidx, char* argv[],
                              const size_t offset,
                              const struct sunKeyIntLongPair* testpairs,
                              int numpairs, sunbooleantype* arg_used,
                              int* failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    *arg_used = SUNFALSE;
    if (strcmp(argv[*argidx] + offset, testpairs[j].key) == 0)
    {
      (*argidx) += 1;
      int iarg = atoi(argv[*argidx]);
      (*argidx) += 1;
      long int large = atol(argv[*argidx]);
      int retval     = testpairs[j].set(mem, iarg, large);
      if (retval != SUN_SUCCESS)
      {
        *failedarg = j;
        return retval;
      }
      *arg_used = SUNTRUE;
      return SUN_SUCCESS;
    }
  }
  return SUN_SUCCESS;
}

int sunCheckAndSetRealArgs(void* mem, int* argidx, char* argv[],
                           const size_t offset,
                           const struct sunKeyRealPair* testpairs, int numpairs,
                           sunbooleantype* arg_used, int* failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    *arg_used = SUNFALSE;
    if (strcmp(argv[*argidx] + offset, testpairs[j].key) == 0)
    {
      (*argidx) += 1;
      sunrealtype rarg = SUNStrToReal(argv[*argidx]);
      int retval       = testpairs[j].set(mem, rarg);
      if (retval != SUN_SUCCESS)
      {
        *failedarg = j;
        return retval;
      }
      *arg_used = SUNTRUE;
      return SUN_SUCCESS;
    }
  }
  return SUN_SUCCESS;
}

int sunCheckAndSetTwoRealArgs(void* mem, int* argidx, char* argv[],
                              const size_t offset,
                              const struct sunKeyTwoRealPair* testpairs,
                              int numpairs, sunbooleantype* arg_used,
                              int* failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    *arg_used = SUNFALSE;
    if (strcmp(argv[*argidx] + offset, testpairs[j].key) == 0)
    {
      (*argidx) += 1;
      sunrealtype rarg1 = SUNStrToReal(argv[*argidx]);
      (*argidx) += 1;
      sunrealtype rarg2 = SUNStrToReal(argv[*argidx]);
      int retval        = testpairs[j].set(mem, rarg1, rarg2);
      if (retval != SUN_SUCCESS)
      {
        *failedarg = j;
        return retval;
      }
      *arg_used = SUNTRUE;
      return SUN_SUCCESS;
    }
  }
  return SUN_SUCCESS;
}

int sunCheckAndSetCharArgs(void* mem, int* argidx, char* argv[],
                           const size_t offset,
                           const struct sunKeyCharPair* testpairs, int numpairs,
                           sunbooleantype* arg_used, int* failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    *arg_used = SUNFALSE;
    if (strcmp(argv[*argidx] + offset, testpairs[j].key) == 0)
    {
      (*argidx) += 1;
      int retval = testpairs[j].set(mem, argv[*argidx]);
      if (retval != SUN_SUCCESS)
      {
        *failedarg = j;
        return retval;
      }
      *arg_used = SUNTRUE;
      return SUN_SUCCESS;
    }
  }
  return SUN_SUCCESS;
}

int sunCheckAndSetTwoCharArgs(void* mem, int* argidx, char* argv[],
                              const size_t offset,
                              const struct sunKeyTwoCharPair* testpairs,
                              int numpairs, sunbooleantype* arg_used,
                              int* failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    *arg_used = SUNFALSE;
    if (strcmp(argv[*argidx] + offset, testpairs[j].key) == 0)
    {
      int retval = testpairs[j].set(mem, argv[*argidx + 1], argv[*argidx + 2]);
      (*argidx) += 2;
      if (retval != SUN_SUCCESS)
      {
        *failedarg = j;
        return retval;
      }
      *arg_used = SUNTRUE;
      return SUN_SUCCESS;
    }
  }
  return SUN_SUCCESS;
}

int sunCheckAndSetActionArgs(void* mem, int* argidx, char* argv[],
                             const size_t offset,
                             const struct sunKeyActionPair* testpairs,
                             int numpairs, sunbooleantype* arg_used,
                             int* failedarg)
{
  for (int j = 0; j < numpairs; j++)
  {
    *arg_used = SUNFALSE;
    if (strcmp(argv[*argidx] + offset, testpairs[j].key) == 0)
    {
      int retval = testpairs[j].set(mem);
      if (retval != SUN_SUCCESS)
      {
        *failedarg = j;
        return retval;
      }
      *arg_used = SUNTRUE;
      return SUN_SUCCESS;
    }
  }
  return SUN_SUCCESS;
}
