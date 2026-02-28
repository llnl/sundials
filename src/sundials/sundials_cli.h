/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds
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
 * SUNDIALS command-line utility definitions.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_CLI_H
#define _SUNDIALS_CLI_H

#include <stdio.h>

#include <sundials/sundials_errors.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* utilities for integer "set" routines */
typedef int (*sunIntSetFn)(void*, int);

struct sunKeyIntPair
{
  const char* key;
  sunIntSetFn set;
};

SUNDIALS_EXPORT
int sunCheckAndSetIntArgs(void* mem, int* argidx, char* argv[],
                          const size_t offset,
                          const struct sunKeyIntPair* testpairs, int numpairs,
                          sunbooleantype* arg_used, int* failedarg);

/* utilities for pair-of-integer "set" routines */
typedef int (*sunTwoIntSetFn)(void*, int, int);

struct sunKeyTwoIntPair
{
  const char* key;
  sunTwoIntSetFn set;
};

SUNDIALS_EXPORT
int sunCheckAndSetTwoIntArgs(void* mem, int* argidx, char* argv[],
                             const size_t offset,
                             const struct sunKeyTwoIntPair* testpairs,
                             int numpairs, sunbooleantype* arg_used,
                             int* failedarg);

/* utilities for long int "set" routines */
typedef int (*sunLongSetFn)(void*, long int);

struct sunKeyLongPair
{
  const char* key;
  sunLongSetFn set;
};

SUNDIALS_EXPORT
int sunCheckAndSetLongArgs(void* mem, int* argidx, char* argv[],
                           const size_t offset,
                           const struct sunKeyLongPair* testpairs, int numpairs,
                           sunbooleantype* arg_used, int* failedarg);

/* utilities for pair int/sunrealtype "set" routines */
typedef int (*sunIntRealSetFn)(void*, int, sunrealtype);

struct sunKeyIntRealPair
{
  const char* key;
  sunIntRealSetFn set;
};

SUNDIALS_EXPORT
int sunCheckAndSetIntRealArgs(void* mem, int* argidx, char* argv[],
                              const size_t offset,
                              const struct sunKeyIntRealPair* testpairs,
                              int numpairs, sunbooleantype* arg_used,
                              int* failedarg);

/* utilities for triplet int/sunrealtype/sunrealtype "set" routines */
typedef int (*sunIntRealRealSetFn)(void*, int, sunrealtype, sunrealtype);

struct sunKeyIntRealRealPair
{
  const char* key;
  sunIntRealRealSetFn set;
};

SUNDIALS_EXPORT
int sunCheckAndSetIntRealRealArgs(void* mem, int* argidx, char* argv[],
                                  const size_t offset,
                                  const struct sunKeyIntRealRealPair* testpairs,
                                  int numpairs, sunbooleantype* arg_used,
                                  int* failedarg);

/* utilities for pair int/long int "set" routines */
typedef int (*sunIntLongSetFn)(void*, int, long int);

struct sunKeyIntLongPair
{
  const char* key;
  sunIntLongSetFn set;
};

SUNDIALS_EXPORT
int sunCheckAndSetIntLongArgs(void* mem, int* argidx, char* argv[],
                              const size_t offset,
                              const struct sunKeyIntLongPair* testpairs,
                              int numpairs, sunbooleantype* arg_used,
                              int* failedarg);

/* utilities for sunrealtype "set" routines */
typedef int (*sunRealSetFn)(void*, sunrealtype);

struct sunKeyRealPair
{
  const char* key;
  sunRealSetFn set;
};

SUNDIALS_EXPORT
int sunCheckAndSetRealArgs(void* mem, int* argidx, char* argv[],
                           const size_t offset,
                           const struct sunKeyRealPair* testpairs, int numpairs,
                           sunbooleantype* arg_used, int* failedarg);

/* utilities for pair-of-sunrealtype "set" routines */
typedef int (*sunTwoRealSetFn)(void*, sunrealtype, sunrealtype);

struct sunKeyTwoRealPair
{
  const char* key;
  sunTwoRealSetFn set;
};

SUNDIALS_EXPORT
int sunCheckAndSetTwoRealArgs(void* mem, int* argidx, char* argv[],
                              const size_t offset,
                              const struct sunKeyTwoRealPair* testpairs,
                              int numpairs, sunbooleantype* arg_used,
                              int* failedarg);

/* utilities for char* "set" routines */
typedef int (*sunCharSetFn)(void*, const char*);

struct sunKeyCharPair
{
  const char* key;
  sunCharSetFn set;
};

SUNDIALS_EXPORT
int sunCheckAndSetCharArgs(void* mem, int* argidx, char* argv[],
                           const size_t offset,
                           const struct sunKeyCharPair* testpairs, int numpairs,
                           sunbooleantype* arg_used, int* failedarg);

/* utilities for pair-of-char* "set" routines */
typedef int (*sunTwoCharSetFn)(void*, const char*, const char*);

struct sunKeyTwoCharPair
{
  const char* key;
  sunTwoCharSetFn set;
};

SUNDIALS_EXPORT
int sunCheckAndSetTwoCharArgs(void* mem, int* argidx, char* argv[],
                              const size_t offset,
                              const struct sunKeyTwoCharPair* testpairs,
                              int numpairs, sunbooleantype* arg_used,
                              int* failedarg);

/* utilities for action "set" routines */
typedef int (*sunActionSetFn)(void*);

struct sunKeyActionPair
{
  const char* key;
  sunActionSetFn set;
};

SUNDIALS_EXPORT
int sunCheckAndSetActionArgs(void* mem, int* argidx, char* argv[],
                             const size_t offset,
                             const struct sunKeyActionPair* testpairs,
                             int numpairs, sunbooleantype* arg_used,
                             int* failedarg);

#ifdef __cplusplus
}
#endif

#endif
