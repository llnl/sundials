/* -----------------------------------------------------------------
 * Programmer: Cody J. Balos @ LLNL
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
 * A simple hashmap implementation for char* keys and
 * void* values. Uses linear probing to resolve collisions.
 * The values can be anything, but will be freed by
 * the hash map upon its destruction.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_HASHMAP_IMPL_H
#define _SUNDIALS_HASHMAP_IMPL_H

#include <stdlib.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SUNHASHMAP_ERROR       -99
#define SUNHASHMAP_KEYNOTFOUND -1
#define SUNHASHMAP_DUPLICATE   -2

struct SUNHashMapKeyValue_
{
  char* key;
  void* value;
};

typedef struct SUNHashMapKeyValue_* SUNHashMapKeyValue;

#define TTYPE SUNHashMapKeyValue
#include "stl/sunstl_vector.h"
#undef TTYPE

struct SUNHashMap_
{
  SUNErrCode (*destroyKeyValue)(SUNHashMapKeyValue*);
  SUNStlVector_SUNHashMapKeyValue buckets;
};

typedef struct SUNHashMap_* SUNHashMap;

SUNErrCode SUNHashMap_New(int64_t capacity,
                          SUNErrCode (*destroyValue)(SUNHashMapKeyValue* value_ptr),
                          SUNHashMap* map);

int64_t SUNHashMap_Capacity(SUNHashMap map);

SUNErrCode SUNHashMap_Destroy(SUNHashMap* map);

int64_t SUNHashMap_Iterate(SUNHashMap map, int64_t start,
                           int64_t (*yieldfn)(int64_t, SUNHashMapKeyValue,
                                              const void*),
                           const void* ctx);

int64_t SUNHashMap_Insert(SUNHashMap map, const char* key, void* value);

int64_t SUNHashMap_GetValue(SUNHashMap map, const char* key, void** value);

int64_t SUNHashMap_Remove(SUNHashMap map, const char* key, void** value);

SUNErrCode SUNHashMap_Sort(SUNHashMap map, SUNHashMapKeyValue** sorted,
                           int (*compar)(const void*, const void*));

SUNErrCode SUNHashMap_Values(SUNHashMap map, void*** values, int64_t value_size);

SUNErrCode SUNHashMap_PrintKeys(SUNHashMap map, FILE* file);

#ifdef __cplusplus
}
#endif

#endif
