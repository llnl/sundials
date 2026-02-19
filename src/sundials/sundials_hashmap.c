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

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <sundials/sundials_errors.h>
#include <sundials/sundials_types.h>

#include "sundials_hashmap_impl.h"
#include "sundials_macros.h"

static const uint64_t HASH_PRIME        = 14695981039346656037U;
static const uint64_t HASH_OFFSET_BASIS = 1099511628211U;

/*
  For a nice discussion on popular hashing algorithms see:
  https://softwareengineering.stackexchange.com/questions/49550/which-hashing-algorithm-is-best-for-uniqueness-and-speed/145633#145633

  This is a 64-bit implementation of the 'a' modification of the
  Fowler-Noll-Vo hash (i.e., FNV1-a).
 */
static uint64_t fnv1a_hash(const char* str)
{
  uint64_t hash = HASH_OFFSET_BASIS;
  char c;
  while ((c = *str++)) { hash = (hash ^ c) * HASH_PRIME; }
  return hash;
}

static inline int64_t sunHashMapIdxFromKey(SUNHashMap map, const char* key)
{
  /* We want the index to be in [0, SUNHashMap_Capacity(map)) */
  int64_t end = SUNHashMap_Capacity(map) - 1;
  int64_t idx = end == 0 ? end : (int64_t)(fnv1a_hash(key) % end);
  return idx;
}

/*
  This function creates a new SUNHashMap object allocated to hold
  up to 'capacity' entries.

  **Arguments:**
    * ``capacity`` -- the initial capactity number of the hashmap
    * ``destroyKeyValue`` -- a callback function that frees both the key and value.
    * ``map`` -- on input, a SUNHasMap pointer, on output the SUNHashMap will be
                 allocated

  **Returns:**
    * A SUNErrCode indicating success or a failure
 */
SUNErrCode SUNHashMap_New(int64_t capacity,
                          SUNErrCode (*destroyKeyValue)(SUNHashMapKeyValue* kv_ptr),
                          SUNHashMap* map)
{
  if (capacity <= 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  *map = NULL;
  *map = (SUNHashMap)malloc(sizeof(**map));

  if (!map) { return SUN_ERR_MALLOC_FAIL; }

  SUNStlVector_SUNHashMapKeyValue buckets =
    SUNStlVector_SUNHashMapKeyValue_New(capacity, destroyKeyValue);
  if (!buckets)
  {
    free(*map);
    return SUN_ERR_MALLOC_FAIL;
  }

  /* Initialize all buckets to NULL */
  for (int64_t i = 0; i < capacity; i++)
  {
    SUNErrCode err = SUNStlVector_SUNHashMapKeyValue_PushBack(buckets, NULL);
    if (err) { return err; };
  }

  (*map)->buckets = buckets;

  return SUN_SUCCESS;
}

/*
  This function returns the capacity of the hashmap.

  **Arguments:**
    * ``map`` -- the SUNHashMap object

  **Returns:**
    * The capacity of the hashmap
 */
int64_t SUNHashMap_Capacity(SUNHashMap map)
{
  return SUNStlVector_SUNHashMapKeyValue_Capacity(map->buckets);
}

/*
  This function frees the SUNHashMap object.

  **Arguments:**
    * ``map`` -- on input, a SUNHasMap pointer, on output the SUNHashMap will be
                 deallocated and set to ``NULL``

  **Returns:**
    * A SUNErrCode indicating success or a failure
 */
SUNErrCode SUNHashMap_Destroy(SUNHashMap* map)
{
  if (map == NULL) { return SUN_SUCCESS; }

  SUNErrCode err = SUNStlVector_SUNHashMapKeyValue_Destroy(&(*map)->buckets);
  if (err) { return err; }
  free(*map);
  *map = NULL;

  return SUN_SUCCESS;
}

/*
  This function iterates the map over the range [start, N]. N is either the
  index at which ``yieldfn`` indicates the iteration should stop, or the max
  entries in the map.

  **Arguments:**
    * ``map`` -- the ``SUNHashMap`` object to operate on
    * ``start`` -- the start of the iteration range
    * ``yieldfn`` -- the callback function to call every iteration
                     this should return SUNHASHMAP_ERROR to continue the iteration, or [0, SUNHASHMAP_KEYNOTFOUND]
                     to stop; the first argument is the current index, the second argument
                     is the current key-value pair, and the final argument is the same
                    pointer ``ctx`` as the final argument  to SUNHashMapIterate.
    * ``ctx`` -- a pointer to pass on to ``yieldfn``

  **Returns:**
    * ``SUNHASHMAP_ERROR`` -- an error occurred
    * ``capacity`` -- iterated the whole map
    * ``>=0`` -- the index at which the iteration stopped
 */
int64_t SUNHashMap_Iterate(SUNHashMap map, int64_t start,
                           int64_t (*yieldfn)(int64_t, SUNHashMapKeyValue,
                                              const void*),
                           const void* ctx)
{
  if (map == NULL || yieldfn == NULL) { return SUNHASHMAP_ERROR; }

  for (int64_t i = start;
       i < SUNStlVector_SUNHashMapKeyValue_Size(map->buckets); i++)
  {
    int64_t retval =
      yieldfn(i, *SUNStlVector_SUNHashMapKeyValue_At(map->buckets, i), ctx);
    if (retval == SUNHASHMAP_ERROR) { continue; /* keep looking */ }
    else { return (retval); /* yieldfn indicates the loop should break */ }
  }

  return SUNHashMap_Capacity(map);
}

static int64_t sunHashMapLinearProbeInsert(int64_t idx, SUNHashMapKeyValue kv,
                                           SUNDIALS_MAYBE_UNUSED const void* ctx)
{
  /* find the next open spot */
  if (kv == NULL) { return (idx); /* open spot found at idx */ }
  return SUNHASHMAP_ERROR; /* keep looking */
}

static SUNErrCode sunHashMapResize(SUNHashMap map)
{
  int64_t old_capacity = SUNHashMap_Capacity(map);
  int64_t new_capacity = old_capacity == 0
                           ? 2
                           : (int64_t)(ceill(((long double)old_capacity) *
                                             SUNSTLVECTOR_GROWTH_FACTOR));

  SUNStlVector_SUNHashMapKeyValue old_buckets = map->buckets;
  map->buckets = SUNStlVector_SUNHashMapKeyValue_New(new_capacity,
                                                     map->buckets->destroyValue);

  /* Set all buckets to NULL */
  for (int64_t i = 0; i < new_capacity; i++)
  {
    SUNErrCode err = SUNStlVector_SUNHashMapKeyValue_PushBack(map->buckets, NULL);
    if (err) { return err; }
  }

  /* Rehash and reinsert */
  for (int64_t i = old_capacity - 1; i >= 0; i--)
  {
    SUNHashMapKeyValue kvp = *SUNStlVector_SUNHashMapKeyValue_At(old_buckets, i);
    if (kvp)
    {
      SUNHashMap_Insert(map, kvp->key, kvp->value);
      free(kvp->key);
      free(kvp);
    }
    SUNErrCode err = SUNStlVector_SUNHashMapKeyValue_PopBack(old_buckets);
    if (err) { return err; }
  }

  return SUNStlVector_SUNHashMapKeyValue_Destroy(&old_buckets);
}

/*
  This function creates a key-value pair and attempts to insert it into the map.
  Will use linear probing if there is a collision.

  **Arguments:**
    * ``map`` -- the ``SUNHashMap`` object to operate on
    * ``key`` -- the key to store (we will make a copy)
    * ``value`` -- the value associated with the key

  **Returns:**
    * ``0`` -- success
    * ``SUNHASHMAP_ERROR`` -- an error occurred
    * ``SUNHASHMAP_DUPLICATE`` -- duplicate key
 */
int64_t SUNHashMap_Insert(SUNHashMap map, const char* key, void* value)
{
  int64_t idx;
  int64_t retval;
  SUNHashMapKeyValue kvp;

  if (map == NULL || key == NULL || value == NULL) { return SUNHASHMAP_ERROR; }

  idx = sunHashMapIdxFromKey(map, key);

  /* Check if the bucket is already filled (i.e., we might have had a collision) */
  kvp = *SUNStlVector_SUNHashMapKeyValue_At(map->buckets, idx);
  if (kvp != NULL)
  {
    /* Determine if key is actually a duplicate (not allowed) */
    if (!strcmp(key, kvp->key)) { return SUNHASHMAP_DUPLICATE; }

    /* OK, it was a real collision, so find the next open spot */
    retval = SUNHashMap_Iterate(map, idx + 1, sunHashMapLinearProbeInsert, NULL);
    if (retval == SUNHASHMAP_ERROR)
    {
      /* an error occurred */
      return retval;
    }
    else if (retval == SUNHashMap_Capacity(map))
    {
      /* the map is out of empty buckets, so we grow it */
      SUNErrCode err = sunHashMapResize(map);
      if (err) { return err; }
      return SUNHashMap_Insert(map, key, value);
    }

    idx = retval;
  }

  /* Create the key-value pair */
  kvp = (SUNHashMapKeyValue)malloc(sizeof(*kvp));

  /* Copy the original_key so that the hashmap owns it */
  size_t len     = strlen(key) + 1;
  char* key_copy = malloc(sizeof(*key) * len);
  if (!key_copy)
  {
    free(kvp);
    return SUNHASHMAP_ERROR;
  }
  strcpy(key_copy, key);

  kvp->key   = key_copy;
  kvp->value = value;

  /* Insert the key-value pair */
  return SUNStlVector_SUNHashMapKeyValue_Set(map->buckets, idx, kvp);
}

static int64_t sunHashMapLinearProbeGet(int64_t idx, SUNHashMapKeyValue kv,
                                        const void* key)
{
  /* target key cannot be NULL */
  if (key == NULL) { return SUNHASHMAP_KEYNOTFOUND; }

  /* find the matching entry */
  if (kv == NULL)
  {
    return SUNHASHMAP_ERROR; /* keep looking since this bucket is empty */
  }
  if (!strcmp(kv->key, (const char*)key))
  {
    return (idx); /* found it at idx */
  }
  return SUNHASHMAP_ERROR; /* keep looking */
}

/*
  This function gets the value for the given key.

  **Arguments:**
    * ``map`` -- the ``SUNHashMap`` object to operate on
    * ``key`` -- the key to look up
    * ``value`` -- the value associated with the key

  **Returns:**
    * ``0`` -- success
    * ``SUNHASHMAP_ERROR`` -- an error occurred
    * ``SUNHASHMAP_KEYNOTFOUND`` -- key not found
 */
int64_t SUNHashMap_GetValue(SUNHashMap map, const char* key, void** value)
{
  int64_t idx;
  int64_t retval;
  sunbooleantype collision;

  if (map == NULL || key == NULL || value == NULL) { return SUNHASHMAP_ERROR; }

  idx = sunHashMapIdxFromKey(map, key);

  SUNHashMapKeyValue kvp = *SUNStlVector_SUNHashMapKeyValue_At(map->buckets, idx);

  /* Check for a collision (NULL kvp means there was a collision at one point, but
     the colliding key has since been removed)*/
  collision = kvp ? strcmp(kvp->key, key) : SUNTRUE;

  /* Resolve a collision via linear probing */
  if (collision)
  {
    retval = SUNHashMap_Iterate(map, idx + 1, sunHashMapLinearProbeGet, key);
    if (retval == SUNHASHMAP_ERROR)
    {
      /* the key was either not found anywhere or an error occurred */
      return retval;
    }
    else { idx = retval; }
  }

  /* Return a reference to the value only */
  SUNHashMapKeyValue* kvp_ptr = SUNStlVector_SUNHashMapKeyValue_At(map->buckets,
                                                                   idx);
  if (kvp_ptr) { *value = (*kvp_ptr)->value; }
  else { return SUNHASHMAP_KEYNOTFOUND; }

  return (0);
}

/*
  This function remove the key-value pair.

  **Arguments:**
    * ``map`` -- the ``SUNHashMap`` object to operate on
    * ``key`` -- the key to remove
    * ``value`` -- the value to remove

  **Returns:**
    * ``0`` -- success
    * ``SUNHASHMAP_ERROR`` -- an error occurred
    * ``SUNHASHMAP_KEYNOTFOUND`` -- key not found
 */
int64_t SUNHashMap_Remove(SUNHashMap map, const char* key, void** value)
{
  int64_t idx;
  int64_t retval;
  sunbooleantype collision;

  if (map == NULL || key == NULL) { return SUNHASHMAP_ERROR; }

  idx = sunHashMapIdxFromKey(map, key);

  SUNHashMapKeyValue kvp = *SUNStlVector_SUNHashMapKeyValue_At(map->buckets, idx);

  /* Check for a collision (NULL kvp means there was a collision at one point, but
     the colliding key has since been removed)*/
  collision = kvp ? strcmp(kvp->key, key) : SUNTRUE;

  /* Check to see if this is a collision */
  if (collision)
  {
    /* Keys did not match, so we have a collision and need to probe */
    retval = SUNHashMap_Iterate(map, idx + 1, sunHashMapLinearProbeGet,
                                (const void*)key);
    if (retval == SUNHASHMAP_ERROR)
    {
      /* an error occurred or the key was not found anywhere */
      return retval;
    }
    else { idx = retval; }
  }

  /* Return a reference to the value only */
  kvp    = *SUNStlVector_SUNHashMapKeyValue_At(map->buckets, idx);
  *value = kvp->value;

  /* Since we are returning the value only, we must free the key and the kvp itself. */
  free(kvp->key);
  free(kvp);

  /* Clear the bucket by setting it to NULL */
  return SUNStlVector_SUNHashMapKeyValue_Set(map->buckets, idx, NULL);
}

/*
  This function allocates a new array the same max_size as the map,
  then it sorts map into a new array of key-value pairs leaving
  the map unchanged.

  **Arguments:**
    * ``map`` -- the ``SUNHashMap`` object to operate on
    * ``sorted`` -- pointer to the sorted array of key-value pairs, this
                    function will allocate the array
    * ``compar`` -- comparator function that is passed to the C standard qsort
                    function

  **Returns:**
    * A SUNErrCode indicating success or a failure
 */
SUNErrCode SUNHashMap_Sort(SUNHashMap map, SUNHashMapKeyValue** sorted,
                           int (*compar)(const void*, const void*))
{
  if (!map || !compar) { return SUN_ERR_ARG_CORRUPT; }

  *sorted =
    (SUNHashMapKeyValue*)malloc(SUNHashMap_Capacity(map) * sizeof(**sorted));
  if (!(*sorted)) { return SUN_ERR_MALLOC_FAIL; }

  /* Copy the buckets into a new array */
  for (int64_t i = 0; i < SUNHashMap_Capacity(map); i++)
  {
    (*sorted)[i] = *SUNStlVector_SUNHashMapKeyValue_At(map->buckets, i);
  }

  qsort(*sorted, SUNHashMap_Capacity(map), sizeof(SUNHashMapKeyValue), compar);

  return SUN_SUCCESS;
}

/*
  This function allocates a new array with just they values of the map.

  **Arguments:**
    * ``map`` -- the ``SUNHashMap`` object to operate on
    * ``values`` -- pointer to the array of keys
    * ``value_size`` -- the size of the values in bytes

  **Returns:**
    * A SUNErrCode indicating success or a failure
 */
SUNDIALS_MAYBE_UNUSED
SUNErrCode SUNHashMap_Values(SUNHashMap map, void*** values, int64_t value_size)
{
  int count = 0;

  if (!map) { return SUN_ERR_ARG_CORRUPT; }

  *values = (void**)malloc(SUNHashMap_Capacity(map) * value_size);
  if (!values) { return SUN_ERR_MALLOC_FAIL; }

  /* Copy the values into a new array */
  for (int64_t i = 0; i < SUNHashMap_Capacity(map); i++)
  {
    SUNHashMapKeyValue kvp = *SUNStlVector_SUNHashMapKeyValue_At(map->buckets, i);
    if (kvp) { (*values)[count++] = kvp->value; }
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNHashMap_PrintKeys(SUNHashMap map, FILE* file)
{
  if (!map) { return SUN_ERR_ARG_CORRUPT; }

  /* Print keys into a new array */
  fprintf(file, "[");
  for (int64_t i = 0; i < SUNHashMap_Capacity(map); i++)
  {
    SUNHashMapKeyValue kvp = *SUNStlVector_SUNHashMapKeyValue_At(map->buckets, i);
    if (kvp) { fprintf(file, "%s, ", kvp->key); }
  }
  fprintf(file, "]\n");

  return SUN_SUCCESS;
}
