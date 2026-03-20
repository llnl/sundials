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
 * Implementation of a resizable container similar to a std::vector.
 * The values can be anything but data must be contiguous.
 *
 * To use the StlVector, first define TTYPE with your data type
 * before including this header. The name of the class for your data
 * type will then be SUNStlVector_TTYPE and functions will be
 * SUNStlVector_TTYPE_<function>. E.g.
 *   #define TTYPE int
 *   #include "sunstl_vector.h'
 *   #undef TTYPE
 *   SUNStlVector_int_New(10, destroyIntFn);
 * If you need StlVectors that hold different types in the same file,
 * then define TTYPE for the first, include this header, undefine
 * TTYPE then repeat.
 * -----------------------------------------------------------------*/

#include <stdlib.h>
#include <sundials/sundials_core.h>

#ifndef TTYPE
#error "Must define template type for SUNStlVector"
#endif

#define CONCAT(a, b)            a##b
#define PASTE(a, b)             CONCAT(a, b)
#define MAKE_NAME(prefix, name) PASTE(prefix, PASTE(_, name))

#define SUNStlVectorTtype_s MAKE_NAME(SUNStlVector, PASTE(TTYPE, _s))
#define SUNStlVectorTtype   MAKE_NAME(SUNStlVector, TTYPE)

typedef struct SUNStlVectorTtype_s* SUNStlVectorTtype;

struct SUNStlVectorTtype_s
{
  int64_t size;
  int64_t capacity;
  TTYPE* values;
  SUNErrCode (*destroyValue)(TTYPE*);
};

// This constant controls how much space will be allocated when a resize is needed.
// The new capacity is SUNSTLVECTOR_GROWTH_FACTOR*current_capacity.
// Some std::vector implementations use 2, but 1.5 will be more conservative in terms
// of the memory usage but yields a larger constant factor in terms of the
// amortized constant time complexity.
#define SUNSTLVECTOR_GROWTH_FACTOR 1.5L

/**
 * Creates a new SUNStlVector with the specified initial capacity.
 *
 * :param init_capacity: Initial capacity of the vector.
 * :param destroyValue: Function pointer to destroy the value.
 * :return: New vector instance or NULL on failure.
 */
static inline SUNStlVectorTtype MAKE_NAME(SUNStlVectorTtype,
                                          New)(int64_t init_capacity,
                                               SUNErrCode (*destroyValue)(TTYPE*))
{
  if (init_capacity < 0 || !destroyValue) { return NULL; }
  SUNStlVectorTtype self = (SUNStlVectorTtype)malloc(sizeof(*self));
  if (!self) { return NULL; }
  self->values = (TTYPE*)malloc(sizeof(TTYPE) * init_capacity);
  if (!(self->values)) { return NULL; }
  self->size         = 0;
  self->capacity     = init_capacity;
  self->destroyValue = destroyValue;
  return self;
}

/**
 * Checks if the vector is empty.
 *
 * :param self: Pointer to the vector.
 * :return: True if the vector is empty, false otherwise.
 */
static inline sunbooleantype MAKE_NAME(SUNStlVectorTtype,
                                       IsEmpty)(SUNStlVectorTtype self)
{
  return self->size == 0;
}

/**
 * Allocates more memory (capacity) for the vector.
 *
 * :param self: Pointer to the vector.
 * :param new_capacity: New capacity to reserve.
 * :return: SUN_SUCCESS on success, SUN_ERR_MALLOC_FAIL on failure.
 */
static inline SUNErrCode MAKE_NAME(SUNStlVectorTtype,
                                   Reserve)(SUNStlVectorTtype self,
                                            int64_t new_capacity)
{
  if (new_capacity <= self->capacity) { return SUN_SUCCESS; }
  TTYPE* new_values = (TTYPE*)realloc(self->values, sizeof(TTYPE) * new_capacity);
  if (!new_values) { return SUN_ERR_MALLOC_FAIL; }
  self->values   = new_values;
  self->capacity = new_capacity;
  return SUN_SUCCESS;
}

/**
 * Grows the vector capacity if needed.
 *
 * :param self: Pointer to the vector.
 * :return: SUN_SUCCESS on success, SUN_ERR_MALLOC_FAIL on failure.
 */
static inline SUNErrCode MAKE_NAME(SUNStlVectorTtype, Grow)(SUNStlVectorTtype self)
{
  if (self->size == self->capacity)
  {
    /* It is possible, although unlikely, that new_capacity overflows a long double.
       We explicitly cast capacity to a long double to silence any implicit
       conversion compiler warning. */
    int64_t new_capacity = self->capacity == 0
                             ? 2
                             : (int64_t)(ceill(((long double)self->capacity) *
                                               SUNSTLVECTOR_GROWTH_FACTOR));

    return MAKE_NAME(SUNStlVectorTtype, Reserve)(self, new_capacity);
  }
  return SUN_SUCCESS;
}

/**
 * Adds an element to the end of the vector.
 *
 * :param self: Pointer to the vector.
 * :param element: Element to add.
 * :return: SUN_SUCCESS on success, SUN_ERR_MALLOC_FAIL on failure.
 */
static inline SUNErrCode MAKE_NAME(SUNStlVectorTtype,
                                   PushBack)(SUNStlVectorTtype self, TTYPE element)
{
  if (self->size == self->capacity)
  {
    SUNErrCode err = MAKE_NAME(SUNStlVectorTtype, Grow)(self);
    if (err != SUN_SUCCESS) { return err; }
  }
  self->values[self->size++] = element;
  return SUN_SUCCESS;
}

/**
 * Returns a pointer to the element at the specified index.
 *
 * :param self: Pointer to the vector.
 * :param index: Index of the element.
 * :return: Pointer to the element or NULL if out of bounds.
 */
static inline TTYPE* MAKE_NAME(SUNStlVectorTtype, At)(SUNStlVectorTtype self,
                                                      int64_t index)
{
  if (index >= self->size || index < 0)
  {
    // Handle index out of bounds
    return NULL;
  }
  return &(self->values[index]);
}

/**
 * Sets the element at the specified index.
 *
 * :param self: Pointer to the vector.
 * :param index: Index of the element.
 * :param element: Element to set.
 * :return: SUN_SUCCESS on success, SUN_ERR_OUTOFRANGE if index is out of bounds.
 */
static inline SUNErrCode MAKE_NAME(SUNStlVectorTtype,
                                   Set)(SUNStlVectorTtype self, int64_t index,
                                        TTYPE element)
{
  if (index >= self->size || index < 0)
  {
    // Handle index out of bounds
    return SUN_ERR_OUTOFRANGE;
  }
  self->values[index] = element;
  return SUN_SUCCESS;
}

/**
 * Removes the last element from the vector.
 *
 * :param self: Pointer to the vector.
 * :return: SUN_SUCCESS on success.
 */
static inline SUNErrCode MAKE_NAME(SUNStlVectorTtype,
                                   PopBack)(SUNStlVectorTtype self)
{
  /* `static` results in implicit empty initialization in C99. */
  static TTYPE nullish;
  if (self->size == 0) return SUN_SUCCESS;
  SUNErrCode err = MAKE_NAME(SUNStlVectorTtype, Set)(self, self->size - 1,
                                                     nullish);
  if (err) { return err; }
  self->size--;
  return SUN_SUCCESS;
}

/**
 * Removes the element at the specified index and shifts all other elements.
 *
 * :param self: Pointer to the vector.
 * :param index: Index of the element to erase.
 * :return: SUN_SUCCESS on success, SUN_ERR_OUTOFRANGE if index is out of bounds.
 */
static inline SUNErrCode MAKE_NAME(SUNStlVectorTtype,
                                   Erase)(SUNStlVectorTtype self, int64_t index)
{
  static TTYPE nullish;
  if (self->size == 0) return SUN_SUCCESS;
  if (index >= self->size || index < 0) return SUN_ERR_OUTOFRANGE;

  SUNErrCode err = MAKE_NAME(SUNStlVectorTtype, Set)(self, index, nullish);
  if (err != SUN_SUCCESS) return err;

  for (int64_t i = index; i < self->size - 1; i++)
  {
    self->values[i] = self->values[i + 1];
  }

  err = MAKE_NAME(SUNStlVectorTtype, Set)(self, self->size - 1, nullish);
  if (err != SUN_SUCCESS) return err;

  self->size -= 1;
  return SUN_SUCCESS;
}

/**
 * Returns the size of the vector.
 *
 * :param self: Pointer to the vector.
 * :return: Size of the vector.
 */
static inline int64_t MAKE_NAME(SUNStlVectorTtype, Size)(SUNStlVectorTtype self)
{
  return self->size;
}

/**
 * Returns the capacity of the vector.
 *
 * :param self: Pointer to the vector.
 * :return: Capacity of the vector.
 */
static inline int64_t MAKE_NAME(SUNStlVectorTtype,
                                Capacity)(SUNStlVectorTtype self)
{
  return self->capacity;
}

/**
 * Destroys the vector and frees its memory.
 *
 * :param self_ptr: Pointer to the vector pointer.
 * :return: SUN_SUCCESS on success.
 */
static inline SUNErrCode MAKE_NAME(SUNStlVectorTtype,
                                   Destroy)(SUNStlVectorTtype* self_ptr)
{
  static TTYPE nullish;

  if (!self_ptr || !(*self_ptr)) return SUN_SUCCESS;

  SUNStlVectorTtype self = *self_ptr;

  for (int64_t i = 0; i < MAKE_NAME(SUNStlVectorTtype, Size)(self); i++)
  {
    SUNErrCode err = self->destroyValue(&(self->values[i]));
    if (err) { return err; }
    self->values[i] = nullish;
  }

  free(self->values);
  free(self);
  *self_ptr = NULL;

  return SUN_SUCCESS;
}
