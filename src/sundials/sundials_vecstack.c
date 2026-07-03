/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * Stack of temporary N_Vectors for reuse.
 * -----------------------------------------------------------------*/

#include <stdlib.h>

#include <sundials/sundials_errors.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_vecstack.h>

#define TTYPE N_Vector
#include "stl/sunstl_vector.h"
#undef TTYPE

struct SUNVecStack_
{
  N_Vector tmpl;
  SUNStlVector_N_Vector vecs;
  int64_t num_owned;
  int64_t num_checked_out;
};

static SUNErrCode SUNVecStack_DestroyValue(N_Vector* vec)
{
  if (vec && *vec)
  {
    N_VDestroy(*vec);
    *vec = NULL;
  }
  return SUN_SUCCESS;
}

SUNErrCode SUNVecStack_Create(N_Vector tmpl, SUNVecStack* stack_out)
{
  if (!tmpl || !stack_out) { return SUN_ERR_ARG_CORRUPT; }

  *stack_out = NULL;

  SUNVecStack stack = (SUNVecStack)malloc(sizeof(*stack));
  if (!stack) { return SUN_ERR_MALLOC_FAIL; }

  stack->tmpl            = tmpl;
  stack->num_owned       = 0;
  stack->num_checked_out = 0;
  stack->vecs            = SUNStlVector_N_Vector_New(0, SUNVecStack_DestroyValue);

  if (!stack->vecs)
  {
    free(stack);
    return SUN_ERR_MALLOC_FAIL;
  }

  *stack_out = stack;
  return SUN_SUCCESS;
}

SUNErrCode SUNVecStack_Destroy(SUNVecStack* stack_in)
{
  if (stack_in == NULL || *stack_in == NULL) { return SUN_SUCCESS; }

  SUNVecStack stack = *stack_in;
  if (stack->num_checked_out < 0 || stack->num_owned < 0)
  {
    return SUN_ERR_CORRUPT;
  }
  if (stack->num_checked_out != 0) { return SUN_ERR_MEM_FAIL; }

  SUNErrCode err = SUNStlVector_N_Vector_Destroy(&stack->vecs);
  if (err != SUN_SUCCESS) { return err; }

  free(stack);
  *stack_in = NULL;
  return SUN_SUCCESS;
}

SUNErrCode SUNVecStack_Pop(SUNVecStack stack, N_Vector* vec_out)
{
  if (!stack || !vec_out) { return SUN_ERR_ARG_CORRUPT; }

  *vec_out = NULL;

  N_Vector* cached = SUNStlVector_N_Vector_Back(stack->vecs);
  if (cached)
  {
    *vec_out = *cached;
    SUNErrCode err = SUNStlVector_N_Vector_PopBack(stack->vecs);
    if (err != SUN_SUCCESS) { return err; }
  }
  else
  {
    *vec_out = N_VClone(stack->tmpl);
    if (!*vec_out) { return SUN_ERR_MALLOC_FAIL; }
    stack->num_owned++;
  }

  stack->num_checked_out++;
  return SUN_SUCCESS;
}

SUNErrCode SUNVecStack_Push(SUNVecStack stack, N_Vector* vec_in)
{
  if (!stack || !vec_in || !*vec_in) { return SUN_ERR_ARG_CORRUPT; }

  int64_t idle = SUNStlVector_N_Vector_Size(stack->vecs);
  if (idle < 0) { return SUN_ERR_CORRUPT; }
  if (stack->num_owned < 0 || stack->num_checked_out < 0)
  {
    return SUN_ERR_CORRUPT;
  }
  if (stack->num_checked_out <= 0) { return SUN_ERR_MEM_FAIL; }
  if (stack->num_checked_out > stack->num_owned) { return SUN_ERR_CORRUPT; }
  if (idle + stack->num_checked_out != stack->num_owned)
  {
    return SUN_ERR_CORRUPT;
  }
  if (idle >= stack->num_owned) { return SUN_ERR_MEM_FAIL; }

  SUNErrCode err = SUNStlVector_N_Vector_PushBack(stack->vecs, *vec_in);
  if (err != SUN_SUCCESS) { return err; }

  *vec_in = NULL;
  stack->num_checked_out--;

  idle = SUNStlVector_N_Vector_Size(stack->vecs);
  if (idle + stack->num_checked_out != stack->num_owned)
  {
    return SUN_ERR_CORRUPT;
  }

  return SUN_SUCCESS;
}
