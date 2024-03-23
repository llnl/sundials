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
  SUNContext sunctx;
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

SUNErrCode SUNVecStack_Create(N_Vector tmpl, int init_size, SUNContext sunctx,
                              SUNVecStack* stack_out)
{
  SUNFunctionBegin(sunctx);
  SUNAssert(tmpl, SUN_ERR_ARG_CORRUPT);
  SUNAssert(init_size, SUN_ERR_ARG_OUTOFRANGE);
  SUNAssert(stack_out, SUN_ERR_ARG_CORRUPT);

  *stack_out = NULL;

  SUNVecStack stack = (SUNVecStack)malloc(sizeof(*stack));
  if (!stack) { return SUN_ERR_MALLOC_FAIL; }

  stack->sunctx          = sunctx;
  stack->tmpl            = tmpl;
  stack->num_owned       = 0;
  stack->num_checked_out = 0;
  stack->vecs            = SUNStlVector_N_Vector_New(init_size, SUNVecStack_DestroyValue);

  if (!stack->vecs)
  {
    SUNVecStack_Destroy(&stack);
    return SUN_ERR_MALLOC_FAIL;
  }

  for (int i = 0; i < init_size; i++)
  {
    N_Vector vec = N_VClone(tmpl);
    if (vec == NULL)
    {
      SUNVecStack_Destroy(&stack);
      return SUN_ERR_MALLOC_FAIL;
    }
    stack->num_owned++;
    stack->num_checked_out++;

    SUNErrCode err = SUNVecStack_Push(stack, &vec);
    if (err != SUN_SUCCESS)
    {
      SUNVecStack_Destroy(&stack);
      return SUN_ERR_MALLOC_FAIL;
    }
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

  SUNCheckCall(SUNStlVector_N_Vector_Destroy(&stack->vecs));

  free(stack);
  *stack_in = NULL;
  return SUN_SUCCESS;
}

SUNErrCode SUNVecStack_Pop(SUNVecStack stack, N_Vector* vec_out)
{
  SUNFunctionBegin(stack->sunctx);
  SUNAssert(vec_out, SUN_ERR_ARG_CORRUPT);

  *vec_out = NULL;

  N_Vector* cached = SUNStlVector_N_Vector_Back(stack->vecs);
  if (cached)
  {
    *vec_out = *cached;
    SUNCheckCall(SUNStlVector_N_Vector_PopBack(stack->vecs));
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
  SUNFunctionBegin(stack->sunctx);
  SUNAssert(vec_in, SUN_ERR_ARG_CORRUPT);
  SUNAssert(*vec_in, SUN_ERR_ARG_CORRUPT);

  SUNAssert(stack->num_owned > 0, SUN_ERR_CORRUPT);
  SUNAssert(stack->num_checked_out > 0, SUN_ERR_CORRUPT);
  SUNAssert(stack->num_checked_out <= stack->num_owned, SUN_ERR_CORRUPT);

  SUNAssert(SUNStlVector_N_Vector_Size(stack->vecs) + stack->num_checked_out ==
              stack->num_owned,
            SUN_ERR_CORRUPT);

  SUNCheckCall(SUNStlVector_N_Vector_PushBack(stack->vecs, *vec_in));
  *vec_in = NULL;
  stack->num_checked_out--;

  SUNAssert(SUNStlVector_N_Vector_Size(stack->vecs) + stack->num_checked_out ==
              stack->num_owned,
            SUN_ERR_CORRUPT);

  return SUN_SUCCESS;
}

SUNErrCode SUNVecStack_GetNumVecs(SUNVecStack stack, int64_t* num_vecs)
{
  SUNFunctionBegin(stack->sunctx);
  SUNAssert(num_vecs, SUN_ERR_ARG_CORRUPT);

  *num_vecs = stack->num_owned;

  return SUN_SUCCESS;
}

SUNErrCode SUNVecStack_GetNumActiveVecs(SUNVecStack stack, int64_t* num_active_vecs)
{
  SUNFunctionBegin(stack->sunctx);
  SUNAssert(num_active_vecs, SUN_ERR_ARG_CORRUPT);

  *num_active_vecs = stack->num_checked_out;

  return SUN_SUCCESS;
}

SUNErrCode SUNVecStack_GetNumIdleVecs(SUNVecStack stack, int64_t* num_idle_vecs)
{
  SUNFunctionBegin(stack->sunctx);
  SUNAssert(num_idle_vecs, SUN_ERR_ARG_CORRUPT);

  *num_idle_vecs = SUNStlVector_N_Vector_Size(stack->vecs);

  return SUN_SUCCESS;
}
