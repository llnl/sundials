/* ---------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ---------------------------------------------------------------------------
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
 * ---------------------------------------------------------------------------
 * This is the implementation file for the dense implementation of the
 * SUNLINEARSOLVER class using the Intel oneAPI Math Kernel Library (oneMKL).
 * ---------------------------------------------------------------------------*/

#include <oneapi/mkl/lapack.hpp>
#include <stdio.h>
#include <stdlib.h>
using namespace oneapi::mkl::lapack;

// SUNDIALS public headers
#include <sunlinsol/sunlinsol_onemkldense.h>
#include <sunmatrix/sunmatrix_onemkldense.h>

// SUNDIALS private headers
#include "sundials_debug.h"

// Check for a valid precision and index size
#if defined(SUNDIALS_EXTENDED_PRECISION)
#error "oneMLK unsupported precision"
#endif

#if defined(SUNDIALS_INT32_T)
#error "oneMLK unsupported index size"
#endif

// Accessor macros

// Content and last error flag
#define LS_CONTENT(S)  ((SUNLinearSolverContent_OneMklDense)(S->content))
#define LS_LASTFLAG(S) (LS_CONTENT(S)->last_flag)

// Pivots array length and memory
#define LS_ROWS(S)    (LS_CONTENT(S)->rows)
#define LS_PIVOTS(S)  (LS_CONTENT(S)->pivots)
#define LS_PIVOTSp(S) ((sunindextype*)LS_CONTENT(S)->pivots->ptr)

// Getrf scratch space size and memory
#define LS_F_SCRATCH_SIZE(S) (LS_CONTENT(S)->f_scratch_size)
#define LS_F_SCRATCH(S)      (LS_CONTENT(S)->f_scratchpad)
#define LS_F_SCRATCHp(S)     ((sunrealtype*)LS_CONTENT(S)->f_scratchpad->ptr)

// Getrs scratch space size and memory
#define LS_S_SCRATCH_SIZE(S) (LS_CONTENT(S)->s_scratch_size)
#define LS_S_SCRATCH(S)      (LS_CONTENT(S)->s_scratchpad)
#define LS_S_SCRATCHp(S)     ((sunrealtype*)LS_CONTENT(S)->s_scratchpad->ptr)

// Memory type, helper, and SYCL queue
#define LS_MEM_TYPE(S)   (LS_CONTENT(S)->mem_type)
#define LS_MEM_HELPER(S) (LS_CONTENT(S)->mem_helper)
#define LS_QUEUE(S)      (LS_CONTENT(S)->queue)

/* --------------------------------------------------------------------------
 * Constructors
 * -------------------------------------------------------------------------- */

SUNLinearSolver SUNLinSol_OneMklDense(N_Vector y, SUNMatrix Amat,
                                      SUNContext sunctx)
{
  int retval = 0;

  // Check inputs
  if (!y || !Amat)
  {
    SUNDIALS_DEBUG_ERROR("Illegal input, y or A is NULL\n");
    return NULL;
  }

  if (!(y->ops) || !(Amat->ops))
  {
    SUNDIALS_DEBUG_ERROR("Illegal input, y->ops or A->ops is NULL\n");
    return NULL;
  }

  if (!(y->ops->nvgetlength) || !(y->ops->nvgetdevicearraypointer) ||
      !(Amat->ops->getid))
  {
    SUNDIALS_DEBUG_ERROR("Illegal input, y or A missing required operations\n");
    return NULL;
  }

  // Check compatibility with supplied SUNMatrix
  if (SUNMatGetID(Amat) != SUNMATRIX_ONEMKLDENSE)
  {
    SUNDIALS_DEBUG_ERROR("Illegal input, SUNMatID != SUNMATRIX_ONEMKLDENSE\n");
    return NULL;
  }

  if (!(Amat->content))
  {
    SUNDIALS_DEBUG_ERROR("Illegal input, SUNMatID != SUNMATRIX_ONEMKLDENSE\n");
    return NULL;
  }

  SUNMatrixContent_OneMklDense A = (SUNMatrixContent_OneMklDense)Amat->content;

  // Check that the matrix is square
  if (A->rows != A->cols)
  {
    SUNDIALS_DEBUG_ERROR("Illegal input, A is not square\n");
    return NULL;
  }

  // Check that the matrix and vector dimensions agree
  if (A->cols != N_VGetLength(y))
  {
    SUNDIALS_DEBUG_ERROR(
      "Illegal input, number of columns in A != length of y\n");
    return NULL;
  }

  // Create the linear solver
  SUNLinearSolver S = SUNLinSolNewEmpty(sunctx);
  if (!S)
  {
    SUNDIALS_DEBUG_ERROR("SUNLinSolNewEmpty returned NULL\n");
    return NULL;
  }

  // Attach operations
  S->ops->gettype    = SUNLinSolGetType_OneMklDense;
  S->ops->getid      = SUNLinSolGetID_OneMklDense;
  S->ops->initialize = SUNLinSolInitialize_OneMklDense;
  S->ops->setup      = SUNLinSolSetup_OneMklDense;
  S->ops->solve      = SUNLinSolSolve_OneMklDense;
  S->ops->lastflag   = SUNLinSolLastFlag_OneMklDense;
  S->ops->space      = SUNLinSolSpace_OneMklDense;
  S->ops->free       = SUNLinSolFree_OneMklDense;

  // Create content
  S->content = (SUNLinearSolverContent_OneMklDense)malloc(
    sizeof(_SUNLinearSolverContent_OneMklDense));
  if (!(S->content))
  {
    SUNDIALS_DEBUG_ERROR("Content allocation failed\n");
    SUNLinSolFree(S);
    return NULL;
  }

  // Fill content
  LS_CONTENT(S)->last_flag      = 0;
  LS_CONTENT(S)->rows           = A->rows;
  LS_CONTENT(S)->pivots         = NULL;
  LS_CONTENT(S)->f_scratch_size = 0;
  LS_CONTENT(S)->f_scratchpad   = NULL;
  LS_CONTENT(S)->s_scratch_size = 0;
  LS_CONTENT(S)->s_scratchpad   = NULL;
  LS_CONTENT(S)->mem_type       = A->mem_type;
  LS_CONTENT(S)->mem_helper     = A->mem_helper;
  LS_CONTENT(S)->queue          = A->queue;

  // Allocate data
  retval = SUNMemoryHelper_Alloc(LS_MEM_HELPER(S), &(LS_PIVOTS(S)),
                                 A->rows * sizeof(sunindextype), LS_MEM_TYPE(S),
                                 A->queue);
  if (retval)
  {
    SUNDIALS_DEBUG_ERROR("Pivots allocation failed\n");
    SUNLinSolFree(S);
    return NULL;
  }

  // Compute scratchpad size for factorization and solve
  ::sycl::queue* queue    = A->queue;
  sunindextype M          = SUNMatrix_OneMklDense_BlockRows(Amat);
  sunindextype N          = SUNMatrix_OneMklDense_BlockColumns(Amat);
  sunindextype num_blocks = SUNMatrix_OneMklDense_NumBlocks(Amat);

  if (num_blocks > 1)
  {
#ifdef SUNDIALS_ONEMKL_USE_GETRF_LOOP
    LS_F_SCRATCH_SIZE(S) =
      getrf_scratchpad_size<sunrealtype>(*queue, // device queue
                                         M,      // rows in A_i
                                         N,      // columns in A_i
                                         M);     // leading dimension
#else
    LS_F_SCRATCH_SIZE(S) =
      getrf_batch_scratchpad_size<sunrealtype>(*queue, // device queue
                                               M,      // rows in A_i
                                               N,      // columns in A_i
                                               M,      // leading dimension
                                               M * N,  // stride between A_i
                                               M,      // stride in P_i
                                               num_blocks); // number of blocks
#endif

#ifdef SUNDIALS_ONEMKL_USE_GETRS_LOOP
    LS_S_SCRATCH_SIZE(S) =
      getrs_scratchpad_size<sunrealtype>(*queue, // device queue
                                         oneapi::mkl::transpose::nontrans,
                                         M,  // number of rows in A
                                         1,  // number of right-hand sizes
                                         M,  // leading dimension of A
                                         M); // leading dimension of B
#else
    LS_S_SCRATCH_SIZE(S) =
      getrs_batch_scratchpad_size<sunrealtype>(*queue, // device queue
                                               oneapi::mkl::transpose::nontrans,
                                               M, // number of rows in A_i
                                               1, // number of right-hand sides
                                               M, // leading dimension of A_i
                                               M * N, // stride between A_i
                                               M,     // stride between pivots
                                               M, // leading dimension of B_i
                                               M, // stride between B_i
                                               num_blocks); // number of blocks
#endif
  }
  else
  {
    LS_F_SCRATCH_SIZE(S) =
      getrf_scratchpad_size<sunrealtype>(*queue, // device queue
                                         M,      // rows in A_i
                                         N,      // columns in A_i
                                         M);     // leading dimension

    LS_S_SCRATCH_SIZE(S) =
      getrs_scratchpad_size<sunrealtype>(*queue, // device queue
                                         oneapi::mkl::transpose::nontrans,
                                         M,  // number of rows in A
                                         1,  // number of right-hand sizes
                                         M,  // leading dimension of A
                                         M); // leading dimension of B
  }

  // Allocate factorization scratchpad if necessary
  retval = SUNMemoryHelper_Alloc(LS_MEM_HELPER(S), &(LS_F_SCRATCH(S)),
                                 LS_F_SCRATCH_SIZE(S) * sizeof(sunrealtype),
                                 LS_MEM_TYPE(S), queue);
  if (retval)
  {
    SUNDIALS_DEBUG_ERROR("Scratchpad allocation failed\n");
    SUNLinSolFree(S);
    return NULL;
  }

  // Allocate solve scratchpad if necessary
  retval = SUNMemoryHelper_Alloc(LS_MEM_HELPER(S), &(LS_S_SCRATCH(S)),
                                 LS_S_SCRATCH_SIZE(S) * sizeof(sunrealtype),
                                 LS_MEM_TYPE(S), queue);
  if (retval)
  {
    SUNDIALS_DEBUG_ERROR("Scratchpad allocation failed\n");
    SUNLinSolFree(S);
    return NULL;
  }

  return S;
}

/* --------------------------------------------------------------------------
 * Implementation of SUNLinearSolver operations
 * -------------------------------------------------------------------------- */

SUNErrCode SUNLinSolInitialize_OneMklDense(SUNLinearSolver S)
{
  // All solver-specific memory has already been allocated
  if (!S)
  {
    SUNDIALS_DEBUG_ERROR("Linear solver is NULL\n");
    return SUN_ERR_ARG_CORRUPT;
  }

  LS_LASTFLAG(S) = SUN_SUCCESS;
  return SUN_SUCCESS;
}

int SUNLinSolSetup_OneMklDense(SUNLinearSolver S, SUNMatrix A)
{
  // Check for valid inputs
  if (!S)
  {
    SUNDIALS_DEBUG_ERROR("Linear solver is NULL\n");
    return SUN_ERR_ARG_CORRUPT;
  }

  if (!A)
  {
    SUNDIALS_DEBUG_ERROR("Matrix is NULL\n");
    LS_LASTFLAG(S) = SUN_ERR_ARG_CORRUPT;
    return SUN_ERR_ARG_CORRUPT;
  }

  // Ensure that A is a oneMKL dense matrix
  if (SUNMatGetID(A) != SUNMATRIX_ONEMKLDENSE)
  {
    SUNDIALS_DEBUG_ERROR("Matrix is not the oneMKL matrix\n");
    LS_LASTFLAG(S) = SUN_ERR_ARG_INCOMPATIBLE;
    return SUN_ERR_ARG_INCOMPATIBLE;
  }

  // Access A matrix data array
  sunrealtype* Adata = SUNMatrix_OneMklDense_Data(A);
  if (!Adata)
  {
    SUNDIALS_DEBUG_ERROR("Matrix data array is NULL\n");
    LS_LASTFLAG(S) = SUN_ERR_MEM_FAIL;
    return SUN_ERR_MEM_FAIL;
  }

  // Access pivots data array
  sunindextype* pivots = LS_PIVOTSp(S);
  if (!pivots)
  {
    SUNDIALS_DEBUG_ERROR("Matrix data array is NULL\n");
    LS_LASTFLAG(S) = SUN_ERR_MEM_FAIL;
    return SUN_ERR_MEM_FAIL;
  }

  // Call oneMKL to do LU factorization of A
  ::sycl::queue* queue      = LS_QUEUE(S);
  sunindextype ier          = 0;
  sunindextype M            = SUNMatrix_OneMklDense_BlockRows(A);
  sunindextype N            = SUNMatrix_OneMklDense_BlockColumns(A);
  sunindextype num_blocks   = SUNMatrix_OneMklDense_NumBlocks(A);
  sunindextype scratch_size = LS_F_SCRATCH_SIZE(S);
  sunrealtype* scratchpad   = LS_F_SCRATCHp(S);

  if (num_blocks > 1)
  {
#ifdef SUNDIALS_ONEMKL_USE_GETRF_LOOP
    try
    {
      for (sunindextype i = 0; i < num_blocks; i++)
      {
        getrf(*queue,            // device queue
              M,                 // number of rows
              N,                 // number of columns
              Adata + i * M * N, // matrix data
              M,                 // leading dimension of A
              pivots + i * M,    // array of pivots
              scratchpad,        // scratchpad memory
              scratch_size);     // scratchpad size
      }
    }
    catch (oneapi::mkl::lapack::exception const& e)
    {
      SUNDIALS_DEBUG_ERROR("An exception occurred in getrf\n");
      if (e.info())
      {
        // An illegal value was providied or the scratch pad is too small
        ier = -1;
      }
      else
      {
        // The diagonal element of some of U_i is zero
        ier = 1;
      }
    }
#else
    try
    {
      getrf_batch(*queue,        // device queue
                  M,             // number of block rows
                  N,             // number of block columns
                  Adata,         // matrix data
                  M,             // leading dimension of A
                  M * N,         // stride between A_i
                  pivots,        // array of pivots
                  M,             // stride between P_i
                  num_blocks,    // number of blocks
                  scratchpad,    // scratchpad memory
                  scratch_size); // scratchpad size
    }
    catch (oneapi::mkl::lapack::exception const& e)
    {
      SUNDIALS_DEBUG_ERROR("An exception occurred in getrf_batch\n");
      if (e.info())
      {
        // An illegal value was providied or the scratch pad is too small
        ier = -1;
      }
      else
      {
        // The diagonal element of some of U_i is zero
        ier = 1;
      }
    }
#endif
  }
  else
  {
    try
    {
      getrf(*queue,        // device queue
            M,             // number of rows
            N,             // number of columns
            Adata,         // matrix data
            M,             // leading dimension of A
            pivots,        // array of pivots
            scratchpad,    // scratchpad memory
            scratch_size); // scratchpad size
    }
    catch (oneapi::mkl::lapack::exception const& e)
    {
      SUNDIALS_DEBUG_ERROR("An exception occurred in getrf\n");
      if (e.info())
      {
        // An illegal value was providied or the scratch pad is too small
        ier = -1;
      }
      else
      {
        // The diagonal element of some of U_i is zero
        ier = 1;
      }
    }
  }

  if (ier > 0)
  {
    LS_LASTFLAG(S) = ier;
    return SUNLS_LUFACT_FAIL;
  }

  if (ier < 0)
  {
    LS_LASTFLAG(S) = ier;
    return SUN_ERR_EXT_FAIL;
  }

  LS_LASTFLAG(S) = SUN_SUCCESS;
  return SUN_SUCCESS;
}

int SUNLinSolSolve_OneMklDense(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                               N_Vector b, sunrealtype tol)
{
  // Check for valid inputs
  if (!S)
  {
    SUNDIALS_DEBUG_ERROR("Linear solver is NULL\n");
    return SUN_ERR_ARG_CORRUPT;
  }

  if (!A || !x || !b)
  {
    SUNDIALS_DEBUG_ERROR("A, x, or b is NULL\n");
    LS_LASTFLAG(S) = SUN_ERR_ARG_CORRUPT;
    return SUN_ERR_ARG_CORRUPT;
  }

  // Ensure that A is a onemkl dense matrix
  if (SUNMatGetID(A) != SUNMATRIX_ONEMKLDENSE)
  {
    SUNDIALS_DEBUG_ERROR("Matrix is not the oneMKL matrix\n");
    LS_LASTFLAG(S) = SUN_ERR_ARG_INCOMPATIBLE;
    return SUN_ERR_ARG_INCOMPATIBLE;
  }

  // Copy b into x
  N_VScale(SUN_RCONST(1.0), b, x);

  // Access x vector data array
  sunrealtype* xdata = N_VGetDeviceArrayPointer(x);
  if (!xdata)
  {
    SUNDIALS_DEBUG_ERROR("Vector data array is NULL\n");
    LS_LASTFLAG(S) = SUN_ERR_MEM_FAIL;
    return SUN_ERR_MEM_FAIL;
  }

  // Access A matrix data array
  sunrealtype* Adata = SUNMatrix_OneMklDense_Data(A);
  if (!Adata)
  {
    SUNDIALS_DEBUG_ERROR("Matrix data array is NULL\n");
    LS_LASTFLAG(S) = SUN_ERR_MEM_FAIL;
    return SUN_ERR_MEM_FAIL;
  }

  // Access pivots data array
  sunindextype* pivots = LS_PIVOTSp(S);
  if (!pivots)
  {
    SUNDIALS_DEBUG_ERROR("Matrix data array is NULL\n");
    LS_LASTFLAG(S) = SUN_ERR_MEM_FAIL;
    return SUN_ERR_MEM_FAIL;
  }

  // Call oneMKL to solve the linear system
  sunindextype ier          = 0;
  ::sycl::queue* queue      = LS_QUEUE(S);
  sunindextype M            = SUNMatrix_OneMklDense_BlockRows(A);
  sunindextype N            = SUNMatrix_OneMklDense_BlockColumns(A);
  sunindextype num_blocks   = SUNMatrix_OneMklDense_NumBlocks(A);
  sunindextype scratch_size = LS_S_SCRATCH_SIZE(S);
  sunrealtype* scratchpad   = LS_S_SCRATCHp(S);

  if (num_blocks > 1)
  {
#ifdef SUNDIALS_ONEMKL_USE_GETRS_LOOP
    try
    {
      for (sunindextype i = 0; i < num_blocks; i++)
      {
        getrs(*queue, // device queue
              oneapi::mkl::transpose::nontrans,
              M,                 // number of rows
              1,                 // number of right-hand sides
              Adata + i * M * N, // factorized matrix data
              M,                 // leading dimension of A
              pivots,            // array of pivots
              xdata + i * M,     // right-hand side data
              M,                 // leading dimension of B_i
              scratchpad,        // scratchpad memory
              scratch_size);     // scratchpad size
      }
    }
    catch (oneapi::mkl::lapack::exception const& e)
    {
      SUNDIALS_DEBUG_ERROR("An exception occurred in getrs\n");
      ier = -1;
    }
#else
    try
    {
      getrs_batch(*queue, // device queue
                  oneapi::mkl::transpose::nontrans,
                  M,             // number of rows
                  1,             // number of right-hand sides
                  Adata,         // factorized matrix data
                  M,             // leading dimension of A_i
                  M * N,         // stride between A_i
                  pivots,        // array of pivots
                  M,             // stride between pivots
                  xdata,         // right-hand side data
                  M,             // leading dimension of B_i
                  M,             // stride between B_i
                  num_blocks,    // number of blocks
                  scratchpad,    // scratchpad memory
                  scratch_size); // scratchpad size
    }
    catch (oneapi::mkl::lapack::exception const& e)
    {
      SUNDIALS_DEBUG_ERROR("An exception occurred in getrs_batch\n");
      ier = -1;
    }
#endif
  }
  else
  {
    try
    {
      getrs(*queue, // device queue
            oneapi::mkl::transpose::nontrans,
            M,             // number of rows
            1,             // number of right-hand sides
            Adata,         // factorized matrix data
            M,             // leading dimension of A
            pivots,        // array of pivots
            xdata,         // right-hand side data
            M,             // leading dimension of B_i
            scratchpad,    // scratchpad memory
            scratch_size); // scratchpad size
    }
    catch (oneapi::mkl::lapack::exception const& e)
    {
      SUNDIALS_DEBUG_ERROR("An exception occurred in getrs\n");
      ier = -1;
    }
  }

  if (ier < 0)
  {
    LS_LASTFLAG(S) = ier;
    return SUN_ERR_EXT_FAIL;
  }

  LS_LASTFLAG(S) = SUN_SUCCESS;
  return SUN_SUCCESS;
}

sunindextype SUNLinSolLastFlag_OneMklDense(SUNLinearSolver S)
{
  return LS_LASTFLAG(S);
}

SUNErrCode SUNLinSolSpace_OneMklDense(SUNLinearSolver S, long int* lenrwLS,
                                      long int* leniwLS)
{
  if (!S)
  {
    SUNDIALS_DEBUG_ERROR("Linear solver is NULL\n");
    return SUN_ERR_ARG_CORRUPT;
  }

  *lenrwLS = 0;
  *leniwLS = 2 + LS_CONTENT(S)->rows;

  LS_LASTFLAG(S) = SUN_SUCCESS;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolFree_OneMklDense(SUNLinearSolver S)
{
  // return if S is already free
  if (!S) { return SUN_SUCCESS; }

  // delete items from contents, then delete generic structure
  if (S->content)
  {
    // Pivots memory
    if (LS_PIVOTS(S))
    {
      SUNMemoryHelper_Dealloc(LS_MEM_HELPER(S), LS_PIVOTS(S), LS_QUEUE(S));
    }

    // Factorization scratch memory
    if (LS_F_SCRATCH(S))
    {
      SUNMemoryHelper_Dealloc(LS_MEM_HELPER(S), LS_F_SCRATCH(S), LS_QUEUE(S));
    }
    LS_F_SCRATCH_SIZE(S) = 0;

    // Solve scratch memory
    if (LS_S_SCRATCH(S))
    {
      SUNMemoryHelper_Dealloc(LS_MEM_HELPER(S), LS_S_SCRATCH(S), LS_QUEUE(S));
    }
    LS_S_SCRATCH_SIZE(S) = 0;
  }

  SUNLinSolFreeEmpty(S);
  S = NULL;

  return SUN_SUCCESS;
}
