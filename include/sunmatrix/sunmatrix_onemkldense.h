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
 * This is the header file for the dense implementation of the SUNMATRIX
 * class using the Intel oneAPI Math Kernel Library (oneMKL).
 * ---------------------------------------------------------------------------*/

#ifndef _SUNMATRIX_ONEMKLDENSE_H
#define _SUNMATRIX_ONEMKLDENSE_H

#include <stdio.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_memory.h>
#include <sundials/sundials_sycl_policies.hpp>
#include <sycl/sycl.hpp>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

struct _SUNMatrixContent_OneMklDense
{
  int last_flag;                  /* last error code returned       */
  sunindextype block_rows;        /* number of rows in a block      */
  sunindextype block_cols;        /* number of columns in a block   */
  sunindextype num_blocks;        /* number of blocks in the matrix */
  sunindextype rows;              /* total number of rows           */
  sunindextype cols;              /* total number of columns        */
  sunindextype ldata;             /* length of data array           */
  SUNMemory data;                 /* matrix data; column-major      */
  SUNMemory blocks;               /* device pointers to blocks of A */
  SUNSyclExecPolicy* exec_policy; /* execution policy               */
  SUNMemoryType mem_type;         /* memory type                    */
  SUNMemoryHelper mem_helper;     /* memory helper                  */
  ::sycl::queue* queue;           /* operation queue                */
};

typedef struct _SUNMatrixContent_OneMklDense* SUNMatrixContent_OneMklDense;

/* ---------------------------------------------------------------------------
 * Implementation specific functions
 * ---------------------------------------------------------------------------*/

/* Constructors */

SUNDIALS_EXPORT
SUNMatrix SUNMatrix_OneMklDense(sunindextype M, sunindextype N,
                                SUNMemoryType mem_type,
                                SUNMemoryHelper mem_helper,
                                ::sycl::queue* queue, SUNContext sunctx);

SUNDIALS_EXPORT
SUNMatrix SUNMatrix_OneMklDenseBlock(sunindextype num_blocks,
                                     sunindextype M_block, sunindextype N_block,
                                     SUNMemoryType mem_type,
                                     SUNMemoryHelper mem_helper,
                                     ::sycl::queue* queue, SUNContext sunctx);

/* Get matrix dimensions */

SUNDIALS_EXPORT
sunindextype SUNMatrix_OneMklDense_Rows(SUNMatrix A);

SUNDIALS_EXPORT
sunindextype SUNMatrix_OneMklDense_Columns(SUNMatrix A);

/* Get matrix block dimensions */

SUNDIALS_EXPORT
sunindextype SUNMatrix_OneMklDense_NumBlocks(SUNMatrix A);

SUNDIALS_EXPORT
sunindextype SUNMatrix_OneMklDense_BlockRows(SUNMatrix A);

SUNDIALS_EXPORT
sunindextype SUNMatrix_OneMklDense_BlockColumns(SUNMatrix A);

/* Get matrix data */

SUNDIALS_EXPORT
sunindextype SUNMatrix_OneMklDense_LData(SUNMatrix A);

SUNDIALS_EXPORT
sunrealtype* SUNMatrix_OneMklDense_Data(SUNMatrix A);

static inline sunrealtype* SUNMatrix_OneMklDense_Column(SUNMatrix Amat,
                                                        sunindextype j)
{
  SUNMatrixContent_OneMklDense A = (SUNMatrixContent_OneMklDense)Amat->content;
  return (((sunrealtype*)A->data->ptr) + j * A->block_rows);
}

/* Get matrix block data */

SUNDIALS_EXPORT
sunindextype SUNMatrix_OneMklDense_BlockLData(SUNMatrix A);

SUNDIALS_EXPORT
sunrealtype** SUNMatrix_OneMklDense_BlockData(SUNMatrix A);

static inline sunrealtype* SUNMatrix_OneMklDense_Block(SUNMatrix Amat,
                                                       sunindextype k)
{
  SUNMatrixContent_OneMklDense A = (SUNMatrixContent_OneMklDense)Amat->content;
  return (((sunrealtype*)A->data->ptr) + k * A->block_rows * A->block_cols);
}

static inline sunrealtype* SUNMatrix_OneMklDense_BlockColumn(SUNMatrix Amat,
                                                             sunindextype k,
                                                             sunindextype j)
{
  SUNMatrixContent_OneMklDense A = (SUNMatrixContent_OneMklDense)Amat->content;
  return (((sunrealtype*)A->data->ptr) + k * A->block_rows * A->block_cols +
          j * A->block_rows);
}

/* Copy data */

SUNDIALS_EXPORT
SUNErrCode SUNMatrix_OneMklDense_CopyToDevice(SUNMatrix A, sunrealtype* h_data);

SUNDIALS_EXPORT
SUNErrCode SUNMatrix_OneMklDense_CopyFromDevice(SUNMatrix A, sunrealtype* h_data);

/* ---------------------------------------------------------------------------
 * SUNMatrix API functions
 * ---------------------------------------------------------------------------*/

static inline SUNMatrix_ID SUNMatGetID_OneMklDense(SUNMatrix A)
{
  return SUNMATRIX_ONEMKLDENSE;
}

SUNDIALS_EXPORT
SUNMatrix SUNMatClone_OneMklDense(SUNMatrix A);

SUNDIALS_EXPORT
void SUNMatDestroy_OneMklDense(SUNMatrix A);

SUNDIALS_EXPORT
SUNErrCode SUNMatZero_OneMklDense(SUNMatrix A);

SUNDIALS_EXPORT
SUNErrCode SUNMatCopy_OneMklDense(SUNMatrix A, SUNMatrix B);

SUNDIALS_EXPORT
SUNErrCode SUNMatScaleAdd_OneMklDense(sunrealtype c, SUNMatrix A, SUNMatrix B);

SUNDIALS_EXPORT
SUNErrCode SUNMatScaleAddI_OneMklDense(sunrealtype c, SUNMatrix A);

SUNDIALS_EXPORT
SUNErrCode SUNMatMatvecSetup_OneMklDense(SUNMatrix A);

SUNDIALS_EXPORT
SUNErrCode SUNMatMatvec_OneMklDense(SUNMatrix A, N_Vector x, N_Vector y);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
SUNErrCode SUNMatSpace_OneMklDense(SUNMatrix A, long int* lenrw, long int* leniw);

#ifdef __cplusplus
}
#endif

#endif
