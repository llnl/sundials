/*
 * -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * This is the header file for the dense implementation of the
 * SUNMATRIX module, SUNMATRIX_MAGMADENSE.
 * -----------------------------------------------------------------
 */

#ifndef _SUNMATRIX_MAGMADENSE_H
#define _SUNMATRIX_MAGMADENSE_H

#include <stdio.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_memory.h>

#if defined(SUNDIALS_MAGMA_BACKENDS_CUDA)
#define HAVE_CUBLAS
#elif defined(SUNDIALS_MAGMA_BACKENDS_HIP)
#define HAVE_HIP
#endif
#include <magma_v2.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

struct _SUNMatrixContent_MagmaDense
{
  int last_flag;           /* last error code returned by magma  */
  int device_id;           /* device ID used by magma            */
  sunindextype M;          /* number of rows in block            */
  sunindextype N;          /* number of columns in block         */
  sunindextype nblocks;    /* number of blocks in matrix         */
  sunindextype ldata;      /* length of data array               */
  SUNMemory data;          /* matrix data; column-major          */
  SUNMemory blocks;        /* device pointers to blocks of A     */
  SUNMemory xblocks;       /* device pointers to blocks of x     */
  SUNMemory yblocks;       /* device pointers to blocks of y     */
  SUNMemoryHelper memhelp; /* memory helper                      */
  magma_queue_t q;         /* operation queue (i.e. stream)      */
};

typedef struct _SUNMatrixContent_MagmaDense* SUNMatrixContent_MagmaDense;

/* ---------------------------------------
 * Implementation specific functions
 * ---------------------------------------*/

SUNDIALS_EXPORT SUNMatrix SUNMatrix_MagmaDense(sunindextype M, sunindextype N,
                                               SUNMemoryType memtype,
                                               SUNMemoryHelper memhelper,
                                               void* queue, SUNContext sunctx);
SUNDIALS_EXPORT SUNMatrix SUNMatrix_MagmaDenseBlock(
  sunindextype nblocks, sunindextype M, sunindextype N, SUNMemoryType memtype,
  SUNMemoryHelper memhelper, void* queue, SUNContext sunctx);
SUNDIALS_EXPORT void SUNMatrix_MagmaDense_Print(SUNMatrix A);
SUNDIALS_EXPORT sunrealtype* SUNMatrix_MagmaDense_Data(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_MagmaDense_LData(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_MagmaDense_Rows(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_MagmaDense_Columns(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_MagmaDense_BlockRows(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_MagmaDense_BlockColumns(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_MagmaDense_BlockLData(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_MagmaDense_NumBlocks(SUNMatrix A);
SUNDIALS_EXPORT sunrealtype** SUNMatrix_MagmaDense_BlockData(SUNMatrix A);
SUNDIALS_EXPORT SUNErrCode SUNMatrix_MagmaDense_CopyToDevice(SUNMatrix A,
                                                             sunrealtype* h_data);
SUNDIALS_EXPORT
SUNErrCode SUNMatrix_MagmaDense_CopyFromDevice(SUNMatrix A, sunrealtype* h_data);

static inline sunrealtype* SUNMatrix_MagmaDense_Block(SUNMatrix Amat,
                                                      sunindextype k)
{
  SUNMatrixContent_MagmaDense A = (SUNMatrixContent_MagmaDense)Amat->content;
  return (((sunrealtype*)A->data->ptr) + k * A->M * A->N);
}

static inline sunrealtype* SUNMatrix_MagmaDense_Column(SUNMatrix Amat,
                                                       sunindextype j)
{
  SUNMatrixContent_MagmaDense A = (SUNMatrixContent_MagmaDense)Amat->content;
  return (((sunrealtype*)A->data->ptr) + j * A->M);
}

static inline sunrealtype* SUNMatrix_MagmaDense_BlockColumn(SUNMatrix Amat,
                                                            sunindextype k,
                                                            sunindextype j)
{
  SUNMatrixContent_MagmaDense A = (SUNMatrixContent_MagmaDense)Amat->content;
  return (((sunrealtype*)A->data->ptr) + k * A->M * A->N + j * A->M);
}

/* ---------------------------------------
 * SUNMatrix API functions
 * ---------------------------------------*/

static inline SUNMatrix_ID SUNMatGetID_MagmaDense(SUNMatrix A)
{
  return SUNMATRIX_MAGMADENSE;
}

SUNDIALS_EXPORT SUNMatrix SUNMatClone_MagmaDense(SUNMatrix A);
SUNDIALS_EXPORT void SUNMatDestroy_MagmaDense(SUNMatrix A);
SUNDIALS_EXPORT SUNErrCode SUNMatZero_MagmaDense(SUNMatrix A);
SUNDIALS_EXPORT SUNErrCode SUNMatCopy_MagmaDense(SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT SUNErrCode SUNMatScaleAdd_MagmaDense(sunrealtype c, SUNMatrix A,
                                                     SUNMatrix B);
SUNDIALS_EXPORT SUNErrCode SUNMatScaleAddI_MagmaDense(sunrealtype c, SUNMatrix A);
SUNDIALS_EXPORT SUNErrCode SUNMatMatvecSetup_MagmaDense(SUNMatrix A);
SUNDIALS_EXPORT SUNErrCode SUNMatMatvec_MagmaDense(SUNMatrix A, N_Vector x,
                                                   N_Vector y);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
SUNErrCode SUNMatSpace_MagmaDense(SUNMatrix A, long int* lenrw, long int* leniw);

#ifdef __cplusplus
}
#endif

#endif
