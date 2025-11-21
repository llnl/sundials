/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ UMBC
 *                David Gardner, Carol Woodward, Slaven Peles,
 *                Cody Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
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
 * This is the header file for a generic matrix package.
 * It defines the SUNMatrix structure (_generic_SUNMatrix) which
 * contains the following fields:
 *   - an implementation-dependent 'content' field which contains
 *     the description and actual data of the matrix
 *   - an 'ops' filed which contains a structure listing operations
 *     acting on such matrices
 * -----------------------------------------------------------------
 * This header file contains:
 *   - enumeration constants for all SUNDIALS-defined matrix types,
 *     as well as a generic type for user-supplied matrix types,
 *   - type declarations for the _generic_SUNMatrix and
 *     _generic_SUNMatrix_Ops structures, as well as references to
 *     pointers to such structures (SUNMatrix), and
 *   - prototypes for the matrix functions which operate on
 *     SUNMatrix objects.
 * -----------------------------------------------------------------
 * At a minimum, a particular implementation of a SUNMatrix must
 * do the following:
 *   - specify the 'content' field of SUNMatrix,
 *   - implement the operations on those SUNMatrix objects,
 *   - provide a constructor routine for new SUNMatrix objects
 *
 * Additionally, a SUNMatrix implementation may provide the following:
 *   - macros to access the underlying SUNMatrix data
 *   - a routine to print the content of a SUNMatrix
 * -----------------------------------------------------------------*/

#ifndef _SUNMATRIX_H
#define _SUNMATRIX_H

#include <sundials/sundials_config.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------------------
 * Implemented SUNMatrix types
 * ----------------------------------------------------------------- */

enum SUNMatrix_ID
{
  SUNMATRIX_DENSE,
  SUNMATRIX_MAGMADENSE,
  SUNMATRIX_ONEMKLDENSE,
  SUNMATRIX_BAND,
  SUNMATRIX_SPARSE,
  SUNMATRIX_SLUNRLOC,
  SUNMATRIX_CUSPARSE,
  SUNMATRIX_GINKGO,
  SUNMATRIX_GINKGOBATCH,
  SUNMATRIX_KOKKOSDENSE,
  SUNMATRIX_CUSTOM
};

#ifndef SWIG
typedef enum SUNMatrix_ID SUNMatrix_ID;
#endif

/* -----------------------------------------------------------------
 * Generic definition of SUNMatrix
 * ----------------------------------------------------------------- */

/* Forward reference for pointer to SUNMatrix_Ops object */
typedef _SUNDIALS_STRUCT_ _generic_SUNMatrix_Ops* SUNMatrix_Ops;

/* Forward reference for pointer to SUNMatrix object */
typedef _SUNDIALS_STRUCT_ _generic_SUNMatrix* SUNMatrix;

/* Structure containing function pointers to matrix operations  */
struct _generic_SUNMatrix_Ops
{
  SUNMatrix_ID (*getid)(SUNMatrix);
  SUNMatrix (*clone)(SUNMatrix);
  void (*destroy)(SUNMatrix);
  SUNErrCode (*zero)(SUNMatrix);
  SUNErrCode (*copy)(SUNMatrix, SUNMatrix);
  SUNErrCode (*scaleadd)(sunrealtype, SUNMatrix, SUNMatrix);
  SUNErrCode (*scaleaddi)(sunrealtype, SUNMatrix);
  SUNErrCode (*matvecsetup)(SUNMatrix);
  SUNErrCode (*matvec)(SUNMatrix, N_Vector, N_Vector);
  SUNErrCode (*mathermitiantransposevec)(SUNMatrix, N_Vector, N_Vector);
  SUNErrCode (*space)(SUNMatrix, long int*, long int*);
};

/* A matrix is a structure with an implementation-dependent
   'content' field, and a pointer to a structure of matrix
   operations corresponding to that implementation.  */
struct _generic_SUNMatrix
{
  void* content;
  SUNMatrix_Ops ops;
  SUNContext sunctx;
};

/* -----------------------------------------------------------------
 * Functions exported by SUNMatrix module
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT
SUNMatrix SUNMatNewEmpty(SUNContext sunctx);

SUNDIALS_EXPORT
void SUNMatFreeEmpty(SUNMatrix A);

SUNDIALS_EXPORT
SUNErrCode SUNMatCopyOps(SUNMatrix A, SUNMatrix B);

SUNDIALS_EXPORT
SUNMatrix_ID SUNMatGetID(SUNMatrix A);

SUNDIALS_EXPORT
SUNMatrix SUNMatClone(SUNMatrix A);

SUNDIALS_EXPORT
void SUNMatDestroy(SUNMatrix A);

SUNDIALS_EXPORT
SUNErrCode SUNMatZero(SUNMatrix A);

SUNDIALS_EXPORT
SUNErrCode SUNMatCopy(SUNMatrix A, SUNMatrix B);

SUNDIALS_EXPORT
SUNErrCode SUNMatScaleAdd(sunrealtype c, SUNMatrix A, SUNMatrix B);

SUNDIALS_EXPORT
SUNErrCode SUNMatScaleAddI(sunrealtype c, SUNMatrix A);

SUNDIALS_EXPORT
SUNErrCode SUNMatMatvecSetup(SUNMatrix A);

SUNDIALS_EXPORT
SUNErrCode SUNMatMatvec(SUNMatrix A, N_Vector x, N_Vector y);

SUNDIALS_EXPORT
SUNErrCode SUNMatHermitianTransposeVec(SUNMatrix A, N_Vector x, N_Vector y);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
SUNErrCode SUNMatSpace(SUNMatrix A, long int* lenrw, long int* leniw);

#ifdef __cplusplus
}
#endif

#endif /* _SUNMATRIX_H */
