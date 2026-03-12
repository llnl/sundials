/*
 * ----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * ----------------------------------------------------------------------------
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
 * ----------------------------------------------------------------------------
 * This is the implementation file for the SuperLU SuperMatrix SLU_NR_loc
 * format compatible SUNMatrix.
 * ----------------------------------------------------------------------------
 */

#include <mpi.h>
#include <stdarg.h>
#include <stdlib.h>

#include <superlu_ddefs.h>

#include <sundials/sundials_math.h>
#include <sundials/sundials_mpi_types.h>
#include <sunmatrix/sunmatrix_slunrloc.h>

#include "sundials_macros.h"

/*
 * ----------------------------------------------------------------------------
 *  Macros for accessing the SUNMatrix_SLUNRloc content structure members
 *  and useful values stored in nested SuperLU structures.
 * ----------------------------------------------------------------------------
 */

#define SM_CONTENT_SLUNRLOC(A) ((SUNMatrixContent_SLUNRloc)(A->content))

#define SM_SUPERMATRIX_SLUNRLOC(A) (SM_CONTENT_SLUNRLOC(A)->A_super)

#define SM_GRID_SLUNRLOC(A) (SM_CONTENT_SLUNRLOC(A)->grid)

#define SM_COMMPATTERN_SLUNRLOC(A) (SM_CONTENT_SLUNRLOC(A)->gsmv_comm)

#define SM_ROWTOPROC_SLUNRLOC(A) (SM_CONTENT_SLUNRLOC(A)->row_to_proc)

#define SM_OWNDATA_SLUNRLOC(A) (SM_CONTENT_SLUNRLOC(A)->own_data)

#define SM_COLSORTED_SLUNRLOC(A) (SM_CONTENT_SLUNRLOC(A)->ACS_super)

#define SM_SUPERSTORE_SLUNRLOC(A) \
  ((NRformat_loc*)(SM_SUPERMATRIX_SLUNRLOC(A)->Store))

#define SM_GLOBALROWS_SLUNRLOC(A) (SM_SUPERMATRIX_SLUNRLOC(A)->nrow)

#define SM_GLOBALCOLS_SLUNRLOC(A) (SM_SUPERMATRIX_SLUNRLOC(A)->ncol)

#define SM_LOCALROWS_SLUNRLOC(A) (SM_SUPERSTORE_SLUNRLOC(A)->m_loc)

#define SM_LOCALNNZ_SLUNRLOC(A) (SM_SUPERSTORE_SLUNRLOC(A)->nnz_loc)

#define SM_FSTROW_SLUNRLOC(A) (SM_SUPERSTORE_SLUNRLOC(A)->fst_row)

/* Private function prototypes */
static sunbooleantype SMCompatible_SLUNRloc(SUNMatrix A, SUNMatrix B);

/*
 * ----------------------------------------------------------------------------
 *  Exported functions
 * ----------------------------------------------------------------------------
 */

SUNMatrix SUNMatrix_SLUNRloc(SuperMatrix* A_super, gridinfo_t* grid,
                             SUNContext sunctx)
{
  SUNMatrix A;
  SUNMatrixContent_SLUNRloc content;

  /* Check for valid inputs */
  if (A_super == NULL || grid == NULL) { return (NULL); }

  if (A_super->Stype != SLU_NR_loc || A_super->Dtype != SLU_D ||
      A_super->Mtype != SLU_GE)
  {
    return (NULL);
  }

  /* Create an empty matrix object */
  A = NULL;
  A = SUNMatNewEmpty(sunctx);
  if (A == NULL) { return (NULL); }

  /* Attach operations */
  A->ops->getid       = SUNMatGetID_SLUNRloc;
  A->ops->clone       = SUNMatClone_SLUNRloc;
  A->ops->destroy     = SUNMatDestroy_SLUNRloc;
  A->ops->zero        = SUNMatZero_SLUNRloc;
  A->ops->copy        = SUNMatCopy_SLUNRloc;
  A->ops->scaleadd    = SUNMatScaleAdd_SLUNRloc;
  A->ops->scaleaddi   = SUNMatScaleAddI_SLUNRloc;
  A->ops->matvecsetup = SUNMatMatvecSetup_SLUNRloc;
  A->ops->matvec      = SUNMatMatvec_SLUNRloc;
  A->ops->space       = SUNMatSpace_SLUNRloc;

  /* Create content */
  content = NULL;
  content =
    (SUNMatrixContent_SLUNRloc)malloc(sizeof(struct _SUNMatrixContent_SLUNRloc));
  if (content == NULL)
  {
    SUNMatDestroy(A);
    return (NULL);
  }

  /* Attach content to SuperMatrix */
  content->A_super     = A_super;
  content->own_data    = SUNFALSE;
  content->grid        = grid;
  content->row_to_proc = NULL;
  content->gsmv_comm   = NULL;
  content->ACS_super   = NULL;

  /* Attach content to SUNMatrix */
  A->content = content;

  return (A);
}

void SUNMatrix_SLUNRloc_Print(SUNMatrix A, FILE* fp)
{
  fprintf(fp, "====== START SUNMatrix_SLUNRloc_Print %p  ======\n", (void*)A);
  fprintf(fp, "A->content->A_super = %p\n", (void*)SM_SUPERMATRIX_SLUNRLOC(A));

  /* Call SuperLU_DIST print routine */
  file_dPrint_CompRowLoc_Matrix_dist(fp, SM_SUPERMATRIX_SLUNRLOC(A));

  fprintf(fp, "======= END SUNMatrix_SLUNRloc_Print %p  =======\n", (void*)A);
}

/*
 * ----------------------------------------------------------------------------
 *  Implementation of accessor functions
 * ----------------------------------------------------------------------------
 */

SuperMatrix* SUNMatrix_SLUNRloc_SuperMatrix(SUNMatrix A)
{
  return (SM_SUPERMATRIX_SLUNRLOC(A));
}

gridinfo_t* SUNMatrix_SLUNRloc_ProcessGrid(SUNMatrix A)
{
  return (SM_GRID_SLUNRLOC(A));
}

sunbooleantype SUNMatrix_SLUNRloc_OwnData(SUNMatrix A)
{
  return (SM_OWNDATA_SLUNRLOC(A));
}

/*
 * ----------------------------------------------------------------------------
 *  Implementation of matrix operations
 * ----------------------------------------------------------------------------
 */

SUNMatrix_ID SUNMatGetID_SLUNRloc(SUNDIALS_MAYBE_UNUSED SUNMatrix A)
{
  return (SUNMATRIX_SLUNRLOC);
}

SUNMatrix SUNMatClone_SLUNRloc(SUNMatrix A)
{
  SUNMatrix B;
  SuperMatrix* B_super;

  /* allocate new SuperMatrix */
  B_super = NULL;
  B_super = malloc(sizeof(SuperMatrix));
  if (B_super == NULL) { return (NULL); }

  /* call the SuperLU-DIST clone function */
  dClone_CompRowLoc_Matrix_dist(SM_SUPERMATRIX_SLUNRLOC(A), B_super);

  /* create the new SUNMatrix */
  B = NULL;
  B = SUNMatrix_SLUNRloc(B_super, SM_GRID_SLUNRLOC(A), A->sunctx);
  if (B == NULL)
  {
    /* call SuperLU-DIST destroy function */
    Destroy_CompRowLoc_Matrix_dist(B_super);
    free(B_super);
    return (NULL);
  }

  /* Allocated the SuperMatrix ourselves, so SUNMatrix now owns the data */
  SM_OWNDATA_SLUNRLOC(B) = SUNTRUE;

  return (B);
}

void SUNMatDestroy_SLUNRloc(SUNMatrix A)
{
  if (!A) { return; }

  if (A->content)
  {
    if (SM_OWNDATA_SLUNRLOC(A))
    {
      /* call SuperLU-DIST destroy function */
      Destroy_CompRowLoc_Matrix_dist(SM_SUPERMATRIX_SLUNRLOC(A));
      free(SM_SUPERMATRIX_SLUNRLOC(A));
    }
    if (SM_COLSORTED_SLUNRLOC(A))
    {
      /* if CS exists, then the Matvec has been initialized and we must finalize
         it by calling the SuperLU DIST pdgsmv_finalize routine to free up
         memory allocated */

      pdgsmv_finalize(SM_COMMPATTERN_SLUNRLOC(A));
      Destroy_CompRowLoc_Matrix_dist(SM_COLSORTED_SLUNRLOC(A));
      free(SM_COLSORTED_SLUNRLOC(A));
      SM_COLSORTED_SLUNRLOC(A) = NULL;
    }
    if (SM_ROWTOPROC_SLUNRLOC(A))
    {
      free(SM_ROWTOPROC_SLUNRLOC(A));
      SM_ROWTOPROC_SLUNRLOC(A) = NULL;
    }
    if (SM_COMMPATTERN_SLUNRLOC(A))
    {
      free(SM_COMMPATTERN_SLUNRLOC(A));
      SM_COMMPATTERN_SLUNRLOC(A) = NULL;
    }
    free(A->content);
    A->content = NULL;
  }

  if (A->ops)
  {
    free(A->ops);
    A->ops = NULL;
  }

  free(A);
  A = NULL;

  return;
}

SUNErrCode SUNMatZero_SLUNRloc(SUNMatrix A)
{
  /* call SuperLU-DIST clone function */
  dZero_CompRowLoc_Matrix_dist(SM_SUPERMATRIX_SLUNRLOC(A));
  return (SUN_SUCCESS);
}

SUNErrCode SUNMatCopy_SLUNRloc(SUNMatrix A, SUNMatrix B)
{
  /* call SuperLU-DIST copy function */
  dCopy_CompRowLoc_Matrix_dist(SM_SUPERMATRIX_SLUNRLOC(A),
                               SM_SUPERMATRIX_SLUNRLOC(B));
  return (SUN_SUCCESS);
}

SUNErrCode SUNMatScaleAdd_SLUNRloc(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  /* check that B can be added into A */
  if (!SMCompatible_SLUNRloc(A, B)) { return (SUN_ERR_ARG_INCOMPATIBLE); }

  /* call SuperLU-DIST ScaleAdd function */
  dScaleAdd_CompRowLoc_Matrix_dist(SM_SUPERMATRIX_SLUNRLOC(A),
                                   SM_SUPERMATRIX_SLUNRLOC(B), c);
  return (SUN_SUCCESS);
}

SUNErrCode SUNMatScaleAddI_SLUNRloc(sunrealtype c, SUNMatrix A)
{
  /* call SuperLU-DIST ScaleAddI function */
  dScaleAddId_CompRowLoc_Matrix_dist(SM_SUPERMATRIX_SLUNRLOC(A), c);
  return (SUN_SUCCESS);
}

SUNErrCode SUNMatMatvec_SLUNRloc(SUNMatrix A, N_Vector x, N_Vector y)
{
  SuperMatrix* ACS;
  sunrealtype *xdata, *ydata;

  /* Extract the column-sorted A */
  ACS = SM_COLSORTED_SLUNRLOC(A);

  /* The column-sorted matrix ACS and the communication pattern must be
     established by calling SUNMatMatvecSetup prior to calling SUNMatMatvec. */
  if (ACS == NULL || SM_ROWTOPROC_SLUNRLOC(A) == NULL ||
      SM_COMMPATTERN_SLUNRLOC(A) == NULL)
  {
    return (SUN_ERR_OP_FAIL);
  }

  xdata = N_VGetArrayPointer(x);
  ydata = N_VGetArrayPointer(y);
  if (xdata == NULL || ydata == NULL) { return (SUN_ERR_MEM_FAIL); }

  /* Call SuperLU-DIST Matvec routine to perform the actual Matvec. */
  pdgsmv(0, ACS, SM_GRID_SLUNRLOC(A), SM_COMMPATTERN_SLUNRLOC(A), xdata, ydata);

  return (SUN_SUCCESS);
}

SUNErrCode SUNMatMatvecSetup_SLUNRloc(SUNMatrix A)
{
  sunindextype* temp;
  sunindextype nprocs;
  sunindextype i, j;

  SuperMatrix* ACS = SM_COLSORTED_SLUNRLOC(A);

  /* If ACS is NULL, then this is the first setup call and things must be
     allocated */
  if (ACS == NULL)
  {
    ACS = (SuperMatrix*)malloc(sizeof(SuperMatrix));

    /* Clone and copy A to create ACS which will be A but with column-sorted
       column indices to [internal, external]. ACS is used with the Matvec
       routine. */
    dClone_CompRowLoc_Matrix_dist(SM_SUPERMATRIX_SLUNRLOC(A), ACS);
    dCopy_CompRowLoc_Matrix_dist(SM_SUPERMATRIX_SLUNRLOC(A), ACS);

    SM_ROWTOPROC_SLUNRLOC(A) =
      (sunindextype*)malloc(SM_GLOBALROWS_SLUNRLOC(A) * sizeof(sunindextype));
    if (SM_ROWTOPROC_SLUNRLOC(A) == NULL)
    {
      Destroy_CompRowLoc_Matrix_dist(ACS);
      free(ACS);
      return (SUN_ERR_MEM_FAIL);
    }

    SM_COMMPATTERN_SLUNRLOC(A) = (pdgsmv_comm_t*)malloc(sizeof(pdgsmv_comm_t));
    if (SM_COMMPATTERN_SLUNRLOC(A) == NULL)
    {
      free(SM_ROWTOPROC_SLUNRLOC(A));
      Destroy_CompRowLoc_Matrix_dist(ACS);
      free(ACS);
      return (SUN_ERR_MEM_FAIL);
    }

    SM_COLSORTED_SLUNRLOC(A) = ACS;
  }
  else
  {
    /* if ACS has already been created, we can reuse it to save allocations,
       but we must finalize the last matvec to avoid leaking memory.*/
    pdgsmv_finalize(SM_COMMPATTERN_SLUNRLOC(A));
  }

  /* calculate the number of processes the matrix is spread across */
  nprocs = SM_GRID_SLUNRLOC(A)->nprow * SM_GRID_SLUNRLOC(A)->npcol;

  /* establish a row number to process mapping */
  temp = (sunindextype*)malloc((nprocs + 1) * sizeof(sunindextype));
  if (temp == NULL)
  {
    SUNMatDestroy(A);
    return (SUN_ERR_MEM_FAIL);
  }

  MPI_Allgather(&SM_FSTROW_SLUNRLOC(A), 1, MPI_SUNINDEXTYPE, temp, 1,
                MPI_SUNINDEXTYPE, SM_GRID_SLUNRLOC(A)->comm);

  temp[nprocs] = SM_GLOBALROWS_SLUNRLOC(A);
  for (i = 0; i < nprocs; i++)
  {
    for (j = temp[i]; j < temp[i + 1]; j++) { SM_ROWTOPROC_SLUNRLOC(A)[j] = i; }
  }
  free(temp);

  /* Call SuperLU-DIST routine to establish communication pattern for the
     Matvec. WARNING: This will overwrite the matrix provided. It is modified
     with colind permuted to [internal, external]. */
  pdgsmv_init(ACS, SM_ROWTOPROC_SLUNRLOC(A), SM_GRID_SLUNRLOC(A),
              SM_COMMPATTERN_SLUNRLOC(A));

  return (SUN_SUCCESS);
}

SUNErrCode SUNMatSpace_SLUNRloc(SUNMatrix A, long int* lenrw, long int* leniw)
{
  /* since the SuperLU_DIST structures are opaque objects, we omit those
     from these results */

  *leniw = SM_GLOBALROWS_SLUNRLOC(A); /* length(row_to_proc) */
  *lenrw = 0;

  return (SUN_SUCCESS);
}

/*
 * ----------------------------------------------------------------------------
 * private functions
 * ----------------------------------------------------------------------------
 */

/* Function to check compatibility of two sparse SUNMatrix objects.
   Checks to make sure that the the matrices are both SLUNRLOC,
   have the same number of rows, cols and nonzeros.*/
static sunbooleantype SMCompatible_SLUNRloc(SUNMatrix A, SUNMatrix B)
{
  /* both matrices must be SLUNRLOC */
  if ((SUNMatGetID(A) != SUNMATRIX_SLUNRLOC) ||
      (SUNMatGetID(B) != SUNMATRIX_SLUNRLOC))
  {
    return SUNFALSE;
  }

  if (SM_GLOBALCOLS_SLUNRLOC(A) != SM_GLOBALCOLS_SLUNRLOC(B))
  {
    return SUNFALSE;
  }

  if (SM_LOCALROWS_SLUNRLOC(A) != SM_LOCALROWS_SLUNRLOC(B)) { return SUNFALSE; }

  if (SM_LOCALNNZ_SLUNRLOC(A) != SM_LOCALNNZ_SLUNRLOC(B)) { return SUNFALSE; }

  return SUNTRUE;
}
