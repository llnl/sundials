/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * Based on N_Vector_Parallel by Scott D. Cohen, Alan C. Hindmarsh,
 * Radu Serban, and Aaron Collier @ LLNL
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
 * This is the implementation file for a PETSc implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_petsc_deprecated.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_math.h>

#include "sundials_macros.h"

#define ZERO   SUN_RCONST(0.0)
#define HALF   SUN_RCONST(0.5)
#define ONE    SUN_RCONST(1.0)
#define ONEPT5 SUN_RCONST(1.5)

/* Error Message */
#define BAD_N1 "N_VNewEmpty_Petsc -- Sum of local vector lengths differs from "
#define BAD_N2 "input global length. \n\n"
#define BAD_N  BAD_N1 BAD_N2

/*
 * -----------------------------------------------------------------
 * Simplifying macros NV_CONTENT_PTC, NV_OWN_DATA_PTC,
 *                    NV_LOCLENGTH_PTC, NV_GLOBLENGTH_PTC,
 *                    NV_COMM_PTC
 * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations
 * are assumed:
 *
 * N_Vector v;
 * sunindextype v_len, s_len, i;
 *
 * (1) NV_CONTENT_PTC
 *
 *     This routines gives access to the contents of the PETSc
 *     vector wrapper N_Vector.
 *
 *     The assignment v_cont = NV_CONTENT_PTC(v) sets v_cont to be
 *     a pointer to the N_Vector (PETSc wrapper) content structure.
 *
 * (2) NV_PVEC_PTC, NV_OWN_DATA_PTC, NV_LOCLENGTH_PTC, NV_GLOBLENGTH_PTC,
 *     and NV_COMM_PTC
 *
 *     These routines give access to the individual parts of
 *     the content structure of a PETSc N_Vector wrapper.
 *
 *     NV_PVEC_PTC(v) returns the PETSc vector (Vec) object.
 *
 *     The assignment v_llen = NV_LOCLENGTH_PTC(v) sets v_llen to
 *     be the length of the local part of the vector v. The call
 *     NV_LOCLENGTH_PTC(v) = llen_v should NOT be used! It will
 *     change the value stored in the N_Vector content structure,
 *     but it will NOT change the length of the actual PETSc vector.
 *
 *     The assignment v_glen = NV_GLOBLENGTH_PTC(v) sets v_glen to
 *     be the global length of the vector v. The call
 *     NV_GLOBLENGTH_PTC(v) = glen_v should NOT be used! It will
 *     change the value stored in the N_Vector content structure,
 *     but it will NOT change the length of the actual PETSc vector.
 *
 *     The assignment v_comm = NV_COMM_PTC(v) sets v_comm to be the
 *     MPI communicator of the vector v. The assignment
 *     NV_COMM_PTC(v) = comm_v should NOT be used! It will change
 *     the value stored in the N_Vector content structure, but it
 *     will NOT change the MPI communicator of the actual PETSc
 *     vector.
 *
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_PTC(v) ((N_VectorContent_Petsc)(v->content))

#define NV_LOCLENGTH_PTC(v) (NV_CONTENT_PTC(v)->local_length)

#define NV_GLOBLENGTH_PTC(v) (NV_CONTENT_PTC(v)->global_length)

#define NV_OWN_DATA_PTC(v) (NV_CONTENT_PTC(v)->own_data)

#define NV_PVEC_PTC(v) (NV_CONTENT_PTC(v)->pvec)

#define NV_COMM_PTC(v) (NV_CONTENT_PTC(v)->comm)

/* Functions attached to the N_Vector */
static void nvAbs_Petsc(N_Vector x, N_Vector z);
static void nvAddConst_Petsc(N_Vector x, sunrealtype b, N_Vector z);
static SUNErrCode nvBufPack_Petsc(N_Vector x, void* buf);
static SUNErrCode nvBufSize_Petsc(N_Vector x, sunindextype* size);
static SUNErrCode nvBufUnpack_Petsc(N_Vector x, void* buf);
static N_Vector nvCloneEmpty_Petsc(N_Vector w);
static N_Vector nvClone_Petsc(N_Vector w);
static void nvCompare_Petsc(sunrealtype c, N_Vector x, N_Vector z);
static SUNErrCode nvConstVectorArray_Petsc(int nvecs, sunrealtype c, N_Vector* Z);
static void nvConst_Petsc(sunrealtype c, N_Vector z);
static sunbooleantype nvConstrMaskLocal_Petsc(N_Vector c, N_Vector x, N_Vector m);
static sunbooleantype nvConstrMask_Petsc(N_Vector c, N_Vector x, N_Vector m);
static void nvDestroy_Petsc(N_Vector v);
static void nvDiv_Petsc(N_Vector x, N_Vector y, N_Vector z);
static sunrealtype nvDotProdLocal_Petsc(N_Vector x, N_Vector y);
static SUNErrCode nvDotProdMultiAllReduce_Petsc(int nvec, N_Vector x,
                                                sunrealtype* sum);
static SUNErrCode nvDotProdMultiLocal_Petsc(int nvec, N_Vector x, N_Vector* Y,
                                            sunrealtype* dotprods);
static SUNErrCode nvDotProdMulti_Petsc(int nvec, N_Vector x, N_Vector* Y,
                                       sunrealtype* dotprods);
static sunrealtype nvDotProd_Petsc(N_Vector x, N_Vector y);
static MPI_Comm nvGetCommunicator_Petsc(N_Vector v);
static sunindextype nvGetLength_Petsc(N_Vector v);
static N_Vector_ID nvGetVectorID_Petsc(N_Vector v);
static sunbooleantype nvInvTestLocal_Petsc(N_Vector x, N_Vector z);
static sunbooleantype nvInvTest_Petsc(N_Vector x, N_Vector z);
static void nvInv_Petsc(N_Vector x, N_Vector z);
static sunrealtype nvL1NormLocal_Petsc(N_Vector x);
static sunrealtype nvL1Norm_Petsc(N_Vector x);
static SUNErrCode nvLinearCombinationVectorArray_Petsc(int nvec, int nsum,
                                                       sunrealtype* c,
                                                       N_Vector** X, N_Vector* Z);
static SUNErrCode nvLinearCombination_Petsc(int nvec, sunrealtype* c,
                                            N_Vector* X, N_Vector z);
static SUNErrCode nvLinearSumVectorArray_Petsc(int nvec, sunrealtype a,
                                               N_Vector* X, sunrealtype b,
                                               N_Vector* Y, N_Vector* Z);
static void nvLinearSum_Petsc(sunrealtype a, N_Vector x, sunrealtype b,
                              N_Vector y, N_Vector z);
static sunrealtype nvMaxNormLocal_Petsc(N_Vector x);
static sunrealtype nvMaxNorm_Petsc(N_Vector x);
static sunrealtype nvMinLocal_Petsc(N_Vector x);
static sunrealtype nvMinQuotientLocal_Petsc(N_Vector num, N_Vector denom);
static sunrealtype nvMinQuotient_Petsc(N_Vector num, N_Vector denom);
static sunrealtype nvMin_Petsc(N_Vector x);
static void nvPrintFile_Petsc(N_Vector x, const char fname[]);
static void nvPrint_Petsc(N_Vector x);
static void nvProd_Petsc(N_Vector x, N_Vector y, N_Vector z);
static SUNErrCode nvScaleAddMultiVectorArray_Petsc(int nvec, int nsum,
                                                   sunrealtype* a, N_Vector* X,
                                                   N_Vector** Y, N_Vector** Z);
static SUNErrCode nvScaleAddMulti_Petsc(int nvec, sunrealtype* a, N_Vector x,
                                        N_Vector* Y, N_Vector* Z);
static SUNErrCode nvScaleVectorArray_Petsc(int nvec, sunrealtype* c,
                                           N_Vector* X, N_Vector* Z);
static void nvScale_Petsc(sunrealtype c, N_Vector x, N_Vector z);
static sunrealtype nvWL2Norm_Petsc(N_Vector x, N_Vector w);
static sunrealtype nvWSqrSumLocal_Petsc(N_Vector x, N_Vector w);
static sunrealtype nvWSqrSumMaskLocal_Petsc(N_Vector x, N_Vector w, N_Vector id);
static SUNErrCode nvWrmsNormMaskVectorArray_Petsc(int nvec, N_Vector* X,
                                                  N_Vector* W, N_Vector id,
                                                  sunrealtype* nrm);
static sunrealtype nvWrmsNormMask_Petsc(N_Vector x, N_Vector w, N_Vector id);
static SUNErrCode nvWrmsNormVectorArray_Petsc(int nvecs, N_Vector* X,
                                              N_Vector* W, sunrealtype* nrm);
static sunrealtype nvWrmsNorm_Petsc(N_Vector x, N_Vector w);

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID nvGetVectorID_Petsc(SUNDIALS_MAYBE_UNUSED N_Vector v)
{
  return SUNDIALS_NVEC_PETSC;
}

/* ----------------------------------------------------------------
 * Function to create a new N_Vector wrapper with an empty (NULL)
 * PETSc vector.
 */

N_Vector N_VNewEmpty_Petsc(MPI_Comm comm, sunindextype local_length,
                           sunindextype global_length, SUNContext sunctx)
{
  N_Vector v;
  N_VectorContent_Petsc content;
  sunindextype n, Nsum;
  PetscErrorCode ierr;

  /* Compute global length as sum of local lengths */
  n    = local_length;
  ierr = MPI_Allreduce(&n, &Nsum, 1, MPI_SUNINDEXTYPE, MPI_SUM, comm);
  CHKERRABORT(comm, ierr);
  if (Nsum != global_length)
  {
    fprintf(stderr, BAD_N);
    return (NULL);
  }

  /* Create an empty vector object */
  v = NULL;
  v = N_VNewEmpty(sunctx);
  if (v == NULL) { return (NULL); }

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = nvGetVectorID_Petsc;
  v->ops->nvclone           = nvClone_Petsc;
  v->ops->nvcloneempty      = nvCloneEmpty_Petsc;
  v->ops->nvdestroy         = nvDestroy_Petsc;
  v->ops->nvgetcommunicator = nvGetCommunicator_Petsc;
  v->ops->nvgetlength       = nvGetLength_Petsc;

  /* standard vector operations */
  v->ops->nvlinearsum    = nvLinearSum_Petsc;
  v->ops->nvconst        = nvConst_Petsc;
  v->ops->nvprod         = nvProd_Petsc;
  v->ops->nvdiv          = nvDiv_Petsc;
  v->ops->nvscale        = nvScale_Petsc;
  v->ops->nvabs          = nvAbs_Petsc;
  v->ops->nvinv          = nvInv_Petsc;
  v->ops->nvaddconst     = nvAddConst_Petsc;
  v->ops->nvdotprod      = nvDotProd_Petsc;
  v->ops->nvmaxnorm      = nvMaxNorm_Petsc;
  v->ops->nvwrmsnormmask = nvWrmsNormMask_Petsc;
  v->ops->nvwrmsnorm     = nvWrmsNorm_Petsc;
  v->ops->nvmin          = nvMin_Petsc;
  v->ops->nvwl2norm      = nvWL2Norm_Petsc;
  v->ops->nvl1norm       = nvL1Norm_Petsc;
  v->ops->nvcompare      = nvCompare_Petsc;
  v->ops->nvinvtest      = nvInvTest_Petsc;
  v->ops->nvconstrmask   = nvConstrMask_Petsc;
  v->ops->nvminquotient  = nvMinQuotient_Petsc;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvdotprodlocal     = nvDotProdLocal_Petsc;
  v->ops->nvmaxnormlocal     = nvMaxNormLocal_Petsc;
  v->ops->nvminlocal         = nvMinLocal_Petsc;
  v->ops->nvl1normlocal      = nvL1NormLocal_Petsc;
  v->ops->nvinvtestlocal     = nvInvTestLocal_Petsc;
  v->ops->nvconstrmasklocal  = nvConstrMaskLocal_Petsc;
  v->ops->nvminquotientlocal = nvMinQuotientLocal_Petsc;
  v->ops->nvwsqrsumlocal     = nvWSqrSumLocal_Petsc;
  v->ops->nvwsqrsummasklocal = nvWSqrSumMaskLocal_Petsc;

  /* single buffer reduction operations */
  v->ops->nvdotprodmultilocal     = nvDotProdMultiLocal_Petsc;
  v->ops->nvdotprodmultiallreduce = nvDotProdMultiAllReduce_Petsc;

  /* XBraid interface operations */
  v->ops->nvbufsize   = nvBufSize_Petsc;
  v->ops->nvbufpack   = nvBufPack_Petsc;
  v->ops->nvbufunpack = nvBufUnpack_Petsc;

  /* debugging functions */
  v->ops->nvprint = nvPrint_Petsc;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Petsc)malloc(sizeof *content);
  if (content == NULL)
  {
    N_VDestroy(v);
    return (NULL);
  }

  /* Attach content */
  v->content = content;

  /* Initialize content */
  content->local_length  = local_length;
  content->global_length = global_length;
  content->comm          = comm;
  content->own_data      = SUNFALSE;
  content->pvec          = NULL;

  return (v);
}

/* ----------------------------------------------------------------
 * Function to create an N_Vector wrapper for a PETSc vector.
 */

N_Vector N_VMake_Petsc(Vec pvec, SUNContext sunctx)
{
  N_Vector v = NULL;
  MPI_Comm comm;
  PetscInt local_length;
  PetscInt global_length;

  VecGetLocalSize(pvec, &local_length);
  VecGetSize(pvec, &global_length);
  PetscObjectGetComm((PetscObject)pvec, &comm);

  v = N_VNewEmpty_Petsc(comm, local_length, global_length, sunctx);
  if (v == NULL) { return (NULL); }

  /* Attach data */
  NV_OWN_DATA_PTC(v) = SUNFALSE;
  NV_PVEC_PTC(v)     = pvec;

  return (v);
}

/* ----------------------------------------------------------------
 * Function to extract PETSc vector
 */

Vec N_VGetVector_Petsc(N_Vector v) { return NV_PVEC_PTC(v); }

/* ----------------------------------------------------------------
 * Function to set the PETSc vector
 */

void N_VSetVector_Petsc(N_Vector v, Vec p) { NV_PVEC_PTC(v) = p; }

/* ----------------------------------------------------------------
 * Function to print the global data in a PETSc vector to stdout
 */

void nvPrint_Petsc(N_Vector x)
{
  Vec xv        = NV_PVEC_PTC(x);
  MPI_Comm comm = NV_COMM_PTC(x);

  VecView(xv, PETSC_VIEWER_STDOUT_(comm));

  return;
}

/* ----------------------------------------------------------------
 * Function to print the global data in a PETSc vector to fname
 */

void nvPrintFile_Petsc(N_Vector x, const char fname[])
{
  Vec xv        = NV_PVEC_PTC(x);
  MPI_Comm comm = NV_COMM_PTC(x);
  PetscViewer viewer;

  PetscViewerASCIIOpen(comm, fname, &viewer);

  VecView(xv, viewer);

  PetscViewerDestroy(&viewer);

  return;
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector nvCloneEmpty_Petsc(N_Vector w)
{
  N_Vector v;
  N_VectorContent_Petsc content;

  if (w == NULL) { return (NULL); }

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty(w->sunctx);
  if (v == NULL) { return (NULL); }

  /* Attach operations */
  if (N_VCopyOps(w, v))
  {
    N_VDestroy(v);
    return (NULL);
  }

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Petsc)malloc(sizeof *content);
  if (content == NULL)
  {
    N_VDestroy(v);
    return (NULL);
  }

  /* Attach content */
  v->content = content;

  /* Initialize content */
  content->local_length  = NV_LOCLENGTH_PTC(w);
  content->global_length = NV_GLOBLENGTH_PTC(w);
  content->comm          = NV_COMM_PTC(w);
  content->own_data      = SUNFALSE;
  content->pvec          = NULL;

  return (v);
}

N_Vector nvClone_Petsc(N_Vector w)
{
  N_Vector v = NULL;
  Vec pvec   = NULL;
  Vec wvec   = NV_PVEC_PTC(w);

  /* PetscErrorCode ierr; */

  v = nvCloneEmpty_Petsc(w);
  if (v == NULL) { return (NULL); }

  /* Duplicate vector */
  VecDuplicate(wvec, &pvec);
  if (pvec == NULL)
  {
    nvDestroy_Petsc(v);
    return (NULL);
  }

  /* Attach data */
  NV_OWN_DATA_PTC(v) = SUNTRUE;
  NV_PVEC_PTC(v)     = pvec;

  return (v);
}

void nvDestroy_Petsc(N_Vector v)
{
  if (v == NULL) { return; }

  /* free content */
  if (v->content != NULL)
  {
    if (NV_OWN_DATA_PTC(v) && NV_PVEC_PTC(v) != NULL)
    {
      VecDestroy(&(NV_PVEC_PTC(v)));
      NV_PVEC_PTC(v) = NULL;
    }
    free(v->content);
    v->content = NULL;
  }

  /* free ops and vector */
  if (v->ops != NULL)
  {
    free(v->ops);
    v->ops = NULL;
  }
  free(v);
  v = NULL;

  return;
}

MPI_Comm nvGetCommunicator_Petsc(N_Vector v) { return (NV_COMM_PTC(v)); }

sunindextype nvGetLength_Petsc(N_Vector v) { return (NV_GLOBLENGTH_PTC(v)); }

void nvLinearSum_Petsc(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                       N_Vector z)
{
  Vec xv = NV_PVEC_PTC(x);
  Vec yv = NV_PVEC_PTC(y);
  Vec zv = NV_PVEC_PTC(z);

  if (x == y)
  {
    nvScale_Petsc(a + b, x, z); /* z <~ ax+bx */
    return;
  }

  if (z == y)
  {
    if (b == ONE)
    {
      VecAXPY(yv, a, xv); /* BLAS usage: axpy  y <- ax+y */
      return;
    }
    VecAXPBY(yv, a, b, xv); /* BLAS usage: axpby y <- ax+by */
    return;
  }

  if (z == x)
  {
    if (a == ONE)
    {
      VecAXPY(xv, b, yv); /* BLAS usage: axpy  x <- by+x */
      return;
    }
    VecAXPBY(xv, b, a, yv); /* BLAS usage: axpby x <- by+ax */
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */

  VecAXPBYPCZ(zv, a, b, 0.0, xv, yv); /* PETSc, probably not optimal */

  return;
}

void nvConst_Petsc(sunrealtype c, N_Vector z)
{
  Vec zv = NV_PVEC_PTC(z);

  VecSet(zv, c);

  return;
}

void nvProd_Petsc(N_Vector x, N_Vector y, N_Vector z)
{
  Vec xv = NV_PVEC_PTC(x);
  Vec yv = NV_PVEC_PTC(y);
  Vec zv = NV_PVEC_PTC(z);

  VecPointwiseMult(zv, xv, yv);

  return;
}

void nvDiv_Petsc(N_Vector x, N_Vector y, N_Vector z)
{
  Vec xv = NV_PVEC_PTC(x);
  Vec yv = NV_PVEC_PTC(y);
  Vec zv = NV_PVEC_PTC(z);

  VecPointwiseDivide(zv, xv, yv); /* z = x/y */

  return;
}

void nvScale_Petsc(sunrealtype c, N_Vector x, N_Vector z)
{
  Vec xv = NV_PVEC_PTC(x);
  Vec zv = NV_PVEC_PTC(z);

  if (z == x)
  { /* BLAS usage: scale x <- cx */
    VecScale(xv, c);
    return;
  }

  VecAXPBY(zv, c, 0.0, xv);

  return;
}

void nvAbs_Petsc(N_Vector x, N_Vector z)
{
  Vec xv = NV_PVEC_PTC(x);
  Vec zv = NV_PVEC_PTC(z);

  if (z != x) VecCopy(xv, zv); /* copy x~>z */
  VecAbs(zv);

  return;
}

void nvInv_Petsc(N_Vector x, N_Vector z)
{
  Vec xv = NV_PVEC_PTC(x);
  Vec zv = NV_PVEC_PTC(z);

  if (z != x) VecCopy(xv, zv); /* copy x~>z */
  VecReciprocal(zv);

  return;
}

void nvAddConst_Petsc(N_Vector x, sunrealtype b, N_Vector z)
{
  Vec xv = NV_PVEC_PTC(x);
  Vec zv = NV_PVEC_PTC(z);

  if (z != x) VecCopy(xv, zv); /* copy x~>z */
  VecShift(zv, b);

  return;
}

sunrealtype nvDotProdLocal_Petsc(N_Vector x, N_Vector y)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec xv         = NV_PVEC_PTC(x);
  Vec yv         = NV_PVEC_PTC(y);
  PetscScalar* xd;
  PetscScalar* yd;
  PetscReal sum = ZERO;

  VecGetArray(xv, &xd);
  VecGetArray(yv, &yd);
  for (i = 0; i < N; i++) sum += xd[i] * yd[i];
  VecRestoreArray(xv, &xd);
  VecRestoreArray(yv, &yd);
  return ((sunrealtype)sum);
}

sunrealtype nvDotProd_Petsc(N_Vector x, N_Vector y)
{
  Vec xv = NV_PVEC_PTC(x);
  Vec yv = NV_PVEC_PTC(y);
  PetscScalar dotprod;

  VecDot(xv, yv, &dotprod);
  return dotprod;
}

sunrealtype nvMaxNormLocal_Petsc(N_Vector x)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec xv         = NV_PVEC_PTC(x);
  PetscScalar* xd;
  PetscReal max = ZERO;

  VecGetArray(xv, &xd);
  for (i = 0; i < N; i++)
  {
    if (PetscAbsScalar(xd[i]) > max) max = PetscAbsScalar(xd[i]);
  }
  VecRestoreArray(xv, &xd);
  return ((sunrealtype)max);
}

sunrealtype nvMaxNorm_Petsc(N_Vector x)
{
  Vec xv = NV_PVEC_PTC(x);
  PetscReal norm;

  VecNorm(xv, NORM_INFINITY, &norm);
  return norm;
}

sunrealtype nvWSqrSumLocal_Petsc(N_Vector x, N_Vector w)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec xv         = NV_PVEC_PTC(x);
  Vec wv         = NV_PVEC_PTC(w);
  PetscScalar* xd;
  PetscScalar* wd;
  PetscReal sum = ZERO;

  VecGetArray(xv, &xd);
  VecGetArray(wv, &wd);
  for (i = 0; i < N; i++) { sum += PetscSqr(PetscAbsScalar(xd[i] * wd[i])); }
  VecRestoreArray(xv, &xd);
  VecRestoreArray(wv, &wd);
  return ((sunrealtype)sum);
}

sunrealtype nvWrmsNorm_Petsc(N_Vector x, N_Vector w)
{
  sunrealtype global_sum;
  sunindextype N_global = NV_GLOBLENGTH_PTC(x);
  sunrealtype sum       = nvWSqrSumLocal_Petsc(x, w);
  (void)MPI_Allreduce(&sum, &global_sum, 1, MPI_SUNREALTYPE, MPI_SUM,
                      NV_COMM_PTC(x));
  return (SUNRsqrt(global_sum / N_global));
}

sunrealtype nvWSqrSumMaskLocal_Petsc(N_Vector x, N_Vector w, N_Vector id)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec xv         = NV_PVEC_PTC(x);
  Vec wv         = NV_PVEC_PTC(w);
  Vec idv        = NV_PVEC_PTC(id);
  PetscScalar* xd;
  PetscScalar* wd;
  PetscScalar* idd;
  PetscReal sum = ZERO;

  VecGetArray(xv, &xd);
  VecGetArray(wv, &wd);
  VecGetArray(idv, &idd);
  for (i = 0; i < N; i++)
  {
    PetscReal tag = (PetscReal)idd[i];
    if (tag > ZERO) { sum += PetscSqr(PetscAbsScalar(xd[i] * wd[i])); }
  }
  VecRestoreArray(xv, &xd);
  VecRestoreArray(wv, &wd);
  VecRestoreArray(idv, &idd);
  return sum;
}

sunrealtype nvWrmsNormMask_Petsc(N_Vector x, N_Vector w, N_Vector id)
{
  sunrealtype global_sum;
  sunindextype N_global = NV_GLOBLENGTH_PTC(x);
  sunrealtype sum       = nvWSqrSumMaskLocal_Petsc(x, w, id);
  (void)MPI_Allreduce(&sum, &global_sum, 1, MPI_SUNREALTYPE, MPI_SUM,
                      NV_COMM_PTC(x));
  return (SUNRsqrt(global_sum / N_global));
}

sunrealtype nvMinLocal_Petsc(N_Vector x)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec xv         = NV_PVEC_PTC(x);
  PetscScalar* xd;
  PetscReal min = SUN_BIG_REAL;

  VecGetArray(xv, &xd);
  for (i = 0; i < N; i++)
  {
    if (xd[i] < min) min = xd[i];
  }
  VecRestoreArray(xv, &xd);
  return ((sunrealtype)min);
}

sunrealtype nvMin_Petsc(N_Vector x)
{
  Vec xv = NV_PVEC_PTC(x);
  PetscReal minval;
  PetscInt i;

  VecMin(xv, &i, &minval);
  return minval;
}

sunrealtype nvWL2Norm_Petsc(N_Vector x, N_Vector w)
{
  sunrealtype global_sum;
  sunrealtype sum = nvWSqrSumLocal_Petsc(x, w);
  (void)MPI_Allreduce(&sum, &global_sum, 1, MPI_SUNREALTYPE, MPI_SUM,
                      NV_COMM_PTC(x));
  return (SUNRsqrt(global_sum));
}

sunrealtype nvL1NormLocal_Petsc(N_Vector x)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec xv         = NV_PVEC_PTC(x);
  PetscScalar* xd;
  PetscReal sum = ZERO;

  VecGetArray(xv, &xd);
  for (i = 0; i < N; i++) { sum += PetscAbsScalar(xd[i]); }
  VecRestoreArray(xv, &xd);
  return ((sunrealtype)sum);
}

sunrealtype nvL1Norm_Petsc(N_Vector x)
{
  Vec xv = NV_PVEC_PTC(x);
  PetscReal norm;

  VecNorm(xv, NORM_1, &norm);
  return norm;
}

void nvCompare_Petsc(sunrealtype c, N_Vector x, N_Vector z)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec xv         = NV_PVEC_PTC(x);
  Vec zv         = NV_PVEC_PTC(z);
  PetscReal cpet = c; /* <~ sunrealtype should typedef to PETScReal */
  PetscScalar* xdata;
  PetscScalar* zdata;

  VecGetArray(xv, &xdata);
  VecGetArray(zv, &zdata);
  for (i = 0; i < N; i++)
  {
    zdata[i] = PetscAbsScalar(xdata[i]) >= cpet ? ONE : ZERO;
  }
  VecRestoreArray(xv, &xdata);
  VecRestoreArray(zv, &zdata);

  return;
}

sunbooleantype nvInvTestLocal_Petsc(N_Vector x, N_Vector z)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec xv         = NV_PVEC_PTC(x);
  Vec zv         = NV_PVEC_PTC(z);
  PetscScalar* xd;
  PetscScalar* zd;
  PetscReal val = ONE;

  VecGetArray(xv, &xd);
  VecGetArray(zv, &zd);
  for (i = 0; i < N; i++)
  {
    if (xd[i] == ZERO) val = ZERO;
    else zd[i] = ONE / xd[i];
  }
  VecRestoreArray(xv, &xd);
  VecRestoreArray(zv, &zd);

  if (val == ZERO) { return (SUNFALSE); }
  else { return (SUNTRUE); }
}

sunbooleantype nvInvTest_Petsc(N_Vector x, N_Vector z)
{
  sunrealtype val2;
  sunrealtype val = (nvInvTestLocal_Petsc(x, z)) ? ONE : ZERO;
  (void)MPI_Allreduce(&val, &val2, 1, MPI_SUNREALTYPE, MPI_MIN, NV_COMM_PTC(x));
  if (val2 == ZERO) { return (SUNFALSE); }
  else { return (SUNTRUE); }
}

sunbooleantype nvConstrMaskLocal_Petsc(N_Vector c, N_Vector x, N_Vector m)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  sunrealtype temp;
  sunbooleantype test;
  Vec xv = NV_PVEC_PTC(x);
  Vec cv = NV_PVEC_PTC(c);
  Vec mv = NV_PVEC_PTC(m);
  PetscScalar* xd;
  PetscScalar* cd;
  PetscScalar* md;

  temp = ZERO;

  VecGetArray(xv, &xd);
  VecGetArray(cv, &cd);
  VecGetArray(mv, &md);
  for (i = 0; i < N; i++)
  {
    PetscReal cc = (PetscReal)cd[i]; /* <~ Drop imaginary parts if any. */
    PetscReal xx = (PetscReal)xd[i]; /* <~ Constraints defined on Re{x} */
    md[i]        = ZERO;

    /* Continue if no constraints were set for the variable */
    if (cc == ZERO) { continue; }

    /* Check if a set constraint has been violated */
    test = (SUNRabs(cc) > ONEPT5 && xx * cc <= ZERO) ||
           (SUNRabs(cc) > HALF && xx * cc < ZERO);
    if (test) { temp = md[i] = ONE; }
  }
  VecRestoreArray(xv, &xd);
  VecRestoreArray(cv, &cd);
  VecRestoreArray(mv, &md);

  /* Return false if any constraint was violated */
  return (temp == ONE) ? SUNFALSE : SUNTRUE;
}

sunbooleantype nvConstrMask_Petsc(N_Vector c, N_Vector x, N_Vector m)
{
  sunrealtype temp2;
  sunrealtype temp = (nvConstrMaskLocal_Petsc(c, x, m)) ? ZERO : ONE;
  (void)MPI_Allreduce(&temp, &temp2, 1, MPI_SUNREALTYPE, MPI_MAX, NV_COMM_PTC(x));
  return (temp2 == ONE) ? SUNFALSE : SUNTRUE;
}

sunrealtype nvMinQuotientLocal_Petsc(N_Vector num, N_Vector denom)
{
  sunbooleantype notEvenOnce = SUNTRUE;
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(num);

  Vec nv = NV_PVEC_PTC(num);
  Vec dv = NV_PVEC_PTC(denom);
  PetscScalar* nd;
  PetscScalar* dd;
  PetscReal minval = SUN_BIG_REAL;

  VecGetArray(nv, &nd);
  VecGetArray(dv, &dd);
  for (i = 0; i < N; i++)
  {
    PetscReal nr = (PetscReal)nd[i];
    PetscReal dr = (PetscReal)dd[i];
    if (dr == ZERO) { continue; }
    else
    {
      if (!notEvenOnce) { minval = SUNMIN(minval, nr / dr); }
      else
      {
        minval      = nr / dr;
        notEvenOnce = SUNFALSE;
      }
    }
  }
  VecRestoreArray(nv, &nd);
  VecRestoreArray(dv, &dd);
  return ((sunrealtype)minval);
}

sunrealtype nvMinQuotient_Petsc(N_Vector num, N_Vector denom)
{
  PetscReal gmin;
  sunrealtype minval = nvMinQuotientLocal_Petsc(num, denom);
  (void)MPI_Allreduce(&minval, &gmin, 1, MPI_SUNREALTYPE, MPI_MIN,
                      NV_COMM_PTC(num));
  return (gmin);
}

/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */

SUNErrCode nvLinearCombination_Petsc(int nvec, sunrealtype* c, N_Vector* X,
                                     N_Vector z)
{
  int i;
  Vec* xv;
  Vec zv;

  /* invalid number of vectors */
  if (nvec < 1) { return SUN_ERR_GENERIC; }

  /* should have called N_VScale */
  if (nvec == 1)
  {
    nvScale_Petsc(c[0], X[0], z);
    return SUN_SUCCESS;
  }

  /* should have called N_VLinearSum */
  if (nvec == 2)
  {
    nvLinearSum_Petsc(c[0], X[0], c[1], X[1], z);
    return SUN_SUCCESS;
  }

  /* get petsc vectors */
  xv = (Vec*)malloc(nvec * sizeof(Vec));
  for (i = 0; i < nvec; i++) xv[i] = NV_PVEC_PTC(X[i]);

  zv = NV_PVEC_PTC(z);

  /*
   * X[0] += c[i]*X[i], i = 1,...,nvec-1
   */
  if ((X[0] == z) && (c[0] == ONE))
  {
    VecMAXPY(zv, nvec - 1, c + 1, xv + 1);
    free(xv);
    return SUN_SUCCESS;
  }

  /*
   * X[0] = c[0] * X[0] + sum{ c[i] * X[i] }, i = 1,...,nvec-1
   */
  if (X[0] == z)
  {
    VecScale(zv, c[0]);
    VecMAXPY(zv, nvec - 1, c + 1, xv + 1);
    free(xv);
    return SUN_SUCCESS;
  }

  /*
   * z = sum{ c[i] * X[i] }, i = 0,...,nvec-1
   */
  VecAXPBY(zv, c[0], 0.0, xv[0]);
  VecMAXPY(zv, nvec - 1, c + 1, xv + 1);
  free(xv);

  return SUN_SUCCESS;
}

SUNErrCode nvScaleAddMulti_Petsc(int nvec, sunrealtype* a, N_Vector x,
                                 N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  PetscScalar *xd, *yd, *zd;

  /* invalid number of vectors */
  if (nvec < 1) { return SUN_ERR_GENERIC; }

  /* should have called N_VLinearSum */
  if (nvec == 1)
  {
    nvLinearSum_Petsc(a[0], x, ONE, Y[0], Z[0]);
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N = NV_LOCLENGTH_PTC(x);
  VecGetArray(NV_PVEC_PTC(x), &xd);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z)
  {
    for (i = 0; i < nvec; i++)
    {
      VecGetArray(NV_PVEC_PTC(Y[i]), &yd);
      for (j = 0; j < N; j++) { yd[j] += a[i] * xd[j]; }
      VecRestoreArray(NV_PVEC_PTC(Y[i]), &yd);
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i = 0; i < nvec; i++)
  {
    VecGetArray(NV_PVEC_PTC(Y[i]), &yd);
    VecGetArray(NV_PVEC_PTC(Z[i]), &zd);
    for (j = 0; j < N; j++) { zd[j] = a[i] * xd[j] + yd[j]; }
    VecRestoreArray(NV_PVEC_PTC(Y[i]), &yd);
    VecRestoreArray(NV_PVEC_PTC(Z[i]), &zd);
  }
  return SUN_SUCCESS;
}

SUNErrCode nvDotProdMulti_Petsc(int nvec, N_Vector x, N_Vector* Y,
                                sunrealtype* dotprods)
{
  int i;
  Vec* yv;
  Vec xv;

  /* invalid number of vectors */
  if (nvec < 1) { return SUN_ERR_GENERIC; }

  /* should have called N_VDotProd */
  if (nvec == 1)
  {
    dotprods[0] = nvDotProd_Petsc(x, Y[0]);
    return SUN_SUCCESS;
  }

  /* get petsc vectors */
  yv = (Vec*)malloc(nvec * sizeof(Vec));
  for (i = 0; i < nvec; i++) yv[i] = NV_PVEC_PTC(Y[i]);

  xv = NV_PVEC_PTC(x);

  VecMDot(xv, nvec, yv, dotprods);
  free(yv);

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------------------
 * single buffer reduction operations
 * -----------------------------------------------------------------------------
 */

SUNErrCode nvDotProdMultiLocal_Petsc(int nvec, N_Vector x, N_Vector* Y,
                                     sunrealtype* dotprods)
{
  int j;
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec xv         = NV_PVEC_PTC(x);
  Vec yv;
  PetscScalar* xd;
  PetscScalar* yd;

  VecGetArray(xv, &xd);
  for (j = 0; j < nvec; j++)
  {
    yv = NV_PVEC_PTC(Y[j]);
    VecGetArray(yv, &yd);
    dotprods[j] = ZERO;
    for (i = 0; i < N; i++) { dotprods[j] += xd[i] * yd[i]; }
    VecRestoreArray(yv, &yd);
  }
  VecRestoreArray(xv, &xd);

  return SUN_SUCCESS;
}

SUNErrCode nvDotProdMultiAllReduce_Petsc(int nvec, N_Vector x,
                                         sunrealtype* dotprods)
{
  int retval;
  retval = MPI_Allreduce(MPI_IN_PLACE, dotprods, nvec, MPI_SUNREALTYPE, MPI_SUM,
                         NV_COMM_PTC(x));
  return retval == MPI_SUCCESS ? SUN_SUCCESS : SUN_ERR_GENERIC;
}

/*
 * -----------------------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------------------
 */

SUNErrCode nvLinearSumVectorArray_Petsc(int nvec, sunrealtype a, N_Vector* X,
                                        sunrealtype b, N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  PetscScalar *xd, *yd, *zd;

  /* invalid number of vectors */
  if (nvec < 1) { return SUN_ERR_GENERIC; }

  /* should have called N_VLinearSum */
  if (nvec == 1)
  {
    nvLinearSum_Petsc(a, X[0], b, Y[0], Z[0]);
    return SUN_SUCCESS;
  }

  /* get vector length */
  N = NV_LOCLENGTH_PTC(Z[0]);

  /* compute linear sum for each vector pair in vector arrays */
  for (i = 0; i < nvec; i++)
  {
    VecGetArray(NV_PVEC_PTC(X[i]), &xd);
    VecGetArray(NV_PVEC_PTC(Y[i]), &yd);
    VecGetArray(NV_PVEC_PTC(Z[i]), &zd);
    for (j = 0; j < N; j++) { zd[j] = a * xd[j] + b * yd[j]; }
    VecRestoreArray(NV_PVEC_PTC(X[i]), &xd);
    VecRestoreArray(NV_PVEC_PTC(Y[i]), &yd);
    VecRestoreArray(NV_PVEC_PTC(Z[i]), &zd);
  }

  return SUN_SUCCESS;
}

SUNErrCode nvScaleVectorArray_Petsc(int nvec, sunrealtype* c, N_Vector* X,
                                    N_Vector* Z)
{
  int i;
  sunindextype j, N;
  PetscScalar *xd, *zd;

  /* invalid number of vectors */
  if (nvec < 1) { return SUN_ERR_GENERIC; }

  /* should have called N_VScale */
  if (nvec == 1)
  {
    nvScale_Petsc(c[0], X[0], Z[0]);
    return SUN_SUCCESS;
  }

  /* get vector length */
  N = NV_LOCLENGTH_PTC(Z[0]);

  /*
   * X[i] *= c[i]
   */
  if (X == Z)
  {
    for (i = 0; i < nvec; i++)
    {
      VecGetArray(NV_PVEC_PTC(X[i]), &xd);
      for (j = 0; j < N; j++) { xd[j] *= c[i]; }
      VecRestoreArray(NV_PVEC_PTC(X[i]), &xd);
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[i] = c[i] * X[i]
   */
  for (i = 0; i < nvec; i++)
  {
    VecGetArray(NV_PVEC_PTC(X[i]), &xd);
    VecGetArray(NV_PVEC_PTC(Z[i]), &zd);
    for (j = 0; j < N; j++) { zd[j] = c[i] * xd[j]; }
    VecRestoreArray(NV_PVEC_PTC(X[i]), &xd);
    VecRestoreArray(NV_PVEC_PTC(Z[i]), &zd);
  }
  return SUN_SUCCESS;
}

SUNErrCode nvConstVectorArray_Petsc(int nvec, sunrealtype c, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  PetscScalar* zd;

  /* invalid number of vectors */
  if (nvec < 1) { return SUN_ERR_GENERIC; }

  /* should have called N_VConst */
  if (nvec == 1)
  {
    nvConst_Petsc(c, Z[0]);
    return SUN_SUCCESS;
  }

  /* get vector length */
  N = NV_LOCLENGTH_PTC(Z[0]);

  /* set each vector in the vector array to a constant */
  for (i = 0; i < nvec; i++)
  {
    VecGetArray(NV_PVEC_PTC(Z[i]), &zd);
    for (j = 0; j < N; j++) { zd[j] = c; }
    VecRestoreArray(NV_PVEC_PTC(Z[i]), &zd);
  }

  return SUN_SUCCESS;
}

SUNErrCode nvWrmsNormVectorArray_Petsc(int nvec, N_Vector* X, N_Vector* W,
                                       sunrealtype* nrm)
{
  int i, retval;
  sunindextype j, Nl, Ng;
  sunrealtype* wd = NULL;
  sunrealtype* xd = NULL;
  MPI_Comm comm;

  /* invalid number of vectors */
  if (nvec < 1) { return SUN_ERR_GENERIC; }

  /* should have called N_VWrmsNorm */
  if (nvec == 1)
  {
    nrm[0] = nvWrmsNorm_Petsc(X[0], W[0]);
    return SUN_SUCCESS;
  }

  /* get vector lengths and communicator */
  Nl   = NV_LOCLENGTH_PTC(X[0]);
  Ng   = NV_GLOBLENGTH_PTC(X[0]);
  comm = NV_COMM_PTC(X[0]);

  /* compute the WRMS norm for each vector in the vector array */
  for (i = 0; i < nvec; i++)
  {
    VecGetArray(NV_PVEC_PTC(X[i]), &xd);
    VecGetArray(NV_PVEC_PTC(W[i]), &wd);
    nrm[i] = ZERO;
    for (j = 0; j < Nl; j++)
    {
      nrm[i] += PetscSqr(PetscAbsScalar(xd[j] * wd[j]));
    }
    VecRestoreArray(NV_PVEC_PTC(X[i]), &xd);
    VecRestoreArray(NV_PVEC_PTC(W[i]), &wd);
  }
  retval = MPI_Allreduce(MPI_IN_PLACE, nrm, nvec, MPI_SUNREALTYPE, MPI_SUM, comm);

  for (i = 0; i < nvec; i++) { nrm[i] = SUNRsqrt(nrm[i] / Ng); }

  return retval == MPI_SUCCESS ? SUN_SUCCESS : SUN_ERR_GENERIC;
}

SUNErrCode nvWrmsNormMaskVectorArray_Petsc(int nvec, N_Vector* X, N_Vector* W,
                                           N_Vector id, sunrealtype* nrm)
{
  int i, retval;
  sunindextype j, Nl, Ng;
  PetscScalar *wd, *xd, *idd;
  MPI_Comm comm;

  /* invalid number of vectors */
  if (nvec < 1) { return SUN_ERR_GENERIC; }

  /* should have called N_VWrmsNorm */
  if (nvec == 1)
  {
    nrm[0] = nvWrmsNormMask_Petsc(X[0], W[0], id);
    return SUN_SUCCESS;
  }

  /* get vector lengths and communicator */
  Nl   = NV_LOCLENGTH_PTC(X[0]);
  Ng   = NV_GLOBLENGTH_PTC(X[0]);
  comm = NV_COMM_PTC(X[0]);

  /* compute the WRMS norm for each vector in the vector array */
  VecGetArray(NV_PVEC_PTC(id), &idd);
  for (i = 0; i < nvec; i++)
  {
    VecGetArray(NV_PVEC_PTC(X[i]), &xd);
    VecGetArray(NV_PVEC_PTC(W[i]), &wd);
    nrm[i] = ZERO;
    for (j = 0; j < Nl; j++)
    {
      if (idd[j] > ZERO) nrm[i] += SUNSQR(xd[j] * wd[j]);
    }
    VecRestoreArray(NV_PVEC_PTC(X[i]), &xd);
    VecRestoreArray(NV_PVEC_PTC(W[i]), &wd);
  }
  VecRestoreArray(NV_PVEC_PTC(id), &idd);

  retval = MPI_Allreduce(MPI_IN_PLACE, nrm, nvec, MPI_SUNREALTYPE, MPI_SUM, comm);

  for (i = 0; i < nvec; i++) { nrm[i] = SUNRsqrt(nrm[i] / Ng); }

  return retval == MPI_SUCCESS ? SUN_SUCCESS : SUN_ERR_GENERIC;
}

SUNErrCode nvScaleAddMultiVectorArray_Petsc(int nvec, int nsum, sunrealtype* a,
                                            N_Vector* X, N_Vector** Y,
                                            N_Vector** Z)
{
  int i, j;
  sunindextype k, N;
  PetscScalar *xd, *yd, *zd;

  int retval;
  N_Vector* YY;
  N_Vector* ZZ;

  /* invalid number of vectors */
  if (nvec < 1) { return SUN_ERR_GENERIC; }
  if (nsum < 1) { return SUN_ERR_GENERIC; }

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1)
  {
    /* should have called N_VLinearSum */
    if (nsum == 1)
    {
      nvLinearSum_Petsc(a[0], X[0], ONE, Y[0][0], Z[0][0]);
      return SUN_SUCCESS;
    }

    /* should have called N_VScaleAddMulti */
    YY = (N_Vector*)malloc(nsum * sizeof(N_Vector));
    ZZ = (N_Vector*)malloc(nsum * sizeof(N_Vector));

    for (j = 0; j < nsum; j++)
    {
      YY[j] = Y[j][0];
      ZZ[j] = Z[j][0];
    }

    retval = nvScaleAddMulti_Petsc(nsum, a, X[0], YY, ZZ);

    free(YY);
    free(ZZ);
    return (retval);
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 1)
  {
    retval = nvLinearSumVectorArray_Petsc(nvec, a[0], X, ONE, Y[0], Z[0]);
    return (retval);
  }

  /* ----------------------------
   * Compute multiple linear sums
   * ---------------------------- */

  /* get vector length */
  N = NV_LOCLENGTH_PTC(X[0]);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z)
  {
    for (i = 0; i < nvec; i++)
    {
      VecGetArray(NV_PVEC_PTC(X[i]), &xd);
      for (j = 0; j < nsum; j++)
      {
        VecGetArray(NV_PVEC_PTC(Y[j][i]), &yd);
        for (k = 0; k < N; k++) { yd[k] += a[j] * xd[k]; }
        VecRestoreArray(NV_PVEC_PTC(Y[j][i]), &yd);
      }
      VecRestoreArray(NV_PVEC_PTC(X[i]), &xd);
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i = 0; i < nvec; i++)
  {
    VecGetArray(NV_PVEC_PTC(X[i]), &xd);
    for (j = 0; j < nsum; j++)
    {
      VecGetArray(NV_PVEC_PTC(Y[j][i]), &yd);
      VecGetArray(NV_PVEC_PTC(Z[j][i]), &zd);
      for (k = 0; k < N; k++) { zd[k] = a[j] * xd[k] + yd[k]; }
      VecRestoreArray(NV_PVEC_PTC(Y[j][i]), &yd);
      VecRestoreArray(NV_PVEC_PTC(Z[j][i]), &zd);
    }
    VecRestoreArray(NV_PVEC_PTC(X[i]), &xd);
  }

  return SUN_SUCCESS;
}

SUNErrCode nvLinearCombinationVectorArray_Petsc(int nvec, int nsum,
                                                sunrealtype* c, N_Vector** X,
                                                N_Vector* Z)
{
  int i;          /* vector arrays index in summation [0,nsum) */
  int j;          /* vector index in vector array     [0,nvec) */
  sunindextype k; /* element index in vector          [0,N)    */
  sunindextype N;
  PetscScalar *zd, *xd;

  sunrealtype* ctmp;
  N_Vector* Y;

  /* invalid number of vectors */
  if (nvec < 1) { return SUN_ERR_GENERIC; }
  if (nsum < 1) { return SUN_ERR_GENERIC; }

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1)
  {
    /* should have called N_VScale */
    if (nsum == 1)
    {
      nvScale_Petsc(c[0], X[0][0], Z[0]);
      return SUN_SUCCESS;
    }

    /* should have called N_VLinearSum */
    if (nsum == 2)
    {
      nvLinearSum_Petsc(c[0], X[0][0], c[1], X[1][0], Z[0]);
      return SUN_SUCCESS;
    }

    /* should have called N_VLinearCombination */
    Y = (N_Vector*)malloc(nsum * sizeof(N_Vector));

    for (i = 0; i < nsum; i++) { Y[i] = X[i][0]; }

    nvLinearCombination_Petsc(nsum, c, Y, Z[0]);

    free(Y);
    return SUN_SUCCESS;
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VScaleVectorArray */
  if (nsum == 1)
  {
    ctmp = (sunrealtype*)malloc(nvec * sizeof(sunrealtype));

    for (j = 0; j < nvec; j++) { ctmp[j] = c[0]; }

    nvScaleVectorArray_Petsc(nvec, ctmp, X[0], Z);

    free(ctmp);
    return SUN_SUCCESS;
  }

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 2)
  {
    nvLinearSumVectorArray_Petsc(nvec, c[0], X[0], c[1], X[1], Z);
    return SUN_SUCCESS;
  }

  /* --------------------------
   * Compute linear combination
   * -------------------------- */

  /* get vector length */
  N = NV_LOCLENGTH_PTC(Z[0]);

  /*
   * X[0][j] += c[i]*X[i][j], i = 1,...,nvec-1
   */
  if ((X[0] == Z) && (c[0] == ONE))
  {
    for (j = 0; j < nvec; j++)
    {
      VecGetArray(NV_PVEC_PTC(Z[j]), &zd);
      for (i = 1; i < nsum; i++)
      {
        VecGetArray(NV_PVEC_PTC(X[i][j]), &xd);
        for (k = 0; k < N; k++) { zd[k] += c[i] * xd[k]; }
        VecRestoreArray(NV_PVEC_PTC(X[i][j]), &xd);
      }
      VecRestoreArray(NV_PVEC_PTC(Z[j]), &zd);
    }
    return SUN_SUCCESS;
  }

  /*
   * X[0][j] = c[0] * X[0][j] + sum{ c[i] * X[i][j] }, i = 1,...,nvec-1
   */
  if (X[0] == Z)
  {
    for (j = 0; j < nvec; j++)
    {
      VecGetArray(NV_PVEC_PTC(Z[j]), &zd);
      for (k = 0; k < N; k++) { zd[k] *= c[0]; }
      for (i = 1; i < nsum; i++)
      {
        VecGetArray(NV_PVEC_PTC(X[i][j]), &xd);
        for (k = 0; k < N; k++) { zd[k] += c[i] * xd[k]; }
        VecRestoreArray(NV_PVEC_PTC(X[i][j]), &xd);
      }
      VecRestoreArray(NV_PVEC_PTC(Z[j]), &zd);
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[j] = sum{ c[i] * X[i][j] }, i = 0,...,nvec-1
   */
  for (j = 0; j < nvec; j++)
  {
    VecGetArray(NV_PVEC_PTC(X[0][j]), &xd);
    VecGetArray(NV_PVEC_PTC(Z[j]), &zd);
    for (k = 0; k < N; k++) { zd[k] = c[0] * xd[k]; }
    VecRestoreArray(NV_PVEC_PTC(X[0][j]), &xd);
    for (i = 1; i < nsum; i++)
    {
      VecGetArray(NV_PVEC_PTC(X[i][j]), &xd);
      for (k = 0; k < N; k++) { zd[k] += c[i] * xd[k]; }
      VecRestoreArray(NV_PVEC_PTC(X[i][j]), &xd);
    }
    VecRestoreArray(NV_PVEC_PTC(Z[j]), &zd);
  }
  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * OPTIONAL XBraid interface operations
 * -----------------------------------------------------------------
 */

SUNErrCode nvBufSize_Petsc(N_Vector x, sunindextype* size)
{
  if (x == NULL) { return SUN_ERR_GENERIC; }
  *size = NV_LOCLENGTH_PTC(x) * ((sunindextype)sizeof(PetscScalar));
  return SUN_SUCCESS;
}

SUNErrCode nvBufPack_Petsc(N_Vector x, void* buf)
{
  Vec xv;
  sunindextype i, N;
  PetscScalar* xd = NULL;
  PetscScalar* bd = NULL;

  if (x == NULL || buf == NULL) { return SUN_ERR_GENERIC; }

  xv = NV_PVEC_PTC(x);
  N  = NV_LOCLENGTH_PTC(x);
  bd = (PetscScalar*)buf;

  VecGetArray(xv, &xd);
  for (i = 0; i < N; i++) bd[i] = xd[i];
  VecRestoreArray(xv, &xd);

  return SUN_SUCCESS;
}

SUNErrCode nvBufUnpack_Petsc(N_Vector x, void* buf)
{
  Vec xv;
  sunindextype i, N;
  PetscScalar* xd = NULL;
  PetscScalar* bd = NULL;

  if (x == NULL || buf == NULL) { return SUN_ERR_GENERIC; }

  xv = NV_PVEC_PTC(x);
  N  = NV_LOCLENGTH_PTC(x);
  bd = (PetscScalar*)buf;

  VecGetArray(xv, &xd);
  for (i = 0; i < N; i++) xd[i] = bd[i];
  VecRestoreArray(xv, &xd);

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

SUNErrCode N_VEnableFusedOps_Petsc(N_Vector v, sunbooleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) { return SUN_ERR_GENERIC; }

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }

  if (tf)
  {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = nvLinearCombination_Petsc;
    v->ops->nvscaleaddmulti     = nvScaleAddMulti_Petsc;
    v->ops->nvdotprodmulti      = nvDotProdMulti_Petsc;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray     = nvLinearSumVectorArray_Petsc;
    v->ops->nvscalevectorarray         = nvScaleVectorArray_Petsc;
    v->ops->nvconstvectorarray         = nvConstVectorArray_Petsc;
    v->ops->nvwrmsnormvectorarray      = nvWrmsNormVectorArray_Petsc;
    v->ops->nvwrmsnormmaskvectorarray  = nvWrmsNormMaskVectorArray_Petsc;
    v->ops->nvscaleaddmultivectorarray = nvScaleAddMultiVectorArray_Petsc;
    v->ops->nvlinearcombinationvectorarray = nvLinearCombinationVectorArray_Petsc;
    /* enable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = nvDotProdMultiLocal_Petsc;
  }
  else
  {
    /* disable all fused vector operations */
    v->ops->nvlinearcombination = NULL;
    v->ops->nvscaleaddmulti     = NULL;
    v->ops->nvdotprodmulti      = NULL;
    /* disable all vector array operations */
    v->ops->nvlinearsumvectorarray         = NULL;
    v->ops->nvscalevectorarray             = NULL;
    v->ops->nvconstvectorarray             = NULL;
    v->ops->nvwrmsnormvectorarray          = NULL;
    v->ops->nvwrmsnormmaskvectorarray      = NULL;
    v->ops->nvscaleaddmultivectorarray     = NULL;
    v->ops->nvlinearcombinationvectorarray = NULL;
    /* disable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = NULL;
  }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearCombination_Petsc(N_Vector v, sunbooleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) { return SUN_ERR_GENERIC; }

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }

  /* enable/disable operation */
  if (tf) { v->ops->nvlinearcombination = nvLinearCombination_Petsc; }
  else { v->ops->nvlinearcombination = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleAddMulti_Petsc(N_Vector v, sunbooleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) { return SUN_ERR_GENERIC; }

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }

  /* enable/disable operation */
  if (tf) { v->ops->nvscaleaddmulti = nvScaleAddMulti_Petsc; }
  else { v->ops->nvscaleaddmulti = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableDotProdMulti_Petsc(N_Vector v, sunbooleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) { return SUN_ERR_GENERIC; }

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }

  /* enable/disable operation */
  if (tf) { v->ops->nvdotprodmulti = nvDotProdMulti_Petsc; }
  else { v->ops->nvdotprodmulti = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearSumVectorArray_Petsc(N_Vector v, sunbooleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) { return SUN_ERR_GENERIC; }

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }

  /* enable/disable operation */
  if (tf) { v->ops->nvlinearsumvectorarray = nvLinearSumVectorArray_Petsc; }
  else { v->ops->nvlinearsumvectorarray = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleVectorArray_Petsc(N_Vector v, sunbooleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) { return SUN_ERR_GENERIC; }

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }

  /* enable/disable operation */
  if (tf) { v->ops->nvscalevectorarray = nvScaleVectorArray_Petsc; }
  else { v->ops->nvscalevectorarray = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableConstVectorArray_Petsc(N_Vector v, sunbooleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) { return SUN_ERR_GENERIC; }

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }

  /* enable/disable operation */
  if (tf) { v->ops->nvconstvectorarray = nvConstVectorArray_Petsc; }
  else { v->ops->nvconstvectorarray = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableWrmsNormVectorArray_Petsc(N_Vector v, sunbooleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) { return SUN_ERR_GENERIC; }

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }

  /* enable/disable operation */
  if (tf) { v->ops->nvwrmsnormvectorarray = nvWrmsNormVectorArray_Petsc; }
  else { v->ops->nvwrmsnormvectorarray = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableWrmsNormMaskVectorArray_Petsc(N_Vector v, sunbooleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) { return SUN_ERR_GENERIC; }

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }

  /* enable/disable operation */
  if (tf)
  {
    v->ops->nvwrmsnormmaskvectorarray = nvWrmsNormMaskVectorArray_Petsc;
  }
  else { v->ops->nvwrmsnormmaskvectorarray = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableScaleAddMultiVectorArray_Petsc(N_Vector v, sunbooleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) { return SUN_ERR_GENERIC; }

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }

  /* enable/disable operation */
  if (tf)
  {
    v->ops->nvscaleaddmultivectorarray = nvScaleAddMultiVectorArray_Petsc;
  }
  else { v->ops->nvscaleaddmultivectorarray = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableLinearCombinationVectorArray_Petsc(N_Vector v,
                                                       sunbooleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) { return SUN_ERR_GENERIC; }

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }

  /* enable/disable operation */
  if (tf)
  {
    v->ops->nvlinearcombinationvectorarray = nvLinearCombinationVectorArray_Petsc;
  }
  else { v->ops->nvlinearcombinationvectorarray = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode N_VEnableDotProdMultiLocal_Petsc(N_Vector v, sunbooleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) { return SUN_ERR_GENERIC; }

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) { return SUN_ERR_GENERIC; }

  /* enable/disable operation */
  if (tf) { v->ops->nvdotprodmultilocal = nvDotProdMultiLocal_Petsc; }
  else { v->ops->nvdotprodmultilocal = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

/* Deprecated concrete operation wrappers */

void N_VAbs_Petsc(N_Vector x, N_Vector z) { nvAbs_Petsc(x, z); }

void N_VAddConst_Petsc(N_Vector x, sunrealtype b, N_Vector z)
{
  nvAddConst_Petsc(x, b, z);
}

SUNErrCode N_VBufPack_Petsc(N_Vector x, void* buf)
{
  return nvBufPack_Petsc(x, buf);
}

SUNErrCode N_VBufSize_Petsc(N_Vector x, sunindextype* size)
{
  return nvBufSize_Petsc(x, size);
}

SUNErrCode N_VBufUnpack_Petsc(N_Vector x, void* buf)
{
  return nvBufUnpack_Petsc(x, buf);
}

N_Vector N_VCloneEmpty_Petsc(N_Vector w) { return nvCloneEmpty_Petsc(w); }

N_Vector N_VClone_Petsc(N_Vector w) { return nvClone_Petsc(w); }

void N_VCompare_Petsc(sunrealtype c, N_Vector x, N_Vector z)
{
  nvCompare_Petsc(c, x, z);
}

SUNErrCode N_VConstVectorArray_Petsc(int nvecs, sunrealtype c, N_Vector* Z)
{
  return nvConstVectorArray_Petsc(nvecs, c, Z);
}

void N_VConst_Petsc(sunrealtype c, N_Vector z) { nvConst_Petsc(c, z); }

sunbooleantype N_VConstrMaskLocal_Petsc(N_Vector c, N_Vector x, N_Vector m)
{
  return nvConstrMaskLocal_Petsc(c, x, m);
}

sunbooleantype N_VConstrMask_Petsc(N_Vector c, N_Vector x, N_Vector m)
{
  return nvConstrMask_Petsc(c, x, m);
}

void N_VDestroy_Petsc(N_Vector v) { nvDestroy_Petsc(v); }

void N_VDiv_Petsc(N_Vector x, N_Vector y, N_Vector z) { nvDiv_Petsc(x, y, z); }

sunrealtype N_VDotProdLocal_Petsc(N_Vector x, N_Vector y)
{
  return nvDotProdLocal_Petsc(x, y);
}

SUNErrCode N_VDotProdMultiAllReduce_Petsc(int nvec, N_Vector x, sunrealtype* sum)
{
  return nvDotProdMultiAllReduce_Petsc(nvec, x, sum);
}

SUNErrCode N_VDotProdMultiLocal_Petsc(int nvec, N_Vector x, N_Vector* Y,
                                      sunrealtype* dotprods)
{
  return nvDotProdMultiLocal_Petsc(nvec, x, Y, dotprods);
}

SUNErrCode N_VDotProdMulti_Petsc(int nvec, N_Vector x, N_Vector* Y,
                                 sunrealtype* dotprods)
{
  return nvDotProdMulti_Petsc(nvec, x, Y, dotprods);
}

sunrealtype N_VDotProd_Petsc(N_Vector x, N_Vector y)
{
  return nvDotProd_Petsc(x, y);
}

MPI_Comm N_VGetCommunicator_Petsc(N_Vector v)
{
  return nvGetCommunicator_Petsc(v);
}

sunindextype N_VGetLength_Petsc(N_Vector v) { return nvGetLength_Petsc(v); }

N_Vector_ID N_VGetVectorID_Petsc(N_Vector v) { return nvGetVectorID_Petsc(v); }

sunbooleantype N_VInvTestLocal_Petsc(N_Vector x, N_Vector z)
{
  return nvInvTestLocal_Petsc(x, z);
}

sunbooleantype N_VInvTest_Petsc(N_Vector x, N_Vector z)
{
  return nvInvTest_Petsc(x, z);
}

void N_VInv_Petsc(N_Vector x, N_Vector z) { nvInv_Petsc(x, z); }

sunrealtype N_VL1NormLocal_Petsc(N_Vector x) { return nvL1NormLocal_Petsc(x); }

sunrealtype N_VL1Norm_Petsc(N_Vector x) { return nvL1Norm_Petsc(x); }

SUNErrCode N_VLinearCombinationVectorArray_Petsc(int nvec, int nsum,
                                                 sunrealtype* c, N_Vector** X,
                                                 N_Vector* Z)
{
  return nvLinearCombinationVectorArray_Petsc(nvec, nsum, c, X, Z);
}

SUNErrCode N_VLinearCombination_Petsc(int nvec, sunrealtype* c, N_Vector* X,
                                      N_Vector z)
{
  return nvLinearCombination_Petsc(nvec, c, X, z);
}

SUNErrCode N_VLinearSumVectorArray_Petsc(int nvec, sunrealtype a, N_Vector* X,
                                         sunrealtype b, N_Vector* Y, N_Vector* Z)
{
  return nvLinearSumVectorArray_Petsc(nvec, a, X, b, Y, Z);
}

void N_VLinearSum_Petsc(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                        N_Vector z)
{
  nvLinearSum_Petsc(a, x, b, y, z);
}

sunrealtype N_VMaxNormLocal_Petsc(N_Vector x)
{
  return nvMaxNormLocal_Petsc(x);
}

sunrealtype N_VMaxNorm_Petsc(N_Vector x) { return nvMaxNorm_Petsc(x); }

sunrealtype N_VMinLocal_Petsc(N_Vector x) { return nvMinLocal_Petsc(x); }

sunrealtype N_VMinQuotientLocal_Petsc(N_Vector num, N_Vector denom)
{
  return nvMinQuotientLocal_Petsc(num, denom);
}

sunrealtype N_VMinQuotient_Petsc(N_Vector num, N_Vector denom)
{
  return nvMinQuotient_Petsc(num, denom);
}

sunrealtype N_VMin_Petsc(N_Vector x) { return nvMin_Petsc(x); }

void N_VProd_Petsc(N_Vector x, N_Vector y, N_Vector z)
{
  nvProd_Petsc(x, y, z);
}

SUNErrCode N_VScaleAddMultiVectorArray_Petsc(int nvec, int nsum, sunrealtype* a,
                                             N_Vector* X, N_Vector** Y,
                                             N_Vector** Z)
{
  return nvScaleAddMultiVectorArray_Petsc(nvec, nsum, a, X, Y, Z);
}

SUNErrCode N_VScaleAddMulti_Petsc(int nvec, sunrealtype* a, N_Vector x,
                                  N_Vector* Y, N_Vector* Z)
{
  return nvScaleAddMulti_Petsc(nvec, a, x, Y, Z);
}

SUNErrCode N_VScaleVectorArray_Petsc(int nvec, sunrealtype* c, N_Vector* X,
                                     N_Vector* Z)
{
  return nvScaleVectorArray_Petsc(nvec, c, X, Z);
}

void N_VScale_Petsc(sunrealtype c, N_Vector x, N_Vector z)
{
  nvScale_Petsc(c, x, z);
}

sunrealtype N_VWL2Norm_Petsc(N_Vector x, N_Vector w)
{
  return nvWL2Norm_Petsc(x, w);
}

sunrealtype N_VWSqrSumLocal_Petsc(N_Vector x, N_Vector w)
{
  return nvWSqrSumLocal_Petsc(x, w);
}

sunrealtype N_VWSqrSumMaskLocal_Petsc(N_Vector x, N_Vector w, N_Vector id)
{
  return nvWSqrSumMaskLocal_Petsc(x, w, id);
}

SUNErrCode N_VWrmsNormMaskVectorArray_Petsc(int nvec, N_Vector* X, N_Vector* W,
                                            N_Vector id, sunrealtype* nrm)
{
  return nvWrmsNormMaskVectorArray_Petsc(nvec, X, W, id, nrm);
}

sunrealtype N_VWrmsNormMask_Petsc(N_Vector x, N_Vector w, N_Vector id)
{
  return nvWrmsNormMask_Petsc(x, w, id);
}

SUNErrCode N_VWrmsNormVectorArray_Petsc(int nvecs, N_Vector* X, N_Vector* W,
                                        sunrealtype* nrm)
{
  return nvWrmsNormVectorArray_Petsc(nvecs, X, W, nrm);
}

sunrealtype N_VWrmsNorm_Petsc(N_Vector x, N_Vector w)
{
  return nvWrmsNorm_Petsc(x, w);
}

void N_VPrint_Petsc(N_Vector v) { nvPrint_Petsc(v); }

void N_VPrintFile_Petsc(N_Vector v, const char fname[])
{
  nvPrintFile_Petsc(v, fname);
}
