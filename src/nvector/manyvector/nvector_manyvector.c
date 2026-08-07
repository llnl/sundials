/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
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
 * This is the implementation file for the ManyVector implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_types.h>
#ifdef MANYVECTOR_BUILD_WITH_MPI
#include <nvector/nvector_mpimanyvector.h>
#include <sundials/priv/sundials_mpi_errors_impl.h>
#else
#include <nvector/nvector_manyvector.h>
#endif
#include <sundials/priv/sundials_context_impl.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>

#include "sundials_macros.h"

/* Macro to handle separate MPI-aware/unaware installations */
#ifdef MANYVECTOR_BUILD_WITH_MPI
#define MVAPPEND(fun) fun##_MPIManyVector
#else
#define MVAPPEND(fun) fun##_ManyVector
#endif

#define ZERO SUN_RCONST(0.0)
#ifdef MANYVECTOR_BUILD_WITH_MPI
#define ONE SUN_RCONST(1.0)
#endif

/* -----------------------------------------------------------------
   ManyVector content accessor macros
   -----------------------------------------------------------------*/
#ifdef MANYVECTOR_BUILD_WITH_MPI
#define MANYVECTOR_CONTENT(v) ((N_VectorContent_MPIManyVector)(v->content))
#define MANYVECTOR_COMM(v)    (MANYVECTOR_CONTENT(v)->comm)
#else
#define MANYVECTOR_CONTENT(v) ((N_VectorContent_ManyVector)(v->content))
#endif
#define MANYVECTOR_NUM_SUBVECS(v) (MANYVECTOR_CONTENT(v)->num_subvectors)
#define MANYVECTOR_GLOBLENGTH(v)  (MANYVECTOR_CONTENT(v)->global_length)
#define MANYVECTOR_SUBVECS(v)     (MANYVECTOR_CONTENT(v)->subvec_array)
#define MANYVECTOR_SUBVEC(v, i)   (MANYVECTOR_SUBVECS(v)[i])
#define MANYVECTOR_OWN_DATA(v)    (MANYVECTOR_CONTENT(v)->own_data)

/* -----------------------------------------------------------------
   Prototypes of utility routines
   -----------------------------------------------------------------*/
static N_Vector ManyVectorClone(N_Vector w, sunbooleantype cloneempty);
#ifdef MANYVECTOR_BUILD_WITH_MPI
static int SubvectorMPIRank(N_Vector w);
#endif

/* -----------------------------------------------------------------
   ManyVector API routines
   -----------------------------------------------------------------*/

#ifdef MANYVECTOR_BUILD_WITH_MPI

/* This function creates an MPIManyVector from a set of existing
   N_Vector objects, along with a user-created MPI (inter/intra)communicator
   that couples all subvectors together. */
N_Vector N_VMake_MPIManyVector(MPI_Comm comm, sunindextype num_subvectors,
                               N_Vector* vec_array, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  N_Vector v;
  N_VectorContent_MPIManyVector content;
  sunindextype i, local_length;
  int rank;

  /* Check that input N_Vectors are non-NULL */
  SUNAssertNull(vec_array, SUN_ERR_ARG_CORRUPT);
  for (i = 0; i < num_subvectors; i++)
  {
    SUNAssertNull(vec_array[i], SUN_ERR_ARG_CORRUPT);
  }

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = N_VGetVectorID_MPIManyVector;
  v->ops->nvcloneempty      = N_VCloneEmpty_MPIManyVector;
  v->ops->nvclone           = N_VClone_MPIManyVector;
  v->ops->nvdestroy         = N_VDestroy_MPIManyVector;
  v->ops->nvspace           = N_VSpace_MPIManyVector;
  v->ops->nvgetcommunicator = N_VGetCommunicator_MPIManyVector;
  v->ops->nvgetlength       = N_VGetLength_MPIManyVector;

  /* standard vector operations */
  v->ops->nvlinearsum    = N_VLinearSum_MPIManyVector;
  v->ops->nvconst        = N_VConst_MPIManyVector;
  v->ops->nvprod         = N_VProd_MPIManyVector;
  v->ops->nvdiv          = N_VDiv_MPIManyVector;
  v->ops->nvscale        = N_VScale_MPIManyVector;
  v->ops->nvabs          = N_VAbs_MPIManyVector;
  v->ops->nvinv          = N_VInv_MPIManyVector;
  v->ops->nvaddconst     = N_VAddConst_MPIManyVector;
  v->ops->nvdotprod      = N_VDotProd_MPIManyVector;
  v->ops->nvmaxnorm      = N_VMaxNorm_MPIManyVector;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_MPIManyVector;
  v->ops->nvwrmsnormmask = N_VWrmsNormMask_MPIManyVector;
  v->ops->nvmin          = N_VMin_MPIManyVector;
  v->ops->nvwl2norm      = N_VWL2Norm_MPIManyVector;
  v->ops->nvl1norm       = N_VL1Norm_MPIManyVector;
  v->ops->nvcompare      = N_VCompare_MPIManyVector;
  v->ops->nvinvtest      = N_VInvTest_MPIManyVector;
  v->ops->nvconstrmask   = N_VConstrMask_MPIManyVector;
  v->ops->nvminquotient  = N_VMinQuotient_MPIManyVector;

  /* fused vector operations */
  v->ops->nvlinearcombination = N_VLinearCombination_MPIManyVector;
  v->ops->nvscaleaddmulti     = N_VScaleAddMulti_MPIManyVector;
  v->ops->nvdotprodmulti      = N_VDotProdMulti_MPIManyVector;

  /* vector array operations */
  v->ops->nvwrmsnormvectorarray     = N_VWrmsNormVectorArray_MPIManyVector;
  v->ops->nvwrmsnormmaskvectorarray = N_VWrmsNormMaskVectorArray_MPIManyVector;

  /* local reduction operations */
  v->ops->nvdotprodlocal     = N_VDotProdLocal_MPIManyVector;
  v->ops->nvmaxnormlocal     = N_VMaxNormLocal_MPIManyVector;
  v->ops->nvminlocal         = N_VMinLocal_MPIManyVector;
  v->ops->nvl1normlocal      = N_VL1NormLocal_MPIManyVector;
  v->ops->nvinvtestlocal     = N_VInvTestLocal_MPIManyVector;
  v->ops->nvconstrmasklocal  = N_VConstrMaskLocal_MPIManyVector;
  v->ops->nvminquotientlocal = N_VMinQuotientLocal_MPIManyVector;
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_MPIManyVector;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_MPIManyVector;

  /* single buffer reduction operations */
  v->ops->nvdotprodmultilocal     = N_VDotProdMultiLocal_MPIManyVector;
  v->ops->nvdotprodmultiallreduce = N_VDotProdMultiAllReduce_MPIManyVector;

  /* XBraid interface operations */
  v->ops->nvbufsize   = N_VBufSize_MPIManyVector;
  v->ops->nvbufpack   = N_VBufPack_MPIManyVector;
  v->ops->nvbufunpack = N_VBufUnpack_MPIManyVector;

  /* debugging functions */
  v->ops->nvprint     = N_VPrint_MPIManyVector;
  v->ops->nvprintfile = N_VPrintFile_MPIManyVector;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_MPIManyVector)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  v->content = content;

  /* Attach content components */

  /* set scalar content entries, and allocate/set subvector array */
  content->comm           = MPI_COMM_NULL;
  content->num_subvectors = num_subvectors;
  content->own_data       = SUNFALSE;
  content->subvec_array   = NULL;
  content->subvec_array = (N_Vector*)malloc(num_subvectors * sizeof(N_Vector));
  SUNAssertNull(content->subvec_array, SUN_ERR_MALLOC_FAIL);

  for (i = 0; i < num_subvectors; i++)
  {
    content->subvec_array[i] = vec_array[i];
  }

  /* duplicate input communicator (if non-NULL) */
  if (comm != MPI_COMM_NULL)
  {
    SUNCheckMPICallNull(MPI_Comm_dup(comm, &(content->comm)));
  }

  /* Determine overall MPIManyVector length: sum contributions from all
     subvectors where this rank is the root, then perform reduction */
  local_length = 0;
  for (i = 0; i < num_subvectors; i++)
  {
    /* determine rank of this task in subvector communicator
       (serial vectors default to rank=0) */
    rank = SubvectorMPIRank(vec_array[i]);
    if (rank < 0)
    {
      N_VDestroy(v);
      SUNCheckLastErrNull();
      return (NULL);
    }

    /* accumulate contribution from root tasks */
    if (vec_array[i]->ops->nvgetlength)
    {
      if (rank == 0)
      {
        local_length += N_VGetLength(vec_array[i]);
        SUNCheckLastErrNull();
      }
    }
    else
    {
      N_VDestroy(v);
      SUNCheckLastErrNull();
      return (NULL);
    }
  }
  if (content->comm != MPI_COMM_NULL)
  {
    SUNCheckMPICallNull(MPI_Allreduce(&local_length, &(content->global_length), 1,
                                      MPI_SUNINDEXTYPE, MPI_SUM, content->comm));
  }
  else { content->global_length = local_length; }

  return (v);
}
#endif

#ifdef MANYVECTOR_BUILD_WITH_MPI
/* This function creates an MPIManyVector from a set of existing
   N_Vector objects, under the requirement that all MPI-aware
   sub-vectors use the same MPI communicator (this is verified
   internally).  If no sub-vector is MPI-aware, then this
   function will return NULL. For single-node partitioning,
   the regular (not MPI aware) manyvector should be used. */
N_Vector N_VNew_MPIManyVector(sunindextype num_subvectors, N_Vector* vec_array,
                              SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  sunbooleantype nocommfound;
  MPI_Comm comm, vcomm;
  int comparison;
  N_Vector v = NULL;

  /* Check that all subvectors have identical MPI communicators (if present) */
  nocommfound = SUNTRUE;
  comm        = MPI_COMM_NULL;
  for (sunindextype i = 0; i < num_subvectors; i++)
  {
    /* access MPI communicator for subvector i (vcomm);
       if none is present then continue to next subvector */
    vcomm = N_VGetCommunicator(vec_array[i]);
    if (vcomm == MPI_COMM_NULL) { continue; }

    /* if this is the first communicator, create a copy */
    if (nocommfound)
    {
      /* set comm to duplicate this first subvector communicator */
      SUNCheckMPICallNull(MPI_Comm_dup(vcomm, &comm));
      nocommfound = SUNFALSE;

      /* otherwise, verify that vcomm matches stored comm */
    }
    else
    {
      SUNCheckMPICallNoRet(MPI_Comm_compare(vcomm, comm, &comparison));
      if ((comparison != MPI_IDENT) && (comparison != MPI_CONGRUENT))
      {
        SUNHandleErrWithMsg(__LINE__, __func__, __FILE__,
                            "All subvectors must have the same communicator, "
                            "i.e., MPI_Comm_compare must return MPI_IDENT or "
                            "MPI_CONGRUENT.",
                            SUN_ERR_ARG_INCOMPATIBLE, SUNCTX_);
        return NULL;
      }
    }
  }

  if (!nocommfound)
  {
    /* Create vector using "Make" routine and shared communicator (if non-NULL) */
    v = N_VMake_MPIManyVector(comm, num_subvectors, vec_array, sunctx);
    SUNCheckLastErrNull();
    if (comm != MPI_COMM_NULL) { SUNCheckMPICallNull(MPI_Comm_free(&comm)); }
  }

  return (v);
}
#else
/* This function creates a ManyVector from a set of existing
   N_Vector objects.  ManyVector objects created with this constructor
   may be used to describe data partitioning within a single node. */
N_Vector N_VNew_ManyVector(sunindextype num_subvectors, N_Vector* vec_array,
                           SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  N_Vector v;
  N_VectorContent_ManyVector content;
  sunindextype i, local_length;

  /* Check that input N_Vectors are non-NULL */
  SUNAssertNull(vec_array, SUN_ERR_ARG_CORRUPT);
  for (i = 0; i < num_subvectors; i++)
  {
    SUNAssertNull(vec_array[i], SUN_ERR_ARG_CORRUPT);
  }

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid = N_VGetVectorID_ManyVector;
  v->ops->nvcloneempty  = N_VCloneEmpty_ManyVector;
  v->ops->nvclone       = N_VClone_ManyVector;
  v->ops->nvdestroy     = N_VDestroy_ManyVector;
  v->ops->nvspace       = N_VSpace_ManyVector;
  v->ops->nvgetlength   = N_VGetLength_ManyVector;

  /* standard vector operations */
  v->ops->nvlinearsum    = N_VLinearSum_ManyVector;
  v->ops->nvconst        = N_VConst_ManyVector;
  v->ops->nvprod         = N_VProd_ManyVector;
  v->ops->nvdiv          = N_VDiv_ManyVector;
  v->ops->nvscale        = N_VScale_ManyVector;
  v->ops->nvabs          = N_VAbs_ManyVector;
  v->ops->nvinv          = N_VInv_ManyVector;
  v->ops->nvaddconst     = N_VAddConst_ManyVector;
  v->ops->nvdotprod      = N_VDotProdLocal_ManyVector;
  v->ops->nvmaxnorm      = N_VMaxNormLocal_ManyVector;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_ManyVector;
  v->ops->nvwrmsnormmask = N_VWrmsNormMask_ManyVector;
  v->ops->nvmin          = N_VMinLocal_ManyVector;
  v->ops->nvwl2norm      = N_VWL2Norm_ManyVector;
  v->ops->nvl1norm       = N_VL1NormLocal_ManyVector;
  v->ops->nvcompare      = N_VCompare_ManyVector;
  v->ops->nvinvtest      = N_VInvTestLocal_ManyVector;
  v->ops->nvconstrmask   = N_VConstrMaskLocal_ManyVector;
  v->ops->nvminquotient  = N_VMinQuotientLocal_ManyVector;

  /* fused vector operations */
  v->ops->nvlinearcombination = N_VLinearCombination_ManyVector;
  v->ops->nvscaleaddmulti     = N_VScaleAddMulti_ManyVector;
  v->ops->nvdotprodmulti      = N_VDotProdMulti_ManyVector;

  /* vector array operations */
  v->ops->nvwrmsnormvectorarray     = N_VWrmsNormVectorArray_ManyVector;
  v->ops->nvwrmsnormmaskvectorarray = N_VWrmsNormMaskVectorArray_ManyVector;

  /* local reduction operations */
  v->ops->nvdotprodlocal     = N_VDotProdLocal_ManyVector;
  v->ops->nvmaxnormlocal     = N_VMaxNormLocal_ManyVector;
  v->ops->nvminlocal         = N_VMinLocal_ManyVector;
  v->ops->nvl1normlocal      = N_VL1NormLocal_ManyVector;
  v->ops->nvinvtestlocal     = N_VInvTestLocal_ManyVector;
  v->ops->nvconstrmasklocal  = N_VConstrMaskLocal_ManyVector;
  v->ops->nvminquotientlocal = N_VMinQuotientLocal_ManyVector;
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_ManyVector;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_ManyVector;

  /* single buffer reduction operations */
  v->ops->nvdotprodmultilocal = N_VDotProdMultiLocal_ManyVector;

  /* XBraid interface operations */
  v->ops->nvbufsize   = N_VBufSize_ManyVector;
  v->ops->nvbufpack   = N_VBufPack_ManyVector;
  v->ops->nvbufunpack = N_VBufUnpack_ManyVector;

  /* debugging functions */
  v->ops->nvprint     = N_VPrint_ManyVector;
  v->ops->nvprintfile = N_VPrintFile_ManyVector;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_ManyVector)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  v->content = content;

  /* Attach content components */

  /* allocate and set subvector array */
  content->num_subvectors = num_subvectors;
  content->own_data       = SUNFALSE;

  content->subvec_array = NULL;
  content->subvec_array = (N_Vector*)malloc(num_subvectors * sizeof(N_Vector));
  SUNAssertNull(content->subvec_array, SUN_ERR_MALLOC_FAIL);

  for (i = 0; i < num_subvectors; i++)
  {
    content->subvec_array[i] = vec_array[i];
  }

  /* Determine overall ManyVector length: sum contributions from all subvectors */
  local_length = 0;
  for (i = 0; i < num_subvectors; i++)
  {
    SUNAssertNull(vec_array[i]->ops->nvgetlength, SUN_ERR_NOT_IMPLEMENTED);
    local_length += N_VGetLength(vec_array[i]);
    SUNCheckLastErrNull();
  }
  content->global_length = local_length;

  return (v);
}
#endif

/* This function returns the vec_num sub-N_Vector from the N_Vector
   array.  If vec_num is outside of applicable bounds, NULL is returned. */
N_Vector MVAPPEND(N_VGetSubvector)(N_Vector v, sunindextype vec_num)
{
  SUNFunctionBegin(v->sunctx);
  SUNAssertNull(vec_num >= 0, SUN_ERR_ARG_OUTOFRANGE);
  SUNAssertNull(vec_num < MANYVECTOR_NUM_SUBVECS(v), SUN_ERR_ARG_OUTOFRANGE);
  return (MANYVECTOR_SUBVEC(v, vec_num));
}

/* This function returns data pointer for the vec_num sub-N_Vector from
   the N_Vector array.  If vec_num is outside of applicable bounds, or if
   the subvector does not support the N_VGetArrayPointer routine, then
   NULL is returned. */
sunrealtype* MVAPPEND(N_VGetSubvectorArrayPointer)(N_Vector v,
                                                   sunindextype vec_num)
{
  SUNFunctionBegin(v->sunctx);
  SUNAssertNull(vec_num >= 0, SUN_ERR_ARG_OUTOFRANGE);
  SUNAssertNull(vec_num < MANYVECTOR_NUM_SUBVECS(v), SUN_ERR_ARG_OUTOFRANGE);
  sunrealtype* arr = NULL;
  if (MANYVECTOR_SUBVEC(v, vec_num)->ops->nvgetarraypointer)
  {
    arr = N_VGetArrayPointer(MANYVECTOR_SUBVEC(v, vec_num));
    SUNCheckLastErrNull();
  }
  return arr;
}

/* This function sets the data pointer for the vec_num sub-N_Vector from
   the N_Vector array.  If vec_num is outside of applicable bounds, or if
   the subvector does not support the N_VSetArrayPointer routine, then
   -1 is returned; otherwise this routine returns 0. */
SUNErrCode MVAPPEND(N_VSetSubvectorArrayPointer)(sunrealtype* v_data, N_Vector v,
                                                 sunindextype vec_num)
{
  SUNFunctionBegin(v->sunctx);
  SUNAssert(vec_num >= 0, SUN_ERR_ARG_OUTOFRANGE);
  SUNAssert(vec_num < MANYVECTOR_NUM_SUBVECS(v), SUN_ERR_ARG_OUTOFRANGE);
  N_VSetArrayPointer(v_data, MANYVECTOR_SUBVEC(v, vec_num));
  SUNCheckLastErr();
  return SUN_SUCCESS;
}

/* This function returns the overall number of sub-vectors.
   It returns a locally stored integer, and is therefore a local call. */
sunindextype MVAPPEND(N_VGetNumSubvectors)(N_Vector v)
{
  return (MANYVECTOR_NUM_SUBVECS(v));
}

/* -----------------------------------------------------------------
   ManyVector implementations of generic NVector routines
   -----------------------------------------------------------------*/

/* Returns vector type ID. Used to identify vector implementation
   from abstract N_Vector interface. */
N_Vector_ID MVAPPEND(N_VGetVectorID)(SUNDIALS_MAYBE_UNUSED N_Vector v)
{
#ifdef MANYVECTOR_BUILD_WITH_MPI
  return (SUNDIALS_NVEC_MPIMANYVECTOR);
#else
  return (SUNDIALS_NVEC_MANYVECTOR);
#endif
}

/* Prints the vector to stdout, calling Print on subvectors. */
void MVAPPEND(N_VPrint)(N_Vector x)
{
  sunindextype i;
  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    N_VPrint(MANYVECTOR_SUBVEC(x, i));
  }
  return;
}

/* Prints the vector to outfile, calling PrintFile on subvectors. */
void MVAPPEND(N_VPrintFile)(N_Vector x, FILE* outfile)
{
  sunindextype i;
  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    N_VPrintFile(MANYVECTOR_SUBVEC(x, i), outfile);
  }
  return;
}

/* Clones a ManyVector, calling CloneEmpty on subvectors. */
N_Vector MVAPPEND(N_VCloneEmpty)(N_Vector w)
{
  return (ManyVectorClone(w, SUNTRUE));
}

/* Clones a ManyVector, calling Clone on subvectors. */
N_Vector MVAPPEND(N_VClone)(N_Vector w)
{
  return (ManyVectorClone(w, SUNFALSE));
}

/* Destroys a ManyVector */
void MVAPPEND(N_VDestroy)(N_Vector v)
{
  SUNFunctionBegin(v->sunctx);

  sunindextype i;

  if (v == NULL) { return; }

  /* free content */
  if (v->content != NULL)
  {
    /* if subvectors are owned by v, then Destroy those */
    if ((MANYVECTOR_OWN_DATA(v) == SUNTRUE) && (MANYVECTOR_SUBVECS(v) != NULL))
    {
      for (i = 0; i < MANYVECTOR_NUM_SUBVECS(v); i++)
      {
        N_VDestroy(MANYVECTOR_SUBVEC(v, i));
        MANYVECTOR_SUBVEC(v, i) = NULL;
      }
    }

    /* free subvector array */
    free(MANYVECTOR_SUBVECS(v));
    MANYVECTOR_SUBVECS(v) = NULL;

#ifdef MANYVECTOR_BUILD_WITH_MPI
    /* free communicator */
    if (MANYVECTOR_COMM(v) != MPI_COMM_NULL)
    {
      SUNCheckMPICallVoid(MPI_Comm_free(&(MANYVECTOR_COMM(v))));
    }
    MANYVECTOR_COMM(v) = MPI_COMM_NULL;
#endif

    /* free content structure */
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

/* Returns the space requirements for the ManyVector, by accumulating this
   information from all subvectors. */
void MVAPPEND(N_VSpace)(N_Vector v, sunindextype* lrw, sunindextype* liw)
{
  SUNFunctionBegin(v->sunctx);
  sunindextype i, lrw1, liw1;
  *lrw = 0;
  *liw = 0;
  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(v); i++)
  {
    /* update space requirements for this subvector (if 'nvspace' is implemented) */
    if ((MANYVECTOR_SUBVEC(v, i))->ops->nvspace != NULL)
    {
      N_VSpace(MANYVECTOR_SUBVEC(v, i), &lrw1, &liw1);
      SUNCheckLastErrVoid();
      *lrw += lrw1;
      *liw += liw1;
    }
  }
  return;
}

#ifdef MANYVECTOR_BUILD_WITH_MPI
/* This function retrieves the MPI Communicator from an MPIManyVector object. */
MPI_Comm N_VGetCommunicator_MPIManyVector(N_Vector v)
{
  return (MANYVECTOR_COMM(v));
}
#endif

/* This function retrieves the global length of a ManyVector object. */
sunindextype MVAPPEND(N_VGetLength)(N_Vector v)
{
  return (MANYVECTOR_GLOBLENGTH(v));
}

sunindextype MVAPPEND(N_VGetSubvectorLocalLength)(N_Vector v, sunindextype vec_num)
{
  SUNFunctionBegin(v->sunctx);
  N_Vector subvector = MVAPPEND(N_VGetSubvector)(v, vec_num);
  SUNCheckLastErrNoRet();
  sunindextype loc_length = N_VGetLocalLength(subvector);
  SUNCheckLastErrNoRet();
  return loc_length;
}

/* Performs the linear sum z = a*x + b*y by calling N_VLinearSum on all subvectors;
   this routine does not check that x, y and z are ManyVectors, if they have the
   same number of subvectors, or if these subvectors are compatible. */
void MVAPPEND(N_VLinearSum)(sunrealtype a, N_Vector x, sunrealtype b,
                            N_Vector y, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i;
  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    N_VLinearSum(a, MANYVECTOR_SUBVEC(x, i), b, MANYVECTOR_SUBVEC(y, i),
                 MANYVECTOR_SUBVEC(z, i));
    SUNCheckLastErrVoid();
  }
  return;
}

/* Performs the operation z = c by calling N_VConst on all subvectors. */
void MVAPPEND(N_VConst)(sunrealtype c, N_Vector z)
{
  SUNFunctionBegin(z->sunctx);
  sunindextype i;
  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(z); i++)
  {
    N_VConst(c, MANYVECTOR_SUBVEC(z, i));
    SUNCheckLastErrVoid();
  }
  return;
}

/* Performs the operation z_j = x_j*y_j by calling N_VProd on all subvectors;
   this routine does not check that x, y and z are ManyVectors, if they have the
   same number of subvectors, or if these subvectors are compatible. */
void MVAPPEND(N_VProd)(N_Vector x, N_Vector y, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i;
  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    N_VProd(MANYVECTOR_SUBVEC(x, i), MANYVECTOR_SUBVEC(y, i),
            MANYVECTOR_SUBVEC(z, i));
    SUNCheckLastErrVoid();
  }
  return;
}

/* Performs the operation z_j = x_j/y_j by calling N_VDiv on all subvectors;
   this routine does not check that x, y and z are ManyVectors, if they have the
   same number of subvectors, or if these subvectors are compatible. */
void MVAPPEND(N_VDiv)(N_Vector x, N_Vector y, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i;
  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    N_VDiv(MANYVECTOR_SUBVEC(x, i), MANYVECTOR_SUBVEC(y, i),
           MANYVECTOR_SUBVEC(z, i));
    SUNCheckLastErrVoid();
  }
  return;
}

/* Performs the operation z_j = c*x_j by calling N_VScale on all subvectors;
   this routine does not check that x and z are ManyVectors, if they have the
   same number of subvectors, or if these subvectors are compatible. */
void MVAPPEND(N_VScale)(sunrealtype c, N_Vector x, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i;
  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    N_VScale(c, MANYVECTOR_SUBVEC(x, i), MANYVECTOR_SUBVEC(z, i));
    SUNCheckLastErrVoid();
  }
  return;
}

/* Performs the operation z_j = |x_j| by calling N_VAbs on all subvectors;
   this routine does not check that x and z are ManyVectors, if they have the
   same number of subvectors, or if these subvectors are compatible. */
void MVAPPEND(N_VAbs)(N_Vector x, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i;
  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    N_VAbs(MANYVECTOR_SUBVEC(x, i), MANYVECTOR_SUBVEC(z, i));
    SUNCheckLastErrVoid();
  }
  return;
}

/* Performs the operation z_j = 1/x_j by calling N_VInv on all subvectors;
   this routine does not check that x and z are ManyVectors, if they have the
   same number of subvectors, or if these subvectors are compatible. */
void MVAPPEND(N_VInv)(N_Vector x, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i;
  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    N_VInv(MANYVECTOR_SUBVEC(x, i), MANYVECTOR_SUBVEC(z, i));
    SUNCheckLastErrVoid();
  }
  return;
}

/* Performs the operation z_j = x_j + b by calling N_VAddConst on all subvectors;
   this routine does not check that x and z are ManyVectors, if they have the
   same number of subvectors, or if these subvectors are compatible. */
void MVAPPEND(N_VAddConst)(N_Vector x, sunrealtype b, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i;
  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    N_VAddConst(MANYVECTOR_SUBVEC(x, i), b, MANYVECTOR_SUBVEC(z, i));
    SUNCheckLastErrVoid();
  }
  return;
}

/* Performs the MPI task-local dot product of two ManyVectors by calling
   N_VDotProdLocal on all subvectors; this routine does not check that x and
   y are ManyVectors, if they have the same number of subvectors, or if these
   subvectors are compatible.

   If any subvector does not implement the N_VDotProdLocal routine (NULL
   function pointer), then this routine will call N_VDotProd, but only
   accumulate the sum if this is the root task for that subvector's
   communicator (note: serial vectors are always root task). */
sunrealtype MVAPPEND(N_VDotProdLocal)(N_Vector x, N_Vector y)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i;
  sunrealtype sum;
#ifdef MANYVECTOR_BUILD_WITH_MPI
  sunrealtype contrib;
  int rank;
#endif

  /* initialize output*/
  sum = ZERO;

  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
#ifdef MANYVECTOR_BUILD_WITH_MPI

    /* check for nvdotprodlocal in subvector */
    if (MANYVECTOR_SUBVEC(x, i)->ops->nvdotprodlocal)
    {
      sum += N_VDotProdLocal(MANYVECTOR_SUBVEC(x, i), MANYVECTOR_SUBVEC(y, i));
      SUNCheckLastErrNoRet();
      /* otherwise, call nvdotprod and root tasks accumulate to overall sum */
    }
    else
    {
      contrib = N_VDotProd(MANYVECTOR_SUBVEC(x, i), MANYVECTOR_SUBVEC(y, i));
      SUNCheckLastErrNoRet();

      /* get this task's rank in subvector communicator (note: serial
         subvectors will result in rank==0) */
      rank = SubvectorMPIRank(MANYVECTOR_SUBVEC(x, i));
      if (rank < 0) { return (ZERO); }
      if (rank == 0) { sum += contrib; }
    }

#else

    /* add subvector contribution */
    sum += N_VDotProd(MANYVECTOR_SUBVEC(x, i), MANYVECTOR_SUBVEC(y, i));
    SUNCheckLastErrNoRet();

#endif
  }

  return (sum);
}

#ifdef MANYVECTOR_BUILD_WITH_MPI
/* Performs the dot product of two ManyVectors by calling N_VDotProdLocal and
   combining the results.  This routine does not check that x and y are
   ManyVectors, if they have the same number of subvectors, or if these
   subvectors are compatible. */
sunrealtype N_VDotProd_MPIManyVector(N_Vector x, N_Vector y)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype lsum, gsum;
  lsum = gsum = N_VDotProdLocal_MPIManyVector(x, y);
  SUNCheckLastErrNoRet();
  if (MANYVECTOR_COMM(x) != MPI_COMM_NULL)
  {
    SUNCheckMPICallNoRet(MPI_Allreduce(&lsum, &gsum, 1, MPI_SUNREALTYPE,
                                       MPI_SUM, MANYVECTOR_COMM(x)));
  }
  return (gsum);
}
#endif

/* Performs the MPI task-local maximum norm of a ManyVector by calling
   N_VMaxNormLocal on all subvectors.

   If any subvector does not implement the N_VMaxNormLocal routine (NULL
   function pointer), then this routine will call N_VMaxNorm instead. */
sunrealtype MVAPPEND(N_VMaxNormLocal)(N_Vector x)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i;
  sunrealtype max, lmax;

  /* initialize output*/
  max = ZERO;

  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    /* check for nvmaxnormlocal in subvector */
    if (MANYVECTOR_SUBVEC(x, i)->ops->nvmaxnormlocal)
    {
      lmax = N_VMaxNormLocal(MANYVECTOR_SUBVEC(x, i));
      SUNCheckLastErrNoRet();
      max = (max > lmax) ? max : lmax;

      /* otherwise, call nvmaxnorm and accumulate to overall max */
    }
    else
    {
      lmax = N_VMaxNorm(MANYVECTOR_SUBVEC(x, i));
      SUNCheckLastErrNoRet();
      max = (max > lmax) ? max : lmax;
    }
  }

  return (max);
}

#ifdef MANYVECTOR_BUILD_WITH_MPI
/* Performs the maximum norm of a ManyVector by calling N_VMaxNormLocal and
   combining the results. */
sunrealtype N_VMaxNorm_MPIManyVector(N_Vector x)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype lmax, gmax;
  lmax = gmax = N_VMaxNormLocal_MPIManyVector(x);
  SUNCheckLastErrNoRet();
  if (MANYVECTOR_COMM(x) != MPI_COMM_NULL)
  {
    SUNCheckMPICallNoRet(MPI_Allreduce(&lmax, &gmax, 1, MPI_SUNREALTYPE,
                                       MPI_MAX, MANYVECTOR_COMM(x)));
  }
  return (gmax);
}
#endif

/* Performs the MPI task-local weighted squared sum of a ManyVector by calling
   N_VWSqrSumLocal on all subvectors; this routine does not check that x and
   w are ManyVectors, if they have the same number of subvectors, or if these
   subvectors are compatible.

   If any subvector does not implement the N_VWSqrSumLocal routine (NULL
   function pointer), then this routine will call N_VWrmsNorm and N_VGetLength
   to unravel the squared sum of the subvector components.  It will then only
   accumulate this to the overall sum if this is the root task for that
   subvector's communicator (note: serial vectors are always root task). */
sunrealtype MVAPPEND(N_VWSqrSumLocal)(N_Vector x, N_Vector w)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i, N;
  sunrealtype sum, contrib;
#ifdef MANYVECTOR_BUILD_WITH_MPI
  int rank;
#endif

  /* initialize output*/
  sum = ZERO;

  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
#ifdef MANYVECTOR_BUILD_WITH_MPI

    /* check for nvwsqrsumlocal in subvector */
    if (MANYVECTOR_SUBVEC(x, i)->ops->nvwsqrsumlocal)
    {
      sum += N_VWSqrSumLocal(MANYVECTOR_SUBVEC(x, i), MANYVECTOR_SUBVEC(w, i));
      SUNCheckLastErrNoRet();

      /* otherwise, call nvwrmsnorm, and accumulate to overall sum on root task */
    }
    else
    {
      contrib = N_VWrmsNorm(MANYVECTOR_SUBVEC(x, i), MANYVECTOR_SUBVEC(w, i));
      SUNCheckLastErrNoRet();

      /* get this task's rank in subvector communicator (note: serial
         subvectors will result in rank==0) */
      rank = SubvectorMPIRank(MANYVECTOR_SUBVEC(x, i));
      if (rank < 0) { return (ZERO); }
      if (rank == 0)
      {
        N = N_VGetLength(MANYVECTOR_SUBVEC(x, i));
        SUNCheckLastErrNoRet();
        sum += (contrib * contrib * N);
      }
    }

#else

    /* accumulate subvector contribution to overall sum */
    contrib = N_VWrmsNorm(MANYVECTOR_SUBVEC(x, i), MANYVECTOR_SUBVEC(w, i));
    SUNCheckLastErrNoRet();
    N = N_VGetLength(MANYVECTOR_SUBVEC(x, i));
    SUNCheckLastErrNoRet();
    sum += (contrib * contrib * N);

#endif
  }

  return (sum);
}

/* Performs the WRMS norm of a ManyVector by calling N_VWSqrSumLocal and
   combining the results; this routine does not check that x and
   w are ManyVectors, if they have the same number of subvectors, or if these
   subvectors are compatible. */
sunrealtype MVAPPEND(N_VWrmsNorm)(N_Vector x, N_Vector w)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype gsum;
#ifdef MANYVECTOR_BUILD_WITH_MPI
  sunrealtype lsum;
  lsum = gsum = N_VWSqrSumLocal_MPIManyVector(x, w);
  SUNCheckLastErrNoRet();
  if (MANYVECTOR_COMM(x) != MPI_COMM_NULL)
  {
    SUNCheckMPICallNoRet(MPI_Allreduce(&lsum, &gsum, 1, MPI_SUNREALTYPE,
                                       MPI_SUM, MANYVECTOR_COMM(x)));
  }
#else
  gsum = N_VWSqrSumLocal_ManyVector(x, w);
  SUNCheckLastErrNoRet();
#endif
  return (SUNRsqrt(gsum / (MANYVECTOR_GLOBLENGTH(x))));
}

/* Performs the MPI task-local masked weighted squared sum of a ManyVector by
   calling N_VWSqrSumMaskLocal on all subvectors; this routine does not check
   that x, w and id are ManyVectors, if they have the same number of
   subvectors, or if these subvectors are compatible.

   If any subvector does not implement the N_VWSqrSumLocal routine (NULL
   function pointer), then this routine will call N_VWrmsNormMask and
   N_VGetLength to unravel the masked squared sum of the subvector components.
   It will then only accumulate this to the overall sum if this is the root
   task for that subvector's communicator (note: serial vectors are always
   root task). */
sunrealtype MVAPPEND(N_VWSqrSumMaskLocal)(N_Vector x, N_Vector w, N_Vector id)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i, N;
  sunrealtype sum, contrib;
#ifdef MANYVECTOR_BUILD_WITH_MPI
  int rank;
#endif

  /* initialize output*/
  sum = ZERO;

  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
#ifdef MANYVECTOR_BUILD_WITH_MPI

    /* check for nvwsqrsummasklocal in subvector */
    if (MANYVECTOR_SUBVEC(x, i)->ops->nvwsqrsummasklocal)
    {
      sum += N_VWSqrSumMaskLocal(MANYVECTOR_SUBVEC(x, i), MANYVECTOR_SUBVEC(w, i),
                                 MANYVECTOR_SUBVEC(id, i));
      SUNCheckLastErrNoRet();

      /* otherwise, call nvwrmsnormmask, and accumulate to overall sum on root task */
    }
    else
    {
      contrib = N_VWrmsNormMask(MANYVECTOR_SUBVEC(x, i), MANYVECTOR_SUBVEC(w, i),
                                MANYVECTOR_SUBVEC(id, i));
      SUNCheckLastErrNoRet();

      /* get this task's rank in subvector communicator (note: serial
         subvectors will result in rank==0) */
      rank = SubvectorMPIRank(MANYVECTOR_SUBVEC(x, i));
      if (rank < 0) { return (ZERO); }
      if (rank == 0)
      {
        N = N_VGetLength(MANYVECTOR_SUBVEC(x, i));
        SUNCheckLastErrNoRet();
        sum += (contrib * contrib * N);
      }
    }

#else

    /* accumulate subvector contribution to overall sum */
    contrib = N_VWrmsNormMask(MANYVECTOR_SUBVEC(x, i), MANYVECTOR_SUBVEC(w, i),
                              MANYVECTOR_SUBVEC(id, i));
    SUNCheckLastErrNoRet();
    N = N_VGetLength(MANYVECTOR_SUBVEC(x, i));
    SUNCheckLastErrNoRet();
    sum += (contrib * contrib * N);

#endif
  }

  return (sum);
}

/* Performs the masked WRMS norm of a ManyVector by calling N_VWSqrSumMaskLocal
   and combining the results; this routine does not check that x, w and id are
   ManyVectors, if they have the same number of subvectors, or if these subvectors
   are compatible. */
sunrealtype MVAPPEND(N_VWrmsNormMask)(N_Vector x, N_Vector w, N_Vector id)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype gsum;
#ifdef MANYVECTOR_BUILD_WITH_MPI
  sunrealtype lsum;
  lsum = gsum = N_VWSqrSumMaskLocal_MPIManyVector(x, w, id);
  SUNCheckLastErrNoRet();
  if (MANYVECTOR_COMM(x) != MPI_COMM_NULL)
  {
    SUNCheckMPICallNoRet(MPI_Allreduce(&lsum, &gsum, 1, MPI_SUNREALTYPE,
                                       MPI_SUM, MANYVECTOR_COMM(x)));
  }
#else
  gsum = N_VWSqrSumMaskLocal_ManyVector(x, w, id);
  SUNCheckLastErrNoRet();
#endif
  return (SUNRsqrt(gsum / (MANYVECTOR_GLOBLENGTH(x))));
}

/* Computes the MPI task-local minimum entry of a ManyVector by calling
   N_VMinLocal on all subvectors.

   If any subvector does not implement the N_VMinLocal routine (NULL
   function pointer), then this routine will call N_VMin instead. */
sunrealtype MVAPPEND(N_VMinLocal)(N_Vector x)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i;
  sunrealtype min, lmin;

  /* initialize output*/
  min = SUN_BIG_REAL;

  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    /* check for nvminlocal in subvector */
    if (MANYVECTOR_SUBVEC(x, i)->ops->nvminlocal)
    {
      lmin = N_VMinLocal(MANYVECTOR_SUBVEC(x, i));
      SUNCheckLastErrNoRet();
      min = (min < lmin) ? min : lmin;

      /* otherwise, call nvmin and accumulate to overall min */
    }
    else
    {
      lmin = N_VMin(MANYVECTOR_SUBVEC(x, i));
      SUNCheckLastErrNoRet();
      min = (min < lmin) ? min : lmin;
    }
  }

  return (min);
}

#ifdef MANYVECTOR_BUILD_WITH_MPI
/* Computes the minimum entry of a ManyVector by calling N_VMinLocal and
   combining the results. */
sunrealtype N_VMin_MPIManyVector(N_Vector x)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype lmin, gmin;
  lmin = gmin = N_VMinLocal_MPIManyVector(x);
  SUNCheckLastErrNoRet();
  if (MANYVECTOR_COMM(x) != MPI_COMM_NULL)
  {
    SUNCheckMPICallNoRet(MPI_Allreduce(&lmin, &gmin, 1, MPI_SUNREALTYPE,
                                       MPI_MIN, MANYVECTOR_COMM(x)));
  }
  return (gmin);
}
#endif

/* Performs the WL2 norm of a ManyVector by calling N_VSqrSumLocal and
   'massaging' the result.  This routine does not check that x and w are
   ManyVectors, if they have the same number of subvectors, or if these
   subvectors are compatible. */
sunrealtype MVAPPEND(N_VWL2Norm)(N_Vector x, N_Vector w)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype gsum;
#ifdef MANYVECTOR_BUILD_WITH_MPI
  sunrealtype lsum;
  lsum = gsum = N_VWSqrSumLocal_MPIManyVector(x, w);
  SUNCheckLastErrNoRet();
  if (MANYVECTOR_COMM(x) != MPI_COMM_NULL)
  {
    SUNCheckMPICallNoRet(MPI_Allreduce(&lsum, &gsum, 1, MPI_SUNREALTYPE,
                                       MPI_SUM, MANYVECTOR_COMM(x)));
  }
#else
  gsum = N_VWSqrSumLocal_ManyVector(x, w);
  SUNCheckLastErrNoRet();
#endif
  return (SUNRsqrt(gsum));
}

/* Performs the MPI task-local L1 norm of a ManyVector by calling N_VL1NormLocal on
   all subvectors.  If any subvector does not implement the N_VL1NormLocal routine
   (NULL function pointer), then this routine will call N_VL1Norm, but only
   accumulate the sum if this is the root task for that subvector's
   communicator (note: serial vectors are always root task). */
sunrealtype MVAPPEND(N_VL1NormLocal)(N_Vector x)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i;
  sunrealtype sum;
#ifdef MANYVECTOR_BUILD_WITH_MPI
  sunrealtype contrib;
  int rank;
#endif

  /* initialize output*/
  sum = ZERO;

  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
#ifdef MANYVECTOR_BUILD_WITH_MPI

    /* check for nvl1normlocal in subvector */
    if (MANYVECTOR_SUBVEC(x, i)->ops->nvl1normlocal)
    {
      sum += N_VL1NormLocal(MANYVECTOR_SUBVEC(x, i));
      SUNCheckLastErrNoRet();

      /* otherwise, call nvl1norm and root tasks accumulate to overall sum */
    }
    else
    {
      contrib = N_VL1Norm(MANYVECTOR_SUBVEC(x, i));
      SUNCheckLastErrNoRet();

      /* get this task's rank in subvector communicator (note: serial
         subvectors will result in rank==0) */
      rank = SubvectorMPIRank(MANYVECTOR_SUBVEC(x, i));
      if (rank < 0) { return (ZERO); }
      if (rank == 0) { sum += contrib; }
    }

#else

    /* accumulate subvector contribution to overall sum */
    sum += N_VL1Norm(MANYVECTOR_SUBVEC(x, i));
    SUNCheckLastErrNoRet();

#endif
  }

  return (sum);
}

#ifdef MANYVECTOR_BUILD_WITH_MPI
/* Performs the L1 norm of a ManyVector by calling N_VL1NormLocal and
   combining the results. */
sunrealtype N_VL1Norm_MPIManyVector(N_Vector x)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype lsum, gsum;
  lsum = gsum = N_VL1NormLocal_MPIManyVector(x);
  SUNCheckLastErrNoRet();
  if (MANYVECTOR_COMM(x) != MPI_COMM_NULL)
  {
    SUNCheckMPICallNoRet(MPI_Allreduce(&lsum, &gsum, 1, MPI_SUNREALTYPE,
                                       MPI_SUM, MANYVECTOR_COMM(x)));
  }
  return (gsum);
}
#endif

/* Performs N_VCompare on all subvectors; this routine does not check that x and z are
   ManyVectors, if they have the same number of subvectors, or if these subvectors are
   compatible. */
void MVAPPEND(N_VCompare)(sunrealtype c, N_Vector x, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i;
  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    N_VCompare(c, MANYVECTOR_SUBVEC(x, i), MANYVECTOR_SUBVEC(z, i));
    SUNCheckLastErrVoid();
  }
  return;
}

/* Performs the MPI task-local InvTest for a ManyVector by calling N_VInvTestLocal
   on all subvectors and combining the results appropriately.  This routine does
   not check that x and z are ManyVectors, if they have the same number of
   subvectors, or if these subvectors are compatible.

   If any subvector does not implement the N_VInvTestLocal routine (NULL
   function pointer), then this routine will call N_VInvTest instead. */
sunbooleantype MVAPPEND(N_VInvTestLocal)(N_Vector x, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i;
  sunbooleantype val, subval;

  /* initialize output*/
  val = SUNTRUE;

  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    /* check for nvinvtestlocal in subvector */
    if (MANYVECTOR_SUBVEC(x, i)->ops->nvinvtestlocal)
    {
      subval = N_VInvTestLocal(MANYVECTOR_SUBVEC(x, i), MANYVECTOR_SUBVEC(z, i));
      SUNCheckLastErrNoRet();
      val = (val && subval);

      /* otherwise, call nvinvtest and accumulate to overall val */
    }
    else
    {
      subval = N_VInvTest(MANYVECTOR_SUBVEC(x, i), MANYVECTOR_SUBVEC(z, i));
      SUNCheckLastErrNoRet();
      val = (val && subval);
    }
  }

  return (val);
}

#ifdef MANYVECTOR_BUILD_WITH_MPI
/* Performs the InvTest for a ManyVector by calling N_VInvTestLocal and
   combining the results. This routine does not check that x and z
   are ManyVectors, if they have the same number of subvectors, or if these
   subvectors are compatible. */
sunbooleantype N_VInvTest_MPIManyVector(N_Vector x, N_Vector z)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype val, gval;
  sunbooleantype invtest = N_VInvTestLocal_MPIManyVector(x, z);
  SUNCheckLastErrNoRet();
  val = gval = (invtest) ? ONE : ZERO;
  if (MANYVECTOR_COMM(x) != MPI_COMM_NULL)
  {
    SUNCheckMPICallNoRet(MPI_Allreduce(&val, &gval, 1, MPI_SUNREALTYPE, MPI_MIN,
                                       MANYVECTOR_COMM(x)));
  }
  return (gval != ZERO);
}
#endif

/* Performs the MPI task-local ConstrMask for a ManyVector by calling N_VConstrMaskLocal
   on all subvectors and combining the results appropriately.  This routine does not
   check that c, x and m are ManyVectors, if they have the same number of subvectors,
   or if these subvectors are compatible.

   If any subvector does not implement the N_VConstrMaskLocal routine (NULL
   function pointer), then this routine will call N_VConstrMask instead. */
sunbooleantype MVAPPEND(N_VConstrMaskLocal)(N_Vector c, N_Vector x, N_Vector m)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i;
  sunbooleantype val, subval;

  /* initialize output*/
  val = SUNTRUE;

  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    /* check for nvconstrmasklocal in subvector */
    if (MANYVECTOR_SUBVEC(x, i)->ops->nvconstrmasklocal)
    {
      subval = N_VConstrMaskLocal(MANYVECTOR_SUBVEC(c, i),
                                  MANYVECTOR_SUBVEC(x, i),
                                  MANYVECTOR_SUBVEC(m, i));
      SUNCheckLastErrNoRet();
      val = (val && subval);

      /* otherwise, call nvconstrmask and accumulate to overall val */
    }
    else
    {
      subval = N_VConstrMask(MANYVECTOR_SUBVEC(c, i), MANYVECTOR_SUBVEC(x, i),
                             MANYVECTOR_SUBVEC(m, i));
      SUNCheckLastErrNoRet();
      val = (val && subval);
    }
  }

  return (val);
}

#ifdef MANYVECTOR_BUILD_WITH_MPI
/* Performs the ConstrMask for a ManyVector by calling N_VConstrMaskLocal and
   combining the results.  This routine does not check that c, x and m
   are ManyVectors, if they have the same number of subvectors, or if these
   subvectors are compatible. */
sunbooleantype N_VConstrMask_MPIManyVector(N_Vector c, N_Vector x, N_Vector m)
{
  SUNFunctionBegin(x->sunctx);
  sunrealtype val, gval;
  sunbooleantype constrmask = N_VConstrMaskLocal_MPIManyVector(c, x, m);
  SUNCheckLastErrNoRet();
  val = gval = (constrmask) ? ONE : ZERO;
  if (MANYVECTOR_COMM(x) != MPI_COMM_NULL)
  {
    SUNCheckMPICallNoRet(MPI_Allreduce(&val, &gval, 1, MPI_SUNREALTYPE, MPI_MIN,
                                       MANYVECTOR_COMM(x)));
  }
  return (gval != ZERO);
}
#endif

/* Performs the MPI task-local MinQuotient for a ManyVector by calling N_VMinQuotientLocal
   on all subvectors and combining the results appropriately.  This routine does not check
   that num and denom are ManyVectors, if they have the same number of subvectors, or if
   these subvectors are compatible.

   If any subvector does not implement the N_VMinQuotientLocal routine (NULL
   function pointer), then this routine will call N_VMinQuotient instead. */
sunrealtype MVAPPEND(N_VMinQuotientLocal)(N_Vector num, N_Vector denom)
{
  SUNFunctionBegin(num->sunctx);
  sunindextype i;
  sunrealtype min, lmin;

  /* initialize output*/
  min = SUN_BIG_REAL;

  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(num); i++)
  {
    /* check for nvminquotientlocal in subvector */
    if (MANYVECTOR_SUBVEC(num, i)->ops->nvminquotientlocal)
    {
      lmin = N_VMinQuotientLocal(MANYVECTOR_SUBVEC(num, i),
                                 MANYVECTOR_SUBVEC(denom, i));
      SUNCheckLastErrNoRet();
      min = (min < lmin) ? min : lmin;

      /* otherwise, call nvmin and accumulate to overall min */
    }
    else
    {
      lmin = N_VMinQuotient(MANYVECTOR_SUBVEC(num, i),
                            MANYVECTOR_SUBVEC(denom, i));
      SUNCheckLastErrNoRet();
      min = (min < lmin) ? min : lmin;
    }
  }

  return (min);
}

#ifdef MANYVECTOR_BUILD_WITH_MPI
/* Performs the MinQuotient for a ManyVector by calling N_VMinQuotientLocal
   and combining the results.  This routine does not check that num and
   denom are ManyVectors, if they have the same number of subvectors, or if
   these subvectors are compatible. */
sunrealtype N_VMinQuotient_MPIManyVector(N_Vector num, N_Vector denom)
{
  SUNFunctionBegin(num->sunctx);
  sunrealtype lmin, gmin;
  lmin = gmin = N_VMinQuotientLocal_MPIManyVector(num, denom);
  SUNCheckLastErrNoRet();
  if (MANYVECTOR_COMM(num) != MPI_COMM_NULL)
  {
    SUNCheckMPICallNoRet(MPI_Allreduce(&lmin, &gmin, 1, MPI_SUNREALTYPE,
                                       MPI_MIN, MANYVECTOR_COMM(num)));
  }
  return (gmin);
}
#endif

/* -----------------------------------------------------------------
   Single buffer reduction operations
   ----------------------------------------------------------------- */

SUNErrCode MVAPPEND(N_VDotProdMultiLocal)(int nvec, N_Vector x, N_Vector* Y,
                                          sunrealtype* dotprods)
{
  SUNFunctionBegin(x->sunctx);
  int j;
  sunindextype i;
  N_Vector* Ysub;
  sunrealtype* contrib;

  /* create temporary workspace arrays */
  Ysub = NULL;
  Ysub = (N_Vector*)malloc(nvec * sizeof(N_Vector));
  SUNAssert(Ysub, SUN_ERR_MALLOC_FAIL);

  contrib = NULL;
  contrib = (sunrealtype*)malloc(nvec * sizeof(sunrealtype));
  SUNAssert(contrib, SUN_ERR_MALLOC_FAIL);

  /* initialize output */
  for (j = 0; j < nvec; j++) { dotprods[j] = ZERO; }

  /* loop over subvectors */
  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    /* extract subvectors from vector array */
    for (j = 0; j < nvec; j++) { Ysub[j] = MANYVECTOR_SUBVEC(Y[j], i); }

    /* compute dot products */
    SUNCheckCall(
      N_VDotProdMultiLocal(nvec, MANYVECTOR_SUBVEC(x, i), Ysub, contrib));

    /* accumulate contributions */
    for (j = 0; j < nvec; j++) { dotprods[j] += contrib[j]; }
  }

  free(Ysub);
  free(contrib);

  /* return with success */
  return SUN_SUCCESS;
}

#ifdef MANYVECTOR_BUILD_WITH_MPI
SUNErrCode N_VDotProdMultiAllReduce_MPIManyVector(int nvec_total, N_Vector x,
                                                  sunrealtype* sum)
{
  SUNFunctionBegin(x->sunctx);

  /* accumulate totals and return */
  if (MANYVECTOR_COMM(x) == MPI_COMM_NULL) { return SUN_ERR_ARG_CORRUPT; }

  SUNCheckMPICall(MPI_Allreduce(MPI_IN_PLACE, sum, nvec_total, MPI_SUNREALTYPE,
                                MPI_SUM, MANYVECTOR_COMM(x)));

  return SUN_SUCCESS;
}
#endif

/* -----------------------------------------------------------------
   Fused vector operations
   ----------------------------------------------------------------- */

/* Performs the linear combination z = sum_j c[j]*X[j] by calling
   N_VLinearCombination on all subvectors; this routine does not check that z
   or the components of X are ManyVectors, if they have the same number of
   subvectors, or if these subvectors are compatible.

   NOTE: implementation of this routine is more challenging, due to the
   array-of-arrays of N_Vectors that comprise X.  This routine will be
   passed an array of ManyVectors, so to call the subvector-specific routines
   we must unravel the subvectors while retaining an array of outer vectors. */
SUNErrCode MVAPPEND(N_VLinearCombination)(int nvec, sunrealtype* c, N_Vector* X,
                                          N_Vector z)
{
  SUNFunctionBegin(z->sunctx);
  sunindextype i, j;
  N_Vector* Xsub;

  /* create array of nvec N_Vector pointers for reuse within loop */
  Xsub = NULL;
  Xsub = (N_Vector*)malloc(nvec * sizeof(N_Vector));
  SUNAssert(Xsub, SUN_ERR_MALLOC_FAIL);

  /* perform operation by calling N_VLinearCombination for each subvector */
  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(z); i++)
  {
    /* for each subvector, create the array of subvectors of X */
    for (j = 0; j < nvec; j++) { Xsub[j] = MANYVECTOR_SUBVEC(X[j], i); }

    /* now call N_VLinearCombination for this array of subvectors */
    SUNCheckCall(N_VLinearCombination(nvec, c, Xsub, MANYVECTOR_SUBVEC(z, i)));
  }

  /* clean up and return */
  free(Xsub);
  return SUN_SUCCESS;
}

/* Performs the ScaleAddMulti operation by calling N_VScaleAddMulti on all
   subvectors; this routine does not check that x, or the components of X and Z are
   ManyVectors, if they have the same number of subvectors, or if these subvectors
   are compatible.

   NOTE: this routine is more challenging, due to the array-of-arrays of
   N_Vectors that comprise Y and Z.  This routine will be passed an array of
   ManyVectors, so to call the subvector-specific routines we must unravel
   the subvectors while retaining an array of outer vectors. */
SUNErrCode MVAPPEND(N_VScaleAddMulti)(int nvec, sunrealtype* a, N_Vector x,
                                      N_Vector* Y, N_Vector* Z)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i, j;
  N_Vector *Ysub, *Zsub;

  /* create arrays of nvec N_Vector pointers for reuse within loop */
  Ysub = Zsub = NULL;
  Ysub        = (N_Vector*)malloc(nvec * sizeof(N_Vector));
  SUNAssert(Ysub, SUN_ERR_MALLOC_FAIL);
  Zsub = (N_Vector*)malloc(nvec * sizeof(N_Vector));
  SUNAssert(Zsub, SUN_ERR_MALLOC_FAIL);

  /* perform operation by calling N_VScaleAddMulti for each subvector */
  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    /* for each subvector, create the array of subvectors of Y and Z */
    for (j = 0; j < nvec; j++)
    {
      Ysub[j] = MANYVECTOR_SUBVEC(Y[j], i);
      Zsub[j] = MANYVECTOR_SUBVEC(Z[j], i);
    }

    /* now call N_VScaleAddMulti for this array of subvectors */
    SUNCheckCall(N_VScaleAddMulti(nvec, a, MANYVECTOR_SUBVEC(x, i), Ysub, Zsub));
  }

  /* clean up and return */
  free(Ysub);
  free(Zsub);
  return SUN_SUCCESS;
}

/* Performs the DotProdMulti operation by calling N_VDotProdLocal and combining results.
   This routine does not check that x, or the components of Y, are ManyVectors, if
   they have the same number of subvectors, or if these subvectors are compatible.

   NOTE: if all subvectors implement the N_VDotProdLocal routine, then this routine
   will require only a single array-valued reduction operation (in contrast to calling
   N_VDotProdMulti on all subvectors, where we would require num_subvectors separate
   reductions). */
SUNErrCode MVAPPEND(N_VDotProdMulti)(int nvec, N_Vector x, N_Vector* Y,
                                     sunrealtype* dotprods)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype i;

  /* call N_VDotProdLocal for each <x,Y[i]> pair */
  for (i = 0; i < nvec; i++)
  {
    dotprods[i] = N_VDotProdLocal(x, Y[i]);
    SUNCheckLastErr();
  }

#ifdef MANYVECTOR_BUILD_WITH_MPI
  /* accumulate totals and return */
  if (MANYVECTOR_COMM(x) != MPI_COMM_NULL)
  {
    SUNCheckMPICall(MPI_Allreduce(MPI_IN_PLACE, dotprods, nvec, MPI_SUNREALTYPE,
                                  MPI_SUM, MANYVECTOR_COMM(x)));
  }
#endif

  /* return with success */
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
   Vector array operations
   ----------------------------------------------------------------- */

/* Performs the LinearSumVectorArray operation by calling N_VLinearSumVectorArray
   on all subvectors; this routine does not check that the components of X, Y or Z are
   ManyVectors, if they have the same number of subvectors, or if these subvectors
   are compatible.

   NOTE: this routine is more challenging, due to the array-of-arrays of
   N_Vectors that comprise X, Y and Z.  This routine will be passed arrays of
   ManyVectors, so to call the subvector-specific routines we must unravel
   the subvectors while retaining arrays of outer vectors. */
SUNErrCode MVAPPEND(N_VLinearSumVectorArray)(int nvec, sunrealtype a,
                                             N_Vector* X, sunrealtype b,
                                             N_Vector* Y, N_Vector* Z)
{
  SUNFunctionBegin(X[0]->sunctx);
  sunindextype i, j;
  N_Vector *Xsub, *Ysub, *Zsub;

  SUNAssert(nvec > 0, SUN_ERR_ARG_OUTOFRANGE);

  /* create arrays of nvec N_Vector pointers for reuse within loop */
  Xsub = NULL;
  Ysub = NULL;
  Zsub = NULL;
  Xsub = (N_Vector*)malloc(nvec * sizeof(N_Vector));
  SUNAssert(Xsub, SUN_ERR_MALLOC_FAIL);
  Ysub = (N_Vector*)malloc(nvec * sizeof(N_Vector));
  SUNAssert(Ysub, SUN_ERR_MALLOC_FAIL);
  Zsub = (N_Vector*)malloc(nvec * sizeof(N_Vector));
  SUNAssert(Zsub, SUN_ERR_MALLOC_FAIL);

  /* perform operation by calling N_VLinearSumVectorArray for each subvector */
  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(X[0]); i++)
  {
    /* for each subvector, create the array of subvectors of X, Y and Z */
    for (j = 0; j < nvec; j++)
    {
      Xsub[j] = MANYVECTOR_SUBVEC(X[j], i);
      Ysub[j] = MANYVECTOR_SUBVEC(Y[j], i);
      Zsub[j] = MANYVECTOR_SUBVEC(Z[j], i);
    }

    /* now call N_VLinearSumVectorArray for this array of subvectors */
    SUNCheckCall(N_VLinearSumVectorArray(nvec, a, Xsub, b, Ysub, Zsub));
  }

  /* clean up and return */
  free(Xsub);
  free(Ysub);
  free(Zsub);
  return SUN_SUCCESS;
}

/* Performs the ScaleVectorArray operation by calling N_VScaleVectorArray
   on all subvectors; this routine does not check that the components of X or Z are
   ManyVectors, if they have the same number of subvectors, or if these subvectors
   are compatible.

   NOTE: this routine is more challenging, due to the array-of-arrays of
   N_Vectors that comprise X and Z.  This routine will be passed arrays of
   ManyVectors, so to call the subvector-specific routines we must unravel
   the subvectors while retaining arrays of outer vectors. */
SUNErrCode MVAPPEND(N_VScaleVectorArray)(int nvec, sunrealtype* c, N_Vector* X,
                                         N_Vector* Z)
{
  SUNFunctionBegin(X[0]->sunctx);
  sunindextype i, j;
  N_Vector *Xsub, *Zsub;

  SUNAssert(nvec > 0, SUN_ERR_ARG_OUTOFRANGE);

  /* create arrays of nvec N_Vector pointers for reuse within loop */
  Xsub = NULL;
  Zsub = NULL;
  Xsub = (N_Vector*)malloc(nvec * sizeof(N_Vector));
  SUNAssert(Xsub, SUN_ERR_MALLOC_FAIL);
  Zsub = (N_Vector*)malloc(nvec * sizeof(N_Vector));
  SUNAssert(Zsub, SUN_ERR_MALLOC_FAIL);

  /* perform operation by calling N_VScaleVectorArray for each subvector */
  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(X[0]); i++)
  {
    /* for each subvector, create the array of subvectors of X, Y and Z */
    for (j = 0; j < nvec; j++)
    {
      Xsub[j] = MANYVECTOR_SUBVEC(X[j], i);
      Zsub[j] = MANYVECTOR_SUBVEC(Z[j], i);
    }

    /* now call N_VScaleVectorArray for this array of subvectors */
    SUNCheckCall(N_VScaleVectorArray(nvec, c, Xsub, Zsub));
  }

  /* clean up and return */
  free(Xsub);
  free(Zsub);
  return SUN_SUCCESS;
}

/* Performs the ConstVectorArray operation by calling N_VConstVectorArray
   on all subvectors.

   NOTE: this routine is more challenging, due to the array-of-arrays of
   N_Vectors that comprise Z.  This routine will be passed an array of
   ManyVectors, so to call the subvector-specific routines we must unravel
   the subvectors while retaining an array of outer vectors. */
SUNErrCode MVAPPEND(N_VConstVectorArray)(int nvec, sunrealtype c, N_Vector* Z)
{
  SUNFunctionBegin(Z[0]->sunctx);
  sunindextype i, j;
  N_Vector* Zsub;

  SUNAssert(nvec > 0, SUN_ERR_ARG_OUTOFRANGE);

  /* create array of N_Vector pointers for reuse within loop */
  Zsub = NULL;
  Zsub = (N_Vector*)malloc(nvec * sizeof(N_Vector));
  SUNAssert(Zsub, SUN_ERR_MALLOC_FAIL);

  /* perform operation by calling N_VConstVectorArray for each subvector */
  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(Z[0]); i++)
  {
    /* for each subvector, create the array of subvectors of X, Y and Z */
    for (j = 0; j < nvec; j++) { Zsub[j] = MANYVECTOR_SUBVEC(Z[j], i); }

    /* now call N_VConstVectorArray for this array of subvectors */
    SUNCheckCall(N_VConstVectorArray(nvec, c, Zsub));
  }

  /* clean up and return */
  free(Zsub);
  return SUN_SUCCESS;
}

/* Performs the WrmsNormVectorArray operation by calling N_VWSqrSumLocal and combining
   results.  This routine does not check that the components of X or W are ManyVectors, if
   they have the same number of subvectors, or if these subvectors are compatible.

   NOTE: if all subvectors implement the N_VWSqrSumLocal routine, then this routine
   will require only a single array-valued reduction operation (in contrast to calling
   N_VWrmsNormVectorArray on all subvectors, where we would require num_subvectors
   separate reductions). */
SUNErrCode MVAPPEND(N_VWrmsNormVectorArray)(int nvec, N_Vector* X, N_Vector* W,
                                            sunrealtype* nrm)
{
  SUNFunctionBegin(X[0]->sunctx);
  sunindextype i;

  SUNAssert(nvec > 0, SUN_ERR_ARG_OUTOFRANGE);

  /* call N_VWSqrSumLocal for each (X[i],W[i]) pair */
  for (i = 0; i < nvec; i++)
  {
    nrm[i] = N_VWSqrSumLocal(X[i], W[i]);
    SUNCheckLastErr();
  }

  /* accumulate totals */
#ifdef MANYVECTOR_BUILD_WITH_MPI
  if (MANYVECTOR_COMM(X[0]) != MPI_COMM_NULL)
  {
    SUNCheckMPICall(MPI_Allreduce(MPI_IN_PLACE, nrm, nvec, MPI_SUNREALTYPE,
                                  MPI_SUM, MANYVECTOR_COMM(X[0])));
  }
#endif

  /* finish off WRMS norms and return */
  for (i = 0; i < nvec; i++)
  {
    nrm[i] = SUNRsqrt(nrm[i] / (MANYVECTOR_GLOBLENGTH(X[i])));
  }

  return SUN_SUCCESS;
}

/* Performs the WrmsNormMaskVectorArray operation by calling N_VWSqrSumMaskLocal and
   combining results.  This routine does not check that id or the components of X and
   W are ManyVectors, if they have the same number of subvectors, or if these
   subvectors are compatible.

   NOTE: if all subvectors implement the N_VWSqrSumMaskLocal routine, then this
   routine will require only a single array-valued reduction operation (in contrast
   to calling N_VWrmsNormMaskVectorArray on all subvectors, where we would require
   num_subvectors separate reductions). */
SUNErrCode MVAPPEND(N_VWrmsNormMaskVectorArray)(int nvec, N_Vector* X,
                                                N_Vector* W, N_Vector id,
                                                sunrealtype* nrm)
{
  SUNFunctionBegin(X[0]->sunctx);

  sunindextype i;

  SUNAssert(nvec > 0, SUN_ERR_ARG_OUTOFRANGE);

  /* call N_VWSqrSumMaskLocal for each (X[i],W[i]) pair */
  for (i = 0; i < nvec; i++)
  {
    nrm[i] = N_VWSqrSumMaskLocal(X[i], W[i], id);
    SUNCheckLastErr();
  }

  /* accumulate totals */
#ifdef MANYVECTOR_BUILD_WITH_MPI
  if (MANYVECTOR_COMM(X[0]) != MPI_COMM_NULL)
  {
    SUNCheckMPICall(MPI_Allreduce(MPI_IN_PLACE, nrm, nvec, MPI_SUNREALTYPE,
                                  MPI_SUM, MANYVECTOR_COMM(X[0])));
  }
#endif

  /* finish off WRMS norms and return */
  for (i = 0; i < nvec; i++)
  {
    nrm[i] = SUNRsqrt(nrm[i] / (MANYVECTOR_GLOBLENGTH(X[i])));
  }

  return SUN_SUCCESS;
}

/* Performs the BufSize operation by calling N_VBufSize for each subvector and
   combining results */
SUNErrCode MVAPPEND(N_VBufSize)(N_Vector x, sunindextype* size)
{
  SUNFunctionBegin(x->sunctx);
  sunindextype subvec_size; /* subvector buffer size */
  sunindextype i;

  /* initialize total size */
  *size = 0;

  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    /* get buffer sized needed for this subvector */
    SUNCheckCall(N_VBufSize(MANYVECTOR_SUBVEC(x, i), &subvec_size));

    /* update total buffer size */
    *size += subvec_size;
  }

  return SUN_SUCCESS;
}

/* Performs the BufPack operation by calling N_VBufPack for each subvector where
   the output buffer is offset by the buffer size used by the the previous
   subvector in the set */
SUNErrCode MVAPPEND(N_VBufPack)(N_Vector x, void* buf)
{
  SUNFunctionBegin(x->sunctx);
  void* loc;           /* location in output buffer */
  sunindextype offset; /* subvector buffer offset   */
  sunindextype i;

  SUNAssert(buf, SUN_ERR_ARG_CORRUPT);

  /* start at the beginning of the output buffer */
  loc = buf;

  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    /* pack the output buffer starting at the given buffer location */
    SUNCheckCall(N_VBufPack(MANYVECTOR_SUBVEC(x, i), loc));

    /* get the offset from this subvector */
    SUNCheckCall(N_VBufSize(MANYVECTOR_SUBVEC(x, i), &offset));

    /* update the buffer location for the next vector */
    loc = (char*)buf + offset;
  }

  return SUN_SUCCESS;
}

/* Performs the BufUnpack operation by calling N_VBufUnpack for each subvector
   where the input buffer is offset by the buffer size used by the the previous
   subvector in the set */
SUNErrCode MVAPPEND(N_VBufUnpack)(N_Vector x, void* buf)
{
  SUNFunctionBegin(x->sunctx);
  void* loc;           /* location in input buffer */
  sunindextype offset; /* subvector buffer offset   */
  sunindextype i;

  SUNAssert(buf, SUN_ERR_ARG_CORRUPT);

  /* start at the beginning of the input buffer */
  loc = buf;

  for (i = 0; i < MANYVECTOR_NUM_SUBVECS(x); i++)
  {
    /* unpack the input buffer starting at the given buffer location */
    SUNCheckCall(N_VBufUnpack(MANYVECTOR_SUBVEC(x, i), loc));

    /* get the offset from this subvector */
    SUNCheckCall(N_VBufSize(MANYVECTOR_SUBVEC(x, i), &offset));

    /* update the buffer location for the next vector */
    loc = (char*)buf + offset;
  }

  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
   Enable / Disable fused and vector array operations
   ----------------------------------------------------------------- */

SUNErrCode MVAPPEND(N_VEnableFusedOps)(N_Vector v, sunbooleantype tf)
{
  if (tf)
  {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = MVAPPEND(N_VLinearCombination);
    v->ops->nvscaleaddmulti     = MVAPPEND(N_VScaleAddMulti);
    v->ops->nvdotprodmulti      = MVAPPEND(N_VDotProdMulti);
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray     = MVAPPEND(N_VLinearSumVectorArray);
    v->ops->nvscalevectorarray         = MVAPPEND(N_VScaleVectorArray);
    v->ops->nvconstvectorarray         = MVAPPEND(N_VConstVectorArray);
    v->ops->nvwrmsnormvectorarray      = MVAPPEND(N_VWrmsNormVectorArray);
    v->ops->nvwrmsnormmaskvectorarray  = MVAPPEND(N_VWrmsNormMaskVectorArray);
    v->ops->nvscaleaddmultivectorarray = NULL;
    v->ops->nvlinearcombinationvectorarray = NULL;
    /* enable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = MVAPPEND(N_VDotProdMultiLocal);
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

SUNErrCode MVAPPEND(N_VEnableLinearCombination)(N_Vector v, sunbooleantype tf)
{
  /* enable/disable operation */
  if (tf) { v->ops->nvlinearcombination = MVAPPEND(N_VLinearCombination); }
  else { v->ops->nvlinearcombination = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode MVAPPEND(N_VEnableScaleAddMulti)(N_Vector v, sunbooleantype tf)
{
  /* enable/disable operation */
  if (tf) { v->ops->nvscaleaddmulti = MVAPPEND(N_VScaleAddMulti); }
  else { v->ops->nvscaleaddmulti = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode MVAPPEND(N_VEnableDotProdMulti)(N_Vector v, sunbooleantype tf)
{
  /* enable/disable operation */
  if (tf) { v->ops->nvdotprodmulti = MVAPPEND(N_VDotProdMulti); }
  else { v->ops->nvdotprodmulti = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode MVAPPEND(N_VEnableLinearSumVectorArray)(N_Vector v, sunbooleantype tf)
{
  /* enable/disable operation */
  if (tf)
  {
    v->ops->nvlinearsumvectorarray = MVAPPEND(N_VLinearSumVectorArray);
  }
  else { v->ops->nvlinearsumvectorarray = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode MVAPPEND(N_VEnableScaleVectorArray)(N_Vector v, sunbooleantype tf)
{
  /* enable/disable operation */
  if (tf) { v->ops->nvscalevectorarray = MVAPPEND(N_VScaleVectorArray); }
  else { v->ops->nvscalevectorarray = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode MVAPPEND(N_VEnableConstVectorArray)(N_Vector v, sunbooleantype tf)
{
  /* enable/disable operation */
  if (tf) { v->ops->nvconstvectorarray = MVAPPEND(N_VConstVectorArray); }
  else { v->ops->nvconstvectorarray = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode MVAPPEND(N_VEnableWrmsNormVectorArray)(N_Vector v, sunbooleantype tf)
{
  /* enable/disable operation */
  if (tf) { v->ops->nvwrmsnormvectorarray = MVAPPEND(N_VWrmsNormVectorArray); }
  else { v->ops->nvwrmsnormvectorarray = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode MVAPPEND(N_VEnableWrmsNormMaskVectorArray)(N_Vector v,
                                                      sunbooleantype tf)
{
  /* enable/disable operation */
  if (tf)
  {
    v->ops->nvwrmsnormmaskvectorarray = MVAPPEND(N_VWrmsNormMaskVectorArray);
  }
  else { v->ops->nvwrmsnormmaskvectorarray = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

SUNErrCode MVAPPEND(N_VEnableDotProdMultiLocal)(N_Vector v, sunbooleantype tf)
{
  /* enable/disable operation */
  if (tf) { v->ops->nvdotprodmultilocal = MVAPPEND(N_VDotProdMultiLocal); }
  else { v->ops->nvdotprodmultilocal = NULL; }

  /* return success */
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
   Implementation of utility routines
   -----------------------------------------------------------------*/

/* This function performs a generic clone operation on an input N_Vector.
   Based on the 'cloneempty' flag it will either call "nvclone" or
   "nvcloneempty" when creating subvectors in the cloned vector. */
static N_Vector ManyVectorClone(N_Vector w, sunbooleantype cloneempty)
{
  SUNFunctionBegin(w->sunctx);
  N_Vector v;
  MVAPPEND(N_VectorContent) content;
  sunindextype i;

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty(w->sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  SUNCheckCallNull(N_VCopyOps(w, v));

  /* Create content */
  content = NULL;
  content = (MVAPPEND(N_VectorContent))malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content and ops to new vector, and return */
  v->content = content;

  /* Attach content components */

  /* Set scalar components */
#ifdef MANYVECTOR_BUILD_WITH_MPI
  content->comm = MPI_COMM_NULL;
#endif
  content->num_subvectors = MANYVECTOR_NUM_SUBVECS(w);
  content->global_length  = MANYVECTOR_GLOBLENGTH(w);
  content->own_data       = SUNTRUE;

  /* Allocate the subvector array */
  content->subvec_array = NULL;
  content->subvec_array =
    (N_Vector*)malloc(content->num_subvectors * sizeof(N_Vector));
  SUNAssertNull(content->subvec_array, SUN_ERR_MALLOC_FAIL);

  /* Initialize the subvector array to NULL */
  for (i = 0; i < content->num_subvectors; i++)
  {
    content->subvec_array[i] = NULL;
  }

  /* Duplicate the input communicator (if applicable) */
#ifdef MANYVECTOR_BUILD_WITH_MPI
  if (MANYVECTOR_COMM(w) != MPI_COMM_NULL)
  {
    SUNCheckMPICallNull(MPI_Comm_dup(MANYVECTOR_COMM(w), &(content->comm)));
  }
#endif

  /* Clone vectors into the subvector array */
  for (i = 0; i < content->num_subvectors; i++)
  {
    if (cloneempty)
    {
      content->subvec_array[i] = N_VCloneEmpty(MANYVECTOR_SUBVEC(w, i));
      SUNCheckLastErrNull();
    }
    else
    {
      content->subvec_array[i] = N_VClone(MANYVECTOR_SUBVEC(w, i));
      SUNCheckLastErrNull();
    }
  }

  return (v);
}

#ifdef MANYVECTOR_BUILD_WITH_MPI
/* This function returns the rank of this task in the MPI communicator
   associated with the input N_Vector.  If the input N_Vector is MPI-unaware, it
   returns 0.  If an error occurs in the call to MPI_Comm_Rank, it returns -1. */
static int SubvectorMPIRank(N_Vector x)
{
  SUNFunctionBegin(x->sunctx);

  int rank      = 0;
  MPI_Comm comm = N_VGetCommunicator(x);
  SUNCheckLastErrNoRet();

  if (comm != MPI_COMM_NULL)
  {
    int status = MPI_Comm_rank(comm, &rank);
    SUNCheckMPICallNoRet(status);
    if (status != MPI_SUCCESS) { return -1; }
  }

  return rank;
}
#endif
