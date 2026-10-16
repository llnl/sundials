/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 *
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
 * This is the implementation file for a Trilinos implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_trilinos_deprecated.h>
#include <nvector/trilinos/SundialsTpetraVectorInterface.hpp>
#include <nvector/trilinos/SundialsTpetraVectorKernels.hpp>

#include "sundials_macros.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * using statements
 * -----------------------------------------------------------------
 */

using Teuchos::Comm;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::REDUCE_SUM;
using Teuchos::reduceAll;
using namespace sundials::trilinos::nvector_tpetra;

/*
 * -----------------------------------------------------------------
 * N_Vector operation declarations
 * -----------------------------------------------------------------
 */
static void nvAbs_Trilinos(N_Vector x, N_Vector z);
static void nvAddConst_Trilinos(N_Vector x, sunrealtype b, N_Vector z);
static N_Vector nvCloneEmpty_Trilinos(N_Vector w);
static N_Vector nvClone_Trilinos(N_Vector w);
static void nvCompare_Trilinos(sunrealtype c, N_Vector x, N_Vector z);
static void nvConst_Trilinos(sunrealtype c, N_Vector z);
static sunbooleantype nvConstrMaskLocal_Trilinos(N_Vector c, N_Vector x,
                                                 N_Vector m);
static sunbooleantype nvConstrMask_Trilinos(N_Vector c, N_Vector x, N_Vector m);
static void nvDestroy_Trilinos(N_Vector v);
static void nvDiv_Trilinos(N_Vector x, N_Vector y, N_Vector z);
static sunrealtype nvDotProdLocal_Trilinos(N_Vector x, N_Vector y);
static sunrealtype nvDotProd_Trilinos(N_Vector x, N_Vector y);
static SUNComm nvGetCommunicator_Trilinos(N_Vector v);
static sunindextype nvGetLength_Trilinos(N_Vector v);
static N_Vector_ID nvGetVectorID_Trilinos(N_Vector v);
static sunbooleantype nvInvTestLocal_Trilinos(N_Vector x, N_Vector z);
static sunbooleantype nvInvTest_Trilinos(N_Vector x, N_Vector z);
static void nvInv_Trilinos(N_Vector x, N_Vector z);
static sunrealtype nvL1NormLocal_Trilinos(N_Vector x);
static sunrealtype nvL1Norm_Trilinos(N_Vector x);
static void nvLinearSum_Trilinos(sunrealtype a, N_Vector x, sunrealtype b,
                                 N_Vector y, N_Vector z);
static sunrealtype nvMaxNormLocal_Trilinos(N_Vector x);
static sunrealtype nvMaxNorm_Trilinos(N_Vector x);
static sunrealtype nvMinLocal_Trilinos(N_Vector x);
static sunrealtype nvMinQuotientLocal_Trilinos(N_Vector num, N_Vector denom);
static sunrealtype nvMinQuotient_Trilinos(N_Vector num, N_Vector denom);
static sunrealtype nvMin_Trilinos(N_Vector x);
static void nvProd_Trilinos(N_Vector x, N_Vector y, N_Vector z);
static void nvScale_Trilinos(sunrealtype c, N_Vector x, N_Vector z);
static sunrealtype nvWL2Norm_Trilinos(N_Vector x, N_Vector w);
static sunrealtype nvWSqrSumLocal_Trilinos(N_Vector x, N_Vector w);
static sunrealtype nvWSqrSumMaskLocal_Trilinos(N_Vector x, N_Vector w,
                                               N_Vector id);
static sunrealtype nvWrmsNormMask_Trilinos(N_Vector x, N_Vector w, N_Vector id);
static sunrealtype nvWrmsNorm_Trilinos(N_Vector x, N_Vector w);

/*
 * -----------------------------------------------------------------
 * type definitions
 * -----------------------------------------------------------------
 */

typedef TpetraVectorInterface::vector_type vector_type;

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID nvGetVectorID_Trilinos(SUNDIALS_MAYBE_UNUSED N_Vector v)
{
  return SUNDIALS_NVEC_TRILINOS;
}

/* ----------------------------------------------------------------
 * Function to create a new Trilinos vector with empty data array
 */

N_Vector N_VNewEmpty_Trilinos(SUNContext sunctx)
{
  N_Vector v;

  /* Create an empty vector object */
  v = NULL;
  v = N_VNewEmpty(sunctx);
  if (v == NULL) { return (NULL); }

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = nvGetVectorID_Trilinos;
  v->ops->nvclone           = nvClone_Trilinos;
  v->ops->nvcloneempty      = nvCloneEmpty_Trilinos;
  v->ops->nvdestroy         = nvDestroy_Trilinos;
  v->ops->nvgetcommunicator = nvGetCommunicator_Trilinos;
  v->ops->nvgetlength       = nvGetLength_Trilinos;

  /* standard vector operations */
  v->ops->nvlinearsum    = nvLinearSum_Trilinos;
  v->ops->nvconst        = nvConst_Trilinos;
  v->ops->nvprod         = nvProd_Trilinos;
  v->ops->nvdiv          = nvDiv_Trilinos;
  v->ops->nvscale        = nvScale_Trilinos;
  v->ops->nvabs          = nvAbs_Trilinos;
  v->ops->nvinv          = nvInv_Trilinos;
  v->ops->nvaddconst     = nvAddConst_Trilinos;
  v->ops->nvdotprod      = nvDotProd_Trilinos;
  v->ops->nvmaxnorm      = nvMaxNorm_Trilinos;
  v->ops->nvwrmsnorm     = nvWrmsNorm_Trilinos;
  v->ops->nvwrmsnormmask = nvWrmsNormMask_Trilinos;
  v->ops->nvmin          = nvMin_Trilinos;
  v->ops->nvwl2norm      = nvWL2Norm_Trilinos;
  v->ops->nvl1norm       = nvL1Norm_Trilinos;
  v->ops->nvcompare      = nvCompare_Trilinos;
  v->ops->nvinvtest      = nvInvTest_Trilinos;
  v->ops->nvconstrmask   = nvConstrMask_Trilinos;
  v->ops->nvminquotient  = nvMinQuotient_Trilinos;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvdotprodlocal     = nvDotProdLocal_Trilinos;
  v->ops->nvmaxnormlocal     = nvMaxNormLocal_Trilinos;
  v->ops->nvminlocal         = nvMinLocal_Trilinos;
  v->ops->nvl1normlocal      = nvL1NormLocal_Trilinos;
  v->ops->nvinvtestlocal     = nvInvTestLocal_Trilinos;
  v->ops->nvconstrmasklocal  = nvConstrMaskLocal_Trilinos;
  v->ops->nvminquotientlocal = nvMinQuotientLocal_Trilinos;
  v->ops->nvwsqrsumlocal     = nvWSqrSumLocal_Trilinos;
  v->ops->nvwsqrsummasklocal = nvWSqrSumMaskLocal_Trilinos;

  return (v);
}

/* ----------------------------------------------------------------
 * Function to create an N_Vector attachment to Tpetra vector.
 * void* argument is to allow for calling this method from C code.
 *
 */

N_Vector N_VMake_Trilinos(Teuchos::RCP<vector_type> vec, SUNContext sunctx)
{
  N_Vector v = NULL;

  // Create an N_Vector with operators attached and empty content
  v = N_VNewEmpty_Trilinos(sunctx);
  if (v == NULL) return (NULL);

  // Create vector content using a pointer to Tpetra vector
  v->content = new TpetraVectorInterface(vec);
  if (v->content == NULL)
  {
    N_VDestroy(v);
    return NULL;
  }

  return (v);
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector nvCloneEmpty_Trilinos(N_Vector w)
{
  N_Vector v;

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

  return (v);
}

N_Vector nvClone_Trilinos(N_Vector w)
{
  N_Vector v = nvCloneEmpty_Trilinos(w);
  if (v == NULL) { return (NULL); }

  // Get raw pointer to Tpetra vector
  Teuchos::RCP<vector_type> wvec = N_VGetVector_Trilinos(w);

  // Clone wvec and get raw pointer to the clone
  Teuchos::RCP<vector_type> tvec =
    Teuchos::rcp(new vector_type(*wvec, Teuchos::Copy));

  // Create vector content using the raw pointer to the cloned Tpetra vector
  v->content = new TpetraVectorInterface(tvec);
  if (v->content == NULL)
  {
    N_VDestroy(v);
    return NULL;
  }

  return (v);
}

void nvDestroy_Trilinos(N_Vector v)
{
  if (v == NULL) { return; }

  if (v->content != NULL)
  {
    TpetraVectorInterface* iface =
      reinterpret_cast<TpetraVectorInterface*>(v->content);

    // iface was created with 'new', so use 'delete' to destroy it.
    delete iface;
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

/*
 * MPI communicator accessor
 */
SUNComm nvGetCommunicator_Trilinos(SUNDIALS_MAYBE_UNUSED N_Vector x)
{
#ifdef SUNDIALS_TRILINOS_HAVE_MPI
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  /* Access Teuchos::Comm* (which is actually a Teuchos::MpiComm*) */
  auto comm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(
    xv->getMap()->getComm());
  return (*(comm->getRawMpiComm().get())); /* extract MPI_Comm */
#else
  return (SUN_COMM_NULL);
#endif
}

/*
 * Global vector length accessor
 */
sunindextype nvGetLength_Trilinos(N_Vector x)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);

  return ((sunindextype)xv->getGlobalLength());
}

/*
 * Linear combination of two vectors: z = a*x + b*y
 */
void nvLinearSum_Trilinos(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                          N_Vector z)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> yv = N_VGetVector_Trilinos(y);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  if (x == z) { zv->update(b, *yv, a); }
  else if (y == z) { zv->update(a, *xv, b); }
  else { zv->update(a, *xv, b, *yv, ZERO); }
}

/*
 * Set all vector elements to a constant: z[i] = c
 */
void nvConst_Trilinos(sunrealtype c, N_Vector z)
{
  Teuchos::RCP<vector_type> zv = N_VGetVector_Trilinos(z);

  zv->putScalar(c);
}

/*
 * Elementwise multiply vectors: z[i] = x[i]*y[i]
 */
void nvProd_Trilinos(N_Vector x, N_Vector y, N_Vector z)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> yv = N_VGetVector_Trilinos(y);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  zv->elementWiseMultiply(ONE, *xv, *yv, ZERO);
}

/*
 * Elementwise divide vectors: z[i] = x[i]/y[i]
 */
void nvDiv_Trilinos(N_Vector x, N_Vector y, N_Vector z)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> yv = N_VGetVector_Trilinos(y);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  elementWiseDivide(*xv, *yv, *zv);
}

/*
 * Scale vector: z = c*x
 */
void nvScale_Trilinos(sunrealtype c, N_Vector x, N_Vector z)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  zv->scale(c, *xv);
}

/*
 * Elementwise absolute value: z[i] = |x[i]|
 */
void nvAbs_Trilinos(N_Vector x, N_Vector z)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  zv->abs(*xv);
}

/*
 * Elementwise inverse: z[i] = 1/x[i]
 */
void nvInv_Trilinos(N_Vector x, N_Vector z)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  zv->reciprocal(*xv);
}

/*
 * Add constant: z = x + b
 */
void nvAddConst_Trilinos(N_Vector x, sunrealtype b, N_Vector z)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  addConst(*xv, b, *zv);
}

/*
 * Scalar product of vectors x and y
 */
sunrealtype nvDotProd_Trilinos(N_Vector x, N_Vector y)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> yv = N_VGetVector_Trilinos(y);

  return xv->dot(*yv);
}

/*
 * Max norm (L infinity) of vector x
 */
sunrealtype nvMaxNorm_Trilinos(N_Vector x)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);

  return xv->normInf();
}

/*
 * Weighted RMS norm
 */
sunrealtype nvWrmsNorm_Trilinos(N_Vector x, N_Vector w)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> wv = N_VGetVector_Trilinos(w);

  return normWrms(*xv, *wv);
}

/*
 * Masked weighted RMS norm
 */
sunrealtype nvWrmsNormMask_Trilinos(N_Vector x, N_Vector w, N_Vector id)
{
  Teuchos::RCP<const vector_type> xv  = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> wv  = N_VGetVector_Trilinos(w);
  Teuchos::RCP<const vector_type> idv = N_VGetVector_Trilinos(id);

  return normWrmsMask(*xv, *wv, *idv);
}

/*
 * Returns minimum vector element
 */
sunrealtype nvMin_Trilinos(N_Vector x)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);

  return minElement(*xv);
}

/*
 * Weighted L2 norm
 */
sunrealtype nvWL2Norm_Trilinos(N_Vector x, N_Vector w)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> wv = N_VGetVector_Trilinos(w);

  return normWL2(*xv, *wv);
}

/*
 * L1 norm
 */
sunrealtype nvL1Norm_Trilinos(N_Vector x)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);

  return xv->norm1();
}

/*
 * Elementwise z[i] = |x[i]| >= c ? 1 : 0
 */
void nvCompare_Trilinos(sunrealtype c, N_Vector x, N_Vector z)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  compare(c, *xv, *zv);
}

/*
 * Elementwise inverse with zero checking: z[i] = 1/x[i], x[i] != 0
 */
sunbooleantype nvInvTest_Trilinos(N_Vector x, N_Vector z)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  return invTest(*xv, *zv) ? SUNTRUE : SUNFALSE;
}

/*
 * Checks constraint violations for vector x. Constraints are defined in
 * vector c, and constraint violation flags are stored in vector m.
 */
sunbooleantype nvConstrMask_Trilinos(N_Vector c, N_Vector x, N_Vector m)
{
  Teuchos::RCP<const vector_type> cv = N_VGetVector_Trilinos(c);
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<vector_type> mv       = N_VGetVector_Trilinos(m);

  return constraintMask(*cv, *xv, *mv) ? SUNTRUE : SUNFALSE;
}

/*
 * Find minimum quotient: minq  = min ( num[i]/denom[i]), denom[i] != 0.
 */
sunrealtype nvMinQuotient_Trilinos(N_Vector num, N_Vector denom)
{
  Teuchos::RCP<const vector_type> numv = N_VGetVector_Trilinos(num);
  Teuchos::RCP<const vector_type> denv = N_VGetVector_Trilinos(denom);

  return minQuotient(*numv, *denv);
}

/*
 * MPI task-local dot product
 */
sunrealtype nvDotProdLocal_Trilinos(N_Vector x, N_Vector y)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> yv = N_VGetVector_Trilinos(y);

  return dotProdLocal(*xv, *yv);
}

/*
 * MPI task-local maximum norm
 */
sunrealtype nvMaxNormLocal_Trilinos(N_Vector x)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);

  return maxNormLocal(*xv);
}

/*
 * MPI task-local minimum element
 */
sunrealtype nvMinLocal_Trilinos(N_Vector x)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);

  return minLocal(*xv);
}

/*
 * MPI task-local L1 norm
 */
sunrealtype nvL1NormLocal_Trilinos(N_Vector x)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);

  return L1NormLocal(*xv);
}

/*
 * MPI task-local weighted squared sum
 */
sunrealtype nvWSqrSumLocal_Trilinos(N_Vector x, N_Vector w)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> wv = N_VGetVector_Trilinos(w);

  return WSqrSumLocal(*xv, *wv);
}

/*
 * MPI task-local weighted masked squared sum
 */
sunrealtype nvWSqrSumMaskLocal_Trilinos(N_Vector x, N_Vector w, N_Vector id)
{
  Teuchos::RCP<const vector_type> xv  = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> wv  = N_VGetVector_Trilinos(w);
  Teuchos::RCP<const vector_type> idv = N_VGetVector_Trilinos(id);

  return WSqrSumMaskLocal(*xv, *wv, *idv);
}

/*
 * MPI task-local elementwise inverse with zero checking: z[i] = 1/x[i], x[i] != 0
 */
sunbooleantype nvInvTestLocal_Trilinos(N_Vector x, N_Vector z)
{
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  return invTestLocal(*xv, *zv) ? SUNTRUE : SUNFALSE;
}

/*
 * MPI task-local constraint checking for vector x. Constraints are defined in
 * vector c, and constraint violation flags are stored in vector m.
 */
sunbooleantype nvConstrMaskLocal_Trilinos(N_Vector c, N_Vector x, N_Vector m)
{
  Teuchos::RCP<const vector_type> cv = N_VGetVector_Trilinos(c);
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<vector_type> mv       = N_VGetVector_Trilinos(m);

  return constraintMaskLocal(*cv, *xv, *mv) ? SUNTRUE : SUNFALSE;
}

/*
 * MPI task-local minimum quotient: minq  = min ( num[i]/denom[i]), denom[i] != 0.
 */
sunrealtype nvMinQuotientLocal_Trilinos(N_Vector num, N_Vector denom)
{
  Teuchos::RCP<const vector_type> numv = N_VGetVector_Trilinos(num);
  Teuchos::RCP<const vector_type> denv = N_VGetVector_Trilinos(denom);

  return minQuotientLocal(*numv, *denv);
}

/* Deprecated concrete operation wrappers */

extern "C" {

void N_VAbs_Trilinos(N_Vector x, N_Vector z) { nvAbs_Trilinos(x, z); }

void N_VAddConst_Trilinos(N_Vector x, sunrealtype b, N_Vector z)
{
  nvAddConst_Trilinos(x, b, z);
}

N_Vector N_VCloneEmpty_Trilinos(N_Vector w) { return nvCloneEmpty_Trilinos(w); }

N_Vector N_VClone_Trilinos(N_Vector w) { return nvClone_Trilinos(w); }

void N_VCompare_Trilinos(sunrealtype c, N_Vector x, N_Vector z)
{
  nvCompare_Trilinos(c, x, z);
}

void N_VConst_Trilinos(sunrealtype c, N_Vector z) { nvConst_Trilinos(c, z); }

sunbooleantype N_VConstrMaskLocal_Trilinos(N_Vector c, N_Vector x, N_Vector m)
{
  return nvConstrMaskLocal_Trilinos(c, x, m);
}

sunbooleantype N_VConstrMask_Trilinos(N_Vector c, N_Vector x, N_Vector m)
{
  return nvConstrMask_Trilinos(c, x, m);
}

void N_VDestroy_Trilinos(N_Vector v) { nvDestroy_Trilinos(v); }

void N_VDiv_Trilinos(N_Vector x, N_Vector y, N_Vector z)
{
  nvDiv_Trilinos(x, y, z);
}

sunrealtype N_VDotProdLocal_Trilinos(N_Vector x, N_Vector y)
{
  return nvDotProdLocal_Trilinos(x, y);
}

sunrealtype N_VDotProd_Trilinos(N_Vector x, N_Vector y)
{
  return nvDotProd_Trilinos(x, y);
}

SUNComm N_VGetCommunicator_Trilinos(N_Vector v)
{
  return nvGetCommunicator_Trilinos(v);
}

sunindextype N_VGetLength_Trilinos(N_Vector v)
{
  return nvGetLength_Trilinos(v);
}

N_Vector_ID N_VGetVectorID_Trilinos(N_Vector v)
{
  return nvGetVectorID_Trilinos(v);
}

sunbooleantype N_VInvTestLocal_Trilinos(N_Vector x, N_Vector z)
{
  return nvInvTestLocal_Trilinos(x, z);
}

sunbooleantype N_VInvTest_Trilinos(N_Vector x, N_Vector z)
{
  return nvInvTest_Trilinos(x, z);
}

void N_VInv_Trilinos(N_Vector x, N_Vector z) { nvInv_Trilinos(x, z); }

sunrealtype N_VL1NormLocal_Trilinos(N_Vector x)
{
  return nvL1NormLocal_Trilinos(x);
}

sunrealtype N_VL1Norm_Trilinos(N_Vector x) { return nvL1Norm_Trilinos(x); }

void N_VLinearSum_Trilinos(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                           N_Vector z)
{
  nvLinearSum_Trilinos(a, x, b, y, z);
}

sunrealtype N_VMaxNormLocal_Trilinos(N_Vector x)
{
  return nvMaxNormLocal_Trilinos(x);
}

sunrealtype N_VMaxNorm_Trilinos(N_Vector x) { return nvMaxNorm_Trilinos(x); }

sunrealtype N_VMinLocal_Trilinos(N_Vector x) { return nvMinLocal_Trilinos(x); }

sunrealtype N_VMinQuotientLocal_Trilinos(N_Vector num, N_Vector denom)
{
  return nvMinQuotientLocal_Trilinos(num, denom);
}

sunrealtype N_VMinQuotient_Trilinos(N_Vector num, N_Vector denom)
{
  return nvMinQuotient_Trilinos(num, denom);
}

sunrealtype N_VMin_Trilinos(N_Vector x) { return nvMin_Trilinos(x); }

void N_VProd_Trilinos(N_Vector x, N_Vector y, N_Vector z)
{
  nvProd_Trilinos(x, y, z);
}

void N_VScale_Trilinos(sunrealtype c, N_Vector x, N_Vector z)
{
  nvScale_Trilinos(c, x, z);
}

sunrealtype N_VWL2Norm_Trilinos(N_Vector x, N_Vector w)
{
  return nvWL2Norm_Trilinos(x, w);
}

sunrealtype N_VWSqrSumLocal_Trilinos(N_Vector x, N_Vector w)
{
  return nvWSqrSumLocal_Trilinos(x, w);
}

sunrealtype N_VWSqrSumMaskLocal_Trilinos(N_Vector x, N_Vector w, N_Vector id)
{
  return nvWSqrSumMaskLocal_Trilinos(x, w, id);
}

sunrealtype N_VWrmsNormMask_Trilinos(N_Vector x, N_Vector w, N_Vector id)
{
  return nvWrmsNormMask_Trilinos(x, w, id);
}

sunrealtype N_VWrmsNorm_Trilinos(N_Vector x, N_Vector w)
{
  return nvWrmsNorm_Trilinos(x, w);
}

} // extern "C"
