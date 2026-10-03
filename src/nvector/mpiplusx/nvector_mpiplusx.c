/* -----------------------------------------------------------------
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
 * This is the implementation file for the MPIPlusX NVECTOR.
 * -----------------------------------------------------------------*/

#include <nvector/nvector_mpiplusx_deprecated.h>
#include <sundials/priv/sundials_context_impl.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_errors.h>

#include "sundials_macros.h"

#define MPIPLUSX_LOCAL_VECTOR(v) (N_VGetSubvector_MPIManyVector(v, 0))

/* Functions attached to the N_Vector */
static sunrealtype* nvGetArrayPointer_MPIPlusX(N_Vector v);
static sunindextype nvGetLocalLength_MPIPlusX(N_Vector v);
static N_Vector_ID nvGetVectorID_MPIPlusX(N_Vector v);
static void nvPrintFile_MPIPlusX(N_Vector x, FILE* outfile);
static void nvPrint_MPIPlusX(N_Vector x);
static void nvSetArrayPointer_MPIPlusX(sunrealtype* vdata, N_Vector v);

N_Vector N_VMake_MPIPlusX(MPI_Comm comm, N_Vector X, SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);
  N_Vector v;

  SUNAssertNull(X, SUN_ERR_ARG_CORRUPT);

  v = NULL;
  v = N_VMake_MPIManyVector(comm, 1, &X, SUNCTX_);
  SUNCheckLastErrNull();

  /* override certain ops */
  v->ops->nvgetvectorid     = nvGetVectorID_MPIPlusX;
  v->ops->nvgetarraypointer = nvGetArrayPointer_MPIPlusX;
  v->ops->nvsetarraypointer = nvSetArrayPointer_MPIPlusX;
  v->ops->nvgetlocallength  = nvGetLocalLength_MPIPlusX;

  /* debugging functions */
  if (X->ops->nvprint) { v->ops->nvprint = nvPrint_MPIPlusX; }

  if (X->ops->nvprintfile) { v->ops->nvprintfile = nvPrintFile_MPIPlusX; }

  return v;
}

N_Vector_ID nvGetVectorID_MPIPlusX(SUNDIALS_MAYBE_UNUSED N_Vector v)
{
  return SUNDIALS_NVEC_MPIPLUSX;
}

sunrealtype* nvGetArrayPointer_MPIPlusX(N_Vector v)
{
  SUNFunctionBegin(v->sunctx);
  sunrealtype* arr = N_VGetSubvectorArrayPointer_MPIManyVector(v, 0);
  SUNCheckLastErrNull();
  return arr;
}

void nvSetArrayPointer_MPIPlusX(sunrealtype* vdata, N_Vector v)
{
  SUNFunctionBegin(v->sunctx);
  N_VSetSubvectorArrayPointer_MPIManyVector(vdata, v, 0);
  SUNCheckLastErrVoid();
}

N_Vector N_VGetLocalVector_MPIPlusX(N_Vector v)
{
  SUNFunctionBegin(v->sunctx);
  N_Vector result = N_VGetSubvector_MPIManyVector(v, 0);
  SUNCheckLastErrNull();
  return result;
}

sunindextype nvGetLocalLength_MPIPlusX(N_Vector v)
{
  SUNFunctionBegin(v->sunctx);
  N_Vector local_vector = N_VGetLocalVector_MPIPlusX(v);
  SUNCheckLastErrNoRet();
  sunindextype len = N_VGetLength(local_vector);
  SUNCheckLastErrNoRet();
  return len;
}

SUNErrCode N_VEnableFusedOps_MPIPlusX(N_Vector v, sunbooleantype tf)
{
  SUNFunctionBegin(v->sunctx);
  SUNCheckCall(N_VEnableFusedOps_MPIManyVector(v, tf));
  return SUN_SUCCESS;
}

void nvPrint_MPIPlusX(N_Vector v)
{
  N_Vector x = MPIPLUSX_LOCAL_VECTOR(v);
  if (x->ops->nvprint) { x->ops->nvprint(x); }
}

void nvPrintFile_MPIPlusX(N_Vector v, FILE* outfile)
{
  N_Vector x = MPIPLUSX_LOCAL_VECTOR(v);
  if (x->ops->nvprintfile) { x->ops->nvprintfile(x, outfile); }
}

/* Deprecated concrete operation wrappers */

sunrealtype* N_VGetArrayPointer_MPIPlusX(N_Vector v)
{
  return nvGetArrayPointer_MPIPlusX(v);
}

sunindextype N_VGetLocalLength_MPIPlusX(N_Vector v)
{
  return nvGetLocalLength_MPIPlusX(v);
}

N_Vector_ID N_VGetVectorID_MPIPlusX(N_Vector v)
{
  return nvGetVectorID_MPIPlusX(v);
}

void N_VPrintFile_MPIPlusX(N_Vector x, FILE* outfile)
{
  nvPrintFile_MPIPlusX(x, outfile);
}

void N_VPrint_MPIPlusX(N_Vector x) { nvPrint_MPIPlusX(x); }

void N_VSetArrayPointer_MPIPlusX(sunrealtype* vdata, N_Vector v)
{
  nvSetArrayPointer_MPIPlusX(vdata, v);
}
