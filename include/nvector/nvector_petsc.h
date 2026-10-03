/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
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
 * This is the main header file for the PETSc vector wrapper
 * for NVECTOR module.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be
 *     found in the header file sundials_nvector.h.
 *
 *   - The definition of the type sunrealtype can be found in the
 *     header file sundials_types.h, and it may be changed (at the
 *     build configuration stage) according to the user's needs.
 *     The sundials_types.h file also contains the definition
 *     for the type sunbooleantype.
 *
 *   - N_Vector arguments to arithmetic vector operations need not
 *     be distinct. For example, the following call:
 *
 *        N_VLinearSum(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_NVECTOR_PETSC_H
#define SUNDIALS_NVECTOR_PETSC_H

#include <mpi.h>
#include <petscvec.h>
#include <sundials/sundials_mpi_types.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * PETSc implementation of N_Vector
 * -----------------------------------------------------------------
 */

struct _N_VectorContent_Petsc
{
  sunindextype local_length;  /* copy of local vector length  */
  sunindextype global_length; /* copy of global vector length */
  sunbooleantype own_data;    /* ownership of data            */
  Vec pvec;                   /* the PETSc Vec object         */
  MPI_Comm comm;              /* copy of MPI communicator     */
};

typedef struct _N_VectorContent_Petsc* N_VectorContent_Petsc;

/*
 * -----------------------------------------------------------------
 * Functions exported by nvector_petsc
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_Petsc(MPI_Comm comm,
                                           sunindextype local_length,
                                           sunindextype global_length,
                                           SUNContext sunctx);

SUNDIALS_EXPORT N_Vector N_VMake_Petsc(Vec v, SUNContext sunctx);

SUNDIALS_EXPORT Vec N_VGetVector_Petsc(N_Vector v);

SUNDIALS_EXPORT void N_VSetVector_Petsc(N_Vector v, Vec p);

/* nvector API functions */

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNErrCode N_VEnableFusedOps_Petsc(N_Vector v, sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableLinearCombination_Petsc(N_Vector v,
                                                            sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleAddMulti_Petsc(N_Vector v,
                                                        sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableDotProdMulti_Petsc(N_Vector v,
                                                       sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableLinearSumVectorArray_Petsc(N_Vector v,
                                                               sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleVectorArray_Petsc(N_Vector v,
                                                           sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableConstVectorArray_Petsc(N_Vector v,
                                                           sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableWrmsNormVectorArray_Petsc(N_Vector v,
                                                              sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormMaskVectorArray_Petsc(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMultiVectorArray_Petsc(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombinationVectorArray_Petsc(N_Vector v,
                                                       sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableDotProdMultiLocal_Petsc(N_Vector v,
                                                            sunbooleantype tf);

#ifdef __cplusplus
}
#endif

#endif /* SUNDIALS_NVECTOR_PETSC_H */
