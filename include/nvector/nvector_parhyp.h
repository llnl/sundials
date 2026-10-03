/* -----------------------------------------------------------------
 * Programmer(s): Jean M. Sexton @ UMBC
 *                Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * Based on work by: Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                   and Aaron Collier @ LLNL
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
 * This is the main header file for the ParHyp implementation
 * of the NVECTOR module.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be
 *     found in the header file sundials_nvector.h.
 *
 *   - The definition of the type sunrealtype can be found in the
 *     header file sundials_types.h, and it may be changed (at the
 *     configuration stage) according to the user's needs.
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

#ifndef SUNDIALS_NVECTOR_PARHYP_H
#define SUNDIALS_NVECTOR_PARHYP_H

#include <mpi.h>
#include <stdio.h>
#include <sundials/sundials_mpi_types.h>
#include <sundials/sundials_nvector.h>

/* hypre header files */
#include <_hypre_parcsr_mv.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * ParHyp implementation of N_Vector
 * -----------------------------------------------------------------
 */

struct N_VectorContent_ParHyp_
{
  sunindextype local_length;    /* local vector length         */
  sunindextype global_length;   /* global vector length        */
  sunbooleantype own_parvector; /* ownership of HYPRE vector   */
  MPI_Comm comm;                /* pointer to MPI communicator */

  HYPRE_ParVector x; /* the actual HYPRE_ParVector object */
};

typedef struct N_VectorContent_ParHyp_* N_VectorContent_ParHyp;

/*
 * -----------------------------------------------------------------
 * Functions exported by nvector_parhyp
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_ParHyp(MPI_Comm comm,
                                            sunindextype local_length,
                                            sunindextype global_length,
                                            SUNContext sunctx);

SUNDIALS_EXPORT N_Vector N_VMake_ParHyp(HYPRE_ParVector x, SUNContext sunctx);

SUNDIALS_EXPORT HYPRE_ParVector N_VGetVector_ParHyp(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNErrCode N_VEnableFusedOps_ParHyp(N_Vector v,
                                                    sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableLinearCombination_ParHyp(N_Vector v,
                                                             sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleAddMulti_ParHyp(N_Vector v,
                                                         sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableDotProdMulti_ParHyp(N_Vector v,
                                                        sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableLinearSumVectorArray_ParHyp(N_Vector v,
                                                                sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableScaleVectorArray_ParHyp(N_Vector v,
                                                            sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableConstVectorArray_ParHyp(N_Vector v,
                                                            sunbooleantype tf);
SUNDIALS_EXPORT SUNErrCode N_VEnableWrmsNormVectorArray_ParHyp(N_Vector v,
                                                               sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableWrmsNormMaskVectorArray_ParHyp(N_Vector v, sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableScaleAddMultiVectorArray_ParHyp(N_Vector v,
                                                    sunbooleantype tf);
SUNDIALS_EXPORT
SUNErrCode N_VEnableLinearCombinationVectorArray_ParHyp(N_Vector v,
                                                        sunbooleantype tf);

SUNDIALS_EXPORT SUNErrCode N_VEnableDotProdMultiLocal_ParHyp(N_Vector v,
                                                             sunbooleantype tf);

#ifdef __cplusplus
}
#endif

#endif
