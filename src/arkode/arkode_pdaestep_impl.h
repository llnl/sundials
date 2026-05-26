/*---------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 *---------------------------------------------------------------
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
 *---------------------------------------------------------------
 * This header defines the step memory for PDAEStep.
 *--------------------------------------------------------------*/

#ifndef ARKODE_PDAESTEP_IMPL_H_
#define ARKODE_PDAESTEP_IMPL_H_

#include "arkode/arkode_pdaestep.h"
#include "arkode_impl.h"

#define DEFAULT_MAX_COR 5
#define DEFAULT_ORDER 1
#define DEFAULT_NLSCOEF SUN_RCONST(0.1)

/* This struct wraps user data provided to the overall PDAEStep integrator and
 * includes additional members needed by an IDA instance
 */
typedef struct
{
  void *ark_mem;
  PDAEStepComponentResFn component_res_fn;
  PDAEStepLsComponentJacFn component_res_jac;

  /* Manyvectors with two componenets (diff/alg). These don't own the subvectors
   * but instead have pointers to the subvectors of PDAE state
   */
  N_Vector y;
  N_Vector yp;
} IDAUserData;

typedef struct
{
  PDAEStepAlgebraicResFn algebraic_res_fn;
  PDAEStepLsAlgebraicJacFn algebraic_res_jac;
  N_Vector yp;
  void **ida_mems;
  IDAUserData *user_datas;
  SUNNonlinearSolver NLS;

  sunrealtype nlscoef;

  long int nsetups;
  long int nls_iters;
  long int nls_fails;

  long int nh;
  
  int partitions;
  int order;
  int maxcor;
  
  sunbooleantype ownNLS;
}* ARKodePDAEStepMem;

/* Interface routines supplied to ARKODE */
int pdaeStep_GetNumRhsEvals(ARKodeMem ark_mem, int partition_index,
                            long int* rhs_evals);

/* Internal utility routines */
int pdaeStep_AccessARKODEStepMem(void* arkode_mem, const char* fname,
                                 ARKodeMem* ark_mem,
                                 ARKodePDAEStepMem* step_mem);
int pdaeStep_AccessStepMem(ARKodeMem ark_mem, const char* fname,
                           ARKodePDAEStepMem* step_mem);

#endif
