/*------------------------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
 *------------------------------------------------------------------------------
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
 *------------------------------------------------------------------------------
 * This is the implementation file for ARKODE's partitioned DAE method
 *----------------------------------------------------------------------------*/

#include <arkode/arkode_pdaestep.h>
#include <nvector/nvector_manyvector.h>
#include <ida/ida.h>

#include "arkode_pdaestep_impl.h"
#include "arkode_impl.h"
#include "sundials_utils.h"

/*------------------------------------------------------------------------------
  Shortcut routine to unpack step_mem structure from ark_mem. If missing it
  returns ARK_MEM_NULL.
  ----------------------------------------------------------------------------*/
int pdaeStep_AccessStepMem(ARKodeMem ark_mem, const char* fname,
                           ARKodePDAEStepMem* step_mem)
{
  if (ark_mem->step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    "Time step module memory is NULL.");
    return ARK_MEM_NULL;
  }
  *step_mem = (ARKodePDAEStepMem)ark_mem->step_mem;
  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  Shortcut routine to unpack ark_mem and step_mem structures from void* pointer.
  If either is missing it returns ARK_MEM_NULL.
  ----------------------------------------------------------------------------*/
int pdaeStep_AccessARKODEStepMem(void* arkode_mem, const char* fname,
                                 ARKodeMem* ark_mem,
                                 ARKodePDAEStepMem* step_mem)
{
  /* access ARKodeMem structure */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_ARK_NO_MEM);
    return ARK_MEM_NULL;
  }
  *ark_mem = (ARKodeMem)arkode_mem;

  return pdaeStep_AccessStepMem(*ark_mem, __func__, step_mem);
}

/*------------------------------------------------------------------------------
  This routine is called just prior to performing internal time steps (after
  all user "set" routines have been called) from within arkInitialSetup.
  ----------------------------------------------------------------------------*/
static int pdaeStep_Init(ARKodeMem ark_mem,
                         SUNDIALS_MAYBE_UNUSED sunrealtype tout, SUNDIALS_MAYBE_UNUSED int init_type)
{
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  /* assume fixed outer step size */
  if (!ark_mem->fixedstep)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Adaptive outer time stepping is not currently supported");
    return ARK_ILL_INPUT;
  }

  if (ark_mem->interp_type == ARK_INTERP_HERMITE)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Hermite interpolation is not supported");
    return ARK_ILL_INPUT;
  }

  //TODO(SBR): Check if anything else needed here, e.g., checking init_type

  return ARK_SUCCESS;
}


/*------------------------------------------------------------------------------
  This routine performs a single step of the partitioned DAE method.
  ----------------------------------------------------------------------------*/
static int pdaeStep_TakeStep(ARKodeMem ark_mem, sunrealtype* dsmPtr,
                                int* nflagPtr)
{
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  *nflagPtr = ARK_SUCCESS; /* No algebraic solver */
  *dsmPtr   = ZERO;        /* No error estimate */

  //TODO(SBR): extend to more stages

  for (int i = 0; i < step_mem->partitions; i++) {
    void *ida_mem = step_mem->ida_mems[i];
    sunrealtype tout = ark_mem->tn + ark_mem->h;
    IDASetStopTime(ida_mem, tout);
    retval = IDASolve(ida_mem, tout, &tout, step_mem->user_datas[i].y,
            step_mem->user_datas[i].yp, IDA_NORMAL);

    if (retval != IDA_SUCCESS) {
      arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_INNERSTEP_FAILED, ark_mem->tcur);
      return ARK_INNERSTEP_FAIL;
    }
  }

  ark_mem->tcur += ark_mem->h;

  //TODO(SBR): perform nonlinear solve

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  Prints integrator statistics
  ----------------------------------------------------------------------------*/
static int pdaeStep_PrintAllStats(ARKodeMem ark_mem, FILE* outfile,
                                     SUNOutputFormat fmt)
{
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  // TODO(SBR): print stats

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  Frees all PDAEStep memory.
  ----------------------------------------------------------------------------*/
static void pdaeStep_Free(ARKodeMem ark_mem)
{
  ARKodePDAEStepMem step_mem = ark_mem->step_mem;
  if (step_mem != NULL) {
    N_VDestroy(step_mem->yp);
    for (int i = 0; i < step_mem->partitions; i++) {
      IDAFree(&step_mem->ida_mems[i]);
      N_VDestroy(step_mem->user_datas[i].y);
      N_VDestroy(step_mem->user_datas[i].yp);
    }
    free(step_mem->ida_mems);
    free(step_mem->user_datas);
    if (step_mem->ownNLS) {
      SUNNonlinSolFree(step_mem->NLS);
    }
    free(step_mem);
  }
  ark_mem->step_mem = NULL;
}

/*------------------------------------------------------------------------------
  This routine outputs the memory from the PDAEStep structure to a specified
  file pointer (useful when debugging).
  ----------------------------------------------------------------------------*/
static void pdaeStep_PrintMem(ARKodeMem ark_mem, FILE* outfile)
{
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return; }

  // TODO(SBR): Print stats and other props
}

static int pdaeStep_Residual(sunrealtype t, N_Vector y, N_Vector yp,
  N_Vector res, void *user_data)
{
  IDAUserData *ida_user_data = (IDAUserData*) user_data;
  void *arkode_user_data = ((ARKodeMem) ida_user_data->ark_mem)->user_data;

  // TODO(SBR): generalize w to polynomial. This is using w=w_n for now
  N_Vector w = PDAEStepGetCouplingSubvector(((ARKodeMem) ida_user_data->ark_mem)->yn);

  return ida_user_data->component_res_fn(t, y, w, yp, res, arkode_user_data);
}

static int pdaeStep_Jac(sunrealtype t, sunrealtype c_j, N_Vector y, N_Vector yp,
                        N_Vector r, SUNMatrix Jac, void *user_data,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  IDAUserData *ida_user_data = (IDAUserData*) user_data;
  void *arkode_user_data = ((ARKodeMem) ida_user_data->ark_mem)->user_data;
  // TODO(SBR): generalize w to polynomial. This is using w=w_n for now
  N_Vector w = PDAEStepGetCouplingSubvector(((ARKodeMem) ida_user_data->ark_mem)->yn);

  return ida_user_data->component_res_jac(t, c_j, y, w, yp, r, Jac, arkode_user_data);
}

/*------------------------------------------------------------------------------
  Creates the PDAEStep integrator
  ----------------------------------------------------------------------------*/
void* PDAEStepCreate(PDAEStepComponentResFn *componenet_res_fns,
                     PDAEStepAlgebraicResFn algebraic_res_fn,
                     sunrealtype t0,
                     N_Vector y0,
                     N_Vector yp0,
                     int partitions, SUNContext sunctx)
{
  if (componenet_res_fns == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "componenet_res_fns must be a non-NULL array");
    return (NULL);
  }
  for (int i = 0; i < partitions; i++) {
    if (componenet_res_fns[i] == NULL) {
      arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "componenet_res_fns[%i] must be a non-NULL", i);
      return (NULL);
    }
  }

  if (algebraic_res_fn == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "algebraic_res_fn must be non-NULL");
    return (NULL);
  }

  /* Check for legal input parameters */
  if (y0 == NULL || yp0 == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Initial conditions must be non-NULL");
    return (NULL);
  }

  if (sunctx == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_SUNCTX);
    return NULL;
  }

  /* Create ark_mem structure and set default values */
  ARKodeMem ark_mem = arkCreate(sunctx);
  if (ark_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return NULL;
  }

  ARKodePDAEStepMem step_mem = malloc(sizeof(*step_mem));
  if (step_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }

  // ark_mem->step_attachlinsol              = pdaeStep_AttachLinsol;
  ark_mem->step_init                      = pdaeStep_Init;
  ark_mem->step                           = pdaeStep_TakeStep;
  // ark_mem->step_setuserdata               = pdaeStep_SetUserData;
  ark_mem->step_printallstats             = pdaeStep_PrintAllStats;
  // ark_mem->step_writeparameters           = pdaeStep_WriteParameters;
  ark_mem->step_setusecompensatedsums     = NULL;
  ark_mem->step_free                      = pdaeStep_Free;
  ark_mem->step_printmem                  = pdaeStep_PrintMem;
  // ark_mem->step_setdefaults               = pdaeStep_SetDefaults;
  // ark_mem->step_setorder                  = pdaeStep_SetOrder;
  // ark_mem->step_setnonlinearsolver        = pdaeStep_SetNonlinearSolver;
  // ark_mem->step_setmaxnonliniters         = pdaeStep_SetMaxNonlinIters;
  // ark_mem->step_setnonlinconvcoef         = pdaeStep_SetNonlinConvCoef;
  // ark_mem->step_getnumrhsevals            = pdaeStep_GetNumRhsEvals;
  // ark_mem->step_getnumlinsolvsetups       = pdaeStep_GetNumLinSolvSetups;
  // ark_mem->step_getnumnonlinsolviters     = pdaeStep_GetNumNonlinSolvIters;
  // ark_mem->step_getnumnonlinsolvconvfails = pdaeStep_GetNumNonlinSolvConvFails;
  // ark_mem->step_getnonlinsolvstats        = pdaeStep_GetNonlinSolvStats;
  ark_mem->step_supports_adaptive         = SUNFALSE;
  /* The existing implicit infrastructure does not map very well to the
   * algebraic solve for PDAE methods. The dimension of the nonlinear system is
   * smaller than the overall problem, and there is effectively a 0 mass matrix.
   * Therefore, we set this to false and manually handle the algebraic solve.
   */
  ark_mem->step_supports_implicit         = SUNFALSE;
  ark_mem->step_supports_massmatrix       = SUNFALSE;
  ark_mem->step_supports_relaxation       = SUNFALSE;
  ark_mem->step_mem                       = (void*)step_mem;
  
  step_mem->algebraic_res_fn = algebraic_res_fn;
  // TODO(SBR): pedantic error checking
  step_mem->yp = N_VClone(yp0);
  step_mem->partitions = partitions;
  step_mem->ownNLS = SUNFALSE;

  int retval;
  //TODO(SBR): uncomment
  // int retval = pdaeStep_SetDefaults((void*)ark_mem);
  // if (retval != ARK_SUCCESS)
  // {
  //   arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
  //                   "Error setting default solver options");
  //   ARKodeFree((void**)&ark_mem);
  //   return NULL;
  // }

  ARKodeSetInterpolantType(ark_mem, ARK_INTERP_LAGRANGE);

  /* Initialize main ARKODE infrastructure */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to initialize main ARKODE infrastructure");
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }

  step_mem->ida_mems = malloc(partitions * sizeof(*step_mem->ida_mems));
  if (step_mem->ida_mems == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Unable to allocate IDA array");
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }
  step_mem->user_datas = malloc(partitions * sizeof(*step_mem->user_datas));
  if (step_mem->user_datas == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Unable to allocate user data array");
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }

  for (int i = 0; i < partitions; i++) {
    step_mem->ida_mems[i] = IDACreate(sunctx);
    if (step_mem->ida_mems[i] == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to allocate IDA for partition %i", i);
      ARKodeFree((void**)&ark_mem);
      return NULL;
    }

    {
      N_Vector yi_sub_vecs[] = {PDAEStepGetDifferentialSubvector(ark_mem->ycur, i),
                                PDAEStepGetAlgebraicSubvector(ark_mem->ycur, i)};
      N_Vector yi = N_VNew_ManyVector(2, yi_sub_vecs, sunctx);
      if (yi == NULL) {
        arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                        "Unable to create y N_Vector for partition %i", i);
        ARKodeFree((void**)&ark_mem);
        return NULL;
      }
      step_mem->user_datas[i].y = yi;
    }

    {
      N_Vector ypi_sub_vecs[] = {PDAEStepGetDifferentialSubvector(step_mem->yp, i),
                                PDAEStepGetAlgebraicSubvector(step_mem->yp, i)};
      N_Vector ypi = N_VNew_ManyVector(2, ypi_sub_vecs, sunctx);
      if (ypi == NULL) {
        arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                        "Unable to create yp N_Vector for partition %i", i);
        ARKodeFree((void**)&ark_mem);
        return NULL;
      }
      step_mem->user_datas[i].yp = ypi;
    }

    step_mem->user_datas[i].ark_mem = ark_mem;
    step_mem->user_datas[i].component_res_fn = componenet_res_fns[i];
    step_mem->user_datas[i].component_res_jac = NULL;

    IDAInit(step_mem->ida_mems[i], pdaeStep_Residual, t0,
      step_mem->user_datas[i].y, step_mem->user_datas[i].yp);
    IDASetUserData(step_mem->ida_mems[i], &step_mem->user_datas[i]);
    IDASetMaxNumSteps(step_mem->ida_mems[i], -1);
  }

  return ark_mem;
}

int PDAEStepGetNumPartitions(void *arkode_mem, int *partitions)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                                 &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (partitions == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "partitions is NULL");
    return ARK_ILL_INPUT;
  }

  *partitions = step_mem->partitions;

  return ARK_SUCCESS;
}

int PDAEStepGetPartitionIntegrator(void *arkode_mem, int partition, void **ida_mem)
{
  ARKodeMem ark_mem               = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                                 &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (partition < 0 || partition >= step_mem->partitions)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The partition index is %i but must be between 0 and %i",
                    partition, step_mem->partitions - 1);
    return ARK_ILL_INPUT;
  }

  if (ida_mem == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "ida_mem is NULL");
    return ARK_ILL_INPUT;
  }

  *ida_mem = step_mem->ida_mems[partition];

  return ARK_SUCCESS;
}

N_Vector PDAEStepManyVector(N_Vector *x, N_Vector *z, N_Vector w, int partitions)
{
  N_Vector *sub_vecs = malloc((2 * partitions + 1) * sizeof (*sub_vecs));
  if (sub_vecs == NULL) { return NULL; }
  for (int i = 0; i < partitions; i++) {
    sub_vecs[i] = x[i];
    sub_vecs[partitions + i] = z[i];
  }
  sub_vecs[2 * partitions] = w;
  N_Vector vec = N_VNew_ManyVector(2 * partitions + 1, sub_vecs, w->sunctx);
  
  free(sub_vecs);
  return vec;
}

N_Vector PDAEStepGetCouplingSubvector(N_Vector y)
{
  return N_VGetSubvector_ManyVector(y, N_VGetNumSubvectors_ManyVector(y) - 1);
}

N_Vector PDAEStepGetDifferentialSubvector(N_Vector y, int partition)
{
  return N_VGetSubvector_ManyVector(y, 2 * partition);
}

N_Vector PDAEStepGetAlgebraicSubvector(N_Vector y, int partition)
{
  return N_VGetSubvector_ManyVector(y, 2 * partition + 1);
}
