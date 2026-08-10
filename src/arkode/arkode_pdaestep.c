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
#include <ida/ida.h>
#include <ida/ida_ls.h>
#include <nvector/nvector_manyvector.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#include <stdlib.h>
#include <string.h>

#include "arkode_impl.h"
#include "arkode_pdaestep_impl.h"
#include "sundials_utils.h"

static int pdaeStep_Jac(sunrealtype t, sunrealtype c_j, N_Vector y, N_Vector yp,
                        N_Vector r, SUNMatrix Jac, void* user_data,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int pdaeStep_JacTimes(sunrealtype t, N_Vector y, N_Vector yp, N_Vector r,
                             N_Vector v, N_Vector Jv, sunrealtype c_j,
                             void* user_data, N_Vector tmp1, N_Vector tmp2);
static int pdaeStep_ReInitPartitions(ARKodeMem ark_mem);

static N_Vector pdaeStep_CreateAlgebraicVector(N_Vector y, int partitions,
                                               SUNContext sunctx)
{
  N_Vector vec       = NULL;
  N_Vector* sub_vecs = NULL;

  sub_vecs = (N_Vector*)malloc((partitions + 1) * sizeof(*sub_vecs));
  if (sub_vecs == NULL) { return NULL; }

  for (int i = 0; i < partitions; i++)
  {
    sub_vecs[i] = PDAEStepGetAlgebraicSubvector(y, i);
  }
  sub_vecs[partitions] = PDAEStepGetCouplingSubvector(y);

  vec = N_VNew_ManyVector(partitions + 1, sub_vecs, sunctx);
  free(sub_vecs);

  return vec;
}

static int pdaeStep_AllocAlgebraicVectors(ARKodeMem ark_mem,
                                          ARKodePDAEStepMem step_mem,
                                          SUNContext sunctx)
{
  step_mem->alg = pdaeStep_CreateAlgebraicVector(ark_mem->ycur,
                                                 step_mem->partitions, sunctx);
  if (step_mem->alg == NULL) { return ARK_MEM_FAIL; }

  step_mem->alg_pred     = N_VClone(step_mem->alg);
  step_mem->alg_cor      = N_VClone(step_mem->alg);
  step_mem->alg_res      = N_VClone(step_mem->alg);
  step_mem->alg_fcur     = N_VClone(step_mem->alg);
  step_mem->alg_tmp      = N_VClone(step_mem->alg);
  step_mem->alg_linsol_x = N_VClone(step_mem->alg);
  step_mem->alg_ewt =
    pdaeStep_CreateAlgebraicVector(ark_mem->ewt, step_mem->partitions, sunctx);

  if ((step_mem->alg_pred == NULL) || (step_mem->alg_cor == NULL) ||
      (step_mem->alg_res == NULL) || (step_mem->alg_fcur == NULL) ||
      (step_mem->alg_tmp == NULL) || (step_mem->alg_linsol_x == NULL) ||
      (step_mem->alg_ewt == NULL))
  {
    return ARK_MEM_FAIL;
  }

  return ARK_SUCCESS;
}

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
                                 ARKodeMem* ark_mem, ARKodePDAEStepMem* step_mem)
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
static int pdaeStep_Init(ARKodeMem ark_mem, SUNDIALS_MAYBE_UNUSED int init_type)
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

  for (int i = 0; i < step_mem->partitions; i++)
  {
    if (step_mem->user_datas[i].component_res_jac != NULL)
    {
      retval = IDASetJacFn(step_mem->ida_mems[i], pdaeStep_Jac);
      if (retval != IDA_SUCCESS)
      {
        arkProcessError(ark_mem, ARK_INNERSTEP_ATTACH_ERR, __LINE__, __func__,
                        __FILE__,
                        "Unable to set IDA Jacobian function for partition %i",
                        i);
        return ARK_INNERSTEP_ATTACH_ERR;
      }
    }

    if (step_mem->user_datas[i].component_res_jtimes != NULL)
    {
      retval = IDASetJacTimes(step_mem->ida_mems[i], NULL, pdaeStep_JacTimes);
      if (retval != IDA_SUCCESS)
      {
        arkProcessError(ark_mem, ARK_INNERSTEP_ATTACH_ERR, __LINE__, __func__,
                        __FILE__,
                        "Unable to set IDA Jacobian-times function for "
                        "partition %i",
                        i);
        return ARK_INNERSTEP_ATTACH_ERR;
      }
    }
  }

  return pdaeStep_NlsInit(ark_mem);
}

/*------------------------------------------------------------------------------
  This routine performs a single step of the partitioned DAE method.
  ----------------------------------------------------------------------------*/
static int pdaeStep_TakeStep(ARKodeMem ark_mem, sunrealtype* dsmPtr, int* nflagPtr)
{
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  *nflagPtr = ARK_SUCCESS;
  *dsmPtr   = ZERO; /* No error estimate */

  //TODO(SBR): extend to more stages

  for (int i = 0; i < step_mem->partitions; i++)
  {
    void* ida_mem    = step_mem->ida_mems[i];
    sunrealtype tout = ark_mem->tn + ark_mem->h;
    IDASetStopTime(ida_mem, tout);
    retval = IDASolve(ida_mem, tout, &tout, step_mem->user_datas[i].y,
                      step_mem->user_datas[i].yp, IDA_NORMAL);

    if (retval < IDA_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_INNERSTEP_FAILED, ark_mem->tcur);
      return ARK_INNERSTEP_FAIL;
    }
  }

  ark_mem->tcur = ark_mem->tn + ark_mem->h;

  *nflagPtr = pdaeStep_Nls(ark_mem, *nflagPtr);
  if (*nflagPtr < 0) { return *nflagPtr; }
  if (*nflagPtr > 0) { return TRY_AGAIN; }

  retval = pdaeStep_ReInitPartitions(ark_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  Prints integrator statistics
  ----------------------------------------------------------------------------*/
static int pdaeStep_PrintAllStats(ARKodeMem ark_mem,
                                  SUNDIALS_MAYBE_UNUSED FILE* outfile,
                                  SUNDIALS_MAYBE_UNUSED SUNOutputFormat fmt)
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
  if (step_mem != NULL)
  {
    N_VDestroy(step_mem->alg_linsol_x);
    N_VDestroy(step_mem->alg_tmp);
    N_VDestroy(step_mem->alg_fcur);
    N_VDestroy(step_mem->alg_res);
    N_VDestroy(step_mem->alg_cor);
    N_VDestroy(step_mem->alg_pred);
    N_VDestroy(step_mem->alg_ewt);
    N_VDestroy(step_mem->alg);
    N_VDestroy(step_mem->yp);
    if (step_mem->ida_mems != NULL)
    {
      for (int i = 0; i < step_mem->partitions; i++)
      {
        if (step_mem->ida_mems[i] != NULL) { IDAFree(&step_mem->ida_mems[i]); }
      }
    }
    if (step_mem->user_datas != NULL)
    {
      for (int i = 0; i < step_mem->partitions; i++)
      {
        N_VDestroy(step_mem->user_datas[i].y);
        N_VDestroy(step_mem->user_datas[i].yp);
      }
    }
    free(step_mem->ida_mems);
    free(step_mem->user_datas);
    if (step_mem->ownNLS) { SUNNonlinSolFree(step_mem->NLS); }
    free(step_mem);
  }
  ark_mem->step_mem = NULL;
}

/*------------------------------------------------------------------------------
  This routine outputs the memory from the PDAEStep structure to a specified
  file pointer (useful when debugging).
  ----------------------------------------------------------------------------*/
static void pdaeStep_PrintMem(ARKodeMem ark_mem,
                              SUNDIALS_MAYBE_UNUSED FILE* outfile)
{
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return; }

  // TODO(SBR): Print stats and other props
}

static int pdaeStep_Residual(sunrealtype t, N_Vector y, N_Vector yp,
                             N_Vector res, void* user_data)
{
  IDAUserData* ida_user_data = (IDAUserData*)user_data;
  void* arkode_user_data     = ((ARKodeMem)ida_user_data->ark_mem)->user_data;

  // TODO(SBR): generalize w to polynomial. This is using w=w_n for now
  N_Vector w =
    PDAEStepGetCouplingSubvector(((ARKodeMem)ida_user_data->ark_mem)->yn);

  return ida_user_data->component_res_fn(ida_user_data->partition, t, y, w, yp,
                                         res, arkode_user_data);
}

static int pdaeStep_Jac(sunrealtype t, sunrealtype c_j, N_Vector y, N_Vector yp,
                        N_Vector r, SUNMatrix Jac, void* user_data,
                        SUNDIALS_MAYBE_UNUSED N_Vector tmp1,
                        SUNDIALS_MAYBE_UNUSED N_Vector tmp2,
                        SUNDIALS_MAYBE_UNUSED N_Vector tmp3)
{
  IDAUserData* ida_user_data = (IDAUserData*)user_data;
  void* arkode_user_data     = ((ARKodeMem)ida_user_data->ark_mem)->user_data;
  // TODO(SBR): generalize w to polynomial. This is using w=w_n for now
  N_Vector w =
    PDAEStepGetCouplingSubvector(((ARKodeMem)ida_user_data->ark_mem)->yn);

  return ida_user_data->component_res_jac(ida_user_data->partition, t, c_j, y,
                                          w, yp, r, Jac, arkode_user_data);
}

static int pdaeStep_JacTimes(sunrealtype t, N_Vector y, N_Vector yp, N_Vector r,
                             N_Vector v, N_Vector Jv, sunrealtype c_j,
                             void* user_data, N_Vector tmp1, N_Vector tmp2)
{
  IDAUserData* ida_user_data = (IDAUserData*)user_data;
  void* arkode_user_data     = ((ARKodeMem)ida_user_data->ark_mem)->user_data;
  // TODO(SBR): generalize w to polynomial. This is using w=w_n for now
  N_Vector w =
    PDAEStepGetCouplingSubvector(((ARKodeMem)ida_user_data->ark_mem)->yn);

  return ida_user_data->component_res_jtimes(ida_user_data->partition, t, y, w,
                                             yp, r, v, Jv, c_j,
                                             arkode_user_data, tmp1, tmp2);
}

static int pdaeStep_ReInitPartitions(ARKodeMem ark_mem)
{
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  for (int i = 0; i < step_mem->partitions; i++)
  {
    // printf("-----------------------\ni = %d\nstep = %ld\n", i, ark_mem->nst);
    // IDAPrintAllStats(step_mem->ida_mems[i], stdout, SUN_OUTPUTFORMAT_TABLE);
    retval = IDAReInit(step_mem->ida_mems[i], ark_mem->tcur,
                       step_mem->user_datas[i].y, step_mem->user_datas[i].yp);
    if (retval != IDA_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to reinitialize IDA for partition %i", i);
      return ARK_INNERSTEP_FAIL;
    }
  }

  return ARK_SUCCESS;
}

/*------------------------------------------------------------------------------
  Creates the PDAEStep integrator
  ----------------------------------------------------------------------------*/
void* PDAEStepCreate(PDAEStepComponentResFn component_res_fn,
                     PDAEStepAlgebraicResFn algebraic_res_fn, sunrealtype t0,
                     N_Vector y0, N_Vector yp0, int partitions, SUNContext sunctx)
{
  if (component_res_fn == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "component_res_fn must be non-NULL");
    return NULL;
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
  memset(step_mem, 0, sizeof(*step_mem));

  // ark_mem->step_attachlinsol              = pdaeStep_AttachLinsol;
  ark_mem->step_init = pdaeStep_Init;
  ark_mem->step      = pdaeStep_TakeStep;
  // ark_mem->step_setuserdata               = pdaeStep_SetUserData;
  ark_mem->step_printallstats = pdaeStep_PrintAllStats;
  // ark_mem->step_writeparameters           = pdaeStep_WriteParameters;
  ark_mem->step_setusecompensatedsums = NULL;
  ark_mem->step_free                  = pdaeStep_Free;
  ark_mem->step_printmem              = pdaeStep_PrintMem;
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
  ark_mem->step_supports_adaptive = SUNFALSE;
  /* The existing implicit infrastructure does not map very well to the
   * algebraic solve for PDAE methods. The dimension of the nonlinear system is
   * smaller than the overall problem, and there is effectively a 0 mass matrix.
   * Therefore, we set this to false and manually handle the algebraic solve.
   */
  ark_mem->step_supports_implicit   = SUNFALSE;
  ark_mem->step_supports_massmatrix = SUNFALSE;
  ark_mem->step_supports_relaxation = SUNFALSE;
  ark_mem->step_mem                 = (void*)step_mem;

  step_mem->algebraic_res_fn     = algebraic_res_fn;
  step_mem->algebraic_res_jac    = NULL;
  step_mem->algebraic_res_jtimes = NULL;
  step_mem->yp                   = N_VClone(yp0);
  if (step_mem->yp == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Unable to clone yp vector");
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }
  step_mem->partitions = partitions;
  step_mem->order      = DEFAULT_ORDER;
  step_mem->maxcor     = DEFAULT_MAX_COR;
  step_mem->nlscoef    = DEFAULT_NLSCOEF;
  step_mem->crdown     = DEFAULT_CRDOWN;
  step_mem->rdiv       = DEFAULT_RDIV;
  step_mem->crate      = ONE;
  step_mem->jcur       = SUNFALSE;
  step_mem->ownNLS     = SUNFALSE;

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

  retval = pdaeStep_AllocAlgebraicVectors(ark_mem, step_mem, sunctx);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to allocate PDAEStep algebraic vectors");
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }

  step_mem->NLS = SUNNonlinSol_Newton(step_mem->alg, sunctx);
  if (step_mem->NLS == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Unable to allocate default PDAEStep nonlinear solver");
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }

  retval = PDAEStepSetNonlinearSolver((void*)ark_mem, step_mem->NLS);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Unable to attach default PDAEStep nonlinear solver");
    SUNNonlinSolFree(step_mem->NLS);
    step_mem->NLS = NULL;
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }
  step_mem->ownNLS = SUNTRUE;

  step_mem->ida_mems = calloc(partitions, sizeof(*step_mem->ida_mems));
  if (step_mem->ida_mems == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Unable to allocate IDA array");
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }
  step_mem->user_datas = calloc(partitions, sizeof(*step_mem->user_datas));
  if (step_mem->user_datas == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Unable to allocate user data array");
    ARKodeFree((void**)&ark_mem);
    return NULL;
  }

  for (int i = 0; i < partitions; i++)
  {
    step_mem->ida_mems[i] = IDACreate(sunctx);
    if (step_mem->ida_mems[i] == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to allocate IDA for partition %i", i);
      ARKodeFree((void**)&ark_mem);
      return NULL;
    }

    {
      N_Vector yi_sub_vecs[] = {PDAEStepGetDifferentialSubvector(ark_mem->ycur, i),
                                PDAEStepGetAlgebraicSubvector(ark_mem->ycur, i)};
      N_Vector yi = N_VNew_ManyVector(2, yi_sub_vecs, sunctx);
      if (yi == NULL)
      {
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
      if (ypi == NULL)
      {
        arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                        "Unable to create yp N_Vector for partition %i", i);
        ARKodeFree((void**)&ark_mem);
        return NULL;
      }
      step_mem->user_datas[i].yp = ypi;
    }

    step_mem->user_datas[i].ark_mem              = ark_mem;
    step_mem->user_datas[i].partition            = i;
    step_mem->user_datas[i].component_res_fn     = component_res_fn;
    step_mem->user_datas[i].component_res_jac    = NULL;
    step_mem->user_datas[i].component_res_jtimes = NULL;

    retval = IDAInit(step_mem->ida_mems[i], pdaeStep_Residual, t0,
                     step_mem->user_datas[i].y, step_mem->user_datas[i].yp);
    if (retval != IDA_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_ATTACH_ERR, __LINE__, __func__,
                      __FILE__, "Unable to initialize IDA for partition %i", i);
      ARKodeFree((void**)&ark_mem);
      return NULL;
    }

    retval = IDASetUserData(step_mem->ida_mems[i], &step_mem->user_datas[i]);
    if (retval != IDA_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_ATTACH_ERR, __LINE__, __func__,
                      __FILE__, "Unable to set IDA user data for partition %i",
                      i);
      ARKodeFree((void**)&ark_mem);
      return NULL;
    }

    retval = IDASetMaxNumSteps(step_mem->ida_mems[i], -1);
    if (retval != IDA_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_ATTACH_ERR, __LINE__, __func__,
                      __FILE__, "Unable to set IDA max steps for partition %i",
                      i);
      ARKodeFree((void**)&ark_mem);
      return NULL;
    }

    retval = IDASStolerances(step_mem->ida_mems[i], ark_mem->reltol,
                             ark_mem->Sabstol);
    if (retval != IDA_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_INNERSTEP_ATTACH_ERR, __LINE__, __func__,
                      __FILE__, "Unable to set IDA tolerances for partition %i",
                      i);
      ARKodeFree((void**)&ark_mem);
      return NULL;
    }
  }

  return ark_mem;
}

int PDAEStepGetNumPartitions(void* arkode_mem, int* partitions)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (partitions == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "partitions is NULL");
    return ARK_ILL_INPUT;
  }

  *partitions = step_mem->partitions;

  return ARK_SUCCESS;
}

int PDAEStepGetPartitionVectorTemplate(void* arkode_mem, int partition,
                                       N_Vector* y)
{
  ARKodeMem ark_mem          = NULL;
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

  if (y == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "y is NULL");
    return ARK_ILL_INPUT;
  }

  *y = step_mem->user_datas[partition].y;

  return ARK_SUCCESS;
}

int PDAEStepGetPartitionIntegrator(void* arkode_mem, int partition, void** ida_mem)
{
  ARKodeMem ark_mem          = NULL;
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

  if (ida_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "ida_mem is NULL");
    return ARK_ILL_INPUT;
  }

  *ida_mem = step_mem->ida_mems[partition];

  return ARK_SUCCESS;
}

int PDAEStepSetPartitionJacobian(void* arkode_mem, PDAEStepLsComponentJacFn jac)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  for (int partition = 0; partition < step_mem->partitions; partition++)
  {
    step_mem->user_datas[partition].component_res_jac = jac;
  }

  return ARK_SUCCESS;
}

int PDAEStepSetPartitionJacTimes(void* arkode_mem,
                                 PDAEStepLsComponentJacTimesVecFn jtimes)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  for (int partition = 0; partition < step_mem->partitions; partition++)
  {
    step_mem->user_datas[partition].component_res_jtimes = jtimes;
  }

  return ARK_SUCCESS;
}

int PDAEStepSetCouplingJacobian(void* arkode_mem, PDAEStepLsAlgebraicJacFn jac)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  step_mem->algebraic_res_jac = jac;

  return ARK_SUCCESS;
}

int PDAEStepGetAlgebraicVectorTemplate(void* arkode_mem, N_Vector* y)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (y == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "y is NULL");
    return ARK_ILL_INPUT;
  }

  *y = step_mem->alg;

  return ARK_SUCCESS;
}

N_Vector PDAEStepManyVector(N_Vector* x, N_Vector* z, N_Vector w, int partitions)
{
  N_Vector* sub_vecs = malloc((2 * partitions + 1) * sizeof(*sub_vecs));
  if (sub_vecs == NULL) { return NULL; }
  for (int i = 0; i < partitions; i++)
  {
    sub_vecs[i]              = x[i];
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
  return N_VGetSubvector_ManyVector(y, partition);
}

N_Vector PDAEStepGetAlgebraicSubvector(N_Vector y, int partition)
{
  sunindextype num_subvectors = N_VGetNumSubvectors_ManyVector(y);
  if (num_subvectors == 2) { return N_VGetSubvector_ManyVector(y, 1); }

  sunindextype partitions = (num_subvectors - 1) / 2;
  return N_VGetSubvector_ManyVector(y, partitions + partition);
}
