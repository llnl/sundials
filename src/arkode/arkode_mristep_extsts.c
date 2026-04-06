/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 * This is the implementation file for MRIStep's ExtSTS utility routines.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <arkode/arkode_lsrkstep.h>
#include "arkode_impl.h"
#include "arkode_mristep_impl.h"

/* TODO:
   * update LSRKStep to internally support external polynomial forcing terms
   * update MRIStep's SetUserData function so that if an inner STS method is used,
     it will also set the user data for the STS integrator
   * add desired ExtSTS MRI tables to arkode_mri_tables.def
*/

/* MRIStep ExtSTS constructor routine */
void* MRIStepCreateExtSTS(ARKRhsFn fd, ARKRhsFn fe, ARKRhsFn fi,
                          sunrealtype t0, N_Vector y0, SUNContext sunctx)
{

  /* Create LSRKStep integrator and configure with ExtSTS defaults */
  void* sts_mem = LSRKStepCreateSTS(fd, t0, y0, sunctx);
  if (sts_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Failed to create LSRKStep memory for ExtSTS method.");
    return NULL;
  }

  /*   RKC method */
  int retval = LSRKStepSetSTSMethod(sts_mem, ARKODE_LSRK_RKC_2);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to set RKC method for ExtSTS method.");
    ARKodeFree(&sts_mem);
    return NULL;
  }

  /*   Increase the maximum number of internal stages */
  retval = LSRKStepSetMaxNumStages(sts_mem, 10000);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to set maximum number of STS stages.");
    ARKodeFree(&sts_mem);
    return NULL;
  }

  /* Disable temporal interpolation for inner STS method */
  retval = ARKodeSetInterpolantType(sts_mem, ARK_INTERP_NONE);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to disable LRSKStep interpolation in ExtSTS method.");
    ARKodeFree(&sts_mem);
    return NULL;
  }

  /* Create the inner stepper wrapper and attach objects and function pointers */
  extSTSInnerStepper inner_content = (extSTSInnerStepper)calloc(1, sizeof(*inner_content));
  if (inner_content == NULL)
  {
    arkProcessError(NULL, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Failed to create extSTSInnerStepper content for ExtSTS method.");
    ARKodeFree(&sts_mem);
    return NULL;
  }
  inner_content->sts_mem = sts_mem;
  inner_content->f_diffusion = fd;
  inner_content->inner_stepper = NULL;
  retval = MRIStepInnerStepper_Create(sunctx, &(inner_content->inner_stepper));
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to create MRIStep inner stepper for ExtSTS method.");
    free(inner_content);
    ARKodeFree(&sts_mem);
    return NULL;
  }
  retval = MRIStepInnerStepper_SetContent(inner_content->inner_stepper, inner_content);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to set MRIStep inner stepper content.");
    MRIStepInnerStepper_Free(&(inner_content->inner_stepper));
    free(inner_content);
    ARKodeFree(&sts_mem);
    return NULL;
  }
  retval = MRIStepInnerStepper_SetEvolveFn(inner_content->inner_stepper, extSTSInnerStepper_Evolve);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to set MRIStep inner stepper evolve function.");
    MRIStepInnerStepper_Free(&(inner_content->inner_stepper));
    free(inner_content);
    ARKodeFree(&sts_mem);
    return NULL;
  }
  retval = MRIStepInnerStepper_SetFullRhsFn(inner_content->inner_stepper, extSTSInnerStepper_FullRhs);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to set MRIStep inner stepper full RHS function.");
    MRIStepInnerStepper_Free(&(inner_content->inner_stepper));
    free(inner_content);
    ARKodeFree(&sts_mem);
    return NULL;
  }
  retval = MRIStepInnerStepper_SetResetFn(inner_content->inner_stepper, extSTSInnerStepper_Reset);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to set MRIStep inner stepper reset function.");
    MRIStepInnerStepper_Free(&(inner_content->inner_stepper));
    free(inner_content);
    ARKodeFree(&sts_mem);
    return NULL;
  }

  /* Create the MRIStep integrator, attaching the inner stepper for diffusion */
  void* arkode_mem = MRIStepCreate(fe, fi, t0, y0, inner_content->inner_stepper, sunctx);
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Failed to create MRIStep integrator.");
    return NULL;
  }

  /* Store a pointer to the inner stepper content directly in MRIStep */
  ARKodeMem ark_mem = NULL;
  ARKodeMRIStepMem step_mem = NULL;
  retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to access stepper memory from MRIStep.");
    ARKodeFree(&arkode_mem);
    MRIStepInnerStepper_Free(&(inner_content->inner_stepper));
    free(inner_content);
    ARKodeFree(&sts_mem);
    return NULL;
  }
  step_mem->extsts_inner_stepper = inner_content;

  /* Configure the MRIStep integrator with ExtSTS defaults */
  /*   Select default ExtSTS method based on provided RHS functions */
  MRIStepCoupling MRIC = NULL;
  if (fe == NULL && fi != NULL)       /* Implicit ExtSTS method */
  {
    MRIC = MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_IRK21a);  /* UPDATE THIS */
  }
  else if (fe != NULL && fi == NULL)  /* Explicit ExtSTS method */
  {
    MRIC = MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_ERK22a);  /* UPDATE THIS */
  }
  else                                /* ImEx ExtSTS method */
  {
    MRIC = MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_SR21);  /* UPDATE THIS */
  }
  if (MRIC == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Failed to create MRIStep coupling table for ExtSTS method.");
    ARKodeFree(&arkode_mem);
    MRIStepInnerStepper_Free(&(inner_content->inner_stepper));
    free(inner_content);
    ARKodeFree(&sts_mem);
    return NULL;
  }
  retval = MRIStepSetCoupling(arkode_mem, MRIC);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to set MRIStep coupling table for ExtSTS method.");
    ARKodeFree(&arkode_mem);
    MRIStepInnerStepper_Free(&(inner_content->inner_stepper));
    free(inner_content);
    ARKodeFree(&sts_mem);
    return NULL;
  }

  /* return with the constructed MRIStep object */
  return arkode_mem;
}

/* MRIStep ExtSTS reinitialization routine */
int MRIStepReInitExtSTS(void* arkode_mem, ARKRhsFn fd, ARKRhsFn fe,
                        ARKRhsFn fi, sunrealtype t0, N_Vector y0)
{
  /* access ARKodeMem structure */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ARKodeMem ark_mem = (ARKodeMem)arkode_mem;
  int retval = MRIStepReInit(arkode_mem, fe, fi, t0, y0);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to set MRIStep coupling table for ExtSTS method.");
    return retval;
  }
  void* sts_mem = MRIStep_GetSTSStepper(arkode_mem);
  if (sts_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Failed to access STS integrator from MRIStep memory.");
    return ARK_MEM_FAIL;
  }
  retval = LSRKStepReInitSTS(sts_mem, fd, t0, y0);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to set reinitialize LSRKStep for ExtSTS method.");
    return retval;
  }
  return retval;
}

/* Accessor routine for the inner LSRKStep solver from an ExtSTS method */
void* MRIStep_GetSTSStepper(void* arkode_mem)
{
  /* access ARKodeMem and ARKodeMRIStepMem structures */
  ARKodeMem ark_mem = NULL;
  ARKodeMRIStepMem step_mem = NULL;
  int retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) { return NULL; }

  /* return pointer to STS integrator */
  return (void*) step_mem->extsts_inner_stepper->sts_mem;
}

/* Inner stepper utility routines */
int extSTSInnerStepper_Evolve(MRIStepInnerStepper sts_mem, sunrealtype t0,
                              sunrealtype tout, N_Vector y)
{
  /* Reset STS integrator to current state */
  int retval = ARKodeReset(EXTSTS_STS(sts_mem), t0, y);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to reset LSRKStep for ExtSTS method.");
    return retval;
  }

  /* Set step size to reach tout in a single step */
  retval = ARKodeSetFixedStep(EXTSTS_STS(sts_mem), tout - t0);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to set LSRKStep step size for ExtSTS method.");
    return retval;
  }

  /* Set stop time */
  retval = ARKodeSetStopTime(EXTSTS_STS(sts_mem), tout);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to set LSRKStep stop time for ExtSTS method.");
    return retval;
  }

  /* Evolve a single time step */
  sunrealtype tret;
  retval = ARKodeEvolve(EXTSTS_STS(sts_mem), tout, y, &tret, ARK_ONE_STEP);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to evolve LSRKStep for ExtSTS method.");
    return retval;
  }

  return 0;
}

int extSTSInnerStepper_FullRhs(MRIStepInnerStepper sts_mem, sunrealtype t,
                               N_Vector y, N_Vector f, int mode)
{
  /* Call diffusion RHS function */
  int retval = EXTSTS_FD(sts_mem)(t, y, f, EXTSTS_UDATA(sts_mem));
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to evaluate diffusion RHS for ExtSTS method.");
    return retval;
  }

  return 0;
}

int extSTSInnerStepper_Reset(MRIStepInnerStepper sts_mem, sunrealtype tR,
                             N_Vector yR)
{
  /* Reset STS integrator to current state */
  int retval = ARKodeReset(EXTSTS_STS(sts_mem), tR, yR);
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to reset LSRKStep for ExtSTS method.");
  }
  return 0;
}

int extSTSInnerStepper_Free(MRIStepInnerStepper* sts_mem)
{
  /* If the input is NULL, do nothing */
  if (sts_mem == NULL || *sts_mem == NULL) { return 0; }

  /* Access the inner stepper content */
  extSTSInnerStepper inner_content = (extSTSInnerStepper)(*sts_mem)->content;

  /* Free the LSRKStep memory */
  ARKodeFree(&(inner_content->sts_mem));

  /* Free the inner stepper content */
  free(inner_content);

  /* Set the input pointer to NULL */
  *sts_mem = NULL;
  return 0;
}

/*===============================================================
  EOF
  ===============================================================*/
