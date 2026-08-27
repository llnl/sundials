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

/* MRIStep ExtSTS constructor routine */
void* MRIStepExtSTSCreate(ARKRhsFn fd, ARKRhsFn fe, ARKRhsFn fi, sunrealtype t0,
                          N_Vector y0, SUNContext sunctx)
{
  /* Create LSRKStep integrator and configure with ExtSTS defaults */
  void* sts_mem = LSRKStepCreateSTS(fd, t0, y0, sunctx);
  if (sts_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Failed to create LSRKStep memory for ExtSTS method.");
    return NULL;
  }

  /* Increase the maximum number of internal stages */
  int retval = LSRKStepSetMaxNumStages(sts_mem, 10000);
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
    arkProcessError(NULL, retval, __LINE__, __func__,
                    __FILE__, "Failed to disable LRSKStep interpolation in ExtSTS method.");
    ARKodeFree(&sts_mem);
    return NULL;
  }

  /* Create the inner stepper and attach objects and function pointers */
  MRIStepInnerStepper inner_stepper = NULL;
  retval = MRIStepInnerStepper_Create(sunctx, &inner_stepper);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__,
                    __FILE__, "Failed to create MRIStep inner stepper for ExtSTS method.");
    ARKodeFree(&sts_mem);
    return NULL;
  }
  retval = MRIStepInnerStepper_SetContent(inner_stepper, sts_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to set MRIStep inner stepper content.");
    ARKodeFree(&sts_mem);
    MRIStepInnerStepper_Free(&inner_stepper);
    return NULL;
  }
  retval = MRIStepInnerStepper_SetEvolveFn(inner_stepper,
                                           extSTSInnerStepper_Evolve);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to set MRIStep inner stepper evolve function.");
    ARKodeFree(&sts_mem);
    MRIStepInnerStepper_Free(&inner_stepper);
    return NULL;
  }
  retval = MRIStepInnerStepper_SetFullRhsFn(inner_stepper,
                                            ark_MRIStepInnerFullRhs);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to set MRIStep inner stepper full RHS function.");
    ARKodeFree(&sts_mem);
    MRIStepInnerStepper_Free(&inner_stepper);
    return NULL;
  }
  retval = MRIStepInnerStepper_SetResetFn(inner_stepper, ark_MRIStepInnerReset);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to set MRIStep inner stepper reset function.");
    ARKodeFree(&sts_mem);
    MRIStepInnerStepper_Free(&inner_stepper);
    return NULL;
  }

  /* Create the MRIStep integrator, attaching the inner stepper for diffusion */
  void* arkode_mem = MRIStepCreate(fe, fi, t0, y0, inner_stepper, sunctx);
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Failed to create MRIStep integrator.");
    return NULL;
  }

  /* Store a pointer to the inner stepper content directly in MRIStep */
  ARKodeMem ark_mem         = NULL;
  ARKodeMRIStepMem step_mem = NULL;
  retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(NULL, retval, __LINE__, __func__, __FILE__,
                    "Failed to access stepper memory from MRIStep.");
    ARKodeFree(&arkode_mem);
    ARKodeFree(&sts_mem);
    MRIStepInnerStepper_Free(&inner_stepper);
    return NULL;
  }
  step_mem->extsts_method = SUNTRUE;

  /* Configure the MRIStep integrator with ExtSTS defaults */
  /*   Select default ExtSTS method based on provided RHS functions */
  MRIStepCoupling MRIC = NULL;
  if (fe == NULL && fi != NULL) /* Implicit ExtSTS method */
  {
    MRIC = MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_GARK_GIRALDO2);
  }
  else if (fe != NULL && fi == NULL) /* Explicit ExtSTS method */
  {
    MRIC = MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_EXP_GIRALDO2);
  }
  else /* ImEx ExtSTS method */
  {
    MRIC = MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_GARK_GIRALDO2);
  }
  if (MRIC == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__,
                    __FILE__, "Failed to create MRIStep coupling table for ExtSTS method.");
    ARKodeFree(&arkode_mem);
    ARKodeFree(&sts_mem);
    MRIStepInnerStepper_Free(&inner_stepper);
    return NULL;
  }
  retval = MRIStepSetCoupling(arkode_mem, MRIC);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to set MRIStep coupling table for ExtSTS method.");
    ARKodeFree(&arkode_mem);
    ARKodeFree(&sts_mem);
    MRIStepInnerStepper_Free(&inner_stepper);
    return NULL;
  }
  MRIStepCoupling_Free(MRIC);

  /* return with the constructed MRIStep object */
  return arkode_mem;
}

/* MRIStep ExtSTS reinitialization routine */
int MRIStepExtSTSReInit(void* arkode_mem, ARKRhsFn fd, ARKRhsFn fe, ARKRhsFn fi,
                        sunrealtype t0, N_Vector y0)
{
  /* access ARKodeMem and ARKodeMRIStepMem structures */
  ARKodeMem ark_mem         = NULL;
  ARKodeMRIStepMem step_mem = NULL;
  int retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                           &step_mem);
  if (retval) { return retval; }

  /* Reinitialize the MRIStep integrator */
  retval = MRIStepReInit(arkode_mem, fe, fi, t0, y0);
  if (retval)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to set MRIStep coupling table for ExtSTS method.");
    return retval;
  }

  /* Reinitialize the LSRKStep integrator */
  retval = LSRKStepReInitSTS(step_mem->stepper->content, fd, t0, y0);
  if (retval)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to set reinitialize LSRKStep for ExtSTS method.");
    return retval;
  }

  return retval;
}

/* Accessor routine for the inner LSRKStep solver from an ExtSTS method */
int MRIStepExtSTSGetSTS(void* arkode_mem, void** sts_mem)
{
  /* access ARKodeMem and ARKodeMRIStepMem structures */
  ARKodeMem ark_mem         = NULL;
  ARKodeMRIStepMem step_mem = NULL;
  int retval = mriStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                           &step_mem);
  if (retval) { return retval; }

  /* return pointer to stored STS integrator */
  *sts_mem = step_mem->stepper->content;
  return ARK_SUCCESS;
}

/* Inner stepper utility routines */
int extSTSInnerStepper_Evolve(MRIStepInnerStepper sts_mem, sunrealtype t0,
                              sunrealtype tout, N_Vector y)
{
  /* Get the forcing data */
  ARKodeMem ark_mem = (ARKodeMem)sts_mem->content;
  sunrealtype tshift, tscale;
  N_Vector* forcing;
  int nforcing;
  int retval;
  retval = MRIStepInnerStepper_GetForcingData(sts_mem, &tshift, &tscale,
                                              &forcing, &nforcing);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to retrieve forcing data for ExtSTS method.");
    return retval;
  }

  /* Set the inner forcing data */
  retval = ark_mem->step_setforcing(ark_mem, tshift, tscale, forcing, nforcing);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to set LSRKStep forcing for ExtSTS method.");
    return retval;
  }

  /* Reset STS integrator to current state */
  retval = ARKodeReset(ark_mem, t0, y);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to reset LSRKStep for ExtSTS method.");
    return retval;
  }

  /* Set step size to reach tout in a single step */
  retval = ARKodeSetFixedStep(ark_mem, tout - t0);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to set LSRKStep step size for ExtSTS method.");
    return retval;
  }

  /* Set stop time */
  retval = ARKodeSetStopTime(ark_mem, tout);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to set LSRKStep stop time for ExtSTS method.");
    return retval;
  }

  /* Evolve a single time step */
  sunrealtype tret;
  retval = ARKodeEvolve(ark_mem, tout, y, &tret, ARK_ONE_STEP);
  if (retval < 0)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to evolve LSRKStep for ExtSTS method.");
    return retval;
  }

  /* Disable inner forcing */
  retval = ark_mem->step_setforcing(ark_mem, ZERO, ONE, NULL, 0);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Failed to reset LSRKStep forcing for ExtSTS method.");
    return retval;
  }

  return ARK_SUCCESS;
}

int extSTSInnerStepper_Free(MRIStepInnerStepper* sts_mem)
{
  /* If the input is NULL, do nothing */
  if (sts_mem == NULL || *sts_mem == NULL) { return ARK_SUCCESS; }

  /* Free the LSRKStep memory */
  ARKodeFree(&((*sts_mem)->content));

  /* Set the input pointer to NULL */
  *sts_mem = NULL;
  return ARK_SUCCESS;
}

/*===============================================================
  EOF
  ===============================================================*/
