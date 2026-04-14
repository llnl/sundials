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
 * This is the implementation file for MRIStep's multirate adaptivity controller
 * layer.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode_impl.h"
#include "arkode_mristep_impl.h"

/*--------------------------------------------
  MRIStep SUNAdaptController wrapper functions
  --------------------------------------------*/

SUNAdaptController SUNAdaptController_MRIStep(ARKodeMem ark_mem,
                                              SUNAdaptController CMRI)
{
  SUNAdaptController C;
  mriStepControlContent content;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* Return with failure if input controller is NULL or has
     unsupported type */
  if (CMRI == NULL) { return (NULL); }
  if (SUNAdaptController_GetType(CMRI) != SUN_ADAPTCONTROLLER_MRI_H_TOL)
  {
    return (NULL);
  }

  /* Return with failure if stepper is inaccessible */
  retval = mriStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) return (NULL);

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(CMRI->sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->gettype      = SUNAdaptController_GetType_MRIStep;
  C->ops->estimatestep = SUNAdaptController_EstimateStep_MRIStep;
  C->ops->reset        = SUNAdaptController_Reset_MRIStep;
  C->ops->write        = SUNAdaptController_Write_MRIStep;
  C->ops->updateh      = SUNAdaptController_UpdateH_MRIStep;
  C->ops->space        = SUNAdaptController_Space_MRIStep;

  /* Create content */
  content = NULL;
  content = (mriStepControlContent)malloc(sizeof *content);
  if (content == NULL)
  {
    SUNAdaptController_Destroy(C);
    return (NULL);
  }

  /* Attach ARKODE memory, MRI stepper memory and MRI controller objects */
  content->ark_mem  = ark_mem;
  content->step_mem = step_mem;
  content->C        = CMRI;

  /* Attach content and return */
  C->content = content;
  return (C);
}

SUNErrCode SUNAdaptController_EstimateStep_MRIStep(SUNAdaptController C,
                                                   sunrealtype H, int P,
                                                   sunrealtype DSM,
                                                   sunrealtype* Hnew)
{
  /* Shortcuts to ARKODE and MRIStep memory */
  ARKodeMem ark_mem         = MRICONTROL_A(C);
  ARKodeMRIStepMem step_mem = MRICONTROL_S(C);
  if ((ark_mem == NULL) || (step_mem == NULL)) { return SUN_ERR_MEM_FAIL; }

  /* Estimate slow stepsize from MRI controller */
  return SUNAdaptController_EstimateStepTol(MRICONTROL_C(C), H,
                                            step_mem->inner_rtol_factor, P, DSM,
                                            step_mem->inner_dsm, Hnew,
                                            &(step_mem->inner_rtol_factor_new));
}

SUNErrCode SUNAdaptController_UpdateH_MRIStep(SUNAdaptController C,
                                              sunrealtype H, sunrealtype DSM)
{
  /* Shortcuts to ARKODE and MRIStep memory */
  ARKodeMem ark_mem         = MRICONTROL_A(C);
  ARKodeMRIStepMem step_mem = MRICONTROL_S(C);
  if ((ark_mem == NULL) || (step_mem == NULL)) { return SUN_ERR_MEM_FAIL; }

  /* Update MRI controller */
  SUNErrCode retval = SUNAdaptController_UpdateMRIHTol(MRICONTROL_C(C), H,
                                                       step_mem->inner_rtol_factor,
                                                       DSM, step_mem->inner_dsm);
  if (retval != SUN_SUCCESS) { return (retval); }

  /* Update inner controller parameter to most-recent prediction */
  step_mem->inner_rtol_factor = step_mem->inner_rtol_factor_new;

  /* return with success*/
  return SUN_SUCCESS;
}

SUNAdaptController_Type SUNAdaptController_GetType_MRIStep(SUNAdaptController C)
{
  return SUNAdaptController_GetType(MRICONTROL_C(C));
}

SUNErrCode SUNAdaptController_Reset_MRIStep(SUNAdaptController C)
{
  return SUNAdaptController_Reset(MRICONTROL_C(C));
}

SUNErrCode SUNAdaptController_Write_MRIStep(SUNAdaptController C, FILE* fptr)
{
  return SUNAdaptController_Write(MRICONTROL_C(C), fptr);
}

SUNErrCode SUNAdaptController_Space_MRIStep(SUNAdaptController C,
                                            long int* lenrw, long int* leniw)
{
  return SUNAdaptController_Space(MRICONTROL_C(C), lenrw, leniw);
}

/*===============================================================
  EOF
  ===============================================================*/
