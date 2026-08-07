/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
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
 * This is the implementation file for the ARKUserControl
 * SUNAdaptController module.
 * -----------------------------------------------------------------*/

#include "arkode_user_controller.h"

#include <stdio.h>
#include <stdlib.h>

/* ---------------
 * Macro accessors
 * --------------- */

#define SC_CONTENT(C) ((ARKUserControlContent)(C->content))
#define SC_HP(C)      (SC_CONTENT(C)->hp)
#define SC_HPP(C)     (SC_CONTENT(C)->hpp)
#define SC_EP(C)      (SC_CONTENT(C)->ep)
#define SC_EPP(C)     (SC_CONTENT(C)->epp)
#define SC_ARKMEM(C)  (SC_CONTENT(C)->ark_mem)
#define SC_HADAPT(C)  (SC_CONTENT(C)->hadapt)
#define SC_DATA(C)    (SC_CONTENT(C)->hadapt_data)

/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new ARKUserControl controller
 * ----------------------------------------------------------------- */

SUNAdaptController ARKUserControl(SUNContext sunctx, void* arkode_mem,
                                  ARKAdaptFn hadapt, void* hadapt_data)
{
  SUNAdaptController C;
  ARKUserControlContent content;

  /* Return with failure if hadapt, arkode_mem, or context are NULL */
  if ((hadapt == NULL) || (arkode_mem == NULL) || (sunctx == NULL))
  {
    return (NULL);
  }

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->gettype      = SUNAdaptController_GetType_ARKUserControl;
  C->ops->estimatestep = SUNAdaptController_EstimateStep_ARKUserControl;
  C->ops->reset        = SUNAdaptController_Reset_ARKUserControl;
  C->ops->write        = SUNAdaptController_Write_ARKUserControl;
  C->ops->updateh      = SUNAdaptController_UpdateH_ARKUserControl;
  C->ops->space        = SUNAdaptController_Space_ARKUserControl;

  /* Create content */
  content = NULL;
  content = (ARKUserControlContent)malloc(sizeof *content);
  if (content == NULL)
  {
    (void)SUNAdaptController_Destroy(C);
    return (NULL);
  }

  /* Attach content */
  C->content = content;

  /* Attach ARKODE memory structure */
  content->ark_mem = (ARKodeMem)arkode_mem;

  /* Attach user-provided adaptivity function and data */
  content->hadapt      = hadapt;
  content->hadapt_data = hadapt_data;

  /* Fill content with default/reset values */
  SUNAdaptController_Reset_ARKUserControl(C);

  return (C);
}

/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptController_GetType_ARKUserControl(
  SUNDIALS_MAYBE_UNUSED SUNAdaptController C)
{
  return SUN_ADAPTCONTROLLER_H;
}

SUNErrCode SUNAdaptController_EstimateStep_ARKUserControl(
  SUNAdaptController C, sunrealtype h, SUNDIALS_MAYBE_UNUSED int p,
  sunrealtype dsm, sunrealtype* hnew)
{
  /* call user-provided function to compute new step */
  sunrealtype ttmp = (dsm <= ONE) ? SC_ARKMEM(C)->tn + SC_ARKMEM(C)->h
                                  : SC_ARKMEM(C)->tn;
  int retval = SC_HADAPT(C)(SC_ARKMEM(C)->ycur, ttmp, h, SC_HP(C), SC_HPP(C),
                            dsm, SC_EP(C), SC_EPP(C), SC_ARKMEM(C)->hadapt_mem->q,
                            SC_ARKMEM(C)->hadapt_mem->p, hnew, SC_DATA(C));
  if (retval != SUN_SUCCESS) { return (SUN_ERR_USER_FCN_FAIL); }
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Reset_ARKUserControl(SUNAdaptController C)
{
  SC_EP(C)  = SUN_RCONST(1.0);
  SC_EPP(C) = SUN_RCONST(1.0);
  SC_HP(C)  = SUN_RCONST(0.0);
  SC_HPP(C) = SUN_RCONST(0.0);
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Write_ARKUserControl(SUNAdaptController C,
                                                   FILE* fptr)
{
  fprintf(fptr, "ARKUserControl module:\n");
  fprintf(fptr, "  hp = " SUN_FORMAT_G "\n", SC_HP(C));
  fprintf(fptr, "  hpp = " SUN_FORMAT_G "\n", SC_HPP(C));
  fprintf(fptr, "  ep = " SUN_FORMAT_G "\n", SC_EP(C));
  fprintf(fptr, "  epp = " SUN_FORMAT_G "\n", SC_EPP(C));
  fprintf(fptr, "  hadapt_data = %p\n", SC_DATA(C));
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_UpdateH_ARKUserControl(SUNAdaptController C,
                                                     sunrealtype h,
                                                     sunrealtype dsm)
{
  SC_HPP(C) = SC_HP(C);
  SC_HP(C)  = h;
  SC_EPP(C) = SC_EP(C);
  SC_EP(C)  = dsm;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Space_ARKUserControl(
  SUNDIALS_MAYBE_UNUSED SUNAdaptController C, long int* lenrw, long int* leniw)
{
  *lenrw = 4;
  *leniw = 2;
  return SUN_SUCCESS;
}
