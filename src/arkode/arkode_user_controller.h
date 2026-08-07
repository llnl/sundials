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
 * This is the header file for a user-provided controller function
 * SUNAdaptController module implementation.  This provides backwards-
 * compatibility for ARKODE's previous "ARKAdaptFn"
 * -----------------------------------------------------------------*/

#ifndef _ARK_USERCONTROL_H
#define _ARK_USERCONTROL_H

#include <stdio.h>
#include <sundials/sundials_adaptcontroller.h>

#include "arkode_impl.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ---------------------------------------------------
 * ARKUserControl implementation of SUNAdaptController
 * --------------------------------------------------- */

struct _ARKUserControlContent
{
  sunrealtype hp;    /* h from previous step */
  sunrealtype hpp;   /* h from 2 steps ago */
  sunrealtype ep;    /* error from previous step */
  sunrealtype epp;   /* error from 2 steps ago */
  ARKodeMem ark_mem; /* main ARKODE memory structure */
  ARKAdaptFn hadapt; /* user-provided adaptivity fn */
  void* hadapt_data; /* user-provided data pointer */
};

typedef struct _ARKUserControlContent* ARKUserControlContent;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController ARKUserControl(SUNContext sunctx, void* arkode_mem,
                                  ARKAdaptFn hadapt, void* hadapt_data);
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptController_GetType_ARKUserControl(
  SUNAdaptController C);
SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_EstimateStep_ARKUserControl(SUNAdaptController C,
                                                          sunrealtype h, int p,
                                                          sunrealtype dsm,
                                                          sunrealtype* hnew);
SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_Reset_ARKUserControl(SUNAdaptController C);
SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_Write_ARKUserControl(SUNAdaptController C,
                                                   FILE* fptr);
SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_UpdateH_ARKUserControl(SUNAdaptController C,
                                                     sunrealtype h,
                                                     sunrealtype dsm);
SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_Space_ARKUserControl(SUNAdaptController C,
                                                   long int* lenrw,
                                                   long int* leniw);

#ifdef __cplusplus
}
#endif

#endif /* _ARK_USERCONTROL_H */
