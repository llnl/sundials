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
 * This is the header file for the SUNAdaptController_MRIHTol module.
 * -----------------------------------------------------------------*/

#ifndef _SUNADAPTCONTROLLER_MRIHTOL_H
#define _SUNADAPTCONTROLLER_MRIHTOL_H

#include <stdio.h>
#include <sundials/sundials_adaptcontroller.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ----------------------------------------------------
 * MRI H+tolerance implementation of SUNAdaptController
 * ---------------------------------------------------- */

struct SUNAdaptControllerContent_MRIHTol_
{
  SUNAdaptController HControl;
  SUNAdaptController TolControl;
  sunrealtype inner_max_relch;
  sunrealtype inner_min_tolfac;
  sunrealtype inner_max_tolfac;
};

typedef struct SUNAdaptControllerContent_MRIHTol_* SUNAdaptControllerContent_MRIHTol;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_MRIHTol(SUNAdaptController HControl,
                                              SUNAdaptController TolControl,
                                              SUNContext sunctx);

SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_SetParams_MRIHTol(SUNAdaptController C,
                                                sunrealtype inner_max_relch,
                                                sunrealtype inner_min_tolfac,
                                                sunrealtype inner_max_tolfac);

SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_GetSlowController_MRIHTol(
  SUNAdaptController C,
  SUNAdaptController* Cslow); // nb::rv_policy::reference

SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_GetFastController_MRIHTol(
  SUNAdaptController C,
  SUNAdaptController* Cfast); // nb::rv_policy::reference

SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptController_GetType_MRIHTol(SUNAdaptController C);

SUNDIALS_EXPORT
int SUNAdaptController_EstimateStepTol_MRIHTol(
  SUNAdaptController C, sunrealtype H, sunrealtype tolfac, int P,
  sunrealtype DSM, sunrealtype dsm, sunrealtype* Hnew, sunrealtype* tolfacnew);

SUNDIALS_EXPORT
int SUNAdaptController_Reset_MRIHTol(SUNAdaptController C);

SUNDIALS_EXPORT
int SUNAdaptController_SetDefaults_MRIHTol(SUNAdaptController C);

SUNDIALS_EXPORT
int SUNAdaptController_Write_MRIHTol(SUNAdaptController C, FILE* fptr);

SUNDIALS_EXPORT
int SUNAdaptController_SetErrorBias_MRIHTol(SUNAdaptController C,
                                            sunrealtype bias);

SUNDIALS_EXPORT
int SUNAdaptController_UpdateMRIHTol_MRIHTol(SUNAdaptController C,
                                             sunrealtype H, sunrealtype tolfac,
                                             sunrealtype DSM, sunrealtype dsm);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
int SUNAdaptController_Space_MRIHTol(SUNAdaptController C, long int* lenrw,
                                     long int* leniw);

#ifdef __cplusplus
}
#endif

#endif /* _SUNADAPTCONTROLLER_MRIHTOL_H */
