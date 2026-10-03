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
 * Temporary header file for deprecated functions.
 * TODO:(SBR) remove in version 8.1.0
 * -----------------------------------------------------------------*/

#ifndef SUNDIALS_SUNADAPTCONTROLLER_SUNADAPTCONTROLLER_MRIHTOL_DEPRECATED_H
#define SUNDIALS_SUNADAPTCONTROLLER_SUNADAPTCONTROLLER_MRIHTOL_DEPRECATED_H

#include <sunadaptcontroller/sunadaptcontroller_mrihtol.h>

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNAdaptController_GetType instead; will be removed in version 8.1.0")
SUNAdaptController_Type SUNAdaptController_GetType_MRIHTol(SUNAdaptController C);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNAdaptController_EstimateStepTol "
                               "instead; will be removed in version 8.1.0")
int SUNAdaptController_EstimateStepTol_MRIHTol(
  SUNAdaptController C, sunrealtype H, sunrealtype tolfac, int P,
  sunrealtype DSM, sunrealtype dsm, sunrealtype* Hnew, sunrealtype* tolfacnew);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNAdaptController_Reset instead; will be removed in version 8.1.0")
int SUNAdaptController_Reset_MRIHTol(SUNAdaptController C);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNAdaptController_SetDefaults instead; "
                               "will be removed in version 8.1.0")
int SUNAdaptController_SetDefaults_MRIHTol(SUNAdaptController C);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNAdaptController_Write instead; will be removed in version 8.1.0")
int SUNAdaptController_Write_MRIHTol(SUNAdaptController C, FILE* fptr);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNAdaptController_SetErrorBias instead; "
                               "will be removed in version 8.1.0")
int SUNAdaptController_SetErrorBias_MRIHTol(SUNAdaptController C,
                                            sunrealtype bias);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNAdaptController_UpdateMRIHTol instead; "
                               "will be removed in version 8.1.0")
int SUNAdaptController_UpdateMRIHTol_MRIHTol(SUNAdaptController C,
                                             sunrealtype H, sunrealtype tolfac,
                                             sunrealtype DSM, sunrealtype dsm);

#ifdef __cplusplus
}
#endif

#endif
