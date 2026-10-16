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

#ifndef SUNDIALS_SUNADAPTCONTROLLER_SUNADAPTCONTROLLER_SODERLIND_DEPRECATED_H
#define SUNDIALS_SUNADAPTCONTROLLER_SUNADAPTCONTROLLER_SODERLIND_DEPRECATED_H

#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>

#ifdef __cplusplus
extern "C" {
#endif

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNAdaptController_GetType instead; will be removed in version 8.1.0")
SUNAdaptController_Type SUNAdaptController_GetType_Soderlind(SUNAdaptController C);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNAdaptController_EstimateStep instead; "
                               "will be removed in version 8.1.0")
SUNErrCode SUNAdaptController_EstimateStep_Soderlind(SUNAdaptController C,
                                                     sunrealtype h, int p,
                                                     sunrealtype dsm,
                                                     sunrealtype* hnew);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNAdaptController_Reset instead; will be removed in version 8.1.0")
SUNErrCode SUNAdaptController_Reset_Soderlind(SUNAdaptController C);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNAdaptController_SetDefaults instead; "
                               "will be removed in version 8.1.0")
SUNErrCode SUNAdaptController_SetDefaults_Soderlind(SUNAdaptController C);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNAdaptController_Write instead; will be removed in version 8.1.0")
SUNErrCode SUNAdaptController_Write_Soderlind(SUNAdaptController C, FILE* fptr);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNAdaptController_SetErrorBias instead; "
                               "will be removed in version 8.1.0")
SUNErrCode SUNAdaptController_SetErrorBias_Soderlind(SUNAdaptController C,
                                                     sunrealtype bias);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use SUNAdaptController_UpdateH instead; will be removed in version 8.1.0")
SUNErrCode SUNAdaptController_UpdateH_Soderlind(SUNAdaptController C,
                                                sunrealtype h, sunrealtype dsm);

#ifdef __cplusplus
}
#endif

#endif
