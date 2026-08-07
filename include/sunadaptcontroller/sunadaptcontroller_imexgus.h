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
 * This is the header file for the SUNAdaptController_ImExGus module.
 * -----------------------------------------------------------------*/

#ifndef _SUNADAPTCONTROLLER_IMEXGUS_H
#define _SUNADAPTCONTROLLER_IMEXGUS_H

#include <stdio.h>
#include <sundials/sundials_adaptcontroller.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ----------------------------------------------------
 * ImEx Gustafsson implementation of SUNAdaptController
 * ---------------------------------------------------- */

struct _SUNAdaptControllerContent_ImExGus
{
  sunrealtype k1i; /* internal controller parameters */
  sunrealtype k2i;
  sunrealtype k1e;
  sunrealtype k2e;
  sunrealtype bias; /* error bias factor */
  sunrealtype ep;   /* error from previous step */
  sunrealtype hp;   /* previous step size */
  /* TODO(SRB): Consider removing firststep. We can use a negative ep to
   * indicate we are on the first step */
  sunbooleantype firststep; /* flag indicating first step */
};

typedef struct _SUNAdaptControllerContent_ImExGus* SUNAdaptControllerContent_ImExGus;

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_ImExGus(SUNContext sunctx);

SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_SetParams_ImExGus(SUNAdaptController C,
                                                sunrealtype k1e, sunrealtype k2e,
                                                sunrealtype k1i, sunrealtype k2i);

SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptController_GetType_ImExGus(SUNAdaptController C);

SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_EstimateStep_ImExGus(SUNAdaptController C,
                                                   sunrealtype h, int p,
                                                   sunrealtype dsm,
                                                   sunrealtype* hnew);

SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_Reset_ImExGus(SUNAdaptController C);

SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_SetDefaults_ImExGus(SUNAdaptController C);

SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_Write_ImExGus(SUNAdaptController C, FILE* fptr);

SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_SetErrorBias_ImExGus(SUNAdaptController C,
                                                   sunrealtype bias);
SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_UpdateH_ImExGus(SUNAdaptController C,
                                              sunrealtype h, sunrealtype dsm);

SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
SUNErrCode SUNAdaptController_Space_ImExGus(SUNAdaptController C,
                                            long int* lenrw, long int* leniw);

#ifdef __cplusplus
}
#endif

#endif /* _SUNADAPTCONTROLLER_IMEXGUS_H */
