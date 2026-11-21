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
 * This is the implementation file for the SUNAdaptController_ImExGus
 * module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sunadaptcontroller/sunadaptcontroller_imexgus_deprecated.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_errors.h>

#include "sundials_cli.h"
#include "sundials_macros.h"

/* ---------------
 * Macro accessors
 * --------------- */

#define SACIMEXGUS_CONTENT(C)   ((SUNAdaptControllerContent_ImExGus)(C->content))
#define SACIMEXGUS_K1E(C)       (SACIMEXGUS_CONTENT(C)->k1e)
#define SACIMEXGUS_K2E(C)       (SACIMEXGUS_CONTENT(C)->k2e)
#define SACIMEXGUS_K1I(C)       (SACIMEXGUS_CONTENT(C)->k1i)
#define SACIMEXGUS_K2I(C)       (SACIMEXGUS_CONTENT(C)->k2i)
#define SACIMEXGUS_BIAS(C)      (SACIMEXGUS_CONTENT(C)->bias)
#define SACIMEXGUS_EP(C)        (SACIMEXGUS_CONTENT(C)->ep)
#define SACIMEXGUS_HP(C)        (SACIMEXGUS_CONTENT(C)->hp)
#define SACIMEXGUS_FIRSTSTEP(C) (SACIMEXGUS_CONTENT(C)->firststep)

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K1E  SUN_RCONST(0.367)
#define DEFAULT_K2E  SUN_RCONST(0.268)
#define DEFAULT_K1I  SUN_RCONST(0.95)
#define DEFAULT_K2I  SUN_RCONST(0.95)
#define DEFAULT_BIAS SUN_RCONST(1.0)

/*
 * ----------------------------------------------------------------------------
 * Un-exported implementation specific routines
 * ----------------------------------------------------------------------------
 */

static SUNErrCode sunAdaptControllerEstimateStep_ImExGus(SUNAdaptController C,
                                                         sunrealtype h, int p,
                                                         sunrealtype dsm,
                                                         sunrealtype* hnew);

static SUNAdaptController_Type sunAdaptControllerGetType_ImExGus(
  SUNAdaptController C);

static SUNErrCode sunAdaptControllerReset_ImExGus(SUNAdaptController C);

static SUNErrCode sunAdaptControllerSetDefaults_ImExGus(SUNAdaptController C);

static SUNErrCode sunAdaptControllerSetErrorBias_ImExGus(SUNAdaptController C,
                                                         sunrealtype bias);

static SUNErrCode sunAdaptControllerUpdateH_ImExGus(SUNAdaptController C,
                                                    sunrealtype h,
                                                    sunrealtype dsm);

static SUNErrCode sunAdaptControllerWrite_ImExGus(SUNAdaptController C,
                                                  FILE* fptr);

static SUNErrCode setFromCommandLine_ImExGus(SUNAdaptController C,
                                             const char* Cid, int argc,
                                             char* argv[]);
SUNErrCode SUNAdaptController_SetOptions_ImExGus(SUNAdaptController C,
                                                 const char* Cid,
                                                 const char* file_name,
                                                 int argc, char* argv[]);

/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new ImExGus controller
 */

SUNAdaptController SUNAdaptController_ImExGus(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C;
  SUNAdaptControllerContent_ImExGus content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  C->ops->gettype      = sunAdaptControllerGetType_ImExGus;
  C->ops->estimatestep = sunAdaptControllerEstimateStep_ImExGus;
  C->ops->reset        = sunAdaptControllerReset_ImExGus;
  C->ops->setoptions   = SUNAdaptController_SetOptions_ImExGus;
  C->ops->setdefaults  = sunAdaptControllerSetDefaults_ImExGus;
  C->ops->write        = sunAdaptControllerWrite_ImExGus;
  C->ops->seterrorbias = sunAdaptControllerSetErrorBias_ImExGus;
  C->ops->updateh      = sunAdaptControllerUpdateH_ImExGus;
  /* Create content */
  content = NULL;
  content = (SUNAdaptControllerContent_ImExGus)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach content */
  C->content = content;

  /* Fill content with default/reset values */
  SUNCheckCallNull(sunAdaptControllerSetDefaults_ImExGus(C));
  SUNCheckCallNull(sunAdaptControllerReset_ImExGus(C));

  return (C);
}

/* ----------------------------------------------------------------------------
 * Function to control set routines via the command line or file
 */

SUNErrCode SUNAdaptController_SetOptions_ImExGus(
  SUNAdaptController C, const char* Cid,
  SUNDIALS_MAYBE_UNUSED const char* file_name, int argc, char* argv[])
{
  SUNFunctionBegin(C->sunctx);

  /* File-based option control is currently unimplemented */
  SUNAssert((file_name == NULL || strlen(file_name) == 0),
            SUN_ERR_ARG_INCOMPATIBLE);

  if (argc > 0 && argv != NULL)
  {
    SUNCheckCall(setFromCommandLine_ImExGus(C, Cid, argc, argv));
  }

  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * Function to control ImExGus parameters from the command line
 */

static SUNErrCode setFromCommandLine_ImExGus(SUNAdaptController C,
                                             const char* Cid, int argc,
                                             char* argv[])
{
  SUNFunctionBegin(C->sunctx);

  /* Prefix for options to set */
  const char* default_id = "sunadaptcontroller";
  size_t offset          = strlen(default_id) + 1;
  if (Cid != NULL && strlen(Cid) > 0) { offset = strlen(Cid) + 1; }
  char* prefix = (char*)malloc(sizeof(char) * (offset + 1));
  if (Cid != NULL && strlen(Cid) > 0) { strcpy(prefix, Cid); }
  else { strcpy(prefix, default_id); }
  strcat(prefix, ".");

  int retval;
  sunbooleantype write_parameters = SUNFALSE;
  for (int idx = 1; idx < argc; idx++)
  {
    /* skip command-line arguments that do not begin with correct prefix */
    if (strncmp(argv[idx], prefix, strlen(prefix)) != 0) { continue; }

    /* control over SetParams function */
    if (strcmp(argv[idx] + offset, "params_imexgus") == 0)
    {
      idx += 1;
      sunrealtype rarg1 = SUNStrToReal(argv[idx]);
      idx += 1;
      sunrealtype rarg2 = SUNStrToReal(argv[idx]);
      idx += 1;
      sunrealtype rarg3 = SUNStrToReal(argv[idx]);
      idx += 1;
      sunrealtype rarg4 = SUNStrToReal(argv[idx]);
      retval = SUNAdaptController_SetParams_ImExGus(C, rarg1, rarg2, rarg3,
                                                    rarg4);
      if (retval != SUN_SUCCESS)
      {
        free(prefix);
        return retval;
      }
      continue;
    }

    /* check whether it was requested that all parameters be printed to screen */
    if (strcmp(argv[idx] + offset, "write_parameters") == 0)
    {
      write_parameters = SUNTRUE;
      continue;
    }
  }

  /* Call SUNAdaptController_Write (if requested) now that all
     command-line options have been set -- WARNING: this knows
     nothing about MPI, so it could be redundantly written by all
     processes if requested. */
  if (write_parameters)
  {
    retval = SUNAdaptController_Write(C, stdout);
    if (retval != SUN_SUCCESS)
    {
      free(prefix);
      return retval;
    }
  }

  free(prefix);
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * Function to set ImExGus parameters
 */

SUNErrCode SUNAdaptController_SetParams_ImExGus(SUNAdaptController C,
                                                sunrealtype k1e, sunrealtype k2e,
                                                sunrealtype k1i, sunrealtype k2i)
{
  SUNFunctionBegin(C->sunctx);
  SACIMEXGUS_K1E(C) = k1e;
  SACIMEXGUS_K2E(C) = k2e;
  SACIMEXGUS_K1I(C) = k1i;
  SACIMEXGUS_K2I(C) = k2i;
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type sunAdaptControllerGetType_ImExGus(
  SUNDIALS_MAYBE_UNUSED SUNAdaptController C)
{
  return SUN_ADAPTCONTROLLER_H;
}

SUNErrCode sunAdaptControllerEstimateStep_ImExGus(SUNAdaptController C,
                                                  sunrealtype h, int p,
                                                  sunrealtype dsm,
                                                  sunrealtype* hnew)
{
  SUNFunctionBegin(C->sunctx);

  /* order parameter to use */
  const int ord = p + 1;

  /* compute estimated time step size, modifying the first step formula */
  if (SACIMEXGUS_FIRSTSTEP(C))
  {
    /* set usable time-step adaptivity parameters -- first step */
    const sunrealtype k = -SUN_RCONST(1.0) / ord;
    const sunrealtype e = SACIMEXGUS_BIAS(C) * dsm;
    *hnew               = h * SUNRpowerR(e, k);
  }
  else
  {
    /* set usable time-step adaptivity parameters -- subsequent steps */
    const sunrealtype k1e  = -SACIMEXGUS_K1E(C) / ord;
    const sunrealtype k2e  = -SACIMEXGUS_K2E(C) / ord;
    const sunrealtype k1i  = -SACIMEXGUS_K1I(C) / ord;
    const sunrealtype k2i  = -SACIMEXGUS_K2I(C) / ord;
    const sunrealtype e1   = SACIMEXGUS_BIAS(C) * dsm;
    const sunrealtype e2   = e1 / SACIMEXGUS_EP(C);
    const sunrealtype hrat = h / SACIMEXGUS_HP(C);
    *hnew = h * SUNMIN(hrat * SUNRpowerR(e1, k1i) * SUNRpowerR(e2, k2i),
                       SUNRpowerR(e1, k1e) * SUNRpowerR(e2, k2e));
  }

  /* return with success */
  return SUN_SUCCESS;
}

SUNErrCode sunAdaptControllerReset_ImExGus(SUNAdaptController C)
{
  SACIMEXGUS_EP(C)        = SUN_RCONST(1.0);
  SACIMEXGUS_FIRSTSTEP(C) = SUNTRUE;
  return SUN_SUCCESS;
}

SUNErrCode sunAdaptControllerSetDefaults_ImExGus(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  SACIMEXGUS_BIAS(C) = DEFAULT_BIAS;
  return SUNAdaptController_SetParams_ImExGus(C, DEFAULT_K1E, DEFAULT_K2E,
                                              DEFAULT_K1I, DEFAULT_K2I);
}

SUNErrCode sunAdaptControllerWrite_ImExGus(SUNAdaptController C, FILE* fptr)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(fptr, SUN_ERR_ARG_CORRUPT);
  fprintf(fptr, "ImEx Gustafsson SUNAdaptController module:\n");
  fprintf(fptr, "  k1e = " SUN_FORMAT_G "\n", SACIMEXGUS_K1E(C));
  fprintf(fptr, "  k2e = " SUN_FORMAT_G "\n", SACIMEXGUS_K2E(C));
  fprintf(fptr, "  k1i = " SUN_FORMAT_G "\n", SACIMEXGUS_K1I(C));
  fprintf(fptr, "  k2i = " SUN_FORMAT_G "\n", SACIMEXGUS_K2I(C));
  fprintf(fptr, "  bias factor = " SUN_FORMAT_G "\n", SACIMEXGUS_BIAS(C));
  fprintf(fptr, "  previous error = " SUN_FORMAT_G "\n", SACIMEXGUS_EP(C));
  fprintf(fptr, "  previous step = " SUN_FORMAT_G "\n", SACIMEXGUS_HP(C));
  return SUN_SUCCESS;
}

SUNErrCode sunAdaptControllerSetErrorBias_ImExGus(SUNAdaptController C,
                                                  sunrealtype bias)
{
  SUNFunctionBegin(C->sunctx);
  /* set allowed value, otherwise set default */
  if (bias <= SUN_RCONST(0.0)) { SACIMEXGUS_BIAS(C) = DEFAULT_BIAS; }
  else { SACIMEXGUS_BIAS(C) = bias; }

  return SUN_SUCCESS;
}

SUNErrCode sunAdaptControllerUpdateH_ImExGus(SUNAdaptController C,
                                             sunrealtype h, sunrealtype dsm)
{
  SUNFunctionBegin(C->sunctx);
  SACIMEXGUS_EP(C)        = SACIMEXGUS_BIAS(C) * dsm;
  SACIMEXGUS_HP(C)        = h;
  SACIMEXGUS_FIRSTSTEP(C) = SUNFALSE;
  return SUN_SUCCESS;
}

/* Deprecated concrete operation wrappers */

SUNErrCode SUNAdaptController_EstimateStep_ImExGus(SUNAdaptController C,
                                                   sunrealtype h, int p,
                                                   sunrealtype dsm,
                                                   sunrealtype* hnew)
{
  return sunAdaptControllerEstimateStep_ImExGus(C, h, p, dsm, hnew);
}

SUNAdaptController_Type SUNAdaptController_GetType_ImExGus(SUNAdaptController C)
{
  return sunAdaptControllerGetType_ImExGus(C);
}

SUNErrCode SUNAdaptController_Reset_ImExGus(SUNAdaptController C)
{
  return sunAdaptControllerReset_ImExGus(C);
}

SUNErrCode SUNAdaptController_SetDefaults_ImExGus(SUNAdaptController C)
{
  return sunAdaptControllerSetDefaults_ImExGus(C);
}

SUNErrCode SUNAdaptController_SetErrorBias_ImExGus(SUNAdaptController C,
                                                   sunrealtype bias)
{
  return sunAdaptControllerSetErrorBias_ImExGus(C, bias);
}

SUNErrCode SUNAdaptController_UpdateH_ImExGus(SUNAdaptController C,
                                              sunrealtype h, sunrealtype dsm)
{
  return sunAdaptControllerUpdateH_ImExGus(C, h, dsm);
}

SUNErrCode SUNAdaptController_Write_ImExGus(SUNAdaptController C, FILE* fptr)
{
  return sunAdaptControllerWrite_ImExGus(C, fptr);
}
