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
 * This is the implementation file for the
 * SUNAdaptController_MRIHTol module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sunadaptcontroller/sunadaptcontroller_mrihtol_deprecated.h>
#include <sundials/sundials_core.h>
#include "sundials/priv/sundials_errors_impl.h"
#include "sundials/sundials_errors.h"

#include "sundials_cli.h"
#include "sundials_macros.h"

/* ------------------
 * Default parameters
 * ------------------ */

/*   maximum relative change for inner tolerance factor */
#define INNER_MAX_RELCH SUN_RCONST(20.0)
/*   minimum tolerance factor for inner solver */
#define INNER_MIN_TOLFAC SUN_RCONST(1.e-5)
/*   maximum tolerance factor for inner solver */
#define INNER_MAX_TOLFAC SUN_RCONST(1.0)

/* ---------------
 * Macro accessors
 * --------------- */

#define MRIHTOL_CONTENT(C)          ((SUNAdaptControllerContent_MRIHTol)(C->content))
#define MRIHTOL_CSLOW(C)            (MRIHTOL_CONTENT(C)->HControl)
#define MRIHTOL_CFAST(C)            (MRIHTOL_CONTENT(C)->TolControl)
#define MRIHTOL_INNER_MAX_RELCH(C)  (MRIHTOL_CONTENT(C)->inner_max_relch)
#define MRIHTOL_INNER_MIN_TOLFAC(C) (MRIHTOL_CONTENT(C)->inner_min_tolfac)
#define MRIHTOL_INNER_MAX_TOLFAC(C) (MRIHTOL_CONTENT(C)->inner_max_tolfac)

/*
 * ----------------------------------------------------------------------------
 * Un-exported implementation specific routines
 * ----------------------------------------------------------------------------
 */

static int sunAdaptControllerEstimateStepTol_MRIHTol(
  SUNAdaptController C, sunrealtype H, sunrealtype tolfac, int P,
  sunrealtype DSM, sunrealtype dsm, sunrealtype* Hnew, sunrealtype* tolfacnew);

static SUNAdaptController_Type sunAdaptControllerGetType_MRIHTol(
  SUNAdaptController C);

static int sunAdaptControllerReset_MRIHTol(SUNAdaptController C);

static int sunAdaptControllerSetDefaults_MRIHTol(SUNAdaptController C);

static int sunAdaptControllerSetErrorBias_MRIHTol(SUNAdaptController C,
                                                  sunrealtype bias);

static int sunAdaptControllerUpdateMRIHTol_MRIHTol(SUNAdaptController C,
                                                   sunrealtype H,
                                                   sunrealtype tolfac,
                                                   sunrealtype DSM,
                                                   sunrealtype dsm);

static int sunAdaptControllerWrite_MRIHTol(SUNAdaptController C, FILE* fptr);

static SUNErrCode setFromCommandLine_MRIHTol(SUNAdaptController C,
                                             const char* Cid, int argc,
                                             char* argv[]);

SUNErrCode SUNAdaptController_SetOptions_MRIHTol(SUNAdaptController C,
                                                 const char* Cid,
                                                 const char* file_name,
                                                 int argc, char* argv[]);

/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new MRIHTol controller
 */

SUNAdaptController SUNAdaptController_MRIHTol(SUNAdaptController HControl,
                                              SUNAdaptController TolControl,
                                              SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C;
  SUNAdaptControllerContent_MRIHTol content;

  /* Verify that input controllers have the appropriate type */
  SUNAssertNull(SUNAdaptController_GetType(HControl) == SUN_ADAPTCONTROLLER_H,
                SUN_ERR_ARG_INCOMPATIBLE);
  SUNAssertNull(SUNAdaptController_GetType(TolControl) == SUN_ADAPTCONTROLLER_H,
                SUN_ERR_ARG_INCOMPATIBLE);

  /* Create an empty controller object */
  C = NULL;
  C = SUNAdaptController_NewEmpty(sunctx);
  SUNCheckLastErrNull();

  /* Attach operations */
  C->ops->gettype         = sunAdaptControllerGetType_MRIHTol;
  C->ops->estimatesteptol = sunAdaptControllerEstimateStepTol_MRIHTol;
  C->ops->reset           = sunAdaptControllerReset_MRIHTol;
  C->ops->setoptions      = SUNAdaptController_SetOptions_MRIHTol;
  C->ops->setdefaults     = sunAdaptControllerSetDefaults_MRIHTol;
  C->ops->write           = sunAdaptControllerWrite_MRIHTol;
  C->ops->seterrorbias    = sunAdaptControllerSetErrorBias_MRIHTol;
  C->ops->updatemrihtol   = sunAdaptControllerUpdateMRIHTol_MRIHTol;
  /* Create content */
  content = NULL;
  content = (SUNAdaptControllerContent_MRIHTol)malloc(sizeof *content);
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);

  /* Attach input controllers */
  content->HControl   = HControl;
  content->TolControl = TolControl;

  /* Set parameters to default values */
  content->inner_max_relch  = INNER_MAX_RELCH;
  content->inner_min_tolfac = INNER_MIN_TOLFAC;
  content->inner_max_tolfac = INNER_MAX_TOLFAC;

  /* Attach content */
  C->content = content;

  return C;
}

/* ----------------------------------------------------------------------------
 * Function to control set routines via the command line or file
 */

SUNErrCode SUNAdaptController_SetOptions_MRIHTol(
  SUNAdaptController C, const char* Cid,
  SUNDIALS_MAYBE_UNUSED const char* file_name, int argc, char* argv[])
{
  SUNFunctionBegin(C->sunctx);

  /* File-based option control is currently unimplemented */
  SUNAssert((file_name == NULL || strlen(file_name) == 0),
            SUN_ERR_ARG_INCOMPATIBLE);

  if (argc > 0 && argv != NULL)
  {
    SUNCheckCall(setFromCommandLine_MRIHTol(C, Cid, argc, argv));
  }

  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * Function to control MRIHTol parameters from the command line
 */

static SUNErrCode setFromCommandLine_MRIHTol(SUNAdaptController C,
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
    if (strcmp(argv[idx] + offset, "params_mrihtol") == 0)
    {
      idx += 1;
      sunrealtype rarg1 = SUNStrToReal(argv[idx]);
      idx += 1;
      sunrealtype rarg2 = SUNStrToReal(argv[idx]);
      idx += 1;
      sunrealtype rarg3 = SUNStrToReal(argv[idx]);
      retval = SUNAdaptController_SetParams_MRIHTol(C, rarg1, rarg2, rarg3);
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
 * Function to set MRIHTol parameters
 */

SUNErrCode SUNAdaptController_SetParams_MRIHTol(SUNAdaptController C,
                                                sunrealtype inner_max_relch,
                                                sunrealtype inner_min_tolfac,
                                                sunrealtype inner_max_tolfac)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(inner_max_tolfac > inner_min_tolfac, SUN_ERR_ARG_OUTOFRANGE);
  if (inner_max_relch < SUN_RCONST(1.0))
  {
    MRIHTOL_INNER_MAX_RELCH(C) = INNER_MAX_RELCH;
  }
  else { MRIHTOL_INNER_MAX_RELCH(C) = inner_max_relch; }
  if (inner_min_tolfac <= SUN_RCONST(0.0))
  {
    MRIHTOL_INNER_MIN_TOLFAC(C) = INNER_MIN_TOLFAC;
  }
  else { MRIHTOL_INNER_MIN_TOLFAC(C) = inner_min_tolfac; }
  if (inner_max_tolfac <= SUN_RCONST(0.0))
  {
    MRIHTOL_INNER_MAX_TOLFAC(C) = INNER_MAX_TOLFAC;
  }
  else { MRIHTOL_INNER_MAX_TOLFAC(C) = inner_max_tolfac; }
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * Function to get slow and fast sub-controllers
 */

SUNErrCode SUNAdaptController_GetSlowController_MRIHTol(SUNAdaptController C,
                                                        SUNAdaptController* Cslow)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(Cslow, SUN_ERR_ARG_CORRUPT);
  *Cslow = MRIHTOL_CSLOW(C);
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_GetFastController_MRIHTol(SUNAdaptController C,
                                                        SUNAdaptController* Cfast)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(Cfast, SUN_ERR_ARG_CORRUPT);
  *Cfast = MRIHTOL_CFAST(C);
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type sunAdaptControllerGetType_MRIHTol(
  SUNDIALS_MAYBE_UNUSED SUNAdaptController C)
{
  return SUN_ADAPTCONTROLLER_MRI_H_TOL;
}

SUNErrCode sunAdaptControllerEstimateStepTol_MRIHTol(
  SUNAdaptController C, sunrealtype H, sunrealtype tolfac, int P,
  sunrealtype DSM, sunrealtype dsm, sunrealtype* Hnew, sunrealtype* tolfacnew)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(Hnew, SUN_ERR_ARG_CORRUPT);
  SUNAssert(tolfacnew, SUN_ERR_ARG_CORRUPT);
  sunrealtype tolfacest;

  /* Call slow time scale sub-controller to fill Hnew -- note that all heuristics
     bounds on Hnew will be enforced by the time integrator itself */
  SUNCheckCall(SUNAdaptController_EstimateStep(MRIHTOL_CSLOW(C), H, P, DSM, Hnew));

  /* Call fast time scale sub-controller with order=1: no matter the integrator
     order, we expect its error to be proportional to the tolerance factor */
  SUNCheckCall(SUNAdaptController_EstimateStep(MRIHTOL_CFAST(C), tolfac, 0, dsm,
                                               &tolfacest));

  /* Enforce bounds on estimated tolerance factor */
  /*     keep relative change within bounds */
  tolfacest = SUNMAX(tolfacest, tolfac / MRIHTOL_INNER_MAX_RELCH(C));
  tolfacest = SUNMIN(tolfacest, tolfac * MRIHTOL_INNER_MAX_RELCH(C));
  /*     enforce absolute min/max bounds */
  tolfacest = SUNMAX(tolfacest, MRIHTOL_INNER_MIN_TOLFAC(C));
  tolfacest = SUNMIN(tolfacest, MRIHTOL_INNER_MAX_TOLFAC(C));

  /* Set result and return */
  *tolfacnew = tolfacest;
  return SUN_SUCCESS;
}

SUNErrCode sunAdaptControllerReset_MRIHTol(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  SUNCheckCall(SUNAdaptController_Reset(MRIHTOL_CSLOW(C)));
  SUNCheckCall(SUNAdaptController_Reset(MRIHTOL_CFAST(C)));
  return SUN_SUCCESS;
}

SUNErrCode sunAdaptControllerSetDefaults_MRIHTol(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  SUNCheckCall(SUNAdaptController_SetDefaults(MRIHTOL_CSLOW(C)));
  SUNCheckCall(SUNAdaptController_SetDefaults(MRIHTOL_CFAST(C)));
  MRIHTOL_INNER_MAX_RELCH(C)  = INNER_MAX_RELCH;
  MRIHTOL_INNER_MIN_TOLFAC(C) = INNER_MIN_TOLFAC;
  MRIHTOL_INNER_MAX_TOLFAC(C) = INNER_MAX_TOLFAC;
  return SUN_SUCCESS;
}

SUNErrCode sunAdaptControllerWrite_MRIHTol(SUNAdaptController C, FILE* fptr)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(fptr, SUN_ERR_ARG_CORRUPT);
  fprintf(fptr, "Multirate H-Tol SUNAdaptController module:\n");
  fprintf(fptr, "  inner_max_relch  = " SUN_FORMAT_G "\n",
          MRIHTOL_INNER_MAX_RELCH(C));
  fprintf(fptr, "  inner_min_tolfac = " SUN_FORMAT_G "\n",
          MRIHTOL_INNER_MIN_TOLFAC(C));
  fprintf(fptr, "  inner_max_tolfac = " SUN_FORMAT_G "\n",
          MRIHTOL_INNER_MAX_TOLFAC(C));
  fprintf(fptr, "\nSlow step controller:\n");
  SUNCheckCall(SUNAdaptController_Write(MRIHTOL_CSLOW(C), fptr));
  fprintf(fptr, "\nFast tolerance controller:\n");
  SUNCheckCall(SUNAdaptController_Write(MRIHTOL_CFAST(C), fptr));
  return SUN_SUCCESS;
}

SUNErrCode sunAdaptControllerSetErrorBias_MRIHTol(SUNAdaptController C,
                                                  sunrealtype bias)
{
  SUNFunctionBegin(C->sunctx);
  SUNCheckCall(SUNAdaptController_SetErrorBias(MRIHTOL_CSLOW(C), bias));
  SUNCheckCall(SUNAdaptController_SetErrorBias(MRIHTOL_CFAST(C), bias));
  return SUN_SUCCESS;
}

SUNErrCode sunAdaptControllerUpdateMRIHTol_MRIHTol(SUNAdaptController C,
                                                   sunrealtype H,
                                                   sunrealtype tolfac,
                                                   sunrealtype DSM,
                                                   sunrealtype dsm)
{
  SUNFunctionBegin(C->sunctx);
  SUNCheckCall(SUNAdaptController_UpdateH(MRIHTOL_CSLOW(C), H, DSM));
  SUNCheckCall(SUNAdaptController_UpdateH(MRIHTOL_CFAST(C), tolfac, dsm));
  return SUN_SUCCESS;
}

/* Deprecated concrete operation wrappers */

int SUNAdaptController_EstimateStepTol_MRIHTol(
  SUNAdaptController C, sunrealtype H, sunrealtype tolfac, int P,
  sunrealtype DSM, sunrealtype dsm, sunrealtype* Hnew, sunrealtype* tolfacnew)
{
  return sunAdaptControllerEstimateStepTol_MRIHTol(C, H, tolfac, P, DSM, dsm,
                                                   Hnew, tolfacnew);
}

SUNAdaptController_Type SUNAdaptController_GetType_MRIHTol(SUNAdaptController C)
{
  return sunAdaptControllerGetType_MRIHTol(C);
}

int SUNAdaptController_Reset_MRIHTol(SUNAdaptController C)
{
  return sunAdaptControllerReset_MRIHTol(C);
}

int SUNAdaptController_SetDefaults_MRIHTol(SUNAdaptController C)
{
  return sunAdaptControllerSetDefaults_MRIHTol(C);
}

int SUNAdaptController_SetErrorBias_MRIHTol(SUNAdaptController C, sunrealtype bias)
{
  return sunAdaptControllerSetErrorBias_MRIHTol(C, bias);
}

int SUNAdaptController_UpdateMRIHTol_MRIHTol(SUNAdaptController C,
                                             sunrealtype H, sunrealtype tolfac,
                                             sunrealtype DSM, sunrealtype dsm)
{
  return sunAdaptControllerUpdateMRIHTol_MRIHTol(C, H, tolfac, DSM, dsm);
}

int SUNAdaptController_Write_MRIHTol(SUNAdaptController C, FILE* fptr)
{
  return sunAdaptControllerWrite_MRIHTol(C, fptr);
}
