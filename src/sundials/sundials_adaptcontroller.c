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
 * This is the implementation file for a generic SUNAdaptController
 * package. It contains the implementation of the SUNAdaptController
 * operations listed in sundials_adaptcontroller.h
 * -----------------------------------------------------------------*/

#include <math.h>
#include <string.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_adaptcontroller.h>
#include <sundials/sundials_core.h>

#include "sundials/sundials_errors.h"

/* internal function prototypes */
SUNErrCode sunadctrlSetFromCommandLine(SUNAdaptController C, const char* Cid,
                                       int argc, char* argv[]);

/* -----------------------------------------------------------------
 * Create a new empty SUNAdaptController object
 * ----------------------------------------------------------------- */

SUNAdaptController SUNAdaptController_NewEmpty(SUNContext sunctx)
{
  SUNAdaptController C;
  SUNAdaptController_Ops ops;

  /* a context is required */
  if (sunctx == NULL) { return (NULL); }

  SUNFunctionBegin(sunctx);

  /* create controller object */
  C = NULL;
  C = (SUNAdaptController)malloc(sizeof *C);
  SUNAssertNull(C, SUN_ERR_MALLOC_FAIL);

  /* create matrix ops structure */
  ops = NULL;
  ops = (SUNAdaptController_Ops)malloc(sizeof *ops);
  SUNAssertNull(ops, SUN_ERR_MALLOC_FAIL);

  /* initialize operations to NULL */
  ops->gettype         = NULL;
  ops->destroy         = NULL;
  ops->reset           = NULL;
  ops->estimatestep    = NULL;
  ops->estimatesteptol = NULL;
  ops->setoptions      = NULL;
  ops->setdefaults     = NULL;
  ops->write           = NULL;
  ops->seterrorbias    = NULL;
  ops->updateh         = NULL;
  ops->updatemrihtol   = NULL;
  ops->space           = NULL;

  /* attach ops and initialize content to NULL */
  C->ops     = ops;
  C->content = NULL;
  C->sunctx  = sunctx;

  return (C);
}

/* -----------------------------------------------------------------
 * Free a generic SUNAdaptController (assumes content is already empty)
 * ----------------------------------------------------------------- */

void SUNAdaptController_DestroyEmpty(SUNAdaptController C)
{
  if (C == NULL) { return; }

  /* free non-NULL ops structure */
  if (C->ops) { free(C->ops); }
  C->ops = NULL;

  /* free overall SUNAdaptController object and return */
  free(C);
  return;
}

/* -----------------------------------------------------------------
 * Required functions in the 'ops' structure for non-NULL controller
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptController_GetType(SUNAdaptController C)
{
  if (C == NULL) { return SUN_ADAPTCONTROLLER_NONE; }
  if (C->ops->gettype) { return C->ops->gettype(C); }
  return SUN_ADAPTCONTROLLER_NONE;
}

/* -----------------------------------------------------------------
 * internal utility routines
 * ----------------------------------------------------------------- */

SUNErrCode sunadctrlSetFromCommandLine(SUNAdaptController C, const char* Cid,
                                       int argc, char* argv[])
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
  for (int idx = 1; idx < argc; idx++)
  {
    /* skip command-line arguments that do not begin with correct prefix */
    if (strncmp(argv[idx], prefix, strlen(prefix)) != 0) { continue; }

    /* control over SetDefaults function */
    if (strcmp(argv[idx] + offset, "defaults") == 0)
    {
      retval = SUNAdaptController_SetDefaults(C);
      if (retval != SUN_SUCCESS)
      {
        free(prefix);
        return retval;
      }
      continue;
    }

    /* control over SetErrorBias function */
    if (strcmp(argv[idx] + offset, "error_bias") == 0)
    {
      idx += 1;
      sunrealtype rarg = SUNStrToReal(argv[idx]);
      retval           = SUNAdaptController_SetErrorBias(C, rarg);
      if (retval != SUN_SUCCESS)
      {
        free(prefix);
        return retval;
      }
      continue;
    }

    /* Note that although SUNAdaptController_Write is part of the base class,
       it should NOT be called until all options have been set, so we process
       that in the implementations instead of here. */
  }

  free(prefix);
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * Optional functions in the 'ops' structure
 * ----------------------------------------------------------------- */

SUNErrCode SUNAdaptController_Destroy(SUNAdaptController C)
{
  if (C == NULL) { return (SUN_SUCCESS); }

  /* if the destroy operation exists use it */
  if (C->ops)
  {
    if (C->ops->destroy) { return (C->ops->destroy(C)); }
  }

  /* if we reach this point, either ops == NULL or destroy == NULL,
     try to cleanup by freeing the content, ops, and matrix */
  if (C->content)
  {
    free(C->content);
    C->content = NULL;
  }
  if (C->ops)
  {
    free(C->ops);
    C->ops = NULL;
  }
  free(C);
  C = NULL;

  return (SUN_SUCCESS);
}

SUNErrCode SUNAdaptController_EstimateStep(SUNAdaptController C, sunrealtype h,
                                           int p, sunrealtype dsm,
                                           sunrealtype* hnew)
{
  SUNErrCode ier = SUN_SUCCESS;
  if (C == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(C->sunctx);
  SUNAssert(isfinite(h), SUN_ERR_ARG_OUTOFRANGE);
  SUNAssert(p >= 0, SUN_ERR_ARG_OUTOFRANGE);
  SUNAssert(dsm >= SUN_RCONST(0.0), SUN_ERR_ARG_OUTOFRANGE);
  SUNAssert(hnew, SUN_ERR_ARG_CORRUPT);
  *hnew = h; /* initialize output with identity */
  if (C->ops->estimatestep) { ier = C->ops->estimatestep(C, h, p, dsm, hnew); }
  return (ier);
}

SUNErrCode SUNAdaptController_EstimateStepTol(SUNAdaptController C,
                                              sunrealtype H, sunrealtype tolfac,
                                              int P, sunrealtype DSM,
                                              sunrealtype dsm, sunrealtype* Hnew,
                                              sunrealtype* tolfacnew)
{
  SUNErrCode ier = SUN_SUCCESS;
  if (C == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(C->sunctx);
  SUNAssert(Hnew, SUN_ERR_ARG_CORRUPT);
  SUNAssert(tolfacnew, SUN_ERR_ARG_CORRUPT);
  *Hnew      = H; /* initialize outputs with identity */
  *tolfacnew = tolfac;
  if (C->ops->estimatesteptol)
  {
    ier = C->ops->estimatesteptol(C, H, tolfac, P, DSM, dsm, Hnew, tolfacnew);
  }
  return (ier);
}

SUNErrCode SUNAdaptController_Reset(SUNAdaptController C)
{
  SUNErrCode ier = SUN_SUCCESS;
  if (C == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(C->sunctx);
  if (C->ops->reset) { ier = C->ops->reset(C); }
  return (ier);
}

SUNErrCode SUNAdaptController_SetOptions(SUNAdaptController C, const char* Cid,
                                         const char* file_name, int argc,
                                         char* argv[])
{
  if (C == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(C->sunctx);

  /* File-based option control is currently unimplemented */
  SUNAssert((file_name == NULL || strlen(file_name) == 0),
            SUN_ERR_ARG_INCOMPATIBLE);

  /* First, process all base-class options */
  if (argc > 0 && argv != NULL)
  {
    SUNCheckCall(sunadctrlSetFromCommandLine(C, Cid, argc, argv));
  }

  /* Second, ask the implementation to process any remaining options */
  if (C->ops->setoptions)
  {
    return (C->ops->setoptions(C, Cid, file_name, argc, argv));
  }
  else { return (SUN_SUCCESS); }
}

SUNErrCode SUNAdaptController_SetDefaults(SUNAdaptController C)
{
  SUNErrCode ier = SUN_SUCCESS;
  if (C == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(C->sunctx);
  if (C->ops->setdefaults) { ier = C->ops->setdefaults(C); }
  return (ier);
}

SUNErrCode SUNAdaptController_Write(SUNAdaptController C, FILE* fptr)
{
  SUNErrCode ier = SUN_SUCCESS;
  if (C == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(C->sunctx);
  SUNAssert(fptr, SUN_ERR_ARG_CORRUPT);
  if (C->ops->write) { ier = C->ops->write(C, fptr); }
  return (ier);
}

SUNErrCode SUNAdaptController_SetErrorBias(SUNAdaptController C, sunrealtype bias)
{
  SUNErrCode ier = SUN_SUCCESS;
  if (C == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(C->sunctx);
  if (C->ops->seterrorbias) { ier = C->ops->seterrorbias(C, bias); }
  return (ier);
}

SUNErrCode SUNAdaptController_UpdateH(SUNAdaptController C, sunrealtype h,
                                      sunrealtype dsm)
{
  SUNErrCode ier = SUN_SUCCESS;
  if (C == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(C->sunctx);
  SUNAssert(isfinite(h), SUN_ERR_ARG_OUTOFRANGE);
  SUNAssert(dsm >= SUN_RCONST(0.0), SUN_ERR_ARG_OUTOFRANGE);
  if (C->ops->updateh) { ier = C->ops->updateh(C, h, dsm); }
  return (ier);
}

SUNErrCode SUNAdaptController_UpdateMRIHTol(SUNAdaptController C, sunrealtype H,
                                            sunrealtype tolfac, sunrealtype DSM,
                                            sunrealtype dsm)
{
  SUNErrCode ier = SUN_SUCCESS;
  if (C == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(C->sunctx);
  if (C->ops->updatemrihtol)
  {
    ier = C->ops->updatemrihtol(C, H, tolfac, DSM, dsm);
  }
  return (ier);
}

SUNErrCode SUNAdaptController_Space(SUNAdaptController C, long int* lenrw,
                                    long int* leniw)
{
  SUNErrCode ier = SUN_SUCCESS;
  if (C == NULL) { return SUN_ERR_ARG_CORRUPT; }
  SUNFunctionBegin(C->sunctx);
  SUNAssert(lenrw, SUN_ERR_ARG_CORRUPT);
  SUNAssert(leniw, SUN_ERR_ARG_CORRUPT);
  *lenrw = 0; /* initialize outputs with identity */
  *leniw = 0;
  if (C->ops->space) { ier = C->ops->space(C, lenrw, leniw); }
  return (ier);
}
