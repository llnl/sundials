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
 * SUNAdaptController_Soderlind module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_core.h>
#include <sundials/sundials_errors.h>

#include "sundials_cli.h"
#include "sundials_macros.h"

/* ---------------
 * Macro accessors
 * --------------- */

#define SODERLIND_CONTENT(C)     ((SUNAdaptControllerContent_Soderlind)(C->content))
#define SODERLIND_K1(C)          (SODERLIND_CONTENT(C)->k1)
#define SODERLIND_K2(C)          (SODERLIND_CONTENT(C)->k2)
#define SODERLIND_K3(C)          (SODERLIND_CONTENT(C)->k3)
#define SODERLIND_K4(C)          (SODERLIND_CONTENT(C)->k4)
#define SODERLIND_K5(C)          (SODERLIND_CONTENT(C)->k5)
#define SODERLIND_BIAS(C)        (SODERLIND_CONTENT(C)->bias)
#define SODERLIND_EP(C)          (SODERLIND_CONTENT(C)->ep)
#define SODERLIND_EPP(C)         (SODERLIND_CONTENT(C)->epp)
#define SODERLIND_HP(C)          (SODERLIND_CONTENT(C)->hp)
#define SODERLIND_HPP(C)         (SODERLIND_CONTENT(C)->hpp)
#define SODERLIND_FIRSTSTEPS(C)  (SODERLIND_CONTENT(C)->firststeps)
#define SODERLIND_HISTORYSIZE(C) (SODERLIND_CONTENT(C)->historysize)

/* ------------------
 * Default parameters
 * ------------------ */

#define DEFAULT_K1        SUN_RCONST(1.25) /* H_{0}321 parameters */
#define DEFAULT_K2        SUN_RCONST(0.5)
#define DEFAULT_K3        SUN_RCONST(-0.75)
#define DEFAULT_K4        SUN_RCONST(0.25)
#define DEFAULT_K5        SUN_RCONST(0.75)
#define DEFAULT_PID_K1    SUN_RCONST(0.58) /* PID parameters */
#define DEFAULT_PID_K2    -SUN_RCONST(0.21)
#define DEFAULT_PID_K3    SUN_RCONST(0.1)
#define DEFAULT_PI_K1     SUN_RCONST(0.8) /* PI parameters */
#define DEFAULT_PI_K2     -SUN_RCONST(0.31)
#define DEFAULT_I_K1      SUN_RCONST(1.0)   /* I parameters */
#define DEFAULT_EXPGUS_K1 SUN_RCONST(0.367) /* Explicit Gustafsson parameters */
#define DEFAULT_EXPGUS_K2 SUN_RCONST(0.268)
#define DEFAULT_IMPGUS_K1 SUN_RCONST(0.98) /* Implicit Gustafsson parameters */
#define DEFAULT_IMPGUS_K2 SUN_RCONST(0.95)
#define DEFAULT_BIAS      SUN_RCONST(1.0)

/*
 * ----------------------------------------------------------------------------
 * Un-exported implementation specific routines
 * ----------------------------------------------------------------------------
 */

static SUNErrCode setFromCommandLine_Soderlind(SUNAdaptController C,
                                               const char* Cid, int argc,
                                               char* argv[]);

SUNErrCode SUNAdaptController_SetOptions_Soderlind(SUNAdaptController C,
                                                   const char* Cid,
                                                   const char* file_name,
                                                   int argc, char* argv[]);

/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new Soderlind controller (a.k.a., H_{0}321)
 */

SUNAdaptController SUNAdaptController_Soderlind(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  /* Create an empty controller object */
  SUNAdaptController C = SUNAdaptController_NewEmpty(sunctx);
  SUNCheckLastErrNull();
  SUNAdaptControllerContent_Soderlind content = NULL;

  /* Attach operations */
  C->ops->gettype      = SUNAdaptController_GetType_Soderlind;
  C->ops->estimatestep = SUNAdaptController_EstimateStep_Soderlind;
  C->ops->reset        = SUNAdaptController_Reset_Soderlind;
  C->ops->setoptions   = SUNAdaptController_SetOptions_Soderlind;
  C->ops->setdefaults  = SUNAdaptController_SetDefaults_Soderlind;
  C->ops->write        = SUNAdaptController_Write_Soderlind;
  C->ops->seterrorbias = SUNAdaptController_SetErrorBias_Soderlind;
  C->ops->updateh      = SUNAdaptController_UpdateH_Soderlind;
  C->ops->space        = SUNAdaptController_Space_Soderlind;

  /* Create content */
  content = (SUNAdaptControllerContent_Soderlind)malloc(sizeof(*content));
  SUNAssertNull(content, SUN_ERR_MALLOC_FAIL);
  C->content = content;

  /* Fill content with default/reset values */
  SUNCheckCallNull(SUNAdaptController_SetDefaults_Soderlind(C));
  SUNCheckCallNull(SUNAdaptController_Reset_Soderlind(C));

  return (C);
}

/* ----------------------------------------------------------------------------
 * Function to control set routines via the command line or file
 */

SUNErrCode SUNAdaptController_SetOptions_Soderlind(
  SUNAdaptController C, const char* Cid,
  SUNDIALS_MAYBE_UNUSED const char* file_name, int argc, char* argv[])
{
  SUNFunctionBegin(C->sunctx);

  /* File-based option control is currently unimplemented */
  SUNAssert((file_name == NULL || strlen(file_name) == 0),
            SUN_ERR_ARG_INCOMPATIBLE);

  if (argc > 0 && argv != NULL)
  {
    SUNCheckCall(setFromCommandLine_Soderlind(C, Cid, argc, argv));
  }

  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * Function to control Soderlind parameters from the command line
 */

static SUNErrCode setFromCommandLine_Soderlind(SUNAdaptController C,
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

    /* control over SetParams_Soderlind function */
    if (strcmp(argv[idx] + offset, "params_soderlind") == 0)
    {
      idx += 1;
      sunrealtype rarg1 = SUNStrToReal(argv[idx]);
      idx += 1;
      sunrealtype rarg2 = SUNStrToReal(argv[idx]);
      idx += 1;
      sunrealtype rarg3 = SUNStrToReal(argv[idx]);
      idx += 1;
      sunrealtype rarg4 = SUNStrToReal(argv[idx]);
      idx += 1;
      sunrealtype rarg5 = SUNStrToReal(argv[idx]);
      retval = SUNAdaptController_SetParams_Soderlind(C, rarg1, rarg2, rarg3,
                                                      rarg4, rarg5);
      if (retval != SUN_SUCCESS)
      {
        free(prefix);
        return retval;
      }
      continue;
    }

    /* control over SetParams_PID function */
    if (strcmp(argv[idx] + offset, "params_pid") == 0)
    {
      idx += 1;
      sunrealtype rarg1 = SUNStrToReal(argv[idx]);
      idx += 1;
      sunrealtype rarg2 = SUNStrToReal(argv[idx]);
      idx += 1;
      sunrealtype rarg3 = SUNStrToReal(argv[idx]);
      retval = SUNAdaptController_SetParams_PID(C, rarg1, rarg2, rarg3);
      if (retval != SUN_SUCCESS)
      {
        free(prefix);
        return retval;
      }
      continue;
    }

    /* control over SetParams_PI function */
    if (strcmp(argv[idx] + offset, "params_pi") == 0)
    {
      idx += 1;
      sunrealtype rarg1 = SUNStrToReal(argv[idx]);
      idx += 1;
      sunrealtype rarg2 = SUNStrToReal(argv[idx]);
      retval            = SUNAdaptController_SetParams_PI(C, rarg1, rarg2);
      if (retval != SUN_SUCCESS)
      {
        free(prefix);
        return retval;
      }
      continue;
    }

    /* control over SetParams_I function */
    if (strcmp(argv[idx] + offset, "params_i") == 0)
    {
      idx += 1;
      sunrealtype rarg1 = SUNStrToReal(argv[idx]);
      retval            = SUNAdaptController_SetParams_I(C, rarg1);
      if (retval != SUN_SUCCESS)
      {
        free(prefix);
        return retval;
      }
      continue;
    }

    /* control over SetParams_ExpGus function */
    if (strcmp(argv[idx] + offset, "params_expgus") == 0)
    {
      idx += 1;
      sunrealtype rarg1 = SUNStrToReal(argv[idx]);
      idx += 1;
      sunrealtype rarg2 = SUNStrToReal(argv[idx]);
      retval            = SUNAdaptController_SetParams_ExpGus(C, rarg1, rarg2);
      if (retval != SUN_SUCCESS)
      {
        free(prefix);
        return retval;
      }
      continue;
    }

    /* control over SetParams_ImpGus function */
    if (strcmp(argv[idx] + offset, "params_impgus") == 0)
    {
      idx += 1;
      sunrealtype rarg1 = SUNStrToReal(argv[idx]);
      idx += 1;
      sunrealtype rarg2 = SUNStrToReal(argv[idx]);
      retval            = SUNAdaptController_SetParams_ImpGus(C, rarg1, rarg2);
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
 * Function to set Soderlind parameters
 */

SUNErrCode SUNAdaptController_SetParams_Soderlind(SUNAdaptController C,
                                                  sunrealtype k1,
                                                  sunrealtype k2, sunrealtype k3,
                                                  sunrealtype k4, sunrealtype k5)
{
  SUNFunctionBegin(C->sunctx);
  SODERLIND_K1(C) = k1;
  SODERLIND_K2(C) = k2;
  SODERLIND_K3(C) = k3;
  SODERLIND_K4(C) = k4;
  SODERLIND_K5(C) = k5;

  if (k5 != SUN_RCONST(0.0) || k3 != SUN_RCONST(0.0))
  {
    SODERLIND_HISTORYSIZE(C) = 2;
  }
  else if (k4 != SUN_RCONST(0.0) || k2 != SUN_RCONST(0.0))
  {
    SODERLIND_HISTORYSIZE(C) = 1;
  }
  else { SODERLIND_HISTORYSIZE(C) = 0; }
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * Function to create a PID controller (subset of Soderlind)
 */

SUNAdaptController SUNAdaptController_PID(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C = SUNAdaptController_Soderlind(sunctx);
  SUNCheckLastErrNull();

  SUNCheckCallNull(SUNAdaptController_SetParams_PID(C, DEFAULT_PID_K1,
                                                    DEFAULT_PID_K2,
                                                    DEFAULT_PID_K3));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set PID parameters
 */

SUNErrCode SUNAdaptController_SetParams_PID(SUNAdaptController C, sunrealtype k1,
                                            sunrealtype k2, sunrealtype k3)
{
  SUNFunctionBegin(C->sunctx);
  return SUNAdaptController_SetParams_Soderlind(C, k1, k2, k3, SUN_RCONST(0.0),
                                                SUN_RCONST(0.0));
}

/* -----------------------------------------------------------------
 * Function to create a PI controller (subset of Soderlind)
 */

SUNAdaptController SUNAdaptController_PI(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C = SUNAdaptController_Soderlind(sunctx);
  SUNCheckLastErrNull();

  SUNCheckCallNull(
    SUNAdaptController_SetParams_PI(C, DEFAULT_PI_K1, DEFAULT_PI_K2));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set PI parameters
 */

SUNErrCode SUNAdaptController_SetParams_PI(SUNAdaptController C, sunrealtype k1,
                                           sunrealtype k2)
{
  SUNFunctionBegin(C->sunctx);
  return SUNAdaptController_SetParams_Soderlind(C, k1, k2, SUN_RCONST(0.0),
                                                SUN_RCONST(0.0), SUN_RCONST(0.0));
}

/* -----------------------------------------------------------------
 * Function to create an I controller (subset of Soderlind)
 */

SUNAdaptController SUNAdaptController_I(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C = SUNAdaptController_Soderlind(sunctx);
  SUNCheckLastErrNull();

  SUNCheckCallNull(SUNAdaptController_SetParams_I(C, DEFAULT_I_K1));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set PI parameters
 */

SUNErrCode SUNAdaptController_SetParams_I(SUNAdaptController C, sunrealtype k1)
{
  SUNFunctionBegin(C->sunctx);
  return SUNAdaptController_SetParams_Soderlind(C, k1, SUN_RCONST(0.0),
                                                SUN_RCONST(0.0), SUN_RCONST(0.0),
                                                SUN_RCONST(0.0));
}

/* -----------------------------------------------------------------
 * Function to create an explicit Gustafsson controller
 */

SUNAdaptController SUNAdaptController_ExpGus(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C = SUNAdaptController_Soderlind(sunctx);
  SUNCheckLastErrNull();

  SUNCheckCallNull(SUNAdaptController_SetParams_ExpGus(C, DEFAULT_EXPGUS_K1,
                                                       DEFAULT_EXPGUS_K2));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set explicit Gustafsson parameters
 */

SUNErrCode SUNAdaptController_SetParams_ExpGus(SUNAdaptController C,
                                               sunrealtype k1, sunrealtype k2)
{
  SUNFunctionBegin(C->sunctx);
  return SUNAdaptController_SetParams_Soderlind(C, k1 + k2, -k2, SUN_RCONST(0.0),
                                                SUN_RCONST(0.0), SUN_RCONST(0.0));
}

/* -----------------------------------------------------------------
 * Function to create an implicit Gustafsson controller
 */

SUNAdaptController SUNAdaptController_ImpGus(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C = SUNAdaptController_Soderlind(sunctx);
  SUNCheckLastErrNull();

  SUNCheckCallNull(SUNAdaptController_SetParams_ImpGus(C, DEFAULT_IMPGUS_K1,
                                                       DEFAULT_IMPGUS_K2));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to set explicit Gustafsson parameters
 */

SUNErrCode SUNAdaptController_SetParams_ImpGus(SUNAdaptController C,
                                               sunrealtype k1, sunrealtype k2)
{
  SUNFunctionBegin(C->sunctx);
  return SUNAdaptController_SetParams_Soderlind(C, k1 + k2, -k2, SUN_RCONST(0.0),
                                                SUN_RCONST(1.0), SUN_RCONST(0.0));
  return SUN_SUCCESS;
}

/* -----------------------------------------------------------------
 * Function to create an H_{0}211 controller (subset of Soderlind)
 */

SUNAdaptController SUNAdaptController_H0211(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C = SUNAdaptController_Soderlind(sunctx);
  SUNCheckLastErrNull();

  SUNCheckCallNull(
    SUNAdaptController_SetParams_Soderlind(C, SUN_RCONST(0.5), SUN_RCONST(0.5),
                                           SUN_RCONST(0.0), -SUN_RCONST(0.5),
                                           SUN_RCONST(0.0)));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to create an H_{0}321 controller (subset of Soderlind)
 */

SUNAdaptController SUNAdaptController_H0321(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C = SUNAdaptController_Soderlind(sunctx);
  SUNCheckLastErrNull();

  SUNCheckCallNull(
    SUNAdaptController_SetParams_Soderlind(C, SUN_RCONST(1.25), SUN_RCONST(0.5),
                                           -SUN_RCONST(0.75), SUN_RCONST(0.25),
                                           SUN_RCONST(0.75)));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to create an H211 controller (subset of Soderlind)
 */

SUNAdaptController SUNAdaptController_H211(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C = SUNAdaptController_Soderlind(sunctx);
  SUNCheckLastErrNull();

  SUNCheckCallNull(
    SUNAdaptController_SetParams_Soderlind(C, SUN_RCONST(0.25),
                                           SUN_RCONST(0.25), SUN_RCONST(0.0),
                                           -SUN_RCONST(0.25), SUN_RCONST(0.0)));

  return (C);
}

/* -----------------------------------------------------------------
 * Function to create an H312 controller (subset of Soderlind)
 */

SUNAdaptController SUNAdaptController_H312(SUNContext sunctx)
{
  SUNFunctionBegin(sunctx);

  SUNAdaptController C = SUNAdaptController_Soderlind(sunctx);
  SUNCheckLastErrNull();

  SUNCheckCallNull(
    SUNAdaptController_SetParams_Soderlind(C, SUN_RCONST(1.0) / SUN_RCONST(8.0),
                                           SUN_RCONST(0.25),
                                           SUN_RCONST(1.0) / SUN_RCONST(8.0),
                                           -SUN_RCONST(3.0) / SUN_RCONST(8.0),
                                           -SUN_RCONST(1.0) / SUN_RCONST(8.0)));

  return (C);
}

/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNAdaptController_Type SUNAdaptController_GetType_Soderlind(
  SUNDIALS_MAYBE_UNUSED SUNAdaptController C)
{
  return SUN_ADAPTCONTROLLER_H;
}

SUNErrCode SUNAdaptController_EstimateStep_Soderlind(SUNAdaptController C,
                                                     sunrealtype h, int p,
                                                     sunrealtype dsm,
                                                     sunrealtype* hnew)
{
  SUNFunctionBegin(C->sunctx);

  /* order parameter to use */
  const int ord        = p + 1;
  const sunrealtype e1 = SODERLIND_BIAS(C) * dsm;

  /* Handle the case of insufficient history */
  if (SODERLIND_FIRSTSTEPS(C) < SODERLIND_HISTORYSIZE(C))
  {
    /* Fall back onto an I controller */
    *hnew = h * SUNRpowerR(e1, -SUN_RCONST(1.0) / ord);
    return SUN_SUCCESS;
  }

  const sunrealtype k1 = -SODERLIND_K1(C) / ord;
  *hnew                = h * SUNRpowerR(e1, k1);

  /* This branching is not ideal, but it's more efficient than computing extra
   * math operations with degenerate k values. */
  /* TODO(SBR): Consider making separate subclasses for I and PI controllers for
   * improved efficiency */
  if (SODERLIND_HISTORYSIZE(C) > 0)
  {
    const sunrealtype k2    = -SODERLIND_K2(C) / ord;
    const sunrealtype hrat1 = h / SODERLIND_HP(C);
    *hnew *= SUNRpowerR(SODERLIND_EP(C), k2) * SUNRpowerR(hrat1, SODERLIND_K4(C));

    if (SODERLIND_HISTORYSIZE(C) > 1)
    {
      const sunrealtype k3    = -SODERLIND_K3(C) / ord;
      const sunrealtype hrat2 = SODERLIND_HP(C) / SODERLIND_HPP(C);
      *hnew *= SUNRpowerR(SODERLIND_EPP(C), k3) *
               SUNRpowerR(hrat2, SODERLIND_K5(C));
    }
  }

  /* return with success */
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Reset_Soderlind(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  SODERLIND_EP(C)         = SUN_RCONST(1.0);
  SODERLIND_EPP(C)        = SUN_RCONST(1.0);
  SODERLIND_HP(C)         = SUN_RCONST(1.0);
  SODERLIND_HPP(C)        = SUN_RCONST(1.0);
  SODERLIND_FIRSTSTEPS(C) = 0;
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_SetDefaults_Soderlind(SUNAdaptController C)
{
  SUNFunctionBegin(C->sunctx);
  SODERLIND_BIAS(C) = DEFAULT_BIAS;
  return SUNAdaptController_SetParams_Soderlind(C, DEFAULT_K1, DEFAULT_K2,
                                                DEFAULT_K3, DEFAULT_K4,
                                                DEFAULT_K5);
}

SUNErrCode SUNAdaptController_Write_Soderlind(SUNAdaptController C, FILE* fptr)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(fptr, SUN_ERR_ARG_CORRUPT);
  fprintf(fptr, "Soderlind SUNAdaptController module:\n");
  fprintf(fptr, "  k1 = " SUN_FORMAT_G "\n", SODERLIND_K1(C));
  fprintf(fptr, "  k2 = " SUN_FORMAT_G "\n", SODERLIND_K2(C));
  fprintf(fptr, "  k3 = " SUN_FORMAT_G "\n", SODERLIND_K3(C));
  fprintf(fptr, "  k4 = " SUN_FORMAT_G "\n", SODERLIND_K4(C));
  fprintf(fptr, "  k5 = " SUN_FORMAT_G "\n", SODERLIND_K5(C));
  fprintf(fptr, "  bias factor = " SUN_FORMAT_G "\n", SODERLIND_BIAS(C));
  fprintf(fptr, "  previous error = " SUN_FORMAT_G "\n", SODERLIND_EP(C));
  fprintf(fptr, "  previous-previous error = " SUN_FORMAT_G "\n",
          SODERLIND_EPP(C));
  fprintf(fptr, "  previous step = " SUN_FORMAT_G "\n", SODERLIND_HP(C));
  fprintf(fptr, "  previous-previous step = " SUN_FORMAT_G "\n",
          SODERLIND_HPP(C));
  fprintf(fptr, "  firststeps = %i\n", SODERLIND_FIRSTSTEPS(C));
  fprintf(fptr, "  historysize = %i\n", SODERLIND_HISTORYSIZE(C));
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_SetErrorBias_Soderlind(SUNAdaptController C,
                                                     sunrealtype bias)
{
  SUNFunctionBegin(C->sunctx);

  /* set allowed value, otherwise set default */
  if (bias <= SUN_RCONST(0.0)) { SODERLIND_BIAS(C) = DEFAULT_BIAS; }
  else { SODERLIND_BIAS(C) = bias; }

  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_UpdateH_Soderlind(SUNAdaptController C,
                                                sunrealtype h, sunrealtype dsm)
{
  SUNFunctionBegin(C->sunctx);
  SODERLIND_EPP(C) = SODERLIND_EP(C);
  SODERLIND_EP(C)  = SODERLIND_BIAS(C) * dsm;
  SODERLIND_HPP(C) = SODERLIND_HP(C);
  SODERLIND_HP(C)  = h;
  if (SODERLIND_FIRSTSTEPS(C) < SODERLIND_HISTORYSIZE(C))
  {
    SODERLIND_FIRSTSTEPS(C) += 1;
  }
  return SUN_SUCCESS;
}

SUNErrCode SUNAdaptController_Space_Soderlind(SUNAdaptController C,
                                              long int* lenrw, long int* leniw)
{
  SUNFunctionBegin(C->sunctx);
  SUNAssert(lenrw, SUN_ERR_ARG_CORRUPT);
  SUNAssert(leniw, SUN_ERR_ARG_CORRUPT);
  *lenrw = 10;
  *leniw = 2;
  return SUN_SUCCESS;
}
