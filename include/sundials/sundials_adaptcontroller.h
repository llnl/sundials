/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
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
 * SUNDIALS accuracy-based adaptivity controller class. These
 * objects estimate step sizes for time integration methods such
 * that the next step solution should satisfy a desired temporal
 * accuracy.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_ADAPTCONTROLLER_H
#define _SUNDIALS_ADAPTCONTROLLER_H

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_context.h>

#include "sundials/sundials_types.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------------------
 * SUNAdaptController types:
 *    NONE - empty controller (does nothing)
 *    H    - controls a single-rate step size
 *    MRI_H_TOL - controls slow step and fast relative tolerances
 * ----------------------------------------------------------------- */

enum SUNAdaptController_Type
{
  SUN_ADAPTCONTROLLER_NONE,
  SUN_ADAPTCONTROLLER_H,
  SUN_ADAPTCONTROLLER_MRI_H_TOL
};

#ifndef SWIG
typedef enum SUNAdaptController_Type SUNAdaptController_Type;
#endif

/* -----------------------------------------------------------------
 * Generic definition of SUNAdaptController
 * ----------------------------------------------------------------- */

/* Forward reference for pointer to SUNAdaptController_Ops object */
typedef _SUNDIALS_STRUCT_ _generic_SUNAdaptController_Ops* SUNAdaptController_Ops;

/* Forward reference for pointer to SUNAdaptController object */
typedef _SUNDIALS_STRUCT_ _generic_SUNAdaptController* SUNAdaptController;

/* Structure containing function pointers to controller operations  */
struct _generic_SUNAdaptController_Ops
{
  /* REQUIRED of all controller implementations. */
  SUNAdaptController_Type (*gettype)(SUNAdaptController C);

  /* REQUIRED for controllers of SUN_ADAPTCONTROLLER_H type. */
  SUNErrCode (*estimatestep)(SUNAdaptController C, sunrealtype h, int p,
                             sunrealtype dsm, sunrealtype* hnew);

  /* REQUIRED for controllers of SUN_ADAPTCONTROLLER_MRI_H_TOL type. */
  SUNErrCode (*estimatesteptol)(SUNAdaptController C, sunrealtype H,
                                sunrealtype tolfac, int P, sunrealtype DSM,
                                sunrealtype dsm, sunrealtype* Hnew,
                                sunrealtype* tolfacnew);

  /* OPTIONAL for all SUNAdaptController implementations. */
  SUNErrCode (*destroy)(SUNAdaptController C);
  SUNErrCode (*reset)(SUNAdaptController C);
  SUNErrCode (*setoptions)(SUNAdaptController C, const char* Cid,
                           const char* file_name, int argc, char* argv[]);
  SUNErrCode (*setdefaults)(SUNAdaptController C);
  SUNErrCode (*write)(SUNAdaptController C, FILE* fptr);
  SUNErrCode (*seterrorbias)(SUNAdaptController C, sunrealtype bias);
  SUNErrCode (*updateh)(SUNAdaptController C, sunrealtype h, sunrealtype dsm);
  SUNErrCode (*updatemrihtol)(SUNAdaptController C, sunrealtype H,
                              sunrealtype tolfac, sunrealtype DSM,
                              sunrealtype dsm);
  SUNErrCode (*space)(SUNAdaptController C, long int* lenrw, long int* leniw);
};

/* A SUNAdaptController is a structure with an implementation-dependent
   'content' field, and a pointer to a structure of
   operations corresponding to that implementation. */
struct _generic_SUNAdaptController
{
  void* content;
  SUNAdaptController_Ops ops;
  SUNContext sunctx;
};

/* -----------------------------------------------------------------
 * Functions exported by SUNAdaptController module
 * ----------------------------------------------------------------- */

/* Function to create an empty SUNAdaptController data structure. */
SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_NewEmpty(SUNContext sunctx);

/* Function to free a generic SUNAdaptController (assumes content is already empty) */
SUNDIALS_EXPORT
void SUNAdaptController_DestroyEmpty(SUNAdaptController C);

/* Function to report the type of a SUNAdaptController object. */
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptController_GetType(SUNAdaptController C);

/* Function to deallocate a SUNAdaptController object.

   Any return value other than SUN_SUCCESS will be treated as
   an unrecoverable failure. */
SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_Destroy(SUNAdaptController C);

/* Main step size controller function.  This is called following
   a time step with size 'h' and local error factor 'dsm', and the
   controller should estimate 'hnew' so that the ensuing step
   will have 'dsm' value JUST BELOW 1.

   Any return value other than SUN_SUCCESS will be treated as
   an unrecoverable failure. */
SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_EstimateStep(SUNAdaptController C, sunrealtype h,
                                           int p, sunrealtype dsm,
                                           sunrealtype* hnew);

/* Combined slow step/fast tolerance multirate controller function.
   This is called following a slow multirate time step with size 'H'
   and fast/slow relative tolerance ratio 'tolfac', and error factors
   'DSM' and 'dsm' (slow and fast, resp.).  The controller should
   estimate slow stepsize 'Hnew' and updated relative tolerance ratio
   'tolfacnew', so that the ensuing step will have 'DSM' and 'dsm'
   values JUST BELOW 1 with minimal computational effort. */
SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_EstimateStepTol(SUNAdaptController C,
                                              sunrealtype H, sunrealtype tolfac,
                                              int P, sunrealtype DSM,
                                              sunrealtype dsm, sunrealtype* Hnew,
                                              sunrealtype* tolfacnew);

/* Function to reset the controller to its initial state, e.g., if
   it stores a small number of previous dsm or step size values. */
SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_Reset(SUNAdaptController C);

/* Function to update internal controller parameters from the command line. */
SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_SetOptions(SUNAdaptController C, const char* Cid,
                                         const char* file_name, int argc,
                                         char* argv[]);

/* Function to set the controller parameters to their default values. */
SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_SetDefaults(SUNAdaptController C);

/* Function to write all controller parameters to the indicated file
   pointer. */
SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_Write(SUNAdaptController C, FILE* fptr);

/* Function to set an error bias factor to use for scaling the local error
   'dsm' factors above. */
SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_SetErrorBias(SUNAdaptController C,
                                           sunrealtype bias);

/* Function to notify a controller of type SUN_ADAPTCONTROLLER_H that
   a successful time step was taken with stepsize h and local error factor
   dsm, indicating that these can be saved for subsequent controller functions. */
SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_UpdateH(SUNAdaptController C, sunrealtype h,
                                      sunrealtype dsm);

/* Function to notify the controller of a successful multirate time step
   with size H and fast tolerance factor tolfac, and local error factors
   DSM and dsm, indicating that the step size, tolerance factor, or local
   error factors can be saved for subsequent controller functions. */
SUNDIALS_EXPORT
SUNErrCode SUNAdaptController_UpdateMRIHTol(SUNAdaptController C, sunrealtype H,
                                            sunrealtype tolfac, sunrealtype DSM,
                                            sunrealtype dsm);

/* Function to return the memory requirements of the controller object. */
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
SUNErrCode SUNAdaptController_Space(SUNAdaptController C, long int* lenrw,
                                    long int* leniw);

#ifdef __cplusplus
}
#endif

#endif /* _SUNDIALS_ADAPTCONTROLLER_H */
