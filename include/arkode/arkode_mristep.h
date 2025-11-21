/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ UMBC
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
 * This is the header file for the ARKODE MRIStep module.
 * -----------------------------------------------------------------*/

#ifndef _MRISTEP_H
#define _MRISTEP_H

#include <arkode/arkode.h>
#include <arkode/arkode_butcher_dirk.h>
#include <arkode/arkode_butcher_erk.h>
#include <arkode/arkode_ls.h>
#include <arkode/arkode_mristep_deprecated.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>
#include <sundials/sundials_stepper.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------
 * MRIStep Constants
 * ----------------- */

/* MRIStep method types */
enum MRISTEP_METHOD_TYPE
{
  MRISTEP_EXPLICIT,
  MRISTEP_IMPLICIT,
  MRISTEP_IMEX,
  MRISTEP_MERK,
  MRISTEP_SR
};

#ifndef SWIG
typedef enum MRISTEP_METHOD_TYPE MRISTEP_METHOD_TYPE;
#endif

/* MRI coupling table IDs */
enum ARKODE_MRITableID
{
  ARKODE_MRI_NONE = -1, /* ensure enum is signed int */
  /* WARNING:  ARKODE_MIN_MRI_NUM must come after the first entry, ARKODE_MIS_KW3,
     because Python enums will only expose the member that is defined first. Due to
     this and how pybind/nanobind handle the enums, if we defined ARKODE_MRI_NUM first,
     then ARKODE_MIS_KW3 would not be usable from the module scope (the MIN/MAX) entries
     will still be usable when accessing through the IntEnum object, but not from module scope. */
  ARKODE_MIS_KW3     = 200,
  ARKODE_MIN_MRI_NUM = 200,
  ARKODE_MRI_GARK_ERK33a,
  ARKODE_MRI_GARK_ERK45a,
  ARKODE_MRI_GARK_IRK21a,
  ARKODE_MRI_GARK_ESDIRK34a,
  ARKODE_MRI_GARK_ESDIRK46a,
  ARKODE_IMEX_MRI_GARK3a,
  ARKODE_IMEX_MRI_GARK3b,
  ARKODE_IMEX_MRI_GARK4,
  ARKODE_MRI_GARK_FORWARD_EULER,
  ARKODE_MRI_GARK_RALSTON2,
  ARKODE_MRI_GARK_ERK22a,
  ARKODE_MRI_GARK_ERK22b,
  ARKODE_MRI_GARK_RALSTON3,
  ARKODE_MRI_GARK_BACKWARD_EULER,
  ARKODE_MRI_GARK_IMPLICIT_MIDPOINT,
  ARKODE_IMEX_MRI_GARK_EULER,
  ARKODE_IMEX_MRI_GARK_TRAPEZOIDAL,
  ARKODE_IMEX_MRI_GARK_MIDPOINT,
  ARKODE_MERK21,
  ARKODE_MERK32,
  ARKODE_MERK43,
  ARKODE_MERK54,
  ARKODE_IMEX_MRI_SR21,
  ARKODE_IMEX_MRI_SR32,
  ARKODE_IMEX_MRI_SR43,
  ARKODE_MAX_MRI_NUM = ARKODE_IMEX_MRI_SR43
};

#ifndef SWIG
typedef enum ARKODE_MRITableID ARKODE_MRITableID;
#endif

/* Default MRI coupling tables for each order and type */
static const int MRISTEP_DEFAULT_EXPL_1 = ARKODE_MRI_GARK_FORWARD_EULER;
static const int MRISTEP_DEFAULT_EXPL_2 = ARKODE_MRI_GARK_ERK22b;
static const int MRISTEP_DEFAULT_EXPL_3 = ARKODE_MIS_KW3;
static const int MRISTEP_DEFAULT_EXPL_4 = ARKODE_MRI_GARK_ERK45a;

static const int MRISTEP_DEFAULT_EXPL_2_AD = ARKODE_MRI_GARK_ERK22b;
static const int MRISTEP_DEFAULT_EXPL_3_AD = ARKODE_MRI_GARK_ERK33a;
static const int MRISTEP_DEFAULT_EXPL_4_AD = ARKODE_MRI_GARK_ERK45a;
static const int MRISTEP_DEFAULT_EXPL_5_AD = ARKODE_MERK54;

static const int MRISTEP_DEFAULT_IMPL_SD_1 = ARKODE_MRI_GARK_BACKWARD_EULER;
static const int MRISTEP_DEFAULT_IMPL_SD_2 = ARKODE_MRI_GARK_IRK21a;
static const int MRISTEP_DEFAULT_IMPL_SD_3 = ARKODE_MRI_GARK_ESDIRK34a;
static const int MRISTEP_DEFAULT_IMPL_SD_4 = ARKODE_MRI_GARK_ESDIRK46a;

static const int MRISTEP_DEFAULT_IMEX_SD_1 = ARKODE_IMEX_MRI_GARK_EULER;
static const int MRISTEP_DEFAULT_IMEX_SD_2 = ARKODE_IMEX_MRI_GARK_TRAPEZOIDAL;
static const int MRISTEP_DEFAULT_IMEX_SD_3 = ARKODE_IMEX_MRI_GARK3b;
static const int MRISTEP_DEFAULT_IMEX_SD_4 = ARKODE_IMEX_MRI_GARK4;

static const int MRISTEP_DEFAULT_IMEX_SD_2_AD = ARKODE_IMEX_MRI_SR21;
static const int MRISTEP_DEFAULT_IMEX_SD_3_AD = ARKODE_IMEX_MRI_SR32;
static const int MRISTEP_DEFAULT_IMEX_SD_4_AD = ARKODE_IMEX_MRI_SR43;

/* ------------------------------------
 * MRIStep Inner Stepper Function Types
 * ------------------------------------ */

typedef int (*MRIStepInnerEvolveFn)(MRIStepInnerStepper stepper, sunrealtype t0,
                                    sunrealtype tout, N_Vector y);

typedef int (*MRIStepInnerFullRhsFn)(MRIStepInnerStepper stepper, sunrealtype t,
                                     N_Vector y, N_Vector f, int mode);

typedef int (*MRIStepInnerResetFn)(MRIStepInnerStepper stepper, sunrealtype tR,
                                   N_Vector yR);

typedef int (*MRIStepInnerGetAccumulatedError)(MRIStepInnerStepper stepper,
                                               sunrealtype* accum_error);

typedef int (*MRIStepInnerResetAccumulatedError)(MRIStepInnerStepper stepper);

typedef int (*MRIStepInnerSetRTol)(MRIStepInnerStepper stepper, sunrealtype rtol);

/*---------------------------------------------------------------
  MRI coupling data structure and associated utility routines
  ---------------------------------------------------------------*/
struct MRIStepCouplingMem
{
  MRISTEP_METHOD_TYPE type; /* flag to encode the MRI method type                  */
  int nmat;         /* number of MRI coupling matrices                     */
  int stages;       /* size of coupling matrices ((stages+1) * stages)     */
  int q;            /* method order of accuracy                            */
  int p;            /* embedding order of accuracy                         */
  sunrealtype* c;   /* stage abscissae                                     */
  sunrealtype*** W; /* explicit coupling matrices [nmat][stages+1][stages] */
  sunrealtype*** G; /* implicit coupling matrices [nmat][stages+1][stages] */

  int ngroup;  /* number of stage groups (MERK-specific)              */
  int** group; /* stages to integrate together (MERK-specific)        */
};

typedef _SUNDIALS_STRUCT_ MRIStepCouplingMem* MRIStepCoupling;

/* Accessor routine to load built-in MRI table */
SUNDIALS_EXPORT MRIStepCoupling MRIStepCoupling_LoadTable(ARKODE_MRITableID method);

/* Accessor routine to load built-in MRI table from string */
SUNDIALS_EXPORT MRIStepCoupling MRIStepCoupling_LoadTableByName(const char* method);

/* Utility routines to allocate/free/output coupling table structures */
SUNDIALS_EXPORT MRIStepCoupling MRIStepCoupling_Alloc(int nmat, int stages,
                                                      MRISTEP_METHOD_TYPE type);
SUNDIALS_EXPORT MRIStepCoupling MRIStepCoupling_Create(int nmat, int stages,
                                                       int q, int p,
                                                       sunrealtype* W_1d,
                                                       sunrealtype* G_1d,
                                                       sunrealtype* c_1d);
SUNDIALS_EXPORT MRIStepCoupling MRIStepCoupling_MIStoMRI(ARKodeButcherTable B,
                                                         int q, int p);
SUNDIALS_EXPORT MRIStepCoupling MRIStepCoupling_Copy(MRIStepCoupling MRIC);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
void MRIStepCoupling_Space(MRIStepCoupling MRIC, sunindextype* liw,
                           sunindextype* lrw);
SUNDIALS_EXPORT void MRIStepCoupling_Free(MRIStepCoupling MRIC);
SUNDIALS_EXPORT void MRIStepCoupling_Write(MRIStepCoupling MRIC, FILE* outfile);

/* ------------------------------
 * User-Supplied Function Types
 * ------------------------------ */

typedef int (*MRIStepPreInnerFn)(sunrealtype t, N_Vector* f_1d, int nvecs,
                                 void* user_data);

typedef int (*MRIStepPostInnerFn)(sunrealtype t, N_Vector y, void* user_data);

/* -------------------
 * Exported Functions
 * ------------------- */

/* Creation and Reinitialization functions */
SUNDIALS_EXPORT void* MRIStepCreate(ARKRhsFn fse, ARKRhsFn fsi, sunrealtype t0,
                                    N_Vector y0, MRIStepInnerStepper stepper,
                                    SUNContext sunctx);
SUNDIALS_EXPORT int MRIStepReInit(void* arkode_mem, ARKRhsFn fse, ARKRhsFn fsi,
                                  sunrealtype t0, N_Vector y0);

/* Optional input functions -- must be called AFTER MRIStepCreate */
SUNDIALS_EXPORT int MRIStepSetCoupling(void* arkode_mem, MRIStepCoupling MRIC);
SUNDIALS_EXPORT int MRIStepSetPreInnerFn(void* arkode_mem,
                                         MRIStepPreInnerFn prefn);
SUNDIALS_EXPORT int MRIStepSetPostInnerFn(void* arkode_mem,
                                          MRIStepPostInnerFn postfn);

/* Optional output functions */
SUNDIALS_EXPORT int MRIStepGetCurrentCoupling(
  void* arkode_mem,
  MRIStepCoupling* MRIC); // nb::rv_policy::reference
SUNDIALS_EXPORT int MRIStepGetLastInnerStepFlag(void* arkode_mem, int* flag);
SUNDIALS_EXPORT int MRIStepGetNumInnerStepperFails(void* arkode_mem,
                                                   long int* inner_fails);

/* Custom inner stepper functions */
SUNDIALS_EXPORT int MRIStepInnerStepper_Create(SUNContext sunctx,
                                               MRIStepInnerStepper* stepper);

SUNDIALS_EXPORT int MRIStepInnerStepper_CreateFromSUNStepper(
  SUNStepper sunstepper, MRIStepInnerStepper* stepper);

SUNDIALS_EXPORT int MRIStepInnerStepper_Free(MRIStepInnerStepper* stepper);
SUNDIALS_EXPORT int MRIStepInnerStepper_SetContent(MRIStepInnerStepper stepper,
                                                   void* content);
SUNDIALS_EXPORT int MRIStepInnerStepper_GetContent(MRIStepInnerStepper stepper,
                                                   void** content);
SUNDIALS_EXPORT int MRIStepInnerStepper_SetEvolveFn(MRIStepInnerStepper stepper,
                                                    MRIStepInnerEvolveFn fn);
SUNDIALS_EXPORT int MRIStepInnerStepper_SetFullRhsFn(MRIStepInnerStepper stepper,
                                                     MRIStepInnerFullRhsFn fn);
SUNDIALS_EXPORT int MRIStepInnerStepper_SetResetFn(MRIStepInnerStepper stepper,
                                                   MRIStepInnerResetFn fn);
SUNDIALS_EXPORT int MRIStepInnerStepper_SetAccumulatedErrorGetFn(
  MRIStepInnerStepper stepper, MRIStepInnerGetAccumulatedError fn);
SUNDIALS_EXPORT int MRIStepInnerStepper_SetAccumulatedErrorResetFn(
  MRIStepInnerStepper stepper, MRIStepInnerResetAccumulatedError fn);
SUNDIALS_EXPORT int MRIStepInnerStepper_SetRTolFn(MRIStepInnerStepper stepper,
                                                  MRIStepInnerSetRTol fn);
SUNDIALS_EXPORT int MRIStepInnerStepper_AddForcing(MRIStepInnerStepper stepper,
                                                   sunrealtype t, N_Vector f);
SUNDIALS_EXPORT int MRIStepInnerStepper_GetForcingData(
  MRIStepInnerStepper stepper, sunrealtype* tshift, sunrealtype* tscale,
  N_Vector** forcing, int* nforcing);

#ifdef __cplusplus
}
#endif

#endif
