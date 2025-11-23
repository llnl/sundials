/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *---------------------------------------------------------------
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
 *---------------------------------------------------------------
 * Implementation header file for ARKODE's root-finding (in time)
 * utility.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_ROOT_IMPL_H
#define _ARKODE_ROOT_IMPL_H

#include <arkode/arkode.h>
#include <stdarg.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  ARKODE Root-finding constants
  ===============================================================*/

#define ARK_ROOT_LRW 5
#define ARK_ROOT_LIW 12

/* Numeric constants */
#define HUND SUN_RCONST(100.0)

/*===============================================================
  ARKODE Root-finding Data Structure
  ===============================================================*/

/*---------------------------------------------------------------
  Types : struct ARKodeRootMemRec, ARKodeRootMem
  -----------------------------------------------------------------
  The type ARKodeRootMem is type pointer to struct
  ARKodeRootMemRec.  This structure contains data pertaining to
  the use of root-finding capabilities in ARKODE.
  ---------------------------------------------------------------*/
typedef struct ARKodeRootMemRec
{
  ARKRootFn gfun;          /* function g for roots sought                  */
  int nrtfn;               /* number of components of g                    */
  int* iroots;             /* array for root information                   */
  int* rootdir;            /* array specifying direction of zero-crossing  */
  sunrealtype tlo;         /* nearest endpoint of interval in root search  */
  sunrealtype thi;         /* farthest endpoint of interval in root search */
  sunrealtype trout;       /* t value returned by rootfinding routine      */
  sunrealtype* glo;        /* saved array of g values at t = tlo           */
  sunrealtype* ghi;        /* saved array of g values at t = thi           */
  sunrealtype* grout;      /* array of g values at t = trout               */
  sunrealtype ttol;        /* tolerance on root location                   */
  int irfnd;               /* flag showing whether last step had a root    */
  long int nge;            /* counter for g evaluations                    */
  sunbooleantype* gactive; /* array with active/inactive event functions   */
  int mxgnull;             /* num. warning messages about possible g==0    */
  void* root_data;         /* pointer to user_data                         */

}* ARKodeRootMem;

/*===============================================================
  ARKODE Root-finding Routines
===============================================================*/

int arkRootFree(void* arkode_mem);
int arkPrintRootMem(void* arkode_mem, FILE* outfile);
int arkRootCheck1(void* arkode_mem);
int arkRootCheck2(void* arkode_mem);
int arkRootCheck3(void* arkode_mem, sunrealtype tout, int itask);
int arkRootfind(void* arkode_mem);

#ifdef __cplusplus
}
#endif

#endif
