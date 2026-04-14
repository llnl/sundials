/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *---------------------------------------------------------------
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
 *---------------------------------------------------------------
 * Implementation header file for ARKODE's time step adaptivity
 * utilities.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_ADAPT_IMPL_H
#define _ARKODE_ADAPT_IMPL_H

#include <arkode/arkode.h>
#include <stdarg.h>
#include <sundials/sundials_adaptcontroller.h>

#include "arkode_types_impl.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  ARKODE Time Step Adaptivity Private Constants
  ===============================================================*/

/* size constants for the adaptivity memory structure */
#define ARK_ADAPT_LRW 10
#define ARK_ADAPT_LIW 7 /* includes function/data pointers */

/* Time step controller default values */
#define CFLFAC    SUN_RCONST(0.5)
#define SAFETY    SUN_RCONST(0.9)  /* CVODE uses 1.0  */
#define GROWTH    SUN_RCONST(20.0) /* CVODE uses 10.0 */
#define HFIXED_LB SUN_RCONST(1.0)  /* CVODE uses 1.0  */
#define HFIXED_UB SUN_RCONST(1.0)  /* CVODE uses 1.5  */

/* maximum step size change on first step */
#define ETAMX1 SUN_RCONST(10000.0)
/* step size reduction factor on multiple error test failures (multiple implies >= SMALL_NEF) */
#define ETAMXF SUN_RCONST(0.3)
/* smallest allowable step size reduction factor on an error test failure */
#define ETAMIN SUN_RCONST(0.1)
/* step size reduction factor on nonlinear convergence failure */
#define ETACF SUN_RCONST(0.25)
/* if an error failure occurs and SMALL_NEF <= nef, then reset eta = MIN(eta, ETAMXF) */
#define SMALL_NEF 2
/* order to use for controller:
     0=embedding,
     1=method,
     otherwise min(method,embedding)
   DEPRECATED, REMOVE AT SAME TIME AS ARKStepSetAdaptivityMethod */
#define PQ 0
/* adjustment to apply within controller to method order of accuracy */
#define ADJUST 0

/*===============================================================
  ARKODE Time Step Adaptivity Data Structure
  ===============================================================*/

/*---------------------------------------------------------------
  Types : struct ARKodeHAdaptMemRec, ARKodeHAdaptMem
  -----------------------------------------------------------------
  The type ARKodeHAdaptMem is type pointer to struct
  ARKodeHAdaptMemRec.  This structure contains fields to
  keep track of temporal adaptivity.
  ---------------------------------------------------------------*/
typedef struct ARKodeHAdaptMemRec
{
  sunrealtype etamax; /* eta <= etamax                              */
  sunrealtype etamx1; /* max step size change on first step         */
  sunrealtype etamxf; /* h reduction factor on multiple error fails */
  sunrealtype etamin; /* eta >= etamin on error test fail           */
  int small_nef;      /* bound to determine 'multiple' above        */
  sunrealtype etacf;  /* h reduction factor on nonlinear conv fail  */
  sunrealtype cfl;    /* cfl safety factor                          */
  sunrealtype safety; /* accuracy safety factor on h                */
  sunrealtype growth; /* maximum step growth safety factor          */
  sunrealtype lbound; /* eta lower bound to leave h unchanged       */
  sunrealtype ubound; /* eta upper bound to leave h unchanged       */
  int p;              /* embedding order                            */
  int q;              /* method order                               */
  int pq;             /* decision flag for controller order         */
  int adjust;         /* controller order adjustment factor         */

  SUNAdaptController hcontroller; /* temporal error controller            */
  sunbooleantype owncontroller;   /* flag indicating hcontroller ownership   */

  ARKExpStabFn expstab; /* step stability function                    */
  void* estab_data;     /* user pointer passed to expstab             */

  long int nst_acc; /* num accuracy-limited internal steps        */
  long int nst_exp; /* num stability-limited internal steps       */

}* ARKodeHAdaptMem;

/*===============================================================
  ARKODE Time Step Adaptivity Routines
  ===============================================================*/

ARKodeHAdaptMem arkAdaptInit(void);
void arkPrintAdaptMem(ARKodeHAdaptMem hadapt_mem, FILE* outfile);
int arkAdapt(ARKodeMem ark_mem, ARKodeHAdaptMem hadapt_mem, N_Vector ycur,
             sunrealtype tcur, sunrealtype hcur, sunrealtype dsm);

#ifdef __cplusplus
}
#endif

#endif
