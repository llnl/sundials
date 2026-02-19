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
 * Implementation header file for ARKODE's temporal interpolation
 * utilities.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_INTERP_IMPL_H
#define _ARKODE_INTERP_IMPL_H

#include <arkode/arkode.h>
#include <stdarg.h>

#include "arkode_impl.h"

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  ARKODE temporal interpolation constants
  ===============================================================*/

/* Numeric constants */
#define FOURTH SUN_RCONST(0.25)
#define THREE  SUN_RCONST(3.0)
#define SIX    SUN_RCONST(6.0)
#define TWELVE SUN_RCONST(12.0)

/*===============================================================
  ARKODE Hermite Temporal Interpolation Data Structure
  ===============================================================*/

/* Hermite interpolation structure */

struct _ARKInterpContent_Hermite
{
  int degree;       /* maximum interpolant degree to use           */
  N_Vector fold;    /* f(t,y) at beginning of last successful step */
  N_Vector yold;    /* y at beginning of last successful step      */
  N_Vector fa;      /* f(t,y) used in higher-order interpolation   */
  N_Vector fb;      /* f(t,y) used in higher-order interpolation   */
  sunrealtype told; /* t at beginning of last successful step      */
  sunrealtype tnew; /* t at end of last successful step            */
  sunrealtype h;    /* last successful step size                   */
};

typedef struct _ARKInterpContent_Hermite* ARKInterpContent_Hermite;

/* Hermite structure accessor macros */

#define HINT_CONTENT(I) ((ARKInterpContent_Hermite)(I->content))
#define HINT_DEGREE(I)  (HINT_CONTENT(I)->degree)
#define HINT_FOLD(I)    (HINT_CONTENT(I)->fold)
#define HINT_YOLD(I)    (HINT_CONTENT(I)->yold)
#define HINT_FA(I)      (HINT_CONTENT(I)->fa)
#define HINT_FB(I)      (HINT_CONTENT(I)->fb)
#define HINT_TOLD(I)    (HINT_CONTENT(I)->told)
#define HINT_TNEW(I)    (HINT_CONTENT(I)->tnew)
#define HINT_H(I)       (HINT_CONTENT(I)->h)

/* Hermite structure operations */

ARKInterp arkInterpCreate_Hermite(ARKodeMem ark_mem, int degree);
int arkInterpResize_Hermite(ARKodeMem ark_mem, ARKInterp interp,
                            ARKVecResizeFn resize, void* resize_data,
                            sunindextype lrw_diff, sunindextype liw_diff,
                            N_Vector tmpl);
void arkInterpFree_Hermite(ARKodeMem ark_mem, ARKInterp interp);
void arkInterpPrintMem_Hermite(ARKInterp interp, FILE* outfile);
int arkInterpSetDegree_Hermite(ARKodeMem ark_mem, ARKInterp interp, int degree);
int arkInterpInit_Hermite(ARKodeMem ark_mem, ARKInterp interp, sunrealtype tnew);
int arkInterpUpdate_Hermite(ARKodeMem ark_mem, ARKInterp interp,
                            sunrealtype tnew);
int arkInterpEvaluate_Hermite(ARKodeMem ark_mem, ARKInterp interp,
                              sunrealtype tau, int d, int order, N_Vector yout);

/*===============================================================
  ARKODE Lagrange Temporal Interpolation Data Structure
  ===============================================================*/

/* Lagrange interpolation structure */

struct _ARKInterpContent_Lagrange
{
  int nmax;           /* number of previous solutions to use      */
  int nmaxalloc;      /* vectors allocated for previous solutions */
  N_Vector* yhist;    /* previous solution vectors                */
  sunrealtype* thist; /* 't' values associated with yhist         */
  int nhist;          /* number of 'active' vectors in yhist      */
  sunrealtype tround; /* unit roundoff for 't' values             */
};

typedef struct _ARKInterpContent_Lagrange* ARKInterpContent_Lagrange;

/* Lagrange structure accessor macros */

#define LINT_CONTENT(I)   ((ARKInterpContent_Lagrange)(I->content))
#define LINT_NMAX(I)      (LINT_CONTENT(I)->nmax)
#define LINT_NMAXALLOC(I) (LINT_CONTENT(I)->nmaxalloc)
#define LINT_YHIST(I)     (LINT_CONTENT(I)->yhist)
#define LINT_THIST(I)     (LINT_CONTENT(I)->thist)
#define LINT_YJ(I, j)     ((LINT_YHIST(I))[j])
#define LINT_TJ(I, j)     ((LINT_THIST(I))[j])
#define LINT_NHIST(I)     (LINT_CONTENT(I)->nhist)
#define LINT_TROUND(I)    (LINT_CONTENT(I)->tround)

/* Lagrange structure operations */

ARKInterp arkInterpCreate_Lagrange(ARKodeMem ark_mem, int degree);
int arkInterpResize_Lagrange(ARKodeMem ark_mem, ARKInterp interp,
                             ARKVecResizeFn resize, void* resize_data,
                             sunindextype lrw_diff, sunindextype liw_diff,
                             N_Vector tmpl);
void arkInterpFree_Lagrange(ARKodeMem ark_mem, ARKInterp interp);
void arkInterpPrintMem_Lagrange(ARKInterp interp, FILE* outfile);
int arkInterpSetDegree_Lagrange(ARKodeMem ark_mem, ARKInterp interp, int degree);
int arkInterpInit_Lagrange(ARKodeMem ark_mem, ARKInterp interp, sunrealtype tnew);
int arkInterpUpdate_Lagrange(ARKodeMem ark_mem, ARKInterp interp,
                             sunrealtype tnew);
int arkInterpEvaluate_Lagrange(ARKodeMem ark_mem, ARKInterp interp,
                               sunrealtype tau, int d, int order, N_Vector yout);

/* Lagrange structure utility routines */
sunrealtype LBasis(ARKInterp interp, int idx, sunrealtype t);
sunrealtype LBasisD(ARKInterp interp, int idx, sunrealtype t);
sunrealtype LBasisD2(ARKInterp interp, int idx, sunrealtype t);
sunrealtype LBasisD3(ARKInterp interp, int idx, sunrealtype t);

#ifdef __cplusplus
}
#endif

#endif
