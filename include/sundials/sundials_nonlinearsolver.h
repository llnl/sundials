/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner, and Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 * This is the header file for a generic nonlinear solver package. It defines
 * the SUNNonlinearSolver structure (_generic_SUNNonlinearSolver) which contains
 * the following fields:
 *   - an implementation-dependent 'content' field which contains any internal
 *     data required by the solver
 *   - an 'ops' filed which contains a structure listing operations acting on/by
 *     such solvers
 *
 * We consider iterative nonlinear solvers for systems in both root finding
 * (F(y) = 0) or fixed-point (G(y) = y) form. As a result, some of the routines
 * are applicable only to one type of nonlinear solver.
 * -----------------------------------------------------------------------------
 * This header file contains:
 *   - function types supplied to a SUNNonlinearSolver,
 *   - enumeration constants for SUNDIALS-defined nonlinear solver types,
 *   - type declarations for the _generic_SUNNonlinearSolver and
 *     _generic_SUNNonlinearSolver_Ops structures, as well as references to
 *     pointers to such structures (SUNNonlinearSolver),
 *   - prototypes for the nonlinear solver functions which operate
 *     on/by SUNNonlinearSolver objects, and
 *   - return codes for SUNNonLinearSolver objects.
 * -----------------------------------------------------------------------------
 * At a minimum, a particular implementation of a SUNNonlinearSolver must do the
 * following:
 *   - specify the 'content' field of a SUNNonlinearSolver,
 *   - implement the operations on/by the SUNNonlinearSovler objects,
 *   - provide a constructor routine for new SUNNonlinearSolver objects
 *
 * Additionally, a SUNNonlinearSolver implementation may provide the following:
 *   - "Set" routines to control solver-specific parameters/options
 *   - "Get" routines to access solver-specific performance metrics
 * ---------------------------------------------------------------------------*/

#ifndef _SUNNONLINEARSOLVER_H
#define _SUNNONLINEARSOLVER_H

#include <sundials/sundials_config.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------------------------------
 *  Forward references for SUNNonlinearSolver types defined below
 * ---------------------------------------------------------------------------*/

/* Forward reference for pointer to SUNNonlinearSolver_Ops object */
typedef _SUNDIALS_STRUCT_ _generic_SUNNonlinearSolver_Ops* SUNNonlinearSolver_Ops;

/* Forward reference for pointer to SUNNonlinearSolver object */
typedef _SUNDIALS_STRUCT_ _generic_SUNNonlinearSolver* SUNNonlinearSolver;

/* -----------------------------------------------------------------------------
 * Integrator supplied function types
 * ---------------------------------------------------------------------------*/

typedef int (*SUNNonlinSolSysFn)(N_Vector y, N_Vector F, void* mem);

typedef int (*SUNNonlinSolLSetupFn)(sunbooleantype jbad, sunbooleantype* jcur,
                                    void* mem);

typedef int (*SUNNonlinSolLSolveFn)(N_Vector b, void* mem);

typedef int (*SUNNonlinSolConvTestFn)(SUNNonlinearSolver NLS, N_Vector y,
                                      N_Vector del, sunrealtype tol,
                                      N_Vector ewt, void* mem);

/* -----------------------------------------------------------------------------
 * SUNNonlinearSolver types
 * ---------------------------------------------------------------------------*/

enum SUNNonlinearSolver_Type
{
  SUNNONLINEARSOLVER_ROOTFIND,
  SUNNONLINEARSOLVER_FIXEDPOINT
};

#ifndef SWIG
typedef enum SUNNonlinearSolver_Type SUNNonlinearSolver_Type;
#endif

/* -----------------------------------------------------------------------------
 * Generic definition of SUNNonlinearSolver
 * ---------------------------------------------------------------------------*/

/* Structure containing function pointers to nonlinear solver operations */
struct _generic_SUNNonlinearSolver_Ops
{
  SUNNonlinearSolver_Type (*gettype)(SUNNonlinearSolver);
  SUNErrCode (*initialize)(SUNNonlinearSolver);
  int (*setup)(SUNNonlinearSolver, N_Vector, void*);
  int (*solve)(SUNNonlinearSolver, N_Vector, N_Vector, N_Vector, sunrealtype,
               sunbooleantype, void*);
  SUNErrCode (*free)(SUNNonlinearSolver);
  SUNErrCode (*setsysfn)(SUNNonlinearSolver, SUNNonlinSolSysFn);
  SUNErrCode (*setlsetupfn)(SUNNonlinearSolver, SUNNonlinSolLSetupFn);
  SUNErrCode (*setlsolvefn)(SUNNonlinearSolver, SUNNonlinSolLSolveFn);
  SUNErrCode (*setctestfn)(SUNNonlinearSolver, SUNNonlinSolConvTestFn, void*);
  SUNErrCode (*setoptions)(SUNNonlinearSolver NLS, const char* NLSid,
                           const char* file_name, int argc, char* argv[]);
  SUNErrCode (*setmaxiters)(SUNNonlinearSolver, int);
  SUNErrCode (*getnumiters)(SUNNonlinearSolver, long int*);
  SUNErrCode (*getcuriter)(SUNNonlinearSolver, int*);
  SUNErrCode (*getnumconvfails)(SUNNonlinearSolver, long int*);
};

/* A nonlinear solver is a structure with an implementation-dependent 'content'
   field, and a pointer to a structure of solver nonlinear solver operations
   corresponding to that implementation. */
struct _generic_SUNNonlinearSolver
{
  void* content;
  void* python;
  SUNNonlinearSolver_Ops ops;
  SUNContext sunctx;
};

/* -----------------------------------------------------------------------------
 * Functions exported by SUNNonlinearSolver module
 * ---------------------------------------------------------------------------*/

/* empty constructor/destructor */
SUNDIALS_EXPORT
SUNNonlinearSolver SUNNonlinSolNewEmpty(SUNContext sunctx);

SUNDIALS_EXPORT
void SUNNonlinSolFreeEmpty(SUNNonlinearSolver NLS);

/* core functions */
SUNDIALS_EXPORT
SUNNonlinearSolver_Type SUNNonlinSolGetType(SUNNonlinearSolver NLS);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolInitialize(SUNNonlinearSolver NLS);

SUNDIALS_EXPORT
int SUNNonlinSolSetup(SUNNonlinearSolver NLS, N_Vector y, void* mem);

SUNDIALS_EXPORT
int SUNNonlinSolSolve(SUNNonlinearSolver NLS, N_Vector y0, N_Vector y, N_Vector w,
                      sunrealtype tol, sunbooleantype callLSetup, void* mem);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolFree(SUNNonlinearSolver NLS);

/* set functions */
SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetSysFn(SUNNonlinearSolver NLS, SUNNonlinSolSysFn SysFn);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetLSetupFn(SUNNonlinearSolver NLS,
                                   SUNNonlinSolLSetupFn SetupFn);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetLSolveFn(SUNNonlinearSolver NLS,
                                   SUNNonlinSolLSolveFn SolveFn);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetConvTestFn(SUNNonlinearSolver NLS,
                                     SUNNonlinSolConvTestFn CTestFn,
                                     void* ctest_data);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetOptions(SUNNonlinearSolver NLS, const char* NLSid,
                                  const char* file_name, int argc, char* argv[]);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetMaxIters(SUNNonlinearSolver NLS, int maxiters);

/* get functions */
SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetNumIters(SUNNonlinearSolver NLS, long int* niters);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetCurIter(SUNNonlinearSolver NLS, int* iter);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetNumConvFails(SUNNonlinearSolver NLS,
                                       long int* nconvfails);

/* -----------------------------------------------------------------------------
 * SUNNonlinearSolver return values
 * ---------------------------------------------------------------------------*/

/* Recoverable */
#define SUN_NLS_CONTINUE   +901 /* not converged, keep iterating      */
#define SUN_NLS_CONV_RECVR +902 /* convergece failure, try to recover */

#ifdef __cplusplus
}
#endif

#endif
