/* ----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 * ----------------------------------------------------------------
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
 * ----------------------------------------------------------------
 * This is the header file for ARKode's linear solver interface.
 * ----------------------------------------------------------------*/

#ifndef _ARKLS_H
#define _ARKLS_H

#include <sundials/sundials_direct.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*=================================================================
  ARKLS Constants
  =================================================================*/

#define ARKLS_SUCCESS          0
#define ARKLS_MEM_NULL         -1
#define ARKLS_LMEM_NULL        -2
#define ARKLS_ILL_INPUT        -3
#define ARKLS_MEM_FAIL         -4
#define ARKLS_PMEM_NULL        -5
#define ARKLS_MASSMEM_NULL     -6
#define ARKLS_JACFUNC_UNRECVR  -7
#define ARKLS_JACFUNC_RECVR    -8
#define ARKLS_MASSFUNC_UNRECVR -9
#define ARKLS_MASSFUNC_RECVR   -10
#define ARKLS_SUNMAT_FAIL      -11
#define ARKLS_SUNLS_FAIL       -12

/*=================================================================
  ARKLS user-supplied function prototypes
  =================================================================*/

typedef int (*ARKLsJacFn)(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac,
                          void* user_data, N_Vector tmp1, N_Vector tmp2,
                          N_Vector tmp3);

typedef int (*ARKLsMassFn)(sunrealtype t, SUNMatrix M, void* user_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

typedef int (*ARKLsPrecSetupFn)(sunrealtype t, N_Vector y, N_Vector fy,
                                sunbooleantype jok, sunbooleantype* jcurPtr,
                                sunrealtype gamma, void* user_data);

typedef int (*ARKLsPrecSolveFn)(sunrealtype t, N_Vector y, N_Vector fy,
                                N_Vector r, N_Vector z, sunrealtype gamma,
                                sunrealtype delta, int lr, void* user_data);

typedef int (*ARKLsJacTimesSetupFn)(sunrealtype t, N_Vector y, N_Vector fy,
                                    void* user_data);

typedef int (*ARKLsJacTimesVecFn)(N_Vector v, N_Vector Jv, sunrealtype t,
                                  N_Vector y, N_Vector fy, void* user_data,
                                  N_Vector tmp);

typedef int (*ARKLsLinSysFn)(sunrealtype t, N_Vector y, N_Vector fy,
                             SUNMatrix A, SUNMatrix M, sunbooleantype jok,
                             sunbooleantype* jcur, sunrealtype gamma,
                             void* user_data, N_Vector tmp1, N_Vector tmp2,
                             N_Vector tmp3);

typedef int (*ARKLsMassTimesSetupFn)(sunrealtype t, void* mtimes_data);

typedef int (*ARKLsMassTimesVecFn)(N_Vector v, N_Vector Mv, sunrealtype t,
                                   void* mtimes_data);

typedef int (*ARKLsMassPrecSetupFn)(sunrealtype t, void* user_data);

typedef int (*ARKLsMassPrecSolveFn)(sunrealtype t, N_Vector r, N_Vector z,
                                    sunrealtype delta, int lr, void* user_data);

/* Linear solver set functions */
SUNDIALS_EXPORT int ARKodeSetLinearSolver(void* arkode_mem, SUNLinearSolver LS,
                                          SUNMatrix A);
SUNDIALS_EXPORT int ARKodeSetMassLinearSolver(void* arkode_mem,
                                              SUNLinearSolver LS, SUNMatrix M,
                                              sunbooleantype time_dep);

/* Linear solver interface optional input functions -- must be called
   AFTER ARKodeSetLinearSolver and/or ARKodeSetMassLinearSolver */
SUNDIALS_EXPORT int ARKodeSetJacFn(void* arkode_mem, ARKLsJacFn jac);
SUNDIALS_EXPORT int ARKodeSetMassFn(void* arkode_mem, ARKLsMassFn mass);
SUNDIALS_EXPORT int ARKodeSetJacEvalFrequency(void* arkode_mem, long int msbj);
SUNDIALS_EXPORT int ARKodeSetLinearSolutionScaling(void* arkode_mem,
                                                   sunbooleantype onoff);
SUNDIALS_EXPORT int ARKodeSetEpsLin(void* arkode_mem, sunrealtype eplifac);
SUNDIALS_EXPORT int ARKodeSetMassEpsLin(void* arkode_mem, sunrealtype eplifac);
SUNDIALS_EXPORT int ARKodeSetLSNormFactor(void* arkode_mem, sunrealtype nrmfac);
SUNDIALS_EXPORT int ARKodeSetMassLSNormFactor(void* arkode_mem,
                                              sunrealtype nrmfac);
SUNDIALS_EXPORT int ARKodeSetPreconditioner(void* arkode_mem,
                                            ARKLsPrecSetupFn psetup,
                                            ARKLsPrecSolveFn psolve);
SUNDIALS_EXPORT int ARKodeSetMassPreconditioner(void* arkode_mem,
                                                ARKLsMassPrecSetupFn psetup,
                                                ARKLsMassPrecSolveFn psolve);
SUNDIALS_EXPORT int ARKodeSetJacTimes(void* arkode_mem,
                                      ARKLsJacTimesSetupFn jtsetup,
                                      ARKLsJacTimesVecFn jtimes);
SUNDIALS_EXPORT int ARKodeSetJacTimesRhsFn(void* arkode_mem,
                                           ARKRhsFn jtimesRhsFn);
SUNDIALS_EXPORT int ARKodeSetMassTimes(void* arkode_mem,
                                       ARKLsMassTimesSetupFn msetup,
                                       ARKLsMassTimesVecFn mtimes,
                                       void* mtimes_data);
SUNDIALS_EXPORT int ARKodeSetLinSysFn(void* arkode_mem, ARKLsLinSysFn linsys);

#ifdef __cplusplus
}
#endif

#endif
