/* ----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *                Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
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
 * This is the header file for CVODE's linear solver interface.
 * ----------------------------------------------------------------*/

#ifndef _CVLS_H
#define _CVLS_H

#include <sundials/sundials_direct.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/*=================================================================
  CVLS Constants
  =================================================================*/

#define CVLS_SUCCESS         0
#define CVLS_MEM_NULL        -1
#define CVLS_LMEM_NULL       -2
#define CVLS_ILL_INPUT       -3
#define CVLS_MEM_FAIL        -4
#define CVLS_PMEM_NULL       -5
#define CVLS_JACFUNC_UNRECVR -6
#define CVLS_JACFUNC_RECVR   -7
#define CVLS_SUNMAT_FAIL     -8
#define CVLS_SUNLS_FAIL      -9

/*=================================================================
  CVLS user-supplied function prototypes
  =================================================================*/

typedef int (*CVLsJacFn)(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac,
                         void* user_data, N_Vector tmp1, N_Vector tmp2,
                         N_Vector tmp3);

typedef int (*CVLsPrecSetupFn)(sunrealtype t, N_Vector y, N_Vector fy,
                               sunbooleantype jok, sunbooleantype* jcurPtr,
                               sunrealtype gamma, void* user_data);

typedef int (*CVLsPrecSolveFn)(sunrealtype t, N_Vector y, N_Vector fy,
                               N_Vector r, N_Vector z, sunrealtype gamma,
                               sunrealtype delta, int lr, void* user_data);

typedef int (*CVLsJacTimesSetupFn)(sunrealtype t, N_Vector y, N_Vector fy,
                                   void* user_data);

typedef int (*CVLsJacTimesVecFn)(N_Vector v, N_Vector Jv, sunrealtype t,
                                 N_Vector y, N_Vector fy, void* user_data,
                                 N_Vector tmp);

typedef int (*CVLsLinSysFn)(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix A,
                            sunbooleantype jok, sunbooleantype* jcur,
                            sunrealtype gamma, void* user_data, N_Vector tmp1,
                            N_Vector tmp2, N_Vector tmp3);

/*=================================================================
  CVLS Exported functions
  =================================================================*/

SUNDIALS_EXPORT int CVodeSetLinearSolver(void* cvode_mem, SUNLinearSolver LS,
                                         SUNMatrix A);

/*-----------------------------------------------------------------
  Optional inputs to the CVLS linear solver interface
  -----------------------------------------------------------------*/

SUNDIALS_EXPORT int CVodeSetJacFn(void* cvode_mem, CVLsJacFn jac);
SUNDIALS_EXPORT int CVodeSetJacEvalFrequency(void* cvode_mem, long int msbj);
SUNDIALS_EXPORT int CVodeSetLinearSolutionScaling(void* cvode_mem,
                                                  sunbooleantype onoff);
SUNDIALS_EXPORT int CVodeSetDeltaGammaMaxBadJac(void* cvode_mem,
                                                sunrealtype dgmax_jbad);
SUNDIALS_EXPORT int CVodeSetEpsLin(void* cvode_mem, sunrealtype eplifac);
SUNDIALS_EXPORT int CVodeSetLSNormFactor(void* arkode_mem, sunrealtype nrmfac);
SUNDIALS_EXPORT int CVodeSetPreconditioner(void* cvode_mem, CVLsPrecSetupFn pset,
                                           CVLsPrecSolveFn psolve);
SUNDIALS_EXPORT int CVodeSetJacTimes(void* cvode_mem, CVLsJacTimesSetupFn jtsetup,
                                     CVLsJacTimesVecFn jtimes);
SUNDIALS_EXPORT int CVodeSetLinSysFn(void* cvode_mem, CVLsLinSysFn linsys);

/*-----------------------------------------------------------------
  Optional outputs from the CVLS linear solver interface
  -----------------------------------------------------------------*/

SUNDIALS_EXPORT int CVodeGetJac(void* cvode_mem,
                                SUNMatrix* J); // nb::rv_policy::reference
SUNDIALS_EXPORT int CVodeGetJacTime(void* cvode_mem, sunrealtype* t_J);
SUNDIALS_EXPORT int CVodeGetJacNumSteps(void* cvode_mem, long int* nst_J);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
int CVodeGetLinWorkSpace(void* cvode_mem, long int* lenrwLS, long int* leniwLS);
SUNDIALS_EXPORT int CVodeGetNumJacEvals(void* cvode_mem, long int* njevals);
SUNDIALS_EXPORT int CVodeGetNumPrecEvals(void* cvode_mem, long int* npevals);
SUNDIALS_EXPORT int CVodeGetNumPrecSolves(void* cvode_mem, long int* npsolves);
SUNDIALS_EXPORT int CVodeGetNumLinIters(void* cvode_mem, long int* nliters);
SUNDIALS_EXPORT int CVodeGetNumLinConvFails(void* cvode_mem, long int* nlcfails);
SUNDIALS_EXPORT int CVodeGetNumJTSetupEvals(void* cvode_mem, long int* njtsetups);
SUNDIALS_EXPORT int CVodeGetNumJtimesEvals(void* cvode_mem, long int* njvevals);
SUNDIALS_EXPORT int CVodeGetNumLinRhsEvals(void* cvode_mem, long int* nfevalsLS);
SUNDIALS_EXPORT int CVodeGetLinSolveStats(void* cvode_mem, long int* njevals,
                                          long int* nfevalsLS,
                                          long int* nliters, long int* nlcfails,
                                          long int* npevals, long int* npsolves,
                                          long int* njtsetups, long int* njtimes);
SUNDIALS_EXPORT int CVodeGetLastLinFlag(void* cvode_mem, long int* flag);
SUNDIALS_EXPORT char* CVodeGetLinReturnFlagName(long int flag);

#ifdef __cplusplus
}
#endif

#endif
