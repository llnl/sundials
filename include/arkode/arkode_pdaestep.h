/*---------------------------------------------------------------
 * Programmer(s): Steven B. Roberts @ LLNL
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
 * This is the header file for the ARKODE PDAEStep module.
 *--------------------------------------------------------------*/

#ifndef ARKODE_PDAESTEP_H_
#define ARKODE_PDAESTEP_H_

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_nonlinearsolver.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_stepper.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* y must be a manyvector since we need to update the subvector z independent of
 * x. We use single residual function rather than separate RHS functions for
 * diff/alg so that the corresponding Jacobian function can create the entire
 * Jacobian. Otherwise, we would have to merge the diff/alg blocks of the
 * Jacobian which is not a currently supported operation and would probably be
 * less efficient.
 */
/* TODO(SBR): How do we split this for MPI?
 * - Multiple processes per partition or multiple partitions per process?
 * - Is w split across ranks?
 */
typedef int (*PDAEStepComponentResFn)(int partition, sunrealtype t, N_Vector y,
                                      N_Vector w, N_Vector yp, N_Vector res,
                                      void* user_data); // [y_d' - f(y,w); g(y,w)]

// TODO(SBR): Could add tmp vectors and f/res if needed
typedef int (*PDAEStepLsComponentJacFn)(int partition, sunrealtype t,
                                        sunrealtype c_j, N_Vector y, N_Vector w,
                                        N_Vector yp, N_Vector r, SUNMatrix Jac,
                                        void* user_data);

typedef int (*PDAEStepLsComponentJacTimesVecFn)(int partition, sunrealtype t,
                                                N_Vector y, N_Vector w,
                                                N_Vector yp, N_Vector r,
                                                N_Vector v, N_Vector Jv,
                                                sunrealtype c_j, void* user_data,
                                                N_Vector tmp1, N_Vector tmp2);

/* y is a manyvector so the residual can use all z^{r}. res is also a manyvector
 * This include h and g functions so we can defer the exploiting of structure to
 * the (non)linear solver
 */
typedef int (*PDAEStepAlgebraicResFn)(sunrealtype t, N_Vector y, N_Vector w,
                                      N_Vector res, void* user_data);

// TODO(SBR): Could add tmp vectors and f/res if needed
typedef int (*PDAEStepLsAlgebraicJacFn)(sunrealtype t, N_Vector y, N_Vector w,
                                        SUNMatrix Jac, void* user_data);

typedef int (*PDAEStepLsAlgebraicJacTimesVecFn)(sunrealtype t, N_Vector y,
                                                N_Vector w, N_Vector v,
                                                N_Vector Jv, void* user_data,
                                                N_Vector tmp);

// y0 and yp0 are manyvectors with 2*n_partitions+1 subvectors
SUNDIALS_EXPORT void* PDAEStepCreate(PDAEStepComponentResFn component_res_fn,
                                     PDAEStepAlgebraicResFn algebraic_res_fn,
                                     sunrealtype t0, N_Vector y0, N_Vector yp0,
                                     int n_partitions, SUNContext sunctx);

SUNDIALS_EXPORT int PDAEStepGetNumPartitions(void* arkode_mem, int* partitions);

SUNDIALS_EXPORT int PDAEStepGetPartitionVectorTemplate(void* arkode_mem,
                                                       int partition,
                                                       N_Vector* y);

SUNDIALS_EXPORT int PDAEStepGetPartitionIntegrator(void* arkode_mem,
                                                   int partition, void** ida_mem);

SUNDIALS_EXPORT int PDAEStepSetPartitionJacobian(void* arkode_mem,
                                                 PDAEStepLsComponentJacFn jac);

SUNDIALS_EXPORT int PDAEStepSetPartitionJacTimes(
  void* arkode_mem, PDAEStepLsComponentJacTimesVecFn jtimes);

SUNDIALS_EXPORT int PDAEStepSetCouplingJacobian(void* arkode_mem,
                                                PDAEStepLsAlgebraicJacFn jac);

SUNDIALS_EXPORT int PDAEStepSetCouplingJacTimes(
  void* arkode_mem, PDAEStepLsAlgebraicJacTimesVecFn jtimes);

/* We provide PDAE-specific (non)linear solver functions because there are
 * customizations we need to use that are not easily supported by the default
 * ARKODE infrastructure. First, the dimension of the nonlinear system is the
 * total number of algebraic variables rather than the total number of
 * variables. Second, the system is purely algebraic or as if there is a 0 mass
 * matrix.
 */
SUNDIALS_EXPORT int PDAEStepSetNonlinearSolver(void* arkode_mem,
                                               SUNNonlinearSolver nls);

SUNDIALS_EXPORT int PDAEStepSetLinearSolver(void* arkode_mem,
                                            SUNLinearSolver ls, SUNMatrix j);

SUNDIALS_EXPORT int PDAEStepSetMaxNonlinIters(void* arkode_mem, int maxcor);

SUNDIALS_EXPORT int PDAEStepSetNonlinConvCoef(void* arkode_mem,
                                              sunrealtype nlscoef);

SUNDIALS_EXPORT int PDAEStepSetNonlinCRDown(void* arkode_mem, sunrealtype crdown);

SUNDIALS_EXPORT int PDAEStepSetNonlinRDiv(void* arkode_mem, sunrealtype rdiv);

SUNDIALS_EXPORT int PDAEStepGetNumLinSolvSetups(void* arkode_mem,
                                                long int* nlinsetups);

SUNDIALS_EXPORT int PDAEStepGetNumNonlinSolvIters(void* arkode_mem,
                                                  long int* nniters);

SUNDIALS_EXPORT int PDAEStepGetNumNonlinSolvConvFails(void* arkode_mem,
                                                      long int* nnfails);

SUNDIALS_EXPORT int PDAEStepGetNonlinSolvStats(void* arkode_mem,
                                               long int* nniters,
                                               long int* nnfails);

SUNDIALS_EXPORT int PDAEStepGetAlgebraicVectorTemplate(void* arkode_mem,
                                                       N_Vector* y);

SUNDIALS_EXPORT N_Vector PDAEStepManyVector(N_Vector* x, N_Vector* z,
                                            N_Vector w, int partitions);
SUNDIALS_EXPORT N_Vector PDAEStepGetCouplingSubvector(N_Vector y);
SUNDIALS_EXPORT N_Vector PDAEStepGetDifferentialSubvector(N_Vector y,
                                                          int partition);
SUNDIALS_EXPORT N_Vector PDAEStepGetAlgebraicSubvector(N_Vector y, int partition);

/* TODO(SBR): For setting tolerances, we probably want ARKodeSSTolerances to
 * only specify tolerances for the nonlinear solver and not the inner IDA
 * instances. This will allow separate tuning of those, possible by partition.
 */

// TODO(SBR): Add reinit function? Probably not needed

#ifdef __cplusplus
}
#endif

#endif
