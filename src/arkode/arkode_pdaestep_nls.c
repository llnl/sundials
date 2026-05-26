// /*---------------------------------------------------------------
//  * Programmer(s): Steven B. Roberts @ LLNL
//  *---------------------------------------------------------------
//  * SUNDIALS Copyright Start
//  * Copyright (c) 2025-2026, Lawrence Livermore National Security,
//  * University of Maryland Baltimore County, and the SUNDIALS contributors.
//  * Copyright (c) 2013-2025, Lawrence Livermore National Security
//  * and Southern Methodist University.
//  * Copyright (c) 2002-2013, Lawrence Livermore National Security.
//  * All rights reserved.
//  *
//  * See the top-level LICENSE and NOTICE files for details.
//  *
//  * SPDX-License-Identifier: BSD-3-Clause
//  * SUNDIALS Copyright End
//  *---------------------------------------------------------------
//  * This is the interface between ARKStep and the
//  * SUNNonlinearSolver object
//  *--------------------------------------------------------------*/

// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <sundials/sundials_math.h>
// #include <sundials/sundials_nonlinearsolver.h>

// #include "arkode_pdaestep_impl.h"
// #include "arkode_impl.h"

// static int pdaeStep_NlsResidual(N_Vector wcor, N_Vector r, void* arkode_mem)
// {
//   /* temporary variables */
//   ARKodeMem ark_mem;
//   ARKodePDAEStepMem step_mem;
//   int retval;

//   /* access ARKodeMem and ARKodeARKStepMem structures */
//   retval = arkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
//   if (retval != ARK_SUCCESS) { return (retval); }
  
//   /* update 'ycur' value as stored predictor + current corrector */
//   // TODO(SBR): ycur has wrong dimension
//   N_VLinearSum(ONE, step_mem->wpred, ONE, wcor, ark_mem->ycur);

//   // TODO(SBR): pass correct args
//   retval = step_mem->hFn(ark_mem->tcur, r,
//                             ???, ark_mem->user_data);

//   if (retval < 0) { return (ARK_RHSFUNC_FAIL); }
//   if (retval > 0) { return (RHSFUNC_RECVR); }
//   step_mem->nhi++;

//   return ARK_SUCCESS;
// }

// static int pdaeStep_SetNlsSysFn(ARKodeMem ark_mem)
// {
//   ARKodePDAEStepMem step_mem;
//   int retval;

//   /* access ARKodePDAEStepMem structure */
//   retval = arkStep_AccessStepMem(ark_mem, __func__, &step_mem);
//   if (retval != ARK_SUCCESS) { return (retval); }

//   retval = SUNNonlinSolSetSysFn(step_mem->NLS, pdaeStep_NlsResidual);
//   if (retval != ARK_SUCCESS)
//   {
//     arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
//                     "Setting nonlinear system function failed");
//     return (ARK_ILL_INPUT);
//   }

//   return ARK_SUCCESS;
// }

// static int pdaeStep_NlsConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del,
//   sunrealtype tol, N_Vector ewt, void* arkode_mem)
// {
//   /* temporary variables */
//   ARKodeMem ark_mem;
//   ARKodePDAEStepMem step_mem;
//   sunrealtype delnrm, dcon;
//   int m, retval;

//   /* access ARKodeMem and ARKodePDAEStepMem structures */
//   retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
//   if (retval != ARK_SUCCESS) { return (retval); }

//   /* compute the norm of the correction */
//   delnrm = N_VWrmsNorm(del, ewt);

//   /* get the current nonlinear solver iteration count */
//   retval = SUNNonlinSolGetCurIter(NLS, &m);
//   if (retval != ARK_SUCCESS) { return (ARK_MEM_NULL); }

//   /* update the stored estimate of the convergence rate (assumes linear convergence) */
//   if (m > 0)
//   {
//     step_mem->crate = SUNMAX(step_mem->crdown * step_mem->crate,
//                              delnrm / step_mem->delp);
//   }

//   /* compute our scaled error norm for testing convergence */
//   dcon = SUNMIN(step_mem->crate, ONE) * delnrm / tol;

//   /* check for convergence; if so return with success */
//   if (dcon <= ONE) { return (SUN_SUCCESS); }

//   /* check for divergence */
//   if ((m >= 1) && (delnrm > step_mem->rdiv * step_mem->delp))
//   {
//     return (SUN_NLS_CONV_RECVR);
//   }

//   /* save norm of correction for next iteration */
//   step_mem->delp = delnrm;

//   /* return with flag that there is more work to do */
//   return (SUN_NLS_CONTINUE);
// }

// /*------------------------------------------------------------------------------
//   TODO(SBR): Description
//   ----------------------------------------------------------------------------*/
// int pdaeStep_SetNonlinearSolver(ARKodeMem ark_mem, SUNNonlinearSolver NLS)
// {
//   ARKodePDAEStepMem step_mem = NULL;
//   int retval = pdaeStep_AccessStepMem(ark_mem, __func__, &step_mem);
//   if (retval != ARK_SUCCESS) { return retval; }

//   /* Return immediately if NLS input is NULL */
//   if (NLS == NULL)
//   {
//     arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
//                     "The NLS input must be non-NULL");
//     return (ARK_ILL_INPUT);
//   }

//   /* check for required nonlinear solver functions */
//   if ((NLS->ops->gettype == NULL) || (NLS->ops->solve == NULL) ||
//       (NLS->ops->setsysfn == NULL))
//   {
//     arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
//                     "NLS does not support required operations");
//     return (ARK_ILL_INPUT);
//   }

//   /* free any existing nonlinear solver */
//   if (step_mem->NLS != NULL && step_mem->ownNLS)
//   {
//     retval = SUNNonlinSolFree(step_mem->NLS);
//     if (retval != ARK_SUCCESS) { return retval; }
//   }

//   /* set SUNNonlinearSolver pointer */
//   step_mem->NLS    = NLS;
//   step_mem->ownNLS = SUNFALSE;

//   retval = SUNNonlinSolSetConvTestFn(step_mem->NLS, pdaeStep_NlsConvTest,
//                                      (void*)ark_mem);
//   if (retval != ARK_SUCCESS)
//   {
//     arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
//                     "Setting convergence test function failed");
//     return (ARK_ILL_INPUT);
//   }

//   /* set default nonlinear iterations */
//   retval = SUNNonlinSolSetMaxIters(step_mem->NLS, step_mem->maxcor);
//   if (retval != ARK_SUCCESS)
//   {
//     arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
//                     "Setting maximum number of nonlinear iterations failed");
//     return (ARK_ILL_INPUT);
//   }

//   return ARK_SUCCESS;
// }

// static int pdaeStep_NlsLSetup(sunbooleantype jbad, sunbooleantype* jcur, void* arkode_mem)
// {
//   ARKodeMem ark_mem;
//   ARKodePDAEStepMem step_mem;
//   int retval;

//   /* access ARKodeMem and ARKodeARKStepMem structures */
//   retval = arkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
//   if (retval != ARK_SUCCESS) { return (retval); }

//   /* update convfail based on jbad flag */
//   if (jbad) { step_mem->convfail = ARK_FAIL_BAD_J; }

//   /* Use ARKODE's tempv1, tempv2 and tempv3 as
//      temporary vectors for the linear solver setup routine */
//   step_mem->nsetups++;
//   retval = step_mem->lsetup(ark_mem, step_mem->convfail, ark_mem->tcur,
//                             ???, ???, &(step_mem->jcur), ???, ???, ???);

//   /* update Jacobian status */
//   *jcur = step_mem->jcur;

//   step_mem->crate = ONE;
//   step_mem->nstlp                    = ark_mem->nst;

//   if (retval < 0) { return (ARK_LSETUP_FAIL); }
//   if (retval > 0) { return (CONV_FAIL); }

//   return (ARK_SUCCESS);
// }

// static int pdaeStep_NlsLSolve(N_Vector b, void* arkode_mem)
// {
//   ARKodeMem ark_mem;
//   ARKodePDAEStepMem step_mem;
//   int retval, nonlin_iter;

//   /* access ARKodeMem and ARKodeARKStepMem structures */
//   retval = arkStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem, &step_mem);
//   if (retval != ARK_SUCCESS) { return (retval); }

//   /* retrieve nonlinear solver iteration from module */
//   retval = SUNNonlinSolGetCurIter(step_mem->NLS, &nonlin_iter);
//   if (retval != SUN_SUCCESS) { return (ARK_NLS_OP_ERR); }

//   /* call linear solver interface, and handle return value */
//   retval = step_mem->lsolve(ark_mem, b, ark_mem->tcur, ???,
//                             ???, step_mem->eRNrm,
//                             nonlin_iter);

//   if (retval < 0) { return (ARK_LSOLVE_FAIL); }
//   if (retval > 0) { return (CONV_FAIL); }

//   return (ARK_SUCCESS);
// }

// int pdaeStep_NlsInit(ARKodeMem ark_mem)
// {
//   ARKodePDAEStepMem step_mem;
//   int retval;

//   /* access ARKodePDAEStepMem structure */
//   if (ark_mem->step_mem == NULL)
//   {
//     arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
//                     MSG_ARKSTEP_NO_MEM);
//     return (ARK_MEM_NULL);
//   }
//   step_mem = (ARKodePDAEStepMem)ark_mem->step_mem;

//   /* reset counters */
//   step_mem->nls_iters = 0;
//   step_mem->nls_fails = 0;

//   /* set the linear solver setup wrapper function */
//   if (step_mem->lsetup)
//   {
//     retval = SUNNonlinSolSetLSetupFn(step_mem->NLS, pdaeStep_NlsLSetup);
//   }
//   else { retval = SUNNonlinSolSetLSetupFn(step_mem->NLS, NULL); }
//   if (retval != ARK_SUCCESS)
//   {
//     arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
//                     "Setting the linear solver setup function failed");
//     return (ARK_NLS_INIT_FAIL);
//   }

//   /* set the linear solver solve wrapper function */
//   if (step_mem->lsolve)
//   {
//     retval = SUNNonlinSolSetLSolveFn(step_mem->NLS, pdaeStep_NlsLSolve);
//   }
//   else { retval = SUNNonlinSolSetLSolveFn(step_mem->NLS, NULL); }
//   if (retval != ARK_SUCCESS)
//   {
//     arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
//                     "Setting linear solver solve function failed");
//     return (ARK_NLS_INIT_FAIL);
//   }

//   retval = arkStep_SetNlsSysFn(ark_mem);
//   if (retval != ARK_SUCCESS)
//   {
//     arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
//                     "Setting nonlinear system function failed");
//     return (ARK_ILL_INPUT);
//   }

//   /* initialize nonlinear solver */
//   retval = SUNNonlinSolInitialize(step_mem->NLS);
//   if (retval != ARK_SUCCESS)
//   {
//     arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
//                     MSG_NLS_INIT_FAIL);
//     return (ARK_NLS_INIT_FAIL);
//   }

//   return (ARK_SUCCESS);
// }

// int pdaeStep_Nls(ARKodeMem ark_mem, int nflag)
// {
//   ARKodePDAEStepMem step_mem;
//   sunbooleantype callLSetup;
//   long int nls_iters_inc = 0;
//   long int nls_fails_inc = 0;
//   int retval;

//   /* access ARKodeARKStepMem structure */
//   if (ark_mem->step_mem == NULL)
//   {
//     arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
//                     MSG_ARKSTEP_NO_MEM);
//     return (ARK_MEM_NULL);
//   }
//   step_mem = (ARKodePDAEStepMem)ark_mem->step_mem;

//   /* If a linear solver 'setup' is supplied, set various flags for
//      determining whether it should be called */
//   if (step_mem->lsetup)
//   {
//     /* Set interface 'convfail' flag for use inside lsetup */
//     step_mem->convfail = ((nflag == FIRST_CALL) || (nflag == PREV_ERR_FAIL))
//                              ? ARK_NO_FAILURES
//                              : ARK_FAIL_OTHER;

//     /* Decide whether to recommend call to lsetup within nonlinear solver */
//     //TODO(SBR): figure out if firststage is needed
//     callLSetup = (ark_mem->firststage) || (step_mem->msbp < 0);
//     callLSetup = callLSetup || (nflag == PREV_CONV_FAIL) ||
//                    (nflag == PREV_ERR_FAIL) ||
//                    (ark_mem->nst >= step_mem->nstlp + abs(step_mem->msbp));
//   }
//   else
//   {
//     step_mem->crate = ONE;
//     callLSetup      = SUNFALSE;
//   }

//   /* set a zero guess for correction */
//   // TODO: predictor here or elsewhere?
//   N_VConst(ZERO, step_mem->wcor);

//   /* Reset the stored residual norm (for iterative linear solvers) */
//   step_mem->eRNrm = step_mem->nlscoef;

//   /* solve the nonlinear system for the actual correction */
//   retval = SUNNonlinSolSolve(step_mem->NLS, step_mem->wpred, step_mem->wcor,
//                              ark_mem->ewt, step_mem->nlscoef, callLSetup,
//                              ark_mem);

//   /* increment counters */
//   SUNNonlinSolGetNumIters(step_mem->NLS, &nls_iters_inc);
//   step_mem->nls_iters += nls_iters_inc;

//   SUNNonlinSolGetNumConvFails(step_mem->NLS, &nls_fails_inc);
//   step_mem->nls_fails += nls_fails_inc;

//   /* successful solve -- reset jcur flag and apply correction */
//   if (retval == SUN_SUCCESS)
//   {
//     step_mem->jcur = SUNFALSE;
//     N_VLinearSum(ONE, step_mem->wcor, ONE, step_mem->wpred, ark_mem->ycur);

//     SUNLogInfo(ARK_LOGGER, "end-nonlinear-solve",
//                "status = success, iters = %li", nls_iters_inc);

//     return (ARK_SUCCESS);
//   }

//   SUNLogInfo(ARK_LOGGER, "end-nonlinear-solve",
//              "status = failed, retval = %i, iters = %li", retval, nls_iters_inc);

//   /* check for recoverable failure, return ARKODE::CONV_FAIL */
//   if (retval == SUN_NLS_CONV_RECVR) { return (CONV_FAIL); }

//   return (retval);
// }
