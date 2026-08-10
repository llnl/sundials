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
 * This is the interface between PDAEStep and SUNNonlinearSolver.
 *--------------------------------------------------------------*/

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nonlinearsolver.h>

#include "arkode_impl.h"
#include "arkode_pdaestep_impl.h"
#include "sundials_utils.h"

/* private functions */
static int pdaeStep_NlsResidual(N_Vector alg_cor, N_Vector r, void* arkode_mem);
static int pdaeStep_NlsLSetup(sunbooleantype jbad, sunbooleantype* jcur,
                              void* arkode_mem);
static int pdaeStep_NlsLSolve(N_Vector b, void* arkode_mem);
static int pdaeStep_AlgJtimes(void* arkode_mem, N_Vector v, N_Vector Jv);
static int pdaeStep_AlgDQJtimes(void* arkode_mem, N_Vector v, N_Vector Jv);
static int pdaeStep_AttachNlsCallbacks(ARKodeMem ark_mem,
                                       ARKodePDAEStepMem step_mem);
static SUNErrCode pdaeStep_NlsNorm(N_Vector del, N_Vector ewt,
                                   sunrealtype* delnrm, void* arkode_mem);
static SUNErrCode pdaeStep_NlsGetUpdateNorm(sunrealtype* delnrm,
                                            void* arkode_mem);

/*---------------------------------------------------------------
  PDAEStepSetNonlinearSolver:

  Attaches a root-finding SUNNonlinearSolver to the PDAEStep
  algebraic solve.
  ---------------------------------------------------------------*/
int PDAEStepSetNonlinearSolver(void* arkode_mem, SUNNonlinearSolver NLS)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (NLS == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The NLS input must be non-NULL");
    return ARK_ILL_INPUT;
  }

  if ((NLS->ops == NULL) || (NLS->ops->gettype == NULL) ||
      (NLS->ops->solve == NULL))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "NLS does not support required root-finding operations");
    return ARK_ILL_INPUT;
  }

  if (SUNNonlinSolGetType(NLS) != SUNNONLINEARSOLVER_ROOTFIND)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "PDAEStep currently supports root-finding nonlinear "
                    "solvers only");
    return ARK_ILL_INPUT;
  }

  if (NLS->ops->setsysfn == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "NLS does not support the required system function "
                    "operation");
    return ARK_ILL_INPUT;
  }

  if ((step_mem->NLS != NULL) && (step_mem->ownNLS))
  {
    retval = SUNNonlinSolFree(step_mem->NLS);
    if (retval != SUN_SUCCESS) { return ARK_NLS_OP_ERR; }
  }

  step_mem->NLS    = NLS;
  step_mem->ownNLS = SUNFALSE;

  return pdaeStep_AttachNlsCallbacks(ark_mem, step_mem);
}

/*---------------------------------------------------------------
  PDAEStepSetLinearSolver:

  Attaches a SUNLinearSolver to the PDAEStep Newton solve.
  ---------------------------------------------------------------*/
int PDAEStepSetLinearSolver(void* arkode_mem, SUNLinearSolver LS, SUNMatrix J)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (LS == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The LS input must be non-NULL");
    return ARK_ILL_INPUT;
  }

  if ((LS->ops == NULL) || (LS->ops->gettype == NULL) || (LS->ops->solve == NULL))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "LS does not support required operations");
    return ARK_ILL_INPUT;
  }

  if ((J == NULL) && (LS->ops->setatimes == NULL))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "A matrix-free PDAEStep linear solver must support "
                    "SUNLinSolSetATimes");
    return ARK_ILL_INPUT;
  }

  if ((J == NULL) && (SUNLinSolGetType(LS) == SUNLINEARSOLVER_DIRECT))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "A direct PDAEStep linear solver requires a SUNMatrix");
    return ARK_ILL_INPUT;
  }

  step_mem->LS = LS;
  step_mem->J  = J;

  if (LS->ops->setatimes)
  {
    retval = SUNLinSolSetATimes(LS, ark_mem,
                                (step_mem->algebraic_res_jtimes == NULL)
                                  ? pdaeStep_AlgDQJtimes
                                  : pdaeStep_AlgJtimes);
    if (retval != SUN_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Error in calling SUNLinSolSetATimes");
      return ARK_ILL_INPUT;
    }
  }

  if (LS->ops->setpreconditioner)
  {
    retval = SUNLinSolSetPreconditioner(LS, ark_mem, NULL, NULL);
    if (retval != SUN_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Error in calling SUNLinSolSetPreconditioner");
      return ARK_ILL_INPUT;
    }
  }

  if (step_mem->NLS != NULL)
  {
    retval = SUNNonlinSolSetLSetupFn(step_mem->NLS, pdaeStep_NlsLSetup);
    if (retval != SUN_SUCCESS) { return ARK_NLS_OP_ERR; }

    retval = SUNNonlinSolSetLSolveFn(step_mem->NLS, pdaeStep_NlsLSolve);
    if (retval != SUN_SUCCESS) { return ARK_NLS_OP_ERR; }
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  PDAEStepSetCouplingJacTimes:

  Attaches a Jacobian-vector product to the PDAEStep Newton solve.
  ---------------------------------------------------------------*/
int PDAEStepSetCouplingJacTimes(void* arkode_mem,
                                PDAEStepLsAlgebraicJacTimesVecFn jtimes)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  step_mem->algebraic_res_jtimes = jtimes;

  if ((step_mem->LS != NULL) && (step_mem->LS->ops != NULL) &&
      (step_mem->LS->ops->setatimes != NULL))
  {
    retval = SUNLinSolSetATimes(step_mem->LS, ark_mem,
                                (jtimes == NULL) ? pdaeStep_AlgDQJtimes
                                                 : pdaeStep_AlgJtimes);
    if (retval != SUN_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Error in calling SUNLinSolSetATimes");
      return ARK_ILL_INPUT;
    }
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  PDAEStep nonlinear solver settings.
  ---------------------------------------------------------------*/
int PDAEStepSetMaxNonlinIters(void* arkode_mem, int maxcor)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (step_mem->NLS == NULL)
  {
    arkProcessError(ark_mem, ARK_NLS_OP_ERR, __LINE__, __func__, __FILE__,
                    "No SUNNonlinearSolver object is present");
    return ARK_ILL_INPUT;
  }

  step_mem->maxcor = (maxcor <= 0) ? DEFAULT_MAX_COR : maxcor;

  retval = SUNNonlinSolSetMaxIters(step_mem->NLS, step_mem->maxcor);
  if (retval != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_NLS_OP_ERR, __LINE__, __func__, __FILE__,
                    "Error setting maxcor in SUNNonlinearSolver object");
    return ARK_NLS_OP_ERR;
  }

  return ARK_SUCCESS;
}

int PDAEStepSetNonlinConvCoef(void* arkode_mem, sunrealtype nlscoef)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  step_mem->nlscoef = (nlscoef <= ZERO) ? DEFAULT_NLSCOEF : nlscoef;

  return ARK_SUCCESS;
}

int PDAEStepSetNonlinCRDown(void* arkode_mem, sunrealtype crdown)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  step_mem->crdown = (crdown <= ZERO) ? DEFAULT_CRDOWN : crdown;

  return ARK_SUCCESS;
}

int PDAEStepSetNonlinRDiv(void* arkode_mem, sunrealtype rdiv)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  step_mem->rdiv = (rdiv <= ZERO) ? DEFAULT_RDIV : rdiv;

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  PDAEStep nonlinear solver statistics.
  ---------------------------------------------------------------*/
int PDAEStepGetNumLinSolvSetups(void* arkode_mem, long int* nlinsetups)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (nlinsetups == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "nlinsetups is NULL");
    return ARK_ILL_INPUT;
  }

  *nlinsetups = step_mem->nsetups;

  return ARK_SUCCESS;
}

int PDAEStepGetNumNonlinSolvIters(void* arkode_mem, long int* nniters)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (nniters == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "nniters is NULL");
    return ARK_ILL_INPUT;
  }

  *nniters = step_mem->nls_iters;

  return ARK_SUCCESS;
}

int PDAEStepGetNumNonlinSolvConvFails(void* arkode_mem, long int* nnfails)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (nnfails == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "nnfails is NULL");
    return ARK_ILL_INPUT;
  }

  *nnfails = step_mem->nls_fails;

  return ARK_SUCCESS;
}

int PDAEStepGetNonlinSolvStats(void* arkode_mem, long int* nniters,
                               long int* nnfails)
{
  int retval = PDAEStepGetNumNonlinSolvIters(arkode_mem, nniters);
  if (retval != ARK_SUCCESS) { return retval; }

  return PDAEStepGetNumNonlinSolvConvFails(arkode_mem, nnfails);
}

/*---------------------------------------------------------------
  pdaeStep_NlsInit:

  Initializes the PDAEStep nonlinear and linear solver objects.
  ---------------------------------------------------------------*/
int pdaeStep_NlsInit(ARKodeMem ark_mem)
{
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (step_mem->NLS == NULL)
  {
    arkProcessError(ark_mem, ARK_NLS_INIT_FAIL, __LINE__, __func__, __FILE__,
                    "No SUNNonlinearSolver object is present");
    return ARK_NLS_INIT_FAIL;
  }

  if (step_mem->LS == NULL)
  {
    arkProcessError(ark_mem, ARK_NLS_INIT_FAIL, __LINE__, __func__, __FILE__,
                    "No PDAEStep SUNLinearSolver object is present. Call "
                    "PDAEStepSetLinearSolver before ARKodeEvolve.");
    return ARK_NLS_INIT_FAIL;
  }

  retval = pdaeStep_AttachNlsCallbacks(ark_mem, step_mem);
  if (retval != ARK_SUCCESS) { return ARK_NLS_INIT_FAIL; }

  retval = SUNLinSolInitialize(step_mem->LS);
  if (retval != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_LINIT_FAIL, __LINE__, __func__, __FILE__,
                    "PDAEStep linear solver initialization failed");
    return ARK_LINIT_FAIL;
  }

  retval = SUNNonlinSolInitialize(step_mem->NLS);
  if (retval != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_NLS_INIT_FAIL, __LINE__, __func__, __FILE__,
                    "PDAEStep nonlinear solver initialization failed");
    return ARK_NLS_INIT_FAIL;
  }

  step_mem->nsetups   = 0;
  step_mem->nls_iters = 0;
  step_mem->nls_fails = 0;
  step_mem->nfeDQ     = 0;

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  pdaeStep_Nls:

  Solves the final coupled algebraic system over [Z, W].
  ---------------------------------------------------------------*/
int pdaeStep_Nls(ARKodeMem ark_mem, int nflag)
{
  ARKodePDAEStepMem step_mem = NULL;
  long int nls_iters_inc     = 0;
  long int nls_fails_inc     = 0;
  int retval = pdaeStep_AccessStepMem(ark_mem, __func__, &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  (void)nflag;

  N_VScale(ONE, step_mem->alg, step_mem->alg_pred);
  N_VConst(ZERO, step_mem->alg_cor);

  step_mem->crate    = ONE;
  step_mem->delnrm   = ZERO;
  step_mem->delnrm_p = ZERO;

  retval = SUNNonlinSolSolve(step_mem->NLS, step_mem->alg_pred,
                             step_mem->alg_cor, step_mem->alg_ewt,
                             step_mem->nlscoef, SUNTRUE, ark_mem);

  (void)SUNNonlinSolGetNumIters(step_mem->NLS, &nls_iters_inc);
  step_mem->nls_iters += nls_iters_inc;

  (void)SUNNonlinSolGetNumConvFails(step_mem->NLS, &nls_fails_inc);
  step_mem->nls_fails += nls_fails_inc;

  if (retval == SUN_SUCCESS)
  {
    step_mem->jcur = SUNFALSE;
    N_VLinearSum(ONE, step_mem->alg_pred, ONE, step_mem->alg_cor, step_mem->alg);
    return ARK_SUCCESS;
  }

  if (retval == SUN_NLS_CONV_RECVR) { return CONV_FAIL; }

  return retval;
}

/*---------------------------------------------------------------
  Interface routines supplied to SUNNonlinearSolver.
  ---------------------------------------------------------------*/
static int pdaeStep_NlsResidual(N_Vector alg_cor, N_Vector r, void* arkode_mem)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  N_VLinearSum(ONE, step_mem->alg_pred, ONE, alg_cor, step_mem->alg);

  retval = step_mem->algebraic_res_fn(ark_mem->tcur, ark_mem->ycur,
                                      PDAEStepGetCouplingSubvector(ark_mem->ycur),
                                      step_mem->alg_res, ark_mem->user_data);

  if (retval < 0) { return ARK_RHSFUNC_FAIL; }
  if (retval > 0) { return RHSFUNC_RECVR; }

  N_VScale(ONE, step_mem->alg_res, r);
  N_VScale(ONE, step_mem->alg_res, step_mem->alg_fcur);
  step_mem->nh++;

  return ARK_SUCCESS;
}

static int pdaeStep_NlsLSetup(SUNDIALS_MAYBE_UNUSED sunbooleantype jbad,
                              sunbooleantype* jcur, void* arkode_mem)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (step_mem->LS == NULL)
  {
    arkProcessError(ark_mem, ARK_NLS_OP_ERR, __LINE__, __func__, __FILE__,
                    "No PDAEStep SUNLinearSolver object is present");
    return ARK_NLS_OP_ERR;
  }

  if ((step_mem->J != NULL) && (step_mem->algebraic_res_jac != NULL))
  {
    if (step_mem->J->ops && step_mem->J->ops->zero)
    {
      retval = SUNMatZero(step_mem->J);
      if (retval != SUN_SUCCESS) { return ARK_LSETUP_FAIL; }
    }

    retval =
      step_mem->algebraic_res_jac(ark_mem->tcur, ark_mem->ycur,
                                  PDAEStepGetCouplingSubvector(ark_mem->ycur),
                                  step_mem->J, ark_mem->user_data);

    if (retval < 0) { return ARK_LSETUP_FAIL; }
    if (retval > 0) { return CONV_FAIL; }
  }

  retval = SUNLinSolSetup(step_mem->LS, step_mem->J);
  if (retval < 0) { return ARK_LSETUP_FAIL; }
  if (retval > 0) { return CONV_FAIL; }

  step_mem->nsetups++;
  step_mem->jcur = SUNTRUE;
  if (jcur != NULL) { *jcur = SUNTRUE; }
  step_mem->crate = ONE;

  return ARK_SUCCESS;
}

static int pdaeStep_NlsLSolve(N_Vector b, void* arkode_mem)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (step_mem->LS == NULL)
  {
    arkProcessError(ark_mem, ARK_NLS_OP_ERR, __LINE__, __func__, __FILE__,
                    "No PDAEStep SUNLinearSolver object is present");
    return ARK_NLS_OP_ERR;
  }

  N_VConst(ZERO, step_mem->alg_linsol_x);

  retval = SUNLinSolSetScalingVectors(step_mem->LS, step_mem->alg_ewt,
                                      step_mem->alg_ewt);
  if (retval != SUN_SUCCESS) { return ARK_LSOLVE_FAIL; }

  retval = SUNLinSolSetZeroGuess(step_mem->LS, SUNTRUE);
  if (retval != SUN_SUCCESS) { return ARK_LSOLVE_FAIL; }

  retval = SUNLinSolSolve(step_mem->LS, step_mem->J, step_mem->alg_linsol_x, b,
                          step_mem->nlscoef);
  if (retval == SUN_SUCCESS)
  {
    N_VScale(ONE, step_mem->alg_linsol_x, b);
    return ARK_SUCCESS;
  }

  if (retval < 0) { return ARK_LSOLVE_FAIL; }

  return CONV_FAIL;
}

static int pdaeStep_AlgJtimes(void* arkode_mem, N_Vector v, N_Vector Jv)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  retval =
    step_mem->algebraic_res_jtimes(ark_mem->tcur, ark_mem->ycur,
                                   PDAEStepGetCouplingSubvector(ark_mem->ycur),
                                   v, Jv, ark_mem->user_data, step_mem->alg_tmp);

  if (retval < 0) { return -1; }
  if (retval > 0) { return 1; }

  return 0;
}

static int pdaeStep_AlgDQJtimes(void* arkode_mem, N_Vector v, N_Vector Jv)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  sunrealtype vnorm = N_VWrmsNorm(v, step_mem->alg_ewt);
  if (vnorm == ZERO)
  {
    N_VConst(ZERO, Jv);
    return ARK_SUCCESS;
  }

  sunrealtype ynorm = N_VWrmsNorm(step_mem->alg, step_mem->alg_ewt);
  sunrealtype sigma = SUNRsqrt(ark_mem->uround) * SUNMAX(ynorm, ONE) / vnorm;

  N_VLinearSum(ONE, step_mem->alg, sigma, v, step_mem->alg_tmp);
  N_VScale(ONE, step_mem->alg_tmp, step_mem->alg);

  retval = step_mem->algebraic_res_fn(ark_mem->tcur, ark_mem->ycur,
                                      PDAEStepGetCouplingSubvector(ark_mem->ycur),
                                      Jv, ark_mem->user_data);

  N_VLinearSum(ONE, step_mem->alg_pred, ONE, step_mem->alg_cor, step_mem->alg);

  if (retval < 0) { return -1; }
  if (retval > 0) { return 1; }

  N_VLinearSum(ONE / sigma, Jv, -ONE / sigma, step_mem->alg_fcur, Jv);
  step_mem->nfeDQ++;

  return 0;
}

static int pdaeStep_NlsConvTest(SUNNonlinearSolver NLS,
                                SUNDIALS_MAYBE_UNUSED N_Vector y, N_Vector del,
                                sunrealtype tol, N_Vector ewt, void* arkode_mem)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  sunrealtype dcon;
  int m, retval;

  retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                        &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  if (pdaeStep_NlsNorm(del, ewt, &step_mem->delnrm, arkode_mem) != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_NLS_OP_ERR, __LINE__, __func__, __FILE__,
                    MSG_ARK_NLS_FAIL, ark_mem->tcur);
    return ARK_NLS_OP_ERR;
  }

  retval = SUNNonlinSolGetCurIter(NLS, &m);
  if (retval != SUN_SUCCESS) { return ARK_NLS_OP_ERR; }

  if (m > 0)
  {
    step_mem->crate = SUNMAX(step_mem->crdown * step_mem->crate,
                             step_mem->delnrm / step_mem->delnrm_p);
  }

  dcon = SUNMIN(step_mem->crate, ONE) * step_mem->delnrm / tol;

  if (dcon <= ONE) { return SUN_SUCCESS; }

  if ((m >= 1) && (step_mem->delnrm > step_mem->rdiv * step_mem->delnrm_p))
  {
    return SUN_NLS_CONV_RECVR;
  }

  step_mem->delnrm_p = step_mem->delnrm;

  return SUN_NLS_CONTINUE;
}

static int pdaeStep_AttachNlsCallbacks(ARKodeMem ark_mem,
                                       ARKodePDAEStepMem step_mem)
{
  int retval;

  retval = SUNNonlinSolSetSysFn(step_mem->NLS, pdaeStep_NlsResidual);
  if (retval != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Setting nonlinear system function failed");
    return ARK_ILL_INPUT;
  }

  retval = SUNNonlinSolSetConvTestFn(step_mem->NLS, pdaeStep_NlsConvTest,
                                     (void*)ark_mem);
  if (retval != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Setting convergence test function failed");
    return ARK_ILL_INPUT;
  }

  retval = SUNNonlinSolSetNormFn(step_mem->NLS, pdaeStep_NlsNorm, ark_mem);
  if (retval != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Setting convergence-test norm function failed");
    return ARK_ILL_INPUT;
  }

  retval = SUNNonlinSolSetGetUpdateNormFn(step_mem->NLS,
                                          pdaeStep_NlsGetUpdateNorm, ark_mem);
  if (retval != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Setting update-norm getter failed");
    return ARK_ILL_INPUT;
  }

  retval = SUNNonlinSolSetLSetupFn(step_mem->NLS, pdaeStep_NlsLSetup);
  if (retval != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Setting linear solver setup function failed");
    return ARK_ILL_INPUT;
  }

  retval = SUNNonlinSolSetLSolveFn(step_mem->NLS, pdaeStep_NlsLSolve);
  if (retval != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Setting linear solver solve function failed");
    return ARK_ILL_INPUT;
  }

  retval = SUNNonlinSolSetMaxIters(step_mem->NLS, step_mem->maxcor);
  if (retval != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Setting maximum number of nonlinear iterations failed");
    return ARK_ILL_INPUT;
  }

  return ARK_SUCCESS;
}

static SUNErrCode pdaeStep_NlsNorm(N_Vector del, N_Vector ewt,
                                   sunrealtype* delnrm,
                                   SUNDIALS_MAYBE_UNUSED void* arkode_mem)
{
  *delnrm = N_VWrmsNorm(del, ewt);
  return SUN_SUCCESS;
}

static SUNErrCode pdaeStep_NlsGetUpdateNorm(sunrealtype* delnrm, void* arkode_mem)
{
  ARKodeMem ark_mem          = NULL;
  ARKodePDAEStepMem step_mem = NULL;
  int retval = pdaeStep_AccessARKODEStepMem(arkode_mem, __func__, &ark_mem,
                                            &step_mem);
  if (retval != ARK_SUCCESS) { return retval; }

  *delnrm = step_mem->delnrm;

  return SUN_SUCCESS;
}
