/* -----------------------------------------------------------------
 * Programmer(s): Alan C. Hindmarsh and Radu Serban @ LLNL
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
 * This is the implementation file for the optional input and output
 * functions for the CVODES solver.
 * ----------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>

#include "cvodes_impl.h"
#include "cvodes_ls_impl.h"
#include "sundials_utils.h"

#define ZERO   SUN_RCONST(0.0)
#define HALF   SUN_RCONST(0.5)
#define ONE    SUN_RCONST(1.0)
#define TWOPT5 SUN_RCONST(2.5)

/*
 * =================================================================
 * CVODES optional input functions
 * =================================================================
 */

/*
 * CVodeSetDeltaGammaMaxLSetup
 *
 * Specifies the gamma ratio threshold to signal for a linear solver setup
 */

int CVodeSetDeltaGammaMaxLSetup(void* cvode_mem, sunrealtype dgmax_lsetup)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  /* Set value or use default */
  if (dgmax_lsetup < ZERO) { cv_mem->cv_dgmax_lsetup = DGMAX_LSETUP_DEFAULT; }
  else { cv_mem->cv_dgmax_lsetup = dgmax_lsetup; }

  return (CV_SUCCESS);
}

/*
 * CVodeSetUserData
 *
 * Specifies the user data pointer for f
 */

int CVodeSetUserData(void* cvode_mem, void* user_data)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  cv_mem->cv_user_data = user_data;

  return (CV_SUCCESS);
}

/*
 * CVodeSetMonitorFn
 *
 * Specifies the user function to call for monitoring
 * the solution and/or integrator statistics.
 */

int CVodeSetMonitorFn(void* cvode_mem, CVMonitorFn fn)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

#ifdef SUNDIALS_BUILD_WITH_MONITORING
  cv_mem->cv_monitorfun = fn;
  return (CV_SUCCESS);
#else
  /* silence warnings when monitoring is disabled */
  ((void)fn);
  cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                 "SUNDIALS was not built with monitoring enabled.");
  return (CV_ILL_INPUT);
#endif
}

/*
 * CVodeSetMonitorFrequency
 *
 * Specifies the frequency with which to call the user function.
 */

int CVodeSetMonitorFrequency(void* cvode_mem, long int nst)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  if (nst < 0)
  {
    cvProcessError(NULL, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   "step interval must be >= 0\n");
    return (CV_ILL_INPUT);
  }

  cv_mem = (CVodeMem)cvode_mem;

#ifdef SUNDIALS_BUILD_WITH_MONITORING
  cv_mem->cv_monitor_interval = nst;
  return (CV_SUCCESS);
#else
  cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                 "SUNDIALS was not built with monitoring enabled.");
  return (CV_ILL_INPUT);
#endif
}

/*
 * CVodeSetMaxOrd
 *
 * Specifies the maximum method order
 */

int CVodeSetMaxOrd(void* cvode_mem, int maxord)
{
  CVodeMem cv_mem;
  int qmax_alloc;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (maxord <= 0)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   MSGCV_NEG_MAXORD);
    return (CV_ILL_INPUT);
  }

  /* Cannot increase maximum order beyond the value that
     was used when allocating memory */
  qmax_alloc = cv_mem->cv_qmax_alloc;
  qmax_alloc = SUNMIN(qmax_alloc, cv_mem->cv_qmax_allocQ);
  qmax_alloc = SUNMIN(qmax_alloc, cv_mem->cv_qmax_allocS);

  if (maxord > qmax_alloc)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   MSGCV_BAD_MAXORD);
    return (CV_ILL_INPUT);
  }

  cv_mem->cv_qmax = maxord;

  return (CV_SUCCESS);
}

/*
 * CVodeSetMaxNumSteps
 *
 * Specifies the maximum number of integration steps
 */

int CVodeSetMaxNumSteps(void* cvode_mem, long int mxsteps)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  /* Passing mxsteps=0 sets the default. Passing mxsteps<0 disables the test. */
  if (mxsteps == 0) { cv_mem->cv_mxstep = MXSTEP_DEFAULT; }
  else { cv_mem->cv_mxstep = mxsteps; }

  return (CV_SUCCESS);
}

/*
 * CVodeSetMaxHnilWarns
 *
 * Specifies the maximum number of warnings for small h
 */

int CVodeSetMaxHnilWarns(void* cvode_mem, int mxhnil)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  cv_mem->cv_mxhnil = mxhnil;

  return (CV_SUCCESS);
}

/*
 *CVodeSetStabLimDet
 *
 * Turns on/off the stability limit detection algorithm
 */

int CVodeSetStabLimDet(void* cvode_mem, sunbooleantype sldet)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (sldet && (cv_mem->cv_lmm != CV_BDF))
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   MSGCV_SET_SLDET);
    return (CV_ILL_INPUT);
  }

  cv_mem->cv_sldeton = sldet;

  return (CV_SUCCESS);
}

/*
 * CVodeSetInitStep
 *
 * Specifies the initial step size
 */

int CVodeSetInitStep(void* cvode_mem, sunrealtype hin)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  cv_mem->cv_hin = hin;

  return (CV_SUCCESS);
}

/*
 * CVodeSetMinStep
 *
 * Specifies the minimum step size
 */

int CVodeSetMinStep(void* cvode_mem, sunrealtype hmin)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (hmin < ZERO)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   MSGCV_NEG_HMIN);
    return (CV_ILL_INPUT);
  }

  /* Passing 0 sets hmin = zero */
  if (hmin == ZERO)
  {
    cv_mem->cv_hmin = HMIN_DEFAULT;
    return (CV_SUCCESS);
  }

  if (hmin * cv_mem->cv_hmax_inv > ONE)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   MSGCV_BAD_HMIN_HMAX);
    return (CV_ILL_INPUT);
  }

  cv_mem->cv_hmin = hmin;

  return (CV_SUCCESS);
}

/*
 * CVodeSetMaxStep
 *
 * Specifies the maximum step size
 */

int CVodeSetMaxStep(void* cvode_mem, sunrealtype hmax)
{
  sunrealtype hmax_inv;
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (hmax < ZERO)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   MSGCV_NEG_HMAX);
    return (CV_ILL_INPUT);
  }

  /* Passing 0 sets hmax = infinity */
  if (hmax == ZERO)
  {
    cv_mem->cv_hmax_inv = HMAX_INV_DEFAULT;
    return (CV_SUCCESS);
  }

  hmax_inv = ONE / hmax;
  if (hmax_inv * cv_mem->cv_hmin > ONE)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   MSGCV_BAD_HMIN_HMAX);
    return (CV_ILL_INPUT);
  }

  cv_mem->cv_hmax_inv = hmax_inv;

  return (CV_SUCCESS);
}

/*
 * CVodeSetEtaFixedStepBounds
 *
 * Specifies the bounds for retaining the current step size
 */

int CVodeSetEtaFixedStepBounds(void* cvode_mem, sunrealtype eta_min_fx,
                               sunrealtype eta_max_fx)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  /* set allowed value or use default */
  if (eta_min_fx >= ZERO && eta_min_fx <= ONE)
  {
    cv_mem->cv_eta_min_fx = eta_min_fx;
  }
  else { cv_mem->cv_eta_min_fx = ETA_MIN_FX_DEFAULT; }

  if (eta_max_fx >= ONE) { cv_mem->cv_eta_max_fx = eta_max_fx; }
  else { cv_mem->cv_eta_max_fx = ETA_MAX_FX_DEFAULT; }

  return (CV_SUCCESS);
}

/*
 * CVodeSetEtaMaxFirstStep
 *
 * Specifies the maximum step size change on the first step
 */

int CVodeSetEtaMaxFirstStep(void* cvode_mem, sunrealtype eta_max_fs)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  /* set allowed value or use default */
  if (eta_max_fs <= ONE) { cv_mem->cv_eta_max_fs = ETA_MAX_FS_DEFAULT; }
  else { cv_mem->cv_eta_max_fs = eta_max_fs; }

  return (CV_SUCCESS);
}

/*
 * CVodeSetEtaMaxEarlyStep
 *
 * Specifies the maximum step size change on steps early in the integration
 * when nst <= small_nst
 */

int CVodeSetEtaMaxEarlyStep(void* cvode_mem, sunrealtype eta_max_es)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  /* set allowed value or use default */
  if (eta_max_es <= ONE) { cv_mem->cv_eta_max_es = ETA_MAX_ES_DEFAULT; }
  else { cv_mem->cv_eta_max_es = eta_max_es; }

  return (CV_SUCCESS);
}

/*
 * CVodeSetNumStepsEtaMaxEarlyStep
 *
 * Specifies the maximum number of steps for using the early integration change
 * factor
 */

int CVodeSetNumStepsEtaMaxEarlyStep(void* cvode_mem, long int small_nst)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  /* set allowed value or use default */
  if (small_nst < 0) { cv_mem->cv_small_nst = SMALL_NST_DEFAULT; }
  else { cv_mem->cv_small_nst = small_nst; }

  return (CV_SUCCESS);
}

/*
 * CVodeSetEtaMax
 *
 * Specifies the maximum step size change on a general steps (nst > small_nst)
 */

int CVodeSetEtaMax(void* cvode_mem, sunrealtype eta_max_gs)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  /* set allowed value or use default */
  if (eta_max_gs <= ONE) { cv_mem->cv_eta_max_gs = ETA_MAX_GS_DEFAULT; }
  else { cv_mem->cv_eta_max_gs = eta_max_gs; }

  return (CV_SUCCESS);
}

/*
 * CVodeSetEtaMin
 *
 * Specifies the minimum change on a general steps
 */

int CVodeSetEtaMin(void* cvode_mem, sunrealtype eta_min)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  /* set allowed value or use default */
  if (eta_min <= ZERO || eta_min >= ONE)
  {
    cv_mem->cv_eta_min = ETA_MIN_DEFAULT;
  }
  else { cv_mem->cv_eta_min = eta_min; }

  return (CV_SUCCESS);
}

/*
 * CVodeSetEtaMinErrFail
 *
 * Specifies the minimum step size change after an error test failure
 */

int CVodeSetEtaMinErrFail(void* cvode_mem, sunrealtype eta_min_ef)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  /* set allowed value or use default */
  if (eta_min_ef <= ZERO || eta_min_ef >= ONE)
  {
    cv_mem->cv_eta_min_ef = ETA_MIN_EF_DEFAULT;
  }
  else { cv_mem->cv_eta_min_ef = eta_min_ef; }

  return (CV_SUCCESS);
}

/*
 * CVodeSetEtaMaxErrFail
 *
 * Specifies the maximum step size change after multiple (>= small_nef) error
 * test failures
 */

int CVodeSetEtaMaxErrFail(void* cvode_mem, sunrealtype eta_max_ef)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  /* set allowed value or use default */
  if (eta_max_ef <= ZERO || eta_max_ef >= ONE)
  {
    cv_mem->cv_eta_max_ef = ETA_MAX_EF_DEFAULT;
  }
  else { cv_mem->cv_eta_max_ef = eta_max_ef; }

  return (CV_SUCCESS);
}

/*
 * CVodeSetNumFailsEtaMaxErrFail
 *
 * Specifies the maximum number of error test failures necessary to enforce
 * eta_max_ef
 */

int CVodeSetNumFailsEtaMaxErrFail(void* cvode_mem, int small_nef)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  /* set allowed value or use default */
  if (small_nef < 0) { cv_mem->cv_small_nef = SMALL_NEF_DEFAULT; }
  else { cv_mem->cv_small_nef = small_nef; }

  return (CV_SUCCESS);
}

/*
 * CVodeSetEtaConvFail
 *
 * Specifies the step size change after a nonlinear solver failure
 */

int CVodeSetEtaConvFail(void* cvode_mem, sunrealtype eta_cf)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  /* set allowed value or use default */
  if (eta_cf <= ZERO || eta_cf >= ONE) { cv_mem->cv_eta_cf = ETA_CF_DEFAULT; }
  else { cv_mem->cv_eta_cf = eta_cf; }

  return (CV_SUCCESS);
}

/*
 * CVodeSetStopTime
 *
 * Specifies the time beyond which the integration is not to proceed.
 */

int CVodeSetStopTime(void* cvode_mem, sunrealtype tstop)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }
  cv_mem = (CVodeMem)cvode_mem;

  /* If CVode was called at least once, test if tstop is legal
   * (i.e. if it was not already passed).
   * If CVodeSetStopTime is called before the first call to CVode,
   * tstop will be checked in CVode. */
  if (cv_mem->cv_nst > 0)
  {
    if ((tstop - cv_mem->cv_tn) * cv_mem->cv_h < ZERO)
    {
      cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                     MSGCV_BAD_TSTOP, tstop, cv_mem->cv_tn);
      return (CV_ILL_INPUT);
    }
  }

  cv_mem->cv_tstop    = tstop;
  cv_mem->cv_tstopset = SUNTRUE;

  return (CV_SUCCESS);
}

/*
 * CVodeSetInterpolateStopTime
 *
 * Specifies to use interpolation to fill the output solution at
 * the stop time (instead of a copy).
 */

int CVodeSetInterpolateStopTime(void* cvode_mem, sunbooleantype interp)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }
  cv_mem                 = (CVodeMem)cvode_mem;
  cv_mem->cv_tstopinterp = interp;

  return (CV_SUCCESS);
}

/*
 * CVodeClearStopTime
 *
 * Disable the stop time.
 */

int CVodeClearStopTime(void* cvode_mem)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }
  cv_mem = (CVodeMem)cvode_mem;

  cv_mem->cv_tstopset = SUNFALSE;

  return (CV_SUCCESS);
}

/*
 * CVodeSetMaxErrTestFails
 *
 * Specifies the maximum number of error test failures during one
 * step try.
 */

int CVodeSetMaxErrTestFails(void* cvode_mem, int maxnef)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  cv_mem->cv_maxnef = maxnef;

  return (CV_SUCCESS);
}

/*
 * CVodeSetMaxConvFails
 *
 * Specifies the maximum number of nonlinear convergence failures
 * during one step try.
 */

int CVodeSetMaxConvFails(void* cvode_mem, int maxncf)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  cv_mem->cv_maxncf = maxncf;

  return (CV_SUCCESS);
}

/*
 * CVodeSetMaxNonlinIters
 *
 * Specifies the maximum number of nonlinear iterations during
 * one solve.
 */

int CVodeSetMaxNonlinIters(void* cvode_mem, int maxcor)
{
  CVodeMem cv_mem;
  sunbooleantype sensi_sim;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  /* Are we computing sensitivities with the simultaneous approach? */
  sensi_sim = (cv_mem->cv_sensi && (cv_mem->cv_ism == CV_SIMULTANEOUS));

  if (sensi_sim)
  {
    /* check that the NLS is non-NULL */
    if (cv_mem->NLSsim == NULL)
    {
      cvProcessError(NULL, CV_MEM_FAIL, __LINE__, __func__, __FILE__,
                     MSGCV_MEM_FAIL);
      return (CV_MEM_FAIL);
    }

    return (SUNNonlinSolSetMaxIters(cv_mem->NLSsim, maxcor));
  }
  else
  {
    /* check that the NLS is non-NULL */
    if (cv_mem->NLS == NULL)
    {
      cvProcessError(NULL, CV_MEM_FAIL, __LINE__, __func__, __FILE__,
                     MSGCV_MEM_FAIL);
      return (CV_MEM_FAIL);
    }

    return (SUNNonlinSolSetMaxIters(cv_mem->NLS, maxcor));
  }

  return (CV_SUCCESS);
}

/*
 * CVodeSetNonlinConvCoef
 *
 * Specifies the coefficient in the nonlinear solver convergence
 * test
 */

int CVodeSetNonlinConvCoef(void* cvode_mem, sunrealtype nlscoef)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  cv_mem->cv_nlscoef = nlscoef;

  return (CV_SUCCESS);
}

/*
 * CVodeSetLSetupFrequency
 *
 * Specifies the frequency for calling the linear solver setup function to
 * recompute the Jacobian matrix and/or preconditioner
 */

int CVodeSetLSetupFrequency(void* cvode_mem, long int msbp)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  /* check for a valid input */
  if (msbp < 0)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   "A negative setup frequency was provided");
    return (CV_ILL_INPUT);
  }

  /* use default or user provided value */
  cv_mem->cv_msbp = (msbp == 0) ? MSBP_DEFAULT : msbp;

  return (CV_SUCCESS);
}

/*
 * CVodeSetRootDirection
 *
 * Specifies the direction of zero-crossings to be monitored.
 * The default is to monitor both crossings.
 */

int CVodeSetRootDirection(void* cvode_mem, int* rootdir)
{
  CVodeMem cv_mem;
  int i, nrt;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  nrt = cv_mem->cv_nrtfn;
  if (nrt == 0)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   MSGCV_NO_ROOT);
    return (CV_ILL_INPUT);
  }

  for (i = 0; i < nrt; i++) { cv_mem->cv_rootdir[i] = rootdir[i]; }

  return (CV_SUCCESS);
}

/*
 * CVodeSetNoInactiveRootWarn
 *
 * Disables issuing a warning if some root function appears
 * to be identically zero at the beginning of the integration
 */

int CVodeSetNoInactiveRootWarn(void* cvode_mem)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  cv_mem->cv_mxgnull = 0;

  return (CV_SUCCESS);
}

/*
 * CVodeSetConstraints
 *
 * Setup for constraint handling feature
 */

int CVodeSetConstraints(void* cvode_mem, N_Vector constraints)
{
  CVodeMem cv_mem;
  sunrealtype temptest;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  /* Disable constraints */
  if (constraints == NULL)
  {
    if (cv_mem->cv_constraints)
    {
      N_VDestroy(cv_mem->cv_constraints);
      cv_mem->cv_lrw -= cv_mem->cv_lrw1;
      cv_mem->cv_liw -= cv_mem->cv_liw1;
    }
    return (CV_SUCCESS);
  }

  /* Test if required vector ops. are defined */
  if (constraints->ops->nvdiv == NULL || constraints->ops->nvmaxnorm == NULL ||
      constraints->ops->nvcompare == NULL ||
      constraints->ops->nvconstrmask == NULL ||
      constraints->ops->nvminquotient == NULL)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   MSGCV_BAD_NVECTOR);
    return (CV_ILL_INPUT);
  }

  /* Check the constraints vector */
  temptest = N_VMaxNorm(constraints);
  if ((temptest > TWOPT5) || (temptest < HALF))
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   MSGCV_BAD_CONSTR);
    return (CV_ILL_INPUT);
  }

  /* Enable constraints */
  if (cv_mem->cv_constraints == NULL)
  {
    cv_mem->cv_constraints = N_VClone(constraints);
    if (cv_mem->cv_constraints == NULL)
    {
      cvProcessError(NULL, CV_MEM_FAIL, __LINE__, __func__, __FILE__,
                     MSGCV_MEM_FAIL);
      return (CV_MEM_FAIL);
    }
    cv_mem->cv_lrw += cv_mem->cv_lrw1;
    cv_mem->cv_liw += cv_mem->cv_liw1;
  }

  /* Load the constraints vector */
  N_VScale(ONE, constraints, cv_mem->cv_constraints);

  return (CV_SUCCESS);
}

/*
 * CVodeSetMaxNumConstraintFails
 *
 * Set the maximum number of constraint failure allowed in a step
 */

int CVodeSetMaxNumConstraintFails(void* cvode_mem, int max_fails)
{
  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }
  CVodeMem cv_mem = (CVodeMem)cvode_mem;

  if (max_fails <= 0) { cv_mem->max_constraint_fails = MAX_CONSTRAINT_FAILS; }
  else { cv_mem->max_constraint_fails = max_fails; }

  return CV_SUCCESS;
}

/*
 * CVodeGetNumConstraintFails
 *
 * Get the number of failed steps due to constraint violation
 */

int CVodeGetNumConstraintFails(void* cvode_mem, long int* num_fails_out)
{
  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }
  CVodeMem cv_mem = (CVodeMem)cvode_mem;

  *num_fails_out = cv_mem->constraint_fails;

  return CV_SUCCESS;
}

/*
 * CVodeGetNumConstraintCorrections
 *
 * Get the number of constraint corrections
 */

int CVodeGetNumConstraintCorrections(void* cvode_mem,
                                     long int* num_corrections_out)
{
  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }
  CVodeMem cv_mem = (CVodeMem)cvode_mem;

  *num_corrections_out = cv_mem->constraint_corrections;

  return CV_SUCCESS;
}

/*
 * =================================================================
 * Quadrature optional input functions
 * =================================================================
 */

int CVodeSetQuadErrCon(void* cvode_mem, sunbooleantype errconQ)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }
  cv_mem = (CVodeMem)cvode_mem;

  cv_mem->cv_errconQ = errconQ;

  return (CV_SUCCESS);
}

/*
 * =================================================================
 * FSA optional input functions
 * =================================================================
 */

int CVodeSetSensDQMethod(void* cvode_mem, int DQtype, sunrealtype DQrhomax)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if ((DQtype != CV_CENTERED) && (DQtype != CV_FORWARD))
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   MSGCV_BAD_DQTYPE);
    return (CV_ILL_INPUT);
  }

  if (DQrhomax < ZERO)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   MSGCV_BAD_DQRHO);
    return (CV_ILL_INPUT);
  }

  cv_mem->cv_DQtype   = DQtype;
  cv_mem->cv_DQrhomax = DQrhomax;

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeSetSensErrCon(void* cvode_mem, sunbooleantype errconS)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }
  cv_mem = (CVodeMem)cvode_mem;

  cv_mem->cv_errconS = errconS;

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeSetSensMaxNonlinIters(void* cvode_mem, int maxcorS)
{
  CVodeMem cv_mem;
  sunbooleantype sensi_stg;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  /* Are we computing sensitivities with a staggered approach? */
  sensi_stg = (cv_mem->cv_sensi && (cv_mem->cv_ism == CV_STAGGERED));

  if (sensi_stg)
  {
    /* check that the NLS is non-NULL */
    if (cv_mem->NLSstg == NULL)
    {
      cvProcessError(NULL, CV_MEM_FAIL, __LINE__, __func__, __FILE__,
                     MSGCV_MEM_FAIL);
      return (CV_MEM_FAIL);
    }

    return (SUNNonlinSolSetMaxIters(cv_mem->NLSstg, maxcorS));
  }
  else
  {
    /* check that the NLS is non-NULL */
    if (cv_mem->NLSstg1 == NULL)
    {
      cvProcessError(NULL, CV_MEM_FAIL, __LINE__, __func__, __FILE__,
                     MSGCV_MEM_FAIL);
      return (CV_MEM_FAIL);
    }

    return (SUNNonlinSolSetMaxIters(cv_mem->NLSstg1, maxcorS));
  }

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeSetSensParams(void* cvode_mem, sunrealtype* p, sunrealtype* pbar,
                       int* plist)
{
  CVodeMem cv_mem;
  int is, Ns;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  /* Was sensitivity initialized? */

  if (cv_mem->cv_SensMallocDone == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_SENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_SENSI);
    return (CV_NO_SENS);
  }

  Ns = cv_mem->cv_Ns;

  /* Parameters */

  cv_mem->cv_p = p;

  /* pbar */

  if (pbar != NULL)
  {
    for (is = 0; is < Ns; is++)
    {
      if (pbar[is] == ZERO)
      {
        cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                       MSGCV_BAD_PBAR);
        return (CV_ILL_INPUT);
      }
      cv_mem->cv_pbar[is] = SUNRabs(pbar[is]);
    }
  }
  else
  {
    for (is = 0; is < Ns; is++) { cv_mem->cv_pbar[is] = ONE; }
  }

  /* plist */

  if (plist != NULL)
  {
    for (is = 0; is < Ns; is++)
    {
      if (plist[is] < 0)
      {
        cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                       MSGCV_BAD_PLIST);
        return (CV_ILL_INPUT);
      }
      cv_mem->cv_plist[is] = plist[is];
    }
  }
  else
  {
    for (is = 0; is < Ns; is++) { cv_mem->cv_plist[is] = is; }
  }

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeSetQuadSensErrCon(void* cvode_mem, sunbooleantype errconQS)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }
  cv_mem = (CVodeMem)cvode_mem;

  /* Was sensitivity initialized? */

  if (cv_mem->cv_SensMallocDone == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_SENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_SENSI);
    return (CV_NO_SENS);
  }

  /* Check if quadrature sensitivity was initialized? */

  if (cv_mem->cv_QuadSensMallocDone == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_QUADSENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_QUADSENSI);
    return (CV_NO_QUAD);
  }

  cv_mem->cv_errconQS = errconQS;

  return (CV_SUCCESS);
}

/*
 * =================================================================
 * CVODES optional output functions
 * =================================================================
 */

/*
 * CVodeGetNumSteps
 *
 * Returns the current number of integration steps
 */

int CVodeGetNumSteps(void* cvode_mem, long int* nsteps)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *nsteps = cv_mem->cv_nst;

  return (CV_SUCCESS);
}

/*
 * CVodeGetNumRhsEvals
 *
 * Returns the current number of calls to f
 */

int CVodeGetNumRhsEvals(void* cvode_mem, long int* nfevals)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *nfevals = cv_mem->cv_nfe;

  return (CV_SUCCESS);
}

/*
 * CVodeGetNumLinSolvSetups
 *
 * Returns the current number of calls to the linear solver setup routine
 */

int CVodeGetNumLinSolvSetups(void* cvode_mem, long int* nlinsetups)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *nlinsetups = cv_mem->cv_nsetups;

  return (CV_SUCCESS);
}

/*
 * CVodeGetNumErrTestFails
 *
 * Returns the current number of error test failures
 */

int CVodeGetNumErrTestFails(void* cvode_mem, long int* netfails)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *netfails = cv_mem->cv_netf;

  return (CV_SUCCESS);
}

/*
 * CVodeGetLastOrder
 *
 * Returns the order on the last successful step
 */

int CVodeGetLastOrder(void* cvode_mem, int* qlast)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *qlast = cv_mem->cv_qu;

  return (CV_SUCCESS);
}

/*
 * CVodeGetCurrentOrder
 *
 * Returns the order to be attempted on the next step
 */

int CVodeGetCurrentOrder(void* cvode_mem, int* qcur)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *qcur = cv_mem->cv_next_q;

  return (CV_SUCCESS);
}

/*
 * CVodeGetCurrentGamma
 *
 * Returns the value of gamma for the current step.
 */

int CVodeGetCurrentGamma(void* cvode_mem, sunrealtype* gamma)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *gamma = cv_mem->cv_gamma;

  return (CV_SUCCESS);
}

/*
 * CVodeGetNumStabLimOrderReds
 *
 * Returns the number of order reductions triggered by the stability
 * limit detection algorithm
 */

int CVodeGetNumStabLimOrderReds(void* cvode_mem, long int* nslred)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_sldeton == SUNFALSE) { *nslred = 0; }
  else { *nslred = cv_mem->cv_nor; }

  return (CV_SUCCESS);
}

/*
 * CVodeGetActualInitStep
 *
 * Returns the step size used on the first step
 */

int CVodeGetActualInitStep(void* cvode_mem, sunrealtype* hinused)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *hinused = cv_mem->cv_h0u;

  return (CV_SUCCESS);
}

/*
 * CVodeGetLastStep
 *
 * Returns the step size used on the last successful step
 */

int CVodeGetLastStep(void* cvode_mem, sunrealtype* hlast)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *hlast = cv_mem->cv_hu;

  return (CV_SUCCESS);
}

/*
 * CVodeGetCurrentStep
 *
 * Returns the step size to be attempted on the next step
 */

int CVodeGetCurrentStep(void* cvode_mem, sunrealtype* hcur)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *hcur = cv_mem->cv_next_h;

  return (CV_SUCCESS);
}

/*
 * CVodeGetCurrentState
 *
 * Returns the current state vector
 */

int CVodeGetCurrentState(void* cvode_mem, N_Vector* y)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *y = cv_mem->cv_y;

  return (CV_SUCCESS);
}

/*
 * CVodeGetCurrentStateSens
 *
 * Returns the current sensitivity state vector array
 */

int CVodeGetCurrentStateSens(void* cvode_mem, N_Vector** yS)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *yS = cv_mem->cv_yS;

  return (CV_SUCCESS);
}

/*
 * CVodeGetCurrentSensSolveIndex
 *
 * Returns the current index of the sensitivity solve when using
 * the staggered1 nonlinear solver.
 */

int CVodeGetCurrentSensSolveIndex(void* cvode_mem, int* index)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *index = cv_mem->sens_solve_idx;

  return (CV_SUCCESS);
}

/*
 * CVodeGetCurrentTime
 *
 * Returns the current value of the independent variable
 */

int CVodeGetCurrentTime(void* cvode_mem, sunrealtype* tcur)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *tcur = cv_mem->cv_tn;

  return (CV_SUCCESS);
}

/*
 * CVodeGetTolScaleFactor
 *
 * Returns a suggested factor for scaling tolerances
 */

int CVodeGetTolScaleFactor(void* cvode_mem, sunrealtype* tolsfact)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *tolsfact = cv_mem->cv_tolsf;

  return (CV_SUCCESS);
}

/*
 * CVodeGetErrWeights
 *
 * This routine returns the current weight vector.
 */

int CVodeGetErrWeights(void* cvode_mem, N_Vector eweight)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  N_VScale(ONE, cv_mem->cv_ewt, eweight);

  return (CV_SUCCESS);
}

/*
 * CVodeGetEstLocalErrors
 *
 * Returns an estimate of the local error
 */

int CVodeGetEstLocalErrors(void* cvode_mem, N_Vector ele)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  N_VScale(ONE, cv_mem->cv_acor, ele);

  return (CV_SUCCESS);
}

/*
 * CVodeGetWorkSpace
 *
 * Returns integrator work space requirements
 */

int CVodeGetWorkSpace(void* cvode_mem, long int* lenrw, long int* leniw)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *leniw = cv_mem->cv_liw;
  *lenrw = cv_mem->cv_lrw;

  return (CV_SUCCESS);
}

/*
 * CVodeGetIntegratorStats
 *
 * Returns integrator statistics
 */

int CVodeGetIntegratorStats(void* cvode_mem, long int* nsteps, long int* nfevals,
                            long int* nlinsetups, long int* netfails, int* qlast,
                            int* qcur, sunrealtype* hinused, sunrealtype* hlast,
                            sunrealtype* hcur, sunrealtype* tcur)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *nsteps     = cv_mem->cv_nst;
  *nfevals    = cv_mem->cv_nfe;
  *nlinsetups = cv_mem->cv_nsetups;
  *netfails   = cv_mem->cv_netf;
  *qlast      = cv_mem->cv_qu;
  *qcur       = cv_mem->cv_next_q;
  *hinused    = cv_mem->cv_h0u;
  *hlast      = cv_mem->cv_hu;
  *hcur       = cv_mem->cv_next_h;
  *tcur       = cv_mem->cv_tn;

  return (CV_SUCCESS);
}

/*
 * CVodeGetNumGEvals
 *
 * Returns the current number of calls to g (for rootfinding)
 */

int CVodeGetNumGEvals(void* cvode_mem, long int* ngevals)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *ngevals = cv_mem->cv_nge;

  return (CV_SUCCESS);
}

/*
 * CVodeGetRootInfo
 *
 * Returns pointer to array rootsfound showing roots found
 */

int CVodeGetRootInfo(void* cvode_mem, int* rootsfound)
{
  CVodeMem cv_mem;
  int i, nrt;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  nrt = cv_mem->cv_nrtfn;

  for (i = 0; i < nrt; i++) { rootsfound[i] = cv_mem->cv_iroots[i]; }

  return (CV_SUCCESS);
}

/*
 * CVodeGetNumNonlinSolvIters
 *
 * Returns the current number of iterations in the nonlinear solver
 */

int CVodeGetNumNonlinSolvIters(void* cvode_mem, long int* nniters)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *nniters = cv_mem->cv_nni;

  return (CV_SUCCESS);
}

/*
 * CVodeGetNumNonlinSolvConvFails
 *
 * Returns the current number of convergence failures in the
 * nonlinear solver
 */

int CVodeGetNumNonlinSolvConvFails(void* cvode_mem, long int* nnfails)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *nnfails = cv_mem->cv_nnf;

  return (CV_SUCCESS);
}

/*
 * CVodeGetNonlinSolvStats
 *
 * Returns nonlinear solver statistics
 */

int CVodeGetNonlinSolvStats(void* cvode_mem, long int* nniters, long int* nnfails)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *nniters = cv_mem->cv_nni;
  *nnfails = cv_mem->cv_nnf;

  return (CV_SUCCESS);
}

/*
 * CVodeGetNumStepSolveFails
 *
 * Returns the current number of failed steps due to a nonlinear solver
 * convergence failure
 */

int CVodeGetNumStepSolveFails(void* cvode_mem, long int* nncfails)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *nncfails = cv_mem->cv_ncfn;

  return (CV_SUCCESS);
}

/*
 * =================================================================
 * Quadrature optional output functions
 * =================================================================
 */

/*-----------------------------------------------------------------*/

int CVodeGetQuadNumRhsEvals(void* cvode_mem, long int* nfQevals)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_quadr == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_QUAD, __LINE__, __func__, __FILE__,
                   MSGCV_NO_QUAD);
    return (CV_NO_QUAD);
  }

  *nfQevals = cv_mem->cv_nfQe;

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetQuadNumErrTestFails(void* cvode_mem, long int* nQetfails)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_quadr == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_QUAD, __LINE__, __func__, __FILE__,
                   MSGCV_NO_QUAD);
    return (CV_NO_QUAD);
  }

  *nQetfails = cv_mem->cv_netfQ;

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetQuadErrWeights(void* cvode_mem, N_Vector eQweight)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_quadr == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_QUAD, __LINE__, __func__, __FILE__,
                   MSGCV_NO_QUAD);
    return (CV_NO_QUAD);
  }

  if (cv_mem->cv_errconQ) { N_VScale(ONE, cv_mem->cv_ewtQ, eQweight); }

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetQuadStats(void* cvode_mem, long int* nfQevals, long int* nQetfails)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_quadr == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_QUAD, __LINE__, __func__, __FILE__,
                   MSGCV_NO_QUAD);
    return (CV_NO_QUAD);
  }

  *nfQevals  = cv_mem->cv_nfQe;
  *nQetfails = cv_mem->cv_netfQ;

  return (CV_SUCCESS);
}

/*
 * =================================================================
 * Quadrature FSA optional output functions
 * =================================================================
 */

/*-----------------------------------------------------------------*/

int CVodeGetQuadSensNumRhsEvals(void* cvode_mem, long int* nfQSevals)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_quadr_sensi == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_QUADSENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_QUADSENSI);
    return (CV_NO_QUADSENS);
  }

  *nfQSevals = cv_mem->cv_nfQSe;

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetQuadSensNumErrTestFails(void* cvode_mem, long int* nQSetfails)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_quadr_sensi == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_QUADSENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_QUADSENSI);
    return (CV_NO_QUADSENS);
  }

  *nQSetfails = cv_mem->cv_netfQS;

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetQuadSensErrWeights(void* cvode_mem, N_Vector* eQSweight)
{
  CVodeMem cv_mem;
  int is, Ns;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_quadr_sensi == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_QUADSENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_QUADSENSI);
    return (CV_NO_QUADSENS);
  }
  Ns = cv_mem->cv_Ns;

  if (cv_mem->cv_errconQS)
  {
    for (is = 0; is < Ns; is++)
    {
      N_VScale(ONE, cv_mem->cv_ewtQS[is], eQSweight[is]);
    }
  }

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetQuadSensStats(void* cvode_mem, long int* nfQSevals,
                          long int* nQSetfails)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_quadr_sensi == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_QUADSENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_QUADSENSI);
    return (CV_NO_QUADSENS);
  }

  *nfQSevals  = cv_mem->cv_nfQSe;
  *nQSetfails = cv_mem->cv_netfQS;

  return (CV_SUCCESS);
}

/*
 * =================================================================
 * FSA optional output functions
 * =================================================================
 */

/*-----------------------------------------------------------------*/

int CVodeGetSensNumRhsEvals(void* cvode_mem, long int* nfSevals)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_sensi == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_SENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_SENSI);
    return (CV_NO_SENS);
  }

  *nfSevals = cv_mem->cv_nfSe;

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetNumRhsEvalsSens(void* cvode_mem, long int* nfevalsS)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_sensi == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_SENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_SENSI);
    return (CV_NO_SENS);
  }

  *nfevalsS = cv_mem->cv_nfeS;

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetSensNumErrTestFails(void* cvode_mem, long int* nSetfails)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_sensi == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_SENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_SENSI);
    return (CV_NO_SENS);
  }

  *nSetfails = cv_mem->cv_netfS;

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetSensNumLinSolvSetups(void* cvode_mem, long int* nlinsetupsS)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_sensi == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_SENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_SENSI);
    return (CV_NO_SENS);
  }

  *nlinsetupsS = cv_mem->cv_nsetupsS;

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetSensErrWeights(void* cvode_mem, N_Vector* eSweight)
{
  CVodeMem cv_mem;
  int is, Ns;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_sensi == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_SENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_SENSI);
    return (CV_NO_SENS);
  }

  Ns = cv_mem->cv_Ns;

  for (is = 0; is < Ns; is++)
  {
    N_VScale(ONE, cv_mem->cv_ewtS[is], eSweight[is]);
  }

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetSensStats(void* cvode_mem, long int* nfSevals, long int* nfevalsS,
                      long int* nSetfails, long int* nlinsetupsS)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_sensi == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_SENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_SENSI);
    return (CV_NO_SENS);
  }

  *nfSevals    = cv_mem->cv_nfSe;
  *nfevalsS    = cv_mem->cv_nfeS;
  *nSetfails   = cv_mem->cv_netfS;
  *nlinsetupsS = cv_mem->cv_nsetupsS;

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetSensNumNonlinSolvIters(void* cvode_mem, long int* nSniters)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_sensi == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_SENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_SENSI);
    return (CV_NO_SENS);
  }

  *nSniters = cv_mem->cv_nniS;

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetSensNumNonlinSolvConvFails(void* cvode_mem, long int* nSnfails)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_sensi == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_SENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_SENSI);
    return (CV_NO_SENS);
  }

  *nSnfails = cv_mem->cv_nnfS;

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetSensNonlinSolvStats(void* cvode_mem, long int* nSniters,
                                long int* nSnfails)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_sensi == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_SENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_SENSI);
    return (CV_NO_SENS);
  }

  *nSniters = cv_mem->cv_nniS;
  *nSnfails = cv_mem->cv_nnfS;

  return (CV_SUCCESS);
}

int CVodeGetNumStepSensSolveFails(void* cvode_mem, long int* nSncfails)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (cv_mem->cv_sensi == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_SENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_SENSI);
    return (CV_NO_SENS);
  }

  *nSncfails = cv_mem->cv_ncfnS;

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetStgrSensNumNonlinSolvIters(void* cvode_mem, long int* nSTGR1niters)
{
  CVodeMem cv_mem;
  int is, Ns;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  Ns = cv_mem->cv_Ns;

  if (cv_mem->cv_sensi == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_SENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_SENSI);
    return (CV_NO_SENS);
  }

  if (cv_mem->cv_ism == CV_STAGGERED1)
  {
    for (is = 0; is < Ns; is++) { nSTGR1niters[is] = cv_mem->cv_nniS1[is]; }
  }

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetStgrSensNumNonlinSolvConvFails(void* cvode_mem, long int* nSTGR1nfails)
{
  CVodeMem cv_mem;
  int is, Ns;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  Ns = cv_mem->cv_Ns;

  if (cv_mem->cv_sensi == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_SENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_SENSI);
    return (CV_NO_SENS);
  }

  if (cv_mem->cv_ism == CV_STAGGERED1)
  {
    for (is = 0; is < Ns; is++) { nSTGR1nfails[is] = cv_mem->cv_nnfS1[is]; }
  }

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetStgrSensNonlinSolvStats(void* cvode_mem, long int* nSTGR1niters,
                                    long int* nSTGR1nfails)
{
  CVodeMem cv_mem;
  int is, Ns;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  Ns = cv_mem->cv_Ns;

  if (cv_mem->cv_sensi == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_SENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_SENSI);
    return (CV_NO_SENS);
  }

  if (cv_mem->cv_ism == CV_STAGGERED1)
  {
    for (is = 0; is < Ns; is++) { nSTGR1niters[is] = cv_mem->cv_nniS1[is]; }
    for (is = 0; is < Ns; is++) { nSTGR1nfails[is] = cv_mem->cv_nnfS1[is]; }
  }

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetNumStepStgrSensSolveFails(void* cvode_mem, long int* nSTGR1ncfails)
{
  CVodeMem cv_mem;
  int is, Ns;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  Ns = cv_mem->cv_Ns;

  if (cv_mem->cv_sensi == SUNFALSE)
  {
    cvProcessError(cv_mem, CV_NO_SENS, __LINE__, __func__, __FILE__,
                   MSGCV_NO_SENSI);
    return (CV_NO_SENS);
  }

  if (cv_mem->cv_ism == CV_STAGGERED1)
  {
    for (is = 0; is < Ns; is++) { nSTGR1ncfails[is] = cv_mem->cv_ncfnS1[is]; }
  }

  return (CV_SUCCESS);
}

/*
 * CVodePrintAllStats
 *
 * Print all integrator statistics
 */

int CVodePrintAllStats(void* cvode_mem, FILE* outfile, SUNOutputFormat fmt)
{
  CVodeMem cv_mem;
  CVLsMem cvls_mem;
  CVodeProjMem cvproj_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  if (fmt != SUN_OUTPUTFORMAT_TABLE && fmt != SUN_OUTPUTFORMAT_CSV)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   "Invalid formatting option.");
    return (CV_ILL_INPUT);
  }

  /* step and method stats */
  sunfprintf_real(outfile, fmt, SUNTRUE, "Current time", cv_mem->cv_tn);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Steps", cv_mem->cv_nst);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Error test fails", cv_mem->cv_netf);
  sunfprintf_long(outfile, fmt, SUNFALSE, "NLS step fails", cv_mem->cv_ncfn);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Constraint fails",
                  cv_mem->constraint_fails);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Constraint corrections",
                  cv_mem->constraint_corrections);
  sunfprintf_real(outfile, fmt, SUNFALSE, "Initial step size", cv_mem->cv_h0u);
  sunfprintf_real(outfile, fmt, SUNFALSE, "Last step size", cv_mem->cv_hu);
  sunfprintf_real(outfile, fmt, SUNFALSE, "Current step size", cv_mem->cv_next_h);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Last method order", cv_mem->cv_qu);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Current method order",
                  cv_mem->cv_next_q);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Stab. lim. order reductions",
                  cv_mem->cv_nor);
  /* function evaluations */
  sunfprintf_long(outfile, fmt, SUNFALSE, "RHS fn evals", cv_mem->cv_nfe);
  /* nonlinear solver stats */
  sunfprintf_long(outfile, fmt, SUNFALSE, "NLS iters", cv_mem->cv_nni);
  sunfprintf_long(outfile, fmt, SUNFALSE, "NLS fails", cv_mem->cv_nnf);
  if (cv_mem->cv_nst > 0)
  {
    sunfprintf_real(outfile, fmt, SUNFALSE, "NLS iters per step",
                    (sunrealtype)cv_mem->cv_nni / (sunrealtype)cv_mem->cv_nst);
  }

  /* linear solver stats */
  sunfprintf_long(outfile, fmt, SUNFALSE, "LS setups", cv_mem->cv_nsetups);
  if (cv_mem->cv_lmem)
  {
    cvls_mem = (CVLsMem)(cv_mem->cv_lmem);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Jac fn evals", cvls_mem->nje);
    sunfprintf_long(outfile, fmt, SUNFALSE, "LS RHS fn evals", cvls_mem->nfeDQ);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Prec setup evals", cvls_mem->npe);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Prec solves", cvls_mem->nps);
    sunfprintf_long(outfile, fmt, SUNFALSE, "LS iters", cvls_mem->nli);
    sunfprintf_long(outfile, fmt, SUNFALSE, "LS fails", cvls_mem->ncfl);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Jac-times setups",
                    cvls_mem->njtsetup);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Jac-times evals", cvls_mem->njtimes);
    if (cv_mem->cv_nni > 0)
    {
      sunfprintf_real(outfile, fmt, SUNFALSE, "LS iters per NLS iter",
                      (sunrealtype)cvls_mem->nli / (sunrealtype)cv_mem->cv_nni);
      sunfprintf_real(outfile, fmt, SUNFALSE, "Jac evals per NLS iter",
                      (sunrealtype)cvls_mem->nje / (sunrealtype)cv_mem->cv_nni);
      sunfprintf_real(outfile, fmt, SUNFALSE, "Prec evals per NLS iter",
                      (sunrealtype)cvls_mem->npe / (sunrealtype)cv_mem->cv_nni);
    }
  }

  /* rootfinding stats */
  sunfprintf_long(outfile, fmt, SUNFALSE, "Root fn evals", cv_mem->cv_nge);

  /* projection stats */
  if (cv_mem->proj_mem)
  {
    cvproj_mem = (CVodeProjMem)(cv_mem->proj_mem);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Projection fn evals",
                    cvproj_mem->nproj);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Projection fails",
                    cvproj_mem->npfails);
  }

  /* quadrature stats */
  if (cv_mem->cv_quadr)
  {
    sunfprintf_long(outfile, fmt, SUNFALSE, "Quad fn evals", cv_mem->cv_nfQe);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Quad error test fails",
                    cv_mem->cv_netfQ);
  }

  /* sensitivity stats */
  if (cv_mem->cv_sensi)
  {
    sunfprintf_long(outfile, fmt, SUNFALSE, "Sens fn evals", cv_mem->cv_nfSe);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Sens RHS fn evals", cv_mem->cv_nfeS);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Sens error test fails",
                    cv_mem->cv_netfS);
    if (cv_mem->cv_ism != CV_SIMULTANEOUS)
    {
      sunfprintf_long(outfile, fmt, SUNFALSE, "Sens NLS iters", cv_mem->cv_nniS);
      sunfprintf_long(outfile, fmt, SUNFALSE, "Sens NLS fails", cv_mem->cv_nnfS);
      sunfprintf_long(outfile, fmt, SUNFALSE, "Sens NLS step fails",
                      cv_mem->cv_ncfnS);
    }
    if (cv_mem->cv_ism == CV_STAGGERED1)
    {
      sunfprintf_long_array(outfile, fmt, SUNFALSE, "Sens stgr1 NLS iters",
                            cv_mem->cv_nniS1, cv_mem->cv_Ns);
      sunfprintf_long_array(outfile, fmt, SUNFALSE, "Sens stgr1 NLS fails",
                            cv_mem->cv_nnfS1, cv_mem->cv_Ns);
      sunfprintf_long_array(outfile, fmt, SUNFALSE, "Sens stgr1 NLS step fails",
                            cv_mem->cv_ncfnS1, cv_mem->cv_Ns);
    }
    sunfprintf_long(outfile, fmt, SUNFALSE, "Sens LS setups",
                    cv_mem->cv_nsetupsS);
  }

  /* quadrature-sensitivity stats */
  if (cv_mem->cv_quadr_sensi)
  {
    sunfprintf_long(outfile, fmt, SUNFALSE, "QuadSens fn evals",
                    cv_mem->cv_nfQSe);
    sunfprintf_long(outfile, fmt, SUNFALSE, "QuadSens error test fails",
                    cv_mem->cv_netfQS);
  }

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

int CVodeGetUserData(void* cvode_mem, void** user_data)
{
  CVodeMem cv_mem;

  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return (CV_MEM_NULL);
  }

  cv_mem = (CVodeMem)cvode_mem;

  *user_data = cv_mem->cv_user_data;

  return (CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

char* CVodeGetReturnFlagName(long int flag)
{
  char* name;

  name = (char*)malloc(24 * sizeof(char));

  switch (flag)
  {
  case CV_SUCCESS: sprintf(name, "CV_SUCCESS"); break;
  case CV_TSTOP_RETURN: sprintf(name, "CV_TSTOP_RETURN"); break;
  case CV_ROOT_RETURN: sprintf(name, "CV_ROOT_RETURN"); break;
  case CV_TOO_MUCH_WORK: sprintf(name, "CV_TOO_MUCH_WORK"); break;
  case CV_TOO_MUCH_ACC: sprintf(name, "CV_TOO_MUCH_ACC"); break;
  case CV_ERR_FAILURE: sprintf(name, "CV_ERR_FAILURE"); break;
  case CV_CONV_FAILURE: sprintf(name, "CV_CONV_FAILURE"); break;
  case CV_LINIT_FAIL: sprintf(name, "CV_LINIT_FAIL"); break;
  case CV_LSETUP_FAIL: sprintf(name, "CV_LSETUP_FAIL"); break;
  case CV_LSOLVE_FAIL: sprintf(name, "CV_LSOLVE_FAIL"); break;
  case CV_RHSFUNC_FAIL: sprintf(name, "CV_RHSFUNC_FAIL"); break;
  case CV_FIRST_RHSFUNC_ERR: sprintf(name, "CV_FIRST_RHSFUNC_ERR"); break;
  case CV_REPTD_RHSFUNC_ERR: sprintf(name, "CV_REPTD_RHSFUNC_ERR"); break;
  case CV_UNREC_RHSFUNC_ERR: sprintf(name, "CV_UNREC_RHSFUNC_ERR"); break;
  case CV_RTFUNC_FAIL: sprintf(name, "CV_RTFUNC_FAIL"); break;
  case CV_MEM_FAIL: sprintf(name, "CV_MEM_FAIL"); break;
  case CV_MEM_NULL: sprintf(name, "CV_MEM_NULL"); break;
  case CV_ILL_INPUT: sprintf(name, "CV_ILL_INPUT"); break;
  case CV_NO_MALLOC: sprintf(name, "CV_NO_MALLOC"); break;
  case CV_BAD_K: sprintf(name, "CV_BAD_K"); break;
  case CV_BAD_T: sprintf(name, "CV_BAD_T"); break;
  case CV_BAD_DKY: sprintf(name, "CV_BAD_DKY"); break;
  case CV_NO_QUAD: sprintf(name, "CV_NO_QUAD"); break;
  case CV_QRHSFUNC_FAIL: sprintf(name, "CV_QRHSFUNC_FAIL"); break;
  case CV_FIRST_QRHSFUNC_ERR: sprintf(name, "CV_FIRST_QRHSFUNC_ERR"); break;
  case CV_REPTD_QRHSFUNC_ERR: sprintf(name, "CV_REPTD_QRHSFUNC_ERR"); break;
  case CV_UNREC_QRHSFUNC_ERR: sprintf(name, "CV_UNREC_QRHSFUNC_ERR"); break;
  case CV_BAD_IS: sprintf(name, "CV_BAD_IS"); break;
  case CV_NO_SENS: sprintf(name, "CV_NO_SENS"); break;
  case CV_SRHSFUNC_FAIL: sprintf(name, "CV_SRHSFUNC_FAIL"); break;
  case CV_FIRST_SRHSFUNC_ERR: sprintf(name, "CV_FIRST_SRHSFUNC_ERR"); break;
  case CV_REPTD_SRHSFUNC_ERR: sprintf(name, "CV_REPTD_SRHSFUNC_ERR"); break;
  case CV_UNREC_SRHSFUNC_ERR: sprintf(name, "CV_UNREC_SRHSFUNC_ERR"); break;
  case CV_TOO_CLOSE: sprintf(name, "CV_TOO_CLOSE"); break;
  case CV_NLS_INIT_FAIL: sprintf(name, "CV_NLS_INIT_FAIL"); break;
  case CV_NLS_SETUP_FAIL: sprintf(name, "CV_NLS_SETUPT_FAIL"); break;
  case CV_NO_ADJ: sprintf(name, "CV_NO_ADJ"); break;
  case CV_NO_FWD: sprintf(name, "CV_NO_FWD"); break;
  case CV_NO_BCK: sprintf(name, "CV_NO_BCK"); break;
  case CV_BAD_TB0: sprintf(name, "CV_BAD_TB0"); break;
  case CV_REIFWD_FAIL: sprintf(name, "CV_REIFWD_FAIL"); break;
  case CV_FWD_FAIL: sprintf(name, "CV_FWD_FAIL"); break;
  case CV_GETY_BADT: sprintf(name, "CV_GETY_BADT"); break;
  case CV_NLS_FAIL: sprintf(name, "CV_NLS_FAIL"); break;
  case CV_PROJ_MEM_NULL: sprintf(name, "CV_PROJ_MEM_NULL"); break;
  case CV_PROJFUNC_FAIL: sprintf(name, "CV_PROJFUNC_FAIL"); break;
  case CV_REPTD_PROJFUNC_ERR: sprintf(name, "CV_REPTD_PROJFUNC_ERR"); break;
  default: sprintf(name, "NONE");
  }

  return (name);
}
