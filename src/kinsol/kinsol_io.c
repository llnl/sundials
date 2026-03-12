/* -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier, Shelby Lockhart @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * This is the implementation file for the optional input and output
 * functions for the KINSOL solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "kinsol_impl.h"
#include "kinsol_ls_impl.h"
#include "sundials_utils.h"

#define ZERO      SUN_RCONST(0.0)
#define POINT1    SUN_RCONST(0.1)
#define ONETHIRD  SUN_RCONST(0.3333333333333333)
#define TWOTHIRDS SUN_RCONST(0.6666666666666667)
#define POINT9    SUN_RCONST(0.9)
#define ONE       SUN_RCONST(1.0)
#define TWO       SUN_RCONST(2.0)
#define TWOPT5    SUN_RCONST(2.5)

/*
 * =================================================================
 * KINSOL optional input functions
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Function : KINSetUserData
 * -----------------------------------------------------------------
 */

int KINSetUserData(void* kinmem, void* user_data)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem                = (KINMem)kinmem;
  kin_mem->kin_user_data = user_data;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetDamping
 * -----------------------------------------------------------------
 */

int KINSetDamping(void* kinmem, sunrealtype beta)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  /* check for illegal input value */
  if (beta <= ZERO)
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "beta <= 0 illegal");
    return (KIN_ILL_INPUT);
  }

  if (beta < ONE)
  {
    /* enable damping */
    kin_mem->kin_beta    = beta;
    kin_mem->kin_damping = SUNTRUE;
  }
  else
  {
    /* disable damping */
    kin_mem->kin_beta    = ONE;
    kin_mem->kin_damping = SUNFALSE;
  }

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetMAA
 * -----------------------------------------------------------------
 */

int KINSetMAA(void* kinmem, long int maa)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  if (maa < 0)
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_BAD_MAA);
    return (KIN_ILL_INPUT);
  }

  // To allow for setting the depth and max number of iterations in any order we
  // do not limit maa here and instead enforce maa < mxiter in the AA
  // initialization function (KINInitAA)
  kin_mem->kin_m_aa = maa;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetDelayAA
 * -----------------------------------------------------------------
 */

int KINSetDelayAA(void* kinmem, long int delay)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  /* check for illegal input value */
  if (delay < 0)
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "delay < 0 illegal");
    return (KIN_ILL_INPUT);
  }

  kin_mem->kin_delay_aa = delay;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetOrthAA
 * -----------------------------------------------------------------
 */

int KINSetOrthAA(void* kinmem, int orthaa)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  if ((orthaa < KIN_ORTH_MGS) || (orthaa > KIN_ORTH_DCGS2))
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_BAD_ORTHAA);
    return (KIN_ILL_INPUT);
  }

  kin_mem->kin_orth_aa = orthaa;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetDampingAA
 * -----------------------------------------------------------------
 */

int KINSetDampingAA(void* kinmem, sunrealtype beta)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  /* check for illegal input value */
  if (beta <= ZERO)
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "beta <= 0 illegal");
    return (KIN_ILL_INPUT);
  }

  if (beta < ONE)
  {
    /* enable damping */
    kin_mem->kin_beta_aa    = beta;
    kin_mem->kin_damping_aa = SUNTRUE;
  }
  else
  {
    /* disable damping */
    kin_mem->kin_beta_aa    = ONE;
    kin_mem->kin_damping_aa = SUNFALSE;
  }

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetDampingFn
 * -----------------------------------------------------------------
 */

int KINSetDampingFn(void* kinmem, KINDampingFn damping_fn)
{
  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return KIN_MEM_NULL;
  }
  KINMem kin_mem = (KINMem)kinmem;

  kin_mem->kin_damping_fn = damping_fn;

  return KIN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetDepthFn
 * -----------------------------------------------------------------
 */

int KINSetDepthFn(void* kinmem, KINDepthFn depth_fn)
{
  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return KIN_MEM_NULL;
  }
  KINMem kin_mem = (KINMem)kinmem;

  kin_mem->kin_depth_fn = depth_fn;

  return KIN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetReturnNewest
 * -----------------------------------------------------------------
 */

int KINSetReturnNewest(void* kinmem, sunbooleantype ret_newest)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  kin_mem->kin_ret_newest = ret_newest;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetNumMaxIters
 * -----------------------------------------------------------------
 */

int KINSetNumMaxIters(void* kinmem, long int mxiter)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  if (mxiter < 0)
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_BAD_MXITER);
    return (KIN_ILL_INPUT);
  }

  if (mxiter == 0) { kin_mem->kin_mxiter = MXITER_DEFAULT; }
  else { kin_mem->kin_mxiter = mxiter; }

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetNoInitSetup
 * -----------------------------------------------------------------
 */

int KINSetNoInitSetup(void* kinmem, sunbooleantype noInitSetup)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem                  = (KINMem)kinmem;
  kin_mem->kin_noInitSetup = noInitSetup;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetNoResMon
 * -----------------------------------------------------------------
 */

int KINSetNoResMon(void* kinmem, sunbooleantype noResMon)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem               = (KINMem)kinmem;
  kin_mem->kin_noResMon = noResMon;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetMaxSetupCalls
 * -----------------------------------------------------------------
 */

int KINSetMaxSetupCalls(void* kinmem, long int msbset)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  if (msbset < 0)
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_BAD_MSBSET);
    return (KIN_ILL_INPUT);
  }

  if (msbset == 0) { kin_mem->kin_msbset = MSBSET_DEFAULT; }
  else { kin_mem->kin_msbset = msbset; }

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetMaxSubSetupCalls
 * -----------------------------------------------------------------
 */

int KINSetMaxSubSetupCalls(void* kinmem, long int msbsetsub)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  if (msbsetsub < 0)
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_BAD_MSBSETSUB);
    return (KIN_ILL_INPUT);
  }

  if (msbsetsub == 0) { kin_mem->kin_msbset_sub = MSBSET_SUB_DEFAULT; }
  else { kin_mem->kin_msbset_sub = msbsetsub; }

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetEtaForm
 * -----------------------------------------------------------------
 */

int KINSetEtaForm(void* kinmem, int etachoice)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  if ((etachoice != KIN_ETACONSTANT) && (etachoice != KIN_ETACHOICE1) &&
      (etachoice != KIN_ETACHOICE2))
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_BAD_ETACHOICE);
    return (KIN_ILL_INPUT);
  }

  kin_mem->kin_etaflag = etachoice;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetEtaConstValue
 * -----------------------------------------------------------------
 */

int KINSetEtaConstValue(void* kinmem, sunrealtype eta)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  if ((eta < ZERO) || (eta > ONE))
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_BAD_ETACONST);
    return (KIN_ILL_INPUT);
  }

  if (eta == ZERO) { kin_mem->kin_eta = POINT1; }
  else { kin_mem->kin_eta = eta; }

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetEtaParams
 * -----------------------------------------------------------------
 */

int KINSetEtaParams(void* kinmem, sunrealtype egamma, sunrealtype ealpha)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  if ((ealpha <= ONE) || (ealpha > TWO))
  {
    if (ealpha != ZERO)
    {
      KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                      MSG_BAD_ALPHA);
      return (KIN_ILL_INPUT);
    }
  }

  if (ealpha == ZERO) { kin_mem->kin_eta_alpha = TWO; }
  else { kin_mem->kin_eta_alpha = ealpha; }

  if ((egamma <= ZERO) || (egamma > ONE))
  {
    if (egamma != ZERO)
    {
      KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                      MSG_BAD_GAMMA);
      return (KIN_ILL_INPUT);
    }
  }

  if (egamma == ZERO) { kin_mem->kin_eta_gamma = POINT9; }
  else { kin_mem->kin_eta_gamma = egamma; }

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetResMonParams
 * -----------------------------------------------------------------
 */

int KINSetResMonParams(void* kinmem, sunrealtype omegamin, sunrealtype omegamax)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  /* check omegamin */

  if (omegamin < ZERO)
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_BAD_OMEGA);
    return (KIN_ILL_INPUT);
  }

  if (omegamin == ZERO) { kin_mem->kin_omega_min = OMEGA_MIN; }
  else { kin_mem->kin_omega_min = omegamin; }

  /* check omegamax */

  if (omegamax < ZERO)
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_BAD_OMEGA);
    return (KIN_ILL_INPUT);
  }

  if (omegamax == ZERO)
  {
    if (kin_mem->kin_omega_min > OMEGA_MAX)
    {
      KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                      MSG_BAD_OMEGA);
      return (KIN_ILL_INPUT);
    }
    else { kin_mem->kin_omega_max = OMEGA_MAX; }
  }
  else
  {
    if (kin_mem->kin_omega_min > omegamax)
    {
      KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                      MSG_BAD_OMEGA);
      return (KIN_ILL_INPUT);
    }
    else { kin_mem->kin_omega_max = omegamax; }
  }

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetResMonConstValue
 * -----------------------------------------------------------------
 */

int KINSetResMonConstValue(void* kinmem, sunrealtype omegaconst)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  /* check omegaconst */

  if (omegaconst < ZERO)
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_BAD_OMEGA);
    return (KIN_ILL_INPUT);
  }

  /* Load omega value. A value of 0 will force using omega_min and omega_max */
  kin_mem->kin_omega = omegaconst;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetNoMinEps
 * -----------------------------------------------------------------
 */

int KINSetNoMinEps(void* kinmem, sunbooleantype noMinEps)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem               = (KINMem)kinmem;
  kin_mem->kin_noMinEps = noMinEps;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetMaxNewtonStep
 * -----------------------------------------------------------------
 */

int KINSetMaxNewtonStep(void* kinmem, sunrealtype mxnewtstep)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  if (mxnewtstep < ZERO)
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_BAD_MXNEWTSTEP);
    return (KIN_ILL_INPUT);
  }

  /* Note: passing a value of 0.0 will use the default
     value (computed in KINSolInit) */

  kin_mem->kin_mxnstepin = mxnewtstep;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetMaxBetaFails
 * -----------------------------------------------------------------
 */

int KINSetMaxBetaFails(void* kinmem, long int mxnbcf)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  if (mxnbcf < 0)
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_BAD_MXNBCF);
    return (KIN_ILL_INPUT);
  }

  if (mxnbcf == 0) { kin_mem->kin_mxnbcf = MXNBCF_DEFAULT; }
  else { kin_mem->kin_mxnbcf = mxnbcf; }

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetRelErrFunc
 * -----------------------------------------------------------------
 */

int KINSetRelErrFunc(void* kinmem, sunrealtype relfunc)
{
  KINMem kin_mem;
  sunrealtype uround;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  if (relfunc < ZERO)
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_BAD_RELFUNC);
    return (KIN_ILL_INPUT);
  }

  if (relfunc == ZERO)
  {
    uround                    = kin_mem->kin_uround;
    kin_mem->kin_sqrt_relfunc = SUNRsqrt(uround);
  }
  else { kin_mem->kin_sqrt_relfunc = SUNRsqrt(relfunc); }

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetFuncNormTol
 * -----------------------------------------------------------------
 */

int KINSetFuncNormTol(void* kinmem, sunrealtype fnormtol)
{
  KINMem kin_mem;
  sunrealtype uround;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  if (fnormtol < ZERO)
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_BAD_FNORMTOL);
    return (KIN_ILL_INPUT);
  }

  if (fnormtol == ZERO)
  {
    uround                = kin_mem->kin_uround;
    kin_mem->kin_fnormtol = SUNRpowerR(uround, ONETHIRD);
  }
  else { kin_mem->kin_fnormtol = fnormtol; }

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetScaledStepTol
 * -----------------------------------------------------------------
 */

int KINSetScaledStepTol(void* kinmem, sunrealtype scsteptol)
{
  KINMem kin_mem;
  sunrealtype uround;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  if (scsteptol < ZERO)
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_BAD_SCSTEPTOL);
    return (KIN_ILL_INPUT);
  }

  if (scsteptol == ZERO)
  {
    uround                 = kin_mem->kin_uround;
    kin_mem->kin_scsteptol = SUNRpowerR(uround, TWOTHIRDS);
  }
  else { kin_mem->kin_scsteptol = scsteptol; }

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetConstraints
 * -----------------------------------------------------------------
 */

int KINSetConstraints(void* kinmem, N_Vector constraints)
{
  KINMem kin_mem;
  sunrealtype temptest;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  if (constraints == NULL)
  {
    if (kin_mem->kin_constraintsSet)
    {
      N_VDestroy(kin_mem->kin_constraints);
      kin_mem->kin_lrw -= kin_mem->kin_lrw1;
      kin_mem->kin_liw -= kin_mem->kin_liw1;
    }
    kin_mem->kin_constraintsSet = SUNFALSE;
    return (KIN_SUCCESS);
  }

  /* Check the constraints vector */

  temptest = N_VMaxNorm(constraints);
  if (temptest > TWOPT5)
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_BAD_CONSTRAINTS);
    return (KIN_ILL_INPUT);
  }

  if (!kin_mem->kin_constraintsSet)
  {
    kin_mem->kin_constraints = N_VClone(constraints);
    kin_mem->kin_lrw += kin_mem->kin_lrw1;
    kin_mem->kin_liw += kin_mem->kin_liw1;
    kin_mem->kin_constraintsSet = SUNTRUE;
  }

  /* Load the constraint vector */

  N_VScale(ONE, constraints, kin_mem->kin_constraints);

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetSysFunc
 * -----------------------------------------------------------------
 */

int KINSetSysFunc(void* kinmem, KINSysFn func)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  if (func == NULL)
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_FUNC_NULL);
    return (KIN_ILL_INPUT);
  }

  kin_mem->kin_func = func;

  return (KIN_SUCCESS);
}

/*
 * =================================================================
 * KINSOL optional output functions
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Function : KINGetWorkSpace
 * -----------------------------------------------------------------
 */

int KINGetWorkSpace(void* kinmem, long int* lenrw, long int* leniw)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  *lenrw = kin_mem->kin_lrw;
  *leniw = kin_mem->kin_liw;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetNumNonlinSolvIters
 * -----------------------------------------------------------------
 */

int KINGetNumNonlinSolvIters(void* kinmem, long int* nniters)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem  = (KINMem)kinmem;
  *nniters = kin_mem->kin_nni;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetNumFuncEvals
 * -----------------------------------------------------------------
 */

int KINGetNumFuncEvals(void* kinmem, long int* nfevals)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem  = (KINMem)kinmem;
  *nfevals = kin_mem->kin_nfe;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetNumBetaCondFails
 * -----------------------------------------------------------------
 */

int KINGetNumBetaCondFails(void* kinmem, long int* nbcfails)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem   = (KINMem)kinmem;
  *nbcfails = kin_mem->kin_nbcf;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetNumBacktrackOps
 * -----------------------------------------------------------------
 */

int KINGetNumBacktrackOps(void* kinmem, long int* nbacktr)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem  = (KINMem)kinmem;
  *nbacktr = kin_mem->kin_nbktrk;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetFuncNorm
 * -----------------------------------------------------------------
 */

int KINGetFuncNorm(void* kinmem, sunrealtype* funcnorm)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem   = (KINMem)kinmem;
  *funcnorm = kin_mem->kin_fnorm;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetStepLength
 * -----------------------------------------------------------------
 */

int KINGetStepLength(void* kinmem, sunrealtype* steplength)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem     = (KINMem)kinmem;
  *steplength = kin_mem->kin_stepl;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetUserData
 * -----------------------------------------------------------------
 */

int KINGetUserData(void* kinmem, void** user_data)
{
  KINMem kin_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  *user_data = kin_mem->kin_user_data;

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINPrintAllStats
 * -----------------------------------------------------------------
 */

int KINPrintAllStats(void* kinmem, FILE* outfile, SUNOutputFormat fmt)
{
  KINMem kin_mem;
  KINLsMem kinls_mem;

  if (kinmem == NULL)
  {
    KINProcessError(NULL, KIN_MEM_NULL, __LINE__, __func__, __FILE__, MSG_NO_MEM);
    return (KIN_MEM_NULL);
  }

  kin_mem = (KINMem)kinmem;

  if (fmt != SUN_OUTPUTFORMAT_TABLE && fmt != SUN_OUTPUTFORMAT_CSV)
  {
    KINProcessError(kin_mem, KIN_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid formatting option.");
    return (KIN_ILL_INPUT);
  }

  sunfprintf_long(outfile, fmt, SUNTRUE, "Nonlinear iters", kin_mem->kin_nni);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Nonlinear fn evals", kin_mem->kin_nfe);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Beta condition fails",
                  kin_mem->kin_nbcf);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Backtrack operations",
                  kin_mem->kin_nbktrk);
  sunfprintf_real(outfile, fmt, SUNFALSE, "Nonlinear fn norm",
                  kin_mem->kin_fnorm);
  sunfprintf_real(outfile, fmt, SUNFALSE, "Step length", kin_mem->kin_stepl);

  /* linear solver stats */
  if (kin_mem->kin_lmem)
  {
    kinls_mem = (KINLsMem)(kin_mem->kin_lmem);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Jac fn evals", kinls_mem->nje);
    sunfprintf_long(outfile, fmt, SUNFALSE, "LS Nonlinear fn evals",
                    kinls_mem->nfeDQ);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Prec setup evals", kinls_mem->npe);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Prec solves", kinls_mem->nps);
    sunfprintf_long(outfile, fmt, SUNFALSE, "LS iters", kinls_mem->nli);
    sunfprintf_long(outfile, fmt, SUNFALSE, "LS fails", kinls_mem->ncfl);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Jac-times evals",
                    kinls_mem->njtimes);
    if (kin_mem->kin_nni > 0)
    {
      sunfprintf_real(outfile, fmt, SUNFALSE, "LS iters per NLS iter",
                      (sunrealtype)kinls_mem->nli / (sunrealtype)kin_mem->kin_nni);
      sunfprintf_real(outfile, fmt, SUNFALSE, "Jac evals per NLS iter",
                      (sunrealtype)kinls_mem->nje / (sunrealtype)kin_mem->kin_nni);
      sunfprintf_real(outfile, fmt, SUNFALSE, "Prec evals per NLS iter",
                      (sunrealtype)kinls_mem->npe / (sunrealtype)kin_mem->kin_nni);
    }
  }

  return (KIN_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetReturnFlagName
 * -----------------------------------------------------------------
 */

char* KINGetReturnFlagName(long int flag)
{
  char* name;

  name = (char*)malloc(24 * sizeof(char));

  switch (flag)
  {
  case KIN_SUCCESS: sprintf(name, "KIN_SUCCESS"); break;
  case KIN_INITIAL_GUESS_OK: sprintf(name, "KIN_INITIAL_GUESS_OK"); break;
  case KIN_STEP_LT_STPTOL: sprintf(name, "KIN_STEP_LT_STPTOL"); break;
  case KIN_WARNING: sprintf(name, "KIN_WARNING"); break;
  case KIN_MEM_NULL: sprintf(name, "KIN_MEM_NULL"); break;
  case KIN_ILL_INPUT: sprintf(name, "KIN_ILL_INPUT"); break;
  case KIN_NO_MALLOC: sprintf(name, "KIN_NO_MALLOC"); break;
  case KIN_MEM_FAIL: sprintf(name, "KIN_MEM_FAIL"); break;
  case KIN_LINESEARCH_NONCONV: sprintf(name, "KIN_LINESEARCH_NONCONV"); break;
  case KIN_MAXITER_REACHED: sprintf(name, "KIN_MAXITER_REACHED"); break;
  case KIN_MXNEWT_5X_EXCEEDED: sprintf(name, "KIN_MXNEWT_5X_EXCEEDED"); break;
  case KIN_LINESEARCH_BCFAIL: sprintf(name, "KIN_LINESEARCH_BCFAIL"); break;
  case KIN_LINSOLV_NO_RECOVERY: sprintf(name, "KIN_LINSOLV_NO_RECOVERY"); break;
  case KIN_LINIT_FAIL: sprintf(name, "KIN_LINIT_FAIL"); break;
  case KIN_LSETUP_FAIL: sprintf(name, "KIN_LSETUP_FAIL"); break;
  case KIN_LSOLVE_FAIL: sprintf(name, "KIN_LSOLVE_FAIL"); break;
  case KIN_SYSFUNC_FAIL: sprintf(name, "KIN_SYSFUNC_FAIL"); break;
  case KIN_FIRST_SYSFUNC_ERR: sprintf(name, "KIN_FIRST_SYSFUNC_ERR"); break;
  case KIN_REPTD_SYSFUNC_ERR: sprintf(name, "KIN_REPTD_SYSFUNC_ERR"); break;
  case KIN_VECTOROP_ERR: sprintf(name, "KIN_VECTOROP_ERR"); break;
  case KIN_CONTEXT_ERR: sprintf(name, "KIN_CONTEXT_ERR"); break;
  case KIN_DAMPING_FN_ERR: sprintf(name, "KIN_DAMPING_FN_ERR"); break;
  default: sprintf(name, "NONE");
  }

  return (name);
}
