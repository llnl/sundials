/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
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
 * -----------------------------------------------------------------------------
 * Build Nordsieck array from solution history
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_types.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#include "cvodes/cvodes.h"
#include "cvodes_impl.h"

#define ZERO SUN_RCONST(0.0) /* real 0.0 */
#define ONE  SUN_RCONST(1.0) /* real 1.0 */

/* -----------------------------------------------------------------------------
 * Build Adams Nordsieck array from f(t,y) history and y(t) value
 * ---------------------------------------------------------------------------*/

static int cvBuildNordsieckArrayAdams(sunrealtype* t, N_Vector y, N_Vector* f,
                                      N_Vector* wrk, int order,
                                      sunrealtype hscale, N_Vector* zn)
{
  /* Check for valid inputs */
  if (!t || !y || !f || !wrk || order < 1 || !zn) { return CV_ILL_INPUT; }

  for (int i = 0; i < order; i++)
  {
    if (!f[i]) { return CV_ILL_INPUT; }
    if (!wrk[i]) { return CV_ILL_INPUT; }
  }

  /* Compute Nordsieck array */
  if (order > 1)
  {
    /* Compute Newton polynomial coefficients interpolating f history */
    for (int i = 0; i < order; i++) { N_VScale(ONE, f[i], wrk[i]); }

    for (int i = 1; i < order; i++)
    {
      for (int j = order - 1; j >= i; j--)
      {
        /* Divided difference */
        sunrealtype delta_t = ONE / (t[j - i] - t[j]);
        N_VLinearSum(delta_t, wrk[j - 1], -delta_t, wrk[j], wrk[j]);
      }
    }

    /* Compute derivatives of Newton polynomial of f history */
    N_VScale(ONE, wrk[order - 1], zn[1]);
    for (int i = 2; i <= order; i++) { N_VConst(ZERO, zn[i]); }

    for (int i = order - 2; i >= 0; i--)
    {
      for (int j = order - 1; j > 0; j--)
      {
        N_VLinearSum(t[0] - t[i], zn[j + 1], j, zn[j], zn[j + 1]);
      }
      N_VLinearSum(t[0] - t[i], zn[1], ONE, wrk[i], zn[1]);
    }
  }

  /* Overwrite first two columns with input values */
  N_VScale(ONE, y, zn[0]);
  N_VScale(ONE, f[0], zn[1]);

  /* Scale entries */
  sunrealtype scale = ONE;
  for (int i = 1; i <= order; i++)
  {
    scale *= hscale / ((sunrealtype)i);
    N_VScale(scale, zn[i], zn[i]);
  }

  return CV_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Build BDF Nordsieck array from y(t) history and f(t,y) value
 * ---------------------------------------------------------------------------*/

static int cvBuildNordsieckArrayBDF(sunrealtype* t, N_Vector* y, N_Vector f,
                                    N_Vector* wrk, int order,
                                    sunrealtype hscale, N_Vector* zn)
{
  /* Check for valid inputs */
  if (!t || !y || !f || !wrk || order < 1 || !zn) { return CV_ILL_INPUT; }

  for (int i = 0; i < order; i++)
  {
    if (!y[i]) { return CV_ILL_INPUT; }
  }

  for (int i = 0; i < order + 1; i++)
  {
    if (!wrk[i]) { return CV_ILL_INPUT; }
  }

  /* Compute Nordsieck array */
  if (order > 1)
  {
    /* Setup extended array of times to incorporate derivative value */
    sunrealtype t_ext[BDF_Q_MAX + 1];

    t_ext[0] = t[0];
    for (int i = 1; i <= order; i++) { t_ext[i] = t[i - 1]; }

    /* Compute Hermite polynomial coefficients interpolating y history and f */
    N_VScale(ONE, y[0], wrk[0]);
    for (int i = 1; i <= order; i++) { N_VScale(ONE, y[i - 1], wrk[i]); }

    for (int i = 1; i <= order; i++)
    {
      for (int j = order; j > i - 1; j--)
      {
        if (i == 1 && j == 1)
        {
          /* Replace with actual derivative value */
          N_VScale(ONE, f, wrk[j]);
        }
        else
        {
          /* Divided difference */
          sunrealtype delta_t = ONE / (t_ext[j - i] - t_ext[j]);
          N_VLinearSum(delta_t, wrk[j - 1], -delta_t, wrk[j], wrk[j]);
        }
      }
    }

    /* Compute derivatives of Hermite polynomial */
    N_VScale(ONE, wrk[order], zn[0]);
    for (int i = 1; i <= order; i++) { N_VConst(ZERO, zn[i]); }

    for (int i = order - 1; i >= 0; i--)
    {
      for (int j = order; j > 0; j--)
      {
        N_VLinearSum(t_ext[0] - t_ext[i], zn[j], j, zn[j - 1], zn[j]);
      }
      N_VLinearSum(t_ext[0] - t_ext[i], zn[0], ONE, wrk[i], zn[0]);
    }
  }

  /* Overwrite first two columns with input values */
  N_VScale(ONE, y[0], zn[0]);
  N_VScale(ONE, f, zn[1]);

  /* Scale entries */
  sunrealtype scale = ONE;
  for (int i = 1; i <= order; i++)
  {
    scale *= hscale / ((sunrealtype)i);
    N_VScale(scale, zn[i], zn[i]);
  }

  return CV_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Compute predicted new state (simplified cvPredict for k = 1 and j = q...1)
 * ---------------------------------------------------------------------------*/

static int cvPredictY(int order, N_Vector* zn, N_Vector ypred)
{
  N_VScale(ONE, zn[0], ypred);
  for (int j = 1; j <= order; j++)
  {
    N_VLinearSum(ONE, zn[j], ONE, ypred, ypred);
  }

  return CV_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Resize CVODE and build new history array
 * ---------------------------------------------------------------------------*/

int CVodeResizeHistory(void* cvode_mem, sunrealtype* t_hist, N_Vector* y_hist,
                       N_Vector* f_hist, int num_y_hist, int num_f_hist)
{
  int retval = 0;

  /* ------------ *
   * Check inputs *
   * ------------ */

  if (!cvode_mem)
  {
    cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
    return CV_MEM_NULL;
  }
  CVodeMem cv_mem = (CVodeMem)cvode_mem;

  if (!t_hist)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   "Time history array is NULL");
    return CV_ILL_INPUT;
  }

  if (!y_hist)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   "State history array is NULL");
    return CV_ILL_INPUT;
  }

  if (!f_hist)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                   "RHS history array is NULL");
    return CV_ILL_INPUT;
  }

  /* Check that the input history is sufficient for the current (next) order */
  int n_hist = SUNMIN(cv_mem->cv_q + 1, cv_mem->cv_qmax);

  if (cv_mem->cv_lmm == CV_ADAMS)
  {
    if (num_y_hist < 2)
    {
      cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Insufficient solution history");
      return CV_ILL_INPUT;
    }

    for (int i = 0; i < n_hist; i++)
    {
      if (!f_hist[i])
      {
        cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                       "Insufficient right-hand side history");
        return CV_ILL_INPUT;
      }
    }
  }
  else
  {
    if (num_f_hist < 2)
    {
      cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                     "Insufficient right-hand side history");
      return CV_ILL_INPUT;
    }

    for (int i = 0; i < n_hist; i++)
    {
      if (!y_hist[i])
      {
        cvProcessError(cv_mem, CV_ILL_INPUT, __LINE__, __func__, __FILE__,
                       "Insufficient solution history");
        return CV_ILL_INPUT;
      }
    }
  }

  /* -------------- *
   * Resize vectors *
   * -------------- */

  N_VDestroy(cv_mem->cv_ewt);
  cv_mem->cv_ewt = N_VClone(y_hist[0]);
  if (!(cv_mem->cv_ewt))
  {
    cvProcessError(cv_mem, CV_MEM_FAIL, __LINE__, __func__, __FILE__,
                   "A vector allocation failed");
    return CV_MEM_FAIL;
  }

  N_VDestroy(cv_mem->cv_acor);
  cv_mem->cv_acor = N_VClone(y_hist[0]);
  if (!(cv_mem->cv_acor))
  {
    cvProcessError(cv_mem, CV_MEM_FAIL, __LINE__, __func__, __FILE__,
                   "A vector allocation failed");
    return CV_MEM_FAIL;
  }

  N_VDestroy(cv_mem->cv_tempv);
  cv_mem->cv_tempv = N_VClone(y_hist[0]);
  if (!(cv_mem->cv_tempv))
  {
    cvProcessError(cv_mem, CV_MEM_FAIL, __LINE__, __func__, __FILE__,
                   "A vector allocation failed");
    return CV_MEM_FAIL;
  }

  N_VDestroy(cv_mem->cv_ftemp);
  cv_mem->cv_ftemp = N_VClone(y_hist[0]);
  if (!(cv_mem->cv_ftemp))
  {
    cvProcessError(cv_mem, CV_MEM_FAIL, __LINE__, __func__, __FILE__,
                   "A vector allocation failed");
    return CV_MEM_FAIL;
  }

  N_VDestroy(cv_mem->cv_vtemp1);
  cv_mem->cv_vtemp1 = N_VClone(y_hist[0]);
  if (!(cv_mem->cv_vtemp1))
  {
    cvProcessError(cv_mem, CV_MEM_FAIL, __LINE__, __func__, __FILE__,
                   "A vector allocation failed");
    return CV_MEM_FAIL;
  }

  N_VDestroy(cv_mem->cv_vtemp2);
  cv_mem->cv_vtemp2 = N_VClone(y_hist[0]);
  if (!(cv_mem->cv_vtemp2))
  {
    cvProcessError(cv_mem, CV_MEM_FAIL, __LINE__, __func__, __FILE__,
                   "A vector allocation failed");
    return CV_MEM_FAIL;
  }

  N_VDestroy(cv_mem->cv_vtemp3);
  cv_mem->cv_vtemp3 = N_VClone(y_hist[0]);
  if (!(cv_mem->cv_vtemp3))
  {
    cvProcessError(cv_mem, CV_MEM_FAIL, __LINE__, __func__, __FILE__,
                   "A vector allocation failed");
    return CV_MEM_FAIL;
  }

  /* User will need to set a new vector of absolute tolerances */
  if (cv_mem->cv_VabstolMallocDone)
  {
    N_VDestroy(cv_mem->cv_Vabstol);
    cv_mem->cv_Vabstol = N_VClone(y_hist[0]);
  }

  /* User will need to set a new constraints vector */
  if (cv_mem->cv_constraints)
  {
    N_VDestroy(cv_mem->cv_constraints);
    cv_mem->cv_constraints = NULL;
  }

  for (int j = 0; j <= cv_mem->cv_qmax_alloc; j++)
  {
    N_VDestroy(cv_mem->cv_zn[j]);
    cv_mem->cv_zn[j] = N_VClone(y_hist[0]);
    if (!(cv_mem->cv_zn[j]))
    {
      cvProcessError(cv_mem, CV_MEM_FAIL, __LINE__, __func__, __FILE__,
                     "A vector allocation failed");
      return CV_MEM_FAIL;
    }
  }

  /* ----------------------- *
   * Resize nonlinear solver *
   * ----------------------- */

  if (cv_mem->NLS && cv_mem->ownNLS)
  {
    retval = SUNNonlinSolFree(cv_mem->NLS);
    if (retval)
    {
      cvProcessError(cv_mem, CV_MEM_FAIL, __LINE__, __func__, __FILE__,
                     "Destroying the Newton solver failed");
      return CV_MEM_FAIL;
    }
    cv_mem->NLS    = NULL;
    cv_mem->ownNLS = SUNFALSE;

    SUNNonlinearSolver NLS = SUNNonlinSol_Newton(y_hist[0], cv_mem->cv_sunctx);
    if (!NLS)
    {
      cvProcessError(cv_mem, CV_MEM_FAIL, __LINE__, __func__, __FILE__,
                     "Error creating the Newton solver");
      return CV_MEM_FAIL;
    }

    retval = CVodeSetNonlinearSolver(cv_mem, NLS);
    if (retval)
    {
      SUNNonlinSolFree(NLS);
      cvProcessError(cv_mem, CV_MEM_FAIL, __LINE__, __func__, __FILE__,
                     "Error attaching default Newton solver");
      return CV_MEM_FAIL;
    }
    cv_mem->ownNLS = SUNTRUE;
  }

  /* ----------------------------- *
   * Create workspace for resizing *
   * ----------------------------- */

  N_Vector resize_wrk[L_MAX];

  int wrk_space_size = SUNMAX(cv_mem->cv_q, cv_mem->cv_qprime);
  if (cv_mem->cv_lmm == CV_BDF) { wrk_space_size++; }

  for (int j = 0; j < wrk_space_size; j++)
  {
    resize_wrk[j] = N_VClone(y_hist[0]);
    if (!resize_wrk[j])
    {
      for (int i = 0; i < j; i++) { N_VDestroy(resize_wrk[i]); }
      cvProcessError(cv_mem, CV_MEM_FAIL, __LINE__, __func__, __FILE__,
                     "A vector allocation failed");
      return CV_MEM_FAIL;
    }
  }

  /* ------------------------------------------------------------------------ *
   * Construct Nordsieck array at the old time but with the new size to
   * compute correction vector at the new state size.
   * ------------------------------------------------------------------------ */

  if (cv_mem->cv_q < cv_mem->cv_qmax)
  {
    /* Compute z_{n-1} with new history size */
    if (cv_mem->cv_lmm == CV_ADAMS)
    {
      retval = cvBuildNordsieckArrayAdams(t_hist + 1, y_hist[1], f_hist + 1,
                                          resize_wrk, cv_mem->cv_q,
                                          cv_mem->cv_hscale, cv_mem->cv_zn);
    }
    else
    {
      retval = cvBuildNordsieckArrayBDF(t_hist + 1, y_hist + 1, f_hist[1],
                                        resize_wrk, cv_mem->cv_q,
                                        cv_mem->cv_hscale, cv_mem->cv_zn);
    }

    if (retval)
    {
      cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
                     "Building the Nordsieck array failed");
      return retval;
    }

    /* Get predicted value */
    retval = cvPredictY(cv_mem->cv_q, cv_mem->cv_zn, cv_mem->cv_vtemp1);

    if (retval)
    {
      cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
                     "Computing the predictor failed");
      return retval;
    }

    /* Resized correction */
    N_VLinearSum(ONE, y_hist[0], -ONE, cv_mem->cv_vtemp1,
                 cv_mem->cv_zn[cv_mem->cv_qmax]);
  }

  /* ----------------------------- *
   * Construct new Nordsieck Array *
   * ----------------------------- */

  if (cv_mem->cv_lmm == CV_ADAMS)
  {
    retval = cvBuildNordsieckArrayAdams(t_hist, y_hist[0], f_hist, resize_wrk,
                                        cv_mem->cv_qprime, cv_mem->cv_hscale,
                                        cv_mem->cv_zn);
  }
  else
  {
    retval = cvBuildNordsieckArrayBDF(t_hist, y_hist, f_hist[0], resize_wrk,
                                      cv_mem->cv_qprime, cv_mem->cv_hscale,
                                      cv_mem->cv_zn);
  }

  if (retval)
  {
    cvProcessError(cv_mem, retval, __LINE__, __func__, __FILE__,
                   "Building the Nordsieck array failed");
    return retval;
  }

  /* ------------------- *
   * Update time history *
   * ------------------- */

  /* Ensure internal time and step history match the input history */
  cv_mem->cv_tn = t_hist[0];

  for (int i = 1; i < n_hist; i++)
  {
    cv_mem->cv_tau[i] = t_hist[i - 1] - t_hist[i];
  }

  /* In the next step, perform initialization needed after a resize */
  cv_mem->first_step_after_resize = SUNTRUE;

  /* ------------------------------ *
   * Destroy workspace for resizing *
   * ------------------------------ */

  for (int i = 0; i < wrk_space_size; i++) { N_VDestroy(resize_wrk[i]); }

  return CV_SUCCESS;
}
