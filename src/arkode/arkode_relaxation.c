/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * This is the implementation file for ARKODE's relaxation (in time)
 * functionality
 *
 * Temporary vectors utilized in the functions below:
 *   tempv2 - holds delta_y, the update direction vector
 *   tempv3 - holds y_relax, the relaxed solution vector
 *   tempv4 - holds J_relax, the Jacobian of the relaxation function
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode_impl.h"
#include "arkode_relaxation_impl.h"
#include "sundials_utils.h"

/* =============================================================================
 * Private Functions
 * ===========================================================================*/

/* Access the ARKODE and relaxation memory structures */
static int arkRelaxAccessMem(void* arkode_mem, const char* fname,
                             ARKodeMem* ark_mem, ARKodeRelaxMem* relax_mem)
{
  if (!arkode_mem)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_ARK_NO_MEM);
    return ARK_MEM_NULL;
  }
  *ark_mem = (ARKodeMem)arkode_mem;

  if (!((*ark_mem)->relax_mem))
  {
    arkProcessError(*ark_mem, ARK_RELAX_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_RELAX_MEM_NULL);
    return ARK_RELAX_MEM_NULL;
  }
  *relax_mem = (ARKodeRelaxMem)((*ark_mem)->relax_mem);

  return ARK_SUCCESS;
}

/* Evaluates the relaxation residual function */
static int arkRelaxResidual(sunrealtype relax_param, sunrealtype* relax_res,
                            ARKodeMem ark_mem)
{
  int retval;
  sunrealtype e_old   = ark_mem->relax_mem->e_old;
  sunrealtype delta_e = ark_mem->relax_mem->delta_e;
  N_Vector delta_y    = ark_mem->tempv2;
  N_Vector y_relax    = ark_mem->tempv3;
  void* user_data     = ark_mem->user_data;

  /* y_relax = y_n + r * delta_y */
  N_VLinearSum(ONE, ark_mem->yn, relax_param, delta_y, y_relax);

  /* Evaluate entropy function */
  retval = ark_mem->relax_mem->relax_fn(y_relax, relax_res, user_data);
  ark_mem->relax_mem->num_relax_fn_evals++;
  if (retval < 0) { return ARK_RELAX_FUNC_FAIL; }
  if (retval > 0) { return ARK_RELAX_FUNC_RECV; }

  /* Compute relaxation residual */
  *relax_res = *relax_res - e_old - relax_param * delta_e;

  return ARK_SUCCESS;
}

/* Evaluates the Jacobian of the relaxation residual function */
static int arkRelaxResidualJacobian(sunrealtype relax_param,
                                    sunrealtype* relax_jac, ARKodeMem ark_mem)
{
  int retval;
  N_Vector delta_y    = ark_mem->tempv2;
  N_Vector y_relax    = ark_mem->tempv3;
  N_Vector J_relax    = ark_mem->tempv4;
  sunrealtype delta_e = ark_mem->relax_mem->delta_e;
  void* user_data     = ark_mem->user_data;

  /* y_relax = y_n + r * delta_y */
  N_VLinearSum(ONE, ark_mem->yn, relax_param, delta_y, y_relax);

  /* Evaluate Jacobian of entropy functions */
  retval = ark_mem->relax_mem->relax_jac_fn(y_relax, J_relax, user_data);
  ark_mem->relax_mem->num_relax_jac_evals++;
  if (retval < 0) { return ARK_RELAX_JAC_FAIL; }
  if (retval > 0) { return ARK_RELAX_JAC_RECV; }

  /* Compute relaxation residual Jacobian */
  *relax_jac = N_VDotProd(delta_y, J_relax);
  *relax_jac -= delta_e;

  return ARK_SUCCESS;
}

/* Solve the relaxation residual equation using Newton's method */
static int arkRelaxNewtonSolve(ARKodeMem ark_mem)
{
  int i, retval;
  sunrealtype tol, delta;
  ARKodeRelaxMem relax_mem = ark_mem->relax_mem;

  for (i = 0; i < ark_mem->relax_mem->max_iters; i++)
  {
    /* Compute the current residual */
    retval = arkRelaxResidual(relax_mem->relax_param, &(relax_mem->res), ark_mem);
    if (retval) { return retval; }

    SUNLogExtraDebug(ARK_LOGGER, "residual",
                     "iter = %i, relax_param = " SUN_FORMAT_G
                     ", residual = " SUN_FORMAT_G,
                     i, relax_mem->relax_param, relax_mem->res);

    /* Check for convergence */
    if (SUNRabs(relax_mem->res) < relax_mem->res_tol) { return ARK_SUCCESS; }

    /* Compute Jacobian */
    retval = arkRelaxResidualJacobian(relax_mem->relax_param, &(relax_mem->jac),
                                      ark_mem);
    if (retval) { return retval; }

    /* Update step length tolerance and solution */
    tol = (relax_mem->rel_tol * SUNRabs(relax_mem->relax_param) +
           relax_mem->abs_tol);

    delta = relax_mem->res / relax_mem->jac;
    relax_mem->relax_param -= delta;

    /* Update cumulative iteration count */
    relax_mem->nls_iters++;

    /* Check for small update */
    if (SUNRabs(delta) < tol) { return ARK_SUCCESS; }
  }

  return ARK_RELAX_SOLVE_RECV;
}

/* Solve the relaxation residual equation using Brent's method */
static int arkRelaxBrentSolve(ARKodeMem ark_mem)
{
  int i, retval;
  sunrealtype xa, fa;         /* previous solution and function value */
  sunrealtype xb, fb;         /* current solution and function value  */
  sunrealtype xc, fc;         /* together brac and curr bracket zero  */
  sunrealtype xm;             /* midpoint between brac and curr       */
  sunrealtype old_update;     /* previous iteration update            */
  sunrealtype new_update;     /* new iteration update                 */
  sunrealtype tol;            /* iteration tolerance                  */
  sunrealtype pt, qt, rt, st; /* temporary values                     */

  ARKodeRelaxMem relax_mem = ark_mem->relax_mem;

  /* Compute interval that brackets the root */
  xa = SUN_RCONST(0.9) * relax_mem->relax_param;
  xb = SUN_RCONST(1.1) * relax_mem->relax_param;

  for (i = 0; i < 10; i++)
  {
    /* Compute relaxation residual */
    retval = arkRelaxResidual(xa, &fa, ark_mem);
    ark_mem->relax_mem->num_relax_fn_evals++;
    if (retval < 0) { return ARK_RELAX_FUNC_FAIL; }
    if (retval > 0) { return ARK_RELAX_FUNC_RECV; }

    /* Check if we got lucky */
    if (SUNRabs(fa) < relax_mem->res_tol)
    {
      relax_mem->res         = fa;
      relax_mem->relax_param = xa;
      return ARK_SUCCESS;
    }

    if (fa < ZERO) { break; }

    fb = fa;
    xb = xa;
    xa *= SUN_RCONST(0.9);
  }
  if (fa > ZERO) { return ARK_RELAX_SOLVE_RECV; }

  for (i = 0; i < 10; i++)
  {
    /* Compute relaxation residual */
    retval = arkRelaxResidual(xb, &fb, ark_mem);
    ark_mem->relax_mem->num_relax_fn_evals++;
    if (retval < 0) { return ARK_RELAX_FUNC_FAIL; }
    if (retval > 0) { return ARK_RELAX_FUNC_RECV; }

    /* Check if we got lucky */
    if (SUNRabs(fb) < relax_mem->res_tol)
    {
      relax_mem->res         = fb;
      relax_mem->relax_param = xb;
      return ARK_SUCCESS;
    }

    if (fb > ZERO) { break; }

    fa = fb;
    xa = xb;
    xb *= SUN_RCONST(1.1);
  }
  if (fb < ZERO) { return ARK_RELAX_SOLVE_RECV; }

  /* Initialize values bracketing values to lower bound and updates */
  xc = xa;
  fc = fa;

  old_update = ZERO;
  new_update = ZERO;

  /* Find root */
  for (i = 0; i < ark_mem->relax_mem->max_iters; i++)
  {
    /* Ensure xc and xb bracket zero */
    if (SUNRsamesign(fc, fb))
    {
      xc         = xa;
      fc         = fa;
      old_update = new_update = xb - xa;
    }

    /* Ensure xb is closer to zero than xc */
    if (SUNRabs(fb) > SUNRabs(fc))
    {
      xa = xb;
      xb = xc;
      xc = xa;

      fa = fb;
      fb = fc;
      fc = fa;
    }

    /* Update tolerance */
    tol = relax_mem->rel_tol * SUNRabs(xb) + HALF * relax_mem->abs_tol;

    /* Compute midpoint for bisection */
    xm = HALF * (xc - xb);

    /* Check for convergence */
    if (SUNRabs(xm) < tol || SUNRabs(fb) < relax_mem->res_tol)
    {
      relax_mem->res         = fb;
      relax_mem->relax_param = xb;
      return ARK_SUCCESS;
    }

    /* Compute iteration update */
    if (SUNRabs(old_update) >= tol && SUNRabs(fb) < SUNRabs(fa))
    {
      /* Converging sufficiently fast, interpolate solution */
      st = fb / fa;

      if (xa == xc)
      {
        /* Two unique values available, try linear interpolant (secant) */
        pt = TWO * xm * st;
        qt = ONE - st;
      }
      else
      {
        /* Three unique values available, try inverse quadratic interpolant */
        qt = fa / fc;
        rt = fb / fc;
        pt = st * (TWO * xm * qt * (qt - rt) - (xb - xa) * (rt - ONE));
        qt = (qt - ONE) * (rt - ONE) * (st - ONE);
      }

      /* Ensure updates produce values within [xc, xb] or [xb, xc] */
      if (pt > ZERO) { qt = -qt; }
      else { pt = -pt; }

      /* Check if interpolant is acceptable, otherwise use bisection */
      st = THREE * xm * qt - SUNRabs(tol * qt);
      rt = SUNRabs(old_update * qt);

      if (TWO * pt < SUNMIN(st, rt))
      {
        old_update = new_update;
        new_update = pt / qt;
      }
      else
      {
        new_update = xm;
        old_update = xm;
      }
    }
    else
    {
      /* Converging too slowly, use bisection */
      new_update = xm;
      old_update = xm;
    }

    /* Update solution */
    xa = xb;
    fa = fb;

    /* If update is small, use tolerance in bisection direction */
    if (SUNRabs(new_update) > tol) { xb += new_update; }
    else { xb += SUNRcopysign(tol, xm); }

    /* Compute relaxation residual */
    retval = arkRelaxResidual(xb, &fb, ark_mem);
    ark_mem->relax_mem->num_relax_fn_evals++;
    if (retval < 0) { return ARK_RELAX_FUNC_FAIL; }
    if (retval > 0) { return ARK_RELAX_FUNC_RECV; }
  }

  return ARK_RELAX_SOLVE_RECV;
}

/* Compute and apply relaxation parameter */
static int arkRelaxSolve(ARKodeMem ark_mem, ARKodeRelaxMem relax_mem,
                         sunrealtype* relax_val_out)
{
  int retval;

  /* Get the change in entropy (uses temp vectors 2 and 3) */
  retval = relax_mem->delta_e_fn(ark_mem, relax_mem->relax_jac_fn,
                                 &(relax_mem->num_relax_jac_evals),
                                 &(relax_mem->delta_e));
  if (retval) { return retval; }

  SUNLogExtraDebug(ARK_LOGGER, "compute delta e", "delta_e = " SUN_FORMAT_G,
                   relax_mem->delta_e);

  /* Get the change in state (delta_y = tempv2) */
  N_VLinearSum(ONE, ark_mem->ycur, -ONE, ark_mem->yn, ark_mem->tempv2);

  SUNLogExtraDebugVec(ARK_LOGGER, "compute delta y", ark_mem->tempv2,
                      "delta_y(:) =");

  /* Store the current relaxation function value */
  retval = relax_mem->relax_fn(ark_mem->yn, &(relax_mem->e_old),
                               ark_mem->user_data);
  relax_mem->num_relax_fn_evals++;
  if (retval < 0) { return ARK_RELAX_FUNC_FAIL; }
  if (retval > 0) { return ARK_RELAX_FUNC_RECV; }

  SUNLogExtraDebug(ARK_LOGGER, "compute old e", "e_old = " SUN_FORMAT_G,
                   relax_mem->e_old);

  /* Initial guess for relaxation parameter */
  relax_mem->relax_param = relax_mem->relax_param_prev;

  switch (relax_mem->solver)
  {
  case (ARK_RELAX_BRENT): retval = arkRelaxBrentSolve(ark_mem); break;
  case (ARK_RELAX_NEWTON): retval = arkRelaxNewtonSolve(ark_mem); break;
  default: return ARK_ILL_INPUT; break;
  }

  /* Check for solver failure */
  if (retval)
  {
    relax_mem->nls_fails++;
    return retval;
  }

  /* Check for bad relaxation value */
  if (ark_mem->relax_mem->relax_param < relax_mem->lower_bound ||
      ark_mem->relax_mem->relax_param > relax_mem->upper_bound)
  {
    relax_mem->bound_fails++;
    return ARK_RELAX_SOLVE_RECV;
  }

  /* Save parameter for next initial guess */
  relax_mem->relax_param_prev = relax_mem->relax_param;

  /* Return relaxation value */
  *relax_val_out = ark_mem->relax_mem->relax_param;

  return ARK_SUCCESS;
}

/* =============================================================================
 * User Functions
 * ===========================================================================*/

/* -----------------------------------------------------------------------------
 * Set functions
 * ---------------------------------------------------------------------------*/

int ARKodeSetRelaxFn(void* arkode_mem, ARKRelaxFn rfn, ARKRelaxJacFn rjac)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Ensure that the current N_Vector supports N_VDotProd */
  if (ark_mem->tempv1->ops->nvdotprod == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "N_VDotProd unimplemented (required for relaxation)");
    return (ARK_ILL_INPUT);
  }

  /* Call stepper-specific routine (if it exists) */
  if (ark_mem->step_setrelaxfn)
  {
    return ark_mem->step_setrelaxfn(arkode_mem, rfn, rjac);
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support relaxation");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

int ARKodeSetRelaxEtaFail(void* arkode_mem, sunrealtype eta_fail)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, __func__, &ark_mem, &relax_mem);
  if (retval) { return retval; }

  /* Guard against use for time steppers that do not allow relaxation */
  if (!ark_mem->step_supports_relaxation)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support relaxation");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  if (eta_fail > ZERO && eta_fail < ONE) { relax_mem->eta_fail = eta_fail; }
  else { relax_mem->eta_fail = ARK_RELAX_DEFAULT_ETA_FAIL; }

  return ARK_SUCCESS;
}

int ARKodeSetRelaxLowerBound(void* arkode_mem, sunrealtype lower)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, __func__, &ark_mem, &relax_mem);
  if (retval) { return retval; }

  /* Guard against use for time steppers that do not allow relaxation */
  if (!ark_mem->step_supports_relaxation)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support relaxation");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  if (lower > ZERO && lower < ONE) { relax_mem->lower_bound = lower; }
  else { relax_mem->lower_bound = ARK_RELAX_DEFAULT_LOWER_BOUND; }

  return ARK_SUCCESS;
}

int ARKodeSetRelaxMaxFails(void* arkode_mem, int max_fails)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, __func__, &ark_mem, &relax_mem);
  if (retval) { return retval; }

  /* Guard against use for time steppers that do not allow relaxation */
  if (!ark_mem->step_supports_relaxation)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support relaxation");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  if (max_fails > 0) { relax_mem->max_fails = max_fails; }
  else { relax_mem->max_fails = ARK_RELAX_DEFAULT_MAX_FAILS; }

  return ARK_SUCCESS;
}

int ARKodeSetRelaxMaxIters(void* arkode_mem, int max_iters)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, __func__, &ark_mem, &relax_mem);
  if (retval) { return retval; }

  /* Guard against use for time steppers that do not allow relaxation */
  if (!ark_mem->step_supports_relaxation)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support relaxation");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  if (max_iters > 0) { relax_mem->max_iters = max_iters; }
  else { relax_mem->max_iters = ARK_RELAX_DEFAULT_MAX_ITERS; }

  return ARK_SUCCESS;
}

int ARKodeSetRelaxSolver(void* arkode_mem, ARKRelaxSolver solver)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, __func__, &ark_mem, &relax_mem);
  if (retval) { return retval; }

  /* Guard against use for time steppers that do not allow relaxation */
  if (!ark_mem->step_supports_relaxation)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support relaxation");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  if (solver != ARK_RELAX_BRENT && solver != ARK_RELAX_NEWTON)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "An invalid relaxation solver option was provided.");
    return ARK_ILL_INPUT;
  }

  relax_mem->solver = solver;

  return ARK_SUCCESS;
}

int ARKodeSetRelaxResTol(void* arkode_mem, sunrealtype res_tol)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, __func__, &ark_mem, &relax_mem);
  if (retval) { return retval; }

  /* Guard against use for time steppers that do not allow relaxation */
  if (!ark_mem->step_supports_relaxation)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support relaxation");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  if (res_tol > ZERO) { relax_mem->res_tol = res_tol; }
  else { relax_mem->res_tol = ARK_RELAX_DEFAULT_RES_TOL; }

  return ARK_SUCCESS;
}

int ARKodeSetRelaxTol(void* arkode_mem, sunrealtype rel_tol, sunrealtype abs_tol)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, __func__, &ark_mem, &relax_mem);
  if (retval) { return retval; }

  /* Guard against use for time steppers that do not allow relaxation */
  if (!ark_mem->step_supports_relaxation)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support relaxation");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  if (rel_tol > ZERO) { relax_mem->rel_tol = rel_tol; }
  else { relax_mem->rel_tol = ARK_RELAX_DEFAULT_REL_TOL; }

  if (abs_tol > ZERO) { relax_mem->abs_tol = abs_tol; }
  else { relax_mem->abs_tol = ARK_RELAX_DEFAULT_ABS_TOL; }

  return ARK_SUCCESS;
}

int ARKodeSetRelaxUpperBound(void* arkode_mem, sunrealtype upper)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, __func__, &ark_mem, &relax_mem);
  if (retval) { return retval; }

  /* Guard against use for time steppers that do not allow relaxation */
  if (!ark_mem->step_supports_relaxation)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support relaxation");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  if (upper > ONE) { relax_mem->upper_bound = upper; }
  else { relax_mem->upper_bound = ARK_RELAX_DEFAULT_UPPER_BOUND; }

  return ARK_SUCCESS;
}

/* -----------------------------------------------------------------------------
 * Get functions
 * ---------------------------------------------------------------------------*/

int ARKodeGetNumRelaxFnEvals(void* arkode_mem, long int* r_evals)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, __func__, &ark_mem, &relax_mem);
  if (retval) { return retval; }

  /* Guard against use for time steppers that do not allow relaxation */
  if (!ark_mem->step_supports_relaxation)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support relaxation");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  *r_evals = relax_mem->num_relax_fn_evals;

  return ARK_SUCCESS;
}

int ARKodeGetNumRelaxJacEvals(void* arkode_mem, long int* J_evals)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, __func__, &ark_mem, &relax_mem);
  if (retval) { return retval; }

  /* Guard against use for time steppers that do not allow relaxation */
  if (!ark_mem->step_supports_relaxation)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support relaxation");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  *J_evals = relax_mem->num_relax_jac_evals;

  return ARK_SUCCESS;
}

int ARKodeGetNumRelaxFails(void* arkode_mem, long int* relax_fails)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, __func__, &ark_mem, &relax_mem);
  if (retval) { return retval; }

  /* Guard against use for time steppers that do not allow relaxation */
  if (!ark_mem->step_supports_relaxation)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support relaxation");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  *relax_fails = relax_mem->num_fails;

  return ARK_SUCCESS;
}

int ARKodeGetNumRelaxSolveFails(void* arkode_mem, long int* fails)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, __func__, &ark_mem, &relax_mem);
  if (retval) { return retval; }

  /* Guard against use for time steppers that do not allow relaxation */
  if (!ark_mem->step_supports_relaxation)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support relaxation");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  *fails = relax_mem->nls_fails;

  return ARK_SUCCESS;
}

int ARKodeGetNumRelaxBoundFails(void* arkode_mem, long int* fails)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, __func__, &ark_mem, &relax_mem);
  if (retval) { return retval; }

  /* Guard against use for time steppers that do not allow relaxation */
  if (!ark_mem->step_supports_relaxation)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support relaxation");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  *fails = relax_mem->bound_fails;

  return ARK_SUCCESS;
}

int ARKodeGetNumRelaxSolveIters(void* arkode_mem, long int* iters)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, __func__, &ark_mem, &relax_mem);
  if (retval) { return retval; }

  /* Guard against use for time steppers that do not allow relaxation */
  if (!ark_mem->step_supports_relaxation)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support relaxation");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  *iters = relax_mem->nls_iters;

  return ARK_SUCCESS;
}

/* =============================================================================
 * Driver and Stepper Functions
 * ===========================================================================*/

/* Constructor called by stepper */
int arkRelaxCreate(ARKodeMem ark_mem, ARKRelaxFn relax_fn,
                   ARKRelaxJacFn relax_jac_fn, ARKRelaxDeltaEFn delta_e_fn,
                   ARKRelaxGetOrderFn get_order_fn)
{
  /* Disable relaxation if both user inputs are NULL */
  if (!relax_fn && !relax_jac_fn)
  {
    ark_mem->relax_enabled = SUNFALSE;
    return ARK_SUCCESS;
  }

  /* Ensure both the relaxation function and Jacobian are provided */
  if (!relax_fn)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The relaxation function is NULL.");
    return ARK_ILL_INPUT;
  }

  if (!relax_jac_fn)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The relaxation Jacobian function is NULL.");
    return ARK_ILL_INPUT;
  }

  /* Ensure stepper supplied inputs are provided */
  if (!delta_e_fn || !get_order_fn)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "The Delta y, Delta e, or get order function is NULL.");
    return ARK_ILL_INPUT;
  }

  /* Allocate and initialize relaxation memory structure */
  if (!(ark_mem->relax_mem))
  {
    ark_mem->relax_mem = (ARKodeRelaxMem)malloc(sizeof(*(ark_mem->relax_mem)));
    if (!(ark_mem->relax_mem)) { return ARK_MEM_FAIL; }
    memset(ark_mem->relax_mem, 0, sizeof(struct ARKodeRelaxMemRec));

    /* Set defaults */
    ark_mem->relax_mem->max_fails   = ARK_RELAX_DEFAULT_MAX_FAILS;
    ark_mem->relax_mem->lower_bound = ARK_RELAX_DEFAULT_LOWER_BOUND;
    ark_mem->relax_mem->upper_bound = ARK_RELAX_DEFAULT_UPPER_BOUND;
    ark_mem->relax_mem->eta_fail    = ARK_RELAX_DEFAULT_ETA_FAIL;
    ark_mem->relax_mem->solver      = ARK_RELAX_NEWTON;
    ark_mem->relax_mem->res_tol     = ARK_RELAX_DEFAULT_RES_TOL;
    ark_mem->relax_mem->rel_tol     = ARK_RELAX_DEFAULT_REL_TOL;
    ark_mem->relax_mem->abs_tol     = ARK_RELAX_DEFAULT_ABS_TOL;
    ark_mem->relax_mem->max_iters   = ARK_RELAX_DEFAULT_MAX_ITERS;

    /* Initialize values */
    ark_mem->relax_mem->relax_param_prev = ONE;

    /* Update workspace sizes */
    ark_mem->lrw += 12;
    ark_mem->liw += 14;
  }

  /* Set function pointers */
  ark_mem->relax_mem->relax_fn     = relax_fn;
  ark_mem->relax_mem->relax_jac_fn = relax_jac_fn;
  ark_mem->relax_mem->delta_e_fn   = delta_e_fn;
  ark_mem->relax_mem->get_order_fn = get_order_fn;

  /* Enable relaxation */
  ark_mem->relax_enabled = SUNTRUE;

  return ARK_SUCCESS;
}

/* Destructor called by driver */
int arkRelaxDestroy(ARKodeRelaxMem relax_mem)
{
  if (!relax_mem) { return ARK_SUCCESS; }

  /* Free structure */
  free(relax_mem);

  return ARK_SUCCESS;
}

/* Compute and apply relaxation, called by driver */
int arkRelax(ARKodeMem ark_mem, int* relax_fails, sunrealtype* dsm_inout)
{
  int retval;
  sunrealtype relax_val;
  ARKodeRelaxMem relax_mem = ark_mem->relax_mem;

  /* Get the relaxation memory structure */
  if (!relax_mem)
  {
    arkProcessError(ark_mem, ARK_RELAX_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_RELAX_MEM_NULL);
    return ARK_RELAX_MEM_NULL;
  }

  /* Compute the relaxation parameter */
  retval = arkRelaxSolve(ark_mem, relax_mem, &relax_val);
  if (retval < 0) { return retval; }
  if (retval > 0)
  {
    /* Update failure counts */
    relax_mem->num_fails++;
    (*relax_fails)++;

    /* Check for max fails in a step */
    if (*relax_fails == relax_mem->max_fails) { return ARK_RELAX_FAIL; }

    /* Return with an error if |h| == hmin */
    if (SUNRabs(ark_mem->h) <= ark_mem->hmin * ONEPSM)
    {
      return ARK_RELAX_FAIL;
    }

    /* Return with error if using fixed step sizes */
    if (ark_mem->fixedstep) { return (ARK_RELAX_FAIL); }

    /* Cut step size and try again */
    ark_mem->eta = relax_mem->eta_fail;

    return TRY_AGAIN;
  }

  /* Update step size and error estimate */
  ark_mem->h *= relax_val;
  *dsm_inout *= SUNRpowerI(relax_val, relax_mem->get_order_fn(ark_mem));

  /* Relax solution */
  N_VLinearSum(relax_val, ark_mem->ycur, (ONE - relax_val), ark_mem->yn,
               ark_mem->ycur);

  SUNLogDebug(ARK_LOGGER, "relaxation",
              "relaxation parameter = " SUN_FORMAT_G
              ", relaxed h = " SUN_FORMAT_G ", relaxed error = " SUN_FORMAT_G,
              relax_val, ark_mem->h, *dsm_inout);

  return ARK_SUCCESS;
}

/* Print relaxation solver statistics, called by ARKODE */
int arkRelaxPrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRelaxMem relax_mem;

  retval = arkRelaxAccessMem(arkode_mem, __func__, &ark_mem, &relax_mem);
  if (retval) { return retval; }

  sunfprintf_long(outfile, fmt, SUNFALSE, "Relax fn evals",
                  relax_mem->num_relax_fn_evals);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Relax Jac evals",
                  relax_mem->num_relax_jac_evals);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Relax fails", relax_mem->num_fails);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Relax bound fails",
                  relax_mem->bound_fails);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Relax NLS iters",
                  relax_mem->nls_iters);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Relax NLS fails",
                  relax_mem->nls_fails);

  return ARK_SUCCESS;
}

/* =============================================================================
 * EOF
 * ===========================================================================*/
