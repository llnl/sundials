/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *---------------------------------------------------------------
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
 *---------------------------------------------------------------
 * This is the implementation file for the optional input and
 * output functions for the ARKODE infrastructure.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sunadaptcontroller/sunadaptcontroller_imexgus.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "arkode/arkode.h"
#include "arkode_impl.h"
#include "arkode_interp_impl.h"
#include "arkode_user_controller.h"
#include "sundials_utils.h"

/*===============================================================
  ARKODE optional input functions
  ===============================================================*/

/*---------------------------------------------------------------
  ARKodeSetDefaults:

  Resets all optional inputs to ARKODE default values.  Does not
  change problem-defining function pointers fe and fi or
  user_data pointer.  Also leaves alone any data
  structures/options related to root-finding (those can be reset
  using ARKodeRootInit) or post-processing a step (ProcessStep).
  ---------------------------------------------------------------*/
int ARKodeSetDefaults(void* arkode_mem)
{
  ARKodeMem ark_mem;
  int retval;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Set default values for integrator optional inputs */
  ark_mem->use_compensated_sums = SUNFALSE;
  ark_mem->fixedstep            = SUNFALSE; /* default to use adaptive steps */
  ark_mem->reltol               = SUN_RCONST(1.e-4); /* relative tolerance */
  ark_mem->itol      = ARK_SS; /* scalar-scalar solution tolerances */
  ark_mem->ritol     = ARK_SS; /* scalar-scalar residual tolerances */
  ark_mem->Sabstol   = SUN_RCONST(1.e-9); /* solution absolute tolerance */
  ark_mem->atolmin0  = SUNFALSE;          /* min(abstol) > 0 */
  ark_mem->SRabstol  = SUN_RCONST(1.e-9); /* residual absolute tolerance */
  ark_mem->Ratolmin0 = SUNFALSE;          /* min(Rabstol) > 0 */
  ark_mem->user_efun = SUNFALSE;          /* no user-supplied ewt function */
  ark_mem->efun      = arkEwtSetSS;    /* built-in scalar-scalar ewt function */
  ark_mem->e_data    = ark_mem;        /* ewt function data */
  ark_mem->user_rfun = SUNFALSE;       /* no user-supplied rwt function */
  ark_mem->rfun      = arkRwtSet;      /* built-in rwt function */
  ark_mem->r_data    = ark_mem;        /* rwt function data */
  ark_mem->mxstep    = MXSTEP_DEFAULT; /* max number of steps */
  ark_mem->mxhnil    = MXHNIL;         /* max warns of t+h==t */
  ark_mem->maxnef    = MAXNEF;         /* max error test fails */
  ark_mem->maxncf    = MAXNCF;         /* max convergence fails */
  ark_mem->maxconstrfails = MAXCONSTRFAILS; /* max number of constraint fails */
  ark_mem->hin            = ZERO;       /* determine initial step on-the-fly */
  ark_mem->hmin           = ZERO;       /* no minimum step size */
  ark_mem->hmax_inv       = ZERO;       /* no maximum step size */
  ark_mem->tstopset       = SUNFALSE;   /* no stop time set */
  ark_mem->tstopinterp    = SUNFALSE;   /* copy at stop time */
  ark_mem->tstop          = ZERO;       /* no fixed stop time */
  ark_mem->hadapt_mem->etamx1 = ETAMX1; /* max change on first step */
  ark_mem->hadapt_mem->etamxf = ETAMXF; /* max change on error-failed step */
  ark_mem->hadapt_mem->etamin = ETAMIN; /* min bound on time step reduction */
  ark_mem->hadapt_mem->small_nef =
    SMALL_NEF; /* num error fails before ETAMXF enforced */
  ark_mem->hadapt_mem->etacf  = ETACF;  /* max change on convergence failure */
  ark_mem->hadapt_mem->cfl    = CFLFAC; /* explicit stability factor */
  ark_mem->hadapt_mem->safety = SAFETY; /* step adaptivity safety factor  */
  ark_mem->hadapt_mem->growth = GROWTH; /* step adaptivity growth factor */
  ark_mem->hadapt_mem->lbound = HFIXED_LB; /* step adaptivity no-change lower bound */
  ark_mem->hadapt_mem->ubound = HFIXED_UB; /* step adaptivity no-change upper bound */
  ark_mem->hadapt_mem->expstab    = NULL;   /* no explicit stability fn */
  ark_mem->hadapt_mem->estab_data = NULL;   /* no explicit stability fn data */
  ark_mem->hadapt_mem->pq         = PQ;     /* embedding order */
  ark_mem->hadapt_mem->p          = 0;      /* no default embedding order */
  ark_mem->hadapt_mem->q          = 0;      /* no default method order */
  ark_mem->hadapt_mem->adjust     = ADJUST; /* controller order adjustment */

  /* Set stepper defaults (if provided) */
  if (ark_mem->step_setdefaults)
  {
    retval = ark_mem->step_setdefaults(arkode_mem);
    if (retval != ARK_SUCCESS) { return retval; }
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetOrder:

  Specifies the method order
  ---------------------------------------------------------------*/
int ARKodeSetOrder(void* arkode_mem, int ord)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Call stepper routine (if provided) */
  if (ark_mem->step_setorder)
  {
    return (ark_mem->step_setorder(arkode_mem, ord));
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

/*---------------------------------------------------------------
  ARKodeSetInterpolantType:

  Specifies use of the Lagrange or Hermite interpolation modules.
    itype == ARK_INTERP_HERMITE specifies the Hermite (nonstiff)
      interpolation module.
    itype == ARK_INTERP_LAGRANGE specifies the Lagrange (stiff)
      interpolation module.
    itype == ARK_INTERP_NONE disables interpolation.

  Return values:
     ARK_SUCCESS on success.
     ARK_MEM_NULL on NULL-valued arkode_mem input.
     ARK_MEM_FAIL if the interpolation module cannot be allocated.
     ARK_ILL_INPUT if the itype argument is not recognized.
  ---------------------------------------------------------------*/
int ARKodeSetInterpolantType(void* arkode_mem, int itype)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* check for legal itype input */
  if ((itype != ARK_INTERP_HERMITE) && (itype != ARK_INTERP_LAGRANGE) &&
      (itype != ARK_INTERP_NONE))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Illegal interpolation type input.");
    return (ARK_ILL_INPUT);
  }

  /* do not change type once the module has been initialized */
  if (ark_mem->initialized)
  {
    arkProcessError(ark_mem, ARK_INTERP_FAIL, __LINE__, __func__, __FILE__,
                    "Type cannot be specified after module initialization.");
    return (ARK_ILL_INPUT);
  }

  /* delete any existing interpolation module */
  if (ark_mem->interp != NULL)
  {
    arkInterpFree(ark_mem, ark_mem->interp);
    ark_mem->interp = NULL;
  }

  /* create requested interpolation module, initially specifying
     the maximum possible interpolant degree. */
  if (itype == ARK_INTERP_HERMITE)
  {
    ark_mem->interp = arkInterpCreate_Hermite(arkode_mem, ark_mem->interp_degree);
    if (ark_mem->interp == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to allocate interpolation structure");
      return ARK_MEM_FAIL;
    }
    ark_mem->interp_type = ARK_INTERP_HERMITE;
  }
  else if (itype == ARK_INTERP_LAGRANGE)
  {
    ark_mem->interp = arkInterpCreate_Lagrange(arkode_mem,
                                               ark_mem->interp_degree);
    if (ark_mem->interp == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to allocate interpolation structure");
    }
    ark_mem->interp_type = ARK_INTERP_LAGRANGE;
  }
  else
  {
    ark_mem->interp      = NULL;
    ark_mem->interp_type = ARK_INTERP_NONE;
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetInterpolantDegree:

  Specifies the polynomial degree for the dense output
  interpolation module.

  Return values:
     ARK_SUCCESS on success.
     ARK_MEM_NULL on NULL-valued arkode_mem input or nonexistent
       interpolation module.
     ARK_INTERP_FAIL if the interpolation module is already
       initialized.
     ARK_ILL_INPUT if the degree is illegal.
  ---------------------------------------------------------------*/
int ARKodeSetInterpolantDegree(void* arkode_mem, int degree)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* do not change degree once the module has been initialized */
  if (ark_mem->initialized)
  {
    arkProcessError(ark_mem, ARK_INTERP_FAIL, __LINE__, __func__, __FILE__,
                    "Degree cannot be specified after module initialization.");
    return (ARK_ILL_INPUT);
  }

  if (degree > ARK_INTERP_MAX_DEGREE)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Illegal degree specified.");
    return ARK_ILL_INPUT;
  }
  else if (degree < 0) { ark_mem->interp_degree = ARK_INTERP_MAX_DEGREE; }
  else { ark_mem->interp_degree = degree; }

  /* Set the degree now if possible otherwise it will be used when creating the
     interpolation module */
  if (ark_mem->interp)
  {
    return arkInterpSetDegree(ark_mem, ark_mem->interp, ark_mem->interp_degree);
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  ARKodeSetNonlinearSolver:

  This routine attaches a SUNNonlinearSolver object to the
  time-stepping module.
  ---------------------------------------------------------------*/
int ARKodeSetNonlinearSolver(void* arkode_mem, SUNNonlinearSolver NLS)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not need an algebraic solver */
  if (!ark_mem->step_supports_implicit)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not require an algebraic solver");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Call stepper routine (if provided) */
  if (ark_mem->step_setnonlinearsolver)
  {
    return (ark_mem->step_setnonlinearsolver(arkode_mem, NLS));
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

/*---------------------------------------------------------------
  ARKodeSetLinear:

  Specifies that the implicit portion of the problem is linear,
  and to tighten the linear solver tolerances while taking only
  one Newton iteration.  DO NOT USE IN COMBINATION WITH THE
  FIXED-POINT SOLVER.  Automatically tightens DeltaGammaMax
  to ensure that step size changes cause Jacobian recomputation.

  The argument should be 1 or 0, where 1 indicates that the
  Jacobian of fi with respect to y depends on time, and
  0 indicates that it is not time dependent.  Alternately, when
  using an iterative linear solver this flag denotes time
  dependence of the preconditioner.
  ---------------------------------------------------------------*/
int ARKodeSetLinear(void* arkode_mem, int timedepend)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not need an algebraic solver */
  if (!ark_mem->step_supports_implicit)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not require an algebraic solver");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Call stepper routine (if provided) */
  if (ark_mem->step_setlinear)
  {
    return (ark_mem->step_setlinear(arkode_mem, timedepend));
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

/*---------------------------------------------------------------
  ARKodeSetNonlinear:

  Specifies that the implicit portion of the problem is nonlinear.
  Used to undo a previous call to ARKodeSetLinear.
  ---------------------------------------------------------------*/
int ARKodeSetNonlinear(void* arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not need an algebraic solver */
  if (!ark_mem->step_supports_implicit)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not require an algebraic solver");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Call stepper routine (if provided) */
  if (ark_mem->step_setnonlinear)
  {
    return (ark_mem->step_setnonlinear(arkode_mem));
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

int ARKodeSetAutonomous(void* arkode_mem, sunbooleantype autonomous)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return ARK_MEM_NULL;
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not need an algebraic solver */
  if (!ark_mem->step_supports_implicit)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not require an algebraic solver");
    return ARK_STEPPER_UNSUPPORTED;
  }

  /* Call stepper routine (if provided) */
  if (ark_mem->step_setautonomous)
  {
    return ark_mem->step_setautonomous(arkode_mem, autonomous);
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return ARK_STEPPER_UNSUPPORTED;
  }
}

/*---------------------------------------------------------------
  ARKodeSetNlsRhsFn:

  This routine sets an alternative user-supplied implicit ODE
  right-hand side function to use in the evaluation of nonlinear
  system functions.
  ---------------------------------------------------------------*/
int ARKodeSetNlsRhsFn(void* arkode_mem, ARKRhsFn nls_fi)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not need an algebraic solver */
  if (!ark_mem->step_supports_implicit)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not require an algebraic solver");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Call stepper routine (if provided) */
  if (ark_mem->step_setnlsrhsfn)
  {
    return (ark_mem->step_setnlsrhsfn(arkode_mem, nls_fi));
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

/*---------------------------------------------------------------
  ARKodeSetDeduceImplicitRhs:

  Specifies if an optimization is used to avoid an evaluation of
  fi after a nonlinear solve for an implicit stage.  If stage
  postprocessecing in enabled, this option is ignored, and the
  RHS is never deduced.

  An argument of SUNTRUE indicates that the RHS should be deduced,
  and SUNFALSE indicates that the RHS should be computed with
  an additional evaluation.
  ---------------------------------------------------------------*/
int ARKodeSetDeduceImplicitRhs(void* arkode_mem, sunbooleantype deduce)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not need an algebraic solver */
  if (!ark_mem->step_supports_implicit)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not require an algebraic solver");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Call stepper routine (if provided) */
  if (ark_mem->step_setdeduceimplicitrhs)
  {
    return (ark_mem->step_setdeduceimplicitrhs(arkode_mem, deduce));
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

/*---------------------------------------------------------------
  ARKodeSetNonlinCRDown:

  Specifies the user-provided nonlinear convergence constant
  crdown.  Legal values are strictly positive; illegal values
  imply a reset to the default.
  ---------------------------------------------------------------*/
int ARKodeSetNonlinCRDown(void* arkode_mem, sunrealtype crdown)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not need an algebraic solver */
  if (!ark_mem->step_supports_implicit)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not require an algebraic solver");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Call stepper routine (if provided) */
  if (ark_mem->step_setnonlincrdown)
  {
    return (ark_mem->step_setnonlincrdown(arkode_mem, crdown));
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

/*---------------------------------------------------------------
  ARKodeSetNonlinRDiv:

  Specifies the user-provided nonlinear convergence constant
  rdiv.  Legal values are strictly positive; illegal values
  imply a reset to the default.
  ---------------------------------------------------------------*/
int ARKodeSetNonlinRDiv(void* arkode_mem, sunrealtype rdiv)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not need an algebraic solver */
  if (!ark_mem->step_supports_implicit)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not require an algebraic solver");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Call stepper routine (if provided) */
  if (ark_mem->step_setnonlinrdiv)
  {
    return (ark_mem->step_setnonlinrdiv(arkode_mem, rdiv));
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

/*---------------------------------------------------------------
  ARKodeSetDeltaGammaMax:

  Specifies the user-provided linear setup decision constant
  dgmax.  Legal values are strictly positive; illegal values imply
  a reset to the default.
  ---------------------------------------------------------------*/
int ARKodeSetDeltaGammaMax(void* arkode_mem, sunrealtype dgmax)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not need an algebraic solver */
  if (!ark_mem->step_supports_implicit)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not require an algebraic solver");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Call stepper routine (if provided) */
  if (ark_mem->step_setdeltagammamax)
  {
    return (ark_mem->step_setdeltagammamax(arkode_mem, dgmax));
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

/*---------------------------------------------------------------
  ARKodeSetLSetupFrequency:

  Specifies the user-provided linear setup decision constant
  msbp.  Positive values give the frequency for calling lsetup;
  negative values imply recomputation of lsetup at each nonlinear
  solve; a zero value implies a reset to the default.
  ---------------------------------------------------------------*/
int ARKodeSetLSetupFrequency(void* arkode_mem, int msbp)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not need an algebraic solver */
  if (!ark_mem->step_supports_implicit)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not require an algebraic solver");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Call stepper routine (if provided) */
  if (ark_mem->step_setlsetupfrequency)
  {
    return (ark_mem->step_setlsetupfrequency(arkode_mem, msbp));
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

/*---------------------------------------------------------------
  ARKodeSetPredictorMethod:

  Specifies the method to use for predicting implicit solutions.
  ---------------------------------------------------------------*/
int ARKodeSetPredictorMethod(void* arkode_mem, int pred_method)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not need an algebraic solver */
  if (!ark_mem->step_supports_implicit)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not require an algebraic solver");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Higher-order predictors require interpolation */
  if (ark_mem->interp_type == ARK_INTERP_NONE && pred_method != 0)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Non-trival predictors require an interpolation module");
    return ARK_ILL_INPUT;
  }

  /* Call stepper routine (if provided) */
  if (ark_mem->step_setpredictormethod)
  {
    return (ark_mem->step_setpredictormethod(arkode_mem, pred_method));
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

/*---------------------------------------------------------------
  ARKodeSetMaxNonlinIters:

  Specifies the maximum number of nonlinear iterations during
  one solve.  A non-positive input implies a reset to the
  default value.
  ---------------------------------------------------------------*/
int ARKodeSetMaxNonlinIters(void* arkode_mem, int maxcor)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not need an algebraic solver */
  if (!ark_mem->step_supports_implicit)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not require an algebraic solver");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Call stepper routine (if provided) */
  if (ark_mem->step_setmaxnonliniters)
  {
    return (ark_mem->step_setmaxnonliniters(arkode_mem, maxcor));
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

/*---------------------------------------------------------------
  ARKodeSetNonlinConvCoef:

  Specifies the coefficient in the nonlinear solver convergence
  test.  A non-positive input implies a reset to the default value.
  ---------------------------------------------------------------*/
int ARKodeSetNonlinConvCoef(void* arkode_mem, sunrealtype nlscoef)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not need an algebraic solver */
  if (!ark_mem->step_supports_implicit)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not require an algebraic solver");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Call stepper routine (if provided) */
  if (ark_mem->step_setnonlinconvcoef)
  {
    return (ark_mem->step_setnonlinconvcoef(arkode_mem, nlscoef));
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

/*---------------------------------------------------------------
  ARKodeSetStagePredictFn:  Specifies a user-provided step
  predictor function having type ARKStagePredictFn.  A
  NULL input function disables calls to this routine.
  ---------------------------------------------------------------*/
int ARKodeSetStagePredictFn(void* arkode_mem, ARKStagePredictFn PredictStage)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not need an algebraic solver */
  if (!ark_mem->step_supports_implicit)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not require an algebraic solver");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Call stepper routine (if provided) */
  if (ark_mem->step_setstagepredictfn)
  {
    return (ark_mem->step_setstagepredictfn(arkode_mem, PredictStage));
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

/*---------------------------------------------------------------
  ARKodeSetUserData:

  Specifies the user data pointer for f
  ---------------------------------------------------------------*/
int ARKodeSetUserData(void* arkode_mem, void* user_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem            = (ARKodeMem)arkode_mem;
  ark_mem->user_data = user_data;

  /* Set data for efun */
  if (ark_mem->user_efun) { ark_mem->e_data = user_data; }

  /* Set data for rfun */
  if (ark_mem->user_rfun) { ark_mem->r_data = user_data; }

  /* Set data for root finding */
  if (ark_mem->root_mem != NULL) { ark_mem->root_mem->root_data = user_data; }

  /* Set data for post-processing a step */
  if (ark_mem->ProcessStep != NULL) { ark_mem->ps_data = user_data; }

  /* Set user data into stepper (if provided) */
  if (ark_mem->step_setuserdata)
  {
    return (ark_mem->step_setuserdata(arkode_mem, user_data));
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetAdaptController:

  Specifies a non-default SUNAdaptController time step controller
  object. If a NULL-valued SUNAdaptController is input, the
  default will be re-enabled.
  ---------------------------------------------------------------*/
int ARKodeSetAdaptController(void* arkode_mem, SUNAdaptController C)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* If the stepper has provided a custom function, then call it and return */
  if (ark_mem->step_setadaptcontroller)
  {
    return (ark_mem->step_setadaptcontroller(ark_mem, C));
  }

  /* Otherwise, call a utility routine to replace the current controller object */
  return (arkReplaceAdaptController(ark_mem, C, SUNFALSE));
}

/*---------------------------------------------------------------
  ARKodeSetAdaptControllerByName:

  Specifies a SUNAdaptController time step controller object by
  its name.
  ---------------------------------------------------------------*/
int ARKodeSetAdaptControllerByName(void* arkode_mem, const char* cname)
{
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Create new controller based on the name */
  SUNAdaptController C = NULL;
  if (strcmp(cname, "Soderlind") == 0)
  {
    C = SUNAdaptController_Soderlind(ark_mem->sunctx);
  }
  else if (strcmp(cname, "PID") == 0)
  {
    C = SUNAdaptController_PID(ark_mem->sunctx);
  }
  else if (strcmp(cname, "PI") == 0)
  {
    C = SUNAdaptController_PI(ark_mem->sunctx);
  }
  else if (strcmp(cname, "I") == 0)
  {
    C = SUNAdaptController_I(ark_mem->sunctx);
  }
  else if (strcmp(cname, "ExpGus") == 0)
  {
    C = SUNAdaptController_ExpGus(ark_mem->sunctx);
  }
  else if (strcmp(cname, "ImpGus") == 0)
  {
    C = SUNAdaptController_ImpGus(ark_mem->sunctx);
  }
  else if (strcmp(cname, "ImExGus") == 0)
  {
    C = SUNAdaptController_ImExGus(ark_mem->sunctx);
  }
  else if (strcmp(cname, "H0211") == 0)
  {
    C = SUNAdaptController_H0211(ark_mem->sunctx);
  }
  else if (strcmp(cname, "H0321") == 0)
  {
    C = SUNAdaptController_H0321(ark_mem->sunctx);
  }
  else if (strcmp(cname, "H211") == 0)
  {
    C = SUNAdaptController_H211(ark_mem->sunctx);
  }
  else if (strcmp(cname, "H312") == 0)
  {
    C = SUNAdaptController_H312(ark_mem->sunctx);
  }
  else
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Unknown controller");
    return ARK_ILL_INPUT;
  }

  if (C == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "SUNAdaptController allocation failure");
    return (ARK_MEM_FAIL);
  }

  /* Send controller to be used by ARKODE */
  retval = ARKodeSetAdaptController(arkode_mem, C);
  if (retval != ARK_SUCCESS) { return retval; }

  /* Update controller ownership flag */
  ark_mem->hadapt_mem->owncontroller = SUNTRUE;

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  ARKodeSetMaxNumSteps:

  Specifies the maximum number of integration steps
  ---------------------------------------------------------------*/
int ARKodeSetMaxNumSteps(void* arkode_mem, long int mxsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Passing mxsteps=0 sets the default. Passing mxsteps<0 disables the test. */
  if (mxsteps == 0) { ark_mem->mxstep = MXSTEP_DEFAULT; }
  else { ark_mem->mxstep = mxsteps; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetMaxHnilWarns:

  Specifies the maximum number of warnings for small h
  ---------------------------------------------------------------*/
int ARKodeSetMaxHnilWarns(void* arkode_mem, int mxhnil)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Passing mxhnil=0 sets the default, otherwise use input. */
  if (mxhnil == 0) { ark_mem->mxhnil = 10; }
  else { ark_mem->mxhnil = mxhnil; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetInitStep:

  Specifies the initial step size
  ---------------------------------------------------------------*/
int ARKodeSetInitStep(void* arkode_mem, sunrealtype hin)
{
  ARKodeMem ark_mem;
  int retval;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against hin==0 for non-adaptive time stepper modules */
  if ((!ark_mem->step_supports_adaptive) && (hin == ZERO))
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Passing hin=0 sets the default, otherwise use input. */
  if (hin == ZERO) { ark_mem->hin = ZERO; }
  else { ark_mem->hin = hin; }

  /* Clear previous initial step */
  ark_mem->h0u = ZERO;

  /* Reset error controller (e.g., error and step size history) */
  if (ark_mem->hadapt_mem->hcontroller != NULL)
  {
    retval = SUNAdaptController_Reset(ark_mem->hadapt_mem->hcontroller);
    if (retval != SUN_SUCCESS) { return (ARK_CONTROLLER_ERR); }
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetMinStep:

  Specifies the minimum step size
  ---------------------------------------------------------------*/
int ARKodeSetMinStep(void* arkode_mem, sunrealtype hmin)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Passing a value <= 0 sets hmin = 0 */
  if (hmin <= ZERO)
  {
    ark_mem->hmin = ZERO;
    return (ARK_SUCCESS);
  }

  /* check that hmin and hmax are agreeable */
  if (hmin * ark_mem->hmax_inv > ONE)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_BAD_HMIN_HMAX);
    return (ARK_ILL_INPUT);
  }

  /* set the value */
  ark_mem->hmin = hmin;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetMaxStep:

  Specifies the maximum step size
  ---------------------------------------------------------------*/
int ARKodeSetMaxStep(void* arkode_mem, sunrealtype hmax)
{
  sunrealtype hmax_inv;
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Passing a value <= 0 sets hmax = infinity */
  if (hmax <= ZERO)
  {
    ark_mem->hmax_inv = ZERO;
    return (ARK_SUCCESS);
  }

  /* check that hmax and hmin are agreeable */
  hmax_inv = ONE / hmax;
  if (hmax_inv * ark_mem->hmin > ONE)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_BAD_HMIN_HMAX);
    return (ARK_ILL_INPUT);
  }

  /* set the value */
  ark_mem->hmax_inv = hmax_inv;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetStopTime:

  Specifies the time beyond which the integration is not to proceed.
  ---------------------------------------------------------------*/
int ARKodeSetStopTime(void* arkode_mem, sunrealtype tstop)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* If ARKODE was called at least once, test if tstop is legal
     (i.e. if it was not already passed).
     If ARKodeSetStopTime is called before the first call to ARKODE,
     tstop will be checked in ARKODE. */
  if (ark_mem->nst > 0)
  {
    if ((tstop - ark_mem->tcur) * ark_mem->h < ZERO)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      MSG_ARK_BAD_TSTOP, tstop, ark_mem->tcur);
      return (ARK_ILL_INPUT);
    }
  }

  ark_mem->tstop    = tstop;
  ark_mem->tstopset = SUNTRUE;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetInterpolateStopTime:

  Specifies to use interpolation to fill the solution output at
  the stop time (instead of a copy).
  ---------------------------------------------------------------*/
int ARKodeSetInterpolateStopTime(void* arkode_mem, sunbooleantype interp)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem              = (ARKodeMem)arkode_mem;
  ark_mem->tstopinterp = interp;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeClearStopTime:

  Disable the stop time.
  ---------------------------------------------------------------*/
int ARKodeClearStopTime(void* arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  ark_mem->tstopset = SUNFALSE;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetFixedStep:

  Specifies to use a fixed time step size instead of performing
  any form of temporal adaptivity.  ARKODE will use this step size
  for all steps (unless tstop is set, in which case it may need to
  modify that last step approaching tstop.  If any solver failure
  occurs in the timestepping module, ARKODE will typically
  immediately return with an error message indicating that the
  selected step size cannot be used.

  Any nonzero argument will result in the use of that fixed step
  size; an argument of 0 will re-enable temporal adaptivity.
  ---------------------------------------------------------------*/
int ARKodeSetFixedStep(void* arkode_mem, sunrealtype hfixed)
{
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* ensure that when hfixed=0, the time step module supports adaptivity */
  if ((hfixed == ZERO) && (!ark_mem->step_supports_adaptive))
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "temporal adaptivity is not supported by this time step module");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* re-attach internal error weight functions if necessary */
  if ((hfixed == ZERO) && (!ark_mem->user_efun))
  {
    if (ark_mem->itol == ARK_SV && ark_mem->Vabstol != NULL)
    {
      retval = ARKodeSVtolerances(ark_mem, ark_mem->reltol, ark_mem->Vabstol);
    }
    else
    {
      retval = ARKodeSStolerances(ark_mem, ark_mem->reltol, ark_mem->Sabstol);
    }
    if (retval != ARK_SUCCESS) { return (retval); }
  }

  /* set ark_mem "fixedstep" entry */
  if (hfixed != ZERO)
  {
    ark_mem->fixedstep = SUNTRUE;
    ark_mem->hin       = hfixed;
  }
  else { ark_mem->fixedstep = SUNFALSE; }

  /* Notify ARKODE to use hfixed as the initial step size, and return */
  return ARKodeSetInitStep(arkode_mem, hfixed);
}

/*---------------------------------------------------------------
  ARKodeSetStepDirection:

  Specifies the direction of integration (forward or backward)
  based on the sign of stepdir. If 0, the direction will remain
  unchanged. Note that if a fixed step size was previously set,
  this function can change the sign of that.

  This should only be called after ARKodeReset, or between
  creating a stepper and ARKodeEvolve.
  ---------------------------------------------------------------*/
int ARKodeSetStepDirection(void* arkode_mem, sunrealtype stepdir)
{
  /* stepdir is a sunrealtype because the direction typically comes from a time
   * step h or tend-tstart which are sunrealtypes. If stepdir was in int,
   * conversions would be required which can cause undefined behavior when
   * greater than MAX_INT */
  int retval;
  ARKodeMem ark_mem;
  sunrealtype h;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return ARK_MEM_NULL;
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* do not change direction once the module has been initialized i.e., after calling
     ARKodeEvolve unless ReInit or Reset are called. */
  if (!ark_mem->initsetup)
  {
    arkProcessError(ark_mem, ARK_STEP_DIRECTION_ERR, __LINE__, __func__,
                    __FILE__, "Step direction cannot be specified after module initialization.");
    return ARK_STEP_DIRECTION_ERR;
  }

  if (stepdir != ZERO)
  {
    retval = ARKodeGetStepDirection(arkode_mem, &h);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                      "Unable to access step direction");
      return retval;
    }

    if (h != SUNRcopysign(h, stepdir))
    {
      /* Reverse the sign of h. If adaptive, h will be overwritten anyway by the
       * initial step estimation since ARKodeReset must be called before this.
       * However, the sign of h will be used to check if the integration
       * direction and stop time are consistent, e.g., in ARKodeSetStopTime, so
       * we should not set h = 0. */
      ark_mem->h = -h;
      /* Clear previous initial step and force an initial step recomputation.
       * Normally, this would not occur after a reset, but it is necessary here
       * because the timestep used in one direction may not be suitable for the
       * other */
      ark_mem->h0u = ZERO;
      /* Reverse the step if in fixed mode. If adaptive, reset to 0 to clear any
       * old value from a call to ARKodeSetInit */
      ark_mem->hin = ark_mem->fixedstep ? -h : ZERO;

      /* Reset error controller (e.g., error and step size history) */
      if (ark_mem->hadapt_mem && ark_mem->hadapt_mem->hcontroller)
      {
        SUNErrCode err =
          SUNAdaptController_Reset(ark_mem->hadapt_mem->hcontroller);
        if (err != SUN_SUCCESS)
        {
          arkProcessError(ark_mem, ARK_CONTROLLER_ERR, __LINE__, __func__,
                          __FILE__, "Unable to reset error controller object");
          return ARK_CONTROLLER_ERR;
        }
      }
    }
  }

  if (ark_mem->step_setstepdirection != NULL)
  {
    return ark_mem->step_setstepdirection(ark_mem, stepdir);
  }

  return ARK_SUCCESS;
}

/*---------------------------------------------------------------
  ARKodeSetRootDirection:

  Specifies the direction of zero-crossings to be monitored.
  The default is to monitor both crossings.
  ---------------------------------------------------------------*/
int ARKodeSetRootDirection(void* arkode_mem, int* rootdir)
{
  ARKodeMem ark_mem;
  ARKodeRootMem ark_root_mem;
  int i;

  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  if (ark_mem->root_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_root_mem = (ARKodeRootMem)ark_mem->root_mem;

  if (ark_root_mem->nrtfn == 0)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_ROOT);
    return (ARK_ILL_INPUT);
  }
  for (i = 0; i < ark_root_mem->nrtfn; i++)
  {
    ark_root_mem->rootdir[i] = rootdir[i];
  }
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetNoInactiveRootWarn:

  Disables issuing a warning if some root function appears
  to be identically zero at the beginning of the integration
  ---------------------------------------------------------------*/
int ARKodeSetNoInactiveRootWarn(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeRootMem ark_root_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  if (ark_mem->root_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_ROOT);
    return (ARK_MEM_NULL);
  }
  ark_root_mem          = (ARKodeRootMem)ark_mem->root_mem;
  ark_root_mem->mxgnull = 0;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetPostprocessStepFn:

  Specifies a user-provided step postprocessing function having
  type ARKPostProcessFn.  A NULL input function disables step
  postprocessing.

  IF THE SUPPLIED FUNCTION MODIFIES ANY OF THE ACTIVE STATE DATA,
  THEN ALL THEORETICAL GUARANTEES OF SOLUTION ACCURACY AND
  STABILITY ARE LOST.
  ---------------------------------------------------------------*/
int ARKodeSetPostprocessStepFn(void* arkode_mem, ARKPostProcessFn ProcessStep)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* NULL argument sets default, otherwise set inputs */
  ark_mem->ProcessStep = ProcessStep;
  ark_mem->ps_data     = ark_mem->user_data;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetPostprocessStageFn:

  Specifies a user-provided stage postprocessing function having
  type ARKPostProcessFn.  A NULL input function disables
  stage postprocessing.

  IF THE SUPPLIED FUNCTION MODIFIES ANY OF THE ACTIVE STATE DATA,
  THEN ALL THEORETICAL GUARANTEES OF SOLUTION ACCURACY AND
  STABILITY ARE LOST.

  While it is possible to perform postprocessing when
  ARKodeSetDeduceImplicitRhs is enabled, this can cause implicit
  RHS evaluations to be inconsistent with the postprocessed stage
  values.  It is strongly recommended to disable
  ARKodeSetDeduceImplicitRhs in order to guarantee
  postprocessing constraints are enforced.
  ---------------------------------------------------------------*/
int ARKodeSetPostprocessStageFn(void* arkode_mem, ARKPostProcessFn ProcessStage)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* NULL argument sets default, otherwise set inputs */
  ark_mem->ProcessStage = ProcessStage;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetConstraints:

  Activates or Deactivates inequality constraint checking.
  ---------------------------------------------------------------*/
int ARKodeSetConstraints(void* arkode_mem, N_Vector constraints)
{
  sunrealtype temptest;
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive && (constraints != NULL))
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* If there are no constraints, destroy data structures */
  if (constraints == NULL)
  {
    arkFreeVec(ark_mem, &ark_mem->constraints);
    return (ARK_SUCCESS);
  }

  /* Test if required vector ops. are defined */
  if (constraints->ops->nvdiv == NULL || constraints->ops->nvmaxnorm == NULL ||
      constraints->ops->nvcompare == NULL || constraints->ops->nvprod == NULL ||
      constraints->ops->nvconstrmask == NULL ||
      constraints->ops->nvminquotient == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_BAD_NVECTOR);
    return (ARK_ILL_INPUT);
  }

  /* Check the constraints vector */
  temptest = N_VMaxNorm(constraints);
  if ((temptest > SUN_RCONST(2.5)) || (temptest < HALF))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_BAD_CONSTR);
    return (ARK_ILL_INPUT);
  }

  /* Allocate the internal constrains vector (if necessary) */
  if (!arkAllocVec(ark_mem, constraints, &ark_mem->constraints))
  {
    return (ARK_MEM_FAIL);
  }

  /* Load the constraints vector */
  N_VScale(ONE, constraints, ark_mem->constraints);

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetMaxNumConstrFails:

  Set max number of allowed constraint failures in a step before
  returning an error
  ---------------------------------------------------------------*/
int ARKodeSetMaxNumConstrFails(void* arkode_mem, int maxfails)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Passing maxfails = 0 sets the default, otherwise set to input */
  if (maxfails <= 0) { ark_mem->maxconstrfails = MAXCONSTRFAILS; }
  else { ark_mem->maxconstrfails = maxfails; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetCFLFraction:

  Specifies the safety factor to use on the maximum explicitly-
  stable step size.  Allowable values must be within the open
  interval (0,1).  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int ARKodeSetCFLFraction(void* arkode_mem, sunrealtype cfl_frac)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, __func__, &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* set positive-valued parameters, otherwise set default */
  if (cfl_frac <= ZERO) { hadapt_mem->cfl = CFLFAC; }
  else { hadapt_mem->cfl = cfl_frac; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetAdaptivityAdjustment:

  Adjusts the method order supplied to the temporal adaptivity
  controller.  For example, if the user expects order reduction
  due to problem stiffness, they may request that the controller
  assume a reduced order of accuracy for the method by specifying
  a value adjust < 0.
  ---------------------------------------------------------------*/
int ARKodeSetAdaptivityAdjustment(void* arkode_mem, int adjust)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, __func__, &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* store requested adjustment */
  hadapt_mem->adjust = adjust;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetSafetyFactor:

  Specifies the safety factor to use on the error-based predicted
  time step size.  Allowable values must be within the open
  interval (0,1).  A non-positive input implies a reset to the
  default value.
  ---------------------------------------------------------------*/
int ARKodeSetSafetyFactor(void* arkode_mem, sunrealtype safety)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, __func__, &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* check for allowable parameters */
  if (safety > ONE)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Illegal safety factor");
    return (ARK_ILL_INPUT);
  }

  /* set positive-valued parameters, otherwise set default */
  if (safety <= ZERO) { hadapt_mem->safety = SAFETY; }
  else { hadapt_mem->safety = safety; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetErrorBias:

  Specifies the error bias to use when performing adaptive-step
  error control.  Allowable values must be >= 1.0.  Any illegal
  value implies a reset to the default value.
  ---------------------------------------------------------------*/
int ARKodeSetErrorBias(void* arkode_mem, sunrealtype bias)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, __func__, &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Return an error if there is not a current SUNAdaptController */
  if (ark_mem->hadapt_mem->hcontroller == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__,
                    __FILE__, "SUNAdaptController NULL -- must be set before setting the error bias");
    return (ARK_MEM_NULL);
  }

  /* set allowed value, otherwise set default */
  if (bias < ONE)
  {
    retval = SUNAdaptController_SetErrorBias(hadapt_mem->hcontroller, -ONE);
  }
  else
  {
    retval = SUNAdaptController_SetErrorBias(hadapt_mem->hcontroller, bias);
  }
  if (retval != SUN_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_CONTROLLER_ERR, __LINE__, __func__, __FILE__,
                    "SUNAdaptController_SetErrorBias failure");
    return (ARK_CONTROLLER_ERR);
  }
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetMaxGrowth:

  Specifies the maximum step size growth factor to be allowed
  between successive integration steps.  Note: the first step uses
  a separate maximum growth factor.  Allowable values must be
  > 1.0.  Any illegal value implies a reset to the default.
  ---------------------------------------------------------------*/
int ARKodeSetMaxGrowth(void* arkode_mem, sunrealtype mx_growth)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, __func__, &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* set allowed value, otherwise set default */
  if (mx_growth <= ONE) { hadapt_mem->growth = GROWTH; }
  else { hadapt_mem->growth = mx_growth; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetMinReduction:

  Specifies the minimum possible step size reduction factor to be
  allowed between successive integration steps. Allowable values
  must be > 0.0 and < 1.0. Any illegal value implies a reset to
  the default.
  ---------------------------------------------------------------*/
int ARKodeSetMinReduction(void* arkode_mem, sunrealtype eta_min)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, __func__, &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* set allowed value, otherwise set default */
  if (eta_min >= ONE || eta_min <= ZERO) { hadapt_mem->etamin = ETAMIN; }
  else { hadapt_mem->etamin = eta_min; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetFixedStepBounds:

  Specifies the step size growth interval within which the step
  size will remain unchanged.  Allowable values must enclose the
  value 1.0.  Any illegal interval implies a reset to the default.
  ---------------------------------------------------------------*/
int ARKodeSetFixedStepBounds(void* arkode_mem, sunrealtype lb, sunrealtype ub)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, __func__, &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* set allowable interval, otherwise set defaults */
  if ((lb <= ONE) && (ub >= ONE))
  {
    hadapt_mem->lbound = lb;
    hadapt_mem->ubound = ub;
  }
  else
  {
    hadapt_mem->lbound = HFIXED_LB;
    hadapt_mem->ubound = HFIXED_UB;
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetMaxFirstGrowth:

  Specifies the user-provided time step adaptivity constant
  etamx1.  Legal values are greater than 1.0.  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int ARKodeSetMaxFirstGrowth(void* arkode_mem, sunrealtype etamx1)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, __func__, &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* if argument legal set it, otherwise set default */
  if (etamx1 <= ONE) { hadapt_mem->etamx1 = ETAMX1; }
  else { hadapt_mem->etamx1 = etamx1; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetMaxEFailGrowth:

  Specifies the user-provided time step adaptivity constant
  etamxf. Legal values are in the interval (0,1].  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int ARKodeSetMaxEFailGrowth(void* arkode_mem, sunrealtype etamxf)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, __func__, &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* if argument legal set it, otherwise set default */
  if ((etamxf <= ZERO) || (etamxf > ONE)) { hadapt_mem->etamxf = ETAMXF; }
  else { hadapt_mem->etamxf = etamxf; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetSmallNumEFails:

  Specifies the user-provided time step adaptivity constant
  small_nef.  Legal values are > 0.  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int ARKodeSetSmallNumEFails(void* arkode_mem, int small_nef)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, __func__, &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* if argument legal set it, otherwise set default */
  if (small_nef <= 0) { hadapt_mem->small_nef = SMALL_NEF; }
  else { hadapt_mem->small_nef = small_nef; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetMaxCFailGrowth:

  Specifies the user-provided time step adaptivity constant
  etacf. Legal values are in the interval (0,1].  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int ARKodeSetMaxCFailGrowth(void* arkode_mem, sunrealtype etacf)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, __func__, &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* if argument legal set it, otherwise set default */
  if ((etacf <= ZERO) || (etacf > ONE)) { hadapt_mem->etacf = ETACF; }
  else { hadapt_mem->etacf = etacf; }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetStabilityFn:

  Specifies the user-provided explicit time step stability
  function to use.  A NULL input function implies a reset to
  the default function (empty).
  ---------------------------------------------------------------*/
int ARKodeSetStabilityFn(void* arkode_mem, ARKExpStabFn EStab, void* estab_data)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, __func__, &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* NULL argument sets default, otherwise set inputs */
  hadapt_mem->expstab    = EStab;
  hadapt_mem->estab_data = estab_data;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetMaxErrTestFails:

  Specifies the maximum number of error test failures during one
  step try.  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int ARKodeSetMaxErrTestFails(void* arkode_mem, int maxnef)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* argument <= 0 sets default, otherwise set input */
  if (maxnef <= 0) { ark_mem->maxnef = MAXNEF; }
  else { ark_mem->maxnef = maxnef; }
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetMaxConvFails:

  Specifies the maximum number of nonlinear convergence failures
  during one step try.  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int ARKodeSetMaxConvFails(void* arkode_mem, int maxncf)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not need an algebraic solver */
  if (!ark_mem->step_supports_implicit)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not require an algebraic solver");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* argument <= 0 sets default, otherwise set input */
  if (maxncf <= 0) { ark_mem->maxncf = MAXNCF; }
  else { ark_mem->maxncf = maxncf; }
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSetAccumulatedErrorType:

  This routine sets the accumulated temporal error estimation
  strategy.
  ---------------------------------------------------------------*/
int ARKodeSetAccumulatedErrorType(void* arkode_mem, ARKAccumError accum_type)
{
  int retval = ARKodeResetAccumulatedError(arkode_mem);
  if (retval != ARK_SUCCESS) { return retval; }
  ((ARKodeMem)arkode_mem)->AccumErrorType = accum_type;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeResetAccumulatedError:

  This routine resets the accumulated temporal error estimate.
  ---------------------------------------------------------------*/
int ARKodeResetAccumulatedError(void* arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for non-adaptive time stepper modules */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support temporal adaptivity");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Reset value and counter, and return */
  ark_mem->AccumErrorStart = ark_mem->tn;
  ark_mem->AccumError      = ZERO;

  return (ARK_SUCCESS);
}

int ARKodeSetAdjointCheckpointScheme(void* arkode_mem,
                                     SUNAdjointCheckpointScheme checkpoint_scheme)

{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  ark_mem->checkpoint_scheme = checkpoint_scheme;

  return (ARK_SUCCESS);
}

int ARKodeSetAdjointCheckpointIndex(void* arkode_mem, suncountertype step_index)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  if (step_index < 0)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "step_index must be >= 0");
    return ARK_ILL_INPUT;
  }

  ark_mem->checkpoint_step_idx = step_index;

  return (ARK_SUCCESS);
}

int ARKodeSetUseCompensatedSums(void* arkode_mem, sunbooleantype onoff)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  ark_mem->use_compensated_sums = onoff;

  /* Call stepper routine (if provided) */
  if (ark_mem->step_setusecompensatedsums)
  {
    return ark_mem->step_setusecompensatedsums(arkode_mem, onoff);
  }

  return (ARK_SUCCESS);
}

/*===============================================================
  ARKODE optional output utility functions
  ===============================================================*/

/*---------------------------------------------------------------
  ARKodeGetNumStepAttempts:

  Returns the current number of RHS evaluations
  ---------------------------------------------------------------*/
int ARKodeGetNumRhsEvals(void* arkode_mem, int partition_index,
                         long int* num_rhs_evals)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return ARK_MEM_NULL;
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Call stepper routine (if provided) */
  if (ark_mem->step_getnumrhsevals)
  {
    return ark_mem->step_getnumrhsevals(arkode_mem, partition_index,
                                        num_rhs_evals);
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return ARK_STEPPER_UNSUPPORTED;
  }
}

/*---------------------------------------------------------------
  ARKodeGetNumStepAttempts:

   Returns the current number of steps attempted by the solver
  ---------------------------------------------------------------*/
int ARKodeGetNumStepAttempts(void* arkode_mem, long int* nstep_attempts)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  *nstep_attempts = ark_mem->nst_attempts;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetNumSteps:

  Returns the current number of integration steps
  ---------------------------------------------------------------*/
int ARKodeGetNumSteps(void* arkode_mem, long int* nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  *nsteps = ark_mem->nst;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetActualInitStep:

  Returns the step size used on the first step
  ---------------------------------------------------------------*/
int ARKodeGetActualInitStep(void* arkode_mem, sunrealtype* hinused)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  *hinused = ark_mem->h0u;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetLastStep:

  Returns the step size used on the last successful step
  ---------------------------------------------------------------*/
int ARKodeGetLastStep(void* arkode_mem, sunrealtype* hlast)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  *hlast = ark_mem->hold;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetCurrentStep:

  Returns the step size to be attempted on the next step
  ---------------------------------------------------------------*/
int ARKodeGetCurrentStep(void* arkode_mem, sunrealtype* hcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  *hcur = ark_mem->next_h;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetStepDirection:

  Gets the direction of integration (forward or backward) based
  on the sign of stepdir. A value of 0 indicates integration can
  proceed in either direction.
  ---------------------------------------------------------------*/
int ARKodeGetStepDirection(void* arkode_mem, sunrealtype* stepdir)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  if (stepdir == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "stepdir cannot be NULL");
  }

  *stepdir = (ark_mem->fixedstep || ark_mem->h == ZERO) ? ark_mem->hin
                                                        : ark_mem->h;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetCurrentState:

  Returns the current solution (before or after as step) or
  stage value (during step solve).
  ---------------------------------------------------------------*/
int ARKodeGetCurrentState(void* arkode_mem, N_Vector* state)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  *state = ark_mem->ycur;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetEstLocalErrors:

  Returns an estimate of the local error
  ---------------------------------------------------------------*/
int ARKodeGetEstLocalErrors(void* arkode_mem, N_Vector ele)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Call stepper-specific routine (if provided); otherwise return an error */
  if (ark_mem->step_getestlocalerrors)
  {
    return (ark_mem->step_getestlocalerrors(ark_mem, ele));
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does provide a temporal error estimate");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

/*---------------------------------------------------------------
  ARKodeGetCurrentTime:

  Returns the current value of the independent variable
  ---------------------------------------------------------------*/
int ARKodeGetCurrentTime(void* arkode_mem, sunrealtype* tcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  *tcur = ark_mem->tcur;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetCurrentGamma: Returns the current value of gamma
  ---------------------------------------------------------------*/
int ARKodeGetCurrentGamma(void* arkode_mem, sunrealtype* gamma)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not need an algebraic solver */
  if (!ark_mem->step_supports_implicit)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not require an algebraic solver");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Call stepper routine to compute the state (if provided) */
  if (ark_mem->step_getcurrentgamma)
  {
    return (ark_mem->step_getcurrentgamma(ark_mem, gamma));
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

/*---------------------------------------------------------------
  ARKodeGetTolScaleFactor:

  Returns a suggested factor for scaling tolerances
  ---------------------------------------------------------------*/
int ARKodeGetTolScaleFactor(void* arkode_mem, sunrealtype* tolsfact)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not use tolerances
     (i.e., neither supports adaptivity nor needs an algebraic solver) */
  if ((!ark_mem->step_supports_implicit) && (!ark_mem->step_supports_adaptive))
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not use tolerances");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  *tolsfact = ark_mem->tolsf;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetErrWeights:

  This routine returns the current error weight vector.
  ---------------------------------------------------------------*/
int ARKodeGetErrWeights(void* arkode_mem, N_Vector eweight)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not use tolerances
     (i.e., neither supports adaptivity nor needs an algebraic solver) */
  if ((!ark_mem->step_supports_implicit) && (!ark_mem->step_supports_adaptive))
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not use tolerances");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  N_VScale(ONE, ark_mem->ewt, eweight);
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetResWeights:

  This routine returns the current residual weight vector.
  ---------------------------------------------------------------*/
int ARKodeGetResWeights(void* arkode_mem, N_Vector rweight)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for time steppers that do not support mass matrices */
  if (!ark_mem->step_supports_massmatrix)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support non-identity mass matrices");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  N_VScale(ONE, ark_mem->rwt, rweight);
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetWorkSpace:

  Returns integrator work space requirements
  ---------------------------------------------------------------*/
int ARKodeGetWorkSpace(void* arkode_mem, long int* lenrw, long int* leniw)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  *leniw = ark_mem->liw;
  *lenrw = ark_mem->lrw;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetNumGEvals:

  Returns the current number of calls to g (for rootfinding)
  ---------------------------------------------------------------*/
int ARKodeGetNumGEvals(void* arkode_mem, long int* ngevals)
{
  ARKodeMem ark_mem;
  ARKodeRootMem ark_root_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  if (ark_mem->root_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_root_mem = (ARKodeRootMem)ark_mem->root_mem;
  *ngevals     = ark_root_mem->nge;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetRootInfo:

  Returns pointer to array rootsfound showing roots found
  ---------------------------------------------------------------*/
int ARKodeGetRootInfo(void* arkode_mem, int* rootsfound)
{
  int i;
  ARKodeMem ark_mem;
  ARKodeRootMem ark_root_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;
  if (ark_mem->root_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_root_mem = (ARKodeRootMem)ark_mem->root_mem;
  for (i = 0; i < ark_root_mem->nrtfn; i++)
  {
    rootsfound[i] = ark_root_mem->iroots[i];
  }
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetStepStats:

  Returns step statistics
  ---------------------------------------------------------------*/
int ARKodeGetStepStats(void* arkode_mem, long int* nsteps, sunrealtype* hinused,
                       sunrealtype* hlast, sunrealtype* hcur, sunrealtype* tcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  *nsteps  = ark_mem->nst;
  *hinused = ark_mem->h0u;
  *hlast   = ark_mem->hold;
  *hcur    = ark_mem->next_h;
  *tcur    = ark_mem->tcur;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetAccumulatedError:

  This routine returns the accumulated temporal error estimate.
  ---------------------------------------------------------------*/
int ARKodeGetAccumulatedError(void* arkode_mem, sunrealtype* accum_error)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Return an error if the stepper cannot accumulate temporal error */
  if (!ark_mem->step_supports_adaptive)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support accumulated error estimation");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Get time since last accumulated error reset */
  sunrealtype time_interval = ark_mem->tcur - ark_mem->AccumErrorStart;

  /* Fill output based on error accumulation type */
  if (ark_mem->AccumErrorType == ARK_ACCUMERROR_MAX)
  {
    *accum_error = ark_mem->AccumError * ark_mem->reltol;
  }
  else if (ark_mem->AccumErrorType == ARK_ACCUMERROR_SUM)
  {
    *accum_error = ark_mem->AccumError * ark_mem->reltol;
  }
  else if (ark_mem->AccumErrorType == ARK_ACCUMERROR_AVG)
  {
    *accum_error = ark_mem->AccumError * ark_mem->reltol / time_interval;
  }
  else
  {
    arkProcessError(ark_mem, ARK_WARNING, __LINE__, __func__, __FILE__,
                    "temporal error accumulation is currently disabled");
    return (ARK_WARNING);
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetNumConstrFails:

  Returns the current number of constraint fails
  ---------------------------------------------------------------*/
int ARKodeGetNumConstrFails(void* arkode_mem, long int* nconstrfails)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  *nconstrfails = ark_mem->nconstrfails;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetNumExpSteps:

  Returns the current number of stability-limited steps
  ---------------------------------------------------------------*/
int ARKodeGetNumExpSteps(void* arkode_mem, long int* nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  *nsteps = ark_mem->hadapt_mem->nst_exp;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetNumAccSteps:

  Returns the current number of accuracy-limited steps
  ---------------------------------------------------------------*/
int ARKodeGetNumAccSteps(void* arkode_mem, long int* nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  *nsteps = ark_mem->hadapt_mem->nst_acc;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetNumErrTestFails:

  Returns the current number of error test failures
  ---------------------------------------------------------------*/
int ARKodeGetNumErrTestFails(void* arkode_mem, long int* netfails)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  *netfails = ark_mem->netf;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeComputeState:

  Computes y based on the current prediction and a given
  correction.
  ---------------------------------------------------------------*/
int ARKodeComputeState(void* arkode_mem, N_Vector zcor, N_Vector z)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for incompatible time stepper modules */
  if (!ark_mem->step_supports_implicit)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support algebraic solvers");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Call stepper routine to compute the state (if provided) */
  if (ark_mem->step_computestate)
  {
    return (ark_mem->step_computestate(ark_mem, zcor, z));
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

/*---------------------------------------------------------------
  ARKodeGetNonlinearSystemData:

  This routine provides access to the relevant data needed to
  compute the nonlinear system function.
  ---------------------------------------------------------------*/
int ARKodeGetNonlinearSystemData(void* arkode_mem, sunrealtype* tcur,
                                 N_Vector* zpred, N_Vector* z, N_Vector* Fi,
                                 sunrealtype* gamma, N_Vector* sdata,
                                 void** user_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Guard against use for incompatible time stepper modules */
  if (!ark_mem->step_supports_implicit)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support algebraic solvers");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  /* Call stepper routine to compute the state (if provided) */
  if (ark_mem->step_getnonlinearsystemdata)
  {
    return (ark_mem->step_getnonlinearsystemdata(ark_mem, tcur, zpred, z, Fi,
                                                 gamma, sdata, user_data));
  }
  else
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__,
                    "time-stepping module does not support this function");
    return (ARK_STEPPER_UNSUPPORTED);
  }
}

/*---------------------------------------------------------------
  ARKodeGetNumNonlinSolvIters:

  Returns the current number of nonlinear solver iterations
  ---------------------------------------------------------------*/
int ARKodeGetNumNonlinSolvIters(void* arkode_mem, long int* nniters)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Call stepper routine to compute the state (if provided) */
  if (ark_mem->step_getnumnonlinsolviters)
  {
    return (ark_mem->step_getnumnonlinsolviters(ark_mem, nniters));
  }
  else
  {
    *nniters = 0;
    return (ARK_SUCCESS);
  }
}

/*---------------------------------------------------------------
  ARKodeGetNumNonlinSolvConvFails:

  Returns the current number of nonlinear solver convergence fails
  ---------------------------------------------------------------*/
int ARKodeGetNumNonlinSolvConvFails(void* arkode_mem, long int* nnfails)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Call stepper routine to compute the state (if provided) */
  if (ark_mem->step_getnumnonlinsolvconvfails)
  {
    return (ark_mem->step_getnumnonlinsolvconvfails(ark_mem, nnfails));
  }
  else
  {
    *nnfails = 0;
    return (ARK_SUCCESS);
  }
}

/*---------------------------------------------------------------
  ARKodeGetNonlinSolvStats:

  Returns nonlinear solver statistics
  ---------------------------------------------------------------*/
int ARKodeGetNonlinSolvStats(void* arkode_mem, long int* nniters,
                             long int* nnfails)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Call stepper routine to compute the state (if provided) */
  if (ark_mem->step_getnonlinsolvstats)
  {
    return (ark_mem->step_getnonlinsolvstats(ark_mem, nniters, nnfails));
  }
  else
  {
    *nniters = *nnfails = 0;
    return (ARK_SUCCESS);
  }
}

/*---------------------------------------------------------------
  ARKodeGetNumStepSolveFails:

  Returns the current number of failed steps due to an algebraic
  solver convergence failure.
  ---------------------------------------------------------------*/
int ARKodeGetNumStepSolveFails(void* arkode_mem, long int* nncfails)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  *nncfails = ark_mem->ncfn;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeGetNumLinSolvSetups:

  Returns the current number of calls to the lsetup routine
  ---------------------------------------------------------------*/
int ARKodeGetNumLinSolvSetups(void* arkode_mem, long int* nlinsetups)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Call stepper routine to compute the state (if provided) */
  if (ark_mem->step_getnumlinsolvsetups)
  {
    return (ark_mem->step_getnumlinsolvsetups(ark_mem, nlinsetups));
  }
  else
  {
    *nlinsetups = 0;
    return (ARK_SUCCESS);
  }
}

/*---------------------------------------------------------------
  ARKodeGetUserData:

  Returns the user data pointer
  ---------------------------------------------------------------*/
int ARKodeGetUserData(void* arkode_mem, void** user_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  *user_data = ark_mem->user_data;

  return (ARK_SUCCESS);
}

/*-----------------------------------------------------------------
  ARKodePrintAllStats

  Prints the current value of all statistics
  ---------------------------------------------------------------*/

int ARKodePrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeRootMem ark_root_mem;

  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  if (fmt != SUN_OUTPUTFORMAT_TABLE && fmt != SUN_OUTPUTFORMAT_CSV)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Invalid formatting option.");
    return (ARK_ILL_INPUT);
  }

  sunfprintf_real(outfile, fmt, SUNTRUE, "Current time", ark_mem->tcur);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Steps", ark_mem->nst);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Step attempts", ark_mem->nst_attempts);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Stability limited steps",
                  ark_mem->hadapt_mem->nst_exp);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Accuracy limited steps",
                  ark_mem->hadapt_mem->nst_acc);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Error test fails", ark_mem->netf);
  sunfprintf_long(outfile, fmt, SUNFALSE, "NLS step fails", ark_mem->ncfn);
  sunfprintf_long(outfile, fmt, SUNFALSE, "Inequality constraint fails",
                  ark_mem->nconstrfails);
  sunfprintf_real(outfile, fmt, SUNFALSE, "Initial step size", ark_mem->h0u);
  sunfprintf_real(outfile, fmt, SUNFALSE, "Last step size", ark_mem->hold);
  sunfprintf_real(outfile, fmt, SUNFALSE, "Current step size", ark_mem->next_h);
  if (ark_mem->root_mem)
  {
    ark_root_mem = (ARKodeRootMem)ark_mem->root_mem;
    sunfprintf_long(outfile, fmt, SUNFALSE, "Root fn evals", ark_root_mem->nge);
  }

  /* Print relaxation stats */
  if (ark_mem->relax_enabled)
  {
    retval = arkRelaxPrintAllStats(ark_mem, outfile, fmt);
    if (retval != ARK_SUCCESS) { return (retval); }
  }

  /* Print stepper stats (if provided) */
  if (ark_mem->step_printallstats)
  {
    return (ark_mem->step_printallstats(ark_mem, outfile, fmt));
  }

  return (ARK_SUCCESS);
}

/*-----------------------------------------------------------------*/

char* ARKodeGetReturnFlagName(long int flag)
{
  char* name;
  name = (char*)malloc(27 * sizeof(char));

  switch (flag)
  {
  case ARK_SUCCESS: sprintf(name, "ARK_SUCCESS"); break;
  case ARK_TSTOP_RETURN: sprintf(name, "ARK_TSTOP_RETURN"); break;
  case ARK_ROOT_RETURN: sprintf(name, "ARK_ROOT_RETURN"); break;
  case ARK_WARNING: sprintf(name, "ARK_WARNING"); break;
  case ARK_TOO_MUCH_WORK: sprintf(name, "ARK_TOO_MUCH_WORK"); break;
  case ARK_TOO_MUCH_ACC: sprintf(name, "ARK_TOO_MUCH_ACC"); break;
  case ARK_ERR_FAILURE: sprintf(name, "ARK_ERR_FAILURE"); break;
  case ARK_CONV_FAILURE: sprintf(name, "ARK_CONV_FAILURE"); break;
  case ARK_LINIT_FAIL: sprintf(name, "ARK_LINIT_FAIL"); break;
  case ARK_LSETUP_FAIL: sprintf(name, "ARK_LSETUP_FAIL"); break;
  case ARK_LSOLVE_FAIL: sprintf(name, "ARK_LSOLVE_FAIL"); break;
  case ARK_RHSFUNC_FAIL: sprintf(name, "ARK_RHSFUNC_FAIL"); break;
  case ARK_FIRST_RHSFUNC_ERR: sprintf(name, "ARK_FIRST_RHSFUNC_ERR"); break;
  case ARK_REPTD_RHSFUNC_ERR: sprintf(name, "ARK_REPTD_RHSFUNC_ERR"); break;
  case ARK_UNREC_RHSFUNC_ERR: sprintf(name, "ARK_UNREC_RHSFUNC_ERR"); break;
  case ARK_RTFUNC_FAIL: sprintf(name, "ARK_RTFUNC_FAIL"); break;
  case ARK_LFREE_FAIL: sprintf(name, "ARK_LFREE_FAIL"); break;
  case ARK_MASSINIT_FAIL: sprintf(name, "ARK_MASSINIT_FAIL"); break;
  case ARK_MASSSETUP_FAIL: sprintf(name, "ARK_MASSSETUP_FAIL"); break;
  case ARK_MASSSOLVE_FAIL: sprintf(name, "ARK_MASSSOLVE_FAIL"); break;
  case ARK_MASSFREE_FAIL: sprintf(name, "ARK_MASSFREE_FAIL"); break;
  case ARK_MASSMULT_FAIL: sprintf(name, "ARK_MASSMULT_FAIL"); break;
  case ARK_CONSTR_FAIL: sprintf(name, "ARK_CONSTR_FAIL"); break;
  case ARK_MEM_FAIL: sprintf(name, "ARK_MEM_FAIL"); break;
  case ARK_MEM_NULL: sprintf(name, "ARK_MEM_NULL"); break;
  case ARK_ILL_INPUT: sprintf(name, "ARK_ILL_INPUT"); break;
  case ARK_NO_MALLOC: sprintf(name, "ARK_NO_MALLOC"); break;
  case ARK_BAD_K: sprintf(name, "ARK_BAD_K"); break;
  case ARK_BAD_T: sprintf(name, "ARK_BAD_T"); break;
  case ARK_BAD_DKY: sprintf(name, "ARK_BAD_DKY"); break;
  case ARK_TOO_CLOSE: sprintf(name, "ARK_TOO_CLOSE"); break;
  case ARK_VECTOROP_ERR: sprintf(name, "ARK_VECTOROP_ERR"); break;
  case ARK_NLS_INIT_FAIL: sprintf(name, "ARK_NLS_INIT_FAIL"); break;
  case ARK_NLS_SETUP_FAIL: sprintf(name, "ARK_NLS_SETUP_FAIL"); break;
  case ARK_NLS_SETUP_RECVR: sprintf(name, "ARK_NLS_SETUP_RECVR"); break;
  case ARK_NLS_OP_ERR: sprintf(name, "ARK_NLS_OP_ERR"); break;
  case ARK_INNERSTEP_ATTACH_ERR:
    sprintf(name, "ARK_INNERSTEP_ATTACH_ERR");
    break;
  case ARK_INNERSTEP_FAIL: sprintf(name, "ARK_INNERSTEP_FAIL"); break;
  case ARK_OUTERTOINNER_FAIL: sprintf(name, "ARK_OUTERTOINNER_FAIL"); break;
  case ARK_INNERTOOUTER_FAIL: sprintf(name, "ARK_INNERTOOUTER_FAIL"); break;
  case ARK_POSTPROCESS_STEP_FAIL:
    sprintf(name, "ARK_POSTPROCESS_STEP_FAIL");
    break;
  case ARK_POSTPROCESS_STAGE_FAIL:
    sprintf(name, "ARK_POSTPROCESS_STAGE_FAIL");
    break;
  case ARK_USER_PREDICT_FAIL: sprintf(name, "ARK_USER_PREDICT_FAIL"); break;
  case ARK_INTERP_FAIL: sprintf(name, "ARK_INTERP_FAIL"); break;
  case ARK_INVALID_TABLE: sprintf(name, "ARK_INVALID_TABLE"); break;
  case ARK_CONTEXT_ERR: sprintf(name, "ARK_CONTEXT_ERR"); break;
  case ARK_RELAX_FAIL: sprintf(name, "ARK_RELAX_FAIL"); break;
  case ARK_RELAX_MEM_NULL: sprintf(name, "ARK_RELAX_MEM_NULL"); break;
  case ARK_RELAX_FUNC_FAIL: sprintf(name, "ARK_RELAX_FUNC_FAIL"); break;
  case ARK_RELAX_JAC_FAIL: sprintf(name, "ARK_RELAX_JAC_FAIL"); break;
  case ARK_CONTROLLER_ERR: sprintf(name, "ARK_CONTROLLER_ERR"); break;
  case ARK_STEPPER_UNSUPPORTED: sprintf(name, "ARK_STEPPER_UNSUPPORTED"); break;
  case ARK_ADJ_RECOMPUTE_FAIL: sprintf(name, "ARK_ADJ_RECOMPUTE_FAIL"); break;
  case ARK_ADJ_CHECKPOINT_FAIL: sprintf(name, "ARK_ADJ_CHECKPOINT_FAIL"); break;
  case ARK_SUNADJSTEPPER_ERR: sprintf(name, "ARK_SUNADJSTEPPER_ERR"); break;
  case ARK_DOMEIG_FAIL: sprintf(name, "ARK_DOMEIG_FAIL"); break;
  case ARK_MAX_STAGE_LIMIT_FAIL:
    sprintf(name, "ARK_MAX_STAGE_LIMIT_FAIL");
    break;
  case ARK_SUNSTEPPER_ERR: sprintf(name, "ARK_SUNSTEPPER_ERR"); break;
  case ARK_STEP_DIRECTION_ERR: sprintf(name, "ARK_STEP_DIRECTION_ERR"); break;
  case ARK_UNRECOGNIZED_ERROR: sprintf(name, "ARK_UNRECOGNIZED_ERROR"); break;
  default: sprintf(name, "NONE");
  }

  return (name);
}

/*===============================================================
  ARKODE parameter output utility routine
  ===============================================================*/

/*---------------------------------------------------------------
  ARKodeWriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int ARKodeWriteParameters(void* arkode_mem, FILE* fp)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* print integrator parameters to file */
  fprintf(fp, "ARKODE solver parameters:\n");
  if (ark_mem->hmin != ZERO)
  {
    fprintf(fp, "  Minimum step size = " SUN_FORMAT_G "\n", ark_mem->hmin);
  }
  if (ark_mem->hmax_inv != ZERO)
  {
    fprintf(fp, "  Maximum step size = " SUN_FORMAT_G "\n",
            ONE / ark_mem->hmax_inv);
  }
  if (ark_mem->fixedstep) { fprintf(fp, "  Fixed time-stepping enabled\n"); }
  if (ark_mem->itol == ARK_WF)
  {
    fprintf(fp, "  User provided error weight function\n");
  }
  else
  {
    fprintf(fp, "  Solver relative tolerance = " SUN_FORMAT_G "\n",
            ark_mem->reltol);
    if (ark_mem->itol == ARK_SS)
    {
      fprintf(fp, "  Solver absolute tolerance = " SUN_FORMAT_G "\n",
              ark_mem->Sabstol);
    }
    else { fprintf(fp, "  Vector-valued solver absolute tolerance\n"); }
  }
  if (!ark_mem->rwt_is_ewt)
  {
    if (ark_mem->ritol == ARK_WF)
    {
      fprintf(fp, "  User provided residual weight function\n");
    }
    else
    {
      if (ark_mem->ritol == ARK_SS)
      {
        fprintf(fp, "  Absolute residual tolerance = " SUN_FORMAT_G "\n",
                ark_mem->SRabstol);
      }
      else { fprintf(fp, "  Vector-valued residual absolute tolerance\n"); }
    }
  }
  if (ark_mem->hin != ZERO)
  {
    fprintf(fp, "  Initial step size = " SUN_FORMAT_G "\n", ark_mem->hin);
  }
  fprintf(fp, "\n");
  fprintf(fp, "  Maximum step increase (first step) = " SUN_FORMAT_G "\n",
          ark_mem->hadapt_mem->etamx1);
  fprintf(fp,
          "  Step reduction factor on multiple error fails = " SUN_FORMAT_G "\n",
          ark_mem->hadapt_mem->etamxf);
  fprintf(fp, "  Minimum error fails before above factor is used = %i\n",
          ark_mem->hadapt_mem->small_nef);
  fprintf(fp, "  Step reduction factor on nonlinear convergence failure = " SUN_FORMAT_G "\n",
          ark_mem->hadapt_mem->etacf);
  fprintf(fp, "  Explicit safety factor = " SUN_FORMAT_G "\n",
          ark_mem->hadapt_mem->cfl);
  fprintf(fp, "  Safety factor = " SUN_FORMAT_G "\n",
          ark_mem->hadapt_mem->safety);
  fprintf(fp, "  Growth factor = " SUN_FORMAT_G "\n",
          ark_mem->hadapt_mem->growth);
  fprintf(fp, "  Step growth lower bound = " SUN_FORMAT_G "\n",
          ark_mem->hadapt_mem->lbound);
  fprintf(fp, "  Step growth upper bound = " SUN_FORMAT_G "\n",
          ark_mem->hadapt_mem->ubound);
  if (ark_mem->hadapt_mem->expstab == NULL)
  {
    fprintf(fp, "  No explicit stability function supplied\n");
  }
  else { fprintf(fp, "  User provided explicit stability function\n"); }
  if (ark_mem->hadapt_mem->hcontroller != NULL)
  {
    (void)SUNAdaptController_Write(ark_mem->hadapt_mem->hcontroller, fp);
  }

  fprintf(fp, "  Maximum number of error test failures = %i\n", ark_mem->maxnef);
  fprintf(fp, "  Maximum number of convergence test failures = %i\n",
          ark_mem->maxncf);

  /* Call stepper routine (if provided) */
  if (ark_mem->step_writeparameters)
  {
    return (ark_mem->step_writeparameters(ark_mem, fp));
  }

  return (ARK_SUCCESS);
}

/*===============================================================
  ARKODE-IO internal utility functions
  ===============================================================*/

/*---------------------------------------------------------------
  arkReplaceAdaptController

  Replaces the current SUNAdaptController time step controller
  object. If a NULL-valued SUNAdaptController is input, the
  default will be re-enabled.
  ---------------------------------------------------------------*/
int arkReplaceAdaptController(ARKodeMem ark_mem, SUNAdaptController C,
                              sunbooleantype take_ownership)
{
  int retval;
  long int lenrw, leniw;

  /* Remove current SUNAdaptController object
     (delete if owned, and then nullify pointer) */
  if (ark_mem->hadapt_mem->owncontroller &&
      (ark_mem->hadapt_mem->hcontroller != NULL))
  {
    retval = SUNAdaptController_Space(ark_mem->hadapt_mem->hcontroller, &lenrw,
                                      &leniw);
    if (retval == SUN_SUCCESS)
    {
      ark_mem->liw -= leniw;
      ark_mem->lrw -= lenrw;
    }

    retval = SUNAdaptController_Destroy(ark_mem->hadapt_mem->hcontroller);
    ark_mem->hadapt_mem->owncontroller = SUNFALSE;
    if (retval != SUN_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "SUNAdaptController_Destroy failure");
      return (ARK_MEM_FAIL);
    }
  }
  ark_mem->hadapt_mem->hcontroller = NULL;

  /* On NULL-valued input, create default SUNAdaptController object */
  if (C == NULL)
  {
    C = SUNAdaptController_I(ark_mem->sunctx);
    if (C == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "SUNAdaptController_I allocation failure");
      return (ARK_MEM_FAIL);
    }
    ark_mem->hadapt_mem->owncontroller = SUNTRUE;
  }
  else { ark_mem->hadapt_mem->owncontroller = take_ownership; }

  /* Attach new SUNAdaptController object */
  retval = SUNAdaptController_Space(C, &lenrw, &leniw);
  if (retval == SUN_SUCCESS)
  {
    ark_mem->liw += leniw;
    ark_mem->lrw += lenrw;
  }
  ark_mem->hadapt_mem->hcontroller = C;

  return (ARK_SUCCESS);
}

/*===============================================================
  ARKODE + XBraid interface utility functions
  ===============================================================*/

/*---------------------------------------------------------------
  arkSetForcePass:

  Ignore the value of kflag after the temporal error test and
  force the step to pass.
  ---------------------------------------------------------------*/
int arkSetForcePass(void* arkode_mem, sunbooleantype force_pass)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  ark_mem->force_pass = force_pass;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkGetLastKFlag:

  The last kflag value returned by the temporal error test.
  ---------------------------------------------------------------*/
int arkGetLastKFlag(void* arkode_mem, int* last_kflag)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  *last_kflag = ark_mem->last_kflag;

  return (ARK_SUCCESS);
}

/*===============================================================
  Deprecated functions
  ===============================================================*/

/*---------------------------------------------------------------
  arkSetAdaptivityMethod:

  Specifies the built-in time step adaptivity algorithm (and
  optionally, its associated parameters) to use.  All parameters
  will be checked for validity when used by the solver.

  Users should transition to constructing non-default SUNAdaptController
  objects directly, and providing those directly to the integrator
  via the time-stepping module *SetController routines.
  ---------------------------------------------------------------*/
int arkSetAdaptivityMethod(void* arkode_mem, int imethod, int idefault, int pq,
                           sunrealtype adapt_params[3])
{
  int retval;
  long int lenrw, leniw;
  sunrealtype k1, k2, k3;
  ARKodeMem ark_mem;
  SUNAdaptController C;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Check for illegal inputs */
  if ((idefault != 1) && (adapt_params == NULL))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "NULL-valued adapt_params provided");
    return (ARK_ILL_INPUT);
  }

  /* Remove current SUNAdaptController object
     (delete if owned, and then nullify pointer) */
  if (ark_mem->hadapt_mem->owncontroller &&
      (ark_mem->hadapt_mem->hcontroller != NULL))
  {
    retval = SUNAdaptController_Space(ark_mem->hadapt_mem->hcontroller, &lenrw,
                                      &leniw);
    if (retval == SUN_SUCCESS)
    {
      ark_mem->liw -= leniw;
      ark_mem->lrw -= lenrw;
    }

    retval = SUNAdaptController_Destroy(ark_mem->hadapt_mem->hcontroller);
    ark_mem->hadapt_mem->owncontroller = SUNFALSE;
    if (retval != SUN_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "SUNAdaptController_Destroy failure");
      return (ARK_MEM_FAIL);
    }
  }
  ark_mem->hadapt_mem->hcontroller = NULL;

  /* set adaptivity parameters from inputs */
  k1 = k2 = k3 = ZERO;
  if (idefault != 1)
  {
    k1 = adapt_params[0];
    k2 = adapt_params[1];
    k3 = adapt_params[2];
  }
  ark_mem->hadapt_mem->pq = pq;

  /* Create new SUNAdaptController object based on "imethod" input, optionally setting
     the specified controller parameters */
  C = NULL;
  switch (imethod)
  {
  case (ARK_ADAPT_PID):
    C = SUNAdaptController_PID(ark_mem->sunctx);
    if (C == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "SUNAdaptController_PID allocation failure");
      return (ARK_MEM_FAIL);
    }
    if (idefault != 1)
    {
      retval = SUNAdaptController_SetParams_PID(C, k1, -k2, k3);
      if (retval != SUN_SUCCESS)
      {
        (void)SUNAdaptController_Destroy(C);
        arkProcessError(ark_mem, ARK_CONTROLLER_ERR, __LINE__, __func__,
                        __FILE__, "SUNAdaptController_SetParams_PID failure");
        return (ARK_CONTROLLER_ERR);
      }
    }
    break;
  case (ARK_ADAPT_PI):
    C = SUNAdaptController_PI(ark_mem->sunctx);
    if (C == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "SUNAdaptController_PI allocation failure");
      return (ARK_MEM_FAIL);
    }
    if (idefault != 1)
    {
      retval = SUNAdaptController_SetParams_PI(C, k1, -k2);
      if (retval != SUN_SUCCESS)
      {
        (void)SUNAdaptController_Destroy(C);
        arkProcessError(ark_mem, ARK_CONTROLLER_ERR, __LINE__, __func__,
                        __FILE__, "SUNAdaptController_SetParams_PI failure");
        return (ARK_CONTROLLER_ERR);
      }
    }
    break;
  case (ARK_ADAPT_I):
    C = SUNAdaptController_I(ark_mem->sunctx);
    if (C == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "SUNAdaptController_I allocation failure");
      return (ARK_MEM_FAIL);
    }
    if (idefault != 1)
    {
      retval = SUNAdaptController_SetParams_I(C, k1);
      if (retval != SUN_SUCCESS)
      {
        (void)SUNAdaptController_Destroy(C);
        arkProcessError(ark_mem, ARK_CONTROLLER_ERR, __LINE__, __func__,
                        __FILE__, "SUNAdaptController_SetParams_I failure");
        return (ARK_CONTROLLER_ERR);
      }
    }
    break;
  case (ARK_ADAPT_EXP_GUS):
    C = SUNAdaptController_ExpGus(ark_mem->sunctx);
    if (C == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "SUNAdaptController_ExpGus allocation failure");
      return (ARK_MEM_FAIL);
    }
    if (idefault != 1)
    {
      retval = SUNAdaptController_SetParams_ExpGus(C, k1, k2);
      if (retval != SUN_SUCCESS)
      {
        (void)SUNAdaptController_Destroy(C);
        arkProcessError(ark_mem, ARK_CONTROLLER_ERR, __LINE__, __func__,
                        __FILE__, "SUNAdaptController_SetParams_ExpGus failure");
        return (ARK_CONTROLLER_ERR);
      }
    }
    break;
  case (ARK_ADAPT_IMP_GUS):
    C = SUNAdaptController_ImpGus(ark_mem->sunctx);
    if (C == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "SUNAdaptController_ImpGus allocation failure");
      return (ARK_MEM_FAIL);
    }
    if (idefault != 1)
    {
      retval = SUNAdaptController_SetParams_ImpGus(C, k1, k2);
      if (retval != SUN_SUCCESS)
      {
        (void)SUNAdaptController_Destroy(C);
        arkProcessError(ark_mem, ARK_CONTROLLER_ERR, __LINE__, __func__,
                        __FILE__, "SUNAdaptController_SetParams_ImpGus failure");
        return (ARK_CONTROLLER_ERR);
      }
    }
    break;
  case (ARK_ADAPT_IMEX_GUS):
    C = SUNAdaptController_ImExGus(ark_mem->sunctx);
    if (C == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "SUNAdaptController_ImExGus allocation failure");
      return (ARK_MEM_FAIL);
    }
    if (idefault != 1)
    {
      retval = SUNAdaptController_SetParams_ImExGus(C, k1, k2, k3, k3);
      if (retval != SUN_SUCCESS)
      {
        (void)SUNAdaptController_Destroy(C);
        arkProcessError(ark_mem, ARK_CONTROLLER_ERR, __LINE__, __func__,
                        __FILE__, "SUNAdaptController_SetParams_ImExGus failure");
        return (ARK_CONTROLLER_ERR);
      }
    }
    break;
  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Illegal imethod");
    return (ARK_ILL_INPUT);
  }

  /* Attach new SUNAdaptController object */
  retval = SUNAdaptController_Space(C, &lenrw, &leniw);
  if (retval == SUN_SUCCESS)
  {
    ark_mem->liw += leniw;
    ark_mem->lrw += lenrw;
  }
  ark_mem->hadapt_mem->hcontroller   = C;
  ark_mem->hadapt_mem->owncontroller = SUNTRUE;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkSetAdaptivityFn:

  Specifies the user-provided time step adaptivity function to use.
  If 'hfun' is NULL-valued, then the default I controller will
  be used instead.

  Users should transition to constructing a custom SUNAdaptController
  object, and providing this directly to the integrator
  via the time-stepping module *SetController routines.
  ---------------------------------------------------------------*/
int arkSetAdaptivityFn(void* arkode_mem, ARKAdaptFn hfun, void* h_data)
{
  int retval;
  long int lenrw, leniw;
  ARKodeMem ark_mem;
  SUNAdaptController C;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Remove current SUNAdaptController object
     (delete if owned, and then nullify pointer) */
  if (ark_mem->hadapt_mem->owncontroller &&
      (ark_mem->hadapt_mem->hcontroller != NULL))
  {
    retval = SUNAdaptController_Space(ark_mem->hadapt_mem->hcontroller, &lenrw,
                                      &leniw);
    if (retval == SUN_SUCCESS)
    {
      ark_mem->liw -= leniw;
      ark_mem->lrw -= lenrw;
    }

    retval = SUNAdaptController_Destroy(ark_mem->hadapt_mem->hcontroller);
    ark_mem->hadapt_mem->owncontroller = SUNFALSE;
    if (retval != SUN_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "SUNAdaptController_Destroy failure");
      return (ARK_MEM_FAIL);
    }
  }
  ark_mem->hadapt_mem->hcontroller = NULL;

  /* Create new SUNAdaptController object depending on NULL-ity of 'hfun' */
  C = NULL;
  if (hfun == NULL)
  {
    C = SUNAdaptController_I(ark_mem->sunctx);
    if (C == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "SUNAdaptController_I allocation failure");
      return (ARK_MEM_FAIL);
    }
  }
  else
  {
    C = ARKUserControl(ark_mem->sunctx, arkode_mem, hfun, h_data);
    if (C == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "ARKUserControl allocation failure");
      return (ARK_MEM_FAIL);
    }
  }

  /* Attach new SUNAdaptController object */
  retval = SUNAdaptController_Space(C, &lenrw, &leniw);
  if (retval == SUN_SUCCESS)
  {
    ark_mem->liw += leniw;
    ark_mem->lrw += lenrw;
  }
  ark_mem->hadapt_mem->hcontroller   = C;
  ark_mem->hadapt_mem->owncontroller = SUNTRUE;

  return (ARK_SUCCESS);
}

/*===============================================================
  EOF
  ===============================================================*/
