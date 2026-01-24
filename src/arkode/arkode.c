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
 * This is the implementation file for the main ARKODE
 * infrastructure.  It is independent of the ARKODE time step
 * module, nonlinear solver, linear solver and vector modules in
 * use.
 *--------------------------------------------------------------*/

#include "arkode/arkode.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include "arkode_impl.h"
#include "arkode_interp_impl.h"
#include "sundials/priv/sundials_errors_impl.h"
#include "sundials/sundials_context.h"
#include "sundials/sundials_logger.h"
#include "sundials_utils.h"

#include "sundials_macros.h"

/*===============================================================
  Exported functions
  ===============================================================*/

/*---------------------------------------------------------------
  ARKodeResize:

  ARKodeResize re-initializes ARKODE's memory for a problem with a
  changing vector size.  It is assumed that the problem dynamics
  before and after the vector resize will be comparable, so that
  all time-stepping heuristics prior to calling ARKodeResize
  remain valid after the call.  If instead the dynamics should be
  re-calibrated, the ARKODE memory structure should be deleted
  with a call to ARKodeFree, and re-created with a call to
  *StepCreate.

  To aid in the vector-resize operation, the user can supply a
  vector resize function, that will take as input an N_Vector with
  the previous size, and return as output a corresponding vector
  of the new size.  If this function (of type ARKVecResizeFn) is
  not supplied (i.e. is set to NULL), then all existing N_Vectors
  will be destroyed and re-cloned from the input vector.

  In the case that the dynamical time scale should be modified
  slightly from the previous time scale, an input "hscale" is
  allowed, that will re-scale the upcoming time step by the
  specified factor.  If a value <= 0 is specified, the default of
  1.0 will be used.

  Other arguments:
  ark_mem          Existing ARKODE memory data structure.
  y0               The newly-sized solution vector, holding
                   the current dependent variable values.
  t0               The current value of the independent
                   variable.
  resize_data      User-supplied data structure that will be
                   passed to the supplied resize function.

  The return value is ARK_SUCCESS = 0 if no errors occurred, or
  a negative value otherwise.
  ---------------------------------------------------------------*/
int ARKodeResize(void* arkode_mem, N_Vector y0, sunrealtype hscale,
                 sunrealtype t0, ARKVecResizeFn resize, void* resize_data)
{
  sunbooleantype resizeOK;
  sunindextype lrw1, liw1, lrw_diff, liw_diff;
  int retval;
  ARKodeMem ark_mem;

  /* Check ark_mem */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE)
  {
    arkProcessError(ark_mem, ARK_NO_MALLOC, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MALLOC);
    return (ARK_NO_MALLOC);
  }

  /* Check for legal input parameters */
  if (y0 == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_Y0);
    return (ARK_ILL_INPUT);
  }

  /* Copy the input parameters into ARKODE state */
  ark_mem->tcur = t0;
  ark_mem->tn   = t0;

  /* Update time-stepping parameters */
  /*   adjust upcoming step size depending on hscale */
  if (hscale <= ZERO) { hscale = ONE; }
  if (hscale != ONE)
  {
    /* Encode hscale into ark_mem structure */
    ark_mem->eta = hscale;
    ark_mem->hprime *= hscale;

    /* If next step would overtake tstop, adjust stepsize */
    if (ark_mem->tstopset)
    {
      if ((ark_mem->tcur + ark_mem->hprime - ark_mem->tstop) * ark_mem->hprime >
          ZERO)
      {
        ark_mem->hprime = (ark_mem->tstop - ark_mem->tcur) *
                          (ONE - FOUR * ark_mem->uround);
        ark_mem->eta = ark_mem->hprime / ark_mem->h;
      }
    }
  }

  /* Determine change in vector sizes */
  lrw1 = liw1 = 0;
  if (y0->ops->nvspace != NULL) { N_VSpace(y0, &lrw1, &liw1); }
  lrw_diff      = lrw1 - ark_mem->lrw1;
  liw_diff      = liw1 - ark_mem->liw1;
  ark_mem->lrw1 = lrw1;
  ark_mem->liw1 = liw1;

  /* Disable constraints, the user will need to set a new constraint vector for
     the updated problem size */
  arkFreeVec(ark_mem, &ark_mem->constraints);

  /* Resize the solver vectors (using y0 as a template) */
  resizeOK = arkResizeVectors(ark_mem, resize, resize_data, lrw_diff, liw_diff,
                              y0);
  if (!resizeOK)
  {
    arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Unable to resize vector");
    return (ARK_MEM_FAIL);
  }

  /* Resize the interpolation structure memory */
  if (ark_mem->interp != NULL)
  {
    retval = arkInterpResize(ark_mem, ark_mem->interp, resize, resize_data,
                             lrw_diff, liw_diff, y0);
    if (retval != ARK_SUCCESS)
    {
      arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                      "Interpolation module resize failure");
      return (retval);
    }
  }

  /* Copy y0 into ark_yn to set the current solution */
  N_VScale(ONE, y0, ark_mem->yn);
  ark_mem->fn_is_current = SUNFALSE;

  /* Indicate that problem needs to be initialized */
  ark_mem->initsetup  = SUNTRUE;
  ark_mem->init_type  = RESIZE_INIT;
  ark_mem->firststage = SUNTRUE;

  /* Call the stepper-specific resize (if provided) */
  if (ark_mem->step_resize)
  {
    return (ark_mem->step_resize(ark_mem, y0, hscale, t0, resize, resize_data));
  }

  /* Problem has been successfully re-sized */
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeReset:

  This routine resets an ARKode module to solve the same
  problem from the given time with the input state (all counter
  values are retained).
  ---------------------------------------------------------------*/
int ARKodeReset(void* arkode_mem, sunrealtype tR, N_Vector yR)
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

  /* Reset main ARKODE infrastructure */
  retval = arkInit(ark_mem, tR, yR, RESET_INIT);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "ARKode reset failure");
    return (retval);
  }

  /* Call stepper routine to perform remaining reset operations (if provided) */
  if (ark_mem->step_reset) { return (ark_mem->step_reset(ark_mem, tR, yR)); }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeSStolerances, ARKodeSVtolerances, ARKodeWFtolerances:

  These functions specify the integration tolerances. One of them
  SHOULD be called before the first call to ARKodeEvolve; otherwise
  default values of reltol=1e-4 and abstol=1e-9 will be used,
  which may be entirely incorrect for a specific problem.

  ARKodeSStolerances specifies scalar relative and absolute
  tolerances.

  ARKodeSVtolerances specifies scalar relative tolerance and a
  vector absolute tolerance (a potentially different absolute
  tolerance for each vector component).

  ARKodeWFtolerances specifies a user-provides function (of type
  ARKEwtFn) which will be called to set the error weight vector.
  ---------------------------------------------------------------*/
int ARKodeSStolerances(void* arkode_mem, sunrealtype reltol, sunrealtype abstol)
{
  /* unpack ark_mem */
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Check inputs */
  if (ark_mem->MallocDone == SUNFALSE)
  {
    arkProcessError(ark_mem, ARK_NO_MALLOC, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MALLOC);
    return (ARK_NO_MALLOC);
  }
  if (reltol < ZERO)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_BAD_RELTOL);
    return (ARK_ILL_INPUT);
  }
  if (abstol < ZERO)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_BAD_ABSTOL);
    return (ARK_ILL_INPUT);
  }

  /* Ensure that vector supports N_VAddConst */
  if (!ark_mem->tempv1->ops->nvaddconst)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "N_VAddConst unimplemented (required for scalar abstol)");
    return (ARK_ILL_INPUT);
  }

  /* Set flag indicating whether abstol == 0 */
  ark_mem->atolmin0 = (abstol == ZERO);

  /* Copy tolerances into memory */
  ark_mem->reltol  = reltol;
  ark_mem->Sabstol = abstol;
  ark_mem->itol    = ARK_SS;

  /* enforce use of arkEwtSetSS */
  ark_mem->user_efun = SUNFALSE;
  ark_mem->efun      = arkEwtSetSS;
  ark_mem->e_data    = ark_mem;

  return (ARK_SUCCESS);
}

int ARKodeSVtolerances(void* arkode_mem, sunrealtype reltol, N_Vector abstol)
{
  /* local variables */
  sunrealtype abstolmin;

  /* unpack ark_mem */
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Check inputs */
  if (ark_mem->MallocDone == SUNFALSE)
  {
    arkProcessError(ark_mem, ARK_NO_MALLOC, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MALLOC);
    return (ARK_NO_MALLOC);
  }
  if (reltol < ZERO)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_BAD_RELTOL);
    return (ARK_ILL_INPUT);
  }
  if (abstol == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_ABSTOL);
    return (ARK_ILL_INPUT);
  }
  if (abstol->ops->nvmin == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Missing N_VMin routine from N_Vector");
    return (ARK_ILL_INPUT);
  }
  abstolmin = N_VMin(abstol);
  if (abstolmin < ZERO)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_BAD_ABSTOL);
    return (ARK_ILL_INPUT);
  }

  /* Set flag indicating whether min(abstol) == 0 */
  ark_mem->atolmin0 = (abstolmin == ZERO);

  /* Copy tolerances into memory */
  if (!(ark_mem->VabstolMallocDone))
  {
    if (!arkAllocVec(ark_mem, ark_mem->ewt, &(ark_mem->Vabstol)))
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_ARKMEM_FAIL);
      return (ARK_ILL_INPUT);
    }
    ark_mem->VabstolMallocDone = SUNTRUE;
  }
  N_VScale(ONE, abstol, ark_mem->Vabstol);
  ark_mem->reltol = reltol;
  ark_mem->itol   = ARK_SV;

  /* enforce use of arkEwtSetSV */
  ark_mem->user_efun = SUNFALSE;
  ark_mem->efun      = arkEwtSetSV;
  ark_mem->e_data    = ark_mem;

  return (ARK_SUCCESS);
}

int ARKodeWFtolerances(void* arkode_mem, ARKEwtFn efun)
{
  /* unpack ark_mem */
  ARKodeMem ark_mem;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  if (ark_mem->MallocDone == SUNFALSE)
  {
    arkProcessError(ark_mem, ARK_NO_MALLOC, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MALLOC);
    return (ARK_NO_MALLOC);
  }

  /* Copy tolerance data into memory */
  ark_mem->itol      = ARK_WF;
  ark_mem->user_efun = SUNTRUE;
  ark_mem->efun      = efun;
  ark_mem->e_data    = ark_mem->user_data;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeResStolerance, ARKodeResVtolerance, ARKodeResFtolerance:

  These functions specify the absolute residual tolerance.
  Specification of the absolute residual tolerance is only
  necessary for problems with non-identity mass matrices in which
  the units of the solution vector y dramatically differ from the
  units of the ODE right-hand side f(t,y).  If this occurs, one
  of these routines SHOULD be called before the first call to
  ARKODE; otherwise the default value of rabstol=1e-9 will be
  used, which may be entirely incorrect for a specific problem.

  ARKodeResStolerances specifies a scalar residual tolerance.

  ARKodeResVtolerances specifies a vector residual tolerance
  (a potentially different absolute residual tolerance for
  each vector component).

  ARKodeResFtolerances specifies a user-provides function (of
  type ARKRwtFn) which will be called to set the residual
  weight vector.
  ---------------------------------------------------------------*/
int ARKodeResStolerance(void* arkode_mem, sunrealtype rabstol)
{
  /* unpack ark_mem */
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

  /* Check inputs */
  if (ark_mem->MallocDone == SUNFALSE)
  {
    arkProcessError(ark_mem, ARK_NO_MALLOC, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MALLOC);
    return (ARK_NO_MALLOC);
  }
  if (rabstol < ZERO)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_BAD_RABSTOL);
    return (ARK_ILL_INPUT);
  }

  /* Ensure that vector supports N_VAddConst */
  if (!ark_mem->tempv1->ops->nvaddconst)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "N_VAddConst unimplemented (required for scalar rabstol)");
    return (ARK_ILL_INPUT);
  }

  /* Set flag indicating whether rabstol == 0 */
  ark_mem->Ratolmin0 = (rabstol == ZERO);

  /* Allocate space for rwt if necessary */
  if (ark_mem->rwt_is_ewt)
  {
    ark_mem->rwt = NULL;
    if (!arkAllocVec(ark_mem, ark_mem->ewt, &(ark_mem->rwt)))
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_ARKMEM_FAIL);
      return (ARK_ILL_INPUT);
    }
    ark_mem->rwt_is_ewt = SUNFALSE;
  }

  /* Copy tolerances into memory */
  ark_mem->SRabstol = rabstol;
  ark_mem->ritol    = ARK_SS;

  /* enforce use of arkRwtSet */
  ark_mem->user_efun = SUNFALSE;
  ark_mem->rfun      = arkRwtSet;
  ark_mem->r_data    = ark_mem;

  return (ARK_SUCCESS);
}

int ARKodeResVtolerance(void* arkode_mem, N_Vector rabstol)
{
  /* local variables */
  sunrealtype rabstolmin;

  /* unpack ark_mem */
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

  /* Check inputs */
  if (ark_mem->MallocDone == SUNFALSE)
  {
    arkProcessError(ark_mem, ARK_NO_MALLOC, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MALLOC);
    return (ARK_NO_MALLOC);
  }
  if (rabstol == NULL)
  {
    arkProcessError(ark_mem, ARK_NO_MALLOC, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_RABSTOL);
    return (ARK_NO_MALLOC);
  }
  if (rabstol->ops->nvmin == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Missing N_VMin routine from N_Vector");
    return (ARK_ILL_INPUT);
  }
  rabstolmin = N_VMin(rabstol);
  if (rabstolmin < ZERO)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_BAD_RABSTOL);
    return (ARK_ILL_INPUT);
  }

  /* Set flag indicating whether min(abstol) == 0 */
  ark_mem->Ratolmin0 = (rabstolmin == ZERO);

  /* Allocate space for rwt if necessary */
  if (ark_mem->rwt_is_ewt)
  {
    ark_mem->rwt = NULL;
    if (!arkAllocVec(ark_mem, ark_mem->ewt, &(ark_mem->rwt)))
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_ARKMEM_FAIL);
      return (ARK_ILL_INPUT);
    }
    ark_mem->rwt_is_ewt = SUNFALSE;
  }

  /* Copy tolerances into memory */
  if (!(ark_mem->VRabstolMallocDone))
  {
    if (!arkAllocVec(ark_mem, ark_mem->rwt, &(ark_mem->VRabstol)))
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_ARKMEM_FAIL);
      return (ARK_ILL_INPUT);
    }
    ark_mem->VRabstolMallocDone = SUNTRUE;
  }
  N_VScale(ONE, rabstol, ark_mem->VRabstol);
  ark_mem->ritol = ARK_SV;

  /* enforce use of arkRwtSet */
  ark_mem->user_efun = SUNFALSE;
  ark_mem->rfun      = arkRwtSet;
  ark_mem->r_data    = ark_mem;

  return (ARK_SUCCESS);
}

int ARKodeResFtolerance(void* arkode_mem, ARKRwtFn rfun)
{
  /* unpack ark_mem */
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

  if (ark_mem->MallocDone == SUNFALSE)
  {
    arkProcessError(ark_mem, ARK_NO_MALLOC, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MALLOC);
    return (ARK_NO_MALLOC);
  }

  /* Allocate space for rwt if necessary */
  if (ark_mem->rwt_is_ewt)
  {
    ark_mem->rwt = NULL;
    if (!arkAllocVec(ark_mem, ark_mem->ewt, &(ark_mem->rwt)))
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_ARKMEM_FAIL);
      return (ARK_ILL_INPUT);
    }
    ark_mem->rwt_is_ewt = SUNFALSE;
  }

  /* Copy tolerance data into memory */
  ark_mem->ritol     = ARK_WF;
  ark_mem->user_rfun = SUNTRUE;
  ark_mem->rfun      = rfun;
  ark_mem->r_data    = ark_mem->user_data;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeEvolve:

  This routine is the main driver of ARKODE-based integrators.

  It integrates over a time interval defined by the user, by
  calling the time step module to do internal time steps.

  The first time that ARKodeEvolve is called for a successfully
  initialized problem, it computes a tentative initial step size.

  ARKodeEvolve supports two modes as specified by itask: ARK_NORMAL and
  ARK_ONE_STEP.  In the ARK_NORMAL mode, the solver steps until
  it reaches or passes tout and then interpolates to obtain
  y(tout).  In the ARK_ONE_STEP mode, it takes one internal step
  and returns.  The behavior of both modes can be over-rided
  through user-specification of ark_tstop (through the
  *StepSetStopTime function), in which case if a solver step
  would pass tstop, the step is shortened so that it stops at
  exactly the specified stop time, and hence interpolation of
  y(tout) is not required.
  ---------------------------------------------------------------*/
int ARKodeEvolve(void* arkode_mem, sunrealtype tout, N_Vector yout,
                 sunrealtype* tret, int itask)
{
  long int nstloc;
  int retval, kflag, istate, ir;
  int ewtsetOK;
  sunrealtype troundoff, nrm;
  sunbooleantype inactive_roots;
  sunrealtype dsm;
  int nflag, ncf, nef, constrfails;
  int relax_fails;
  ARKodeMem ark_mem;

  /* used only with debugging logging */
  SUNDIALS_MAYBE_UNUSED int attempts;

  /* Check and process inputs */

  /* Check if ark_mem exists */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE)
  {
    arkProcessError(ark_mem, ARK_NO_MALLOC, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MALLOC);
    return (ARK_NO_MALLOC);
  }

  /* Check for yout != NULL */
  if ((ark_mem->ycur = yout) == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_YOUT_NULL);
    return (ARK_ILL_INPUT);
  }

  /* Check for tret != NULL */
  if (tret == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_TRET_NULL);
    return (ARK_ILL_INPUT);
  }

  /* Check for valid itask */
  if ((itask != ARK_NORMAL) && (itask != ARK_ONE_STEP))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_BAD_ITASK);
    return (ARK_ILL_INPUT);
  }

  /* start profiler */
  SUNDIALS_MARK_FUNCTION_BEGIN(ARK_PROFILER);

  /* perform first-step-specific initializations:
     - initialize tret values to initialization time
     - perform initial integrator setup  */
  if (ark_mem->initsetup)
  {
    ark_mem->tretlast = *tret = ark_mem->tcur;
    retval                    = arkInitialSetup(ark_mem, tout);
    if (retval != ARK_SUCCESS)
    {
      SUNDIALS_MARK_FUNCTION_END(ARK_PROFILER);
      return (retval);
    }
  }

  /* perform stopping tests */
  if (!ark_mem->initsetup)
  {
    if (arkStopTests(ark_mem, tout, yout, tret, itask, &retval))
    {
      SUNDIALS_MARK_FUNCTION_END(ARK_PROFILER);
      return (retval);
    }
  }

  /*--------------------------------------------------
    Looping point for successful internal steps

    - update the ewt/rwt vectors for upcoming step
    - check for errors (too many steps, too much
      accuracy requested, step size too small)
    - loop over attempts at a new step:
      * try to take step (via time stepper module),
        handle solver convergence or other failures
      * if the stepper requests ARK_RETRY_STEP, we
        retry the step without accumulating failures.
        A stepper should never request this multiple
        times in a row.
      * perform constraint-handling (if selected)
      * check temporal error
      * if all of the above pass, complete step by
        updating current time, solution, error &
        stepsize history arrays.
    - perform stop tests:
      * check for root in last step taken
      * check if tout was passed
      * check if close to tstop
      * check if in ONE_STEP mode (must return)
    --------------------------------------------------*/
  nstloc = 0;
  for (;;)
  {
    ark_mem->next_h = ark_mem->h;

    /* Reset and check ewt and rwt */
    if (!ark_mem->initsetup)
    {
      ewtsetOK = ark_mem->efun(ark_mem->yn, ark_mem->ewt, ark_mem->e_data);
      if (ewtsetOK != 0)
      {
        if (ark_mem->itol == ARK_WF)
        {
          arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                          MSG_ARK_EWT_NOW_FAIL, ark_mem->tcur);
        }
        else
        {
          arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                          MSG_ARK_EWT_NOW_BAD, ark_mem->tcur);
        }

        istate            = ARK_ILL_INPUT;
        ark_mem->tretlast = *tret = ark_mem->tcur;
        N_VScale(ONE, ark_mem->yn, yout);
        break;
      }

      if (!ark_mem->rwt_is_ewt)
      {
        ewtsetOK = ark_mem->rfun(ark_mem->yn, ark_mem->rwt, ark_mem->r_data);
        if (ewtsetOK != 0)
        {
          if (ark_mem->itol == ARK_WF)
          {
            arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                            __FILE__, MSG_ARK_RWT_NOW_FAIL, ark_mem->tcur);
          }
          else
          {
            arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                            __FILE__, MSG_ARK_RWT_NOW_BAD, ark_mem->tcur);
          }

          istate            = ARK_ILL_INPUT;
          ark_mem->tretlast = *tret = ark_mem->tcur;
          N_VScale(ONE, ark_mem->yn, yout);
          break;
        }
      }
    }

    /* Check for too many steps */
    if ((ark_mem->mxstep > 0) && (nstloc >= ark_mem->mxstep))
    {
      arkProcessError(ark_mem, ARK_TOO_MUCH_WORK, __LINE__, __func__, __FILE__,
                      MSG_ARK_MAX_STEPS, ark_mem->tcur);
      istate            = ARK_TOO_MUCH_WORK;
      ark_mem->tretlast = *tret = ark_mem->tcur;
      N_VScale(ONE, ark_mem->yn, yout);
      break;
    }

    /* Check for too much accuracy requested */
    nrm            = N_VWrmsNorm(ark_mem->yn, ark_mem->ewt);
    ark_mem->tolsf = ark_mem->uround * nrm;
    if (ark_mem->tolsf > ONE && !ark_mem->fixedstep)
    {
      arkProcessError(ark_mem, ARK_TOO_MUCH_ACC, __LINE__, __func__, __FILE__,
                      MSG_ARK_TOO_MUCH_ACC, ark_mem->tcur);
      istate            = ARK_TOO_MUCH_ACC;
      ark_mem->tretlast = *tret = ark_mem->tcur;
      N_VScale(ONE, ark_mem->yn, yout);
      ark_mem->tolsf *= TWO;
      break;
    }
    else { ark_mem->tolsf = ONE; }

    /* Check for h below roundoff level in tn */
    if (ark_mem->tcur + ark_mem->h == ark_mem->tcur)
    {
      ark_mem->nhnil++;
      if (ark_mem->nhnil <= ark_mem->mxhnil)
      {
        arkProcessError(ark_mem, ARK_WARNING, __LINE__, __func__, __FILE__,
                        MSG_ARK_HNIL, ark_mem->tcur, ark_mem->h);
      }
      if (ark_mem->nhnil == ark_mem->mxhnil)
      {
        arkProcessError(ark_mem, ARK_WARNING, __LINE__, __func__, __FILE__,
                        MSG_ARK_HNIL_DONE);
      }
    }

    /* Update parameter for upcoming step size */
    if (ark_mem->hprime != ark_mem->h)
    {
      ark_mem->h      = ark_mem->h * ark_mem->eta;
      ark_mem->next_h = ark_mem->h;
    }
    if (ark_mem->fixedstep)
    {
      ark_mem->h      = ark_mem->hin;
      ark_mem->next_h = ark_mem->h;

      /* patch for 'fixedstep' + 'tstop' use case:
         limit fixed step size if step would overtake tstop */
      if (ark_mem->tstopset)
      {
        if ((ark_mem->tcur + ark_mem->h - ark_mem->tstop) * ark_mem->h > ZERO)
        {
          ark_mem->h = (ark_mem->tstop - ark_mem->tcur) *
                       (ONE - FOUR * ark_mem->uround);
        }
      }
    }

    /* Looping point for step attempts */
    dsm      = ZERO;
    kflag    = ARK_SUCCESS;
    attempts = ncf = nef = constrfails = ark_mem->last_kflag = 0;
    relax_fails                                              = 0;
    nflag                                                    = FIRST_CALL;
    for (;;)
    {
      /* increment attempt counters
         Note: kflag can only equal ARK_RETRY_STEP if the stepper rejected
         the current step size before performing calculations. Thus, we do
         not include those when keeping track of step "attempts". */
      if (kflag != ARK_RETRY_STEP)
      {
        attempts++;
        ark_mem->nst_attempts++;
      }

      SUNLogInfo(ARK_LOGGER, "begin-step-attempt",
                 "step = %li, tn = " SUN_FORMAT_G ", h = " SUN_FORMAT_G,
                 ark_mem->nst + 1, ark_mem->tn, ark_mem->h);

      /* Call time stepper module to attempt a step:
            0 => step completed successfully
           >0 => step encountered recoverable failure; reduce step if possible
           <0 => step encountered unrecoverable failure */
      kflag = ark_mem->step((void*)ark_mem, &dsm, &nflag);
      if (kflag < 0)
      {
        /* Log fatal errors here, other returns handled below */
        SUNLogInfo(ARK_LOGGER, "end-step-attempt",
                   "status = failed step, kflag = %i", kflag);
        break;
      }

      /* handle solver convergence failures */
      kflag = arkCheckConvergence(ark_mem, &nflag, &ncf);

      SUNLogInfoIf(kflag != ARK_SUCCESS, ARK_LOGGER, "end-step-attempt",
                   "status = failed step, kflag = %i", kflag);

      if (kflag < 0) { break; }

      /* Perform relaxation:
           - computes relaxation parameter
           - on success, updates ycur, h, and dsm
           - on recoverable failure, updates eta and signals to retry step
           - on fatal error, returns negative error flag */
      if (ark_mem->relax_enabled && (kflag == ARK_SUCCESS))
      {
        kflag = arkRelax(ark_mem, &relax_fails, &dsm);

        SUNLogInfoIf(kflag != ARK_SUCCESS, ARK_LOGGER, "end-step-attempt",
                     "status = failed relaxtion, kflag = %i", kflag);

        if (kflag < 0) { break; }
      }

      /* perform constraint-handling (if selected, and if solver check passed) */
      if (ark_mem->constraints && (kflag == ARK_SUCCESS))
      {
        kflag = arkCheckConstraints(ark_mem, &constrfails, &nflag);

        SUNLogInfoIf(kflag != ARK_SUCCESS, ARK_LOGGER, "end-step-attempt",
                     "status = failed constraints, kflag = %i", kflag);

        if (kflag < 0) { break; }
      }

      /* when fixed time-stepping is enabled, 'success' == successful stage solves
         (checked in previous block), so just enforce no step size change */
      if (ark_mem->fixedstep)
      {
        ark_mem->eta = ONE;
        SUNLogInfo(ARK_LOGGER, "end-step-attempt", "status = success");
        break;
      }

      /* check temporal error (if checks above passed) */
      if (kflag == ARK_SUCCESS)
      {
        kflag = arkCheckTemporalError(ark_mem, &nflag, &nef, dsm);

        SUNLogInfoIf(kflag != ARK_SUCCESS, ARK_LOGGER, "end-step-attempt",
                     "status = failed error test, dsm = " SUN_FORMAT_G
                     ", kflag = %i",
                     dsm, kflag);

        if (kflag < 0) { break; }
      }

      /* if ignoring temporal error test result (XBraid) force step to pass */
      if (ark_mem->force_pass)
      {
        ark_mem->last_kflag = kflag;
        kflag               = ARK_SUCCESS;
        SUNLogInfo(ARK_LOGGER, "end-step-attempt", "status = success");
        break;
      }

      /* break attempt loop on successful step */
      if (kflag == ARK_SUCCESS)
      {
        SUNLogInfo(ARK_LOGGER, "end-step-attempt",
                   "status = success, dsm = " SUN_FORMAT_G, dsm);
        break;
      }

      /* unsuccessful step, if |h| = hmin, return ARK_ERR_FAILURE */
      if (SUNRabs(ark_mem->h) <= ark_mem->hmin * ONEPSM)
      {
        SUNDIALS_MARK_FUNCTION_END(ARK_PROFILER);
        return (ARK_ERR_FAILURE);
      }

      /* update h, hprime and next_h for next iteration */
      ark_mem->h *= ark_mem->eta;
      ark_mem->next_h = ark_mem->hprime = ark_mem->h;

    } /* end looping for step attempts */

    /* If step attempt loop succeeded, complete step (update current time, solution,
       error stepsize history arrays; call user-supplied step postprocessing function)
       (added stuff from arkStep_PrepareNextStep -- revisit) */
    if (kflag == ARK_SUCCESS) { kflag = arkCompleteStep(ark_mem, dsm); }

    /* If step attempt loop failed, process flag and return to user */
    if (kflag != ARK_SUCCESS)
    {
      istate            = arkHandleFailure(ark_mem, kflag);
      ark_mem->tretlast = *tret = ark_mem->tcur;
      N_VScale(ONE, ark_mem->yn, yout);
      break;
    }

    nstloc++;

    /* Check for root in last step taken. */
    if (ark_mem->root_mem != NULL)
    {
      if (ark_mem->root_mem->nrtfn > 0)
      {
        retval = arkRootCheck3((void*)ark_mem, tout, itask);
        if (retval == RTFOUND)
        { /* A new root was found */
          ark_mem->root_mem->irfnd = 1;
          istate                   = ARK_ROOT_RETURN;
          ark_mem->tretlast = *tret = ark_mem->root_mem->tlo;
          break;
        }
        else if (retval == ARK_RTFUNC_FAIL)
        { /* g failed */
          arkProcessError(ark_mem, ARK_RTFUNC_FAIL, __LINE__, __func__, __FILE__,
                          MSG_ARK_RTFUNC_FAILED, ark_mem->root_mem->tlo);
          istate = ARK_RTFUNC_FAIL;
          break;
        }

        /* If we are at the end of the first step and we still have
           some event functions that are inactive, issue a warning
           as this may indicate a user error in the implementation
           of the root function. */
        if (ark_mem->nst == 1)
        {
          inactive_roots = SUNFALSE;
          for (ir = 0; ir < ark_mem->root_mem->nrtfn; ir++)
          {
            if (!ark_mem->root_mem->gactive[ir])
            {
              inactive_roots = SUNTRUE;
              break;
            }
          }
          if ((ark_mem->root_mem->mxgnull > 0) && inactive_roots)
          {
            arkProcessError(ark_mem, ARK_WARNING, __LINE__, __func__, __FILE__,
                            MSG_ARK_INACTIVE_ROOTS);
          }
        }
      }
    }

    /* Check if tn is at tstop or near tstop */
    if (ark_mem->tstopset)
    {
      troundoff = FUZZ_FACTOR * ark_mem->uround *
                  (SUNRabs(ark_mem->tcur) + SUNRabs(ark_mem->h));

      if (SUNRabs(ark_mem->tcur - ark_mem->tstop) <= troundoff)
      {
        /* Ensure tout >= tstop, otherwise check for tout return below */
        if ((tout - ark_mem->tstop) * ark_mem->h >= ZERO ||
            SUNRabs(tout - ark_mem->tstop) <= troundoff)
        {
          if (ark_mem->tstopinterp && ark_mem->interp)
          {
            retval = ARKodeGetDky(ark_mem, ark_mem->tstop, 0, yout);
            if (retval != ARK_SUCCESS)
            {
              arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                              MSG_ARK_INTERPOLATION_FAIL, ark_mem->tstop);
              istate = retval;
              break;
            }
          }
          else { N_VScale(ONE, ark_mem->yn, yout); }
          ark_mem->tretlast = *tret = ark_mem->tstop;
          ark_mem->tstopset         = SUNFALSE;
          istate                    = ARK_TSTOP_RETURN;
          break;
        }
      }
      /* limit upcoming step if it will overcome tstop */
      else if ((ark_mem->tcur + ark_mem->hprime - ark_mem->tstop) * ark_mem->h >
               ZERO)
      {
        ark_mem->hprime = (ark_mem->tstop - ark_mem->tcur) *
                          (ONE - FOUR * ark_mem->uround);
        ark_mem->eta = ark_mem->hprime / ark_mem->h;
      }
    }

    /* In NORMAL mode, check if tout reached */
    if ((itask == ARK_NORMAL) && (ark_mem->tcur - tout) * ark_mem->h >= ZERO)
    {
      if (ark_mem->interp)
      {
        retval = ARKodeGetDky(ark_mem, tout, 0, yout);
        if (retval != ARK_SUCCESS)
        {
          arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                          MSG_ARK_INTERPOLATION_FAIL, tout);
          istate = retval;
          break;
        }
        ark_mem->tretlast = *tret = tout;
      }
      else
      {
        N_VScale(ONE, ark_mem->yn, yout);
        ark_mem->tretlast = *tret = ark_mem->tcur;
      }
      ark_mem->next_h = ark_mem->hprime;
      istate          = ARK_SUCCESS;
      break;
    }

    /* In ONE_STEP mode, copy y and exit loop */
    if (itask == ARK_ONE_STEP)
    {
      istate            = ARK_SUCCESS;
      ark_mem->tretlast = *tret = ark_mem->tcur;
      N_VScale(ONE, ark_mem->yn, yout);
      ark_mem->next_h = ark_mem->hprime;
      break;
    }

  } /* end looping for internal steps */

  /* stop profiler and return */
  SUNDIALS_MARK_FUNCTION_END(ARK_PROFILER);
  return (istate);
}

/*---------------------------------------------------------------
  ARKodeGetDky:

  This routine computes the k-th derivative of the interpolating
  polynomial at the time t and stores the result in the vector
  dky. This routine internally calls arkInterpEvaluate to perform
  the interpolation.  We have the restriction that 0 <= k <= 3.
  This routine uses an interpolating polynomial of degree
  max(deg, k), i.e. it will form a polynomial of the degree
  available by the interpolation module and/or requested by
  the user through deg, unless higher-order derivatives are
  requested.

  This function is called by ARKodeEvolve with k=0 and t=tout to
  perform interpolation of outputs, but may also be called
  indirectly by the user via time step module *StepGetDky calls.
  Note: in all cases it will be called after ark_tcur has been
  updated to correspond with the end time of the last successful
  step.
  ---------------------------------------------------------------*/
int ARKodeGetDky(void* arkode_mem, sunrealtype t, int k, N_Vector dky)
{
  sunrealtype s, tfuzz, tp, tn1;
  int retval;
  ARKodeMem ark_mem;

  /* Check if ark_mem exists */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Check all inputs for legality */
  if (dky == NULL)
  {
    arkProcessError(ark_mem, ARK_BAD_DKY, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_DKY);
    return (ARK_BAD_DKY);
  }
  if (ark_mem->interp == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "Missing interpolation structure");
    return (ARK_MEM_NULL);
  }

  /* Allow for some slack */
  tfuzz = FUZZ_FACTOR * ark_mem->uround *
          (SUNRabs(ark_mem->tcur) + SUNRabs(ark_mem->hold));
  if (ark_mem->hold < ZERO) { tfuzz = -tfuzz; }
  tp  = ark_mem->tcur - ark_mem->hold - tfuzz;
  tn1 = ark_mem->tcur + tfuzz;
  if ((t - tp) * (t - tn1) > ZERO)
  {
    arkProcessError(ark_mem, ARK_BAD_T, __LINE__, __func__, __FILE__,
                    MSG_ARK_BAD_T, t, ark_mem->tcur - ark_mem->hold,
                    ark_mem->tcur);
    return (ARK_BAD_T);
  }

  /* call arkInterpEvaluate to evaluate result */
  s      = (t - ark_mem->tcur) / ark_mem->h;
  retval = arkInterpEvaluate(ark_mem, ark_mem->interp, s, k,
                             ARK_INTERP_MAX_DEGREE, dky);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Error calling arkInterpEvaluate");
    return (retval);
  }
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  ARKodeFree:

  This routine frees the ARKODE infrastructure memory.
  ---------------------------------------------------------------*/
void ARKodeFree(void** arkode_mem)
{
  ARKodeMem ark_mem;

  if (*arkode_mem == NULL) { return; }

  ark_mem = (ARKodeMem)(*arkode_mem);

  /* free the time-stepper module memory (if provided) */
  if (ark_mem->step_free) { ark_mem->step_free(ark_mem); }

  /* free vector storage */
  arkFreeVectors(ark_mem);

  /* free the time step adaptivity module */
  if (ark_mem->hadapt_mem != NULL)
  {
    if (ark_mem->hadapt_mem->owncontroller)
    {
      (void)SUNAdaptController_Destroy(ark_mem->hadapt_mem->hcontroller);
      ark_mem->hadapt_mem->owncontroller = SUNFALSE;
    }
    free(ark_mem->hadapt_mem);
    ark_mem->hadapt_mem = NULL;
  }

  /* free the interpolation module */
  if (ark_mem->interp != NULL)
  {
    arkInterpFree(ark_mem, ark_mem->interp);
    ark_mem->interp = NULL;
  }

  /* free the root-finding module */
  if (ark_mem->root_mem != NULL)
  {
    (void)arkRootFree(ark_mem);
    ark_mem->root_mem = NULL;
  }

  /* free the relaxation module */
  if (ark_mem->relax_mem)
  {
    (void)arkRelaxDestroy(ark_mem->relax_mem);
    ark_mem->relax_mem = NULL;
  }

#if defined(SUNDIALS_ENABLE_PYTHON)
  arkode_user_supplied_fn_table_destroy(ark_mem->python);
#endif
  ark_mem->python = NULL;

  free(*arkode_mem);
  *arkode_mem = NULL;
}

/*---------------------------------------------------------------
  ARKodePrintMem:

  This routine outputs the ark_mem structure to a specified file
  pointer.
  ---------------------------------------------------------------*/
void ARKodePrintMem(void* arkode_mem, FILE* outfile)
{
  ARKodeMem ark_mem;

  /* Check if ark_mem exists */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return;
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* if outfile==NULL, set it to stdout */
  if (outfile == NULL) { outfile = stdout; }

  /* output general values */
  fprintf(outfile, "itol = %i\n", ark_mem->itol);
  fprintf(outfile, "ritol = %i\n", ark_mem->ritol);
  fprintf(outfile, "mxhnil = %i\n", ark_mem->mxhnil);
  fprintf(outfile, "mxstep = %li\n", ark_mem->mxstep);
  fprintf(outfile, "lrw1 = %li\n", (long int)ark_mem->lrw1);
  fprintf(outfile, "liw1 = %li\n", (long int)ark_mem->liw1);
  fprintf(outfile, "lrw = %li\n", (long int)ark_mem->lrw);
  fprintf(outfile, "liw = %li\n", (long int)ark_mem->liw);
  fprintf(outfile, "user_efun = %i\n", ark_mem->user_efun);
  fprintf(outfile, "tstopset = %i\n", ark_mem->tstopset);
  fprintf(outfile, "tstopinterp = %i\n", ark_mem->tstopinterp);
  fprintf(outfile, "tstop = " SUN_FORMAT_G "\n", ark_mem->tstop);
  fprintf(outfile, "VabstolMallocDone = %i\n", ark_mem->VabstolMallocDone);
  fprintf(outfile, "MallocDone = %i\n", ark_mem->MallocDone);
  fprintf(outfile, "initsetup = %i\n", ark_mem->initsetup);
  fprintf(outfile, "init_type = %i\n", ark_mem->init_type);
  fprintf(outfile, "firststage = %i\n", ark_mem->firststage);
  fprintf(outfile, "uround = " SUN_FORMAT_G "\n", ark_mem->uround);
  fprintf(outfile, "reltol = " SUN_FORMAT_G "\n", ark_mem->reltol);
  fprintf(outfile, "Sabstol = " SUN_FORMAT_G "\n", ark_mem->Sabstol);
  fprintf(outfile, "fixedstep = %i\n", ark_mem->fixedstep);
  fprintf(outfile, "tolsf = " SUN_FORMAT_G "\n", ark_mem->tolsf);
  fprintf(outfile, "call_fullrhs = %i\n", ark_mem->call_fullrhs);
  fprintf(outfile, "do_adjoint = %i\n", ark_mem->do_adjoint);

  /* output counters */
  fprintf(outfile, "nhnil = %i\n", ark_mem->nhnil);
  fprintf(outfile, "nst_attempts = %li\n", ark_mem->nst_attempts);
  fprintf(outfile, "nst = %li\n", ark_mem->nst);
  fprintf(outfile, "ncfn = %li\n", ark_mem->ncfn);
  fprintf(outfile, "netf = %li\n", ark_mem->netf);

  /* output time-stepping values */
  fprintf(outfile, "hin = " SUN_FORMAT_G "\n", ark_mem->hin);
  fprintf(outfile, "h = " SUN_FORMAT_G "\n", ark_mem->h);
  fprintf(outfile, "hprime = " SUN_FORMAT_G "\n", ark_mem->hprime);
  fprintf(outfile, "next_h = " SUN_FORMAT_G "\n", ark_mem->next_h);
  fprintf(outfile, "eta = " SUN_FORMAT_G "\n", ark_mem->eta);
  fprintf(outfile, "tcur = " SUN_FORMAT_G "\n", ark_mem->tcur);
  fprintf(outfile, "tretlast = " SUN_FORMAT_G "\n", ark_mem->tretlast);
  fprintf(outfile, "hmin = " SUN_FORMAT_G "\n", ark_mem->hmin);
  fprintf(outfile, "hmax_inv = " SUN_FORMAT_G "\n", ark_mem->hmax_inv);
  fprintf(outfile, "h0u = " SUN_FORMAT_G "\n", ark_mem->h0u);
  fprintf(outfile, "tn = " SUN_FORMAT_G "\n", ark_mem->tn);
  fprintf(outfile, "hold = " SUN_FORMAT_G "\n", ark_mem->hold);
  fprintf(outfile, "maxnef = %i\n", ark_mem->maxnef);
  fprintf(outfile, "maxncf = %i\n", ark_mem->maxncf);

  /* output time-stepping adaptivity structure */
  fprintf(outfile, "timestep adaptivity structure:\n");
  arkPrintAdaptMem(ark_mem->hadapt_mem, outfile);

  /* output inequality constraints quantities */
  fprintf(outfile, "maxconstrfails = %i\n", ark_mem->maxconstrfails);

  /* output root-finding quantities */
  if (ark_mem->root_mem != NULL)
  {
    (void)arkPrintRootMem((void*)ark_mem, outfile);
  }

  /* output interpolation quantities */
  if (ark_mem->interp) { arkInterpPrintMem(ark_mem->interp, outfile); }
  else { fprintf(outfile, "interpolation = NULL\n"); }

#ifdef SUNDIALS_DEBUG_PRINTVEC
  /* output vector quantities */
  fprintf(outfile, "Vapbsol:\n");
  N_VPrintFile(ark_mem->Vabstol, outfile);
  fprintf(outfile, "ewt:\n");
  N_VPrintFile(ark_mem->ewt, outfile);
  if (!ark_mem->rwt_is_ewt)
  {
    fprintf(outfile, "rwt:\n");
    N_VPrintFile(ark_mem->rwt, outfile);
  }
  fprintf(outfile, "ycur:\n");
  N_VPrintFile(ark_mem->ycur, outfile);
  fprintf(outfile, "yn:\n");
  N_VPrintFile(ark_mem->yn, outfile);
  fprintf(outfile, "fn:\n");
  if (ark_mem->fn) { N_VPrintFile(ark_mem->fn, outfile); }
  fprintf(outfile, "tempv1:\n");
  N_VPrintFile(ark_mem->tempv1, outfile);
  fprintf(outfile, "tempv2:\n");
  N_VPrintFile(ark_mem->tempv2, outfile);
  fprintf(outfile, "tempv3:\n");
  N_VPrintFile(ark_mem->tempv3, outfile);
  fprintf(outfile, "tempv4:\n");
  N_VPrintFile(ark_mem->tempv4, outfile);
  fprintf(outfile, "tempv5:\n");
  N_VPrintFile(ark_mem->tempv5, outfile);
  fprintf(outfile, "constraints:\n");
  N_VPrintFile(ark_mem->constraints, outfile);
#endif

  /* Call stepper PrintMem function (if provided) */
  if (ark_mem->step_printmem) { ark_mem->step_printmem(ark_mem, outfile); }
}

/*------------------------------------------------------------------------------
  ARKodeCreateMRIStepInnerStepper

  Wraps an ARKODE integrator as an MRIStep inner stepper.
  ----------------------------------------------------------------------------*/

int ARKodeCreateMRIStepInnerStepper(void* inner_arkode_mem,
                                    MRIStepInnerStepper* stepper)
{
  ARKodeMem ark_mem;
  int retval;

  /* Check if ark_mem exists */
  if (inner_arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)inner_arkode_mem;

  /* return with an error if the ARKODE solver does not support forcing */
  if (ark_mem->step_setforcing == NULL)
  {
    arkProcessError(ark_mem, ARK_STEPPER_UNSUPPORTED, __LINE__, __func__,
                    __FILE__, "time-stepping module does not support forcing");
    return (ARK_STEPPER_UNSUPPORTED);
  }

  retval = MRIStepInnerStepper_Create(ark_mem->sunctx, stepper);
  if (retval != ARK_SUCCESS) { return (retval); }

  retval = MRIStepInnerStepper_SetContent(*stepper, inner_arkode_mem);
  if (retval != ARK_SUCCESS) { return (retval); }

  retval = MRIStepInnerStepper_SetEvolveFn(*stepper, ark_MRIStepInnerEvolve);
  if (retval != ARK_SUCCESS) { return (retval); }

  retval = MRIStepInnerStepper_SetFullRhsFn(*stepper, ark_MRIStepInnerFullRhs);
  if (retval != ARK_SUCCESS) { return (retval); }

  retval = MRIStepInnerStepper_SetResetFn(*stepper, ark_MRIStepInnerReset);
  if (retval != ARK_SUCCESS) { return (retval); }

  retval =
    MRIStepInnerStepper_SetAccumulatedErrorGetFn(*stepper,
                                                 ark_MRIStepInnerGetAccumulatedError);
  if (retval != ARK_SUCCESS) { return (retval); }

  retval =
    MRIStepInnerStepper_SetAccumulatedErrorResetFn(*stepper,
                                                   ark_MRIStepInnerResetAccumulatedError);
  if (retval != ARK_SUCCESS) { return (retval); }

  retval = MRIStepInnerStepper_SetRTolFn(*stepper, ark_MRIStepInnerSetRTol);
  if (retval != ARK_SUCCESS) { return (retval); }

  return (ARK_SUCCESS);
}

/*===============================================================
  Private internal functions
  ===============================================================*/

/*---------------------------------------------------------------
  arkCreate:

  arkCreate creates an internal memory block for a problem to
  be solved by a time step module built on ARKODE.  If successful,
  arkCreate returns a pointer to the problem memory. If an
  initialization error occurs, arkCreate prints an error message
  to standard err and returns NULL.
  ---------------------------------------------------------------*/
ARKodeMem arkCreate(SUNContext sunctx)
{
  int iret;
  ARKodeMem ark_mem;

  if (!sunctx)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_SUNCTX);
    return (NULL);
  }

  ark_mem = NULL;
  ark_mem = (ARKodeMem)malloc(sizeof(struct ARKodeMemRec));
  if (ark_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_ARKMEM_FAIL);
    return (NULL);
  }

  /* Zero out ark_mem */
  memset(ark_mem, 0, sizeof(struct ARKodeMemRec));

  /* Set the context */
  ark_mem->sunctx = sunctx;

  /* Set the Python context to NULL */
  ark_mem->python = NULL;

  /* Set uround */
  ark_mem->uround = SUN_UNIT_ROUNDOFF;

  /* Initialize time step module to NULL */
  ark_mem->step_attachlinsol              = NULL;
  ark_mem->step_attachmasssol             = NULL;
  ark_mem->step_disablelsetup             = NULL;
  ark_mem->step_disablemsetup             = NULL;
  ark_mem->step_getlinmem                 = NULL;
  ark_mem->step_getmassmem                = NULL;
  ark_mem->step_getimplicitrhs            = NULL;
  ark_mem->step_mmult                     = NULL;
  ark_mem->step_getgammas                 = NULL;
  ark_mem->step_init                      = NULL;
  ark_mem->step_fullrhs                   = NULL;
  ark_mem->step                           = NULL;
  ark_mem->step_setuserdata               = NULL;
  ark_mem->step_printallstats             = NULL;
  ark_mem->step_writeparameters           = NULL;
  ark_mem->step_resize                    = NULL;
  ark_mem->step_reset                     = NULL;
  ark_mem->step_free                      = NULL;
  ark_mem->step_printmem                  = NULL;
  ark_mem->step_setdefaults               = NULL;
  ark_mem->step_computestate              = NULL;
  ark_mem->step_setrelaxfn                = NULL;
  ark_mem->step_setorder                  = NULL;
  ark_mem->step_setnonlinearsolver        = NULL;
  ark_mem->step_setlinear                 = NULL;
  ark_mem->step_setnonlinear              = NULL;
  ark_mem->step_setautonomous             = NULL;
  ark_mem->step_setnlsrhsfn               = NULL;
  ark_mem->step_setdeduceimplicitrhs      = NULL;
  ark_mem->step_setnonlincrdown           = NULL;
  ark_mem->step_setnonlinrdiv             = NULL;
  ark_mem->step_setdeltagammamax          = NULL;
  ark_mem->step_setlsetupfrequency        = NULL;
  ark_mem->step_setpredictormethod        = NULL;
  ark_mem->step_setmaxnonliniters         = NULL;
  ark_mem->step_setnonlinconvcoef         = NULL;
  ark_mem->step_setstagepredictfn         = NULL;
  ark_mem->step_getnumrhsevals            = NULL;
  ark_mem->step_setstepdirection          = NULL;
  ark_mem->step_setoptions                = NULL;
  ark_mem->step_getnumlinsolvsetups       = NULL;
  ark_mem->step_setadaptcontroller        = NULL;
  ark_mem->step_getestlocalerrors         = NULL;
  ark_mem->step_getcurrentgamma           = NULL;
  ark_mem->step_getnonlinearsystemdata    = NULL;
  ark_mem->step_getnumnonlinsolviters     = NULL;
  ark_mem->step_getnumnonlinsolvconvfails = NULL;
  ark_mem->step_getnonlinsolvstats        = NULL;
  ark_mem->step_setforcing                = NULL;
  ark_mem->step_mem                       = NULL;
  ark_mem->step_supports_adaptive         = SUNFALSE;
  ark_mem->step_supports_implicit         = SUNFALSE;
  ark_mem->step_supports_massmatrix       = SUNFALSE;
  ark_mem->step_supports_relaxation       = SUNFALSE;

  /* Initialize root finding variables */
  ark_mem->root_mem = NULL;

  /* Initialize inequality constraints variables */
  ark_mem->constraints = NULL;

  /* Initialize relaxation variables */
  ark_mem->relax_enabled = SUNFALSE;
  ark_mem->relax_mem     = NULL;

  /* Initialize lrw and liw */
  ark_mem->lrw = 18;
  ark_mem->liw = 53; /* fcn/data ptr, int, long int, sunindextype, sunbooleantype */

  /* No mallocs have been done yet */
  ark_mem->VabstolMallocDone  = SUNFALSE;
  ark_mem->VRabstolMallocDone = SUNFALSE;
  ark_mem->MallocDone         = SUNFALSE;

  /* No user-supplied step postprocessing function yet */
  ark_mem->ProcessStep = NULL;
  ark_mem->ps_data     = NULL;

  /* No user-supplied stage postprocessing function yet */
  ark_mem->ProcessStage = NULL;

  /* No user_data pointer yet */
  ark_mem->user_data = NULL;

  /* Allocate step adaptivity structure and note storage */
  ark_mem->hadapt_mem = arkAdaptInit();
  if (ark_mem->hadapt_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                    "Allocation of step adaptivity structure failed");
    ARKodeFree((void**)&ark_mem);
    return (NULL);
  }
  ark_mem->lrw += ARK_ADAPT_LRW;
  ark_mem->liw += ARK_ADAPT_LIW;

  /* Initialize the interpolation structure to NULL */
  ark_mem->interp        = NULL;
  ark_mem->interp_type   = ARK_INTERP_HERMITE;
  ark_mem->interp_degree = ARK_INTERP_MAX_DEGREE;

  /* Initially, rwt should point to ewt */
  ark_mem->rwt_is_ewt = SUNTRUE;

  /* Indicate that calling the full RHS function is not required, this flag is
     updated to SUNTRUE by the interpolation module initialization function
     and/or the stepper initialization function in arkInitialSetup */
  ark_mem->call_fullrhs = SUNFALSE;

  /* Indicate that the problem needs to be initialized */
  ark_mem->initsetup   = SUNTRUE;
  ark_mem->init_type   = FIRST_INIT;
  ark_mem->firststage  = SUNTRUE;
  ark_mem->initialized = SUNFALSE;

  /* Initial step size has not been determined yet */
  ark_mem->h   = ZERO;
  ark_mem->h0u = ZERO;

  /* Accumulated error estimation strategy */
  ark_mem->AccumErrorType = ARK_ACCUMERROR_NONE;
  ark_mem->AccumError     = ZERO;

  /* Set default values for integrator and stepper optional inputs */
  iret = ARKodeSetDefaults(ark_mem);
  if (iret != ARK_SUCCESS)
  {
    arkProcessError(NULL, 0, __LINE__, __func__, __FILE__,
                    "Error setting default solver options");
    ARKodeFree((void**)&ark_mem);
    return (NULL);
  }

  ark_mem->load_checkpoint_fail = SUNFALSE;
  ark_mem->do_adjoint           = SUNFALSE;

  /* Return pointer to ARKODE memory block */
  return (ark_mem);
}

/*---------------------------------------------------------------
  arkRwtSet

  This routine is responsible for setting the residual weight
  vector rwt, according to tol_type, as follows:

  (1) rwt[i] = 1 / (reltol * SUNRabs(M*ycur[i]) + rabstol), i=0,...,neq-1
      if tol_type = ARK_SS
  (2) rwt[i] = 1 / (reltol * SUNRabs(M*ycur[i]) + rabstol[i]), i=0,...,neq-1
      if tol_type = ARK_SV
  (3) unset if tol_type is any other value (occurs rwt=ewt)

  arkRwtSet returns 0 if rwt is successfully set as above to a
  positive vector and -1 otherwise. In the latter case, rwt is
  considered undefined.

  All the real work is done in the routines arkRwtSetSS, arkRwtSetSV.
  ---------------------------------------------------------------*/
int arkRwtSet(N_Vector y, N_Vector weight, void* data)
{
  ARKodeMem ark_mem;
  N_Vector My;
  int flag = 0;

  /* data points to ark_mem here */
  ark_mem = (ARKodeMem)data;

  /* return if rwt is just ewt */
  if (ark_mem->rwt_is_ewt) { return (0); }

  /* put M*y into ark_tempv1 */
  My = ark_mem->tempv1;
  if (ark_mem->step_mmult != NULL)
  {
    flag = ark_mem->step_mmult((void*)ark_mem, y, My);
    if (flag != ARK_SUCCESS) { return (ARK_MASSMULT_FAIL); }
  }
  else
  { /* this condition should not apply, but just in case */
    N_VScale(ONE, y, My);
  }

  /* call appropriate routine to fill rwt */
  switch (ark_mem->ritol)
  {
  case ARK_SS: flag = arkRwtSetSS(ark_mem, My, weight); break;
  case ARK_SV: flag = arkRwtSetSV(ark_mem, My, weight); break;
  }

  return (flag);
}

/*---------------------------------------------------------------
  arkInit:

  arkInit allocates and initializes memory for a problem. All
  inputs are checked for errors. If any error occurs during
  initialization, an error flag is returned. Otherwise, it returns
  ARK_SUCCESS.  This routine should be called by an ARKODE
  timestepper module (not by the user).  This routine must be
  called prior to calling ARKodeEvolve to evolve the problem. The
  initialization type indicates if the values of internal counters
  should be reinitialized (FIRST_INIT) or retained (RESET_INIT).
  ---------------------------------------------------------------*/
int arkInit(ARKodeMem ark_mem, sunrealtype t0, N_Vector y0, int init_type)
{
  sunbooleantype stepperOK, nvectorOK, allocOK;
  int retval;
  sunindextype lrw1, liw1;

  /* Check ark_mem */
  if (ark_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }

  /* Check for legal input parameters */
  if (y0 == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_NULL_Y0);
    return (ARK_ILL_INPUT);
  }

  /* Check if reset was called before the first Evolve call */
  if (init_type == RESET_INIT && !(ark_mem->initialized))
  {
    init_type = FIRST_INIT;
  }

  /* Check if allocations have been done i.e., is this first init call */
  if (ark_mem->MallocDone == SUNFALSE)
  {
    /* Test if all required time stepper operations are implemented */
    stepperOK = arkCheckTimestepper(ark_mem);
    if (!stepperOK)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Time stepper module is missing required functionality");
      return (ARK_ILL_INPUT);
    }

    /* Test if all required vector operations are implemented */
    nvectorOK = arkCheckNvectorRequired(y0);
    if (!nvectorOK)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      MSG_ARK_BAD_NVECTOR);
      return (ARK_ILL_INPUT);
    }

    /* Set space requirements for one N_Vector */
    if (y0->ops->nvspace != NULL) { N_VSpace(y0, &lrw1, &liw1); }
    else
    {
      lrw1 = 0;
      liw1 = 0;
    }
    ark_mem->lrw1 = lrw1;
    ark_mem->liw1 = liw1;

    /* Allocate the solver vectors (using y0 as a template) */
    allocOK = arkAllocVectors(ark_mem, y0);
    if (!allocOK)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_MEM_FAIL);
      return (ARK_MEM_FAIL);
    }

    /* All allocations are complete */
    ark_mem->MallocDone = SUNTRUE;
  }

  /* All allocation and error checking is complete at this point */

  /* Copy the input parameters into ARKODE state */
  ark_mem->tcur = t0;
  ark_mem->tn   = t0;

  /* Initialize yn */
  N_VScale(ONE, y0, ark_mem->yn);
  ark_mem->fn_is_current = SUNFALSE;

  /* Clear any previous 'tstop' */
  ark_mem->tstopset = SUNFALSE;

  /* Initializations on (re-)initialization call, skip on reset */
  if (init_type == FIRST_INIT)
  {
    /* Counters */
    ark_mem->nst_attempts = 0;
    ark_mem->nst          = 0;
    ark_mem->nhnil        = 0;
    ark_mem->ncfn         = 0;
    ark_mem->netf         = 0;
    ark_mem->nconstrfails = 0;

    /* Initial, old, and next step sizes */
    ark_mem->h0u    = ZERO;
    ark_mem->hold   = ZERO;
    ark_mem->next_h = ZERO;

    /* Tolerance scale factor */
    ark_mem->tolsf = ONE;

    /* Reset error controller object */
    if (ark_mem->hadapt_mem->hcontroller)
    {
      retval = SUNAdaptController_Reset(ark_mem->hadapt_mem->hcontroller);
      if (retval != SUN_SUCCESS)
      {
        arkProcessError(ark_mem, ARK_CONTROLLER_ERR, __LINE__, __func__,
                        __FILE__, "Unable to reset error controller object");
        return (ARK_CONTROLLER_ERR);
      }
    }

    /* Adaptivity counters */
    ark_mem->hadapt_mem->nst_acc = 0;
    ark_mem->hadapt_mem->nst_exp = 0;

    /* Accumulated error estimate */
    ark_mem->AccumError = ZERO;

    /* Indicate that calling the full RHS function is not required, this flag is
       updated to SUNTRUE by the interpolation module initialization function
       and/or the stepper initialization function in arkInitialSetup */
    ark_mem->call_fullrhs = SUNFALSE;

    /* Adjoint related */
    ark_mem->checkpoint_step_idx = 0;

    /* Indicate that initialization has not been done before */
    ark_mem->initialized = SUNFALSE;
  }

  /* Indicate initialization is needed */
  ark_mem->initsetup  = SUNTRUE;
  ark_mem->init_type  = init_type;
  ark_mem->firststage = SUNTRUE;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkCheckTimestepper:

  This routine checks if all required time stepper function
  pointers have been supplied.  If any of them is missing it
  returns SUNFALSE.
  ---------------------------------------------------------------*/
sunbooleantype arkCheckTimestepper(ARKodeMem ark_mem)
{
  if ((ark_mem->step_init == NULL) || (ark_mem->step == NULL) ||
      (ark_mem->step_mem == NULL))
  {
    return (SUNFALSE);
  }
  return (SUNTRUE);
}

/*---------------------------------------------------------------
  arkCheckNvectorRequired:

  This routine checks if all absolutely-required vector
  operations are present.  If any of them is missing it returns
  SUNFALSE.
  ---------------------------------------------------------------*/
sunbooleantype arkCheckNvectorRequired(N_Vector tmpl)
{
  if ((tmpl->ops->nvclone == NULL) || (tmpl->ops->nvdestroy == NULL) ||
      (tmpl->ops->nvlinearsum == NULL) || (tmpl->ops->nvconst == NULL) ||
      (tmpl->ops->nvdiv == NULL) || (tmpl->ops->nvscale == NULL) ||
      (tmpl->ops->nvabs == NULL) || (tmpl->ops->nvinv == NULL) ||
      (tmpl->ops->nvwrmsnorm == NULL))
  {
    return (SUNFALSE);
  }
  else { return (SUNTRUE); }
}

/*---------------------------------------------------------------
  arkCheckNvectorOptional:

  This routine perform conditional checks on required vector
  operations are present (i.e., if the current ARKODE
  configuration requires additional N_Vector routines).  If any
  of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
sunbooleantype arkCheckNvectorOptional(ARKodeMem ark_mem)
{
  /* If using a built-in routine for error/residual weights with abstol==0,
     ensure that N_VMin is available */
  if ((!ark_mem->user_efun) && (ark_mem->atolmin0) &&
      (!ark_mem->tempv1->ops->nvmin))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "N_VMin unimplemented (required by error-weight function)");
    return (SUNFALSE);
  }
  if ((!ark_mem->user_rfun) && (!ark_mem->rwt_is_ewt) && (ark_mem->Ratolmin0) &&
      (!ark_mem->tempv1->ops->nvmin))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "N_VMin unimplemented (required by residual-weight function)");
    return (SUNFALSE);
  }

  /* If the user has not specified a step size (and it will be estimated
     internally), ensure that N_VDiv and N_VMaxNorm are available */
  if ((ark_mem->h0u == ZERO) && (ark_mem->hin == ZERO) &&
      (!ark_mem->tempv1->ops->nvdiv))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "N_VDiv unimplemented (required for initial step estimation)");
    return (SUNFALSE);
  }
  if ((ark_mem->h0u == ZERO) && (ark_mem->hin == ZERO) &&
      (!ark_mem->tempv1->ops->nvmaxnorm))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__,
                    __FILE__, "N_VMaxNorm unimplemented (required for initial step estimation)");
    return (SUNFALSE);
  }

  /* If using a scalar-valued absolute tolerance (for either the state or
     residual), then ensure that N_VAddConst is available */
  if ((ark_mem->itol == ARK_SS) && (!ark_mem->tempv1->ops->nvaddconst))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "N_VAddConst unimplemented (required for scalar abstol)");
    return (SUNFALSE);
  }
  if ((!ark_mem->rwt_is_ewt) && (ark_mem->ritol == ARK_SS) &&
      (!ark_mem->tempv1->ops->nvaddconst))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "N_VAddConst unimplemented (required for scalar rabstol)");
    return (SUNFALSE);
  }

  /* If we made it here, then the vector is sufficient */
  return (SUNTRUE);
}

/*---------------------------------------------------------------
  arkInitialSetup

  This routine performs all necessary items to prepare ARKODE for
  the first internal step after initialization, reinitialization,
  a reset() call, or a resize() call, including:
  - input consistency checks
  - (re)initializes the stepper
  - computes error and residual weights
  - (re)initialize the interpolation structure
  - checks for valid initial step input or estimates first step
  - checks for approach to tstop
  - checks for root near t0
  ---------------------------------------------------------------*/
int arkInitialSetup(ARKodeMem ark_mem, sunrealtype tout)
{
  int retval, hflag, istate;
  sunrealtype tout_hin, rh, htmp;
  sunbooleantype conOK;

  /* Is tout too close to tn? */
  sunrealtype tdist  = SUNRabs(tout - ark_mem->tcur);
  sunrealtype tround = ark_mem->uround *
                       SUNMAX(SUNRabs(ark_mem->tcur), SUNRabs(tout));

  if (tdist == ZERO || tdist < TWO * tround)
  {
    arkProcessError(ark_mem, ARK_TOO_CLOSE, __LINE__, __func__, __FILE__,
                    MSG_ARK_TOO_CLOSE);
    return (ARK_TOO_CLOSE);
  }

  /* Check that user has supplied an initial step size if fixedstep mode is on */
  if ((ark_mem->fixedstep) && (ark_mem->hin == ZERO))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Fixed step mode enabled, but no step size set");
    return (ARK_ILL_INPUT);
  }

  /* Perform additional N_Vector checks here, now that ARKODE has been
     fully configured by the user */
  if (!arkCheckNvectorOptional(ark_mem))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    MSG_ARK_BAD_NVECTOR);
    return (ARK_ILL_INPUT);
  }

  /* Test input tstop for legality (correct direction of integration) */
  if (ark_mem->tstopset)
  {
    htmp = (ark_mem->h == ZERO) ? tout - ark_mem->tcur : ark_mem->h;
    if ((ark_mem->tstop - ark_mem->tcur) * htmp <= ZERO)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      MSG_ARK_BAD_TSTOP, ark_mem->tstop, ark_mem->tcur);
      return (ARK_ILL_INPUT);
    }
  }

  /* Check to see if y0 satisfies constraints */
  if (ark_mem->constraints)
  {
    conOK = N_VConstrMask(ark_mem->constraints, ark_mem->yn, ark_mem->tempv1);
    if (!conOK)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      MSG_ARK_Y0_FAIL_CONSTR);
      return (ARK_ILL_INPUT);
    }
  }

  /* Load initial error weights */
  retval = ark_mem->efun(ark_mem->yn, ark_mem->ewt, ark_mem->e_data);
  if (retval != 0)
  {
    if (ark_mem->itol == ARK_WF)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      MSG_ARK_EWT_FAIL);
    }
    else
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      MSG_ARK_BAD_EWT);
    }
    return (ARK_ILL_INPUT);
  }

  /* Set up the time stepper module */
  if (ark_mem->step_init == NULL)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Time stepper module is missing");
    return (ARK_ILL_INPUT);
  }
  retval = ark_mem->step_init(ark_mem, tout, ark_mem->init_type);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, retval, __LINE__, __func__, __FILE__,
                    "Error in initialization of time stepper module");
    return (retval);
  }

  /* Load initial residual weights */
  if (ark_mem->rwt_is_ewt)
  { /* update pointer to ewt */
    ark_mem->rwt = ark_mem->ewt;
  }
  else
  {
    retval = ark_mem->rfun(ark_mem->yn, ark_mem->rwt, ark_mem->r_data);
    if (retval != 0)
    {
      if (ark_mem->itol == ARK_WF)
      {
        arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                        MSG_ARK_RWT_FAIL);
      }
      else
      {
        arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                        MSG_ARK_BAD_RWT);
      }
      return (ARK_ILL_INPUT);
    }
  }

  /* Create default interpolation module (if needed) */
  if (ark_mem->interp_type != ARK_INTERP_NONE && !(ark_mem->interp))
  {
    if (ark_mem->interp_type == ARK_INTERP_LAGRANGE)
    {
      ark_mem->interp = arkInterpCreate_Lagrange(ark_mem, ark_mem->interp_degree);
    }
    else
    {
      ark_mem->interp = arkInterpCreate_Hermite(ark_mem, ark_mem->interp_degree);
    }
    if (ark_mem->interp == NULL)
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      "Unable to allocate interpolation module");
      return ARK_MEM_FAIL;
    }
  }

  /* Fill initial interpolation data (if needed) */
  if (ark_mem->interp != NULL)
  {
    /* Stepper init may have limited the interpolation degree */
    if (arkInterpSetDegree(ark_mem, ark_mem->interp, ark_mem->interp_degree))
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Unable to update interpolation polynomial degree");
      return ARK_ILL_INPUT;
    }

    if (arkInterpInit(ark_mem, ark_mem->interp, ark_mem->tcur))
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      "Unable to initialize interpolation module");
      return ARK_ILL_INPUT;
    }
  }

  /* Check if the configuration requires interpolation */
  if (ark_mem->root_mem && !(ark_mem->interp))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Rootfinding requires an interpolation module");
    return ARK_ILL_INPUT;
  }

  if (ark_mem->tstopinterp && !(ark_mem->interp))
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Stop time interpolation requires an interpolation module");
    return ARK_ILL_INPUT;
  }

  /* If fullrhs will be called (to estimate initial step, explicit steppers, Hermite
     interpolation module, and possibly (but not always) arkRootCheck1), then
     ensure that it is provided, and space is allocated for fn.  Otherwise,
     we should free ark_mem->fn if it is allocated. */
  if (ark_mem->call_fullrhs || (ark_mem->h0u == ZERO && ark_mem->hin == ZERO) ||
      ark_mem->root_mem)
  {
    if (!ark_mem->step_fullrhs)
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      MSG_ARK_MISSING_FULLRHS);
      return ARK_ILL_INPUT;
    }

    if (!arkAllocVec(ark_mem, ark_mem->yn, &ark_mem->fn))
    {
      arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                      MSG_ARK_MEM_FAIL);
      return (ARK_MEM_FAIL);
    }
  }
  else
  {
    if (ark_mem->fn != NULL) { arkFreeVec(ark_mem, &ark_mem->fn); }
  }

  /* initialization complete */
  ark_mem->initialized = SUNTRUE;

  /* Set initial step size */
  if (ark_mem->h0u == ZERO)
  {
    /* Check input h for validity */
    ark_mem->h = ark_mem->hin;
    if ((ark_mem->h != ZERO) && ((tout - ark_mem->tcur) * ark_mem->h < ZERO))
    {
      arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                      MSG_ARK_BAD_H0);
      return (ARK_ILL_INPUT);
    }

    /* Estimate initial h if not set */
    if (ark_mem->h == ZERO)
    {
      /* If necessary, temporarily set h as it is used to compute the tolerance
         in a potential mass matrix solve when computing the full rhs */
      ark_mem->h = SUNRabs(tout - ark_mem->tcur);
      if (ark_mem->h == ZERO) { ark_mem->h = ONE; }

      /* Estimate the first step size */
      tout_hin = tout;
      if (ark_mem->tstopset &&
          (tout - ark_mem->tcur) * (tout - ark_mem->tstop) > ZERO)
      {
        tout_hin = ark_mem->tstop;
      }
      hflag = arkHin(ark_mem, tout_hin);
      if (hflag != ARK_SUCCESS)
      {
        istate = arkHandleFailure(ark_mem, hflag);
        return (istate);
      }

      /* Use first step growth factor for estimated h */
      ark_mem->hadapt_mem->etamax = ark_mem->hadapt_mem->etamx1;
    }
    else if (ark_mem->nst == 0)
    {
      /* Use first step growth factor for user defined h */
      ark_mem->hadapt_mem->etamax = ark_mem->hadapt_mem->etamx1;
    }
    else
    {
      /* Use standard growth factor (e.g., for reset) */
      ark_mem->hadapt_mem->etamax = ark_mem->hadapt_mem->growth;
    }

    /* Enforce step size bounds */
    rh = SUNRabs(ark_mem->h) * ark_mem->hmax_inv;
    if (rh > ONE) { ark_mem->h /= rh; }
    if (SUNRabs(ark_mem->h) < ark_mem->hmin)
    {
      ark_mem->h *= ark_mem->hmin / SUNRabs(ark_mem->h);
    }

    /* Check for approach to tstop */
    if (ark_mem->tstopset)
    {
      if ((ark_mem->tcur + ark_mem->h - ark_mem->tstop) * ark_mem->h > ZERO)
      {
        ark_mem->h = (ark_mem->tstop - ark_mem->tcur) *
                     (ONE - FOUR * ark_mem->uround);
      }
    }

    /* Set initial time step factors */
    ark_mem->h0u    = ark_mem->h;
    ark_mem->eta    = ONE;
    ark_mem->hprime = ark_mem->h;
  }
  else
  {
    /* If next step would overtake tstop, adjust stepsize */
    if (ark_mem->tstopset)
    {
      if ((ark_mem->tcur + ark_mem->hprime - ark_mem->tstop) * ark_mem->h > ZERO)
      {
        ark_mem->hprime = (ark_mem->tstop - ark_mem->tcur) *
                          (ONE - FOUR * ark_mem->uround);
        ark_mem->eta = ark_mem->hprime / ark_mem->h;
      }
    }
  }

  /* Check for zeros of root function g at and near t0. */
  if (ark_mem->root_mem != NULL)
  {
    if (ark_mem->root_mem->nrtfn > 0)
    {
      retval = arkRootCheck1((void*)ark_mem);
      if (retval != ARK_SUCCESS) { return retval; }
    }
  }

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkStopTests

  This routine performs relevant stopping tests:
  - check for root in last step
  - check if we passed tstop
  - check if we passed tout (NORMAL mode)
  - check if current tn was returned (ONE_STEP mode)
  - check if we are close to tstop
  (adjust step size if needed)
  ---------------------------------------------------------------*/
int arkStopTests(ARKodeMem ark_mem, sunrealtype tout, N_Vector yout,
                 sunrealtype* tret, int itask, int* ier)
{
  int irfndp, retval;
  sunrealtype troundoff;

  /* Estimate an infinitesimal time interval to be used as
     a roundoff for time quantities (based on current time
     and step size) */
  troundoff = FUZZ_FACTOR * ark_mem->uround *
              (SUNRabs(ark_mem->tcur) + SUNRabs(ark_mem->h));

  /* First, check for a root in the last step taken, other than the
     last root found, if any.  If itask = ARK_ONE_STEP and y(tn) was not
     returned because of an intervening root, return y(tn) now.     */
  if (ark_mem->root_mem != NULL)
  {
    if (ark_mem->root_mem->nrtfn > 0)
    {
      /* Shortcut to roots found in previous step */
      irfndp = ark_mem->root_mem->irfnd;

      /* If the full RHS was not computed in the last call to arkCompleteStep
         and roots were found in the previous step, then compute the full rhs
         for possible use in arkRootCheck2 (not always necessary) */
      if (!(ark_mem->fn_is_current) && irfndp != 0)
      {
        retval = ark_mem->step_fullrhs(ark_mem, ark_mem->tn, ark_mem->yn,
                                       ark_mem->fn, ARK_FULLRHS_END);
        if (retval != 0)
        {
          arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__,
                          __FILE__, MSG_ARK_RHSFUNC_FAILED);
          *ier = ARK_RHSFUNC_FAIL;
          return (1);
        }
        ark_mem->fn_is_current = SUNTRUE;
      }

      retval = arkRootCheck2((void*)ark_mem);

      if (retval == CLOSERT)
      {
        arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                        MSG_ARK_CLOSE_ROOTS, ark_mem->root_mem->tlo);
        *ier = ARK_ILL_INPUT;
        return (1);
      }
      else if (retval == ARK_RTFUNC_FAIL)
      {
        arkProcessError(ark_mem, ARK_RTFUNC_FAIL, __LINE__, __func__, __FILE__,
                        MSG_ARK_RTFUNC_FAILED, ark_mem->root_mem->tlo);
        *ier = ARK_RTFUNC_FAIL;
        return (1);
      }
      else if (retval == RTFOUND)
      {
        ark_mem->tretlast = *tret = ark_mem->root_mem->tlo;
        *ier                      = ARK_ROOT_RETURN;
        return (1);
      }

      /* If tn is distinct from tretlast (within roundoff),
         check remaining interval for roots */
      if (SUNRabs(ark_mem->tcur - ark_mem->tretlast) > troundoff)
      {
        retval = arkRootCheck3((void*)ark_mem, tout, itask);

        if (retval == ARK_SUCCESS)
        { /* no root found */
          ark_mem->root_mem->irfnd = 0;
          if ((irfndp == 1) && (itask == ARK_ONE_STEP))
          {
            ark_mem->tretlast = *tret = ark_mem->tcur;
            N_VScale(ONE, ark_mem->yn, yout);
            *ier = ARK_SUCCESS;
            return (1);
          }
        }
        else if (retval == RTFOUND)
        { /* a new root was found */
          ark_mem->root_mem->irfnd = 1;
          ark_mem->tretlast = *tret = ark_mem->root_mem->tlo;
          *ier                      = ARK_ROOT_RETURN;
          return (1);
        }
        else if (retval == ARK_RTFUNC_FAIL)
        { /* g failed */
          arkProcessError(ark_mem, ARK_RTFUNC_FAIL, __LINE__, __func__, __FILE__,
                          MSG_ARK_RTFUNC_FAILED, ark_mem->root_mem->tlo);
          *ier = ARK_RTFUNC_FAIL;
          return (1);
        }
      }

    } /* end of root stop check */
  }

  /* Test for tn at tstop or near tstop */
  if (ark_mem->tstopset)
  {
    if (SUNRabs(ark_mem->tcur - ark_mem->tstop) <= troundoff)
    {
      /* Ensure tout >= tstop, otherwise check for tout return below */
      if ((tout - ark_mem->tstop) * ark_mem->h >= ZERO ||
          SUNRabs(tout - ark_mem->tstop) <= troundoff)
      {
        if (ark_mem->tstopinterp && ark_mem->interp)
        {
          *ier = ARKodeGetDky(ark_mem, ark_mem->tstop, 0, yout);
          if (*ier != ARK_SUCCESS)
          {
            arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                            MSG_ARK_BAD_TSTOP, ark_mem->tstop, ark_mem->tcur);
            *ier = ARK_ILL_INPUT;
            return (1);
          }
        }
        else { N_VScale(ONE, ark_mem->yn, yout); }
        ark_mem->tretlast = *tret = ark_mem->tstop;
        ark_mem->tstopset         = SUNFALSE;
        *ier                      = ARK_TSTOP_RETURN;
        return (1);
      }
    }
    /* If next step would overtake tstop, adjust stepsize */
    else if ((ark_mem->tcur + ark_mem->hprime - ark_mem->tstop) * ark_mem->h >
             ZERO)
    {
      ark_mem->hprime = (ark_mem->tstop - ark_mem->tcur) *
                        (ONE - FOUR * ark_mem->uround);
      ark_mem->eta = ark_mem->hprime / ark_mem->h;
    }
  }

  /* In ARK_NORMAL mode, test if tout was reached */
  if ((itask == ARK_NORMAL) && ((ark_mem->tcur - tout) * ark_mem->h >= ZERO))
  {
    if (ark_mem->interp)
    {
      *ier = ARKodeGetDky(ark_mem, tout, 0, yout);
      if (*ier != ARK_SUCCESS)
      {
        arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                        MSG_ARK_BAD_TOUT, tout);
        *ier = ARK_ILL_INPUT;
        return (1);
      }
      ark_mem->tretlast = *tret = tout;
    }
    else
    {
      N_VScale(ONE, ark_mem->yn, yout);
      ark_mem->tretlast = *tret = ark_mem->tcur;
    }
    *ier = ARK_SUCCESS;
    return (1);
  }

  /* In ARK_ONE_STEP mode, test if tn was returned */
  if (itask == ARK_ONE_STEP &&
      SUNRabs(ark_mem->tcur - ark_mem->tretlast) > troundoff)
  {
    ark_mem->tretlast = *tret = ark_mem->tcur;
    N_VScale(ONE, ark_mem->yn, yout);
    *ier = ARK_SUCCESS;
    return (1);
  }

  return (0);
}

/*---------------------------------------------------------------
  arkHin

  This routine computes a tentative initial step size h0.
  Note that here tout is either the value passed to ARKodeEvolve
  at the first call or the value of tstop (if tstop is enabled and
  it is closer to t0=tn than tout). If the RHS function fails
  unrecoverably, arkHin returns ARK_RHSFUNC_FAIL. If the RHS
  function fails recoverably too many times and recovery is not
  possible, arkHin returns ARK_REPTD_RHSFUNC_ERR. Otherwise, arkHin
  sets h to the chosen value h0 and returns ARK_SUCCESS.

  The algorithm used seeks to find h0 as a solution of
  (WRMS norm of (h0^2 ydd / 2)) = 1,
  where ydd = estimated second derivative of y.  Although this
  choice is based on an error expansion of the Backward Euler
  method, and hence results in an overly-conservative time step
  for our higher-order ARK methods, it does find an order-of-
  magnitude estimate of the initial time scale of the solution.
  Since this method is only used on the first time step, the
  additional caution will not overly hinder solver efficiency.

  We start with an initial estimate equal to the geometric mean
  of the lower and upper bounds on the step size.

  Loop up to H0_ITERS times to find h0.
  Stop if new and previous values differ by a factor < 2.
  Stop if hnew/hg > 2 after one iteration, as this probably
  means that the ydd value is bad because of cancellation error.

  For each new proposed hg, we allow H0_ITERS attempts to
  resolve a possible recoverable failure from f() by reducing
  the proposed stepsize by a factor of 0.2. If a legal stepsize
  still cannot be found, fall back on a previous value if
  possible, or else return ARK_REPTD_RHSFUNC_ERR.

  Finally, we apply a bias (0.5) and verify that h0 is within
  bounds.
  ---------------------------------------------------------------*/
int arkHin(ARKodeMem ark_mem, sunrealtype tout)
{
  int retval, sign, count1, count2;
  sunrealtype tdiff, tdist, tround, hlb, hub;
  sunrealtype hg, hgs, hs, hnew, hrat, h0, yddnrm;
  sunbooleantype hgOK;

  /* arkInitialSetup checks for tdiff = 0 or < 2 * troundoff */
  tdiff  = tout - ark_mem->tcur;
  sign   = (tdiff > ZERO) ? 1 : -1;
  tdist  = SUNRabs(tdiff);
  tround = ark_mem->uround * SUNMAX(SUNRabs(ark_mem->tcur), SUNRabs(tout));

  /* call full RHS if needed */
  if (!(ark_mem->fn_is_current))
  {
    /* NOTE: The step size (h) is used in setting the tolerance in a potential
       mass matrix solve when computing the full RHS. Before calling arkHin, h
       is set to |tout - tcur| or 1 and so we do not need to guard against
       h == 0 here before calling the full RHS. */
    retval = ark_mem->step_fullrhs(ark_mem, ark_mem->tn, ark_mem->yn,
                                   ark_mem->fn, ARK_FULLRHS_START);
    if (retval) { return ARK_RHSFUNC_FAIL; }
    ark_mem->fn_is_current = SUNTRUE;
  }

  /* Set lower and upper bounds on h0, and take geometric mean
     as first trial value.
     Exit with this value if the bounds cross each other. */
  hlb = H0_LBFACTOR * tround;
  hub = arkUpperBoundH0(ark_mem, tdist);

  hg = SUNRsqrt(hlb * hub);

  if (hub < hlb)
  {
    if (sign == -1) { ark_mem->h = -hg; }
    else { ark_mem->h = hg; }
    return (ARK_SUCCESS);
  }

  /* Outer loop */
  hs = hg; /* safeguard against 'uninitialized variable' warning */
  for (count1 = 1; count1 <= H0_ITERS; count1++)
  {
    /* Attempts to estimate ydd */
    hgOK = SUNFALSE;

    for (count2 = 1; count2 <= H0_ITERS; count2++)
    {
      hgs    = hg * sign;
      retval = arkYddNorm(ark_mem, hgs, &yddnrm);
      /* If f() failed unrecoverably, give up */
      if (retval < 0) { return (ARK_RHSFUNC_FAIL); }
      /* If successful, we can use ydd */
      if (retval == ARK_SUCCESS)
      {
        hgOK = SUNTRUE;
        break;
      }
      /* f() failed recoverably; cut step size and test it again */
      hg *= SUN_RCONST(0.2);
    }

    /* If f() failed recoverably H0_ITERS times */
    if (!hgOK)
    {
      /* Exit if this is the first or second pass. No recovery possible */
      if (count1 <= 2) { return (ARK_REPTD_RHSFUNC_ERR); }
      /* We have a fall-back option. The value hs is a previous hnew which
         passed through f(). Use it and break */
      hnew = hs;
      break;
    }

    /* The proposed step size is feasible. Save it. */
    hs = hg;

    /* Propose new step size */
    hnew = (yddnrm * hub * hub > TWO) ? SUNRsqrt(TWO / yddnrm)
                                      : SUNRsqrt(hg * hub);

    /* If last pass, stop now with hnew */
    if (count1 == H0_ITERS) { break; }

    hrat = hnew / hg;

    /* Accept hnew if it does not differ from hg by more than a factor of 2 */
    if ((hrat > HALF) && (hrat < TWO)) { break; }

    /* After one pass, if ydd seems to be bad, use fall-back value. */
    if ((count1 > 1) && (hrat > TWO))
    {
      hnew = hg;
      break;
    }

    /* Send this value back through f() */
    hg = hnew;
  }

  /* Apply bounds, bias factor, and attach sign */
  h0 = H0_BIAS * hnew;
  if (h0 < hlb) { h0 = hlb; }
  if (h0 > hub) { h0 = hub; }
  if (sign == -1) { h0 = -h0; }
  ark_mem->h = h0;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkUpperBoundH0

  This routine sets an upper bound on abs(h0) based on
  tdist = tn - t0 and the values of y[i]/y'[i].
  ---------------------------------------------------------------*/
sunrealtype arkUpperBoundH0(ARKodeMem ark_mem, sunrealtype tdist)
{
  sunrealtype hub_inv, hub;
  N_Vector temp1, temp2;

  /* Bound based on |y0|/|y0'| -- allow at most an increase of
   * H0_UBFACTOR in y0 (based on a forward Euler step). The weight
   * factor is used as a safeguard against zero components in y0. */
  temp1 = ark_mem->tempv1;
  temp2 = ark_mem->tempv2;

  N_VAbs(ark_mem->yn, temp2);
  ark_mem->efun(ark_mem->yn, temp1, ark_mem->e_data);
  N_VInv(temp1, temp1);
  N_VLinearSum(H0_UBFACTOR, temp2, ONE, temp1, temp1);

  N_VAbs(ark_mem->fn, temp2);

  N_VDiv(temp2, temp1, temp1);
  hub_inv = N_VMaxNorm(temp1);

  /* bound based on tdist -- allow at most a step of magnitude
   * H0_UBFACTOR * tdist */
  hub = H0_UBFACTOR * tdist;

  /* Use the smaller of the two */
  if (hub * hub_inv > ONE) { hub = ONE / hub_inv; }

  return (hub);
}

/*---------------------------------------------------------------
  arkYddNorm

  This routine computes an estimate of the second derivative of y
  using a difference quotient, and returns its WRMS norm.
  ---------------------------------------------------------------*/
int arkYddNorm(ARKodeMem ark_mem, sunrealtype hg, sunrealtype* yddnrm)
{
  int retval;

  /* increment y with a multiple of f */
  N_VLinearSum(hg, ark_mem->fn, ONE, ark_mem->yn, ark_mem->ycur);

  /* compute y', via the ODE RHS routine */
  retval = ark_mem->step_fullrhs(ark_mem, ark_mem->tcur + hg, ark_mem->ycur,
                                 ark_mem->tempv1, ARK_FULLRHS_OTHER);
  if (retval != 0) { return (ARK_RHSFUNC_FAIL); }

  /* difference new f and original f to estimate y'' */
  N_VLinearSum(ONE / hg, ark_mem->tempv1, -ONE / hg, ark_mem->fn,
               ark_mem->tempv1);

  /* reset ycur to equal yn (unnecessary?) */
  N_VScale(ONE, ark_mem->yn, ark_mem->ycur);

  /* compute norm of y'' */
  *yddnrm = N_VWrmsNorm(ark_mem->tempv1, ark_mem->ewt);

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkCompleteStep

  This routine performs various update operations when the step
  solution is complete.  It is assumed that the timestepper
  module has stored the time-evolved solution in ark_mem->ycur,
  and the step that gave rise to this solution in ark_mem->h.
  We update the current time (tn), the current solution (yn),
  increment the overall step counter nst, record the values hold
  and tnew, allow for user-provided postprocessing, and update
  the interpolation structure.
  ---------------------------------------------------------------*/
int arkCompleteStep(ARKodeMem ark_mem, sunrealtype dsm)
{
  int retval;
  sunrealtype troundoff;

  /* Set current time to the end of the step (in case the last
     stage time does not coincide with the step solution time).
     If tstop is enabled, it is possible for tn + h to be past
     tstop by roundoff, and in that case, we reset tn (after
     incrementing by h) to tstop. */

  /* During long-time integration, roundoff can creep into tcur.
     Compensated summation fixes this but with increased cost, so it is optional. */
  if (ark_mem->use_compensated_sums)
  {
    sunCompensatedSum(ark_mem->tn, ark_mem->h, &ark_mem->tcur, &ark_mem->terr);
  }
  else { ark_mem->tcur = ark_mem->tn + ark_mem->h; }

  if (ark_mem->tstopset)
  {
    troundoff = FUZZ_FACTOR * ark_mem->uround *
                (SUNRabs(ark_mem->tcur) + SUNRabs(ark_mem->h));
    if (SUNRabs(ark_mem->tcur - ark_mem->tstop) <= troundoff)
    {
      ark_mem->tcur = ark_mem->tstop;
    }
  }

  /* store this step's contribution to accumulated temporal error */
  if (ark_mem->AccumErrorType != ARK_ACCUMERROR_NONE)
  {
    if (ark_mem->AccumErrorType == ARK_ACCUMERROR_MAX)
    {
      ark_mem->AccumError = SUNMAX(dsm, ark_mem->AccumError);
    }
    else if (ark_mem->AccumErrorType == ARK_ACCUMERROR_SUM)
    {
      ark_mem->AccumError += dsm;
    }
    else /* ARK_ACCUMERROR_AVG */ { ark_mem->AccumError += (dsm * ark_mem->h); }
  }

  /* apply user-supplied step postprocessing function (if supplied) */
  if (ark_mem->ProcessStep != NULL)
  {
    retval = ark_mem->ProcessStep(ark_mem->tcur, ark_mem->ycur, ark_mem->ps_data);
    if (retval != 0) { return (ARK_POSTPROCESS_STEP_FAIL); }
  }

  /* update interpolation structure

     NOTE: This must be called before updating yn with ycur as the interpolation
     module may need to save tn, yn from the start of this step. */
  if (ark_mem->interp != NULL)
  {
    retval = arkInterpUpdate(ark_mem, ark_mem->interp, ark_mem->tcur);
    if (retval != ARK_SUCCESS) { return (retval); }
  }

  /* update yn to current solution */
  N_VScale(ONE, ark_mem->ycur, ark_mem->yn);
  ark_mem->fn_is_current = SUNFALSE;

  /* Notify time step controller object of successful step */
  if (ark_mem->hadapt_mem->hcontroller)
  {
    retval = SUNAdaptController_UpdateH(ark_mem->hadapt_mem->hcontroller,
                                        ark_mem->h, dsm);
    if (retval != SUN_SUCCESS)
    {
      arkProcessError(ark_mem, ARK_CONTROLLER_ERR, __LINE__, __func__, __FILE__,
                      "Failure updating controller object");
      return (ARK_CONTROLLER_ERR);
    }
  }

  /* update scalar quantities */
  ark_mem->nst++;
  ark_mem->checkpoint_step_idx++;
  ark_mem->hold   = ark_mem->h;
  ark_mem->tn     = ark_mem->tcur;
  ark_mem->hprime = ark_mem->h * ark_mem->eta;

  /* Reset growth factor for subsequent time step */
  ark_mem->hadapt_mem->etamax = ark_mem->hadapt_mem->growth;

  /* Turn off flag indicating initial step and first stage */
  ark_mem->initsetup  = SUNFALSE;
  ark_mem->firststage = SUNFALSE;

  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkHandleFailure

  This routine prints error messages for all cases of failure by
  arkHin and ark_step. It returns to ARKODE the value that ARKODE
  is to return to the user.
  ---------------------------------------------------------------*/
int arkHandleFailure(ARKodeMem ark_mem, int flag)
{
  /* Depending on flag, print error message and return error flag */
  switch (flag)
  {
  case ARK_ERR_FAILURE:
    arkProcessError(ark_mem, ARK_ERR_FAILURE, __LINE__, __func__, __FILE__,
                    MSG_ARK_ERR_FAILS, ark_mem->tcur, ark_mem->h);
    break;
  case ARK_CONV_FAILURE:
    arkProcessError(ark_mem, ARK_CONV_FAILURE, __LINE__, __func__, __FILE__,
                    MSG_ARK_CONV_FAILS, ark_mem->tcur, ark_mem->h);
    break;
  case ARK_LSETUP_FAIL:
    arkProcessError(ark_mem, ARK_LSETUP_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_SETUP_FAILED, ark_mem->tcur);
    break;
  case ARK_LSOLVE_FAIL:
    arkProcessError(ark_mem, ARK_LSOLVE_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_SOLVE_FAILED, ark_mem->tcur);
    break;
  case ARK_RHSFUNC_FAIL:
    arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_RHSFUNC_FAILED, ark_mem->tcur);
    break;
  case ARK_UNREC_RHSFUNC_ERR:
    arkProcessError(ark_mem, ARK_UNREC_RHSFUNC_ERR, __LINE__, __func__,
                    __FILE__, MSG_ARK_RHSFUNC_UNREC, ark_mem->tcur);
    break;
  case ARK_REPTD_RHSFUNC_ERR:
    arkProcessError(ark_mem, ARK_REPTD_RHSFUNC_ERR, __LINE__, __func__,
                    __FILE__, MSG_ARK_RHSFUNC_REPTD, ark_mem->tcur);
    break;
  case ARK_RTFUNC_FAIL:
    arkProcessError(ark_mem, ARK_RTFUNC_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_RTFUNC_FAILED, ark_mem->tcur);
    break;
  case ARK_TOO_CLOSE:
    arkProcessError(ark_mem, ARK_TOO_CLOSE, __LINE__, __func__, __FILE__,
                    MSG_ARK_TOO_CLOSE);
    break;
  case ARK_CONSTR_FAIL:
    arkProcessError(ark_mem, ARK_CONSTR_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_FAILED_CONSTR, ark_mem->tcur);
    break;
  case ARK_MASSSOLVE_FAIL:
    arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_MASSSOLVE_FAIL);
    break;
  case ARK_NLS_SETUP_FAIL:
    arkProcessError(ark_mem, ARK_NLS_SETUP_FAIL, __LINE__, __func__, __FILE__,
                    "At t = " SUN_FORMAT_G
                    " the nonlinear solver setup failed unrecoverably",
                    ark_mem->tcur);
    break;
  case ARK_VECTOROP_ERR:
    arkProcessError(ark_mem, ARK_VECTOROP_ERR, __LINE__, __func__, __FILE__,
                    MSG_ARK_VECTOROP_ERR, ark_mem->tcur);
    break;
  case ARK_INNERSTEP_FAIL:
    arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, __LINE__, __func__, __FILE__,
                    MSG_ARK_INNERSTEP_FAILED, ark_mem->tcur);
    break;
  case ARK_NLS_OP_ERR:
    arkProcessError(ark_mem, ARK_NLS_OP_ERR, __LINE__, __func__, __FILE__,
                    MSG_ARK_NLS_FAIL, ark_mem->tcur);
    break;
  case ARK_USER_PREDICT_FAIL:
    arkProcessError(ark_mem, ARK_USER_PREDICT_FAIL, __LINE__, __func__,
                    __FILE__, MSG_ARK_USER_PREDICT_FAIL, ark_mem->tcur);
    break;
  case ARK_POSTPROCESS_STEP_FAIL:
    arkProcessError(ark_mem, ARK_POSTPROCESS_STEP_FAIL, __LINE__, __func__,
                    __FILE__, MSG_ARK_POSTPROCESS_STEP_FAIL, ark_mem->tcur);
    break;
  case ARK_POSTPROCESS_STAGE_FAIL:
    arkProcessError(ark_mem, ARK_POSTPROCESS_STAGE_FAIL, __LINE__, __func__,
                    __FILE__, MSG_ARK_POSTPROCESS_STAGE_FAIL, ark_mem->tcur);
    break;
  case ARK_INTERP_FAIL:
    arkProcessError(ark_mem, ARK_INTERP_FAIL, __LINE__, __func__, __FILE__,
                    "At t = " SUN_FORMAT_G
                    " the interpolation module failed unrecoverably",
                    ark_mem->tcur);
    break;
  case ARK_INVALID_TABLE:
    arkProcessError(ark_mem, ARK_INVALID_TABLE, __LINE__, __func__, __FILE__,
                    "ARKODE was provided an invalid method table");
    break;
  case ARK_RELAX_FAIL:
    arkProcessError(ark_mem, ARK_RELAX_FAIL, __LINE__, __func__, __FILE__,
                    "At t = " SUN_FORMAT_G " the relaxation module failed",
                    ark_mem->tcur);
    break;
  case ARK_RELAX_MEM_NULL:
    arkProcessError(ark_mem, ARK_RELAX_MEM_NULL, __LINE__, __func__, __FILE__,
                    "The ARKODE relaxation module memory is NULL");
    break;
  case ARK_RELAX_FUNC_FAIL:
    arkProcessError(ark_mem, ARK_RELAX_FUNC_FAIL, __LINE__, __func__, __FILE__,
                    "The relaxation function failed unrecoverably");
    break;
  case ARK_RELAX_JAC_FAIL:
    arkProcessError(ark_mem, ARK_RELAX_JAC_FAIL, __LINE__, __func__, __FILE__,
                    "The relaxation Jacobian failed unrecoverably");
    break;
  case ARK_ADJ_RECOMPUTE_FAIL:
    arkProcessError(ark_mem, ARK_ADJ_RECOMPUTE_FAIL, __LINE__, __func__, __FILE__,
                    "The forward recomputation of step failed unrecoverably");
    break;
  case ARK_ADJ_CHECKPOINT_FAIL:
    arkProcessError(ark_mem, ARK_ADJ_CHECKPOINT_FAIL, __LINE__, __func__,
                    __FILE__, "A checkpoint operation failed unrecoverably");
    break;
  case ARK_SUNADJSTEPPER_ERR:
    arkProcessError(ark_mem, ARK_SUNADJSTEPPER_ERR, __LINE__, __func__,
                    __FILE__, "A SUNAdjStepper operation failed unrecoverably");
    break;
  case ARK_DOMEIG_FAIL:
    arkProcessError(ark_mem, ARK_DOMEIG_FAIL, __LINE__, __func__, __FILE__,
                    "The dominant eigenvalue function failed unrecoverably");
    break;
  case ARK_MAX_STAGE_LIMIT_FAIL:
    arkProcessError(ark_mem, ARK_MAX_STAGE_LIMIT_FAIL, __LINE__, __func__,
                    __FILE__, "The max stage limit failed unrecoverably");
    break;
  case ARK_SUNSTEPPER_ERR:
    arkProcessError(ark_mem, ARK_SUNSTEPPER_ERR, __LINE__, __func__, __FILE__,
                    "An inner SUNStepper error occurred");
    break;
  default:
    /* This return should never happen */
    arkProcessError(ark_mem, ARK_UNRECOGNIZED_ERROR, __LINE__, __func__, __FILE__,
                    "ARKODE encountered an unrecognized error. Please report "
                    "this to the Sundials developers at "
                    "sundials-users@llnl.gov");
    return (ARK_UNRECOGNIZED_ERROR);
  }

  return (flag);
}

/*---------------------------------------------------------------
  arkEwtSetSS

  This routine is responsible for setting the error weight vector
  ewt as follows:

  ewt[i] = 1 / (reltol * SUNRabs(ycur[i]) + abstol), i=0,...,neq-1

  When the absolute tolerance is zero, it tests for non-positive
  components before inverting. arkEwtSetSS returns 0 if ewt is
  successfully set to a positive vector and -1 otherwise. In the
  latter case, ewt is considered undefined.
  ---------------------------------------------------------------*/
int arkEwtSetSS(N_Vector ycur, N_Vector weight, void* arkode_mem)
{
  ARKodeMem ark_mem = (ARKodeMem)arkode_mem;
  N_VAbs(ycur, ark_mem->tempv1);
  N_VScale(ark_mem->reltol, ark_mem->tempv1, ark_mem->tempv1);
  N_VAddConst(ark_mem->tempv1, ark_mem->Sabstol, ark_mem->tempv1);
  if (ark_mem->atolmin0)
  {
    if (N_VMin(ark_mem->tempv1) <= ZERO) { return (-1); }
  }
  N_VInv(ark_mem->tempv1, weight);
  return (0);
}

/*---------------------------------------------------------------
  arkEwtSetSV

  This routine is responsible for setting the error weight vector
  ewt as follows:

  ewt[i] = 1 / (reltol * SUNRabs(ycur[i]) + abstol[i]), i=0,...,neq-1

  When any absolute tolerance is zero, it tests for non-positive
  components before inverting. arkEwtSetSV returns 0 if ewt is
  successfully set to a positive vector and -1 otherwise. In the
  latter case, ewt is considered undefined.
  ---------------------------------------------------------------*/
int arkEwtSetSV(N_Vector ycur, N_Vector weight, void* arkode_mem)
{
  ARKodeMem ark_mem = (ARKodeMem)arkode_mem;
  N_VAbs(ycur, ark_mem->tempv1);
  N_VLinearSum(ark_mem->reltol, ark_mem->tempv1, ONE, ark_mem->Vabstol,
               ark_mem->tempv1);
  if (ark_mem->atolmin0)
  {
    if (N_VMin(ark_mem->tempv1) <= ZERO) { return (-1); }
  }
  N_VInv(ark_mem->tempv1, weight);
  return (0);
}

/*---------------------------------------------------------------
  arkEwtSetSmallReal

  This routine is responsible for setting the error weight vector
  ewt as follows:

  ewt[i] = SUN_SMALL_REAL

  This is routine is only used with explicit time stepping with
  a fixed step size to avoid a potential too much error return
  to the user.
  ---------------------------------------------------------------*/
int arkEwtSetSmallReal(SUNDIALS_MAYBE_UNUSED N_Vector ycur, N_Vector weight,
                       SUNDIALS_MAYBE_UNUSED void* arkode_mem)
{
  N_VConst(SUN_SMALL_REAL, weight);
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkRwtSetSS

  This routine sets rwt as described above in the case tol_type = ARK_SS.
  When the absolute tolerance is zero, it tests for non-positive
  components before inverting. arkRwtSetSS returns 0 if rwt is
  successfully set to a positive vector and -1 otherwise. In the
  latter case, rwt is considered undefined.
  ---------------------------------------------------------------*/
int arkRwtSetSS(ARKodeMem ark_mem, N_Vector My, N_Vector weight)
{
  N_VAbs(My, ark_mem->tempv1);
  N_VScale(ark_mem->reltol, ark_mem->tempv1, ark_mem->tempv1);
  N_VAddConst(ark_mem->tempv1, ark_mem->SRabstol, ark_mem->tempv1);
  if (ark_mem->Ratolmin0)
  {
    if (N_VMin(ark_mem->tempv1) <= ZERO) { return (-1); }
  }
  N_VInv(ark_mem->tempv1, weight);
  return (0);
}

/*---------------------------------------------------------------
  arkRwtSetSV

  This routine sets rwt as described above in the case tol_type = ARK_SV.
  When any absolute tolerance is zero, it tests for non-positive
  components before inverting. arkRwtSetSV returns 0 if rwt is
  successfully set to a positive vector and -1 otherwise. In the
  latter case, rwt is considered undefined.
  ---------------------------------------------------------------*/
int arkRwtSetSV(ARKodeMem ark_mem, N_Vector My, N_Vector weight)
{
  N_VAbs(My, ark_mem->tempv1);
  N_VLinearSum(ark_mem->reltol, ark_mem->tempv1, ONE, ark_mem->VRabstol,
               ark_mem->tempv1);
  if (ark_mem->Ratolmin0)
  {
    if (N_VMin(ark_mem->tempv1) <= ZERO) { return (-1); }
  }
  N_VInv(ark_mem->tempv1, weight);
  return (0);
}

/*---------------------------------------------------------------
  arkPredict_MaximumOrder

  This routine predicts the nonlinear implicit stage solution
  using the ARKode interpolation module.  This uses the
  highest-degree interpolant supported by the module (stored
  in the interpolation module).
  ---------------------------------------------------------------*/
int arkPredict_MaximumOrder(ARKodeMem ark_mem, sunrealtype tau, N_Vector yguess)
{
  /* verify that ark_mem and interpolation structure are provided */
  if (ark_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "ARKodeMem structure is NULL");
    return (ARK_MEM_NULL);
  }
  if (ark_mem->interp == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "ARKodeInterpMem structure is NULL");
    return (ARK_MEM_NULL);
  }

  /* call the interpolation module to do the work */
  return (arkInterpEvaluate(ark_mem, ark_mem->interp, tau, 0,
                            ARK_INTERP_MAX_DEGREE, yguess));
}

/*---------------------------------------------------------------
  arkPredict_VariableOrder

  This routine predicts the nonlinear implicit stage solution
  using the ARKODE interpolation module.  The degree of the
  interpolant is based on the level of extrapolation outside the
  preceding time step.
  ---------------------------------------------------------------*/
int arkPredict_VariableOrder(ARKodeMem ark_mem, sunrealtype tau, N_Vector yguess)
{
  int ord;
  sunrealtype tau_tol  = HALF;
  sunrealtype tau_tol2 = SUN_RCONST(0.75);

  /* verify that ark_mem and interpolation structure are provided */
  if (ark_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "ARKodeMem structure is NULL");
    return (ARK_MEM_NULL);
  }
  if (ark_mem->interp == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "ARKodeInterpMem structure is NULL");
    return (ARK_MEM_NULL);
  }

  /* set the polynomial order based on tau input */
  if (tau <= tau_tol) { ord = 3; }
  else if (tau <= tau_tol2) { ord = 2; }
  else { ord = 1; }

  /* call the interpolation module to do the work */
  return (arkInterpEvaluate(ark_mem, ark_mem->interp, tau, 0, ord, yguess));
}

/*---------------------------------------------------------------
  arkPredict_CutoffOrder

  This routine predicts the nonlinear implicit stage solution
  using the ARKODE interpolation module.  If the level of
  extrapolation is small enough, it uses the maximum degree
  polynomial available (stored in the interpolation module
  structure); otherwise it uses a linear polynomial.
  ---------------------------------------------------------------*/
int arkPredict_CutoffOrder(ARKodeMem ark_mem, sunrealtype tau, N_Vector yguess)
{
  int ord;
  sunrealtype tau_tol = HALF;

  /* verify that ark_mem and interpolation structure are provided */
  if (ark_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "ARKodeMem structure is NULL");
    return (ARK_MEM_NULL);
  }
  if (ark_mem->interp == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "ARKodeInterpMem structure is NULL");
    return (ARK_MEM_NULL);
  }

  /* set the polynomial order based on tau input */
  if (tau <= tau_tol) { ord = ARK_INTERP_MAX_DEGREE; }
  else { ord = 1; }

  /* call the interpolation module to do the work */
  return (arkInterpEvaluate(ark_mem, ark_mem->interp, tau, 0, ord, yguess));
}

/*---------------------------------------------------------------
  arkPredict_Bootstrap

  This routine predicts the nonlinear implicit stage solution
  using a quadratic Hermite interpolating polynomial, based on
  the data {y_n, f(t_n,y_n), f(t_n+hj,z_j)}.

  Note: we assume that ftemp = f(t_n+hj,z_j) can be computed via
     N_VLinearCombination(nvec, cvals, Xvecs, ftemp),
  i.e. the inputs cvals[0:nvec-1] and Xvecs[0:nvec-1] may be
  combined to form f(t_n+hj,z_j).
  ---------------------------------------------------------------*/
int arkPredict_Bootstrap(ARKodeMem ark_mem, sunrealtype hj, sunrealtype tau,
                         int nvec, sunrealtype* cvals, N_Vector* Xvecs,
                         N_Vector yguess)
{
  sunrealtype a0, a1, a2;
  int i, retval;

  /* verify that ark_mem and interpolation structure are provided */
  if (ark_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "ARKodeMem structure is NULL");
    return (ARK_MEM_NULL);
  }
  if (ark_mem->interp == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    "ARKodeInterpMem structure is NULL");
    return (ARK_MEM_NULL);
  }

  /* set coefficients for Hermite interpolant */
  a0 = ONE;
  a2 = tau * tau / TWO / hj;
  a1 = tau - a2;

  /* set arrays for fused vector operation; shift inputs for
     f(t_n+hj,z_j) to end of queue */
  for (i = 0; i < nvec; i++)
  {
    cvals[2 + i] = a2 * cvals[i];
    Xvecs[2 + i] = Xvecs[i];
  }
  cvals[0] = a0;
  Xvecs[0] = ark_mem->yn;
  cvals[1] = a1;
  Xvecs[1] = ark_mem->fn;

  /* call fused vector operation to compute prediction */
  retval = N_VLinearCombination(nvec + 2, cvals, Xvecs, yguess);
  if (retval != 0) { return (ARK_VECTOROP_ERR); }
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkCheckConvergence

  This routine checks the return flag from the time-stepper's
  "step" routine for algebraic solver convergence issues.

  Returns ARK_SUCCESS (0) if successful, PREDICT_AGAIN (>0)
  on a recoverable convergence failure, or a relevant
  nonrecoverable failure flag (<0).
  --------------------------------------------------------------*/
int arkCheckConvergence(ARKodeMem ark_mem, int* nflagPtr, int* ncfPtr)
{
  ARKodeHAdaptMem hadapt_mem;

  /* If nonlinear solver succeeded, return with ARK_SUCCESS */
  if (*nflagPtr == ARK_SUCCESS) { return (ARK_SUCCESS); }
  /* Returns with an ARK_RETRY_STEP flag occur at a stage well before
  any algebraic solvers are involved. On the other hand,
  the arkCheckConvergence function handles the results from algebraic
  solvers, which never take place with an ARK_RETRY_STEP flag.
  Therefore, we immediately return from arkCheckConvergence,
  as it is irrelevant in the case of an ARK_RETRY_STEP */
  if (*nflagPtr == ARK_RETRY_STEP) { return (ARK_RETRY_STEP); }

  /* The nonlinear soln. failed; increment ncfn */
  ark_mem->ncfn++;

  /* If fixed time stepping, then return with convergence failure */
  if (ark_mem->fixedstep) { return (ARK_CONV_FAILURE); }

  /* Otherwise, access adaptivity structure */
  if (ark_mem->hadapt_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARKADAPT_NO_MEM);
    return (ARK_MEM_NULL);
  }
  hadapt_mem = ark_mem->hadapt_mem;

  /* Return if lsetup, lsolve, or rhs failed unrecoverably */
  if (*nflagPtr < 0)
  {
    if (*nflagPtr == ARK_LSETUP_FAIL) { return (ARK_LSETUP_FAIL); }
    else if (*nflagPtr == ARK_LSOLVE_FAIL) { return (ARK_LSOLVE_FAIL); }
    else if (*nflagPtr == ARK_RHSFUNC_FAIL) { return (ARK_RHSFUNC_FAIL); }
    else { return (ARK_NLS_OP_ERR); }
  }

  /* At this point, nflag = CONV_FAIL or RHSFUNC_RECVR; increment ncf */
  (*ncfPtr)++;
  hadapt_mem->etamax = ONE;

  /* If we had maxncf failures, or if |h| = hmin,
     return ARK_CONV_FAILURE or ARK_REPTD_RHSFUNC_ERR. */
  if ((*ncfPtr == ark_mem->maxncf) ||
      (SUNRabs(ark_mem->h) <= ark_mem->hmin * ONEPSM))
  {
    if (*nflagPtr == CONV_FAIL) { return (ARK_CONV_FAILURE); }
    if (*nflagPtr == RHSFUNC_RECVR) { return (ARK_REPTD_RHSFUNC_ERR); }
  }

  /* Reduce step size due to convergence failure */
  ark_mem->eta = hadapt_mem->etacf;

  /* Signal for Jacobian/preconditioner setup */
  *nflagPtr = PREV_CONV_FAIL;

  /* Return to reattempt the step */
  return (PREDICT_AGAIN);
}

/*---------------------------------------------------------------
  arkCheckConstraints

  This routine determines if the constraints of the problem
  are satisfied by the proposed step

  Returns ARK_SUCCESS if successful, otherwise CONSTR_RECVR
  --------------------------------------------------------------*/
int arkCheckConstraints(ARKodeMem ark_mem, int* constrfails, int* nflag)
{
  SUNLogInfo(ARK_LOGGER, "begin-constraint-check", "");

  sunbooleantype constraintsPassed;
  N_Vector mm  = ark_mem->tempv4;
  N_Vector tmp = ark_mem->tempv3;

  /* Check constraints and get mask vector mm for where constraints failed */
  constraintsPassed = N_VConstrMask(ark_mem->constraints, ark_mem->ycur, mm);
  if (constraintsPassed)
  {
    SUNLogInfo(ARK_LOGGER, "end-constraint-check", "status = success");
    return (ARK_SUCCESS);
  }

  /* Constraints not met */

  /* Update total fails and fails in current step */
  ark_mem->nconstrfails++;
  (*constrfails)++;

  /* Return with error if reached max fails in a step */
  if (*constrfails == ark_mem->maxconstrfails)
  {
    SUNLogInfo(ARK_LOGGER, "end-constraint-check",
               "status = failed max attempts");
    return (ARK_CONSTR_FAIL);
  }

  /* Return with error if using fixed step sizes */
  if (ark_mem->fixedstep)
  {
    SUNLogInfo(ARK_LOGGER, "end-constraint-check", "status = failed fixed step");
    return (ARK_CONSTR_FAIL);
  }

  /* Return with error if |h| == hmin */
  if (SUNRabs(ark_mem->h) <= ark_mem->hmin * ONEPSM)
  {
    SUNLogInfo(ARK_LOGGER, "end-constraint-check", "status = failed min step");
    return (ARK_CONSTR_FAIL);
  }

  /* Reduce h by computing eta = h'/h */
  N_VLinearSum(ONE, ark_mem->yn, -ONE, ark_mem->ycur, tmp);
  N_VProd(mm, tmp, tmp);
  ark_mem->eta = SUN_RCONST(0.9) * N_VMinQuotient(ark_mem->yn, tmp);
  ark_mem->eta = SUNMAX(ark_mem->eta, TENTH);

  /* Signal for Jacobian/preconditioner setup */
  *nflag = PREV_CONV_FAIL;

  SUNLogInfo(ARK_LOGGER, "end-constraint-check",
             "status = failed, eta = " SUN_FORMAT_G, ark_mem->eta);

  /* Return to reattempt the step */
  return (CONSTR_RECVR);
}

/*---------------------------------------------------------------
  arkCheckTemporalError

  This routine performs the local error test for the method.
  The weighted local error norm dsm is passed in.  This value is
  used to predict the next step to attempt based on dsm.
  The test dsm <= 1 is made, and if this fails then additional
  checks are performed based on the number of successive error
  test failures.

  Returns ARK_SUCCESS if the test passes.

  If the test fails:
    - if maxnef error test failures have occurred or if
      SUNRabs(h) = hmin, we return ARK_ERR_FAILURE.
    - otherwise: set *nflagPtr to PREV_ERR_FAIL, and
      return TRY_AGAIN.
  --------------------------------------------------------------*/
int arkCheckTemporalError(ARKodeMem ark_mem, int* nflagPtr, int* nefPtr,
                          sunrealtype dsm)
{
  int retval;
  sunrealtype ttmp;
  ARKodeHAdaptMem hadapt_mem;

  /* Access hadapt_mem structure */
  if (ark_mem->hadapt_mem == NULL)
  {
    arkProcessError(ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARKADAPT_NO_MEM);
    return (ARK_MEM_NULL);
  }
  hadapt_mem = ark_mem->hadapt_mem;

  /* consider change of step size for next step attempt (may be
     larger/smaller than current step, depending on dsm) */
  ttmp   = (dsm <= ONE) ? ark_mem->tn + ark_mem->h : ark_mem->tn;
  retval = arkAdapt(ark_mem, hadapt_mem, ark_mem->ycur, ttmp, ark_mem->h, dsm);
  if (retval != ARK_SUCCESS) { return (ARK_ERR_FAILURE); }

  /* if we've made it here then no nonrecoverable failures occurred; someone above
     has recommended an 'eta' value for the next step -- enforce bounds on that value
     and set upcoming step size */
  ark_mem->eta = SUNMIN(ark_mem->eta, ark_mem->hadapt_mem->etamax);
  ark_mem->eta = SUNMAX(ark_mem->eta, ark_mem->hmin / SUNRabs(ark_mem->h));
  ark_mem->eta /=
    SUNMAX(ONE, SUNRabs(ark_mem->h) * ark_mem->hmax_inv * ark_mem->eta);

  /* If est. local error norm dsm passes test, return ARK_SUCCESS */
  if (dsm <= ONE) { return (ARK_SUCCESS); }

  /* Test failed; increment counters, set nflag */
  (*nefPtr)++;
  ark_mem->netf++;
  *nflagPtr = PREV_ERR_FAIL;

  /* At maxnef failures, return ARK_ERR_FAILURE */
  if (*nefPtr == ark_mem->maxnef) { return (ARK_ERR_FAILURE); }

  /* Set etamax=1 to prevent step size increase at end of this step */
  hadapt_mem->etamax = ONE;

  /* Enforce failure bounds on eta */
  if (*nefPtr >= hadapt_mem->small_nef)
  {
    ark_mem->eta = SUNMIN(ark_mem->eta, hadapt_mem->etamxf);
  }

  /* Enforce min/max step bounds once again due to adjustments above */
  ark_mem->eta = SUNMIN(ark_mem->eta, ark_mem->hadapt_mem->etamax);
  ark_mem->eta = SUNMAX(ark_mem->eta, ark_mem->hmin / SUNRabs(ark_mem->h));
  ark_mem->eta /=
    SUNMAX(ONE, SUNRabs(ark_mem->h) * ark_mem->hmax_inv * ark_mem->eta);

  return (TRY_AGAIN);
}

/*---------------------------------------------------------------
  arkAllocVec and arkAllocVecArray:

  These routines allocate (respectively) single vector or a vector
  array based on a template vector.  If the target vector or vector
  array already exists it is left alone; otherwise it is allocated
  by cloning the input vector.

  This routine also updates the optional outputs lrw and liw, which
  are (respectively) the lengths of the overall ARKODE real and
  integer work spaces.

  SUNTRUE is returned if the allocation is successful (or if the
  target vector or vector array already exists) otherwise SUNFALSE
  is returned.
  ---------------------------------------------------------------*/
sunbooleantype arkAllocVec(ARKodeMem ark_mem, N_Vector tmpl, N_Vector* v)
{
  /* return failure if N_VClone or N_VDestroy is not implemented */
  if ((!tmpl->ops->nvclone) || (!tmpl->ops->nvdestroy)) { return SUNFALSE; }

  /* allocate the new vector if necessary */
  if (*v == NULL)
  {
    *v = N_VClone(tmpl);
    if (*v == NULL)
    {
      arkFreeVectors(ark_mem);
      return (SUNFALSE);
    }
    else
    {
      ark_mem->lrw += ark_mem->lrw1;
      ark_mem->liw += ark_mem->liw1;
    }
  }
  return (SUNTRUE);
}

sunbooleantype arkAllocVecArray(int count, N_Vector tmpl, N_Vector** v,
                                sunindextype lrw1, long int* lrw,
                                sunindextype liw1, long int* liw)
{
  /* allocate the new vector array if necessary */
  if (*v == NULL)
  {
    *v = N_VCloneVectorArray(count, tmpl);
    if (*v == NULL) { return (SUNFALSE); }
    *lrw += count * lrw1;
    *liw += count * liw1;
  }
  return (SUNTRUE);
}

/*---------------------------------------------------------------
  arkFreeVec and arkFreeVecArray:

  These routines (respectively) free a single vector or a vector
  array. If the target vector or vector array is already NULL it
  is left alone; otherwise it is freed and the optional outputs
  lrw and liw are updated accordingly.
  ---------------------------------------------------------------*/
void arkFreeVec(ARKodeMem ark_mem, N_Vector* v)
{
  if (*v != NULL)
  {
    N_VDestroy(*v);
    *v = NULL;
    ark_mem->lrw -= ark_mem->lrw1;
    ark_mem->liw -= ark_mem->liw1;
  }
}

void arkFreeVecArray(int count, N_Vector** v, sunindextype lrw1, long int* lrw,
                     sunindextype liw1, long int* liw)
{
  if (*v != NULL)
  {
    N_VDestroyVectorArray(*v, count);
    *v = NULL;
    *lrw -= count * lrw1;
    *liw -= count * liw1;
  }
}

/*---------------------------------------------------------------
  arkResizeVec and arkResizeVecArray:

  This routines (respectively) resize a single vector or a vector
  array based on a template vector. If the ARKVecResizeFn function
  is non-NULL, then it calls that routine to perform the resize;
  otherwise it deallocates and reallocates the target vector or
  vector array based on the template vector. These routines also
  updates the optional outputs lrw and liw, which are
  (respectively) the lengths of the overall ARKODE real and
  integer work spaces.

  SUNTRUE is returned if the resize is successful otherwise
  SUNFALSE is returned.
  ---------------------------------------------------------------*/
sunbooleantype arkResizeVec(ARKodeMem ark_mem, ARKVecResizeFn resize,
                            void* resize_data, sunindextype lrw_diff,
                            sunindextype liw_diff, N_Vector tmpl, N_Vector* v)
{
  if (*v != NULL)
  {
    if (resize == NULL)
    {
      N_VDestroy(*v);
      *v = NULL;
      *v = N_VClone(tmpl);
      if (*v == NULL)
      {
        arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                        "Unable to clone vector");
        return (SUNFALSE);
      }
    }
    else
    {
      if (resize(*v, tmpl, resize_data))
      {
        arkProcessError(ark_mem, ARK_MEM_FAIL, __LINE__, __func__, __FILE__,
                        MSG_ARK_RESIZE_FAIL);
        return (SUNFALSE);
      }
    }
    ark_mem->lrw += lrw_diff;
    ark_mem->liw += liw_diff;
  }
  return (SUNTRUE);
}

sunbooleantype arkResizeVecArray(ARKVecResizeFn resize, void* resize_data,
                                 int count, N_Vector tmpl, N_Vector** v,
                                 sunindextype lrw_diff, long int* lrw,
                                 sunindextype liw_diff, long int* liw)
{
  int i;

  if (*v != NULL)
  {
    if (resize == NULL)
    {
      N_VDestroyVectorArray(*v, count);
      *v = NULL;
      *v = N_VCloneVectorArray(count, tmpl);
      if (*v == NULL) { return (SUNFALSE); }
    }
    else
    {
      for (i = 0; i < count; i++)
      {
        if (resize((*v)[i], tmpl, resize_data)) { return (SUNFALSE); }
      }
    }
    *lrw += count * lrw_diff;
    *liw += count * liw_diff;
  }
  return (SUNTRUE);
}

/*---------------------------------------------------------------
  arkAllocVectors:

  This routine allocates the ARKODE vectors ewt, yn, tempv* and
  ftemp. If any of these vectors already exist, they are left
  alone. Otherwise, it will allocate each vector by cloning the
  input vector. This routine also updates the optional outputs
  lrw and liw, which are (respectively) the lengths of the real
  and integer work spaces.

  If all memory allocations are successful, arkAllocVectors
  returns SUNTRUE, otherwise it returns SUNFALSE.
  ---------------------------------------------------------------*/
sunbooleantype arkAllocVectors(ARKodeMem ark_mem, N_Vector tmpl)
{
  /* Allocate ewt if needed */
  if (!arkAllocVec(ark_mem, tmpl, &ark_mem->ewt)) { return (SUNFALSE); }

  /* Set rwt to point at ewt */
  if (ark_mem->rwt_is_ewt) { ark_mem->rwt = ark_mem->ewt; }

  /* Allocate yn if needed */
  if (!arkAllocVec(ark_mem, tmpl, &ark_mem->yn)) { return (SUNFALSE); }

  /* Allocate tempv1 if needed */
  if (!arkAllocVec(ark_mem, tmpl, &ark_mem->tempv1)) { return (SUNFALSE); }

  /* Allocate tempv2 if needed */
  if (!arkAllocVec(ark_mem, tmpl, &ark_mem->tempv2)) { return (SUNFALSE); }

  /* Allocate tempv3 if needed */
  if (!arkAllocVec(ark_mem, tmpl, &ark_mem->tempv3)) { return (SUNFALSE); }

  /* Allocate tempv4 if needed */
  if (!arkAllocVec(ark_mem, tmpl, &ark_mem->tempv4)) { return (SUNFALSE); }

  return (SUNTRUE);
}

/*---------------------------------------------------------------
  arkResizeVectors:

  This routine resizes all ARKODE vectors if they exist,
  otherwise they are left alone. If a resize function is provided
  it is called to resize the vectors otherwise the vector is
  freed and a new vector is created by cloning in input vector.
  This routine also updates the optional outputs lrw and liw,
  which are (respectively) the lengths of the real and integer
  work spaces.

  If all memory allocations are successful, arkResizeVectors
  returns SUNTRUE, otherwise it returns SUNFALSE.
  ---------------------------------------------------------------*/
sunbooleantype arkResizeVectors(ARKodeMem ark_mem, ARKVecResizeFn resize,
                                void* resize_data, sunindextype lrw_diff,
                                sunindextype liw_diff, N_Vector tmpl)
{
  /* Vabstol */
  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl,
                    &ark_mem->Vabstol))
  {
    return (SUNFALSE);
  }

  /* VRabstol */
  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl,
                    &ark_mem->VRabstol))
  {
    return (SUNFALSE);
  }

  /* ewt */
  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl,
                    &ark_mem->ewt))
  {
    return (SUNFALSE);
  }

  /* rwt  */
  if (ark_mem->rwt_is_ewt)
  { /* update pointer to ewt */
    ark_mem->rwt = ark_mem->ewt;
  }
  else
  { /* resize if distinct from ewt */
    if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl,
                      &ark_mem->rwt))
    {
      return (SUNFALSE);
    }
  }

  /* yn */
  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl,
                    &ark_mem->yn))
  {
    return (SUNFALSE);
  }

  /* fn */
  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl,
                    &ark_mem->fn))
  {
    return (SUNFALSE);
  }

  /* tempv* */
  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl,
                    &ark_mem->tempv1))
  {
    return (SUNFALSE);
  }

  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl,
                    &ark_mem->tempv2))
  {
    return (SUNFALSE);
  }

  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl,
                    &ark_mem->tempv3))
  {
    return (SUNFALSE);
  }

  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl,
                    &ark_mem->tempv4))
  {
    return (SUNFALSE);
  }

  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl,
                    &ark_mem->tempv5))
  {
    return (SUNFALSE);
  }

  return (SUNTRUE);
}

/*---------------------------------------------------------------
  arkFreeVectors

  This routine frees the ARKODE vectors allocated in both
  arkAllocVectors and arkAllocRKVectors.
  ---------------------------------------------------------------*/
void arkFreeVectors(ARKodeMem ark_mem)
{
  arkFreeVec(ark_mem, &ark_mem->ewt);
  if (!ark_mem->rwt_is_ewt) { arkFreeVec(ark_mem, &ark_mem->rwt); }
  arkFreeVec(ark_mem, &ark_mem->tempv1);
  arkFreeVec(ark_mem, &ark_mem->tempv2);
  arkFreeVec(ark_mem, &ark_mem->tempv3);
  arkFreeVec(ark_mem, &ark_mem->tempv4);
  arkFreeVec(ark_mem, &ark_mem->tempv5);
  arkFreeVec(ark_mem, &ark_mem->yn);
  arkFreeVec(ark_mem, &ark_mem->fn);
  arkFreeVec(ark_mem, &ark_mem->Vabstol);
  arkFreeVec(ark_mem, &ark_mem->constraints);
}

/*---------------------------------------------------------------
  arkAccessHAdaptMem:

  Shortcut routine to unpack ark_mem and hadapt_mem structures from
  void* pointer.  If either is missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int arkAccessHAdaptMem(void* arkode_mem, const char* fname, ARKodeMem* ark_mem,
                       ARKodeHAdaptMem* hadapt_mem)
{
  /* access ARKodeMem structure */
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, fname, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *ark_mem = (ARKodeMem)arkode_mem;
  if ((*ark_mem)->hadapt_mem == NULL)
  {
    arkProcessError(*ark_mem, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARKADAPT_NO_MEM);
    return (ARK_MEM_NULL);
  }
  *hadapt_mem = (ARKodeHAdaptMem)(*ark_mem)->hadapt_mem;
  return (ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkProcessError is a high level error handling function that
  calls the appropriate SUNDIALS error handler
  ---------------------------------------------------------------*/
void arkProcessError(ARKodeMem ark_mem, int error_code, int line,
                     const char* func, const char* file, const char* msgfmt, ...)
{
  /* We initialize the argument pointer variable before each vsnprintf call to avoid undefined behavior
     (msgfmt is the last required argument to arkProcessError) */
  va_list ap;

  /* Compose the message */
  va_start(ap, msgfmt);
  size_t msglen = 1;
  if (msgfmt) { msglen += vsnprintf(NULL, 0, msgfmt, ap); }
  va_end(ap);

  char* msg = (char*)malloc(msglen);

  va_start(ap, msgfmt);
  vsnprintf(msg, msglen, msgfmt, ap);
  va_end(ap);

  do {
    if (ark_mem == NULL)
    {
      SUNGlobalFallbackErrHandler(line, func, file, msg, error_code);
      break;
    }

    if (error_code == ARK_WARNING)
    {
#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_WARNING
      char* file_and_line = sunCombineFileAndLine(line, file);
      SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_WARNING, file_and_line, func,
                         msg);
      free(file_and_line);
#endif
      break;
    }

    /* Call the SUNDIALS main error handler */
    SUNHandleErrWithMsg(line, func, file, msg, error_code, ark_mem->sunctx);

    /* Clear the error now */
    (void)SUNContext_GetLastError(ark_mem->sunctx);
  }
  while (0);

  free(msg);

  return;
}

/*---------------------------------------------------------------
  Utility routines for ARKODE to serve as an MRIStepInnerStepper
  ---------------------------------------------------------------*/

/*------------------------------------------------------------------------------
  ark_MRIStepInnerEvolve

  Implementation of MRIStepInnerStepperEvolveFn to advance the inner (fast)
  ODE IVP.  Since the raw return value from an MRIStepInnerStepper is
  meaningless, aside from whether it is 0 (success), >0 (recoverable failure),
  and <0 (unrecoverable failure), we map various ARKODE return values
  accordingly.
  ----------------------------------------------------------------------------*/

int ark_MRIStepInnerEvolve(MRIStepInnerStepper stepper,
                           SUNDIALS_MAYBE_UNUSED sunrealtype t0,
                           sunrealtype tout, N_Vector y)
{
  void* arkode_mem; /* arkode memory             */
  ARKodeMem ark_mem;
  sunrealtype tret;           /* return time               */
  sunrealtype tshift, tscale; /* time normalization values */
  N_Vector* forcing;          /* forcing vectors           */
  int nforcing;               /* number of forcing vectors */
  int retval;                 /* return value              */

  /* extract the ARKODE memory struct */
  retval = MRIStepInnerStepper_GetContent(stepper, &arkode_mem);
  if (retval != ARK_SUCCESS) { return -1; }
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return -1;
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* get the forcing data */
  retval = MRIStepInnerStepper_GetForcingData(stepper, &tshift, &tscale,
                                              &forcing, &nforcing);
  if (retval != ARK_SUCCESS) { return -1; }

  /* set the inner forcing data */
  retval = ark_mem->step_setforcing(ark_mem, tshift, tscale, forcing, nforcing);
  if (retval != ARK_SUCCESS) { return -1; }

  /* set the stop time */
  retval = ARKodeSetStopTime(arkode_mem, tout);
  if (retval != ARK_SUCCESS) { return -1; }

  /* evolve inner ODE, consider all positive return values as 'success' */
  retval = ARKodeEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);
  if (retval > 0) { retval = 0; }

  /* set a recoverable failure for a few ARKODE failure modes;
     on other ARKODE errors return with an unrecoverable failure */
  if (retval < 0)
  {
    if ((retval == ARK_TOO_MUCH_WORK) || (retval == ARK_CONV_FAILURE) ||
        (retval == ARK_ERR_FAILURE))
    {
      retval = 1;
    }
    else { return -1; }
  }

  /* disable inner forcing */
  if (ark_mem->step_setforcing(ark_mem, ZERO, ONE, NULL, 0) != ARK_SUCCESS)
  {
    return -1;
  }

  return retval;
}

/*------------------------------------------------------------------------------
  ark_MRIStepInnerFullRhs

  Implementation of MRIStepInnerStepperFullRhsFn to compute the full inner
  (fast) ODE IVP RHS.
  ----------------------------------------------------------------------------*/

int ark_MRIStepInnerFullRhs(MRIStepInnerStepper stepper, sunrealtype t,
                            N_Vector y, N_Vector f, int mode)
{
  void* arkode_mem; /* arkode memory */
  ARKodeMem ark_mem;
  int retval = MRIStepInnerStepper_GetContent(stepper, &arkode_mem);
  if (retval != ARK_SUCCESS) { return -1; }
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return -1;
  }
  ark_mem = (ARKodeMem)arkode_mem;
  retval  = ark_mem->step_fullrhs(arkode_mem, t, y, f, mode);
  if (retval == ARK_SUCCESS) { return 0; }
  return -1;
}

/*------------------------------------------------------------------------------
  ark_MRIStepInnerReset

  Implementation of MRIStepInnerStepperResetFn to reset the inner (fast) stepper
  state.
  ----------------------------------------------------------------------------*/

int ark_MRIStepInnerReset(MRIStepInnerStepper stepper, sunrealtype tR, N_Vector yR)
{
  void* arkode_mem;
  int retval = MRIStepInnerStepper_GetContent(stepper, &arkode_mem);
  if (retval != ARK_SUCCESS) { return -1; }
  retval = ARKodeReset(arkode_mem, tR, yR);
  if (retval == ARK_SUCCESS) { return 0; }
  return -1;
}

/*------------------------------------------------------------------------------
  ark_MRIStepInnerGetAccumulatedError

  Implementation of MRIStepInnerGetAccumulatedError to retrieve the accumulated
  temporal error estimate from the inner (fast) stepper.
  ----------------------------------------------------------------------------*/

int ark_MRIStepInnerGetAccumulatedError(MRIStepInnerStepper stepper,
                                        sunrealtype* accum_error)
{
  void* arkode_mem;
  int retval = MRIStepInnerStepper_GetContent(stepper, &arkode_mem);
  if (retval != ARK_SUCCESS) { return -1; }
  retval = ARKodeGetAccumulatedError(arkode_mem, accum_error);
  if (retval == ARK_SUCCESS) { return 0; }
  if (retval > 0) { return 1; }
  return -1;
}

/*------------------------------------------------------------------------------
  ark_MRIStepInnerResetAccumulatedError

  Implementation of MRIStepInnerResetAccumulatedError to reset the accumulated
  temporal error estimator in the inner (fast) stepper.
  ----------------------------------------------------------------------------*/

int ark_MRIStepInnerResetAccumulatedError(MRIStepInnerStepper stepper)
{
  void* arkode_mem;
  int retval = MRIStepInnerStepper_GetContent(stepper, &arkode_mem);
  if (retval != ARK_SUCCESS) { return -1; }
  retval = ARKodeResetAccumulatedError(arkode_mem);
  if (retval == ARK_SUCCESS) { return 0; }
  return -1;
}

/*------------------------------------------------------------------------------
  ark_MRIStepInnerSetRTol

  Implementation of MRIStepInnerSetRTol to set a relative tolerance for the
  upcoming evolution using the inner (fast) stepper.
  ----------------------------------------------------------------------------*/

int ark_MRIStepInnerSetRTol(MRIStepInnerStepper stepper, sunrealtype rtol)
{
  void* arkode_mem;
  ARKodeMem ark_mem;
  int retval = MRIStepInnerStepper_GetContent(stepper, &arkode_mem);
  if (retval != ARK_SUCCESS) { return -1; }
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return -1;
  }
  ark_mem = (ARKodeMem)arkode_mem;
  if (rtol > ZERO)
  {
    ark_mem->reltol = rtol;
    return 0;
  }
  else { return -1; }
}

/*===============================================================
  EOF
  ===============================================================*/
