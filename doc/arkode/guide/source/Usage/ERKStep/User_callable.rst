.. ----------------------------------------------------------------
   Programmer(s): Daniel R. Reynolds @ UMBC
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025-2026, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage.ERKStep.UserCallable:

ERKStep User-callable functions
==================================

This section describes the ERKStep-specific functions that may be called
by the user to setup and then solve an IVP using the ERKStep time-stepping
module.  The large majority of these routines merely wrap :ref:`underlying
ARKODE functions <ARKODE.Usage.UserCallable>`, and are now deprecated
-- each of these are clearly marked.  However, some
of these user-callable functions are specific to ERKStep, as explained
below.

As discussed in the main :ref:`ARKODE user-callable function introduction
<ARKODE.Usage.UserCallable>`, each of ARKODE's time-stepping modules
clarifies the categories of user-callable functions that it supports.
ERKStep supports the following categories:

* temporal adaptivity
* relaxation Runge--Kutta methods

ERKStep also has forcing function support when converted to a
:c:type:`SUNStepper` or :c:type:`MRIStepInnerStepper`. See
:c:func:`ARKodeCreateSUNStepper` and :c:func:`ARKStepCreateMRIStepInnerStepper`
for additional details.


.. _ARKODE.Usage.ERKStep.Initialization:

ERKStep initialization and deallocation functions
------------------------------------------------------


.. c:function:: void* ERKStepCreate(ARKRhsFn f, sunrealtype t0, N_Vector y0, SUNContext sunctx)

   This function allocates and initializes memory for a problem to
   be solved using the ERKStep time-stepping module in ARKODE.

   **Arguments:**
      * *f* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the right-hand side function in
        :math:`\dot{y} = f(t,y)`.
      * *t0* -- the initial value of :math:`t`.
      * *y0* -- the initial condition vector :math:`y(t_0)`.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a pointer to initialized problem memory
      of type ``void*``, to be passed to all user-facing ERKStep routines
      listed below.  If unsuccessful, a ``NULL`` pointer will be
      returned, and an error message will be printed to ``stderr``.


.. c:function:: void ERKStepFree(void** arkode_mem)

   This function frees the problem memory *arkode_mem* created by
   :c:func:`ERKStepCreate`.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.

   **Return value:**  None

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeFree` instead.



.. _ARKODE.Usage.ERKStep.Tolerances:

ERKStep tolerance specification functions
------------------------------------------------------

.. c:function:: int ERKStepSStolerances(void* arkode_mem, sunrealtype reltol, sunrealtype abstol)

   This function specifies scalar relative and absolute tolerances.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *reltol* -- scalar relative tolerance.
      * *abstol* -- scalar absolute tolerance.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ERKStep memory was ``NULL``
      * *ARK_NO_MALLOC*  if the ERKStep memory was not allocated by the time-stepping module
      * *ARK_ILL_INPUT* if an argument had an illegal value (e.g. a negative tolerance).

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSStolerances` instead.



.. c:function:: int ERKStepSVtolerances(void* arkode_mem, sunrealtype reltol, N_Vector abstol)

   This function specifies a scalar relative tolerance and a vector
   absolute tolerance (a potentially different absolute tolerance for
   each vector component).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *reltol* -- scalar relative tolerance.
      * *abstol* -- vector containing the absolute tolerances for each
        solution component.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ERKStep memory was ``NULL``
      * *ARK_NO_MALLOC*  if the ERKStep memory was not allocated by the time-stepping module
      * *ARK_ILL_INPUT* if an argument had an illegal value (e.g. a negative tolerance).

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSVtolerances` instead.



.. c:function:: int ERKStepWFtolerances(void* arkode_mem, ARKEwtFn efun)

   This function specifies a user-supplied function *efun* to compute
   the error weight vector ``ewt``.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *efun* -- the name of the function (of type :c:func:`ARKEwtFn()`)
        that implements the error weight vector computation.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ERKStep memory was ``NULL``
      * *ARK_NO_MALLOC*  if the ERKStep memory was not allocated by the time-stepping module

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeWFtolerances` instead.



.. _ARKODE.Usage.ERKStep.RootFinding:

Rootfinding initialization function
--------------------------------------

.. c:function:: int ERKStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g)

   Initializes a rootfinding problem to be solved during the
   integration of the ODE system.  It must be called after
   :c:func:`ERKStepCreate`, and before :c:func:`ERKStepEvolve()`.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *nrtfn* -- number of functions :math:`g_i`, an integer :math:`\ge` 0.
      * *g* -- name of user-supplied function, of type :c:func:`ARKRootFn()`,
        defining the functions :math:`g_i` whose roots are sought.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ERKStep memory was ``NULL``
      * *ARK_MEM_FAIL*  if there was a memory allocation failure
      * *ARK_ILL_INPUT* if *nrtfn* is greater than zero but *g* = ``NULL``.

   **Notes:**
      To disable the rootfinding feature after it has already
      been initialized, or to free memory associated with ERKStep's
      rootfinding module, call *ERKStepRootInit* with *nrtfn = 0*.

      Similarly, if a new IVP is to be solved with a call to
      :c:func:`ERKStepReInit()`, where the new IVP has no rootfinding
      problem but the prior one did, then call *ERKStepRootInit* with
      *nrtfn = 0*.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeRootInit` instead.




.. _ARKODE.Usage.ERKStep.Integration:

ERKStep solver function
-------------------------

.. c:function:: int ERKStepEvolve(void* arkode_mem, sunrealtype tout, N_Vector yout, sunrealtype *tret, int itask)

   Integrates the ODE over an interval in :math:`t`.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *tout* -- the next time at which a computed solution is desired.
      * *yout* -- the computed solution vector.
      * *tret* -- the time corresponding to *yout* (output).
      * *itask* -- a flag indicating the job of the solver for the next
        user step.

        The *ARK_NORMAL* option causes the solver to take internal
        steps until it has just overtaken a user-specified output
        time, *tout*, in the direction of integration,
        i.e. :math:`t_{n-1} <` *tout* :math:`\le t_{n}` for forward
        integration, or :math:`t_{n} \le` *tout* :math:`< t_{n-1}` for
        backward integration.  It will then compute an approximation
        to the solution :math:`y(tout)` by interpolation (using one
        of the dense output routines described in
        :numref:`ARKODE.Mathematics.Interpolation`).

        The *ARK_ONE_STEP* option tells the solver to only take a
        single internal step, :math:`y_{n-1} \to y_{n}`, and return the solution
        at that point, :math:`y_{n}`, in the vector *yout*.

   **Return value:**
      * *ARK_SUCCESS* if successful.
      * *ARK_ROOT_RETURN* if :c:func:`ERKStepEvolve()` succeeded, and
        found one or more roots.  If the number of root functions,
        *nrtfn*, is greater than 1, call
        :c:func:`ERKStepGetRootInfo()` to see which :math:`g_i` were
        found to have a root at (*\*tret*).
      * *ARK_TSTOP_RETURN* if :c:func:`ERKStepEvolve()` succeeded and
        returned at *tstop*.
      * *ARK_MEM_NULL* if the *arkode_mem* argument was ``NULL``.
      * *ARK_NO_MALLOC* if *arkode_mem* was not allocated.
      * *ARK_ILL_INPUT* if one of the inputs to
        :c:func:`ERKStepEvolve()` is illegal, or some other input to
        the solver was either illegal or missing.  Details will be
        provided in the error message.  Typical causes of this failure:

        (a) A component of the error weight vector became zero during
            internal time-stepping.

        (b) A root of one of the root functions was found both at a
            point :math:`t` and also very near :math:`t`.

        (c) The initial condition violates the inequality constraints.

      * *ARK_TOO_MUCH_WORK* if the solver took *mxstep* internal steps
        but could not reach *tout*.  The default value for *mxstep* is
        *MXSTEP_DEFAULT = 500*.
      * *ARK_TOO_MUCH_ACC* if the solver could not satisfy the accuracy
        demanded by the user for some internal step.
      * *ARK_ERR_FAILURE* if error test failures occurred either too many
        times (*ark_maxnef*) during one internal time step or occurred
        with :math:`|h| = h_{min}`.
      * *ARK_VECTOROP_ERR* a vector operation error occurred.

   **Notes:**
      The input vector *yout* can use the same memory as the
      vector *y0* of initial conditions that was passed to
      :c:func:`ERKStepCreate`.

      In *ARK_ONE_STEP* mode, *tout* is used only on the first call, and
      only to get the direction and a rough scale of the independent
      variable.

      All failure return values are negative and so testing the
      return argument for negative values will trap all
      :c:func:`ERKStepEvolve()` failures.

      Since interpolation may reduce the accuracy in the reported
      solution, if full method accuracy is desired the user should issue
      a call to :c:func:`ERKStepSetStopTime()` before the call to
      :c:func:`ERKStepEvolve()` to specify a fixed stop time to
      end the time step and return to the user.  Upon return from
      :c:func:`ERKStepEvolve()`, a copy of the internal solution
      :math:`y_{n}` will be returned in the vector *yout*.  Once the
      integrator returns at a *tstop* time, any future testing for
      *tstop* is disabled (and can be re-enabled only though a new call
      to :c:func:`ERKStepSetStopTime()`).

      On any error return in which one or more internal steps were taken
      by :c:func:`ERKStepEvolve()`, the returned values of *tret* and
      *yout* correspond to the farthest point reached in the integration.
      On all other error returns, *tret* and *yout* are left unchanged
      from those provided to the routine.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeEvolve` instead.




.. _ARKODE.Usage.ERKStep.OptionalInputs:

Optional input functions
-------------------------


.. _ARKODE.Usage.ERKStep.ERKStepInput:

Optional inputs for ERKStep
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: int ERKStepSetDefaults(void* arkode_mem)

   Resets all optional input parameters to ERKStep's original
   default values.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      Does not change problem-defining function pointer *f*
      or the *user_data* pointer.

      Also leaves alone any data structures or options related to
      root-finding (those can be reset using :c:func:`ERKStepRootInit()`).

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetDefaults` instead.



.. c:function:: int ERKStepSetInterpolantType(void* arkode_mem, int itype)

   .. deprecated:: 6.1.0

      This function is now a wrapper to :c:func:`ARKodeSetInterpolantType`, see
      the documentation for that function instead.



.. c:function:: int ERKStepSetInterpolantDegree(void* arkode_mem, int degree)

   Specifies the degree of the polynomial interpolant
   used for dense output (i.e. interpolation of solution output values
   and implicit method predictors).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *degree* -- requested polynomial degree.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory or interpolation module are ``NULL``
      * *ARK_INTERP_FAIL* if this is called after :c:func:`ERKStepEvolve()`
      * *ARK_ILL_INPUT* if an argument had an illegal value or the
        interpolation module has already been initialized

   **Notes:**
      Allowed values are between 0 and 5.

      This routine should be called *after* :c:func:`ERKStepCreate` and *before*
      :c:func:`ERKStepEvolve()`. After the first call to :c:func:`ERKStepEvolve()`
      the interpolation degree may not be changed without first calling
      :c:func:`ERKStepReInit()`.

      If a user calls both this routine and :c:func:`ERKStepSetInterpolantType()`, then
      :c:func:`ERKStepSetInterpolantType()` must be called first.

      Since the accuracy of any polynomial interpolant is limited by the
      accuracy of the time-step solutions on which it is based, the *actual*
      polynomial degree that is used by ERKStep will be the minimum of
      :math:`q-1` and the input *degree*, for :math:`q > 1` where :math:`q` is
      the order of accuracy for the time integration method.

      .. versionchanged:: 5.5.1

         When :math:`q=1`, a linear interpolant is the default to ensure values
         obtained by the integrator are returned at the ends of the time
         interval.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetInterpolantDegree` instead.



.. c:function:: int ERKStepSetDenseOrder(void* arkode_mem, int dord)

   .. deprecated:: 5.2.0

      Use :c:func:`ARKodeSetInterpolantDegree` instead.



.. c:function:: int ERKStepSetDiagnostics(void* arkode_mem, FILE* diagfp)

   Specifies the file pointer for a diagnostics file where
   all ERKStep step adaptivity and solver information is written.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *diagfp* -- pointer to the diagnostics output file.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      This parameter can be ``stdout`` or ``stderr``, although the
      suggested approach is to specify a pointer to a unique file opened
      by the user and returned by ``fopen``.  If not called, or if called
      with a ``NULL`` file pointer, all diagnostics output is disabled.

      When run in parallel, only one process should set a non-NULL value
      for this pointer, since statistics from all processes would be
      identical.

   .. deprecated:: 5.2.0

      Use :c:func:`SUNLogger_SetInfoFilename` instead.


.. c:function:: int ERKStepSetFixedStep(void* arkode_mem, sunrealtype hfixed)

   Disabled time step adaptivity within ERKStep, and specifies the
   fixed time step size to use for the following internal step(s).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *hfixed* -- value of the fixed step size to use.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      Pass 0.0 to return ERKStep to the default (adaptive-step) mode.

      Use of this function is not generally recommended, since we it gives no
      assurance of the validity of the computed solutions.  It is
      primarily provided for code-to-code verification testing purposes.

      When using :c:func:`ERKStepSetFixedStep()`, any values provided to
      the functions
      :c:func:`ERKStepSetInitStep()`,
      :c:func:`ERKStepSetAdaptivityFn()`,
      :c:func:`ERKStepSetMaxErrTestFails()`,
      :c:func:`ERKStepSetAdaptivityMethod()`,
      :c:func:`ERKStepSetCFLFraction()`,
      :c:func:`ERKStepSetErrorBias()`,
      :c:func:`ERKStepSetFixedStepBounds()`,
      :c:func:`ERKStepSetMaxEFailGrowth()`,
      :c:func:`ERKStepSetMaxFirstGrowth()`,
      :c:func:`ERKStepSetMaxGrowth()`,
      :c:func:`ERKStepSetMinReduction()`,
      :c:func:`ERKStepSetSafetyFactor()`,
      :c:func:`ERKStepSetSmallNumEFails()`,
      :c:func:`ERKStepSetStabilityFn()`, and
      :c:func:`ERKStepSetAdaptController()`
      will be ignored, since temporal adaptivity is disabled.

      If both :c:func:`ERKStepSetFixedStep()` and
      :c:func:`ERKStepSetStopTime()` are used, then the fixed step size
      will be used for all steps until the final step preceding the
      provided stop time (which may be shorter).  To resume use of the
      previous fixed step size, another call to
      :c:func:`ERKStepSetFixedStep()` must be made prior to calling
      :c:func:`ERKStepEvolve()` to resume integration.

      It is *not* recommended that :c:func:`ERKStepSetFixedStep()` be used
      in concert with :c:func:`ERKStepSetMaxStep()` or
      :c:func:`ERKStepSetMinStep()`, since at best those latter two
      routines will provide no useful information to the solver, and at
      worst they may interfere with the desired fixed step size.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetFixedStep` instead.



.. c:function:: int ERKStepSetInitStep(void* arkode_mem, sunrealtype hin)

   Specifies the initial time step size ERKStep should use after
   initialization, re-initialization, or resetting.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *hin* -- value of the initial step to be attempted :math:`(\ne 0)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      Pass 0.0 to use the default value.

      By default, ERKStep estimates the initial step size to be
      :math:`h = \sqrt{\dfrac{2}{\left\| \ddot{y} \right\|}}`, where
      :math:`\ddot{y}` is an estimate of the second derivative of the
      solution at :math:`t_0`.

      This routine will also reset the step size and error history.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetInitStep` instead.



.. c:function:: int ERKStepSetMaxHnilWarns(void* arkode_mem, int mxhnil)

   Specifies the maximum number of messages issued by the
   solver to warn that :math:`t+h=t` on the next internal step, before
   ERKStep will instead return with an error.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *mxhnil* -- maximum allowed number of warning messages :math:`(>0)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      The default value is 10; set *mxhnil* to zero to specify
      this default.

      A negative value indicates that no warning messages should be issued.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetMaxHnilWarns` instead.



.. c:function:: int ERKStepSetMaxNumSteps(void* arkode_mem, long int mxsteps)

   Specifies the maximum number of steps to be taken by the
   solver in its attempt to reach the next output time, before ERKStep
   will return with an error.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *mxsteps* -- maximum allowed number of internal steps.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      Passing *mxsteps* = 0 results in ERKStep using the
      default value (500).

      Passing *mxsteps* < 0 disables the test (not recommended).

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetMaxNumSteps` instead.



.. c:function:: int ERKStepSetMaxStep(void* arkode_mem, sunrealtype hmax)

   Specifies the upper bound on the magnitude of the time step size.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *hmax* -- maximum absolute value of the time step size :math:`(\ge 0)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      Pass *hmax* :math:`\le 0.0` to set the default value of :math:`\infty`.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetMaxStep` instead.



.. c:function:: int ERKStepSetMinStep(void* arkode_mem, sunrealtype hmin)

   Specifies the lower bound on the magnitude of the time step size.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *hmin* -- minimum absolute value of the time step size :math:`(\ge 0)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      Pass *hmin* :math:`\le 0.0` to set the default value of 0.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetMinStep` instead.



.. c:function:: int ERKStepSetStopTime(void* arkode_mem, sunrealtype tstop)

   Specifies the value of the independent variable
   :math:`t` past which the solution is not to proceed.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *tstop* -- stopping time for the integrator.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      The default is that no stop time is imposed.

      Once the integrator returns at a stop time, any future testing for
      ``tstop`` is disabled (and can be re-enabled only though a new call to
      :c:func:`ERKStepSetStopTime`).

      A stop time not reached before a call to :c:func:`ERKStepReInit` or
      :c:func:`ERKStepReset` will remain active but can be disabled by calling
      :c:func:`ERKStepClearStopTime`.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetStopTime` instead.



.. c:function:: int ERKStepSetInterpolateStopTime(void* arkode_mem, sunbooleantype interp)

   Specifies that the output solution should be interpolated when the current
   :math:`t` equals the specified ``tstop`` (instead of merely copying the
   internal solution :math:`y_n`).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *interp* -- flag indicating to use interpolation (1) or copy (0).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``

   .. versionadded:: 5.6.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetInterpolateStopTime` instead.



.. c:function:: int ERKStepClearStopTime(void* arkode_mem)

   Disables the stop time set with :c:func:`ERKStepSetStopTime`.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``

   **Notes:**
      The stop time can be re-enabled though a new call to
      :c:func:`ERKStepSetStopTime`.

   .. versionadded:: 5.5.1

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeClearStopTime` instead.



.. c:function:: int ERKStepSetUserData(void* arkode_mem, void* user_data)

   Specifies the user data block *user_data* and
   attaches it to the main ERKStep memory block.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *user_data* -- pointer to the user data.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      If specified, the pointer to *user_data* is passed to all
      user-supplied functions for which it is an argument; otherwise
      ``NULL`` is passed.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetUserData` instead.



.. c:function:: int ERKStepSetMaxErrTestFails(void* arkode_mem, int maxnef)

   Specifies the maximum number of error test failures
   permitted in attempting one step, before returning with an error.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *maxnef* -- maximum allowed number of error test failures :math:`(>0)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      The default value is 7; set *maxnef* :math:`\le 0`
      to specify this default.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetMaxErrTestFails` instead.



.. c:function:: int ERKStepSetConstraints(void* arkode_mem, N_Vector constraints)

   Specifies a vector defining inequality constraints for each component of the
   solution vector :math:`y`.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *constraints* -- vector of constraint flags. Each component specifies
        the type of solution constraint:

        .. math::

           \texttt{constraints[i]} = \left\{ \begin{array}{rcl}
           0.0  &\Rightarrow\;& \text{no constraint is imposed on}\; y_i,\\
           1.0  &\Rightarrow\;& y_i \geq 0,\\
           -1.0  &\Rightarrow\;& y_i \leq 0,\\
           2.0  &\Rightarrow\;& y_i > 0,\\
           -2.0  &\Rightarrow\;& y_i < 0.\\
           \end{array}\right.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if the constraints vector contains illegal values

   **Notes:**
      The presence of a non-``NULL`` constraints vector that is not 0.0
      in all components will cause constraint checking to be performed. However, a
      call with 0.0 in all components of ``constraints`` will result in an illegal
      input return. A ``NULL`` constraints vector will disable constraint checking.

      After a call to :c:func:`ERKStepResize()` inequality constraint checking
      will be disabled and a call to :c:func:`ERKStepSetConstraints()` is
      required to re-enable constraint checking.

      Since constraint-handling is performed through cutting time steps that would
      violate the constraints, it is possible that this feature will cause some
      problems to fail due to an inability to enforce constraints even at the
      minimum time step size.  Additionally, the features :c:func:`ERKStepSetConstraints()`
      and :c:func:`ERKStepSetFixedStep()` are incompatible, and should not be used
      simultaneously.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetConstraints` instead.



.. c:function:: int ERKStepSetMaxNumConstrFails(void* arkode_mem, int maxfails)

   Specifies the maximum number of constraint failures in a step before ERKStep
   will return with an error.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *maxfails* -- maximum allowed number of constrain failures.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``

   **Notes:**
      Passing *maxfails* <= 0 results in ERKStep using the
      default value (10).

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetMaxNumConstrFails` instead.



.. _ARKODE.Usage.ERKStep.ERKStepMethodInput:

Optional inputs for IVP method selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _ARKODE.Usage.ERKStep.ERKStepMethodInputTable:
.. table:: Optional inputs for IVP method selection

   +--------------------------------------+---------------------------------+------------------+
   | Optional input                       | Function name                   | Default          |
   +--------------------------------------+---------------------------------+------------------+
   | Set integrator method order          | :c:func:`ERKStepSetOrder()`     | 4                |
   +--------------------------------------+---------------------------------+------------------+
   | Set explicit RK table                | :c:func:`ERKStepSetTable()`     | internal         |
   +--------------------------------------+---------------------------------+------------------+
   | Set explicit RK table via its number | :c:func:`ERKStepSetTableNum()`  | internal         |
   +--------------------------------------+---------------------------------+------------------+
   | Set explicit RK table via its name   | :c:func:`ERKStepSetTableName()` | internal         |
   +--------------------------------------+---------------------------------+------------------+



.. c:function:: int ERKStepSetOrder(void* arkode_mem, int ord)

   Specifies the order of accuracy for the ERK integration method.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *ord* -- requested order of accuracy.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      The allowed values are :math:`2 \le` *ord* :math:`\le
      8`.  Any illegal input will result in the default value of 4.

      Since *ord* affects the memory requirements for the internal
      ERKStep memory block, it cannot be changed after the first call to
      :c:func:`ERKStepEvolve()`, unless :c:func:`ERKStepReInit()` is called.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetOrder` instead.



.. c:function:: int ERKStepSetTable(void* arkode_mem, ARKodeButcherTable B)

   Specifies a customized Butcher table for the ERK method.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *B* -- the Butcher table for the explicit RK method.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**

      For a description of the :c:type:`ARKodeButcherTable` type and related
      functions for creating Butcher tables, see :numref:`ARKodeButcherTable`.

      No error checking is performed to ensure that either the method order *p* or
      the embedding order *q* specified in the Butcher table structure correctly
      describe the coefficients in the Butcher table.

      Error checking is performed to ensure that the Butcher table is strictly
      lower-triangular (i.e. that it specifies an ERK method).

      If the Butcher table does not contain an embedding, the user *must* call
      :c:func:`ERKStepSetFixedStep()` to enable fixed-step mode and set the desired
      time step size.

   **Warning:**
      This should not be used with :c:func:`ARKodeSetOrder`.


.. c:function:: int ERKStepSetTableNum(void* arkode_mem, ARKODE_ERKTableID etable)

   Indicates to use a specific built-in Butcher table for the ERK method.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *etable* -- index of the Butcher table.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      *etable* should match an existing explicit method from
      :numref:`Butcher.explicit`.  Error-checking is performed
      to ensure that the table exists, and is not implicit.

   **Warning:**
      This should not be used with :c:func:`ARKodeSetOrder`.



.. c:function:: int ERKStepSetTableName(void* arkode_mem, const char *etable)

   Indicates to use a specific built-in Butcher table for the ERK method.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *etable* -- name of the Butcher table.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      *etable* should match an existing explicit method from
      :numref:`Butcher.explicit`.  Error-checking is performed
      to ensure that the table exists, and is not implicit.
      This function is case sensitive.

   **Warning:**
      This should not be used with :c:func:`ARKodeSetOrder`.

   .. note::

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.table_name".


.. _ARKODE.Usage.ERKStep.ERKStepAdaptivityInput:

Optional inputs for time step adaptivity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mathematical explanation of ARKODE's time step adaptivity
algorithm, including how each of the parameters below is used within
the code, is provided in :numref:`ARKODE.Mathematics.Adaptivity`.


.. c:function:: int ERKStepSetAdaptController(void* arkode_mem, SUNAdaptController C)

   Sets a user-supplied time-step controller object.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *C* -- user-supplied time adaptivity controller.  If ``NULL`` then the I controller will be created (see :numref:`SUNAdaptController.Soderlind`).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_MEM_FAIL* if *C* was ``NULL`` and the I controller could not be allocated.

   .. versionadded:: 5.7.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetAdaptController` instead.

   .. versionchanged:: 6.3.0

      The default controller was changed from PI to I. Additionally, in prior
      versions, passing ``NULL`` to this function would attach the PID
      controller.



.. c:function:: int ERKStepSetAdaptivityFn(void* arkode_mem, ARKAdaptFn hfun, void* h_data)

   Sets a user-supplied time-step adaptivity function.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *hfun* -- name of user-supplied adaptivity function.
      * *h_data* -- pointer to user data passed to *hfun* every time
        it is called.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      This function should focus on accuracy-based time step
      estimation; for stability based time steps the function
      :c:func:`ERKStepSetStabilityFn()` should be used instead.


   .. deprecated:: 5.7.0

      Use the SUNAdaptController infrastructure instead (see :numref:`SUNAdaptController.Description`).


.. c:function:: int ERKStepSetAdaptivityMethod(void* arkode_mem, int imethod, int idefault, int pq, sunrealtype* adapt_params)

   Specifies the method (and associated parameters) used for time step adaptivity.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *imethod* -- accuracy-based adaptivity method choice
        (0 :math:`\le` `imethod` :math:`\le` 5):
        0 is PID, 1 is PI, 2 is I, 3 is explicit Gustafsson, 4 is
        implicit Gustafsson, and 5 is the ImEx Gustafsson.
      * *idefault* -- flag denoting whether to use default adaptivity
        parameters (1), or that they will be supplied in the
        *adapt_params* argument (0).
      * *pq* -- flag denoting whether to use the embedding order of
        accuracy *p* (0), the method order of accuracy *q* (1), or the
        minimum of the two (any input not equal to 0 or 1)
        within the adaptivity algorithm.  *p* is the default.
      * *adapt_params[0]* -- :math:`k_1` parameter within accuracy-based adaptivity algorithms.
      * *adapt_params[1]* -- :math:`k_2` parameter within accuracy-based adaptivity algorithms.
      * *adapt_params[2]* -- :math:`k_3` parameter within accuracy-based adaptivity algorithms.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      If custom parameters are supplied, they will be checked
      for validity against published stability intervals.  If other
      parameter values are desired, it is recommended to instead provide
      a custom function through a call to :c:func:`ERKStepSetAdaptivityFn()`.

      .. versionchanged:: 5.7.0

         Prior to version 5.7.0, any nonzero value for *pq* would result in use of the
         embedding order of accuracy.

   .. deprecated:: 5.7.0

      Use the SUNAdaptController infrastructure instead (see :numref:`SUNAdaptController.Description`).


.. c:function:: int ERKStepSetAdaptivityAdjustment(void* arkode_mem, int adjust)

   Called by a user to adjust the method order supplied to the temporal adaptivity
   controller.  For example, if the user expects order reduction due to problem stiffness,
   they may request that the controller assume a reduced order of accuracy for the method
   by specifying a value :math:`adjust < 0`.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *adjust* -- adjustment factor (default is 0).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      This should be called prior to calling :c:func:`ERKStepEvolve()`, and can only be
      reset following a call to :c:func:`ERKStepReInit()`.

   .. versionadded:: 5.7.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetAdaptivityAdjustment` instead.

   .. versionchanged:: 6.3.0

      The default value was changed from -1 to 0



.. c:function:: int ERKStepSetCFLFraction(void* arkode_mem, sunrealtype cfl_frac)

   Specifies the fraction of the estimated explicitly stable step to use.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *cfl_frac* -- maximum allowed fraction of explicitly stable step (default is 0.5).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      Any non-positive parameter will imply a reset to the default
      value.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetCFLFraction` instead.



.. c:function:: int ERKStepSetErrorBias(void* arkode_mem, sunrealtype bias)

   Specifies the bias to be applied to the error estimates within
   accuracy-based adaptivity strategies.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *bias* -- bias applied to error in accuracy-based time
        step estimation (default is 1.0).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      Any value below 1.0 will imply a reset to the default value.

      If both this and one of :c:func:`ERKStepSetAdaptivityMethod` or
      :c:func:`ERKStepSetAdaptController` will be called, then this routine must be called
      *second*.

   .. deprecated:: 5.7.0

      Use the SUNAdaptController infrastructure instead (see :numref:`SUNAdaptController.Description`).

   .. versionchanged:: 6.3.0

      The default value was changed from 1.5 to 1.0



.. c:function:: int ERKStepSetFixedStepBounds(void* arkode_mem, sunrealtype lb, sunrealtype ub)

   Specifies the step growth interval in which the step size will remain unchanged.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *lb* -- lower bound on window to leave step size fixed (default is 1.0).
      * *ub* -- upper bound on window to leave step size fixed (default is 1.0).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      Any interval *not* containing 1.0 will imply a reset to the default values.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetFixedStepBounds` instead.

   .. versionchanged:: 6.3.0

      The default upper bound was changed from 1.5 to 1.0



.. c:function:: int ERKStepSetMaxEFailGrowth(void* arkode_mem, sunrealtype etamxf)

   Specifies the maximum step size growth factor upon multiple successive
   accuracy-based error failures in the solver.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *etamxf* -- time step reduction factor on multiple error fails (default is 0.3).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      Any value outside the interval :math:`(0,1]` will imply a reset to the default value.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetMaxEFailGrowth` instead.



.. c:function:: int ERKStepSetMaxFirstGrowth(void* arkode_mem, sunrealtype etamx1)

   Specifies the maximum allowed growth factor in step size following the very
   first integration step.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *etamx1* -- maximum allowed growth factor after the first time
        step (default is 10000.0).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      Any value :math:`\le 1.0` will imply a reset to the default value.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetMaxFirstGrowth` instead.



.. c:function:: int ERKStepSetMaxGrowth(void* arkode_mem, sunrealtype mx_growth)

   Specifies the maximum allowed growth factor in step size between
   consecutive steps in the integration process.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *mx_growth* -- maximum allowed growth factor between consecutive time steps (default is 20.0).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      Any value :math:`\le 1.0` will imply a reset to the default
      value.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetMaxGrowth` instead.



.. c:function:: int ERKStepSetMinReduction(void* arkode_mem, sunrealtype eta_min)

   Specifies the minimum allowed reduction factor in step size between
   step attempts, resulting from a temporal error failure in the integration
   process.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *eta_min* -- minimum allowed reduction factor time step after an error
        test failure (default is 0.1).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      Any value :math:`\ge 1.0` or :math:`\le 0.0` will imply a reset to
      the default value.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetMinReduction` instead.



.. c:function:: int ERKStepSetSafetyFactor(void* arkode_mem, sunrealtype safety)

   Specifies the safety factor to be applied to the accuracy-based
   estimated step.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *safety* -- safety factor applied to accuracy-based time step (default is 0.9).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      Any non-positive parameter will imply a reset to the default
      value.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetSafetyFactor` instead.

   .. versionchanged:: 6.3.0

      The default default was changed from 0.96 to 0.9. The maximum value is now
      exactly 1.0 rather than strictly less than 1.0.



.. c:function:: int ERKStepSetSmallNumEFails(void* arkode_mem, int small_nef)

   Specifies the threshold for "multiple" successive error failures
   before the *etamxf* parameter from
   :c:func:`ERKStepSetMaxEFailGrowth()` is applied.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *small_nef* -- bound to determine "multiple" for *etamxf* (default is 2).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      Any non-positive parameter will imply a reset to the default value.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetSmallNumEFails` instead.



.. c:function:: int ERKStepSetStabilityFn(void* arkode_mem, ARKExpStabFn EStab, void* estab_data)

   Sets the problem-dependent function to estimate a stable
   time step size for the explicit portion of the ODE system.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *EStab* -- name of user-supplied stability function.
      * *estab_data* -- pointer to user data passed to *EStab* every time
        it is called.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      This function should return an estimate of the absolute
      value of the maximum stable time step for the the ODE system.  It
      is not required, since accuracy-based adaptivity may be sufficient
      for retaining stability, but this can be quite useful for problems
      where the right-hand side function :math:`f(t,y)` contains stiff
      terms.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetStabilityFn` instead.




.. _ARKODE.Usage.ERKStep.ERKStepRootfindingInput:


Rootfinding optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: int ERKStepSetRootDirection(void* arkode_mem, int* rootdir)

   Specifies the direction of zero-crossings to be located and returned.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *rootdir* -- state array of length *nrtfn*, the number of root
        functions :math:`g_i`  (the value of *nrtfn* was supplied in
        the call to :c:func:`ERKStepRootInit()`).  If ``rootdir[i] ==
        0`` then crossing in either direction for :math:`g_i` should be
        reported.  A value of +1 or -1 indicates that the solver
        should report only zero-crossings where :math:`g_i` is
        increasing or decreasing, respectively.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      The default behavior is to monitor for both zero-crossing directions.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetRootDirection` instead.



.. c:function:: int ERKStepSetNoInactiveRootWarn(void* arkode_mem)

   Disables issuing a warning if some root function appears
   to be identically zero at the beginning of the integration.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``

   **Notes:**
      ERKStep will not report the initial conditions as a
      possible zero-crossing (assuming that one or more components
      :math:`g_i` are zero at the initial time).  However, if it appears
      that some :math:`g_i` is identically zero at the initial time
      (i.e., :math:`g_i` is zero at the initial time *and* after the
      first step), ERKStep will issue a warning which can be disabled with
      this optional input function.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetNoInactiveRootWarn` instead.




.. _ARKODE.Usage.ERKStep.InterpolatedOutput:

Interpolated output function
--------------------------------

.. c:function:: int ERKStepGetDky(void* arkode_mem, sunrealtype t, int k, N_Vector dky)

   Computes the *k*-th derivative of the function
   :math:`y` at the time *t*,
   i.e., :math:`y^{(k)}(t)`, for values of the
   independent variable satisfying :math:`t_n-h_n \le t \le t_n`, with
   :math:`t_n` as current internal time reached, and :math:`h_n` is
   the last internal step size successfully used by the solver.  This
   routine uses an interpolating polynomial of degree *min(degree, 5)*,
   where *degree* is the argument provided to
   :c:func:`ERKStepSetInterpolantDegree()`.  The user may request *k* in the
   range {0,..., *min(degree, kmax)*} where *kmax* depends on the choice of
   interpolation module. For Hermite interpolants *kmax = 5* and for Lagrange
   interpolants *kmax = 3*.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *t* -- the value of the independent variable at which the
        derivative is to be evaluated.
      * *k* -- the derivative order requested.
      * *dky* -- output vector (must be allocated by the user).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_BAD_K* if *k* is not in the range {0,..., *min(degree, kmax)*}.
      * *ARK_BAD_T* if *t* is not in the interval :math:`[t_n-h_n, t_n]`
      * *ARK_BAD_DKY* if the *dky* vector was ``NULL``
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``

   **Notes:**
      It is only legal to call this function after a successful
      return from :c:func:`ERKStepEvolve()`.

      A user may access the values :math:`t_n` and :math:`h_n` via the
      functions :c:func:`ERKStepGetCurrentTime()` and
      :c:func:`ERKStepGetLastStep()`, respectively.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetDky` instead.



.. _ARKODE.Usage.ERKStep.OptionalOutputs:

Optional output functions
------------------------------


.. _ARKODE.Usage.ERKStep.ERKStepMainOutputs:

Main solver optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. c:function:: int ERKStepGetWorkSpace(void* arkode_mem, long int* lenrw, long int* leniw)

   Returns the ERKStep real and integer workspace sizes.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *lenrw* -- the number of ``sunrealtype`` values in the ERKStep workspace.
      * *leniw* -- the number of integer values in the ERKStep workspace.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetWorkSpace` instead.



.. c:function:: int ERKStepGetNumSteps(void* arkode_mem, long int* nsteps)

   Returns the cumulative number of internal steps taken by
   the solver (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *nsteps* -- number of steps taken in the solver.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetNumSteps` instead.



.. c:function:: int ERKStepGetActualInitStep(void* arkode_mem, sunrealtype* hinused)

   Returns the value of the integration step size used on the first step.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *hinused* -- actual value of initial step size.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   **Notes:**
      Even if the value of the initial integration step was
      specified by the user through a call to
      :c:func:`ERKStepSetInitStep()`, this value may have been changed by
      ERKStep to ensure that the step size fell within the prescribed
      bounds :math:`(h_{min} \le h_0 \le h_{max})`, or to satisfy the
      local error test condition.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetActualInitStep` instead.



.. c:function:: int ERKStepGetLastStep(void* arkode_mem, sunrealtype* hlast)

   Returns the integration step size taken on the last successful
   internal step.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *hlast* -- step size taken on the last internal step.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetLastStep` instead.



.. c:function:: int ERKStepGetCurrentStep(void* arkode_mem, sunrealtype* hcur)

   Returns the integration step size to be attempted on the next internal step.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *hcur* -- step size to be attempted on the next internal step.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetCurrentStep` instead.



.. c:function:: int ERKStepGetCurrentTime(void* arkode_mem, sunrealtype* tcur)

   Returns the current internal time reached by the solver.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *tcur* -- current internal time reached.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetCurrentTime` instead.



.. c:function:: int ERKStepGetTolScaleFactor(void* arkode_mem, sunrealtype* tolsfac)

   Returns a suggested factor by which the user's
   tolerances should be scaled when too much accuracy has been
   requested for some internal step.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *tolsfac* -- suggested scaling factor for user-supplied tolerances.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetTolScaleFactor` instead.



.. c:function:: int ERKStepGetErrWeights(void* arkode_mem, N_Vector eweight)

   Returns the current error weight vector.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *eweight* -- solution error weights at the current time.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   **Notes:**
      The user must allocate space for *eweight*, that will be
      filled in by this function.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetErrWeights` instead.



.. c:function:: int ERKStepGetStepStats(void* arkode_mem, long int* nsteps, sunrealtype* hinused, sunrealtype* hlast, sunrealtype* hcur, sunrealtype* tcur)

   Returns many of the most useful optional outputs in a single call.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *nsteps* -- number of steps taken in the solver.
      * *hinused* -- actual value of initial step size.
      * *hlast* -- step size taken on the last internal step.
      * *hcur* -- step size to be attempted on the next internal step.
      * *tcur* -- current internal time reached.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetStepStats` instead.



.. c:function:: int ERKStepPrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt)

   Outputs all of the integrator and other statistics.

   **Arguments:**
     * *arkode_mem* -- pointer to the ERKStep memory block.
     * *outfile* -- pointer to output file.
     * *fmt* -- the output format:

       * :c:enumerator:`SUN_OUTPUTFORMAT_TABLE` -- prints a table of values
       * :c:enumerator:`SUN_OUTPUTFORMAT_CSV` -- prints a comma-separated list
         of key and value pairs e.g., ``key1,value1,key2,value2,...``

   **Return value:**
     * *ARK_SUCCESS* -- if the output was successfully.
     * *CV_MEM_NULL* -- if the ERKStep memory was ``NULL``.
     * *CV_ILL_INPUT* -- if an invalid formatting option was provided.

   .. note::

      The Python module ``tools/suntools`` provides utilities to read and output
      the data from a SUNDIALS CSV output file using the key and value pair
      format.

   .. versionadded:: 5.2.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodePrintAllStats` instead.



.. c:function:: char* ERKStepGetReturnFlagName(long int flag)

   Returns the name of the ERKStep constant corresponding to *flag*.
   See :ref:`ARKODE.Constants`.

   **Arguments:**
      * *flag* -- a return flag from an ERKStep function.

   **Return value:**
      The return value is a string containing the name of
      the corresponding constant.

   .. warning::

      The user is responsible for freeing the returned string.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetReturnFlagName` instead.



.. c:function:: int ERKStepGetNumExpSteps(void* arkode_mem, long int* expsteps)

   Returns the cumulative number of stability-limited steps
   taken by the solver (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *expsteps* -- number of stability-limited steps taken in the solver.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetNumExpSteps` instead.



.. c:function:: int ERKStepGetNumAccSteps(void* arkode_mem, long int* accsteps)

   Returns the cumulative number of accuracy-limited steps
   taken by the solver (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *accsteps* -- number of accuracy-limited steps taken in the solver.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetNumAccSteps` instead.



.. c:function:: int ERKStepGetNumStepAttempts(void* arkode_mem, long int* step_attempts)

   Returns the cumulative number of steps attempted by the solver (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *step_attempts* -- number of steps attempted by solver.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetNumStepAttempts` instead.



.. c:function:: int ERKStepGetNumRhsEvals(void* arkode_mem, long int* nf_evals)

   Returns the number of calls to the user's right-hand
   side function, :math:`f` (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *nf_evals* -- number of calls to the user's :math:`f(t,y)` function.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   .. deprecated:: 6.2.0

      Use :c:func:`ARKodeGetNumRhsEvals` instead.


.. c:function:: int ERKStepGetNumErrTestFails(void* arkode_mem, long int* netfails)

   Returns the number of local error test failures that
   have occurred (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *netfails* -- number of error test failures.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetNumErrTestFails` instead.



.. c:function:: int ERKStepGetCurrentButcherTable(void* arkode_mem, ARKodeButcherTable *B)

   Returns the Butcher table currently in use by the solver.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *B* -- pointer to the Butcher table structure.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   **Notes:**
      The :c:type:`ARKodeButcherTable` data structure is defined as a
      pointer to the following C structure:

      .. code-block:: c

         typedef struct ARKodeButcherTableMem {

           int q;              /* method order of accuracy       */
           int p;              /* embedding order of accuracy    */
           int stages;         /* number of stages               */
           sunrealtype **A;    /* Butcher table coefficients     */
           sunrealtype *c;     /* canopy node coefficients       */
           sunrealtype *b;     /* root node coefficients         */
           sunrealtype *d;     /* embedding coefficients         */

         } *ARKodeButcherTable;

      For more details see :numref:`ARKodeButcherTable`.

.. c:function:: int ERKStepGetEstLocalErrors(void* arkode_mem, N_Vector ele)

   Returns the vector of estimated local truncation errors
   for the current step.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *ele* -- vector of estimated local truncation errors.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   **Notes:**
      The user must allocate space for *ele*, that will be
      filled in by this function.

      The values returned in *ele* are valid only after a successful call
      to :c:func:`ERKStepEvolve()` (i.e., it returned a non-negative value).

      The *ele* vector, together with the *eweight* vector from
      :c:func:`ERKStepGetErrWeights()`, can be used to determine how the
      various components of the system contributed to the estimated local
      error test.  Specifically, that error test uses the WRMS norm of a
      vector whose components are the products of the components of these
      two vectors.  Thus, for example, if there were recent error test
      failures, the components causing the failures are those with largest
      values for the products, denoted loosely as ``eweight[i]*ele[i]``.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetEstLocalErrors` instead.



.. c:function:: int ERKStepGetTimestepperStats(void* arkode_mem, long int* expsteps, long int* accsteps, long int* step_attempts, long int* nf_evals, long int* netfails)

   Returns many of the most useful time-stepper statistics in a single call.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *expsteps* -- number of stability-limited steps taken in the solver.
      * *accsteps* -- number of accuracy-limited steps taken in the solver.
      * *step_attempts* -- number of steps attempted by the solver.
      * *nf_evals* -- number of calls to the user's :math:`f(t,y)` function.
      * *netfails* -- number of error test failures.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``



.. c:function:: int ERKStepGetNumConstrFails(void* arkode_mem, long int* nconstrfails)

   Returns the cumulative number of constraint test failures (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *nconstrfails* -- number of constraint test failures.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetNumConstrFails` instead.



.. c:function:: int ERKStepGetUserData(void* arkode_mem, void** user_data)

   Returns the user data pointer previously set with
   :c:func:`ERKStepSetUserData`.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *user_data* -- memory reference to a user data pointer

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   .. versionadded:: 5.3.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetUserData` instead.



.. _ARKODE.Usage.ERKStep.ERKStepRootOutputs:

Rootfinding optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: int ERKStepGetRootInfo(void* arkode_mem, int* rootsfound)

   Returns an array showing which functions were found to
   have a root.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *rootsfound* -- array of length *nrtfn* with the indices of the
        user functions :math:`g_i` found to have a root (the value of
        *nrtfn* was supplied in the call to
        :c:func:`ERKStepRootInit()`).  For :math:`i = 0 \ldots`
        *nrtfn*-1, ``rootsfound[i]`` is nonzero if :math:`g_i` has a
        root, and 0 if not.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   **Notes:**
      The user must allocate space for *rootsfound* prior to
      calling this function.

      For the components of :math:`g_i` for which a root was found, the
      sign of ``rootsfound[i]`` indicates the direction of
      zero-crossing.  A value of +1 indicates that :math:`g_i` is
      increasing, while a value of -1 indicates a decreasing :math:`g_i`.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetRootInfo` instead.



.. c:function:: int ERKStepGetNumGEvals(void* arkode_mem, long int* ngevals)

   Returns the cumulative number of calls made to the
   user's root function :math:`g`.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *ngevals* -- number of calls made to :math:`g` so far.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetNumGEvals` instead.



.. _ARKODE.Usage.ERKStep.ERKStepExtraOutputs:

General usability functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. c:function:: int ERKStepWriteParameters(void* arkode_mem, FILE *fp)

   Outputs all ERKStep solver parameters to the provided file pointer.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *fp* -- pointer to use for printing the solver parameters.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   **Notes:**
      The *fp* argument can be ``stdout`` or ``stderr``, or it
      may point to a specific file created using ``fopen``.

      When run in parallel, only one process should set a non-NULL value
      for this pointer, since parameters for all processes would be
      identical.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeWriteParameters` instead.



.. c:function:: int ERKStepWriteButcher(void* arkode_mem, FILE *fp)

   Outputs the current Butcher table to the provided file pointer.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *fp* -- pointer to use for printing the Butcher table.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   **Notes:**
      The *fp* argument can be ``stdout`` or ``stderr``, or it
      may point to a specific file created using ``fopen``.

      When run in parallel, only one process should set a non-NULL value
      for this pointer, since tables for all processes would be
      identical.

   .. deprecated:: 6.1.0

      Use :c:func:`ERKStepGetCurrentButcherTable` and :c:func:`ARKodeButcherTable_Write`
      instead.





.. _ARKODE.Usage.ERKStep.Reinitialization:

ERKStep re-initialization function
-------------------------------------

To reinitialize the ERKStep module for the solution of a new problem,
where a prior call to :c:func:`ERKStepCreate` has been made, the
user must call the function :c:func:`ERKStepReInit()`.  The new
problem must have the same size as the previous one.  This routine
retains the current settings for all ERKstep module options and
performs the same input checking and initializations that are done in
:c:func:`ERKStepCreate`, but it performs no memory allocation as it
assumes that the existing internal memory is sufficient for the new
problem.  A call to this re-initialization routine deletes the
solution history that was stored internally during the previous
integration, and deletes any previously-set *tstop* value specified via a
call to :c:func:`ERKStepSetStopTime()`.  Following a successful call to
:c:func:`ERKStepReInit()`, call :c:func:`ERKStepEvolve()` again for the
solution of the new problem.

The use of :c:func:`ERKStepReInit()` requires that the number of
Runge--Kutta stages, denoted by *s*, be no larger for the new problem than
for the previous problem.  This condition is automatically fulfilled
if the method order *q* is left unchanged.

One important use of the :c:func:`ERKStepReInit()` function is in the
treating of jump discontinuities in the RHS function.  Except in cases
of fairly small jumps, it is usually more efficient to stop at each
point of discontinuity and restart the integrator with a readjusted
ODE model, using a call to this routine.  To stop when the location
of the discontinuity is known, simply make that location a value of
``tout``.  To stop when the location of the discontinuity is
determined by the solution, use the rootfinding feature.  In either
case, it is critical that the RHS function *not* incorporate the
discontinuity, but rather have a smooth extension over the
discontinuity, so that the step across it (and subsequent rootfinding,
if used) can be done efficiently.  Then use a switch within the RHS
function (communicated through ``user_data``) that can be flipped
between the stopping of the integration and the restart, so that the
restarted problem uses the new values (which have jumped).  Similar
comments apply if there is to be a jump in the dependent variable
vector.


.. c:function:: int ERKStepReInit(void* arkode_mem, ARKRhsFn f, sunrealtype t0, N_Vector y0)

   Provides required problem specifications and re-initializes the
   ERKStep time-stepper module.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *f* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the right-hand side function in :math:`\dot{y} = f(t,y)`.
      * *t0* -- the initial value of :math:`t`.
      * *y0* -- the initial condition vector :math:`y(t_0)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ERKStep memory was ``NULL``
      * *ARK_MEM_FAIL*  if a memory allocation failed
      * *ARK_ILL_INPUT* if an argument had an illegal value.

   **Notes:**
      All previously set options are retained but may be updated by calling
      the appropriate "Set" functions.

      If an error occurred, :c:func:`ERKStepReInit()` also
      sends an error message to the error handler function.




.. _ARKODE.Usage.ERKStep.Reset:

ERKStep reset function
----------------------


.. c:function:: int ERKStepReset(void* arkode_mem, sunrealtype tR, N_Vector yR)

   Resets the current ERKStep time-stepper module state to the provided
   independent variable value and dependent variable vector.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *tR* -- the value of the independent variable :math:`t`.
      * *yR* -- the value of the dependent variable vector :math:`y(t_R)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ERKStep memory was ``NULL``
      * *ARK_MEM_FAIL*  if a memory allocation failed
      * *ARK_ILL_INPUT* if an argument had an illegal value.

   **Notes:**
      By default the next call to :c:func:`ERKStepEvolve()` will use the step size
      computed by ERKStep prior to calling :c:func:`ERKStepReset()`. To set a
      different step size or have ERKStep estimate a new step size use
      :c:func:`ERKStepSetInitStep()`.

      All previously set options are retained but may be updated by calling
      the appropriate "Set" functions.

      If an error occurred, :c:func:`ERKStepReset()` also sends an error message to
      the error handler function.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeReset` instead.




.. _ARKODE.Usage.ERKStep.Resizing:

ERKStep system resize function
-------------------------------------


.. c:function:: int ERKStepResize(void* arkode_mem, N_Vector yR, sunrealtype hscale, sunrealtype tR, ARKVecResizeFn resize, void* resize_data)

   Re-sizes ERKStep with a different state vector but with comparable
   dynamical time scale.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *yR* -- the newly-sized solution vector, holding the current
        dependent variable values :math:`y(t_R)`.
      * *hscale* -- the desired time step scaling factor (i.e. the next
        step will be of size *h\*hscale*).
      * *tR* -- the current value of the independent variable
        :math:`t_R` (this must be consistent with *yR*).
      * *resize* -- the user-supplied vector resize function (of type
        :c:func:`ARKVecResizeFn()`.
      * *resize_data* -- the user-supplied data structure to be passed
        to *resize* when modifying internal ERKStep vectors.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ERKStep memory was ``NULL``
      * *ARK_NO_MALLOC* if *arkode_mem* was not allocated.
      * *ARK_ILL_INPUT* if an argument had an illegal value.

   **Notes:**
      If an error occurred, :c:func:`ERKStepResize()` also sends an error
      message to the error handler function.

      If inequality constraint checking is enabled a call to
      :c:func:`ERKStepResize()` will disable constraint checking. A call
      to :c:func:`ERKStepSetConstraints()` is required to re-enable constraint
      checking.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeResize` instead.
