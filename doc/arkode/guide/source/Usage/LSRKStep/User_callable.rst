.. ----------------------------------------------------------------
   Programmer(s): Mustafa Aggul @ UMBC
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

.. _ARKODE.Usage.LSRKStep.UserCallable:

LSRKStep User-callable functions
==================================

This section describes the LSRKStep-specific functions that may be called
by the user to setup and then solve an IVP using the LSRKStep time-stepping
module.  As mentioned in Section :numref:`ARKODE.Usage.UserCallable`,
shared ARKODE-level routines may be used for the large majority of LSRKStep
configuration and use.  In this section, we describe only those routines
that are specific to LSRKStep.

As discussed in the main :ref:`ARKODE user-callable function introduction
<ARKODE.Usage.UserCallable>`, each of ARKODE's time-stepping modules
clarifies the categories of user-callable functions that it supports.
LSRKStep supports the following categories:

* temporal adaptivity

LSRKStep does not have forcing function support when converted to a
:c:type:`SUNStepper` or :c:type:`MRIStepInnerStepper`. See
:c:func:`ARKodeCreateSUNStepper` and :c:func:`ARKStepCreateMRIStepInnerStepper`
for additional details.


.. _ARKODE.Usage.LSRKStep.Initialization:

LSRKStep initialization functions
---------------------------------


.. c:function:: void* LSRKStepCreateSTS(ARKRhsFn rhs, sunrealtype t0, N_Vector y0, SUNContext sunctx);

   This function allocates and initializes memory for a problem to
   be solved using STS methods from the LSRKStep time-stepping module in ARKODE.

   **Arguments:**
      * *rhs* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the right-hand side function.
      * *t0* -- the initial value of :math:`t`.
      * *y0* -- the initial condition vector :math:`y(t_0)`.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a pointer to initialized problem memory
      of type ``void*``, to be passed to all user-facing LSRKStep routines
      listed below.  If unsuccessful, a ``NULL`` pointer will be
      returned, and an error message will be printed to ``stderr``.


.. c:function:: void* LSRKStepCreateSSP(ARKRhsFn rhs, sunrealtype t0, N_Vector y0, SUNContext sunctx);

   This function allocates and initializes memory for a problem to
   be solved using SSP methods from the LSRKStep time-stepping module in ARKODE.

   **Arguments:**
      * *rhs* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the right-hand side function.
      * *t0* -- the initial value of :math:`t`.
      * *y0* -- the initial condition vector :math:`y(t_0)`.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a pointer to initialized problem memory
      of type ``void*``, to be passed to all user-facing LSRKStep routines
      listed below.  If unsuccessful, a ``NULL`` pointer will be
      returned, and an error message will be printed to ``stderr``.


.. _ARKODE.Usage.LSRKStep.OptionalInputs:

Optional input functions
-------------------------


.. c:function:: int LSRKStepSetSTSMethod(void* arkode_mem, ARKODE_LSRKMethodType method);

   This function selects the LSRK STS method that should be used.  The list of allowable
   values for this input is below. :c:func:`LSRKStepCreateSTS` defaults to using
   :c:enumerator:`ARKODE_LSRK_RKC_2`.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *method* -- Type of the method.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_ILL_INPUT* if an argument had an illegal value (e.g. typo in the method type).


.. c:function:: int LSRKStepSetSSPMethod(void* arkode_mem, ARKODE_LSRKMethodType method);

   This function selects the LSRK SSP method that should be used.  The list of allowable
   values for this input is below. :c:func:`LSRKStepCreateSSP` defaults to using
   :c:enumerator:`ARKODE_LSRK_SSP_S_2`.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *method* -- Type of the method.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_ILL_INPUT* if an argument had an illegal value (e.g. typo in the method type).


Allowable Method Families

.. c:enum:: ARKODE_LSRKMethodType

   .. c:enumerator:: ARKODE_LSRK_RKC_2

      Second order Runge--Kutta--Chebyshev method

   .. c:enumerator:: ARKODE_LSRK_RKL_2

      Second order Runge--Kutta--Legendre method

   .. c:enumerator:: ARKODE_LSRK_SSP_S_2

      Second order, s-stage SSP(s,2) method

   .. c:enumerator:: ARKODE_LSRK_SSP_S_3

      Third order, s-stage SSP(s,3) method

   .. c:enumerator:: ARKODE_LSRK_SSP_10_4

      Fourth order, 10-stage SSP(10,4) method


.. c:function:: int LSRKStepSetSTSMethodByName(void* arkode_mem, const char* emethod);

   This function selects the LSRK STS method by name. The list of allowable values
   for this input is above. :c:func:`LSRKStepCreateSTS` defaults to using
   :c:enumerator:`ARKODE_LSRK_RKC_2`.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *emethod* -- the method name.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_ILL_INPUT* if an argument had an illegal value (e.g. typo in the method name).

   .. note::

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.sts_method_name".

.. c:function:: int LSRKStepSetSSPMethodByName(void* arkode_mem, const char* emethod);

   This function selects the LSRK SSP method by name. The list of allowable values
   for this input is above. :c:func:`LSRKStepCreateSSP` defaults to using
   :c:enumerator:`ARKODE_LSRK_SSP_S_2`.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *emethod* -- the method name.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_ILL_INPUT* if an argument had an illegal value (e.g. typo in the method name).

   .. note::

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.ssp_method_name".

.. c:function:: int LSRKStepSetDomEigFn(void* arkode_mem, ARKDomEigFn dom_eig);

   Specifies the user-supplied dominant eigenvalue approximation routine to
   be used for determining the number of stages that will be used by either the
   RKC or RKL methods.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *dom_eig* -- name of user-supplied dominant eigenvalue approximation function (of type :c:func:`ARKDomEigFn()`).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if ``arkode_mem`` was ``NULL``.

   .. note::

      When using RKC or RKL methods, users must supply a :c:type:`ARKDomEigFn` function
      or attach a dominant eigenvalue estimator with :c:func:`LSRKStepSetDomEigEstimator`.


.. c:function:: int LSRKStepSetDomEigEstimator(void* arkode_mem, SUNDomEigEstimator DEE);

   Specifies the dominant eigenvalue estimator (DEE) used to determine the number of
   stages in an RKC or RKL method. This function is an alternative to supplying a
   dominant eigenvalue function with :c:func:`LSRKStepSetDomEigFn`.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *DEE* -- the dominant eigenvalue estimator to use.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if ``arkode_mem`` was ``NULL``.
      * *ARK_ILL_INPUT* if an argument had an illegal value (e.g., ``DEE`` does
        not implement the required operations)
      * *ARK_DEE_FAIL* if the call to :c:func:`SUNDomEigEstimator_SetATimes`
        failed

   .. versionadded:: 6.5.0

   .. note::

      When using RKC or RKL methods, users must supply a :c:type:`ARKDomEigFn`
      function or attach a dominant eigenvalue estimator with
      :c:func:`LSRKStepSetDomEigEstimator`.  If both are provided then the
      estimator ``DEE`` will be used and the function ignored.

      ARKODE will supply the :c:type:`SUNDomEigEstimator` with an internal
      Jacobian-vector product approximation function. Users may supply their own
      Jacobian-vector product function by calling
      :c:func:`SUNDomEigEstimator_SetATimes` after attaching the estimator with
      :c:func:`LSRKStepSetDomEigEstimator`.


.. c:function:: int LSRKStepSetDomEigFrequency(void* arkode_mem, long int nsteps);

   Specifies the number of steps after which the dominant eigenvalue information is
   considered out-of-date, and should be recomputed. This only applies to RKL and RKC methods.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *nsteps* -- the dominant eigenvalue re-computation update frequency.  A value  ``nsteps = 0``  indicates that the dominant eigenvalue will not change throughout the simulation.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if ``arkode_mem`` was ``NULL``.

   .. note::

      If LSRKStepSetDomEigFrequency routine is not called, then the default ``nsteps``
      is set to :math:`25` as recommended in :cite:p:`VSH:04`. Calling this function
      with ``nsteps < 0`` resets the default value while ``nsteps = 0`` refers to
      constant dominant eigenvalue.


      Calling this function with ``nsteps < 0`` resets the default value while
      ``nsteps = 0`` refers to constant dominant eigenvalue.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.dom_eig_frequency".

.. c:function:: int LSRKStepSetMaxNumStages(void* arkode_mem, int stage_max_limit);

   Specifies the maximum number of stages allowed within each time step.  This bound only applies to
   RKL and RKC methods.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *stage_max_limit* -- maximum allowed number of stages :math:`(>=2)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if ``arkode_mem`` was ``NULL``.

   .. note::

      If :c:func:`LSRKStepSetMaxNumStages` is not called, the default ``stage_max_limit``
      is set to :math:`200`. Calling this function with ``stage_max_limit < 2`` resets the
      default value.

      This limit should be chosen with consideration of the following proportionality:
      :math:`s^2 \sim - h\lambda`,  where :math:`s` is the number of stages used, :math:`h`
      is the current step size and :math:`\lambda` is the dominant eigenvalue.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.max_num_stages".


.. c:function:: int LSRKStepSetDomEigSafetyFactor(void* arkode_mem, sunrealtype dom_eig_safety);

   Specifies a safety factor to use for the result of the dominant eigenvalue estimation function.
   This value is used to scale the magnitude of the dominant eigenvalue, in the hope of ensuring
   a sufficient number of stages for the method to be stable.  This input is only used for RKC
   and RKL methods.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *dom_eig_safety* -- safety factor :math:`(\ge 1)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if ``arkode_mem`` was ``NULL``.

   .. note::

      If :c:func:`LSRKStepSetDomEigSafetyFactor` is not called, then the default
      ``dom_eig_safety`` is set to :math:`1.01`. Calling this function with
      ``dom_eig_safety < 1`` resets the default value.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.dom_eig_safety_factor".


.. c:function:: int LSRKStepSetNumDomEigEstInitPreprocessIters(void* arkode_mem, int num_iters);

   Specifies the number of the preprocessing iterations before the very first estimate call.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *num_iters* -- the number of iterations.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if ``arkode_mem`` was ``NULL``.

   .. versionadded:: 6.5.0

   .. note::

      If LSRKStepSetNumDomEigEstInitPreprocessIters routine is not called, then the
      default value of the estimator is used. Calling this function with ``num_iters < 0``
      resets the default.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.num_dom_eig_est_init_preprocess_iters".

.. c:function:: int LSRKStepSetNumDomEigEstPreprocessIters(void* arkode_mem, int num_iters);

   Specifies the number of the preprocessing iterations before each estimate call after the very first estimate.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *num_iters* -- the number of iterations.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if ``arkode_mem`` was ``NULL``.
      * *ARK_DEE_FAIL* if the call to
        :c:func:`SUNDomEigEstimator_SetNumPreprocessIters` failed.

   .. versionadded:: 6.5.0

   .. note::

      If LSRKStepSetNumDomEigEstPreprocessIters routine is not called, then the
      default value of 0 is used. Calling this function with ``num_iters < 0`` resets
      the default.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.num_dom_eig_est_preprocess_iters".


.. c:function:: int LSRKStepSetNumSSPStages(void* arkode_mem, int num_of_stages);

   Sets the number of stages, ``s`` in ``SSP(s, p)`` methods. This input is only utilized by SSPRK methods.

   * :c:enumerator:`ARKODE_LSRK_SSP_S_2`  -- ``num_of_stages`` must be greater than or equal to 2
   * :c:enumerator:`ARKODE_LSRK_SSP_S_3`  -- ``num_of_stages`` must be a perfect-square greater than or equal to 4
   * :c:enumerator:`ARKODE_LSRK_SSP_10_4` -- ``num_of_stages`` cannot be modified from 10, so this function should not be called.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *num_of_stages* -- number of stages :math:`(>1)` for ``SSP(s,2)`` and :math:`(n^2 = s \geq 4)` for ``SSP(s,3)``.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if ``arkode_mem`` was ``NULL``.
      * *ARK_ILL_INPUT* if an argument had an illegal value (e.g. SSP method is not declared)

   .. note::

      If :c:func:`LSRKStepSetNumSSPStages` is not called, the default ``num_of_stages``
      is set. Calling this function with ``num_of_stages <= 0`` resets the default values:

      * ``num_of_stages = 10`` for :c:enumerator:`ARKODE_LSRK_SSP_S_2`
      * ``num_of_stages = 9`` for :c:enumerator:`ARKODE_LSRK_SSP_S_3`
      * ``num_of_stages = 10`` for :c:enumerator:`ARKODE_LSRK_SSP_10_4`

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.num_ssp_stages".

.. _ARKODE.Usage.LSRKStep.OptionalOutputs:

Optional output functions
------------------------------

.. c:function:: int LSRKStepGetNumDomEigUpdates(void* arkode_mem, long int* dom_eig_num_evals);

   Returns the number of dominant eigenvalue evaluations (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *dom_eig_num_evals* -- number of calls to the user's ``dom_eig`` function.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the LSRKStep memory was ``NULL``


.. c:function:: int LSRKStepGetMaxNumStages(void* arkode_mem, int* stage_max);

   Returns the max number of stages used in any single step (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *stage_max* -- max number of stages used.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the LSRKStep memory was ``NULL``
      * *ARK_ILL_INPUT* if ``stage_max`` is illegal


.. c:function:: int LSRKStepGetNumDomEigEstRhsEvals(void* arkode_mem, long int* nfeDQ);

   Returns the number of RHS function evaluations used in the difference quotient
   Jacobian approximations (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *nfeDQ* -- number of rhs calls.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the LSRKStep memory was ``NULL``
      * *ARK_ILL_INPUT* if ``nfeDQ`` is illegal

   .. versionadded:: 6.5.0

   .. note::

      The number of RHS evaluations is non-zero only when using a dominant eigenvalue
      estimator and the internal Jacobian-vector product approximation.


.. c:function:: int LSRKStepGetNumDomEigEstIters(void* arkode_mem, long int* num_iters);

   Returns the number of iterations used in the dominant eigenvalue estimator (DEE) (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *num_iters* -- number of iterations.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the LSRKStep memory was ``NULL``
      * *ARK_ILL_INPUT* if ``num_iters`` is illegal

   .. versionadded:: 6.5.0


.. _ARKODE.Usage.LSRKStep.Reinitialization:

LSRKStep re-initialization function
-------------------------------------

To reinitialize the LSRKStep module for the solution of a new problem,
where a prior call to :c:func:`LSRKStepCreateSTS` or :c:func:`LSRKStepCreateSSP`
has been made, the user must call the function :c:func:`LSRKStepReInitSTS()`
or :c:func:`LSRKStepReInitSSP()`, accordingly.  The new problem must have
the same size as the previous one.  This routine retains the current settings
for all LSRKstep module options and performs the same input checking and
initializations that are done in :c:func:`LSRKStepCreateSTS` or
:c:func:`LSRKStepCreateSSP`, but it performs no memory allocation as it
assumes that the existing internal memory is sufficient for the new problem.
A call to this re-initialization routine deletes the solution history that
was stored internally during the previous integration, and deletes any
previously-set *tstop* value specified via a call to
:c:func:`ARKodeSetStopTime()`.  Following a successful call to
:c:func:`LSRKStepReInitSTS()` or :c:func:`LSRKStepReInitSSP()`,
call :c:func:`ARKodeEvolve()` again for the solution of the new problem.

One important use of the :c:func:`LSRKStepReInitSTS()` and
:c:func:`LSRKStepReInitSSP()` function is in the treating of jump
discontinuities in the RHS function.  Except in cases of fairly small
jumps, it is usually more efficient to stop at each point of discontinuity
and restart the integrator with a readjusted ODE model, using a call to this
routine.  To stop when the location of the discontinuity is known, simply
make that location a value of ``tout``.  To stop when the location of
the discontinuity is determined by the solution, use the rootfinding feature.
In either case, it is critical that the RHS function *not* incorporate the
discontinuity, but rather have a smooth extension over the discontinuity,
so that the step across it (and subsequent rootfinding, if used) can be done
efficiently.  Then use a switch within the RHS function (communicated through
``user_data``) that can be flipped between the stopping of the integration
and the restart, so that the restarted problem uses the new values (which
have jumped).  Similar comments apply if there is to be a jump in the
dependent variable vector.


.. c:function:: int LSRKStepReInitSTS(void* arkode_mem, ARKRhsFn rhs, sunrealtype t0, N_Vector y0);

   Provides required problem specifications and re-initializes the
   LSRKStep time-stepper module when using STS methods.

   All previously set options are retained but may be updated by calling
   the appropriate "Set" functions.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *rhs* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the right-hand side function.
      * *t0* -- the initial value of :math:`t`.
      * *y0* -- the initial condition vector :math:`y(t_0)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the LSRKStep memory was ``NULL``
      * *ARK_MEM_FAIL* if memory allocation failed
      * *ARK_NO_MALLOC* if memory allocation failed
      * *ARK_CONTROLLER_ERR* if unable to reset error controller object
      * *ARK_ILL_INPUT* if an argument had an illegal value.

   .. note::

      If using a :c:type:`SUNDomEigEstimator`, the initial guess for the
      dominant eigenvalue should be reinitialized with
      :c:func:`SUNDomEigEstimator_SetInitialGuess`.

.. c:function:: int LSRKStepReInitSSP(void* arkode_mem, ARKRhsFn rhs, sunrealtype t0, N_Vector y0);

   Provides required problem specifications and re-initializes the
   LSRKStep time-stepper module when using SSP methods.

   All previously set options are retained but may be updated by calling
   the appropriate "Set" functions.

   **Arguments:**
      * *arkode_mem* -- pointer to the LSRKStep memory block.
      * *rhs* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the right-hand side function.
      * *t0* -- the initial value of :math:`t`.
      * *y0* -- the initial condition vector :math:`y(t_0)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the LSRKStep memory was ``NULL``
      * *ARK_MEM_FAIL* if memory allocation failed
      * *ARK_NO_MALLOC* if memory allocation failed
      * *ARK_CONTROLLER_ERR* if unable to reset error controller object
      * *ARK_ILL_INPUT* if an argument had an illegal value.
