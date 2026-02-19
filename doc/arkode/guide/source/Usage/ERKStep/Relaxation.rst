.. -----------------------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
   -----------------------------------------------------------------------------
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
   -----------------------------------------------------------------------------

.. _ARKODE.Usage.ERKStep.Relaxation:

Relaxation Methods
==================

This section describes ERKStep-specific user-callable functions for applying 
relaxation methods with ERKStep.  All of these routines have been deprecated in 
favor of :ref:`shared ARKODE-level routines <ARKODE.Usage.Relaxation>`, but 
this documentation will be retained for as long as these functions are present 


Enabling or Disabling Relaxation
--------------------------------

.. c:function:: int ERKStepSetRelaxFn(void* arkode_mem, ARKRelaxFn rfn, ARKRelaxJacFn rjac)

   Attaches the user supplied functions for evaluating the relaxation function
   (``rfn``) and its Jacobian (``rjac``).

   Both ``rfn`` and ``rjac`` are required and an error will be returned if only
   one of the functions is ``NULL``. If both ``rfn`` and ``rjac`` are ``NULL``,
   relaxation is disabled.

   :param arkode_mem: the ERKStep memory structure
   :param rfn: the user-defined function to compute the relaxation function
               :math:`\xi(y)`
   :param rjac: the user-defined function to compute the relaxation Jacobian
                :math:`\xi'(y)`

   :retval ARK_SUCCESS: the function exited successfully
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_ILL_INPUT: an invalid input combination was provided (see the
                          output error message for more details)
   :retval ARK_MEM_FAIL: a memory allocation failed

   .. warning::

      Applying relaxation requires using a method of at least second order with
      :math:`b_i \geq 0`. If these conditions are not satisfied,
      :c:func:`ERKStepEvolve` will return with an error during initialization.

   .. note::

      When combined with fixed time step sizes, ERKStep will attempt each step
      using the specified step size. If the step is successful, relaxation will
      be applied, effectively modifying the step size for the current step. If
      the step fails or applying relaxation fails, :c:func:`ERKStepEvolve` will
      return with an error.

   .. versionadded:: 5.6.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetRelaxFn` instead.


Optional Input Functions
------------------------

This section describes optional input functions used to control applying
relaxation.

.. c:function:: int ERKStepSetRelaxEtaFail(void* arkode_mem, sunrealtype eta_rf)

   Sets the step size reduction factor applied after a failed relaxation
   application.

   The default value is 0.25. Input values :math:`\leq 0` or :math:`\geq 1` will
   result in the default value being used.

   :param arkode_mem: the ERKStep memory structure
   :param eta_rf: the step size reduction factor

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 5.6.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetRelaxEtaFail` instead.


.. c:function:: int ERKStepSetRelaxLowerBound(void* arkode_mem, sunrealtype lower)

   Sets the smallest acceptable value for the relaxation parameter.

   Values smaller than the lower bound will result in a failed relaxation
   application and the step will be repeated with a smaller step size
   (determined by :c:func:`ERKStepSetRelaxEtaFail`).

   The default value is 0.8. Input values :math:`\leq 0` or :math:`\geq 1` will
   result in the default value being used.

   :param arkode_mem: the ERKStep memory structure
   :param lower: the relaxation parameter lower bound

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 5.6.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetRelaxLowerBound` instead.


.. c:function:: int ERKStepSetRelaxUpperBound(void* arkode_mem, sunrealtype upper)

   Sets the largest acceptable value for the relaxation parameter.

   Values larger than the upper bound will result in a failed relaxation
   application and the step will be repeated with a smaller step size
   (determined by :c:func:`ERKStepSetRelaxEtaFail`).

   The default value is 1.2. Input values :math:`\leq 1` will result in the
   default value being used.

   :param arkode_mem: the ERKStep memory structure
   :param upper: the relaxation parameter upper bound

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 5.6.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetRelaxUpperBound` instead.


.. c:function:: int ERKStepSetRelaxMaxFails(void* arkode_mem, int max_fails)

   Sets the maximum number of times applying relaxation can fail within a step
   attempt before the integration is halted with an error.

   The default value is 10. Input values :math:`\leq 0` will result in the
   default value being used.

   :param arkode_mem: the ERKStep memory structure
   :param max_fails: the maximum number of failed relaxation applications
                     allowed in a step

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 5.6.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetRelaxMaxFails` instead.


.. c:function:: int ERKStepSetRelaxMaxIters(void* arkode_mem, int max_iters)

   Sets the maximum number of nonlinear iterations allowed when solving for the
   relaxation parameter.

   If the maximum number of iterations is reached before meeting the solve
   tolerance (determined by :c:func:`ERKStepSetRelaxResTol` and
   :c:func:`ERKStepSetRelaxTol`), the step will be repeated with a smaller
   step size (determined by :c:func:`ERKStepSetRelaxEtaFail`).

   The default value is 10. Input values :math:`\leq 0` will result in the
   default value being used.

   :param arkode_mem: the ERKStep memory structure
   :param max_iters: the maximum number of solver iterations allowed

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 5.6.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetRelaxMaxIters` instead.


.. c:function:: int ERKStepSetRelaxSolver(void* arkode_mem, ARKRelaxSolver solver)

   Sets the nonlinear solver method used to compute the relaxation parameter.

   The default value is :c:enumerator:`ARK_RELAX_NEWTON`

   :param arkode_mem: the ERKStep memory structure
   :param solver: the nonlinear solver to use
   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``
   :retval ARK_ILL_INPUT: an invalid solver option was provided

   .. versionadded:: 5.6.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetRelaxSolver` instead.


.. c:function:: int ERKStepSetRelaxResTol(void* arkode_mem, sunrealtype res_tol)

   Sets the nonlinear solver residual tolerance to use when solving
   :eq:`ARKODE_RELAX_NLS`.

   If the residual or solution tolerance (see :c:func:`ERKStepSetRelaxMaxIters`)
   is not reached within the maximum number of  iterations (determined by
   :c:func:`ERKStepSetRelaxMaxIters`), the step will be repeated with a smaller
   step size (determined by :c:func:`ERKStepSetRelaxEtaFail`).

   The default value is :math:`4 \epsilon` where :math:`\epsilon` is
   floating-point precision. Input values :math:`\leq 0` will result in the
   default value being used.

   :param arkode_mem: the ERKStep memory structure
   :param res_tol: the nonlinear solver residual tolerance to use

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 5.6.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetRelaxResTol` instead.


.. c:function:: int ERKStepSetRelaxTol(void* arkode_mem, sunrealtype rel_tol, sunrealtype abs_tol)

   Sets the nonlinear solver relative and absolute tolerance on changes in
   :math:`r` when solving :eq:`ARKODE_RELAX_NLS`.


   If the residual (see :c:func:`ERKStepSetRelaxResTol`) or solution tolerance
   is not reached within the maximum number of iterations (determined by
   :c:func:`ERKStepSetRelaxMaxIters`), the step will be repeated with a smaller
   step size (determined by :c:func:`ERKStepSetRelaxEtaFail`).

   The default relative and absolute tolerances are :math:`4 \epsilon` and
   :math:`10^{-14}`, respectively, where :math:`\epsilon` is floating-point
   precision. Input values :math:`\leq 0` will result in the default value being
   used.

   :param arkode_mem: the ERKStep memory structure
   :param rel_tol: the nonlinear solver relative solution tolerance to use
   :param abs_tol: the nonlinear solver absolute solution tolerance to use

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 5.6.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeSetRelaxTol` instead.


Optional Output Functions
-------------------------

This section describes optional output functions used to retrieve information
about the performance of the relaxation method.

.. c:function:: int ERKStepGetNumRelaxFnEvals(void* arkode_mem, long int* r_evals)

   Get the number of times the user's relaxation function was evaluated.

   :param arkode_mem: the ERKStep memory structure
   :param r_evals: the number of relaxation function evaluations

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 5.6.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetNumRelaxFnEvals` instead.


.. c:function:: int ERKStepGetNumRelaxJacEvals(void* arkode_mem, long int* J_evals)

   Get the number of times the user's relaxation Jacobian was evaluated.

   :param arkode_mem: the ERKStep memory structure
   :param J_evals: the number of relaxation Jacobian evaluations

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 5.6.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetNumRelaxJacEvals` instead.


.. c:function:: int ERKStepGetNumRelaxFails(void* arkode_mem, long int* fails)

   Get the total number of times applying relaxation failed.

   The counter includes the sum of the number of nonlinear solver failures
   (see :c:func:`ERKStepGetNumRelaxSolveFails`) and the number of failures due
   an unacceptable relaxation value (see :c:func:`ERKStepSetRelaxLowerBound` and
   :c:func:`ERKStepSetRelaxUpperBound`).

   :param arkode_mem: the ERKStep memory structure
   :param fails: the total number of failed relaxation attempts

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 5.6.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetNumRelaxFails` instead.



.. c:function:: int ERKStepGetNumRelaxBoundFails(void* arkode_mem, long int* fails)

   Get the number of times the relaxation parameter was deemed unacceptable.

   :param arkode_mem: the ERKStep memory structure
   :param fails: the number of failures due to an unacceptable relaxation
                 parameter value

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 5.6.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetNumRelaxBoundFails` instead.


.. c:function:: int ERKStepGetNumRelaxSolveFails(void* arkode_mem, long int* fails)

   Get the number of times the relaxation parameter nonlinear solver failed.

   :param arkode_mem: the ERKStep memory structure
   :param fails: the number of relaxation nonlinear solver failures

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 5.6.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetNumRelaxSolveFails` instead.


.. c:function:: int ERKStepGetNumRelaxSolveIters(void* arkode_mem, long int* iters)

   Get the number of relaxation parameter nonlinear solver iterations.

   :param arkode_mem: the ERKStep memory structure
   :param iters: the number of relaxation nonlinear solver iterations

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 5.6.0

   .. deprecated:: 6.1.0

      Use :c:func:`ARKodeGetNumRelaxSolveIters` instead.

