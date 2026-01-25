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

.. _ARKODE.Usage.Relaxation:

Relaxation Methods
==================

This section describes user-callable functions for applying relaxation methods
with ARKODE. For more information on relaxation Runge--Kutta methods see
:numref:`ARKODE.Mathematics.Relaxation`.

.. warning::

   Relaxation support as not been evaluated with non-identity mass matrices.
   While this usage mode is supported, feedback from users who explore this
   combination would be appreciated.

Enabling or Disabling Relaxation
--------------------------------

.. c:function:: int ARKodeSetRelaxFn(void* arkode_mem, ARKRelaxFn rfn, ARKRelaxJacFn rjac)

   Attaches the user supplied functions for evaluating the relaxation function
   (``rfn``) and its Jacobian (``rjac``).

   Both ``rfn`` and ``rjac`` are required and an error will be returned if only
   one of the functions is ``NULL``. If both ``rfn`` and ``rjac`` are ``NULL``,
   relaxation is disabled.

   With DIRK and IMEX-ARK methods or when a fixed mass matrix is present,
   applying relaxation requires allocating :math:`s` additional state vectors
   (where :math:`s` is the number of stages in the method).

   :param arkode_mem: the ARKODE memory structure
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
      :math:`b^E_i \geq 0` and :math:`b^I_i \geq 0`. If these conditions are not
      satisfied, :c:func:`ARKodeEvolve` will return with an error during
      initialization.

   .. note::

      When combined with fixed time step sizes, ARKODE will attempt each step
      using the specified step size. If the step is successful, relaxation will
      be applied, effectively modifying the step size for the current step. If
      the step fails or applying relaxation fails, :c:func:`ARKodeEvolve` will
      return with an error.

   .. versionadded:: 6.1.0

Optional Input Functions
------------------------

This section describes optional input functions used to control applying
relaxation.

.. c:function:: int ARKodeSetRelaxEtaFail(void* arkode_mem, sunrealtype eta_rf)

   Sets the step size reduction factor applied after a failed relaxation
   application.

   The default value is 0.25. Input values :math:`\leq 0` or :math:`\geq 1` will
   result in the default value being used.

   :param arkode_mem: the ARKODE memory structure
   :param eta_rf: the step size reduction factor

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetRelaxLowerBound(void* arkode_mem, sunrealtype lower)

   Sets the smallest acceptable value for the relaxation parameter.

   Values smaller than the lower bound will result in a failed relaxation
   application and the step will be repeated with a smaller step size
   (determined by :c:func:`ARKodeSetRelaxEtaFail`).

   The default value is 0.8. Input values :math:`\leq 0` or :math:`\geq 1` will
   result in the default value being used.

   :param arkode_mem: the ARKODE memory structure
   :param lower: the relaxation parameter lower bound

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetRelaxUpperBound(void* arkode_mem, sunrealtype upper)

   Sets the largest acceptable value for the relaxation parameter.

   Values larger than the upper bound will result in a failed relaxation
   application and the step will be repeated with a smaller step size
   (determined by :c:func:`ARKodeSetRelaxEtaFail`).

   The default value is 1.2. Input values :math:`\leq 1` will result in the
   default value being used.

   :param arkode_mem: the ARKODE memory structure
   :param upper: the relaxation parameter upper bound

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetRelaxMaxFails(void* arkode_mem, int max_fails)

   Sets the maximum number of times applying relaxation can fail within a step
   attempt before the integration is halted with an error.

   The default value is 10. Input values :math:`\leq 0` will result in the
   default value being used.

   :param arkode_mem: the ARKODE memory structure
   :param max_fails: the maximum number of failed relaxation applications
                     allowed in a step

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetRelaxMaxIters(void* arkode_mem, int max_iters)

   Sets the maximum number of nonlinear iterations allowed when solving for the
   relaxation parameter.

   If the maximum number of iterations is reached before meeting the solve
   tolerance (determined by :c:func:`ARKodeSetRelaxResTol` and
   :c:func:`ARKodeSetRelaxTol`), the step will be repeated with a smaller
   step size (determined by :c:func:`ARKodeSetRelaxEtaFail`).

   The default value is 10. Input values :math:`\leq 0` will result in the
   default value being used.

   :param arkode_mem: the ARKODE memory structure
   :param max_iters: the maximum number of solver iterations allowed

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetRelaxSolver(void* arkode_mem, ARKRelaxSolver solver)

   Sets the nonlinear solver method used to compute the relaxation parameter.

   The default value is ``ARK_RELAX_NEWTON``.

   :param arkode_mem: the ARKODE memory structure
   :param solver: the nonlinear solver to use: ``ARK_RELAX_BRENT`` or
                  ``ARK_RELAX_NEWTON``

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``
   :retval ARK_ILL_INPUT: an invalid solver option was provided

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetRelaxResTol(void* arkode_mem, sunrealtype res_tol)

   Sets the nonlinear solver residual tolerance to use when solving
   :eq:`ARKODE_RELAX_NLS`.

   If the residual or iteration update tolerance (see
   :c:func:`ARKodeSetRelaxMaxIters`) is not reached within the maximum number of
   iterations (determined by :c:func:`ARKodeSetRelaxMaxIters`), the step will
   be repeated with a smaller step size (determined by
   :c:func:`ARKodeSetRelaxEtaFail`).

   The default value is :math:`4 \epsilon` where :math:`\epsilon` is
   floating-point precision. Input values :math:`\leq 0` will result in the
   default value being used.

   :param arkode_mem: the ARKODE memory structure
   :param res_tol: the nonlinear solver residual tolerance to use

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeSetRelaxTol(void* arkode_mem, sunrealtype rel_tol, sunrealtype abs_tol)

   Sets the nonlinear solver relative and absolute tolerance on changes in
   :math:`r` iterates when solving :eq:`ARKODE_RELAX_NLS`.

   If the residual (see :c:func:`ARKodeSetRelaxResTol`) or iterate update
   tolerance is not reached within the maximum number of iterations (determined
   by :c:func:`ARKodeSetRelaxMaxIters`), the step will be repeated with a
   smaller step size (determined by :c:func:`ARKodeSetRelaxEtaFail`).

   The default relative and absolute tolerances are :math:`4 \epsilon` and
   :math:`10^{-14}`, respectively, where :math:`\epsilon` is floating-point
   precision. Input values :math:`\leq 0` will result in the default value being
   used.

   :param arkode_mem: the ARKODE memory structure
   :param rel_tol: the nonlinear solver relative solution tolerance to use
   :param abs_tol: the nonlinear solver absolute solution tolerance to use

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 6.1.0


Optional Output Functions
-------------------------

This section describes optional output functions used to retrieve information
about the performance of the relaxation method.

.. c:function:: int ARKodeGetNumRelaxFnEvals(void* arkode_mem, long int* r_evals)

   Get the number of times the user's relaxation function was evaluated.

   :param arkode_mem: the ARKODE memory structure
   :param r_evals: the number of relaxation function evaluations

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumRelaxJacEvals(void* arkode_mem, long int* J_evals)

   Get the number of times the user's relaxation Jacobian was evaluated.

   :param arkode_mem: the ARKODE memory structure
   :param J_evals: the number of relaxation Jacobian evaluations

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumRelaxFails(void* arkode_mem, long int* fails)

   Get the total number of times applying relaxation failed.

   The counter includes the sum of the number of nonlinear solver failures
   (see :c:func:`ARKodeGetNumRelaxSolveFails`) and the number of failures due
   an unacceptable relaxation value (see :c:func:`ARKodeSetRelaxLowerBound` and
   :c:func:`ARKodeSetRelaxUpperBound`).

   :param arkode_mem: the ARKODE memory structure
   :param fails: the total number of failed relaxation attempts

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumRelaxBoundFails(void* arkode_mem, long int* fails)

   Get the number of times the relaxation parameter was deemed unacceptable.

   :param arkode_mem: the ARKODE memory structure
   :param fails: the number of failures due to an unacceptable relaxation
                 parameter value

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumRelaxSolveFails(void* arkode_mem, long int* fails)

   Get the number of times the relaxation parameter nonlinear solver failed.

   :param arkode_mem: the ARKODE memory structure
   :param fails: the number of relaxation nonlinear solver failures

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 6.1.0


.. c:function:: int ARKodeGetNumRelaxSolveIters(void* arkode_mem, long int* iters)

   Get the number of relaxation parameter nonlinear solver iterations.

   :param arkode_mem: the ARKODE memory structure
   :param iters: the number of relaxation nonlinear solver iterations

   :retval ARK_SUCCESS: the value was successfully set
   :retval ARK_MEM_NULL: ``arkode_mem`` was ``NULL``
   :retval ARK_RELAX_MEM_NULL: the internal relaxation memory structure was
                               ``NULL``

   .. versionadded:: 6.1.0
