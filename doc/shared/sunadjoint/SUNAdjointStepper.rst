..
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

.. _SUNAdjoint.Stepper:

The SUNAdjointStepper Class
===========================

.. versionadded:: 7.3.0

.. c:type:: SUNAdjointStepper

   The :c:type:`SUNAdjointStepper` class provides a package-agnostic interface to
   SUNDIALS ASA capabilities. It currently only supports the discrete ASA
   capabilities in the ARKODE package, but in the future this support may be
   expanded.

Class Methods
-------------

The :c:type:`SUNAdjointStepper` class has the following methods:

.. c:function:: SUNErrCode SUNAdjointStepper_Create(SUNStepper fwd_sunstepper, sunbooleantype own_fwd, \
   SUNStepper adj_sunstepper, sunbooleantype own_adj, suncountertype final_step_idx, sunrealtype tf, N_Vector sf, \
   SUNAdjointCheckpointScheme checkpoint_scheme, SUNContext sunctx, SUNAdjointStepper* adj_stepper)

   Creates the ``SUNAdjointStepper`` object needed to solve the adjoint problem.

   :param fwd_sunstepper: The :c:type:`SUNStepper` to be used for forward computations of the original ODE.
   :param own_fwd: Should `fwd_sunstepper` be owned (and destroyed) by the `SUNAdjointStepper` or not.
   :param adj_sunstepper: The :c:type:`SUNStepper` to be used for the backward integration of the adjoint ODE.
   :param own_adj: Should `adj_sunstepper` be owned (and destroyed) by the `SUNAdjointStepper` or not.
   :param final_step_idx: The index (step number) of the step corresponding to ``t_f`` for the forward ODE.
   :param tf: The terminal time for the forward ODE (the initial time for the adjoint ODE).
   :param sf: The terminal condition for the adjoint ODE.
   :param checkpoint_scheme: The :c:type:`SUNAdjointCheckpointScheme` object that determines the checkpointing strategy to use. This should be the same object provided to the forward integrator/stepper.
   :param sunctx: The :c:type:`SUNContext` for the simulation.
   :param adj_stepper: The :c:type:`SUNAdjointStepper` to construct (will be ``NULL`` on failure).

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_ReInit(SUNAdjointStepper self, sunrealtype t0, \
                                                    N_Vector y0,  sunrealtype tf, N_Vector sf)

   Reinitializes the adjoint stepper to solve a new problem of the same size.

   :param adj_stepper: The adjoint solver object.
   :param t0: The new initial time.
   :param y0: The new initial condition.
   :param tf: The time to start integrating the adjoint system from.
   :param sf: The terminal condition vector of sensitivity solutions :math:`\partial g/\partial y_0` and :math:`\partial g/\partial p`.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_Evolve(SUNAdjointStepper adj_stepper, sunrealtype tout,\
   N_Vector sens, sunrealtype* tret)

   Integrates the adjoint system.

   :param adj_stepper: The adjoint solver object.
   :param tout: The time at which the adjoint solution is desired.
   :param sens: The vector of sensitivity solutions :math:`\partial g/\partial y_0` and :math:`\partial g/\partial p`.
   :param tret: On return, the time reached by the adjoint solver.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_OneStep(SUNAdjointStepper adj_stepper, sunrealtype tout,\
   N_Vector sens, sunrealtype* tret)

   Evolves the adjoint system backwards one step.

   :param adj_stepper: The adjoint solver object.
   :param tout: The time at which the adjoint solution is desired.
   :param sens: The vector of sensitivity solutions :math:`\partial g/\partial y_0` and :math:`\partial g/\partial p`.
   :param tret: On return, the time reached by the adjoint solver.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_RecomputeFwd(SUNAdjointStepper adj_stepper, suncountertype start_idx,\
                                                          sunrealtype t0, N_Vector y0, sunrealtype tf)

   Evolves the forward system in time from (``start_idx``, ``t0``) to (``stop_idx``, ``tf``) with dense checkpointing.

   :param adj_stepper: The SUNAdjointStepper object.
   :param start_idx: the index of the step, w.r.t. the original forward integration, to begin forward integration from.
   :param t0: the initial time, w.r.t. the original forward integration, to start forward integration from.
   :param y0: the initial state, w.r.t. the original forward integration, to start forward integration from.
   :param tf: the final time, w.r.t. the original forward integration, to stop forward integration at.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_SetUserData(SUNAdjointStepper adj_stepper, void* user_data)

   Sets the user data pointer.

   :param adj_stepper: The SUNAdjointStepper object.
   :param user_data: the user data pointer that will be passed back to user-supplied callback functions.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_GetNumSteps(SUNAdjointStepper adj_stepper, suncountertype* num_steps)

   Retrieves the number of steps taken by the adjoint stepper.

   :param adj_stepper: The SUNAdjointStepper object.
   :param num_steps: Pointer to store the number of steps.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_GetNumRecompute(SUNAdjointStepper adj_stepper, suncountertype* num_recompute)

   Retrieves the number of recomputation steps (in the forward direction) performed by the adjoint stepper.

   :param adj_stepper: The SUNAdjointStepper object.
   :param num_recompute: Pointer to store the number of recomputations.

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. c:function:: SUNErrCode SUNAdjointStepper_PrintAllStats(SUNAdjointStepper adj_stepper, \
                                                           FILE* outfile, SUNOutputFormat fmt)

   Prints the adjoint stepper statistics/counters in a human-readable table format or CSV format.

   :param adj_stepper: The SUNAdjointStepper object.
   :param outfile: A file to write the output to.
   :param fmt: the format to write in (:c:type:`SUN_OUTPUTFORMAT_TABLE` or :c:type:`SUN_OUTPUTFORMAT_CSV`).

   :return: A :c:type:`SUNErrCode` indicating failure or success.


.. _SUNAdjoint.Stepper.UserSupplied:

User-Supplied Functions
-----------------------

.. c:type:: int (*SUNAdjRhsFn)(sunrealtype t, N_Vector y, N_Vector sens, N_Vector sens_dot, void* user_data)

   These functions compute the adjoint ODE right-hand side.

   For :ref:`ARKODE <ARKODE.Mathematics.ASA>`, this is

   .. math::
      \Lambda &= f_y^*(t, y, p) \lambda, \quad \text{and if the systems has parameters}, \\
          \nu &= f_p^*(t, y, p) \lambda.

   and corresponds to :eq:`ARKODE_ERK_ADJOINT` for explicit Runge--Kutta methods.

   **Parameters:**

   * **t** -- the current value of the independent variable.
   * **y** -- the current value of the forward solution vector.
   * **sens** -- a :ref:`NVECTOR_MANYVECTOR <NVectors.ManyVector>` object with two
     subvectors, the first subvector holds :math:`\lambda` and the second holds
     :math:`\mu` and is unused in this function.
   * **sens_dot** -- a :ref:`NVECTOR_MANYVECTOR <NVectors.ManyVector>` object with
     two subvectors, the first subvector holds :math:`\Lambda` and the second holds
     :math:`\nu`.
   * **user_data** -- the `user_data` pointer that was passed to
     :c:func:`SUNAdjointStepper_SetUserData`.

   **Returns:**

     A :c:type:`SUNAdjRhsFn` should return 0 if successful, a positive value if a
     recoverable error occurred (in which case the integrator may attempt to
     correct), or a negative value if it failed unrecoverably (in which
     case the integration is halted and an error is raised).

   .. note::

      Allocation of memory for ``y`` is handled within the integrator.

      The vector ``sens_dot`` may be uninitialized on input; it is the user's
      responsibility to fill this entire vector with meaningful values.
