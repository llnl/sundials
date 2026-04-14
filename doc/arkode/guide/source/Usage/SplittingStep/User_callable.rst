.. ----------------------------------------------------------------
   Programmer(s): Steven B. Roberts @ LLNL
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

.. _ARKODE.Usage.SplittingStep.UserCallable:

SplittingStep User-callable functions
=====================================

This section describes the SplittingStep-specific functions that may be called
by the user to setup and then solve an IVP using the SplittingStep time-stepping
module.

As discussed in the main :ref:`ARKODE user-callable function introduction
<ARKODE.Usage.UserCallable>`, each of ARKODE's time-stepping modules
clarifies the categories of user-callable functions that it supports.
SplittingStep does not support any of the categories beyond the functions that
apply for all time-stepping modules.


.. _ARKODE.Usage.SplittingStep.Initialization:

SplittingStep initialization functions
--------------------------------------

.. c:function:: void* SplittingStepCreate(SUNStepper* steppers, int partitions, sunrealtype t0, N_Vector y0, SUNContext sunctx)

   This function allocates and initializes memory for a problem to be solved
   using the SplittingStep time-stepping module in ARKODE.

   :param steppers: an array of :c:type:`SUNStepper` objects with one for each
      partition of the IVP. At minimum, they must implement the
      :c:func:`SUNStepper_Evolve`, :c:func:`SUNStepper_Reset`,
      :c:func:`SUNStepper_SetStopTime`, and
      :c:func:`SUNStepper_SetStepDirection` operations.
   :param partitions: the number of partitions, :math:`P > 1`, in the IVP.
   :param t0: the initial value of :math:`t`.
   :param y0: the initial condition vector :math:`y(t_0)`.
   :param sunctx: the :c:type:`SUNContext` object (see
      :numref:`SUNDIALS.SUNContext`)

   :return: If successful, a pointer to initialized problem memory of type
      ``void*``, to be passed to all user-facing SplittingStep routines listed
      below. If unsuccessful, a ``NULL`` pointer will be returned, and an error
      message will be printed to ``stderr``.

   **Example usage:**

      .. code-block:: C

         /* ARKODE objects for integrating individual partitions */
         void *partition_mem[] = {NULL, NULL};

         /* SUNSteppers to wrap the ARKODE objects */
         SUNStepper steppers[] = {NULL, NULL};

         /* create ARKODE objects, setting right-hand side functions and the
            initial condition */
         partition_mem[0] = ERKStepCreate(f1, t0, y0, sunctx);
         partition_mem[1] = ARKStepCreate(fe2, fi2, t0, y0, sunctx);

         /* setup ARKODE objects */
         . . .

         /* create SUNStepper wrappers for the ARKODE memory blocks */
         flag = ARKodeCreateSUNStepper(partition_mem[0], &stepper[0]);
         flag = ARKodeCreateSUNStepper(partition_mem[1], &stepper[1]);

         /* create a SplittingStep object with two partitions */
         arkode_mem = SplittingStepCreate(steppers, 2, t0, y0, sunctx);

   **Example codes:**
      * ``examples/arkode/C_serial/ark_advection_diffusion_reaction_splitting.c``
      * ``examples/arkode/C_serial/ark_analytic_partitioned.c``

   .. versionadded:: 6.2.0


Optional inputs for IVP method selection
----------------------------------------

.. c:function:: int SplittingStepSetCoefficients(void* arkode_mem, SplittingStepCoefficients coefficients)

   Specifies a set of coefficients for the operator splitting method.

   :param arkode_mem: pointer to the SplittingStep memory block.
   :param coefficients: the splitting coefficients for the method.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SplittingStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument has an illegal value

   .. note::

      For a description of the :c:type:`SplittingStepCoefficients` type and
      related functions for creating splitting coefficients see
      :numref:`ARKODE.Usage.SplittingStep.SplittingStepCoefficients`.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.splitting_coefficients_name", where ``coefficients``
      is itself constructed by passing the command-line option to
      :c:func:`SplittingStepCoefficients_LoadCoefficientsByName`.

   .. warning::

      This should not be used with :c:func:`ARKodeSetOrder`.

   .. versionadded:: 6.2.0


.. _ARKODE.Usage.SplittingStep.OptionalOutputs:


Optional output functions
------------------------------

.. c:function:: int SplittingStepGetNumEvolves(void* arkode_mem, int partition, long int *evolves)

   Returns the number of times the :c:type:`SUNStepper` for the given partition
   index has been evolved (so far).

   :param arkode_mem: pointer to the SplittingStep memory block.
   :param partition: index of the partition between 0 and :math:`P - 1` or a
      negative number to indicate the total number across all
      partitions.
   :param evolves: number of :c:type:`SUNStepper` evolves.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SplittingStep memory was ``NULL``
   :retval ARK_ILL_INPUT: if *partition* was out of bounds

   .. versionadded:: 6.2.0


SplittingStep re-initialization function
----------------------------------------

To reinitialize the SplittingStep module for the solution of a new problem,
where a prior call to :c:func:`SplittingStepCreate` has been made, the user must
call the function :c:func:`SplittingStepReInit` and re-initialize each
:c:type:`SUNStepper`.  The new problem must have the same size as the previous
one.  This routine retains the current settings for all SplittingStep module
options and performs the same input checking and initializations that are done
in :c:func:`SplittingStepCreate`, but it performs no memory allocation as it
assumes that the existing internal memory is sufficient for the new problem.  A
call to this re-initialization routine deletes the solution history that was
stored internally during the previous integration, and deletes any
previously-set *tstop* value specified via a call to
:c:func:`ARKodeSetStopTime`.  Following a successful call to
:c:func:`SplittingStepReInit`, call :c:func:`ARKodeEvolve` again for
the solution of the new problem.

One important use of the :c:func:`SplittingStepReInit` function is in the
treating of jump discontinuities in the RHS function.  Except in cases of fairly
small jumps, it is usually more efficient to stop at each point of discontinuity
and restart the integrator with a readjusted ODE model, using a call to this
routine.  To stop when the location of the discontinuity is known, simply make
that location a value of ``tout``.  To stop when the location of the
discontinuity is determined by the solution, use the rootfinding feature.  In
either case, it is critical that the RHS function *not* incorporate the
discontinuity, but rather have a smooth extension over the discontinuity, so
that the step across it (and subsequent rootfinding, if used) can be done
efficiently.  Then use a switch within the RHS function (communicated through
``user_data``) that can be flipped between the stopping of the integration and
the restart, so that the restarted problem uses the new values (which have
jumped).  Similar comments apply if there is to be a jump in the dependent
variable vector.

Another use of :c:func:`SplittingStepReInit` is changing the partitioning of
the ODE and the :c:type:`SUNStepper` objects used to evolve each partition.


.. c:function:: int SplittingStepReInit(void* arkode_mem, SUNStepper* steppers, int partitions, sunrealtype t0, N_Vector y0)

   Provides required problem specifications and re-initializes the SplittingStep
   time-stepper module.

   :param arkode_mem: pointer to the SplittingStep memory block.
   :param steppers: an array of :c:type:`SUNStepper` objects with one for each
      partition of the IVP. At minimum, they must implement the
      :c:func:`SUNStepper_Evolve`, :c:func:`SUNStepper_Reset`,
      :c:func:`SUNStepper_SetStopTime`, and
      :c:func:`SUNStepper_SetStepDirection` operations.
   :param partitions: the number of partitions, :math:`P > 1`, in the IVP.
   :param t0: the initial value of :math:`t`.
   :param y0: the initial condition vector :math:`y(t_0)`.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SplittingStep memory was ``NULL``
   :retval ARK_MEM_FAIL: if a memory allocation failed
   :retval ARK_ILL_INPUT: if an argument has an illegal value

   .. warning::

      This function does not perform any re-initialization of the
      :c:type:`SUNStepper` objects. It is up to the user to do this, if
      necessary.

   .. warning::

      If the number of partitions changes and coefficients were previously
      specified with :c:func:`SplittingStepSetCoefficients`, the coefficients
      will be reset since they are no longer compatible. Otherwise, all
      previously set options are retained but may be updated by calling the
      appropriate "Set" functions.

   .. versionadded:: 6.2.0
