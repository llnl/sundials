.. _ARKODE.Usage.ASA:

Adjoint Sensitivity Analysis
============================

.. versionadded:: 7.3.0 (ARKODE 6.3.0)

The previous sections discuss using ARKODE for the integration of forward ODE models. This section
discusses how to use ARKODE for adjoint sensitivity analysis as introduced in
:numref:`ARKODE.Mathematics.ASA`. To use ARKStep for adjoint sensitivity analysis (ASA), users
simply setup the forward integration as usual (following :numref:`ARKODE.Usage.Skeleton`) with a few
differences. Below we provide an updated version of the ARKODE usage in section
:numref:`ARKODE.Usage.Skeleton` where steps that are unchanged are *italicized*. The example
code ``examples/arkode/C_serial/ark_lotka_volterra_asa.c`` demonstrates these steps in detail.

.. index:: Adjoint Sensitivity Analysis user main program

#. *Initialize parallel or multi-threaded environment, if appropriate.*

#. *Create the SUNDIALS simulation context object*

#. *Set problem dimensions, etc.*

#. *Set vector of initial values*

#. *Create ARKODE object*

#. Specify a fixed time step size.

   Currently the discrete ASA capability only allows a fixed time step size
   to be used. Call :c:func:`ARKodeSetFixedStep` to set the time step.

#. *Set optional inputs*

#. *Specify rootfinding problem*

#. Create a :c:type:`SUNAdjointCheckpointScheme` object

   Create the :c:type:`SUNAdjointCheckpointScheme` object by calling ``SUNAdjointCheckpointScheme_Create_*``.
   Available :c:type:`SUNAdjointCheckpointScheme` implementations are found in
   section :numref:`SUNAdjoint.CheckpointScheme`.

#. Attach the checkpoint scheme object to ARKODE

   Call :c:func:`ARKodeSetAdjointCheckpointScheme`.

#. *Advance solution in time*

#. *Get optional outputs*

#. Create the sensitivities vector with the terminal condition

   The sensitivities vector must be an instance of the :ref:`ManyVector N_Vector implementation <NVectors.ManyVector>`.
   You will have one subvector for the initial condition sensitivities and
   an additional subvector if you want sensitivities with respect to parameters. The vectors should
   contain the terminal conditions for the adjoint problem. The first subvector should contain
   :math:`dg(t_f,y(t_f),p)/dy(t_f)` and the second subvector should contain
   :math:`dg(t_f,y(t_f),p)/dp`.
   The subvectors can be any implementation of the :ref:`N_Vector class <NVectors>`.

   For example, in a problem with 10 state variables and 4 parameters using serial
   computations, the ManyVector can be constructed as follows:

   .. code-block:: C

      sunindextype num_equations = 10;
      sunindextype num_params    = 4;
      N_Vector sensu0            = N_VNew_Serial(num_equations, sunctx);
      N_Vector sensp             = N_VNew_Serial(num_params, sunctx);
      N_Vector sens[2]           = {sensu0, sensp};
      N_Vector sf                = N_VNew_ManyVector(2, sens, sunctx);
      // Set the terminal condition for the adjoint system, which
      // should be the the gradient of our cost function at tf.
      dgdu(u, sensu0, params);
      dgdp(u, sensp, params);

#. Create the :c:type:`SUNAdjointStepper` object

   Call :c:func:`ERKStepCreateAdjointStepper` or :c:func:`ARKStepCreateAdjointStepper`.

#. Set optional ASA input

   Refer to :numref:`SUNAdjoint.Stepper` for options.

#. Advance the adjoint sensitivity analysis ODE

   Call :c:func:`SUNAdjointStepper_Evolve` or :c:func:`SUNAdjointStepper_OneStep`.

#. Get optional ASA outputs

   Refer to :numref:`SUNAdjoint.Stepper` for options.

#. Deallocate memory for ASA objects

   Deallocate the sensitivities vector, :c:type:`SUNAdjointStepper`,
   and :c:type:`SUNAdjointCheckpointScheme` objects.

#. *Deallocate memory for solution vector*

#. Free solver memory

   Call :c:func:`SUNStepper_Destroy` and :c:func:`ARKodeFree` to free the memory
   allocated for the SUNStepper and ARKODE integrator objects.

#. *Free the SUNContext object*

#. *Finalize MPI, if used*



User Callable Functions
-----------------------

This section describes user-callable functions for performing
adjoint sensitivity analysis with methods with ERKStep and ARKStep.

.. c:function:: int ERKStepCreateAdjointStepper(void* arkode_mem, SUNAdjRhsFn f, sunrealtype tf, N_Vector sf, SUNContext sunctx, SUNAdjointStepper* adj_stepper_ptr)

   Creates a :c:type:`SUNAdjointStepper` object compatible with the provided ERKStep instance for
   integrating the adjoint sensitivity system :eq:`ARKODE_DISCRETE_ADJOINT`.

   :param arkode_mem: a pointer to the ERKStep memory block.
   :param f: the adjoint right hand side function which implements
             :math:`\Lambda = f_y^*(t, y, p) \lambda` and, if sensitivities
             with respect to parameters should be computed,
             :math:`\nu = f_p^*(t, y, p) \lambda`.
   :param tf: the terminal time for the adjoint sensitivity system.
   :param sf: the sensitivity vector holding the adjoint system terminal
              condition. This must be an :ref:`NVECTOR_MANYVECTOR
              <NVectors.ManyVector>` instance. The first subvector must be
              :math:`g^*_y(t_f, y(t_f), p) \in \mathbb{R}^N`. If sensitivities
              to parameters should be computed, then the second subvector must
              be :math:`g^*_p(t_f, y(t_f), p) \in \mathbb{R}^{N_s}`, otherwise
              only one subvector should be provided.
   :param sunctx: the SUNDIALS simulation context object.
   :param adj_stepper_ptr: the newly created :c:type:`SUNAdjointStepper` object.

   :retval ARK_SUCCESS: if successful.
   :retval ARK_MEM_FAIL: if a memory allocation failed.
   :retval ARK_ILL_INPUT: if an argument has an illegal value.

   .. versionadded:: 7.3.0 (ARKODE 6.3.0)

   .. note::

      Currently fixed time steps must be used.
      Furthermore, the explicit stability function, inequality constraints, and relaxation
      features are not yet compatible as they require adaptive time steps.


.. c:function:: int ARKStepCreateAdjointStepper(void* arkode_mem, SUNAdjRhsFn fe, SUNAdjRhsFn fi, sunrealtype tf, N_Vector sf, SUNContext sunctx, SUNAdjointStepper* adj_stepper_ptr)

   Creates a :c:type:`SUNAdjointStepper` object compatible with the provided ARKStep instance for
   integrating the adjoint sensitivity system :eq:`ARKODE_DISCRETE_ADJOINT`.

   :param arkode_mem: a pointer to the ARKStep memory block.
   :param fe: the adjoint right hand side function which implements
              :math:`\Lambda = f_y^{E,*}(t, y, p) \lambda` and, if sensitivities
              with respect to parameters should be computed,
              :math:`\nu = f_p^*(t, y, p) \lambda`.
   :param fi: not yet support, the user should pass ``NULL``.
   :param tf: the terminal time for the adjoint sensitivity system.
   :param sf: the sensitivity vector holding the adjoint system terminal
              condition. This must be a :ref:`NVECTOR_MANYVECTOR
              <NVectors.ManyVector>` instance. The first subvector must be
              :math:`g^*_y(t_f, y(t_f), p) \in \mathbb{R}^N`. If sensitivities
              to parameters should be computed, then the second subvector must
              be :math:`g^*_p(t_f, y(t_f), p) \in \mathbb{R}^{N_s}`, otherwise
              only one subvector should be provided.
   :param sunctx: the SUNDIALS simulation context object.
   :param adj_stepper_ptr: the newly created :c:type:`SUNAdjointStepper` object.

   :retval ARK_SUCCESS: if successful.
   :retval ARK_MEM_FAIL: if a memory allocation failed.
   :retval ARK_ILL_INPUT: if an argument has an illegal value.

   .. versionadded:: 7.3.0 (ARKODE 6.3.0)

   .. note::

      Currently only explicit methods with identity mass matrices are supported for ASA,
      and fixed time steps must be used.
      Furthermore, the explicit stability function, inequality constraints, and relaxation
      features are not yet compatible as they require adaptive time steps.
