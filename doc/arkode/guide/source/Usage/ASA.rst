.. _ARKODE.Usage.ASA:

Adjoint Sensitivity Analysis
============================

.. versionadded:: 7.3.0 (ARKODE 6.3.0)

The previous sections discuss using ARKODE for the integration of forward ODE
models. This section discusses how to use ARKODE for adjoint sensitivity
analysis as introduced in :numref:`ARKODE.Mathematics.ASA`. To use ARKStep,
ERKStep, or MRIStep for adjoint sensitivity analysis (ASA), users simply setup
the forward integration as usual (following :numref:`ARKODE.Usage.Skeleton`)
with a few differences. Below we provide an updated version of the ARKODE usage
in section :numref:`ARKODE.Usage.Skeleton` where steps that are unchanged are
*italicized*. The example code
``examples/arkode/C_serial/ark_lotka_volterra_asa.c`` demonstrates these steps
in detail.

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

   Call :c:func:`ERKStepCreateAdjointStepper`,
   :c:func:`ARKStepCreateAdjointStepper`, or
   :c:func:`MRIStepCreateAdjointStepper`. MRIStep additionally requires an
   inner adjoint stepper and an :c:type:`MRIStepInnerAdjointProblem`, as
   described in :numref:`ARKODE.Usage.ASA.MRIStep`.

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
adjoint sensitivity analysis with ERKStep, ARKStep, and MRIStep.

.. c:function:: int ERKStepCreateAdjointStepper(void* arkode_mem, SUNAdjRhsFn f, sunrealtype tf, N_Vector sf, SUNContext sunctx, SUNAdjointStepper* adj_stepper_ptr)

   Creates a :c:type:`SUNAdjointStepper` object compatible with the provided ERKStep instance for
   integrating the adjoint sensitivity system :eq:`ARKODE_DISCRETE_ADJOINT`.

   :param arkode_mem: a pointer to the ERKStep memory block.
   :param f: the adjoint right hand side function which implements
             :math:`\Lambda = f_y(t, y, p)^* \lambda` and, if sensitivities
             with respect to parameters should be computed,
             :math:`\nu = f_p(t, y, p)^* \lambda`.
   :param tf: the terminal time for the adjoint sensitivity system.
   :param sf: the sensitivity vector holding the adjoint system terminal
              condition. This must be an :ref:`NVECTOR_MANYVECTOR
              <NVectors.ManyVector>` instance. The first subvector must be
              :math:`g_y(t_f, y(t_f), p)^* \in \mathbb{R}^N`. If sensitivities
              to parameters should be computed, then the second subvector must
              be :math:`g_p(t_f, y(t_f), p)^* \in \mathbb{R}^{N_s}`, otherwise
              only one subvector should be provided.
   :param sunctx: the SUNDIALS simulation context object.
   :param adj_stepper_ptr: the newly created :c:type:`SUNAdjointStepper` object.

   :retval ARK_SUCCESS: if successful.
   :retval ARK_MEM_FAIL: if a memory allocation failed.
   :retval ARK_ILL_INPUT: if an argument has an illegal value.
   :retval ARK_SUNSTEPPER_ERR: if creation or configuration of an underlying
                               stepper fails.

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
              :math:`\Lambda = f_y^{E}(t, y, p)^* \lambda` and, if sensitivities
              with respect to parameters should be computed,
              :math:`\nu = f_p^{E}(t, y, p)^* \lambda`.
   :param fi: not yet supported, the user should pass ``NULL``.
   :param tf: the terminal time for the adjoint sensitivity system.
   :param sf: the sensitivity vector holding the adjoint system terminal
              condition. This must be a :ref:`NVECTOR_MANYVECTOR
              <NVectors.ManyVector>` instance. The first subvector must be
              :math:`g_y(t_f, y(t_f), p)^* \in \mathbb{R}^N`. If sensitivities
              to parameters should be computed, then the second subvector must
              be :math:`g_p(t_f, y(t_f), p)^* \in \mathbb{R}^{N_s}`, otherwise
              only one subvector should be provided.
   :param sunctx: the SUNDIALS simulation context object.
   :param adj_stepper_ptr: the newly created :c:type:`SUNAdjointStepper` object.

   :retval ARK_SUCCESS: if successful.
   :retval ARK_MEM_FAIL: if a memory allocation failed.
   :retval ARK_ILL_INPUT: if an argument has an illegal value.
   :retval ARK_SUNSTEPPER_ERR: if creation or configuration of an underlying
                               stepper fails.

   .. versionadded:: 7.3.0 (ARKODE 6.3.0)

   .. note::

      Currently only explicit methods with identity mass matrices are supported for ASA,
      and fixed time steps must be used.
      Furthermore, the explicit stability function, inequality constraints, and relaxation
      features are not yet compatible as they require adaptive time steps.

.. c:function:: int MRIStepCreateAdjointStepper(void* arkode_mem, SUNAdjointStepper inner_stepper, SUNAdjRhsFn fse, SUNAdjRhsFn fsi, MRIStepInnerAdjointProblem inner_problem, sunrealtype tf, N_Vector sf, SUNContext sunctx, SUNAdjointStepper* adj_stepper_ptr)

   Creates a :c:type:`SUNAdjointStepper` that reverses an explicit MRI-GARK
   integration. The returned stepper uses ``inner_stepper`` to reverse each
   fast stage IVP and applies ``fse`` to accumulate the slow contributions in
   the reverse stage recursion described in
   :numref:`ARKODE.Mathematics.ASA.MRIGARK`.

   :param arkode_mem: a pointer to the forward MRIStep memory block.
   :param inner_stepper: an adjoint stepper for the inner integrator, built
                         using the wrapped callback and terminal state from
                         ``inner_problem``. Its forward stepper must wrap the
                         same inner integrator used by ``arkode_mem`` and must
                         have a checkpoint scheme attached.
   :param fse: the adjoint right hand side function which implements
               :math:`\Lambda = f_y^{E}(t, y, p)^* \lambda` and, if
               sensitivities with respect to parameters should be computed,
               :math:`\nu = f_p^{E}(t, y, p)^* \lambda`.
   :param fsi: not yet supported, the user should pass ``NULL``.
   :param inner_problem: the inner adjoint problem used to construct
                         ``inner_stepper``.
   :param tf: the terminal time for the adjoint integration.
   :param sf: the sensitivity vector holding the adjoint system terminal
              condition. This must be a :ref:`NVECTOR_MANYVECTOR
              <NVectors.ManyVector>` instance. The first subvector must be
              :math:`g_y(t_f, y(t_f), p)^* \in \mathbb{R}^N`. If sensitivities
              to parameters should be computed, then the second subvector must
              be :math:`g_p(t_f, y(t_f), p)^* \in \mathbb{R}^{N_s}`, otherwise
              only one subvector should be provided.
   :param sunctx: the SUNDIALS simulation context object.
   :param adj_stepper_ptr: the newly created :c:type:`SUNAdjointStepper` object.

   :retval ARK_SUCCESS: if successful.
   :retval ARK_MEM_FAIL: if a memory allocation failed.
   :retval ARK_ILL_INPUT: if an argument, method, or vector layout is invalid.
   :retval ARK_SUNSTEPPER_ERR: if creation or configuration of an underlying
                               stepper fails.

   The returned stepper stores the user-data pointer from ``arkode_mem`` for
   calls to ``fse`` and installs ``inner_problem`` as the user data for
   ``inner_stepper``. It does not take ownership of either ``inner_stepper`` or
   ``inner_problem``.

   .. note::

      Currently, only fixed-step explicit MRI-GARK methods are supported.
      Inequality constraints and relaxation are not compatible with the
      MRIStep adjoint stepper.

   .. versionadded:: X.Y.Z


.. _ARKODE.Usage.ASA.MRIStep:

MRIStep adjoint setup
---------------------

As discussed in :numref:`ARKODE.Mathematics.ASA.MRIGARK`, an MRIStep adjoint
calculation reverses both the outer MRI-GARK method and the inner method used to
integrate each fast stage IVP. Consequently, it requires an outer
:c:type:`SUNAdjointStepper`, an inner :c:type:`SUNAdjointStepper`, and separate
checkpoint schemes for the two integration levels. The outer checkpoint scheme
must be attached before evolving the forward MRIStep problem, but the inner
checkpoint scheme must be attached after.

First, construct the forward inner integrator and wrap it as a
:c:type:`SUNStepper` for use by MRIStep. Configure the outer MRIStep integrator
with a fixed step size and an explicit MRI-GARK coupling table. Attach an outer
checkpoint scheme before evolving the forward MRIStep problem. The inner
checkpoint scheme may be attached after the outer forward evolution because
MRIStep recomputes each fast stage IVP before reversing it.

After the forward solution reaches :math:`t_f`, create the public terminal
sensitivity vector ``sf`` as described above. It has the layout
:math:`[\lambda]` when only initial-state sensitivities are requested and
:math:`[\lambda,\mu]` when parameter sensitivities are also requested. Then:

#. Call :c:func:`MRIStepInnerAdjointProblem_Create` with the fast adjoint
   right-hand side and ``sf``. This creates the additional forcing-moment
   variables required by the inner adjoint problems.

#. Obtain the wrapped inner adjoint right-hand side and augmented terminal
   state with :c:func:`MRIStepInnerAdjointProblem_GetAdjRhsFn` and
   :c:func:`MRIStepInnerAdjointProblem_GetTerminalState`.

#. Create the inner :c:type:`SUNAdjointStepper` using the wrapped right-hand
   side and augmented terminal state.

#. Call :c:func:`MRIStepCreateAdjointStepper` with the inner adjoint stepper,
   the slow adjoint right-hand side, the inner adjoint problem, and the public
   terminal state.

#. Evolve the outer adjoint stepper with :c:func:`SUNAdjointStepper_Evolve` or
   :c:func:`SUNAdjointStepper_OneStep`.

A condensed setup is shown below. Error checks are omitted for readability.

.. code-block:: C

   /* Evolve the configured outer MRIStep problem to tf first. */
   ARKodeSetAdjointCheckpointScheme(outer_mem, outer_checkpoint_scheme);
   ARKodeEvolve(outer_mem, tf, y, &tret, ARK_NORMAL);

   /* sf is [lambda] or [lambda, mu]. */
   MRIStepInnerAdjointProblem inner_problem = NULL;
   MRIStepInnerAdjointProblem_Create(outer_mem, fast_adj_rhs, sf,
                                     fast_user_data, &inner_problem);

   SUNAdjRhsFn inner_adj_rhs = NULL;
   N_Vector inner_sf         = NULL;
   MRIStepInnerAdjointProblem_GetAdjRhsFn(inner_problem, &inner_adj_rhs);
   MRIStepInnerAdjointProblem_GetTerminalState(inner_problem, &inner_sf);

   ARKodeSetAdjointCheckpointScheme(inner_mem, inner_checkpoint_scheme);

   SUNAdjointStepper inner_adjoint = NULL;
   ERKStepCreateAdjointStepper(inner_mem, inner_adj_rhs, tf, inner_sf,
                               sunctx, &inner_adjoint);

   SUNAdjointStepper outer_adjoint = NULL;
   MRIStepCreateAdjointStepper(outer_mem, inner_adjoint, slow_adj_rhs, NULL,
                               inner_problem, tf, sf, sunctx, &outer_adjoint);

   SUNAdjointStepper_Evolve(outer_adjoint, t0, sf, &tret);

The outer adjoint stepper does not own ``inner_adjoint`` or ``inner_problem``.
Both must remain valid throughout the outer adjoint integration. Destroy the
outer adjoint stepper first, then the inner adjoint stepper, and finally free
the inner adjoint problem. Values returned by the inner-problem getter
functions are borrowed and must not be destroyed separately.


MRIStep inner adjoint problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:type:: MRIStepInnerAdjointProblem

   An opaque object that augments the adjoint problem for an MRIStep inner
   integrator with extra parameters needed to differentiate forcing terms.

   .. versionadded:: X.Y.Z


.. c:function:: int MRIStepInnerAdjointProblem_Create(void* arkode_mem, SUNAdjRhsFn adj_f, N_Vector sf, void* user_data, MRIStepInnerAdjointProblem* problem)

   Creates an :c:type:`MRIStepInnerAdjointProblem` for the fast component of an
   MRIStep problem. The number of forcing-moment vectors is determined from the
   coupling table installed in ``arkode_mem``.

   The supplied ``sf`` is copied into a new augmented terminal state. If
   ``sf`` has the layout :math:`[\lambda]`, the augmented state is

   .. math::
      [\lambda,[\omega_0,\ldots,\omega_{n_f-1}]],

   and if ``sf`` has the layout :math:`[\lambda,\mu]`, it is

   .. math::
      [\lambda,[\omega_0,\ldots,\omega_{n_f-1},\mu]],

   where :math:`n_f` is the number of MRI forcing polynomials. Each
   :math:`\omega_k` is initialized to zero. The wrapped callback presents only
   :math:`[\lambda]` or :math:`[\lambda,\mu]` to ``adj_f`` and computes the
   forcing-moment components internally.

   :param arkode_mem: a pointer to the forward MRIStep memory block.
   :param adj_f: the fast adjoint right-hand side function. It computes
                  :math:`f_y^{F,*}(t,y,p)\lambda` and, when ``sf`` includes
                  :math:`\mu`, :math:`f_p^{F,*}(t,y,p)\lambda`.
   :param sf: the public terminal sensitivity vector. This must be a
              :ref:`ManyVector <NVectors.ManyVector>` containing either
              :math:`[\lambda]` or :math:`[\lambda,\mu]`.
   :param user_data: the pointer passed to ``adj_f`` whenever the wrapped
                     callback is evaluated.
   :param problem: on success, the newly allocated inner adjoint problem; set
                   to ``NULL`` on failure.

   :retval ARK_SUCCESS: if successful.
   :retval ARK_MEM_FAIL: if a memory allocation failed.
   :retval ARK_ILL_INPUT: if an argument or terminal-state layout is invalid.

   .. versionadded:: X.Y.Z


.. c:function:: int MRIStepInnerAdjointProblem_GetAdjRhsFn(MRIStepInnerAdjointProblem problem, SUNAdjRhsFn* adj_f)

   Returns the wrapped fast adjoint right-hand side used to construct the inner
   :c:type:`SUNAdjointStepper`.

   :param problem: the inner adjoint problem.
   :param adj_f: on return, the wrapped adjoint right-hand side function.

   :retval ARK_SUCCESS: if successful.
   :retval ARK_ILL_INPUT: if an argument is ``NULL``.

   .. versionadded:: X.Y.Z


.. c:function:: int MRIStepInnerAdjointProblem_GetTerminalState(MRIStepInnerAdjointProblem problem, N_Vector* sf)

   Returns the augmented terminal state owned by ``problem``. The returned
   vector is borrowed and must not be destroyed or retained after
   :c:func:`MRIStepInnerAdjointProblem_Free` is called.

   :param problem: the inner adjoint problem.
   :param sf: on return, the augmented terminal state.

   :retval ARK_SUCCESS: if successful.
   :retval ARK_ILL_INPUT: if an argument is ``NULL``.

   .. versionadded:: X.Y.Z


.. c:function:: int MRIStepInnerAdjointProblem_GetUserData(MRIStepInnerAdjointProblem problem, void** user_data)

   Returns the callback data required by the wrapped adjoint right-hand side.
   This is the :c:type:`MRIStepInnerAdjointProblem` itself, not the original
   user-data pointer supplied at construction. The wrapped callback passes that
   original pointer to the user's ``adj_f`` function internally.

   :param problem: the inner adjoint problem.
   :param user_data: on return, the borrowed callback-data pointer.

   :retval ARK_SUCCESS: if successful.
   :retval ARK_ILL_INPUT: if an argument is ``NULL``.

   .. note::

      :c:func:`MRIStepCreateAdjointStepper` installs this callback data on the
      inner adjoint stepper automatically. This getter is primarily useful when
      assembling or configuring an inner adjoint stepper manually.

   .. versionadded:: X.Y.Z


.. c:function:: void MRIStepInnerAdjointProblem_Free(MRIStepInnerAdjointProblem problem)

   Frees an inner adjoint problem and its augmented vectors. Calling this
   function with ``NULL`` has no effect. The problem must not be freed while an
   inner or outer adjoint stepper may still use it.

   :param problem: the inner adjoint problem to free.

   .. versionadded:: X.Y.Z
