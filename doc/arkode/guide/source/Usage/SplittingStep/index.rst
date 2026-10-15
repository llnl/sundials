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

.. _ARKODE.Usage.SplittingStep:

============================================
Using the SplittingStep time-stepping module
============================================

This section is concerned with the use of the SplittingStep time-stepping module
for the solution of initial value problems (IVPs) in a C or C++ language
setting.  Usage of SplittingStep follows that of the rest of ARKODE, and so in
this section we primarily focus on those usage aspects that are specific to
SplittingStep.

As discussed in the main :ref:`ARKODE user-callable function introduction
<ARKODE.Usage.UserCallable>`, each of ARKODE's time-stepping modules clarifies
the categories of user-callable functions that it supports.  SplittingStep does
not support any of the categories beyond the functions that apply for all
time-stepping modules.


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

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)


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

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)


.. _ARKODE.Usage.SplittingStep.OptionalOutputs:


Optional output functions
-------------------------

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

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)


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

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)


.. _ARKODE.Usage.SplittingStep.SplittingStepCoefficients:

Operator Splitting Coefficients Data Structure
----------------------------------------------

SplittingStep supplies several functions to construct operator splitting
coefficients of various orders and partitions. There are also a number of
built-in methods of fixed orders and partitions (see
:numref:`ARKODE.Usage.SplittingStep.SplittingStepCoefficients.Coefficients`).
Finally, a user may construct a custom set of coefficients and attach it with
:c:func:`SplittingStepSetCoefficients`. The operator splitting coefficients are
stored in a :c:type:`SplittingStepCoefficients` object which is a pointer to a
:c:struct:`SplittingStepCoefficientsMem` structure:

.. c:type:: SplittingStepCoefficientsMem *SplittingStepCoefficients

.. c:struct::  SplittingStepCoefficientsMem

   Structure for storing the coefficients defining an operator splitting method.

   As described in :numref:`ARKODE.Mathematics.SplittingStep`, an operator
   splitting method is defined by a vector :math:`\alpha \in \mathbb{R}^{r}` and
   a tensor :math:`\beta \in \mathbb{R}^{r \times (s + 1) \times P}` where :math:`r` is
   the number of sequential methods, :math:`s` is the number of stages, and :math:`P`
   is the number of partitions.

   .. c:member:: sunrealtype *alpha

      An array containing the weight of each sequential method used to produce the
      overall operator splitting solution. The array is of length
      ``[sequential_methods]``.

   .. c:member:: sunrealtype ***beta

      A three-dimensional array containing the time nodes of the partition integrations.
      The array has dimensions ``[sequential_methods][stages+1][partitions]``.

   .. c:member:: int sequential_methods

      The number of sequential methods, :math:`r`, combined to produce the overall
      operator splitting solution

   .. c:member:: int stages

      The number of stages, :math:`s`

   .. c:member:: int partitions

      The number of partitions, :math:`P`, in the IVP

   .. c:member:: int order

      The method order of accuracy


.. _ARKODE.Usage.SplittingStep.SplittingStepCoefficients.Functions:

SplittingStepCoefficients Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


This section describes the functions for creating and interacting with operator
splitting coefficients. The function prototypes and as well as the relevant
integer constants are defined ``arkode/arkode_splittingstep.h``.

.. _ARKODE.Usage.SplittingStep.SplittingStepCoefficients.Functions.Table:
.. table:: SplittingStepCoefficients functions

   +--------------------------------------------------------------+--------------------------------------------------------------+
   | Function name                                                | Description                                                  |
   +--------------------------------------------------------------+--------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_LoadCoefficients()`       | Load a pre-defined SplittingStepCoefficients by ID           |
   +--------------------------------------------------------------+--------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_LoadCoefficientsByName()` | Load a pre-defined SplittingStepCoefficients by name         |
   +--------------------------------------------------------------+--------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_IDToName()`               | Convert a pre-defined SplittingStepCoefficients to its name  |
   +--------------------------------------------------------------+--------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_LieTrotter()`             | Create a Lie--Trotter splitting method                       |
   +--------------------------------------------------------------+--------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_Strang()`                 | Create a Strang splitting method                             |
   +--------------------------------------------------------------+--------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_SymmetricParallel()`      | Create a symmetrization of the Lie--Trotter splitting method |
   +--------------------------------------------------------------+--------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_ThirdOrderSuzuki()`       | Create a third order composition method of Suzuki            |
   +--------------------------------------------------------------+--------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_TripleJump()`             | Create an arbitrary order, three-jump composition method     |
   +--------------------------------------------------------------+--------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_SuzukiFractal()`          | Create an arbitrary order, five-jump composition method      |
   +--------------------------------------------------------------+--------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_Alloc()`                  | Allocate an empty :c:type:`SplittingStepCoefficients`        |
   |                                                              | object                                                       |
   +--------------------------------------------------------------+--------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_Create()`                 | Create a new :c:type:`SplittingStepCoefficients` object      |
   |                                                              | from coefficient arrays                                      |
   +--------------------------------------------------------------+--------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_Copy()`                   | Create a copy of a :c:type:`SplittingStepCoefficients`       |
   |                                                              | object                                                       |
   +--------------------------------------------------------------+--------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_Destroy()`                | Deallocate a :c:type:`SplittingStepCoefficients` object      |
   +--------------------------------------------------------------+--------------------------------------------------------------+
   | :c:func:`SplittingStepCoefficients_Write()`                  | Write the :c:type:`SplittingStepCoefficients` object to an   |
   |                                                              | output file                                                  |
   +--------------------------------------------------------------+--------------------------------------------------------------+


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_LoadCoefficients(ARKODE_SplittingCoefficientsID method)

   Retrieves specified splitting coefficients. For further information on the
   current set of splitting coefficients and their corresponding identifiers,
   see
   :numref:`ARKODE.Usage.SplittingStep.SplittingStepCoefficients.Coefficients`.

   :param method: the splitting coefficients identifier.
   :return: A :c:type:`SplittingStepCoefficients` structure if successful or a
      ``NULL`` pointer if ``method`` was invalid or an allocation error occurred.

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)



.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_LoadCoefficientsByName(const char *method)

   Retrieves specified splitting coefficients. For further information on the
   current set of splitting coefficients and their corresponding name, see
   :numref:`ARKODE.Usage.SplittingStep.SplittingStepCoefficients.Coefficients`.

   :param method: the splitting coefficients identifier.
   :return: A :c:type:`SplittingStepCoefficients` structure if successful or a
      ``NULL`` pointer if ``method`` was invalid or an allocation error occurred.

   .. note::

      This function is case sensitive.

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)


.. c:function:: const char* SplittingStepCoefficients_IDToName(ARKODE_SplittingCoefficientsID method)

   Converts specified splitting coefficients ID to a string of the same name.
   For further information on the current set of splitting coefficients and
   their corresponding name, see
   :numref:`ARKODE.Usage.SplittingStep.SplittingStepCoefficients.Coefficients`.

   :param method: the splitting coefficients identifier.
   :return: A :c:type:`SplittingStepCoefficients` structure if successful or a
      ``NULL`` pointer if ``method`` was invalid or an allocation error occurred.

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_LieTrotter(int partitions)

   Create the coefficients for the first order Lie--Trotter splitting,
   see :eq:`ARKODE_Lie-Trotter`.

   :param partitions: The number of partitions, :math:`P > 1`, in the IVP.
   :return: A :c:type:`SplittingStepCoefficients` structure if successful or a
      ``NULL`` pointer if ``partitions`` was invalid or an allocation error
      occurred.

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_Strang(int partitions)

   Create the coefficients for the second order Strang splitting
   :cite:p:`Strang:68`, see :eq:`ARKODE_Strang`.

   :param partitions: The number of partitions, :math:`P > 1`, in the IVP.
   :return: A :c:type:`SplittingStepCoefficients` structure if successful or a
      ``NULL`` pointer if ``partitions`` was invalid or an allocation error
      occurred.

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_Parallel(int partitions)

   Create the coefficients for the first order parallel splitting method

   .. math::
      y_n = \phi^1_{h_n}(y_{n-1}) + \phi^2_{h_n}(y_{n-1}) + \dots
      + \phi^P_{h_n}(y_{n-1}) + (1 - P) y_{n-1}.

   :param partitions: The number of partitions, :math:`P > 1`, in the IVP.
   :return: A :c:type:`SplittingStepCoefficients` structure if successful or a
      ``NULL`` pointer if ``partitions`` was invalid or an allocation error
      occurred.

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_SymmetricParallel(int partitions)

   Create the coefficients for the second order, symmetrized Lie--Trotter
   splitting :cite:p:`Strang:63`

   .. math::
      y_n = \frac{1}{2} \left( L_{h_n}(y_{n-1}) + L^*_{h_n}(y_{n-1}) \right),

   where :math:`L_{h_n}` is the Lie--Trotter splitting :eq:`ARKODE_Lie-Trotter`
   and :math:`L^*_{h_n}` is its adjoint :eq:`ARKODE_Lie-Trotter_adjoint`.

   :param partitions: The number of partitions, :math:`P > 1`, in the IVP.
   :return: A :c:type:`SplittingStepCoefficients` structure if successful or a
      ``NULL`` pointer if ``partitions`` was invalid or an allocation error
      occurred.

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_ThirdOrderSuzuki(int partitions)

   Create the coefficients for a splitting method of Suzuki :cite:p:`Suzuki:92`

   .. math::
      y_n = \left( L_{p_1 h_n} \circ L^*_{p_2 h_n} \circ L_{p_3 h_n}
      \circ L^*_{p_4 h_n} \circ L_{p_5 h_n} \right) (y_{n-1}),

   where :math:`L_{h_n}` is the Lie--Trotter splitting :eq:`ARKODE_Lie-Trotter`
   and :math:`L^*_{h_n}` is its adjoint :eq:`ARKODE_Lie-Trotter_adjoint`. The
   parameters :math:`p_1, \dots, p_5` are selected to give third order.

   :param partitions: The number of partitions, :math:`P > 1`, in the IVP.
   :return: A :c:type:`SplittingStepCoefficients` structure if successful or a
      ``NULL`` pointer if ``partitions`` was invalid or an allocation error
      occurred.

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_TripleJump(int partitions, int order)

   Create the coefficients for the triple jump splitting method
   :cite:p:`CrGo:89`

   .. math::
      T_{h_n}^{[2]} &= S_{h_n}, \\
      T_{h_n}^{[i+2]} &= T_{\gamma_1 h_n}^{[i]} \circ
      T_{(1 - 2 \gamma_1) h_n}^{[i]} \circ T_{\gamma_1 h_n}^{[i]}, \\
      y_n &= T_{h_n}^{[\mathrm{order}]}(y_{n-1}),

   where :math:`S` is the Strang splitting :eq:`ARKODE_Strang` and
   :math:`\gamma_1` is selected to increase the order by two each recursion.

   :param partitions: The number of partitions, :math:`P > 1`, in the IVP.
   :param order: A positive even number for the method order of accuracy.
   :return: A :c:type:`SplittingStepCoefficients` structure if successful or a
      ``NULL`` pointer if an argument was invalid or an allocation error
      occurred.

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_SuzukiFractal(int partitions, int order)

   Create the coefficients for the quintuple jump splitting method
   :cite:p:`Suzuki:90`

   .. math::
      Q_{h_n}^{[2]} &= S_{h_n}, \\
      Q_{h_n}^{[i+2]} &= Q_{\gamma_1 h_n}^{[i]} \circ
      Q_{\gamma_1 h_n}^{[i]} \circ Q_{(1 - 4 \gamma_1) h_n}^i \circ
      Q_{\gamma_1 h_n}^{[i]} \circ Q_{\gamma_1 h_n}^{[i]}, \\
      y_n &= Q_{h_n}^{[\mathrm{order}]}(y_{n-1}),

   where :math:`S` is the Strang splitting :eq:`ARKODE_Strang` and
   :math:`\gamma_1` is selected to increase the order by two each recursion.

   :param partitions: The number of partitions, :math:`P > 1`, in the IVP.
   :param order: A positive even number for the method order of accuracy.
   :return: A :c:type:`SplittingStepCoefficients` structure if successful or a
      ``NULL`` pointer if an argument was invalid or an allocation error
      occurred.

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_Alloc(int sequential_methods, int stages, int partitions)

   Allocates an empty :c:type:`SplittingStepCoefficients` object.

   :param sequential_methods: The number of sequential methods, :math:`r \geq 1`,
      combined to produce the overall operator splitting solution.
   :param stages: The number of stages, :math:`s \geq 1`.
   :param partitions: The number of partitions, :math:`P > 1`, in the IVP.
   :return: A :c:type:`SplittingStepCoefficients` structure if successful or a
      ``NULL`` pointer if an argument was invalid or an allocation error
      occurred.

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_Create(int sequential_methods, int stages, int partitions, int order, sunrealtype* alpha, sunrealtype* beta)

   Allocates a :c:type:`SplittingStepCoefficients` object and fills it with the
   given values.

   :param sequential_methods: The number of sequential methods, :math:`r \geq 1`
      combined to produce the overall operator splitting solution.
   :param stages: The number of stages, :math:`s \geq 1`.
   :param partitions: The number of partitions, :math:`P > 1` in the IVP.
   :param order: The method order of accuracy.
   :param alpha: An array of length ``sequential_methods`` containing the weight
      of each sequential method used to produce the overall operator splitting
      solution.
   :param beta: An array of length
      ``sequential_methods * (stages+1) * partitions`` containing the time nodes
      of the partition integrations in the C order

      .. math::
         & \beta_{1,1,1}, \dots, \beta_{1,1,P},
         \dots,
         \beta_{1,s+1,1}, \dots, \beta_{1,s+1,P},
         \dots, \\
         & \beta_{2,1,1}, \dots, \beta_{2,1,P},
         \dots,
         \beta_{2,s+1,1}, \dots, \beta_{2,s+1,P},
         \dots, \\
         & \beta_{r,1,1}, \dots, \beta_{r,1,P},
         \dots,
         \beta_{r,s+1,1}, \dots, \beta_{r,s+1,P}.

   :return: A :c:type:`SplittingStepCoefficients` structure if successful or a
      ``NULL`` pointer if an argument was invalid or an allocation error
      occurred.

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)


.. c:function:: SplittingStepCoefficients SplittingStepCoefficients_Copy(SplittingStepCoefficients coefficients)

   Creates copy of the given splitting coefficients.

   :param coefficients: The splitting coefficients to copy.
   :return: A :c:type:`SplittingStepCoefficients` structure if successful or a
      ``NULL`` pointer if an allocation error occurred.

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)


.. c:function:: void SplittingStepCoefficients_Destroy(SplittingStepCoefficients* coefficients)

   Deallocate the splitting coefficients memory.

   :param coefficients: A pointer to the splitting coefficients.

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)


.. c:function:: void SplittingStepCoefficients_Write(SplittingStepCoefficients coefficients, FILE* outfile)

   Write the splitting coefficients to the provided file pointer.

   :param coefficients: The splitting coefficients.
   :param outfile: Pointer to use for printing the splitting coefficients. It
      can be ``stdout`` or ``stderr``, or it may point to a specific file
      created using ``fopen``.

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)


.. _ARKODE.Usage.SplittingStep.SplittingStepCoefficients.Coefficients:

Operator Splitting Coefficients
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SplittingStep currently provides several pre-defined coefficients for problems
with two partitions.  We list the identifiers, order of accuracy, and relevant
references for each in the table below. We use the naming convention

.. code-block:: text

   ARKODE_SPLITTING_<name>_<stages>_<order>_<partitions>

Each of the splitting coefficients that are packaged with SplittingStep are
specified by a unique ID having type:

.. c:enum:: ARKODE_SplittingCoefficientsID

with values specified for each method below (e.g.,
``ARKODE_SPLITTING_LIE_TROTTER_1_1_2``).

.. table:: Operator splitting coefficients.

   ======================================  =====  =====================
   Table name                              Order        Reference
   ======================================  =====  =====================
   ``ARKODE_SPLITTING_LIE_TROTTER_1_1_2``  1
   ``ARKODE_SPLITTING_STRANG_2_2_2``       2      :cite:p:`Strang:68`
   ``ARKODE_SPLITTING_BEST_2_2_2``         2      :cite:p:`AuHoKeKo:16`
   ``ARKODE_SPLITTING_SUZUKI_3_3_2``       3      :cite:p:`Suzuki:92`
   ``ARKODE_SPLITTING_RUTH_3_3_2``         3      :cite:p:`Ruth:93`
   ``ARKODE_SPLITTING_YOSHIDA_4_4_2``      4      :cite:p:`Yoshida:90`
   ``ARKODE_SPLITTING_YOSHIDA_8_6_2``      6      :cite:p:`Yoshida:90`
   ======================================  =====  =====================


.. _ARKODE.Usage.SplittingStep.SplittingStepCoefficients.Default:

Default Operator Splitting Coefficients
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default SplittingStep coefficients are Lie--Trotter. If a particular order
is requested with :c:func:`ARKodeSetOrder`, the following are the default for
each order

.. table:: Default operator splitting coefficients by order.

   ============  ==========================================================
   Order         Default operator splitting coefficients
   ============  ==========================================================
   1             :c:func:`SplittingStepCoefficients_LieTrotter`
   2             :c:func:`SplittingStepCoefficients_Strang`
   3             :c:func:`SplittingStepCoefficients_ThirdOrderSuzuki`
   4, 6, 8, ...  :c:func:`SplittingStepCoefficients_TripleJump`
   5, 7, 9, ...  Warning: this will select a triple jump method of the next
                 even order
   ============  ==========================================================
