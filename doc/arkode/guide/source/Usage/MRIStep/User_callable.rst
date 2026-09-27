.. ----------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
                  Daniel R. Reynolds @ UMBC
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

.. _ARKODE.Usage.MRIStep.UserCallable:

MRIStep User-callable functions
==================================

This section describes the currently supported MRIStep-specific functions
that may be called by the user to setup and then solve an IVP using the
MRIStep time-stepping module.  Legacy compatibility wrappers were removed;
use the underlying :ref:`ARKODE functions <ARKODE.Usage.UserCallable>` for
the deprecated functionality.

As discussed in the main :ref:`ARKODE user-callable function introduction
<ARKODE.Usage.UserCallable>`, each of ARKODE's time-stepping modules
clarifies the categories of user-callable functions that it supports.
MRIStep supports the following categories:

* temporal adaptivity
* implicit nonlinear and/or linear solvers

MRIStep also has forcing function support when converted to a
:c:type:`SUNStepper` or :c:type:`MRIStepInnerStepper`. See
:c:func:`ARKodeCreateSUNStepper` and
:c:func:`ARKodeCreateMRIStepInnerStepper` for additional details.


.. _ARKODE.Usage.MRIStep.Initialization:

MRIStep initialization and reinitialization functions
------------------------------------------------------

.. c:function:: void* MRIStepCreate(ARKRhsFn fse, ARKRhsFn fsi, sunrealtype t0, N_Vector y0, MRIStepInnerStepper stepper, SUNContext sunctx)

   This function allocates and initializes memory for a problem to
   be solved using the MRIStep time-stepping module in ARKODE.

   :param fse: the name of the function (of type :c:func:`ARKRhsFn()`)
               defining the explicit slow portion of the right-hand side function in
               :math:`\dot{y} = f^E(t,y) + f^I(t,y) + f^F(t,y)`.
   :param fsi: the name of the function (of type :c:func:`ARKRhsFn()`)
               defining the implicit slow portion of the right-hand side function in
               :math:`\dot{y} = f^E(t,y) + f^I(t,y) + f^F(t,y)`.
   :param t0: the initial value of :math:`t`.
   :param y0: the initial condition vector :math:`y(t_0)`.
   :param stepper: an :c:type:`MRIStepInnerStepper` for integrating the fast
                   time scale.
   :param sunctx: the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   :returns: If successful, a pointer to initialized problem memory of type ``void*``, to
             be passed to all user-facing MRIStep routines listed below.  If unsuccessful,
             a ``NULL`` pointer will be returned, and an error message will be printed to
             ``stderr``.

   **Example usage:**

      .. code-block:: C

         /* fast (inner) and slow (outer) ARKODE objects */
         void *inner_arkode_mem = NULL;
         void *outer_arkode_mem = NULL;

         /* MRIStepInnerStepper to wrap the inner (fast) object */
         MRIStepInnerStepper stepper = NULL;

         /* create an ARKODE object, setting fast (inner) right-hand side
            functions and the initial condition */
         inner_arkode_mem = *StepCreate(...);

         /* configure the inner integrator */
         retval = ARKodeSet*(inner_arkode_mem, ...);

         /* create MRIStepInnerStepper wrapper for the ARKODE integrator */
         flag = ARKodeCreateMRIStepInnerStepper(inner_arkode_mem, &stepper);

         /* create an MRIStep object, setting the slow (outer) right-hand side
            functions and the initial condition */
         outer_arkode_mem = MRIStepCreate(fse, fsi, t0, y0, stepper, sunctx)

   **Example codes:**
      * ``examples/arkode/C_serial/ark_brusselator_mri.c``
      * ``examples/arkode/C_serial/ark_twowaycouple_mri.c``
      * ``examples/arkode/C_serial/ark_brusselator_1D_mri.c``
      * ``examples/arkode/C_serial/ark_onewaycouple_mri.c``
      * ``examples/arkode/C_serial/ark_reaction_diffusion_mri.c``
      * ``examples/arkode/C_serial/ark_kpr_mri.c``
      * ``examples/arkode/CXX_parallel/ark_diffusion_reaction_p.cpp``
      * ``examples/arkode/CXX_serial/ark_test_kpr_nestedmri.cpp``
        (uses MRIStep within itself)



.. c:function:: void* MRIStepCreateExtSTS(ARKRhsFn fd, ARKRhsFn fe, ARKRhsFn fi, sunrealtype t0, N_Vector y0, SUNContext sunctx)

   This function allocates and initializes memory for a problem to be solved
   using an ExtSTS time-stepping method.

   :param fd: the user-defined function for the diffusive dynamics (required).
   :param fe: the user-defined function for the explicit dynamics
      (may be ``NULL`` if no explicit dynamics are present).
   :param fi: the user-defined function for the implicit dynamics
      (may be ``NULL`` if no implicit dynamics are present).
   :param t0: the initial value of :math:`t`.
   :param y0: the initial condition vector :math:`y(t_0)`.
   :param sunctx: the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`).

   :return: If successful, a pointer to initialized problem memory of type
      ``void*``, to be passed to all user-facing MRIStep or ARKODE routines.
      If unsuccessful, a ``NULL`` pointer will be returned, and an error
      message will be printed to ``stderr``.

   .. warning::

      Although ExtSTS methods are packaged within MRIStep, the inner STS solver
      is not subcycled as with other MRI methods.  By default, the inner STS
      solver is configured to be of Runge--Kutta--Chebyshev type, with a maximum
      of 10000 stages per step.  Although users may request the STS solver to
      change these defaults or configure other relevant options, they should never
      change any LSRKStep settings related to time step adaptivity or
      interpolated output, nor should they configure options that would
      request the STS solver to stop a time step prematurely (e.g., through
      calls to :c:func:`ARKodeRootInit` or :c:func:`ARKodeSetConstraints`).
      Similarly, users should not employ "H-Tol" multirate time step
      controllers on the object returned using :c:func:`MRIStepCreateExtSTS`.

   .. versionadded:: 7.9.0 (ARKODE 6.9.0)

   **Example usage:**

      .. code-block:: C

         /* create ExtSTS instantiation of MRIStep object */
         void* extsts_mem = MRIStepCreateExtSTS(fd, fe, fi, t0, y0, sunctx);

         /* access the inner STS stepper object */
         void* sts_mem = NULL;
         retval = MRIStepGetSTS(extsts_mem, &sts_mem);

         /* configure ExtSTS integrator */
         retval = MRIStepSet*(extsts_mem, ...);
         retval = ARKodeSet*(extsts_mem, ...);

         /* configure inner STS integrator */
         retval = LSRKStepSet*(sts_mem, ...);

   **Example codes:**
      * ``examples/arkode/CXX_serial/ark_adr1d_extsts.cpp``


.. c:function:: int MRIStepSetPreInnerFn(void* arkode_mem, MRIStepPreInnerFn prefn)

   Specifies the function called *before* each inner integration.

   :param arkode_mem: pointer to the MRIStep memory block.
   :param prefn: the name of the C function (of type :c:func:`MRIStepPreInnerFn()`)
                 defining pre inner integration function.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the MRIStep memory is ``NULL``



.. c:function:: int MRIStepSetPostInnerFn(void* arkode_mem, MRIStepPostInnerFn postfn)

   Specifies the function called *after* each inner integration.

   :param arkode_mem: pointer to the MRIStep memory block.
   :param postfn: the name of the C function (of type :c:func:`MRIStepPostInnerFn()`)
                  defining post inner integration function.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the MRIStep memory is ``NULL``




.. _ARKODE.Usage.MRIStep.MRIStepMethodInput:

Optional inputs for IVP method selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _ARKODE.Usage.MRIStep.MRIStepMethodInputTable:
.. table:: Optional inputs for IVP method selection

   +--------------------------------+-------------------------------------+----------+
   | Optional input                 | Function name                       | Default  |
   +================================+=====================================+==========+
   | Select the default MRI method  | :c:func:`ARKodeSetOrder()`          | 3        |
   | of a given order               |                                     |          |
   +--------------------------------+-------------------------------------+----------+
   | Set MRI coupling coefficients  | :c:func:`MRIStepSetCoupling()`      | internal |
   +--------------------------------+-------------------------------------+----------+


.. c:function:: int MRIStepSetCoupling(void* arkode_mem, MRIStepCoupling C)

   Specifies a customized set of slow-to-fast coupling coefficients for the MRI method.

   :param arkode_mem: pointer to the MRIStep memory block.
   :param C: the table of coupling coefficients for the MRI method.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the MRIStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument has an illegal value

   .. note::

      For a description of the :c:type:`MRIStepCoupling` type and related
      functions for creating Butcher tables see :numref:`ARKODE.Usage.MRIStep.MRIStepCoupling`.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.coupling_table_name", where ``C``
      is itself constructed by passing the command-line option to
      :c:func:`MRIStepCoupling_LoadTableByName`.

   .. warning::

      This should not be used with :c:func:`ARKodeSetOrder`.



.. _ARKODE.Usage.MRIStep.MRIStepSolverInput:

Optional inputs for implicit stage solves
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: int MRIStepGetNumInnerStepperFails(void* arkode_mem, long int* inner_fails)

   Returns the number of recoverable failures reported by the inner stepper (so far).

   :param arkode_mem: pointer to the MRIStep memory block.
   :param inner_fails: number of failed fast (inner) integrations.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the MRIStep memory was ``NULL``

   .. versionadded:: 7.2.0 (ARKODE 6.2.0)



.. c:function:: int MRIStepGetSTS(void* arkode_mem, void** stsptr)

   Returns a pointer to the LSRKStep super-time-stepping object
   used within the MRIStep module when performing ExtSTS time-stepping.

   :param arkode_mem: pointer to the MRIStep memory block.
   :param stsptr: pointer to the LSRKStep super-time-stepping object.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the MRIStep memory was ``NULL``

   .. warning::

      If the user wishes to set different ``user_data`` pointers for the
      MRIStep and LSRKStep components of the ExtSTS solver, they should first
      call :c:func:`ARKodeSetUserData` on the MRIStep object, and then call
      :c:func:`ARKodeSetUserData` on the LSRKStep object returned from
      :c:func:`MRIStepCreateExtSTS`.  By default, the LSRKStep object will
      inherit the ``user_data`` pointer from the MRIStep object.

   .. versionadded:: 7.9.0 (ARKODE 6.9.0)



.. c:function:: int MRIStepGetCurrentCoupling(void* arkode_mem, MRIStepCoupling *C)

   Returns the MRI coupling table currently in use by the solver.

   :param arkode_mem: pointer to the MRIStep memory block.
   :param C: pointer to slow-to-fast MRI coupling structure.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the MRIStep memory was ``NULL``

   .. note::

      The *MRIStepCoupling* data structure is defined in
      the header file ``arkode/arkode_mristep.h``.  For more details
      see :numref:`ARKODE.Usage.MRIStep.MRIStepCoupling`.


.. c:function:: int MRIStepGetLastInnerStepFlag(void* arkode_mem, int* flag)

   Returns the last return value from the inner stepper.

   :param arkode_mem: pointer to the MRIStep memory block.
   :param flag: inner stepper return value.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the MRIStep memory was ``NULL``



.. c:function:: int MRIStepReInit(void* arkode_mem, ARKRhsFn fse, ARKRhsFn fsi, sunrealtype t0, N_Vector y0)

   Provides required problem specifications and re-initializes the
   MRIStep outer (slow) stepper.

   :param arkode_mem: pointer to the MRIStep memory block.
   :param fse: the name of the function (of type :c:func:`ARKRhsFn()`)
               defining the explicit slow portion of the right-hand side function in
               :math:`\dot{y} = f^E(t,y) + f^I(t,y) + f^F(t,y)`.
   :param fsi: the name of the function (of type :c:func:`ARKRhsFn()`)
               defining the implicit slow portion of the right-hand side function in
               :math:`\dot{y} = f^E(t,y) + f^I(t,y) + f^F(t,y)`.
   :param t0: the initial value of :math:`t`.
   :param y0: the initial condition vector :math:`y(t_0)`.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL:  if the MRIStep memory was ``NULL``
   :retval ARK_MEM_FAIL:  if a memory allocation failed
   :retval ARK_ILL_INPUT: if an argument has an illegal value.

   .. note::

      If the inner (fast) stepper also needs to be reinitialized, its
      reinitialization function should be called before calling
      :c:func:`MRIStepReInit()` to reinitialize the outer stepper.

      All previously set options are retained but may be updated by calling
      the appropriate "Set" functions.

      If an error occurred, :c:func:`MRIStepReInit()` also
      sends an error message to the error handler function.




.. c:function:: int MRIStepReInitExtSTS(void* arkode_mem, ARKRhsFn fd, ARKRhsFn fe, ARKRhsFn fi, sunrealtype t0, N_Vector y0)

   Provides required problem specifications and re-initializes
   MRIStep for an Extended Super Time Stepping (ExtSTS) method.

   :param arkode_mem: pointer to the MRIStep memory block.
   :param fd: the user-defined function for the diffusive dynamics
      of the problem-defining right-hand side function
      :math:`\dot{y} = f^D(t,y) + f^E(t,y) + f^I(t,y)` (required).
   :param fe: the user-defined function for the explicit dynamics
      of the problem-defining right-hand side function
      :math:`\dot{y} = f^D(t,y) + f^E(t,y) + f^I(t,y)`
      (may be ``NULL`` if no explicit dynamics are present).
   :param fi: the user-defined function for the implicit dynamics
      of the problem-defining right-hand side function
      :math:`\dot{y} = f^D(t,y) + f^E(t,y) + f^I(t,y)`
      (may be ``NULL`` if no implicit dynamics are present).
   :param t0: the initial value of :math:`t`.
   :param y0: the initial condition vector :math:`y(t_0)`.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL:  if the MRIStep memory was ``NULL``
   :retval ARK_MEM_FAIL:  if a memory allocation failed
   :retval ARK_ILL_INPUT: if an argument has an illegal value.

   .. note::

      All previously set options are retained but may be updated by calling
      the appropriate "Set" functions.

      If an error occurred, :c:func:`MRIStepReInitExtSTS()` also
      sends an error message to the error handler function.

   .. versionadded:: 7.9.0 (ARKODE 6.9.0)


.. _ARKODE.Usage.MRIStep.Reset:

MRIStep reset function
----------------------
