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

This section describes the currently supported ERKStep-specific functions
that may be called by the user to setup and then solve an IVP using the
ERKStep time-stepping module.  Legacy compatibility wrappers were removed;
use the underlying :ref:`ARKODE functions <ARKODE.Usage.UserCallable>` for
the deprecated functionality.

As discussed in the main :ref:`ARKODE user-callable function introduction
<ARKODE.Usage.UserCallable>`, each of ARKODE's time-stepping modules
clarifies the categories of user-callable functions that it supports.
ERKStep supports the following categories:

* temporal adaptivity
* relaxation Runge--Kutta methods

ERKStep also has forcing function support when converted to a
:c:type:`SUNStepper` or :c:type:`MRIStepInnerStepper`. See
:c:func:`ARKodeCreateSUNStepper` and
:c:func:`ARKodeCreateMRIStepInnerStepper` for additional details.


.. _ARKODE.Usage.ERKStep.Initialization:

ERKStep initialization and reinitialization functions
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
       :c:func:`ARKodeSetFixedStep()` to enable fixed-step mode and set the desired
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

