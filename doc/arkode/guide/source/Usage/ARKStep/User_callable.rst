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

.. _ARKODE.Usage.ARKStep.UserCallable:

ARKStep User-callable functions
================================

This section describes ARKStep-specific functions that may be called by the user
to setup and then solve an IVP using the ARKStep time-stepping module.

As discussed in the main :ref:`ARKODE user-callable function introduction
<ARKODE.Usage.UserCallable>`, each of ARKODE's time-stepping modules
clarifies the categories of user-callable functions that it supports.
ARKStep supports *all categories*:

* temporal adaptivity
* implicit nonlinear and/or linear solvers
* non-identity mass matrices
* relaxation Runge--Kutta methods

ARKStep also has forcing function support when converted to a
:c:type:`SUNStepper` or :c:type:`MRIStepInnerStepper`. See
:c:func:`ARKodeCreateSUNStepper` and
:c:func:`ARKodeCreateMRIStepInnerStepper` for additional details.


.. _ARKODE.Usage.ARKStep.Initialization:

ARKStep initialization function
-------------------------------


.. c:function:: void* ARKStepCreate(ARKRhsFn fe, ARKRhsFn fi, sunrealtype t0, N_Vector y0, SUNContext sunctx)

   This function creates an internal memory block for a problem to be
   solved using the ARKStep time-stepping module in ARKODE.

   **Arguments:**
      * *fe* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the explicit portion of the right-hand side function in
        :math:`M(t)\, y'(t) = f^E(t,y) + f^I(t,y)`.
      * *fi* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the implicit portion of the right-hand side function in
        :math:`M(t)\, y'(t) = f^E(t,y) + f^I(t,y)`.
      * *t0* -- the initial value of :math:`t`.
      * *y0* -- the initial condition vector :math:`y(t_0)`.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**  If successful, a pointer to initialized problem memory
   of type ``void*``, to be passed to all user-facing ARKStep routines
   listed below.  If unsuccessful, a ``NULL`` pointer will be
   returned, and an error message will be printed to ``stderr``.

.. _ARKODE.Usage.ARKStep.OptionalInputs:

Optional input functions
------------------------

.. _ARKODE.Usage.ARKStep.ARKStepMethodInputTable:

Optional inputs for IVP method selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: int ARKStepSetImEx(void* arkode_mem)

   Specifies that both the implicit and explicit portions
   of problem are enabled, and to use an additive Runge--Kutta method.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      This is automatically deduced when neither of the function
      pointers *fe* or *fi* passed to :c:func:`ARKStepCreate` are
      ``NULL``, but may be set directly by the user if desired.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.imex".


.. c:function:: int ARKStepSetExplicit(void* arkode_mem)

   Specifies that the implicit portion of problem is disabled,
   and to use an explicit RK method.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      This is automatically deduced when the function pointer `fi`
      passed to :c:func:`ARKStepCreate` is ``NULL``, but may be set
      directly by the user if desired.

      If the problem is posed in explicit form, i.e. :math:`\dot{y} =
      f(t,y)`, then we recommend that the ERKStep time-stepper module be
      used instead.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.explicit".

.. c:function:: int ARKStepSetImplicit(void* arkode_mem)

   Specifies that the explicit portion of problem is disabled,
   and to use a diagonally implicit RK method.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      This is automatically deduced when the function pointer `fe`
      passed to :c:func:`ARKStepCreate` is ``NULL``, but may be set
      directly by the user if desired.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.implicit".



.. c:function:: int ARKStepSetTables(void* arkode_mem, int q, int p, ARKodeButcherTable Bi, ARKodeButcherTable Be)

   Specifies a customized Butcher table (or pair) for the ERK, DIRK, or ARK method.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.
      * *q* -- global order of accuracy for the ARK method.
      * *p* -- global order of accuracy for the embedded ARK method.
      * *Bi* -- the Butcher table for the implicit RK method.
      * *Be* -- the Butcher table for the explicit RK method.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      For a description of the :c:type:`ARKodeButcherTable` type and related
      functions for creating Butcher tables, see :numref:`ARKodeButcherTable`.

      To set an explicit table, *Bi* must be ``NULL``.  This automatically calls
      :c:func:`ARKStepSetExplicit()`.  However, if the problem is posed
      in explicit form, i.e. :math:`\dot{y} = f(t,y)`, then we recommend
      that the ERKStep time-stepper module be used instead of ARKStep.

      To set an implicit table, *Be* must be ``NULL``.  This automatically calls
      :c:func:`ARKStepSetImplicit()`.

      If both *Bi* and *Be* are provided, this routine automatically calls
      :c:func:`ARKStepSetImEx()`.

      When only one table is provided (i.e., *Bi* or *Be* is ``NULL``) then the
      input values of *q* and *p* are ignored and the global order of the method
      and embedding (if applicable) are obtained from the Butcher table
      structures. If both *Bi* and *Be* are non-NULL (e.g, an ImEx method is
      provided) then the input values of *q* and *p* are used as the order of the
      ARK method may be less than the orders of the individual tables. No error
      checking is performed to ensure that either *p* or *q* correctly describe the
      coefficients that were input.

      Error checking is subsequently performed at ARKStep initialization to ensure
      that *Bi* and *Be* (if non-NULL) specify DIRK and ERK methods, respectively.
      Specifically, the *A* member of *Bi* must be lower triangular with at least
      one nonzero value on the diagonal, and the *A* member of *Be* must be strictly
      lower triangular.  When both *Bi* and *Be* are non-NULL, they must agree on
      the number of internal stages, i.e., the *stages* members of both structures
      must match.

       If the inputs *Bi* or *Be* do not contain an embedding (when the
       corresponding explicit or implicit table is non-NULL), the user *must* call
       :c:func:`ARKodeSetFixedStep()` to enable fixed-step mode and set the
       desired time step size.

   **Warning:**
      This should not be used with :c:func:`ARKodeSetOrder`.



.. c:function:: int ARKStepSetTableNum(void* arkode_mem, ARKODE_DIRKTableID itable, ARKODE_ERKTableID etable)

   Indicates to use specific built-in Butcher tables for the ERK, DIRK
   or ARK method.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.
      * *itable* -- index of the DIRK Butcher table.
      * *etable* -- index of the ERK Butcher table.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      The allowable values for both the *itable* and *etable* arguments
      corresponding to built-in tables may be found in :numref:`Butcher`.

      To choose an explicit table, set *itable* to a negative value.  This
      automatically calls :c:func:`ARKStepSetExplicit()`.  However, if
      the problem is posed in explicit form, i.e. :math:`\dot{y} =
      f(t,y)`, then we recommend that the ERKStep time-stepper module be
      used instead of ARKStep.

      To select an implicit table, set *etable* to a negative value.
      This automatically calls :c:func:`ARKStepSetImplicit()`.

      If both *itable* and *etable* are non-negative, then these should
      match an existing implicit/explicit pair, listed in
      :numref:`Butcher.additive`.  This automatically calls
      :c:func:`ARKStepSetImEx()`.

      In all cases, error-checking is performed to ensure that the tables
      exist.

   **Warning:**
      This should not be used with :c:func:`ARKodeSetOrder`.




.. c:function:: int ARKStepSetTableName(void* arkode_mem, const char *itable, const char *etable)

   Indicates to use specific built-in Butcher tables for the ERK, DIRK
   or ARK method.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.
      * *itable* -- name of the DIRK Butcher table.
      * *etable* -- name of the ERK Butcher table.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument had an illegal value

   **Notes:**
      The allowable values for both the *itable* and *etable* arguments
      corresponding to built-in tables may be found in :numref:`Butcher`.
      This function is case sensitive.

      To choose an explicit table, set *itable* to ``"ARKODE_DIRK_NONE"``.
      This automatically calls :c:func:`ARKStepSetExplicit()`.  However,
      if the problem is posed in explicit form, i.e. :math:`\dot{y} =
      f(t,y)`, then we recommend that the ERKStep time-stepper module be
      used instead of ARKStep.

      To select an implicit table, set *etable* to ``"ARKODE_ERK_NONE"``.
      This automatically calls :c:func:`ARKStepSetImplicit()`.

      If both *itable* and *etable* are not none, then these should match
      an existing implicit/explicit pair, listed in
      :numref:`Butcher.additive`.  This automatically calls
      :c:func:`ARKStepSetImEx()`.

      In all cases, error-checking is performed to ensure that the tables
      exist.

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.table_names".


   **Warning:**
      This should not be used with :c:func:`ARKodeSetOrder`.


.. _ARKODE.Usage.ARKStep.OptionalOutputs:

Optional output functions
-------------------------

.. _ARKODE.Usage.ARKStep.ARKStepMainOutputs:

Main solver optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: int ARKStepGetCurrentButcherTables(void* arkode_mem, ARKodeButcherTable *Bi, ARKodeButcherTable *Be)

   Returns the explicit and implicit Butcher tables
   currently in use by the solver.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.
      * *Bi* -- pointer to the implicit Butcher table structure.
      * *Be* -- pointer to the explicit Butcher table structure.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKStep memory was ``NULL``

   **Note:**  The :c:type:`ARKodeButcherTable` data structure is defined as a
   pointer to the following C structure:

   .. code-block:: c

      typedef struct ARKStepButcherTableMem {

        int q;           /* method order of accuracy       */
        int p;           /* embedding order of accuracy    */
        int stages;      /* number of stages               */
        sunrealtype **A;    /* Butcher table coefficients     */
        sunrealtype *c;     /* canopy node coefficients       */
        sunrealtype *b;     /* root node coefficients         */
        sunrealtype *d;     /* embedding coefficients         */

      } *ARKStepButcherTable;

   For more details see :numref:`ARKodeButcherTable`.


.. c:function:: int ARKStepGetTimestepperStats(void* arkode_mem, long int* expsteps, long int* accsteps, long int* step_attempts, long int* nfe_evals, long int* nfi_evals, long int* nlinsetups, long int* netfails)

   Returns many of the most useful time-stepper statistics in a single call.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.
      * *expsteps* -- number of stability-limited steps taken in the solver.
      * *accsteps* -- number of accuracy-limited steps taken in the solver.
      * *step_attempts* -- number of steps attempted by the solver.
      * *nfe_evals* -- number of calls to the user's :math:`f^E(t,y)` function.
      * *nfi_evals* -- number of calls to the user's :math:`f^I(t,y)` function.
      * *nlinsetups* -- number of linear solver setup calls made.
      * *netfails* -- number of error test failures.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKStep memory was ``NULL``

.. _ARKODE.Usage.ARKStep.Reinitialization:

ARKStep re-initialization function
----------------------------------

.. c:function:: int ARKStepReInit(void* arkode_mem, ARKRhsFn fe, ARKRhsFn fi, sunrealtype t0, N_Vector y0)

   Provides required problem specifications and re-initializes the
   ARKStep time-stepper module.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.
      * *fe* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the explicit portion of the right-hand side function in
        :math:`M\, \dot{y} = f^E(t,y) + f^I(t,y)`.
      * *fi* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the implicit portion of the right-hand side function in
        :math:`M\, \dot{y} = f^E(t,y) + f^I(t,y)`.
      * *t0* -- the initial value of :math:`t`.
      * *y0* -- the initial condition vector :math:`y(t_0)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ARKStep memory was ``NULL``
      * *ARK_MEM_FAIL*  if a memory allocation failed
      * *ARK_ILL_INPUT* if an argument had an illegal value.

   **Notes:**
      All previously set options are retained but may be updated by calling
      the appropriate "Set" functions.

      If an error occurred, :c:func:`ARKStepReInit()` also
      sends an error message to the error handler function.
