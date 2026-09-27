.. ----------------------------------------------------------------
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

.. _ARKODE.Usage.SPRKStep.UserCallable:

SPRKStep User-callable functions
==================================

This section describes the SPRKStep-specific functions that may be called
by the user to setup and then solve an IVP using the SPRKStep time-stepping
module. Some of these user-callable functions are specific to SPRKStep,
as explained below.

As discussed in the main :ref:`ARKODE user-callable function introduction
<ARKODE.Usage.UserCallable>`, each of ARKODE's time-stepping modules
clarifies the categories of user-callable functions that it supports.
SPRKStep supports only the basic set of user-callable functions, and
does not support any of the restricted groups (time adaptivity, implicit
solvers, etc.).

SPRKStep does not have forcing function support when converted to a
:c:type:`SUNStepper` or :c:type:`MRIStepInnerStepper`. See
:c:func:`ARKodeCreateSUNStepper` and :c:func:`ARKodeCreateMRIStepInnerStepper`
for additional details.


.. _ARKODE.Usage.SPRKStep.Initialization:

SPRKStep initialization and deallocation functions
------------------------------------------------------


.. c:function:: void* SPRKStepCreate(ARKRhsFn f1, ARKRhsFn f2, sunrealtype t0,\
                                     N_Vector y0, SUNContext sunctx)

   This function allocates and initializes memory for a problem to
   be solved using the SPRKStep time-stepping module in ARKODE.

   :param f1: the name of the C function (of type :c:func:`ARKRhsFn()`) defining :math:`f_1(t,q) = -\frac{\partial V(t,q)}{\partial q}`
   :param f2: the name of the C function (of type :c:func:`ARKRhsFn()`) defining :math:`f_2(t,p) = \frac{\partial T(t,p)}{\partial p}`
   :param t0: the initial value of :math:`t`
   :param y0: the initial condition vector :math:`y(t_0)`
   :param sunctx: the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   :returns: If successful, a pointer to initialized problem memory of type
             ``void*``, to be passed to all user-facing SPRKStep routines listed
             below.  If unsuccessful, a ``NULL`` pointer will be returned, and
             an error message will be printed to ``stderr``.

   .. warning::

      SPRKStep requires a partitioned problem where ``f1`` should only modify
      the q variables and ``f2`` should only modify the p variables (or vice
      versa). However, the vector passed to these functions is the full vector
      with both p and q. The ordering of the variables is determined implicitly
      by the user when they set the initial conditions.

.. c:function:: int SPRKStepSetMethod(void* arkode_mem, ARKodeSPRKTable sprk_table)

   Specifies the SPRK method.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param sprk_table: the SPRK method table.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument had an illegal value

   .. note::

      No error checking is performed on the coefficients contained in the
      table to ensure its declared order of accuracy.

   .. warning::

      This should not be used with :c:func:`ARKodeSetOrder`.


.. c:function:: int SPRKStepSetMethodName(void* arkode_mem, const char* method)

   Specifies the SPRK method by its name.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param method: the SPRK method name.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory is ``NULL``
   :retval ARK_ILL_INPUT: if an argument had an illegal value

   .. note::

      This routine will be called by :c:func:`ARKodeSetOptions`
      when using the key "arkid.method_name".

   .. warning::

      This should not be used with :c:func:`ARKodeSetOrder`.

.. c:function:: int SPRKStepGetCurrentMethod(void* arkode_mem, ARKodeSPRKTable *sprk_table)

   Returns the SPRK method coefficient table currently in use by the solver.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param sprk_table: pointer to the SPRK method table.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL: if the SPRKStep memory was ``NULL``

.. c:function:: int SPRKStepReInit(void* arkode_mem, ARKRhsFn f1, ARKRhsFn f2, sunrealtype t0, N_Vector y0)

   Provides required problem specifications and re-initializes the SPRKStep
   time-stepper module.

   All previously set options are retained but may be updated by calling the
   appropriate "Set" functions.

   If an error occurred, :c:func:`SPRKStepReInit()` also sends an error message
   to the error handler function.

   :param arkode_mem: pointer to the SPRKStep memory block.
   :param f1: the name of the C function (of type :c:func:`ARKRhsFn()`) defining :math:`f1(t,q) = \frac{\partial V(t,q)}{\partial q}`
   :param f2: the name of the C function (of type :c:func:`ARKRhsFn()`) defining :math:`f2(t,p) = \frac{\partial T(t,p)}{\partial p}`
   :param t0: the initial value of :math:`t`.
   :param y0: the initial condition vector :math:`y(t_0)`.

   :retval ARK_SUCCESS: if successful
   :retval ARK_MEM_NULL:  if the SPRKStep memory was ``NULL``
   :retval ARK_MEM_FAIL:  if a memory allocation failed
   :retval ARK_ILL_INPUT: if an argument had an illegal value.
