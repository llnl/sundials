..
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

.. _SUNNonlinSol.API:

===============================
The SUNNonlinearSolver API
===============================

The SUNNonlinSol API defines several nonlinear solver operations that enable
SUNDIALS integrators to utilize any SUNNonlinSol implementation that
provides the required functions. These functions can be divided into three
categories. The first are the core nonlinear solver functions. The second
consists of "set" routines to supply the nonlinear solver with
functions provided by the SUNDIALS time integrators and to modify solver
parameters. The final group consists of "get" routines for retrieving nonlinear
solver statistics. All of these functions are defined in the header file
``sundials/sundials_nonlinearsolver.h``.



.. _SUNNonlinSol.API.CoreFn:

SUNNonlinearSolver core functions
-----------------------------------------------------

The core nonlinear solver functions consist of two required functions to get the
nonlinear solver type (:c:func:`SUNNonlinSolGetType`) and solve the nonlinear system
(:c:func:`SUNNonlinSolSolve`). The remaining three functions for nonlinear solver
initialization (:c:func:`SUNNonlinSolInitialize`), setup
(:c:func:`SUNNonlinSolSetup`), and destruction (:c:func:`SUNNonlinSolFree`) are optional.

.. c:enum:: SUNNonlinearSolver_Type

   An identifier indicating the form of the nonlinear system expected by the
   nonlinear solver.

   .. c:enumerator:: SUNNONLINEARSOLVER_ROOTFIND

      The nonlinear solver expects systems in rootfinding form :math:`F(y) = 0`

   .. c:enumerator:: SUNNONLINEARSOLVER_FIXEDPOINT

      The nonlinear solver expects systems in fixed-point form :math:`G(y) = y`.


.. c:function:: SUNNonlinearSolver_Type SUNNonlinSolGetType(SUNNonlinearSolver NLS)

   This *required* function returns the nonlinear solver type.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.

   **Return value:**
      The :c:enum:`SUNNonlinearSolver_Type` type identifier for the nonlinear
      solver.


.. c:function:: SUNErrCode SUNNonlinSolInitialize(SUNNonlinearSolver NLS)

   This *optional* function handles nonlinear solver initialization
   and may perform any necessary memory allocations.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.

   **Return value:**
      A :c:type:`SUNErrCode`.

   **Notes:**
      It is assumed all solver-specific options have been set
      prior to calling :c:func:`SUNNonlinSolInitialize`. SUNNonlinSol
      implementations that do not require initialization may set this
      operation to ``NULL``.


.. c:function:: SUNErrCode SUNNonlinSolSetup(SUNNonlinearSolver NLS, N_Vector y, void* mem)

   This *optional* function performs any solver setup needed for a nonlinear solve.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *y* -- the initial guess passed to the nonlinear solver.
      * *mem* -- the SUNDIALS integrator memory structure.

   **Return value:**
      A :c:type:`SUNErrCode`.

   **Notes:**
      SUNDIALS integrators call :c:func:`SUNNonlinSolSetup` before
      each step attempt. SUNNonlinSol implementations that do not
      require setup may set this operation to ``NULL``.


.. c:function:: int SUNNonlinSolSolve(SUNNonlinearSolver NLS, N_Vector y0, N_Vector ycor, N_Vector w, sunrealtype tol, sunbooleantype callLSetup, void *mem)

   This *required* function solves the nonlinear system
   :math:`F(y)=0` or :math:`G(y)=y`.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *y0* -- the predicted value for the new solution state. This
        *must* remain unchanged throughout the solution process.
      * *ycor* -- on input the initial guess for the correction to the predicted
        state (zero) and on output the final correction to the predicted
        state.
      * *w* -- the solution error weight vector used for computing weighted error norms.
      * *tol* -- the requested solution tolerance in the weighted root-mean-squared norm.
      * *callLSetup* -- a flag indicating that the integrator
        recommends for the linear solver setup function to be called.
      * *mem* -- the SUNDIALS integrator memory structure.

   **Return value:**
      The return value is zero for a successful solve, a positive value
      for a recoverable error (i.e., the solve failed and the integrator
      should reduce the step size and reattempt the step), and a negative
      value for an unrecoverable error (i.e., the solve failed the and
      the integrator should halt and return an error to the user).


.. c:function:: SUNErrCode SUNNonlinSolFree(SUNNonlinearSolver NLS)

   This *optional* function frees any memory allocated by the
   nonlinear solver.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.

   **Return value:**
      * A :c:type:`SUNErrCode`




.. _SUNNonlinSol.API.SetFn:

SUNNonlinearSolver "set" functions
-------------------------------------

The following functions are used to supply nonlinear solver modules with
functions defined by the SUNDIALS integrators and to modify solver
parameters. Only the routine for setting the nonlinear system defining function
(:c:func:`SUNNonlinSolSetSysFn`) is required. All other set functions are optional.


.. c:function:: SUNErrCode SUNNonlinSolSetOptions(SUNNonlinearSolver NLS, const char* NLSid, const char* file_name, int argc, char* argv[])

   This *optional* routine sets SUNNonlinearSolver options from an array of strings or a file.

   :param NLS: the :c:type:`SUNNonlinearSolver` object.
   :param NLSid: the prefix for options to read. The default is "sunnonlinearsolver".
   :param file_name: the name of a file containing options to read. If this is
                     ``NULL`` or an empty string, ``""``, then no file is read.
   :param argc: length of the ``argv`` array.
   :param argv: an array of strings containing the options to set and their values.

   :return: :c:type:`SUNErrCode` indicating success or failure.

   .. note::

      The ``argc`` and ``argv`` arguments are typically those supplied to the user's
      ``main`` routine however, this is not required. The inputs are left unchanged by
      :c:func:`SUNNonlinSolSetOptions`.

      If the ``NLSid`` argument is ``NULL``, then the default prefix, ``sunnonlinearsolver``, must
      be used for all SUNNonlinearSolver options.  Whether ``NLSid`` is supplied or not, a ``"."``
      must be used to separate an option key from the prefix.  For example, when
      using the default ``NLSid``, the option ``sunnonlinearsolver.max_iters`` followed by the value
      can be used to set the maximum number of nonlinear solver iterations.
      When using a combination of SUNNonlinearSolver objects (e.g., when using MRIStep), it is recommended
      that users call :c:func:`SUNNonlinSolSetOptions` for each nonlinear solver using distinct
      *NLSid* inputs, so that each solver object can be configured separately.

      SUNNonlinearSolver options set via command-line arguments to
      :c:func:`SUNNonlinSolSetOptions` will overwrite any previously-set values.
      Options are set in the order they are given in ``argv`` and, if an
      option with the same prefix appears multiple times in ``argv``, the value of the
      last occurrence will used.

      The supported option names are noted within the documentation for the
      corresponding "set" function. For options that take a :c:type:`sunbooleantype` as
      input, use ``1`` to indicate ``true`` and ``0`` for ``false``.

   .. warning::

      This function is not available in the Fortran interface.

      File-based options are not yet supported, so the ``file_name`` argument
      should be set to either ``NULL`` or the empty string ``""``.

   .. versionadded:: 7.5.0


.. c:function:: SUNErrCode SUNNonlinSolSetSysFn(SUNNonlinearSolver NLS, SUNNonlinSolSysFn SysFn)

   This *required* function is used to provide the nonlinear solver
   with the function defining the nonlinear system. This is the function
   :math:`F(y)` in :math:`F(y)=0` for ``SUNNONLINEARSOLVER_ROOTFIND`` modules or
   :math:`G(y)` in :math:`G(y)=y` for ``SUNNONLINEARSOLVER_FIXEDPOINT`` modules.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *SysFn* -- the function defining the nonlinear system. See
        :numref:`SUNNonlinSol.API.SUNSuppliedFn` for the definition of
        :c:type:`SUNNonlinSolSysFn`.

   **Return value:**
      * A :c:type:`SUNErrCode`


.. c:function:: SUNErrCode SUNNonlinSolSetLSetupFn(SUNNonlinearSolver NLS, SUNNonlinSolLSetupFn SetupFn)

   This *optional* function is called by SUNDIALS integrators to provide
   the nonlinear solver with access to its linear solver setup function.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *SetupFn* -- a wrapper function to the SUNDIALS integrator's linear solver setup
        function. See :numref:`SUNNonlinSol.API.SUNSuppliedFn`  for the
        definition of :c:type:`SUNNonlinSolLSetupFn`.

   **Return value:**
      * A :c:type:`SUNErrCode`

   **Notes:**
      The :c:type:`SUNNonlinSolLSetupFn` function sets up the
      linear system :math:`Ax=b` where :math:`A = \frac{\partial
      F}{\partial y}` is the linearization of the nonlinear residual
      function :math:`F(y) = 0` (when using SUNLinSol direct linear
      solvers) or calls the user-defined preconditioner setup function
      (when using SUNLinSol iterative linear solvers). SUNNonlinSol
      implementations that do not require solving this system, do not
      utilize SUNLinSol linear solvers, or use SUNLinSol linear solvers
      that do not require setup may set this operation to ``NULL``.



.. c:function:: SUNErrCode SUNNonlinSolSetLSolveFn(SUNNonlinearSolver NLS, SUNNonlinSolLSolveFn SolveFn)

   This *optional* function is called by SUNDIALS integrators to provide
   the nonlinear solver with access to its linear solver solve function.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *SolveFn* -- a wrapper function to the SUNDIALS integrator's
        linear solver solve function. See
        :numref:`SUNNonlinSol.API.SUNSuppliedFn` for the definition of
        :c:type:`SUNNonlinSolLSolveFn`.

   **Return value:**
      * A :c:type:`SUNErrCode`

   **Notes:**
      The :c:type:`SUNNonlinSolLSolveFn` function solves the
      linear system :math:`Ax=b` where :math:`A = \frac{\partial
      F}{\partial y}` is the linearization of the nonlinear residual
      function :math:`F(y) = 0`.  SUNNonlinSol implementations that do
      not require solving this system or do not use SUNLinSol linear
      solvers may set this operation to ``NULL``.



.. c:function:: SUNErrCode SUNNonlinSolSetConvTestFn(SUNNonlinearSolver NLS, SUNNonlinSolConvTestFn CTestFn, void* ctest_data)

   This *optional* function is used to provide the nonlinear solver
   with a function for determining if the nonlinear solver iteration
   has converged. This is typically called by SUNDIALS integrators to
   define their nonlinear convergence criteria, but may be replaced by
   the user.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *CTestFn* -- a SUNDIALS integrator's nonlinear solver
        convergence test function. See
        :numref:`SUNNonlinSol.API.SUNSuppliedFn` for the definition of
        :c:type:`SUNNonlinSolConvTestFn`.
      * *ctest_data* -- is a data pointer passed to *CTestFn* every time it is
        called.

   **Return value:**
      * A :c:type:`SUNErrCode`

   **Notes:**
      SUNNonlinSol implementations utilizing their own convergence test
      criteria may set this function to ``NULL``.



.. c:function:: SUNErrCode SUNNonlinSolSetMaxIters(SUNNonlinearSolver NLS, int maxiters)

   This *optional* function sets the maximum number of nonlinear solver
   iterations. This is typically called by SUNDIALS integrators to
   define their default iteration limit, but may be adjusted by the user.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *maxiters* -- the maximum number of nonlinear iterations.

   **Return value:**
      * A :c:type:`SUNErrCode`

   **Notes:**

      If supported by the SUNNonlinearSolver implementation, this routine will be called
      by :c:func:`SUNNonlinSolSetOptions` when using the key
      "NLSid.max_iters".



.. _SUNNonlinSol.API.GetFn:

SUNNonlinearSolver "get" functions
----------------------------------

The following functions allow SUNDIALS integrators to retrieve nonlinear
solver statistics. The routines to get the number of iterations in the most
recent solve (:c:func:`SUNNonlinSolGetNumIters`) and number of convergence failures
are optional. The routine to get the current nonlinear solver iteration
(:c:func:`SUNNonlinSolGetCurIter`) is required when using the convergence test
provided by the SUNDIALS integrator or when using an iterative SUNLinSol
linear solver module; otherwise :c:func:`SUNNonlinSolGetCurIter` is optional.


.. c:function:: SUNErrCode SUNNonlinSolGetNumIters(SUNNonlinearSolver NLS, long int *niters)

   This *optional* function returns the number of nonlinear solver iterations
   in the most recent solve. This is typically called by the SUNDIALS
   integrator to store the nonlinear solver statistics, but may also be
   called by the user.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *niters* -- the total number of nonlinear solver iterations.

   **Return value:**
      * A :c:type:`SUNErrCode`


.. c:function:: SUNErrCode SUNNonlinSolGetCurIter(SUNNonlinearSolver NLS, int *iter)

   This function returns the iteration index of the current nonlinear
   solve. This function is *required* when using SUNDIALS
   integrator-provided convergence tests or when using an iterative
   SUNLinSol linear solver module; otherwise it is *optional*.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *iter* -- the nonlinear solver iteration in the current solve
        starting from zero.

   **Return value:**
      * A :c:type:`SUNErrCode`


.. c:function:: SUNErrCode SUNNonlinSolGetNumConvFails(SUNNonlinearSolver NLS, long int *nconvfails)

   This *optional* function returns the number of nonlinear solver convergence
   failures in the most recent solve. This is typically called by the SUNDIALS
   integrator to store the nonlinear solver statistics, but may also be called
   by the user.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *nconvfails* -- the total number of nonlinear solver convergence failures.

   **Return value:**
      * A :c:type:`SUNErrCode`


.. _SUNNonlinSol.API.SUNSuppliedFn:

Functions provided by SUNDIALS integrators
--------------------------------------------

To interface with SUNNonlinSol modules, the SUNDIALS integrators
supply a variety of routines for evaluating the nonlinear system,
calling the SUNLinSol setup and solve functions, and testing the
nonlinear iteration for convergence.  These integrator-provided routines
translate between the user-supplied ODE or DAE systems and the generic
interfaces to the nonlinear or linear systems of equations that result
in their solution. The functions provided to a SUNNonlinSol
module have types defined in the header file
``sundials/sundials_nonlinearsolver.h``; these are also described below.


.. c:type:: int (*SUNNonlinSolSysFn)(N_Vector ycor, N_Vector F, void* mem)

   These functions evaluate the nonlinear system :math:`F(y)`
   for ``SUNNONLINEARSOLVER_ROOTFIND`` type modules or :math:`G(y)`
   for ``SUNNONLINEARSOLVER_FIXEDPOINT`` type modules. Memory
   for *F* must by be allocated prior to calling this function. The
   vector *ycor* will be left unchanged.

   **Arguments:**
      * *ycor* -- is the current correction to the predicted state at which the
        nonlinear system should be evaluated.
      * *F* -- is the output vector containing :math:`F(y)` or
        :math:`G(y)`, depending on the solver type.
      * *mem* -- is the SUNDIALS integrator memory structure.

   **Return value:**
      The return value is zero for a successful solve, a positive value for
      a recoverable error, and a negative value for an unrecoverable error.

   **Notes:**
      SUNDIALS integrators formulate nonlinear systems as a function of the
      correction to the predicted solution. On each call to the nonlinear system
      function the integrator will compute and store the current solution based on
      the input correction. Additionally, the residual will store the value of the
      ODE right-hand side function or DAE residual used in computing the nonlinear
      system. These stored values are then directly used in the integrator-supplied
      linear solver setup and solve functions as applicable.


.. c:type:: int (*SUNNonlinSolLSetupFn)(sunbooleantype jbad, sunbooleantype* jcur, void* mem)

   These functions are wrappers to the SUNDIALS integrator's function
   for setting up linear solves with SUNLinSol modules.

   **Arguments:**
      * *jbad* -- is an input indicating whether the nonlinear solver
        believes that :math:`A` has gone stale (``SUNTRUE``) or not (``SUNFALSE``).
      * *jcur* -- is an output indicating whether the routine has updated the
        Jacobian :math:`A` (``SUNTRUE``) or not (``SUNFALSE``).
      * *mem* -- is the SUNDIALS integrator memory structure.

   **Return value:**
      The return value is zero for a successful solve, a positive value for
      a recoverable error, and a negative value for an unrecoverable error.

   **Notes:**
      The :c:type:`SUNNonlinSolLSetupFn` function sets up the linear
      system :math:`Ax=b` where :math:`A = \frac{\partial F}{\partial y}`
      is the linearization of the nonlinear residual function
      :math:`F(y) = 0` (when using SUNLinSol direct linear solvers) or
      calls the user-defined preconditioner setup function (when using
      SUNLinSol iterative linear solvers). SUNNonlinSol implementations
      that do not require solving this system, do not utilize SUNLinSol
      linear solvers, or use SUNLinSol linear solvers that do not
      require setup may ignore these functions.

      As discussed in the description of :c:type:`SUNNonlinSolSysFn`, the linear
      solver setup function assumes that the nonlinear system function has been
      called prior to the linear solver setup function as the setup will utilize
      saved values from the nonlinear system evaluation (e.g., the updated
      solution).


.. c:type:: int (*SUNNonlinSolLSolveFn)(N_Vector b, void* mem)

   These functions are wrappers to the SUNDIALS integrator's function
   for solving linear systems with SUNLinSol modules.

   **Arguments:**
      * *b* -- contains the right-hand side vector for the linear
        solve on input and the solution to the linear system on output.
      * *mem* -- is the SUNDIALS integrator memory structure.

   **Return value:**
      The return value is zero for a successful solve, a positive value for
      a recoverable error, and a negative value for an unrecoverable error.

   **Notes:**
      The :c:type:`SUNNonlinSolLSolveFn` function solves the linear
      system :math:`Ax=b` where :math:`A = \frac{\partial F}{\partial y}`
      is the linearization of the nonlinear residual function
      :math:`F(y) = 0`. SUNNonlinSol implementations that do not
      require solving this system or do not use SUNLinSol linear solvers
      may ignore these functions.

      As discussed in the description of :c:type:`SUNNonlinSolSysFn`, the linear
      solver solve function assumes that the nonlinear system function has been
      called prior to the linear solver solve function as the setup may utilize
      saved values from the nonlinear system evaluation (e.g., the updated
      solution).


.. c:type:: int (*SUNNonlinSolConvTestFn)(SUNNonlinearSolver NLS, N_Vector ycor, N_Vector del, sunrealtype tol, N_Vector ewt, void* ctest_data)

   These functions are SUNDIALS integrator-specific convergence tests for
   nonlinear solvers and are typically supplied by each SUNDIALS integrator,
   but users may supply custom problem-specific versions as desired.

   **Arguments:**
      * *NLS* -- is the SUNNonlinSol object.
      * *ycor* -- is the current correction (nonlinear iterate).
      * *del* -- is the difference between the current and prior nonlinear iterates.
      * *tol* -- is the nonlinear solver tolerance.
      * *ewt* -- is the weight vector used in computing weighted norms.
      * *ctest_data* -- is the data pointer provided to
        :c:func:`SUNNonlinSolSetConvTestFn()`.

   **Return value:**
      The return value of this routine will be a negative value if an
      unrecoverable error occurred or one of the following:

      * ``SUN_SUCCESS`` -- the iteration is converged.

      * ``SUN_NLS_CONTINUE`` -- the iteration has not converged, keep
        iterating.

      * ``SUN_NLS_CONV_RECVR`` -- the iteration appears to be
        diverging, try to recover.

   **Notes:**
      The tolerance passed to this routine by SUNDIALS integrators is
      the tolerance in a weighted root-mean-squared norm with error
      weight vector ``ewt``.  SUNNonlinSol modules utilizing their
      own convergence criteria may ignore these functions.



.. _SUNNonlinSol.API.ReturnCodes:

SUNNonlinearSolver return codes
---------------------------------

The functions provided to SUNNonlinSol modules by each SUNDIALS
integrator, and functions within the SUNDIALS-provided SUNNonlinSol
implementations, utilize a common set of return codes shown in
:numref:`SUNNonlinSol.API.CodeTable`.  Here, negative values correspond to non-recoverable
failures, positive values to recoverable failures, and zero to a
successful call.

.. _SUNNonlinSol.API.CodeTable:
.. table:: Description of the ``SUNNonlinearSolver`` return codes.
   :align: center

   +-----------------------+---------+---------------------------------------------------------------+
   | Name                  | Value   | Description                                                   |
   +=======================+=========+===============================================================+
   | SUN_SUCCESS           |    0    | successful call or converged solve                            |
   +-----------------------+---------+---------------------------------------------------------------+
   | SUN_NLS_CONTINUE      |  901    | the nonlinear solver is not converged, keep iterating         |
   +-----------------------+---------+---------------------------------------------------------------+
   | SUN_NLS_CONV_RECVR    |  902    | the nonlinear solver appears to be diverging, try to recover  |
   +-----------------------+---------+---------------------------------------------------------------+



.. _SUNNonlinSol.API.Generic:

The generic SUNNonlinearSolver module
-----------------------------------------

SUNDIALS integrators interact with specific SUNNonlinSol
implementations through the generic SUNNonlinSol module on which all
other SUNNonlinSol implementations are built. The
``SUNNonlinearSolver`` type is a pointer to a structure containing an
implementation-dependent *content* field and an *ops*
field.

A :c:type:`SUNNonlinearSolver` is a pointer to the
:c:struct:`_generic_SUNNonlinearSolver` structure:

.. c:type:: struct _generic_SUNNonlinearSolver *SUNNonlinearSolver

.. c:struct:: _generic_SUNNonlinearSolver

   The structure defining the SUNDIALS nonlinear solver class.

   .. c:member:: void* content

      Pointer to nonlinear solver-specific member data

   .. c:member:: SUNNonlinearSolver_Ops ops

      A virtual table of nonlinear solver operations provided by a specific
      implementation

   .. c:member:: SUNContext sunctx

      The SUNDIALS simulation context

The virtual table structure is defined as

.. c:type:: struct _generic_SUNNonlinearSolver_Ops *SUNNonlinearSolver_Ops

.. c:struct:: _generic_SUNNonlinearSolver_Ops

   The structure defining :c:type:`SUNNonlinearSolver` operations.

   .. c:member:: SUNNonlinearSolver_Type (*gettype)(SUNNonlinearSolver)

      The function implementing :c:func:`SUNNonlinSolGetType`

   .. c:member:: int (*initialize)(SUNNonlinearSolver)

      The function implementing :c:func:`SUNNonlinSolInitialize`

   .. c:member:: int (*setup)(SUNNonlinearSolver, N_Vector, void*)

      The function implementing :c:func:`SUNNonlinSolSetup`

   .. c:member:: int (*solve)(SUNNonlinearSolver, N_Vector, N_Vector, N_Vector, sunrealtype, sunbooleantype, void*)

      The function implementing :c:func:`SUNNonlinSolSolve`

   .. c:member:: int (*free)(SUNNonlinearSolver)

      The function implementing :c:func:`SUNNonlinSolFree`

   .. c:member:: int (*setsysfn)(SUNNonlinearSolver, SUNNonlinSolSysFn)

      The function implementing :c:func:`SUNNonlinSolSetSysFn`

   .. c:member:: int (*setlsetupfn)(SUNNonlinearSolver, SUNNonlinSolLSetupFn)

      The function implementing :c:func:`SUNNonlinSolSetLSetupFn`

   .. c:member:: int (*setlsolvefn)(SUNNonlinearSolver, SUNNonlinSolLSolveFn)

      The function implementing :c:func:`SUNNonlinSolSetLSolveFn`

   .. c:member:: int (*setctestfn)(SUNNonlinearSolver, SUNNonlinSolConvTestFn, void*)

      The function implementing :c:func:`SUNNonlinSolSetConvTestFn`

   .. c:member:: int (*setmaxiters)(SUNNonlinearSolver, int)

      The function implementing :c:func:`SUNNonlinSolSetMaxIters`

   .. c:member:: int (*getnumiters)(SUNNonlinearSolver, long int*)

      The function implementing :c:func:`SUNNonlinSolGetNumIters`

   .. c:member:: int (*getcuriter)(SUNNonlinearSolver, int*)

      The function implementing :c:func:`SUNNonlinSolGetCurIter`

   .. c:member:: int (*getnumconvfails)(SUNNonlinearSolver, long int*)

      The function implementing :c:func:`SUNNonlinSolGetNumConvFails`


The generic SUNNonlinSol module defines and implements the nonlinear
solver operations defined in
:numref:`SUNNonlinSol.API.CoreFn`--:numref:`SUNNonlinSol.API.GetFn`.
These routines are in fact only wrappers to the nonlinear solver
operations provided by a particular SUNNonlinSol implementation,
which are accessed through the ops field of the ``SUNNonlinearSolver``
structure. To illustrate this point we show below the implementation
of a typical nonlinear solver operation from the generic SUNNonlinSol
module, namely :c:func:`SUNNonlinSolSolve`, which solves the nonlinear
system and returns a flag denoting a successful or failed solve:

.. code-block:: c

   int SUNNonlinSolSolve(SUNNonlinearSolver NLS,
                         N_Vector y0, N_Vector y,
                         N_Vector w, sunrealtype tol,
                         sunbooleantype callLSetup, void* mem)
   {
     return((int) NLS->ops->solve(NLS, y0, y, w, tol, callLSetup, mem));
   }



.. _SUNNonlinSol.API.Custom:

Implementing a Custom SUNNonlinearSolver Module
--------------------------------------------------

A SUNNonlinSol implementation *must* do the following:

* Specify the content of the SUNNonlinSol module.

* Define and implement the required nonlinear solver operations defined
  in :numref:`SUNNonlinSol.API.CoreFn`--:numref:`SUNNonlinSol.API.GetFn`.
  Note that the names of the module routines should be unique to that
  implementation in order to permit using more than one SUNNonlinSol
  module (each with different ``SUNNonlinearSolver`` internal data
  representations) in the same code.

* Define and implement a user-callable constructor to create a
  ``SUNNonlinearSolver`` object.

To aid in the creation of custom ``SUNNonlinearSolver`` modules, the generic
``SUNNonlinearSolver`` module provides the utility functions
:c:func:`SUNNonlinSolNewEmpty` and :c:func:`SUNNonlinSolFreeEmpty`. When used
in custom ``SUNNonlinearSolver`` constructors these functions will ease the
introduction of any new optional nonlinear solver operations to the
``SUNNonlinearSolver`` API by ensuring that only required operations need to
be set.

.. c:function:: SUNNonlinearSolver SUNNonlinSolNewEmpty(SUNContext sunctx)

  This function allocates a new generic ``SUNNonlinearSolver`` object and
  initializes its content pointer and the function pointers in the operations
  structure to ``NULL``.

  **Return value:**
     If successful, this function returns a ``SUNNonlinearSolver`` object.
     If an error occurs when allocating the object, then this routine will
     return ``NULL``.

.. c:function:: void SUNNonlinSolFreeEmpty(SUNNonlinearSolver NLS)

  This routine frees the generic ``SUNNonlinearSolver`` object, under the assumption that any
  implementation-specific data that was allocated within the underlying content structure
  has already been freed. It will additionally test whether the ops pointer is ``NULL``,
  and, if it is not, it will free it as well.

   **Arguments:**
      * *NLS* -- a SUNNonlinearSolver object


Additionally, a ``SUNNonlinearSolver`` implementation *may* do
the following:

* Define and implement additional user-callable "set" routines
  acting on the ``SUNNonlinearSolver`` object, e.g., for setting
  various configuration options to tune the performance of the
  nonlinear solve algorithm.

* Provide additional user-callable "get" routines acting on the
  ``SUNNonlinearSolver`` object, e.g., for returning various solve
  statistics.
