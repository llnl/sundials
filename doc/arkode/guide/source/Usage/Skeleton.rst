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

.. _ARKODE.Usage.Skeleton:

A skeleton of the user's main program
============================================

The following is a skeleton of the user's main program (or calling
program) for the integration of an ODE IVP using ARKODE.
Most of the steps are independent of the NVECTOR, SUNMATRIX, SUNLINSOL
and SUNNONLINSOL implementations used.  For the steps that are not,
refer to :numref:`NVectors`, :numref:`SUNMatrix`,
:numref:`SUNLinSol`, and  :numref:`SUNNonlinSol` for the specific name of
the function to be called or macro to be referenced.

.. index:: User main program

#. Initialize parallel or multi-threaded environment, if appropriate.

   For example, call ``MPI_Init`` to initialize MPI if used, or set
   ``num_threads``, the number of threads to use within the threaded
   vector functions, if used.

#. Create the SUNDIALS simulation context object.

   Call :c:func:`SUNContext_Create` to allocate the ``SUNContext`` object.

#. Set problem dimensions, etc.

   This generally includes the problem size, ``N``, and may include
   the local vector length ``Nlocal``.

   .. note::

      The variables ``N`` and ``Nlocal`` should be of type
      ``sunindextype``.

#. Set vector of initial values

   To set the vector ``y0`` of initial values, use the appropriate
   functions defined by the particular NVECTOR implementation.

   For native SUNDIALS vector implementations (except the CUDA and
   RAJA based ones), use a call of the form

   .. code-block:: c

      y0 = N_VMake_***(..., ydata);

   if the ``sunrealtype`` array ``ydata`` containing the initial values of
   :math:`y` already exists.  Otherwise, create a new vector by making
   a call of the form

   .. code-block:: c

      y0 = N_VNew_***(...);

   and then set its elements by accessing the underlying data where it
   is located with a call of the form

   .. code-block:: c

      ydata = N_VGetArrayPointer_***(y0);

   For details on each of SUNDIALS' provided vector implementations, see
   the corresponding sections in :numref:`NVectors` for details.

#. Create ARKODE object

   Call a stepper-specific constructor, ``arkode_mem = *StepCreate(...)``, to
   create the ARKODE memory block. These routines return a ``void*`` pointer to
   this memory structure. See :numref:`ARKODE.Usage.ARKStep.Initialization`,
   :numref:`ARKODE.Usage.ERKStep.Initialization`,
   :numref:`ARKODE.Usage.MRIStep.Initialization`, or
   :numref:`ARKODE.Usage.SPRKStep.Initialization` for details.

#. Specify integration tolerances

   Call :c:func:`ARKodeSStolerances()` or
   :c:func:`ARKodeSVtolerances()` to specify either a scalar relative
   tolerance and scalar absolute tolerance, or a scalar relative
   tolerance and a vector of absolute tolerances,
   respectively.  Alternatively, call :c:func:`ARKodeWFtolerances()`
   to specify a function which sets directly the weights used in
   evaluating WRMS vector norms. See
   :numref:`ARKODE.Usage.Tolerances` for details.

   If a problem with non-identity mass matrix is used, and the
   solution units differ considerably from the equation units,
   absolute tolerances for the equation residuals (nonlinear and
   linear) may be specified separately through calls to
   :c:func:`ARKodeResStolerance()`, :c:func:`ARKodeResVtolerance()`, or
   :c:func:`ARKodeResFtolerance()`.

#. Create matrix object

   If a nonlinear solver requiring a linear solver will be used (e.g.,
   a Newton iteration) and the linear solver will be a matrix-based linear
   solver, then a template Jacobian matrix must be created by using the
   appropriate functions defined by the particular SUNMATRIX
   implementation.

   For the SUNDIALS-supplied SUNMATRIX implementations, the
   matrix object may be created using a call of the form

   .. code-block:: c

      SUNMatrix A = SUNBandMatrix(..., sunctx);

   or similar for the other matrix modules (see :numref:`SUNMatrix` for
   further information).

   Similarly, if the problem involves a non-identity mass matrix, and
   the mass-matrix linear systems will be solved using a direct linear
   solver, then a template mass matrix must be created by using the
   appropriate functions defined by the particular SUNMATRIX
   implementation.

#. Create linear solver object

   If a nonlinear solver requiring a linear solver will be used (e.g.,
   a Newton iteration), or if the problem involves a non-identity mass
   matrix, then the desired linear solver object(s) must be created by
   using the appropriate functions defined by the particular SUNLINSOL
   implementation.

   For any of the SUNDIALS-supplied SUNLINSOL implementations, the
   linear solver object may be created using a call of the form

   .. code-block:: c

      SUNLinearSolver LS = SUNLinSol_*(...);

   where ``*`` can be replaced with "Dense", "SPGMR", or other
   options, as discussed in :numref:`SUNLinSol`.

#. Set linear solver optional inputs

   Call ``*Set*`` functions from the selected linear solver module
   to change optional inputs specific to that linear solver.  See the
   documentation for each SUNLINSOL module in
   :numref:`SUNLinSol` for details.

#. Attach linear solver module

   If a linear solver was created above for implicit stage solves,
   initialize the ARKLS linear solver interface by attaching the
   linear solver object (and Jacobian matrix object, if applicable)
   with the call (for details see :numref:`ARKODE.Usage.LinearSolvers`):

   .. code-block:: c

      ier = ARKodeSetLinearSolver(...);

   Similarly, if the problem involves a non-identity mass matrix,
   initialize the ARKLS mass matrix linear solver interface by
   attaching the mass linear solver object (and mass matrix object,
   if applicable) with the call (for details see
   :numref:`ARKODE.Usage.LinearSolvers`):

   .. code-block:: c

      ier = ARKodeSetMassLinearSolver(...);

#. Create nonlinear solver object

   If the problem involves an implicit component, and if a non-default
   nonlinear solver object will be used for implicit stage solves
   (see :numref:`ARKODE.Usage.NonlinearSolvers`),
   then the desired nonlinear solver object must be created by using
   the appropriate functions defined by the particular SUNNONLINSOL
   implementation (e.g., ``NLS = SUNNonlinSol_***(...);`` where
   ``***`` is the name of the nonlinear solver (see
   :numref:`SUNNonlinSol` for details).

   For the SUNDIALS-supplied SUNNONLINSOL implementations, the
   nonlinear solver object may be created using a call of the form

   .. code-block:: c

      SUNNonlinearSolver NLS = SUNNonlinSol_*(...);

   where ``*`` can be replaced with "Newton", "FixedPoint", or other
   options, as discussed in :numref:`SUNNonlinSol`.

#. Attach nonlinear solver module

   If a nonlinear solver object was created above, then it must be
   attached to ARKODE using the call (for details see
   :numref:`ARKODE.Usage.NonlinearSolvers`):

   .. code-block:: c

      ier = ARKodeSetNonlinearSolver(...);

#. Set nonlinear solver optional inputs

   Call the appropriate set functions for the selected nonlinear
   solver module to change optional inputs specific to that nonlinear
   solver.  These *must* be called after attaching the nonlinear
   solver to ARKODE, otherwise the optional inputs will be
   overridden by ARKODE defaults.  See
   :numref:`SUNNonlinSol` for more information on optional inputs.

#. Set optional inputs

   Call ``ARKodeSet*`` functions to change any optional inputs that
   control the behavior of ARKODE from their default values. See
   :numref:`ARKODE.Usage.OptionalInputs` for details.

   Additionally, call ``*StepSet*`` routines to change any
   stepper-specific optional inputs from their default values.  See
   :numref:`ARKODE.Usage.ARKStep.OptionalInputs`,
   :numref:`ARKODE.Usage.ERKStep.OptionalInputs`,
   :numref:`ARKODE.Usage.MRIStep.OptionalInputs`, or
   :numref:`ARKODE.Usage.SPRKStep.OptionalInputs` for details.

#. Specify rootfinding problem

   Optionally, call :c:func:`ARKodeRootInit()` to initialize a rootfinding
   problem to be solved during the integration of the ODE system. See
   :numref:`ARKODE.Usage.RootFinding` for general details, and
   :numref:`ARKODE.Usage.OptionalInputs` for relevant optional
   input calls.

#. Advance solution in time

   For each point at which output is desired, call

   .. code-block:: c

      ier = ARKodeEvolve(arkode_mem, tout, yout, &tret, itask);

   Here, ``itask`` specifies the return mode. The vector ``yout``
   (which can be the same as the vector ``y0`` above) will contain
   :math:`y(t_\text{out})`. See
   :numref:`ARKODE.Usage.Integration` for details.

#. Get optional outputs

   Call ``ARKodeGet*`` functions to obtain optional output. See
   :numref:`ARKODE.Usage.OptionalOutputs` for details.

   Additionally, call ``*StepGet*`` routines to retrieve any
   stepper-specific optional outputs.  See
   :numref:`ARKODE.Usage.ARKStep.OptionalOutputs`,
   :numref:`ARKODE.Usage.ERKStep.OptionalOutputs`,
   :numref:`ARKODE.Usage.MRIStep.OptionalOutputs`, or
   :numref:`ARKODE.Usage.SPRKStep.OptionalOutputs` for details.

#. Deallocate memory for solution vector

   Upon completion of the integration, deallocate memory for the
   vector ``y`` (or ``yout``) by calling the destructor function:

   .. code-block:: c

      N_VDestroy(y);

#. Free solver memory

   Call :c:func:`ARKodeFree()` to free the memory allocated for
   the ARKODE module (and any nonlinear solver module).

#. Free linear solver and matrix memory

   Call :c:func:`SUNLinSolFree()` and (possibly)
   :c:func:`SUNMatDestroy()` to free any memory allocated for the
   linear solver and matrix objects created above.

#. Free nonlinear solver memory

   If a user-supplied ``SUNNonlinearSolver`` was provided to ARKODE,
   then call :c:func:`SUNNonlinSolFree()` to free any memory allocated
   for the nonlinear solver object created above.

#. Free the SUNContext object

   Call :c:func:`SUNContext_Free` to free the memory allocated for the ``SUNContext`` object.

#. Finalize MPI, if used

   Call ``MPI_Finalize`` to terminate MPI.
