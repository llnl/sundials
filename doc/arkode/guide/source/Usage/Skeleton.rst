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

.. _ARKODE.Usage.ARKStep.Skeleton:

ARKStep, ERKStep, LSRKStep, and SPRKStep
----------------------------------------

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


.. _ARKODE.Usage.MRIStep.Skeleton:

MRIStep
-------

While MRIStep usage generally follows the same pattern as the rest of ARKODE,
since it involves the solution of both MRIStep for the slow time scale and
another time integrator for the fast time scale, we summarize the differences in
using MRIStep here.  Steps that are unchanged from the skeleton program
presented in :numref:`ARKODE.Usage.Skeleton` are *italicized*.

.. index:: MRIStep user main program

#. *Initialize parallel or multi-threaded environment, if appropriate.*

#. *Create the SUNDIALS simulation context object*

#. *Set problem dimensions, etc.*

#. *Set vector of initial values*

#. Create an inner stepper object to solve the fast (inner) IVP

   * If using an ARKODE stepper module for the fast integrator, create and configure
     the stepper as normal following the steps detailed in the section for the desired
     stepper.

     Once the ARKODE stepper object is setup, create an ``MRIStepInnerStepper`` object
     with :c:func:`ARKodeCreateMRIStepInnerStepper`.

   * If supplying a user-defined fast (inner) integrator, create the
     ``MRIStepInnerStepper`` object as described in section
     :numref:`ARKODE.Usage.MRIStep.CustomInnerStepper`.

   .. note::

      When using ARKStep as a fast (inner) integrator it is the user's
      responsibility to create, configure, and attach the integrator to the
      MRIStep module. User-specified options regarding how this fast integration
      should be performed (e.g., adaptive vs. fixed time step,
      explicit/implicit/ImEx partitioning, algebraic solvers, etc.) will be
      respected during evolution of the fast time scale during MRIStep
      integration.

      Due to the algorithms supported in MRIStep, the ARKStep module used for
      the fast time scale must be configured with an identity mass matrix.

      If a *user_data* pointer needs to be passed to user functions called by
      the fast (inner) integrator then it should be attached here by calling
      :c:func:`ARKodeSetUserData()`. This *user_data* pointer will only be
      passed to user-supplied functions that are attached to the fast (inner)
      integrator. To supply a *user_data* pointer to user-supplied functions
      called by the slow (outer) integrator the desired pointer should be
      attached by calling :c:func:`ARKodeSetUserData()` after creating the
      MRIStep memory below. The *user_data* pointers attached to the inner and
      outer integrators may be the same or different depending on what is
      required by the user code.

      Specifying a rootfinding problem for the fast integration is not
      supported. Rootfinding problems should be created and initialized with
      the slow integrator. See the steps below and :c:func:`ARKodeRootInit()`
      for more details.

#. Create an MRIStep object for the slow (outer) integration

   Create the MRIStep object by calling  :c:func:`MRIStepCreate`. One of the
   inputs to :c:func:`MRIStepCreate` is the ``MRIStepInnerStepper`` object for
   solving the fast (inner) IVP created in the previous step.

#. If using fixed step sizes, then set the slow step size by calling
   :c:func:`ARKodeSetFixedStep()` on the MRIStep object to specify the
   slow time step size.

   If using adaptive slow steps, then specify the desired integration tolerances
   as normal.  By default, MRIStep will use a "decoupled" (see
   :numref:`ARKODE.Mathematics.MultirateControllers`) I controller (see
   :numref:`SUNAdaptController.Soderlind`),  Alternately, create and attach a
   multirate temporal controller (see :numref:`SUNAdaptController.MRIHTol`).

#. Create and configure implicit solvers (*as appropriate*)

   Specifically, if MRIStep is configured with an implicit slow right-hand side
   function in the prior step, then the following steps are recommended:

   #. *Specify integration tolerances*

   #. *Create matrix object*

   #. *Create linear solver object*

   #. *Set linear solver optional inputs*

   #. *Attach linear solver module*

   #. *Create nonlinear solver object*

   #. *Attach nonlinear solver module*

   #. *Set nonlinear solver optional inputs*

#. *Set optional inputs*

#. *Specify rootfinding problem*

#. *Advance solution in time*

#. *Get optional outputs*

#. *Deallocate memory for solution vector*

#. Free solver memory

   * If ARKStep was used as the fast (inner) IVP integrator, call
     :c:func:`MRIStepInnerStepper_Free` and :c:func:`ARKodeFree` to free the
     memory allocated for the fast (inner) integrator.

   * If a user-defined fast (inner) integrator was supplied, free the integrator
     content and call :c:func:`MRIStepInnerStepper_Free` to free the
     ``MRIStepInnerStepper`` object.

   * Call :c:func:`ARKodeFree` to free the memory allocated for the MRIStep
     slow integration object.

#. *Free linear solver and matrix memory (as appropriate)*

#. *Free nonlinear solver memory (as appropriate)*

#. *Free the SUNContext object*

#. *Finalize MPI, if used*


.. _ARKODE.Usage.MRIStep.Skeleton-ExtSTS:

MRIStep ExtSTS
--------------

When using MRIStep to construct an extended super time stepping (ExtSTS) method,
its usage differs somewhat from the standard MRIStep usage.  Namely, since
ExtSTS methods are a special case of MRIStep methods, but where the inner
integrator consists of a single step of a STS method, some of the "setup" for
MRIStep is simplified.  Steps that are unchanged from the skeleton program
presented above are *italicized*.

.. index:: MRIStep ExtSTS method user main program

#. *Initialize parallel or multi-threaded environment, if appropriate.*

#. *Create the SUNDIALS simulation context object*

#. *Set problem dimensions, etc.*

#. *Set vector of initial values*

#. Create the MRIStep object for ExtSTS methods by calling
   :c:func:`MRIStepCreateExtSTS`.

#. Call :c:func:`MRIStepGetSTS` to retrieve the LSRKStep component
   from inside the ExtSTS solver object, and provide either a dominant
   eigenvalue estimation function, or a SUNDomEigEstimator object, to the
   LSRKStep component by calling either :c:func:`LSRKStepSetDomEigFn` or
   :c:func:`LSRKStepSetDomEigEstimator`.  This is required for the LSRKStep
   component to determine the number of stages to use for each inner step.

#. Set optional inputs

   * If an implicit partition is present, or if the user wishes to adjust
     the configuration of the MRI aspects of the ExtSTS solver, they may
     configure these by calling the relevant MRIStep configuration
     functions directly on the MRIStep object, e.g.,
     :c:func:`MRIStepSetCoupling` or :c:func:`ARKodeSetNonlinearSolver`.

   * If the user wishes to adjust other aspects of the configuration of
     the LSRKStep component inside the ExtSTS solver object, call the
     relevant ``ARKodeSet`` or ``LSRKStepSet`` functions on the object returned
     by :c:func:`MRIStepGetSTS` above.

#. *Specify rootfinding problem*

#. *Advance solution in time*

#. *Get optional outputs*

#. *Deallocate memory for solution vector*

#. Free solver memory

   * Call :c:func:`ARKodeFree` to free the memory allocated for the MRIStep
     slow integration object.  This will automatically free the inner LSRKStep
     object as well.

#. *Free the SUNContext object*

#. *Finalize MPI, if used*

.. _ARKODE.Usage.SplittingStep.Skeleton:

SplittingStep
-------------

While SplittingStep usage generally follows the same pattern as the rest of
ARKODE, since it is the composition of other steppers, we summarize the
differences in using SplittingStep here. Steps that are unchanged from the
skeleton program presented in :numref:`ARKODE.Usage.Skeleton` are *italicized*.

.. index:: SplittingStep user main program

#. *Initialize parallel or multi-threaded environment, if appropriate.*

#. *Create the SUNDIALS simulation context object*

#. *Set problem dimensions, etc.*

#. *Set vector of initial values*

#. Create a stepper object for each problem partition

   * If using an ARKODE stepper module as an partition integrator, create and
     configure the stepper as normal following the steps detailed in the section
     for the desired stepper.

     Once the ARKODE stepper object is setup, create a :c:type:`SUNStepper`
     object with :c:func:`ARKodeCreateSUNStepper`.

   * If supplying a user-defined partition integrator, create the
     :c:type:`SUNStepper` object as described in section
     :numref:`SUNStepper.Implementing`.

   .. note::

      When using ARKODE for partition integrators, it is the user's
      responsibility to create and configure the integrator. User-specified
      options regarding how the integration should be performed (e.g., adaptive
      vs. fixed time step, explicit/implicit/ImEx partitioning, algebraic
      solvers, etc.) will be respected during evolution of a partition during
      SplittingStep integration.

      If a *user_data* pointer needs to be passed to user functions called by
      a partition integrator then it should be attached to the partition integrator here by calling
      :c:func:`ARKodeSetUserData`. This *user_data* pointer will only be
      passed to user-supplied functions that are attached to a partition
      integrator. To supply a *user_data* pointer to user-supplied functions
      called by the SplittingStep integrator, the desired pointer should be
      attached by calling :c:func:`ARKodeSetUserData` after creating the
      SplittingStep memory below. The *user_data* pointers attached to the
      partition and SplittingStep integrators may be the same or different
      depending on what is required by the user code.

      Specifying a rootfinding problem for a partition integrator is not
      supported. Rootfinding problems should be created and initialized with
      SplittingStep. See the steps below and :c:func:`ARKodeRootInit` for more
      details.

#. Create a SplittingStep object

   Create the SplittingStep object by calling :c:func:`SplittingStepCreate`. One
   of the inputs to :c:func:`SplittingStepCreate` is an array of
   :c:type:`SUNStepper` objects with one to evolve each partition.

#. Set the SplittingStep step size

   Call :c:func:`ARKodeSetFixedStep()` on the SplittingStep object to specify
   the overall time step size.

#. *Set optional inputs*

#. *Specify rootfinding problem*

#. *Advance solution in time*

#. *Get optional outputs*

#. *Deallocate memory for solution vector*

#. Free solver memory

   * If an ARKODE stepper module was used as a partition IVP integrator, call
     :c:func:`SUNStepper_Destroy` and :c:func:`ARKodeFree` to free the memory
     allocated for that integrator.

   * If a user-defined partition integrator was supplied, free the integrator
     content and call :c:func:`SUNStepper_Destroy` to free the :c:type:`SUNStepper`
     object.

   * Call :c:func:`ARKodeFree` to free the memory allocated for the
     SplittingStep integration object.

#. *Free the SUNContext object*

#. *Finalize MPI, if used*
