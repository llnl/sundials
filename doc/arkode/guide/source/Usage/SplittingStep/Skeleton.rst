.. ----------------------------------------------------------------
   Programmer(s): Steven B. Roberts @ LLNL
   ----------------------------------------------------------------
   Based on MRIStep by David J. Gardner @ LLNL
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

.. _ARKODE.Usage.SplittingStep.Skeleton:

A skeleton of the user's main program
============================================

While SplittingStep usage generally follows the same pattern as the rest of
ARKODE, since it is the composition of other steppers, we summarize the
differences in using SplittingStep here.  Steps that are unchanged from the
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
