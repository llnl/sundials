.. ----------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
                  Steven B. Roberts @ LLNL
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

.. _SUNStepper.Implementing:

Implementing a SUNStepper
=========================

To create a SUNStepper implementation:

#. Define the stepper-specific content.

   This is typically a user-defined structure in C codes, a user-defined class
   or structure in C++ codes, or a user-defined module in Fortran codes. This
   content should hold any data necessary to perform the operations defined by
   the :c:type:`SUNStepper` member functions.

#. Define implementations of the required member functions (see
   :numref:`SUNStepper.Description.ImplMethods`).

   These are typically user-defined functions in C, member functions of the
   user-defined structure or class in C++, or functions contained in the
   user-defined module in Fortran.

   Note that all member functions are passed the :c:type:`SUNStepper` object and
   the stepper-specific content can, if necessary, be retrieved using
   :c:func:`SUNStepper_GetContent`. Stepper-specific warnings and errors can be
   recorded with :c:func:`SUNStepper_SetLastFlag`.

#. In the user code, before creating the outer memory structure that uses the
   :c:type:`SUNStepper`, e.g., with :c:func:`SplittingStepCreate` or
   :c:func:`ForcingStepCreate`, do the following:

   #. Create a :c:type:`SUNStepper` object with :c:func:`SUNStepper_Create`.

   #. Attach a pointer to the stepper content to the :c:type:`SUNStepper` object
      with :c:func:`SUNStepper_SetContent` if necessary, e.g., when the content
      is a C structure.

   #. Attach the member function implementations using the functions described
      in :numref:`SUNStepper.Description.BaseMethods.AttachFunctions`.

#. Attach the :c:type:`SUNStepper` object to the outer memory structure, e.g.,
   with :c:func:`SplittingStepCreate` or :c:func:`ForcingStepCreate`.
