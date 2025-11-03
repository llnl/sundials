.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage:

************
Using ARKODE
************

This chapter discusses usage of ARKODE for the solution of initial value
problems (IVPs) in C, C++ and Fortran applications.  The chapter builds upon
:numref:`SUNDIALS`.  Unlike other packages in SUNDIALS, ARKODE provides an
infrastructure for one-step methods.  However, ARKODE's individual
time-stepping methods, including definition of the IVP itself, are handled
by time-stepping modules that sit on top of ARKODE.  While most of the
routines to use ARKODE generally apply to all of its time-stepping modules,
some of these apply to only a subset of these "steppers," while others are
specific to a given stepper.

Thus, we organize this chapter as follows.  We first discuss commonalities
to each of ARKODE's time-stepping modules.  These commonalities include the
locations and naming conventions for the library and header files, data types
in SUNDIALS, the layout of the user's main program, and a variety of
user-callable and user-supplied functions.  For these user-callable routines,
we distinguish those that apply for only a subset of ARKODE's time-stepping
modules.  We then describe shared utilities that are supported by some of
ARKODE's time stepping modules, including "relaxation" methods and
preconitioners.  Following our discussion of these commonalities, we
separately discuss the usage details that that are specific to each of ARKODE's
time stepping modules: :ref:`ARKStep <ARKODE.Usage.ARKStep>`,
:ref:`ERKStep <ARKODE.Usage.ERKStep>`,
:ref:`ForcingStep <ARKODE.Usage.ForcingStep>`,
:ref:`LSRKStep <ARKODE.Usage.LSRKStep>`,
:ref:`MRIStep <ARKODE.Usage.MRIStep>`,
:ref:`SplittingStep <ARKODE.Usage.SplittingStep>`, and
:ref:`SPRKStep <ARKODE.Usage.SPRKStep>`.

ARKODE also uses various input and output constants; these are defined as
needed throughout this chapter, but for convenience the full list is provided
separately in :numref:`ARKODE.Constants`.

The example programs for ARKODE are located in the source code ``examples/arkode``
folder.  We note that these may be helpful as templates for new codes.  Users
with applications written in Fortran should see the chapter
:numref:`Fortran`, which describes the Fortran interfaces for
SUNDIALS, and we additionally include multiple Fortran example programs
in the ARKODE ``examples`` directory.

When solving problems with an implicit component, we note that not all
SUNLINSOL, SUNMATRIX, and preconditioning modules are compatible with
all NVECTOR implementations.  Details on compatibility are given in the
documentation for each SUNMATRIX (see :numref:`SUNMatrix`) and each
SUNLINSOL module (see :numref:`SUNLinSol`). For example, NVECTOR_PARALLEL
is not compatible with the dense, banded, or sparse SUNMATRIX types,
or with the corresponding dense, banded, or sparse SUNLINSOL modules.
Please check :numref:`SUNMatrix` and :numref:`SUNLinSol` to
verify compatibility between these modules.  In addition to that
documentation, we note that the ARKBANDPRE preconditioning module is
only compatible with the NVECTOR_SERIAL, NVECTOR_OPENMP or
NVECTOR_PTHREADS vector implementations, and the preconditioner module
ARKBBDPRE can only be used with NVECTOR_PARALLEL.


.. toctree::
   :maxdepth: 1

   General
   Skeleton
   User_callable
   User_supplied
   Relaxation
   Preconditioners
   ARKStep/index.rst
   ERKStep/index.rst
   ForcingStep/index.rst
   LSRKStep/index.rst
   MRIStep/index.rst
   SplittingStep/index.rst
   SPRKStep/index.rst
   ASA.rst
