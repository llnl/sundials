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

.. _ARKODE.Usage.Headers:

Access to library and header files
==================================

At this point, it is assumed that the installation of ARKODE, following the
procedure described in :numref:`Installation`, has been completed successfully.
In the proceeding text, the directories ``libdir`` and ``incdir`` are the
installation library and include directories, respectively. For a default
installation, these are ``instdir/lib`` and ``instdir/include``, respectively,
where ``instdir`` is the directory where SUNDIALS was installed.

Regardless of where the user's application program resides, its
associated compilation and load commands must make reference to the
appropriate locations for the library and header files required by
ARKODE. ARKODE symbols are found in ``libdir/libsundials_arkode.lib``. 
Thus, in addition to linking to ``libdir/libsundials_core.lib``, ARKODE
users need to link to the ARKODE library. Symbols for additional SUNDIALS
modules, vectors and algebraic solvers, are found in

.. code-block::

  <libdir>/libsundials_nvec*.lib
  <libdir>/libsundials_sunmat*.lib
  <libdir>/libsundials_sunlinsol*.lib
  <libdir>/libsundials_sunnonlinsol*.lib
  <libdir>/libsundials_sunmem*.lib

The file extension ``.lib`` is typically ``.so`` for shared libraries 
and ``.a`` for static libraries.  

The relevant header files for ARKODE are located in the subdirectories
``incdir/include/arkode``. To use ARKODE the application needs to include 
the header file(s) for the ARKODE time-stepper(s) of choice in addition
to the SUNDIALS core header file. 

.. code:: c

  #include <sundials/sundials_core.h>      // Provides core SUNDIALS types
  #include <arkode/arkode_arkstep.h>       // ARKStep provides explicit, implicit, IMEX additive RK methods.
  #include <arkode/arkode_erkstep.h>       // ERKStep provides explicit RK methods.
  #include <arkode/arkode_forcingstep.h>   // ForcingStep provides a forcing method.
  #include <arkode/arkode_lsrkstep.h>      // LSRKStep provides low-storage RK methods.
  #include <arkode/arkode_mristep.h>       // MRIStep provides multirate RK methods.
  #include <arkode/arkode_splittingstep.h> // SplittingStep provides operator splitting methods.
  #include <arkode/arkode_sprkstep.h>      // SPRKStep provides symplectic partitioned RK methods.

Each of these define several types and various constants, include function
prototypes, and include the shared ``arkode/arkode.h`` and
``arkode/arkode_ls.h`` header files. No other header files are required to be
included, but there are optional ones that can be included as necessary.
Information on optional headers is given in the relevant documentation section.

The calling program must also include an :c:type:`N_Vector` implementation header file,  
of the form ``nvector/nvector_*.h``. See :numref:`NVectors` for the appropriate name.  

If the user includes a non-trivial implicit component to their ODE system in
ARKStep, or if the slow time scale for MRIStep should be treated implicitly,
then each implicit stage will require a nonlinear solver for the resulting
system of algebraic equations -- the default for this is a modified or inexact
Newton iteration, depending on the user's choice of linear solver.  If choosing
to set which nonlinear solver module, or when interacting with a
:c:type:`SUNNonlinearSolver` module directly, the calling program must also include a
:c:type:`SUNNonlinearSolver` header file, of the form ``sunnonlinsol/sunnonlinsol_***.h``
where ``***`` is the name of the nonlinear solver module 
(see :numref:`SUNNonlinSol` for more information). 


If using a nonlinear solver that requires the solution of a linear system of the
form :math:`\mathcal{A}x=b` (e.g., the default Newton iteration), then a linear
solver module header file will also be required.  Similarly, if the ODE system
in ARKStep involves a non-identity mass matrix :math:`M \ne I`, then each time
step will require a linear solver for systems of the form :math:`Mx=b`. In this
case it will be necessary to include the header file for a
:c:type:`SUNLinearSolver` solver, which is of the form
``sunlinsol/sunlinsol_***.h`` (see :numref:`SUNLinSol` for more information). 

If the linear solver is matrix-based, the linear solver header will also include a  
header file of the from ``sunmatrix/sunmatrix_*.h`` where ``*`` is the name of the  
matrix implementation compatible with the linear solver (see :numref:`SUNMatrix` for  
more information). 

Other headers may be needed, according to the choice of preconditioner, etc.
For example, if preconditioning for an iterative linear solver were performed
using the ARKBBDPRE module, the header ``arkode/arkode_bbdpre.h`` is needed to
access the preconditioner initialization routines.


