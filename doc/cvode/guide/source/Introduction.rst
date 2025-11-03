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

.. _CVODE.Introduction:

************
Introduction
************

CVODE is part of a software family called SUNDIALS: SUite of
Nonlinear and DIfferential/ALgebraic equation
Solvers :cite:p:`HBGLSSW:05`. This suite consists of CVODE,
ARKODE, KINSOL, and IDA, and variants of these with sensitivity
analysis capabilities.

.. _CVODE.Introduction.history:

Historical Background
=====================

Fortran solvers for ODE initial value problems are widespread and heavily
used. Two solvers that have been written at LLNL in the past are
VODE :cite:p:`BBH:89` and
VODPK :cite:p:`Byr:92`. VODE is a general purpose
solver that includes methods for both stiff and nonstiff systems, and in
the stiff case uses direct methods (full or banded) for the solution of
the linear systems that arise at each implicit step. Externally,
VODE is very similar to the well known solver
LSODE :cite:p:`RaHi:94`. VODPK is a variant of VODE
that uses a preconditioned Krylov (iterative) method, namely GMRES, for
the solution of the linear systems. VODPK is a powerful tool for
large stiff systems because it combines established methods for stiff
integration, nonlinear iteration, and Krylov (linear) iteration with a
problem-specific treatment of the dominant source of stiffness, in the
form of the user-supplied preconditioner
matrix :cite:p:`BrHi:89`. The capabilities of both VODE
and VODPK have been combined in the C-language package
CVODE :cite:p:`CoHi:96`.

At present, CVODE may utilize a variety of Krylov methods provided
in SUNDIALS that can be used in conjunction with Newton iteration:
these include the GMRES (Generalized Minimal
RESidual) :cite:p:`SaSc:86`, FGMRES (Flexible Generalized
Minimum RESidual) :cite:p:`Saa:93`, Bi-CGStab (Bi-Conjugate
Gradient Stabilized) :cite:p:`Van:92`, TFQMR (Transpose-Free
Quasi-Minimal Residual) :cite:p:`Fre:93`, and PCG
(Preconditioned Conjugate Gradient) :cite:p:`HeSt:52` linear
iterative methods. As Krylov methods, these require almost no matrix
storage for solving the Newton equations as compared to direct methods.
However, the algorithms allow for a user-supplied preconditioner matrix,
and for most problems preconditioning is essential for an efficient
solution. For very large stiff ODE systems, the Krylov methods are
preferable over direct linear solver methods, and are often the only
feasible choice. Among the Krylov methods in SUNDIALS, we recommend
GMRES as the best overall choice. However, users are encouraged to
compare all options, especially if encountering convergence failures
with GMRES. Bi-CGStab and TFQMR have an advantage in storage
requirements, in that the number of workspace vectors they require is
fixed, while that number for GMRES depends on the desired Krylov
subspace size. FGMRES has an advantage in that it is designed to support
preconditioners that vary between iterations (e.g. iterative methods).
PCG exhibits rapid convergence and minimal workspace vectors, but only
works for symmetric linear systems.

In the process of translating the VODE and VODPK algorithms into
C, the overall CVODE organization has been changed considerably.
One key feature of the CVODE organization is that the linear system
solvers comprise a layer of code modules that is separated from the
integration algorithm, allowing for easy modification and expansion of
the linear solver array. A second key feature is a separate module
devoted to vector operations; this facilitated the extension to
multiprosessor environments with minimal impacts on the rest of the
solver, resulting in PVODE :cite:p:`ByHi:99`, the parallel
variant of CVODE.

Around 2002, the functionality of CVODE and PVODE were combined
into one single code, simply called CVODE. Development of this
version of CVODE was concurrent with a redesign of the vector
operations module across the SUNDIALS suite. The key feature of the
``N_Vector`` module is that it is written in terms of abstract vector
operations with the actual vector kernels attached by a particular
implementation (such as serial or parallel) of ``N_Vector``. This allows
writing the SUNDIALS solvers in a manner independent of the actual
``N_Vector`` implementation (which can be user-supplied), as well as
allowing more than one ``N_Vector`` module linked into an executable
file. SUNDIALS (and thus CVODE) is supplied with a wide range of different
``N_Vector`` implementations, including: serial, MPI-parallel, both OpenMP and
Pthreads thread-parallel ``N_Vector`` implementations, a Hypre parallel
implementation, a PETSc implementation, and various GPU-enabled
implementations.

.. There are several motivations for choosing the C language for
.. CVODE. First, a general movement away from Fortran and toward C in
.. scientific computing was apparent. Second, the pointer, structure, and
.. dynamic memory allocation features in C are extremely useful in software
.. of this complexity, with the great variety of method options offered.
.. Finally, we prefer C over |CPP| for CVODE because of the
.. wider availability of C compilers, the potentially greater
.. efficiency of C, and the greater ease of interfacing the solver to
.. applications written in extended Fortran.

Changes to SUNDIALS in release X.Y.Z
====================================

.. include:: ../../../shared/RecentChanges.rst

For changes in prior versions of SUNDIALS see :numref:`Changelog`.

.. _CVODE.Introduction.reading:

Reading this User Guide
=======================

This user guide is a combination of general usage instructions. Specific
example programs are provided as a separate document. We expect that
some readers will want to concentrate on the general instructions, while
others will refer mostly to the examples, and the organization is
intended to accommodate both styles.

There are different possible levels of usage of CVODE. The most
casual user, with a small IVP problem only, can get by with reading
:numref:`CVODE.Mathematics.ivp_sol`, then :numref:`CVODE.Usage.CC` through
:numref:`CVODE.Usage.CC.cvode` only, and looking at examples
in :cite:p:`cvode_ex`.

In a different direction, a more expert user with an IVP problem may
want to (a) use a package preconditioner
(:numref:`CVODE.Usage.CC.precond`), (b) supply his/her own Jacobian
or preconditioner routines
(:numref:`CVODE.Usage.CC.user_fct_sim.jacFn`), (c) do multiple runs of
problems of the same size (:numref:`CVODE.Usage.CC.reinit`), (d)
supply a new ``N_Vector`` module (:numref:`NVectors`), (e)
supply new ``SUNLinearSolver`` and/or ``SUNMatrix`` modules
(:numref:`SUNMatrix` and :numref:`SUNLinSol`),
or even (f) supply new ``SUNNonlinearSolver`` modules
(:numref:`SUNNonlinSol`).

The structure of this document is as follows:

-  In :numref:`CVODE.Mathematics`, we give short descriptions of the
   numerical methods implemented by CVODE for the solution of
   initial value problems for systems of ODEs, and continue with short
   descriptions of preconditioning
   (:numref:`CVODE.Mathematics.preconditioning`), stability limit
   detection (:numref:`CVODE.Mathematics.stablimit`), and rootfinding
   (:numref:`CVODE.Mathematics.rootfinding`).

-  The following chapter describes the software organization of the CVODE
   solver (:numref:`CVODE.Organization`).

-  :numref:`CVODE.Usage.CC` is the main usage document
   for CVODE for C applications. It includes a complete
   description of the user interface for the integration of ODE initial
   value problems.

-  In :numref:`Fortran`, we describe the use of
   CVODE with Fortran applications.

-  :numref:`NVectors` gives a brief overview of the
   generic ``N_Vector`` module shared among the various components of
   SUNDIALS, and details on the ``N_Vector`` implementations provided
   with SUNDIALS.

-  :numref:`SUNMatrix` gives a brief overview of
   the generic ``SUNMatrix`` module shared among the various components
   of SUNDIALS, and details on the ``SUNMatrix`` implementations
   provided with SUNDIALS: a dense implementation
   (:numref:`SUNMatrix.Dense`), a banded
   implementation (:numref:`SUNMatrix.Band`) and a
   sparse implementation
   (:numref:`SUNMatrix.Sparse`).

-  :numref:`SUNLinSol` gives a brief overview of
   the generic ``SUNLinearSolver`` module shared among the various components
   of SUNDIALS. This chapter contains details on the ``SUNLinearSolver``
   implementations provided with SUNDIALS. The chapter also contains
   details on the ``SUNLinearSolver`` implementations provided with
   SUNDIALS that interface with external linear solver libraries.

-  :numref:`SUNNonlinSol` describes the
   ``SUNNonlinearSolver`` API and nonlinear solver implementations shared
   among the various components of SUNDIALS.

-  Finally, in the appendices, we provide detailed instructions for the
   installation of CVODE, within the structure of SUNDIALS
   (:numref:`Installation`), as well as a list of all the
   constants used for input to and output from CVODE functions
   (:numref:`CVODE.Constants`).

Finally, the reader should be aware of the following notational
conventions in this user guide: program listings and identifiers (such
as :c:func:`CVodeInit`) within textual explanations are hyperlinked
to their definitions directly; fields
in C structures (such as *content*) appear in italics; and packages
or modules, such as CVLS, are written in all capitals.


SUNDIALS License and Notices
============================

.. ifconfig:: package_name != 'super'

   .. include:: ../../../shared/LicenseReleaseNumbers.rst

.. ifconfig:: package_name == 'super'

   All SUNDIALS packages are released open source, under the BSD 3-Clause
   license for more details see the LICENSE and NOTICE files provided with all
   SUNDIALS packages.


.. _CVODE.Introduction.Ack:

Acknowledgments
===============

We wish to acknowledge the contributions to previous versions of the
CVODE and PVODE codes and their user guides by Scott D.
Cohen :cite:p:`CoHi:94` and George D.
Byrne :cite:p:`ByHi:98`.
