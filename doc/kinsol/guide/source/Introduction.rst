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

.. _KINSOL.Introduction:

************
Introduction
************

KINSOL is part of a software family called SUNDIALS: SUite of Nonlinear and DIfferential/ALgebraic equation
Solvers :cite:p:`HBGLSSW:05`. This suite consists of CVODE, ARKODE, KINSOL, and IDA,
and variants of these with sensitivity analysis capabilities.

KINSOL is a general-purpose nonlinear system solver based on Newton-Krylov solver technology. A fixed point
iteration is also included with the release of KINSOL v.2.8.0 and higher.

.. _KINSOL.Introduction.Historical:

Historical Background
=====================

The first nonlinear solver packages based on Newton-Krylov methods were written in Fortran. In particular, the NKSOL
package, written at LLNL, was the first Newton-Krylov solver package written for solution of systems arising in the
solution of partial differential equations :cite:p:`BrSa:90`. This Fortran code made use of Newton's method to
solve the discrete nonlinear systems and applied a preconditioned Krylov linear solver for solution of the Jacobian
system at each nonlinear iteration. The key to the Newton-Krylov method was that the matrix-vector multiplies required
by the Krylov method could effectively be approximated by a finite difference of the nonlinear system-defining function,
avoiding a requirement for the formation of the actual Jacobian matrix. Significantly less memory was required for the
solver as a result.

In the late 1990s, there was a push at LLNL to rewrite the nonlinear solver in C and port it to distributed
memory parallel machines. Both Newton and Krylov methods are easily implemented in parallel, and this effort gave rise
to the KINSOL package. KINSOL is similar to NKSOL in functionality, except that it provides for more options
in the choice of linear system methods and tolerances, and has a more modular design to provide flexibility for future
enhancements.

At present, KINSOL may utilize a variety of Krylov methods provided in SUNDIALS. These methods include the
GMRES (Generalized Minimal RESidual) :cite:p:`SaSc:86`, FGMRES (Flexible Generalized Minimum
RESidual) :cite:p:`Saa:93`, Bi-CGStab (Bi-Conjugate Gradient Stabilized) :cite:p:`Van:92`, TFQMR
(Transpose-Free Quasi-Minimal Residual) :cite:p:`Fre:93`, and PCG (Preconditioned Conjugate
Gradient) :cite:p:`HeSt:52` linear iterative methods. As Krylov methods, these require little matrix storage
for solving the Newton equations as compared to direct methods. However, the algorithms allow for a user-supplied
preconditioner, and, for most problems, preconditioning is essential for an efficient solution. For very large
nonlinear algebraic systems, the Krylov methods are preferable over direct linear solver methods, and are often the only
feasible choice. Among the Krylov methods in SUNDIALS, we recommend GMRES as the best overall choice. However,
users are encouraged to compare all options, especially if encountering convergence failures with GMRES. Bi-CGStab and
TFQMR have an advantage in storage requirements, in that the number of workspace vectors they require is fixed, while
that number for GMRES depends on the desired Krylov subspace size. FGMRES has an advantage in that it is designed to
support preconditioners that vary between iterations (e.g., iterative methods). PCG exhibits rapid convergence and
minimal workspace vectors, but only works for symmetric linear systems.

For the sake of completeness in functionality, direct linear system solvers are included in KINSOL. These include
methods for both dense and banded linear systems, with Jacobians that are either user-supplied or generated internally
by difference quotients. KINSOL also includes interfaces to sparse direct solvers, including
KLU :cite:p:`DaPa:10,KLU_site` and the threaded sparse direct solver,
SuperLU_MT :cite:p:`Li:05,DGL:99,SuperLUMT_site`, among others (see Chapter :numref:`SUNLinSol` for further details).

In the process of translating NKSOL into C, the overall KINSOL organization has been changed considerably.
One key feature of the KINSOL organization is that a separate module devoted to vector operations was created.
This module facilitated extension to multiprosessor environments with minimal impact on the rest of the solver. The
vector module design is shared across the SUNDIALS suite. This :c:type:`N_Vector` module is written in terms of
abstract vector operations with the actual routines attached by a particular implementation (such as serial or parallel)
of ``N_Vector``. This abstraction allows writing the SUNDIALS solvers in a manner independent of the actual
``N_Vector`` implementation (which can be user-supplied), as well as allowing more than one ``N_Vector`` module linked
into an executable file. SUNDIALS (and thus KINSOL) is supplied with serial, MPI-parallel, OpenMP
and Pthreads thread-parallel ``N_Vector`` implementations, as well as multiple ``N_Vector``
implementations designed to leverage GPU architectures (see Chapter :numref:`NVectors` for
further details).

There are several motivations for choosing the C language for KINSOL. First, a general movement away from
Fortran and toward C in scientific computing was apparent. Second, the pointer, structure, and dynamic memory
allocation features in C are extremely useful in software of this complexity, with the great variety of method options
offered. Finally, we prefer C over C++ for KINSOL because of the wider availability of C
compilers, the potentially greater efficiency of C, and the greater ease of interfacing the solver to
applications written in Fortran.

.. _KINSOL.Introduction.Changes:

Changes to SUNDIALS in release X.Y.Z
====================================

.. include:: ../../../shared/RecentChanges.rst

For changes in prior versions of SUNDIALS see :numref:`Changelog`.

.. _KINSOL.Introduction.reading:

Reading this User Guide
=======================

This user guide is a combination of general usage instructions and specific examples. We expect that some readers will
want to concentrate on the general instructions, while others will refer mostly to the examples, and the organization is
intended to accommodate both styles.

There are different possible levels of usage of KINSOL. The most casual user, with a small nonlinear system, can
get by with reading all of Chapter :numref:`KINSOL.Mathematics`, then Chapter :numref:KINSOL.Usage.CC through
:numref:`KINSOL.Usage.CC` only, and looking at examples in :cite:p:`kinsol_ex`. In a different
direction, a more expert user with a nonlinear system may want to (a) use a package preconditioner
(:numref:`KINSOL.Usage.CC.kin_bbdpre`), (b) supply his/her own Jacobian or preconditioner routines
(:numref:`KINSOL.Usage.CC.user_fct_sim`), (c) supply a new ``N_Vector`` module
(Chapter :numref:`NVectors`), or even (d) supply a different linear solver module
(:numref:`KINSOL.Usage.CC.callable_fct_sim.lin_solv_init` and Chapter :numref:`SUNLinSol`).

The structure of this document is as follows:

-  In Chapter :numref:`KINSOL.Mathematics`, we provide short descriptions of the numerical methods implemented by KINSOL
   for the solution of nonlinear systems.

-  The following chapter describes the software organization of the KINSOL
   solver (:numref:`KINSOL.Organization`).

-  Chapter :numref:KINSOL.Usage.CC is the main usage document for KINSOL for C applications. It includes a
   complete description of the user interface for the solution of nonlinear algebraic systems.

-  Chapter :numref:`NVectors` gives a brief overview of the generic ``N_Vector`` module shared among the
   various components of SUNDIALS, and details on the four ``N_Vector`` implementations provided with
   SUNDIALS.

-  Chapter :numref:`SUNMatrix` gives a brief overview of the generic ``SUNMatrix`` module shared among
   the various components of SUNDIALS, and details on the ``SUNMatrix`` implementations provided with
   SUNDIALS.

-  Chapter :numref:`SUNLinSol` gives a brief overview of the generic ``SUNLinearSolver`` module shared among
   the various components of SUNDIALS. This chapter contains details on the ``SUNLinearSolver`` implementations
   provided with SUNDIALS. The chapter also contains details on the ``SUNLinearSolver`` implementations provided with
   SUNDIALS that interface with external linear solver libraries.

-  Finally, in the appendices, we provide detailed instructions for the installation of KINSOL, within the
   structure of SUNDIALS (Appendix :numref:`Installation`), as well as a list of all the constants used for
   input to and output from KINSOL functions (Appendix :numref:`KINSOL.Constants`).

Finally, the reader should be aware of the following notational conventions in this user guide: program listings and
identifiers (such as ``KINInit``) within textual explanations appear in typewriter type style; fields in C
structures (such as *content*) appear in italics; and packages or modules are written in all capitals. Usage and


SUNDIALS License and Notices
============================

.. ifconfig:: package_name != 'super'

   .. include:: ../../../shared/LicenseReleaseNumbers.rst

.. ifconfig:: package_name == 'super'

   All SUNDIALS packages are released open source, under the BSD 3-Clause
   license for more details see the LICENSE and NOTICE files provided with all
   SUNDIALS packages.


Acknowledgments
===============

We wish to acknowledge the contributions to previous versions of the KINSOL code and user guide by Allan G.
Taylor.
