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

.. _CVODES.Introduction:

************
Introduction
************

CVODES :cite:p:`SeHi:05` is part of a software family called SUNDIALS: SUite of
Nonlinear and DIfferential/ALgebraic equation Solvers :cite:p:`HBGLSSW:05`. This
suite consists of CVODE, ARKODE, KINSOL, and IDA, and variants of these with
sensitivity analysis capabilities. CVODES is a solver for stiff and nonstiff
initial value problems (IVPs) for systems of ordinary differential equation
(ODEs). In addition to solving stiff and nonstiff ODE systems, CVODES has
sensitivity analysis capabilities, using either the forward or the adjoint
methods.

.. _CVODES.Introduction.History:

Historical Background
=====================

Fortran solvers for ODE initial value problems are widespread and heavily used.
Two solvers that have been written at LLNL in the past are VODE :cite:p:`BBH:89`
and VODPK :cite:p:`Byr:92`. VODE is a general purpose solver that includes
methods for both stiff and nonstiff systems, and in the stiff case uses direct
methods (full or banded) for the solution of the linear systems that arise at
each implicit step. Externally, VODE is very similar to the well known solver
LSODE :cite:p:`RaHi:94`. VODPK is a variant of VODE that uses a preconditioned
Krylov (iterative) method, namely GMRES, for the solution of the linear systems.
VODPK is a powerful tool for large stiff systems because it combines established
methods for stiff integration, nonlinear iteration, and Krylov (linear)
iteration with a problem-specific treatment of the dominant source of stiffness,
in the form of the user-supplied preconditioner matrix :cite:p:`BrHi:89`. The
capabilities of both VODE and VODPK have been combined in the C-language package
CVODE :cite:p:`CoHi:96`.

At present, CVODE may utilize a variety of Krylov methods provided in SUNDIALS
that can be used in conjunction with Newton iteration: these include the GMRES
(Generalized Minimal RESidual) :cite:p:`SaSc:86`, FGMRES (Flexible Generalized
Minimum RESidual) :cite:p:`Saa:93`, Bi-CGStab (Bi-Conjugate Gradient Stabilized)
:cite:p:`Van:92`, TFQMR (Transpose-Free Quasi-Minimal Residual)
:cite:p:`Fre:93`, and PCG (Preconditioned Conjugate Gradient) :cite:p:`HeSt:52`
linear iterative methods. As Krylov methods, these require almost no matrix
storage for solving the Newton equations as compared to direct methods. However,
the algorithms allow for a user-supplied preconditioner matrix, and for most
problems preconditioning is essential for an efficient solution. For very large
stiff ODE systems, the Krylov methods are preferable over direct linear solver
methods, and are often the only feasible choice. Among the Krylov methods in
SUNDIALS, we recommend GMRES as the best overall choice. However, users are
encouraged to compare all options, especially if encountering convergence
failures with GMRES. Bi-CGStab and TFQMR have an advantage in storage
requirements, in that the number of workspace vectors they require is fixed,
while that number for GMRES depends on the desired Krylov subspace size. FGMRES
has an advantage in that it is designed to support preconditioners that vary
between iterations (e.g., iterative methods). PCG exhibits rapid convergence and
minimal workspace vectors, but only works for symmetric linear systems.

In the process of translating the VODE and VODPK algorithms into C, the overall
CVODE organization has been changed considerably. One key feature of the CVODE
organization is that the linear system solvers comprise a layer of code modules
that is separated from the integration algorithm, allowing for easy modification
and expansion of the linear solver array. A second key feature is a separate
module devoted to vector operations; this facilitated the extension to
multiprosessor environments with minimal impacts on the rest of the solver,
resulting in PVODE :cite:p:`ByHi:99`, the parallel variant of CVODE.

CVODES is written with a functionality that is a superset of that of the pair
CVODE/PVODE. Sensitivity analysis capabilities, both forward and adjoint, have
been added to the main integrator. Enabling forward sensitivity computations
in CVODES will result in the code integrating the so-called *sensitivity
equations* simultaneously with the original IVP, yielding both the solution and
its sensitivity with respect to parameters in the model. Adjoint sensitivity
analysis, most useful when the gradients of relatively few functionals of the
solution with respect to many parameters are sought, involves integration of the
original IVP forward in time followed by the integration of the so-called
*adjoint equations* backward in time. CVODES provides the infrastructure needed
to integrate any final-condition ODE dependent on the solution of the original
IVP (in particular the adjoint system).

Development of CVODES was concurrent with a redesign of the vector operations
module across the SUNDIALS suite. The key feature of the ``N_Vector`` module is
that it is written in terms of abstract vector operations with the actual vector
functions attached by a particular implementation (such as serial or parallel)
of ``N_Vector``. This allows writing the SUNDIALS solvers in a manner
independent of the actual ``N_Vector`` implementation (which can be
user-supplied), as well as allowing more than one ``N_Vector`` module to be
linked into an executable file. SUNDIALS (and thus CVODES) is supplied with
serial, MPI-parallel, and both OpenMP and Pthreads thread-parallel ``N_Vector``
implementations.

There were several motivations for choosing the C language for CVODE, and later
for CVODES. First, a general movement away from Fortran and toward C in
scientific computing was apparent. Second, the pointer, structure, and dynamic
memory allocation features in C are extremely useful in software of this
complexity. Finally, we prefer C over C++ for CVODES because of the wider
availability of C compilers, the potentially greater efficiency of C, and the
greater ease of interfacing the solver to applications written in extended
Fortran.


Changes to SUNDIALS in release 7.6.0
====================================

.. include:: ../../../shared/RecentChanges.rst

For changes in prior versions of SUNDIALS see :numref:`Changelog`.


.. _CVODES.Introduction.Reading:

Reading this User Guide
=======================

This user guide is a combination of general usage instructions. Specific example
programs are provided as a separate document. We expect that some readers will
want to concentrate on the general instructions, while others will refer mostly
to the examples, and the organization is intended to accommodate both styles.

There are different possible levels of usage of CVODES. The most casual user,
with a small IVP problem only, can get by with reading
:numref:`CVODES.Mathematics.ivp_sol`, then Chapter :numref:`CVODES.Usage.SIM` up
to :numref:`CVODES.Usage.Purequad` only, and looking at examples in
:cite:p:`cvodes_ex`. In addition, to solve a forward sensitivity problem the
user should read :numref:`CVODES.Mathematics.FSA`, followed by Chapter
:numref:`CVODES.Usage.FSA` and look at examples in :cite:p:`cvodes_ex`.

In a different direction, a more expert user with an IVP problem may want to (a)
use a package preconditioner (:numref:`CVODES.Usage.SIM.precond`), (b) supply
his/her own Jacobian or preconditioner routines
(:numref:`CVODES.Usage.SIM.user_supplied`), (c) do multiple runs of problems of
the same size (:c:func:`CVodeReInit`), (d) supply a new ``N_Vector`` module
(:numref:`NVectors`), or even (e) supply new ``SUNLinearSolver`` and/or
``SUNMatrix`` modules (Chapters :numref:`SUNMatrix` and :numref:`SUNLinSol`). An
advanced user with a forward sensitivity problem may also want to (a) provide
his/her own sensitivity equations right-hand side routine
:numref:`CVODES.Usage.FSA.user_supplied`, (b) perform multiple runs with the
same number of sensitivity parameters
(:numref:`CVODES.Usage.FSA.user_callable.sensi_malloc`, or (c) extract additional
diagnostic information
(:numref:`CVODES.Usage.FSA.user_callable.optional_output`). A user with an
adjoint sensitivity problem needs to understand the IVP solution approach at the
desired level and also go through :numref:`CVODES.Mathematics.ASA` for a short
mathematical description of the adjoint approach, Chapter
:numref:`CVODES.Usage.ADJ` for the usage of the adjoint module in CVODES, and
the examples in :cite:p:`cvodes_ex`.

The structure of this document is as follows:

-  In Chapter :numref:`CVODES.Mathematics`, we give short descriptions of the
   numerical methods implemented by CVODES for the solution of initial value
   problems for systems of ODEs, continue with short descriptions of
   preconditioning :numref:`CVODES.Mathematics.preconditioning`, stability limit
   detection (:numref:`CVODES.Mathematics.stablimit`), and rootfinding
   (:numref:`CVODES.Mathematics.rootfinding`), and conclude with an overview of
   the mathematical aspects of sensitivity analysis, both forward
   (:numref:`CVODES.Mathematics.FSA`) and adjoint
   (:numref:`CVODES.Mathematics.ASA`).

-  The following chapter describes the software organization of the CVODES
   solver (:numref:`CVODES.Organization`).

-  Chapter :numref:`CVODES.Usage.SIM` is the main usage document for CVODES for simulation applications.
   It includes a complete description of the user interface for the integration
   of ODE initial value problems. Readers that are not interested in using
   CVODES for sensitivity analysis can then skip the next two chapters.

-  Chapter :numref:`CVODES.Usage.FSA` describes the usage of CVODES for forward sensitivity analysis as an
   extension of its IVP integration capabilities. We begin with a skeleton of
   the user main program, with emphasis on the steps that are required in
   addition to those already described in Chapter :numref:`CVODES.Usage.SIM`.
   Following that we provide detailed descriptions of the
   user-callable interface routines specific to forward sensitivity analysis and
   of the additional optional user-defined routines.

-  Chapter :numref:`CVODES.Usage.ADJ` describes the usage of CVODES for adjoint sensitivity analysis. We begin
   by describing the CVODES checkpointing implementation for interpolation of
   the original IVP solution during integration of the adjoint system backward
   in time, and with an overview of a user's main program. Following that we
   provide complete descriptions of the user-callable interface routines for
   adjoint sensitivity analysis as well as descriptions of the required
   additional user-defined routines.

-  Chapter :numref:`NVectors` gives a brief overview of the generic ``N_Vector`` module shared among the
   various components of SUNDIALS, and details on the ``N_Vector`` implementations provided with SUNDIALS.

-  Chapter :numref:`SUNMatrix` gives a brief overview of the generic
   ``SUNMatrix`` module shared among the various components of SUNDIALS, and
   details on the ``SUNMatrix`` implementations provided with SUNDIALS: a dense
   implementation (:numref:`SUNMatrix.Dense`), a banded implementation
   (:numref:`SUNMatrix.Band`) and a sparse implementation
   (:numref:`SUNMatrix.Sparse`).

-  Chapter :numref:`SUNLinSol` gives a brief overview of the generic ``SUNLinearSolver`` module shared among
   the various components of SUNDIALS. This chapter contains details on the
   ``SUNLinearSolver`` implementations provided with SUNDIALS. The chapter also
   contains details on the ``SUNLinearSolver`` implementations provided with
   SUNDIALS that interface with external linear solver libraries.

-  Finally, in the appendices, we provide detailed instructions for the installation of CVODES, within the
   structure of SUNDIALS (Appendix :numref:`Installation`), as well as a
   list of all the constants used for input to and output from CVODES functions
   (Appendix :numref:`CVODES.Constants`).

Finally, the reader should be aware of the following notational conventions in
this user guide: program listings and identifiers (such as ``CVodeInit``) within
textual explanations appear in typewriter type style; fields in C structures
(such as *content*) appear in italics; and packages or modules, such as CVDLS,
are written in all capitals.

.. warning::

   Usage and installation instructions that constitute important warnings are
   marked in yellow boxes like this one.


SUNDIALS License and Notices
============================

.. ifconfig:: package_name != 'super'

   .. include:: ../../../shared/LicenseReleaseNumbers.rst

.. ifconfig:: package_name == 'super'

   All SUNDIALS packages are released open source, under the BSD 3-Clause
   license for more details see the LICENSE and NOTICE files provided with all
   SUNDIALS packages.
