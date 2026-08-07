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

.. _IDA.Introduction:

************
Introduction
************

IDA is part of a software family called SUNDIALS: SUite of Nonlinear and
DIfferential/ALgebraic equation Solvers :cite:p:`HBGLSSW:05`.  This suite
consists of CVODE, ARKODE, KINSOL, and IDA, and variants of these with
sensitivity analysis capabilities, CVODES and IDAS.

IDA is a general purpose solver for the initial value problem (IVP) for systems
of differential-algebraic equations (DAEs). The name IDA stands for Implicit
Differential-Algebraic solver. IDA is based on DASPK :cite:p:`BHP:94,BHP:98`,
but is written in ANSI-standard C rather than Fortran77.  Its most notable
features are that,
(1) in the solution of the underlying nonlinear system at each time step, it
offers a choice of Newton/direct methods and a choice of Inexact Newton/Krylov
(iterative) methods; and
(2) it is written in a *data-independent* manner in that it acts on generic
vectors and matrices without any assumptions on the underlying organization of
the data.  Thus IDA shares significant modules previously written within CASC at
LLNL to support the ordinary differential equation (ODE) solvers CVODE
:cite:p:`cvode_ug,CoHi:96` and PVODE :cite:p:`ByHi:98,ByHi:99`, and also the
nonlinear system solver KINSOL :cite:p:`kinsol_ug`.

At present, IDA may utilize a variety of Krylov methods provided in SUNDIALS
that can be used in conjunction with Newton iteration: these include the GMRES
(Generalized Minimal RESidual) :cite:p:`SaSc:86`, FGMRES (Flexible Generalized
Minimum RESidual) :cite:p:`Saa:93`, Bi-CGStab (Bi-Conjugate Gradient Stabilized)
:cite:p:`Van:92`, TFQMR (Transpose-Free Quasi-Minimal Residual)
:cite:p:`Fre:93`, and PCG (Preconditioned Conjugate Gradient) :cite:p:`HeSt:52`
linear iterative methods. As Krylov methods, these require little matrix storage
for solving the Newton equations as compared to direct methods. However, the
algorithms allow for a user-supplied preconditioner, and, for most
problems, preconditioning is essential for an efficient solution.

For very large DAE systems, the Krylov methods are preferable over direct linear
solver methods, and are often the only feasible choice.  Among the Krylov
methods in SUNDIALS, we recommend GMRES as the best overall choice. However,
users are encouraged to compare all options, especially if encountering
convergence failures with GMRES.  Bi-CGFStab and TFQMR have an advantage in
storage requirements, in that the number of workspace vectors they require is
fixed, while that number for GMRES depends on the desired Krylov subspace
size. FGMRES has an advantage in that it is designed to support preconditioners
that vary between iterations (e.g. iterative methods). PCG exhibits rapid
convergence and minimal workspace vectors, but only works for symmetric linear
systems.

..
   There are several motivations for choosing the C language for IDA.  First, a
   general movement away from Fortran and toward C in scientific computing was
   apparent. Second, the pointer, structure, and dynamic memory allocation features
   in C are extremely useful in software of this complexity, with the great variety
   of method options offered.  Finally, we prefer C over C++ for IDA because of the
   wider availability of C compilers, the potentially greater efficiency of C, and
   the greater ease of interfacing the solver to applications written in extended
   Fortran.

Changes to SUNDIALS in release X.Y.Z
====================================

.. include:: ../../../shared/RecentChanges.rst

For changes in prior versions of SUNDIALS see :numref:`Changelog`.


.. _IDA.Introduction.Reading:

Reading this User Guide
=======================

The structure of this document is as follows:

* In Chapter :numref:`IDA.Mathematics`, we give short descriptions of the numerical
  methods implemented by IDA for the solution of initial value problems for
  systems of DAEs, along with short descriptions of preconditioning
  (:numref:`IDA.Mathematics.Preconditioning`) and rootfinding
  (:numref:`IDA.Mathematics.rootfinding`).

* The following chapter describes the software organization of the IDA
  solver (:numref:`IDA.Organization`).

* Chapter :numref:`IDA.Usage.CC` is the main usage document for IDA for C and C++
  applications. It includes a complete description of the user interface for the
  integration of DAE initial value problems. This is followed by documentation
  for using IDA with Fortran applications and on GPU accelerated systems.

* Chapter :numref:`NVectors` gives a brief overview of the generic ``N_Vector``
  module shared among the various components of SUNDIALS, as well as details on
  the ``N_Vector`` implementations provided with SUNDIALS.

* Chapter :numref:`SUNMatrix` gives a brief overview of the generic
  ``SUNMatrix`` module shared among the various components of SUNDIALS, and
  details on the ``SUNMatrix`` implementations provided with SUNDIALS.

* Chapter :numref:`SUNLinSol` gives a brief overview of the generic
  ``SUNLinearSolver`` module shared among the various components of
  SUNDIALS. This chapter contains details on the ``SUNLinearSolver``
  implementations provided with SUNDIALS.  The chapter also contains details on
  the ``SUNLinearSolver`` implementations provided with SUNDIALS that interface
  with external linear solver libraries.

* Chapter :numref:`SUNNonlinSol` describes the ``SUNNonlinearSolver`` API and
  nonlinear solver implementations shared among the various components of
  SUNDIALS.

* Finally, in the appendices, we provide detailed instructions for the
  installation of IDA, within the structure of SUNDIALS (Appendix
  :numref:`Installation`), as well as a list of all the constants used for input
  to and output from IDA functions (Appendix :numref:`IDA.Constants`).

..
   Finally, the reader should be aware of the following notational conventions in
   this user guide: program listings and identifiers (such as :c:func:`IDAInit`)
   within textual explanations appear in typewriter type style; fields in C
   structures (such as *content*) appear in italics; and packages or modules, such
   as IDADLS, are written in all capitals.


SUNDIALS License and Notices
============================

.. ifconfig:: package_name != 'super'

   .. include:: ../../../shared/LicenseReleaseNumbers.rst

.. ifconfig:: package_name == 'super'

   All SUNDIALS packages are released open source, under the BSD 3-Clause
   license for more details see the LICENSE and NOTICE files provided with all
   SUNDIALS packages.
