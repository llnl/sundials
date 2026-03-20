..
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

.. ARKODE_example documentation master file, created by
   sphinx-quickstart on Sat Dec 22 20:38:03 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Introduction
============

The following tables summarize the salient features of each of the
example problems in this document.  Each example is designed to be
relatively self-contained, so that you need only study and/or emulate
the problem that is most closely related to your own.  We group these
examples according to programming language (C, C++, Fortran).

ARKODE example problems written in C are summarized in the table
below, and are further described in the chapters :ref:`serial_c`,
:ref:`openmp_c`, :ref:`parallel_c` and :ref:`parhyp_c`.

.. cssclass:: table-bordered

================================  ==========  ===========  ==========  =============  =====================================================
Problem                           Integrator  Nonlinear    Linear      Size           Extras
================================  ==========  ===========  ==========  =============  =====================================================
:ref:`ark_analytic`               DIRK        Newton       Dense       1
:ref:`ark_analytic_nonlin`        ERK         N.A.         N.A.        1              ERKStep timestepping module
:ref:`ark_brusselator`            DIRK        Newton       Dense       3
:ref:`ark_brusselator_fp`         ARK         Fixed-point  N.A.        3
:ref:`ark_robertson`              DIRK        Newton       Dense       3
:ref:`ark_robertson_root`         DIRK        Newton       Dense       3              rootfinding
:ref:`ark_brusselator1D`          DIRK        Newton       Band        3N
:ref:`ark_brusselator1D_omp`      DIRK        Newton       Band        3N             OpenMP-enabled
:ref:`ark_brusselator1D_klu`      DIRK        Newton       KLU         3N             sparse matrices
:ref:`ark_brusselator1D_FEM_slu`  DIRK        Newton       SuperLU_MT  3N             finite-element, :math:`M\ne I`, sparse matrices
:ref:`ark_heat1D`                 DIRK        Newton       PCG         N
:ref:`ark_heat1D_adapt`           DIRK        Newton       PCG         (dynamic)      adaptive vector resizing
:ref:`ark_KrylovDemo_prec`        DIRK        Newton       SPGMR       216            multiple preconditioners
:ref:`ark_diurnal_kry_bbd_p`      DIRK        Newton       SPGMR       200            parallel, BBD preconditioner
:ref:`ark_diurnal_kry_p`          DIRK        Newton       SPGMR       200            parallel, block-diagonal precond.
:ref:`ark_diurnal_kry_ph`         DIRK        Newton       SPGMR       200            HYPRE parallel vector
================================  ==========  ===========  ==========  =============  =====================================================


ARKODE example problems written in C++ are summarized in the table
below, and are further described in the chapters :ref:`serial_cpp` and
:ref:`parallel_cpp`.

.. cssclass:: table-bordered

=======================  ==========  ===========  ======  =============  =================================
Problem                  Integrator  Nonlinear    Linear  Size           Extras
=======================  ==========  ===========  ======  =============  =================================
:ref:`ark_analytic_sys`  DIRK        Newton       Dense   3
:ref:`ark_heat2D`        DIRK        Newton       PCG     :math:`nx*ny`  parallel
=======================  ==========  ===========  ======  =============  =================================


..
   ARKODE example problems written in Fortran 77 are summarized in the table
   below, and are further described in the chapters :ref:`serial_f77` and
   :ref:`parallel_f77`.

   .. cssclass:: table-bordered

   ==========================   ==========  ===========  ======  =============  =================================
   Problem                      Integrator  Nonlinear    Linear  Size           Extras
   ==========================   ==========  ===========  ======  =============  =================================
   :ref:`fark_diurnal_kry_bp`   DIRK        Newton       SPGMR   10             banded preconditioner
   :ref:`fark_roberts_dnsL`     DIRK        Newton       Dense   3              LAPACK dense solver, rootfinding
   :ref:`fark_diag_kry_bbd_p`   DIRK        Newton       SPGMR   10*NProcs      parallel BBD preconditioner
   :ref:`fark_diag_non_p`       ERK         N.A.         N.A.    10*NProcs      parallel
   ==========================   ==========  ===========  ======  =============  =================================


   ARKODE example problems written in Fortran 90 are summarized in the table
   below, and are further described in the chapters :ref:`serial_f90` and
   :ref:`parallel_f90`.

   .. cssclass:: table-bordered

   ==========================  ==========  =========  ======  =============  ===============================================
   Problem                     Integrator  Nonlinear  Linear  Size           Extras
   ==========================  ==========  =========  ======  =============  ===============================================
   :ref:`ark_bruss`            ARK         Newton     Dense   3
   :ref:`ark_bruss1D_FEM_klu`  DIRK        Newton     KLU     3N             finite-element, :math:`M\ne I`, sparse matrices
   :ref:`fark_heat2D`          DIRK        Newton     PCG     :math:`nx*ny`  parallel
   ==========================  ==========  =========  ======  =============  ===============================================
