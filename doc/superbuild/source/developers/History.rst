..
   Author(s): David J. Gardner and Carol S. Woodward @ LLNL
   -----------------------------------------------------------------------------
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
   -----------------------------------------------------------------------------

.. _History:

History
=======

SUNDIALS originated in 1993 as a rewrite in C of two general purpose Fortran
solvers for Ordinary Differential Equations, VODE and VODPK. These packages had
evolved from some of the first and most heavily used ODE integrator packages ever
developed, first available through the NETLIB FTP site as part of ODEPACK (which
had more than 250,000 downloads). Based on linear multistep methods for both
non-stiff and stiff systems, these solvers used adaptive time-stepping and
variable-order schemes to deliver solutions whose error within each time step
was bounded by user-supplied tolerances. The VODE line of codes represented a
breakthrough in general purpose ODE integrators in that they used
variable-coefficient multistep methods instead of fixed-step-interpolated
methods. The variable coefficient form allowed for efficient solvers tuned
through integrator heuristics. The goal of the 1993 rewrite was a version of
these heavily used solvers that could be easily extended from serial to
distributed memory parallel computers. The resulting serial code, called CVODE,
was followed by a parallel version, called PVODE. The ease of making this
conversion was achieved by the use of two sets of routines for basic vector
operations, one serial and one parallel (based on MPI).

The success of the rewrite prompted a similar rewrite of a nonlinear algebraic
systems solver, NKSOL, which was the first Newton-Krylov-based nonlinear
software package, into the KINSOL package. (CNKSOL was voted out as a name due
to a desire of LLNL management for a pronounceable package name, hence KINSOL,
for Krylov Inexact Newton SOLver.) Further development of the common
infrastructure led to a C rewrite of a pair of differential-algebraic system
solvers, DASSL and DASPK, resulting in the IDA package. The original name for
this group of three codes was "The PVODE Trio." With careful encapsulation of
data vectors, the codes were engineered so that the core integrator and
nonlinear solvers were the same whether run in serial or parallel. With this
advance, the "P" was replaced with "C" in the CVODE name making the PVODE Trio
name obsolete. The group of packages were then given the acronym SUNDIALS, for
SUite of Nonlinear DIfferential and ALgebraic equation Solvers, and first
released under that name in 2002.

Extensions to the time integrators to handle forward and adjoint sensitivity
problems, called CVODES and IDAS, were added later. In 2015 a new ODE
integration package, ARKODE, was added to SUNDIALS. ARKODE includes multistage
methods which are better suited to time evolution systems posed within spatially
adaptive codes than the multistep methods in CVODE. The six packages in SUNDIALS
all share a collection of vector operation modules, which originally included
serial and MPI-parallel versions. In 2015, versions supporting threaded
environments were added, and in 2017 GPU support through CUDA and RAJA was
added. Recent extensions to SUNDIALS have included wrappers for PETSc and hypre
vectors. Richer interoperability layers with those and other DOE solver packages
are planned through the Exascale Computing Project.

.. image:: ../figs/sundials-timeline-2020.png

A timeline of SUNDIALS development:

 * 1982 - First release of ODEPACK
 * 1986 -- 1991 - Development of VODE and VODPK
 * 1993 - Initial rewrite of VODE/VODPK into CVODE
 * 1994 - CVODE first released
 * 1995 - Initial version of PVODE written
 * 1998 - PVODE first released
 * 1998 - KINSOL first released
 * 1999 - IDA first released
 * 2002 - First release of SUNDIALS under a BSD license
 * 2004 - CVODES first released
 * 2009 - IDAS first released
 * 2015 - ARKODE first released as part of SUNDIALS
 * 2016 - Changed procedures to allow for unregistered downloads
 * 2017 - Created GitHub mirror of released software for download
